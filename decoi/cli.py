"""Portable command-line interface for study generation and inspection."""
import argparse
from datetime import datetime, timezone
import fcntl
from functools import partial
import http.server
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import uuid
import webbrowser

import yaml
from . import __version__


def resources():
    source = Path(__file__).resolve().parents[1]
    if (source/'main.nf').is_file():
        return source
    installed = Path(sys.prefix)/'share/decoi'
    if not (installed/'main.nf').is_file():
        raise ValueError('DECOI runtime resources are missing; reinstall the package in this environment.')
    return installed


def positive(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError('must be at least 1')
    return value


def memory(value):
    if not re.fullmatch(r'\d+(?:\.\d+)?\s*(?:MB|GB|TB)', value, re.I) or float(re.match(r'[\d.]+', value)[0]) <= 0:
        raise argparse.ArgumentTypeError('use a positive size such as "4 GB" or "32 GB"')
    return value.upper()


def parser():
    p = argparse.ArgumentParser(prog='decoi', description='Simulate microbial studies with known ground truth.',
                                epilog='Start small: decoi check, then decoi test -o test-output')
    p.add_argument('--version', action='version', version=f'DECOI {__version__}')
    sub = p.add_subparsers(dest='command', required=True)
    sub.add_parser('check', help='Check the supporting environment and bundled resources')
    for name, help_text in [('run', 'Generate a configured study through Nextflow'),
                            ('test', 'Create and run the tiny paired-assay demo')]:
        q = sub.add_parser(name, help=help_text)
        q.add_argument('-o', '--output', required=True, type=Path, help='Output directory (use --resume to reuse a run)')
        q.add_argument('--threads', type=positive, default=None, help='Simulation CPUs (run: 8, test: 1; capped by local CPU budget)')
        q.add_argument('--memory', type=memory, default=None, help='Memory requested per task (run: 32 GB, test: 4 GB)')
        q.add_argument('--max_cpus', type=positive, help='Total local CPU budget; defaults to available CPUs')
        q.add_argument('--max_memory', type=memory, help='Total local memory budget')
        q.add_argument('--max_tasks', type=positive, default=1, help='Maximum submitted/running tasks (default: 1)')
        q.add_argument('--executor', choices=['local', 'slurm'], default='local')
        q.add_argument('--account', help='Slurm account')
        q.add_argument('--partition', help='Slurm partition; otherwise use the cluster default')
        q.add_argument('--time_limit', default='24h', help='Slurm time limit (default: 24h)')
        q.add_argument('--work_dir', type=Path, help='Work directory; default: OUTPUT/.decoi/work')
        q.add_argument('--resume', action='store_true', help='Reuse matching completed Nextflow tasks')
        if name == 'run':
            q.add_argument('-c', '--config', type=Path, required=True, help='Simulation YAML')
            q.add_argument('--study', type=Path, help='Override study-design YAML')
            r = q.add_mutually_exclusive_group()
            r.add_argument('--genome_dir', type=Path, help='Raw genome FASTAs plus taxonomy.tsv')
            r.add_argument('--silva_fasta', type=Path, help='Local SILVA reference FASTA')
            q.add_argument('--wgs_reference', type=Path, help='Prelinked genomes.tsv bundle for SILVA-based runs')
            q.add_argument('--run_dada2', action='store_true', help='Run optional amplicon validation after simulation')
    q = sub.add_parser('prepare-reference', help='Prepare references without simulating a study')
    q.add_argument('-c', '--config', required=True, type=Path)
    q.add_argument('-o', '--output', required=True, type=Path)
    q = sub.add_parser('add-wgs', help='Add linked WGS to a compatible existing study')
    q.add_argument('-c', '--config', required=True, type=Path)
    q.add_argument('-o', '--output', required=True, type=Path, help='Existing direct-run dataset directory')
    q.add_argument('--threads', type=positive, default=8)
    q = sub.add_parser('report', help='Find or serve a completed study report')
    q.add_argument('-o', '--output', required=True, type=Path)
    q.add_argument('--serve', action='store_true')
    q.add_argument('--port', type=positive, default=8765)
    q.add_argument('--no-browser', action='store_true')
    return p


def read_config(path, root):
    path = path.expanduser().resolve()
    cfg = yaml.safe_load(path.read_text())
    if not isinstance(cfg, dict) or not isinstance(cfg.get('reference'), dict):
        raise ValueError('Configuration must contain a reference mapping.')
    # Bundled examples historically use paths relative to the resource root.
    bundled = path.parent.name == 'config' and (path.parent.parent/'main.nf').is_file()
    base = path.parent.parent if bundled else path.parent
    def resolve(value):
        p = Path(value).expanduser()
        return str((p if p.is_absolute() else base/p).resolve())
    if cfg.get('study_design_file'):
        cfg['study_design_file'] = resolve(cfg['study_design_file'])
    for section, keys in [('reference', ['genome_dir', 'output_dir', 'dada2_silva_fasta']), ('wgs', ['reference_dir'])]:
        for key in keys:
            if cfg.get(section, {}).get(key):
                cfg[section][key] = resolve(cfg[section][key])
    fixture = cfg.get('reference_fixtures', {}).get('mitochondria', {})
    if fixture.get('source_fasta'):
        fixture['source_fasta'] = resolve(fixture['source_fasta'])
    return cfg


def save_stable(path, text):
    if not path.exists() or path.read_text() != text:
        temp = path.with_suffix(path.suffix+'.tmp')
        temp.write_text(text)
        temp.replace(path)


def tool_checks(include_nextflow=True, dada2=False):
    names = ['Rscript', 'iss', 'cutadapt'] + (['nextflow', 'java'] if include_nextflow else [])
    failures = []
    for module in ['numpy', 'pandas', 'scipy', 'yaml', 'Bio']:
        found = importlib.util.find_spec(module) is not None
        print(f'  Python {module}: {"OK" if found else "MISSING"}', flush=True)
        if not found:
            failures.append(module)
    for name in names:
        location = shutil.which(name)
        print(f'  {name}: {location or "MISSING"}', flush=True)
        if not location:
            failures.append(name)
    if shutil.which('Rscript'):
        packages = ['SparseDOSSA2'] + (['dada2'] if dada2 else [])
        for package in packages:
            result = subprocess.run(['Rscript', '-e', f'quit(status=if(requireNamespace("{package}",quietly=TRUE)) 0 else 1)'], capture_output=True)
            print(f'  {package}: {"OK" if result.returncode == 0 else "MISSING"}', flush=True)
            if result.returncode:
                failures.append(package)
    if failures:
        raise ValueError('Missing dependencies: ' + ', '.join(failures) + '. Activate the DECOI mamba environment; see the installation guide.')


def stream(command, cwd, log):
    with log.open('w') as out:
        out.write('Command: '+shlex.join(map(str, command))+'\n'); out.flush()
        process = subprocess.Popen(list(map(str, command)), cwd=cwd, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, errors='replace', start_new_session=True)
        try:
            for line in process.stdout:
                print(line, end='', flush=True)
                out.write(line); out.flush()
            return process.wait()
        except KeyboardInterrupt:
            os.killpg(process.pid, signal.SIGTERM)
            process.wait()
            raise


def launch(args, root):
    out = args.output.expanduser().resolve()
    for path in (root, out, args.work_dir):
        if path and any(c in str(path) for c in ["'", '\n', '\r', ',']):
            raise ValueError('Workflow paths cannot contain quotes, commas or newlines.')
    if args.executor == 'local' and (args.account or args.partition):
        raise ValueError('--account and --partition require --executor slurm')
    if args.executor == 'slurm' and (args.max_cpus or args.max_memory):
        raise ValueError('For Slurm use --threads/--memory per job and --max_tasks for queue limits; local maxima do not apply.')
    if out.exists() and any(out.iterdir()) and not args.resume:
        raise ValueError(f'{out} is not empty. Use --resume for this run or choose a new output directory.')
    if args.resume and not (out/'.decoi/run.json').exists():
        raise ValueError('No controller run found to resume in this output directory.')
    out.mkdir(parents=True, exist_ok=True)
    state = out/'.decoi'; state.mkdir(exist_ok=True)
    with (state/'controller.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Another DECOI controller is using this output directory.')
        if not (state/'run.json').exists():
            work = args.work_dir.expanduser().resolve() if args.work_dir else state/'work'
            (state/'run.json').write_text(json.dumps(dict(command=args.command, work_dir=str(work), status='PREPARING')))
        return launch_locked(args, root, out, state)


def launch_locked(args, root, out, state):
    test = args.command == 'test'
    if test:
        inputs = out/'inputs'
        if not inputs.exists():
            subprocess.run([sys.executable, str(root/'scripts/create_paired_smoke.py'), str(inputs), '--genome-first'], check=True)
        cfg = read_config(inputs/'config.yaml', root)
    else:
        cfg = read_config(args.config, root)
        if args.study:
            cfg['study_design_file'] = str(args.study.expanduser().resolve())
        if args.genome_dir:
            cfg['reference']['genome_dir'] = str(args.genome_dir.expanduser().resolve())
            cfg['reference']['source'] = 'genomes'
        if args.silva_fasta:
            cfg['reference'].pop('genome_dir', None)
            cfg['reference']['source'] = 'silva'
            cfg['reference']['dada2_silva_fasta'] = str(args.silva_fasta.expanduser().resolve())
        if args.wgs_reference:
            cfg.setdefault('wgs', {}).update(enabled=True, reference_dir=str(args.wgs_reference.expanduser().resolve()))
    ref = cfg['reference']
    genome = ref.get('genome_dir') if ref.get('source') == 'genomes' else None
    linked = cfg.get('wgs', {}).get('reference_dir') if cfg.get('wgs', {}).get('enabled') else None
    if genome and linked:
        raise ValueError('Use raw genomes or a prelinked WGS bundle, not both.')
    if cfg.get('wgs', {}).get('enabled') and not (genome or linked):
        raise ValueError('WGS requires raw genomes or a prelinked reference bundle.')
    study = cfg.get('study_design_file')
    if not study:
        raise ValueError('Provide study_design_file in the configuration or --study.')
    prepared = ref.get('output_dir') if ref.get('source') == 'prepared' else None
    if prepared and (genome or linked):
        raise ValueError('Prepared references cannot be combined with genome overrides')
    silva = None if genome or prepared else ref.get('dada2_silva_fasta')
    fixture = cfg.get('reference_fixtures', {}).get('mitochondria', {})
    for name, path in [('prepared reference', prepared), ('study', study), ('genomes', genome), ('WGS reference', linked), ('SILVA', silva),
                       ('mitochondrial fixture', fixture.get('source_fasta') if fixture.get('enabled') else None)]:
        if path and any(c in str(path) for c in ["'", '\n', '\r']):
            raise ValueError(f'{name} path cannot contain quotes or newlines.')
        if path and not Path(path).exists():
            raise ValueError(f'Missing {name}: {path}')
    available = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    budget = args.max_cpus or available
    threads = args.threads or (1 if test else 8)
    if args.executor == 'local':
        threads = min(threads, budget)
    mem = args.memory or ('4 GB' if test else '32 GB')
    if args.max_memory:
        def size(value):
            number, unit = re.fullmatch(r'([\d.]+)\s*(MB|GB|TB)', value).groups()
            return float(number) * {'MB': 1, 'GB': 1024, 'TB': 1024**2}[unit]
        if size(mem) > size(args.max_memory):
            raise ValueError('--memory exceeds --max_memory; lower the per-task request or raise the local budget.')
    cfg.setdefault('fastq', {})['cpus'] = threads
    cfg.setdefault('wgs', {})['cpus'] = threads
    dada2 = bool(getattr(args, 'run_dada2', False) or (not test and cfg.get('validation', {}).get('run_dada2')))
    tool_checks(dada2=dada2)
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8]
    logs = out/'logs'/run_id; logs.mkdir(parents=True)
    resolved = state/'config.yaml'; save_stable(resolved, yaml.safe_dump(cfg, sort_keys=False))
    params = dict(config=str(resolved), study=study, genome_dir=genome, silva_fasta=silva,
                  wgs_reference=linked, prepared_reference=prepared, outdir=str(out), threads=threads, dada2_threads=threads, run_dada2=dada2)
    validation = cfg.get('validation', {})
    for key, default in [('forward_primer', 'GTGYCAGCMGCCGCGGTAA'), ('reverse_primer', 'GGACTACNVGGGTWTCTAAT')]:
        value = ref.get(key, default)
        if not re.fullmatch('[ACGTRYSWKMBDHVN]+', value, re.I):
            raise ValueError(f'Invalid {key}: use an IUPAC DNA sequence')
        params[key] = value
    for key, source, default in [('trunc_f', 'trunc_len_f', 240), ('trunc_r', 'trunc_len_r', 200),
                                  ('max_ee_f', 'max_ee_f', 2), ('max_ee_r', 'max_ee_r', 2)]:
        value = validation.get(source, default)
        if not isinstance(value, (int, float)) or value < 0:
            raise ValueError(f'Invalid validation setting: {source}')
        if key in ('trunc_f', 'trunc_r') and not isinstance(value, int):
            raise ValueError(f'{source} must be an integer')
        params[key] = value
    params_path = state/'params.json'; save_stable(params_path, json.dumps(params, indent=2)+'\n')
    # JSON quoting is also valid for these Groovy string literals; reject interpolation.
    def literal(value):
        if '$' in value or '\n' in value or '\r' in value:
            raise ValueError('Resource settings cannot contain dollar signs or newlines.')
        return json.dumps(value)
    config = [f'process.executor = {literal(args.executor)}', f'process.memory = {literal(mem)}',
              f'process {{ withName: SIMULATE_STUDY {{ memory = {literal(mem)} }} }}',
              f'process {{ withName: VALIDATE_DADA2 {{ memory = {literal(mem)} }} }}',
              f'executor.queueSize = {args.max_tasks}']
    if args.executor == 'local':
        config += [f'executor.cpus = {budget}']
        if args.max_memory:
            config += [f'executor.memory = {literal(args.max_memory)}']
    else:
        options = '--nodes=1' + (f' --account={args.account}' if args.account else '')
        if args.account and not re.fullmatch(r'[\w.-]+', args.account):
            raise ValueError('Invalid Slurm account')
        config += [f'process.clusterOptions = {literal(options)}', f'process.time = {literal(args.time_limit)}']
        if args.partition:
            config += [f'process.queue = {literal(args.partition)}']
    runtime = state/'resources.config'; save_stable(runtime, '\n'.join(config)+'\n')
    for source in [resolved, params_path, runtime]:
        shutil.copy2(source, logs/source.name)
    work = args.work_dir.expanduser().resolve() if args.work_dir else state/'work'
    previous = json.loads((state/'run.json').read_text()) if (state/'run.json').exists() else None
    if previous and previous['command'] != args.command:
        raise ValueError('Resume must use the original subcommand.')
    if previous and previous['work_dir'] != str(work):
        raise ValueError('Resume must use the original work directory.')
    record = dict(version=__version__, command=args.command, executor=args.executor, threads=threads,
                  memory=mem, work_dir=str(work), run_id=run_id, status='RUNNING')
    save_stable(state/'run.json', json.dumps(record, indent=2)+'\n')
    command = ['nextflow', '-C', str(root/'nextflow.config')+','+str(runtime), '-log', str(logs/'nextflow.log'),
               'run', str(root/'main.nf'), '-ansi-log', 'false', '-params-file', str(params_path),
               '-work-dir', str(work), '-with-trace', str(logs/'trace.tsv'),
               '-with-report', str(logs/'execution.html'), '-with-timeline', str(logs/'timeline.html')]
    if args.resume:
        command.append('-resume')
    print(f'DECOI {__version__} | {args.command} | {args.executor}\nOutput: {out}\n'
          f'Simulation: {threads} CPU(s), {mem}; task limit: {args.max_tasks}\nLogs: {logs}', flush=True)
    try:
        code = stream(command, state, logs/'console.log')
        record['status'] = 'SUCCESS' if code == 0 else 'FAILED'
        record['exit_code'] = code
    except OSError:
        record['status'] = 'FAILED'
        raise
    except KeyboardInterrupt:
        record['status'] = 'INTERRUPTED'
        raise
    finally:
        save_stable(state/'run.json', json.dumps(record, indent=2)+'\n')
        (logs/'run.json').write_text(json.dumps(record, indent=2)+'\n')
    if code:
        print(f'DECOI stopped. Inspect {logs}/console.log and nextflow.log. Reuse this output with --resume.', file=sys.stderr)
    else:
        print(f'Finished. Report: {out}/dataset/mock_dataset/report.html\nView: decoi report -o {shlex.quote(str(out))} --serve', flush=True)
    return code


def standalone(args, root):
    out = args.output.expanduser().resolve()
    cfg = read_config(args.config, root)
    if args.command == 'prepare-reference':
        out.parent.mkdir(parents=True, exist_ok=True)
        cfg['reference']['output_dir'] = str(out)
        if out.exists():
            raise ValueError('Reference output already exists; choose a new directory.')
        state = out.parent/('.'+out.name+'-logs'); state.mkdir(parents=True, exist_ok=True)
    else:
        if not (out/'manifest.json').is_file():
            raise ValueError('add-wgs requires an existing dataset containing manifest.json.')
        state = out/'logs'; state.mkdir(exist_ok=True)
    run_id = uuid.uuid4().hex
    config = state/(run_id+'.yaml'); config.write_text(yaml.safe_dump(cfg))
    command = [sys.executable, str(root/'mock16s_chem.py'), '--config', str(config)]
    if hasattr(args, 'threads'):
        command += ['--threads', str(args.threads)]
    command += [args.command]
    if args.command == 'add-wgs':
        command += ['--output', str(out)]
    print(f'DECOI {__version__} | {args.command}\nOutput: {out}\nLogs: {state}', flush=True)
    return stream(command, out.parent, state/(run_id+'.log'))


def report(args):
    out = args.output.expanduser().resolve()
    candidates = [out/'report.html', out/'dataset/mock_dataset/report.html']
    path = next((p for p in candidates if p.is_file()), None)
    if not path:
        raise ValueError(f'No completed report found under {out}')
    print(f'Report: {path}', flush=True)
    if args.serve:
        if args.port > 65535:
            raise ValueError('Port must be between 1 and 65535')
        handler = partial(http.server.SimpleHTTPRequestHandler, directory=str(path.parent))
        with http.server.ThreadingHTTPServer(('127.0.0.1', args.port), handler) as server:
            url = f'http://localhost:{args.port}/report.html'
            print(f'Open {url} (Ctrl-C to stop; use SSH forwarding for a remote server)', flush=True)
            if not args.no_browser:
                webbrowser.open(url)
            server.serve_forever()
    return 0


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == 'report':
            return report(args)
        root = resources()
        if args.command == 'check':
            print(f'DECOI {__version__}\nPython: {sys.executable}\nResources: {root}', flush=True)
            tool_checks(dada2=True)
            return 0
        if args.command in {'run', 'test'}:
            return launch(args, root)
        return standalone(args, root)
    except KeyboardInterrupt:
        message = 'Report server stopped.' if args.command == 'report' else 'DECOI interrupted; completed Nextflow tasks remain available for --resume.'
        print('\n'+message, file=sys.stderr)
        return 130
    except (ValueError, OSError, subprocess.CalledProcessError, yaml.YAMLError) as exc:
        print(f'DECOI: {exc}', file=sys.stderr)
        return 1
