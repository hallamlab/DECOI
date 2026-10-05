import json
from pathlib import Path
import pytest
import yaml
from decoi import cli


@pytest.fixture
def study(tmp_path, monkeypatch):
    root = cli.resources()
    (tmp_path/'genomes').mkdir()
    (tmp_path/'study.yaml').write_text('study_name: test\n')
    config = tmp_path/'settings.yaml'
    config.write_text(yaml.safe_dump({'reference': {'source': 'genomes', 'genome_dir': 'genomes'},
                                     'study_design_file': 'study.yaml', 'wgs': {'enabled': True}}))
    monkeypatch.setattr(cli, 'tool_checks', lambda **kw: None)
    monkeypatch.setattr(cli, 'stream', lambda *args: 0)
    return root, config, tmp_path/'output'


def test_run_resolves_yaml_paths_caps_threads_and_preserves_resume_inputs(study):
    root, config, out = study
    argv = ['run', '-c', str(config), '-o', str(out), '--threads', '8', '--max_cpus', '2']
    assert cli.main(argv) == 0
    params = json.loads((out/'.decoi/params.json').read_text())
    assert params['threads'] == 2
    assert params['genome_dir'] == str(config.parent/'genomes')
    assert params['study'] == str(config.parent/'study.yaml')
    stamp = (out/'.decoi/config.yaml').stat().st_mtime_ns
    assert cli.main(argv) == 1  # Do not overwrite an existing run silently.
    assert cli.main(argv+['--resume']) == 0
    assert (out/'.decoi/config.yaml').stat().st_mtime_ns == stamp
    assert len(list((out/'logs').glob('*/config.yaml'))) == 2


def test_slurm_requests_single_node_and_bounded_queue(study):
    _, config, out = study
    assert cli.main(['run','-c',str(config),'-o',str(out),'--executor','slurm','--account','project-1',
                     '--partition','compute','--max_tasks','3','--threads','8','--memory','16 GB']) == 0
    text = (out/'.decoi/resources.config').read_text()
    assert '--nodes=1 --account=project-1' in text
    assert 'executor.queueSize = 3' in text
    assert 'executor.cpus' not in text
    assert 'process.queue = "compute"' in text


def test_failed_preflight_can_be_resumed(study, monkeypatch):
    _, config, out = study
    def fail(**kw): raise ValueError('Missing dependency')
    monkeypatch.setattr(cli, 'tool_checks', fail)
    argv=['run','-c',str(config),'-o',str(out)]
    assert cli.main(argv)==1
    monkeypatch.setattr(cli, 'tool_checks', lambda **kw: None)
    assert cli.main(argv+['--resume'])==0


def test_memory_budget_and_conflicting_references_rejected(study):
    _, config, out = study
    assert cli.main(['run','-c',str(config),'-o',str(out),'--max_memory','4 GB']) == 1
    cfg=yaml.safe_load(config.read_text());cfg['wgs']['reference_dir']='genomes';config.write_text(yaml.safe_dump(cfg))
    assert cli.main(['run','-c',str(config),'-o',str(out),'--resume'])==1


def test_changed_work_directory_is_not_silently_resumed(study):
    _, config, out = study
    argv=['run','-c',str(config),'-o',str(out)]
    assert cli.main(argv)==0
    assert cli.main(argv+['--resume','--work_dir',str(out/'other')])==1


def test_packaged_resources_and_report_location(tmp_path, monkeypatch):
    fake=tmp_path/'site-packages/decoi/cli.py';fake.parent.mkdir(parents=True);fake.touch()
    runtime=tmp_path/'share/decoi';runtime.mkdir(parents=True);(runtime/'main.nf').touch()
    monkeypatch.setattr(cli,'__file__',str(fake));monkeypatch.setattr(cli.sys,'prefix',str(tmp_path))
    assert cli.resources()==runtime
    report=tmp_path/'result/dataset/mock_dataset/report.html';report.parent.mkdir(parents=True);report.write_text('report')
    assert cli.main(['report','-o',str(tmp_path/'result')])==0
