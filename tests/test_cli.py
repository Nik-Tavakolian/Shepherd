"""Tests for the ``shepherd`` command line interface."""
import subprocess
import sys

import pytest

from shepherd import __version__
from shepherd.cli import main


def test_python_m_shepherd_runs_the_cli():
    result = subprocess.run([sys.executable, '-m', 'shepherd', '--help'], capture_output=True, text=True)
    assert result.returncode == 0
    assert 'cluster' in result.stdout and 'track' in result.stdout


def test_version(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(['--version'])
    assert exit_info.value.code == 0
    assert __version__ in capsys.readouterr().out


@pytest.mark.parametrize('argv', [[], ['cluster', '-f', 'reads.txt'], ['track', '-f0', 'reads_t0.txt']])
def test_missing_required_arguments_exit_with_usage_error(argv):
    with pytest.raises(SystemExit) as exit_info:
        main(argv)
    assert exit_info.value.code == 2
