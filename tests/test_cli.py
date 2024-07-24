import numpy as np
import importlib
import rombus
from os.path import isfile
from click.testing import CliRunner

import warnings

warnings.simplefilter("ignore", np.ComplexWarning)


def test_cli_help(tmp_path):
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        result = runner.invoke(
            rombus.cli.cli,
            [
                "--help",
            ],
        )
        assert result.exit_code == 0


def test_cli_version(tmp_path):
    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):
        result = runner.invoke(
            rombus.cli.cli,
            [
                "--version",
            ],
        )
        assert result.exit_code == 0
        assert result.output == f"cli, version {rombus.__version__}\n"


def test_cli_quickstart(tmp_path):
    runner = CliRunner()
    test_project_name = "test_project_name"
    with runner.isolated_filesystem(temp_dir=tmp_path):
        result = runner.invoke(
            rombus.cli.cli,
            [
                "quickstart",
                f"{test_project_name}",
            ],
        )
        assert result.exit_code == 0
        assert isfile(f"{test_project_name}.py")
        assert isfile(f"{test_project_name}_samples.csv")


def test_cli_end_to_end(tmp_path):

    atol = 1e-6
    rtol = 1e-10

    runner = CliRunner()
    test_model = "sinc"
    with runner.isolated_filesystem(temp_dir=tmp_path):

        greedy_filename = str(
            importlib.resources.files("rombus.models").joinpath(
                f"{test_model}_samples.csv"
            )
        )
        result = runner.invoke(
            rombus.cli.cli,
            [
                "build",
                f"rombus.models.{test_model}:Model",
                greedy_filename,
            ],
        )
        assert result.exit_code == 0
        assert isfile(f"{test_model}.hdf5")

        result = runner.invoke(
            rombus.cli.cli,
            [
                "evaluate",
                f"{test_model}.hdf5",
                "A=2",
            ],
        )

        assert result.exit_code == 0
        assert isfile(f"{test_model}_comparison.pdf")
        assert isfile(f"{test_model}_ROM_diff.npy")
        delta_diff = np.load(f"{test_model}_ROM_diff.npy")
        assert np.allclose(delta_diff, 0.0, rtol=rtol, atol=atol)

        result = runner.invoke(
            rombus.cli.cli,
            [
                "evaluate",
                f"{test_model}.hdf5",
                "A=3.5",
            ],
        )
        assert result.exit_code == 0
        assert isfile(f"{test_model}_ROM_diff.npy")
        delta_diff = np.load(f"{test_model}_ROM_diff.npy")
        assert not np.allclose(delta_diff, 0.0, atol=atol)

        result = runner.invoke(
            rombus.cli.cli,
            [
                "refine",
                f"{test_model}.hdf5",
            ],
        )
        assert result.exit_code == 0
        assert isfile(f"{test_model}_refined.hdf5")

        result = runner.invoke(
            rombus.cli.cli,
            [
                "evaluate",
                f"{test_model}_refined.hdf5",
                "A=3.5",
            ],
        )
        assert result.exit_code == 0
        delta_diff = np.load(f"{test_model}_refined_ROM_diff.npy")
        assert np.allclose(delta_diff, 0.0, atol=atol)
