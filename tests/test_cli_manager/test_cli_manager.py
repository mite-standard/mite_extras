from mite_extras.cli.cli_manager import CliManager


def test_init_valid():
    assert isinstance(CliManager(), CliManager)
