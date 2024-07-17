from mite_extras.cli.cli_manager import CliManager


def test_init_valid():
    instance = CliManager()
    assert instance.version_program is not None
    assert instance.version_schema is not None
