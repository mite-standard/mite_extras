import pytest
from mite_extras.main import main_cli


def test_main_cli():
    with pytest.raises(SystemExit) as sample:
        main_cli()
    assert sample.value.code == 0
