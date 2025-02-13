import pytest
from mite_extras.processing.validation_manager import MoleculeValidator


@pytest.fixture
def _molecule_validator():
    return MoleculeValidator()


def test_unescape_string(_molecule_validator):
    """Test canonicalizing a valid SMARTS string"""
    string = r"This\\is\a\\test\\string\\with\\double\backslashes"
    unescaped_string = _molecule_validator._clean_string(string)
    truth = r"This\is\a\test\string\with\double\backslashes"
    assert unescaped_string == truth
