import pytest
from mite_extras.processing.validation_manager import ValidationManager


@pytest.fixture
def validation_manager():
    return ValidationManager()


def test_unescape_string(validation_manager):
    """Test canonicalizing a valid SMARTS string"""
    string = r"This\\is\a\\test\\string\\with\\double\backslashes"
    unescaped_string = validation_manager.unescape_string(string)
    truth = r"This\is\a\test\string\with\double\backslashes"
    assert unescaped_string == truth
