import pytest
from mite_extras.processing.validation_manager import ValidationManager


@pytest.fixture
def validation_manager():
    return ValidationManager()


def test_cleanup_reaction_smarts_valid(validation_manager):
    """Test cleanup of a valid reaction SMARTS string"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    result = validation_manager.cleanup_reaction_smarts(reaction_smarts)
    assert result == reaction_smarts
