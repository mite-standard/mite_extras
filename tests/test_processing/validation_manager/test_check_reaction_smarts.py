import pytest
from mite_extras.processing.validation_manager import ValidationManager
from rdkit.Chem.rdChemReactions import ReactionFromSmarts


@pytest.fixture
def validation_manager():
    return ValidationManager()


def test_check_reaction_smarts_valid(validation_manager):
    """Test checking a valid reaction SMARTS string"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"  # Example valid reaction SMARTS
    result = validation_manager.check_reaction_smarts(reaction_smarts)
    assert result == reaction_smarts
