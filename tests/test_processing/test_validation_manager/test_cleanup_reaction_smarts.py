import pytest
from mite_extras.processing.validation_manager import (
    MoleculeValidator,
    ReactionValidator,
)


@pytest.fixture
def _molecule_validator():
    return MoleculeValidator()


@pytest.fixture
def _reaction_validator():
    return ReactionValidator()


def test_cleanup_reaction_smarts_valid(_molecule_validator, _reaction_validator):
    """Test cleanup of a valid reaction SMARTS string"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    result = _reaction_validator.reaction_cleaner.clean_ketcher_format(
        _reaction_validator.molecule_validator._clean_string(reaction_smarts)
    )
    assert result == reaction_smarts
