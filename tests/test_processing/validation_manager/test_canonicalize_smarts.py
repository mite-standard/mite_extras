import pytest
from mite_extras.processing.validation_manager import ValidationManager
from rdkit.Chem import (
    CanonSmiles,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmarts,
    MolToSmiles,
)


@pytest.fixture
def validation_manager():
    return ValidationManager()


def test_canonicalize_smarts_valid(validation_manager):
    """Test canonicalizing a valid SMARTS string"""
    smarts = "[C;H2][O;H1]"
    canonical_smarts = validation_manager.canonicalize_smarts(smarts)
    expected_canonical_smarts = MolToSmarts(
        MolFromSmiles(CanonSmiles(MolToSmiles(MolFromSmarts(smarts))))
    )
    assert canonical_smarts == expected_canonical_smarts


# TODO (AR 2024-08-09): @MMZ these are overly dumb but well...
def test_canonicalize_smarts_complex(validation_manager):
    """Test canonicalizing a complex SMARTS string"""
    smarts = "[#6]1[#6][#6][#6][#6][#6]1"  # Benzene ring
    canonical_smarts = validation_manager.canonicalize_smarts(smarts)
    expected_canonical_smarts = MolToSmarts(
        MolFromSmiles(CanonSmiles(MolToSmiles(MolFromSmarts(smarts))))
    )
    assert canonical_smarts == expected_canonical_smarts


def test_canonicalize_smarts_with_aromatic(validation_manager):
    """Test canonicalizing a SMARTS string with aromatic atoms"""
    smarts = "c1ccccc1"  # Benzene ring in aromatic form
    canonical_smarts = validation_manager.canonicalize_smarts(smarts)
    expected_canonical_smarts = MolToSmarts(
        MolFromSmiles(CanonSmiles(MolToSmiles(MolFromSmarts(smarts))))
    )
    assert canonical_smarts == expected_canonical_smarts
