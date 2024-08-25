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


def test_canonicalize_smiles_valid(validation_manager):
    """Test canonicalizing a valid SMILES string"""
    smiles = "CCO"  # Ethanol
    canonical_smiles = validation_manager.canonicalize_smiles(smiles)
    expected_canonical_smiles = CanonSmiles(MolToSmiles(MolFromSmiles(smiles)))
    assert canonical_smiles == expected_canonical_smiles


def test_canonicalize_smiles_invalid(validation_manager):
    """Test canonicalizing an invalid SMILES string"""
    smiles = "CCO@"  # Invalid SMILES (unexpected character '@')
    with pytest.raises(ValueError, match=f"Could not read SMILES string '{smiles}'"):
        validation_manager.canonicalize_smiles(smiles)


# TODO (AR 2024-08-09): @MMZ these are overly dumb but well...
def test_canonicalize_smiles_complex(validation_manager):
    """Test canonicalizing a complex SMILES string"""
    smiles = "C1=CC=CC=C1"  # Benzene ring
    canonical_smiles = validation_manager.canonicalize_smiles(smiles)
    expected_canonical_smiles = CanonSmiles(MolToSmiles(MolFromSmiles(smiles)))
    assert canonical_smiles == expected_canonical_smiles


def test_canonicalize_smiles_with_aromatic(validation_manager):
    """Test canonicalizing a SMILES string with aromatic atoms"""
    smiles = "c1ccccc1"  # Benzene ring in aromatic form
    canonical_smiles = validation_manager.canonicalize_smiles(smiles)
    expected_canonical_smiles = CanonSmiles(MolToSmiles(MolFromSmiles(smiles)))
    assert canonical_smiles == expected_canonical_smiles
