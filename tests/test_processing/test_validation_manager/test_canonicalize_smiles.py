import pytest
from mite_extras.processing.validation_manager import MoleculeValidator
from rdkit.Chem import (
    CanonSmiles,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmarts,
    MolToSmiles,
)


@pytest.fixture
def _molecule_validator():
    return MoleculeValidator()


def test_canonicalize_smiles_valid(_molecule_validator):
    """Test canonicalizing a valid SMILES string"""
    smiles = "CCO"  # Ethanol
    canonical_smiles = _molecule_validator.canonicalize_smiles(smiles)
    expected_canonical_smiles = CanonSmiles(MolToSmiles(MolFromSmiles(smiles)))
    assert canonical_smiles == expected_canonical_smiles


def test_canonicalize_smiles_invalid(_molecule_validator):
    """Test canonicalizing an invalid SMILES string"""
    smiles = "CCO@"  # Invalid SMILES (unexpected character '@')
    with pytest.raises(ValueError):
        _molecule_validator.canonicalize_smiles(smiles)


def test_canonicalize_smiles_complex(_molecule_validator):
    """Test canonicalizing a complex SMILES string"""
    smiles = "C1=CC=CC=C1"  # Benzene ring
    canonical_smiles = _molecule_validator.canonicalize_smiles(smiles)
    expected_canonical_smiles = CanonSmiles(MolToSmiles(MolFromSmiles(smiles)))
    assert canonical_smiles == expected_canonical_smiles


def test_canonicalize_smiles_with_aromatic(_molecule_validator):
    """Test canonicalizing a SMILES string with aromatic atoms"""
    smiles = "c1ccccc1"  # Benzene ring in aromatic form
    canonical_smiles = _molecule_validator.canonicalize_smiles(smiles)
    expected_canonical_smiles = CanonSmiles(MolToSmiles(MolFromSmiles(smiles)))
    assert canonical_smiles == expected_canonical_smiles
