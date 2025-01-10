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


def test_canonicalize_smarts_valid(_molecule_validator):
    """Test canonicalizing a valid SMARTS string"""
    smarts = "[C;H2][O;H1]"
    canonical_smarts = _molecule_validator.canonicalize_smarts(smarts)
    mol = MolFromSmarts(smarts)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(i)
    expected_canonical_smarts = MolToSmarts(
        MolFromSmiles(CanonSmiles(MolToSmiles(mol)))
    )
    assert canonical_smarts == expected_canonical_smarts


# TODO (AR 2024-08-09): @MMZ these are overly dumb but well...
def test_canonicalize_smarts_complex(_molecule_validator):
    """Test canonicalizing a complex SMARTS string"""
    smarts = "[#6]1[#6][#6][#6][#6][#6]1"  # Benzene ring
    canonical_smarts = _molecule_validator.canonicalize_smarts(smarts)
    mol = MolFromSmarts(smarts)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(i)
    expected_canonical_smarts = MolToSmarts(
        MolFromSmiles(CanonSmiles(MolToSmiles(mol)))
    )
    assert canonical_smarts == expected_canonical_smarts


def test_canonicalize_smarts_with_aromatic(_molecule_validator):
    """Test canonicalizing a SMARTS string with aromatic atoms"""
    smarts = "c1ccccc1"  # Benzene ring in aromatic form
    canonical_smarts = _molecule_validator.canonicalize_smarts(smarts)
    mol = MolFromSmarts(smarts)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(i)
    expected_canonical_smarts = MolToSmarts(
        MolFromSmiles(CanonSmiles(MolToSmiles(mol)))
    )
    assert canonical_smarts == expected_canonical_smarts
