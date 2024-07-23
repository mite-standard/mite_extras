import pytest
from mite_extras.processing.validation_manager import ValidationManager


def test_canonicalize_smiles_valid():
    assert ValidationManager().canonicalize_smiles("CCCCC") == "CCCCC"


def test_canonicalize_smiles_invalid():
    with pytest.raises(ValueError):
        ValidationManager().canonicalize_smiles("üäö")
