import pytest
from mite_extras.processing.validation_manager import ValidationManager
from rdkit.Chem import CanonSmiles, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdChemReactions import ReactionFromSmarts


@pytest.fixture
def validation_manager():
    return ValidationManager()


def test_validate_reaction_smarts_valid(validation_manager):
    """Test validating a reaction SMARTS string with valid expected and forbidden products"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = ["C=O"]
    forbidden_products = []
    result = validation_manager.validate_reaction_smarts(
        reaction_smarts, substrate_smiles, expected_products, forbidden_products
    )
    assert result == None


def test_validate_reaction_smarts_forbidden_in_expected(validation_manager):
    """Test when forbidden products are mistakenly included in expected products"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = ["C=O"]
    forbidden_products = ["C=O"]  # Forbidden product is also in expected products
    with pytest.raises(
        ValueError, match="Some expected products are listed as forbidden products."
    ):
        validation_manager.validate_reaction_smarts(
            reaction_smarts, substrate_smiles, expected_products, forbidden_products
        )


def test_validate_reaction_smarts_unexpected_products(validation_manager):
    """Test when the predicted products do not match expected products"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "OCO |LN:1:1.2|"
    expected_products = ["C=C"]
    forbidden_products = []
    predicted_smiles_set = {
        "O=CO"
    }  # Assume this is what the reaction predicts
    with pytest.raises(
        ValueError,
        match=f"Products '{predicted_smiles_set}' do not meet expectations '{set(expected_products)}'.",
    ):
        validation_manager.validate_reaction_smarts(
            reaction_smarts, substrate_smiles, expected_products, forbidden_products
        )

# TODO (AR 2024-08-15): Proper enumeration not implemented for now
# def test_validate_reaction_smarts_forbidden_products(validation_manager):
#     """Test when forbidden products are found in the reaction output"""
#     reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
#     substrate_smiles = "OCO |LN:1:1.2|"
#     expected_products = ["O=CO"]
#     forbidden_products = ["O=CCO"]  # Forbidden product is also expected
#     with pytest.raises(
#         ValueError, match="Forbidden products were found in the reaction output."
#     ):
#         validation_manager.validate_reaction_smarts(
#             reaction_smarts, substrate_smiles, expected_products, forbidden_products
#         )


def test_validate_reaction_smarts_empty_expected_products(validation_manager):
    """Test when the expected products list is empty"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = []  # Empty expected products list
    forbidden_products = []
    with pytest.raises(ValueError, match="Expected products list cannot be empty."):
        validation_manager.validate_reaction_smarts(
            reaction_smarts, substrate_smiles, expected_products, forbidden_products
        )


def test_validate_reaction_smarts_invalid_expected_product(validation_manager):
    """Test when one of the expected products cannot be parsed"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = ["C=O", "INVALID_SMILES"]  # Invalid product
    forbidden_products = []
    with pytest.raises(
        ValueError, match=f"Could not read SMILES string 'INVALID_SMILES'"
    ):
        validation_manager.validate_reaction_smarts(
            reaction_smarts, substrate_smiles, expected_products, forbidden_products
        )


def test_validate_reaction_smarts_invalid_forbidden_product(validation_manager):
    """Test when one of the forbidden products cannot be parsed"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = ["C=O"]
    forbidden_products = ["INVALID_SMILES"]  # Invalid forbidden product
    with pytest.raises(
        ValueError, match="Could not read SMILES string 'INVALID_SMILES'"
    ):
        validation_manager.validate_reaction_smarts(
            reaction_smarts, substrate_smiles, expected_products, forbidden_products
        )


def test_validate_reaction_smarts_valid_composite(validation_manager):
    """Test when composite SMILES (Rhea)"""
    # TODO (AR 2024-08-15): @MMZ not sure the reaction gives two products on this one
    reaction_smarts = "([cH:9]1[cH:8][c:7]2[cH:6][cH:5][cH:4][cH:3][c:2]2[nH:1]1.[ClH:10])>>[Cl:10][c:9]1[cH:8][cH:7][cH:6][c:5]2[nH:4][cH:3][cH:2][c:1]12.[H:11][H:12]"
    substrate_smiles = "c1cc2ccccc2[nH]1.Cl"
    expected_products = ["Clc1cccc2[nH]ccc12","[H][H]"]
    forbidden_products = []
    result = validation_manager.validate_reaction_smarts(
        reaction_smarts, substrate_smiles, expected_products, forbidden_products
    )
    assert result == None

# TODO (AR 2024-08-21): CXSMILES parsing error
# TODO (AR 2024-08-21): Proper enumeration not implemented for now
# def test_validate_reaction_smarts_position_variation(validation_manager):
#     """Test when position variation"""
#     reaction_smarts = "([cH:2]1[cH:1][c:5]2[cH:6][cH:7][cH:8][cH:9][c:4]2[nH:3]1.[ClH:10])>>[Cl:10]*.[cH:1]1[cH:2][c:3]2[cH:4][cH:5][cH:6][cH:7][c:8]2[nH:9]1 |m:1:5.6|"  # Example SMARTS with position variation
#     substrate_smiles = "c1cc2ccccc2[nH]1.Cl"
#     expected_products = ["c1c2c(Cl)cccc2[nH]c1", "c1c2cc(Cl)ccc2[nH]c1", "c1c2c(Cl)c(Cl)ccc2[nH]c1"]
#     forbidden_products = []
#     result = validation_manager.validate_reaction_smarts(
#         reaction_smarts, substrate_smiles, expected_products, forbidden_products
#     )
#     assert result == None

# TODO (AR 2024-08-21): Proper enumeration not implemented for now
# def test_validate_reaction_smarts_frequency_variation(validation_manager):
#     """Test when position variation"""
#     reaction_smarts = "([NH2:1][C@@H:2]([CH2:3][c:4]1[cH:5][n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]12)[C:13]([OH:14])=[O:15].[ClH:16])>>[NH2:1][C@@H:2]([CH2:3][c:4]1[cH:5][nH:6][c:7]2[cH:8][cH:9][cH:10][c:11]([Cl:12])[c:13]12)[C:14]([OH:15])=[O:16] |r,Sg:n:2:1-2:ht|"  # Example SMARTS with frequency variation
#     substrate_smiles = "Cl.N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O"
#     expected_products = ["N[C@@H](Cc1c[nH]c2cccc(Cl)c12)C(O)=O", "N[C@@H](CCc1c[nH]c2cccc(Cl)c12)C(O)=O"]
#     forbidden_products = []
#     result = validation_manager.validate_reaction_smarts(
#         reaction_smarts, substrate_smiles, expected_products, forbidden_products
#     )
#     assert result == None
