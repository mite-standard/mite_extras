import pytest

from mite_extras.processing.validation_manager import (
    ReactionValidator,
)


@pytest.fixture
def _reaction_validator():
    return ReactionValidator()


def test_validate_reaction_valid(_reaction_validator):
    """Test validating a reaction SMARTS string with valid expected and forbidden products"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = ["C=O"]
    forbidden_products = []
    result = _reaction_validator.validate_reaction(
        reaction_smarts,
        substrate_smiles,
        expected_products,
        forbidden_products,
    )
    assert result == None


def test_validate_reaction_forbidden_in_expected(_reaction_validator):
    """Test when forbidden products are mistakenly included in expected products"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = ["C=O"]
    forbidden_products = ["C=O"]  # Forbidden product is also in expected products
    with pytest.raises(ValueError):
        _reaction_validator.validate_reaction(
            reaction_smarts,
            substrate_smiles,
            expected_products,
            forbidden_products,
        )


def test_validate_reaction_unexpected_products(_reaction_validator):
    """Test when the predicted products do not match expected products"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "OCO"
    expected_products = ["C=C"]
    forbidden_products = []
    predicted_smiles_set = {"O=CO"}  # Assume this is what the reaction predicts
    with pytest.raises(ValueError):
        _reaction_validator.validate_reaction(
            reaction_smarts,
            substrate_smiles,
            expected_products,
            forbidden_products,
        )


def test_validate_reaction_forbidden_products(_reaction_validator):
    """Test when forbidden products are found in reaction output"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "OCO |LN:1:1.2|"
    expected_products = ["O=CO"]
    forbidden_products = ["O=CCO"]  # Forbidden product is also expected
    with pytest.raises(ValueError):
        _reaction_validator.validate_reaction(
            reaction_smarts,
            substrate_smiles,
            expected_products,
            forbidden_products,
        )


def test_validate_reaction_empty_expected_products(_reaction_validator):
    """Test when the expected products list is empty"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = []  # Empty expected products list
    forbidden_products = []
    with pytest.raises(ValueError):
        _reaction_validator.validate_reaction(
            reaction_smarts,
            substrate_smiles,
            expected_products,
            forbidden_products,
        )


def test_validate_reaction_invalid_expected_product(_reaction_validator):
    """Test when one of the expected products cannot be parsed"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = ["C=O", "INVALID_SMILES"]  # Invalid product
    forbidden_products = []
    with pytest.raises(ValueError):
        _reaction_validator.validate_reaction(
            reaction_smarts,
            substrate_smiles,
            expected_products,
            forbidden_products,
        )


def test_validate_reaction_invalid_forbidden_product(_reaction_validator):
    """Test when one of the forbidden products cannot be parsed"""
    reaction_smarts = "[C:1][O:2]>>[C:1]=[O:2]"
    substrate_smiles = "CO"
    expected_products = ["C=O"]
    forbidden_products = ["INVALID_SMILES"]  # Invalid forbidden product
    with pytest.raises(ValueError):
        _reaction_validator.validate_reaction(
            reaction_smarts,
            substrate_smiles,
            expected_products,
            forbidden_products,
        )


def test_validate_reaction_valid_composite_pass(_reaction_validator):
    """Test when composite SMILES (Rhea)"""
    # TODO (AR 2024-08-15): @MMZ not sure the reaction gives two products on this one
    reaction_smarts = "[cH:9]1[cH:8][c:7]2[cH:6][cH:5][cH:4][cH:3][c:2]2[nH:1]1.[ClH:10]>>[Cl:10][c:9]1[cH:8][cH:7][cH:6][c:5]2[nH:4][cH:3][cH:2][c:1]12"
    substrate_smiles = "c1cc2ccccc2[nH]1.Cl"
    expected_products = ["Clc1cccc2[nH]ccc12"]
    forbidden_products = []
    result = _reaction_validator.validate_reaction(
        reaction_smarts,
        substrate_smiles,
        expected_products,
        forbidden_products,
    )
    assert result == None


def test_validate_reaction_halogenation(_reaction_validator):
    """Test when multiple atoms allowed (halogenation)"""
    reaction_smarts = "[#7:1]-[#6:2](-[#6:3]-[c:4]1[c:5][nH1:6][c:7]2[c:8][c:9][c:10][c:11][c:12]12)-[#6:13](-[#8:14])=[O:15]>>[#7:1]-[#6:2](-[#6:3]-[c:4]1[c:5][nH1:6][c:7]2[c:8][c:9](-[#17,#35:10])[c:11][c:12][c:13]12)-[#6:14](-[#8:15])=[O:16]"  # Example SMARTS with Cl OR Br
    substrate_smiles = "NC(Cc1c[nH]c2ccccc12)C(=O)O"
    expected_products = [
        "NC(Cc1c[nH]c2cc(Br)ccc12)C(=O)O",
        "NC(Cc1c[nH]c2cc(Cl)ccc12)C(=O)O",
    ]
    forbidden_products = []
    result = _reaction_validator.validate_reaction(
        reaction_smarts,
        substrate_smiles,
        expected_products,
        forbidden_products,
    )
    assert result == None


def test_validate_intramolecular_macrolactam(_reaction_validator):
    reaction_smarts = "([#7&H2:1]-[#6&H2:2]-[#6:3](-[#7&H1:5])=[#8:4].[#7&H1:6]-[#6:7](-[#6:13](-[#7&H1:15])=[#8:14])-[#6:8]-[#6:9]-[#6:10](-[#8&H1:12])=[#8:11])>>[#7:1](-[#8:12]-[#6:10](=[#8:11])-[#6:9]-[#6:8]-[#6:7](-[#6:13](-[#7:15])=[#8:14])-[#7:6])-[#6:2]-[#6:3](-[#7:5])=[#8:4]"
    substrate_smiles = (
        "C[C@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)CN)C(=O)O"
    )
    expected_products = [
        "C[C@H](NC(=O)[C@@H]1CCC(=O)ONCC(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N1)C(=O)O",
    ]
    forbidden_products = []
    result = _reaction_validator.validate_reaction(
        reaction_smarts,
        substrate_smiles,
        expected_products,
        forbidden_products,
        intramolecular=True,
    )
    assert result == None


def test_validate_cheatoglobosin_1(_reaction_validator):
    r_smarts = "[#6:1]1(:[#6:36]:[#7:35;h1]:[#6:30]2:[#6:31]:[#6:32]:[#6:33]:[#6:34]:[#6:29]:1:2)-[#6:2]\\[#6@:3]1-[#6@:13]2-[#6@@:7]3(-[#6:27](=[#8:28])-[#6:26]=[#6:25]-[#6:24]-[#6:23]-[#6:21](-[#6:22])=[#6:20]-[#6@@:18](/[#6:19])-[#6:17]-[#6:16]=[#6:15]/[#6@:8]-3-[#6:9]=[#6:10](-[#6:14])-[#6@:11]/2/[#6:12])/[#6:5](=[#8:6])-[#7:4]-1>>[#6:1]1(:[#6:36]:[#7:35;h1]:[#6:30]2:[#6:31]:[#6:32]:[#6:33]:[#6:34]:[#6:29]:1:2)-[#6:2]\\[#6@:3]1-[#6@:13]2-[#6@@:7]3(-[#6:27](=[#8:28])-[#6:26]=[#6:25]-[#6:24]-[#6:23]-[#6:21](-[#6:22])=[#6:20]-[#6@@:18](/[#6:19])-[#6:17]-[#6:16]=[#6:15]/[#6@:8]-3-[#6@@:9]3\\[#8]-[#6@:10]-3(/[#6:14])-[#6@:11]/2/[#6:12])/[#6:5](=[#8:6])-[#7:4]-1"
    substr = "c1(C[C@@H]2NC(=O)[C@@]34[C@H](C=C([C@@H](C)C23)C)C=CC[C@H](C)C=C(C)CCC=CC4=O)c2c(cccc2)[nH]c1 |c:21,t:16,26|"
    prods = [
        "c1(C[C@H]2C3[C@H](C)[C@]4(C)[C@@H](O4)[C@H]4[C@@]3(C(N2)=O)C(C=CCCC(C)=C[C@H](CC=C4)C)=O)c2c(cccc2)[nH]c1",
    ]
    result = _reaction_validator.validate_reaction(
        r_smarts,
        substr,
        prods,
        [],
        intramolecular=False,
    )
    assert result == None
    substr = "c1(C[C@@H]2NC(=O)[C@@]34[C@H](C=C([C@@H](C)C23)C)C=CC[C@H](C)C=C(C)[C@@H](O)[C@@H](O)C=CC4=O)c2c(cccc2)[nH]c1 |c:21,t:16,28|"
    prods = [
        "c1(C[C@@H]2NC(=O)[C@@]34[C@H]([C@H]5[C@]([C@@H](C)C23)(C)O5)C=CC[C@H](C)C=C(C)[C@@H](O)[C@@H](O)C=CC4=O)c2c(cccc2)[nH]c1",
    ]
    result = _reaction_validator.validate_reaction(
        r_smarts,
        substr,
        prods,
        [],
        intramolecular=False,
    )
    assert result == None
    substr = "c1(C[C@@H]2NC(=O)[C@@]34[C@H](C=C([C@@H](C)C23)C)C=CC[C@H](C)C=C(C)[C@@H](O)C(=O)C=CC4=O)c2c(cccc2)[nH]c1 |c:21,t:16,28|"
    prods = [
        "c1(C[C@@H]2NC(=O)[C@@]34[C@H]([C@H]5[C@]([C@@H](C)C23)(C)O5)C=CC[C@H](C)C=C(C)[C@@H](O)C(=O)C=CC4=O)c2c(cccc2)[nH]c1",
    ]
    result = _reaction_validator.validate_reaction(
        r_smarts,
        substr,
        prods,
        [],
        intramolecular=False,
    )
    assert result == None


def test_validate_macrolacton(_reaction_validator):
    r_smarts = "[#8:1]1-[#6:8]-[#6@@:7](/[#8:11])-[#6@@:6](/[#8:10])-[#6:5]-[#6:4]-[#6:3]-[#6:2]-1=[#8:9]>>[#8:1]1-[#6:8]-[#6@@:7](/[#8:11])-[#6@:6](\\[#8:10])-[#6:5]-[#6:4]-[#6:3]-[#6:2]-1=[#8:9]"
    substr = "O1C[C@@H](O)[C@@H](O)CCCC1=O"
    prods = ["O1C[C@@H](O)[C@H](O)CCCC1=O"]
    result = _reaction_validator.validate_reaction(
        r_smarts,
        substr,
        prods,
        [],
        intramolecular=False,
    )
    assert result == None


def test_validate_aromatic_n(_reaction_validator):
    """Checks removal of hydrogen on aromatic N in tryptophan when it becomes substituted with another bond

    Ketcher exports a malformed SMILES string C[C@@H](C(N[C@@H]1C(N[C@H](C(O)=O)Cc2c3c(cccc3)[nH](C1c1c3c(cccc3)[nH]c1)c2)=O)=O)N
    The manually corrected SMILES string is C[C@@H](C(N[C@@H]1C(N[C@H](C(O)=O)Cc2c3c(cccc3)[n](C1c1c3c(cccc3)[nH]c1)c2)=O)=O)N
    """
    r_smarts = "[#6:1]\\[#6@@:2](-[#6:4](-[#7:6]\\[#6@:7](-[#6:18](-[#7:20]/[#6@:21](-[#6:32](-[#8:34])=[#8:33])-[#6:22]-[#6:23]1:[#6:31]2:[#6:26](:[#6:27]:[#6:28]:[#6:29]:[#6:30]:2):[#7:25;h1]:[#6:24]:1)=[#8:19])-[#6:8]-[#6:9]1:[#6:17]2:[#6:12](:[#6:13]:[#6:14]:[#6:15]:[#6:16]:2):[#7:11;h1]:[#6:10]:1)=[#8:5])-[#7:3]>>[#6:1]\\[#6@@:2](-[#6:4](-[#7:6]\\[#6@@:7]1-[#6:18](-[#7:20]/[#6@:21](-[#6:32](-[#8:34])=[#8:33])-[#6:22]-[#6:23]2:[#6:31]3:[#6:26](:[#6:27]:[#6:28]:[#6:29]:[#6:30]:3):[#7:25;h1](-[#6:8]-1-[#6:9]1:[#6:17]3:[#6:12](:[#6:13]:[#6:14]:[#6:15]:[#6:16]:3):[#7:11;h1]:[#6:10]:1):[#6:24]:2)=[#8:19])=[#8:5])-[#7:3]"
    substr = (
        "C[C@H](N)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"
    )
    prods = [
        "C[C@@H](C(N[C@@H]1C(N[C@H](C(O)=O)Cc2c3c(cccc3)[nH](C1c1c3c(cccc3)[nH]c1)c2)=O)=O)N",
    ]
    result = _reaction_validator.validate_reaction(
        r_smarts,
        substr,
        prods,
        [],
        intramolecular=False,
    )
    assert result == None


def test_valid_stereochem_cyclization(_reaction_validator):
    """This reaction with flattened substrates from MITE0000378 succeeds in the current setup"""
    r_smarts = "[#6:1]-[#6:2](-[#6:17]-[#6:18]-[#6:19]=[#6:20](-[#6:22]-[#6:23]\\[#6:24]1-[#6:26](-[#6:28])(-[#6:27])-[#8:25]-1)-[#6:21])=[#6:3]-[#6:4]-[#6:5]1:[#6:15](-[#8:16]):[#6:12](-[#6:13]=[#8:14]):[#6:10](-[#6:11]):[#6:8](-[Cl:9]):[#6:6]:1-[#8:7]>>[#6:1]-[#6:2](-[#6:17]-[#6:18]-[#6:19]1(\\[#6:27])-[#6:26](\\[#6:28])-[#6&H0:24](=[#8:25])-[#6:23]-[#6:22]-[#6:20]-1\\[#6:21])=[#6:3]-[#6:4]-[#6:5]1:[#6:15](-[#8:16]):[#6:12](-[#6:13]=[#8:14]):[#6:10](-[#6:11]):[#6:8](-[Cl:9]):[#6:6]:1-[#8:7]"
    substr = "CC(=CCc1c(O)c(Cl)c(C)c(C=O)c1O)CCC=C(C)CCC1OC1(C)C"
    prods = [
        "CC(=CCc1c(O)c(Cl)c(C)c(C=O)c1O)CCC1(C)C(C)CCC(=O)C1C",
    ]
    result = _reaction_validator.validate_reaction(
        r_smarts,
        substr,
        prods,
        [],
        intramolecular=False,
    )
    assert result == None


def test_failing_w_rdkit_2026_3_3(_reaction_validator):
    """MITE0000346: rdkit produces enantiomer of expected entry in 2026.3.3"""
    r_smarts = "[#6:1]-[#6:2](-[#6:4]-[#6:5]-[#6:6]-[#6:7]-[#6:8]-[#6:9]-[#6:10]-[#6:11]-[#6:12]-[#6:13]1:[#7:30]:[#6:16](-[#6:17]=[#6:18]2-[#7:29]=[#6:23](-[#6:24]3:[#7:28]:[#6:27]:[#6:26]:[#6:25]:3)-[#6:22]=[#6:19]-2-[#8:20]-[#6:21]):[#6:15]:[#6:14]:1)-[#8:3]>>[#6:1]/[#6@:2]1-[#8:3]-[#6:18]2(-[#7:29]=[#6:23](-[#6:24]3-[#7:28]-[#6:27]=[#6:26]-[#6:25]=3)-[#6:22]=[#6:19]-2-[#8:20]-[#6:21])-[#6@:17]2/[#6:16]3=[#6:15]-[#6:14]=[#6:13](-[#6:12]-[#6:11]-[#6:10]-[#6:9]-[#6:8]-[#6:7]-[#6:6]\\[#6@@:5]-2-[#6:4]-1)-[#7:30]-3"
    substr = "COC1=CC(c2ccc[nH]2)=N/C1=C/c1ccc(CCCCCCCCCC(C)O)[nH]1"
    prods = [
        "COC1=CC(c2ccc[nH]2)=NC12O[C@H](C)C[C@H]1CCCCCCCc3ccc([nH]3)[C@@H]12",
    ]
    result = _reaction_validator.validate_reaction(
        r_smarts,
        substr,
        prods,
        [],
        intramolecular=False,
    )
    assert result == None
