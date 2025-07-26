import pytest
from mite_extras.processing.validation_manager import ReactionCleaner


def test_chloride_w_indexing():
    assert ReactionCleaner.clean_ketcher_format("[#6:1]-Cl:2>>[#6:1](-Cl)-Cl:2") == (
        "[#6:1]-[Cl:2]>>[#6:1](-[Cl])-[Cl:2]"
    )


def test_chloride_wo_indexing():
    assert ReactionCleaner.clean_ketcher_format("[#6:1]>>[#6:1]-Cl") == (
        "[#6:1]>>[#6:1]-[Cl]"
    )


def test_chloride_substitute_w_indexing():
    assert ReactionCleaner.clean_ketcher_format(
        "[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1-Cl:7>>[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2](-Cl)-[#6:3]=1-Cl:7"
    ) == (
        "[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1-[Cl:7]>>[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2](-[Cl])-[#6:3]=1-[Cl:7]"
    )


def test_chloride_substitute_wo_indexing():
    assert ReactionCleaner.clean_ketcher_format(
        "[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1>>[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1-Cl"
    ) == (
        "[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1>>[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1-[Cl]"
    )


def test_nitrogen_heterocycle_1():
    assert ReactionCleaner.clean_ketcher_format(
        "[#7:1;h1]1:[#6:3]:[#6:2]:[#6:4](-[#7]-[#6]):[#6:5]:1>>[#7:1;h1]1:[#6:5]:[#6:4]:[#6:2](-[#7]-[#6]):[#6:3]:1"
    ) == (
        "[nH1:1]1:[#6:3]:[#6:2]:[#6:4](-[#7]-[#6]):[#6:5]:1>>[nH1:1]1:[#6:5]:[#6:4]:[#6:2](-[#7]-[#6]):[#6:3]:1"
    )


def test_nitrogen_heterocycle_2():
    assert ReactionCleaner.clean_ketcher_format(
        "[#7;h1]1:[#6:3]:[#6:2]:[#6:4](-[#7]-[#6]):[#6:5]:1>>[#7;h1]1:[#6:5]:[#6:4]:[#6:2](-[#7]-[#6]):[#6:3]:1"
    ) == (
        "[nH1]1:[#6:3]:[#6:2]:[#6:4](-[#7]-[#6]):[#6:5]:1>>[nH1]1:[#6:5]:[#6:4]:[#6:2](-[#7]-[#6]):[#6:3]:1"
    )


def test_numbered_charges_1():
    assert ReactionCleaner.clean_ketcher_format(
        "[#6:1]1:[#6:5]:[#6:6]:[#6:4]:[#6:2](-[#6:7]):[#6:3]:1-[#7:8;+](=[#8:10])-[#8:9;-]>>[#6:1]1:[#6:5]:[#6:6]:[#6:4]:[#6:2](-[#6:7]-[#8]):[#6:3]:1-[#7:8;+](=[#8:10])-[#8:9;-]"
    ) == (
        "[#6:1]1:[#6:5]:[#6:6]:[#6:4]:[#6:2](-[#6:7]):[#6:3]:1-[#7;+:8](=[#8:10])-[#8;-:9]>>[#6:1]1:[#6:5]:[#6:6]:[#6:4]:[#6:2](-[#6:7]-[#8]):[#6:3]:1-[#7;+:8](=[#8:10])-[#8;-:9]"
    )


def test_expl_hydro_valid():
    assert (
        ReactionCleaner.detect_undesired_smarts(
            "[#6:1]1-[#6:2]-[#6:3]-[#6:4]-[#6:5]-[#6:6]-1>>[#6:1]1-[#6:6]-[#6:5](-[#8])-[#6:4]-[#6:3]-[#6:2]-1"
        )
        is None
    )


def test_expl_hydro_invalid():
    with pytest.raises(ValueError):
        assert ReactionCleaner.detect_undesired_smarts(
            "[#6:1]-[#6:2]-[#6:3]-[#7:4]>>[#6:2](-[#6@:3](-[#7:4])(/[H])\\[#8])-[#6:1]"
        )


def test_cxsmarts_valid():
    assert (
        ReactionCleaner.detect_undesired_smarts(
            "[#6:1]1-[#6:2]-[#6:3]-[#6:4]-[#6:5]-[#6:6]-1>>[#6:1]1-[#6:6]-[#6:5](-[#8])-[#6:4]-[#6:3]-[#6:2]-1"
        )
        is None
    )


def test_cxsmarts_invalid():
    with pytest.raises(ValueError):
        assert ReactionCleaner.detect_undesired_smarts(
            "[NH2:1][C@@H:2]([CH2:3][c:4]1[cH:5][n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]12)[C:13]([OH:14])=[O:15].[ClH:16]>>[NH2:1][C@@H:2]([CH2:3][c:4]1[cH:5][nH:6][c:7]2[cH:8][cH:9][cH:10][c:11]([Cl:12])[c:13]12)[C:14]([OH:15])=[O:16] |r,Sg:n:2:1-2:ht|"
        )
