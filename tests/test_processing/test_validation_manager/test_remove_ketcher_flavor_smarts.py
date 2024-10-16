from mite_extras.processing.validation_manager import ValidationManager


def test_chloride_w_indexing():
    assert ValidationManager.remove_ketcher_flavor_smarts(
        "[#6:1]-Cl:2>>[#6:1](-Cl)-Cl:2"
    ) == ("[#6:1]-[Cl:2]>>[#6:1](-[Cl])-[Cl:2]")


def test_chloride_wo_indexing():
    assert ValidationManager.remove_ketcher_flavor_smarts("[#6:1]>>[#6:1]-Cl") == (
        "[#6:1]>>[#6:1]-[Cl]"
    )


def test_chloride_substitute_w_indexing():
    assert ValidationManager.remove_ketcher_flavor_smarts(
        "[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1-Cl:7>>[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2](-Cl)-[#6:3]=1-Cl:7"
    ) == (
        "[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1-[Cl:7]>>[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2](-[Cl])-[#6:3]=1-[Cl:7]"
    )


def test_chloride_substitute_wo_indexing():
    assert ValidationManager.remove_ketcher_flavor_smarts(
        "[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1>>[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1-Cl"
    ) == (
        "[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1>>[#6:1]1-[#6:5]=[#6:6]-[#6:4]=[#6:2]-[#6:3]=1-[Cl]"
    )


def test_nitrogen_heterocycle_1():
    assert ValidationManager.remove_ketcher_flavor_smarts(
        "[#7:1;h1]1:[#6:3]:[#6:2]:[#6:4](-[#7]-[#6]):[#6:5]:1>>[#7:1;h1]1:[#6:5]:[#6:4]:[#6:2](-[#7]-[#6]):[#6:3]:1"
    ) == (
        "[nH1:1]1:[#6:3]:[#6:2]:[#6:4](-[#7]-[#6]):[#6:5]:1>>[nH1:1]1:[#6:5]:[#6:4]:[#6:2](-[#7]-[#6]):[#6:3]:1"
    )


def test_nitrogen_heterocycle_2():
    assert ValidationManager.remove_ketcher_flavor_smarts(
        "[#7;h1]1:[#6:3]:[#6:2]:[#6:4](-[#7]-[#6]):[#6:5]:1>>[#7;h1]1:[#6:5]:[#6:4]:[#6:2](-[#7]-[#6]):[#6:3]:1"
    ) == (
        "[nH1]1:[#6:3]:[#6:2]:[#6:4](-[#7]-[#6]):[#6:5]:1>>[nH1]1:[#6:5]:[#6:4]:[#6:2](-[#7]-[#6]):[#6:3]:1"
    )
