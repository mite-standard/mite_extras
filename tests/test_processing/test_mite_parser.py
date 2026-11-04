import json

import pytest

from mite_extras.processing.data_classes import Cofactors, EnzymeAux
from mite_extras.processing.mite_parser import MiteParser


@pytest.fixture
def mite_json():
    with open("tests/test_processing/example_indir_mite/example_valid.json") as infile:
        return json.load(infile)


def test_get_cofactors(mite_json):
    parser = MiteParser()
    log = parser.get_cofactors(mite_json.get("enzyme").get("cofactors"))
    assert isinstance(log, Cofactors)
    assert log.inorganic == ["Fe"]


def test_get_auxenzymes_valid(mite_json):
    parser = MiteParser()
    log = parser.get_auxenzymes(
        auxenzymes=mite_json.get("enzyme").get("auxiliaryEnzymes"),
    )
    assert len(log) == 1
    assert isinstance(log[0], EnzymeAux)


def test_get_reaction_messages(mite_json):
    parser = MiteParser()
    parser.parse_mite_json(mite_json)
    messages = []
    for r in parser.entry.reactions:
        messages.extend(r.warnings)
    assert len(messages) == 0


def test_parse_mite_json_with_taxid(mite_json):
    mite_json["enzyme"]["taxid"] = "1883"
    parser = MiteParser()
    parser.parse_mite_json(mite_json)

    assert parser.entry.enzyme.taxid == "1883"
    assert parser.entry.to_json()["enzyme"]["taxid"] == "1883" 

def test_parse_mite_json_with_taxid_invalid(mite_json):
    with pytest.raises(ValueError):
        mite_json["enzyme"]["taxid"] = 1883
        parser = MiteParser()
        parser.parse_mite_json(mite_json)

def test_parse_mite_json_without_taxid(mite_json):
    mite_json["enzyme"].pop("taxid", None)
    parser = MiteParser()
    parser.parse_mite_json(mite_json)

    assert parser.entry.enzyme.taxid is None
    assert "taxid" not in parser.entry.to_json().get("enzyme", {})