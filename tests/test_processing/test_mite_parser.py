import json

import pytest
from mite_extras.processing.data_classes import Cofactors, EnzymeAux, Reaction
from mite_extras.processing.mite_parser import MiteParser
from mite_schema import SchemaManager


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
        auxenzymes=mite_json.get("enzyme").get("auxiliaryEnzymes")
    )
    assert len(log) == 1
    assert isinstance(log[0], EnzymeAux)


def test_get_databaseids_reaction_valid(mite_json):
    parser = MiteParser()
    log = parser.get_databaseids_reaction(
        data=mite_json.get("reactions")[0].get("databaseIds")
    )
    assert log.ec == "1.2.3.4"


def test_get_reactions_valid(mite_json):
    parser = MiteParser()
    log = parser.get_reactions(reactions=mite_json.get("reactions"))
    assert len(log) == 1
    assert isinstance(log[0], Reaction)


def test_parse_raw_json_valid(mite_json):
    parser = MiteParser()
    parser.parse_mite_json(data=mite_json)
    assert SchemaManager().validate_mite(instance=parser.to_json()) is None
