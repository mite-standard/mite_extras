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
