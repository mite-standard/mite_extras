import json

import pytest
from mite_extras.processing.data_classes import (
    Changelog,
    ChangelogEntry,
    EnzymeAux,
    Evidence,
    Reaction,
    ReactionEx,
)
from mite_extras.processing.mite_parser import MiteParser
from mite_schema import SchemaManager


@pytest.fixture
def mite_json():
    with open("tests/test_processing/example_indir_mite/example_valid.json") as infile:
        return json.load(infile)


def test_get_changelog_entries_valid(mite_json):
    parser = MiteParser()
    log = parser.get_changelog_entries(
        entries=mite_json.get("changelog").get("releases")[0].get("entries")
    )
    assert len(log) == 2
    assert isinstance(log[0], ChangelogEntry)


def test_get_changelog_valid(mite_json):
    parser = MiteParser()
    log = parser.get_changelog(releases=mite_json.get("changelog").get("releases"))
    assert len(log) == 2
    assert isinstance(log[0], Changelog)


def test_get_databaseids_enzyme_valid(mite_json):
    parser = MiteParser()
    log = parser.get_databaseids_enzyme(data=mite_json.get("enzyme").get("databaseIds"))
    assert log.mibig == "BGC0000026"


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
    assert log.mite == "MITE0000000"


def test_get_reactionex_valid(mite_json):
    parser = MiteParser()
    log = parser.get_reactionex(
        reactions=mite_json.get("reactions")[0].get("reactions")
    )
    assert len(log) == 1
    assert isinstance(log[0], ReactionEx)


def test_get_evidence_valid(mite_json):
    parser = MiteParser()
    log = parser.get_evidence(evidences=mite_json.get("reactions")[0].get("evidence"))
    assert len(log) == 1
    assert isinstance(log[0], Evidence)


def test_get_reactions_valid(mite_json):
    parser = MiteParser()
    log = parser.get_reactions(reactions=mite_json.get("reactions"))
    assert len(log) == 1
    assert isinstance(log[0], Reaction)


def test_parse_raw_json_valid(mite_json):
    parser = MiteParser()
    parser.parse_mite_json(data=mite_json)
    assert SchemaManager().validate_mite(instance=parser.to_json()) is None
