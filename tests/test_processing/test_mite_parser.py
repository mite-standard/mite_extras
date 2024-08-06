import json

import pytest
from mite_extras.processing.data_classes import (
    Changelog,
    ChangelogEntry,
    Entry,
    EnzymeAux,
    Evidence,
    Reaction,
    ReactionEx,
)
from mite_extras.processing.mite_parser import MiteParser


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


def test_get_auxenzymes_valid(mite_json):
    parser = MiteParser()
    log = parser.get_auxenzymes(
        auxenzymes=mite_json.get("enzyme").get("auxiliaryEnzymes")
    )
    assert len(log) == 1
    assert isinstance(log[0], EnzymeAux)


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
    assert isinstance(parser.entry, Entry)
