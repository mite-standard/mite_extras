from builtins import AssertionError

import pytest
from mite_extras.processing.data_classes import Entry, Evidence
from pydantic import ValidationError


def test_entry_init_valid():
    assert isinstance(Entry(), Entry)


def test_entry_init_invalid():
    with pytest.raises(ValidationError):
        isinstance(Entry(accession=12), Entry)


def test_evidence_init_valid():
    assert isinstance(Evidence(evidenceCode=["lorem"], references=["ipsum"]), Evidence)


def test_evidence_init_invalid():
    with pytest.raises(ValidationError):
        isinstance(Evidence(evidenceCode=12, references=12), Evidence)


def test_entry_to_json_valid():
    # TODO: rework data
    data = {
        "accession": "MITE1234567",
        "quality": "high",
        "status": "pending",
        "retirementReasons": ["lorem ipsum", "dolor sit amet"],
        "comment": "No comment",
        "attachments": {"attachment1": "message1"},
    }
    instance = Entry(**data)
    assert instance.to_json() is not None
    # TODO(MMZ 17.07.24): expand data for missing keys


def test_evidence_to_json_valid():
    instance = Evidence(evidenceCode=["lorem"], references=["ipsum"])
    assert instance.to_json() is not None
