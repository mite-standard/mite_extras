from builtins import AssertionError

import pytest
from mite_extras.processing.data_classes import Entry
from pydantic import ValidationError


def test_entry_init_valid():
    assert isinstance(Entry(), Entry)


def test_entry_init_invalid():
    with pytest.raises(ValidationError):
        isinstance(Entry(accession=12), Entry)


def test_entry_to_json_valid():
    data = {
        "accession": "MITE1234567",
        "quality": "high",
        "status": "pending",
        "retirementReasons": ["lorem ipsum", "dolor sit amet"],
        "comment": "No comment",
        "attachments": {"attachment1": "message1"},
    }
    instance = Entry(**data)
    json_dict = instance.to_json()
    assert json_dict == data
    # TODO(MMZ 17.07.24): expand data for missing keys
