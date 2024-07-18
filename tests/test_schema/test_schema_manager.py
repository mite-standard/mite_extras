import json

import pytest
from mite_extras.schema.schema_manager import SchemaManager


def test_init_valid():
    instance = SchemaManager()
    assert instance.main.exists()
    assert instance.enzyme.exists()
    assert instance.reactions.exists()
    assert instance.version_schema is not None


def test_example_valid():
    instance = SchemaManager()
    with open("tests/test_schema/example_valid.json") as infile:
        example = json.load(infile)
    assert instance.validate_against_schema(example) is None


def test_example_invalid():
    instance = SchemaManager()
    with open("tests/test_schema/example_invalid.json") as infile:
        example = json.load(infile)
    with pytest.raises(ValueError):
        instance.validate_against_schema(example)
