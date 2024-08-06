import json

import pytest
from mite_extras.processing.raw_parser import RawParser


@pytest.fixture
def raw_json():
    with open("tests/test_processing/example_indir_raw/example_valid.json") as infile:
        return json.load(infile)


def test_parse_raw_json_valid(raw_json):
    parser = RawParser()
    parser.parse_raw_json(name="newfile", input_data=raw_json)
    out_dict = parser.to_json()
    assert isinstance(out_dict, dict)
