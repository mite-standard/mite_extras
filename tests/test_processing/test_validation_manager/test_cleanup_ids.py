from unittest.mock import patch

import pytest
from mite_extras.processing.validation_manager import ValidationManager


def mock_requests_get(*args, **kwargs):
    class MockResponse:
        def __init__(self, json_data, status_code):
            self.json_data = json_data
            self.status_code = status_code
            self.ok = self.status_code == 200

        def json(self):
            return self.json_data

    if "embl-cds/AAM70353.1" in kwargs["params"]["query"]:
        return MockResponse(
            {
                "results": {
                    "bindings": [
                        {"protein": {"value": "http://purl.uniprot.org/uniprot/Q8KND5"}}
                    ]
                }
            },
            200,
        )
    elif "embl-cds/CAA71118.1" in kwargs["params"]["query"]:
        return MockResponse(
            {
                "results": {
                    "bindings": [
                        {"protein": {"value": "http://purl.uniprot.org/uniprot/Q8KND4"}}
                    ]
                }
            },
            200,
        )
    elif "uniprot/Q8KND5" in kwargs["params"]["query"]:
        return MockResponse(
            {
                "results": {
                    "bindings": [
                        {
                            "protein": {
                                "value": "http://purl.uniprot.org/embl-cds/AAM70353.1"
                            }
                        }
                    ]
                }
            },
            200,
        )
    elif "uniprot/Q8KND4" in kwargs["params"]["query"]:
        return MockResponse(
            {
                "results": {
                    "bindings": [
                        {
                            "protein": {
                                "value": "http://purl.uniprot.org/embl-cds/AAM70354.1"
                            }
                        }
                    ]
                }
            },
            200,
        )
    elif "uniprot/WrongID" in kwargs["params"]["query"]:
        return MockResponse({"results": {"bindings": []}}, 200)

    return MockResponse(None, 404)


@patch("requests.get", side_effect=mock_requests_get)
def test_cleanup_ids_genpept(mock_get):
    result = ValidationManager().cleanup_ids(genpept="AAM70353.1")
    assert result == {"genpept": "AAM70353.1", "uniprot": "Q8KND5"}


@patch("requests.get", side_effect=mock_requests_get)
def test_cleanup_ids_uniprot(mock_get):
    result = ValidationManager().cleanup_ids(uniprot="Q8KND5")
    assert result == {"genpept": "AAM70353.1", "uniprot": "Q8KND5"}


@patch("requests.get", side_effect=mock_requests_get)
def test_cleanup_ids_both_ok(mock_get):
    result = ValidationManager().cleanup_ids(genpept="AAM70353.1", uniprot="Q8KND5")
    assert result == {"genpept": "AAM70353.1", "uniprot": "Q8KND5"}


@patch("requests.get", side_effect=mock_requests_get)
def test_cleanup_ids_both_fail_1(mock_get):
    with pytest.raises(
        ValueError,
        match="The provided genpept ID 'CAA71118.1' and uniprot ID 'Q8KND5' do not match",
    ):
        ValidationManager().cleanup_ids(genpept="CAA71118.1", uniprot="Q8KND5")


@patch("requests.get", side_effect=mock_requests_get)
def test_cleanup_ids_both_fail_2(mock_get):
    with pytest.raises(
        ValueError,
        match="The provided genpept ID 'AAM70353.1' and uniprot ID 'Q8KND4' do not match",
    ):
        ValidationManager().cleanup_ids(genpept="AAM70353.1", uniprot="Q8KND4")


@patch("requests.get", side_effect=mock_requests_get)
def test_cleanup_ids_none(mock_get):
    with pytest.raises(
        ValueError, match="Please provide one of 'genpept' or 'uniprot'"
    ):
        ValidationManager().cleanup_ids()


@patch("requests.get", side_effect=mock_requests_get)
def test_cleanup_ids_invalid(mock_get):
    with pytest.raises(
        ValueError, match="Failed to fetch data from UniProt SPARQL endpoint"
    ):
        ValidationManager().cleanup_ids(genpept="invalidID")
