import json

import jsonschema
from referencing import Registry, Resource


def test_schema_valid():
    with open("src/mite_extras/schema/entry.json") as infile:
        main_schema = json.load(infile)
    with open("src/mite_extras/schema/commons/changelog.json") as infile:
        changelog = json.load(infile)
    with open("src/mite_extras/schema/commons/citation.json") as infile:
        citation = json.load(infile)
    with open("src/mite_extras/schema/definitions/enzyme.json") as infile:
        enzyme = json.load(infile)
    with open("src/mite_extras/schema/definitions/reactions.json") as infile:
        reactions = json.load(infile)

    registry = Registry()
    registry = registry.with_resource(
        resource=Resource.from_contents(changelog), uri="commons/changelog.json"
    )
    registry = registry.with_resource(
        resource=Resource.from_contents(citation), uri="commons/citation.json"
    )
    registry = registry.with_resource(
        resource=Resource.from_contents(enzyme), uri="definitions/enzyme.json"
    )
    registry = registry.with_resource(
        resource=Resource.from_contents(reactions), uri="definitions/reactions.json"
    )

    with open("tests/test_schema/example.json") as infile:
        example = json.load(infile)

    assert (
        jsonschema.validate(instance=example, schema=main_schema, registry=registry)
        is None
    )
