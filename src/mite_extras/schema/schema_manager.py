"""Interface between json schema and calling functionality.

Copyright (c) 2024 to present Mitja Maximilian Zdouc, PhD and individual contributors.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import json
import logging
from pathlib import Path
from typing import Self

import jsonschema
import requests
from pydantic import BaseModel, model_validator
from referencing import Registry, Resource

logger = logging.getLogger("mite_extras")


class SchemaManager(BaseModel):
    """Pydantic-based class to manage validation against json schema

    Attributes:
        main: path to base specs of the MITE json schema
        enzyme: path to specs of the enzyme (meta)data
        reactions: path to specs of reaction data
        changelog_url: changelog-definitions from MIBiG commons
        citation_url: citation-definitions from MIBiG commons
        version_schema: the version of the MITE json schema
    """

    main: Path = Path(__file__).parent.joinpath("entry.json")
    enzyme: Path = Path(__file__).parent.joinpath("definitions/enzyme.json")
    reactions: Path = Path(__file__).parent.joinpath("definitions/reactions.json")
    changelog_url: str = (
        "https://meta.secondarymetabolites.org/schemas/common/changelog.json"
    )
    citation_url: str = (
        "https://meta.secondarymetabolites.org/schemas/common/citation.json"
    )
    version_schema: str | None = None

    @model_validator(mode="after")
    def get_schema_version(self):
        with open(self.main) as infile:
            main = json.load(infile)
        self.version_schema = main.get("$id")
        return self

    def validate_against_schema(self: Self, instance: dict):
        """Validate a dictionary against the MITE JSON schema

        Arguments:
            instance: a dictionary representing a json file

        Raises:
            ValueError: validation of instance against schema led to an error
        """
        logger.debug(
            f"SchemaManager: started validation against MITE schema v{self.version_schema}."
        )

        with open(self.main) as infile:
            main = json.load(infile)
        with open(self.enzyme) as infile:
            enzyme = json.load(infile)
        with open(self.reactions) as infile:
            reactions = json.load(infile)
        changelog = requests.get(self.changelog_url).json()
        citation = requests.get(self.citation_url).json()

        registry = Registry()
        registry = registry.with_resource(
            resource=Resource.from_contents(changelog), uri=self.changelog_url
        )
        registry = registry.with_resource(
            resource=Resource.from_contents(citation), uri=self.citation_url
        )
        registry = registry.with_resource(
            resource=Resource.from_contents(enzyme), uri="definitions/enzyme.json"
        )
        registry = registry.with_resource(
            resource=Resource.from_contents(reactions), uri="definitions/reactions.json"
        )

        try:
            jsonschema.validate(instance=instance, schema=main, registry=registry)
            logger.debug(
                f"SchemaManager: completed validation against MITE schema v{self.version_schema}."
            )
        except jsonschema.exceptions.ValidationError as e:
            raise ValueError(
                f"Validation of instance against MITE schema v{self.version_schema} led to an error: '{e!s}'"
            ) from e
