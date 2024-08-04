"""Validation and internal representation of MITE schema.

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

import logging
from typing import Any, Self

from pydantic import BaseModel, model_validator

from mite_extras.processing.validation_manager import ValidationManager

logger = logging.getLogger("mite_extras")

# TODO(MMZ 18.07.24): add validation functions like for SMILES in ReactionEx class @adafede


class Entry(BaseModel):
    """Pydantic-based class to represent a MITE entry

    Attributes:
        accession: MITE accession number/identifier.
        quality: quality of entry
        status: status of entry
        retirementReasons: a list of retirement reasons
        changelog:  a list of Changelog objects
        enzyme: an Enzyme object
        reactions: a list of Reaction objects
        comment: additional information
        attachments: a dict for further information
    """

    accession: str | None = None
    quality: str | None = None
    status: str | None = None
    retirementReasons: list[str] | None = None
    changelog: list[Any] | None = None
    enzyme: Any | None = None
    reactions: list[Any] | None = None
    comment: str | None = None
    attachments: dict | None = None

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in [
            "accession",
            "quality",
            "status",
            "retirementReasons",
            "comment",
            "attachments",
        ]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val

        if self.changelog is not None:
            json_dict["changelog"] = {
                "releases": [changelog.to_json() for changelog in self.changelog]
            }

        if self.enzyme is not None:
            json_dict["enzyme"] = self.enzyme.to_json()

        if self.reactions is not None:
            json_dict["reactions"] = [reaction.to_json() for reaction in self.reactions]

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        for attr in [
            "accession",
            "quality",
            "status",
            "retirementReasons",
            "comment",
            "attachments",
        ]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        if self.changelog is not None:
            html_dict["changelog"] = [
                changelog.to_html() for changelog in self.changelog
            ]

        if self.enzyme is not None:
            html_dict["enzyme"] = self.enzyme.to_html()

        if self.reactions is not None:
            html_dict["reactions"] = [reaction.to_html() for reaction in self.reactions]

        return html_dict


class Changelog(BaseModel):
    """Pydantic-based class to represent changelog

    Attributes:
        version: the MITE version
        date: the release date in YYYY-MM-DD
        entries: a list of ChangelogEntry objects
    """

    version: str
    date: str
    entries: list

    def to_json(self: Self) -> dict:
        return {
            "version": self.version,
            "date": self.date,
            "entries": [entry.to_json() for entry in self.entries],
        }

    def to_html(self: Self) -> dict:
        return {
            "version": self.version,
            "date": self.date,
            "entries": [entry.to_html() for entry in self.entries],
        }


class ChangelogEntry(BaseModel):
    """Pydantic-based class to represent changelog entries

    Attributes:
        contributors: a list of contributors
        reviewers: a list of reviewers
        date: the date (YYYY-MM-DD) of the last edit (e.g. review)
        comment: comment indicating changes
    """

    contributors: list
    reviewers: list
    date: str
    comment: str

    def to_json(self: Self) -> dict:
        return {
            "contributors": self.contributors,
            "reviewers": self.reviewers,
            "date": self.date,
            "comment": self.comment,
        }

    def to_html(self: Self) -> dict:
        return {
            "contributors": self.contributors,
            "reviewers": self.reviewers,
            "date": self.date,
            "comment": self.comment,
        }


class Enzyme(BaseModel):
    """Pydantic-based class to represent enzyme information

    Attributes:
        name: the protein name
        description: an optional description
        databaseIds: a list of database IDs
        auxiliaryEnzymes: a list of EnzymeAux objects
        references: a list of references
    """

    name: str
    description: str | None = None
    databaseIds: list
    auxiliaryEnzymes: list | None = None
    references: list

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["name", "description", "databaseIds", "references"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val

        if self.auxiliaryEnzymes is not None:
            json_dict["auxiliaryEnzymes"] = [
                entry.to_json() for entry in self.auxiliaryEnzymes
            ]

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        for attr in ["name", "description", "databaseIds", "references"]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        if self.auxiliaryEnzymes is not None:
            html_dict["auxiliaryEnzymes"] = [
                entry.to_html() for entry in self.auxiliaryEnzymes
            ]

        return html_dict


class EnzymeAux(BaseModel):
    """Pydantic-based class to represent auxiliary enzyme information

    Attributes:
        name: name of the auxiliary enzyme
        description: an optional description
        databaseIds: a list of database IDs
    """

    name: str
    description: str | None = None
    databaseIds: list

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["name", "description", "databaseIds"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val
        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        for attr in ["name", "description", "databaseIds"]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val
        return html_dict


class Reaction(BaseModel):
    """Pydantic-based class to represent reactions

    Attributes:
        tailoring: a list of tailoring reaction terms
        description: an optional human-readable description
        reactionSMARTS: a ReactionSmarts object
        reactions: a list of ReactionEx objects
        evidence: list of Evidence objects
        databaseIds: a list of database IDs
    """

    tailoring: list
    description: str | None = None
    reactionSMARTS: Any
    reactions: list
    evidence: list
    databaseIds: list | None = None

    def to_json(self: Self) -> dict:
        json_dict = {}

        for attr in ["tailoring", "description", "databaseIds"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val

        json_dict["reactionSMARTS"] = self.reactionSMARTS.to_json()
        json_dict["reactions"] = [entry.to_json() for entry in self.reactions]
        json_dict["evidence"] = [entry.to_json() for entry in self.evidence]

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}

        for attr in ["tailoring", "description", "databaseIds"]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        html_dict["reactionSMARTS"] = self.reactionSMARTS.to_html()
        html_dict["reactions"] = [entry.to_html() for entry in self.reactions]
        html_dict["evidence"] = [entry.to_html() for entry in self.evidence]

        return html_dict


class ReactionSmarts(BaseModel):
    """Pydantic-based class to represent reactions

    Attributes:
        reactionSMARTS: a reaction SMARTS
        isIterative: flag indicating if reaction SMARTS is applied more than once
    """

    reactionSMARTS: str
    isIterative: bool

    def to_json(self: Self) -> dict:
        json_dict = {}

        for attr in ["reactionSMARTS", "isIterative"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}

        for attr in ["reactionSMARTS", "isIterative"]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        return html_dict


class ReactionEx(BaseModel):
    """Pydantic-based class to represent data of experimentally verified reaction

    Attributes:
        substrate: SMILES string
        products: list of SMILES strings (result from substrate if reactionSMARTS appl)
        forbidden_products: list of SMILES string of forbidden products (must not result)
        isBalanced: is reaction balanced
        isIntermediate: is the reaction an intermediate or a stable product
        description: an optional string
    """

    substrate: str
    products: list
    forbidden_products: list | None = None
    isBalanced: bool
    isIntermediate: bool
    description: str | None = None

    @model_validator(mode="after")
    def validate_smiles(self):
        self.substrate = ValidationManager().canonicalize_smiles(self.substrate)
        self.products = [
            ValidationManager().canonicalize_smiles(prod) for prod in self.products
        ]
        if self.forbidden_products is not None:
            self.forbidden_products = [
                ValidationManager().canonicalize_smiles(prod)
                for prod in self.forbidden_products
            ]

        return self

    def to_json(self: Self) -> dict:
        json_dict = {}

        for attr in [
            "substrate",
            "products",
            "forbidden_products",
            "isBalanced",
            "isIntermediate",
            "description",
        ]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}

        for attr in [
            "substrate",
            "products",
            "forbidden_products",
            "isBalanced",
            "isIntermediate",
            "description",
        ]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        return html_dict


class Evidence(BaseModel):
    """Pydantic-based class to represent evidence information

    Attributes:
        evidenceCode: a list of evidence code strings
        references: a list of references
    """

    evidenceCode: list
    references: list

    def to_json(self: Self) -> dict:
        return {"evidenceCode": self.evidenceCode, "references": self.references}

    def to_html(self: Self) -> dict:
        return {"evidenceCode": self.evidenceCode, "references": self.references}
