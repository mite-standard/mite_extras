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

import base64
import logging
from typing import Any, Self

from pydantic import BaseModel, model_validator
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Draw import rdMolDraw2D

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
        databaseIds: an EnzymeDatabaseIds object
        auxiliaryEnzymes: a list of EnzymeAux objects
        references: a list of references
    """

    name: str
    description: str | None = None
    databaseIds: Any
    auxiliaryEnzymes: list | None = None
    references: list

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["name", "description", "references"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val

        json_dict["databaseIds"] = self.databaseIds.to_json()

        if self.auxiliaryEnzymes is not None:
            json_dict["auxiliaryEnzymes"] = [
                entry.to_json() for entry in self.auxiliaryEnzymes
            ]

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        for attr in ["name", "description", "references"]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        html_dict["databaseIds"] = self.databaseIds.to_html()

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
        databaseIds: an EnyzmeDatabaseIds instance
    """

    name: str
    description: str | None = None
    databaseIds: Any

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["name", "description"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val

        json_dict["databaseIds"] = self.databaseIds.to_json()

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        for attr in ["name", "description"]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        html_dict["databaseIds"] = self.databaseIds.to_html()

        return html_dict


class EnyzmeDatabaseIds(BaseModel):
    """Pydantic-based class to represent enzyme-related database ids

    Attributes:
        uniprot: an uniprot ID
        genpept: an NCBI GenPept ID
        mibig: a MIBiG ID
    """

    uniprot: str | None = None
    genpept: str | None = None
    mibig: str | None = None

    @model_validator(mode="after")
    def populate_ids(self):
        # TODO(MMZ 06.08.24): implement the population function (e.g. call from validationmanager)
        return self

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["uniprot", "genpept", "mibig"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val
        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        # TODO(MMZ 06.08.24): implement URLs to bioregistry/mibig
        for attr in ["uniprot", "genpept", "mibig"]:
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
        databaseIds: a ReactionDatabaseIds object
    """

    tailoring: list
    description: str | None = None
    reactionSMARTS: Any
    reactions: list
    evidence: list
    databaseIds: Any | None = None

    def to_json(self: Self) -> dict:
        json_dict = {}

        for attr in ["tailoring", "description"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val

        json_dict["reactionSMARTS"] = self.reactionSMARTS.to_json()
        json_dict["reactions"] = [entry.to_json() for entry in self.reactions]
        json_dict["evidence"] = [entry.to_json() for entry in self.evidence]

        if self.databaseIds is not None:
            json_dict["databaseIds"] = self.databaseIds.to_json()

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}

        for attr in ["tailoring", "description"]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        html_dict["reactionSMARTS"] = self.reactionSMARTS.to_html()
        html_dict["reactions"] = [entry.to_html() for entry in self.reactions]
        html_dict["evidence"] = [entry.to_html() for entry in self.evidence]

        if self.databaseIds is not None:
            html_dict["databaseIds"] = self.databaseIds.to_html()

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
        self.substrate = ValidationManager().cleanup_smiles(self.substrate)
        self.products = [
            ValidationManager().cleanup_smiles(prod) for prod in self.products
        ]
        if self.forbidden_products is not None:
            self.forbidden_products = [
                ValidationManager().cleanup_smiles(prod)
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
        def _smiles_to_svg(smiles: str) -> str:
            """Generates a base64 encoded SVG strin"""
            m = MolFromSmiles(smiles)

            for atom in m.GetAtoms():
                atom.SetAtomMapNum(0)

            m = rdMolDraw2D.PrepareMolForDrawing(m)

            drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
            drawer.DrawMolecule(m)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            return base64.b64encode(svg.encode("utf-8")).decode("utf-8")

        html_dict = {}

        for attr in [
            "isBalanced",
            "isIntermediate",
            "description",
        ]:
            if (val := getattr(self, attr)) is not None:
                html_dict[attr] = val

        html_dict["substrate"] = (self.substrate, _smiles_to_svg(self.substrate))

        html_dict["products"] = [(mol, _smiles_to_svg(mol)) for mol in self.products]

        if self.forbidden_products:
            html_dict["forbidden_products"] = [
                (mol, _smiles_to_svg(mol)) for mol in self.forbidden_products
            ]

        return html_dict


class ReactionDatabaseIds(BaseModel):
    """Pydantic-based class to represent reaction-related database ids

    Attributes:
        rhea: a RHEA ID
        ec: an EC (Enzyme Commission) number
        mite: a MITE ID
    """

    rhea: str | None = None
    ec: str | None = None
    mite: str | None = None

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["rhea", "ec", "mite"]:
            if (val := getattr(self, attr)) is not None:
                json_dict[attr] = val
        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        # TODO(MMZ 06.08.24): implement URLs to bioregistry/mite
        for attr in ["rhea", "ec", "mite"]:
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
