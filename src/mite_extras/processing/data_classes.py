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
import re
from typing import Any, Self

from pydantic import BaseModel, ValidationError, model_validator
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

from mite_extras.processing.validation_manager import (
    IdValidator,
    MoleculeValidator,
    ReactionValidator,
)

logger = logging.getLogger("mite_extras")


class Entry(BaseModel):
    """Pydantic-based class to represent a MITE entry

    Attributes:
        accession: MITE accession number/identifier.
        status: status of entry
        retirementReasons: a list of retirement reasons
        changelog:  a list of Changelog objects
        enzyme: an Enzyme object
        reactions: a list of Reaction objects
        comment: additional information
    """

    accession: str | None = None
    status: str | None = None
    retirementReasons: list[str] | None = None
    changelog: list | None = None
    enzyme: Any | None = None
    reactions: list[Any] | None = None
    comment: str | None = None

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["accession", "status", "retirementReasons", "comment"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                json_dict[attr] = val

        if self.changelog is not None:
            json_dict["changelog"] = [
                changelog.to_json() for changelog in self.changelog
            ]

        if self.enzyme:
            json_dict["enzyme"] = self.enzyme.to_json()

        if self.reactions:
            json_dict["reactions"] = [reaction.to_json() for reaction in self.reactions]

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        for attr in ["accession", "status", "retirementReasons", "comment"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                html_dict[attr] = val

        if self.changelog is not None:
            html_dict["changelog"] = [
                changelog.to_html() for changelog in self.changelog
            ]

        if self.enzyme:
            html_dict["enzyme"] = self.enzyme.to_html()

        if self.reactions:
            html_dict["reactions"] = [reaction.to_html() for reaction in self.reactions]

        return html_dict


class Changelog(BaseModel):
    """Pydantic-based class to represent changelog

    Attributes:
        version: the MITE version
        contributors: a list of contributors
        reviewers: a list of reviewers
        date: the release date in YYYY-MM-DD
        comment: comment indicating changes
    """

    version: str
    contributors: list
    reviewers: list
    date: str
    comment: str

    @model_validator(mode="after")
    def check_empty(self):
        if self.comment == "":
            raise ValidationError("Changelog: comment must not be an empty string.")
        return self

    def to_json(self: Self) -> dict:
        return {
            attr: getattr(self, attr)
            for attr in ["version", "date", "contributors", "reviewers", "comment"]
        }

    def to_html(self: Self) -> dict:
        return self.to_json()


class Enzyme(BaseModel):
    """Pydantic-based class to represent enzyme information

    Attributes:
        name: the protein name
        description: a brief description of the enzyme
        databaseIds: an EnzymeDatabaseIds object
        auxiliaryEnzymes: a list of EnzymeAux objects
        references: a list of references
        cofactors: optional object containing cofactor lists
    """

    name: str
    description: str | None = None
    databaseIds: Any
    auxiliaryEnzymes: list | None = None
    references: list
    cofactors: Any | None = None

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["name", "description", "references"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                json_dict[attr] = val

        if self.auxiliaryEnzymes:
            json_dict["auxiliaryEnzymes"] = [
                entry.to_json() for entry in self.auxiliaryEnzymes
            ]

        if self.databaseIds.to_json() == {}:
            raise RuntimeError(
                "At least one of 'uniprot' and 'genpept' IDs must be provided."
            )
        else:
            json_dict["databaseIds"] = self.databaseIds.to_json()

        if self.cofactors and self.cofactors.to_json():
            json_dict["cofactors"] = self.cofactors.to_json()

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        for attr in ["name", "description", "references"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                html_dict[attr] = val

        if self.auxiliaryEnzymes:
            html_dict["auxiliaryEnzymes"] = [
                entry.to_html() for entry in self.auxiliaryEnzymes
            ]

        html_dict["databaseIds"] = self.databaseIds.to_html()

        if self.cofactors and self.cofactors.to_html():
            html_dict["cofactors"] = self.cofactors.to_html()

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
            val = getattr(self, attr)
            if val not in (None, ""):
                json_dict[attr] = val

        if self.databaseIds.to_json() != {}:
            json_dict["databaseIds"] = self.databaseIds.to_json()

        return json_dict

    def to_html(self: Self) -> dict:
        html_dict = {}
        for attr in ["name", "description"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                html_dict[attr] = val

        html_dict["databaseIds"] = self.databaseIds.to_html()

        return html_dict


class EnyzmeDatabaseIds(BaseModel):
    """Pydantic-based class to represent enzyme-related database ids

    Class-level instance of IdValidator avoids repeated instantiation

    Attributes:
        uniprot: a UniProt ID
        genpept: an NCBI GenPept ID
        mibig: a MIBiG ID
        wikidata: a wikidata qid reference
    """

    uniprot: str | None = None
    genpept: str | None = None
    mibig: str | None = None
    wikidata: str | None = None

    _id_validator = IdValidator()

    @model_validator(mode="after")
    def validate_ids(self):
        """Cross-reference and populate database IDs using IdValidator.

        Checks if genpept and uniprot IDs correspond
        If only one ID is provided, fetches the other
        Checks if wikidata QID exists

        Raises:
            ValueError: did not provide one of uniprot or genpept IDs
        """
        if not self.genpept and not self.uniprot:
            raise ValueError("Neither 'UniProt' not 'GenPept' IDs provided.")

        try:
            if self.uniprot and self.genpept:
                self._id_validator.cleanup_ids(
                    genpept=self.genpept, uniprot=self.uniprot
                )
            elif self.uniprot:
                data = self._id_validator.cleanup_ids(uniprot=self.uniprot)
                self.genpept = data["genpept"]
            elif self.genpept:
                data = self._id_validator.cleanup_ids(genpept=self.genpept)
                self.uniprot = data["uniprot"]
        except Exception as e:
            logger.warning(f"EnyzmeDatabaseIds: error during ID validation: {e!s}")

        try:
            if self.wikidata:
                self._id_validator.validate_wikidata_qid(self.wikidata)
        except Exception as e:
            logger.error(f"EnyzmeDatabaseIds: error during ID validation: {e!s}")
            raise e from Exception

        return self

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["uniprot", "genpept", "mibig", "wikidata"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                json_dict[attr] = val
        return json_dict

    def to_html(self: Self) -> dict:
        return self.to_json()


class Reaction(BaseModel):
    """Pydantic-based class to represent reactions

    Class-level instance of ReactionValidator avoids repeated instantiation

    Attributes:
        tailoring: a list of tailoring reaction terms
        description: an optional human-readable description
        reactionSMARTS: a ReactionSmarts str
        reactions: a list of ReactionEx objects
        evidence: an evidence object
        databaseIds: a ReactionDatabaseIds object
    """

    tailoring: list
    description: str | None = None
    reactionSMARTS: str
    reactions: list
    evidence: Any
    databaseIds: Any | None = None

    _reaction_validator = ReactionValidator()

    @model_validator(mode="after")
    def cleanup_smarts(self):
        """Clean up reaction SMARTS using the new ReactionCleaner

        Raises:
            ValueError: reaction SMARTS could not be read by RDKit -> invalid
        """
        try:
            cleaned_smarts = (
                self._reaction_validator.reaction_cleaner.clean_ketcher_format(
                    self._reaction_validator.molecule_validator._clean_string(
                        self.reactionSMARTS
                    )
                )
            )

            if ReactionFromSmarts(cleaned_smarts) is None:
                raise ValueError(f"Invalid reaction SMARTS: {self.reactionSMARTS}")

            self.reactionSMARTS = cleaned_smarts
            return self

        except Exception as e:
            raise ValueError(f"Error cleaning reaction SMARTS: {e!s}") from e

    @model_validator(mode="after")
    def validate_reactions(self):
        """Validate all reactions using the new ReactionValidator

        Raises:
            ValueError: Reaction validation failed
        """
        for reaction_id, reaction in enumerate(self.reactions, start=1):
            intramolecular = False
            if re.match(r"^\(.+\)>>|>>\(.+\)$", self.reactionSMARTS):
                intramolecular = True

            try:
                self._reaction_validator.validate_reaction(
                    reaction_smarts=self.reactionSMARTS,
                    substrate_smiles=reaction.substrate,
                    expected_products=reaction.products,
                    forbidden_products=reaction.forbidden_products,
                    intramolecular=intramolecular,
                )
            except Exception as e:
                raise ValueError(
                    f"Reaction SMARTS #{reaction_id} - reaction example #{reaction}: Validation failed for substrate {reaction.substrate}: {e!s}"
                ) from e
        return self

    def to_json(self: Self) -> dict:
        json_dict = {}

        for attr in ["tailoring", "description"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                json_dict[attr] = val

        json_dict["reactionSMARTS"] = self.reactionSMARTS
        json_dict["reactions"] = [entry.to_json() for entry in self.reactions]
        json_dict["evidence"] = self.evidence.to_json()

        if self.databaseIds and self.databaseIds.to_json() != {}:
            json_dict["databaseIds"] = self.databaseIds.to_json()

        return json_dict

    def to_html(self: Self) -> dict:
        def _smarts_to_svg(smarts: str) -> str:
            """Generates a base64 encoded SVG string of the reaction SMARTS"""
            rxn = ReactionFromSmarts(smarts)
            drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
            dopts = drawer.drawOptions()
            dopts.padding = 1e-5
            dopts.clearBackground = False
            drawer.DrawReaction(
                rxn,
                highlightByReactant=True,
                highlightColorsReactants=[(0.69, 0.863, 0.949)],  # RGB blue
            )
            drawer.FinishDrawing()

            svg = drawer.GetDrawingText()
            return base64.b64encode(svg.encode("utf-8")).decode("utf-8")

        html_dict = {}

        for attr in ["tailoring", "description"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                html_dict[attr] = val

        html_dict["reactionSMARTS"] = (
            self.reactionSMARTS,
            _smarts_to_svg(self.reactionSMARTS),
        )
        html_dict["reactions"] = [entry.to_html() for entry in self.reactions]

        evidences = self.evidence.to_html()
        html_dict["evidenceCode"] = evidences["evidenceCode"]
        html_dict["references"] = evidences["references"]

        if self.databaseIds:
            html_dict["databaseIds"] = self.databaseIds.to_html()

        return html_dict


class ReactionEx(BaseModel):
    """Pydantic-based class to represent data of experimentally verified reaction

    Attributes:
        substrate: SMILES string
        products: list of SMILES strings (result from substrate if reactionSMARTS appl)
        forbidden_products: list of SMILES string of forbidden products (must not result)
        isIntermediate: is the reaction an intermediate or a stable product
        description: an optional string
    """

    substrate: str
    products: list
    forbidden_products: list | None = None
    isIntermediate: bool
    description: str | None = None

    # Create validator instance at class level
    _molecule_validator = MoleculeValidator()

    @model_validator(mode="after")
    def validate_smiles(self):
        """Validate and clean SMILES using the new MoleculeValidator"""
        try:
            # Clean substrate SMILES
            self.substrate = self._molecule_validator.canonicalize_smiles(
                self._molecule_validator._clean_string(self.substrate)
            )

            # Clean product SMILES
            self.products = [
                self._molecule_validator.canonicalize_smiles(
                    self._molecule_validator._clean_string(prod)
                )
                for prod in self.products
            ]

            # Clean forbidden product SMILES if present
            if self.forbidden_products:
                cleaned_forbidden = []
                for prod in self.forbidden_products:
                    cleaned = self._molecule_validator.canonicalize_smiles(
                        self._molecule_validator._clean_string(prod)
                    )
                    # Split composite SMILES into individual molecules
                    cleaned_forbidden.extend(cleaned.split("."))
                self.forbidden_products = cleaned_forbidden

            return self

        except Exception as e:
            raise ValueError(f"Error during SMILES cleanup: {e!s}") from e

    def to_json(self: Self) -> dict:
        json_dict = {}

        for attr in [
            "substrate",
            "products",
            "forbidden_products",
            "isIntermediate",
            "description",
        ]:
            val = getattr(self, attr)
            if val not in (None, ""):
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
            dopts = drawer.drawOptions()
            dopts.clearBackground = False

            drawer.DrawMolecule(m)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            return base64.b64encode(svg.encode("utf-8")).decode("utf-8")

        html_dict = {}

        for attr in [
            "isIntermediate",
            "description",
        ]:
            val = getattr(self, attr)
            if val not in (None, ""):
                html_dict[attr] = val

        try:
            html_dict["substrate"] = [
                (mol, _smiles_to_svg(mol)) for mol in self.substrate.split(".")
            ]
        except Exception as e:
            logger.warning(
                f"Error during splitting substrate SMILES for rendering: {e!s}"
            )
            html_dict["substrate"] = [(self.substrate, _smiles_to_svg(self.substrate))]

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
    """

    rhea: str | None = None
    ec: str | None = None

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["rhea", "ec"]:
            val = getattr(self, attr)
            if val not in (None, ""):
                json_dict[attr] = val
        return json_dict

    def to_html(self: Self) -> dict:
        return self.to_json()


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
        return self.to_json()


class Cofactors(BaseModel):
    """Pydantic-based class representing cofactor information

    Attributes:
        organic: a list of organic cofactors
        inorganic: a list of inorganic cofactors
    """

    organic: list | None = None
    inorganic: list | None = None

    def to_json(self: Self) -> dict:
        json_dict = {}
        for attr in ["organic", "inorganic"]:
            val = getattr(self, attr)
            if val not in (None, []):
                json_dict[attr] = val
        return json_dict

    def to_html(self: Self) -> dict:
        return self.to_json()
