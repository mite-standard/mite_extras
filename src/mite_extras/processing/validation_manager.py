"""Validation functionality.

Copyright (c) 2024 to present Mitja Maximilian Zdouc, PhD and Adriano
Rutz (0000-0003-0443-9902) and individual contributors.

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

from __future__ import annotations

import logging
import re
from itertools import permutations, product
from math import pi

import requests
from pydantic import BaseModel
from rdkit.Chem import (
    CanonSmiles,
    GetMolFrags,
    Mol,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmarts,
    MolToSmiles,
    rdMolEnumerator,
)
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts

logger = logging.getLogger("mite_extras")


class IdValidator(BaseModel):
    """Handles validation and cross-referencing of database identifiers."""

    def cleanup_ids(
        self, genpept: str | None = None, uniprot: str | None = None
    ) -> dict[str, str]:
        """Cleans up IDs using the UniProt SPARQL endpoint.

        Args:
            genpept: An EMBL ID to be converted to UniProt ID.
            uniprot: A UniProt ID to be converted to EMBL ID.

        Returns:
            A dictionary with the original and corresponding ID.

        Raises:
            ValueError: If neither genpept nor uniprot is provided, or if the IDs do not match.
        """
        if not genpept and not uniprot:
            raise ValueError("Please provide one of NCBI Genpept or Uniprot IDs.")

        def fetch_result(query: str) -> str:
            response = requests.get(
                "https://sparql.uniprot.org/sparql",
                params={"query": query, "format": "srj"},
                headers={"Accept": "application/sparql-results+json"},
                timeout=pi,
            )
            if not response.ok:
                raise ValueError(
                    f"HTTP Error while querying Uniprot: {response.status_code}"
                )

            response_json = response.json()
            bindings = response_json.get("results", {}).get("bindings", [])

            if not bindings:
                raise ValueError("No results found in Uniprot response")

            protein_data = bindings[0].get("protein")
            if not protein_data or "value" not in protein_data:
                raise ValueError(
                    "'protein' key or its 'value' is missing in the Uniprot response"
                )

            protein_uri = protein_data["value"]
            return protein_uri.rsplit("/", 1)[-1]

        def build_genpept_query(genpept: str) -> str:
            return f"""
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            PREFIX up: <http://purl.uniprot.org/core/>
            SELECT ?protein
            WHERE {{
                VALUES ?target {{<http://purl.uniprot.org/embl-cds/{genpept}>}} 
                ?protein a up:Protein .
                ?protein rdfs:seeAlso ?target.
            }}
            """

        def build_uniprot_query(uniprot: str) -> str:
            return f"""
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            PREFIX up: <http://purl.uniprot.org/core/>
            SELECT ?protein
            WHERE {{
                VALUES ?prefix {{<http://purl.uniprot.org/database/EMBL>}}
                <http://purl.uniprot.org/uniprot/{uniprot}> rdfs:seeAlso ?protein .
                ?protein up:database ?prefix .
            }}
            """

        try:
            genpept_result = None
            uniprot_result = None

            if genpept:
                genpept_result = fetch_result(build_genpept_query(genpept))
            if uniprot:
                uniprot_result = fetch_result(build_uniprot_query(uniprot))

            if genpept and uniprot and genpept_result != uniprot:
                raise ValueError(
                    f"The provided Genpept ID '{genpept}' and Uniprot ID '{uniprot}' do not correspond to each other!"
                )

            return {
                "genpept": genpept if genpept else uniprot_result,
                "uniprot": uniprot if uniprot else genpept_result,
            }

        except Exception as e:
            raise ValueError(f"Error during enzyme ID validation: {e!s}") from e

    @staticmethod
    def validate_wikidata_qid(qid: str) -> None:
        """Checks if Wikidata QID exists

        Args:
            qid: a valid Wikidata QID
        """

        def build_query(qid: str) -> str:
            return f"""
            ASK {{
                wd:{qid} ?p ?o
            }}
            """

        def fetch_result(query: str) -> str | bool:
            response = requests.get(
                "https://query.wikidata.org/sparql",
                params={"query": query},
                headers={"Accept": "application/sparql-results+json"},
                timeout=pi,
            )
            if not response.ok:
                logger.warning(
                    f"HTTP Error while querying Wikidata: {response.status_code}"
                )
                return False

            data = response.json()
            return data.get("boolean")

        if not fetch_result(query=build_query(qid=qid)):
            logger.warning(f"Wikidata QID '{qid}' does not exist or has no statements.")


class MoleculeValidator(BaseModel):
    """Handles basic molecule validation and canonicalization."""

    @staticmethod
    def _clean_string(string: str) -> str:
        """Remove superfluous backslashes and H's from string."""
        string = string.replace("\\\\", "\\")
        string = re.sub(r";h\d", "", string)
        return string

    def canonicalize_smiles(self, smiles: str) -> str:
        """Canonicalize a SMILES string."""
        mol = MolFromSmiles(self._clean_string(smiles))
        if mol is None:
            raise ValueError(
                f"RDKit rejected SMILES string - is it a valid SMILES?\n" f"{smiles}"
            )
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return CanonSmiles(MolToSmiles(mol))

    def canonicalize_smarts(self, smarts: str) -> str:
        """Canonicalize a SMARTS pattern."""
        mol = MolFromSmarts(self._clean_string(smarts))
        if mol is None:
            raise ValueError(
                f"RDKit rejected SMARTS string - is it a valid pattern?\n" f"{smarts}"
            )
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetAtomMapNum(i)
        return MolToSmarts(MolFromSmiles(CanonSmiles(MolToSmiles(mol))))


class ReactionCleaner(BaseModel):
    """Handles cleaning and formatting of reaction SMARTS."""

    @staticmethod
    def clean_ketcher_format(smarts: str) -> str:
        """Clean Ketcher-specific formatting from reaction SMARTS."""
        replacements = {
            # Halogens with indices
            r"-([FCBI]l?):(\d+)": r"-[\1:\2]",
            # Standalone halogens
            r"-([FCBI]l)(?!\[)": r"-[\1]",
            # Substituent halogens with indices
            r"\(-([FCBI]l?):(\d+)\)": r"(-[\1:\2])",
            # Substituent halogens without indices
            r"\(-([FCBI]l)\)": r"(-[\1])",
            # Nitrogen hydrogens in heterocycles
            r"\[#7:(\d+);h(\d+)\]": r"[nH\2:\1]",
            r"\[#7;h(\d+)\]": r"[nH\1]",
            # erroneous specification of charges in indexed atoms
            r"\[(#\d+):(\d+);([+-])\]": r"[\1;\3:\2]",
        }

        for pattern, replacement in replacements.items():
            smarts = re.sub(pattern, replacement, smarts)
        return smarts

    @staticmethod
    def detect_undesired_smarts(smarts: str):
        """Check for undesired pattern in reaction SMARTS

        negative for absence (must not detect)

        Args:
            smarts: reaction SMARTS string

        Raises:
            ValueError: lack of positive pattern/presence of negative pattern
        """

        negative = {
            # explicit hydrogens
            r"\[H\]": "Explicit hydrogen atoms detected in reaction SMARTS (e.g. '[H]'), which is not allowed.\n Depicting chirality? Please specify the stereochemistry using one of the heavy (non-hydrogen) atoms connected to the stereocenter.\n",
            # CXSMARTS
            r"\|": "Reaction SMARTS with CXSMARTS (Chemaxon SMARTS) elements detected, which are not supported by MITE.\nCXSMARTS (aka 'Extended SMARTS') can be recognized with a suffix starting with a pipe character ('|').\n Please export as a Daylight SMARTS or remove the CXSMARTS suffix manually and try again.\n",
        }

        for key, val in negative.items():
            if re.search(key, smarts):
                raise ValueError(f"{val}")


class ReactionEnumerator(BaseModel):
    """Handles enumeration of reactions and molecules."""

    def enumerate_molecule(self, mol: Mol) -> set[Mol]:
        """Enumerate all possible forms of a molecule."""
        if mol is None:
            return set()
        results = rdMolEnumerator.Enumerate(mol)
        return set(results) if results else {mol}

    def enumerate_smarts(self, smarts: str) -> set[str]:
        """Enumerate SMARTS or skip if intramolecular reaction

        Args:
            smarts: a SMARTS string

        Returns:
            A set of SMARTS strings
        """
        if re.match(r"^\(.+\)|\(.+\)$", smarts):
            return {smarts}

        mol = MolFromSmarts(smarts)
        if mol is None:
            raise ValueError(
                f"RDKit rejected SMARTS string - is it a valid pattern?\n" f"{smarts}"
            )
        enumerated_mols = self.enumerate_molecule(mol)
        return {MolToSmarts(m) for m in enumerated_mols if m is not None}

    def generate_smarts_variants(self, smarts: str) -> list[str]:
        """Generate all possible variants of a SMARTS pattern with comma-separated elements."""
        pattern = re.compile(r"\[([^\]:]+(?:,[^\]:]+)*)\:(\d+)\]")
        matches = list(pattern.finditer(smarts))

        if not matches:
            return [smarts]

        options = []
        for m in matches:
            options_list = [f"[{opt}:{m.group(2)}]" for opt in m.group(1).split(",")]
            options.append(options_list)

        variants = []
        for combination in product(*options):
            new_smarts = smarts
            for m, replacement in zip(matches, combination, strict=True):
                new_smarts = new_smarts.replace(m.group(0), replacement, 1)
            variants.append(new_smarts)

        return variants


class ReactionValidator(BaseModel):
    """Main class for validating chemical reactions."""

    molecule_validator: MoleculeValidator = MoleculeValidator()
    reaction_cleaner: ReactionCleaner = ReactionCleaner()
    enumerator: ReactionEnumerator = ReactionEnumerator()

    def validate_reaction(
        self,
        reaction_smarts: str,
        substrate_smiles: str,
        expected_products: list[str],
        forbidden_products: list[str] | None = None,
        intramolecular: bool = False,
    ) -> None:
        """
        Validate a reaction SMARTS against expected and forbidden products.

        Args:
            reaction_smarts: The reaction SMARTS pattern
            substrate_smiles: The substrate SMILES
            expected_products: List of expected product SMILES
            forbidden_products: Optional list of forbidden product SMILES
            intramolecular: Whether the reaction is intramolecular

        Raises:
            ValueError: missing products or overlap expected & forbidden products or not all expected products generated
        """
        if not expected_products:
            raise ValueError("At least one product must be specified")

        self.reaction_cleaner.detect_undesired_smarts(reaction_smarts)

        reaction_smarts = self.reaction_cleaner.clean_ketcher_format(reaction_smarts)

        substrate_smiles = self.molecule_validator._clean_string(substrate_smiles)

        expected_smiles = {
            self.molecule_validator.canonicalize_smiles(p) for p in expected_products
        }
        forbidden_smiles = {
            self.molecule_validator.canonicalize_smiles(p)
            for p in (forbidden_products or [])
        }

        if overlap := forbidden_smiles & expected_smiles:
            raise ValueError(
                f"Overlap between expected and forbidden products:\n"
                f"{'\n'.join(overlap)}\n"
            )

        # Generate reaction variants and validate
        reactions = self._get_reaction_variants(reaction_smarts)
        predicted_products = self._run_reactions(
            reactions, substrate_smiles, intramolecular
        )

        # Validate predictions
        predicted_smiles = {
            self.molecule_validator.canonicalize_smiles(MolToSmiles(p))
            for p in predicted_products
            if p is not None
        }

        if not expected_smiles.issubset(predicted_smiles):
            raise ValueError(
                f"Reaction did not lead to all expected products.\n"
                f"Expected products:\n"
                f"{'\n'.join(expected_smiles)}\n"
                f"Generated products:\n"
                f"{'\n'.join(predicted_smiles)}\n"
            )

        if overlap := forbidden_smiles & predicted_smiles:
            raise ValueError(
                f"Reaction product(s) belong(s) to the specified forbidden product(s):\n"
                f"{'\n'.join(overlap)}\n"
            )

        logger.debug("Successfully validated reaction SMARTS")

    def _get_reaction_variants(self, reaction_smarts: str) -> set[ChemicalReaction]:
        """Get all possible variants of a reaction."""
        try:
            reactants, products = reaction_smarts.split(">>")
        except ValueError as e:
            raise ValueError("Invalid reaction SMARTS format") from e

        reactant_variants = self.enumerator.generate_smarts_variants(reactants)
        product_variants = self.enumerator.generate_smarts_variants(products)

        reactants_list = [
            self.enumerator.enumerate_smarts(reactant) for reactant in reactant_variants
        ]
        products_list = [
            self.enumerator.enumerate_smarts(product) for product in product_variants
        ]

        reactions = set()
        for r in reactants_list:
            for p in products_list:
                reaction_smarts_combinations = [
                    f"{r_smarts}>>{p_smarts}" for r_smarts in r for p_smarts in p
                ]
                for reaction_smarts_combined in reaction_smarts_combinations:
                    reaction = ReactionFromSmarts(reaction_smarts_combined)
                    if reaction is not None:
                        reactions.add(reaction)

        return reactions

    def _run_reactions(
        self,
        reactions: set[ChemicalReaction],
        substrate_smiles: str,
        intramolecular: bool = False,
    ) -> set[Mol]:
        """Run all reaction variants on the substrate."""
        substrate = MolFromSmiles(substrate_smiles)
        if substrate is None:
            raise ValueError(
                f"RDKit rejected SMILES string - is it a valid SMILES?\n"
                f"{substrate_smiles}"
            )

        substrate_variants = self.enumerator.enumerate_molecule(substrate)
        products = set()

        for mol in substrate_variants:
            reactants = GetMolFrags(mol, asMols=True)

            if intramolecular and len(reactants) > 1:
                continue

            reactant_perms = (
                [reactants] if intramolecular else list(permutations(reactants))
            )

            for reaction in reactions:
                for reactant_combo in reactant_perms:
                    try:
                        reaction_products = reaction.RunReactants(reactant_combo)
                        for product_set in reaction_products:
                            products.update(product_set)
                    except Exception as e:
                        logger.debug(f"Error during running of reaction: {e}")
                        continue

        return {p for p in products if p is not None}
