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

import logging
import re
from itertools import permutations, product
from math import pi
from typing import Self

import requests
from pydantic import BaseModel
from rdkit.Chem import (
    CanonSmiles,
    Descriptors,
    GetMolFrags,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmarts,
    MolToSmiles,
    rdMolEnumerator,
)
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

logger = logging.getLogger("mite_extras")


class ValidationManager(BaseModel):
    """Pydantic-based class to manage validation functions"""

    @staticmethod
    def remove_ketcher_flavor_smarts(string: str) -> str:
        """
        Remove error-inducing ketcher-originating features from (reaction) SMARTS.

        Args:
            string: A (reaction) SMARTS string.

        Returns:
            The modified (reaction) SMARTS string
        """
        # missing square brackets for halogens with indexing
        string = re.sub(r"-Cl:(\d+)", r"-[Cl:\1]", string)
        string = re.sub(r"-F:(\d+)", r"-[F:\1]", string)
        string = re.sub(r"-Br:(\d+)", r"-[Br:\1]", string)
        string = re.sub(r"-I:(\d+)", r"-[I:\1]", string)

        # missing square brackets for halogens w/o indexing
        string = re.sub(r"-Cl", r"-[Cl]", string)
        string = re.sub(r"-F", r"-[F]", string)
        string = re.sub(r"-Br", r"-[Br]", string)
        string = re.sub(r"-I", r"-[I]", string)

        # missing square brackets for substituent halogens with indexing
        string = re.sub(r"\(-Cl:(\d+)\)", r"(-[Cl\1])", string)
        string = re.sub(r"\(-F:(\d+)\)", r"(-[F\1])", string)
        string = re.sub(r"\(-Br:(\d+)\)", r"(-[Br\1])", string)
        string = re.sub(r"\(-I:(\d+)\)", r"(-[I\1])", string)

        # missing square brackets for substituent halogens w/o indexing
        string = re.sub(r"\(-Cl\)", r"(-[Cl])", string)
        string = re.sub(r"\(-F\)", r"(-[F])", string)
        string = re.sub(r"\(-Br\)", r"(-[Br])", string)
        string = re.sub(r"\(-I\)", r"(-[I])", string)

        # erroneous specification of nitrogen hydrogens in heterocycles
        string = re.sub(r"\[#7:(\d+);h(\d)+\]", r"[nH\2:\1]", string)
        string = re.sub(r"\[#7;h(\d)+\]", r"[nH\1]", string)

        # erroneous specification of charges in indexed atoms
        string = re.sub(r"\[(#\d+):(\d+);([+-])\]", r"[\1;\3:\2]", string)

        return string

    @staticmethod
    def has_cx_layer(string: str) -> bool:
        """
        Detects if a given SMILES string is a CXSMILES (Canonical Extended SMILES).

        Args:
            string: A string.

        Returns:
            True if the string has a cx layer, False otherwise.
        """
        if " " in string and "|" in string.split(" ", 1)[1]:
            return True
        return False

    @staticmethod
    def unescape_string(string: str) -> str:
        """Remove superfluous backslashes from string

        Args:
            smiles: a user-submitted string

        Returns:
            A cleaned-up string
        """
        return string.replace("\\\\", "\\")

    @staticmethod
    def remove_hs(string: str) -> str:
        """Remove superfluous H's from string

        Args:
            string: a user-submitted string

        Returns:
            A cleaned-up string
        """
        return re.sub(r";h\d", "", string)

    @staticmethod
    def canonicalize_smiles(smiles: str) -> str:
        """Canonicalizes a SMILES

        Args:
            smiles: a user-submitted SMILES string

        Returns:
            A canonical SMILES string

        Raises:
            ValueError: RDKit could not read SMILES
        """
        mol = MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not read SMILES string '{smiles}'")
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return CanonSmiles(MolToSmiles(mol))

    @staticmethod
    def generate_variants(smarts: str) -> list:
        """
        Generates all possible variants of a SMARTS pattern based on comma-separated elements,
        excluding parts enclosed within '| |' (CX layer).

        Args:
            smarts (str): The original SMARTS pattern containing comma-separated elements.

        Returns:
            list: A list of SMARTS patterns with each possible variant substituted.
        """
        # Regular expression to find parts of the pattern not enclosed within '| |'
        pattern = re.compile(r"(?<!\|)\[([^\]:]+(?:,[^\]:]+)*)\:(\d+)\](?!\|)")
        matches = list(pattern.finditer(smarts))

        if not matches:
            return [smarts]

        options = []
        for m in matches:
            options_list = [f"[{opt}:{m.group(2)}]" for opt in m.group(1).split(",")]
            options.append(options_list)

        all_combinations = product(*options)

        all_variants = []
        for combination in all_combinations:
            new_smarts = smarts
            for m, replacement in zip(matches, combination, strict=False):
                new_smarts = new_smarts.replace(m.group(0), replacement, 1)
            all_variants.append(new_smarts)

        return all_variants

    def canonicalize_smarts(self: Self, smarts: str):
        """Canonicalizes a SMARTS

        Args:
            smarts: a user-submitted SMARTS string

        Returns:
            A canonical SMARTS string

        Raises:
            ValueError: RDKit could not read SMARTS
        """
        mol = MolFromSmarts(smarts)
        if mol is None:
            raise ValueError(f"Could not read SMARTS string '{smarts}'")
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetAtomMapNum(i)
        return MolToSmarts(MolFromSmiles(CanonSmiles(MolToSmiles(mol))))

    def cleanup_reaction_smarts(self: Self, reaction_smarts: str) -> str:
        """Checks a reaction SMARTS

        Args:
            reaction_smarts: a reaction SMARTS string

        Returns:
            The same reaction SMARTS string if valid

        Raises:
            ValueError: RDKit could not read reaction SMARTS
        """
        try:
            reaction_smarts = self.remove_ketcher_flavor_smarts(reaction_smarts)
            reaction_smarts_checked = self.remove_hs(
                self.unescape_string(reaction_smarts)
            )
            reaction = ReactionFromSmarts(reaction_smarts_checked)
            if reaction is None:
                raise ValueError(f"Invalid reaction SMARTS string '{reaction_smarts}'")
            return reaction_smarts_checked
        except Exception as e:
            raise ValueError(f"Error parsing reaction SMARTS: {e!s}") from None

    def cleanup_smarts(self: Self, smarts: str) -> str:
        """Cleans up an input SMARTS string

        Args:
            smarts: a SMARTS string

        Returns:
            The SMARTS in RDKit-canonical format

        Raises:
            ValueError: RDKit could not read SMARTS
        """
        unhed = self.remove_hs(self.unescape_string(smarts))
        if self.has_cx_layer(smarts):
            return unhed
        else:
            return self.canonicalize_smarts(unhed)

    def cleanup_smiles(self: Self, smiles: str) -> str:
        """Cleans up an input SMILES string

        Args:
            smiles: a SMILES string

        Returns:
            The SMILES in RDKit-canonical format

        Raises:
            ValueError: RDKit could not read SMILES
        """
        unhed = self.remove_hs(self.unescape_string(smiles))
        if self.has_cx_layer(smiles):
            return unhed
        else:
            return self.canonicalize_smiles(unhed)

    def split_smiles(self: Self, smiles: str) -> list:
        """Split composite SMILES into a list

        Args:
            smiles: a SMILES string

        Returns:
            A list where composite SMILES have been split
        """
        return smiles.split(".")

    def enumerate(self: Self, mol) -> set:
        mols = set()
        if mol is not None:
            res = rdMolEnumerator.Enumerate(mol)
            if len(res) != 0:
                for m in rdMolEnumerator.Enumerate(mol):
                    mols.add(m)
            else:
                mols.add(mol)
        return mols

    def enumerate_reaction_smarts(self: Self, reaction_smarts: str) -> set:
        """Enumerates a reaction SMARTS string

        Args:
            reaction_smarts: a reaction SMARTS string

        Returns:
            A set of enumerated reaction SMARTS

        Raises:
            ValueError: RDKit could not read reaction SMARTS
        """
        try:
            reactants_smarts, products_smarts = reaction_smarts.split(">>")
        except ValueError as e:
            raise ValueError(
                "Invalid reaction SMARTS format. Ensure it contains '>>' separating reactants and products."
            ) from e

        # Generate all variants of reactants and products
        reactants_variants = self.generate_variants(reactants_smarts)
        products_variants = self.generate_variants(products_smarts)

        # Set to hold all possible enumerated reaction SMARTS
        enumerated_reactions = set()

        # Enumerate the reactants and products
        for r_variant in reactants_variants:
            reactant = MolFromSmarts(r_variant)
            reactants_enumerated = self.enumerate(reactant)
            enumerated_reactants_smarts = {
                MolToSmarts(r_mol) for r_mol in reactants_enumerated
            }

            for p_variant in products_variants:
                product = MolFromSmarts(p_variant)
                products_enumerated = self.enumerate(product)
                enumerated_products_smarts = {
                    MolToSmarts(p_mol) for p_mol in products_enumerated
                }

                # Combine all enumerated reactants and products
                for r_smarts in enumerated_reactants_smarts:
                    for p_smarts in enumerated_products_smarts:
                        enumerated_reactions.add(f"{r_smarts}>>{p_smarts}")
        return enumerated_reactions

    @staticmethod
    def cleanup_ids(
        genpept: str | None = None, uniprot: str | None = None
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
            raise ValueError("Please provide one of 'genpept' or 'uniprot'.")

        def fetch_result(query: str) -> str:
            response = requests.get(
                "https://sparql.uniprot.org/sparql",
                params={"query": query, "format": "srj"},
                headers={"Accept": "application/sparql-results+json"},
                timeout=pi,
            )
            if not response.ok:
                raise ValueError(f"HTTP Error: {response.status_code}")
            response_json = response.json()
            bindings = response_json.get("results", {}).get("bindings", [])
            if not bindings:
                raise ValueError("No results found in the response")
            protein_data = bindings[0].get("protein")
            if not protein_data or "value" not in protein_data:
                raise ValueError(
                    "'protein' key or its 'value' is missing in the response"
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

        if genpept:
            genpept_result = fetch_result(build_genpept_query(genpept))
        if uniprot:
            uniprot_result = fetch_result(build_uniprot_query(uniprot))

        if genpept and uniprot and genpept_result != uniprot:
            raise ValueError(
                f"The provided genpept ID '{genpept}' and uniprot ID '{uniprot}' do not match"
            )

        return {
            "genpept": genpept if genpept else uniprot_result,
            "uniprot": uniprot if uniprot else genpept_result,
        }

    def validate_reaction_smarts(
        self: Self,
        reaction_smarts: str,
        substrate_smiles: str,
        expected_products: list[str],
        forbidden_products: list[str],
    ) -> None:
        """Validates the reaction SMARTS

        Args:
            reaction_smarts: a reaction SMARTS string
            substrate_smiles: a SMILES string representing the substrate
            expected_products: a list of expected product SMILES strings
            forbidden_products: a list of forbidden product SMILES strings

        Raises:
            ValueError: If the reaction does not meet the expectations
        """
        # Check that expected products are not empty
        if not expected_products:
            raise ValueError("Expected products list cannot be empty.")

        expected_smiles_set = set()
        for product in expected_products:
            expected_smiles_set.add(self.cleanup_smiles(product))

        # TODO add test
        if not expected_smiles_set:
            raise ValueError("Expected products list is invalid.")

        # Check for forbidden products in expected products
        forbidden_smiles_set = set()
        if forbidden_products:
            for product in forbidden_products:
                forbidden_smiles_set.add(self.cleanup_smiles(product))

        if forbidden_smiles_set.intersection(expected_smiles_set):
            raise ValueError("Some expected products are listed as forbidden products.")

        reactions = self.enumerate_reaction_smarts(
            self.cleanup_reaction_smarts(reaction_smarts)
        )

        expected_mols = set()
        for smiles in expected_smiles_set:
            for m in self.enumerate(MolFromSmiles(self.cleanup_smiles(smiles))):
                expected_mols.add(m)

        forbidden_mols = set()
        for smiles in forbidden_smiles_set:
            for m in self.enumerate(MolFromSmiles(self.cleanup_smiles(smiles))):
                forbidden_mols.add(m)

        if None in expected_mols or None in forbidden_mols:
            raise ValueError(
                "One or more expected/forbidden products could not be parsed."
            )

        reactant_mol = MolFromSmiles(self.cleanup_smiles(substrate_smiles))

        if reactant_mol is None:
            raise ValueError(f"Invalid substrate SMILES '{substrate_smiles}'")

        reactant_mol_enumerated = self.enumerate(reactant_mol)

        predicted_products = set()
        for mol in reactant_mol_enumerated:
            reactants = GetMolFrags(mol, asMols=True)
            # Important to allow reactants to be in whatever order
            reactants_permutations = list(permutations(reactants))
            for reaction in reactions:
                reaction_instance = ReactionFromSmarts(reaction)
                for reactant_combination in reactants_permutations:
                    products = reaction_instance.RunReactants(reactant_combination)
                    for product_tuple in products:
                        for product in product_tuple:
                            predicted_products.add(product)

        predicted_smiles_set = set()
        for product in predicted_products:
            predicted_smiles_set.add(self.cleanup_smiles(MolToSmiles(product)))

        predicted_mols = set()
        for smiles in predicted_smiles_set:
            predicted_mols.update(
                self.enumerate(MolFromSmiles(self.cleanup_smiles(smiles)))
            )

        if None in predicted_mols:
            raise ValueError("One or more predicted products could not be parsed.")

        if not expected_smiles_set.issubset(predicted_smiles_set):
            raise ValueError(
                f"Products '{predicted_smiles_set}' do not meet expectations '{expected_smiles_set}'."
            )

        # Check if any forbidden products are present in the predicted products
        if forbidden_smiles_set.intersection(predicted_smiles_set):
            raise ValueError("Forbidden products were found in the reaction output.")

        logger.debug("ValidationManager: successfully validated reaction SMARTS")
