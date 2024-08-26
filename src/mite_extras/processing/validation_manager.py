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
from itertools import permutations
from typing import (
    Optional,
    Self,
)

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
            reaction_smarts_checked = self.remove_hs(self.unescape_string(reaction_smarts))
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
            A list of reaction SMARTS

        Raises:
            ValueError: RDKit could not read reaction SMARTS
        """
        try:
            reactants_smarts, products_smarts = reaction_smarts.split(">>")
        except ValueError:
            raise ValueError("Invalid reaction SMARTS format. Ensure it contains '>>' separating reactants and products.")
        reactant = MolFromSmarts(reactants_smarts)
        product = MolFromSmarts(products_smarts)

        # Enumerate the reactants and products
        reactants_enumerated = self.enumerate(reactant)
        products_enumerated = self.enumerate(product)

        # Generate all possible combinations of enumerated reactants and products
        enumerated_reactions = set()
        enumerated_reactants_smarts = set()
        enumerated_products_smarts = set()
        for r_mol in reactants_enumerated:
            enumerated_reactants_smarts.add(MolToSmarts(r_mol))
        for p_mol in products_enumerated:
            enumerated_products_smarts.add(MolToSmarts(p_mol))
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
            )
            if not response.ok:
                raise ValueError("Failed to fetch data from UniProt SPARQL endpoint")
            results = response.json().get("results", {}).get("bindings", [])
            if not results:
                raise ValueError("No matching results found")
            return results[0]["protein"]["value"].rsplit("/", 1)[-1]

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

        reactions = self.enumerate_reaction_smarts(self.cleanup_reaction_smarts(reaction_smarts))

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
            predicted_mols.update(self.enumerate(MolFromSmiles(self.cleanup_smiles(smiles))))

        if None in predicted_mols:
            raise ValueError("One or more predicted products could not be parsed.")

        # Check products meet expectations
        # TODO add tests
        # TODO @MMZ tell me if you want to externalize it
        # See MITE0000095 for an example
        for expected_smiles in expected_smiles_set:
            expected_mol = MolFromSmiles(expected_smiles)
            if not any(MolFromSmiles(predicted_smiles).HasSubstructMatch(expected_mol)
                       for predicted_smiles in predicted_smiles_set):
                raise ValueError(
                    f"Products '{predicted_smiles_set}' do not meet expectations '{expected_smiles_set}'."
                )
        # if not expected_smiles_set.issubset(predicted_smiles_set):
        #     raise ValueError(
        #         f"Products '{predicted_smiles_set}' do not meet expectations '{expected_smiles_set}'."
        #     )

        # Check if any forbidden products are present in the predicted products
        if forbidden_smiles_set.intersection(predicted_smiles_set):
            raise ValueError("Forbidden products were found in the reaction output.")

        logger.debug("ValidationManager: successfully validated reaction SMARTS")
