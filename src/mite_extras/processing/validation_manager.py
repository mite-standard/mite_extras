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
from typing import (
    Optional,
    Self,
)

import requests
from pydantic import BaseModel
from rdkit.Chem import (
    CanonSmiles,
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
    def unescape_smiles(smiles: str) -> str:
        """Remove superfluous backslashes from SMILES string

        Args:
            smiles: a user-submitted SMILES string

        Returns:
            A cleaned-up SMILES string
        """
        return smiles.replace("\\\\", "\\")

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
        return CanonSmiles(MolToSmiles(mol))

    @staticmethod
    def canonicalize_smarts(smarts: str) -> str:
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
        return MolToSmarts(MolFromSmiles(CanonSmiles(MolToSmiles(mol))))

    @staticmethod
    def check_reaction_smarts(reaction_smarts: str) -> str:
        """Checks a reaction SMARTS

        Args:
            reaction_smarts: a reaction SMARTS string

        Returns:
            The same reaction SMARTS string if valid

        Raises:
            ValueError: RDKit could not read reaction SMARTS
        """
        try:
            reaction = ReactionFromSmarts(reaction_smarts)
            if reaction is None:
                raise ValueError(f"Invalid reaction SMARTS string '{reaction_smarts}'")
            return reaction_smarts
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
        return self.canonicalize_smarts(smarts)

    def cleanup_smiles(self: Self, smiles: str) -> str:
        """Cleans up an input SMILES string

        Args:
            smiles: a SMILES string

        Returns:
            The SMILES in RDKit-canonical format

        Raises:
            ValueError: RDKit could not read SMILES
        """
        return self.canonicalize_smiles(self.unescape_smiles(smiles))

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
        # Check for forbidden products in expected products
        forbidden_smiles_set = set()
        if forbidden_products:
            forbidden_smiles_set = set(forbidden_smiles_set)

        expected_smiles_set = set(expected_products)
        if forbidden_smiles_set.intersection(expected_smiles_set):
            raise ValueError("Some expected products are listed as forbidden products.")

        # Check that expected products are not empty
        if not expected_products:
            raise ValueError("Expected products list cannot be empty.")

        # Attempt to parse the reaction SMARTS
        reaction = ReactionFromSmarts(reaction_smarts)
        if reaction is None:
            raise ValueError(f"Invalid reaction SMARTS '{reaction_smarts}'")

        # Convert expected and forbidden products to RDKit molecule objects
        expected_mols = {
            MolFromSmiles(self.unescape_smiles(smiles)) for smiles in expected_products
        }
        forbidden_mols = {
            MolFromSmiles(self.unescape_smiles(smiles))
            for smiles in forbidden_smiles_set
        }

        if None in expected_mols or None in forbidden_mols:
            raise ValueError(
                "One or more expected/forbidden products could not be parsed."
            )

        # Convert substrate to RDKit molecule object
        substrate_mol = MolFromSmiles(self.unescape_smiles(substrate_smiles))
        if substrate_mol is None:
            raise ValueError(f"Invalid substrate SMILES '{substrate_smiles}'")

        # Enumerate possible conformations of the substrate molecule
        substrate_mol_enumerated = rdMolEnumerator.Enumerate(substrate_mol)
        if not substrate_mol_enumerated:
            substrate_mol_enumerated = [substrate_mol]

        # Generate products from the reaction and substrate
        predicted_products = []
        for substrate in substrate_mol_enumerated:
            products = reaction.RunReactants([substrate])
            predicted_products.extend(products)

        predicted_smiles = set()
        for product_set in predicted_products:
            for product in product_set:
                predicted_smiles.add(CanonSmiles(MolToSmiles(product)))

        predicted_mols = {MolFromSmiles(smiles) for smiles in predicted_smiles}

        if None in predicted_mols:
            raise ValueError("One or more predicted products could not be parsed.")

        # Check products meet expectations
        if not expected_smiles_set.issubset(predicted_smiles):
            raise ValueError(
                f"Products '{predicted_smiles}' do not meet expectations '{expected_smiles_set}'."
            )

        # Check if any forbidden products are present in the predicted products
        if forbidden_smiles_set.intersection(predicted_smiles):
            raise ValueError("Forbidden products were found in the reaction output.")

        logger.debug("ValidationManager: successfully validated reaction SMARTS")
