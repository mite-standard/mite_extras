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
    AllChem,
    CanonSmiles,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmarts,
    MolToSmiles,
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

    # TODO (AR 2024-08-06): implement tests

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
            reaction = AllChem.ReactionFromSmarts(reaction_smarts)
            if reaction is None:
                raise ValueError(f"Invalid reaction SMARTS string '{reaction_smarts}'")
            return reaction_smarts
        except Exception as e:
            raise ValueError(f"Error parsing reaction SMARTS: {e!s}") from None

    # TODO (AR 2024-08-06): implement tests

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
        """Cleans up IDs using the UniProt SPARQL endpoint

        Args:
            genpept: an EMBL ID to be converted to UniProt ID
            uniprot: a UniProt ID to be converted to EMBL ID

        Returns:
            A dictionary with the original and corresponding ID

        Raises:
            ValueError: if neither embl_id nor uniprot_id is provided, or if the IDs do not match
        """
        if not genpept and not uniprot:
            raise ValueError("Please provide one of 'genpept' or 'uniprot'.")

        genpept_query = (
            f"""
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        SELECT ?protein
        WHERE {{
            VALUES ?target {{<http://purl.uniprot.org/embl-cds/{genpept}>}} 
            ?protein a up:Protein .
            ?protein rdfs:seeAlso ?target.
        }}
        """
            if genpept
            else None
        )

        uniprot_query = (
            f"""
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX up: <http://purl.uniprot.org/core/>
        SELECT ?protein
        WHERE {{
            VALUES ?prefix {{<http://purl.uniprot.org/database/EMBL>}}
            <http://purl.uniprot.org/uniprot/{uniprot}> rdfs:seeAlso ?protein .
            ?protein up:database ?prefix .
        }}
        """
            if uniprot
            else None
        )

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
            return results[0]["protein"]["value"]

        genpept_result = fetch_result(genpept_query) if genpept_query else None
        uniprot_result = fetch_result(uniprot_query) if uniprot_query else None

        if genpept and uniprot:
            if genpept_result != f"http://purl.uniprot.org/uniprot/{uniprot}":
                raise ValueError("The provided genpept ID and uniprot ID do not match")
            return {"genpept": genpept, "uniprot": uniprot}

        if genpept:
            return {"genpept": genpept, "uniprot": uniprot_result}
        elif uniprot:
            return {"uniprot_id": uniprot, "embl_id": genpept_result}

    # TODO (AR 2024-08-06): implement tests

    def validate_reaction_smarts(self: Self, reaction_smarts: str) -> str:
        """Validate an input reaction SMARTS

        Args:
            reaction_smarts: a reactions SMARTS string

        Returns:
            A validated reaction SMARTS

        Raises:
            ValueError: RDKit could not read SMILES
        """
        # TODO(MMZ 19.07.24): stub for implementation of e.g. a SMARTS validation function @adafede
