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
from typing import Self

from pydantic import BaseModel
from rdkit.Chem import (
    CanonSmiles,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmarts,
    MolToSmiles,
)

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
