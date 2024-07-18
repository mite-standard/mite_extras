"""Validation functionality.

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
from typing import Self

from pydantic import BaseModel
from rdkit.Chem import CanonSmiles, MolFromSmiles, MolToSmiles

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

    def canonicalize_smiles(self: Self, smiles: str) -> str:
        """Validate and canonicalize an input SMILES string

        Args:
            smiles: a SMILES string

        Returns:
            The SMILES in RDKit-canonicalized format

        Raises:
            ValueError: RDKit could not read SMILES
        """
        canon = MolFromSmiles(self.unescape_smiles(smiles))

        if canon is None:
            raise ValueError(f"Could not read SMILES string '{smiles}'")

        return self.unescape_smiles(CanonSmiles(MolToSmiles(canon)))
