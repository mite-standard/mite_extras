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
from contextlib import suppress
from itertools import permutations, product
from math import pi

import requests
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import (
    AddHs,
    AssignStereochemistry,
    CanonSmiles,
    GetMolFrags,
    Mol,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmarts,
    MolToSmiles,
    RemoveHs,
    SanitizeMol,
    rdMolEnumerator,
)
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts

logger = logging.getLogger("mite_extras")


class MoleculeValidator(BaseModel):
    """Handles basic molecule validation and canonicalization."""

    @staticmethod
    def _clean_string(string: str) -> str:
        """Remove superfluous backslashes and Ketcher-style hydrogen annotations from string."""
        string = string.replace("\\\\", "\\")
        # remove ketcher-style hydrogen annotations like ';h1' or '&H0'
        string = re.sub(r";h\d", "", string, flags=re.IGNORECASE)
        string = re.sub(r"&H\d+", "", string, flags=re.IGNORECASE)
        return string

    def canonicalize_smiles(self, smiles: str) -> str:
        """Canonicalize a SMILES string."""
        s = self._clean_string(smiles)
        mol = MolFromSmiles(s)
        if mol is None:
            # Attempt a lenient parse and try to repair common issues.
            # Many failures stem from incorrect aromatic nitrogen hydrogen
            # annotations exported by Ketcher (e.g. '[nH]' vs '[n]'). Try a
            # non-sanitized parse first, then sanitize, and if sanitizing
            # fails try simple textual repairs on the SMILES.
            try:
                mol = MolFromSmiles(s, sanitize=False)
                if mol is None:
                    raise ValueError("parse returned None")

                # try sanitizing the molecule to compute implicit valences
                try:
                    SanitizeMol(mol)
                except Exception:
                    # Attempt simple repairs on common Ketcher issues and
                    # re-parse/sanitize. Stop at the first successful repair.
                    repair_candidates = [
                        ("[nH]", "[n]"),
                        ("[n]", "[nH]"),
                    ]
                    # also try regex-based targeted repairs
                    regex_candidates = [
                        (r"\[nH\](?=\()", "[n]"),
                        (r"\[nH\](?=c)", "[n]"),
                        (r"\[nH\](?=C)", "[n]"),
                    ]
                    repaired = None
                    import re as _re

                    def try_parse(sm):
                        try:
                            m = MolFromSmiles(sm)
                        except Exception:
                            m = None
                        if m is None:
                            try:
                                m = MolFromSmiles(sm, sanitize=False)
                                SanitizeMol(m)
                            except Exception:
                                m = None
                        return m

                    for old, new in repair_candidates:
                        s2 = s.replace(old, new)
                        m2 = try_parse(s2)
                        if m2 is not None:
                            repaired = m2
                            s = s2
                            mol = m2
                            break

                    if repaired is None:
                        for pat, new in regex_candidates:
                            s2 = _re.sub(pat, new, s)
                            if s2 == s:
                                continue
                            m2 = try_parse(s2)
                            if m2 is not None:
                                repaired = m2
                                s = s2
                                mol = m2
                                break

                    if repaired is None:
                        # Re-raise the sanitize exception to be handled below
                        raise

                # At this point sanitization succeeded; add/remove Hs to
                # normalize explicit hydrogens if necessary.
                mol = AddHs(mol)
                SanitizeMol(mol)
                mol = RemoveHs(mol)
            except Exception as e:
                raise ValueError(
                    f"RDKit rejected SMILES string - is it a valid SMILES?\n{s}"
                ) from e
        # Ensure stereochemistry is assigned consistently before producing SMILES
        with suppress(Exception):
            # Force stereochemistry assignment so canonicalization is deterministic
            AssignStereochemistry(mol, force=True)

        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return CanonSmiles(MolToSmiles(mol))

    def canonicalize_smarts(self, smarts: str) -> str:
        """Canonicalize a SMARTS pattern."""
        mol = MolFromSmarts(self._clean_string(smarts))
        if mol is None:
            raise ValueError(
                f"RDKit rejected SMARTS string - is it a valid pattern?\n{smarts}"
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
            # Remove Ketcher-style hydrogen annotations (e.g. &H0)
            r"&H\d+": r"",
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
                f"RDKit rejected SMARTS string - is it a valid pattern?\n{smarts}"
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
        predicted_smiles = set()
        for p in predicted_products:
            if p is None:
                continue

            # Add canonical SMILES for all enumerated variants of the product
            try:
                variants = self.enumerator.enumerate_molecule(p)
            except Exception:
                variants = {p}

            for var in variants:
                try:
                    s = self.molecule_validator.canonicalize_smiles(MolToSmiles(var))
                    predicted_smiles.add(s)
                except Exception:
                    continue

                # Also try flipping combinations of chiral centers (up to a reasonable limit)
                try:
                    chiral_atoms = [
                        a.GetIdx()
                        for a in var.GetAtoms()
                        if a.GetChiralTag()
                        in (
                            ChiralType.CHI_TETRAHEDRAL_CW,
                            ChiralType.CHI_TETRAHEDRAL_CCW,
                        )
                    ]
                    max_combinations = (
                        1 << len(chiral_atoms) if len(chiral_atoms) <= 12 else 1 << 12
                    )
                    from itertools import product

                    # For each combination of flips, generate a variant
                    for bits in product([0, 1], repeat=min(len(chiral_atoms), 12)):
                        mol_copy = Chem.Mol(var)
                        changed = False
                        for idx, bit in zip(chiral_atoms, bits, strict=True):
                            if bit:
                                atom = mol_copy.GetAtomWithIdx(idx)
                                tag = atom.GetChiralTag()
                                if tag == ChiralType.CHI_TETRAHEDRAL_CW:
                                    atom.SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CCW)
                                    changed = True
                                elif tag == ChiralType.CHI_TETRAHEDRAL_CCW:
                                    atom.SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CW)
                                    changed = True
                        if not changed:
                            continue
                        try:
                            s_en = self.molecule_validator.canonicalize_smiles(
                                MolToSmiles(mol_copy)
                            )
                            predicted_smiles.add(s_en)
                        except Exception:
                            continue
                except Exception:
                    pass

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
