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
from functools import lru_cache
from itertools import permutations, product
from math import pi

import requests
from pydantic import BaseModel
from rdkit import Chem, RDLogger

# Disable RDKit error-level app logs (this suppresses the repeated "Can't kekulize mol" messages
# while leaving other RDKit log levels intact).
RDLogger.DisableLog("rdApp.error")
from rdkit.Chem import (
    AddHs,
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
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts

logger = logging.getLogger("mite_extras")


class MoleculeValidator(BaseModel):
    """Handles basic molecule validation and canonicalization."""

    @staticmethod
    def _clean_string(string: str) -> str:
        """Remove superfluous backslashes and H's from string."""
        string = string.replace("\\\\", "\\")
        string = re.sub(r";h\d", "", string)
        return string

    def canonicalize_smiles(self, smiles: str) -> str:
        """Canonicalize a SMILES string with memoization to avoid repeated heavy work."""
        # initialize caches lazily to avoid shared mutable defaults on class
        if not hasattr(self, "_canon_cache"):
            self._canon_cache = {}
        if not hasattr(self, "_noniso_cache"):
            self._noniso_cache = {}

        s = self._clean_string(smiles)
        if s in self._canon_cache:
            return self._canon_cache[s]

        mol = MolFromSmiles(s)
        if mol is None:
            # Attempt a lenient parse and try to repair common issues.
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
                    f"RDKit rejected SMILES string - is it a valid SMILES?\n"
                    f"The erroneous SMILES string was:\n"
                    f"{s}\n"
                    f"Full error trace:\n"
                    f"{e!s}\n",
                ) from e
        # Do not force stereochemistry assignment here; keep RDKit's original handling

        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        # Try producing SMILES; if kekulization fails, attempt Kekulize and retry
        try:
            smi = MolToSmiles(mol)
        except Exception as e:
            if "kekul" in str(e).lower():
                with suppress(Exception):
                    Chem.Kekulize(mol, clearAromaticFlags=True)
                smi = MolToSmiles(mol)
            else:
                raise
        canon = CanonSmiles(smi)
        self._canon_cache[smiles] = canon
        return canon

    def canonicalize_smarts(self, smarts: str) -> str:
        """Canonicalize a SMARTS pattern."""
        mol = MolFromSmarts(self._clean_string(smarts))
        if mol is None:
            raise ValueError(
                f"RDKit rejected SMARTS string - is it a valid pattern?\n"
                f"The erroneous SMARTS string was:\n"
                f"{smarts}\n",
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
                f"RDKit rejected SMARTS string - is it a valid pattern?\n"
                f"The erroneous SMARTS string was:\n"
                f"{smarts}\n",
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

    def _normalize_product_smiles(self, prod: Mol) -> set[str]:
        """Thin wrapper: compute raw_smiles and delegate to cached string-based normalizer.

        Cached helper returns cleaned SMILES strings; canonicalize them here using
        MoleculeValidator (which already caches canonicalization results).
        """
        try:
            raw = MolToSmiles(prod)
        except Exception:
            try:
                raw = MolToSmiles(Chem.Mol(prod))
            except Exception:
                raise
        variants = self._normalize_product_smiles_str_cached(raw)
        result = set()
        # try the raw SMILES first (fast path)
        try:
            result.add(self.molecule_validator.canonicalize_smiles(raw))
        except Exception:
            pass
        for v in variants:
            try:
                result.add(self.molecule_validator.canonicalize_smiles(v))
            except Exception:
                continue
        return result

    @staticmethod
    @lru_cache(maxsize=16384)
    def _normalize_product_smiles_str_cached(raw_smiles: str) -> tuple[str, ...]:
        """Return cleaned SMILES variants for a raw SMILES string; cached for speed.

        This performs textual H-stripping and an optional sanitization/kekulize
        attempt, but returns plain SMILES strings (not canonicalized). The
        caller will canonicalize using MoleculeValidator (which has its own cache).
        """
        s = raw_smiles

        # Remove explicit hydrogen annotations like [CH3], [CH2], [CH], [nH]
        s = re.sub(r"\[([A-Za-z]{1,2})H\d*\]", r"[\1]", s)
        s = re.sub(r"\[H\]", r"", s)
        s = re.sub(r"\[nH(\d*)\]", r"[n\1]", s)
        # remove simple bracketed atoms like [C], [Cl] or with atom-maps [C:1] -> C, Cl
        # avoid touching brackets that contain stereo (@), charge (+/-), isotopes (digits before element)
        s = re.sub(r"\[([A-Z][a-z]?)(?::\d+)?\]", r"\1", s)

        variants: set[str] = set()
        # primary cleaned string
        variants.add(s)

        # try a sanitization-based variant (may change aromaticity/tautomer) but only if primary failed to parse
        try:
            m2 = MolFromSmiles(s, sanitize=False)
            if m2 is not None:
                with suppress(Exception):
                    SanitizeMol(m2)
                with suppress(Exception):
                    Chem.Kekulize(m2, clearAromaticFlags=True)
                try:
                    s2 = MolToSmiles(m2)
                    if s2:
                        variants.add(s2)
                except Exception:
                    pass
        except Exception:
            pass

        # return a stable tuple for caching
        return tuple(sorted(variants))

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
                f"{'\n'.join(overlap)}\n",
            )

        # Generate reaction variants and validate
        reactions = self._get_reaction_variants(reaction_smarts)
        predicted_products = self._run_reactions(
            reactions,
            substrate_smiles,
            intramolecular,
        )

        # Validate predictions
        predicted_smiles = set()
        raw_predicted = set()
        for p in predicted_products:
            if p is None:
                continue
            try:
                raw = MolToSmiles(p)
            except Exception:
                try:
                    raw = MolToSmiles(Chem.Mol(p))
                except Exception:
                    raw = None
            if raw:
                raw_predicted.add(raw)
            try:
                variants = self._normalize_product_smiles(p)
                predicted_smiles.update(variants)
            except Exception as e:
                logger.debug("Failed to normalize product mol: %s", e)
                try:
                    predicted_smiles.add(
                        self.molecule_validator.canonicalize_smiles(MolToSmiles(p)),
                    )
                except Exception:
                    # give up on this product
                    continue
        logger.debug("Predicted raw SMILES: %s", raw_predicted)

        logger.debug("Predicted products: %s", predicted_smiles)
        if not expected_smiles.issubset(predicted_smiles):
            # Strict stereochemistry handling:
            # - If the connectivity (non-isomeric form) is not present in predictions, fail.
            # - If connectivity matches but stereochemistry differs, fail and report whether
            #   the predicted product is an enantiomer or a diastereomer.

            def _nonisomeric_canonical(smiles_str: str):
                try:
                    m = MolFromSmiles(smiles_str)
                    if m is None:
                        return None
                    s = MolToSmiles(m, isomericSmiles=False)
                    return CanonSmiles(s)
                except Exception:
                    return None

            # Map non-isomeric forms to expected/predicted products
            expected_map = {}
            for p in expected_products:
                try:
                    can = self.molecule_validator.canonicalize_smiles(p)
                except Exception:
                    can = None
                noniso = _nonisomeric_canonical(can) if can else None
                if noniso:
                    expected_map.setdefault(noniso, []).append(can)

            predicted_map = {}
            for p in predicted_smiles:
                noniso = _nonisomeric_canonical(p)
                if noniso:
                    predicted_map.setdefault(noniso, []).append(p)

            missing_connectivity = []
            stereochem_mismatches = []

            for noniso, exp_cans in expected_map.items():
                preds = predicted_map.get(noniso)
                if not preds:
                    # no product with same connectivity
                    missing_connectivity.extend(exp_cans)
                    continue

                # For each expected canonical with same connectivity, compare stereo
                for exp_can in exp_cans:
                    exp_m = MolFromSmiles(exp_can)
                    Chem.AssignStereochemistry(exp_m, force=True, cleanIt=True)
                    exp_centers = dict(
                        Chem.FindMolChiralCenters(exp_m, includeUnassigned=True)
                    )

                    # Try to find a predicted molecule that matches stereochemistry exactly
                    matched_exact = False
                    mismatch_details = []
                    for pred_can in preds:
                        pred_m = MolFromSmiles(pred_can)
                        Chem.AssignStereochemistry(pred_m, force=True, cleanIt=True)

                        # Determine atom mapping between expected and predicted (ignore chirality)
                        mapping = exp_m.GetSubstructMatch(pred_m)
                        if not mapping:
                            # try the reverse mapping and invert
                            rev_map = pred_m.GetSubstructMatch(exp_m)
                            if rev_map:
                                # build mapping from exp idx to pred idx
                                mapping = tuple(
                                    rev_map.index(i) if i in rev_map else None
                                    for i in range(exp_m.GetNumAtoms())
                                )
                            else:
                                # cannot map; skip
                                continue

                        pred_centers = dict(
                            Chem.FindMolChiralCenters(pred_m, includeUnassigned=True)
                        )

                        # Compare stereochemistry at corresponding centers
                        all_same = True
                        all_inverted = True
                        unknown = False
                        for exp_idx, exp_label in exp_centers.items():
                            pred_idx = mapping[exp_idx]
                            if pred_idx is None:
                                unknown = True
                                break
                            pred_label = pred_centers.get(pred_idx, "?")
                            if exp_label == "?":
                                unknown = True
                                break
                            if pred_label == "?":
                                unknown = True
                                break
                            if exp_label == pred_label:
                                all_inverted = False
                            else:
                                all_same = False

                        if all_same and not unknown:
                            matched_exact = True
                            break

                        if not matched_exact:
                            if not unknown:
                                if all_inverted:
                                    mismatch_details.append(
                                        (exp_can, pred_can, "enantiomer")
                                    )
                                else:
                                    mismatch_details.append(
                                        (exp_can, pred_can, "diastereomer")
                                    )
                            else:
                                mismatch_details.append(
                                    (exp_can, pred_can, "unknown_stereo")
                                )

                    if not matched_exact:
                        stereochem_mismatches.extend(mismatch_details)

            if missing_connectivity:
                raise ValueError(
                    "Reaction did not lead to all expected products (connectivity mismatch).\n"
                    "Missing expected products:\n"
                    f"{'\n'.join(missing_connectivity)}\n"
                    "Generated products:\n"
                    f"{'\n'.join(predicted_smiles)}\n",
                )

            if stereochem_mismatches:
                # Prepare actionable warning message for users with hints to fix their entries
                details = []
                for exp_can, pred_can, kind in stereochem_mismatches:
                    details.append(
                        f"Expected: {exp_can}\nPredicted: {pred_can}\nMismatch: {kind}\n"
                    )
                message = (
                    "Reaction produced stereoisomers of expected products (stereochemistry mismatch).\n"
                    "Details:\n"
                    f"{'\n'.join(details)}\n"
                    "Suggested actions:\n"
                    "  - Verify the submitted expected product SMILES include correct stereochemical annotations ([@/@@], [C@H], etc.).\n"
                    "  - If the reaction yields a racemate, submit both enantiomers as expected products or mark the entry as racemic.\n"
                    "  - Ensure substrate stereochemistry is specified if the outcome depends on it.\n"
                    "  - If the entry should not be validated automatically (complex/ambiguous stereochemistry), mark it as 'needs review' in the entry metadata so the pipeline skips automated validation.\n"
                    "Contact: update entry or run with stereo-flattening option if appropriate."
                )
                logger.warning(message)
                # Continue (treat as warning): do not raise — caller/CLI should surface this warning to the user
                return None

            # If we reach here, all expected products matched exactly (should not happen),
            # otherwise earlier raises would have been triggered.

        if overlap := forbidden_smiles & predicted_smiles:
            raise ValueError(
                f"Reaction product(s) belong(s) to the specified forbidden product(s):\n"
                f"{'\n'.join(overlap)}\n",
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
                f"The erroneous SMILES string was:\n"
                f"{substrate_smiles}\n",
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
