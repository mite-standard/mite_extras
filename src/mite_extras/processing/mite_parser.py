"""Parsing of data from mite json input format

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
from typing import Any, Self

from pydantic import BaseModel

from mite_extras.processing.data_classes import (
    Changelog,
    Cofactors,
    Entry,
    EnyzmeDatabaseIds,
    Enzyme,
    EnzymeAux,
    Evidence,
    Reaction,
    ReactionDatabaseIds,
    ReactionEx,
)

logger = logging.getLogger("mite_extras")


class MiteParser(BaseModel):
    """Assign data from mite input files to internal data structure.

    Attributes:
        entry: an Entry object to represent MITE data
    """

    entry: Any | None = None

    def to_json(self: Self) -> dict:
        """Prepare for export to json

        Returns:
            A dict formatted to be written to json

        Raises:
            RuntimeError: entry has not been assigned yet.
        """
        if self.entry:
            return self.entry.to_json()

    def to_html(self: Self) -> dict:
        """Prepare for export to html

        Returns:
            A dict formatted to be written to json
        """
        if self.entry:
            return self.entry.to_html()

    @staticmethod
    def get_cofactors(cofactors: dict | None) -> Cofactors | None:
        """Parse enzyme-related cofactors info

        Args:
            cofactors: a dict with cofactor info

        Returns:
            An Cofactors object or None
        """
        if isinstance(cofactors, dict):
            return Cofactors(**cofactors)

    @staticmethod
    def get_auxenzymes(auxenzymes: list | None) -> list | None:
        """Extract auxiliary enzyme info and converts into internal data structure

        Args:
            auxenzymes: a list of auxiliary enzyme dicts

        Returns:
            A list of EnzymeAux objects or None
        """
        if not auxenzymes:
            return

        return [
            EnzymeAux(
                name=a.get("name"),
                description=a.get("description"),
                databaseIds=EnyzmeDatabaseIds(**a.get("databaseIds")),
            )
            for a in auxenzymes
        ]

    @staticmethod
    def get_databaseids_reaction(
        data: dict | None,
    ) -> ReactionDatabaseIds | None:
        """Parse reaction-related databaseids info

        Args:
            data: a dict with database IDs

        Returns:
            An ReactionDatabaseIds object
        """
        if isinstance(data, dict):
            return ReactionDatabaseIds(**data)

    def get_reactions(self: Self, reactions: list) -> list:
        """Extract reactions infor and converts into internal data structure

        Args:
            reactions: list with reaction data

        Returns:
            A list of Reaction objects
        """
        return [
            Reaction(
                tailoring=rxn.get("tailoring"),
                description=rxn.get("description"),
                reactionSMARTS=rxn.get("reactionSMARTS"),
                reactions=[ReactionEx(**r) for r in rxn.get("reactions")],
                evidence=Evidence(
                    evidenceCode=rxn.get("evidence", {}).get("evidenceCode"),
                    references=rxn.get("evidence", {}).get("references"),
                ),
                databaseIds=self.get_databaseids_reaction(rxn.get("databaseIds")),
            )
            for rxn in reactions
        ]

    def parse_mite_json(self: Self, data: dict):
        """Parse data from mite-formatted json dict

        Args:
            data: a dict following the mite json schema
        """

        logger.debug("MiteParser: started creating Entry object.")

        try:
            self.entry = Entry(
                accession=data.get("accession"),
                status=data.get("status"),
                retirementReasons=data.get("retirementReasons"),
                changelog=[Changelog(**val) for val in data["changelog"]],
                enzyme=Enzyme(
                    name=data["enzyme"]["name"],
                    description=data.get("enzyme", {}).get("description"),
                    databaseIds=EnyzmeDatabaseIds(**data["enzyme"]["databaseIds"]),
                    auxiliaryEnzymes=self.get_auxenzymes(
                        auxenzymes=data.get("enzyme", {}).get("auxiliaryEnzymes")
                    ),
                    references=data.get("enzyme", {}).get("references"),
                    cofactors=self.get_cofactors(
                        cofactors=data.get("enzyme", {}).get("cofactors")
                    ),
                ),
                reactions=self.get_reactions(reactions=data.get("reactions")),
                comment=data.get("comment"),
            )
        except Exception as e:
            msg = str(e).split("For further information visit")[0].rstrip()
            msg = f"Error in parsing entry {data.get("accession")}: {msg}"
            raise ValueError(msg) from e

        logger.debug("MiteParser: completed creating Entry object.")
