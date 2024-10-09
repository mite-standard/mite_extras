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
    ChangelogEntry,
    Entry,
    EnyzmeDatabaseIds,
    Enzyme,
    EnzymeAux,
    Evidence,
    Reaction,
    ReactionDatabaseIds,
    ReactionEx,
    ReactionSmarts,
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
        else:
            raise RuntimeError(
                "'MiteParser': function 'to_json()' called out of order."
            )

    def to_html(self: Self) -> dict:
        """Prepare for export to html

        Returns:
            A dict formatted to be written to json
        """
        if self.entry:
            return self.entry.to_html()
        else:
            raise RuntimeError(
                "'MiteParser': function 'to_html()' called out of order."
            )

    @staticmethod
    def get_changelog_entries(entries: list) -> list:
        """Extract changlog entries and converts in internal data structure

        Args:
            entries: list of changlog entry data

        Returns:
            A list of ChangelogEntry objects
        """
        logger.debug("MiteParser: started creating ChangelogEntry object(s).")

        log = []

        for entry in entries:
            log.append(
                ChangelogEntry(
                    contributors=entry.get("contributors"),
                    reviewers=entry.get("reviewers"),
                    date=entry.get("date"),
                    comment=entry.get("comment"),
                )
            )

        logger.debug("MiteParser: completed creating ChangelogEntry object(s).")

        return log

    def get_changelog(self: Self, releases: list) -> list:
        """Extract changelog and converts into internal data structure

        Args:
            releases: a list of changelog data

        Returns:
            A list of Changelog objects
        """
        logger.debug("MiteParser: started creating Changelog object(s).")

        log = []

        for release in releases:
            log.append(
                Changelog(
                    version=release.get("version"),
                    date=release.get("date"),
                    entries=self.get_changelog_entries(entries=release.get("entries")),
                )
            )

        logger.debug("MiteParser: completed creating Changelog object(s).")

        return log

    @staticmethod
    def get_databaseids_enzyme(data: list | dict | None) -> EnyzmeDatabaseIds | None:
        """Parse enzyme-related databaseids info

        Args:
            data: a dict with database ids (list for MITE schema < 1.1)

        Returns:
            An EnyzmeDatabaseIds object
        """
        if isinstance(data, list):
            data_dict = {}
            for entry in data:
                if entry.startswith("genpept:"):
                    data_dict["genpept"] = entry.split(":")[1]
                elif entry.startswith("uniprot:"):
                    data_dict["uniprot"] = entry.split(":")[1]
                elif entry.startswith("mibig:"):
                    data_dict["mibig"] = entry.split(":")[1]
                else:
                    continue
            return EnyzmeDatabaseIds(**data_dict)
        elif isinstance(data, dict):
            return EnyzmeDatabaseIds(
                mibig=data.get("mibig"),
                uniprot=data.get("uniprot"),
                genpept=data.get("genpept"),
            )
        else:
            return None

    def get_auxenzymes(self: Self, auxenzymes: list | None) -> list | None:
        """Extract auxiliary enzyme info and converts into internal data structure

        Args:
            auxenzymes: a list of auxiliary enzyme dicts

        Returns:
            A list of EnzymeAux objects
        """
        if auxenzymes is None:
            return None

        logger.debug("MiteParser: started creating EnzymeAux object(s).")

        log = []

        for auxenz in auxenzymes:
            log.append(
                EnzymeAux(
                    name=auxenz.get("name"),
                    description=auxenz.get("description"),
                    databaseIds=self.get_databaseids_enzyme(
                        data=auxenz.get("databaseIds")
                    ),
                )
            )

        logger.debug("MiteParser: complete creating EnzymeAux object(s).")

        return log

    @staticmethod
    def get_databaseids_reaction(
        data: list | dict | None,
    ) -> ReactionDatabaseIds | None:
        """Parse reaction-related databaseids info

        Args:
            data: a dict with database IDs (list for MITE schema < 1.1)

        Returns:
            An ReactionDatabaseIds object
        """
        if isinstance(data, list):
            data_dict = {"rhea": [], "ec": [], "mite": []}
            for entry in data:
                if entry.startswith("rhea:"):
                    data_dict["rhea"].append(entry.split(":")[1])
                elif entry.startswith("EC"):
                    data_dict["ec"].append(entry)
                elif entry.startswith("MITE"):
                    data_dict["mite"].append(entry)
                else:
                    continue
            return ReactionDatabaseIds(**data_dict)
        elif isinstance(data, dict):
            return ReactionDatabaseIds(
                rhea=data.get("rhea"), ec=data.get("ec"), mite=data.get("mite")
            )
        else:
            return None

    @staticmethod
    def get_reactionex(reactions: list) -> list:
        """Extract experimental reaction info, converts to internal data structure

        Args:
            reactions: a list of experimental reaction data

        Returns:
            A list of ReactionEx objects
        """
        logger.debug("MiteParser: started creating ReactionEx object(s).")

        log = []

        for reaction in reactions:
            log.append(
                ReactionEx(
                    substrate=reaction.get("substrate"),
                    products=reaction.get("products"),
                    forbidden_products=reaction.get("forbidden_products"),
                    isIntermediate=reaction.get("isIntermediate"),
                    description=reaction.get("description"),
                )
            )

        logger.debug("MiteParser: completed creating ReactionEx object(s).")

        return log

    @staticmethod
    def get_evidence(evidences: list) -> list:
        """Extract evidence information, convert into internal data structure

        Args:
            evidences: a list with evidence information

        Returns:
            A list of Evidence objects
        """
        logger.debug("MiteParser: started creating Evidence object(s).")

        log = []

        for evidence in evidences:
            log.append(
                Evidence(
                    evidenceCode=evidence.get("evidenceCode"),
                    references=evidence.get("references"),
                )
            )

        logger.debug("MiteParser: started creating Evidence object(s).")

        return log

    def get_reactions(self: Self, reactions: list) -> list:
        """Extract reactions infor and converts into internal data structure

        Args:
            reactions: list with reaction data

        Returns:
            A list of Reaction objects
        """
        logger.debug("MiteParser: started creating Reaction object(s).")

        log = []

        for reaction in reactions:
            log.append(
                Reaction(
                    tailoring=reaction.get("tailoring"),
                    description=reaction.get("description"),
                    reactionSMARTS=ReactionSmarts(
                        reactionSMARTS=reaction.get("reactionSMARTS").get(
                            "reactionSMARTS"
                        )
                    ),
                    reactions=self.get_reactionex(reactions=reaction.get("reactions")),
                    evidence=self.get_evidence(evidences=reaction.get("evidence")),
                    databaseIds=self.get_databaseids_reaction(
                        reaction.get("databaseIds")
                    ),
                )
            )

        logger.debug("MiteParser: completed creating Reaction object(s).")

        return log

    def parse_mite_json(self: Self, data: dict):
        """Parse data from mite-formatted json dict

        Args:
            data: a dict following the mite json schema
        """

        logger.debug("MiteParser: started creating Entry object.")

        self.entry = Entry(
            accession=data.get("accession"),
            status=data.get("status"),
            retirementReasons=data.get("retirementReasons"),
            changelog=self.get_changelog(
                releases=data.get("changelog").get("releases")
            ),
            enzyme=Enzyme(
                name=data.get("enzyme").get("name"),
                description=data.get("enzyme").get("description"),
                databaseIds=self.get_databaseids_enzyme(
                    data=data.get("enzyme").get("databaseIds")
                ),
                auxiliaryEnzymes=self.get_auxenzymes(
                    auxenzymes=data.get("enzyme").get("auxiliaryEnzymes")
                ),
                references=data.get("enzyme").get("references"),
            ),
            reactions=self.get_reactions(reactions=data.get("reactions")),
            comment=data.get("comment"),
            attachments=data.get("attachments"),
        )

        logger.debug("MiteParser: completed creating Entry object.")
