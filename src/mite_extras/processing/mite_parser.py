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
import re
from pathlib import Path
from typing import Any, Self

import polars as pl
from pydantic import BaseModel

from mite_extras.processing.data_classes import (
    Changelog,
    ChangelogEntry,
    Entry,
    Enzyme,
    EnzymeAux,
    Evidence,
    Reaction,
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
        """Extracts changlog entries and converts in internal data structure

        Args:
            entries: list of changlog entry data

        Returns:
            A list of ChangelogEntry objects
        """
        log = []

        for entry in entries:
            log.append(
                ChangelogEntry(
                    contributors=list(entry.get("contributors")),
                    reviewers=list(entry.get("reviewers")),
                    date=entry.get("date"),
                    comment=entry.get("comment"),
                )
            )

        return log

    def get_changelog(self: Self, releases: list) -> list:
        """Extracts changelog and converts into internal data structure

        Args:
            releases: a list of changelog data

        Returns:
            A list of Changelog objects
        """
        log = []

        for release in releases:
            log.append(
                Changelog(
                    version=release.get("version"),
                    date=release.get("date"),
                    entries=self.get_changelog_entries(entries=release.get("entries")),
                )
            )

        return log

    def parse_mite_json(self: Self, data: dict):
        """Parse data from mite-formatted json dict

        Args:
            data: a dict following the mite json schema
        """

        self.entry = Entry(
            accession=data.get("accession"),
            quality=data.get("quality"),
            status=data.get("status"),
            retirementReasons=data.get("retirementReasons"),
            changelog=self.get_changelog(
                releases=data.get("changelog").get("releases")
            ),
        )

        # counting

        # changelog: n releases

        # how many reactions etc -count before
