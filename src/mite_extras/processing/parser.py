"""Parsing of data from input format

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
from typing import Self

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


class Parser(BaseModel):
    """Collection of parsers to assign input data to internal data structure"""

    @staticmethod
    def remove_escaped_quote(instring: str) -> str:
        return instring.replace('"', "")

    @staticmethod
    def parse_mite_json(input_data: dict) -> dict:
        """Read in data from a json file formatted after mite schema

        Args:
            input_data: a dict coming from a mite json file

        Returns:
            A validated and mite-formatted dict for json export
        """
        # TODO(MMZ 18.07.24): Needs to be implemented

        logger.fatal("MITE JSON PARSER CURRENTLY NOT IMPLEMENTED")

        return {}

    def parse_raw_json(self: Self, name: str, input_data: dict) -> dict:
        """Read in data from a json file resulting from the 2024 MIBiG annotathons

        Args:
            name: the name of the file
            input_data: a dict coming from the MIBiG 2024 submission portal

        Returns:
            A validated and mite-formatted dict for json export

        Raises:
            ValueError: data cannot be found
        """

        tailoring_list = input_data.get("Tailoring", [])

        if len(tailoring_list) == 0:
            raise ValueError(f"Parser: file '{name}' does not contain MITE data.")

        tailoring_dict = {}
        for entry in tailoring_list:
            tailoring_dict[entry[0]] = [1]

        reactions = set()
        for key, _value in tailoring_dict.items():
            match = re.match(r"^(enzymes-0-reactions-\d+)", key)
            if match:
                reactions.add(match.group(1))

        # figure out how many reactions included to be able to build key
        # then loop through whole dict and pull out via key, value and regexp

        return reactions

        # figure out how many entries in the list to get the keys/n to loop through

        entry = Entry()

        return entry.to_json()
