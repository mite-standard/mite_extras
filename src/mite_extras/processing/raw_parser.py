"""Parsing of data from input raw json input format

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


class RawParser(BaseModel):
    """Assign data from raw input files to internal data structure

    Attributes:
        entry: an Entry object to represent MITE data

    # TODO(MMZ 19.07.24): Retire this class once the MIBiG Submission portal produces valid MITE jsons
    """

    entry: Any | None = None

    @staticmethod
    def remove_quote(instring: str) -> str:
        return instring.replace('"', "")

    @staticmethod
    def remove_empty_string(instring: str) -> str | None:
        if instring == "":
            return None
        else:
            return instring

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
            raise RuntimeError("'RawParser': function 'to_json()' called out of order.")

    def parse_raw_json(self: Self, name: str, input_data: dict):
        """Read in data from a json file resulting from the 2024 MIBiG annotathons

        Args:
            name: the name of the file
            input_data: a dict coming from the MIBiG 2024 submission portal

        Returns:
            A validated and mite-formatted dict for json export

        Raises:
            ValueError: "Tailoring" data cannot be found
            RuntimeError: "No reaction data found"
        """

        def _compile_aux_enzymes(data: dict) -> list | None:
            """Compile a list of auxiliary enzymes

            Args:
                data: a dict coming from the MIBiG 2024 submission portal

            Returns:
                A list of strings or None (if none were specified
            """
            logger.debug("Parser: started compiling AuxEnzyme info.")

            aux_enz = set()
            for key, _value in data.items():
                match1 = re.match(r"^enzymes-0-enzyme-0-auxiliary_enzymes-(\d+)", key)
                if match1:
                    aux_enz.add(match1.group(1))

            aux_list = []
            for enz in aux_enz:
                aux_list.append(
                    EnzymeAux(
                        name=data.get(
                            f"enzymes-0-enzyme-0-auxiliary_enzymes-{enz}-name"
                        ),
                        description=data.get(
                            f"enzymes-0-enzyme-0-auxiliary_enzymes-{enz}-description"
                        ),
                        databaseIds=[
                            self.remove_quote(string)
                            for string in data.get(
                                f"enzymes-0-enzyme-0-auxiliary_enzymes-{enz}-databaseIds"
                            ).split(", ")
                        ],
                    )
                )

            logger.debug("Parser: completed compiling AuxEnzyme info.")

            if len(aux_list) > 0:
                return aux_list
            else:
                return None

        def _compile_changelog_data(data: dict) -> tuple:
            """Compile the changelog from scratch

            Several params are hardcoded since they are swapped in downstream
             (contributor ID, version, date)

            Args:
                data: a dict coming from the MIBiG 2024 submission portal

            Returns:
                Tuple of (date, contributor list)

            Raises:
                ValueError: file has no changelog information
            """
            logger.debug("Parser: started compiling changelog.")

            changelog_list = data.get("Changelog", [])

            df = pl.read_csv(
                Path(__file__).parent.parent.joinpath("schema/mibig_id_mappings.csv")
            )

            if len(changelog_list) == 0:
                raise ValueError(f"Parser: file '{name}' does not contain a Changelog.")

            contributors = set()
            for i in changelog_list:
                editor_int = int(i.get("Edited_by"))
                editor_hash = df.filter(pl.col("id_int") == editor_int)["id_hash"][0]
                contributors.add(editor_hash)

            date = changelog_list[-1].get("Edited_at").split(",")[0]
            date = date.split("/")
            date = f"{date[2]}-{date[0]}-{date[1]}"

            logger.debug("Parser: completed compiling changelog.")

            return date, list(contributors)

        def _compile_tailoring_function(data: dict, count: str) -> list:
            """Compile list of tailoring functions

            Args:
                data: a dict coming from the MIBiG 2024 submission portal
                count: index of reaction in reactions

            Returns:
                A list of tailoring functions
            """
            logger.debug("Parser: started compiling tailoring functions.")

            number_t_functions = set()
            for key, _value in data.items():
                match = re.match(rf"^enzymes-0-reactions-{count}-tailoring-(\d+)", key)
                if match:
                    number_t_functions.add(match.group(1))

            tailoring_funct = set()
            for funct in number_t_functions:
                tailoring_funct.add(
                    data.get(f"enzymes-0-reactions-{count}-tailoring-{funct}-function")
                )

            logger.debug("Parser: completed compiling tailoring functions.")

            return list(tailoring_funct)

        def _compile_reaction_dbids(data: dict, count: str) -> list | None:
            """Compile the databaseIDs of a reaction

            databaseIDs are parsed from both 'validated_reactions' and 'reaction_smarts'
            due to simplification of schema from previous version.

            Args:
                data: a dict coming from the MIBiG 2024 submission portal
                count: index of reaction in reactions

            Returns:
                A list of databaseIDs
            """
            logger.debug("Parser: started compiling reaction databaseIDs.")

            dbids = set()
            for key, value in data.items():
                match = re.match(
                    rf"^enzymes-0-reactions-{count}-validated_reactions-\d+-databaseIds|enzymes-0-reactions-{count}-reaction_smarts-\d+-databaseIds",
                    key,
                )
                if match:
                    cleaned_val = self.remove_empty_string(value)
                    if cleaned_val is not None:
                        dbids.add(value)

            if len(dbids) > 0:
                db_entries = set()
                for db in dbids:
                    db_entries.update(
                        [self.remove_quote(string) for string in db.split(", ")]
                    )
                logger.debug("Parser: completed compiling reaction databaseIDs.")
                return list(db_entries)
            else:
                logger.debug("Parser: completed compiling reaction databaseIDs.")
                return None

        def _compile_evidence(data: dict, count: str) -> Evidence:
            """Compile the evidence of a reaction

            Args:
                data: a dict coming from the MIBiG 2024 submission portal
                count: index of reaction in reactions

            Returns:
                A single Evidence object
            """
            logger.debug("Parser: started compiling evidence data.")

            codes = set()
            for key, value in data.items():
                match = re.match(
                    rf"^enzymes-0-reactions-{count}-validated_reactions-\d+-evidence_val-\d+-evidenceCode",
                    key,
                )
                if match:
                    codes.update(
                        [
                            self.remove_quote(string)
                            for string in re.split(", |,", value)
                        ]
                    )

            refs = set()
            for key, value in data.items():
                match = re.match(
                    rf"^enzymes-0-reactions-{count}-validated_reactions-\d+-evidence_val-\d+-references",
                    key,
                )
                if match:
                    refs.update(
                        [
                            self.remove_quote(string)
                            for string in re.split(", |,", value)
                        ]
                    )

            logger.debug("Parser: completed compiling evidence data.")

            return Evidence(
                evidenceCode=[
                    s for s in codes if self.remove_empty_string(s) is not None
                ],
                references=[s for s in refs if self.remove_empty_string(s) is not None],
            )

        def _compile_ex_reactions(data: dict, count: str, count_ex: set) -> list:
            """Compile the experimentally verified reaction examples

            Args:
                data: a dict coming from the MIBiG 2024 submission portal
                count: index of reaction in reactions
                count_ex: index of the example reaction of reaction in reactions

            Returns:
                A list of ReactionEx instances
            """
            logger.debug("Parser: started compiling experimental reactions.")

            ex_reactions = []
            for c_ex in count_ex:
                products = set()
                for key, value in data.items():
                    match = re.match(
                        rf"^enzymes-0-reactions-{count}-validated_reactions-{c_ex}-product_substructure-\d+",
                        key,
                    )
                    if match:
                        products.add(value)

                balanced = False
                if data.get(
                    f"enzymes-0-reactions-{count}-validated_reactions-{c_ex}-isBalanced"
                ) in {"yes", "y", "Y", "Yes", "YES"}:
                    balanced = True

                intermediate = False
                if data.get(
                    f"enzymes-0-reactions-{count}-validated_reactions-{c_ex}-isIntermediate"
                ) in {"yes", "y", "Y", "Yes", "YES"}:
                    intermediate = True

                ex_reactions.append(
                    ReactionEx(
                        substrate=data.get(
                            f"enzymes-0-reactions-{count}-validated_reactions-{c_ex}-substrate_substructure"
                        ),
                        products=list(products),
                        isBalanced=balanced,
                        isIntermediate=intermediate,
                        description=self.remove_empty_string(
                            data.get(
                                f"enzymes-0-reactions-{count}-validated_reactions-{c_ex}-description"
                            )
                        ),
                    )
                )

            logger.debug("Parser: completed compiling experimental reactions.")

            return ex_reactions

        def _compile_reactions(data: dict) -> list:
            """Compile the reactions in the MITE entry

            Args:
                data: a dict coming from the MIBiG 2024 submission portal

            Returns:
                A list of Reaction instances
            """
            logger.debug("Parser: started compiling reactions.")

            reactions = set()
            for key, _value in data.items():
                match = re.match(r"^enzymes-0-reactions-(\d+)", key)
                if match:
                    reactions.add(match.group(1))

            reaction_list = []
            for reaction in reactions:
                ex_react = set()
                for key, _value in data.items():
                    match = re.match(
                        rf"^enzymes-0-reactions-{reaction}-validated_reactions-(\d+)",
                        key,
                    )
                    if match:
                        ex_react.add(match.group(1))

                reaction_list.append(
                    Reaction(
                        tailoring=_compile_tailoring_function(
                            data=data, count=reaction
                        ),
                        databaseIds=_compile_reaction_dbids(data=data, count=reaction),
                        evidence=[_compile_evidence(data=data, count=reaction)],
                        reactionSMARTS=ReactionSmarts(
                            reactionSMARTS=data.get(
                                f"enzymes-0-reactions-{reaction}-reaction_smarts-0-reactionSMARTS"
                            ),
                            isIterative=False,
                        ),
                        reactions=_compile_ex_reactions(
                            data=data, count=reaction, count_ex=ex_react
                        ),
                        description=self.remove_empty_string(
                            data.get(f"enzymes-0-reactions-{reaction}-description")
                        ),
                    )
                )

            logger.debug("Parser: completed compiling reactions.")

            return reaction_list

        tailoring_list = input_data.get("Tailoring", [])

        if len(tailoring_list) == 0:
            raise ValueError(f"Parser: file '{name}' does not contain MITE data.")

        tailoring_dict = {}
        for entry in tailoring_list:
            tailoring_dict[entry[0]] = entry[1]

        entry = Entry(
            accession=name,
            quality="medium",
            status="pending",
            comment=self.remove_empty_string(tailoring_dict.get("enzymes-0-comment")),
            changelog=[
                Changelog(
                    date="2024-07-30",
                    version="1.0",
                    entries=[
                        ChangelogEntry(
                            contributors=_compile_changelog_data(input_data)[1],
                            reviewers=["AAAAAAAAAAAAAAAAAAAAAAAA"],
                            date=_compile_changelog_data(input_data)[0],
                            comment="Initial entry.",
                        )
                    ],
                )
            ],
            enzyme=Enzyme(
                name=tailoring_dict.get("enzymes-0-enzyme-0-name"),
                description=tailoring_dict.get("enzymes-0-enzyme-0-description"),
                databaseIds=[
                    self.remove_quote(string)
                    for string in re.split(
                        ", |,", tailoring_dict.get("enzymes-0-enzyme-0-databaseIds")
                    )
                    if self.remove_empty_string(string) is not None
                ],
                auxiliaryEnzymes=_compile_aux_enzymes(tailoring_dict),
                references=[
                    (self.remove_quote(string))
                    for string in re.split(
                        ", |,", tailoring_dict.get("enzymes-0-enzyme-0-references")
                    )
                    if self.remove_empty_string(string) is not None
                ],
            ),
            reactions=_compile_reactions(tailoring_dict),
        )

        entry.enzyme.databaseIds.append(
            f'mibig:{input_data.get("Metadata").get("mibig_id")}'
        )

        if len(entry.reactions) == 0:
            raise RuntimeError("No reaction data found - SKIP.")

        self.entry = entry
