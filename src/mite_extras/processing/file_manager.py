"""Validates and stores input and output file references

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
from pathlib import Path
from typing import Self

from pydantic import BaseModel, model_validator

logger = logging.getLogger("mite_extras")


class FileManager(BaseModel):
    """Pydantic-based class to validate and store input/output file references

    Attributes:
        indir: input directory path
        infiles: input files paths
        outdir: output directory path
    """

    indir: Path
    infiles: list = []
    outdir: Path

    @model_validator(mode="after")
    def validate_exists(self):
        """Validate existence of input/output dirs

        Raises:
            FileNotFoundError: one or more dirs could not be found
        """
        for path in [self.indir, self.outdir]:
            if not path.exists():
                logger.fatal(f"FileManager: could not find directory '{path.name}'")
                raise FileNotFoundError
        return self

    def read_files_indir(self: Self) -> None:
        """Read files in the indir, store to self"""
        for infile in self.indir.iterdir():
            if infile.suffix != ".json":
                continue
            self.infiles.append(infile)


# TODO(MMZ 17.07.24): Add method to export a single file to given dir
