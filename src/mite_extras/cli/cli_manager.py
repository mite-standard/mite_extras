"""Command line interface manager.

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

import argparse
from importlib import metadata
from typing import Self

from pydantic import BaseModel


class CliManager(BaseModel):
    """Manage command line interface and parsing."""

    def run(self: Self, args: list) -> argparse.Namespace:
        """Run command line interface using argparse.

        Arguments:
            args: specified arguments

        Returns:
            argparse.Namespace object with command line parameters
        """
        parser = self.define_cli_args()
        return parser.parse_args(args)

    @staticmethod
    def define_cli_args() -> argparse.ArgumentParser:
        """Define command line interface options.

        Returns:
            argparse object containing command line interface options.
        """
        parser = argparse.ArgumentParser(
            description=(f"'mite_extras' CLI v{metadata.version('mite_extras')}."),
            formatter_class=argparse.RawTextHelpFormatter,
        )

        parser.add_argument(
            "-i",
            "--input_dir",
            type=str,
            required=True,
            help="Specifies a directory containing input files for processing.",
        )

        parser.add_argument(
            "-o",
            "--output_dir",
            type=str,
            required=True,
            help="Specifies an output directory.",
        )

        parser.add_argument(
            "-v",
            "--verboseness",
            type=str,
            default="INFO",
            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            required=False,
            help="Specifies the verboseness of logging (default: 'INFO').",
        )

        return parser
