"""Main entry point to mite_extras.

Copyright (c) 2024 to present Mitja M. Zdouc and individual contributors.

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

import json
import logging
import sys
from importlib import metadata

import coloredlogs

from mite_extras import CliManager, FileManager, Parser, SchemaManager


def config_logger() -> logging.Logger:
    """Set up a named logger with nice formatting

    Returns:
        A Logger object
    """
    logger = logging.getLogger("mite_extras")
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(
        coloredlogs.ColoredFormatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
    )
    logger.addHandler(console_handler)
    return logger


def main_cli() -> None:
    """Entry point for CLI"""
    logger = config_logger()
    logger.debug(f"Started 'mite_extras' v{metadata.version('mite_extras')} as CLI.")

    args = CliManager().run(sys.argv[1:])
    file_manager = FileManager(indir=args.input_dir, outdir=args.output_dir)
    file_manager.read_files_indir()

    for entry in file_manager.infiles:
        with open(entry) as infile:
            input_data = json.load(infile)

        output_data = {}

        try:
            match args.format:
                case "raw":
                    output_data = Parser().parse_raw_json(entry.name, input_data)
                case "mite":
                    output_data = Parser().parse_mite_json(input_data)
                case _:
                    raise RuntimeError(
                        f"Input format '{args.format}' could not be parsed."
                    )
            SchemaManager().validate_against_schema(output_data)
        except Exception as e:
            logger.fatal(f"Could not process file '{entry.name}': {e!s}")
            continue

        file_manager.write_to_outdir(outfile_name=entry.name, payload=output_data)

    logger.debug("Completed 'mite_extras' as CLI.")


if __name__ == "__main__":
    main_cli()
