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

import logging
import sys
from importlib import metadata

import coloredlogs
from mite_schema import SchemaManager

from mite_extras import CliManager, FileManager, MiteParser


def config_logger(verboseness: str) -> logging.Logger:
    """Set up a named logger with nice formatting

    Args:
        verboseness: sets the logging verboseness

    Returns:
        A Logger object
    """
    logger = logging.getLogger("mite_extras")
    logger.setLevel(getattr(logging, verboseness))
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
    args = CliManager().run(sys.argv[1:])

    logger = config_logger(args.verboseness)
    logger.debug(f"Started 'mite_extras' v{metadata.version('mite_extras')} as CLI.")

    schema_manager = SchemaManager()

    file_manager = FileManager(indir=args.input_dir, outdir=args.output_dir)
    file_manager.read_files_indir()

    for entry in file_manager.infiles:
        logger.info(f"CLI: started parsing of file '{entry.name}'.")

        input_data = schema_manager.read_json(infile=entry)

        try:
            parser = MiteParser()
            parser.parse_mite_json(data=input_data)

            schema_manager.validate_mite(instance=parser.to_json())

            file_manager.write_json(outfile_name=entry.stem, payload=parser.to_json())

            logger.info(f"CLI: completed parsing of file '{entry.name}'.")
        except Exception as e:
            logger.fatal(f"Could not process file '{entry.name}': {e!s}")
            continue

    logger.info("Completed 'mite_extras' as CLI.")


if __name__ == "__main__":
    main_cli()
