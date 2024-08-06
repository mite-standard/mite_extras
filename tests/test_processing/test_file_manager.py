import os
from pathlib import Path

import pytest
from mite_extras.processing.file_manager import FileManager


def test_init_valid():
    assert isinstance(
        FileManager(
            indir=Path(__file__).parent.joinpath("example_indir_mite"),
            outdir=Path(__file__).parent.joinpath("example_outdir"),
        ),
        FileManager,
    )


def test_init_invalid():
    with pytest.raises(FileNotFoundError):
        FileManager(
            indir=Path(__file__).parent.joinpath("nonexistent"),
            outdir=Path(__file__).parent.joinpath("nonexistent"),
        )


def test_read_files_indir_valid():
    instance = FileManager(
        indir=Path(__file__).parent.joinpath("example_indir_mite"),
        outdir=Path(__file__).parent.joinpath("example_outdir"),
    )
    instance.read_files_indir()
    assert len(instance.infiles) == 1


def test_read_files_indir_invalid():
    instance = FileManager(
        outdir=Path(__file__).parent.joinpath("example_indir_mite"),
        indir=Path(__file__).parent.joinpath("example_outdir"),
    )
    instance.read_files_indir()
    assert len(instance.infiles) == 0


def test_write_json_valid():
    instance = FileManager(
        indir=Path(__file__).parent.joinpath("example_indir_mite"),
        outdir=Path(__file__).parent.joinpath("example_outdir"),
    )
    instance.write_json(outfile_name="testfile", payload={})
    assert instance.outdir.joinpath("testfile.json").exists()
    os.remove(path=instance.outdir.joinpath("testfile.json"))
