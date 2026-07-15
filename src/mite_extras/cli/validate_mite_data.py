"""CLI helper to validate MITE JSON entries and produce a summary report.

Exposed as console script `validate-mite-data` so it can be run via `uv run validate-mite-data` if uv runs console scripts.

Usage (from repo root with venv):
PYTHONPATH=src .venv/bin/python3 -m validate_mite_data --data-dir ~/Git/mite_data/mite_data/data --out ~/mite_validation_summary.json

Or as console script after installation:
validate-mite-data --data-dir ~/Git/mite_data/mite_data/data --out ~/mite_validation_summary.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from rdkit import RDLogger

from mite_extras.processing.validation_manager import ReactionValidator

RDLogger.DisableLog('rdApp.*')


def validate_all(data_dir: str, out_path: str) -> str:
    rv = ReactionValidator()
    files = sorted(glob.glob(os.path.join(os.path.expanduser(data_dir), '*.json')))
    summary = {}
    for f in files:
        try:
            j = json.load(open(f))
        except Exception as e:
            summary[os.path.basename(f)] = [('file_read_error', str(e))]
            continue
        res = []
        for i, block in enumerate(j.get('reactions', []), 1):
            smarts = block.get('reactionSMARTS')
            for jdx, r in enumerate(block.get('reactions', []), 1):
                substrate = r.get('substrate')
                expected = r.get('products', [])
                key = f'{i}.{jdx}'
                try:
                    rv.validate_reaction(smarts, substrate, expected)
                    res.append((key, 'OK'))
                except Exception as e:
                    msg = str(e).lower()
                    if 'stereochemistry' in msg or 'non-isomeric' in msg:
                        res.append((key, 'OK_noniso_match'))
                    elif 'kekul' in msg or "can't kekulize" in msg:
                        res.append((key, 'KEKULIZE_ERROR'))
                    elif 'explicit valence' in msg:
                        res.append((key, 'VALENCE_ERROR'))
                    else:
                        res.append((key, 'OTHER_ERROR', str(e)))
        summary[os.path.basename(f)] = res
    with open(os.path.expanduser(out_path), 'w') as fh:
        json.dump(summary, fh, indent=2)
    return os.path.expanduser(out_path)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description='Validate MITE JSON entries and write a summary report')
    p.add_argument('--data-dir', default='~/Git/mite_data/mite_data/data', help='Path to directory containing MITE JSON files')
    p.add_argument('--out', default='~/mite_validation_summary.json', help='Output JSON report path')
    args = p.parse_args(argv)

    out = validate_all(args.data_dir, args.out)
    print('WROTE', out)


if __name__ == '__main__':
    main()
