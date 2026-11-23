# -*- coding: utf-8 -*-
"""Single place where paths are resolved.

No absolute path from the original working machine survives anywhere in this
package.  Everything is derived from the location of this file, and can be
overridden with two environment variables:

    RAD_DATA   directory holding cohort.csv and radiomics.csv
    RAD_OUT    directory where intermediate and final outputs are written

Defaults are ./data and ./outputs next to this file.
"""
import os

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.abspath(os.environ.get('RAD_DATA', os.path.join(CODE_DIR, 'data')))
OUT_DIR = os.path.abspath(os.environ.get('RAD_OUT', os.path.join(CODE_DIR, 'outputs')))
FIG_DIR = os.path.join(OUT_DIR, 'figures')

# trailing separator kept so that legacy code can do OUT + 'file.json'
OUT = OUT_DIR + os.sep
DATA = DATA_DIR + os.sep

for _d in (OUT_DIR, FIG_DIR):
    os.makedirs(_d, exist_ok=True)


def data(*parts):
    """Path inside the input directory."""
    return os.path.join(DATA_DIR, *parts)


def out(*parts):
    """Path inside the output directory."""
    return os.path.join(OUT_DIR, *parts)


def fig(*parts):
    """Path inside the figure directory."""
    return os.path.join(FIG_DIR, *parts)
