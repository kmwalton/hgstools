# hgstools

A collection of tools that operate adjacent to Aquanty Inc.'s HydroGeoSphere.

This collection is not affiliated or endorsed by Aquanty. It is provided as-is
with no warranty or claims of fitness for purpose. See LICENSE.

This is the author's collection of scripts, tools, and supporting code to simplify
the HGS work flow, create input files, parse/modify output files, and perform other
analyses. The workflow is mainly centred on HGS' capability as a discrete fracture-
matrix, sub-surface domain, flow and contaminant transport simulator; other aspects
of HGS' functionality are not specifically considered or tested for these tools.

USAGE:
  1. clone the repository to a directory of your choosing, like `hgstools`
  2. add `hgstools/bin` to your system's `$PATH`, `$env:PATH`, or `%PATH%` to use scripts/executables
  3. add `hgstools` to your system's `$PYTHONPATH` for use as a python module

WHAT'S HERE:

Scripts in `bin/`
  1. `hgs-copy-inputs.py` - Copy inputs and adjacent to a different directory.
  2. `hgs-nnnn-to-en.py` - Convert an arbitrary binary output file to a restart
     file.

Most of the above scripts can be run with argument `--help` for full
documentation.

Used/tested with Python3.8 and on Mac and Windows computers.
