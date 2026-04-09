#!/usr/bin/env python
"""Extract the head and tail of the given text file.

This is particularly useful for capturing parts of a Hydrogeosphere .lst file
with the software version (head) and the simulation time and exit status (tail).
"""

import argparse
import sys

from hgstools.pyhgs.utils import excerpt_large_file

# Example usage:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read the first and last lines of a large file.')
    parser.add_argument('input_file', type=str, help='Path to the input file.')
    parser.add_argument('-o', '--output', type=str, default='-',
                        help='Path to the output file. Use "-" for standard output (default).')
    parser.add_argument('-n', '--head-lines', type=int, default=20,
                        help='Number of lines to read from the start of the file (default: 20).')
    parser.add_argument('-m', '--tail-lines', type=int, default=21,
                        help='Number of lines to read from the end of the file (default: 21).')

    args = parser.parse_args()
    excerpt_large_file(args.input_file, args.output, args.head_lines, args.tail_lines)

