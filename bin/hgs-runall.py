#!/usr/bin/env python
"""Launches the extended Hydrogeosphere toolchain.

Basic Toolchain:
1) grok
2) phgs
3) hgs2vtu (use --use-hsplot to swap to hsplot)

Extensions include:
0) preprocess* | sort
1..3)
4) postprocess* | sort

Note: this program captures the standard output and errors streams. This may
interfere with a debugger.
"""

import sys
import argparse

from hgstools.pyhgs.runner import HGSToolChainRun


def main():
    """Main function to parse arguments and execute the toolchain."""
    argp = argparse.ArgumentParser(
        description='Launches the extended Hydrogeosphere toolchain (preprocess, grok, phgs, [hgs2vtu|hsplot], postprocess). hgs2vtu is used by default.',
        formatter_class=argparse.RawTextHelpFormatter  # Allows custom formatting in help text
    )

    argp.add_argument(
        '-d', '--sim-dir',
        default='.',
        type=str,
        help='Directory containing simulation inputs (e.g., input files for grok and phgs).',
    )

    # NEW ARGUMENT to allow swapping the plotting tool
    argp.add_argument(
        '--use-hsplot',
        action='store_true',
        help="Use 'hsplot' for Step 3 instead of the default 'hgs2vtu'."
    )

    argp.add_argument(
        '--dry-run',
        action='store_true',
        help=(
            'List the steps in the toolchain and verify executability of each tool.\n'
            'Also checks file read and write permissions in the simulation directory.'
        )
    )

    argp.add_argument(
        '-sw', '--start-with',
        type=str,
        default=None,
        help=(
            "Start the toolchain with tool N (e.g., '2' or 'phgs').\n"
            "The chain begins with the matching tool with the lowest index."
        )
    )

    argp.add_argument(
        '-ew', '--end-with',
        type=str,
        default=None,
        help=(
            "End the toolchain after tool N finishes (e.g., 'hgs2vtu' or '4').\n"
            "The chain ends after the matching tool with the highest index."
        )
    )

    args = argp.parse_args()

    # The logic relies on the imported class HGSToolChainRun
    tc = HGSToolChainRun(args.sim_dir)

    # Logic to override the default tool for step 3 based on the new flag
    if args.use_hsplot:
        # NOTE: This assumes HGSToolChainRun has a method to override a tool name
        # with a new one (e.g., swapping 'hgs2vtu' with 'hsplot').
        tc.set_tool_override('hgs2vtu', 'hsplot')

    tc.set_start(args.start_with)
    tc.set_end(args.end_with)

    if args.dry_run:
        # Print the toolchain structure (assuming tc.__str__ is implemented)
        print(tc)
        (returncode, msgs) = tc.check_tools()

        if returncode != 0:
            # Use f-string and standard print arguments for cleaner error output
            print(
                f'Error with one or more executables in the toolchain:',
                '\n'.join(msgs),
                file=sys.stderr
            )
        sys.exit(returncode)

    returncode = tc.run()

    # The final exit is clean and simple.
    sys.exit(returncode)

if __name__ == '__main__':
    # No need for os or glob in this specific entry point
    main()

