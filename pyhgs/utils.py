"""General-purpose utility functions for the pyhgs package."""

__docformat__ = 'numpy'

import os
import sys
import io


def excerpt_large_file(input_filename, output_filename, num_head_lines, num_tail_lines, block_size=4096):
    """Read the first and last lines of a large file and write them to output.

    Reads efficiently without loading the entire file into memory. The tail is
    found by reading backwards from the end of the file in blocks.

    Parameters
    ----------
    input_filename : str
        Path to the source file.
    output_filename : str or file-like object
        Path to the output file, '-' for stdout, or an open file-like object.
    num_head_lines : int
        Number of lines to read from the start of the file.
    num_tail_lines : int
        Number of lines to read from the end of the file.
    block_size : int, optional
        Size in bytes of chunks read from the end of the file, by default 4096.

    Notes
    -----
    Files short enough that the head and tail would overlap are reproduced in
    full; no line is ever emitted twice.
    """
    num_head_lines = max(0, num_head_lines)
    num_tail_lines = max(0, num_tail_lines)

    try:
        # Step 1: Read from the beginning, one line further than the excerpt
        # needs. If the file ends within that many lines, the head and the tail
        # would overlap (or abut), so the whole file is the excerpt.
        first_lines = []
        with open(input_filename, 'r', encoding='utf-8') as infile:
            for _ in range(num_head_lines + num_tail_lines + 1):
                line = infile.readline()
                if not line:
                    break
                first_lines.append(line)

        if len(first_lines) <= num_head_lines + num_tail_lines:
            lines_to_write = first_lines

        else:
            del first_lines[num_head_lines:]

            # Step 2: Read the last num_tail_lines from the end, working
            # backwards in blocks until one more line break than needed is
            # buffered -- the extra break guarantees the earliest line kept is
            # complete rather than a fragment of a longer line.
            buf = b''
            with open(input_filename, 'rb') as infile:
                infile.seek(0, os.SEEK_END)
                current_pos = infile.tell()

                while current_pos > 0 and buf.count(b'\n') <= num_tail_lines:
                    read_size = min(block_size, current_pos)
                    current_pos -= read_size
                    infile.seek(current_pos)
                    buf = infile.read(read_size) + buf

            # split on b'\n' only: str.splitlines() would also break on form
            # feeds and other characters that occur in HGS listing files
            tail = buf.decode('utf-8', errors='replace').split('\n')
            if tail and tail[-1] == '':
                tail.pop() # trailing newline; not a line of its own
            last_lines = [ l.removesuffix('\r') + '\n'
                for l in (tail[-num_tail_lines:] if num_tail_lines else []) ]

            # Step 3: Combine head and tail.
            lines_to_write = first_lines
            lines_to_write.append('\n...Content snipped...\n\n')
            lines_to_write.extend(last_lines)

        # Step 4: Write to the specified output.
        combined = ''.join(lines_to_write)
        if output_filename == '-':
            print(combined, end='', file=sys.stdout)
        elif isinstance(output_filename, io.IOBase):
            print(combined, end='', file=output_filename)
        else:
            with open(output_filename, 'w', encoding='utf-8') as outfile:
                print(combined, end='', file=outfile)

    except FileNotFoundError:
        print(f"Error: The file '{input_filename}' was not found.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
