"""General-purpose utility functions for the pyhgs package."""

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
    """
    try:
        # Step 1: Read the first num_head_lines from the beginning.
        first_lines = []
        with open(input_filename, 'r', encoding='utf-8') as infile:
            for i in range(num_head_lines):
                line = infile.readline()
                if not line:
                    break
                first_lines.append(line)

        # Step 2: Read the last num_tail_lines from the end.
        last_lines = []
        with open(input_filename, 'rb') as infile:
            infile.seek(0, os.SEEK_END)
            file_size = infile.tell()

            current_pos = file_size
            lines_found = 0

            while current_pos > 0 and lines_found < num_tail_lines:
                read_size = min(block_size, current_pos)
                current_pos -= read_size
                infile.seek(current_pos)

                chunk = infile.read(read_size)
                lines_in_chunk = chunk.split(b'\n')

                for line in reversed(lines_in_chunk):
                    if lines_found < num_tail_lines:
                        decoded_line = line.decode('utf-8').removesuffix('\r') + '\n'
                        last_lines.insert(0, decoded_line)
                        lines_found += 1
                    else:
                        break

                if current_pos == 0:
                    break

        # Step 3: Combine head and tail.
        lines_to_write = first_lines
        if len(first_lines) >= num_head_lines and len(last_lines) >= num_tail_lines:
            lines_to_write.append('\n...Content snipped...\n\n')
            lines_to_write.extend(last_lines)
        else:
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
