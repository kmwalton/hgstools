#!/usr/bin/env python
"""Extract the head and tail of the given text file.

This is particularly useful for capturing parts of a Hydrogeosphere .lst file
with the software version (head) and the simulation time and exit status (tail).
"""

import io
import os
import argparse
import sys

def excerpt_large_file(input_filename, output_filename, num_head_lines, num_tail_lines, block_size=4096):
    """
    Reads the first and last specified number of lines of a large file efficiently and
    writes them to a new file or standard output, without loading the
    entire file into memory. The last lines are found by reading from
    the end of the file in blocks.

    Args:
        input_filename (str): The path to the source file.
        output_filename (str or a file-like object): The path to the new output file, or '-' for stdout.
        num_head_lines (int): The number of lines to read from the start.
        num_tail_lines (int): The number of lines to read from the end.
        block_size (int): The size of the chunks to read from the end, in bytes.
    """
    try:
        # Step 1: Read the first 'num_head_lines' from the beginning of the file.
        first_lines = []
        with open(input_filename, 'r', encoding='utf-8') as infile:
            for i in range(num_head_lines):
                line = infile.readline()
                if not line:
                    break  # Stop if the file is smaller than num_head_lines
                first_lines.append(line)

        # Step 2: Read the last 'num_tail_lines' from the end of the file.
        last_lines = []
        with open(input_filename, 'rb') as infile:
            # Find the total size of the file
            infile.seek(0, os.SEEK_END)
            file_size = infile.tell()

            # Start position for reading backwards
            current_pos = file_size
            lines_found = 0

            # Loop backwards in blocks
            while current_pos > 0 and lines_found < num_tail_lines:
                # Determine the size of the block to read
                read_size = min(block_size, current_pos)
                current_pos -= read_size
                infile.seek(current_pos)

                # Read the block of bytes and split by newlines
                chunk = infile.read(read_size)
                lines_in_chunk = chunk.split(b'\n')

                # Process lines from the current chunk, from end to start
                # The first line is often an incomplete line segment
                for line in reversed(lines_in_chunk):
                    if lines_found < num_tail_lines:
                        # Decode the line, remove the trailing '\r' if it exists, and add a '\n'
                        decoded_line = line.decode('utf-8').removesuffix('\r') + '\n'
                        last_lines.insert(0, decoded_line)
                        lines_found += 1
                    else:
                        break

                # If we've reached the beginning of the file, we're done
                if current_pos == 0:
                    break
        
        # Step 3: Combine and prepare the lines for output.
        lines_to_write = first_lines
        # Only add ellipses if we read the full number of lines for both head and tail
        if len(first_lines) >= num_head_lines and len(last_lines) >= num_tail_lines:
            lines_to_write.append('\n...Content snipped...\n\n')
            lines_to_write.extend(last_lines)
        else:
            # If the file is small, combine the lists without the separator
            lines_to_write.extend(last_lines)

        # Step 4: Write to the specified output.
        # This conditional logic correctly handles writing to a file or stdout,
        # avoiding the use of 'with' on sys.stdout.
        if output_filename == '-':
            print(''.join(lines_to_write), end='', file=sys.stdout)
            # The success message is also written to stderr to not mix with output.
            #print(f"Successfully processed '{input_filename}' and output to stdout.", file=sys.stderr)
        elif isinstance(output_filename, io.IOBase):
            print(''.join(lines_to_write), end='', file=output_filename)
        else:
            with open(output_filename, 'w', encoding='utf-8') as outfile:
                print(''.join(lines_to_write), end='', file=outfile)
            # The success message is written to stderr to avoid mixing with stdout.
            #print(f"Successfully processed '{input_filename}' and saved output to '{output_filename}'.", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: The file '{input_filename}' was not found.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

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
    process_large_file(args.input_file, args.output, args.head_lines, args.tail_lines)

