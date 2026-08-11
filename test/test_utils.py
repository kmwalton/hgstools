"""Tests for hgstools.pyhgs.utils"""

import io
import os
import tempfile
import unittest
from pathlib import Path

from hgstools.pyhgs.utils import excerpt_large_file

SNIP = '\n...Content snipped...\n\n'


def _excerpt(text, head, tail, **kwargs):
    """Return the excerpt of `text` as a string."""
    with tempfile.TemporaryDirectory() as d:
        fn = Path(d) / 'in.txt'
        fn.write_text(text, encoding='utf-8', newline='')
        out = io.StringIO()
        excerpt_large_file(str(fn), out, head, tail, **kwargs)
        return out.getvalue()


def _numbered(n, eol='\n'):
    return ''.join(f'line {i}{eol}' for i in range(1, n + 1))


class TestExcerptLargeFile(unittest.TestCase):

    def test_short_file_reproduced_in_full(self):
        """A file shorter than head+tail must not repeat any line."""
        txt = _numbered(30)
        self.assertEqual(_excerpt(txt, 20, 21), txt)

    def test_no_duplicated_lines_for_any_short_length(self):
        for n in range(0, 45):
            with self.subTest(n=n):
                got = _excerpt(_numbered(n), 20, 21)
                lines = [l for l in got.splitlines() if l.startswith('line ')]
                self.assertEqual(len(lines), len(set(lines)),
                    f'duplicate lines in excerpt of a {n}-line file')

    def test_exactly_head_plus_tail_lines_is_full_file(self):
        txt = _numbered(41)
        self.assertEqual(_excerpt(txt, 20, 21), txt)
        self.assertNotIn(SNIP, _excerpt(txt, 20, 21))

    def test_one_line_beyond_gets_snipped(self):
        got = _excerpt(_numbered(42), 20, 21)
        self.assertIn(SNIP, got)
        self.assertEqual(got,
            _numbered(20) + SNIP
            + ''.join(f'line {i}\n' for i in range(22, 43)))

    def test_long_file(self):
        got = _excerpt(_numbered(1000), 3, 2)
        self.assertEqual(got,
            'line 1\nline 2\nline 3\n' + SNIP + 'line 999\nline 1000\n')

    def test_tail_spans_multiple_blocks(self):
        """Small block_size forces the backwards reader to loop."""
        got = _excerpt(_numbered(1000), 3, 20, block_size=16)
        self.assertEqual(got,
            'line 1\nline 2\nline 3\n' + SNIP
            + ''.join(f'line {i}\n' for i in range(981, 1001)))

    def test_no_trailing_newline(self):
        txt = _numbered(1000).rstrip('\n')
        got = _excerpt(txt, 2, 2)
        self.assertEqual(got,
            'line 1\nline 2\n' + SNIP + 'line 999\nline 1000\n')

    def test_crlf_line_endings(self):
        got = _excerpt(_numbered(1000, eol='\r\n'), 2, 2)
        self.assertEqual(got,
            'line 1\nline 2\n' + SNIP + 'line 999\nline 1000\n')

    def test_form_feed_is_not_a_line_break(self):
        """HGS listing files contain form feeds; they must not split lines."""
        txt = _numbered(500) + ''.join(f'\x0cpage {i}\n' for i in range(1, 6))
        got = _excerpt(txt, 2, 3)
        self.assertEqual(got,
            'line 1\nline 2\n' + SNIP + '\x0cpage 3\n\x0cpage 4\n\x0cpage 5\n')

    def test_empty_file(self):
        self.assertEqual(_excerpt('', 20, 21), '')

    def test_zero_tail_lines(self):
        got = _excerpt(_numbered(1000), 3, 0)
        self.assertEqual(got, 'line 1\nline 2\nline 3\n' + SNIP)


if __name__ == '__main__':
    unittest.main()
