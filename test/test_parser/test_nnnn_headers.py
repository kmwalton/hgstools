"""Tests for 80-character magic string header parsing in HGS binary .NNNN files.

Covers both "modern" HGS output (r2853+, July 2025 — header contains a
timestamp token plus a structured version/magic token starting with ``#!HGS``)
and "legacy" HGS output (header contains only a timestamp token).

Header byte layout (new format, 1-based within the 80-char record):
  Bytes  1–55  simulation timestamp + spaces
  Bytes 56–60  magic ``#!HGS``
  Bytes 61–63  file format version
  Bytes 64–65  endianness (LE/BE/UK)
  Bytes 66–68  data type (I04/R04/R08/UNK)
  Bytes 69–70  number of components (zero-padded)
  Bytes 71–80  number of values (zero-padded)

Test data:
  modern: test_sims/04b_Saturated_Fracture_Transport   (committed to repo)
  legacy: test_sims/04b_Saturated_Fracture_Transport_lega  (NOT committed)

Recreating the legacy directory
--------------------------------
The ``_lega`` directory must be generated manually and is not version-controlled.
To recreate it, run the 04b_Saturated_Fracture_Transport simulation with an HGS
executable older than r2853 (pre-July 2025) and place (or symlink) the output
directory alongside the modern one, named ``04b_Saturated_Fracture_Transport_lega``.

Any test that depends solely on the legacy data is decorated with
``@unittest.skipUnless(_has_legacy, ...)`` or ``@unittest.skipUnless(_has_legacy_head, ...)``
and will be skipped automatically when the directory is absent.
"""

import struct
import unittest
from pathlib import Path

from hgstools.pyhgs.parser import parse as hgs_parse, peek_NNNN_time
from hgstools.pyhgs.test import sims_join

_MODERN_DIR  = Path(sims_join('04b_Saturated_Fracture_Transport'))
_LEGACY_DIR  = Path(sims_join('04b_Saturated_Fracture_Transport_lega'))

_MODERN_CONC = _MODERN_DIR / 'module4bo.conc_pm.salt.0001'
_LEGACY_CONC = _LEGACY_DIR / 'module4bo.conc_pm.salt.0001'

_MODERN_HEAD = _MODERN_DIR / 'module4bo.head_pm.0001'
_LEGACY_HEAD = _LEGACY_DIR / 'module4bo.head_pm.0001'

_has_modern      = _MODERN_CONC.exists()
_has_legacy      = _LEGACY_CONC.exists()
_has_modern_head = _MODERN_HEAD.exists()
_has_legacy_head = _LEGACY_HEAD.exists()

# Output step 1 in both datasets is t=86400 s (1 day)
_EXPECTED_TIME_STR        = '86400.0000000000'
_EXPECTED_TIME_STR_LEGACY = '8.6400000000E+04'

_MODERN_VER_PREFIX = '#!HGS001'


def _raw_header_tokens(path):
    """Return the stripped split tokens from the 80-byte Fortran header record."""
    with open(path, 'rb') as f:
        raw = f.read(88)
    rec_len = struct.unpack('<i', raw[:4])[0]
    header = raw[4:4 + rec_len].decode('utf-8', errors='replace')
    return header.strip().split()


class TestRawHeaderStructure(unittest.TestCase):
    """Low-level checks: correct Fortran record length and token count."""

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_record_length(self):
        with open(_MODERN_CONC, 'rb') as f:
            rec_len = struct.unpack('<i', f.read(4))[0]
        self.assertEqual(rec_len, 80)

    @unittest.skipUnless(_has_legacy, 'Legacy test data not found')
    def test_legacy_record_length(self):
        with open(_LEGACY_CONC, 'rb') as f:
            rec_len = struct.unpack('<i', f.read(4))[0]
        self.assertEqual(rec_len, 80)

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_header_has_two_tokens(self):
        tokens = _raw_header_tokens(_MODERN_CONC)
        self.assertEqual(len(tokens), 2,
            f'Expected 2 header tokens (timestamp + version), got {tokens!r}')

    @unittest.skipUnless(_has_legacy, 'Legacy test data not found')
    def test_legacy_header_has_one_token(self):
        tokens = _raw_header_tokens(_LEGACY_CONC)
        self.assertEqual(len(tokens), 1,
            f'Expected 1 header token (timestamp only), got {tokens!r}')

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_version_token_prefix(self):
        tokens = _raw_header_tokens(_MODERN_CONC)
        self.assertTrue(tokens[1].startswith(_MODERN_VER_PREFIX),
            f'Version token {tokens[1]!r} does not start with {_MODERN_VER_PREFIX!r}')


class TestPeekNNNNTime(unittest.TestCase):
    """Tests for the ``peek_NNNN_time`` helper."""

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_conc_time(self):
        t = peek_NNNN_time(str(_MODERN_CONC))
        self.assertEqual(t, _EXPECTED_TIME_STR)

    @unittest.skipUnless(_has_legacy, 'Legacy test data not found')
    def test_legacy_conc_time(self):
        t = peek_NNNN_time(str(_LEGACY_CONC))
        self.assertEqual(t, _EXPECTED_TIME_STR_LEGACY)

    @unittest.skipUnless(_has_modern_head, 'Modern head_pm test data not found')
    def test_modern_head_time(self):
        t = peek_NNNN_time(str(_MODERN_HEAD))
        self.assertIsNotNone(t)
        self.assertGreater(len(t), 0)

    @unittest.skipUnless(_has_legacy_head, 'Legacy head_pm test data not found')
    def test_legacy_head_time(self):
        t = peek_NNNN_time(str(_LEGACY_HEAD))
        self.assertIsNotNone(t)
        self.assertGreater(len(t), 0)


class TestParseHeaderFields(unittest.TestCase):
    """Tests that ``parse()`` extracts ``ts`` and ``ver`` correctly from headers."""

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_conc_ts(self):
        d = hgs_parse(str(_MODERN_CONC))
        self.assertEqual(d['ts'], _EXPECTED_TIME_STR)

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_conc_ver_present(self):
        d = hgs_parse(str(_MODERN_CONC))
        self.assertIsNotNone(d['ver'],
            'Modern files should carry a non-None version field')
        self.assertTrue(d['ver'].startswith(_MODERN_VER_PREFIX),
            f"ver={d['ver']!r} should start with {_MODERN_VER_PREFIX!r}")

    @unittest.skipUnless(_has_legacy, 'Legacy test data not found')
    def test_legacy_conc_ts(self):
        d = hgs_parse(str(_LEGACY_CONC))
        self.assertEqual(d['ts'], _EXPECTED_TIME_STR_LEGACY)

    @unittest.skipUnless(_has_legacy, 'Legacy test data not found')
    def test_legacy_conc_ver_is_none(self):
        d = hgs_parse(str(_LEGACY_CONC))
        self.assertIsNone(d['ver'],
            'Legacy files carry no version token; ver should be None')

    @unittest.skipUnless(_has_modern_head, 'Modern head_pm test data not found')
    def test_modern_head_ts(self):
        d = hgs_parse(str(_MODERN_HEAD))
        self.assertIsNotNone(d.get('ts'))

    @unittest.skipUnless(_has_legacy_head, 'Legacy head_pm test data not found')
    def test_legacy_head_ts(self):
        d = hgs_parse(str(_LEGACY_HEAD))
        self.assertIsNotNone(d.get('ts'))

    @unittest.skipUnless(_has_modern and _has_legacy, 'Both datasets required')
    def test_ts_represents_same_time(self):
        """Both formats should encode the same simulation time as a float."""
        modern_t = float(hgs_parse(str(_MODERN_CONC))['ts'])
        legacy_t = float(hgs_parse(str(_LEGACY_CONC))['ts'])
        self.assertAlmostEqual(modern_t, legacy_t, places=6)


class TestDecodedHeaderFields(unittest.TestCase):
    """Tests for the structured sub-fields decoded from new-format (r2853+) headers."""

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_conc_hgs_file_ver(self):
        d = hgs_parse(str(_MODERN_CONC))
        self.assertEqual(d['hgs_file_ver'], '001')

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_conc_endian(self):
        d = hgs_parse(str(_MODERN_CONC))
        self.assertIn(d['endian'], ('LE', 'BE', 'UK'),
            f"endian={d['endian']!r} is not a recognised code")
        self.assertEqual(d['endian'], 'LE')

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_conc_data_type(self):
        # conc_pm is real-8 (float64)
        d = hgs_parse(str(_MODERN_CONC))
        self.assertEqual(d['data_type'], 'R08')

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_conc_n_components(self):
        # concentration is a scalar field
        d = hgs_parse(str(_MODERN_CONC))
        self.assertEqual(d['n_components'], 1)

    @unittest.skipUnless(_has_modern, 'Modern test data not found')
    def test_modern_conc_n_values_matches_data(self):
        d = hgs_parse(str(_MODERN_CONC))
        self.assertEqual(d['n_values'], len(d['data']),
            f"Header n_values={d['n_values']} != actual array length {len(d['data'])}")

    @unittest.skipUnless(_has_modern_head, 'Modern head_pm test data not found')
    def test_modern_head_data_type(self):
        # head_pm is also real-8
        d = hgs_parse(str(_MODERN_HEAD))
        self.assertEqual(d['data_type'], 'R08')

    @unittest.skipUnless(_has_modern_head, 'Modern head_pm test data not found')
    def test_modern_head_n_values_matches_data(self):
        d = hgs_parse(str(_MODERN_HEAD))
        self.assertEqual(d['n_values'], len(d['data']))

    @unittest.skipUnless(_has_legacy, 'Legacy test data not found')
    def test_legacy_decoded_fields_absent(self):
        """Legacy files must not expose the new-format sub-field keys."""
        d = hgs_parse(str(_LEGACY_CONC))
        for key in ('hgs_file_ver', 'endian', 'data_type', 'n_components', 'n_values'):
            self.assertNotIn(key, d,
                f"Legacy parse result should not contain key {key!r}")


if __name__ == '__main__':
    unittest.main()
