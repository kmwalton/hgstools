#!/usr/bin/env python
"""Testing of pyghs.cli module"""

import os
import tempfile
import unittest

from pyhgs import cli


class TestCLI(unittest.TestCase):




    def test_path_to_prefix_strings(self):
        """Test parsing of strings"""
        
        with self.subTest('prefixonly'):
            pth, pfx, oth = cli.PathToPrefix.split('foo')
            self.assertEqual(pth,'.')
            self.assertEqual(pfx,'foo')
            self.assertEqual(oth,'')

        with self.subTest('path/prefix'):
            pth, pfx, oth = cli.PathToPrefix.split('bar/foo')
            self.assertEqual(pth,'bar')
            self.assertEqual(pfx,'foo')
            self.assertEqual(oth,'')

        with self.subTest('path/prefixsuffix'):
            pth, pfx, oth = cli.PathToPrefix.split('bar/fooo.lst')
            self.assertEqual(pth,'bar')
            self.assertEqual(pfx,'foo')
            self.assertEqual(oth,'o.lst')

        with self.subTest('path/prefixsuffix'):
            pth, pfx, oth = cli.PathToPrefix.split('bar/foo.lst')
            self.assertEqual(pth,'bar')
            self.assertEqual(pfx,'fo')
            self.assertEqual(oth,'o.lst')

        with self.subTest('prefixsuffix'):
            pth, pfx, oth = cli.PathToPrefix.split('fooo.lst')
            self.assertEqual(pth,'.')
            self.assertEqual(pfx,'foo')
            self.assertEqual(oth,'o.lst')

        with self.subTest('suffix'):
            pth, pfx, oth = cli.PathToPrefix.split('bar.txt')
            self.assertEqual(pth,'.')
            self.assertEqual(pfx,'')
            self.assertEqual(oth,'bar.txt')


    def test_path_to_prefix_batchpfx(self):
        """Test parsing of strings when batch.pfx is present"""
        owd = os.getcwd()

        with tempfile.TemporaryDirectory(dir='.', prefix='unittest') as tdir:
            os.chdir(tdir)

            with open('batch.pfx','w') as batchpfxfout:
                print('foo', file=batchpfxfout)


            with self.subTest('nothing'):
                pth, pfx, oth = cli.PathToPrefix.split('')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'foo')
                self.assertEqual(oth,'')

            with self.subTest('dot'):
                pth, pfx, oth = cli.PathToPrefix.split('.')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'foo')
                self.assertEqual(oth,'')

            with self.subTest('prefix'):
                pth, pfx, oth = cli.PathToPrefix.split('foo')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'foo')
                self.assertEqual(oth,'')

            with self.subTest('allother'):
                pth, pfx, oth = cli.PathToPrefix.split('bar')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'foo')
                self.assertEqual(oth,'bar')

            with self.subTest('prefixother'):
                pth, pfx, oth = cli.PathToPrefix.split('foo.blahblah')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'foo')
                self.assertEqual(oth,'.blahblah')


        os.chdir(owd)


    def test_argparse_action(self):
        """Test use of this as an argument action"""
        pass
        

if __name__ == '__main__':
    unittest.main()
