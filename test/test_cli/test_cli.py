#!/usr/bin/env python
"""Testing of pyghs.cli module"""

import os
import tempfile
import argparse
import unittest

from hgstools.pyhgs import cli


class TestCLI(unittest.TestCase):

    def test_path_to_prefix_alternative_methods(self):
        """Test aliases to main parsing method"""
        
        pth, pfx, oth = cli.PathToPrefix.split('foo')
        pth, pfx, oth = cli.parse_path_to_prefix('foo')
 
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
                # Change in functionality here -- If the user passes a prefix
                # that does not match what's in batch.pfx, then assume the user
                # intended this to be the prefix
                pth, pfx, oth = cli.PathToPrefix.split('bar')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'bar')
                self.assertEqual(oth,'')

            with self.subTest('prefixother'):
                pth, pfx, oth = cli.PathToPrefix.split('foo.blahblah')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'foo')
                self.assertEqual(oth,'.blahblah')

                pth, pfx, oth = cli.PathToPrefix.split('foo.bar.blahblah')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'foo')
                self.assertEqual(oth,'.bar.blahblah')

                pth, pfx, oth = cli.PathToPrefix.split('fooo.bar.blahblah')
                self.assertEqual(pth,'.')
                self.assertEqual(pfx,'foo')
                self.assertEqual(oth,'o.bar.blahblah')

            os.chdir(owd)


    def test_argparse_action(self):
        """Test use of this as an argparse.Action derived class"""
        argp = argparse.ArgumentParser()
        argp.add_argument('PATH_TO_PREFIX',
            metavar='path/to/prefix',
            action=cli.PathToPrefix,
            help='path/to/prefix of the HGS simulation',
            )
        args = argp.parse_args(['some_path/prefixo.lst',])

        pth, filename = os.path.split(args.PATH_TO_PREFIX)
        self.assertEqual(pth,'some_path')
        self.assertEqual(filename,'prefixo.lst')

        pth, pfx, oth = cli.parse_path_to_prefix(args.PATH_TO_PREFIX)
        self.assertEqual(pth,'some_path')
        self.assertEqual(pfx,'prefix')
        self.assertEqual(oth,'o.lst')
        

if __name__ == '__main__':
    unittest.main()
