#!/usr/bin/env python
"""A unit test for pyghs.parser.grok"""

import unittest
import os
import tempfile
from pathlib import Path

from hgstools.pyhgs.parser.grok import parse as grok_parse


_TEST_GENERAL='''

use domain type
    porous media

! first file
properties file
    bogus1.mprops


Initial head from file
    bogus.hen
Initial concentration from file
    bogus.cen



initial head from file some bad syntax

    don_t_find_this_file.txt



Initial head from output file                  ! a comment!!

!
! more comments about this file
! 
    boguso.pm.hen ! another comment

use domain type
    fracture
properties file


    bogus2.fprops


'''

def _write_test_file(file_path, text):
    # Write the sample text to the file.
    with open(file_path, 'w') as f:
        f.write(text)



class Test_Solute(unittest.TestCase):

    def test_solute_props(self):
        grok = grok_parse(
                os.path.dirname(os.path.abspath(__file__))+os.path.sep+
                'grok_parser_1/iel.grok')
        self.assertEqual(len(grok['solute']),1)
        self.assertTrue(abs(grok['solute']['TCE']['koc']-0.094) < 0.001)



class Test_General_Parser(unittest.TestCase):

    def test_general(self):

        with tempfile.TemporaryDirectory() as temp_dir:
            file_path = Path(temp_dir)/'test1.grok'
            _write_test_file(file_path, _TEST_GENERAL)

            grok = grok_parse(file_path)

        self.assertEqual(grok['files_mprops'], ['bogus1.mprops',])
        self.assertEqual(grok['files_fprops'], ['bogus2.fprops',])
        self.assertEqual(grok['files_initial'], ['bogus.hen', 'bogus.cen', 'boguso.pm.hen',])
    

if __name__ == '__main__':
    unittest.main()
