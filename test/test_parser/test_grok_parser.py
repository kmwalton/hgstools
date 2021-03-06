#!/usr/bin/env python
"""A unit test for pyghs.parser.grok"""

import unittest
import os

import pyhgs.parser.grok


class Test_Solute(unittest.TestCase):

    def test_solute_props(self):
        grok = pyhgs.parser.grok.parse(
                os.path.dirname(os.path.abspath(__file__))+os.path.sep+
                'grok_parser_1/iel.grok')
        self.assertEqual(len(grok['solute']),1)
        self.assertTrue(abs(grok['solute']['TCE']['koc']-0.094) < 0.001)


if __name__ == '__main__':
    unittest.main()
