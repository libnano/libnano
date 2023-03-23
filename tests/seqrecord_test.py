# Copyright (C) 2014-2018. Nick Conway & Ben Pruitt; Wyss Institute
# Copyright (C) 2023 Nick Conway & Ben Pruitt;
# See LICENSE.TXT for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
tests.seqrecord_test
~~~~~~~~~~~~~~~~~~~~

'''
# import json
import unittest
from os.path import (
    abspath,
    dirname,
    join,
)

import conftest

from libnano.seqrecord import (  # Feature,; Location,; SeqRecord,
    fromGenbankLike,
    locationStr2Feature,
)

LOCAL_DIR: str = abspath(dirname(__file__))


class TestSeqRecord(unittest.TestCase):

    def setUp(self):
        self.normal_files = [
            'sample.gb',
            'mds42_full.gb',
            'mds42_recode.gb',
            'failed.gb',
            'sample.gb.json',

        ]
        self.exception_files = [
            'sample_complex.gb',  # fails for compound location
        ]

    def checkFile(self, fn):
        fn_gb = join(LOCAL_DIR, 'test_data', fn)
        a = fromGenbankLike(fn_gb)
        self.assertIsNotNone(a)
        return a

    def test_checkNormalFiles(self):
        for fn in self.normal_files:
            self.checkFile(fn)

    def test_checkExceptionFiles(self):
        def doexc(fn):
            fn_gb = join(LOCAL_DIR, 'test_data', fn)
            self.assertRaises(
                NotImplementedError,
                fromGenbankLike,
                fn_gb,
            )

        for fn in self.exception_files:
            doexc(fn)

    def test_addFeature(self):
        sr = self.checkFile(self.normal_files[0])
        feature_name = 'myfeature'
        ft = locationStr2Feature(
            feature_name,
            '1615..1636',
        )
        sr.addFeature(ft)
        self.assertIn(feature_name, sr)

    def test_removeFeature(self):
        sr = self.checkFile(self.normal_files[0])
        feature_name = 'myfeature'
        ft = locationStr2Feature(
            feature_name,
            '1615..1636',
        )
        sr.addFeature(ft)
        self.assertIn(feature_name, sr)
        sr.removeFeature(ft)
        self.assertNotIn(feature_name, sr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
