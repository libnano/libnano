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
tests.gb_reader_test
~~~~~~~~~~~~~~~~~~~~

'''
import io
import json
import os.path as op
import unittest

import conftest

from libnano.fileio import gb_reader

LOCAL_DIR: str = op.abspath(op.dirname(__file__))


class TestGBReader(unittest.TestCase):

    def setUp(self):
        self.files = [
            'sample.gb',
            'mds42_full.gb',
            'mds42_recode.gb',
            'sample_complex.gb',
            'failed.gb',

        ]

    def check_file(self, fn, should_fail=False):
        fn_gb = op.join(LOCAL_DIR, 'test_data', fn)
        d_gb = gb_reader.parse(fn_gb)
        fn_json = fn_gb + '.json'
        with io.open(fn_json, 'r', encoding='utf-8') as fd_json:
            d_json = json.load(fd_json)
        if should_fail:
            self.assertNotEqual(d_gb['info'], d_json['info'])
        else:
            self.assertEqual(d_gb['info'], d_json['info'])

    def test_check_files(self):
        for fn in self.files:
            self.check_file(fn, should_fail=('fail' in fn))
# end class


if __name__ == '__main__':
    unittest.main(verbosity=2)
