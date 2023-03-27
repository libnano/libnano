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
tests.test_gb_writer
~~~~~~~~~~~~~~~~~~~~

'''
import filecmp
# import json
import os
import os.path as op
import unittest
from tempfile import NamedTemporaryFile

from libnano.fileio import (
    gb_reader,
    gb_writer,
)

# import conftest


LOCAL_DIR: str = op.abspath(op.dirname(__file__))


class TestGBWriter(unittest.TestCase):

    def setUp(self):
        self.good_files = [
            'sample.gb',
            'mds42_full.gb',
            'mds42_recode.gb',
            'sample_complex.gb',
        ]
        self.bad_files = [
            'failed.gb',  # fails due to ApE double COMMENT key
            # biopython has 80 character lines and we do 79
            'mds42_recode_biopython.gb',
        ]
        self.maxDiff = None     # assertEquals will print out whole diff

    def check_file(self, fn, should_be_true):
        fn_gb = op.join(LOCAL_DIR, 'test_data', fn)
        d_gb = gb_reader.parse(fn_gb, is_ordered=True)

        f_temp = NamedTemporaryFile(
            mode='w',
            delete=False,
            encoding='utf-8',
        )

        gb_writer.write(f_temp, d_gb)
        f_temp.close()
        if should_be_true:
            self.assertTrue(
                filecmp.cmp(fn_gb, f_temp.name),
            )
        else:
            self.assertFalse(
                filecmp.cmp(fn_gb, f_temp.name),
            )
        os.unlink(f_temp.name)
        self.assertFalse(
            op.exists(f_temp.name),
        )

    def test_good_files(self):
        for fn in self.good_files:
            self.check_file(fn, True)

    def test_bad_files(self):
        for fn in self.bad_files:
            self.check_file(fn, False)


if __name__ == '__main__':
    unittest.main(verbosity=2)
