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
tests.test_seedfinder
~~~~~~~~~~~~~~~~~~~~~

'''
from __future__ import print_function

import random
# import re
# import time
import unittest

from libnano.search import seedfinder  # type: ignore

# import conftest


'''TODO: Fix this file for the submer module
This also tests
'''


class TestSeedFinder(unittest.TestCase):

    def test_checkSeed(self):
        self.assertEqual(
            seedfinder.checkSeed('#-##--#-##', 15),
            2,
        )
        self.assertEqual(
            seedfinder.checkSeed('#-##-##-##', 15),
            1,
        )
        self.assertNotEqual(
            seedfinder.checkSeed('#######-##', 15),
            2,
        )
        # self.assertEqual(seedfinder.checkSeeds(
        #     ["#-#-#---#-----#-#-#---#-----#-#-#---#", "###-#--###-#--###-#"],
        #      50, k_min=4, k_max=6), 5)

    def test_findSeed(self):
        for x in range(10):
            m = random.randint(5, 20)
            k = random.randint(1, 3)
            sd_list = seedfinder.findSeed(m, k)
            # sds = seedfinder.findSeed(m, k)
            for sd in sd_list:
                the_seed = sd[2]
                self.assertEqual(
                    seedfinder.checkSeed(the_seed, m),
                    k,
                    msg=('m: %d, k: %d, sd: %s' % (m, k, sd)),
                )
                # self.assertEqual(seedfinder.checkSeeds(sds, m), k)


if __name__ == '__main__':
    unittest.main(verbosity=2)
