# -*- coding: utf-8 -*-
import pytest
import os
import _setup
from libnano.scripts import geneprober

if os.environ.get('IS_TRAVIS') is None:
    def testCLI():
        with pytest.raises(SystemExit):
            # need this because click calls sys.exit(0) and pytest fails for that
            geneprober.cli()
    # end def

    def test_listDetails():
        with geneprober.disable_cache():
            geneprober.listDetails('mouse', ['DAXX'])
    # end def

    def test_listDetailsFail():
        a_fake_gene_symbol: str = 'DOOF'
        with geneprober.disable_cache():
            with pytest.raises(KeyError):
                geneprober.listDetails('mouse', [a_fake_gene_symbol])
    # end def



