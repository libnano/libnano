from pygments.formatters import (
    TerminalTrueColorFormatter,
)
from pygments.style import Style
from pygments.token import (
    Text,
    String
)
from pygments import highlight
from pygments.lexer import RegexLexer
from pygments.token import *

TRANTAB = str.maketrans('acgt', 'ACGT')

def align_complement(fwd: str, rev: str):
    rc_rev = reverseComplement(rev)
    alignment = force_align(rc_rev, fwd)
    # print(alignment)
    format_force_align(reverse(rev), fwd, alignment, do_print=True)
# end def

class DNALex(RegexLexer):
    name = 'DNALex'
    aliases = ['dna']
    filenames = ['*.dna']
    flags = 0

    tokens = {
        'root': [
            (r'[acgt]', Text)
        ]
    }

class MisMatchStyle(Style):
    default_style = ""
    styles = {
        Text:                   'bold nounderline #ff0000',
    }

x = highlight('aGcTT', DNALex(), TerminalTrueColorFormatter(style=MisMatchStyle))
print(x.strip().translate(TRANTAB))