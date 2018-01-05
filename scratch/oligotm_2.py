
# cal/mol
DH = [  -7900, # AA
        -8400, # AC
        -7800, # AG
        -7200, # AT
        -8500, # CA
        -8000, # CC
        -10600, # CG
        -7800, # CT
        -8200, # GA
        -9800, # GC
        -8000, # GG
        -8400, # GT
        -7200, # TA
        -8200, # TC
        -8500,# TG
        -7900, # TT
]
# cal / (K * mol)
DS = [  -22.2, # AA
        -22.4, # AC
        -21.0, # AG
        -20.4, # AT
        -22.7, # CA
        -19.9, # CC
        -27.2, # CG
        -21.0, # CT
        -22.2, # GA
        -24.4, # GC
        -19.9, # GG
        -22.4, # GT
        -21.3, # TA
        -22.2, # TC
        -22.7, # TG
        -22.2, # TT
]

# Init. w/term. G·C   
DH_INIT_GC = 100
DS_INIT_GC = -2.8
# Init. w/term. A·T   2.3 4.1
DH_INIT_AT = 2300
DS_INIT_AT = 4.1

DH_SYMMETRY  = 0
DS_SYMMETRY =  -1.4



from libnano import seqint
import math

def dNTPAndMg2PlusToNaPlus(dv_conc, dntp_conc):
    # Ahsen et al., 2001 
    # http://www.clinchem.org/content/47/11/1956.full
    # concentration should be in mmol/L
    # also has a formula for accounting for DMSO
    # Ahsen uses dS = dS0 + 0.847 * (N-1) * math.log(mvc)
    # must be a units thing
    if dv_conc == 0:
        dntp = 0
    if dv_conc < dntp_conc:
        dv_conc = dntp_conc
    return 120*math.sqrt(dv_conc - dntp_conc)

def saltCorrectDS(dS, seq_len, mvc, dvc, dntpc):
    # SantaLucia 1998 eq 7
    # http://bib.oxfordjournals.org/content/early/2010/12/21/bib.bbq081.full.pdf
    # above arg concentrations are in mM/L
    N = seq_len - 1
    d2m = dNTPAndMg2PlusToNaPlus(dvc, dntpc)
    mvc = mvc + d2m
    return dS + 0.368 * N * math.log(mvc/1000.)
# end def

def getTM(seq):
    R = 1.987 # cal/K⋅mol
    C_T = 50e-9 # mol
    mv_conc = 50 # mmol
    dv_conc = 1.5
    dntp_conc = 0.2
    n = len(seq)
    dS = 0
    dH = 0
    if n > 1:
        for i in range(n-1):
            si = seqint.seq2Int(seq[i:i+2])
            dH += DH[si]
            print(seq[i:i+2], si, DH[si], DH[si] - (37+273)*DS[si])
            dS += DS[si]
        if seq[0] == "C" or seq[0] == "G":
            dH += DH_INIT_GC
            dS += DS_INIT_GC
            print("init GC 0")
        else:
            dH += DH_INIT_AT
            dS += DS_INIT_AT
        if seq[-1] == "C" or seq[-1] == "G":
            dH += DH_INIT_GC
            dS += DS_INIT_GC
        else:
            dH += DH_INIT_AT
            dS += DS_INIT_AT

        dS = saltCorrectDS(dS, n, mv_conc, dv_conc, dntp_conc)
        print("ct", C_T, dS)
        TM = dH/(dS + R*math.log(C_T/4))
        dG = dH - (37+273)*dS
        return TM - 273, dG
