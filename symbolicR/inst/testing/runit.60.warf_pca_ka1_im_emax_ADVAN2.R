test.warf_pca_ka1_to_emax1_ADVAN2 = function(){
pks = '
   IF (NEWIND.LE.1) LN2=LOG(2)

   FSZV=WT/70
   FSZCL=FSZV**0.75
; those PPVs should be ETAs, so this control file is not very correct
   CL=FSZCL*POP_CL*EXP(PPV_CL)
   V=FSZV*POP_V*EXP(PPV_V)
   TABS=POP_TABS*EXP(PPV_TABS)
   TLAG=POP_LAG*EXP(PPV_LAG)

   KA=LN2/TABS
   ALAG1=TLAG
   S2=V

   S0=POP_S0*EXP(PPV_S0)
   EMAX=POP_EMAX*EXP(PPV_EMAX)
   C50=POP_C50*EXP(PPV_C50)

'
    str2strs = function(x){
        y = unlist(strsplit(x, '\n'))
        y = y[ !grepl('^\\s*;',y) ]
        y = y[ !grepl('^\\s*$',y) ]
        y = gsub(';.*$','', y)
        y
    }
    eqns = str2strs(pks)
    eqns1 = nonmem.eqns.to.r.eqns(eqns)
    re = analyse.0.0(eqns1, socp.filter=socp.filter.known.pattern)
    checkEquals( re$other.rels[[1]],  parse1(' EMAX = POP_EMAX * exp(PPV_EMAX)') )
    checkEquals( re$other.rels[[2]],  parse1(' C50 = POP_C50 * exp(PPV_C50)  ' ))
}
