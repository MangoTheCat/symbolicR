test.warf_pca_ka1_to_emax1_ADVAN6 = function(){
pks = '
   IF (NEWIND.LE.1) LN2=LOG(2)

   FSZV=WT/70
   FSZCL=FSZV**0.75

   CL=FSZCL*POP_CL*EXP(PPV_CL)
   V=FSZV*POP_V*EXP(PPV_V)
   TABS=POP_TABS*EXP(PPV_TABS)
   TLAG=POP_LAG*EXP(PPV_LAG)

   S0=POP_S0*EXP(PPV_S0)
   EMAX=POP_EMAX*EXP(PPV_EMAX)
   C50=POP_C50*EXP(PPV_C50)
   TEQ=POP_TEQ*EXP(PPV_TEQ)
;  Currently not able to handle the following EXIT statement
;  IF (EMAX.LT.-1) EXIT 1 101

   KA=LN2/TABS
   ALAG1=TLAG
   S2=V
   A_0(3)=S0
   KPCA=LN2/TEQ
   RPCA=S0*KPCA
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
    checkEquals( re$analysed.defns[[1]][[2]]$critical.varname, quote(V))
    checkEquals( re$analysed.defns[[1]][[2]]$defn, 
        quote( FSZV * POP_V * exp(PPV_V) ))
}
