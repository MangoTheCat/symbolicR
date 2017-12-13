test.example.2 = function(){
    pks = '
;  The MUs should always be unconditionally defined, that is, they should never be
; defined in IF?THEN blocks
; THETA(1) models the Volume of sub-population 1
MU_1=THETA(1)
; THETA(2) models the clearance of sub-population 1
MU_2=THETA(2)
; THETA(3) models the Volume of sub-population 2
MU_3=THETA(3)
; THETA(4) models the clearance of sub-population 2
MU_4=THETA(4)
VCM=DEXP(MU_1+ETA(1))
K10M=DEXP(MU_2+ETA(2))
VCF=DEXP(MU_3+ETA(3))
K10F=DEXP(MU_4+ETA(4))
Q=1
IF(MIXNUM.EQ.2) Q=0
V=Q*VCM+(1.0-Q)*VCF
K=Q*K10M+(1.0-Q)*K10F
S1=V
'  
    str2strs = function(x){
        y = unlist(strsplit(x, '\n'))
        y = y[ !grepl('^\\s*;',y) ]
        y = y[ !grepl('^\\s*$',y) ]
        y
    }
    eqns = str2strs(pks)
    eqns1 = nonmem.eqns.to.r.eqns(eqns)
    re = analyse.0.0(eqns1, socp.filter=socp.filter.known.pattern)
    checkEquals( re$analysed.defns[[1]][[9]]$critical.varname, quote(Q))
    checkEquals( re$analysed.defns[[1]][[9]]$defn, quote( ifelse(MIXNUM==2, 0, 1)))
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^K10M', ooo) ], 'K10M=EXP(MU_2+ETA(2))')
}
