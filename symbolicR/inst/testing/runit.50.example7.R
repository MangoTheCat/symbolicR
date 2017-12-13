test.example.07 = function(){
pks = '
MU_1=THETA(1)
MU_2=THETA(2)
V=DEXP(MU_1+ETA(1))
CLB=DEXP(MU_2+ETA(2))
DCL1=DEXP(ETA(3))
DCL2=DEXP(ETA(4))
DCL3=DEXP(ETA(5))
S1=V
DCL=DCL1
IF(TIME.GE.5.0) DCL=DCL2
IF(TIME.GE.10.0) DCL=DCL3
CL=CLB*DCL
VC=V

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
    checkEquals( re$analysed.defns[[2]][[6]]$critical.varname, quote(CL))
    checkEquals( re$analysed.defns[[2]][[6]]$defn, 
        quote(CLB * ifelse(TIME >= 10, DCL3, ifelse(TIME >= 5, DCL2, DCL1))) )
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^CL=', ooo) ], "CL=CLB*DCL1")
}
