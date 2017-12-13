test.example.10 = function(){
pks = '
   CALLFL=1
   MU_1=DLOG(THETA(1))
   KA=DEXP(MU_1+ETA(1))
   MU_2=DLOG(THETA(2))
   V=DEXP(MU_2+ETA(2))
   MU_3=DLOG(THETA(3))
   CL=DEXP(MU_3+ETA(3))
   SC=V/1000
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
    checkEquals( re$analysed.defns[[1]][[1]]$critical.varname, quote(MU_1))
    checkEquals( re$analysed.defns[[1]][[1]]$defn$type, 'fixed')
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^CL=', ooo) ],"CL=EXP(MU_3+ETA(3))") 
    # the CALLFL=1 is dropped in pkObjects
    # This is by design
    # PKObject should only define how PK parameters are constructed
    # The control flow, or say, the dynamic system,
    # should be extracted by Dynamic System extracter,
    # hence we should not capture this here
    checkEquals( grep('CALLFL',ooo) , integer(0))
}
