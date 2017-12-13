test.control51 = function(){
pks = '
   CALLFL=1
   TMP1 = THETA(1)
   TMP2 = TMP1
   TMP3 = TMP2
   KA=TMP3+ETA(1)
   K=THETA(2)+ETA(2)
   CL=THETA(3)*WT+ETA(3)
   SC=CL/K/WT
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
    checkEquals( re$analysed.defns[[1]][[1]]$critical.varname, quote(KA))
    checkEquals( re$analysed.defns[[1]][[1]]$defn$type, 'normal')
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^CL=', ooo) ], "CL=TVCL+ETA(3)")
}
