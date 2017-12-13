
test.example.35 = function(){
    pks = '
;;; CLHCTZ-DEFINITION START
IF(HCTZ.EQ.1) CLHCTZ = 1  ; Most common
IF(HCTZ.EQ.0) CLHCTZ = ( 1 + THETA(7))
;;; CLHCTZ-DEFINITION END


;;; CLSEX-DEFINITION START
IF(SEX.EQ.1) CLSEX = 1  ; Most common
IF(SEX.EQ.2) CLSEX = ( 1 + THETA(6))
;;; CLSEX-DEFINITION END


;;; CLRACE-DEFINITION START
IF(RACE.EQ.1) CLRACE = 1  ; Most common
IF(RACE.EQ.2) CLRACE = ( 1 + THETA(4))
IF(RACE.EQ.3) CLRACE = ( 1 + THETA(5))
;;; CLRACE-DEFINITION END


;;; CL-RELATION START
CLCOV=CLRACE*CLSEX*CLHCTZ
;;; CL-RELATION END

 
    TVCL = THETA(1)
    TVCL = CLCOV * TVCL
    
    TVV  = THETA(2)
    TVKA = THETA(3)
    CL   = EXP(TVCL + ETA(1))
    V    = EXP(TVV  + ETA(2))
    KA   = EXP(TVKA)
    S2   = V

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
    checkEquals( re$analysed.defns[[1]][[1]]$critical.varname, quote(CL))
    checkEquals( re$analysed.defns[[1]][[1]]$defn$var, quote( Var(ETA[1]) ) )
    checkEquals( re$analysed.defns[[3]][[1]]$defn$type, 'fixed')
    checkEquals( re$analysed.defns[[3]][[1]]$defn$mean,  quote(exp(THETA[3])))
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^CL=', ooo) ], 'CL=EXP(TVCL+ETA(1))')
}
