test.warf_base = function(){
    pks = '
   ; COVARIATE MODEL
   TVCL=THETA(1)
   TVV=THETA(2)
   TVKA=THETA(3)
   TVLG=THETA(4)

   ; MODEL FOR RANDOM BETWEEN SUBJECT VARIABILITY
   CL=TVCL*EXP(ETA(1))
   V=TVV*EXP(ETA(2))
   KA=TVKA*EXP(ETA(3))
   ALAG1=TVLG*EXP(ETA(4))


   ; SCALE CONCENTRATIONS
   S2=V

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
    checkEquals( re$analysed.defns[[4]][[1]]$defn$type, 'log normal')
    checkEquals( re$analysed.defns[[4]][[1]]$defn$mean,  quote(log(THETA[4])))
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^CL=', ooo) ], 'CL=TVCL*EXP(ETA(1))')
}

test.warf_allomCL = function(){
    pks = '
   ; COVARIATE MODEL
   TVCL=THETA(1)*(WT/70)**0.75
   TVV=THETA(2)
   TVKA=THETA(3)
   TVLG=THETA(4)

   ; MODEL FOR RANDOM BETWEEN SUBJECT VARIABILITY
   CL=TVCL*EXP(ETA(1))
   V=TVV*EXP(ETA(2))
   KA=TVKA*EXP(ETA(3))
   ALAG1=TVLG*EXP(ETA(4))


   ; SCALE CONCENTRATIONS
   S2=V
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
    checkEquals( re$analysed.defns[[4]][[1]]$defn$type, 'log normal')
    checkEquals( re$analysed.defns[[4]][[1]]$defn$mean,  quote(log(THETA[4])))
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^TVCL=', ooo) ], 'TVCL=THETA(1)*0.0413215372645583*WT**0.75')
}

test.warf_wt_on_CL = function(){
    pks = '
   ; COVARIATE MODEL
   TVCL=THETA(1)*WT/70
   TVV=THETA(2)
   TVKA=THETA(3)
   TVLG=THETA(4)

   ; MODEL FOR RANDOM BETWEEN SUBJECT VARIABILITY
   CL=TVCL*EXP(ETA(1))
   V=TVV*EXP(ETA(2))
   KA=TVKA*EXP(ETA(3))
   ALAG1=TVLG*EXP(ETA(4))


   ; SCALE CONCENTRATIONS
   S2=V

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
    checkEquals( re$analysed.defns[[4]][[1]]$defn$type, 'log normal')
    checkEquals( re$analysed.defns[[4]][[1]]$defn$mean,  quote(log(THETA[4])))
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^TVCL=', ooo) ], 'TVCL=THETA(1)*0.0142857142857143*WT')
}

test.warf_wt_on_V = function(){
    pks = '

   ; COVARIATE MODEL
   TVCL=THETA(1)
   TVV=THETA(2)*WT/70
   TVKA=THETA(3)
   TVLG=THETA(4)

   ; MODEL FOR RANDOM BETWEEN SUBJECT VARIABILITY
   CL=TVCL*EXP(ETA(1))
   V=TVV*EXP(ETA(2))
   KA=TVKA*EXP(ETA(3))
   ALAG1=TVLG*EXP(ETA(4))


   ; SCALE CONCENTRATIONS
   S2=V

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
    checkEquals( re$analysed.defns[[4]][[1]]$defn$type, 'log normal')
    checkEquals( re$analysed.defns[[4]][[1]]$defn$mean,  quote(log(THETA[4])))
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^TVV=', ooo) ], 'TVV=THETA(2)*0.0142857142857143*WT')
}

test.warf_wt_on_Vpwr = function(){
    pks = '
   ; COVARIATE MODEL
   TVCL=THETA(1)
   TVKA=THETA(3)
   TVLG=THETA(4)
   PWR=THETA(5)
   TVV=THETA(2)*(WT/70)**PWR

   ; MODEL FOR RANDOM BETWEEN SUBJECT VARIABILITY
   CL=TVCL*EXP(ETA(1))
   V=TVV*EXP(ETA(2))
   KA=TVKA*EXP(ETA(3))
   ALAG1=TVLG*EXP(ETA(4))


   ; SCALE CONCENTRATIONS
   S2=V
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
    checkEquals( re$analysed.defns[[4]][[1]]$defn$type, 'log normal')
    checkEquals( re$analysed.defns[[4]][[1]]$defn$mean,  quote( log(THETA[2]) + THETA[5] * log(0.0142857142857143 * WT)))
    pkO = pattern.pk.default(L=re)
    txt = export.pk.default(pkO)
    ooo = capture.output( cat(txt$txt))
    checkEquals( ooo[ grep('^TVV=', ooo) ], 'TVV=THETA(4)*(0.0142857142857143*WT)**THETA(5)')
}

