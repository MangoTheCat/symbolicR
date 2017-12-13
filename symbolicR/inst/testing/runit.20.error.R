test.error.block = function(){
    pks = '
CALLFL=0
ETYPE=1
IF(CMT.NE.1) ETYPE=0
IPRED=F
Y = F + F*ETYPE*EPS(1) + F*(1.0-ETYPE)*EPS(2)
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
    re = pattern.error.default(eqns1)
    checkEquals( names(re), c('IPRED','Y'))
}
