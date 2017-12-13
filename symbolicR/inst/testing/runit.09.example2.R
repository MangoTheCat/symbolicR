test.example.2 = function(){
    pks = '
; LCLM=log transformed clearance, male
LCLM=THETA(1)
;LCLF=log transformed clearance, female.
LCLF=THETA(2)
; CLAM=CL age slope, male
CLAM=THETA(3)
; CLAF=CL age slope, female
CLAF=THETA(4)
; LV1M=log transformed V1, male
LV1M=THETA(5)
; LV1F=log transformed V1, female
LV1F=THETA(6)
; V1AM=V1 age slope, male
V1AM=THETA(7)
; V1AF=V1 age slope, female
V1AF=THETA(8)
; LAGE=log transformed age
LAGE=DLOG(AGE)
;Mean of ETA1, the inter-subject deviation of Clearance, is ultimately modeled as linear function
;of THETA(1) to THETA(4).  Relating thetas to Mus by linear functions is not essential for ITS,
;IMP, or IMPMAP methods, but is very helpful for MCMC methods such as SAEM and BAYES.
MU_1=(1.0-GNDR)*(LCLM+LAGE*CLAM) + GNDR*(LCLF+LAGE*CLAF)
;Mean of ETA2, the inter-subject deviation of V1, is ultimately modeled as linear function of
; THETA(5) to THETA(8)
MU_2=(1.0-GNDR)*(LV1M+LAGE*V1AM) + GNDR*(LV1F+LAGE*V1AF)
MU_3=THETA(9)
MU_4=THETA(10)
CL=DEXP(MU_1+ETA(1))
V1=DEXP(MU_2+ETA(2))
Q=DEXP(MU_3+ETA(3))
V2=DEXP(MU_4+ETA(4))
S1=V1
'  
    str2strs = function(x){
        y = unlist(strsplit(x, '\n'))
        y = y[ !grepl('^\\s*;',y) ]
        y = y[ !grepl('^\\s*$',y) ]
        y
    }
    eqns = str2strs(pks)
    eqns1 = nonmem.eqns.to.r.eqns(eqns)
    L = analyse.0.0(eqns1, socp.filter=socp.filter.known.pattern)
    checkEquals( L$analysed.defns[[3]][[2]]$critical.varname, quote(V2))
    checkEquals( L$analysed.defns[[3]][[2]]$defn$type,'log normal')
    checkEquals( L$analysed.defns[[3]][[2]]$defn$mean,quote(MU_4))
    checkEquals( L$analysed.defns[[3]][[2]]$defn$var,quote(Var(ETA[4])))
}
