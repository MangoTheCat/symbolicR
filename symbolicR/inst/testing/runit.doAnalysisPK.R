# git revision: $Rev: $
# Date of last change: $LastChangedDate: 11/10/2012 $
# Last changed by: $LastChangedBy: ccampbell $
# 
# Original author: ccampbell
# Copyright Mango Solutions, Chippenham, UK
###############################################################################

test.analyse.0.0 <- function()
{
    pks <- nonmem.eqns.to.r.eqns(c(
        'TVCL=THETA(1)*WT/70',
        'TVV =THETA(2)',
        'TVKA=THETA(3)',
        'TVLG=THETA(4)',
        'CL=TVCL*EXP(ETA(1))',
        'V =TVV *EXP(ETA(2))',
        'KA=TVKA*EXP(ETA(3))',
        'ALAG1=TVLG*EXP(ETA(4))',
        'S2=V'))
  
    test2 <- analyse.0.0(eqns = pks)
    checkEquals( test2$other.rels , list())
    test2 = test2$analysed.defns
    checkEquals( test2[[2]][[2]]$defn, quote(V) )

    #
    pks <- nonmem.eqns.to.r.eqns(c(
        "KA=THETA(1)+ETA(1)",
        "K=THETA(2)+ETA(2)",
        "CL=THETA(3)*WT+ETA(3)",
        "SC=CL/K/WT")   )
    
    test4 <- analyse.0.0(eqns = pks)$analysed.defns
    checkEquals( test4[[2]][[2]]$defn$mean, quote( THETA[3] * WT))
    
    # TEST 6 : sequence mapping
    
    pkt <- nonmem.eqns.to.r.eqns(c(
        "TMP1 = THETA(1)",
        "TMP2 = TMP1",
        "TMP3 = TMP2",
        "V = TMP3",
        "SC = V + TMP2"
        ))
     #   "V = TMP3 * EXP(ETA(2))")
    
    test6 <- analyse.0.0(eqns = pkt)$analysed.defns
    checkEquals(test6[[1]][[1]]$defn$mean, quote(THETA[1]))

}
