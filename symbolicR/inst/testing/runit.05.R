
test.pattern.pk.default <- function()
{
    test2 <- pattern.pk.default("CL = THETA(1) * EXP(ETA(1))")
    
    checkEquals(test2$CL$pkvar, 'CL')
    checkEquals(test2$CL$distrib, 'log normal')
    checkEquals(instantiate.morphism(test2$CL$mean), quote(log(THETA[1])))
    # TEST 3 : Warf WT ON CL
    
    pks <- c(
        'TVCL=THETA(1)*WT/70',
        'TVV =THETA(2)',
        'TVKA=THETA(3)',
        'TVLG=THETA(4)',
        'CL=TVCL*EXP(ETA(1))',
        'V =TVV *EXP(ETA(2))',
        'KA=TVKA*EXP(ETA(3))',
        'ALAG1=TVLG*EXP(ETA(4))',
        'S2=V'
    )
    test3 <- pattern.pk.default(pks)
    checkEquals(test3$CL$mean$defn$defn[[2]]$defn$defn$grptype, 'multiplicative')
    checkEquals(test3$CL$distrib , 'log normal', msg = "identify relationship")
    
    # TEST 4 : Warf PWR ON WT ON VOL
    
    pks <- c("TVCL=THETA(1)",
    "TVKA=THETA(3)",
    "TVLG=THETA(4)",
    "PWR=THETA(5)",
    "TVV=THETA(2)*(WT/70)**PWR",
    "CL=TVCL*EXP(ETA(1))",
    "V=TVV*EXP(ETA(2))",
    "KA=TVKA*EXP(ETA(3))",
    "ALAG1=TVLG*EXP(ETA(4))",
    "S2=V")
    
    test4 <- pattern.pk.default(pks)
    checkEquals(test4$S2$type , 'related.on.previous.variables')
    checkEquals(test4$V$mean$defn$defn[[2]]$defn$defn$grptype , 'exponential')
}
