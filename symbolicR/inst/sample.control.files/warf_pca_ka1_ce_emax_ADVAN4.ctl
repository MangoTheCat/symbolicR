$PROB ka1_ce_emax_ADVAN4

;O'Reilly RA, Aggeler PM, Leong LS. Studies of the coumarin anticoagulant
;drugs: The pharmacodynamics of warfarin in man.
;Journal of Clinical Investigation 1963;42(10):1542-1551

;O'Reilly RA, Aggeler PM. Studies on coumarin anticoagulant drugs
;Initiation of warfarin therapy without a loading dose.
;Circulation 1968;38:169-177

$INPUT ID TIME WT AGE SEX AMT DVID DV MDV
$DATA warfarin_conc_pca.csv
$EST METHOD=COND INTER 
MAX=9990 SIG=3 PRINT=20 NOABORT
$COV

$THETA
(0.01,0.134,1) ; POP_CL L/h/70kg
(0.01,8.11,20) ; POP_V L/70kg
(0.01,0.523,24) ; POP_TABS h
(0.01,0.823,24) ; POP_LAG h
(0.01,75.6,200) ; POP_S0
(-INF,-237.,0) ; POP_EMAX
(0.01,10.,100) ; POP_C50
(0.01,37.4,100) ; POP_TEQ

$OMEGA
0.0713 ; PPV_CL
0.0181 ; PPV_V
0.696 ; PPV_TABS
0.156 ; PPV_LAG
0.0627 ; PPV_S0
0.05 ; PPV_EMAX
0.0564 ; PPV_C50
0.0731 ; PPV_TEQ

$SIGMA 
0.00752 ; RUV_CV
0.0661 ; RUV_SD mg/L
14.6 ; RUV_FX

$SUBR ADVAN4 TRAN4

$PK
    POP_CL=THETA(1)
    POP_V=THETA(2)
    POP_TABS=THETA(3)
    POP_LAG=THETA(4)
    POP_S0=THETA(5)
    POP_EMAX=THETA(6)
    POP_C50=THETA(7)
    POP_TEQ=THETA(8)

    PPV_CL=ETA(1)
    PPV_V=ETA(2)
    PPV_TABS=ETA(3)
    PPV_LAG=ETA(4)
    PPV_S0=ETA(5)
    PPV_EMAX=ETA(6)
    PPV_C50=ETA(7)
    PPV_TEQ=ETA(8)

   IF (NEWIND.LE.1) LN2=LOG(2)

   FSZV=WT/70
   FSZCL=FSZV**0.75

   CL=FSZCL*POP_CL*EXP(PPV_CL)
   V=FSZV*POP_V*EXP(PPV_V)
   TABS=POP_TABS*EXP(PPV_TABS)
   TLAG=POP_LAG*EXP(PPV_LAG)

   S0=POP_S0*EXP(PPV_S0)
   EMAX=POP_EMAX*EXP(PPV_EMAX)
   C50=POP_C50*EXP(PPV_C50)
   TEQ=POP_TEQ*EXP(PPV_TEQ)

   KEQ=LN2/TEQ
   KA=LN2/TABS
   ALAG1=TLAG
   S2=V
   V2=V
   V3=V2*0.0001
   Q=V3*KEQ

$ERROR

    RUV_CV=EPS(1)
    RUV_SD=EPS(2)
    RUV_FX=EPS(3)
    
   CP=F
   CE=A(3)/V3
   PCA=S0 + EMAX*CE/(C50+CE)
   IF (DVID.LE.1) THEN
      Y=CP*(1+RUV_CV) + RUV_SD
   ENDIF
   IF (DVID.EQ.2) THEN
      Y=PCA + RUV_FX
   ENDIF

$TABLE ID TIME DVID CP CE PCA Y
ONEHEADER NOPRINT FILE=ka1.fit



