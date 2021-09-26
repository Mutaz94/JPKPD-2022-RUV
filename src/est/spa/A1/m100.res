Sat Sep 25 08:22:59 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/A1/dat100.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m100.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1283.88790279514        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4531E+02 -2.0276E+01 -2.5702E-01 -7.6418E+00  8.1710E+01 -4.4317E+00 -2.6413E+01 -1.0456E+01 -1.5909E+01 -2.9788E+01
            -6.2765E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1459.50761566060        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7815E-01  1.0534E+00  1.1112E+00  9.7847E-01  1.0131E+00  9.8421E-01  1.1111E+00  9.3181E-01  1.0144E+00  8.4452E-01
             2.0057E+00
 PARAMETER:  7.7907E-02  1.5206E-01  2.0541E-01  7.8236E-02  1.1298E-01  8.4089E-02  2.0534E-01  2.9378E-02  1.1425E-01 -6.8985E-02
             7.9597E-01
 GRADIENT:   4.9039E+01 -1.1026E+01  4.1866E+00 -2.0910E+01  1.0395E+01 -3.5648E+00  1.5137E+00  7.2321E-01  7.1935E-01 -8.5370E+00
            -7.1104E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1462.12947019989        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.6777E-01  8.2674E-01  9.1100E-01  1.1297E+00  8.3486E-01  9.9552E-01  1.4231E+00  3.0865E-01  9.1396E-01  8.7829E-01
             1.9921E+00
 PARAMETER:  6.7244E-02 -9.0267E-02  6.7885E-03  2.2199E-01 -8.0494E-02  9.5508E-02  4.5282E-01 -1.0756E+00  1.0027E-02 -2.9776E-02
             7.8921E-01
 GRADIENT:   2.2491E+01  6.5068E+00 -7.2118E+00  1.2430E+01  7.4322E+00  7.7670E-01  6.2930E+00  5.0764E-01  7.2453E-01  3.1381E+00
             6.9049E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1462.74174286255        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.5363E-01  7.6695E-01  9.6329E-01  1.1596E+00  8.2889E-01  9.9268E-01  1.3353E+00  1.7554E-01  9.2171E-01  8.8166E-01
             1.9774E+00
 PARAMETER:  5.2522E-02 -1.6534E-01  6.2594E-02  2.4810E-01 -8.7668E-02  9.2649E-02  3.8913E-01 -1.6399E+00  1.8474E-02 -2.5946E-02
             7.8179E-01
 GRADIENT:  -6.4381E+00  2.2391E+00  2.6228E+00  1.3214E+00 -5.6773E+00  1.9853E-01 -1.3753E+00  1.4763E-01  2.2232E-01  1.4438E-01
            -2.3926E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1463.31687538164        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.5680E-01  6.2063E-01  9.8876E-01  1.2487E+00  7.9385E-01  9.9058E-01  1.5620E+00  2.8955E-02  8.7390E-01  8.9270E-01
             1.9661E+00
 PARAMETER:  5.5838E-02 -3.7702E-01  8.8692E-02  3.2210E-01 -1.3086E-01  9.0537E-02  5.4594E-01 -3.4420E+00 -3.4786E-02 -1.3501E-02
             7.7607E-01
 GRADIENT:   4.1589E+00  1.6153E+00  2.3979E-01  5.3572E+00 -7.1450E-01 -4.2259E-02 -6.3865E-01  4.6160E-03 -1.0993E+00  4.1221E-01
            -1.2758E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1463.31964051673        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.5565E-01  6.1277E-01  9.7609E-01  1.2497E+00  7.8545E-01  9.9058E-01  1.5901E+00  2.4786E-02  8.7363E-01  8.8207E-01
             1.9672E+00
 PARAMETER:  5.4634E-02 -3.8976E-01  7.5800E-02  3.2291E-01 -1.4149E-01  9.0531E-02  5.6379E-01 -3.5975E+00 -3.5104E-02 -2.5480E-02
             7.7661E-01
 GRADIENT:   1.5708E+00  5.8528E-01  4.2401E-02  1.8844E+00 -1.8114E-01 -1.7115E-02 -2.1613E-01  3.4780E-03 -3.9272E-01  1.5382E-01
            -4.6521E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1463.55339303131        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      493
 NPARAMETR:  9.5462E-01  4.7747E-01  9.7107E-01  1.3306E+00  7.3813E-01  9.9163E-01  1.9153E+00  1.0000E-02  8.4700E-01  8.6905E-01
             1.9654E+00
 PARAMETER:  5.3553E-02 -6.3924E-01  7.0641E-02  3.8566E-01 -2.0363E-01  9.1597E-02  7.4987E-01 -6.3206E+00 -6.6052E-02 -4.0349E-02
             7.7571E-01
 GRADIENT:  -7.9950E+00  1.5534E+00  3.4823E+00 -2.4405E+00 -6.6942E+00 -2.9828E-01  6.2208E-01  0.0000E+00  7.1049E-01 -5.3605E-01
             2.2621E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1463.67043283324        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      668
 NPARAMETR:  9.5600E-01  4.0832E-01  1.0593E+00  1.3822E+00  7.6279E-01  9.9083E-01  2.0026E+00  1.0000E-02  8.3529E-01  9.2318E-01
             1.9702E+00
 PARAMETER:  5.5005E-02 -7.9571E-01  1.5758E-01  4.2369E-01 -1.7077E-01  9.0789E-02  7.9446E-01 -8.5289E+00 -7.9980E-02  2.0064E-02
             7.7813E-01
 GRADIENT:  -1.5072E+00  3.9592E-01  2.7735E-01  1.0497E+00 -4.6680E-01 -1.2636E-01  2.4384E-02  0.0000E+00 -6.9377E-02 -2.0374E-02
             1.4115E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1463.77068706897        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      848
 NPARAMETR:  9.5711E-01  2.8532E-01  1.1127E+00  1.4569E+00  7.5264E-01  9.9096E-01  2.2086E+00  1.0000E-02  8.2302E-01  9.5977E-01
             1.9687E+00
 PARAMETER:  5.6160E-02 -1.1542E+00  2.0678E-01  4.7631E-01 -1.8417E-01  9.0922E-02  8.9235E-01 -1.4048E+01 -9.4776E-02  5.8940E-02
             7.7738E-01
 GRADIENT:   6.0692E+00 -1.9745E-01  8.6503E-01 -3.1862E+00 -1.2119E+00  4.9531E-01 -6.2666E-01  0.0000E+00  2.0140E-01  4.2365E-01
            -2.8236E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1463.91544461445        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1023
 NPARAMETR:  9.5159E-01  1.4971E-01  1.1375E+00  1.5447E+00  7.2823E-01  9.8781E-01  3.1024E+00  1.0000E-02  7.9995E-01  9.6807E-01
             1.9699E+00
 PARAMETER:  5.0380E-02 -1.7990E+00  2.2885E-01  5.3483E-01 -2.1714E-01  8.7733E-02  1.2322E+00 -2.4950E+01 -1.2321E-01  6.7551E-02
             7.7798E-01
 GRADIENT:  -1.5113E+00  1.1663E+00  1.9379E+00  1.1400E+01 -4.0353E+00 -1.4623E-01  3.6447E-01  0.0000E+00 -3.4429E-01 -1.5781E-02
             4.3079E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1464.00862377001        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1199
 NPARAMETR:  9.4931E-01  7.0362E-02  1.1313E+00  1.5868E+00  7.0812E-01  9.8647E-01  4.3537E+00  1.0000E-02  7.8980E-01  9.7235E-01
             1.9643E+00
 PARAMETER:  4.7978E-02 -2.5541E+00  2.2340E-01  5.6174E-01 -2.4514E-01  8.6380E-02  1.5710E+00 -3.8980E+01 -1.3597E-01  7.1956E-02
             7.7515E-01
 GRADIENT:  -3.1571E+00  2.3915E-01  1.2795E+00  8.2039E+00 -3.1018E+00 -2.9714E-01 -2.2935E-01  0.0000E+00 -1.5485E-01  2.5229E-01
             4.1229E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1464.05917104132        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1374
 NPARAMETR:  9.5041E-01  3.2042E-02  1.1257E+00  1.6048E+00  6.9811E-01  9.8702E-01  6.4827E+00  1.0000E-02  7.8456E-01  9.7478E-01
             1.9564E+00
 PARAMETER:  4.9136E-02 -3.3407E+00  2.1840E-01  5.7299E-01 -2.5938E-01  8.6939E-02  1.9691E+00 -5.4118E+01 -1.4264E-01  7.4455E-02
             7.7112E-01
 GRADIENT:   1.5260E+00  6.5066E-02  7.1810E-01  2.3898E+00 -1.5785E+00  1.1118E-01 -1.9572E-02  0.0000E+00 -2.6497E-01  1.4418E-01
            -8.7327E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1464.08159531762        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1549
 NPARAMETR:  9.4945E-01  1.0833E-02  1.1338E+00  1.6172E+00  6.9717E-01  9.8646E-01  1.1218E+01  1.0000E-02  7.8159E-01  9.7579E-01
             1.9595E+00
 PARAMETER:  4.8128E-02 -4.4251E+00  2.2561E-01  5.8068E-01 -2.6072E-01  8.6369E-02  2.5175E+00 -7.5340E+01 -1.4642E-01  7.5488E-02
             7.7270E-01
 GRADIENT:   2.9007E-01  4.9411E-02  2.7993E-01  1.3882E+00 -5.4022E-01  3.0145E-02  6.3366E-02  0.0000E+00 -1.1068E-01 -2.6289E-02
            -1.3224E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1464.08267434650        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1737
 NPARAMETR:  9.4924E-01  1.0000E-02  1.1328E+00  1.6167E+00  6.9735E-01  9.8637E-01  1.1563E+01  1.0000E-02  7.8205E-01  9.7592E-01
             1.9599E+00
 PARAMETER:  4.7901E-02 -4.5102E+00  2.2467E-01  5.8038E-01 -2.6047E-01  8.6272E-02  2.5478E+00 -7.6876E+01 -1.4583E-01  7.5628E-02
             7.7290E-01
 GRADIENT:  -1.1692E-01  0.0000E+00 -5.6684E-01 -4.2154E-01  9.8334E-01  7.9606E-03 -5.7983E-03  0.0000E+00  1.1779E-01 -2.1746E-02
             1.2212E-01

0ITERATION NO.:   69    OBJECTIVE VALUE:  -1464.08307789733        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1864
 NPARAMETR:  9.4930E-01  1.0000E-02  1.1333E+00  1.6169E+00  6.9703E-01  9.8635E-01  1.1579E+01  1.0000E-02  7.8174E-01  9.7583E-01
             1.9597E+00
 PARAMETER:  4.7968E-02 -4.5102E+00  2.2510E-01  5.8051E-01 -2.6093E-01  8.6259E-02  2.5492E+00 -7.6876E+01 -1.4623E-01  7.5536E-02
             7.7277E-01
 GRADIENT:   1.2944E-02  0.0000E+00  2.5185E-02  2.1873E-02  1.4102E-02  3.2407E-03 -5.4100E-04  0.0000E+00  4.1975E-03 -1.8720E-03
            -3.5933E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1864
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.5777E-05  8.6999E-04  5.2631E-06 -1.1031E-02 -2.5321E-02
 SE:             2.9401E-02  2.0212E-03  1.4158E-04  2.7809E-02  2.2309E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9903E-01  6.6688E-01  9.7035E-01  6.9160E-01  2.5637E-01

 ETASHRINKSD(%)  1.5043E+00  9.3229E+01  9.9526E+01  6.8374E+00  2.5261E+01
 ETASHRINKVR(%)  2.9860E+00  9.9541E+01  9.9998E+01  1.3207E+01  4.4142E+01
 EBVSHRINKSD(%)  1.5072E+00  9.3983E+01  9.9499E+01  6.5249E+00  2.4507E+01
 EBVSHRINKVR(%)  2.9916E+00  9.9638E+01  9.9997E+01  1.2624E+01  4.3009E+01
 RELATIVEINF(%)  9.1705E+01  1.0344E-02  1.7042E-04  3.3363E+00  3.0599E+00
 EPSSHRINKSD(%)  3.5547E+01
 EPSSHRINKVR(%)  5.8458E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1464.0830778973336     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -728.93225133359545     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1464.083       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  1.00E-02  1.13E+00  1.62E+00  6.97E-01  9.86E-01  1.16E+01  1.00E-02  7.82E-01  9.76E-01  1.96E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.24E+03
 
 TH 2
+        0.00E+00  2.94E+03
 
 TH 3
+        3.37E+00  0.00E+00  2.21E+02
 
 TH 4
+       -3.09E+01  0.00E+00 -3.61E+01  6.09E+02
 
 TH 5
+        3.77E+00  0.00E+00 -5.10E+02 -8.66E+01  1.37E+03
 
 TH 6
+        1.64E+01  0.00E+00  3.66E+00 -9.70E+00 -8.89E+00  1.99E+02
 
 TH 7
+       -2.10E-02  0.00E+00 -4.72E-02 -8.78E-02  1.01E-01  1.68E-02  1.79E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.73E+01  0.00E+00  2.08E+01 -9.55E+00 -1.08E+01 -1.11E+00 -2.80E-03  0.00E+00  2.50E+02
 
 TH10
+       -1.22E+01  0.00E+00  5.73E+00 -1.92E+00 -6.46E+01  7.89E+00 -3.87E-02  0.00E+00  1.28E+01  6.69E+01
 
 TH11
+       -1.30E+01  0.00E+00 -1.23E+01 -9.82E+00  1.31E+00  1.74E+00 -4.91E-03  0.00E+00  9.94E+00  2.24E+01  7.03E+01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       28.888
Stop Time:
Sat Sep 25 08:23:29 CDT 2021
