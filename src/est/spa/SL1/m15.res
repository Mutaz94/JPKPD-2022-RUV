Sat Sep 25 10:20:45 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat15.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1697.70505663353        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.6415E+01 -1.0980E+02 -4.9250E+01 -9.3802E+01  1.1635E+02  3.7832E+01 -1.0776E+01 -3.8323E+00 -3.0520E+01 -8.7636E-01
             3.2828E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1715.54310990257        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0425E+00  1.0690E+00  9.7675E-01  1.0266E+00  9.2359E-01  8.9110E-01  1.0143E+00  1.1171E+00  1.1570E+00  8.7337E-01
             9.2351E-01
 PARAMETER:  1.4160E-01  1.6670E-01  7.6472E-02  1.2627E-01  2.0512E-02 -1.5302E-02  1.1419E-01  2.1076E-01  2.4582E-01 -3.5393E-02
             2.0426E-02
 GRADIENT:   6.0985E+01 -8.2747E+00 -7.1791E+00  9.7987E+00  1.3149E+01  4.9529E-01  1.1900E+00 -2.1630E+00  1.0841E+01 -9.9090E-01
             3.1751E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1716.46309288062        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0368E+00  1.2143E+00  1.1183E+00  9.5378E-01  1.0339E+00  8.8083E-01  8.9881E-01  1.6079E+00  1.1907E+00  9.1836E-01
             9.1456E-01
 PARAMETER:  1.3615E-01  2.9418E-01  2.1177E-01  5.2678E-02  1.3332E-01 -2.6895E-02 -6.6829E-03  5.7496E-01  2.7455E-01  1.4834E-02
             1.0691E-02
 GRADIENT:   4.5319E+01  1.5844E+01 -6.9299E+00  1.7450E+01  2.0392E+01 -3.9818E+00  2.7003E+00  2.9372E+00 -1.4311E+00 -3.3599E+00
            -1.6155E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1717.51365857546        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      337
 NPARAMETR:  1.0493E+00  1.3533E+00  9.9369E-01  8.5913E-01  1.0338E+00  9.0096E-01  7.9810E-01  1.4955E+00  1.3116E+00  9.4734E-01
             9.1733E-01
 PARAMETER:  1.4810E-01  4.0255E-01  9.3673E-02 -5.1830E-02  1.3326E-01 -4.2952E-03 -1.2552E-01  5.0249E-01  3.7125E-01  4.5899E-02
             1.3707E-02
 GRADIENT:   6.4542E+00  8.5746E+00  3.2692E+00  8.3466E+00 -7.9124E+00  1.1654E+00 -1.2314E+00 -1.2664E+00 -6.0683E-01  3.0026E-01
            -1.8429E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1718.81932534018        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  1.0449E+00  1.7469E+00  7.4619E-01  5.9123E-01  1.1305E+00  8.9487E-01  7.4445E-01  1.7701E+00  1.6653E+00  9.8440E-01
             9.1858E-01
 PARAMETER:  1.4390E-01  6.5783E-01 -1.9278E-01 -4.2555E-01  2.2262E-01 -1.1082E-02 -1.9511E-01  6.7105E-01  6.0998E-01  8.4280E-02
             1.5075E-02
 GRADIENT:  -7.1067E+00  1.3147E+01  5.2592E+00  1.7835E+00 -9.0517E+00 -1.7738E+00  1.0985E+00 -1.0940E-01 -4.7424E-01  4.1371E-01
            -1.2101E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1718.92049891073        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  1.0448E+00  1.7576E+00  7.2795E-01  5.8296E-01  1.1396E+00  8.9479E-01  7.4303E-01  1.7513E+00  1.6789E+00  9.8429E-01
             9.1888E-01
 PARAMETER:  1.4385E-01  6.6396E-01 -2.1752E-01 -4.3964E-01  2.3067E-01 -1.1170E-02 -1.9702E-01  6.6037E-01  6.1815E-01  8.4163E-02
             1.5401E-02
 GRADIENT:   6.7160E+01  8.9082E+01  2.8291E+00  1.3895E+01  3.9266E+00  1.8741E+00  2.1604E+00  1.5613E-01  3.1646E+00 -5.7346E-01
            -2.0549E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1718.93776736807        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      817
 NPARAMETR:  1.0437E+00  1.7537E+00  7.2618E-01  5.8197E-01  1.1377E+00  8.9860E-01  7.4206E-01  1.7508E+00  1.6778E+00  9.8448E-01
             9.1888E-01
 PARAMETER:  1.4278E-01  6.6172E-01 -2.1995E-01 -4.4134E-01  2.2899E-01 -6.9134E-03 -1.9832E-01  6.6006E-01  6.1748E-01  8.4362E-02
             1.5401E-02
 GRADIENT:   6.3831E+01  8.3284E+01  2.6184E+00  1.1575E+01  3.4376E+00  3.5546E+00  2.0926E+00  3.5179E-01  3.1602E+00 -1.9849E-01
             1.5004E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1718.94222073952        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      989
 NPARAMETR:  1.0439E+00  1.7524E+00  7.2636E-01  5.8198E-01  1.1373E+00  8.9863E-01  7.3996E-01  1.7495E+00  1.6790E+00  9.8459E-01
             9.1898E-01
 PARAMETER:  1.4294E-01  6.6098E-01 -2.1971E-01 -4.4132E-01  2.2862E-01 -6.8829E-03 -2.0116E-01  6.5932E-01  6.1817E-01  8.4474E-02
             1.5513E-02
 GRADIENT:   6.4353E+01  8.1867E+01  2.5905E+00  1.1083E+01  3.3558E+00  3.5720E+00  1.8195E+00  3.3866E-01  3.2686E+00 -1.6985E-01
             2.0481E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1718.95538561543        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1147
 NPARAMETR:  1.0445E+00  1.7522E+00  7.2411E-01  5.8195E-01  1.1370E+00  8.9897E-01  7.3939E-01  1.7477E+00  1.6796E+00  9.8490E-01
             9.1893E-01
 PARAMETER:  1.4350E-01  6.6086E-01 -2.2281E-01 -4.4137E-01  2.2840E-01 -6.5035E-03 -2.0193E-01  6.5830E-01  6.1853E-01  8.4780E-02
             1.5456E-02
 GRADIENT:   6.6169E+01  8.1037E+01  2.2347E+00  1.1232E+01  3.9240E+00  3.7203E+00  1.7417E+00  4.3890E-01  3.3746E+00 -6.3356E-02
             2.0809E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1719.08277333368        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1219
 NPARAMETR:  1.0478E+00  1.7413E+00  6.2773E-01  5.7817E-01  1.0989E+00  8.9749E-01  7.4423E-01  1.4688E+00  1.6681E+00  9.5952E-01
             9.1687E-01
 PARAMETER:  1.4672E-01  6.5465E-01 -3.6564E-01 -4.4788E-01  1.9431E-01 -8.1559E-03 -1.9540E-01  4.8441E-01  6.1169E-01  5.8675E-02
             1.3214E-02
 GRADIENT:   7.6050E+01  6.2026E+01 -2.3476E+00  9.4570E+00  6.0375E+00  2.8573E+00 -3.0142E-01  7.4771E-01  4.4836E+00  1.6087E+00
            -1.0480E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1719.15227473829        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1379
 NPARAMETR:  1.0454E+00  1.7452E+00  6.3422E-01  5.7988E-01  1.0971E+00  8.9884E-01  7.5743E-01  1.4525E+00  1.6552E+00  9.4925E-01
             9.1810E-01
 PARAMETER:  1.4443E-01  6.5689E-01 -3.5535E-01 -4.4493E-01  1.9263E-01 -6.6494E-03 -1.7783E-01  4.7326E-01  6.0389E-01  4.7913E-02
             1.4552E-02
 GRADIENT:  -6.4241E+00 -6.3516E+00  4.9610E-01 -1.2114E+00 -8.2243E-01 -2.5931E-01  4.9515E-03  8.7887E-03 -3.8021E-01  1.7015E-01
             7.0726E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1719.29019112794        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1558
 NPARAMETR:  1.0512E+00  1.8116E+00  5.7496E-01  5.4266E-01  1.1025E+00  9.0055E-01  7.5784E-01  1.4229E+00  1.7236E+00  9.4467E-01
             9.2057E-01
 PARAMETER:  1.4994E-01  6.9422E-01 -4.5345E-01 -5.1127E-01  1.9758E-01 -4.7462E-03 -1.7729E-01  4.5273E-01  6.4439E-01  4.3079E-02
             1.7238E-02
 GRADIENT:   8.2222E+00  9.0349E+00  8.9034E-01  3.9794E+00 -5.0991E+00  3.9681E-01  1.4345E+00  1.4756E-01  5.4993E-02  2.3796E-02
             8.3211E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1719.35465940412        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1745
 NPARAMETR:  1.0452E+00  1.8239E+00  5.6712E-01  5.2492E-01  1.1144E+00  8.9863E-01  7.4211E-01  1.4343E+00  1.7415E+00  9.4722E-01
             9.1693E-01
 PARAMETER:  1.4424E-01  7.0099E-01 -4.6718E-01 -5.4452E-01  2.0828E-01 -6.8848E-03 -1.9826E-01  4.6067E-01  6.5473E-01  4.5777E-02
             1.3277E-02
 GRADIENT:  -7.0745E+00 -5.8959E+00  5.8318E-01 -2.1302E+00  1.7924E+00 -4.0730E-01 -1.2444E+00 -2.3444E-01 -2.1584E+00 -1.3639E+00
            -9.8189E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1719.38997233541        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1921            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0480E+00  1.8270E+00  5.6698E-01  5.2696E-01  1.1134E+00  8.9954E-01  7.4604E-01  1.4295E+00  1.7552E+00  9.5371E-01
             9.1857E-01
 PARAMETER:  1.4686E-01  7.0266E-01 -4.6743E-01 -5.4063E-01  2.0746E-01 -5.8720E-03 -1.9297E-01  4.5735E-01  6.6260E-01  5.2605E-02
             1.5064E-02
 GRADIENT:   7.6088E+01  9.3803E+01  9.0469E-01  1.2415E+01  5.4144E-01  3.6864E+00  1.3905E+00 -2.3247E-02  3.4745E+00 -3.0668E-02
            -8.6960E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1719.39584803768        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2083             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0480E+00  1.8283E+00  5.6485E-01  5.2503E-01  1.1143E+00  8.9957E-01  7.4597E-01  1.4334E+00  1.7580E+00  9.5432E-01
             9.1869E-01
 PARAMETER:  1.4691E-01  7.0339E-01 -4.7120E-01 -5.4430E-01  2.0826E-01 -5.8332E-03 -1.9307E-01  4.6003E-01  6.6417E-01  5.3239E-02
             1.5197E-02
 GRADIENT:   7.6249E+01  9.1950E+01  5.8870E-01  1.1966E+01  1.3274E+00  3.7026E+00  1.4941E+00  7.2501E-02  3.4831E+00  5.7603E-02
             8.1978E-02

0ITERATION NO.:   72    OBJECTIVE VALUE:  -1719.39584803768        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:     2163
 NPARAMETR:  1.0480E+00  1.8281E+00  5.6489E-01  5.2498E-01  1.1145E+00  8.9957E-01  7.4588E-01  1.4333E+00  1.7592E+00  9.5433E-01
             9.1871E-01
 PARAMETER:  1.4691E-01  7.0339E-01 -4.7120E-01 -5.4430E-01  2.0826E-01 -5.8332E-03 -1.9307E-01  4.6003E-01  6.6417E-01  5.3239E-02
             1.5197E-02
 GRADIENT:   1.1292E-01  7.8019E+04 -1.1646E+05  5.0411E+04 -2.8920E-01  1.0284E-02  2.8801E-02  1.1923E+05 -1.8786E-01 -2.7442E+05
            -5.4884E+05
 NUMSIGDIG:         3.6         3.3         3.3         3.3         2.7         3.6         2.7         3.3         2.5         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2163
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.9396E-04 -3.9937E-02 -2.6691E-02  3.2706E-02 -4.8136E-02
 SE:             2.9860E-02  2.2423E-02  1.0877E-02  2.3416E-02  2.1462E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8680E-01  7.4896E-02  1.4128E-02  1.6250E-01  2.4908E-02

 ETASHRINKSD(%)  1.0000E-10  2.4881E+01  6.3562E+01  2.1552E+01  2.8098E+01
 ETASHRINKVR(%)  1.0000E-10  4.3571E+01  8.6723E+01  3.8459E+01  4.8302E+01
 EBVSHRINKSD(%)  4.6892E-01  2.4376E+01  6.5870E+01  2.2105E+01  2.5779E+01
 EBVSHRINKVR(%)  9.3565E-01  4.2811E+01  8.8351E+01  3.9324E+01  4.4913E+01
 RELATIVEINF(%)  9.9046E+01  4.2763E+00  1.3065E+00  4.9230E+00  1.5933E+01
 EPSSHRINKSD(%)  4.6372E+01
 EPSSHRINKVR(%)  7.1240E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1719.3958480376762     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -984.24502147393798     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.54
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1719.396       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.83E+00  5.65E-01  5.25E-01  1.11E+00  9.00E-01  7.46E-01  1.43E+00  1.76E+00  9.54E-01  9.19E-01
 


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
+        5.79E+08
 
 TH 2
+        4.99E+02  8.30E+06
 
 TH 3
+       -2.43E+03 -4.01E+07  1.94E+08
 
 TH 4
+        2.26E+03 -7.64E+03  3.81E+04  1.68E+08
 
 TH 5
+       -3.57E+00  3.48E+02 -2.76E+03  2.54E+03  5.61E+02
 
 TH 6
+        3.10E-01  5.85E+02 -2.83E+03  2.63E+03  2.43E-01  2.44E+02
 
 TH 7
+        3.87E-01  4.01E+02 -1.94E+03  1.80E+03 -1.55E+01 -3.39E+00  1.38E+02
 
 TH 8
+        9.83E+02 -2.12E+03  2.60E+04 -1.55E+04  1.02E+03  1.14E+03  7.88E+02  3.15E+07
 
 TH 9
+        1.00E+00 -6.00E+03  2.89E+04 -2.69E+04 -7.22E+00  2.61E-01  1.24E+01 -1.16E+04  3.11E+01
 
 TH10
+        9.34E+08  3.39E+02 -1.70E+03  1.56E+03  6.20E+08  1.60E+09  9.98E+08  6.84E+02 -1.23E+08  1.51E+09
 
 TH11
+       -7.08E+03 -5.61E+03  2.71E+04 -2.52E+04 -7.34E+03 -8.19E+03 -5.61E+03 -1.09E+04  8.37E+04 -4.87E+03  1.63E+09
 
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
 #CPUT: Total CPU Time in Seconds,       33.629
Stop Time:
Sat Sep 25 10:21:20 CDT 2021
