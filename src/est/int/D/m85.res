Sat Sep 25 06:25:50 CDT 2021
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
$DATA ../../../../data/int/D/dat85.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   42594.8873007965        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.8643E+02  6.6788E+02 -1.2042E+01  5.2440E+02  5.2673E+01 -3.2735E+03 -1.7873E+03 -4.6471E+01 -2.6488E+03 -7.7688E+02
            -8.3590E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -592.567555061799        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1063E+00  2.1050E+00  8.8398E-01  1.8745E+00  1.0564E+00  4.7538E+00  5.0859E+00  9.7111E-01  2.1231E+00  1.3507E+00
             1.2603E+01
 PARAMETER:  2.0105E-01  8.4430E-01 -2.3325E-02  7.2832E-01  1.5491E-01  1.6590E+00  1.7265E+00  7.0681E-02  8.5289E-01  4.0059E-01
             2.6339E+00
 GRADIENT:  -1.2424E+01  3.5040E+01 -5.6572E+01  1.4390E+02  1.2349E+01  1.3675E+02  5.0839E+01  3.9291E+00  1.1824E+01  2.5509E+01
             1.7666E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -647.493328168253        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  7.8556E-01  3.4231E+00  4.0543E+01  2.8612E+00  3.2294E+00  3.0495E+00  1.3020E+01  6.4878E-01  2.3000E+00  1.4602E+00
             1.2646E+01
 PARAMETER: -1.4135E-01  1.3305E+00  3.8024E+00  1.1512E+00  1.2723E+00  1.2150E+00  2.6665E+00 -3.3266E-01  9.3293E-01  4.7854E-01
             2.6374E+00
 GRADIENT:  -8.2290E+01  3.2713E+01 -7.0996E+00  7.7681E+01  3.7565E+01  6.3340E+01  1.3149E+01  7.3151E-03  1.2159E+01  2.6422E+01
             2.1414E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -776.356879960447        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.3124E+00  9.7564E-01  9.9980E+00  1.3717E+00  1.9832E+00  2.6703E+00  5.6571E+00  1.1051E+00  1.2625E+00  6.5416E-01
             1.2584E+01
 PARAMETER:  3.7184E-01  7.5340E-02  2.4024E+00  4.1605E-01  7.8472E-01  1.0822E+00  1.8329E+00  1.9993E-01  3.3308E-01 -3.2440E-01
             2.6324E+00
 GRADIENT:   2.2886E+01 -3.1128E+00  8.9093E-02 -2.7302E+01 -1.6205E+01  4.5232E+01  1.5448E+01  4.2643E-01  1.4484E+01  7.5797E+00
             2.5403E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -809.019462913057        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      295
 NPARAMETR:  1.0676E+00  7.1316E-01  2.3473E+01  1.4400E+00  2.4612E+00  2.1179E+00  6.0506E+00  5.7172E-01  1.1348E+00  2.9603E-01
             1.0695E+01
 PARAMETER:  1.6542E-01 -2.3805E-01  3.2559E+00  4.6463E-01  1.0006E+00  8.5041E-01  1.9002E+00 -4.5911E-01  2.2649E-01 -1.1173E+00
             2.4697E+00
 GRADIENT:  -8.8909E+00 -4.9417E+00 -2.1368E+00  5.5033E+00  1.9485E+01  6.7845E-01  8.7537E-01  2.4704E-02 -8.5796E-01  1.4929E+00
             2.2747E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -809.991241007698        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  1.0912E+00  9.7756E-01  1.6146E+01  1.2652E+00  2.2608E+00  2.1092E+00  5.5042E+00  4.8205E-01  1.0019E+00  2.4957E-01
             1.0672E+01
 PARAMETER:  1.8727E-01  7.7306E-02  2.8817E+00  3.3521E-01  9.1572E-01  8.4630E-01  1.8055E+00 -6.2971E-01  1.0192E-01 -1.2880E+00
             2.4677E+00
 GRADIENT:   1.2723E+00 -1.5171E+00 -2.0862E-01 -2.7156E+00 -4.6825E+00 -9.2211E-01  4.6007E+00  2.6148E-02 -2.6368E-01  1.0645E+00
            -1.3240E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -809.995537541829        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      436
 NPARAMETR:  1.0897E+00  1.0422E+00  1.5943E+01  1.2388E+00  2.2802E+00  2.1134E+00  5.3489E+00  4.6200E-01  9.8829E-01  2.3214E-01
             1.0662E+01
 PARAMETER:  1.8586E-01  1.4138E-01  2.8690E+00  3.1418E-01  9.2425E-01  8.4830E-01  1.7769E+00 -6.7219E-01  8.8217E-02 -1.3604E+00
             2.4667E+00
 GRADIENT:   7.1696E-01 -7.6216E-01 -3.2682E-01 -1.5792E+00 -2.4796E+00 -5.0772E-01  2.8914E+00  2.3277E-02 -2.0939E-01  9.2615E-01
            -2.9186E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -809.998820494547        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      506
 NPARAMETR:  1.0881E+00  1.0876E+00  1.6422E+01  1.2222E+00  2.3035E+00  2.1165E+00  5.2403E+00  4.2981E-01  9.8216E-01  2.1158E-01
             1.0659E+01
 PARAMETER:  1.8442E-01  1.8394E-01  2.8986E+00  3.0063E-01  9.3444E-01  8.4975E-01  1.7564E+00 -7.4441E-01  8.2001E-02 -1.4532E+00
             2.4664E+00
 GRADIENT:   7.4450E-02 -2.2071E-01 -3.7648E-01 -6.3045E-01 -6.0038E-01 -1.7884E-01  1.3200E+00  1.8169E-02 -1.7283E-01  7.7151E-01
            -3.4379E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -810.019458272906        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      576
 NPARAMETR:  1.0865E+00  1.1319E+00  1.8240E+01  1.2085E+00  2.3414E+00  2.1196E+00  5.1356E+00  3.6335E-01  9.8233E-01  1.7547E-01
             1.0659E+01
 PARAMETER:  1.8293E-01  2.2392E-01  3.0036E+00  2.8941E-01  9.5075E-01  8.5121E-01  1.7362E+00 -9.1240E-01  8.2173E-02 -1.6403E+00
             2.4664E+00
 GRADIENT:  -6.3564E-01  3.6879E-01 -3.4333E-01  3.9899E-01  1.4901E+00  1.7520E-01 -5.5575E-01  9.9818E-03 -1.1335E-01  5.3105E-01
            -3.3935E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -810.532946211160        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  1.0892E+00  1.0625E+00  4.4901E+01  1.2533E+00  2.4162E+00  2.1167E+00  5.2939E+00  7.7453E-02  1.0347E+00  4.7442E-02
             1.0674E+01
 PARAMETER:  1.8542E-01  1.6064E-01  3.9045E+00  3.2575E-01  9.8220E-01  8.4985E-01  1.7665E+00 -2.4581E+00  1.3408E-01 -2.9483E+00
             2.4678E+00
 GRADIENT:   4.3249E-01  1.7285E-01  1.8964E-01  1.4256E-01 -2.9337E+00 -2.4780E-02  2.4219E-01  7.2995E-05 -1.1238E-01  3.8128E-02
            -1.2462E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -811.009634801816        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      829
 NPARAMETR:  1.0919E+00  8.4985E-01  9.8082E+01  1.3858E+00  2.4756E+00  2.1269E+00  5.9783E+00  2.4222E-02  1.1258E+00  2.0265E-02
             1.0715E+01
 PARAMETER:  1.8792E-01 -6.2695E-02  4.6858E+00  4.2625E-01  1.0065E+00  8.5467E-01  1.8881E+00 -3.6205E+00  2.1850E-01 -3.7989E+00
             2.4717E+00
 GRADIENT:   1.5994E-03  1.7304E-02  3.8108E-02 -5.6154E-02 -1.3266E-01 -8.1980E-02 -5.6433E-02  1.8529E-06 -1.0880E-01  6.8628E-03
             6.1176E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -811.013681120592        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1004
 NPARAMETR:  1.0920E+00  8.4931E-01  7.5720E+01  1.3857E+00  2.4632E+00  2.1277E+00  5.9791E+00  3.9507E-02  1.1285E+00  3.0653E-02
             1.0713E+01
 PARAMETER:  1.8805E-01 -6.3327E-02  4.4270E+00  4.2619E-01  1.0015E+00  8.5502E-01  1.8883E+00 -3.1313E+00  2.2088E-01 -3.3850E+00
             2.4714E+00
 GRADIENT:   8.0750E-02  2.8379E-02  2.8897E-02 -5.2733E-02 -9.0998E-02  1.7569E-03  1.6380E-02  8.6438E-06  3.2510E-02  1.5716E-02
            -2.2637E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -811.015914028378        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1160
 NPARAMETR:  1.0919E+00  8.4886E-01  6.7800E+01  1.3857E+00  2.4631E+00  2.1276E+00  5.9809E+00  4.0356E-02  1.1271E+00  2.9485E-02
             1.0714E+01
 PARAMETER:  1.8795E-01 -6.3864E-02  4.3166E+00  4.2619E-01  1.0014E+00  8.5500E-01  1.8886E+00 -3.1100E+00  2.1967E-01 -3.4239E+00
             2.4715E+00
 GRADIENT:  -2.0858E-03 -1.4951E-02 -2.5501E-03 -1.2090E-01  8.9385E-01 -1.0483E-03  7.5538E-02  1.1308E-05  5.5815E-02  1.4554E-02
            -1.0622E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -811.023505082651        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1339
 NPARAMETR:  1.0914E+00  8.4769E-01  5.4476E+01  1.3824E+00  2.4417E+00  2.1284E+00  5.9717E+00  4.0256E-02  1.1234E+00  1.5700E-02
             1.0708E+01
 PARAMETER:  1.8747E-01 -6.5238E-02  4.0978E+00  4.2383E-01  9.9271E-01  8.5535E-01  1.8870E+00 -3.1125E+00  2.1636E-01 -4.0541E+00
             2.4709E+00
 GRADIENT:  -9.3278E-02 -1.6184E-01 -1.2563E-03 -5.4542E-02 -9.1050E-02  7.5580E-02 -1.6055E-01  1.7509E-05 -5.0397E-03  4.1211E-03
            -9.1199E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -811.025484310764        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1515
 NPARAMETR:  1.0919E+00  8.5475E-01  5.5117E+01  1.3801E+00  2.4437E+00  2.1276E+00  5.9614E+00  4.0170E-02  1.1205E+00  1.0000E-02
             1.0715E+01
 PARAMETER:  1.8795E-01 -5.6941E-02  4.1095E+00  4.2218E-01  9.9353E-01  8.5498E-01  1.8853E+00 -3.1146E+00  2.1381E-01 -4.5180E+00
             2.4716E+00
 GRADIENT:  -5.0236E-03 -3.5251E-05 -9.8283E-05 -1.0383E-02  1.8990E-04  2.9061E-03 -2.3986E-04  1.6830E-05 -1.4632E-03  0.0000E+00
             5.1043E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -811.025500378474        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1691
 NPARAMETR:  1.0920E+00  8.5428E-01  5.5202E+01  1.3804E+00  2.4438E+00  2.1276E+00  5.9632E+00  4.0168E-02  1.1206E+00  1.0000E-02
             1.0715E+01
 PARAMETER:  1.8798E-01 -5.7500E-02  4.1110E+00  4.2240E-01  9.9356E-01  8.5499E-01  1.8856E+00 -3.1147E+00  2.1389E-01 -4.5104E+00
             2.4716E+00
 GRADIENT:   5.0405E-03  8.9242E-04  3.1197E-05  2.4649E-03 -4.9658E-03  7.6663E-03  1.2964E-02  1.6933E-05 -7.6707E-03  0.0000E+00
             1.1704E-02

0ITERATION NO.:   77    OBJECTIVE VALUE:  -811.025500942717        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1748
 NPARAMETR:  1.0920E+00  8.5419E-01  5.5208E+01  1.3805E+00  2.4439E+00  2.1275E+00  5.9632E+00  4.0154E-02  1.1207E+00  1.0000E-02
             1.0715E+01
 PARAMETER:  1.8797E-01 -5.7600E-02  4.1111E+00  4.2243E-01  9.9358E-01  8.5497E-01  1.8856E+00 -3.1150E+00  2.1399E-01 -4.5106E+00
             2.4716E+00
 GRADIENT:   9.1105E-04 -1.8437E-03 -1.8091E-04 -2.0266E-03  1.2994E-03 -2.3985E-04  9.0418E-03  1.6728E-05 -4.2934E-03  0.0000E+00
             1.3226E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1748
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5340E-02  3.8319E-02 -1.6032E-06 -7.5378E-02  1.2318E-06
 SE:             2.8349E-02  2.3595E-02  9.6640E-06  1.4121E-02  9.7243E-05
 N:                     100         100         100         100         100

 P VAL.:         5.8843E-01  1.0437E-01  8.6824E-01  9.4261E-08  9.8989E-01

 ETASHRINKSD(%)  5.0261E+00  2.0953E+01  9.9968E+01  5.2691E+01  9.9674E+01
 ETASHRINKVR(%)  9.7996E+00  3.7515E+01  1.0000E+02  7.7619E+01  9.9999E+01
 EBVSHRINKSD(%)  5.4822E+00  1.6136E+01  9.9953E+01  5.4926E+01  9.9567E+01
 EBVSHRINKVR(%)  1.0664E+01  2.9667E+01  1.0000E+02  7.9683E+01  9.9998E+01
 RELATIVEINF(%)  8.9143E+01  3.3482E+01  4.2126E-06  9.6799E+00  3.4245E-04
 EPSSHRINKSD(%)  5.2219E+00
 EPSSHRINKVR(%)  1.0171E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -811.02550094271703     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       843.06385882569373     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    45.22
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -811.026       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.09E+00  8.54E-01  5.52E+01  1.38E+00  2.44E+00  2.13E+00  5.96E+00  4.02E-02  1.12E+00  1.00E-02  1.07E+01
 


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
+        1.76E+02
 
 TH 2
+        2.43E+01  3.19E+01
 
 TH 3
+       -4.31E-03  1.54E-02  8.73E-05
 
 TH 4
+       -7.56E+00  3.17E+01  5.03E-03  1.24E+02
 
 TH 5
+       -4.10E-01 -5.81E+00 -4.39E-02 -1.19E+01  3.31E+01
 
 TH 6
+        5.78E-01 -8.44E-01  4.30E-04  1.63E+00 -5.15E-01  3.39E+01
 
 TH 7
+        4.55E-01  4.05E+00 -9.00E-04 -7.79E+00  5.35E-01 -1.56E-01  2.98E+00
 
 TH 8
+       -6.32E+00  1.11E+01 -1.18E-02 -1.18E+00  8.99E-01 -1.51E+00 -1.80E-03  1.40E+01
 
 TH 9
+        2.77E-01  2.28E+01 -8.14E-03 -3.13E+01  3.60E+00 -7.83E-01  2.03E+00 -2.66E+00  2.90E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  5.33E+02
 
 TH11
+       -6.64E+00 -2.80E+00 -3.93E-04 -9.70E+00  7.43E-01  1.76E+00  4.10E-01 -1.43E-01  3.64E+00  0.00E+00  8.13E+00
 
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
 #CPUT: Total CPU Time in Seconds,       60.641
Stop Time:
Sat Sep 25 06:26:52 CDT 2021
