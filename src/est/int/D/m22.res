Sat Sep 25 05:30:36 CDT 2021
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
$DATA ../../../../data/int/D/dat22.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m22.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1359.28658096358        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.1730E+02 -2.1467E+01 -9.7525E+01 -8.7225E+01  3.6394E+02 -1.5105E+03 -3.8667E+02 -7.5878E+01 -6.9213E+02 -6.2510E+02
            -7.5877E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2355.41963396723        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.2427E+00  1.0190E+00  1.0161E+00  1.2812E+00  6.9960E-01  4.0728E+00  1.6064E+00  1.0037E+00  3.0738E+00  3.0711E+00
             2.9310E+00
 PARAMETER:  3.1726E-01  1.1885E-01  1.1596E-01  3.4779E-01 -2.5724E-01  1.5043E+00  5.7399E-01  1.0372E-01  1.2229E+00  1.2220E+00
             1.1754E+00
 GRADIENT:   2.6783E+01  1.4046E+00 -2.1497E+00  4.9033E+01 -3.9763E+01  1.4399E+02  3.0210E+01  1.0528E+01  6.5872E+01  7.3021E+01
             6.7949E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2373.89987575395        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0781E+00  1.4777E+00  3.6813E+00  1.1580E+00  1.1697E+00  3.4028E+00  4.6488E+00  5.4131E-01  2.3843E+00  3.5101E+00
             2.8780E+00
 PARAMETER:  1.7516E-01  4.9051E-01  1.4033E+00  2.4667E-01  2.5674E-01  1.3246E+00  1.6366E+00 -5.1377E-01  9.6889E-01  1.3556E+00
             1.1571E+00
 GRADIENT:  -2.2773E+00  2.7068E+01 -1.9660E+01  3.0037E+01 -5.2669E+01  9.5094E+01  4.4482E+01  9.8321E-01  5.6326E+01  1.3252E+02
             3.2982E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2453.77944468983        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0799E+00  1.2528E+00  1.7779E+00  7.7423E-01  1.3654E+00  2.6429E+00  3.5462E+00  3.4575E-01  2.2014E+00  1.6238E+00
             2.9421E+00
 PARAMETER:  1.7685E-01  3.2537E-01  6.7543E-01 -1.5588E-01  4.1147E-01  1.0719E+00  1.3659E+00 -9.6205E-01  8.8910E-01  5.8478E-01
             1.1791E+00
 GRADIENT:  -1.0588E+01 -3.7020E+01 -3.7070E+00 -1.6959E+01 -4.4797E+01 -7.3735E+00 -5.9668E+00  2.3879E-01  5.5886E+01  6.0360E+01
             9.1275E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2484.05281983408        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0875E+00  1.5734E+00  1.3309E+00  7.0474E-01  1.4285E+00  2.6864E+00  3.4097E+00  8.2917E-01  1.5278E+00  1.1244E+00
             2.8063E+00
 PARAMETER:  1.8385E-01  5.5323E-01  3.8587E-01 -2.4993E-01  4.5662E-01  1.0882E+00  1.3266E+00 -8.7329E-02  5.2382E-01  2.1723E-01
             1.1319E+00
 GRADIENT:  -6.2603E+00 -7.0673E+00 -5.7685E+00 -8.0584E+00 -4.2545E+00  3.4628E+00  2.9897E+00  4.1217E-01  8.0775E+00  4.5646E+00
            -1.2305E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2486.89644490189        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0945E+00  1.5364E+00  1.6380E+00  7.7261E-01  1.4375E+00  2.6724E+00  3.5330E+00  7.4347E-01  1.3323E+00  1.1543E+00
             2.8270E+00
 PARAMETER:  1.9029E-01  5.2941E-01  5.9349E-01 -1.5798E-01  4.6293E-01  1.0830E+00  1.3622E+00 -1.9643E-01  3.8687E-01  2.4351E-01
             1.1392E+00
 GRADIENT:  -4.2968E+00  1.2130E+00 -1.5799E+00 -5.4457E+00 -6.9460E+00  6.4933E-01  1.2835E+00 -6.6408E-01 -5.5550E-02  3.7750E+00
            -3.7632E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2491.85347308113        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  1.1182E+00  1.2841E+00  2.7266E+00  9.4166E-01  1.4918E+00  2.6762E+00  3.9357E+00  7.9208E-01  1.1855E+00  1.1572E+00
             2.8384E+00
 PARAMETER:  2.1173E-01  3.5002E-01  1.1031E+00  3.9888E-02  5.0001E-01  1.0844E+00  1.4701E+00 -1.3309E-01  2.7016E-01  2.4601E-01
             1.1432E+00
 GRADIENT:  -8.0499E+00 -1.2286E+00 -3.0616E-01  5.6234E-01 -8.8816E+00 -3.2418E+01 -3.9484E+01 -9.1175E-02 -5.0855E+00 -5.0561E+00
            -1.4885E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2499.46951245816        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:      668             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1539E+00  1.0377E+00  4.6657E+00  1.1604E+00  1.5586E+00  2.8327E+00  4.8006E+00  2.1657E-01  1.1571E+00  1.2113E+00
             2.8587E+00
 PARAMETER:  2.4311E-01  1.3696E-01  1.6402E+00  2.4880E-01  5.4378E-01  1.1412E+00  1.6687E+00 -1.4298E+00  2.4595E-01  2.9170E-01
             1.1504E+00
 GRADIENT:   1.4848E+01  1.2386E+01  6.2480E+00  1.9572E+01 -1.7996E+01  2.6262E+01  2.5883E+01  5.4967E-02 -1.3933E+01 -1.1246E+01
            -4.4546E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2499.49113437549        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  1.1539E+00  1.0377E+00  4.6653E+00  1.1604E+00  1.5586E+00  2.8325E+00  4.8000E+00  6.2454E-02  1.1572E+00  1.2113E+00
             2.8587E+00
 PARAMETER:  2.4311E-01  1.3696E-01  1.6401E+00  2.4880E-01  5.4381E-01  1.1411E+00  1.6686E+00 -2.6733E+00  2.4596E-01  2.9170E-01
             1.1504E+00
 GRADIENT:   1.5987E+00  1.1170E+01  5.8679E+00  1.4675E+01 -2.0417E+01 -7.7694E+00 -2.4159E+01  4.3153E-03 -1.4515E+01 -1.1653E+01
            -6.5449E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2502.56860978003        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  1.1487E+00  9.3113E-01  4.7046E+00  1.1593E+00  1.5734E+00  2.8944E+00  5.1969E+00  2.3502E-02  1.2875E+00  1.2623E+00
             2.8555E+00
 PARAMETER:  2.3864E-01  2.8643E-02  1.6485E+00  2.4782E-01  5.5326E-01  1.1628E+00  1.7481E+00 -3.6507E+00  3.5269E-01  3.3291E-01
             1.1492E+00
 GRADIENT:   4.5753E-01  4.3277E+00 -1.8871E+00 -1.1138E+01 -1.0620E+00  6.2085E-01 -2.4524E+00  1.1043E-03 -7.5744E-01 -4.5209E-01
             3.4113E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2502.62383141041        NO. OF FUNC. EVALS.: 215
 CUMULATIVE NO. OF FUNC. EVALS.:     1158             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1489E+00  9.2723E-01  4.7246E+00  1.1619E+00  1.5731E+00  2.8952E+00  5.2061E+00  2.2935E-02  1.2875E+00  1.2620E+00
             2.8552E+00
 PARAMETER:  2.3876E-01  2.4444E-02  1.6528E+00  2.5008E-01  5.5302E-01  1.1631E+00  1.7498E+00 -3.6751E+00  3.5268E-01  3.3273E-01
             1.1491E+00
 GRADIENT:   1.3548E+01  5.2528E+00 -1.2482E+00 -5.9478E+00  9.2135E-01  3.5603E+01  5.4452E+01  1.0863E-03 -5.9514E-02 -3.3020E-01
             2.1327E+00

0ITERATION NO.:   53    OBJECTIVE VALUE:  -2502.62426773359        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.1489E+00  9.2714E-01  4.7241E+00  1.1619E+00  1.5731E+00  2.8955E+00  5.2056E+00  2.2850E-02  1.2875E+00  1.2620E+00
             2.8550E+00
 PARAMETER:  2.3876E-01  2.4344E-02  1.6528E+00  2.5008E-01  5.5302E-01  1.1631E+00  1.7498E+00 -3.6751E+00  3.5268E-01  3.3273E-01
             1.1491E+00
 GRADIENT:  -4.1600E+04 -9.9337E+04  6.0088E+03  3.9713E+04 -8.9844E+03 -8.5129E+03  5.6731E+03  1.2275E-03 -1.4085E+04  2.9852E+04
             8.6008E+03
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         3.3         3.3         2.1         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1254
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.0327E-04  2.1344E-02 -3.8196E-04 -5.0303E-02  1.1610E-02
 SE:             2.9801E-02  2.5190E-02  1.2135E-04  1.8112E-02  2.5160E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8920E-01  3.9681E-01  1.6466E-03  5.4815E-03  6.4448E-01

 ETASHRINKSD(%)  1.6315E-01  1.5612E+01  9.9593E+01  3.9321E+01  1.5710E+01
 ETASHRINKVR(%)  3.2603E-01  2.8786E+01  9.9998E+01  6.3181E+01  2.8953E+01
 EBVSHRINKSD(%)  2.2691E-01  1.2407E+01  9.9554E+01  4.0064E+01  1.6426E+01
 EBVSHRINKVR(%)  4.5330E-01  2.3275E+01  9.9998E+01  6.4077E+01  3.0155E+01
 RELATIVEINF(%)  9.9542E+01  3.1590E+01  5.3368E-04  1.1133E+01  3.0632E+01
 EPSSHRINKSD(%)  1.5946E+01
 EPSSHRINKVR(%)  2.9349E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2502.6242677335886     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -848.53490796517781     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.79
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2502.624       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.15E+00  9.27E-01  4.72E+00  1.16E+00  1.57E+00  2.90E+00  5.21E+00  2.29E-02  1.29E+00  1.26E+00  2.86E+00
 


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
+        3.30E+07
 
 TH 2
+        2.07E+02  2.89E+08
 
 TH 3
+       -1.96E+00 -5.74E+00  4.07E+04
 
 TH 4
+        6.24E+00 -1.60E+02 -1.16E+01  2.94E+07
 
 TH 5
+        3.05E+03  5.07E+01 -1.23E+02  7.04E+01  3.28E+06
 
 TH 6
+       -7.67E+02  1.72E+01 -1.41E-01  6.55E-01  2.49E+02  2.18E+05
 
 TH 7
+        2.83E+02 -3.00E+00 -3.66E-01 -9.95E+00 -8.85E+01  3.88E+01  2.99E+04
 
 TH 8
+       -1.18E+00 -1.52E+01  4.68E-03  5.96E+00 -9.62E-01 -5.07E-01 -1.12E-01 -2.16E+01
 
 TH 9
+        2.38E-01  1.24E+02 -2.11E+00 -1.39E+01  1.05E+01 -1.32E-01  2.16E+00 -2.85E+00  1.20E+07
 
 TH10
+        5.22E+03 -1.40E+02 -1.85E+02 -4.93E+03  1.63E+03  4.25E+02 -1.57E+02 -3.59E+00  3.16E+03  1.41E+07
 
 TH11
+        7.84E+02 -1.93E+01 -3.59E-02 -7.12E+00 -2.55E+02  2.30E+03 -3.93E+01  3.42E-01  4.08E+00 -4.28E+02  2.28E+05
 
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
 #CPUT: Total CPU Time in Seconds,       50.030
Stop Time:
Sat Sep 25 05:31:27 CDT 2021
