Sat Sep 18 10:38:32 CDT 2021
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
$DATA ../../../../data/spa/A3/dat69.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   542.447253187372        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3077E+01  9.9466E+01  1.0504E+02  2.8340E+00  1.9139E+02  4.6517E+00 -6.1423E+01 -4.2462E+01 -1.2820E+02 -2.2820E+02
            -3.9536E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1199.97493330508        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0753E+00  9.5677E-01  9.2214E-01  1.1188E+00  9.3559E-01  8.5186E-01  9.6170E-01  9.8941E-01  9.9970E-01  1.0622E+00
             5.3720E+00
 PARAMETER:  1.7256E-01  5.5809E-02  1.8946E-02  2.1223E-01  3.3419E-02 -6.0338E-02  6.0944E-02  8.9357E-02  9.9703E-02  1.6032E-01
             1.7812E+00
 GRADIENT:   8.4782E+00 -1.0216E+01 -1.7223E+01  1.0992E+00  1.2994E+00 -9.7177E+00  1.0893E+01  5.9484E+00  2.0494E+01  2.2919E+01
             1.4376E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1217.01286850040        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0696E+00  6.9797E-01  3.6684E-01  1.1884E+00  4.3799E-01  9.2415E-01  5.1260E-01  2.2157E-01  9.1117E-01  5.2897E-01
             4.9287E+00
 PARAMETER:  1.6731E-01 -2.5958E-01 -9.0282E-01  2.7262E-01 -7.2556E-01  2.1118E-02 -5.6826E-01 -1.4070E+00  6.9698E-03 -5.3683E-01
             1.6951E+00
 GRADIENT:  -2.4871E+01  4.6631E+01  5.6447E-01  8.8758E+01 -4.0868E+01 -3.5983E+00  8.3253E-01  7.2034E-01 -5.5074E-01  1.1061E+01
             9.3730E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1233.26761071955        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.0504E+00  5.1365E-01  4.4192E-01  1.2330E+00  4.4859E-01  9.2490E-01  6.6900E-01  6.8753E-02  9.3184E-01  2.6912E-01
             4.2210E+00
 PARAMETER:  1.4915E-01 -5.6621E-01 -7.1663E-01  3.0943E-01 -7.0165E-01  2.1931E-02 -3.0197E-01 -2.5772E+00  2.9404E-02 -1.2126E+00
             1.5401E+00
 GRADIENT:  -7.5015E+00  1.5007E+01  1.0457E+01  2.4478E+01 -2.1186E+01 -1.3709E+00 -4.8533E-01  5.7797E-02 -2.4272E+00  3.4140E-01
            -1.2145E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1234.65627016083        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0477E+00  2.8977E-01  3.7359E-01  1.2791E+00  3.6033E-01  9.3104E-01  1.3446E+00  1.3607E-02  9.1333E-01  5.6216E-02
             4.2471E+00
 PARAMETER:  1.4659E-01 -1.1387E+00 -8.8459E-01  3.4619E-01 -9.2074E-01  2.8548E-02  3.9612E-01 -4.1972E+00  9.3406E-03 -2.7786E+00
             1.5462E+00
 GRADIENT:  -3.6685E+00  2.8714E+00 -2.6412E-01  1.2075E+01 -4.2813E+00 -4.0182E-01 -2.0183E-01  1.6470E-03 -2.2222E+00 -6.1067E-02
            -1.6062E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1234.74815790959        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.0467E+00  2.3943E-01  3.6202E-01  1.2795E+00  3.4440E-01  9.3136E-01  1.6967E+00  1.0000E-02  9.2449E-01  3.0302E-02
             4.2428E+00
 PARAMETER:  1.4568E-01 -1.3295E+00 -9.1607E-01  3.4643E-01 -9.6595E-01  2.8895E-02  6.2867E-01 -4.8043E+00  2.1491E-02 -3.3966E+00
             1.5452E+00
 GRADIENT:   1.8614E-01  1.6721E+00  3.5342E+00 -6.2034E-01 -6.1134E+00  8.3025E-02  2.3946E-01  0.0000E+00 -6.2691E-02 -2.0924E-02
             2.9904E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1234.92595513088        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  1.0403E+00  1.3740E-01  3.6159E-01  1.3145E+00  3.3318E-01  9.2932E-01  2.7885E+00  1.0000E-02  9.0796E-01  1.0000E-02
             4.2402E+00
 PARAMETER:  1.3952E-01 -1.8849E+00 -9.1725E-01  3.7344E-01 -9.9906E-01  2.6701E-02  1.1255E+00 -6.6688E+00  3.4412E-03 -5.1657E+00
             1.5446E+00
 GRADIENT:   4.4302E-01 -3.4592E-02  4.2367E-01 -1.8530E+00 -5.9203E-01  1.7630E-01 -7.0240E-02  0.0000E+00  5.7019E-01  0.0000E+00
             8.5424E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1235.13269598841        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  1.0412E+00  1.1697E-01  4.1882E-01  1.3676E+00  3.6734E-01  9.2328E-01  3.3409E+00  1.0000E-02  8.6542E-01  1.0000E-02
             4.2796E+00
 PARAMETER:  1.4033E-01 -2.0459E+00 -7.7032E-01  4.1305E-01 -9.0146E-01  2.0180E-02  1.3062E+00 -7.3111E+00 -4.4542E-02 -5.6360E+00
             1.5539E+00
 GRADIENT:   4.1516E+00 -2.1232E-01  3.1759E+00 -3.4426E+00 -3.1901E+00  3.8755E-01  7.2304E-02  0.0000E+00  1.6108E+00  0.0000E+00
             2.0719E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1235.20278072892        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      761
 NPARAMETR:  1.0417E+00  1.6520E-01  4.4137E-01  1.3666E+00  3.8833E-01  9.2136E-01  2.5530E+00  1.0000E-02  8.4982E-01  1.0212E-02
             4.2846E+00
 PARAMETER:  1.4084E-01 -1.7006E+00 -7.1787E-01  4.1230E-01 -8.4589E-01  1.8091E-02  1.0373E+00 -6.1750E+00 -6.2733E-02 -4.4842E+00
             1.5550E+00
 GRADIENT:  -1.6556E-01  5.5650E-02  1.3119E-01  5.4990E-01 -4.1253E-01 -2.1586E-02  1.5804E-02  0.0000E+00 -8.6569E-02 -4.2174E-04
            -9.2273E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1235.20313721043        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  1.0418E+00  1.6871E-01  4.4436E-01  1.3666E+00  3.9078E-01  9.2115E-01  2.5015E+00  1.0000E-02  8.4905E-01  1.0955E-02
             4.2862E+00
 PARAMETER:  1.4099E-01 -1.6796E+00 -7.1112E-01  4.1234E-01 -8.3961E-01  1.7872E-02  1.0169E+00 -6.1067E+00 -6.3633E-02 -4.4140E+00
             1.5554E+00
 GRADIENT:   2.7081E-02 -8.2131E-03 -3.1337E-02 -6.2848E-02  6.5397E-02 -1.7348E-03 -4.0961E-03  0.0000E+00 -2.9044E-03 -4.0518E-04
             2.9129E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1235.22040456029        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1031
 NPARAMETR:  1.0408E+00  1.6831E-01  4.4433E-01  1.3650E+00  3.8996E-01  9.2092E-01  2.4775E+00  1.0000E-02  8.4943E-01  1.5346E-01
             4.2705E+00
 PARAMETER:  1.3997E-01 -1.6820E+00 -7.1120E-01  4.1119E-01 -8.4171E-01  1.7616E-02  1.0073E+00 -6.1067E+00 -6.3193E-02 -1.7743E+00
             1.5517E+00
 GRADIENT:   4.5261E-01  9.2163E-03  4.6381E-01  6.9014E-01  2.7348E+00 -3.7777E-02  5.2243E-02  0.0000E+00 -1.6032E-01  1.4424E-02
             1.8347E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1235.22171577815        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  1.0400E+00  1.6725E-01  4.4099E-01  1.3622E+00  3.8704E-01  9.2147E-01  2.4399E+00  1.0000E-02  8.5438E-01  1.8266E-01
             4.2478E+00
 PARAMETER:  1.3922E-01 -1.6883E+00 -7.1872E-01  4.0907E-01 -8.4922E-01  1.8211E-02  9.9197E-01 -6.1067E+00 -5.7382E-02 -1.6002E+00
             1.5464E+00
 GRADIENT:  -7.7368E-01 -4.3637E-03  1.8870E+00 -1.9226E-01  1.6507E+00 -7.8228E-02 -2.1573E-02  0.0000E+00 -1.5242E-01  3.3366E-03
            -6.1326E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1235.22258196383        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1173
 NPARAMETR:  1.0391E+00  1.6413E-01  4.2847E-01  1.3538E+00  3.7728E-01  9.2336E-01  2.3638E+00  1.0000E-02  8.6766E-01  2.4771E-01
             4.2018E+00
 PARAMETER:  1.3836E-01 -1.7071E+00 -7.4753E-01  4.0291E-01 -8.7476E-01  2.0267E-02  9.6028E-01 -6.1067E+00 -4.1953E-02 -1.2955E+00
             1.5355E+00
 GRADIENT:  -2.5739E+00 -8.4777E-03  3.8612E+00 -1.4115E+00  5.9417E-01 -1.4718E-01 -1.1381E-01  0.0000E+00 -1.1806E-01  2.1993E-03
            -3.7864E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1235.27119535304        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1354
 NPARAMETR:  1.0453E+00  1.5907E-01  4.0223E-01  1.3463E+00  3.5819E-01  9.3147E-01  2.1792E+00  1.0000E-02  8.9625E-01  2.9594E-01
             4.1595E+00
 PARAMETER:  1.4427E-01 -1.7384E+00 -8.1073E-01  3.9734E-01 -9.2668E-01  2.9004E-02  8.7894E-01 -6.1067E+00 -9.5346E-03 -1.1176E+00
             1.5254E+00
 GRADIENT:   5.0363E+00  3.4917E-01  5.2013E+00  6.0751E+00 -5.6247E+00  7.8700E-01 -2.8825E-01  0.0000E+00  9.1706E-01 -2.6341E-01
            -6.9964E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1235.44386229157        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1534
 NPARAMETR:  1.0433E+00  1.5064E-01  3.6491E-01  1.3157E+00  3.3367E-01  9.3166E-01  2.5243E+00  1.0000E-02  9.0973E-01  3.1510E-01
             4.1617E+00
 PARAMETER:  1.4239E-01 -1.7929E+00 -9.0809E-01  3.7434E-01 -9.9760E-01  2.9215E-02  1.0260E+00 -6.1067E+00  5.3886E-03 -1.0549E+00
             1.5259E+00
 GRADIENT:  -8.9210E-01  6.3183E-02 -3.9416E-02 -1.4854E+00  3.1322E-01 -8.2320E-02  3.6890E-02  0.0000E+00 -1.6855E-01 -6.5312E-02
             4.0901E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1235.45120650245        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1711
 NPARAMETR:  1.0427E+00  1.3646E-01  3.6431E-01  1.3212E+00  3.3161E-01  9.3179E-01  2.7420E+00  1.0000E-02  9.0860E-01  3.2540E-01
             4.1546E+00
 PARAMETER:  1.4185E-01 -1.8918E+00 -9.0974E-01  3.7856E-01 -1.0038E+00  2.9351E-02  1.1087E+00 -6.1067E+00  4.1541E-03 -1.0227E+00
             1.5242E+00
 GRADIENT:  -3.9093E-01  2.6570E-02  1.0286E-01  1.3242E-03 -2.0730E-01 -4.7634E-02  6.1350E-03  0.0000E+00 -6.9015E-02 -8.6595E-03
             8.9218E-02

0ITERATION NO.:   77    OBJECTIVE VALUE:  -1235.45132460081        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1768
 NPARAMETR:  1.0428E+00  1.3528E-01  3.6432E-01  1.3217E+00  3.3150E-01  9.3194E-01  2.7596E+00  1.0000E-02  9.0885E-01  3.2666E-01
             4.1535E+00
 PARAMETER:  1.4195E-01 -1.9004E+00 -9.0974E-01  3.7891E-01 -1.0041E+00  2.9514E-02  1.1151E+00 -6.1067E+00  4.4224E-03 -1.0188E+00
             1.5240E+00
 GRADIENT:   3.3482E-02 -8.2686E-04  5.6931E-03 -6.9151E-03  2.3618E-03  3.1558E-03 -1.0053E-03  0.0000E+00  5.7146E-03  9.2032E-04
             5.7582E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1768
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2324E-03 -1.4449E-03  1.1946E-04 -1.6621E-02  5.4976E-04
 SE:             2.8058E-02  5.2005E-03  1.9892E-04  2.4415E-02  9.5931E-03
 N:                     100         100         100         100         100

 P VAL.:         9.6497E-01  7.8113E-01  5.4815E-01  4.9600E-01  9.5430E-01

 ETASHRINKSD(%)  6.0009E+00  8.2578E+01  9.9334E+01  1.8207E+01  6.7862E+01
 ETASHRINKVR(%)  1.1642E+01  9.6965E+01  9.9996E+01  3.3099E+01  8.9671E+01
 EBVSHRINKSD(%)  5.6670E+00  8.3611E+01  9.9284E+01  1.7871E+01  6.7746E+01
 EBVSHRINKVR(%)  1.1013E+01  9.7314E+01  9.9995E+01  3.2549E+01  8.9597E+01
 RELATIVEINF(%)  6.7048E+01  1.8518E-01  1.5806E-04  7.8034E+00  2.7710E-01
 EPSSHRINKSD(%)  2.1448E+01
 EPSSHRINKVR(%)  3.8295E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1235.4513246008119     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -500.30049803707368     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1235.451       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.35E-01  3.64E-01  1.32E+00  3.32E-01  9.32E-01  2.76E+00  1.00E-02  9.09E-01  3.27E-01  4.15E+00
 


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
+        1.05E+03
 
 TH 2
+       -1.31E+02  3.99E+02
 
 TH 3
+       -5.06E+01  4.40E+02  4.01E+03
 
 TH 4
+       -9.16E+01  2.06E+02 -2.20E+02  5.16E+02
 
 TH 5
+        2.75E+02 -1.02E+03 -5.73E+03 -3.43E+02  9.45E+03
 
 TH 6
+       -2.67E+00 -9.96E+00  3.43E+01 -1.95E+01  9.20E+00  1.79E+02
 
 TH 7
+        6.25E-02  1.08E+01 -2.54E+00 -9.41E-01  3.49E+00  1.84E-01  6.58E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.89E+00 -2.75E+01  3.07E+01 -1.25E+01  8.31E+01 -9.90E-01  7.00E-01  0.00E+00  1.02E+02
 
 TH10
+       -1.09E+01 -1.34E+01 -1.03E+02 -9.92E+00  1.84E+02 -1.71E+00  3.09E-01  0.00E+00  1.49E+00  2.57E+01
 
 TH11
+       -1.83E+01 -6.84E+00 -1.22E+01 -8.36E+00  1.70E+01  3.65E+00 -6.52E-02  0.00E+00  9.39E+00  1.20E+01  2.40E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.001
Stop Time:
Sat Sep 18 10:39:02 CDT 2021
