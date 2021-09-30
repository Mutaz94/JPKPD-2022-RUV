Wed Sep 29 20:12:09 CDT 2021
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
$DATA ../../../../data/spa/D/dat61.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12991.1713322449        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.4061E+02  2.7939E+02 -3.1618E+01  2.0629E+02  4.1318E+02 -2.3173E+03 -6.6224E+02 -2.6199E+01 -1.1174E+03 -7.9398E+02
            -2.4062E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -619.761212653929        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4213E+00  1.1168E+00  8.3587E-01  1.8886E+00  1.3586E+00  2.1740E+00  1.2186E+00  9.6207E-01  1.1926E+00  1.3247E+00
             1.4322E+01
 PARAMETER:  4.5155E-01  2.1051E-01 -7.9288E-02  7.3585E-01  4.0643E-01  8.7655E-01  2.9769E-01  6.1332E-02  2.7615E-01  3.8120E-01
             2.7618E+00
 GRADIENT:   9.2240E+00  3.6582E+01 -7.6037E+00  8.0783E+01 -8.8899E+00  4.1036E+01 -9.9307E-01  6.6689E+00  2.9983E+00  2.2113E+00
             8.3569E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -637.940084160482        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.3145E+00  1.2153E+00  1.0148E+00  1.5139E+00  3.1327E+00  1.6840E+00  2.9965E+00  8.9620E-01  9.1508E-01  7.6330E+00
             1.2798E+01
 PARAMETER:  3.7349E-01  2.9502E-01  1.1469E-01  5.1470E-01  1.2419E+00  6.2118E-01  1.1975E+00 -9.5914E-03  1.1252E-02  2.1325E+00
             2.6493E+00
 GRADIENT:   2.4901E+01  3.8993E+01 -1.6049E+00  2.8815E+01 -1.3824E+01 -1.5175E+01  1.8317E+01  1.6007E+00  6.9905E+00  1.0156E+01
             5.0293E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -664.209507419070        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.1222E+00  4.3424E-01  1.2192E+00  1.5863E+00  3.6095E+00  1.6414E+00  2.6506E+00  2.7141E-01  5.0702E-01  8.2556E+00
             1.2276E+01
 PARAMETER:  2.1530E-01 -7.3415E-01  2.9823E-01  5.6142E-01  1.3836E+00  5.9555E-01  1.0748E+00 -1.2041E+00 -5.7921E-01  2.2109E+00
             2.6076E+00
 GRADIENT:  -2.6299E+01  9.4512E+00  2.1531E+01  2.3993E+00 -7.1365E+00  1.4658E+01  4.0286E+00  2.0338E-02  4.8391E+00  1.5795E+00
             5.3877E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -726.901275510832        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  9.1302E-01  7.2735E-02  2.0862E-01  1.1177E+00  1.1246E+01  1.6587E+00  6.4827E-01  1.0000E-02  9.4166E-02  5.8470E+00
             1.1092E+01
 PARAMETER:  9.0048E-03 -2.5209E+00 -1.4672E+00  2.1123E-01  2.5200E+00  6.0604E-01 -3.3345E-01 -8.2825E+00 -2.2627E+00  1.8659E+00
             2.5062E+00
 GRADIENT:  -1.0639E+01 -4.4429E-01  5.3530E+01 -7.0686E+00  1.4861E+01 -2.5599E+01  1.8344E-02  0.0000E+00  5.5687E-01  4.3973E+00
             3.1081E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -740.828785346808        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  9.1370E-01  7.0655E-02  2.0780E-01  1.1167E+00  1.2168E+01  1.7576E+00  2.6432E-01  1.0000E-02  3.9939E-02  5.7351E+00
             1.1165E+01
 PARAMETER:  9.7473E-03 -2.5499E+00 -1.4712E+00  2.1037E-01  2.5988E+00  6.6394E-01 -1.2306E+00 -8.2825E+00 -3.1204E+00  1.8466E+00
             2.5127E+00
 GRADIENT:  -8.0360E-01  6.5618E-01  4.4794E+01  3.1902E+01  7.2163E+00 -2.2170E+00  1.7513E-03  0.0000E+00  2.5622E-02 -3.4814E+00
             1.0575E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -742.378183065622        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:      552             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3668E-01  6.9328E-02  2.0699E-01  1.1147E+00  1.2575E+01  1.7617E+00  1.5125E-02  1.0000E-02  1.0642E-02  5.6711E+00
             1.1308E+01
 PARAMETER:  3.4587E-02 -2.5689E+00 -1.4751E+00  2.0856E-01  2.6317E+00  6.6629E-01 -4.0914E+00 -8.2825E+00 -4.4429E+00  1.8354E+00
             2.5255E+00
 GRADIENT:   1.5288E+01  7.2672E-01  4.1755E+01  3.6449E+01  5.5956E+00 -1.4179E+00  6.3447E-06  0.0000E+00  1.2087E-03 -4.2209E+00
             1.1280E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -743.107626621261        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      675             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2287E-01  6.8258E-02  2.0392E-01  1.1119E+00  1.2575E+01  1.7712E+00  1.1867E-02  1.0000E-02  1.0000E-02  5.6447E+00
             1.1319E+01
 PARAMETER:  1.9730E-02 -2.5845E+00 -1.4900E+00  2.0607E-01  2.6317E+00  6.7165E-01 -4.3340E+00 -8.2825E+00 -4.7510E+00  1.8307E+00
             2.5265E+00
 GRADIENT:   6.6609E+00  8.2624E-01  4.0718E+01  4.2934E+01  5.6012E+00 -6.7194E-01  3.8623E-06  0.0000E+00  0.0000E+00 -4.0840E+00
             1.4209E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -749.292134874301        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  8.9327E-01  6.1143E-02  1.8015E-01  1.0251E+00  1.2992E+01  1.8341E+00  1.0000E-02  1.0000E-02  1.0000E-02  5.5902E+00
             1.1140E+01
 PARAMETER: -1.2861E-02 -2.6945E+00 -1.6140E+00  1.2478E-01  2.6643E+00  7.0655E-01 -4.6841E+00 -8.2825E+00 -4.7510E+00  1.8210E+00
             2.5105E+00
 GRADIENT:   1.0956E+01  2.8209E-01  6.3748E+01 -7.8356E+00  4.1370E+00  1.5417E+01  0.0000E+00  0.0000E+00  0.0000E+00 -2.5367E+00
             1.6417E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -780.921900632024        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      918
 NPARAMETR:  6.6676E-01  1.3301E-02  3.7872E-02  3.9193E-01  2.7387E+01  1.8473E+00  1.0000E-02  1.0000E-02  1.0000E-02  5.1725E+00
             1.2131E+01
 PARAMETER: -3.0532E-01 -4.2200E+00 -3.1735E+00 -8.3668E-01  3.4101E+00  7.1373E-01 -4.6841E+00 -8.2825E+00 -4.7510E+00  1.7434E+00
             2.5958E+00
 GRADIENT:   9.2703E+01 -7.3015E-01  4.0682E+01 -9.4399E+01  1.2775E-01 -1.5482E+01  0.0000E+00  0.0000E+00  0.0000E+00 -1.8387E-03
             1.0571E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -803.893546216451        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1095
 NPARAMETR:  4.6619E-01  1.0000E-02  2.4452E-02  2.8377E-01  2.7883E+01  1.8601E+00  1.0000E-02  1.0000E-02  1.0000E-02  5.0507E+00
             9.9548E+00
 PARAMETER: -6.6316E-01 -4.8932E+00 -3.6111E+00 -1.1596E+00  3.4280E+00  7.2064E-01 -4.6841E+00 -8.2825E+00 -4.7510E+00  1.7195E+00
             2.3981E+00
 GRADIENT:   2.1094E+00  0.0000E+00  1.6611E+01 -2.5854E+01  1.1971E-03 -3.9258E-01  0.0000E+00  0.0000E+00  0.0000E+00 -1.1092E-03
             5.6097E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -803.963332776988        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1283             RESET HESSIAN, TYPE I
 NPARAMETR:  4.7924E-01  1.0000E-02  2.6426E-02  3.0521E-01  2.4490E+01  1.8773E+00  1.0000E-02  1.0000E-02  1.0000E-02  5.2513E+00
             9.9164E+00
 PARAMETER: -6.3555E-01 -4.8226E+00 -3.5334E+00 -1.0868E+00  3.2983E+00  7.2985E-01 -4.6841E+00 -8.2825E+00 -4.7510E+00  1.7585E+00
             2.3942E+00
 GRADIENT:   4.6032E+01  0.0000E+00  5.2746E+01  3.5187E+01  2.9364E-02  1.7020E+01  0.0000E+00  0.0000E+00  0.0000E+00 -7.4585E-04
             1.8867E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -803.984222933675        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1460
 NPARAMETR:  4.7929E-01  1.0000E-02  2.6556E-02  3.0488E-01  1.4857E+01  1.8724E+00  1.0000E-02  1.0000E-02  1.0000E-02  7.1633E+00
             9.9089E+00
 PARAMETER: -6.3546E-01 -4.8226E+00 -3.5285E+00 -1.0878E+00  2.7984E+00  7.2724E-01 -4.6841E+00 -8.2825E+00 -4.7510E+00  2.0690E+00
             2.3934E+00
 GRADIENT:  -9.8706E-01  0.0000E+00  3.0332E+00 -4.6771E+00  2.4107E-02  7.1554E-02  0.0000E+00  0.0000E+00  0.0000E+00 -7.2853E-03
             2.2260E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -803.996232277941        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1636
 NPARAMETR:  4.8129E-01  1.0000E-02  2.6578E-02  3.0592E-01  1.1127E+01  1.8728E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.3160E+00
             9.9133E+00
 PARAMETER: -6.3129E-01 -4.8226E+00 -3.5277E+00 -1.0844E+00  2.5094E+00  7.2745E-01 -4.6841E+00 -8.2825E+00 -4.7510E+00  1.5623E+00
             2.3939E+00
 GRADIENT:   3.7898E-01  0.0000E+00 -2.0077E+00  1.1207E+00  1.2728E-02 -9.8711E-02  0.0000E+00  0.0000E+00  0.0000E+00  6.1479E-03
            -1.1373E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -803.997082308535        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1813
 NPARAMETR:  4.8147E-01  1.0000E-02  2.6620E-02  3.0615E-01  1.0368E+01  1.8733E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.4170E+00
             9.9129E+00
 PARAMETER: -6.3091E-01 -4.8226E+00 -3.5261E+00 -1.0837E+00  2.4387E+00  7.2772E-01 -4.6841E+00 -8.2825E+00 -4.7510E+00  1.3288E+00
             2.3938E+00
 GRADIENT:   1.3007E-01  0.0000E+00 -7.7165E-01 -4.8914E-01 -6.9745E-01 -3.7482E-02  0.0000E+00  0.0000E+00  0.0000E+00  1.9964E-01
            -6.2950E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -803.997562372755        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1989
 NPARAMETR:  4.8256E-01  1.0000E-02  2.6487E-02  3.0632E-01  9.1435E+00  1.8753E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.3485E+00
             9.9365E+00
 PARAMETER: -6.2989E-01 -4.8226E+00 -3.5234E+00 -1.0818E+00  2.3364E+00  7.2801E-01 -4.6841E+00 -8.2825E+00 -4.7510E+00  9.5869E-01
             2.3938E+00
 GRADIENT:  -7.4981E-01  0.0000E+00  9.7094E+00  5.5843E+00  9.6921E+00 -2.8039E-01  0.0000E+00  0.0000E+00  0.0000E+00  5.1106E+00
            -2.7280E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1989
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3714E-03  3.4929E-06  9.8853E-05 -2.1537E-04 -1.8734E-03
 SE:             2.9268E-02  1.0719E-06  3.2531E-04  4.0607E-04  1.1036E-03
 N:                     100         100         100         100         100

 P VAL.:         9.3542E-01  1.1197E-03  7.6123E-01  5.9585E-01  8.9597E-02

 ETASHRINKSD(%)  1.9501E+00  9.9996E+01  9.8910E+01  9.8640E+01  9.6303E+01
 ETASHRINKVR(%)  3.8622E+00  1.0000E+02  9.9988E+01  9.9981E+01  9.9863E+01
 EBVSHRINKSD(%)  2.2005E+00  9.9995E+01  9.8917E+01  9.8648E+01  9.6909E+01
 EBVSHRINKVR(%)  4.3525E+00  1.0000E+02  9.9988E+01  9.9982E+01  9.9904E+01
 RELATIVEINF(%)  1.3370E+01  4.8982E-08  1.2004E-04  1.7898E-04  2.0915E-02
 EPSSHRINKSD(%)  7.7948E+00
 EPSSHRINKVR(%)  1.4982E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -803.99756237275517     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -68.846735809016991     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -803.998       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.82E-01  1.00E-02  2.67E-02  3.07E-01  9.36E+00  1.87E+00  1.00E-02  1.00E-02  1.00E-02  2.36E+00  9.91E+00
 


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
+        1.29E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -4.55E+03  0.00E+00  1.23E+06
 
 TH 4
+       -4.75E+02  0.00E+00 -9.63E+04  1.56E+04
 
 TH 5
+        9.61E+00  0.00E+00  8.82E+02 -8.62E+01  8.01E+00
 
 TH 6
+        2.38E+00  0.00E+00  4.50E+02 -3.55E+01  9.49E-01  5.11E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.59E+01  0.00E+00  3.44E+03  2.43E+02 -1.97E+01  3.55E+00  0.00E+00  0.00E+00  0.00E+00  4.34E+01
 
 TH11
+       -1.17E+01  0.00E+00 -1.85E+02 -1.13E+01  7.47E-01  3.44E-01  0.00E+00  0.00E+00  0.00E+00 -3.02E+00  3.16E+00
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       33.883
Stop Time:
Wed Sep 29 20:12:44 CDT 2021
