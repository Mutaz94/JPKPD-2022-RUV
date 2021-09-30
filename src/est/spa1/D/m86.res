Thu Sep 30 03:42:36 CDT 2021
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
$DATA ../../../../data/spa1/D/dat86.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25145.9418211803        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.9560E+02  6.6039E+02 -7.2115E+01  5.3268E+02  1.8505E+02 -2.6998E+03 -1.1967E+03 -2.5244E+01 -1.9862E+03 -6.1197E+02
            -4.7319E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -475.063727491675        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1973E+00  9.5688E-01  8.8724E-01  1.7716E+00  1.4768E+00  3.1713E+00  1.2870E+00  9.0098E-01  1.7144E+00  8.9807E-01
             1.3423E+01
 PARAMETER:  2.8004E-01  5.5925E-02 -1.9639E-02  6.7188E-01  4.8990E-01  1.2541E+00  3.5229E-01 -4.2717E-03  6.3906E-01 -7.5050E-03
             2.6970E+00
 GRADIENT:  -8.1535E+00  4.2328E+01 -1.4981E+01  6.5937E+01  8.8484E-01  1.1477E+02 -1.5294E+00  3.9041E+00 -5.6779E+01  8.5209E-01
             3.7990E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -512.393363341843        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2588E+00  1.1389E+00  1.1846E+00  1.6841E+00  4.7718E+00  2.4070E+00  9.1049E-01  2.3528E-01  2.8444E+00  3.1223E-01
             1.4112E+01
 PARAMETER:  3.3016E-01  2.3007E-01  2.6943E-01  6.2125E-01  1.6627E+00  9.7837E-01  6.2298E-03 -1.3470E+00  1.1453E+00 -1.0640E+00
             2.7471E+00
 GRADIENT:  -1.2667E+01  2.3176E+01  4.1303E+00  3.8681E+01 -1.0210E+01  1.2736E+01  4.1081E+00  1.6583E-01  1.4487E+01  5.7513E-02
             1.2591E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -530.255561261360        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      258
 NPARAMETR:  1.1706E+00  9.3015E-01  9.0970E-01  1.4792E+00  7.4225E+00  2.2029E+00  6.2143E-01  1.0349E-01  2.6895E+00  8.9387E+00
             1.2613E+01
 PARAMETER:  2.5750E-01  2.7586E-02  5.3590E-03  4.9147E-01  2.1045E+00  8.8977E-01 -3.7574E-01 -2.1683E+00  1.0894E+00  2.2904E+00
             2.6347E+00
 GRADIENT:  -2.5081E+01  1.2684E+01  4.6663E+00  1.9302E+01 -8.3589E+00 -3.6954E+01  2.3620E+00  1.4420E-02  1.5744E+01  1.3177E+00
             4.8555E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -538.188375851956        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      358
 NPARAMETR:  1.1208E+00  6.4176E-01  6.6207E-01  1.3381E+00  7.7907E+00  2.3196E+00  5.4660E-02  1.0000E-02  2.6675E+00  8.7007E+00
             1.1426E+01
 PARAMETER:  2.1400E-01 -3.4353E-01 -3.1238E-01  3.9127E-01  2.1529E+00  9.4140E-01 -2.8066E+00 -4.8216E+00  1.0812E+00  2.2634E+00
             2.5359E+00
 GRADIENT:   6.8010E+00  1.5312E+00 -3.3044E+00  1.8660E+00 -1.5973E+01  1.2214E+01  2.0120E-02  0.0000E+00  5.1677E+01  7.3883E+00
             1.8096E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -538.424625247867        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  9.8904E-01  3.4948E-01  4.5052E-01  1.3069E+00  8.9919E+00  2.1479E+00  1.0000E-02  1.0000E-02  2.5362E+00  8.0093E+00
             1.0425E+01
 PARAMETER:  8.8980E-02 -9.5130E-01 -6.9734E-01  3.6765E-01  2.2963E+00  8.6448E-01 -5.1394E+00 -7.2359E+00  1.0307E+00  2.1806E+00
             2.4442E+00
 GRADIENT:  -7.8993E+00  5.0070E+00  3.3524E+00  9.6857E+00 -2.1160E+01 -8.9982E+00  0.0000E+00  0.0000E+00  8.2517E+01  1.4117E+01
            -6.6371E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -544.043954621585        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.0109E+00  3.3569E-01  4.8567E-01  1.4262E+00  1.0296E+01  2.3270E+00  1.0000E-02  1.0000E-02  2.3850E+00  7.4604E+00
             1.0269E+01
 PARAMETER:  1.1088E-01 -9.9157E-01 -6.2223E-01  4.5501E-01  2.4317E+00  9.4458E-01 -6.2260E+00 -8.4726E+00  9.6918E-01  2.1096E+00
             2.4291E+00
 GRADIENT:  -6.8152E+00 -1.8236E+00  1.1854E+01  1.9058E+01 -7.9674E+00 -5.1024E+00  0.0000E+00  0.0000E+00  6.4156E+01  6.0531E+00
            -1.1132E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -577.737662003422        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      766
 NPARAMETR:  7.0018E-01  1.1994E-01  1.4935E-01  9.0578E-01  7.4823E+01  1.7071E+00  1.0000E-02  1.0000E-02  1.0640E+00  2.6285E+00
             9.5394E+00
 PARAMETER: -2.5641E-01 -2.0208E+00 -1.8015E+00  1.0385E-03  4.4151E+00  6.3482E-01 -2.4013E+01 -3.0784E+01  1.6202E-01  1.0664E+00
             2.3554E+00
 GRADIENT:  -8.9581E+01  3.8431E+01  2.5494E+01  3.1378E+01 -4.8116E-01 -9.1203E+01  0.0000E+00  0.0000E+00 -2.0271E+01  1.0489E-02
            -2.1004E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -648.202285026008        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      946
 NPARAMETR:  5.2349E-01  1.2673E-02  3.4743E-02  3.7077E-01  7.7322E+02  1.7938E+00  1.0000E-02  1.0000E-02  7.5707E-01  6.8101E-01
             1.1531E+01
 PARAMETER: -5.4724E-01 -4.2682E+00 -3.2598E+00 -8.9217E-01  6.7506E+00  6.8434E-01 -2.8737E+01 -3.2044E+01 -1.7830E-01 -2.8418E-01
             2.5450E+00
 GRADIENT:   9.6222E+00  1.8964E-02 -3.7054E+01  5.2491E+01  3.9868E-03  7.7171E+00  0.0000E+00  0.0000E+00 -8.8577E+00  1.9202E-09
            -5.1213E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -651.314221373562        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1133
 NPARAMETR:  4.5040E-01  1.0000E-02  2.4154E-02  2.7280E-01  6.6238E+02  1.7165E+00  1.0000E-02  1.0000E-02  8.1510E-01  6.7485E-01
             1.1296E+01
 PARAMETER: -6.9762E-01 -4.7646E+00 -3.6233E+00 -1.1990E+00  6.5958E+00  6.4028E-01 -2.7429E+01 -2.9106E+01 -1.0445E-01 -2.9327E-01
             2.5245E+00
 GRADIENT:   4.7585E+00  0.0000E+00 -5.7962E+00  3.9918E+00  1.7417E-03  4.6709E-01  0.0000E+00  0.0000E+00 -4.3936E-01  2.9036E-08
            -4.1321E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -651.334736104446        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1329             RESET HESSIAN, TYPE I
 NPARAMETR:  4.4831E-01  1.0000E-02  2.4174E-02  2.7235E-01  4.1395E+02  1.7143E+00  1.0000E-02  1.0000E-02  8.1643E-01  6.6870E-01
             1.1333E+01
 PARAMETER: -7.0228E-01 -4.7646E+00 -3.6225E+00 -1.2007E+00  6.1257E+00  6.3898E-01 -2.7429E+01 -2.9106E+01 -1.0282E-01 -3.0243E-01
             2.5277E+00
 GRADIENT:   5.2117E+01  0.0000E+00  7.1666E+01  2.6029E+01  2.0223E-03  1.2124E+01  0.0000E+00  0.0000E+00  2.3288E-01  8.8867E-08
             2.5324E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -651.339516604519        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1526
 NPARAMETR:  4.4798E-01  1.0000E-02  2.4129E-02  2.7192E-01  3.2479E+02  1.7140E+00  1.0000E-02  1.0000E-02  8.1653E-01  6.5978E-01
             1.1334E+01
 PARAMETER: -7.0300E-01 -4.7646E+00 -3.6244E+00 -1.2023E+00  5.8832E+00  6.3881E-01 -2.7429E+01 -2.9106E+01 -1.0269E-01 -3.1586E-01
             2.5278E+00
 GRADIENT:   5.7860E-01  0.0000E+00 -6.0634E-01 -6.1770E-01  2.6144E-03  2.0883E-01  0.0000E+00  0.0000E+00  1.0639E-01  1.1935E-07
             2.0007E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -651.344612056100        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1724
 NPARAMETR:  4.4787E-01  1.0000E-02  2.4064E-02  2.7173E-01  2.4342E+02  1.7137E+00  1.0000E-02  1.0000E-02  8.1574E-01  6.5387E-01
             1.1328E+01
 PARAMETER: -7.0326E-01 -4.7646E+00 -3.6270E+00 -1.2029E+00  5.5948E+00  6.3863E-01 -2.7429E+01 -2.9106E+01 -1.0365E-01 -3.2484E-01
             2.5273E+00
 GRADIENT:   1.1584E+00  0.0000E+00 -2.8225E+00  1.8752E+00  3.7603E-03  1.5197E-01  0.0000E+00  0.0000E+00 -6.3633E-02  2.0698E-07
            -6.6454E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -651.351800165999        NO. OF FUNC. EVALS.: 209
 CUMULATIVE NO. OF FUNC. EVALS.:     1933             RESET HESSIAN, TYPE I
 NPARAMETR:  4.4722E-01  1.0000E-02  2.4019E-02  2.7092E-01  1.7518E+02  1.7135E+00  1.0000E-02  1.0000E-02  8.1665E-01  6.4598E-01
             1.1334E+01
 PARAMETER: -7.0470E-01 -4.7646E+00 -3.6289E+00 -1.2059E+00  5.2658E+00  6.3851E-01 -2.7429E+01 -2.9106E+01 -1.0255E-01 -3.3698E-01
             2.5278E+00
 GRADIENT:   5.2229E+01  0.0000E+00  7.2276E+01  2.5684E+01  4.6401E-03  1.2104E+01  0.0000E+00  0.0000E+00  3.0820E-01  4.7893E-07
             2.5564E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -651.371029305240        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     2133             RESET HESSIAN, TYPE I
 NPARAMETR:  4.4713E-01  1.0000E-02  2.3769E-02  2.6953E-01  6.6972E+01  1.7147E+00  1.0000E-02  1.0000E-02  8.1688E-01  6.8582E-01
             1.1326E+01
 PARAMETER: -7.0490E-01 -4.7646E+00 -3.6394E+00 -1.2111E+00  4.3043E+00  6.3922E-01 -2.7429E+01 -2.9106E+01 -1.0226E-01 -2.7714E-01
             2.5271E+00
 GRADIENT:   5.4773E+01  0.0000E+00  6.8081E+01  3.0637E+01  1.4800E-02  1.2482E+01  0.0000E+00  0.0000E+00  2.1911E-01  3.6588E-06
             2.4067E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -651.385327577128        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2310
 NPARAMETR:  4.4506E-01  1.0000E-02  2.3906E-02  2.6960E-01  1.7557E+01  1.7096E+00  1.0000E-02  1.0000E-02  8.1350E-01  6.8336E-01
             1.1327E+01
 PARAMETER: -7.0954E-01 -4.7646E+00 -3.6336E+00 -1.2108E+00  2.9655E+00  6.3627E-01 -2.7429E+01 -2.9106E+01 -1.0640E-01 -2.8074E-01
             2.5272E+00
 GRADIENT:  -3.2887E+00  0.0000E+00  3.2997E+00 -4.0076E+00  1.8199E-02 -6.1762E-01  0.0000E+00  0.0000E+00 -7.3036E-02  5.7578E-05
             9.2086E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -651.406573173046        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2499             RESET HESSIAN, TYPE I
 NPARAMETR:  4.4821E-01  1.0000E-02  2.3899E-02  2.7054E-01  1.3333E+01  1.7139E+00  1.0000E-02  1.0000E-02  8.1485E-01  6.3724E-01
             1.1326E+01
 PARAMETER: -7.0250E-01 -4.7646E+00 -3.6339E+00 -1.2073E+00  2.6903E+00  6.3880E-01 -2.7429E+01 -2.9106E+01 -1.0475E-01 -3.5061E-01
             2.5271E+00
 GRADIENT:   5.2603E+01  0.0000E+00  7.0800E+01  2.8062E+01  1.3260E-02  1.2076E+01  0.0000E+00  0.0000E+00  1.1590E-01  1.4979E-04
             2.4625E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -651.409465392781        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2693
 NPARAMETR:  4.4800E-01  1.0000E-02  2.3858E-02  2.7021E-01  1.3311E+01  1.7136E+00  1.0000E-02  1.0000E-02  8.1446E-01  5.6965E-01
             1.1326E+01
 PARAMETER: -7.0297E-01 -4.7646E+00 -3.6356E+00 -1.2085E+00  2.6886E+00  6.3861E-01 -2.7429E+01 -2.9106E+01 -1.0523E-01 -4.6274E-01
             2.5271E+00
 GRADIENT:   1.1908E+00  0.0000E+00 -3.1824E+00  1.8101E+00  4.8246E-03  1.4847E-01  0.0000E+00  0.0000E+00 -9.0027E-02  1.0932E-04
            -6.9756E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -651.459541824548        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     2861
 NPARAMETR:  4.4246E-01  1.0000E-02  2.3179E-02  2.6405E-01  1.3380E+01  1.7107E+00  1.0000E-02  1.0000E-02  8.1647E-01  3.7796E-01
             1.1318E+01
 PARAMETER: -7.1540E-01 -4.7646E+00 -3.6645E+00 -1.2316E+00  2.6937E+00  6.3690E-01 -2.7429E+01 -2.9106E+01 -1.0276E-01 -8.7297E-01
             2.5264E+00
 GRADIENT:   3.8837E-01  0.0000E+00 -2.9308E+00  1.2474E+00 -1.0958E-02  1.8756E-01  0.0000E+00  0.0000E+00  2.4848E-01  5.3240E-05
            -5.2132E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -651.463840454494        NO. OF FUNC. EVALS.: 206
 CUMULATIVE NO. OF FUNC. EVALS.:     3067             RESET HESSIAN, TYPE I
 NPARAMETR:  4.4234E-01  1.0000E-02  2.3136E-02  2.6333E-01  1.3761E+01  1.7100E+00  1.0000E-02  1.0000E-02  8.1251E-01  1.9383E-01
             1.1325E+01
 PARAMETER: -7.1568E-01 -4.7646E+00 -3.6664E+00 -1.2344E+00  2.7218E+00  6.3651E-01 -2.7429E+01 -2.9106E+01 -1.0763E-01 -1.5408E+00
             2.5270E+00
 GRADIENT:   5.3519E+01  0.0000E+00  7.3782E+01  2.6254E+01  4.1486E-03  1.2077E+01  0.0000E+00  0.0000E+00 -1.6945E-02  1.7922E-05
             2.4990E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -651.467435708356        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     3263             RESET HESSIAN, TYPE I
 NPARAMETR:  4.4220E-01  1.0000E-02  2.3052E-02  2.6281E-01  1.3919E+01  1.7104E+00  1.0000E-02  1.0000E-02  8.1113E-01  1.5277E-01
             1.1319E+01
 PARAMETER: -7.1600E-01 -4.7646E+00 -3.6700E+00 -1.2363E+00  2.7333E+00  6.3673E-01 -2.7429E+01 -2.9106E+01 -1.0932E-01 -1.7788E+00
             2.5264E+00
 GRADIENT:   5.4488E+01  0.0000E+00  7.2349E+01  2.7907E+01  1.1063E-02  1.2216E+01  0.0000E+00  0.0000E+00 -2.5077E-01  1.1274E-05
             2.4118E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -651.470378816509        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3458             RESET HESSIAN, TYPE I
 NPARAMETR:  4.4150E-01  1.0000E-02  2.3034E-02  2.6233E-01  1.3835E+01  1.7095E+00  1.0000E-02  1.0000E-02  8.1295E-01  1.4028E-01
             1.1325E+01
 PARAMETER: -7.1759E-01 -4.7646E+00 -3.6708E+00 -1.2381E+00  2.7272E+00  6.3620E-01 -2.7429E+01 -2.9106E+01 -1.0709E-01 -1.8641E+00
             2.5270E+00
 GRADIENT:   5.3559E+01  0.0000E+00  7.4361E+01  2.5842E+01  2.4229E-03  1.2071E+01  0.0000E+00  0.0000E+00  6.4307E-02  9.9548E-06
             2.5174E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -651.471930267205        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3652
 NPARAMETR:  4.4123E-01  1.0000E-02  2.3006E-02  2.6209E-01  1.3735E+01  1.7093E+00  1.0000E-02  1.0000E-02  8.1425E-01  1.3098E-01
             1.1327E+01
 PARAMETER: -7.1818E-01 -4.7646E+00 -3.6720E+00 -1.2391E+00  2.7199E+00  6.3608E-01 -2.7429E+01 -2.9106E+01 -1.0549E-01 -1.9327E+00
             2.5272E+00
 GRADIENT:   4.2625E-01  0.0000E+00 -9.2937E-01 -1.4429E+00 -8.3715E-03  2.0105E-01  0.0000E+00  0.0000E+00  8.7237E-02  5.8316E-06
             2.4927E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -651.473584537905        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3847
 NPARAMETR:  4.4128E-01  1.0000E-02  2.2963E-02  2.6212E-01  1.4122E+01  1.7091E+00  1.0000E-02  1.0000E-02  8.1362E-01  1.2539E-01
             1.1321E+01
 PARAMETER: -7.1808E-01 -4.7646E+00 -3.6739E+00 -1.2389E+00  2.7478E+00  6.3597E-01 -2.7429E+01 -2.9106E+01 -1.0627E-01 -1.9763E+00
             2.5267E+00
 GRADIENT:   1.1648E+00  0.0000E+00 -3.5523E+00  1.5024E+00  5.1961E-03  1.4560E-01  0.0000E+00  0.0000E+00 -8.0462E-02  4.7089E-06
            -7.2237E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -651.475221483024        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     4042
 NPARAMETR:  4.4084E-01  1.0000E-02  2.2959E-02  2.6180E-01  1.4034E+01  1.7087E+00  1.0000E-02  1.0000E-02  8.1393E-01  1.2183E-01
             1.1325E+01
 PARAMETER: -7.1907E-01 -4.7646E+00 -3.6740E+00 -1.2402E+00  2.7415E+00  6.3572E-01 -2.7429E+01 -2.9106E+01 -1.0588E-01 -2.0052E+00
             2.5270E+00
 GRADIENT:   4.4901E-01  0.0000E+00 -1.7414E+00 -4.4481E-01 -1.6203E-03  1.0010E-01  0.0000E+00  0.0000E+00  2.2452E-02  4.5933E-06
            -7.3656E-02

0ITERATION NO.:  121    OBJECTIVE VALUE:  -651.475221483024        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     4064
 NPARAMETR:  4.4084E-01  1.0000E-02  2.2959E-02  2.6180E-01  1.4034E+01  1.7087E+00  1.0000E-02  1.0000E-02  8.1393E-01  1.2183E-01
             1.1325E+01
 PARAMETER: -7.1907E-01 -4.7646E+00 -3.6740E+00 -1.2402E+00  2.7415E+00  6.3572E-01 -2.7429E+01 -2.9106E+01 -1.0588E-01 -2.0052E+00
             2.5270E+00
 GRADIENT:   4.4901E-01  0.0000E+00 -1.7414E+00 -4.4481E-01 -1.6203E-03  1.0010E-01  0.0000E+00  0.0000E+00  2.2452E-02  4.5933E-06
            -7.3656E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4064
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4439E-03  8.8793E-07  1.2973E-04 -2.1784E-02 -3.5227E-06
 SE:             2.8978E-02  9.8004E-07  2.5035E-04  2.2849E-02  8.4200E-06
 N:                     100         100         100         100         100

 P VAL.:         9.3279E-01  3.6493E-01  6.0433E-01  3.4040E-01  6.7567E-01

 ETASHRINKSD(%)  2.9199E+00  9.9997E+01  9.9161E+01  2.3453E+01  9.9972E+01
 ETASHRINKVR(%)  5.7545E+00  1.0000E+02  9.9993E+01  4.1405E+01  1.0000E+02
 EBVSHRINKSD(%)  2.7113E+00  9.9995E+01  9.9243E+01  2.4968E+01  9.9969E+01
 EBVSHRINKVR(%)  5.3490E+00  1.0000E+02  9.9994E+01  4.3702E+01  1.0000E+02
 RELATIVEINF(%)  3.2097E+00  6.1139E-08  4.1860E-05  3.6883E-01  4.6965E-07
 EPSSHRINKSD(%)  8.9893E+00
 EPSSHRINKVR(%)  1.7171E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -651.47522148302403     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       267.46331172164867     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    88.23
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -651.475       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.41E-01  1.00E-02  2.30E-02  2.62E-01  1.40E+01  1.71E+00  1.00E-02  1.00E-02  8.14E-01  1.22E-01  1.13E+01
 


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
+        1.84E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.26E+04  0.00E+00  1.11E+06
 
 TH 4
+       -1.09E+02  0.00E+00 -1.12E+05  1.25E+04
 
 TH 5
+        2.08E-01  0.00E+00 -3.41E+00  2.83E-01  7.12E-04
 
 TH 6
+        5.62E+00  0.00E+00 -6.10E+01 -3.15E+01  3.85E-03  5.95E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.18E+01  0.00E+00  5.62E+02 -6.31E+01 -7.56E-03 -3.20E+00  0.00E+00  0.00E+00  7.77E+01
 
 TH10
+        1.35E-04  0.00E+00 -7.40E-03 -2.79E-03  9.17E-06  7.56E-04  0.00E+00  0.00E+00 -2.25E-03 -3.39E-03
 
 TH11
+       -1.96E+01  0.00E+00  3.70E+02 -2.78E+01 -3.01E-03  1.49E+00  0.00E+00  0.00E+00  3.69E+00 -6.87E-05  3.37E+00
 
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
 #CPUT: Total CPU Time in Seconds,       98.210
Stop Time:
Thu Sep 30 03:44:15 CDT 2021
