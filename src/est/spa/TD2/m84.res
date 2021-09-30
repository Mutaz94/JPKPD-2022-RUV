Wed Sep 29 19:21:37 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat84.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.53954972235        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6061E+02 -1.2639E+01 -4.3997E+01  4.0129E+01  3.5424E+01  8.5394E+01 -8.1694E+00  1.5053E+01 -3.8006E+00  9.3283E+00
            -2.3265E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1711.67048072361        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0346E+00  1.0817E+00  1.2136E+00  1.0020E+00  1.1025E+00  7.9877E-01  1.0701E+00  8.9321E-01  1.0674E+00  9.6043E-01
             1.0903E+00
 PARAMETER:  1.3397E-01  1.7853E-01  2.9356E-01  1.0204E-01  1.9754E-01 -1.2469E-01  1.6778E-01 -1.2932E-02  1.6525E-01  5.9627E-02
             1.8646E-01
 GRADIENT:  -5.1872E+00  1.1489E+01 -7.6753E+00  2.3269E+01  1.5193E+01 -1.9856E+01  7.3238E-01  3.0054E+00  1.0495E+00 -1.3735E+01
             5.6921E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1712.72104802666        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0342E+00  1.0398E+00  1.3797E+00  1.0145E+00  1.1455E+00  8.1279E-01  1.0147E+00  7.5032E-01  1.1095E+00  1.1240E+00
             1.0843E+00
 PARAMETER:  1.3358E-01  1.3907E-01  4.2186E-01  1.1439E-01  2.3581E-01 -1.0728E-01  1.1463E-01 -1.8725E-01  2.0390E-01  2.1687E-01
             1.8093E-01
 GRADIENT:  -2.4209E+00  6.6362E-01  9.2860E+00 -2.8233E+00 -3.7005E+00 -1.1908E+01  1.8589E+00 -3.4241E+00  6.8929E+00  5.1506E-01
             5.0494E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1713.58171225938        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0339E+00  8.9646E-01  1.4294E+00  1.1144E+00  1.1029E+00  8.3720E-01  1.1222E+00  9.2069E-01  9.9826E-01  1.0962E+00
             1.0658E+00
 PARAMETER:  1.3336E-01 -9.3054E-03  4.5724E-01  2.0834E-01  1.9791E-01 -7.7688E-02  2.1526E-01  1.7370E-02  9.8261E-02  1.9183E-01
             1.6374E-01
 GRADIENT:   4.7665E-01  4.7991E+00 -6.3210E-02  5.8477E+00 -2.4915E+00  2.4218E-01 -1.2699E-01  4.6334E-01 -5.5270E-01  3.0248E-01
             3.1506E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1713.76515780597        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0327E+00  6.6922E-01  1.6220E+00  1.2618E+00  1.0817E+00  8.3070E-01  1.2095E+00  1.0053E+00  9.3093E-01  1.1205E+00
             1.0654E+00
 PARAMETER:  1.3218E-01 -3.0164E-01  5.8367E-01  3.3252E-01  1.7853E-01 -8.5490E-02  2.9018E-01  1.0531E-01  2.8425E-02  2.1381E-01
             1.6338E-01
 GRADIENT:   4.2232E+00  2.8414E+00  1.6489E+00  3.7211E+00 -3.8384E+00 -1.6458E+00 -7.8719E-01 -7.6616E-01  3.0364E-01  1.1647E+00
             7.1698E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1713.82771171633        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0296E+00  5.1117E-01  1.8407E+00  1.3748E+00  1.0887E+00  8.3052E-01  1.3048E+00  1.1966E+00  8.7981E-01  1.1393E+00
             1.0626E+00
 PARAMETER:  1.2914E-01 -5.7105E-01  7.1015E-01  4.1830E-01  1.8503E-01 -8.5709E-02  3.6607E-01  2.7952E-01 -2.8045E-02  2.3042E-01
             1.6069E-01
 GRADIENT:   1.0855E-01  5.7980E+00  1.3213E+00  1.5092E+01 -5.6070E+00 -8.3952E-01 -8.8017E-01  1.5930E-01 -1.0743E+00  9.6390E-01
            -2.7988E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1713.87670621994        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.0267E+00  3.9948E-01  2.0045E+00  1.4522E+00  1.0960E+00  8.3307E-01  1.4608E+00  1.3400E+00  8.4156E-01  1.1420E+00
             1.0605E+00
 PARAMETER:  1.2632E-01 -8.1759E-01  7.9540E-01  4.7311E-01  1.9170E-01 -8.2637E-02  4.7901E-01  3.9264E-01 -7.2500E-02  2.3276E-01
             1.5872E-01
 GRADIENT:  -5.0293E+00  5.6088E+00 -6.1334E-02  1.8465E+01 -3.5442E+00  9.6654E-01 -7.5638E-01  1.0392E+00 -2.6247E+00 -2.6092E-01
            -1.3102E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1714.01551173228        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  1.0246E+00  2.8692E-01  2.1776E+00  1.5280E+00  1.1037E+00  8.3441E-01  1.7746E+00  1.4692E+00  8.0574E-01  1.1444E+00
             1.0597E+00
 PARAMETER:  1.2428E-01 -1.1486E+00  8.7823E-01  5.2394E-01  1.9868E-01 -8.1035E-02  6.7359E-01  4.8474E-01 -1.1600E-01  2.3484E-01
             1.5800E-01
 GRADIENT:  -7.7017E+00  4.2916E+00 -7.5032E-01  1.6375E+01 -9.4542E-01  2.2060E+00 -4.4759E-01  1.1637E+00 -3.4166E+00 -1.3677E+00
            -1.8247E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1714.22393612178        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1427
 NPARAMETR:  1.0239E+00  1.6299E-01  2.4202E+00  1.6152E+00  1.1172E+00  8.3246E-01  2.5713E+00  1.6204E+00  7.7305E-01  1.1578E+00
             1.0609E+00
 PARAMETER:  1.2359E-01 -1.7141E+00  9.8386E-01  5.7946E-01  2.1086E-01 -8.3372E-02  1.0444E+00  5.8267E-01 -1.5741E-01  2.4653E-01
             1.5911E-01
 GRADIENT:  -6.2764E+00  3.0759E+00 -2.3656E-02  1.8168E+01 -9.9469E-01  1.9900E+00 -2.5204E-02  4.7245E-01 -2.4240E+00 -1.3917E+00
            -1.4896E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1714.42982790567        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1602
 NPARAMETR:  1.0244E+00  9.0294E-02  2.5869E+00  1.6636E+00  1.1276E+00  8.2911E-01  3.8878E+00  1.7161E+00  7.5695E-01  1.1722E+00
             1.0624E+00
 PARAMETER:  1.2411E-01 -2.3047E+00  1.0505E+00  6.0898E-01  2.2005E-01 -8.7399E-02  1.4578E+00  6.4006E-01 -1.7845E-01  2.5885E-01
             1.6048E-01
 GRADIENT:  -2.6636E+00  1.8048E+00  5.7677E-01  1.1859E+01 -1.2967E+00  8.4277E-01  2.7301E-01  4.0511E-02 -1.0607E+00 -6.6845E-01
            -7.3615E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1714.53425022969        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1778
 NPARAMETR:  1.0245E+00  4.8605E-02  2.6990E+00  1.6913E+00  1.1358E+00  8.2763E-01  5.8913E+00  1.7850E+00  7.5004E-01  1.1817E+00
             1.0630E+00
 PARAMETER:  1.2424E-01 -2.9240E+00  1.0929E+00  6.2547E-01  2.2734E-01 -8.9186E-02  1.8735E+00  6.7940E-01 -1.8763E-01  2.6698E-01
             1.6114E-01
 GRADIENT:  -1.2045E+00  1.2090E+00  5.5995E-01  6.8557E+00 -1.0687E+00  3.7247E-01  5.7998E-01  7.8803E-02  1.0698E-01 -2.6611E-01
            -3.0974E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1714.57447177464        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1953
 NPARAMETR:  1.0246E+00  2.4860E-02  2.7781E+00  1.7085E+00  1.1421E+00  8.2676E-01  8.6412E+00  1.8328E+00  7.4737E-01  1.1897E+00
             1.0634E+00
 PARAMETER:  1.2430E-01 -3.5945E+00  1.1218E+00  6.3562E-01  2.3288E-01 -9.0242E-02  2.2565E+00  7.0585E-01 -1.9119E-01  2.7367E-01
             1.6145E-01
 GRADIENT:  -6.1739E-01  1.0252E+00  5.4907E-01  6.0052E+00 -1.3561E+00  6.4562E-02  9.8800E-01  1.2053E-01  6.5498E-01  6.7864E-02
            -1.0810E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1714.72356985027        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2134             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0254E+00  1.0000E-02  2.8140E+00  1.6989E+00  1.1485E+00  8.2668E-01  1.3322E+01  1.8592E+00  7.4219E-01  1.1905E+00
             1.0637E+00
 PARAMETER:  1.2511E-01 -4.5219E+00  1.1346E+00  6.2997E-01  2.3847E-01 -9.0332E-02  2.6894E+00  7.2013E-01 -1.9815E-01  2.7435E-01
             1.6174E-01
 GRADIENT:   4.9398E+02  9.9844E-01  5.9447E+00  1.1119E+03  1.4487E+01  2.3678E+01  1.5582E+00  1.9673E+00  2.1736E+01  1.6993E+00
             1.7481E+00

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1714.75637570053        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     2229
 NPARAMETR:  1.0252E+00  1.0000E-02  2.8127E+00  1.7081E+00  1.1458E+00  8.2642E-01  1.3321E+01  1.8571E+00  7.4243E-01  1.1923E+00
             1.0634E+00
 PARAMETER:  1.2485E-01 -4.5222E+00  1.1341E+00  6.3538E-01  2.3609E-01 -9.0648E-02  2.6893E+00  7.1904E-01 -1.9783E-01  2.7587E-01
             1.6146E-01
 GRADIENT:  -3.9482E-01  4.6203E-02  4.8788E-01  5.7786E+00  9.2943E-01  1.5825E-02  1.6116E+01  1.6195E-01 -9.5479E-02  1.1426E-01
             1.5552E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2229
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.7404E-04  5.8576E-04 -3.1238E-02 -1.0477E-02 -4.7750E-02
 SE:             2.9753E-02  2.1278E-03  1.5699E-02  2.9098E-02  2.0964E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8729E-01  7.8310E-01  4.6610E-02  7.1881E-01  2.2741E-02

 ETASHRINKSD(%)  3.2359E-01  9.2872E+01  4.7407E+01  2.5192E+00  2.9769E+01
 ETASHRINKVR(%)  6.4614E-01  9.9492E+01  7.2340E+01  4.9750E+00  5.0676E+01
 EBVSHRINKSD(%)  6.6820E-01  9.3019E+01  5.2190E+01  2.7383E+00  2.5436E+01
 EBVSHRINKVR(%)  1.3319E+00  9.9513E+01  7.7142E+01  5.4016E+00  4.4402E+01
 RELATIVEINF(%)  9.5826E+01  1.5227E-02  7.0728E+00  3.1643E+00  1.3784E+01
 EPSSHRINKSD(%)  4.3332E+01
 EPSSHRINKVR(%)  6.7887E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.7563757005296     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -979.60554913679141     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.62
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.756       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.00E-02  2.81E+00  1.71E+00  1.15E+00  8.26E-01  1.33E+01  1.86E+00  7.42E-01  1.19E+00  1.06E+00
 


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
+        1.53E+03
 
 TH 2
+       -1.01E+00  2.24E+06
 
 TH 3
+       -2.69E+00 -1.50E+04  1.07E+01
 
 TH 4
+       -1.62E+01 -4.28E+04 -9.30E+00  2.37E+03
 
 TH 5
+       -4.24E+00  1.77E+05 -4.02E+01 -3.81E+01  3.47E+02
 
 TH 6
+        2.55E-01 -1.04E+01 -2.42E-01 -3.06E+00 -9.84E-01  2.84E+02
 
 TH 7
+       -3.69E-03  1.45E+03  1.22E-01 -5.57E+01 -7.11E-01  5.12E-03  3.49E+00
 
 TH 8
+        4.41E-01 -3.58E+04 -6.77E+00 -1.32E+00 -7.15E+00  3.20E-01  6.74E-02  1.27E+01
 
 TH 9
+        4.50E+00 -3.26E+05  3.10E+00 -1.13E+00 -1.83E-01 -1.78E+00  1.56E+00  5.90E-01  3.31E+02
 
 TH10
+        4.03E-01 -1.45E+05  6.42E-01 -8.41E-01 -5.29E+01  8.95E-02  5.93E-01  3.62E+00  3.66E-01  5.29E+01
 
 TH11
+       -1.40E+01 -1.39E+01 -3.27E+00 -1.28E+01 -1.58E-01  2.90E+00 -8.36E-02  2.06E+00  1.10E+01  1.53E+01  1.95E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.432
Stop Time:
Wed Sep 29 19:22:14 CDT 2021
