Sat Sep 25 10:14:35 CDT 2021
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
$DATA ../../../../data/spa/S1/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1697.88944382532        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0665E+01 -7.3161E+01 -3.3897E+01 -4.1617E+01  4.8496E+01  3.3690E+01 -1.6099E+01 -2.5697E+00  6.9089E+00  1.8333E+01
             3.0369E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1708.26483060536        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9311E-01  1.1794E+00  9.9854E-01  9.2096E-01  1.0142E+00  9.0058E-01  1.2670E+00  1.1608E+00  9.3583E-01  7.6536E-01
             9.4816E-01
 PARAMETER:  9.3082E-02  2.6504E-01  9.8539E-02  1.7659E-02  1.1411E-01 -4.7203E-03  3.3666E-01  2.4915E-01  3.3673E-02 -1.6740E-01
             4.6764E-02
 GRADIENT:   1.3062E+01  2.0617E+01  1.2747E+01 -9.1428E+00 -1.4637E+00 -5.6326E+00  1.4867E+01 -7.5162E+00 -2.0109E+00 -3.2409E-01
             5.0359E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1708.78545030296        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.9430E-01  1.1143E+00  9.2370E-01  9.5955E-01  9.5067E-01  9.1693E-01  1.2548E+00  1.3948E+00  9.0057E-01  6.0870E-01
             9.2565E-01
 PARAMETER:  9.4282E-02  2.0826E-01  2.0630E-02  5.8710E-02  4.9409E-02  1.3277E-02  3.2695E-01  4.3272E-01 -4.7228E-03 -3.9642E-01
             2.2742E-02
 GRADIENT:   1.8865E+01  4.1158E+00 -9.7606E+00  1.3873E+01  1.5106E+01  1.6230E+00  5.0547E+00  4.8330E+00  9.3481E-01 -4.8972E+00
            -5.2292E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1709.43401663831        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.8698E-01  1.0794E+00  1.0901E+00  9.8511E-01  1.0079E+00  9.1458E-01  1.2242E+00  1.4105E+00  9.1230E-01  7.5215E-01
             9.3332E-01
 PARAMETER:  8.6894E-02  1.7643E-01  1.8627E-01  8.4998E-02  1.0790E-01  1.0707E-02  3.0229E-01  4.4396E-01  8.2167E-03 -1.8481E-01
             3.0990E-02
 GRADIENT:  -4.7372E+01 -9.4162E+00 -2.1112E+00 -9.0331E+00  3.1734E+00 -4.7914E+00 -4.6654E-03 -2.0497E-01 -4.6136E-01  2.1326E-01
            -9.1844E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1710.04168901761        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      435
 NPARAMETR:  1.0060E+00  1.1215E+00  1.1427E+00  9.6975E-01  1.0398E+00  9.2443E-01  1.1774E+00  1.5288E+00  9.3517E-01  7.7565E-01
             9.3495E-01
 PARAMETER:  1.0599E-01  2.1469E-01  2.3341E-01  6.9285E-02  1.3903E-01  2.1417E-02  2.6332E-01  5.2448E-01  3.2973E-02 -1.5406E-01
             3.2737E-02
 GRADIENT:   2.7117E+00 -1.4276E+00 -6.4553E-01  3.9366E-01  1.2922E+00  2.3963E-01  2.3068E-01 -1.8913E-01 -1.6216E-01  1.2111E-01
             2.5213E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1710.18620071202        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      617
 NPARAMETR:  1.0060E+00  1.2880E+00  1.0596E+00  8.6585E-01  1.0764E+00  9.2443E-01  1.0476E+00  1.6529E+00  1.0195E+00  7.8566E-01
             9.3256E-01
 PARAMETER:  1.0594E-01  3.5308E-01  1.5792E-01 -4.4046E-02  1.7362E-01  2.1420E-02  1.4647E-01  6.0255E-01  1.1931E-01 -1.4123E-01
             3.0175E-02
 GRADIENT:   2.8794E-01  3.6771E+00  9.4168E-01  4.0272E+00 -2.3506E+00 -1.4461E-01 -6.6811E-01  3.3863E-01 -4.4316E-02 -4.4901E-01
            -9.4813E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1710.21510147859        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      805
 NPARAMETR:  1.0064E+00  1.3195E+00  1.0349E+00  8.4405E-01  1.0816E+00  9.2453E-01  1.0271E+00  1.6661E+00  1.0379E+00  7.8798E-01
             9.3247E-01
 PARAMETER:  1.0634E-01  3.7727E-01  1.3430E-01 -6.9546E-02  1.7841E-01  2.1535E-02  1.2674E-01  6.1047E-01  1.3719E-01 -1.3828E-01
             3.0081E-02
 GRADIENT:   9.2236E-01  2.8816E+00  1.2517E+00  2.6968E+00 -3.1469E+00 -1.5586E-01 -7.9995E-01  3.7741E-01  5.6323E-03 -3.3424E-01
            -9.4474E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1710.22864068534        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      971            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0060E+00  1.3165E+00  1.0298E+00  8.4221E-01  1.0836E+00  9.2490E-01  1.0342E+00  1.6553E+00  1.0379E+00  7.9179E-01
             9.3473E-01
 PARAMETER:  1.0599E-01  3.7501E-01  1.2933E-01 -7.1728E-02  1.8031E-01  2.1933E-02  1.3366E-01  6.0397E-01  1.3721E-01 -1.3345E-01
             3.2500E-02
 GRADIENT:   5.1908E+01  2.1909E+01 -2.5058E-01  5.2058E+00  3.0427E+00  5.2510E+00  9.2009E-01  8.1346E-01  1.1210E+00 -1.7779E-02
             2.6561E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1710.25808284269        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1149
 NPARAMETR:  1.0063E+00  1.3317E+00  9.7645E-01  8.3010E-01  1.0723E+00  9.2515E-01  1.0324E+00  1.5415E+00  1.0447E+00  7.9037E-01
             9.3370E-01
 PARAMETER:  1.0626E-01  3.8649E-01  7.6163E-02 -8.6210E-02  1.6978E-01  2.2204E-02  1.3188E-01  5.3277E-01  1.4371E-01 -1.3525E-01
             3.1403E-02
 GRADIENT:   1.4082E-01 -1.6738E+00  9.3841E-01 -1.3623E+00  3.0072E-01  3.4985E-02  1.9997E-01 -7.2378E-01  6.9888E-02  4.5290E-01
            -1.6371E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1710.34431536982        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1325
 NPARAMETR:  1.0080E+00  1.4646E+00  6.4653E-01  7.2912E-01  9.8952E-01  9.2648E-01  9.9401E-01  1.1471E+00  1.0885E+00  6.8827E-01
             9.3413E-01
 PARAMETER:  1.0794E-01  4.8156E-01 -3.3614E-01 -2.1592E-01  8.9463E-02  2.3633E-02  9.3993E-02  2.3720E-01  1.8482E-01 -2.7357E-01
             3.1865E-02
 GRADIENT:   2.6715E-01  1.7541E+00  7.6975E-02  2.0691E+00 -4.2500E+00  2.8801E-02  5.3324E-01  5.1518E-01 -3.9849E-01 -3.0855E-02
             1.2350E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1710.66917845839        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1505
 NPARAMETR:  1.0083E+00  1.6672E+00  4.0191E-01  5.8301E-01  9.5892E-01  9.2659E-01  8.9606E-01  7.5620E-01  1.2177E+00  6.2871E-01
             9.3414E-01
 PARAMETER:  1.0822E-01  6.1116E-01 -8.1153E-01 -4.3954E-01  5.8050E-02  2.3759E-02 -9.7479E-03 -1.7945E-01  2.9693E-01 -3.6409E-01
             3.1866E-02
 GRADIENT:   4.0744E-01  1.0297E+01  2.4876E+00  2.0682E+00 -1.3124E+01  8.2293E-02 -2.4398E+00  5.0262E-01  2.4466E-01 -6.4675E-01
            -1.6768E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1710.75997888323        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1704             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0081E+00  1.6865E+00  3.9306E-01  5.6983E-01  9.6722E-01  9.2654E-01  8.8769E-01  7.2841E-01  1.2380E+00  6.3692E-01
             9.3428E-01
 PARAMETER:  1.0809E-01  6.2265E-01 -8.3380E-01 -4.6242E-01  6.6666E-02  2.3706E-02 -1.9137E-02 -2.1690E-01  3.1349E-01 -3.5111E-01
             3.2019E-02
 GRADIENT:   5.2532E+01  7.5869E+01  3.0922E+00  1.4929E+01 -1.0325E+01  5.2256E+00 -1.6468E+00  4.3482E-01  1.3424E+00 -5.0498E-01
             3.7058E-03

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1710.76143820142        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1860
 NPARAMETR:  1.0081E+00  1.6865E+00  3.9306E-01  5.6983E-01  9.6721E-01  9.2654E-01  8.8822E-01  7.2841E-01  1.2380E+00  6.3692E-01
             9.3428E-01
 PARAMETER:  1.0809E-01  6.2264E-01 -8.3379E-01 -4.6243E-01  6.6665E-02  2.3697E-02 -1.8537E-02 -2.1689E-01  3.1349E-01 -3.5111E-01
             3.2021E-02
 GRADIENT:  -1.1608E+06  1.0076E+05 -1.5043E+05  1.3566E+05  1.2547E+06 -1.2547E+06  1.2547E+06 -5.7864E+05  2.0007E+05  2.4090E-01
            -6.2737E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1860
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7976E-04 -2.0436E-02 -1.7736E-02  1.8454E-02 -2.9406E-02
 SE:             2.9851E-02  2.6631E-02  7.6258E-03  2.2915E-02  1.8338E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9520E-01  4.4285E-01  2.0031E-02  4.2063E-01  1.0881E-01

 ETASHRINKSD(%)  1.0000E-10  1.0784E+01  7.4453E+01  2.3231E+01  3.8565E+01
 ETASHRINKVR(%)  1.0000E-10  2.0405E+01  9.3473E+01  4.1066E+01  6.2257E+01
 EBVSHRINKSD(%)  4.3082E-01  1.1949E+01  7.5727E+01  2.2992E+01  3.8557E+01
 EBVSHRINKVR(%)  8.5978E-01  2.2470E+01  9.4108E+01  4.0698E+01  6.2247E+01
 RELATIVEINF(%)  9.8946E+01  5.9444E+00  3.1570E-01  3.1895E+00  5.6336E+00
 EPSSHRINKSD(%)  4.5139E+01
 EPSSHRINKVR(%)  6.9902E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1710.7614382014240     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -975.61061163768579     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.47
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1710.761       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.69E+00  3.93E-01  5.70E-01  9.67E-01  9.27E-01  8.88E-01  7.28E-01  1.24E+00  6.37E-01  9.34E-01
 


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
+        2.64E+09
 
 TH 2
+        9.11E+02  2.84E+07
 
 TH 3
+       -2.94E+03 -9.11E+07  2.92E+08
 
 TH 4
+        3.64E+03 -2.22E+04  2.60E+05  4.52E+08
 
 TH 5
+        9.95E+03  5.03E+02  9.90E+08  3.82E+03  3.35E+09
 
 TH 6
+       -1.04E+04  1.73E+03 -5.55E+03  6.91E+03  1.88E+04  3.65E+09
 
 TH 7
+       -3.24E+09 -4.82E+02  1.54E+03 -2.00E+03 -5.37E+03  5.63E+03  3.98E+09
 
 TH 8
+       -6.10E+03 -9.07E+04 -4.34E+05 -3.61E+05 -5.30E+03 -1.15E+04  3.30E+03  1.26E+09
 
 TH 9
+        2.49E+03 -7.70E+07  1.77E+05 -6.08E+04  2.13E+03  4.69E+03 -1.33E+03 -2.45E+05  2.08E+08
 
 TH10
+        1.29E+09  3.08E+02 -2.37E+01  1.26E+03  3.64E+03 -1.51E+09  1.58E+09  8.88E+08  8.83E+02  6.27E+08
 
 TH11
+       -1.03E+04 -6.19E+03 -1.02E+09 -2.46E+04 -8.99E+03 -1.95E+04  5.59E+03  4.11E+04 -1.67E+04 -3.60E+03  3.59E+09
 
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
 #CPUT: Total CPU Time in Seconds,       29.148
Stop Time:
Sat Sep 25 10:15:06 CDT 2021
