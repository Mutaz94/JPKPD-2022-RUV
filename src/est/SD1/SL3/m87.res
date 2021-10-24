Sat Oct 23 15:37:29 CDT 2021
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
$DATA ../../../../data/SD1/SL3/dat87.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      986
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

 TOT. NO. OF OBS RECS:      886
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1203.44108992165        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1645E+02  1.7232E+01  3.0701E+02  1.3153E+02  2.8546E+02  3.9598E+01 -1.4370E+02 -7.5877E+02 -2.4135E+02 -5.0604E+01
            -8.8842E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2349.93401081534        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1124E+00  1.4009E+00  9.4502E-01  9.5309E-01  1.0233E+00  9.4863E-01  1.0640E+00  1.0119E+00  9.8215E-01  8.4414E-01
             5.2792E+00
 PARAMETER:  2.0655E-01  4.3710E-01  4.3447E-02  5.1950E-02  1.2299E-01  4.7268E-02  1.6200E-01  1.1185E-01  8.1989E-02 -6.9439E-02
             1.7638E+00
 GRADIENT:   9.0494E+01  8.9822E+01 -6.1920E-01  4.8744E+01 -4.8926E+01 -5.1563E+00  1.8433E+01  6.4390E+00  1.0761E+01  9.2918E+00
             7.5904E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2430.83409584840        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0270E+00  6.8564E-01  8.5414E-01  1.3825E+00  6.5299E-01  9.7988E-01  2.1556E+00  5.7399E-01  8.0418E-01  1.1296E-01
             4.4830E+00
 PARAMETER:  1.2666E-01 -2.7741E-01 -5.7662E-02  4.2391E-01 -3.2619E-01  7.9673E-02  8.6806E-01 -4.5514E-01 -1.1793E-01 -2.0807E+00
             1.6003E+00
 GRADIENT:  -5.5618E+01  3.5585E+01  1.5562E+01  2.8061E+02 -2.0520E+01 -1.4577E+00  4.8085E+01  4.1292E+00 -9.0867E+00  1.2501E-01
             5.6022E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2582.25046078995        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0129E+00  6.9517E-01  5.3021E-01  1.1235E+00  5.5347E-01  9.9638E-01  1.5381E+00  1.5188E-01  9.0038E-01  2.2762E-01
             3.0713E+00
 PARAMETER:  1.1282E-01 -2.6360E-01 -5.3449E-01  2.1644E-01 -4.9156E-01  9.6369E-02  5.3056E-01 -1.7847E+00 -4.9358E-03 -1.3801E+00
             1.2221E+00
 GRADIENT:   5.9732E+00  1.4061E+01  3.9218E+01 -2.4928E+00 -3.1235E+01  5.9030E+00 -1.1594E+01  2.4564E-01 -5.4154E+00 -4.5088E-01
            -2.0247E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2587.36346387844        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  1.0400E+00  8.9848E-01  6.7318E-01  1.0746E+00  7.1609E-01  9.9600E-01  1.4474E+00  7.6428E-02  9.8025E-01  4.1463E-01
             3.1246E+00
 PARAMETER:  1.3925E-01 -7.0469E-03 -2.9574E-01  1.7198E-01 -2.3395E-01  9.5993E-02  4.6975E-01 -2.4714E+00  8.0048E-02 -7.8038E-01
             1.2393E+00
 GRADIENT:   5.8400E+00  6.1274E+00 -3.3451E+00  6.5117E+00 -1.4569E+00  9.3275E-01  1.3192E+00  1.1088E-01  1.4200E+00 -1.6851E-01
             2.9697E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2591.28089821521        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  1.0379E+00  1.3016E+00  1.0456E+00  9.3986E-01  1.0886E+00  9.8854E-01  1.0043E+00  1.4528E-02  1.0268E+00  1.0205E+00
             3.1094E+00
 PARAMETER:  1.3721E-01  3.6357E-01  1.4462E-01  3.7973E-02  1.8493E-01  8.8477E-02  1.0430E-01 -4.1317E+00  1.2648E-01  1.2029E-01
             1.2344E+00
 GRADIENT:  -1.9601E+00  3.2457E+01 -5.9547E+00  7.8799E+01  7.5810E+00 -2.2144E+00  6.9397E+00  1.5460E-03 -6.5683E+00  4.0228E+00
            -3.1610E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2604.03712208230        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  1.0274E+00  1.9268E+00  1.6265E+00  5.4111E-01  1.6564E+00  9.9347E-01  6.1226E-01  1.0000E-02  1.6594E+00  1.5014E+00
             3.0282E+00
 PARAMETER:  1.2699E-01  7.5585E-01  5.8643E-01 -5.1413E-01  6.0464E-01  9.3450E-02 -3.9059E-01 -5.1395E+00  6.0644E-01  5.0639E-01
             1.2080E+00
 GRADIENT:  -2.1227E+01  6.3515E+01 -2.4189E+00  4.5228E+01  8.0488E+00 -1.3791E+00 -1.6561E+01  0.0000E+00  7.3365E-01  3.6892E+00
            -2.3204E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2610.18839664669        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  1.0366E+00  2.2349E+00  1.8428E+00  2.9769E-01  1.8381E+00  9.9632E-01  6.9071E-01  1.0000E-02  2.0325E+00  1.6157E+00
             3.0342E+00
 PARAMETER:  1.3598E-01  9.0418E-01  7.1129E-01 -1.1117E+00  7.0871E-01  9.6313E-02 -2.7003E-01 -5.2198E+00  8.0926E-01  5.7978E-01
             1.2099E+00
 GRADIENT:   2.7271E-01  5.8224E+00  2.7403E+00  2.0323E+00 -2.0996E+00  5.2793E-01  3.4197E+00  0.0000E+00 -2.9462E+00 -3.0218E+00
            -3.0533E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2611.26508200473        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  1.0399E+00  2.4978E+00  1.7636E+00  1.3536E-01  1.9065E+00  9.9424E-01  6.6642E-01  1.0000E-02  3.1678E+00  1.7083E+00
             3.0373E+00
 PARAMETER:  1.3916E-01  1.0154E+00  6.6737E-01 -1.8998E+00  7.4525E-01  9.4226E-02 -3.0584E-01 -5.3086E+00  1.2530E+00  6.3551E-01
             1.2110E+00
 GRADIENT:   5.7956E+00  4.2612E+01  4.1554E+00  4.6113E+00 -7.9388E+00 -6.4049E-01 -6.5760E-01  0.0000E+00 -2.8466E-01  2.3562E+00
             2.2204E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2613.26645626014        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1288
 NPARAMETR:  1.0346E+00  2.5954E+00  1.4154E+00  4.7832E-02  1.9228E+00  9.9599E-01  6.7513E-01  1.0000E-02  4.8119E+00  1.6655E+00
             3.0261E+00
 PARAMETER:  1.3403E-01  1.0537E+00  4.4740E-01 -2.9401E+00  7.5376E-01  9.5978E-02 -2.9284E-01 -5.6629E+00  1.6711E+00  6.1014E-01
             1.2073E+00
 GRADIENT:  -3.8272E+00  9.9319E+00  4.5204E+00 -2.3874E+00 -4.4514E+00  5.3619E-02  1.0971E+01  0.0000E+00 -8.1420E+00 -4.5054E+00
            -8.3841E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2615.72244308236        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1476
 NPARAMETR:  1.0405E+00  2.6108E+00  5.5511E-01  3.7378E-02  1.9434E+00  9.9779E-01  6.6892E-01  1.0000E-02  5.3071E+00  1.7013E+00
             3.0150E+00
 PARAMETER:  1.3975E-01  1.0597E+00 -4.8858E-01 -3.1867E+00  7.6444E-01  9.7784E-02 -3.0209E-01 -5.7704E+00  1.7690E+00  6.3138E-01
             1.2036E+00
 GRADIENT:   9.0802E+00  3.9260E+00  2.7805E-02 -2.0863E+00  2.3466E+00  3.5241E-01  1.0006E+01  0.0000E+00 -7.6506E+00  7.9087E-01
            -5.6966E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2616.03168140587        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1652
 NPARAMETR:  1.0342E+00  2.5940E+00  5.8258E-01  3.8618E-02  1.9289E+00  9.9453E-01  6.6788E-01  1.0000E-02  5.5009E+00  1.6837E+00
             3.0275E+00
 PARAMETER:  1.3362E-01  1.0532E+00 -4.4029E-01 -3.1540E+00  7.5696E-01  9.4517E-02 -3.0364E-01 -5.7704E+00  1.8049E+00  6.2100E-01
             1.2077E+00
 GRADIENT:  -3.8238E+00 -3.0500E+01  4.3141E-01  1.5530E+00 -1.9151E+00 -5.5756E-01  3.0928E+00  0.0000E+00  2.2093E+00 -9.5332E-01
             3.4242E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2616.16430820136        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1829
 NPARAMETR:  1.0349E+00  2.6058E+00  5.3143E-01  3.7868E-02  1.9355E+00  9.9593E-01  6.5371E-01  1.0000E-02  5.5070E+00  1.6915E+00
             3.0237E+00
 PARAMETER:  1.3429E-01  1.0578E+00 -5.3219E-01 -3.1737E+00  7.6034E-01  9.5918E-02 -3.2509E-01 -5.7704E+00  1.8060E+00  6.2559E-01
             1.2065E+00
 GRADIENT:  -2.8457E+00 -9.9654E+00 -5.9473E-02  6.8947E-01 -1.5815E-01 -2.3787E-01 -6.9063E-01  0.0000E+00 -8.5517E-01 -9.7993E-02
            -2.4855E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2617.29592468907        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2009
 NPARAMETR:  1.0405E+00  2.6446E+00  5.4762E-01  1.2682E-02  1.9349E+00  9.9961E-01  6.4624E-01  1.0000E-02  9.0162E+00  1.6964E+00
             3.0250E+00
 PARAMETER:  1.3968E-01  1.0725E+00 -5.0218E-01 -4.2676E+00  7.6004E-01  9.9612E-02 -3.3658E-01 -5.7704E+00  2.2990E+00  6.2853E-01
             1.2069E+00
 GRADIENT:   8.5906E+00  3.0772E+00 -2.2382E-02 -1.0653E+00 -1.2222E+00  1.0237E+00  1.8830E+00  0.0000E+00 -2.3216E+00  5.1931E-01
             1.1502E+00

0ITERATION NO.:   69    OBJECTIVE VALUE:  -2617.38134668081        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     2143
 NPARAMETR:  1.0338E+00  2.6431E+00  5.4847E-01  1.4758E-02  1.9413E+00  9.9525E-01  6.4516E-01  1.0000E-02  8.8010E+00  1.6904E+00
             3.0228E+00
 PARAMETER:  1.3192E-01  1.0712E+00 -5.0497E-01 -4.1124E+00  7.6386E-01  9.4238E-02 -3.3855E-01 -5.7704E+00  2.2768E+00  6.2447E-01
             1.2062E+00
 GRADIENT:  -7.8273E+00 -1.2671E+02 -7.0387E-02  3.4341E+01  1.7257E+00 -1.0147E+00 -2.2007E+02  0.0000E+00  6.1203E+01 -6.5369E-01
            -1.8783E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2143
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.9140E-03 -1.3698E-02  2.1983E-06  7.1528E-03 -1.7309E-02
 SE:             2.9245E-02  2.7460E-02  2.1930E-06  6.1807E-03  2.6602E-02
 N:                     100         100         100         100         100

 P VAL.:         8.6656E-01  6.1791E-01  3.1614E-01  2.4716E-01  5.1528E-01

 ETASHRINKSD(%)  2.0248E+00  8.0056E+00  9.9993E+01  7.9294E+01  1.0878E+01
 ETASHRINKVR(%)  4.0085E+00  1.5370E+01  1.0000E+02  9.5713E+01  2.0573E+01
 EBVSHRINKSD(%)  2.0883E+00  7.6520E+00  9.9950E+01  8.5478E+01  9.7699E+00
 EBVSHRINKVR(%)  4.1331E+00  1.4718E+01  1.0000E+02  9.7891E+01  1.8585E+01
 RELATIVEINF(%)  9.5756E+01  4.0218E+01  2.4958E-05  1.0058E+00  8.0510E+01
 EPSSHRINKSD(%)  1.5135E+01
 EPSSHRINKVR(%)  2.7980E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          886
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1628.3590808386800     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2617.3813466808074     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -989.02226584212735     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2617.381       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.64E+00  5.46E-01  1.48E-02  1.94E+00  9.94E-01  6.45E-01  1.00E-02  8.82E+00  1.69E+00  3.02E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,      152.881
Stop Time:
Sat Oct 23 15:37:54 CDT 2021
