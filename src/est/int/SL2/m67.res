Sat Sep 18 03:57:55 CDT 2021
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
$DATA ../../../../data/int/SL2/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1659.57961376552        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.2450E+01 -5.9197E+01 -1.0274E+01 -5.6757E+01  8.7791E+01 -1.2405E+01 -1.1722E+02 -1.6066E+02 -1.3416E+02 -9.5381E+00
            -4.0809E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2953.98295237543        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9927E-01  1.3709E+00  1.3088E+00  8.8196E-01  1.2812E+00  1.0786E+00  1.1553E+00  9.4106E-01  1.1507E+00  9.5617E-01
             2.0676E+00
 PARAMETER:  9.9268E-02  4.1544E-01  3.6911E-01 -2.5612E-02  3.4777E-01  1.7570E-01  2.4439E-01  3.9247E-02  2.4040E-01  5.5181E-02
             8.2638E-01
 GRADIENT:   8.3853E+00  1.0785E+01 -4.5334E+00  5.5909E+00  2.9441E+01  1.7814E+01  1.6257E+01 -6.2361E+00 -1.2971E+01 -3.3605E+01
            -2.2204E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2962.40479436731        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9910E-01  1.4844E+00  1.5733E+00  8.6103E-01  1.4350E+00  1.0709E+00  9.0926E-01  5.9898E-01  1.1356E+00  1.2594E+00
             2.1432E+00
 PARAMETER:  9.9104E-02  4.9504E-01  5.5320E-01 -4.9625E-02  4.6114E-01  1.6850E-01  4.8802E-03 -4.1253E-01  2.2718E-01  3.3065E-01
             8.6232E-01
 GRADIENT:   5.5431E+00  3.5929E+01  8.4721E+00  5.5418E+01  1.8348E+01  1.5877E+01 -8.9703E+00 -3.7147E+00 -1.2231E+01 -9.3396E+00
            -1.4781E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2973.85855631201        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0034E+00  1.6862E+00  1.2671E+00  7.0751E-01  1.4839E+00  1.0411E+00  9.0507E-01  1.6528E+00  1.3718E+00  1.2867E+00
             2.2022E+00
 PARAMETER:  1.0338E-01  6.2248E-01  3.3675E-01 -2.4600E-01  4.9470E-01  1.4026E-01  2.5413E-04  6.0245E-01  4.1611E-01  3.5209E-01
             8.8944E-01
 GRADIENT:   1.1963E+01  2.9950E+01 -8.8126E+00  1.7476E+01 -5.1253E+00  4.5066E+00  9.9869E+00  1.6060E+00  7.3692E+00 -8.0113E+00
            -3.9238E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2978.22297815258        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.9333E-01  1.7043E+00  1.9356E+00  6.7917E-01  1.6736E+00  1.0177E+00  8.8647E-01  2.8592E+00  1.2952E+00  1.3888E+00
             2.2270E+00
 PARAMETER:  9.3304E-02  6.3316E-01  7.6041E-01 -2.8689E-01  6.1495E-01  1.1753E-01 -2.0508E-02  1.1505E+00  3.5866E-01  4.2841E-01
             9.0063E-01
 GRADIENT:  -1.0061E+01 -1.1293E+01 -8.8689E+00 -4.2149E+00  1.6690E+01 -4.3425E+00  3.0790E+00 -1.1327E+00  7.9166E+00 -1.4324E+00
             1.1095E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2979.78478291134        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0026E+00  1.7046E+00  2.2796E+00  6.7948E-01  1.6724E+00  1.0353E+00  8.5872E-01  2.8889E+00  1.2933E+00  1.3889E+00
             2.2275E+00
 PARAMETER:  1.0259E-01  6.3333E-01  9.2398E-01 -2.8642E-01  6.1428E-01  1.3472E-01 -5.2315E-02  1.1609E+00  3.5718E-01  4.2848E-01
             9.0089E-01
 GRADIENT:   1.0002E+01 -6.8145E-01  7.9702E-01 -1.4973E+01 -2.9165E+00  2.7351E+00 -1.1879E+00 -6.2919E+00  5.4343E+00 -3.3150E+00
             1.3533E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2979.81669110349        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      560
 NPARAMETR:  1.0009E+00  1.7046E+00  2.2828E+00  6.7948E-01  1.6807E+00  1.0311E+00  8.6786E-01  2.8889E+00  1.2933E+00  1.3889E+00
             2.2275E+00
 PARAMETER:  1.0093E-01  6.3333E-01  9.2542E-01 -2.8642E-01  6.1922E-01  1.3058E-01 -4.1722E-02  1.1609E+00  3.5718E-01  4.2848E-01
             9.0089E-01
 GRADIENT:  -2.3923E+00 -2.1099E+01  1.0879E-01 -1.6845E+01 -2.9785E+00 -4.5497E-01  5.2150E-02 -6.6766E+00  5.5898E+00 -4.1467E+00
             1.2161E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2980.74121853857        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      716
 NPARAMETR:  1.0014E+00  1.7145E+00  2.3046E+00  6.9004E-01  1.6866E+00  1.0317E+00  8.6834E-01  3.1419E+00  1.2460E+00  1.4068E+00
             2.2192E+00
 PARAMETER:  1.0140E-01  6.3913E-01  9.3491E-01 -2.7101E-01  6.2273E-01  1.3120E-01 -4.1178E-02  1.2448E+00  3.1994E-01  4.4132E-01
             8.9713E-01
 GRADIENT:   7.2866E+00  1.7428E+01 -3.4766E+00 -6.6273E-01  5.6365E-01  1.1663E+00 -1.2697E+00 -2.7726E+00  4.4797E+00 -1.5928E+00
             7.2763E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2981.10692995778        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:      849             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0035E+00  1.7111E+00  2.3486E+00  6.9967E-01  1.6888E+00  1.0343E+00  8.9512E-01  3.2152E+00  1.1448E+00  1.4321E+00
             2.2175E+00
 PARAMETER:  1.0345E-01  6.3711E-01  9.5382E-01 -2.5715E-01  6.2399E-01  1.3374E-01 -1.0801E-02  1.2679E+00  2.3526E-01  4.5916E-01
             8.9638E-01
 GRADIENT:   1.1476E+01  2.6682E+01 -3.9932E+00  2.8065E+00 -8.1283E-01  2.1439E+00 -1.3567E+00 -3.3981E+00 -8.1770E-01  1.7136E+00
             5.9535E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2981.12198636216        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  1.0023E+00  1.7111E+00  2.3486E+00  6.9926E-01  1.6888E+00  1.0329E+00  9.0088E-01  3.2152E+00  1.1559E+00  1.4252E+00
             2.2175E+00
 PARAMETER:  1.0229E-01  6.3711E-01  9.5382E-01 -2.5774E-01  6.2399E-01  1.3236E-01 -4.3885E-03  1.2679E+00  2.4485E-01  4.5429E-01
             8.9638E-01
 GRADIENT:   6.4176E-02  7.1667E+00 -4.3036E+00  8.2320E-02 -4.6622E+00  1.9141E-02  8.2471E-05 -3.5559E+00  3.2713E-02  1.1549E-01
             4.5817E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2981.17706519706        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1172             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0024E+00  1.7088E+00  2.3629E+00  7.0163E-01  1.6899E+00  1.0329E+00  9.0117E-01  3.2319E+00  1.1530E+00  1.4247E+00
             2.2179E+00
 PARAMETER:  1.0238E-01  6.3577E-01  9.5988E-01 -2.5435E-01  6.2468E-01  1.3240E-01 -4.0605E-03  1.2731E+00  2.4238E-01  4.5399E-01
             8.9656E-01
 GRADIENT:   9.2002E+00  2.6901E+01 -4.0522E+00  3.1890E+00  2.9990E-01  1.5962E+00  7.7031E-02 -3.0386E+00  3.4900E-01  7.4310E-01
             6.8798E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2981.29293509625        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1328            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0023E+00  1.7053E+00  2.4261E+00  7.0081E-01  1.6941E+00  1.0328E+00  9.0221E-01  3.2309E+00  1.1518E+00  1.4243E+00
             2.2153E+00
 PARAMETER:  1.0227E-01  6.3373E-01  9.8627E-01 -2.5552E-01  6.2714E-01  1.3232E-01 -2.9135E-03  1.2727E+00  2.4136E-01  4.5366E-01
             8.9539E-01
 GRADIENT:   9.1277E+00  2.3508E+01 -2.1994E+00 -9.9937E-01  4.7921E-01  1.5800E+00  1.7620E-01 -4.1389E+00  2.8008E-01  4.4007E-01
             4.6717E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2981.32505881105        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:     1454
 NPARAMETR:  9.9983E-01  1.6930E+00  2.4625E+00  7.0164E-01  1.6939E+00  1.0319E+00  9.0457E-01  3.2298E+00  1.1499E+00  1.4226E+00
             2.2133E+00
 PARAMETER:  9.9831E-02  6.2649E-01  1.0012E+00 -2.5434E-01  6.2704E-01  1.3140E-01 -2.9461E-04  1.2724E+00  2.3968E-01  4.5246E-01
             8.9447E-01
 GRADIENT:   4.3186E+00  1.3238E+01  6.0158E+02 -7.6482E+00  1.1596E+00  1.2539E+00  3.0476E-01 -4.4913E+00  3.4036E-01  4.3245E-01
             3.5247E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2981.33966533589        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1635            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9989E-01  1.6931E+00  2.4625E+00  7.0192E-01  1.6941E+00  1.0331E+00  9.0456E-01  3.2362E+00  1.1491E+00  1.4226E+00
             2.2133E+00
 PARAMETER:  9.9893E-02  6.2658E-01  1.0012E+00 -2.5394E-01  6.2715E-01  1.3252E-01 -3.0231E-04  1.2744E+00  2.3894E-01  4.5248E-01
             8.9448E-01
 GRADIENT:   4.2713E+00  1.3629E+01 -1.6337E+00 -7.2852E+00  1.2256E+00  1.6935E+00  2.6181E-01 -4.3689E+00  3.2894E-01  4.2775E-01
             3.6022E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2981.35551507414        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1776
 NPARAMETR:  1.0027E+00  1.6975E+00  2.4612E+00  7.0201E-01  1.7019E+00  1.0316E+00  9.0446E-01  3.2341E+00  1.1471E+00  1.4223E+00
             2.2143E+00
 PARAMETER:  1.0269E-01  6.2918E-01  1.0007E+00 -2.5381E-01  6.3176E-01  1.3110E-01 -4.1306E-04  1.2737E+00  2.3722E-01  4.5225E-01
             8.9494E-01
 GRADIENT:   1.0493E+00 -1.8309E+00 -2.2650E+00 -6.7733E+00  1.8262E-01 -4.3866E-01 -8.8165E-03 -4.9576E+00 -4.1379E-02 -6.4057E-01
             2.4681E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2981.35736845976        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1954
 NPARAMETR:  1.0022E+00  1.6993E+00  2.4613E+00  7.0201E-01  1.7022E+00  1.0332E+00  9.0459E-01  3.2342E+00  1.1486E+00  1.4223E+00
             2.2143E+00
 PARAMETER:  1.0219E-01  6.3019E-01  1.0007E+00 -2.5381E-01  6.3195E-01  1.3266E-01 -2.6991E-04  1.2738E+00  2.3856E-01  4.5225E-01
             8.9492E-01
 GRADIENT:   7.1423E-03 -1.8434E-01 -2.2276E+00 -5.8168E+00  1.0258E-01  1.3831E-01  1.0847E-01 -4.9724E+00  7.1539E-02 -6.8040E-01
             2.4100E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2981.36303938057        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2142
 NPARAMETR:  1.0021E+00  1.6991E+00  2.4616E+00  7.0213E-01  1.7029E+00  1.0328E+00  9.0343E-01  3.2373E+00  1.1452E+00  1.4221E+00
             2.2141E+00
 PARAMETER:  1.0214E-01  6.3010E-01  1.0008E+00 -2.5363E-01  6.3235E-01  1.3227E-01 -1.5554E-03  1.2747E+00  2.3557E-01  4.5217E-01
             8.9484E-01
 GRADIENT:  -1.0062E-01 -2.8688E-01 -2.3029E+00 -5.8067E+00  4.1355E-01 -7.1358E-03 -2.6997E-01 -4.9730E+00 -2.3215E-01 -7.5970E-01
             2.1591E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2981.37445007943        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2318
 NPARAMETR:  1.0022E+00  1.6994E+00  2.4612E+00  7.0231E-01  1.7026E+00  1.0328E+00  9.0426E-01  3.2432E+00  1.1473E+00  1.4222E+00
             2.2143E+00
 PARAMETER:  1.0221E-01  6.3033E-01  1.0008E+00 -2.5342E-01  6.3234E-01  1.3232E-01 -6.3763E-04  1.2767E+00  2.3765E-01  4.5226E-01
             8.9481E-01
 GRADIENT:   3.5481E-02  2.0864E-01  1.4613E+03 -2.8991E+03  2.4673E-01  3.0458E-03 -2.4266E-03  1.1457E+03  3.9587E-02  3.2430E+03
            -1.6432E+03
 NUMSIGDIG:         3.8         3.7         3.3         3.3         3.0         4.2         4.0         3.3         2.4         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2318
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3643E-03 -2.5212E-02 -2.4478E-02  2.4769E-02 -2.9575E-02
 SE:             2.9560E-02  2.3003E-02  1.6291E-02  2.0871E-02  2.4394E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6319E-01  2.7308E-01  1.3297E-01  2.3532E-01  2.2535E-01

 ETASHRINKSD(%)  9.7039E-01  2.2936E+01  4.5422E+01  3.0081E+01  1.8279E+01
 ETASHRINKVR(%)  1.9314E+00  4.0611E+01  7.0212E+01  5.1113E+01  3.3216E+01
 EBVSHRINKSD(%)  1.1603E+00  2.2869E+01  5.0530E+01  3.5091E+01  1.4385E+01
 EBVSHRINKVR(%)  2.3071E+00  4.0508E+01  7.5527E+01  5.7868E+01  2.6701E+01
 RELATIVEINF(%)  9.7668E+01  7.4588E+00  1.0979E+01  5.1294E+00  3.8152E+01
 EPSSHRINKSD(%)  1.7918E+01
 EPSSHRINKVR(%)  3.2625E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2981.3744500794346     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1327.2850903110239     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    69.44
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2981.374       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.70E+00  2.46E+00  7.02E-01  1.70E+00  1.03E+00  9.04E-01  3.24E+00  1.15E+00  1.42E+00  2.21E+00
 


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
+        1.02E+03
 
 TH 2
+       -6.42E+00  2.91E+02
 
 TH 3
+        6.14E+01 -1.49E+02  6.03E+04
 
 TH 4
+       -8.59E+02  2.53E+03  5.24E+02  1.16E+07
 
 TH 5
+       -1.95E+00 -4.34E+01 -3.03E+02  4.11E+03  3.17E+05
 
 TH 6
+        2.50E+00 -2.38E+00  7.58E+01 -1.06E+03 -1.38E+00  1.80E+02
 
 TH 7
+       -5.44E+00  7.85E+00 -1.41E+02  1.94E+03 -7.86E-01  1.70E+00  1.03E+02
 
 TH 8
+        3.63E+01 -9.64E+01  8.10E+01  3.27E+02 -1.76E+02  4.50E+01 -8.36E+01  2.14E+04
 
 TH 9
+        1.79E+00 -4.43E+00  7.76E+01 -1.06E+03  5.00E+00 -2.22E+00  2.53E+01  4.88E+01  3.41E+01
 
 TH10
+        5.57E+06 -5.33E+05 -5.36E+01 -3.20E+06 -5.30E+05  2.91E+02 -5.41E+02 -3.19E+01  3.03E+02  8.86E+05
 
 TH11
+       -8.91E+01  1.81E+02 -5.28E+02 -6.82E+02  3.63E+02 -9.21E+01  1.81E+02 -4.49E+04 -9.30E+01  7.17E+01  9.43E+04
 
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
 #CPUT: Total CPU Time in Seconds,       84.888
Stop Time:
Sat Sep 18 03:59:21 CDT 2021
