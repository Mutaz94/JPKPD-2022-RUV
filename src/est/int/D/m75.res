Sat Sep 18 07:28:43 CDT 2021
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
$DATA ../../../../data/int/D/dat75.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   29777.9100204960        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1774E+02  3.3511E+02 -7.3632E+01 -5.0600E+00  2.6502E+02 -2.5543E+03 -1.1954E+03 -1.0253E+02 -2.3661E+03 -7.0578E+02
            -6.0301E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -968.451705559780        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.9178E+00  1.6824E+00  9.7334E-01  2.5344E+00  9.2790E-01  5.1811E+00  4.1759E+00  9.9523E-01  3.9355E+00  1.6993E+00
             1.1456E+01
 PARAMETER:  7.5115E-01  6.2024E-01  7.2978E-02  1.0299E+00  2.5168E-02  1.7450E+00  1.5293E+00  9.5221E-02  1.4700E+00  6.3023E-01
             2.5385E+00
 GRADIENT:   2.7111E+01  1.8918E+00 -3.6916E+01  7.7251E+01 -4.1137E+01  1.3143E+02  4.9766E+00  4.4722E+00  5.7028E+01  3.9729E+01
             3.6201E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1045.30197968304        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0643E+00  1.4651E+00  1.3309E+02  7.7441E+00  2.6781E+00  4.4411E+00  7.1408E+00  8.4311E-01  6.0752E+00  1.5264E+00
             1.0736E+01
 PARAMETER:  1.6229E-01  4.8195E-01  4.9910E+00  2.1469E+00  1.0851E+00  1.5909E+00  2.0658E+00 -7.0655E-02  1.9042E+00  5.2294E-01
             2.4736E+00
 GRADIENT:  -2.2461E+01  1.1758E+01 -1.5147E+00  8.5178E+01  1.1219E+01  1.2449E+02  1.3536E+01  2.8894E-01  2.3452E+01  3.6301E+01
             3.3828E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1174.00355992868        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.7144E+00  9.5507E-01  1.1888E+01  1.4690E+00  1.6455E+00  2.9026E+00  5.7421E+00  2.7583E+00  2.6171E+00  1.0684E+00
             1.0706E+01
 PARAMETER:  6.3909E-01  5.4033E-02  2.5755E+00  4.8457E-01  5.9807E-01  1.1656E+00  1.8478E+00  1.1146E+00  1.0621E+00  1.6616E-01
             2.4708E+00
 GRADIENT:   6.1690E+01 -9.7049E+00  5.3878E+00 -2.2556E+01 -8.7921E+01  9.7421E+00 -4.3059E+00  4.7012E+00  3.4164E+01  1.9509E+01
             3.5379E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1264.28652685856        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1658E+00  7.5209E-01  3.0483E+01  1.4632E+00  2.3576E+00  2.3662E+00  6.3658E+00  1.3004E+00  1.8216E+00  4.6483E-01
             8.5582E+00
 PARAMETER:  2.5340E-01 -1.8489E-01  3.5172E+00  4.8064E-01  9.5766E-01  9.6128E-01  1.9509E+00  3.6264E-01  6.9971E-01 -6.6608E-01
             2.2469E+00
 GRADIENT:  -5.5651E+00 -9.3040E+00 -6.6344E-01 -8.9854E+00  5.5170E+00 -9.8878E+00  1.4951E+00  4.1497E-02  1.5710E+01  4.0836E+00
             2.4595E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1269.78364984506        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.1775E+00  1.4530E+00  1.1131E+01  1.0420E+00  2.2878E+00  2.4368E+00  5.2642E+00  7.6171E-01  9.6256E-01  2.3779E-01
             8.4635E+00
 PARAMETER:  2.6339E-01  4.7362E-01  2.5097E+00  1.4110E-01  9.2758E-01  9.9069E-01  1.7609E+00 -1.7219E-01  6.1842E-02 -1.3364E+00
             2.2358E+00
 GRADIENT:  -8.5460E-01 -1.2718E+00 -2.1650E+00 -8.3872E-01  5.6488E-01 -3.0416E-01  2.9257E+00  8.1116E-02  2.1214E+00  9.4705E-01
            -1.1308E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1270.49070286124        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  1.1798E+00  1.6894E+00  1.6620E+01  9.0621E-01  2.4039E+00  2.4377E+00  4.9438E+00  7.5540E-01  5.4578E-01  1.4265E-01
             8.5367E+00
 PARAMETER:  2.6533E-01  6.2435E-01  2.9106E+00  1.5127E-03  9.7710E-01  9.9104E-01  1.6981E+00 -1.8051E-01 -5.0554E-01 -1.8474E+00
             2.2444E+00
 GRADIENT:  -8.0060E-01 -1.1236E+00 -6.6090E-01 -7.5532E-01  1.9238E+00 -8.0441E-02  1.0513E+00  2.6539E-02 -5.5014E-01  3.3240E-01
             3.1408E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1270.51729753404        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      514
 NPARAMETR:  1.1788E+00  1.7558E+00  2.2783E+01  8.8364E-01  2.4350E+00  2.4416E+00  4.8680E+00  7.6795E-01  4.6984E-01  1.2059E-01
             8.5299E+00
 PARAMETER:  2.6447E-01  6.6290E-01  3.2260E+00 -2.3703E-02  9.8995E-01  9.9263E-01  1.6827E+00 -1.6403E-01 -6.5536E-01 -2.0153E+00
             2.2436E+00
 GRADIENT:  -1.0079E+00  7.9102E-02 -2.7794E-01  1.5001E+00  3.6472E-01  3.2766E-01 -5.0265E-01  1.2989E-02 -7.2810E-01  2.3407E-01
             2.8858E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1271.42162219080        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      618
 NPARAMETR:  1.1808E+00  1.5234E+00  1.9632E+03  1.0024E+00  2.5060E+00  2.4364E+00  5.1543E+00  1.7304E+00  8.3675E-01  7.7563E-02
             8.5061E+00
 PARAMETER:  2.6623E-01  5.2094E-01  7.6823E+00  1.0238E-01  1.0187E+00  9.9053E-01  1.7398E+00  6.4835E-01 -7.8232E-02 -2.4567E+00
             2.2408E+00
 GRADIENT:  -2.0246E+00 -2.6481E+00 -1.1250E-03 -3.0550E+00 -3.9005E-01 -3.0409E+00 -9.6032E+00  1.1804E-05  6.1023E-01  9.9437E-02
            -6.2397E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1271.95074621448        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  1.1893E+00  1.3674E+00  5.4373E+03  1.1112E+00  2.5167E+00  2.4554E+00  5.7012E+00  2.1625E+00  9.2952E-01  7.5173E-02
             8.5546E+00
 PARAMETER:  2.7337E-01  4.1291E-01  8.7010E+00  2.0542E-01  1.0230E+00  9.9830E-01  1.8407E+00  8.7128E-01  2.6917E-02 -2.4880E+00
             2.2465E+00
 GRADIENT:  -1.3764E-01  3.7579E-01 -8.6139E-04  3.7267E-01  2.9006E-01  1.9103E-01 -1.6986E-01  4.7985E-06 -5.2346E-02  9.3282E-02
             1.6234E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1271.97706255225        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      968
 NPARAMETR:  1.1973E+00  1.3368E+00  7.7865E+04  1.1185E+00  2.5155E+00  2.4483E+00  5.7447E+00  3.0404E+00  9.4032E-01  4.3715E-02
             8.5505E+00
 PARAMETER:  2.8009E-01  3.9027E-01  1.1363E+01  2.1198E-01  1.0225E+00  9.9539E-01  1.8483E+00  1.2120E+00  3.8469E-02 -3.0301E+00
             2.2460E+00
 GRADIENT:   4.0400E-02 -1.6092E-01 -6.3671E-05 -7.1851E-02  1.4780E-01  2.5890E-01  1.1481E-01  3.5888E-07 -1.1578E-01  3.1513E-02
            -6.1880E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1271.98904083657        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1167             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1923E+00  1.3387E+00  7.6956E+04  1.1215E+00  2.5150E+00  2.4534E+00  5.7446E+00  3.0454E+00  9.4039E-01  3.7508E-02
             8.5541E+00
 PARAMETER:  2.7591E-01  3.9170E-01  1.1351E+01  2.1463E-01  1.0223E+00  9.9748E-01  1.8483E+00  1.2136E+00  3.8534E-02 -3.1832E+00
             2.2464E+00
 GRADIENT:   7.6426E-01  6.9271E-01 -4.5224E-05  1.5588E+00  2.4789E-01  3.8055E+00  1.0297E+01  1.5328E-05 -2.5679E-01  2.3379E-02
             4.4399E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1271.98994288893        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1355             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1903E+00  1.3427E+00  7.5556E+04  1.1209E+00  2.5154E+00  2.4537E+00  5.7463E+00  3.0379E+00  9.4013E-01  3.7087E-02
             8.5538E+00
 PARAMETER:  2.7423E-01  3.9467E-01  1.1333E+01  2.1417E-01  1.0224E+00  9.9759E-01  1.8486E+00  1.2112E+00  3.8263E-02 -3.1945E+00
             2.2464E+00
 GRADIENT:   5.4050E-01  8.2210E-01 -5.2514E-05  1.2386E+00  3.5223E-01  4.0421E+00  1.0608E+01  1.6045E-05 -1.0866E-01  2.2870E-02
             4.2628E+00

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1271.98998142200        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     1450
 NPARAMETR:  1.1907E+00  1.3423E+00  7.5558E+04  1.1209E+00  2.5155E+00  2.4527E+00  5.7464E+00  3.0379E+00  9.4014E-01  3.7068E-02
             8.5537E+00
 PARAMETER:  2.7451E-01  3.9440E-01  1.1333E+01  2.1409E-01  1.0225E+00  9.9721E-01  1.8486E+00  1.2112E+00  3.8275E-02 -3.1950E+00
             2.2464E+00
 GRADIENT:   2.4051E-01  9.2586E-02 -7.3931E-04  1.7017E-01 -5.8617E-02 -1.4885E-01  4.1426E-02 -3.2322E-03  1.2377E-01  6.0325E-03
            -2.6512E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1450
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1362E-02  2.4115E-02  2.7818E-08 -6.2503E-02 -3.3283E-05
 SE:             2.8916E-02  2.5525E-02  4.7679E-07  1.1956E-02  4.5630E-04
 N:                     100         100         100         100         100

 P VAL.:         6.9437E-01  3.4479E-01  9.5348E-01  1.7204E-07  9.4185E-01

 ETASHRINKSD(%)  3.1286E+00  1.4488E+01  9.9998E+01  5.9945E+01  9.8471E+01
 ETASHRINKVR(%)  6.1593E+00  2.6877E+01  1.0000E+02  8.3956E+01  9.9977E+01
 EBVSHRINKSD(%)  3.6515E+00  8.7056E+00  9.9998E+01  6.7139E+01  9.8150E+01
 EBVSHRINKVR(%)  7.1697E+00  1.6653E+01  1.0000E+02  8.9202E+01  9.9966E+01
 RELATIVEINF(%)  9.2741E+01  4.3936E+01  1.1213E-08  5.6719E+00  8.2696E-03
 EPSSHRINKSD(%)  6.6949E+00
 EPSSHRINKVR(%)  1.2942E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1271.9899814220030     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       382.09937834640778     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.78
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1271.990       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  1.34E+00  7.56E+04  1.12E+00  2.52E+00  2.45E+00  5.75E+00  3.04E+00  9.40E-01  3.71E-02  8.55E+00
 


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
+        2.54E+02
 
 TH 2
+        2.79E+01  6.83E+01
 
 TH 3
+        3.75E-05  4.43E-07  3.35E-11
 
 TH 4
+        2.90E+02  1.11E+02 -4.75E-05  4.18E+02
 
 TH 5
+       -1.31E+01 -7.44E+00 -5.03E-06 -2.23E+01  4.23E+01
 
 TH 6
+       -4.24E+00  5.90E+00  5.76E-06 -6.19E+00  6.91E-01  2.98E+01
 
 TH 7
+       -3.49E+00 -1.44E-01 -1.01E-06 -1.15E+01  5.49E-01  4.40E-01  4.03E+00
 
 TH 8
+       -1.19E+01  3.39E-01  5.78E-06  1.21E+01 -5.55E-01  1.60E+00 -3.31E-02  1.59E+00
 
 TH 9
+       -4.79E+01  2.42E+02  3.82E-05 -6.81E+02  5.09E+01  3.54E+01  6.78E+00 -5.50E+01  2.31E+03
 
 TH10
+        9.86E+01  2.05E+02 -9.60E-05  1.80E+02  5.41E+01  3.50E+01 -1.12E+01 -2.12E+01 -1.12E+02  1.19E+03
 
 TH11
+       -5.65E+00 -2.89E+00 -6.13E-07 -1.03E+01  6.54E-01  1.55E+00  5.48E-01  2.55E-01  9.00E+00 -9.00E-01  1.40E+01
 
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
 #CPUT: Total CPU Time in Seconds,       55.388
Stop Time:
Sat Sep 18 07:29:40 CDT 2021
