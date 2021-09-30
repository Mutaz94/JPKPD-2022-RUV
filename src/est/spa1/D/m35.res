Thu Sep 30 02:57:26 CDT 2021
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
$DATA ../../../../data/spa1/D/dat35.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16611.6141924709        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3035E+02  2.9299E+02 -1.1772E+02  1.2092E+02  1.5344E+02 -1.1497E+03 -6.4988E+02 -3.3808E+01 -1.3343E+03 -3.5924E+02
            -3.3665E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -727.691931424391        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3916E+00  1.1459E+00  1.1204E+00  2.0125E+00  1.2252E+00  2.0839E+00  1.5252E+00  9.2381E-01  2.5266E+00  1.0296E+00
             1.3435E+01
 PARAMETER:  4.3048E-01  2.3617E-01  2.1367E-01  7.9937E-01  3.0312E-01  8.3424E-01  5.2211E-01  2.0751E-02  1.0269E+00  1.2919E-01
             2.6979E+00
 GRADIENT:  -3.3532E+01  2.3451E+01 -2.3725E+01  6.5188E+01 -2.7696E+00  4.1132E+01  2.9040E+00  5.4807E+00  1.6019E+01  3.4771E+00
             2.4134E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -753.049458173734        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3788E+00  1.2611E+00  3.7597E+00  1.9585E+00  4.9542E+00  1.8665E+00  1.7005E+00  5.1162E-01  3.0432E+00  4.8056E+00
             1.2477E+01
 PARAMETER:  4.2123E-01  3.3201E-01  1.4243E+00  7.7218E-01  1.7002E+00  7.2409E-01  6.3091E-01 -5.7017E-01  1.2129E+00  1.6698E+00
             2.6239E+00
 GRADIENT:  -2.4350E+01  1.2811E+01  5.8495E+00  4.9728E+01 -2.4602E-01 -1.7035E+00  8.0445E+00 -1.0039E-01  2.5266E+01 -2.5562E-01
             2.3345E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -793.682896335084        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.2584E+00  1.4249E+00  2.5853E+00  1.0891E+00  3.9146E+00  1.9167E+00  5.3977E-01  4.1975E-01  3.0396E+00  5.7708E+00
             1.0249E+01
 PARAMETER:  3.2985E-01  4.5412E-01  1.0498E+00  1.8538E-01  1.4647E+00  7.5062E-01 -5.1662E-01 -7.6809E-01  1.2117E+00  1.8528E+00
             2.4271E+00
 GRADIENT:  -1.8889E+01  1.5619E+01  9.3756E+00 -1.0095E+01 -2.1216E+01  1.5203E+01 -3.6008E-01 -4.7364E-02 -1.3361E+01  2.2664E+01
             1.2271E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -816.177752050391        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.2207E+00  7.9423E-01  9.4832E-01  1.4075E+00  1.0875E+01  1.7841E+00  6.6978E-01  4.1496E-02  2.2703E+00  9.9448E+00
             8.6974E+00
 PARAMETER:  2.9944E-01 -1.3038E-01  4.6933E-02  4.4182E-01  2.4865E+00  6.7894E-01 -3.0081E-01 -3.0822E+00  9.1993E-01  2.3970E+00
             2.2630E+00
 GRADIENT:   2.5315E+01  2.2763E+01 -2.9719E-01  2.2476E+01 -7.4548E+00 -2.6412E+00  1.9167E+00  2.7245E-03 -2.9361E+00  1.8249E+01
            -3.1730E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -856.302281160074        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.4923E-01  1.4363E-01  1.7280E-01  9.9331E-01  3.5716E+01  1.3160E+00  1.0000E-02  1.0000E-02  1.0335E+00  6.0076E+00
             8.2348E+00
 PARAMETER:  4.7894E-02 -1.8405E+00 -1.6556E+00  9.3286E-02  3.6756E+00  3.7458E-01 -4.7985E+00 -9.3094E+00  1.3296E-01  1.8930E+00
             2.2084E+00
 GRADIENT:   7.6097E+01  6.0321E+01 -3.1947E+01  5.6824E+01 -1.6248E+00 -6.0219E+01  0.0000E+00  0.0000E+00 -4.2905E+01  3.0364E-01
            -1.3895E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -874.208933370043        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  7.0991E-01  3.9827E-02  5.1137E-02  4.9012E-01  1.1120E+02  1.2306E+00  1.0000E-02  1.0000E-02  6.3101E-01  4.1540E+00
             8.3754E+00
 PARAMETER: -2.4262E-01 -3.1232E+00 -2.8732E+00 -6.1311E-01  4.8114E+00  3.0753E-01 -1.0021E+01 -1.5038E+01 -3.6043E-01  1.5241E+00
             2.2253E+00
 GRADIENT:   2.0308E+02  7.4386E+00 -1.5716E+02  1.6019E+02  8.5980E-02 -6.7128E+01  0.0000E+00  0.0000E+00 -7.2036E+01  4.2183E-03
            -1.8647E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -928.846033197543        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  4.1578E-01  1.0000E-02  1.8478E-02  2.1953E-01  6.6294E+02  1.2861E+00  1.0000E-02  1.0000E-02  8.6873E-01  2.9836E+00
             8.7975E+00
 PARAMETER: -7.7759E-01 -4.5181E+00 -3.8911E+00 -1.4163E+00  6.5967E+00  3.5160E-01 -1.7798E+01 -2.2299E+01 -4.0719E-02  1.1931E+00
             2.2745E+00
 GRADIENT:   3.9626E+00  4.5862E-02 -3.3844E+01  4.0556E+01  3.5052E-03 -1.2877E+00  0.0000E+00  0.0000E+00 -4.7873E-01  3.9243E-07
            -1.3000E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -930.076281782559        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      816
 NPARAMETR:  3.9255E-01  1.0000E-02  1.6304E-02  1.9484E-01  6.8637E+02  1.2706E+00  1.0000E-02  1.0000E-02  8.5093E-01  2.8648E+00
             8.8638E+00
 PARAMETER: -8.3510E-01 -4.6495E+00 -4.0163E+00 -1.5356E+00  6.6314E+00  3.3947E-01 -1.8588E+01 -2.3020E+01 -6.1431E-02  1.1525E+00
             2.2820E+00
 GRADIENT:   8.5985E+01  0.0000E+00  1.2122E+02  5.7863E+01  1.3469E-03  3.1941E+00  0.0000E+00  0.0000E+00 -6.4953E-01  1.2078E-06
             2.0499E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -930.144291948520        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:      986
 NPARAMETR:  3.9530E-01  1.0000E-02  1.6236E-02  1.9311E-01  5.9630E+02  1.2832E+00  1.0000E-02  1.0000E-02  8.5971E-01  2.8725E+00
             8.8788E+00
 PARAMETER: -8.2811E-01 -4.6495E+00 -4.0205E+00 -1.5445E+00  6.4908E+00  3.4937E-01 -1.8588E+01 -2.3020E+01 -5.1164E-02  1.1552E+00
             2.2837E+00
 GRADIENT:   4.7657E+00  0.0000E+00 -1.1071E+00 -2.8538E+00  2.1166E-03  1.4660E+00  0.0000E+00  0.0000E+00  1.2135E+00  1.0617E-06
            -1.3344E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -930.176299651184        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1184             RESET HESSIAN, TYPE I
 NPARAMETR:  3.9430E-01  1.0000E-02  1.6135E-02  1.9288E-01  3.4820E+02  1.2769E+00  1.0000E-02  1.0000E-02  8.5289E-01  2.8017E+00
             8.8873E+00
 PARAMETER: -8.3064E-01 -4.6495E+00 -4.0268E+00 -1.5457E+00  5.9528E+00  3.4447E-01 -1.8588E+01 -2.3020E+01 -5.9127E-02  1.1302E+00
             2.2846E+00
 GRADIENT:   9.5811E+01  0.0000E+00  1.2075E+02  5.4327E+01  3.7896E-03  5.1529E+00  0.0000E+00  0.0000E+00  2.1220E-02  4.7264E-06
             2.1515E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -930.309824730919        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:     1335
 NPARAMETR:  3.8605E-01  1.0000E-02  1.5455E-02  1.8755E-01  4.5632E+01  1.2724E+00  1.0000E-02  1.0000E-02  8.5281E-01  2.2607E+00
             8.8673E+00
 PARAMETER: -8.5179E-01 -4.6495E+00 -4.0698E+00 -1.5737E+00  3.9206E+00  3.4094E-01 -1.8588E+01 -2.3020E+01 -5.9213E-02  9.1566E-01
             2.2824E+00
 GRADIENT:   9.3085E+01  0.0000E+00  1.1485E+02  6.8515E+01  1.9795E-02  4.7542E+00  0.0000E+00  0.0000E+00 -3.1210E-01  2.2734E-04
             1.9958E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -930.417345768982        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1510
 NPARAMETR:  3.8619E-01  1.0000E-02  1.5464E-02  1.8589E-01  4.2019E+01  1.2729E+00  1.0000E-02  1.0000E-02  8.4756E-01  1.9114E+00
             8.8771E+00
 PARAMETER: -8.5143E-01 -4.6495E+00 -4.0692E+00 -1.5826E+00  3.8381E+00  3.4130E-01 -1.8588E+01 -2.3020E+01 -6.5395E-02  7.4783E-01
             2.2835E+00
 GRADIENT:  -9.1392E-01  0.0000E+00 -2.9729E+00  1.2986E+00  1.3798E-02  1.1671E-02  0.0000E+00  0.0000E+00 -8.5760E-01  1.4901E-04
            -1.2957E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -930.431440766053        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1688
 NPARAMETR:  3.8658E-01  1.0000E-02  1.5460E-02  1.8577E-01  2.1155E+01  1.2726E+00  1.0000E-02  1.0000E-02  8.5160E-01  6.3208E-01
             8.8853E+00
 PARAMETER: -8.5041E-01 -4.6495E+00 -4.0695E+00 -1.5832E+00  3.1519E+00  3.4103E-01 -1.8588E+01 -2.3020E+01 -6.0640E-02 -3.5875E-01
             2.2844E+00
 GRADIENT:  -1.2302E+00  0.0000E+00 -1.4412E+00 -6.7541E-01 -3.7772E-03 -7.9938E-02  0.0000E+00  0.0000E+00  1.0549E-01  1.0755E-04
             3.2396E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -930.462589415027        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1885             RESET HESSIAN, TYPE I
 NPARAMETR:  3.8711E-01  1.0000E-02  1.5318E-02  1.8481E-01  1.7202E+01  1.2727E+00  1.0000E-02  1.0000E-02  8.5080E-01  5.5134E-01
             8.8947E+00
 PARAMETER: -8.4904E-01 -4.6495E+00 -4.0787E+00 -1.5885E+00  2.9450E+00  3.4115E-01 -1.8588E+01 -2.3020E+01 -6.1574E-02 -4.9541E-01
             2.2855E+00
 GRADIENT:   9.8090E+01  0.0000E+00  1.2545E+02  5.3723E+01 -1.6432E-02  5.1737E+00  0.0000E+00  0.0000E+00  1.8678E-01  2.0193E-04
             2.2583E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -930.509787874912        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     2016             RESET HESSIAN, TYPE I
 NPARAMETR:  3.8354E-01  1.0000E-02  1.5146E-02  1.8192E-01  2.4204E+01  1.2715E+00  1.0000E-02  1.0000E-02  8.5504E-01  4.8528E-01
             8.8840E+00
 PARAMETER: -8.5832E-01 -4.6495E+00 -4.0900E+00 -1.6042E+00  3.2865E+00  3.4018E-01 -1.8588E+01 -2.3020E+01 -5.6606E-02 -6.2302E-01
             2.2843E+00
 GRADIENT:   9.5233E+01  0.0000E+00  1.3732E+02  4.3004E+01 -3.1975E-03  5.1585E+00  0.0000E+00  0.0000E+00  1.3585E+00  5.4832E-05
             2.3286E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -930.675149905765        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:     2096
 NPARAMETR:  3.7894E-01  1.0000E-02  1.4433E-02  1.7496E-01  5.8718E+01  1.2732E+00  1.0000E-02  1.0000E-02  8.5421E-01  4.1310E-01
             8.9029E+00
 PARAMETER: -8.7038E-01 -4.6495E+00 -4.1382E+00 -1.6432E+00  4.1728E+00  3.4155E-01 -1.8588E+01 -2.3020E+01 -5.7577E-02 -7.8407E-01
             2.2864E+00
 GRADIENT:   1.0780E+02  0.0000E+00  1.3422E+02  4.5947E+01  1.1323E-02  6.5532E+00  0.0000E+00  0.0000E+00  1.4431E+00  5.8603E-06
             2.3346E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -930.724825708055        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     2265
 NPARAMETR:  3.7468E-01  1.0000E-02  1.4372E-02  1.7514E-01  5.5293E+01  1.2637E+00  1.0000E-02  1.0000E-02  8.4376E-01  4.1600E-01
             8.8993E+00
 PARAMETER: -8.8169E-01 -4.6495E+00 -4.1425E+00 -1.6422E+00  4.1126E+00  3.3404E-01 -1.8588E+01 -2.3020E+01 -6.9888E-02 -7.7706E-01
             2.2860E+00
 GRADIENT:  -1.5476E+00  0.0000E+00 -5.6256E+00  3.3325E+00  2.9359E-03 -7.1494E-01  0.0000E+00  0.0000E+00 -1.0289E+00  5.5131E-06
             1.0966E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -930.740508445060        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     2382
 NPARAMETR:  3.7229E-01  1.0000E-02  1.4199E-02  1.7400E-01  4.7994E+01  1.2648E+00  1.0000E-02  1.0000E-02  8.4913E-01  4.1492E-01
             8.8615E+00
 PARAMETER: -8.8808E-01 -4.6495E+00 -4.1546E+00 -1.6487E+00  3.9711E+00  3.3491E-01 -1.8588E+01 -2.3020E+01 -6.3538E-02 -7.7967E-01
             2.2817E+00
 GRADIENT:   9.6797E+01  0.0000E+00  1.2749E+02  6.2169E+01  1.4818E-03  4.6213E+00  0.0000E+00  0.0000E+00 -1.6927E-01  9.9532E-06
             2.0304E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -930.781982418334        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     2554
 NPARAMETR:  3.7333E-01  1.0000E-02  1.4209E-02  1.7312E-01  6.9928E+01  1.2659E+00  1.0000E-02  1.0000E-02  8.4737E-01  4.0625E-01
             8.8795E+00
 PARAMETER: -8.8529E-01 -4.6495E+00 -4.1539E+00 -1.6538E+00  4.3475E+00  3.3581E-01 -1.8588E+01 -2.3020E+01 -6.5623E-02 -8.0079E-01
             2.2837E+00
 GRADIENT:   7.0020E-01  0.0000E+00 -3.3243E+00 -8.7390E-01  2.5908E-03  5.0398E-02  0.0000E+00  0.0000E+00 -3.0580E-01  3.3757E-06
            -8.4507E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -930.790384526163        NO. OF FUNC. EVALS.: 208
 CUMULATIVE NO. OF FUNC. EVALS.:     2762             RESET HESSIAN, TYPE I
 NPARAMETR:  3.7339E-01  1.0000E-02  1.4167E-02  1.7271E-01  2.6338E+01  1.2658E+00  1.0000E-02  1.0000E-02  8.4899E-01  3.8865E-01
             8.8889E+00
 PARAMETER: -8.8512E-01 -4.6495E+00 -4.1568E+00 -1.6561E+00  3.3710E+00  3.3568E-01 -1.8588E+01 -2.3020E+01 -6.3702E-02 -8.4508E-01
             2.2848E+00
 GRADIENT:   9.9597E+01  0.0000E+00  1.3530E+02  5.0957E+01 -1.9177E-02  5.0669E+00  0.0000E+00  0.0000E+00  4.1226E-01  4.1610E-05
             2.3161E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -930.807520118775        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     2929
 NPARAMETR:  3.7293E-01  1.0000E-02  1.4116E-02  1.7223E-01  5.1977E+01  1.2654E+00  1.0000E-02  1.0000E-02  8.4902E-01  3.7067E-01
             8.8841E+00
 PARAMETER: -8.8636E-01 -4.6495E+00 -4.1604E+00 -1.6589E+00  4.0508E+00  3.3543E-01 -1.8588E+01 -2.3020E+01 -6.3669E-02 -8.9244E-01
             2.2843E+00
 GRADIENT:   2.2946E+00  0.0000E+00 -4.3249E+00 -6.8418E-01  1.7980E-03  4.2808E-02  0.0000E+00  0.0000E+00  8.5897E-02  5.7661E-06
            -5.3026E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -930.837710336541        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     3085             RESET HESSIAN, TYPE I
 NPARAMETR:  3.7122E-01  1.0000E-02  1.3995E-02  1.7071E-01  6.8597E+01  1.2648E+00  1.0000E-02  1.0000E-02  8.5012E-01  3.6448E-01
             8.8807E+00
 PARAMETER: -8.9096E-01 -4.6495E+00 -4.1690E+00 -1.6678E+00  4.3283E+00  3.3491E-01 -1.8588E+01 -2.3020E+01 -6.2373E-02 -9.0929E-01
             2.2839E+00
 GRADIENT:   1.0137E+02  0.0000E+00  1.3757E+02  4.9373E+01  6.8669E-04  5.0481E+00  0.0000E+00  0.0000E+00  6.6707E-01  3.8361E-06
             2.2420E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -930.965437076000        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:     3203
 NPARAMETR:  3.6560E-01  1.0000E-02  1.3480E-02  1.6542E-01  4.9167E+02  1.2634E+00  1.0000E-02  1.0000E-02  8.4725E-01  3.4412E-01
             8.8887E+00
 PARAMETER: -9.0622E-01 -4.6495E+00 -4.2065E+00 -1.6993E+00  6.2978E+00  3.3382E-01 -1.8588E+01 -2.3020E+01 -6.5764E-02 -9.6677E-01
             2.2848E+00
 GRADIENT:   2.2507E+00  0.0000E+00 -2.2770E+00 -4.0381E+00 -1.6709E-05  5.5519E-01  0.0000E+00  0.0000E+00  8.7325E-02  5.5669E-08
             2.2571E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -930.976624173803        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:     3310
 NPARAMETR:  3.6189E-01  1.0000E-02  1.3284E-02  1.6443E-01  1.4650E+04  1.2595E+00  1.0000E-02  1.0000E-02  8.4813E-01  3.3752E-01
             8.8641E+00
 PARAMETER: -9.1641E-01 -4.6495E+00 -4.2212E+00 -1.7053E+00  9.6922E+00  3.3071E-01 -1.8588E+01 -2.3020E+01 -6.4725E-02 -9.8612E-01
             2.2820E+00
 GRADIENT:   1.0100E+02  0.0000E+00  1.3234E+02  6.2343E+01 -9.4338E-06  4.6339E+00  0.0000E+00  0.0000E+00  4.5469E-02  9.2230E-11
             2.0836E+01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -931.023054620227        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     3447
 NPARAMETR:  3.5785E-01  1.0000E-02  1.3062E-02  1.6134E-01  7.3167E+07  1.2578E+00  1.0000E-02  1.0000E-02  8.4999E-01  3.4497E-01
             8.8538E+00
 PARAMETER: -9.2765E-01 -4.6495E+00 -4.2380E+00 -1.7243E+00  1.8208E+01  3.2934E-01 -1.8588E+01 -2.3020E+01 -6.2536E-02 -9.6429E-01
             2.2808E+00
 GRADIENT:   9.7647E+01  0.0000E+00  1.4288E+02  5.4093E+01 -9.1508E-09  4.5549E+00  0.0000E+00  0.0000E+00  7.6385E-01  0.0000E+00
             2.1359E+01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -931.058044053628        NO. OF FUNC. EVALS.: 152
 CUMULATIVE NO. OF FUNC. EVALS.:     3599
 NPARAMETR:  3.5911E-01  1.0000E-02  1.2944E-02  1.6087E-01  5.0007E+07  1.2585E+00  1.0000E-02  1.0000E-02  8.4585E-01  6.4376E-02
             8.8607E+00
 PARAMETER: -9.2412E-01 -4.6495E+00 -4.2471E+00 -1.7272E+00  1.7828E+01  3.2991E-01 -1.8588E+01 -2.3020E+01 -6.7408E-02 -2.6430E+00
             2.2816E+00
 GRADIENT:   1.0639E+02  0.0000E+00  1.3261E+02  6.2688E+01 -3.8810E-09  4.9218E+00  0.0000E+00  0.0000E+00 -3.0142E-01  0.0000E+00
             1.9661E+01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -931.058860217624        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     3741
 NPARAMETR:  3.5905E-01  1.0000E-02  1.2952E-02  1.6091E-01  4.9715E+07  1.2585E+00  1.0000E-02  1.0000E-02  8.4579E-01  6.6581E-02
             8.8616E+00
 PARAMETER: -9.2430E-01 -4.6495E+00 -4.2465E+00 -1.7269E+00  1.7822E+01  3.2993E-01 -1.8588E+01 -2.3020E+01 -6.7487E-02 -2.6093E+00
             2.2817E+00
 GRADIENT:   9.3326E-01  0.0000E+00 -1.2309E+01  7.6154E+00 -4.9062E-09  1.4330E-02  0.0000E+00  0.0000E+00 -5.3710E-01  0.0000E+00
            -2.9978E+00

0ITERATION NO.:  139    OBJECTIVE VALUE:  -931.076782858350        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     3880
 NPARAMETR:  3.5930E-01  1.0000E-02  1.2988E-02  1.6041E-01  2.7029E+07  1.2584E+00  1.0000E-02  1.0000E-02  8.4618E-01  6.4944E-02
             8.8816E+00
 PARAMETER: -9.2359E-01 -4.6495E+00 -4.2437E+00 -1.7300E+00  1.7212E+01  3.2985E-01 -1.8588E+01 -2.3020E+01 -6.7029E-02 -2.6342E+00
             2.2840E+00
 GRADIENT:   3.0690E-02  0.0000E+00 -2.8441E+00 -3.1280E+00 -1.6214E-08  5.4730E-02  0.0000E+00  0.0000E+00  1.8545E-02  0.0000E+00
            -1.2472E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3880
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9518E-03 -9.5197E-07  1.1073E-04 -1.8863E-02 -1.0922E-12
 SE:             2.8951E-02  2.8989E-06  2.8427E-04  2.4716E-02  2.8996E-12
 N:                     100         100         100         100         100

 P VAL.:         9.4625E-01  7.4262E-01  6.9690E-01  4.4534E-01  7.0642E-01

 ETASHRINKSD(%)  3.0092E+00  9.9990E+01  9.9048E+01  1.7198E+01  1.0000E+02
 ETASHRINKVR(%)  5.9278E+00  1.0000E+02  9.9991E+01  3.1438E+01  1.0000E+02
 EBVSHRINKSD(%)  3.2005E+00  9.9977E+01  9.9134E+01  1.8020E+01  1.0000E+02
 EBVSHRINKVR(%)  6.2985E+00  1.0000E+02  9.9992E+01  3.2793E+01  1.0000E+02
 RELATIVEINF(%)  1.7352E+00  3.2237E-06  2.0692E-05  2.0758E-01  0.0000E+00
 EPSSHRINKSD(%)  1.1919E+01
 EPSSHRINKVR(%)  2.2417E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -931.07678285835027     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -12.138249653677576     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    77.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -931.077       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.59E-01  1.00E-02  1.30E-02  1.60E-01  2.70E+07  1.26E+00  1.00E-02  1.00E-02  8.46E-01  6.49E-02  8.88E+00
 


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
+        5.02E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.34E+04  0.00E+00  4.62E+06
 
 TH 4
+       -3.70E+02  0.00E+00 -3.84E+05  3.76E+04
 
 TH 5
+       -5.55E-14  0.00E+00  6.38E-12  4.63E-13 -5.49E-22
 
 TH 6
+       -8.49E+00  0.00E+00 -1.77E+02 -3.60E+01 -1.89E-12  1.03E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -8.01E+00  0.00E+00  1.97E+03 -2.03E+02  4.84E-12 -3.59E+00  0.00E+00  0.00E+00  1.17E+02
 
 TH10
+        5.48E-04  0.00E+00  2.02E-02  6.49E-04  1.43E-12  3.08E-03  0.00E+00  0.00E+00  1.80E-02 -1.25E-02
 
 TH11
+       -3.92E+01  0.00E+00  7.56E+02 -3.15E+01 -2.05E-14  9.55E-01  0.00E+00  0.00E+00  4.44E+00  9.00E-05  6.30E+00
 
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
 #CPUT: Total CPU Time in Seconds,       86.954
Stop Time:
Thu Sep 30 02:58:55 CDT 2021
