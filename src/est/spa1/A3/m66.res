Thu Sep 30 00:27:19 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat66.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -775.863934600419        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2030E+02 -3.3141E+01  1.2947E+02 -1.2185E+02  9.6495E+01  3.2593E+01 -8.3969E+00 -8.9706E+01 -1.8057E+01 -6.4198E+01
            -2.3468E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1629.71286288057        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.4925E-01  1.0951E+00  1.0114E+00  1.1089E+00  9.9995E-01  8.7021E-01  9.3807E-01  9.2589E-01  8.9128E-01  8.9855E-01
             2.6486E+00
 PARAMETER:  4.7919E-02  1.9083E-01  1.1133E-01  2.0341E-01  9.9950E-02 -3.9022E-02  3.6065E-02  2.2999E-02 -1.5102E-02 -6.9774E-03
             1.0740E+00
 GRADIENT:   8.7764E+01  4.5000E+01  8.9185E+00  5.2738E+01 -2.1833E+01 -1.6359E+01 -7.2874E-01  4.8923E+00 -4.6969E+00  6.7238E+00
            -6.4544E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1635.78222084038        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.4149E-01  1.0458E+00  3.8728E-01  1.0919E+00  5.9371E-01  9.2608E-01  8.6570E-01  2.9400E-01  1.0281E+00  3.9704E-01
             2.6005E+00
 PARAMETER:  3.9710E-02  1.4480E-01 -8.4862E-01  1.8789E-01 -4.2136E-01  2.3207E-02 -4.4216E-02 -1.1242E+00  1.2773E-01 -8.2371E-01
             1.0557E+00
 GRADIENT:   5.0968E+01  2.4355E+01 -1.4557E+01  9.2484E+01  3.5764E+01  2.8436E+00 -1.3493E+01 -4.1886E-01  1.1772E+01  1.4662E+00
             9.6922E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1638.24890504926        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      336
 NPARAMETR:  9.3054E-01  1.0803E+00  3.9304E-01  1.0306E+00  6.2308E-01  9.3190E-01  9.1069E-01  2.8877E-01  9.9151E-01  3.8447E-01
             2.6295E+00
 PARAMETER:  2.8005E-02  1.7720E-01 -8.3384E-01  1.3015E-01 -3.7308E-01  2.9469E-02  6.4518E-03 -1.1421E+00  9.1473E-02 -8.5588E-01
             1.0668E+00
 GRADIENT:  -2.1730E+01 -1.7449E+01 -2.2484E+01  2.2435E+01  3.8180E+01  1.9166E+00 -8.7436E+00 -1.3003E-01  7.0173E+00  2.4712E+00
             1.2940E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1644.38096521111        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      517
 NPARAMETR:  9.4096E-01  7.6664E-01  3.9639E-01  1.1929E+00  4.7792E-01  9.2284E-01  1.2944E+00  7.6617E-01  8.2029E-01  1.2844E-01
             2.5250E+00
 PARAMETER:  3.9141E-02 -1.6573E-01 -8.2536E-01  2.7641E-01 -6.3831E-01  1.9700E-02  3.5807E-01 -1.6635E-01 -9.8100E-02 -1.9523E+00
             1.0263E+00
 GRADIENT:   6.0802E+00  2.9742E+01  1.2604E+01  2.4689E+01 -2.2236E+01 -2.4384E+00  1.8175E+00  1.7422E+00 -7.4295E+00  5.1943E-01
             1.6961E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1648.61356566975        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      692
 NPARAMETR:  9.3597E-01  4.8901E-01  2.7553E-01  1.2308E+00  3.2012E-01  9.3964E-01  1.4167E+00  8.7983E-01  9.2298E-01  3.9313E-02
             2.3753E+00
 PARAMETER:  3.3824E-02 -6.1538E-01 -1.1891E+00  3.0763E-01 -1.0391E+00  3.7746E-02  4.4832E-01 -2.8029E-02  1.9850E-02 -3.1362E+00
             9.6511E-01
 GRADIENT:  -2.4509E+00  6.1917E+00  8.8236E+00  1.0052E+01 -1.5298E+01  9.4397E-01 -9.9987E-01 -1.1311E+00 -2.1057E+00  3.8699E-02
             3.1123E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1657.87470570107        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      871
 NPARAMETR:  9.3382E-01  3.2301E-01  1.7018E-01  1.1427E+00  2.1726E-01  9.2774E-01  1.2967E+00  1.4201E+00  1.2496E+00  1.0839E-02
             2.0625E+00
 PARAMETER:  3.1526E-02 -1.0301E+00 -1.6709E+00  2.3336E-01 -1.4267E+00  2.4999E-02  3.5982E-01  4.5075E-01  3.2279E-01 -4.4246E+00
             8.2391E-01
 GRADIENT:   9.8728E+00 -6.5430E+00  2.8280E+00 -1.6334E+01  1.8485E+01 -6.2649E+00  1.3938E+01  7.4059E+00  2.1768E+00  4.3219E-03
             5.6754E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1662.64902304148        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1046
 NPARAMETR:  9.3062E-01  3.2509E-01  1.6858E-01  1.1586E+00  2.1637E-01  9.4262E-01  5.3053E-01  1.4430E+00  1.2813E+00  4.4636E-02
             2.0836E+00
 PARAMETER:  2.8100E-02 -1.0236E+00 -1.6803E+00  2.4719E-01 -1.4308E+00  4.0906E-02 -5.3388E-01  4.6672E-01  3.4789E-01 -3.0092E+00
             8.3410E-01
 GRADIENT:   2.5121E+00 -8.1657E-01  4.0200E+00  3.2113E+00 -4.7423E+00 -8.6896E-01  6.2660E-01  1.4385E+00  4.8033E-01 -1.0704E-03
            -2.4589E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1662.82474469807        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1221
 NPARAMETR:  9.2990E-01  3.2166E-01  1.6709E-01  1.1521E+00  2.1566E-01  9.4409E-01  2.3018E-01  1.4305E+00  1.2762E+00  1.6390E-01
             2.1157E+00
 PARAMETER:  2.7321E-02 -1.0343E+00 -1.6892E+00  2.4157E-01 -1.4341E+00  4.2469E-02 -1.3689E+00  4.5801E-01  3.4388E-01 -1.7085E+00
             8.4939E-01
 GRADIENT:  -2.5242E-01 -5.9358E+00 -1.3530E+00 -2.7909E+00  1.3024E+01 -7.5030E-02  3.4831E-02  1.6262E+00  1.8526E+00  2.3690E-02
             4.4936E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1663.04626638710        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1397
 NPARAMETR:  9.2976E-01  3.1548E-01  1.5846E-01  1.1406E+00  2.0798E-01  9.4511E-01  1.5998E-01  1.3879E+00  1.2804E+00  2.7810E-01
             2.1002E+00
 PARAMETER:  2.7172E-02 -1.0537E+00 -1.7423E+00  2.3154E-01 -1.4703E+00  4.3547E-02 -1.7327E+00  4.2780E-01  3.4720E-01 -1.1798E+00
             8.4201E-01
 GRADIENT:  -1.3186E-02 -1.8937E-01  2.8863E-01  5.7473E-01  4.0633E-01  7.7804E-02  3.3870E-02  5.4473E-01  4.6044E-01  1.7260E-01
            -3.7231E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1663.05864051825        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1574
 NPARAMETR:  9.2972E-01  3.1508E-01  1.5844E-01  1.1403E+00  2.0791E-01  9.4490E-01  7.3755E-02  1.3866E+00  1.2792E+00  2.7685E-01
             2.1038E+00
 PARAMETER:  2.7131E-02 -1.0549E+00 -1.7424E+00  2.3128E-01 -1.4706E+00  4.3325E-02 -2.5070E+00  4.2684E-01  3.4626E-01 -1.1843E+00
             8.4374E-01
 GRADIENT:  -2.4241E-01  1.9696E-01  3.9577E-01  4.3211E-01 -6.2610E-01  2.8695E-03  5.1540E-03  1.8503E-01  6.6703E-02  9.7279E-03
             2.2811E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1663.06354569453        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1752
 NPARAMETR:  9.2997E-01  3.1488E-01  1.5811E-01  1.1391E+00  2.0761E-01  9.4503E-01  1.0000E-02  1.3831E+00  1.2798E+00  2.8207E-01
             2.1031E+00
 PARAMETER:  2.7401E-02 -1.0556E+00 -1.7444E+00  2.3022E-01 -1.4721E+00  4.3465E-02 -4.7604E+00  4.2433E-01  3.4667E-01 -1.1656E+00
             8.4343E-01
 GRADIENT:   4.3215E-01  4.9316E-01  6.4489E-01 -1.5818E-01 -9.8945E-01  5.1915E-02  0.0000E+00 -5.5623E-02  1.9506E-01  1.2056E-02
            -1.9621E-01

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1663.06354569453        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1774
 NPARAMETR:  9.2997E-01  3.1488E-01  1.5811E-01  1.1391E+00  2.0761E-01  9.4503E-01  1.0000E-02  1.3831E+00  1.2798E+00  2.8207E-01
             2.1031E+00
 PARAMETER:  2.7401E-02 -1.0556E+00 -1.7444E+00  2.3022E-01 -1.4721E+00  4.3465E-02 -4.7604E+00  4.2433E-01  3.4667E-01 -1.1656E+00
             8.4343E-01
 GRADIENT:   4.3215E-01  4.9316E-01  6.4489E-01 -1.5818E-01 -9.8945E-01  5.1915E-02  0.0000E+00 -5.5623E-02  1.9506E-01  1.2056E-02
            -1.9621E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1774
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0473E-03 -3.2739E-04  1.2487E-02 -4.6967E-03  2.0270E-02
 SE:             2.9310E-02  1.6294E-04  2.5400E-02  2.7448E-02  1.1517E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4431E-01  4.4511E-02  6.2299E-01  8.6414E-01  7.8410E-02

 ETASHRINKSD(%)  1.8089E+00  9.9454E+01  1.4907E+01  8.0456E+00  6.1416E+01
 ETASHRINKVR(%)  3.5852E+00  9.9997E+01  2.7592E+01  1.5444E+01  8.5113E+01
 EBVSHRINKSD(%)  1.8364E+00  9.9457E+01  1.4454E+01  6.8520E+00  6.3091E+01
 EBVSHRINKVR(%)  3.6390E+00  9.9997E+01  2.6819E+01  1.3235E+01  8.6377E+01
 RELATIVEINF(%)  9.6275E+01  7.0050E-04  1.4390E+01  6.2579E+01  1.2171E+00
 EPSSHRINKSD(%)  3.0707E+01
 EPSSHRINKVR(%)  5.1985E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1663.0635456945345     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -744.12501248986177     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.56
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1663.064       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.30E-01  3.15E-01  1.58E-01  1.14E+00  2.08E-01  9.45E-01  1.00E-02  1.38E+00  1.28E+00  2.82E-01  2.10E+00
 


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
+        1.39E+03
 
 TH 2
+       -2.16E+01  3.50E+03
 
 TH 3
+       -8.39E+01  3.31E+03  1.68E+04
 
 TH 4
+        1.35E+00  1.55E+02 -1.08E+02  4.51E+02
 
 TH 5
+        9.02E+01 -9.97E+03 -2.20E+04 -9.79E+02  4.56E+04
 
 TH 6
+        4.16E+00 -2.20E+01  9.67E+00 -5.65E+00  4.18E+01  2.08E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        2.12E+00 -1.30E+01 -7.12E-01  2.61E-01 -4.05E+00  1.14E+00  0.00E+00  5.77E+01
 
 TH 9
+        1.32E+01 -7.19E+01  8.56E+01 -1.21E+01  4.31E+02 -3.71E-02  0.00E+00 -3.27E+00  8.49E+01
 
 TH10
+       -1.74E+00 -1.89E+02 -2.53E+01 -3.56E+00  6.37E+02  1.72E+00  0.00E+00  2.60E+01  2.74E+01  5.67E+01
 
 TH11
+       -2.13E+01 -4.00E+01 -2.81E+01  9.71E-01  8.04E+01  1.75E+00  0.00E+00  9.04E+00  8.79E+00  9.31E+00  9.38E+01
 
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
 #CPUT: Total CPU Time in Seconds,       35.104
Stop Time:
Thu Sep 30 00:27:56 CDT 2021
