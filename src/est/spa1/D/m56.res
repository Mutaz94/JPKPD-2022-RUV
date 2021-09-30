Thu Sep 30 03:16:13 CDT 2021
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
$DATA ../../../../data/spa1/D/dat56.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22088.9534689376        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8888E+02  4.3824E+02 -7.2141E+01  2.3036E+02  2.6339E+02 -2.4789E+03 -7.7478E+02 -6.1557E+01 -1.5978E+03 -6.6014E+02
            -4.2360E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -585.573985425951        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3444E+00  8.2408E-01  8.4520E-01  1.7800E+00  1.7465E+00  3.9841E+00  1.3690E+00  9.2107E-01  1.9906E+00  1.0632E+00
             1.2882E+01
 PARAMETER:  3.9594E-01 -9.3490E-02 -6.8188E-02  6.7659E-01  6.5763E-01  1.4823E+00  4.1407E-01  1.7781E-02  7.8844E-01  1.6125E-01
             2.6559E+00
 GRADIENT:   1.1889E+01  1.5986E+01 -2.7336E+01  6.4752E+00 -1.0770E+01  1.6445E+02  3.2541E+00  4.6223E+00 -2.2855E+01  2.2218E+00
             1.2326E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -640.606593533202        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3526E+00  7.4241E-01  8.0861E-01  1.9951E+00  5.9435E+00  2.4791E+00  2.1346E+00  4.1511E-01  3.0979E+00  2.2141E+00
             1.2364E+01
 PARAMETER:  4.0202E-01 -1.9785E-01 -1.1244E-01  7.9068E-01  1.8823E+00  1.0079E+00  8.5826E-01 -7.7921E-01  1.2307E+00  8.9483E-01
             2.6148E+00
 GRADIENT:   2.5777E+01  1.4950E+01 -1.8990E+01  4.6344E+01 -1.4064E+01  1.5768E+01  1.3212E+01  8.6557E-02  7.4139E+01  1.0250E+00
             1.2298E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -715.115704484642        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0098E+00  1.7945E-01  8.2120E-01  1.8341E+00  1.1216E+01  1.8582E+00  2.2474E+00  3.8124E+00  1.1833E+00  9.1960E-01
             1.1607E+01
 PARAMETER:  1.0979E-01 -1.6179E+00 -9.6990E-02  7.0657E-01  2.5173E+00  7.1958E-01  9.0975E-01  1.4383E+00  2.6831E-01  1.6187E-02
             2.5516E+00
 GRADIENT:  -7.9194E+01  2.6461E+01  6.7142E+00  1.3241E+02 -1.4180E+00 -3.2436E+01  1.7396E+00  1.0514E+01  2.0922E+00  2.9977E-02
             8.8801E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -783.928481538980        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  8.1953E-01  2.5358E-02  1.2540E-01  8.3260E-01  1.4445E+02  1.5348E+00  2.0539E-01  2.2768E+00  3.6211E-01  4.6362E-01
             1.0207E+01
 PARAMETER: -9.9021E-02 -3.5747E+00 -1.9763E+00 -8.3206E-02  5.0729E+00  5.2839E-01 -1.4829E+00  9.2277E-01 -9.1580E-01 -6.6869E-01
             2.4231E+00
 GRADIENT:   5.1553E+01 -7.9530E-01  2.2741E+01 -1.1387E+01  1.3607E-02 -4.4610E+01  2.1133E-04  1.1707E+01  5.8540E+00 -3.8617E-07
             2.0024E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -784.630328020068        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  7.7386E-01  1.9272E-02  9.3295E-02  7.4749E-01  2.0801E+02  1.5284E+00  1.4091E-01  2.0534E+00  3.0344E-01  4.2395E-01
             1.0018E+01
 PARAMETER: -1.5637E-01 -3.8491E+00 -2.2720E+00 -1.9103E-01  5.4376E+00  5.2423E-01 -1.8596E+00  8.1951E-01 -1.0926E+00 -7.5814E-01
             2.4044E+00
 GRADIENT:   6.0751E+01 -5.1979E-01 -2.8030E+00  4.0842E+01  1.7829E-02 -5.1062E+01  5.9046E-05  1.1698E+01  3.9814E+00 -1.6942E-07
            -4.9311E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -784.638199639112        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:      457
 NPARAMETR:  7.6267E-01  1.8346E-02  8.8189E-02  7.2821E-01  2.2364E+02  1.5282E+00  1.3045E-01  2.0146E+00  2.9200E-01  4.1929E-01
             9.9775E+00
 PARAMETER: -1.7094E-01 -3.8983E+00 -2.3283E+00 -2.1717E-01  5.5100E+00  5.2410E-01 -1.9368E+00  8.0040E-01 -1.1310E+00 -7.6919E-01
             2.4003E+00
 GRADIENT:   6.1704E+01 -4.9863E-01 -6.1483E+00  4.7920E+01  1.7094E-02 -5.1757E+01  4.6252E-05  1.1797E+01  3.6687E+00 -1.3530E-07
            -8.9126E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -793.572790770948        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  7.4544E-01  1.4676E-02  9.3644E-02  7.4330E-01  2.1565E+02  1.8877E+00  1.6034E-01  1.8066E+00  2.1862E-01  3.4016E-01
             9.8194E+00
 PARAMETER: -1.9378E-01 -4.1215E+00 -2.2683E+00 -1.9666E-01  5.4737E+00  7.3535E-01 -1.7305E+00  6.9145E-01 -1.4204E+00 -9.7833E-01
             2.3844E+00
 GRADIENT:   1.5073E+01  1.8494E-02 -2.7760E+01  5.6392E+01  1.8242E-02  1.4211E+00  3.7749E-05  1.0965E+00  1.4556E+00 -2.8993E-08
            -3.8705E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -798.325487207584        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  6.3189E-01  1.0000E-02  5.9040E-02  5.4337E-01  3.0850E+02  1.8656E+00  2.1743E-01  1.4428E+00  4.9097E-02  8.3009E-02
             1.0048E+01
 PARAMETER: -3.5903E-01 -5.0471E+00 -2.7295E+00 -5.0997E-01  5.8317E+00  7.2361E-01 -1.4259E+00  4.6656E-01 -2.9140E+00 -2.3888E+00
             2.4074E+00
 GRADIENT:   1.4684E+00  0.0000E+00 -3.8648E+00  5.0023E+00  4.0903E-03  6.6229E-01  3.2397E-05 -2.2437E-01  7.2669E-02  2.0322E-09
            -1.9199E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -798.361337484044        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  6.3144E-01  1.0000E-02  5.9507E-02  5.4410E-01  3.0059E+02  1.8620E+00  2.0636E-01  1.4485E+00  1.6864E-02  7.9994E-02
             1.0056E+01
 PARAMETER: -3.5976E-01 -5.0514E+00 -2.7217E+00 -5.0862E-01  5.8057E+00  7.2168E-01 -1.4781E+00  4.7050E-01 -3.9826E+00 -2.4258E+00
             2.4082E+00
 GRADIENT:   1.6526E-01  0.0000E+00 -1.0561E-01 -1.8313E-01  3.4704E-03  1.4932E-01  2.9231E-05 -3.4348E-03  8.7899E-03  2.0855E-09
            -1.9980E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -798.362004667639        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1169
 NPARAMETR:  6.3124E-01  1.0000E-02  5.9597E-02  5.4454E-01  2.8987E+02  1.8603E+00  1.5336E-01  1.4494E+00  1.0000E-02  7.9991E-02
             1.0057E+01
 PARAMETER: -3.6007E-01 -5.0514E+00 -2.7201E+00 -5.0782E-01  5.7694E+00  7.2073E-01 -1.7750E+00  4.7116E-01 -4.5080E+00 -2.4258E+00
             2.4082E+00
 GRADIENT:  -2.8385E-01  0.0000E+00  1.1053E-01 -2.5412E-01  3.5019E-03 -1.7630E-01  1.6184E-05 -1.4851E-02  2.0098E-03  2.2355E-09
             1.3662E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -798.378363591359        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  6.3010E-01  1.0000E-02  5.9854E-02  5.4442E-01  1.8237E+01  1.8511E+00  1.3416E-01  1.4454E+00  1.0000E-02  7.9893E-02
             1.0063E+01
 PARAMETER: -3.6188E-01 -5.0514E+00 -2.7158E+00 -5.0804E-01  3.0034E+00  7.1578E-01 -1.9087E+00  4.6841E-01 -4.6404E+00 -2.4271E+00
             2.4088E+00
 GRADIENT:  -3.1033E+00  0.0000E+00  4.5822E+00 -5.7119E+00  3.8515E-02 -1.9881E+00  1.2246E-05 -2.0817E-01  0.0000E+00  5.8437E-07
             2.1621E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -798.453791045007        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1534
 NPARAMETR:  6.3402E-01  1.0000E-02  5.9553E-02  5.4574E-01  9.6491E+00  1.8655E+00  1.3013E-01  1.4470E+00  1.0000E-02  7.9180E-02
             1.0044E+01
 PARAMETER: -3.5567E-01 -5.0514E+00 -2.7209E+00 -5.0562E-01  2.3669E+00  7.2352E-01 -1.9392E+00  4.6951E-01 -4.8647E+00 -2.4360E+00
             2.4070E+00
 GRADIENT:  -5.1112E-01  0.0000E+00 -4.6742E-02 -3.5890E-01 -6.8709E-02  3.9988E-01  1.1388E-05  3.9282E-01  0.0000E+00  2.8705E-05
            -7.1630E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -798.455006302503        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1591
 NPARAMETR:  6.3415E-01  1.0000E-02  5.9573E-02  5.4575E-01  1.0270E+01  1.8635E+00  1.2970E-01  1.4449E+00  1.0000E-02  7.8546E-02
             1.0052E+01
 PARAMETER: -3.5547E-01 -5.0514E+00 -2.7205E+00 -5.0560E-01  2.4293E+00  7.2247E-01 -1.9425E+00  4.6804E-01 -4.8461E+00 -2.4441E+00
             2.4078E+00
 GRADIENT:  -3.7462E-01  0.0000E+00  4.8296E-02 -3.9453E-01  1.5688E-02  4.3830E-02  1.1345E-05  1.5568E-01  0.0000E+00  7.7348E-06
            -1.7577E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1591
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4805E-03 -1.4309E-05  7.3701E-03 -4.3288E-04 -7.2559E-06
 SE:             2.9060E-02  7.9618E-06  2.3471E-02  2.4027E-04  9.7998E-06
 N:                     100         100         100         100         100

 P VAL.:         7.7042E-01  7.2311E-02  7.5351E-01  7.1598E-02  4.5905E-01

 ETASHRINKSD(%)  2.6456E+00  9.9973E+01  2.1370E+01  9.9195E+01  9.9967E+01
 ETASHRINKVR(%)  5.2213E+00  1.0000E+02  3.8174E+01  9.9994E+01  1.0000E+02
 EBVSHRINKSD(%)  2.9447E+00  9.9964E+01  2.0849E+01  9.9138E+01  9.9966E+01
 EBVSHRINKVR(%)  5.8026E+00  1.0000E+02  3.7351E+01  9.9993E+01  1.0000E+02
 RELATIVEINF(%)  1.4476E+01  6.7972E-07  8.3482E-01  6.5676E-05  3.5155E-06
 EPSSHRINKSD(%)  9.1496E+00
 EPSSHRINKVR(%)  1.7462E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -798.45500630250274     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       120.48352690216996     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.22
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -798.455       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.34E-01  1.00E-02  5.96E-02  5.46E-01  1.03E+01  1.86E+00  1.30E-01  1.44E+00  1.00E-02  7.85E-02  1.01E+01
 


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
+        7.37E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.40E+02  0.00E+00  8.94E+04
 
 TH 4
+       -2.76E+02  0.00E+00 -1.42E+04  2.64E+03
 
 TH 5
+        2.00E-01  0.00E+00 -2.58E+00  3.77E-01  4.50E-03
 
 TH 6
+        3.38E+00  0.00E+00  7.09E+01 -2.40E+01  6.69E-03  4.90E+01
 
 TH 7
+        3.24E-05  0.00E+00  1.50E-03 -6.22E-04 -2.67E-05  2.30E-03 -1.06E-02
 
 TH 8
+        2.19E+00  0.00E+00 -1.36E+02 -3.59E+01 -1.69E-02  3.98E+00  2.72E-04  3.65E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -1.51E-04  0.00E+00  2.89E-02  6.99E-03 -2.04E-04 -3.16E-03  4.13E-03  1.08E-04  0.00E+00 -5.53E-04
 
 TH11
+       -1.13E+01  0.00E+00  1.18E+02 -1.73E+01 -6.86E-03  4.53E-01  9.16E-05  2.03E+00  0.00E+00 -7.79E-05  4.75E+00
 
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
 #CPUT: Total CPU Time in Seconds,       39.262
Stop Time:
Thu Sep 30 03:16:54 CDT 2021
