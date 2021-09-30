Wed Sep 29 02:58:36 CDT 2021
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
$DATA ../../../../data/int/SL2/dat27.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      998
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

 TOT. NO. OF OBS RECS:      898
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
 RAW OUTPUT FILE (FILE): m27.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1352.65664288594        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6300E+02  8.4038E+01  1.7548E+02  1.3441E+02  2.1708E+02  4.8529E+01 -7.7703E+01 -2.8192E+02 -6.5749E+01 -1.2884E+01
            -4.5788E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2852.21120741267        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0475E+00  1.2409E+00  9.7963E-01  9.1735E-01  1.0335E+00  9.2551E-01  1.0706E+00  9.2794E-01  1.0718E+00  8.1779E-01
             2.0453E+00
 PARAMETER:  1.4638E-01  3.1587E-01  7.9415E-02  1.3730E-02  1.3297E-01  2.2585E-02  1.6818E-01  2.5209E-02  1.6937E-01 -1.0115E-01
             8.1555E-01
 GRADIENT:   2.8084E+02  1.4862E+02 -2.3868E+00  6.2933E+01 -1.5465E+01 -1.5029E+01  1.0925E+00 -2.9436E-01 -1.8188E+01 -1.2594E+01
            -3.9882E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2857.20203887968        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0362E+00  1.0464E+00  7.7519E-01  1.0141E+00  8.5876E-01  1.0484E+00  1.3137E+00  3.6418E-01  1.0153E+00  5.1532E-01
             2.1031E+00
 PARAMETER:  1.3553E-01  1.4534E-01 -1.5464E-01  1.1398E-01 -5.2266E-02  1.4728E-01  3.7288E-01 -9.1011E-01  1.1516E-01 -5.6297E-01
             8.4342E-01
 GRADIENT:   2.0988E+02  8.1054E+01 -3.7146E+01  8.8279E+01  5.7800E+01  3.7590E+01  5.2678E+00 -1.8740E+00 -2.8699E+01 -1.1494E+01
            -3.3781E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2886.23337979974        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.7184E-01  1.0545E+00  8.7614E-01  9.6449E-01  9.1963E-01  9.5100E-01  1.0896E+00  1.1744E-01  1.0688E+00  7.8412E-01
             2.3676E+00
 PARAMETER:  7.1431E-02  1.5307E-01 -3.2234E-02  6.3842E-02  1.6219E-02  4.9758E-02  1.8579E-01 -2.0418E+00  1.6650E-01 -1.4320E-01
             9.6188E-01
 GRADIENT:   1.0029E+01  6.1777E+00  2.6759E+00  4.2790E+00  4.6810E+00  5.9152E-01 -5.8826E+00 -3.4823E-02  7.7670E-01 -2.0559E+00
            -1.6002E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2890.25189474661        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      406
 NPARAMETR:  1.0164E+00  1.4406E+00  1.0673E+00  7.9475E-01  1.2467E+00  9.6528E-01  8.6869E-01  1.4424E-01  1.2090E+00  1.1293E+00
             2.4158E+00
 PARAMETER:  1.1623E-01  4.6509E-01  1.6510E-01 -1.2972E-01  3.2054E-01  6.4661E-02 -4.0768E-02 -1.8363E+00  2.8975E-01  2.2162E-01
             9.8204E-01
 GRADIENT:   3.2876E+01  2.6127E+01  3.9052E+00  4.1334E+01  6.5337E+00 -1.5354E-01  8.5550E-01 -1.3520E-01 -1.9325E+00  5.7792E+00
             1.3057E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2893.68527719220        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0008E+00  1.7520E+00  9.3173E-01  5.5860E-01  1.4502E+00  9.6725E-01  7.3852E-01  2.3665E-01  1.5612E+00  1.2319E+00
             2.3844E+00
 PARAMETER:  1.0076E-01  6.6073E-01  2.9289E-02 -4.8233E-01  4.7172E-01  6.6700E-02 -2.0311E-01 -1.3412E+00  5.4547E-01  3.0855E-01
             9.6896E-01
 GRADIENT:  -1.3073E+00 -2.3332E+00  3.5417E+00 -4.0841E-01 -5.6663E+00  9.5202E-01 -1.7966E+00 -2.4070E-01  1.3348E+00 -3.7584E-01
            -4.8350E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2894.42085163364        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      760
 NPARAMETR:  1.0029E+00  2.0570E+00  6.6806E-01  3.6598E-01  1.6497E+00  9.6418E-01  7.0688E-01  3.2298E-01  1.9918E+00  1.3617E+00
             2.3782E+00
 PARAMETER:  1.0290E-01  8.2124E-01 -3.0337E-01 -9.0519E-01  6.0058E-01  6.3524E-02 -2.4689E-01 -1.0302E+00  7.8905E-01  4.0872E-01
             9.6636E-01
 GRADIENT:   2.8303E+00  8.7201E+00 -7.8251E-01  4.3822E+00  2.0956E+00 -6.2439E-01  9.7121E-01 -6.8882E-02 -9.6499E-01  2.0094E+00
            -2.6615E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2894.48241401932        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      936
 NPARAMETR:  1.0015E+00  2.1213E+00  5.9882E-01  3.2157E-01  1.6912E+00  9.6722E-01  6.9295E-01  4.0411E-01  2.1654E+00  1.3702E+00
             2.3777E+00
 PARAMETER:  1.0146E-01  8.5205E-01 -4.1279E-01 -1.0345E+00  6.2542E-01  6.6667E-02 -2.6679E-01 -8.0606E-01  8.7259E-01  4.1495E-01
             9.6614E-01
 GRADIENT:  -5.6195E-01  5.8189E+00 -1.8713E+00  3.9020E+00  2.8888E+00  4.7248E-01 -3.2753E-01  2.9687E-02 -2.3306E-01 -8.2754E-01
            -1.4694E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2894.49449383052        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1111
 NPARAMETR:  1.0013E+00  2.1547E+00  5.6873E-01  2.9691E-01  1.7143E+00  9.6701E-01  6.8857E-01  4.4793E-01  2.2684E+00  1.3854E+00
             2.3784E+00
 PARAMETER:  1.0126E-01  8.6764E-01 -4.6435E-01 -1.1143E+00  6.3899E-01  6.6449E-02 -2.7314E-01 -7.0311E-01  9.1908E-01  4.2599E-01
             9.6641E-01
 GRADIENT:  -1.0064E+00  1.9834E+00 -1.6112E+00  2.3701E+00  2.1730E+00  3.8600E-01 -3.1208E-01  9.3735E-02 -1.2400E-01 -5.8019E-01
             4.2340E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2894.52923992325        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  1.0016E+00  2.1639E+00  5.7389E-01  2.8713E-01  1.7198E+00  9.6630E-01  6.8796E-01  3.2273E-01  2.3061E+00  1.3923E+00
             2.3780E+00
 PARAMETER:  1.0158E-01  8.7191E-01 -4.5532E-01 -1.1478E+00  6.4221E-01  6.5719E-02 -2.7403E-01 -1.0309E+00  9.3557E-01  4.3092E-01
             9.6625E-01
 GRADIENT:  -1.5513E-01 -3.9917E-01  6.1118E-01 -7.4071E-01 -8.9998E-01  1.6561E-01 -1.3860E-01  1.7574E-02 -1.0894E+00 -3.9538E-02
             6.1664E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2894.54736532082        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1466
 NPARAMETR:  1.0016E+00  2.1682E+00  5.6099E-01  2.8357E-01  1.7212E+00  9.6588E-01  6.8797E-01  1.2933E-01  2.3425E+00  1.3926E+00
             2.3766E+00
 PARAMETER:  1.0157E-01  8.7389E-01 -4.7806E-01 -1.1603E+00  6.4302E-01  6.5287E-02 -2.7402E-01 -1.9454E+00  9.5123E-01  4.3116E-01
             9.6566E-01
 GRADIENT:  -8.2456E-02 -2.5891E+00  1.4908E-01 -1.1796E-01 -9.2015E-01 -1.2846E-02  1.0846E-01  4.6555E-03  2.5072E-03 -1.1474E-01
            -4.9952E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2894.55006528558        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1649             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0018E+00  2.1686E+00  5.5863E-01  2.8300E-01  1.7226E+00  9.6597E-01  6.8762E-01  1.1348E-02  2.3462E+00  1.3938E+00
             2.3770E+00
 PARAMETER:  1.0175E-01  8.7410E-01 -4.8227E-01 -1.1623E+00  6.4381E-01  6.5381E-02 -2.7451E-01 -4.3787E+00  9.5278E-01  4.3205E-01
             9.6584E-01
 GRADIENT:   8.1424E+01  3.9545E+02  5.0214E-01  1.9535E+01  3.8016E+01  7.8343E+00  5.9512E+00  8.2819E-05  9.3680E+00  4.5391E+00
             1.7322E+01

0ITERATION NO.:   59    OBJECTIVE VALUE:  -2894.55012418704        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1779
 NPARAMETR:  1.0018E+00  2.1686E+00  5.5835E-01  2.8324E-01  1.7226E+00  9.6596E-01  6.8763E-01  1.0000E-02  2.3466E+00  1.3937E+00
             2.3770E+00
 PARAMETER:  1.0175E-01  8.7408E-01 -4.8277E-01 -1.1615E+00  6.4381E-01  6.5367E-02 -2.7450E-01 -4.8997E+00  9.5298E-01  4.3198E-01
             9.6583E-01
 GRADIENT:   3.1064E-01 -3.4098E+00 -3.7651E-02  1.2606E-01 -1.5876E-01  1.8183E-02  2.5083E-02  0.0000E+00  1.6965E-01  3.9641E-02
            -5.5223E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1779
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3292E-03 -3.4004E-02 -2.9286E-05  3.5780E-02 -2.3729E-02
 SE:             2.9486E-02  2.3799E-02  4.7900E-05  2.0144E-02  2.5823E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6404E-01  1.5305E-01  5.4094E-01  7.5701E-02  3.5814E-01

 ETASHRINKSD(%)  1.2175E+00  2.0272E+01  9.9840E+01  3.2515E+01  1.3491E+01
 ETASHRINKVR(%)  2.4201E+00  3.6434E+01  1.0000E+02  5.4458E+01  2.5162E+01
 EBVSHRINKSD(%)  1.4107E+00  1.8790E+01  9.9853E+01  3.7044E+01  1.0994E+01
 EBVSHRINKVR(%)  2.8014E+00  3.4050E+01  1.0000E+02  6.0365E+01  2.0780E+01
 RELATIVEINF(%)  9.7150E+01  1.2293E+01  1.2250E-04  7.2416E+00  3.7703E+01
 EPSSHRINKSD(%)  1.6592E+01
 EPSSHRINKVR(%)  3.0431E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2894.5501241870356     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1244.1365185514435     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    45.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2894.550       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.17E+00  5.58E-01  2.83E-01  1.72E+00  9.66E-01  6.88E-01  1.00E-02  2.35E+00  1.39E+00  2.38E+00
 


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
+        1.15E+03
 
 TH 2
+       -1.15E+01  3.48E+02
 
 TH 3
+       -1.20E-01  4.57E+01  8.85E+01
 
 TH 4
+       -1.51E+01  3.56E+02 -1.16E+02  1.09E+03
 
 TH 5
+       -3.56E+00 -5.39E+01 -2.77E+01  1.08E+02  1.53E+02
 
 TH 6
+        4.67E+00 -3.00E+00  5.40E-01 -5.93E+00 -1.11E+00  2.01E+02
 
 TH 7
+        2.47E+00 -9.30E+00  4.39E-01 -1.66E+01 -6.34E+00 -5.96E-01  2.10E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.52E-01 -7.14E+00 -6.72E+00  5.33E+01  2.01E-01 -1.56E-01  7.69E+00  0.00E+00  1.14E+01
 
 TH10
+        5.43E-01 -1.11E+01  1.20E+00  1.83E+01 -1.18E+01  5.15E-01  1.03E+00  0.00E+00  2.37E+00  6.15E+01
 
 TH11
+       -1.53E+01 -1.37E+01 -6.32E-01 -1.69E+01  1.18E+00  2.39E+00  9.57E+00  0.00E+00  1.81E+00  5.96E+00  2.07E+02
 
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
 #CPUT: Total CPU Time in Seconds,       58.230
Stop Time:
Wed Sep 29 02:59:36 CDT 2021
