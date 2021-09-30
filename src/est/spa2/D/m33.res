Thu Sep 30 08:56:03 CDT 2021
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
$DATA ../../../../data/spa2/D/dat33.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17462.8158684529        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9997E+02  2.3674E+02  2.2080E+01  9.0122E+01  1.8698E+02 -1.5758E+03 -7.1716E+02 -3.1430E+01 -1.1774E+03 -4.8384E+02
            -3.5464E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -681.000321859199        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3853E+00  1.2549E+00  9.2435E-01  1.6399E+00  1.0494E+00  2.4372E+00  1.6916E+00  9.7448E-01  1.5434E+00  1.1453E+00
             1.3721E+01
 PARAMETER:  4.2595E-01  3.2707E-01  2.1336E-02  5.9462E-01  1.4819E-01  9.9085E-01  6.2565E-01  7.4147E-02  5.3400E-01  2.3568E-01
             2.7189E+00
 GRADIENT:   1.4884E+01 -2.7941E+01 -2.9729E+01  1.4997E+01  4.0809E+01  7.4613E+01 -8.6566E+00  4.3955E+00 -1.1028E+01  1.6745E+01
             3.2785E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -733.276493133434        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.3888E+00  1.2077E+00  4.7661E+00  2.0141E+00  2.1084E+00  2.3243E+00  5.7613E+00  5.7300E-01  3.3028E+00  1.9543E+00
             1.2137E+01
 PARAMETER:  4.2845E-01  2.8875E-01  1.6615E+00  8.0016E-01  8.4595E-01  9.4340E-01  1.8512E+00 -4.5686E-01  1.2948E+00  7.7006E-01
             2.5962E+00
 GRADIENT:   2.8350E+01  7.2294E+00 -1.6981E+01  1.9969E+01 -9.2222E+00  2.9464E+01  6.5266E+01  1.6172E-01  7.9842E+01  2.7655E+01
             3.1457E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -837.164842154787        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0790E+00  1.7336E+00  1.1529E+01  9.4193E-01  3.1364E+00  2.2941E+00  3.1750E+00  1.6421E-01  1.3077E+00  1.3190E+00
             7.7782E+00
 PARAMETER:  1.7599E-01  6.5017E-01  2.5448E+00  4.0178E-02  1.2431E+00  9.3032E-01  1.2553E+00 -1.7066E+00  3.6827E-01  3.7685E-01
             2.1513E+00
 GRADIENT:  -2.6480E+01 -3.4138E+00 -7.7143E-01 -8.8824E+00  2.2500E+01  1.2995E+01  2.7739E+00  2.5444E-04  6.2145E+00  1.0907E+01
            -2.4019E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -850.167552082612        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      323
 NPARAMETR:  1.1403E+00  1.6279E+00  6.1390E+00  9.6783E-01  2.1169E+00  2.1239E+00  3.3575E+00  1.2405E-01  1.0223E+00  4.6615E-01
             8.0344E+00
 PARAMETER:  2.3128E-01  5.8726E-01  1.9147E+00  6.7305E-02  8.4998E-01  8.5326E-01  1.3112E+00 -1.9871E+00  1.2207E-01 -6.6324E-01
             2.1837E+00
 GRADIENT:  -1.0744E+01 -8.7220E+00  2.3070E-01 -5.8713E+00  1.5641E+00 -3.2460E+01 -2.3717E+01  1.2126E-03  2.5719E+00  2.3759E+00
            -2.1551E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -859.175466910975        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  1.1917E+00  1.2888E+00  6.5389E+00  1.2021E+00  2.0181E+00  2.4833E+00  4.2638E+00  1.5004E-01  9.7712E-01  1.9384E-01
             8.0304E+00
 PARAMETER:  2.7542E-01  3.5371E-01  1.9778E+00  2.8409E-01  8.0217E-01  1.0096E+00  1.5501E+00 -1.7968E+00  7.6851E-02 -1.5407E+00
             2.1832E+00
 GRADIENT:  -1.1243E+00  1.7465E+00  5.2906E-01  2.5557E+00 -6.4594E-01 -8.9783E-01 -1.1051E+00  3.3814E-03 -1.7125E+00  4.2697E-01
            -2.5323E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -859.583228935793        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      676
 NPARAMETR:  1.2159E+00  7.9476E-01  7.2182E+00  1.5428E+00  1.9276E+00  2.4977E+00  5.2207E+00  3.9685E-01  1.3510E+00  4.7189E-02
             8.1147E+00
 PARAMETER:  2.9551E-01 -1.2971E-01  2.0766E+00  5.3358E-01  7.5627E-01  1.0154E+00  1.7526E+00 -8.2420E-01  4.0082E-01 -2.9536E+00
             2.1937E+00
 GRADIENT:   5.4037E+00  4.0792E-01 -4.5866E-01 -3.3186E+00 -1.6355E-01  2.1382E+00  2.6620E+00  4.3309E-02  1.1547E+00  2.5485E-02
             9.9915E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -859.698585668094        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      854
 NPARAMETR:  1.1867E+00  6.7243E-01  8.7428E+00  1.6412E+00  1.9729E+00  2.4809E+00  5.4642E+00  4.0292E-01  1.4041E+00  3.2090E-02
             8.0356E+00
 PARAMETER:  2.7121E-01 -2.9686E-01  2.2682E+00  5.9540E-01  7.7950E-01  1.0086E+00  1.7982E+00 -8.0902E-01  4.3936E-01 -3.3392E+00
             2.1839E+00
 GRADIENT:  -1.7956E+00  1.9542E-01 -3.7936E-01  1.7606E+00  6.5518E-01  1.0712E-01 -7.9979E-01  3.5627E-02 -3.2482E-01  1.1465E-02
            -2.1194E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -859.740566746243        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      994
 NPARAMETR:  1.1954E+00  6.6644E-01  9.6877E+00  1.6371E+00  1.9753E+00  2.4984E+00  5.4097E+00  1.4622E-01  1.4137E+00  1.0055E-02
             8.0407E+00
 PARAMETER:  2.7851E-01 -3.0581E-01  2.3709E+00  5.9294E-01  7.8070E-01  1.0156E+00  1.7882E+00 -1.8227E+00  4.4619E-01 -4.4997E+00
             2.1845E+00
 GRADIENT:   7.8693E-01 -8.0999E-01  2.6194E-01 -2.0997E-01 -1.7307E+00  2.8440E+00 -2.9506E+00  3.7246E-03 -7.4592E-03  6.8368E-04
            -6.3826E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -859.757917655780        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1133
 NPARAMETR:  1.1922E+00  6.8553E-01  9.2434E+00  1.6365E+00  1.9883E+00  2.4943E+00  5.4083E+00  6.2350E-02  1.4107E+00  1.0000E-02
             8.0416E+00
 PARAMETER:  2.7582E-01 -2.7756E-01  2.3239E+00  5.9254E-01  7.8727E-01  1.0140E+00  1.7879E+00 -2.6750E+00  4.4410E-01 -4.5662E+00
             2.1846E+00
 GRADIENT:  -3.2602E-01  1.6316E-01 -1.3963E-01  1.7848E+00  2.9217E-01  2.0502E+00 -1.4977E+00  7.4045E-04 -2.5093E-02  0.0000E+00
            -1.3191E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -859.768028118451        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1321
 NPARAMETR:  1.1941E+00  6.8865E-01  9.4119E+00  1.6299E+00  1.9898E+00  2.5030E+00  5.3944E+00  1.0000E-02  1.4040E+00  1.0000E-02
             8.0441E+00
 PARAMETER:  2.7742E-01 -2.7302E-01  2.3420E+00  5.8852E-01  7.8804E-01  1.0175E+00  1.7854E+00 -5.8960E+00  4.3932E-01 -4.5662E+00
             2.1849E+00
 GRADIENT:   2.1199E-01 -9.0127E-02 -1.9628E-02  9.5975E-01 -7.0819E-02  3.2970E+00 -1.6600E+00  0.0000E+00 -1.6427E-01  0.0000E+00
            -5.7848E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -859.771253707425        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1519             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1942E+00  6.9394E-01  9.4062E+00  1.6212E+00  1.9935E+00  2.5033E+00  5.3725E+00  1.0000E-02  1.4056E+00  1.0000E-02
             8.0471E+00
 PARAMETER:  2.7747E-01 -2.6537E-01  2.3414E+00  5.8315E-01  7.8988E-01  1.0176E+00  1.7813E+00 -5.8960E+00  4.4044E-01 -4.5662E+00
             2.1853E+00
 GRADIENT:   2.1550E+01  1.5121E+00  2.7962E-02  1.7597E+01  1.6884E+00  5.3038E+01  5.4239E+01  0.0000E+00  2.5118E+00  0.0000E+00
             2.6792E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -859.773534911013        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1684
 NPARAMETR:  1.1943E+00  6.9936E-01  9.4777E+00  1.6225E+00  1.9902E+00  2.5034E+00  5.3653E+00  1.0000E-02  1.4003E+00  1.0000E-02
             8.0438E+00
 PARAMETER:  2.7753E-01 -2.5759E-01  2.3489E+00  5.8397E-01  7.8826E-01  1.0176E+00  1.7800E+00 -5.8960E+00  4.3666E-01 -4.5662E+00
             2.1849E+00
 GRADIENT:   2.3650E-01 -1.6977E-02  5.5982E-02  9.3924E-01 -4.1268E-01  3.3206E+00 -1.6634E+00  0.0000E+00 -1.8886E-01  0.0000E+00
            -6.0838E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1684
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.8726E-03  4.6505E-02  3.0461E-06 -7.6883E-02 -5.7981E-06
 SE:             2.9261E-02  2.2368E-02  6.1958E-06  1.5898E-02  6.0927E-05
 N:                     100         100         100         100         100

 P VAL.:         7.6172E-01  3.7611E-02  6.2298E-01  1.3261E-06  9.2418E-01

 ETASHRINKSD(%)  1.9724E+00  2.5064E+01  9.9979E+01  4.6740E+01  9.9796E+01
 ETASHRINKVR(%)  3.9060E+00  4.3846E+01  1.0000E+02  7.1634E+01  1.0000E+02
 EBVSHRINKSD(%)  2.1991E+00  2.5179E+01  9.9962E+01  4.2523E+01  9.9688E+01
 EBVSHRINKVR(%)  4.3499E+00  4.4019E+01  1.0000E+02  6.6964E+01  9.9999E+01
 RELATIVEINF(%)  9.4966E+01  2.1904E+01  1.5599E-06  1.2514E+01  1.0214E-04
 EPSSHRINKSD(%)  1.1181E+01
 EPSSHRINKVR(%)  2.1112E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -859.77353491101326     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       242.95270493459384     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    38.17
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    11.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -859.774       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  6.99E-01  9.48E+00  1.62E+00  1.99E+00  2.50E+00  5.37E+00  1.00E-02  1.40E+00  1.00E-02  8.04E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        4.46E+01
 
 TH 2
+        4.92E+00  1.38E+01
 
 TH 3
+        3.45E-02  4.70E-02  1.69E-04
 
 TH 4
+       -3.34E+00  2.99E+01  9.58E-02  6.93E+01
 
 TH 5
+       -1.44E+00 -4.61E+00 -1.57E-02 -1.01E+01  1.55E+00
 
 TH 6
+        8.52E+00 -3.07E+00 -6.19E-03 -9.82E+00  1.07E+00  2.89E+00
 
 TH 7
+        2.02E+00 -1.50E+00 -4.03E-03 -4.10E+00  5.15E-01  9.15E-01  3.18E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.52E+00 -8.95E+00 -2.76E-02 -2.16E+01  3.03E+00  3.51E+00  1.38E+00  0.00E+00  6.87E+00
 
 TH10
+        1.01E-73 -7.41E-74 -1.99E-76 -2.02E-73  2.54E-74  4.51E-74  1.57E-74  0.00E+00  6.80E-74  7.8E-128
 
 TH11
+       -9.60E-01 -2.63E+00 -8.20E-03 -5.79E+00  8.63E-01  7.21E-01  3.03E-01  0.00E+00  1.75E+00  1.46E-54  8.88E-01
 
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
+        1.17E+02
 
 TH 2
+       -2.58E+00  4.53E+01
 
 TH 3
+        2.34E-02  1.18E-01  3.10E-02
 
 TH 4
+       -2.13E+00  3.74E+01  2.04E-04  8.53E+01
 
 TH 5
+       -1.50E+00 -8.18E+00 -7.15E-01 -9.08E+00  2.57E+01
 
 TH 6
+        4.44E-01 -1.40E+00  1.81E-02  9.92E-01 -7.24E-01  2.75E+01
 
 TH 7
+       -1.82E-01  5.90E+00 -1.72E-02 -4.90E+00  5.21E-01 -4.83E-01  3.22E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.00E-01 -2.65E+00 -1.07E-01 -2.24E+01  3.58E+00 -4.01E-01  2.02E+00  0.00E+00  2.17E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.33E+01
 
 TH11
+       -4.65E+00 -3.27E+00  7.27E-03 -6.99E+00  4.77E-01  1.05E+00  2.95E-01  0.00E+00  2.07E+00  0.00E+00  1.01E+01
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.19E+02
 
 TH 2
+        6.49E+01  4.69E+01
 
 TH 3
+        1.44E-01  1.45E-01  7.40E-03
 
 TH 4
+        5.40E+01  3.74E+01  1.81E-01  9.00E+01
 
 TH 5
+       -1.18E+01 -9.35E+00 -2.71E-01 -1.46E+01  1.18E+01
 
 TH 6
+        3.73E+01  2.53E+01  2.65E-02 -4.65E+00 -2.03E+00  4.49E+01
 
 TH 7
+        8.31E+00  5.96E+00 -8.83E-03 -5.58E+00  4.77E-01  9.25E+00  3.04E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -3.26E+00 -2.85E+00 -2.16E-02 -1.95E+01  2.68E+00  3.05E+00  1.68E+00  0.00E+00  1.47E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.70E+01 -1.65E+01 -8.68E-02 -1.61E+01  4.90E+00  6.56E-01 -1.74E+00  0.00E+00 -1.19E+00  0.00E+00  2.45E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       48.750
Stop Time:
Thu Sep 30 08:57:06 CDT 2021
