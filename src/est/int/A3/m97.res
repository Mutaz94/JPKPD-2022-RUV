Fri Sep 24 22:54:40 CDT 2021
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
$DATA ../../../../data/int/A3/dat97.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1536.96759783356        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5491E+02  1.0350E+02  2.2952E+02 -1.3904E+02  1.6787E+02  1.2979E+00 -1.4588E+02 -2.8835E+02 -5.1817E+00 -1.3067E+02
            -3.9131E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2811.18729294742        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  8.7350E-01  1.0158E+00  9.2351E-01  1.1104E+00  9.0568E-01  9.1173E-01  9.6136E-01  1.0543E+00  7.5782E-01  7.6354E-01
             2.9952E+00
 PARAMETER: -3.5247E-02  1.1572E-01  2.0430E-02  2.0472E-01  9.3321E-04  7.5858E-03  6.0596E-02  1.5285E-01 -1.7731E-01 -1.6979E-01
             1.1970E+00
 GRADIENT:  -1.4076E+02  4.5391E+01 -5.4966E-01  2.6731E+01 -5.5905E+01 -8.9138E+00 -8.6203E+00  1.7733E+01 -1.2750E+00  1.3095E+01
             4.4007E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2821.29681481083        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  8.8274E-01  7.2456E-01  6.7248E-01  1.3154E+00  6.4180E-01  9.4304E-01  1.0868E+00  8.6850E-01  8.2200E-01  4.9956E-01
             2.8767E+00
 PARAMETER: -2.4730E-02 -2.2219E-01 -2.9679E-01  3.7412E-01 -3.4349E-01  4.1358E-02  1.8328E-01 -4.0982E-02 -9.6013E-02 -5.9402E-01
             1.1566E+00
 GRADIENT:  -1.0612E+02  5.4651E+01 -6.1776E+00  1.8157E+02 -2.5864E+01  3.1538E+00 -1.2062E+01  1.4611E+01  9.3620E+00 -1.6999E+00
             3.8641E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2881.09269114124        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.1459E-01  7.0643E-01  6.8455E-01  1.2176E+00  6.5626E-01  9.2122E-01  1.1029E+00  2.5888E-01  7.8439E-01  7.6337E-01
             2.3295E+00
 PARAMETER:  1.0717E-02 -2.4753E-01 -2.7900E-01  2.9687E-01 -3.2119E-01  1.7939E-02  1.9796E-01 -1.2514E+00 -1.4285E-01 -1.7001E-01
             9.4567E-01
 GRADIENT:   2.0187E+00  4.5200E+00 -3.1101E+00 -4.6258E+00  1.2739E+01 -4.6786E-01  4.8079E+00  1.5125E+00 -2.3103E+00  4.3645E+00
            -4.5737E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2890.93980530325        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.3687E-01  4.1723E-01  3.7010E-01  1.3129E+00  3.5934E-01  9.1770E-01  7.4213E-01  1.0000E-02  9.2670E-01  9.8214E-01
             2.2450E+00
 PARAMETER:  3.4789E-02 -7.7411E-01 -8.9398E-01  3.7221E-01 -9.2349E-01  1.4119E-02 -1.9823E-01 -4.8974E+00  2.3872E-02  8.1976E-02
             9.0871E-01
 GRADIENT:   6.1856E+01  6.0384E+01  1.8731E+01  1.2606E+02 -2.6927E+01 -3.8714E+00 -1.9708E+01  0.0000E+00  1.0768E+01  4.4532E+01
             5.0490E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2898.58525316527        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.2700E-01  2.8419E-01  2.2766E-01  1.2629E+00  2.3916E-01  9.3152E-01  6.9874E-01  1.0000E-02  1.0284E+00  8.0980E-01
             2.2144E+00
 PARAMETER:  2.4197E-02 -1.1581E+00 -1.3799E+00  3.3339E-01 -1.3306E+00  2.9065E-02 -2.5848E-01 -1.0918E+01  1.2800E-01 -1.1097E-01
             8.9499E-01
 GRADIENT:   3.0775E+01  4.5418E+01  3.6340E+01  1.6652E+02 -7.8399E+01  4.2427E-01 -4.2875E+01  0.0000E+00  3.8050E+00  5.6680E+00
             2.0023E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2910.30273982885        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  8.9633E-01  1.8307E-01  1.1776E-01  1.0713E+00  1.5522E-01  9.6077E-01  1.1062E+00  1.0000E-02  1.2697E+00  6.7047E-01
             2.0704E+00
 PARAMETER: -9.4422E-03 -1.5979E+00 -2.0391E+00  1.6883E-01 -1.7629E+00  5.9980E-02  2.0097E-01 -2.0717E+01  3.3880E-01 -2.9978E-01
             8.2776E-01
 GRADIENT:  -4.9352E+01  5.5508E-01  4.4772E+01  1.5799E+02 -9.9951E+01  6.6052E+00 -2.5961E+01  0.0000E+00 -2.4046E+01 -7.1729E+01
             2.2592E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2914.93262899660        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  8.8486E-01  1.5639E-01  9.1791E-02  9.6555E-01  1.3521E-01  9.6524E-01  1.2464E+00  1.0000E-02  1.4980E+00  6.9315E-01
             2.0009E+00
 PARAMETER: -2.2325E-02 -1.7554E+00 -2.2882E+00  6.4943E-02 -1.9009E+00  6.4620E-02  3.2023E-01 -2.4595E+01  5.0410E-01 -2.6650E-01
             7.9360E-01
 GRADIENT:  -7.2996E+01 -2.1051E+01  3.6455E+01  1.1804E+02 -7.2613E+01  6.8441E+00 -9.8915E+00  0.0000E+00 -2.6108E+01 -7.8396E+01
             3.2736E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2918.54597214336        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  8.8626E-01  1.4543E-01  8.0705E-02  8.9606E-01  1.2643E-01  9.6292E-01  1.3152E+00  1.0000E-02  1.6219E+00  7.3340E-01
             1.9746E+00
 PARAMETER: -2.0748E-02 -1.8281E+00 -2.4170E+00 -9.7442E-03 -1.9681E+00  6.2219E-02  3.7398E-01 -2.6658E+01  5.8359E-01 -2.1006E-01
             7.8034E-01
 GRADIENT:  -7.3778E+01 -3.7818E+01  7.7902E+00  9.3169E+01 -1.5655E+02  5.4616E+00 -2.8650E+00  0.0000E+00 -2.7826E+01 -6.6718E+01
            -5.4235E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2932.29906974654        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  8.9975E-01  1.5088E-01  8.3880E-02  8.4567E-01  1.3218E-01  9.4837E-01  1.3311E+00  1.0000E-02  1.6108E+00  7.8929E-01
             1.9725E+00
 PARAMETER: -5.6399E-03 -1.7913E+00 -2.3784E+00 -6.7629E-02 -1.9236E+00  4.6990E-02  3.8602E-01 -2.6809E+01  5.7676E-01 -1.3662E-01
             7.7929E-01
 GRADIENT:  -4.2526E+01 -2.1095E+01 -8.5163E-01  3.5672E+01 -5.2748E+01  2.6083E+00  2.3294E-01  0.0000E+00 -6.8016E+00 -3.4580E+01
            -5.5352E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2933.29237785224        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:      925
 NPARAMETR:  8.9970E-01  1.5102E-01  8.3778E-02  8.0229E-01  1.3230E-01  9.3970E-01  1.3314E+00  1.0000E-02  1.6113E+00  7.8923E-01
             1.9733E+00
 PARAMETER: -5.6905E-03 -1.7904E+00 -2.3796E+00 -1.2028E-01 -1.9226E+00  3.7809E-02  3.8622E-01 -2.6809E+01  5.7705E-01 -1.3669E-01
             7.7969E-01
 GRADIENT:  -4.4862E+01 -1.5866E+01  1.0017E+01 -8.3074E-01 -4.1198E+01  2.9584E-02  9.5014E-01  0.0000E+00 -4.2574E+00 -3.4777E+01
             1.0435E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2933.43546621402        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1108
 NPARAMETR:  8.9971E-01  1.5138E-01  8.3670E-02  8.1849E-01  1.3324E-01  9.4329E-01  1.3314E+00  1.0000E-02  1.6114E+00  7.8925E-01
             1.9733E+00
 PARAMETER: -5.6816E-03 -1.7880E+00 -2.3809E+00 -1.0029E-01 -1.9156E+00  4.1613E-02  3.8622E-01 -2.6809E+01  5.7713E-01 -1.3668E-01
             7.7970E-01
 GRADIENT:  -4.4061E+01 -1.6575E+01 -6.8832E+00  1.2259E+01 -2.0734E+01  1.0644E+00  5.7398E-01  0.0000E+00 -4.2298E+00 -3.4526E+01
             9.1298E-01

0ITERATION NO.:   58    OBJECTIVE VALUE:  -2933.53106865220        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:     1214
 NPARAMETR:  8.9980E-01  1.5122E-01  8.3835E-02  8.1494E-01  1.3329E-01  9.4251E-01  1.3309E+00  1.0000E-02  1.6105E+00  7.8936E-01
             1.9717E+00
 PARAMETER: -5.6785E-03 -1.7872E+00 -2.3813E+00 -1.0475E-01 -1.9133E+00  4.0694E-02  3.8622E-01 -2.6809E+01  5.7715E-01 -1.3667E-01
             7.7970E-01
 GRADIENT:  -1.2603E+05  1.4082E+04 -5.2865E+03 -2.4056E+05  1.3159E+04 -2.5196E+05  3.2620E+04  0.0000E+00  4.3594E+04 -9.2211E+04
             1.6158E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1214
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1522E-02  2.0781E-02  4.0403E-04 -1.2168E-02  1.2828E-02
 SE:             2.9279E-02  2.6362E-02  3.2753E-04  2.5620E-02  3.1045E-02
 N:                     100         100         100         100         100

 P VAL.:         4.6230E-01  4.3052E-01  2.1737E-01  6.3484E-01  6.7945E-01

 ETASHRINKSD(%)  1.9126E+00  1.1684E+01  9.8903E+01  1.4170E+01  1.0000E-10
 ETASHRINKVR(%)  3.7887E+00  2.2003E+01  9.9988E+01  2.6332E+01  1.0000E-10
 EBVSHRINKSD(%)  1.5001E+00  1.1650E+01  9.9230E+01  9.7135E+00  6.0048E+00
 EBVSHRINKVR(%)  2.9777E+00  2.1943E+01  9.9994E+01  1.8483E+01  1.1649E+01
 RELATIVEINF(%)  9.6959E+01  3.7229E+01  1.1061E-03  4.2371E+01  2.0008E+01
 EPSSHRINKSD(%)  2.0446E+01
 EPSSHRINKVR(%)  3.6711E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2933.5310686522025     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1279.4417088837918     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.16
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0T MATRIX UNOBTAINABLE
 Elapsed covariance  time in seconds:    14.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2933.531       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.00E-01  1.51E-01  8.36E-02  8.15E-01  1.34E-01  9.42E-01  1.33E+00  1.00E-02  1.61E+00  7.89E-01  1.97E+00
 


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
+        7.78E+08
 
 TH 2
+        9.08E+03  8.60E+07
 
 TH 3
+       -1.23E+04 -1.23E+04  1.58E+08
 
 TH 4
+        1.21E+05 -4.07E+04  5.29E+04  8.65E+08
 
 TH 5
+        9.57E+03 -5.71E+03  5.40E+05 -4.34E+04  9.66E+07
 
 TH 6
+       -3.47E+04  1.15E+04 -1.56E+04  7.83E+08  1.22E+04  7.09E+08
 
 TH 7
+       -5.77E+02  3.11E+02 -1.49E+02 -2.13E+04  1.27E+02  6.08E+03  2.38E+07
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.04E+05 -6.77E+04  9.24E+04 -1.18E+04 -7.16E+04  3.36E+03  6.05E+01  0.00E+00  7.26E+06
 
 TH10
+       -1.68E+04  5.53E+03 -7.22E+03  1.01E+05  5.91E+03 -2.90E+04 -4.92E+02  0.00E+00  1.63E+03  5.41E+08
 
 TH11
+        4.84E+03 -1.63E+03  2.06E+03 -7.10E+03 -1.62E+03  2.04E+03  4.66E+01  0.00E+00 -1.19E+04  1.01E+03  2.66E+06
 
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
 #CPUT: Total CPU Time in Seconds,       42.819
Stop Time:
Fri Sep 24 22:55:40 CDT 2021
