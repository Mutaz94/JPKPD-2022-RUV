Sat Sep 25 02:02:52 CDT 2021
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
$DATA ../../../../data/int/SL3/dat23.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      980
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

 TOT. NO. OF OBS RECS:      880
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -236.141181810878        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0365E+01  1.7002E+02  1.0373E+02  5.8353E+01  7.3531E+01  1.7435E+01 -1.3254E+02 -2.4312E+02 -2.9481E+01 -3.7167E+01
            -6.6617E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2651.17033228991        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0357E+00  1.1288E+00  1.0414E+00  9.2172E-01  1.1014E+00  8.8573E-01  1.1326E+00  9.9502E-01  7.6231E-01  9.3450E-01
             2.9444E+00
 PARAMETER:  1.3505E-01  2.2111E-01  1.4058E-01  1.8482E-02  1.9655E-01 -2.1340E-02  2.2455E-01  9.5003E-02 -1.7140E-01  3.2260E-02
             1.1799E+00
 GRADIENT:   7.2668E-01  6.3396E+00 -7.2194E+00 -2.4387E+00 -3.3143E+00 -2.4255E+01  1.5356E+00  2.8632E+00  5.6655E-01 -6.3424E+00
             1.0792E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2653.02040747093        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0372E+00  1.2188E+00  1.3422E+00  8.8921E-01  1.2764E+00  9.2921E-01  1.1074E+00  9.9449E-01  7.4357E-01  1.0914E+00
             2.9413E+00
 PARAMETER:  1.3653E-01  2.9786E-01  3.9430E-01 -1.7421E-02  3.4403E-01  2.6574E-02  2.0203E-01  9.4475E-02 -1.9629E-01  1.8745E-01
             1.1789E+00
 GRADIENT:   5.7252E+00  1.5248E+01  1.8816E-01  2.1554E+01  7.6939E+00 -5.0926E+00  5.9759E+00 -4.4933E-02  6.5520E-01 -2.1740E+00
             1.0254E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2658.57360589703        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0299E+00  1.3146E+00  1.8387E+00  8.2257E-01  1.4812E+00  9.4391E-01  1.0039E+00  2.5696E+00  7.3078E-01  1.1904E+00
             2.7359E+00
 PARAMETER:  1.2942E-01  3.7353E-01  7.0903E-01 -9.5316E-02  4.9286E-01  4.2279E-02  1.0387E-01  1.0438E+00 -2.1364E-01  2.7432E-01
             1.1065E+00
 GRADIENT:  -3.2779E+00 -7.8882E+00 -6.6771E+00  4.4038E+00  1.3710E+01 -2.7443E-01 -2.7478E+00 -3.2514E-01  1.7360E+00 -4.5207E-01
            -1.8742E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2661.59883663539        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  1.0308E+00  1.2901E+00  2.6895E+00  8.4313E-01  1.5613E+00  9.4385E-01  1.0770E+00  3.2439E+00  6.1420E-01  1.2246E+00
             2.7401E+00
 PARAMETER:  1.3034E-01  3.5469E-01  1.0894E+00 -7.0636E-02  5.4554E-01  4.2216E-02  1.7420E-01  1.2768E+00 -3.8744E-01  3.0265E-01
             1.1080E+00
 GRADIENT:  -4.6449E-03 -2.1915E+00  9.3733E-02 -9.1561E+00 -5.4131E+00  3.6628E-01  8.5619E-01 -5.8085E-01  4.4119E-01 -8.6862E-01
            -3.7928E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2661.64691628983        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      433
 NPARAMETR:  1.0309E+00  1.2885E+00  2.6982E+00  8.4918E-01  1.5615E+00  9.4305E-01  1.0791E+00  3.2604E+00  6.1078E-01  1.2247E+00
             2.7413E+00
 PARAMETER:  1.3039E-01  3.5352E-01  1.0926E+00 -6.3482E-02  5.4563E-01  4.1359E-02  1.7612E-01  1.2818E+00 -3.9301E-01  3.0272E-01
             1.1084E+00
 GRADIENT:  -3.0498E-01  4.5422E+00 -6.1522E-01 -6.7470E-01 -4.3220E+00 -1.8963E-02  4.4095E-01  7.2515E-03  2.0544E-01 -5.5220E-01
            -2.9696E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2661.64976443257        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      610
 NPARAMETR:  1.0309E+00  1.2886E+00  2.6982E+00  8.5012E-01  1.5615E+00  9.4508E-01  1.0788E+00  3.2736E+00  6.1047E-01  1.2247E+00
             2.7413E+00
 PARAMETER:  1.3039E-01  3.5352E-01  1.0926E+00 -6.2375E-02  5.4564E-01  4.3509E-02  1.7588E-01  1.2859E+00 -3.9353E-01  3.0272E-01
             1.1084E+00
 GRADIENT:  -9.0298E+00  1.0300E+00 -1.2881E+00 -1.2487E-01 -6.3746E+00  3.2014E-02  4.0712E-02 -1.7119E-03  1.6634E-02 -6.6824E-01
            -4.7707E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2661.70221155497        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      796             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0349E+00  1.2896E+00  2.6960E+00  8.4855E-01  1.5734E+00  9.4482E-01  1.0791E+00  3.2691E+00  6.0651E-01  1.2347E+00
             2.7482E+00
 PARAMETER:  1.3429E-01  3.5435E-01  1.0918E+00 -6.4226E-02  5.5321E-01  4.3235E-02  1.7610E-01  1.2845E+00 -4.0003E-01  3.1082E-01
             1.1110E+00
 GRADIENT:   8.9848E+00  2.6459E+00 -1.7744E+00  9.6601E-01  1.2231E+00  7.3649E-01  3.9672E-01  2.7169E-01  2.2997E-01  1.3443E-01
             2.5330E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2661.70364866236        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      952            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0349E+00  1.2915E+00  2.6945E+00  8.4859E-01  1.5752E+00  9.4478E-01  1.0789E+00  3.2745E+00  6.0492E-01  1.2352E+00
             2.7476E+00
 PARAMETER:  1.3427E-01  3.5579E-01  1.0912E+00 -6.4174E-02  5.5438E-01  4.3197E-02  1.7595E-01  1.2862E+00 -4.0266E-01  3.1125E-01
             1.1107E+00
 GRADIENT:   8.8670E+00  4.4015E+00 -2.0171E+00  3.1476E+00  1.9438E+00  6.9601E-01  3.3961E-01  3.7980E-01  1.4255E-01  7.7015E-02
             1.8668E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2661.70379961807        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  1.0347E+00  1.2912E+00  2.6938E+00  8.4741E-01  1.5748E+00  9.4464E-01  1.0783E+00  3.2729E+00  6.0460E-01  1.2351E+00
             2.7482E+00
 PARAMETER:  1.3414E-01  3.5556E-01  1.0909E+00 -6.5574E-02  5.5415E-01  4.3045E-02  1.7538E-01  1.2857E+00 -4.0318E-01  3.1116E-01
             1.1110E+00
 GRADIENT:  -3.1881E-01 -2.1877E+00 -2.2340E+00 -1.9649E-01 -6.3788E-01 -5.9031E-02  5.1066E-02 -9.4426E-02  8.6509E-03 -1.5332E-01
             5.1379E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2661.70427735138        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1246
 NPARAMETR:  1.0347E+00  1.2930E+00  2.6923E+00  8.4744E-01  1.5758E+00  9.4471E-01  1.0779E+00  3.2754E+00  6.0464E-01  1.2362E+00
             2.7498E+00
 PARAMETER:  1.3411E-01  3.5693E-01  1.0904E+00 -6.5540E-02  5.5476E-01  4.3123E-02  1.7503E-01  1.2864E+00 -4.0312E-01  3.1205E-01
             1.1115E+00
 GRADIENT:  -5.3128E-01 -5.6119E-01 -2.3694E+00  1.6586E+00 -4.2127E-01 -4.1820E-02  9.2207E-02 -1.0944E-02 -7.9231E-03 -4.5942E-02
             1.6539E+00

0ITERATION NO.:   53    OBJECTIVE VALUE:  -2661.70651683433        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     1351
 NPARAMETR:  1.0351E+00  1.2938E+00  2.6922E+00  8.4651E-01  1.5767E+00  9.4486E-01  1.0774E+00  3.2767E+00  6.0481E-01  1.2366E+00
             2.7498E+00
 PARAMETER:  1.3445E-01  3.5759E-01  1.0904E+00 -6.6531E-02  5.5532E-01  4.3289E-02  1.7464E-01  1.2868E+00 -4.0276E-01  3.1241E-01
             1.1115E+00
 GRADIENT:  -8.6557E+03  3.2541E+03  1.0625E+03  1.1183E+00  2.0925E+03  1.4075E-02  1.4218E-01 -2.9871E-02  1.4476E-02  3.7246E+03
            -1.0512E+03
 NUMSIGDIG:         3.3         3.3         3.3         1.8         3.3         3.2         2.1         3.6         2.6         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1351
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5339E-03 -5.4636E-03 -3.2180E-02 -2.6852E-03 -2.9946E-02
 SE:             2.9390E-02  2.4867E-02  1.7623E-02  1.5047E-02  2.1886E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5838E-01  8.2610E-01  6.7843E-02  8.5837E-01  1.7123E-01

 ETASHRINKSD(%)  1.5414E+00  1.6691E+01  4.0962E+01  4.9590E+01  2.6678E+01
 ETASHRINKVR(%)  3.0590E+00  3.0597E+01  6.5145E+01  7.4588E+01  4.6239E+01
 EBVSHRINKSD(%)  1.7529E+00  1.7205E+01  4.3562E+01  5.1522E+01  2.4166E+01
 EBVSHRINKVR(%)  3.4751E+00  3.1449E+01  6.8148E+01  7.6499E+01  4.2492E+01
 RELATIVEINF(%)  9.6414E+01  7.0760E+00  1.1077E+01  2.2790E+00  2.5237E+01
 EPSSHRINKSD(%)  1.6718E+01
 EPSSHRINKVR(%)  3.0642E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          880
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1617.3318184402240     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2661.7065168343283     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1044.3746983941044     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2661.707       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.29E+00  2.69E+00  8.47E-01  1.58E+00  9.45E-01  1.08E+00  3.28E+00  6.05E-01  1.24E+00  2.75E+00
 


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
+        1.50E+07
 
 TH 2
+       -4.52E+06  1.36E+06
 
 TH 3
+        3.96E+01  6.79E+01  3.36E+04
 
 TH 4
+       -1.06E+04  3.71E+03  4.78E+02  9.74E+02
 
 TH 5
+       -2.38E+06  7.17E+05 -3.41E+02  1.74E+03  3.79E+05
 
 TH 6
+        1.98E+02 -6.24E+01 -8.70E+00 -7.46E+00 -3.29E+01  2.09E+02
 
 TH 7
+       -1.01E+03  3.27E+02  4.83E+01 -3.98E+01  1.63E+02 -1.41E+00  8.62E+01
 
 TH 8
+       -2.82E+03  8.47E+02  1.32E+02 -8.17E+05  4.48E+02 -4.58E-03 -3.68E+05  1.65E+04
 
 TH 9
+       -1.78E+03  5.21E+02  8.39E+01 -3.57E+01  2.88E+02  4.75E-01  3.39E+01  1.51E+00  3.60E+01
 
 TH10
+       -5.41E+06  1.63E+06 -9.04E+01  3.83E+03  8.59E+05 -6.88E+01  3.66E+02  1.02E+03  6.45E+02  1.95E+06
 
 TH11
+       -5.34E+01 -7.72E+01 -3.54E+02 -4.95E+02  3.10E+02  1.19E+01 -4.10E+01 -1.26E+02 -7.54E+01  9.29E+01  3.16E+04
 
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
 #CPUT: Total CPU Time in Seconds,       49.264
Stop Time:
Sat Sep 25 02:03:43 CDT 2021
