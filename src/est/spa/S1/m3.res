Sat Sep 18 10:51:41 CDT 2021
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
$DATA ../../../../data/spa/S1/dat3.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1695.84063129620        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -8.4420E+00 -7.4817E+01 -2.7349E+01 -7.4290E+01  3.3108E+01  5.1981E+00  3.6701E+00  3.5500E+00  2.3622E+00  2.6112E+01
            -2.0616E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1702.67000821671        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0332E+00  1.0908E+00  1.0197E+00  1.0195E+00  1.0064E+00  1.0062E+00  8.9129E-01  1.0799E+00  1.0133E+00  6.4609E-01
             1.0528E+00
 PARAMETER:  1.3266E-01  1.8692E-01  1.1954E-01  1.1930E-01  1.0635E-01  1.0618E-01 -1.5087E-02  1.7686E-01  1.1325E-01 -3.3681E-01
             1.5148E-01
 GRADIENT:   6.8867E+01  2.8901E+00 -8.8486E+00  3.8170E+01  5.4583E+01  7.0941E+00 -6.0732E+00 -9.5493E+00 -3.3156E-01 -1.7583E+01
            -1.3283E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1703.65632911655        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0436E+00  1.1967E+00  9.4683E-01  9.5473E-01  9.9458E-01  1.0032E+00  9.4690E-01  1.3355E+00  1.0174E+00  5.9143E-01
             1.0601E+00
 PARAMETER:  1.4267E-01  2.7954E-01  4.5362E-02  5.3668E-02  9.4568E-02  1.0320E-01  4.5433E-02  3.8929E-01  1.1721E-01 -4.2522E-01
             1.5832E-01
 GRADIENT:   9.2340E+01  3.0089E+01 -3.1412E+00  3.6378E+01  1.2068E+01  4.5898E+00  2.9814E+00  1.5935E+00 -2.2968E+00 -9.5819E+00
            -7.6718E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1705.14800827282        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      238
 NPARAMETR:  1.0124E+00  1.1428E+00  1.0826E+00  9.7440E-01  1.0373E+00  9.9382E-01  8.6493E-01  1.3760E+00  1.0392E+00  7.2536E-01
             1.0628E+00
 PARAMETER:  1.1236E-01  2.3345E-01  1.7934E-01  7.4068E-02  1.3662E-01  9.3797E-02 -4.5106E-02  4.1919E-01  1.3846E-01 -2.2108E-01
             1.6092E-01
 GRADIENT:   1.5757E+01  3.3969E+00 -1.6517E+00  7.6181E+00  7.3195E+00  1.8421E+00 -3.2315E-01  8.2218E-01  4.3266E-01 -3.6305E+00
            -1.8513E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1705.14818840681        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      313
 NPARAMETR:  1.0121E+00  1.1424E+00  1.0805E+00  9.7424E-01  1.0364E+00  9.9357E-01  8.6547E-01  1.3699E+00  1.0389E+00  7.2715E-01
             1.0628E+00
 PARAMETER:  1.1201E-01  2.3315E-01  1.7746E-01  7.3905E-02  1.3579E-01  9.3547E-02 -4.4477E-02  4.1474E-01  1.3820E-01 -2.1862E-01
             1.6093E-01
 GRADIENT:   1.4801E+01  3.1946E+00 -1.5444E+00  7.1576E+00  6.8650E+00  1.7297E+00 -3.0536E-01  7.7242E-01  4.0281E-01 -3.4106E+00
            -1.7402E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1706.12264748371        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  1.0321E+00  1.3391E+00  9.3152E-01  8.5700E-01  1.0563E+00  1.0102E+00  7.8209E-01  1.2696E+00  1.1425E+00  7.8052E-01
             1.0653E+00
 PARAMETER:  1.3162E-01  3.9199E-01  2.9060E-02 -5.4313E-02  1.5477E-01  1.1010E-01 -1.4579E-01  3.3870E-01  2.3324E-01 -1.4780E-01
             1.6329E-01
 GRADIENT:   1.0365E+01  2.0851E+01  7.3403E+00  1.4444E+01 -1.7110E+01  1.6170E+00 -1.6575E+00 -1.2725E+00 -2.4819E+00  1.5641E+00
             6.3109E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1709.33377959159        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  1.0200E+00  1.8462E+00  4.2774E-01  5.0195E-01  1.1097E+00  1.0038E+00  6.8375E-01  8.4159E-01  1.6129E+00  7.2045E-01
             1.0650E+00
 PARAMETER:  1.1984E-01  7.1314E-01 -7.4923E-01 -5.8926E-01  2.0413E-01  1.0377E-01 -2.8017E-01 -7.2459E-02  5.7801E-01 -2.2788E-01
             1.6300E-01
 GRADIENT:  -1.8679E+01  1.3872E+01 -1.1475E+00  1.1874E+01  1.3992E+00 -1.9245E+00 -2.4174E+00  7.2889E-01  1.1570E+00 -2.2817E+00
            -9.1927E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1709.46621845534        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:      846
 NPARAMETR:  1.0208E+00  1.8698E+00  4.0883E-01  4.8689E-01  1.1127E+00  1.0046E+00  6.8267E-01  8.1296E-01  1.6364E+00  7.2003E-01
             1.0655E+00
 PARAMETER:  1.2055E-01  7.2584E-01 -7.9445E-01 -6.1973E-01  2.0677E-01  1.0458E-01 -2.8175E-01 -1.0707E-01  5.9252E-01 -2.2846E-01
             1.6340E-01
 GRADIENT:  -1.7074E+01  1.8347E+01 -9.2927E-01  1.2081E+01 -2.7439E-01 -1.5661E+00 -1.6654E+00  7.4358E-01  8.2011E-01 -2.0457E+00
            -7.5454E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1709.60686876334        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      972
 NPARAMETR:  1.0280E+00  1.8578E+00  4.0875E-01  4.8681E-01  1.1120E+00  1.0082E+00  6.8850E-01  8.1299E-01  1.6362E+00  7.3630E-01
             1.0671E+00
 PARAMETER:  1.2762E-01  7.1937E-01 -7.9465E-01 -6.1988E-01  2.0614E-01  1.0813E-01 -2.7324E-01 -1.0704E-01  5.9237E-01 -2.0612E-01
             1.6490E-01
 GRADIENT:  -1.4669E+00 -3.0620E+00 -3.3355E+00  7.7852E+00  2.7355E+00 -6.5583E-04  5.1400E-01  1.0692E+00  2.2170E+00  6.1682E-01
             1.2285E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1709.61619926589        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1156            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0287E+00  1.8579E+00  4.0876E-01  4.8680E-01  1.1096E+00  1.0082E+00  6.8737E-01  8.1197E-01  1.6362E+00  7.3036E-01
             1.0646E+00
 PARAMETER:  1.2826E-01  7.1947E-01 -7.9463E-01 -6.1990E-01  2.0401E-01  1.0814E-01 -2.7489E-01 -1.0830E-01  5.9235E-01 -2.1422E-01
             1.6261E-01
 GRADIENT:   5.1692E+01  7.9219E+01 -1.5463E+00  1.8856E+01  1.1879E+00  6.7509E+00  1.7763E+00  1.0045E+00  4.6609E+00  5.9670E-02
             1.3375E-01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1709.61726113891        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  1.0286E+00  1.8578E+00  4.0875E-01  4.8678E-01  1.1096E+00  1.0081E+00  6.8737E-01  8.1142E-01  1.6361E+00  7.3039E-01
             1.0646E+00
 PARAMETER:  1.2821E-01  7.1942E-01 -7.9465E-01 -6.1995E-01  2.0401E-01  1.0810E-01 -2.7488E-01 -1.0897E-01  5.9232E-01 -2.1418E-01
             1.6261E-01
 GRADIENT:  -7.9014E-02 -2.3575E-02  2.3666E+05  1.5167E+05  9.5973E-03 -1.4731E-02 -8.5209E-03 -1.7258E+06  1.5873E+05 -5.1707E-03
             7.7158E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1272
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6796E-04 -3.0989E-02 -1.4925E-02  2.3632E-02 -4.0720E-02
 SE:             2.9864E-02  2.4255E-02  6.8633E-03  2.3006E-02  1.9460E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8483E-01  2.0138E-01  2.9663E-02  3.0431E-01  3.6391E-02

 ETASHRINKSD(%)  1.0000E-10  1.8743E+01  7.7007E+01  2.2929E+01  3.4807E+01
 ETASHRINKVR(%)  1.0000E-10  3.3973E+01  9.4713E+01  4.0600E+01  5.7499E+01
 EBVSHRINKSD(%)  4.6506E-01  1.8866E+01  7.7847E+01  2.2078E+01  3.4775E+01
 EBVSHRINKVR(%)  9.2796E-01  3.4172E+01  9.5092E+01  3.9281E+01  5.7457E+01
 RELATIVEINF(%)  9.9006E+01  4.9061E+00  4.2012E-01  4.6519E+00  8.7559E+00
 EPSSHRINKSD(%)  4.3884E+01
 EPSSHRINKVR(%)  6.8510E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1709.6172611389147     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -974.46643457517655     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.53
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1709.617       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.86E+00  4.09E-01  4.87E-01  1.11E+00  1.01E+00  6.87E-01  8.11E-01  1.64E+00  7.30E-01  1.06E+00
 


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
+       -7.38E+00  4.70E+02
 
 TH 3
+        3.83E+03 -6.10E+03  4.46E+08
 
 TH 4
+        4.36E+03 -6.61E+03 -2.54E+04  5.16E+08
 
 TH 5
+       -5.54E+00 -2.77E+02  2.97E+03  3.70E+03  7.34E+02
 
 TH 6
+       -5.40E-02 -1.04E+00  7.61E+03  8.65E+03 -1.73E-01  1.92E+02
 
 TH 7
+        1.12E+00  8.91E+00  2.23E+03  2.35E+03 -1.20E+01  1.43E+00  2.00E+02
 
 TH 8
+       -1.62E+04 -3.98E+08 -8.59E+04 -1.76E+09 -9.70E+03 -3.19E+04 -7.85E+03  6.01E+09
 
 TH 9
+        1.24E+03 -2.12E+03  3.18E+04 -7.83E+03  1.20E+03  2.46E+03  7.75E+02 -2.98E+04  5.01E+07
 
 TH10
+       -6.54E-01 -1.13E+01  1.69E+03  1.80E+03 -8.08E+01 -1.46E+00  2.67E+01 -6.05E+03  5.88E+02  8.83E+01
 
 TH11
+       -8.08E+00 -1.63E+01 -4.50E+05 -4.85E+05 -1.88E+01  1.67E+00  1.03E+01  1.66E+06 -1.51E+05  2.24E+01  1.89E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.018
Stop Time:
Sat Sep 18 10:52:05 CDT 2021
