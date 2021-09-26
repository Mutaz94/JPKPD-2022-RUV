Sat Sep 25 07:49:54 CDT 2021
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
$DATA ../../../../data/spa/A1/dat3.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1126.71285051310        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4207E+00 -1.7454E+01 -1.5576E+01 -3.1372E+00  1.3309E+02  1.2134E+01 -3.0973E+01 -6.5951E+00 -4.4416E+01 -7.6088E+01
            -9.6959E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1440.99210570221        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0259E+00  9.6706E-01  1.0278E+00  1.0567E+00  9.0127E-01  9.3671E-01  1.0538E+00  9.6035E-01  1.0862E+00  1.0886E+00
             2.2453E+00
 PARAMETER:  1.2556E-01  6.6507E-02  1.2745E-01  1.5518E-01 -3.9519E-03  3.4620E-02  1.5237E-01  5.9538E-02  1.8270E-01  1.8485E-01
             9.0882E-01
 GRADIENT:   1.2533E+00  1.2825E+01  3.4517E+00  7.0852E+00 -1.0056E+01 -1.0529E+01 -9.4371E-01  3.2397E+00  2.7915E+00  1.9527E-01
            -5.9530E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1443.06574531900        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0262E+00  7.3125E-01  8.6989E-01  1.2311E+00  7.5350E-01  9.5090E-01  1.3581E+00  5.2625E-01  9.7100E-01  1.0116E+00
             2.2246E+00
 PARAMETER:  1.2582E-01 -2.1300E-01 -3.9389E-02  3.0790E-01 -1.8303E-01  4.9654E-02  4.0607E-01 -5.4197E-01  7.0572E-02  1.1152E-01
             8.9957E-01
 GRADIENT:  -1.8788E+00  2.8287E+01 -1.1142E+01  6.7810E+01  8.6403E+00 -5.2531E+00  3.2155E+00  1.6185E+00 -1.9715E+00  5.5444E+00
            -1.1189E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1447.98680326648        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0323E+00  5.0008E-01  5.3514E-01  1.2515E+00  4.8613E-01  9.7057E-01  1.3609E+00  8.1989E-02  9.1903E-01  7.4563E-01
             2.2609E+00
 PARAMETER:  1.3178E-01 -5.9299E-01 -5.2523E-01  3.2436E-01 -6.2128E-01  7.0123E-02  4.0815E-01 -2.4012E+00  1.5563E-02 -1.9352E-01
             9.1575E-01
 GRADIENT:   3.7036E+00  8.9756E+00 -6.1384E+00  1.3308E+01 -3.3901E-03 -3.5716E-01  2.8241E-01  9.1305E-02 -4.1761E+00  2.4103E+00
            -1.9740E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1450.95424406593        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0243E+00  2.9053E-01  6.1212E-01  1.3710E+00  4.8656E-01  9.6185E-01  1.7557E+00  2.3131E-02  8.9155E-01  7.6774E-01
             2.2835E+00
 PARAMETER:  1.2399E-01 -1.1360E+00 -3.9083E-01  4.1551E-01 -6.2040E-01  6.1101E-02  6.6288E-01 -3.6666E+00 -1.4788E-02 -1.6430E-01
             9.2571E-01
 GRADIENT:  -1.9855E-01  2.9542E+00  9.0078E-01  7.4924E+00 -1.9477E+00 -2.6627E-01 -7.3430E-01  5.7663E-03 -4.5625E-01 -1.6484E+00
            -9.0666E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1451.91624667088        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0179E+00  1.0075E-01  6.2390E-01  1.4566E+00  4.6332E-01  9.5409E-01  3.6485E+00  1.0000E-02  8.5171E-01  7.7097E-01
             2.2934E+00
 PARAMETER:  1.1774E-01 -2.1951E+00 -3.7176E-01  4.7613E-01 -6.6934E-01  5.3002E-02  1.3943E+00 -7.4359E+00 -6.0508E-02 -1.6010E-01
             9.3006E-01
 GRADIENT:   5.8642E-01  2.5181E-01  4.5103E-01 -2.5360E+00  7.6654E-01 -9.7973E-01  5.6833E-01  0.0000E+00  1.9719E-01 -1.8267E+00
            -4.0854E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1452.39171340338        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  1.0202E+00  7.0180E-02  7.1344E-01  1.5125E+00  5.0554E-01  9.5703E-01  4.4913E+00  1.0000E-02  8.3775E-01  8.2476E-01
             2.3038E+00
 PARAMETER:  1.2003E-01 -2.5567E+00 -2.3766E-01  5.1376E-01 -5.8212E-01  5.6080E-02  1.6021E+00 -8.5314E+00 -7.7039E-02 -9.2668E-02
             9.3455E-01
 GRADIENT:   6.0943E-01  6.1195E-01  4.5786E-01  5.3177E+00 -1.4932E+00  3.6385E-01  4.9202E-01  0.0000E+00 -1.3435E-01  2.1571E-01
             4.0390E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1452.52788587474        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0179E+00  2.1838E-02  7.2116E-01  1.5343E+00  5.0193E-01  9.5421E-01  6.5536E+00  1.0000E-02  8.3168E-01  8.2640E-01
             2.3047E+00
 PARAMETER:  1.1773E-01 -3.7241E+00 -2.2689E-01  5.2809E-01 -5.8930E-01  5.3126E-02  1.9800E+00 -1.3006E+01 -8.4306E-02 -9.0680E-02
             9.3494E-01
 GRADIENT:  -3.7579E-01  2.0955E-02  2.9663E-01 -9.4375E-01 -2.6214E-01 -2.3084E-01 -2.5551E-02  0.0000E+00  4.6477E-01 -1.2222E-01
             2.7400E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1452.54764659116        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0176E+00  1.0000E-02  7.2549E-01  1.5418E+00  5.0221E-01  9.5431E-01  8.6781E+00  1.0000E-02  8.2800E-01  8.2873E-01
             2.3056E+00
 PARAMETER:  1.1743E-01 -4.6886E+00 -2.2090E-01  5.3297E-01 -5.8873E-01  5.3234E-02  2.2608E+00 -1.6763E+01 -8.8737E-02 -8.7861E-02
             9.3532E-01
 GRADIENT:  -1.3580E-01  0.0000E+00  1.9178E-01  2.6715E-02 -2.9395E-01 -7.4298E-02 -1.9287E-02  0.0000E+00  5.1988E-02 -3.8001E-02
            -4.1553E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1452.54992267069        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1043
 NPARAMETR:  1.0176E+00  1.0000E-02  7.2591E-01  1.5420E+00  5.0251E-01  9.5448E-01  1.0378E+01  1.0000E-02  8.2756E-01  8.2914E-01
             2.3055E+00
 PARAMETER:  1.1749E-01 -5.3229E+00 -2.2032E-01  5.3306E-01 -5.8814E-01  5.3409E-02  2.4397E+00 -1.9258E+01 -8.9276E-02 -8.7370E-02
             9.3530E-01
 GRADIENT:   6.2390E-04  0.0000E+00 -1.1480E-03  1.8644E-03  1.6036E-03  3.2927E-04 -1.8462E-05  0.0000E+00  7.8307E-04  2.8372E-04
            -3.3320E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1043
 NO. OF SIG. DIGITS IN FINAL EST.:  4.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.1698E-04 -1.1536E-04  7.7940E-06 -9.9104E-03 -1.5281E-02
 SE:             2.9250E-02  1.5355E-03  1.8602E-04  2.7336E-02  2.1795E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9135E-01  9.4011E-01  9.6658E-01  7.1695E-01  4.8322E-01

 ETASHRINKSD(%)  2.0095E+00  9.4856E+01  9.9377E+01  8.4202E+00  2.6984E+01
 ETASHRINKVR(%)  3.9786E+00  9.9735E+01  9.9996E+01  1.6131E+01  4.6687E+01
 EBVSHRINKSD(%)  2.0902E+00  9.4996E+01  9.9358E+01  7.7803E+00  2.6030E+01
 EBVSHRINKVR(%)  4.1367E+00  9.9750E+01  9.9996E+01  1.4955E+01  4.5285E+01
 RELATIVEINF(%)  8.1437E+01  7.4496E-03  2.0411E-04  4.2021E+00  1.8029E+00
 EPSSHRINKSD(%)  3.3903E+01
 EPSSHRINKVR(%)  5.6312E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1452.5499226706893     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -717.39909610695111     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1452.550       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E-02  7.26E-01  1.54E+00  5.03E-01  9.54E-01  1.04E+01  1.00E-02  8.28E-01  8.29E-01  2.31E+00
 


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
+        1.13E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.38E+00  0.00E+00  8.80E+02
 
 TH 4
+       -3.84E+01  0.00E+00 -1.02E+02  5.82E+02
 
 TH 5
+        4.06E+01  0.00E+00 -1.64E+03 -1.27E+02  3.49E+03
 
 TH 6
+        6.94E+00  0.00E+00  1.08E+01 -1.03E+01 -4.71E+00  2.00E+02
 
 TH 7
+       -2.45E-02  0.00E+00  8.41E-03 -1.02E-02 -7.66E-03 -5.69E-02  1.07E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.48E+00  0.00E+00  1.98E+01 -6.83E+00 -9.63E+00  1.79E+01  3.86E-02  0.00E+00  2.18E+02
 
 TH10
+        1.25E+01  0.00E+00 -1.42E+01  1.13E+00 -5.94E+01  1.09E+01 -1.55E-02  0.00E+00 -1.10E+00  1.15E+02
 
 TH11
+       -1.29E+01  0.00E+00 -1.12E+01 -9.34E+00  5.49E+00  1.66E+00  5.57E-03  0.00E+00  1.13E+01  2.04E+01  5.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       17.995
Stop Time:
Sat Sep 25 07:50:13 CDT 2021
