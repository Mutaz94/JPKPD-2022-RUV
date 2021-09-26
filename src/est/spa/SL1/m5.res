Sat Sep 25 10:17:12 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat5.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1659.14827932835        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6379E+02 -3.6401E+01 -1.9486E+01 -2.4064E+01 -2.9235E+00  1.4745E+01 -4.0418E+00  1.2252E+01  1.9721E+01  3.8039E+00
             3.4158E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.69076190278        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.5432E-01  1.1321E+00  1.1671E+00  9.3999E-01  1.1588E+00  9.3831E-01  1.1682E+00  8.1976E-01  8.0270E-01  1.1498E+00
             9.3670E-01
 PARAMETER:  5.3248E-02  2.2409E-01  2.5449E-01  3.8116E-02  2.4735E-01  3.6322E-02  2.5547E-01 -9.8742E-02 -1.1978E-01  2.3961E-01
             3.4612E-02
 GRADIENT:   7.2179E+01  1.0273E+01  1.1734E+01 -1.1227E+01 -4.4011E+00 -2.0909E+00  1.8408E+00  1.3905E+00 -6.0673E+00 -7.9080E+00
             1.5077E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1671.15395955396        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.5043E-01  1.0174E+00  1.1333E+00  1.0228E+00  1.1021E+00  9.3507E-01  1.2246E+00  4.4509E-01  8.8654E-01  1.1726E+00
             9.2354E-01
 PARAMETER:  4.9157E-02  1.1729E-01  2.2511E-01  1.2252E-01  1.9723E-01  3.2861E-02  3.0259E-01 -7.0948E-01 -2.0425E-02  2.5919E-01
             2.0456E-02
 GRADIENT:   6.4256E+01  9.3379E+00 -4.4631E-01  1.7940E+01  7.8520E+00 -3.5437E+00  5.7862E+00  3.6396E-01  8.7898E+00 -1.8223E-01
            -1.5842E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1671.63281385112        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.3732E-01  1.0148E+00  9.8631E-01  1.0131E+00  1.0259E+00  9.4584E-01  1.2192E+00  3.3854E-01  8.3521E-01  1.0930E+00
             9.2123E-01
 PARAMETER:  3.5268E-02  1.1467E-01  8.6218E-02  1.1298E-01  1.2559E-01  4.4314E-02  2.9816E-01 -9.8312E-01 -8.0073E-02  1.8891E-01
             1.7960E-02
 GRADIENT:   2.6687E+01  1.6103E+00 -9.2234E+00  1.4779E+01  9.7886E+00  1.0151E+00  6.8934E-01  6.3798E-01  2.3061E+00  1.8911E+00
             1.0827E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1671.63362049430        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.3597E-01  1.0151E+00  9.8671E-01  1.0120E+00  1.0256E+00  9.4569E-01  1.2186E+00  3.3040E-01  8.3430E-01  1.0924E+00
             9.2125E-01
 PARAMETER:  3.3832E-02  1.1496E-01  8.6618E-02  1.1197E-01  1.2526E-01  4.4161E-02  2.9771E-01 -1.0075E+00 -8.1167E-02  1.8835E-01
             1.7975E-02
 GRADIENT:   2.3177E+01  1.3107E+00 -8.1294E+00  1.2841E+01  8.5742E+00  9.2168E-01  5.9630E-01  5.8263E-01  2.0481E+00  1.6846E+00
             9.6813E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1671.63386883778        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.3535E-01  1.0153E+00  9.8663E-01  1.0115E+00  1.0253E+00  9.4561E-01  1.2184E+00  3.2441E-01  8.3378E-01  1.0920E+00
             9.2126E-01
 PARAMETER:  3.3166E-02  1.1515E-01  8.6543E-02  1.1145E-01  1.2500E-01  4.4070E-02  2.9756E-01 -1.0258E+00 -8.1780E-02  1.8804E-01
             1.7983E-02
 GRADIENT:   2.1541E+01  1.1976E+00 -7.5836E+00  1.1936E+01  7.9858E+00  8.6598E-01  5.5371E-01  5.4836E-01  1.9140E+00  1.5758E+00
             9.0639E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1671.63393330276        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  9.3504E-01  1.0154E+00  9.8653E-01  1.0112E+00  1.0252E+00  9.4556E-01  1.2183E+00  3.2087E-01  8.3351E-01  1.0918E+00
             9.2126E-01
 PARAMETER:  3.2837E-02  1.1526E-01  8.6436E-02  1.1118E-01  1.2484E-01  4.4022E-02  2.9750E-01 -1.0367E+00 -8.2105E-02  1.8786E-01
             1.7988E-02
 GRADIENT:   2.0733E+01  1.1475E+00 -7.3062E+00  1.1489E+01  7.6903E+00  8.3575E-01  5.3280E-01  5.2956E-01  1.8448E+00  1.5193E+00
             8.7410E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1671.94406065612        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      612
 NPARAMETR:  9.4493E-01  9.5002E-01  1.0368E+00  1.0556E+00  1.0154E+00  9.5166E-01  1.2857E+00  1.9470E-01  8.1039E-01  1.1071E+00
             9.2490E-01
 PARAMETER:  4.3353E-02  4.8729E-02  1.3614E-01  1.5408E-01  1.1525E-01  5.0450E-02  3.5129E-01 -1.5363E+00 -1.1024E-01  2.0174E-01
             2.1928E-02
 GRADIENT:   3.7320E+00  1.9973E+00 -1.0684E-01  3.2293E+00 -1.0116E+00  1.7404E-01 -1.3149E-01  8.4873E-02  3.8340E-02  6.6331E-01
             5.6770E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1671.97689240749        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  9.4324E-01  9.1390E-01  1.0452E+00  1.0753E+00  1.0048E+00  9.5140E-01  1.3238E+00  8.5152E-02  8.0114E-01  1.1050E+00
             9.2467E-01
 PARAMETER:  4.1570E-02  9.9681E-03  1.4417E-01  1.7263E-01  1.0483E-01  5.0181E-02  3.8051E-01 -2.3633E+00 -1.2172E-01  1.9988E-01
             2.1685E-02
 GRADIENT:   7.2774E-02 -6.9896E-01 -2.1597E-01 -7.6020E-01  5.5520E-02  1.4315E-01  1.2194E-01  1.2503E-02  9.9118E-02  3.7116E-01
             2.0285E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1671.98895957663        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      965
 NPARAMETR:  9.4358E-01  9.5308E-01  1.0219E+00  1.0508E+00  1.0106E+00  9.5129E-01  1.2872E+00  1.4031E-02  8.0936E-01  1.0997E+00
             9.2388E-01
 PARAMETER:  4.1921E-02  5.1946E-02  1.2163E-01  1.4954E-01  1.1056E-01  5.0062E-02  3.5248E-01 -4.1665E+00 -1.1151E-01  1.9504E-01
             2.0825E-02
 GRADIENT:  -3.0351E-04  7.9154E-02  5.0985E-03  1.0250E-01 -8.1163E-03 -2.2263E-02 -6.4397E-03  3.6351E-04 -1.1050E-02 -1.2562E-02
            -3.5864E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1671.98907035450        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  9.4365E-01  9.5105E-01  1.0229E+00  1.0520E+00  1.0102E+00  9.5131E-01  1.2890E+00  1.0000E-02  8.0890E-01  1.0998E+00
             9.2385E-01
 PARAMETER:  4.1904E-02  4.9814E-02  1.2266E-01  1.5071E-01  1.1019E-01  5.0106E-02  3.5391E-01 -4.7865E+00 -1.1201E-01  1.9521E-01
             2.0858E-02
 GRADIENT:  -1.2057E-02  5.3576E-05  1.0662E-04  6.7320E-04 -2.7728E-04  3.3279E-04  6.1283E-05  0.0000E+00  5.6642E-04  1.1355E-03
             1.3826E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1134
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.4615E-04 -3.7975E-04 -3.7832E-04 -6.7752E-03 -2.0082E-02
 SE:             2.9823E-02  2.1158E-02  1.4974E-04  2.2737E-02  2.4021E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9341E-01  9.8568E-01  1.1521E-02  7.6572E-01  4.0316E-01

 ETASHRINKSD(%)  9.0775E-02  2.9117E+01  9.9498E+01  2.3828E+01  1.9526E+01
 ETASHRINKVR(%)  1.8147E-01  4.9756E+01  9.9997E+01  4.1978E+01  3.5239E+01
 EBVSHRINKSD(%)  4.0073E-01  2.8924E+01  9.9533E+01  2.4381E+01  1.6398E+01
 EBVSHRINKVR(%)  7.9986E-01  4.9482E+01  9.9998E+01  4.2818E+01  3.0107E+01
 RELATIVEINF(%)  9.8491E+01  2.1445E+00  2.5187E-04  2.5298E+00  9.0632E+00
 EPSSHRINKSD(%)  4.2830E+01
 EPSSHRINKVR(%)  6.7316E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1671.9890703545029     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -936.83824379076475     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1671.989       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.44E-01  9.51E-01  1.02E+00  1.05E+00  1.01E+00  9.51E-01  1.29E+00  1.00E-02  8.09E-01  1.10E+00  9.24E-01
 


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
+        1.37E+03
 
 TH 2
+       -9.70E+00  3.93E+02
 
 TH 3
+        1.41E+01  9.73E+01  2.49E+02
 
 TH 4
+       -6.99E+00  4.27E+02 -1.63E+02  8.81E+02
 
 TH 5
+       -4.21E+00 -1.86E+02 -3.32E+02  1.99E+02  6.30E+02
 
 TH 6
+       -2.49E+00  1.78E-01  4.06E+00 -2.02E+00 -3.69E-02  2.18E+02
 
 TH 7
+        1.81E+00  2.57E+01  6.95E+00 -9.78E+00 -7.96E+00  2.26E-01  3.58E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.03E+00 -1.27E+01 -2.79E+01  3.96E+00  1.77E+01  1.95E+00  2.86E+01  0.00E+00  1.12E+02
 
 TH10
+       -8.98E-01 -7.54E+00 -3.20E+01 -7.36E+00 -4.94E+01 -2.33E-02  1.53E+00  0.00E+00  1.06E+01  7.94E+01
 
 TH11
+       -8.44E+00 -1.75E+01 -4.48E+01 -7.29E-01  1.20E+01 -7.40E-01  5.53E+00  0.00E+00  1.12E+01  2.21E+01  2.63E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.471
Stop Time:
Sat Sep 25 10:17:41 CDT 2021
