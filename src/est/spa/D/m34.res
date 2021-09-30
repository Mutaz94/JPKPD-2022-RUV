Wed Sep 29 19:58:16 CDT 2021
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
$DATA ../../../../data/spa/D/dat34.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13715.7004183556        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7801E+02  1.2935E+02 -3.8992E+01 -4.4911E+00  9.3466E+01 -1.3868E+03 -5.3372E+02 -5.7546E+01 -8.4275E+02 -4.2466E+02
            -2.7214E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -654.602457517617        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4069E+00  1.2957E+00  1.0358E+00  1.6113E+00  1.0124E+00  1.9063E+00  1.1090E+00  1.0150E+00  1.0859E+00  1.0424E+00
             1.4567E+01
 PARAMETER:  4.4136E-01  3.5907E-01  1.3513E-01  5.7705E-01  1.1230E-01  7.4515E-01  2.0344E-01  1.1490E-01  1.8245E-01  1.4149E-01
             2.7788E+00
 GRADIENT:   2.5267E+01  1.1019E+01  1.0601E+00  1.3569E+01 -1.1248E+01  4.4356E+01  2.5813E+00  3.8510E+00  1.4261E+01  5.7911E+00
             1.6717E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -667.769805872413        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.3446E+00  1.0647E+00  1.7321E+00  1.8990E+00  2.7977E+00  1.7238E+00  2.7432E+00  8.7224E-01  7.6843E-01  3.2935E+00
             1.3650E+01
 PARAMETER:  3.9611E-01  1.6270E-01  6.4935E-01  7.4134E-01  1.1288E+00  6.4454E-01  1.1091E+00 -3.6693E-02 -1.6341E-01  1.2920E+00
             2.7137E+00
 GRADIENT:   1.9070E+01  2.9979E+01  3.9009E+00  5.8696E+01 -1.0522E+00  1.1067E+01  1.0952E+01  4.7321E-01  4.6858E+00  8.2355E-01
             1.2283E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -693.479233614722        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0743E+00  5.2016E-01  2.9784E+00  1.6070E+00  4.8107E+00  1.4924E+00  1.7577E+00  4.3701E-01  5.5695E-01  8.4674E+00
             1.1144E+01
 PARAMETER:  1.7166E-01 -5.5363E-01  1.1914E+00  5.7435E-01  1.6708E+00  5.0037E-01  6.6402E-01 -7.2779E-01 -4.8529E-01  2.2362E+00
             2.5109E+00
 GRADIENT:  -3.4259E+01 -5.4754E+00  1.0190E+01 -4.2231E+01 -3.7066E+00  9.1806E+00  3.6996E+00  1.0725E-02  7.6280E+00  3.6597E-02
             7.1380E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -718.609401249510        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.9750E-01  2.8673E-01  7.1575E-01  1.4976E+00  9.3205E+00  1.4353E+00  9.9630E-01  1.0000E-02  1.7642E-01  7.8518E+00
             9.6615E+00
 PARAMETER:  9.7493E-02 -1.1492E+00 -2.3442E-01  5.0384E-01  2.3322E+00  4.6140E-01  9.6293E-02 -6.7234E+00 -1.6349E+00  2.1607E+00
             2.3681E+00
 GRADIENT:  -1.4034E+01  6.0124E+00  2.1911E+01  4.2090E+00 -1.4163E+00 -1.5953E+01  2.0377E-01  0.0000E+00  1.2855E+00  9.2571E+00
            -7.1574E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -766.869465798337        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  8.4939E-01  8.8399E-02  1.3032E-01  1.0059E+00  1.2044E+01  1.6245E+00  6.8856E-01  1.0000E-02  1.0000E-02  1.6737E+00
             1.0808E+01
 PARAMETER: -6.3234E-02 -2.3259E+00 -1.9377E+00  1.0592E-01  2.5886E+00  5.8518E-01 -2.7315E-01 -1.5004E+01 -5.8383E+00  6.1506E-01
             2.4803E+00
 GRADIENT:   3.1363E+01  5.4614E+00 -3.5976E+01  9.9844E+01  4.7529E-02 -3.9540E+00  1.7604E+00  0.0000E+00  0.0000E+00  1.8807E-01
             2.6875E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -778.586574523070        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      487
 NPARAMETR:  6.6060E-01  3.4537E-02  7.0211E-02  6.5353E-01  1.8375E+01  1.5724E+00  2.8907E-01  1.0000E-02  1.0000E-02  9.1950E-01
             9.8350E+00
 PARAMETER: -3.1461E-01 -3.2657E+00 -2.5562E+00 -3.2537E-01  3.0110E+00  5.5260E-01 -1.1411E+00 -1.9571E+01 -8.4290E+00  1.6074E-02
             2.3859E+00
 GRADIENT:   6.9926E+00  3.7982E-01 -3.0375E+01  5.1703E+01  9.3191E-02 -9.7488E+00  5.5623E-03  0.0000E+00  0.0000E+00  1.8394E-05
            -3.5686E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -780.277389867942        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  6.1327E-01  2.3315E-02  5.6866E-02  5.5945E-01  1.9357E+01  1.5991E+00  2.3153E-01  1.0000E-02  1.0000E-02  6.4221E-01
             1.0084E+01
 PARAMETER: -3.8896E-01 -3.6587E+00 -2.7671E+00 -4.8080E-01  3.0630E+00  5.6946E-01 -1.3630E+00 -2.1305E+01 -9.6888E+00 -3.4283E-01
             2.4110E+00
 GRADIENT:   4.1004E-01  2.6757E-01 -8.1052E-01 -4.5304E-01 -4.5790E-03  2.6114E-01  2.8480E-04  0.0000E+00  0.0000E+00 -1.4651E-05
            -2.3360E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -780.298541653672        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  6.1308E-01  1.9920E-02  5.6870E-02  5.6014E-01  1.3075E+01  1.5984E+00  1.1013E-01  1.0000E-02  1.0000E-02  6.5948E-01
             1.0108E+01
 PARAMETER: -3.8927E-01 -3.8160E+00 -2.7670E+00 -4.7956E-01  2.6707E+00  5.6902E-01 -2.1061E+00 -2.1305E+01 -9.6888E+00 -3.1631E-01
             2.4133E+00
 GRADIENT:  -2.4783E-01  2.5995E-02 -7.5340E-01  3.9550E-02  2.5784E-03 -2.0802E-01  1.4451E-05  0.0000E+00  0.0000E+00 -3.8614E-05
            -2.1532E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -780.299317569155        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1033
 NPARAMETR:  6.1332E-01  1.9163E-02  5.6873E-02  5.6028E-01  1.2123E+01  1.6000E+00  8.9258E-02  1.0000E-02  1.0000E-02  6.9010E-01
             1.0108E+01
 PARAMETER: -3.8887E-01 -3.8548E+00 -2.7669E+00 -4.7933E-01  2.5951E+00  5.6999E-01 -2.3162E+00 -2.1305E+01 -9.6888E+00 -2.7091E-01
             2.4133E+00
 GRADIENT:   5.3283E-02 -1.7175E-03 -7.0337E-01 -1.7015E-01  1.4252E-03  4.9641E-02  7.1065E-06  0.0000E+00  0.0000E+00 -3.6348E-05
            -2.4542E-01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -780.299337893643        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:     1140
 NPARAMETR:  6.1348E-01  1.9172E-02  5.6873E-02  5.6032E-01  1.2036E+01  1.6002E+00  7.0855E-02  1.0000E-02  1.0000E-02  6.9616E-01
             1.0110E+01
 PARAMETER: -3.8872E-01 -3.8525E+00 -2.7670E+00 -4.7929E-01  2.5877E+00  5.7011E-01 -2.5219E+00 -2.1305E+01 -9.6888E+00 -2.6343E-01
             2.4135E+00
 GRADIENT:  -6.2861E-02  1.3734E-03 -2.4819E-02 -5.7873E-02 -2.8389E-05 -4.9633E-03  3.9078E-06  0.0000E+00  0.0000E+00 -3.1868E-05
            -2.7745E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1140
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1452E-03  8.4467E-06  1.3202E-04 -2.9739E-04  3.1829E-05
 SE:             2.8845E-02  9.8337E-06  2.5073E-04  3.4936E-04  8.1789E-05
 N:                     100         100         100         100         100

 P VAL.:         9.1317E-01  3.9036E-01  5.9852E-01  3.9463E-01  6.9716E-01

 ETASHRINKSD(%)  3.3649E+00  9.9967E+01  9.9160E+01  9.8830E+01  9.9726E+01
 ETASHRINKVR(%)  6.6165E+00  1.0000E+02  9.9993E+01  9.9986E+01  9.9999E+01
 EBVSHRINKSD(%)  3.5050E+00  9.9962E+01  9.9103E+01  9.8744E+01  9.9735E+01
 EBVSHRINKVR(%)  6.8871E+00  1.0000E+02  9.9992E+01  9.9984E+01  9.9999E+01
 RELATIVEINF(%)  3.7904E+00  1.1868E-06  7.3899E-05  8.3223E-05  2.5990E-05
 EPSSHRINKSD(%)  7.1538E+00
 EPSSHRINKVR(%)  1.3796E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -780.29933789364350     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -45.148511329905318     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -780.299       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.13E-01  1.92E-02  5.69E-02  5.60E-01  1.20E+01  1.60E+00  7.27E-02  1.00E-02  1.00E-02  6.95E-01  1.01E+01
 


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
+        1.06E+03
 
 TH 2
+       -2.68E+02  1.46E+03
 
 TH 3
+       -2.49E+03 -6.15E+02  1.60E+05
 
 TH 4
+       -1.76E+02  1.77E+02 -2.29E+04  3.63E+03
 
 TH 5
+        2.12E-01 -5.57E-01 -2.95E+00  3.76E-01  7.66E-04
 
 TH 6
+       -6.94E-01  1.15E+01  3.98E+02 -6.69E+01  5.67E-03  6.57E+01
 
 TH 7
+        1.03E-03  1.56E-02 -2.36E-03  2.45E-03 -2.33E-05  9.35E-04  9.22E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        4.45E-05  2.23E-02  1.47E-02 -3.14E-03 -4.12E-05  9.51E-04 -4.03E-04  0.00E+00  0.00E+00  1.10E-02
 
 TH11
+       -1.57E+01  2.81E+00  1.75E+02 -2.03E+01 -5.98E-03  4.13E-01 -1.32E-05  0.00E+00  0.00E+00  8.74E-05  4.41E+00
 
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
 #CPUT: Total CPU Time in Seconds,       21.655
Stop Time:
Wed Sep 29 19:58:40 CDT 2021
