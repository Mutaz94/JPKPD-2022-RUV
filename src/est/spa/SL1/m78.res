Sat Sep 25 10:44:32 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1649.39119171851        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.7158E+01 -1.0238E+02 -5.2465E+01 -9.2802E+01  3.0568E+01 -2.9331E+00 -1.4979E+01  1.6054E+01 -7.4831E+00  5.7876E+00
            -6.6561E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1666.97175512141        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7726E-01  1.1627E+00  1.4830E+00  9.8797E-01  1.2268E+00  1.0017E+00  1.1096E+00  7.8971E-01  1.0603E+00  1.1870E+00
             1.0368E+00
 PARAMETER:  7.7002E-02  2.5078E-01  4.9408E-01  8.7894E-02  3.0438E-01  1.0170E-01  2.0397E-01 -1.3609E-01  1.5855E-01  2.7139E-01
             1.3609E-01
 GRADIENT:   4.6252E+01  1.3225E+01  4.3785E+00 -1.5232E+00 -1.1316E+01 -9.9980E-01  8.8456E+00  2.3935E+00  1.9688E+00 -3.0992E+00
             1.7904E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1668.09485657490        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.7358E-01  1.0967E+00  1.6303E+00  1.0406E+00  1.2566E+00  1.0106E+00  9.4669E-01  6.8566E-01  1.1154E+00  1.2512E+00
             1.0459E+00
 PARAMETER:  7.3226E-02  1.9232E-01  5.8874E-01  1.3982E-01  3.2842E-01  1.1052E-01  4.5214E-02 -2.7738E-01  2.0921E-01  3.2409E-01
             1.4487E-01
 GRADIENT:   3.8184E+01  1.1182E+01  2.2737E+00  1.6240E+01 -2.8462E-01  2.8319E+00  4.0991E+00  6.6061E-01  6.2302E+00 -3.4313E-01
             4.8239E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1668.52966164078        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.6927E-01  1.1650E+00  1.4951E+00  9.8743E-01  1.2529E+00  9.9341E-01  8.4654E-01  5.5119E-01  1.1649E+00  1.2536E+00
             1.0351E+00
 PARAMETER:  6.8785E-02  2.5272E-01  5.0219E-01  8.7354E-02  3.2550E-01  9.3386E-02 -6.6596E-02 -4.9567E-01  2.5260E-01  3.2603E-01
             1.3445E-01
 GRADIENT:   2.8049E+01  8.4802E+00 -2.9553E+00  8.8785E+00  4.4932E+00 -4.4397E+00  1.7985E+00  7.5734E-01  3.7690E+00  7.1248E-01
             1.2508E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1669.16683858427        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      392
 NPARAMETR:  9.7641E-01  1.3469E+00  1.4827E+00  8.6681E-01  1.3193E+00  1.0183E+00  5.8746E-01  4.3192E-01  1.3667E+00  1.3035E+00
             1.0389E+00
 PARAMETER:  7.6122E-02  3.9777E-01  4.9386E-01 -4.2931E-02  3.7710E-01  1.1815E-01 -4.3195E-01 -7.3951E-01  4.1240E-01  3.6503E-01
             1.3814E-01
 GRADIENT:   4.2164E+00  1.1711E+01  1.7255E+00  1.1397E+01 -1.1308E+00  1.7058E+00 -5.7504E-01 -2.7148E-01 -1.0077E+00 -4.5841E-01
            -1.2344E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1669.38043792690        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      567
 NPARAMETR:  9.7571E-01  1.4831E+00  1.4111E+00  7.6743E-01  1.3601E+00  1.0174E+00  5.7767E-01  5.1691E-01  1.4955E+00  1.3190E+00
             1.0395E+00
 PARAMETER:  7.5414E-02  4.9412E-01  4.4437E-01 -1.6471E-01  4.0757E-01  1.1730E-01 -4.4875E-01 -5.5988E-01  5.0244E-01  3.7685E-01
             1.3876E-01
 GRADIENT:   1.6228E+00  3.8574E+00 -7.3854E-01  3.3962E+00 -5.7537E-01  1.2559E+00 -8.6828E-01  3.3358E-01 -5.9599E-01  3.0712E-03
            -9.3612E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1669.79182606114        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      744
 NPARAMETR:  9.7468E-01  1.6997E+00  1.1735E+00  6.2648E-01  1.4067E+00  1.0119E+00  6.7459E-01  2.0378E-01  1.6180E+00  1.3099E+00
             1.0434E+00
 PARAMETER:  7.4357E-02  6.3043E-01  2.6001E-01 -3.6764E-01  4.4127E-01  1.1186E-01 -2.9365E-01 -1.4907E+00  5.8122E-01  3.6997E-01
             1.4244E-01
 GRADIENT:  -2.6282E+00  2.3884E+00 -1.3974E+00  1.6906E+00  2.6199E+00 -1.1535E+00 -1.5465E+00  1.4450E-01 -1.8408E+00 -5.8529E-01
             1.0336E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1669.87224282152        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      919
 NPARAMETR:  9.7599E-01  1.7608E+00  1.1377E+00  5.8235E-01  1.4251E+00  1.0149E+00  6.8079E-01  1.0508E-01  1.7061E+00  1.3208E+00
             1.0420E+00
 PARAMETER:  7.5696E-02  6.6576E-01  2.2901E-01 -4.4069E-01  4.5422E-01  1.1480E-01 -2.8450E-01 -2.1530E+00  6.3419E-01  3.7822E-01
             1.4119E-01
 GRADIENT:   3.4986E-02 -1.4252E+00 -1.2504E-01 -7.6569E-01 -4.4070E-02 -3.3742E-02  1.4480E-01  4.0670E-02  2.1064E-01  3.6910E-03
             1.5279E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1669.89553328290        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1096
 NPARAMETR:  9.7590E-01  1.7319E+00  1.1503E+00  6.0363E-01  1.4139E+00  1.0150E+00  6.8525E-01  2.0144E-02  1.6618E+00  1.3155E+00
             1.0409E+00
 PARAMETER:  7.5602E-02  6.4921E-01  2.4002E-01 -4.0479E-01  4.4633E-01  1.1487E-01 -2.7797E-01 -3.8049E+00  6.0787E-01  3.7423E-01
             1.4011E-01
 GRADIENT:  -1.9106E-02  9.3453E-01  5.9264E-02  4.3738E-01 -5.6586E-02  1.5929E-03 -1.0574E-01  1.4046E-03 -1.4861E-01 -7.9158E-03
            -4.5074E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1669.89649519856        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  9.7592E-01  1.7382E+00  1.1464E+00  5.9868E-01  1.4161E+00  1.0150E+00  6.8499E-01  1.0000E-02  1.6711E+00  1.3164E+00
             1.0411E+00
 PARAMETER:  7.5622E-02  6.5287E-01  2.3659E-01 -4.1303E-01  4.4792E-01  1.1491E-01 -2.7836E-01 -4.6989E+00  6.1347E-01  3.7491E-01
             1.4028E-01
 GRADIENT:  -4.4641E-03 -5.0925E-02  1.4189E-02 -5.6995E-02 -5.5463E-02  1.5091E-02 -2.3630E-03  0.0000E+00 -3.9740E-02  2.9762E-03
             3.6970E-03

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1669.89649519856        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1301
 NPARAMETR:  9.7594E-01  1.7382E+00  1.1461E+00  5.9885E-01  1.4163E+00  1.0149E+00  6.8503E-01  1.0000E-02  1.6721E+00  1.3164E+00
             1.0411E+00
 PARAMETER:  7.5622E-02  6.5287E-01  2.3659E-01 -4.1303E-01  4.4792E-01  1.1491E-01 -2.7836E-01 -4.6989E+00  6.1347E-01  3.7491E-01
             1.4028E-01
 GRADIENT:  -1.9948E-02  1.8990E-02  1.4025E-02 -5.6180E-02 -5.4231E-02  1.4243E-02 -1.7549E-03  0.0000E+00 -4.1248E-02  2.5757E-03
             4.1418E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1301
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.5521E-04 -3.8956E-02 -1.3854E-04  1.9892E-02 -4.1109E-02
 SE:             2.9766E-02  1.8097E-02  6.9500E-05  2.3800E-02  2.3717E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9316E-01  3.1351E-02  4.6217E-02  4.0327E-01  8.3039E-02

 ETASHRINKSD(%)  2.7864E-01  3.9371E+01  9.9767E+01  2.0268E+01  2.0545E+01
 ETASHRINKVR(%)  5.5651E-01  6.3241E+01  9.9999E+01  3.6428E+01  3.6870E+01
 EBVSHRINKSD(%)  4.7564E-01  3.6803E+01  9.9761E+01  2.1846E+01  1.7831E+01
 EBVSHRINKVR(%)  9.4901E-01  6.0062E+01  9.9999E+01  3.8919E+01  3.2483E+01
 RELATIVEINF(%)  9.8969E+01  3.0646E+00  2.2963E-04  5.2024E+00  2.7949E+01
 EPSSHRINKSD(%)  4.1370E+01
 EPSSHRINKVR(%)  6.5625E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1669.8964951985554     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -934.74566863481721     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1669.896       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.74E+00  1.15E+00  5.99E-01  1.42E+00  1.02E+00  6.85E-01  1.00E-02  1.67E+00  1.32E+00  1.04E+00
 


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
+        1.12E+03
 
 TH 2
+       -7.94E+00  3.32E+02
 
 TH 3
+        3.94E+00  4.08E+01  5.74E+01
 
 TH 4
+       -7.99E+00  3.83E+02 -2.98E+01  6.86E+02
 
 TH 5
+       -2.04E+00 -8.22E+01 -5.48E+01  3.75E+01  2.19E+02
 
 TH 6
+       -4.51E+00 -1.84E+00 -2.62E-01 -2.43E+00 -2.11E+00  1.86E+02
 
 TH 7
+        2.08E-01 -2.88E+01  1.55E+01 -2.08E+01 -1.74E+01 -8.32E-01  7.02E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.68E+00 -1.59E+01 -6.57E+00  3.76E+01  3.01E+00 -1.92E-01  2.88E+01  0.00E+00  3.12E+01
 
 TH10
+        7.05E-02 -3.99E+00 -7.77E+00 -2.63E+00 -3.38E+01  2.82E-01  5.85E+00  0.00E+00  7.92E-01  5.65E+01
 
 TH11
+       -1.34E+01 -2.30E+01 -1.75E+01 -1.29E+00  4.89E+00 -1.57E+00  8.35E+00  0.00E+00  3.07E+00  1.73E+01  2.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.572
Stop Time:
Sat Sep 25 10:44:58 CDT 2021
