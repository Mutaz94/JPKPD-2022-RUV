Sun Oct 24 03:57:55 CDT 2021
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
$DATA ../../../../data/SD4/TD2/dat6.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1661.37343991064        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3867E+02 -5.6186E+01 -1.0242E+02  7.4542E+01  1.9406E+02  5.4300E+01 -5.9974E+00  9.9354E+00  1.8430E+00 -1.1529E+01
             8.1497E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1673.86542610150        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  1.0291E+00  1.1095E+00  1.2157E+00  9.7318E-01  9.2264E-01  9.9752E-01  1.0295E+00  9.3541E-01  1.0527E+00  9.7220E-01
             9.6422E-01
 PARAMETER:  1.2865E-01  2.0394E-01  2.9531E-01  7.2811E-02  1.9488E-02  9.7519E-02  1.2908E-01  3.3228E-02  1.5139E-01  7.1809E-02
             6.3567E-02
 GRADIENT:   7.1249E+01  4.5009E+01  3.6370E+01  4.4411E+00 -9.0258E+01  1.1029E+00  4.1589E+00  8.7255E-01 -2.7235E+00 -4.7051E+00
            -1.2903E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1675.66818497648        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.0177E+00  1.0982E+00  1.6207E+00  9.7513E-01  1.0424E+00  9.7339E-01  6.5061E-01  8.5578E-01  1.2121E+00  1.1569E+00
             1.0230E+00
 PARAMETER:  1.1756E-01  1.9369E-01  5.8288E-01  7.4814E-02  1.4157E-01  7.3028E-02 -3.2984E-01 -5.5744E-02  2.9238E-01  2.4576E-01
             1.2278E-01
 GRADIENT:   4.9515E+01  2.3235E+01  3.1062E+01  1.0070E+01 -4.5194E+01 -6.6106E+00  3.6958E+00 -6.8135E+00  6.4866E+00  5.1307E-01
             7.7438E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1679.60166509513        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.9699E-01  1.2579E+00  1.2605E+00  8.4920E-01  1.0449E+00  9.8695E-01  5.2334E-01  9.0854E-01  1.3440E+00  1.1277E+00
             9.8573E-01
 PARAMETER:  9.6984E-02  3.2945E-01  3.3152E-01 -6.3465E-02  1.4397E-01  8.6863E-02 -5.4752E-01  4.0821E-03  3.9563E-01  2.2022E-01
             8.5628E-02
 GRADIENT:  -1.1829E+00  1.4812E+00  1.8596E+00 -8.0104E-01 -9.9864E+00 -9.1417E-01  4.9789E-01  6.3292E-01 -7.6064E-01  3.4162E-01
            -1.8184E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1679.72977104022        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.9859E-01  1.3324E+00  1.2545E+00  7.9552E-01  1.0785E+00  9.8979E-01  4.5224E-01  9.0717E-01  1.4552E+00  1.1557E+00
             9.9095E-01
 PARAMETER:  9.8587E-02  3.8702E-01  3.2673E-01 -1.2876E-01  1.7553E-01  8.9735E-02 -6.9355E-01  2.5772E-03  4.7516E-01  2.4470E-01
             9.0914E-02
 GRADIENT:   1.9135E+00 -3.0430E+00  4.4047E-01 -5.8465E-01 -6.8674E-01  2.0692E-01  1.4123E-01 -5.8649E-02  1.3574E+00  4.4799E-01
             1.1112E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1679.73325406023        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.9728E-01  1.3360E+00  1.2486E+00  7.9458E-01  1.0790E+00  9.8896E-01  4.5960E-01  9.0956E-01  1.4473E+00  1.1523E+00
             9.9096E-01
 PARAMETER:  9.7275E-02  3.8965E-01  3.2199E-01 -1.2995E-01  1.7602E-01  8.8900E-02 -6.7739E-01  5.2079E-03  4.6969E-01  2.4172E-01
             9.0922E-02
 GRADIENT:  -1.1446E+00 -1.1111E+00 -1.9225E-02  4.1229E-01  3.9658E-01 -1.3191E-01 -6.9598E-02  4.7756E-02 -1.0795E-01  1.5105E-02
             1.1911E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1679.73659654066        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1037
 NPARAMETR:  9.9835E-01  1.3360E+00  1.2466E+00  7.9409E-01  1.0783E+00  9.8943E-01  4.6483E-01  9.0263E-01  1.4471E+00  1.1518E+00
             9.9089E-01
 PARAMETER:  9.8354E-02  3.8967E-01  3.2046E-01 -1.3056E-01  1.7539E-01  8.9378E-02 -6.6609E-01 -2.4451E-03  4.6958E-01  2.4132E-01
             9.0845E-02
 GRADIENT:   1.3024E+00 -2.1162E+00  2.6358E-01 -3.2528E-01  2.1010E-02  5.7198E-02  5.9639E-02 -1.2163E-02  2.7589E-01  4.2480E-02
            -9.6015E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1679.73715113634        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1205
 NPARAMETR:  9.9859E-01  1.3362E+00  1.2449E+00  7.9421E-01  1.0780E+00  9.8975E-01  4.6622E-01  9.0047E-01  1.4465E+00  1.1514E+00
             9.9088E-01
 PARAMETER:  9.8591E-02  3.8984E-01  3.1907E-01 -1.3041E-01  1.7514E-01  8.9693E-02 -6.6310E-01 -4.8353E-03  4.6914E-01  2.4098E-01
             9.0836E-02
 GRADIENT:   7.1133E-04 -6.9963E-02  2.1430E-01 -1.4006E-02  1.2799E-01 -9.5006E-04 -1.3983E-02 -5.4495E-03  9.4419E-02  1.7017E-02
            -1.1248E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1205
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9219E-04 -3.5689E-02 -1.6314E-02  9.2519E-03 -2.7246E-02
 SE:             2.9843E-02  1.1716E-02  8.6697E-03  2.7160E-02  2.4840E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9486E-01  2.3178E-03  5.9875E-02  7.3337E-01  2.7269E-01

 ETASHRINKSD(%)  2.1841E-02  6.0750E+01  7.0956E+01  9.0119E+00  1.6784E+01
 ETASHRINKVR(%)  4.3677E-02  8.4594E+01  9.1564E+01  1.7212E+01  3.0750E+01
 EBVSHRINKSD(%)  4.3768E-01  6.2509E+01  7.3975E+01  8.4252E+00  1.3939E+01
 EBVSHRINKVR(%)  8.7344E-01  8.5944E+01  9.3227E+01  1.6140E+01  2.5934E+01
 RELATIVEINF(%)  9.9016E+01  1.2724E+00  2.3434E+00  8.9681E+00  2.4951E+01
 EPSSHRINKSD(%)  4.3859E+01
 EPSSHRINKVR(%)  6.8482E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1679.7371511363406     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -944.58632457260239     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1679.737       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  1.34E+00  1.24E+00  7.94E-01  1.08E+00  9.90E-01  4.66E-01  9.00E-01  1.45E+00  1.15E+00  9.91E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       33.882
Stop Time:
Sun Oct 24 03:58:03 CDT 2021
