Wed Sep 29 12:21:38 CDT 2021
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
$DATA ../../../../data/spa/A1/dat78.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1465.68828302923        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7953E+02 -4.5667E+01 -5.0663E+01  2.7051E+00  9.4737E+01  3.2786E+01 -3.1010E+01  1.5415E+01 -1.9122E+01 -1.4918E+01
            -3.1317E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1530.08298886180        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.6287E-01  1.0442E+00  1.2127E+00  1.0258E+00  1.0467E+00  1.0057E+00  1.1816E+00  8.6181E-01  1.0951E+00  9.3187E-01
             1.9143E+00
 PARAMETER:  6.2165E-02  1.4329E-01  2.9281E-01  1.2547E-01  1.4567E-01  1.0568E-01  2.6685E-01 -4.8718E-02  1.9087E-01  2.9443E-02
             7.4933E-01
 GRADIENT:   8.2826E+01 -6.1161E+00 -1.8198E+01  8.5622E+00  2.1082E+01  2.0424E+01  4.5895E+00  5.9822E+00  1.2570E+01  1.1845E+00
             9.0701E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1534.45246484099        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.5296E-01  1.0077E+00  1.3760E+00  1.1017E+00  1.0687E+00  9.9457E-01  1.0869E+00  2.9795E-01  9.9634E-01  1.2060E+00
             1.8089E+00
 PARAMETER:  5.1822E-02  1.0765E-01  4.1915E-01  1.9681E-01  1.6643E-01  9.4559E-02  1.8329E-01 -1.1108E+00  9.6334E-02  2.8728E-01
             6.9272E-01
 GRADIENT:   7.2528E+01  3.7743E+01 -5.8481E+00  7.9152E+01 -1.3245E+01  1.6263E+01 -6.8767E+00  5.7174E-01 -1.2168E+01  2.1960E+01
             7.7170E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1538.21584332927        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.2204E-01  9.7033E-01  1.1047E+00  1.0692E+00  9.4837E-01  9.6546E-01  1.1879E+00  1.7046E-01  1.0540E+00  1.0347E+00
             1.6121E+00
 PARAMETER:  1.8837E-02  6.9883E-02  1.9959E-01  1.6687E-01  4.6988E-02  6.4847E-02  2.7219E-01 -1.6692E+00  1.5261E-01  1.3416E-01
             5.7755E-01
 GRADIENT:   3.0948E+01  1.2071E+01 -3.5530E+00  3.3050E+01 -1.6791E+01  3.3063E+00 -6.5589E-01  3.4220E-01  6.3778E+00  1.1488E+01
             3.0811E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1544.58245739396        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      383
 NPARAMETR:  9.7624E-01  8.7625E-01  1.4762E+00  1.1572E+00  1.0578E+00  9.8621E-01  1.3201E+00  2.1991E-01  1.0283E+00  1.1266E+00
             1.5184E+00
 PARAMETER:  7.5955E-02 -3.2107E-02  4.8949E-01  2.4603E-01  1.5623E-01  8.6114E-02  3.7774E-01 -1.4146E+00  1.2790E-01  2.1920E-01
             5.1764E-01
 GRADIENT:   1.2606E+01  8.0898E+00  3.7453E+00  5.7067E-01 -1.0087E+01  2.5395E+00 -1.2961E-01  2.7307E-01  1.2735E+00  1.4869E+00
            -6.1051E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1547.04156685703        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      558
 NPARAMETR:  9.6352E-01  4.4643E-01  1.4989E+00  1.4411E+00  9.2238E-01  9.5806E-01  1.9899E+00  4.8892E-02  8.9071E-01  1.0607E+00
             1.5288E+00
 PARAMETER:  6.2842E-02 -7.0647E-01  5.0475E-01  4.6542E-01  1.9199E-02  5.7155E-02  7.8810E-01 -2.9181E+00 -1.5735E-02  1.5893E-01
             5.2452E-01
 GRADIENT:  -9.7076E+00  9.6981E+00  5.6380E+00  2.6422E+01 -1.2028E+01 -7.4123E+00 -6.7411E-01  1.7237E-02 -2.8009E+00 -3.9665E-01
             9.9855E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1548.40753757497        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      733
 NPARAMETR:  9.6500E-01  2.0522E-01  1.4296E+00  1.5731E+00  8.3987E-01  9.8336E-01  3.1746E+00  1.0000E-02  8.5322E-01  1.0375E+00
             1.5075E+00
 PARAMETER:  6.4374E-02 -1.4837E+00  4.5739E-01  5.5306E-01 -7.4506E-02  8.3218E-02  1.2552E+00 -5.4085E+00 -5.8737E-02  1.3679E-01
             5.1043E-01
 GRADIENT:   1.9852E+00  2.6227E+00  2.4995E+00  1.3577E+01 -6.6886E+00  3.6047E+00  7.3390E-02  0.0000E+00  1.9615E-01  1.9616E-01
            -9.8246E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1548.70224283633        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  9.6199E-01  9.3175E-02  1.4775E+00  1.6382E+00  8.3048E-01  9.6961E-01  4.6690E+00  1.0000E-02  8.3434E-01  1.0498E+00
             1.5071E+00
 PARAMETER:  6.1249E-02 -2.2733E+00  4.9038E-01  5.9361E-01 -8.5754E-02  6.9134E-02  1.6409E+00 -8.0859E+00 -8.1110E-02  1.4865E-01
             5.1020E-01
 GRADIENT:  -1.3383E+00  8.4138E-01  2.1244E+00  8.2867E+00 -4.6901E+00 -1.1084E+00 -1.8835E-01  0.0000E+00  3.3885E-01  4.1016E-01
            -1.2356E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1548.99931010356        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1086
 NPARAMETR:  9.6181E-01  4.1131E-02  1.5359E+00  1.6675E+00  8.4099E-01  9.7050E-01  6.7145E+00  1.0000E-02  8.1620E-01  1.0711E+00
             1.5069E+00
 PARAMETER:  6.1064E-02 -3.0910E+00  5.2910E-01  6.1130E-01 -7.3177E-02  7.0060E-02  2.0043E+00 -1.1002E+01 -1.0309E-01  1.6869E-01
             5.1008E-01
 GRADIENT:  -2.4888E-01 -1.7138E-01  7.7146E-01  3.1373E+00 -1.2596E+00 -2.6449E-01 -6.2033E-01  0.0000E+00 -3.6719E-01  9.7247E-01
             1.0844E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1549.02074167658        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  9.6201E-01  3.9405E-02  1.5320E+00  1.6670E+00  8.3909E-01  9.7035E-01  6.7841E+00  1.0000E-02  8.1630E-01  1.0682E+00
             1.5075E+00
 PARAMETER:  6.1115E-02 -3.1285E+00  5.2550E-01  6.1035E-01 -7.4440E-02  7.0140E-02  2.0172E+00 -1.1153E+01 -1.0231E-01  1.6527E-01
             5.1017E-01
 GRADIENT:  -8.8695E-01  4.8792E+00 -1.5734E+01 -2.3792E+01  3.5900E+00  1.2624E-01  7.7370E+00  0.0000E+00  3.0850E-01 -5.5490E+00
            -4.4278E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1150
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0723E-04  1.1685E-02 -5.0514E-05 -1.7271E-02 -2.7243E-02
 SE:             2.9596E-02  8.4397E-03  1.0693E-04  2.8207E-02  2.2682E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9711E-01  1.6621E-01  6.3666E-01  5.4034E-01  2.2972E-01

 ETASHRINKSD(%)  8.4807E-01  7.1726E+01  9.9642E+01  5.5024E+00  2.4012E+01
 ETASHRINKVR(%)  1.6889E+00  9.2006E+01  9.9999E+01  1.0702E+01  4.2258E+01
 EBVSHRINKSD(%)  1.0067E+00  8.0497E+01  9.9571E+01  5.2387E+00  2.0220E+01
 EBVSHRINKVR(%)  2.0032E+00  9.6196E+01  9.9998E+01  1.0203E+01  3.6352E+01
 RELATIVEINF(%)  9.7855E+01  1.9246E+00  1.7553E-04  4.3252E+01  6.0862E+00
 EPSSHRINKSD(%)  3.8402E+01
 EPSSHRINKVR(%)  6.2057E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1549.0207416765793     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -813.86991511284111     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.57
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1549.021       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.62E-01  3.96E-02  1.53E+00  1.67E+00  8.40E-01  9.71E-01  6.80E+00  1.00E-02  8.17E-01  1.07E+00  1.51E+00
 


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
+        1.25E+03
 
 TH 2
+       -1.55E+02  1.37E+05
 
 TH 3
+        3.18E+01 -7.82E+02  7.81E+02
 
 TH 4
+       -1.59E+00 -2.70E+03  6.06E+01  1.11E+03
 
 TH 5
+       -7.21E+01  1.62E+03 -7.19E+02 -2.31E+02  2.19E+03
 
 TH 6
+       -5.30E+00  2.31E+02 -5.84E+01 -3.14E+01  1.31E+02  2.16E+02
 
 TH 7
+       -7.11E-01  1.00E+03 -4.80E+00 -1.29E+02  1.09E+01  1.41E+00  1.25E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.40E+00  5.06E+02 -1.53E+02 -8.18E+01  1.47E+02  1.53E+01  2.70E+00  0.00E+00  2.71E+02
 
 TH10
+        2.60E+02 -2.24E+03  2.97E+03  2.28E+02 -2.54E+04 -4.84E+02 -1.31E+01  0.00E+00 -6.00E+02  1.31E+04
 
 TH11
+       -5.53E+00 -1.31E+02  2.02E+01  2.40E+00 -1.17E+02 -1.11E+01 -6.49E-01  0.00E+00 -5.68E+00  2.74E+03  1.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.728
Stop Time:
Wed Sep 29 12:22:03 CDT 2021
