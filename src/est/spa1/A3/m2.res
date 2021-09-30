Wed Sep 29 23:53:39 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat2.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   218.850061168028        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6957E+02  1.1373E+02  1.9932E+02  5.7312E+01  2.2312E+02  3.7165E+01 -9.8407E+01 -1.7087E+02 -2.5003E+01 -2.0182E+02
            -4.0697E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1425.92566654204        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.8246E-01  9.8918E-01  9.0513E-01  1.1047E+00  8.9789E-01  8.5516E-01  1.0046E+00  9.9738E-01  9.2313E-01  1.0922E+00
             5.3324E+00
 PARAMETER:  8.2308E-02  8.9119E-02  3.2148E-04  1.9959E-01 -7.7077E-03 -5.6470E-02  1.0460E-01  9.7379E-02  2.0015E-02  1.8820E-01
             1.7738E+00
 GRADIENT:  -1.1220E+02 -8.9394E+00 -9.8417E+00 -9.2792E+00 -2.3197E+01 -1.4514E+01  1.4758E+01  7.8117E+00  2.9702E+01  3.0305E+01
             3.4865E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1472.06397999176        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      166
 NPARAMETR:  9.8919E-01  6.7893E-01  3.2011E-01  1.2252E+00  3.9882E-01  9.0117E-01  5.8483E-01  1.2497E-01  8.4409E-01  4.1593E-01
             4.7033E+00
 PARAMETER:  8.9133E-02 -2.8724E-01 -1.0391E+00  3.0310E-01 -8.1924E-01 -4.0661E-03 -4.3644E-01 -1.9797E+00 -6.9500E-02 -7.7725E-01
             1.6483E+00
 GRADIENT:  -9.6698E+01  7.0168E+01 -2.1416E+00  1.4693E+02 -4.8984E+01 -1.7489E+01  3.5926E-01  3.1826E-01 -7.3451E+00  8.7258E+00
             2.5719E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1530.57039462207        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      259
 NPARAMETR:  1.0001E+00  5.4106E-01  2.0303E-01  1.0900E+00  2.8566E-01  9.2790E-01  2.0316E-01  1.0000E-02  1.3810E+00  3.7594E-01
             2.9879E+00
 PARAMETER:  1.0008E-01 -5.1423E-01 -1.4944E+00  1.8614E-01 -1.1529E+00  2.5173E-02 -1.4937E+00 -4.9644E+00  4.2277E-01 -8.7832E-01
             1.1946E+00
 GRADIENT:   1.5178E+01  1.2222E+02  1.1209E+01  4.5331E+01 -1.0017E+02 -1.5958E+01 -2.4646E+00  0.0000E+00  1.6181E+01 -3.0222E+01
            -8.5900E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1551.54615766226        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      434
 NPARAMETR:  9.8775E-01  3.9164E-01  2.0536E-01  1.1321E+00  2.5711E-01  9.7211E-01  1.7605E-01  1.0000E-02  1.2051E+00  5.8254E-01
             2.9622E+00
 PARAMETER:  8.7676E-02 -8.3740E-01 -1.4830E+00  2.2411E-01 -1.2583E+00  7.1717E-02 -1.6370E+00 -6.1865E+00  2.8655E-01 -4.4035E-01
             1.1859E+00
 GRADIENT:  -1.7789E+01  2.1133E+01 -7.6819E+00  6.5587E+01 -7.9925E+00  3.4365E+00 -9.6115E-03  0.0000E+00 -9.2559E+00  3.9012E+00
            -2.2796E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1555.02042535775        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  9.9504E-01  3.3304E-01  1.9005E-01  1.0627E+00  2.3570E-01  9.5695E-01  6.3082E-02  1.0000E-02  1.2454E+00  5.5767E-01
             3.0176E+00
 PARAMETER:  9.5028E-02 -9.9950E-01 -1.5605E+00  1.6085E-01 -1.3452E+00  5.6000E-02 -2.6633E+00 -7.3052E+00  3.1946E-01 -4.8398E-01
             1.2045E+00
 GRADIENT:   3.3358E+01  6.2826E+00  2.2497E+01  1.1869E+01  9.5804E+01  2.5640E+00  1.0681E-02  0.0000E+00  4.4308E+00  9.8472E-01
             1.2398E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1555.02260138406        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.9514E-01  3.3254E-01  1.8961E-01  1.0621E+00  2.3535E-01  9.5659E-01  3.2390E-02  1.0000E-02  1.2455E+00  5.5826E-01
             3.0169E+00
 PARAMETER:  9.5124E-02 -1.0010E+00 -1.5628E+00  1.6025E-01 -1.3467E+00  5.5625E-02 -3.3299E+00 -7.3052E+00  3.1957E-01 -4.8292E-01
             1.2042E+00
 GRADIENT:   3.0500E-01 -1.5988E-01  2.3194E-01 -4.0408E-01  2.7810E-02 -1.2041E-01  1.1220E-03  0.0000E+00 -5.3856E-03  8.9654E-02
             2.9203E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1555.02379713581        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.9501E-01  3.3271E-01  1.8917E-01  1.0618E+00  2.3510E-01  9.5707E-01  1.1570E-02  1.0000E-02  1.2469E+00  5.5828E-01
             3.0156E+00
 PARAMETER:  9.4999E-02 -1.0005E+00 -1.5651E+00  1.5993E-01 -1.3477E+00  5.6118E-02 -4.3593E+00 -7.3052E+00  3.2064E-01 -4.8289E-01
             1.2038E+00
 GRADIENT:   4.7336E-02 -1.0330E-01 -7.3412E-02 -1.1074E-02  8.2418E-02  5.6307E-03  1.5125E-04  0.0000E+00  1.8990E-02  7.3891E-03
             8.3408E-04

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1555.02384132508        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      940
 NPARAMETR:  9.9497E-01  3.3284E-01  1.8917E-01  1.0618E+00  2.3511E-01  9.5704E-01  1.0000E-02  1.0000E-02  1.2467E+00  5.5824E-01
             3.0156E+00
 PARAMETER:  9.4962E-02 -1.0001E+00 -1.5651E+00  1.5994E-01 -1.3477E+00  5.6085E-02 -4.6971E+00 -7.3052E+00  3.2049E-01 -4.8297E-01
             1.2038E+00
 GRADIENT:  -5.6254E-02  1.6595E-02  5.4424E-02  4.2842E-02 -2.3395E-01 -1.0488E-02  0.0000E+00  0.0000E+00 -3.2583E-02 -9.3698E-03
            -2.2853E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      940
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0941E-03 -8.1548E-05  1.9538E-04 -1.0849E-02  2.7111E-03
 SE:             2.8894E-02  1.1555E-04  2.2670E-04  2.6630E-02  2.1624E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6979E-01  4.8036E-01  3.8877E-01  6.8371E-01  9.0023E-01

 ETASHRINKSD(%)  3.2024E+00  9.9613E+01  9.9241E+01  1.0785E+01  2.7556E+01
 ETASHRINKVR(%)  6.3023E+00  9.9999E+01  9.9994E+01  2.0406E+01  4.7518E+01
 EBVSHRINKSD(%)  2.9305E+00  9.9613E+01  9.9261E+01  8.8929E+00  2.7851E+01
 EBVSHRINKVR(%)  5.7750E+00  9.9999E+01  9.9995E+01  1.6995E+01  4.7945E+01
 RELATIVEINF(%)  9.4023E+01  2.5054E-04  3.3282E-04  3.6521E+01  1.8073E+00
 EPSSHRINKSD(%)  2.4379E+01
 EPSSHRINKVR(%)  4.2814E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1555.0238413250834     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -636.08530812041067     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.94
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1555.024       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  3.33E-01  1.89E-01  1.06E+00  2.35E-01  9.57E-01  1.00E-02  1.00E-02  1.25E+00  5.58E-01  3.02E+00
 


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
+        1.16E+03
 
 TH 2
+       -6.26E+01  1.50E+03
 
 TH 3
+       -9.41E+01  3.65E+03  1.67E+04
 
 TH 4
+       -1.93E+01  1.43E+02 -5.23E+02  5.28E+02
 
 TH 5
+        2.49E+02 -6.13E+03 -2.16E+04 -4.71E+02  3.40E+04
 
 TH 6
+        3.76E+00 -1.47E+01  5.44E+01 -1.18E+01 -7.03E+00  1.91E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.37E+00 -3.50E+01  1.78E+02 -7.70E+00  5.77E+01 -1.57E-01  0.00E+00  0.00E+00  8.51E+01
 
 TH10
+       -7.85E+00 -5.49E+01 -3.30E+01  6.06E+00  2.00E+02  3.89E+00  0.00E+00  0.00E+00  6.65E+00  1.73E+02
 
 TH11
+       -1.83E+01 -7.31E+00 -6.36E+01 -6.65E+00  5.14E+01  2.61E+00  0.00E+00  0.00E+00  4.32E+00  2.42E+01  5.43E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       22.471
Stop Time:
Wed Sep 29 23:54:17 CDT 2021
