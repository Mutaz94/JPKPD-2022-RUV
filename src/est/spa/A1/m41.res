Sat Sep 25 08:02:47 CDT 2021
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
$DATA ../../../../data/spa/A1/dat41.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1380.53691913981        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.5711E+01 -6.4363E+01 -3.5514E+01 -5.0948E+01  1.1580E+02  1.6980E+01 -1.3833E+01 -7.5248E-02 -3.1545E+01 -3.0806E+01
            -5.4024E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1519.41638991034        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0474E+00  9.8341E-01  1.0354E+00  1.0278E+00  9.7543E-01  9.2416E-01  1.0398E+00  9.7032E-01  1.0594E+00  1.0071E+00
             1.9533E+00
 PARAMETER:  1.4626E-01  8.3271E-02  1.3476E-01  1.2739E-01  7.5120E-02  2.1135E-02  1.3903E-01  6.9875E-02  1.5772E-01  1.0708E-01
             7.6951E-01
 GRADIENT:   4.9129E+01 -5.0725E+01 -2.5814E+01 -4.7292E+01  5.2522E+01 -9.6275E+00 -6.6028E-01  4.5199E+00  4.1729E+00 -9.3205E-01
             3.0778E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1523.98277711153        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0378E+00  7.2294E-01  8.9570E-01  1.2339E+00  7.8218E-01  9.6330E-01  1.5085E+00  2.6981E-01  9.0426E-01  9.2889E-01
             1.9178E+00
 PARAMETER:  1.3706E-01 -2.2444E-01 -1.0152E-02  3.1015E-01 -1.4567E-01  6.2608E-02  5.1115E-01 -1.2100E+00 -6.4208E-04  2.6233E-02
             7.5120E-01
 GRADIENT:   2.1761E+01  8.1900E+00 -2.2740E+01  3.2973E+01  3.6994E+01  5.6138E+00  1.1095E+01  8.6431E-02 -1.4138E+00  4.0214E+00
             1.7446E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1526.51281638522        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.0296E+00  6.2403E-01  6.5353E-01  1.2382E+00  5.9080E-01  9.5436E-01  1.3683E+00  1.3700E-01  8.9263E-01  7.3867E-01
             1.8342E+00
 PARAMETER:  1.2915E-01 -3.7156E-01 -3.2537E-01  3.1366E-01 -4.2628E-01  5.3284E-02  4.1354E-01 -1.8878E+00 -1.3579E-02 -2.0290E-01
             7.0661E-01
 GRADIENT:  -1.0083E+00  7.6919E+00  9.0127E-01  1.3806E+01 -5.4049E+00  8.7000E-01 -1.0499E+00  1.1085E-01 -2.2626E+00  2.2673E-01
            -3.3803E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1527.22280872851        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  1.0279E+00  4.4348E-01  6.7487E-01  1.3309E+00  5.5424E-01  9.4848E-01  1.6847E+00  3.8800E-02  8.5808E-01  7.5805E-01
             1.8646E+00
 PARAMETER:  1.2751E-01 -7.1310E-01 -2.9324E-01  3.8583E-01 -4.9016E-01  4.7104E-02  6.2161E-01 -3.1493E+00 -5.3056E-02 -1.7700E-01
             7.2307E-01
 GRADIENT:   1.3215E-02  2.7941E+00  4.7044E+00  4.2121E+00 -7.7905E+00 -2.4195E-01  3.5137E-01  3.2346E-03 -2.4229E-02 -1.0531E-01
             2.2053E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1527.72918536644        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  1.0302E+00  3.2958E-01  7.7995E-01  1.4162E+00  5.8811E-01  9.4924E-01  1.8978E+00  1.2471E-02  8.3860E-01  8.4321E-01
             1.8634E+00
 PARAMETER:  1.2973E-01 -1.0099E+00 -1.4852E-01  4.4796E-01 -4.3085E-01  4.7907E-02  7.4072E-01 -4.2844E+00 -7.6023E-02 -7.0538E-02
             7.2238E-01
 GRADIENT:   6.6775E-01  8.9905E-01  3.5274E+00 -2.1970E+00 -4.6682E+00  4.9346E-01  4.6362E-01 -3.5155E-04  5.2392E-01  7.3434E-01
            -1.1918E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1527.81079384392        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  1.0279E+00  2.5572E-01  7.7517E-01  1.4559E+00  5.7124E-01  9.4659E-01  2.0679E+00  1.0000E-02  8.2145E-01  8.3949E-01
             1.8725E+00
 PARAMETER:  1.2755E-01 -1.2637E+00 -1.5468E-01  4.7565E-01 -4.5995E-01  4.5109E-02  8.2651E-01 -5.5144E+00 -9.6679E-02 -7.4959E-02
             7.2728E-01
 GRADIENT:  -8.0038E-01  2.7703E-01  8.6746E-01  2.5247E+00 -1.6962E+00 -1.0264E-01 -5.5780E-01  0.0000E+00 -4.2355E-01 -2.9041E-02
            -2.4374E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1527.83422379330        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      795
 NPARAMETR:  1.0271E+00  2.1470E-01  7.7851E-01  1.4765E+00  5.6554E-01  9.4597E-01  2.3412E+00  1.0000E-02  8.1287E-01  8.4125E-01
             1.8759E+00
 PARAMETER:  1.2673E-01 -1.4385E+00 -1.5037E-01  4.8967E-01 -4.6998E-01  4.4454E-02  9.5064E-01 -6.3847E+00 -1.0719E-01 -7.2862E-02
             7.2909E-01
 GRADIENT:  -1.8675E-02  1.1539E-02  3.7348E-02  3.7717E-02 -6.4345E-02 -3.9813E-03  3.3447E-03  0.0000E+00 -3.6556E-03  2.8780E-03
             2.1409E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1527.83422379330        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      817
 NPARAMETR:  1.0271E+00  2.1470E-01  7.7851E-01  1.4765E+00  5.6554E-01  9.4597E-01  2.3412E+00  1.0000E-02  8.1287E-01  8.4125E-01
             1.8759E+00
 PARAMETER:  1.2673E-01 -1.4385E+00 -1.5037E-01  4.8967E-01 -4.6998E-01  4.4454E-02  9.5064E-01 -6.3847E+00 -1.0719E-01 -7.2862E-02
             7.2909E-01
 GRADIENT:  -1.8675E-02  1.1539E-02  3.7348E-02  3.7717E-02 -6.4345E-02 -3.9813E-03  3.3447E-03  0.0000E+00 -3.6556E-03  2.8780E-03
             2.1409E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      817
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9799E-04  2.7544E-03 -7.8083E-05 -9.4768E-03 -1.5654E-02
 SE:             2.9445E-02  8.4545E-03  2.0521E-04  2.7129E-02  2.2986E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9463E-01  7.4458E-01  7.0357E-01  7.2684E-01  4.9586E-01

 ETASHRINKSD(%)  1.3569E+00  7.1676E+01  9.9313E+01  9.1156E+00  2.2993E+01
 ETASHRINKVR(%)  2.6955E+00  9.1978E+01  9.9995E+01  1.7400E+01  4.0700E+01
 EBVSHRINKSD(%)  1.4956E+00  7.3492E+01  9.9322E+01  8.6528E+00  2.1689E+01
 EBVSHRINKVR(%)  2.9688E+00  9.2973E+01  9.9995E+01  1.6557E+01  3.8674E+01
 RELATIVEINF(%)  8.9953E+01  2.4861E-01  2.5642E-04  4.9648E+00  2.0493E+00
 EPSSHRINKSD(%)  3.6872E+01
 EPSSHRINKVR(%)  6.0149E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1527.8342237933000     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -792.68339722956182     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1527.834       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.15E-01  7.79E-01  1.48E+00  5.66E-01  9.46E-01  2.34E+00  1.00E-02  8.13E-01  8.41E-01  1.88E+00
 


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
+        1.14E+03
 
 TH 2
+       -4.95E+01  3.61E+02
 
 TH 3
+        1.65E+01  2.31E+02  8.58E+02
 
 TH 4
+       -2.87E+01  3.33E+02 -1.06E+02  6.47E+02
 
 TH 5
+        1.47E+01 -5.98E+02 -1.53E+03 -5.02E+01  3.02E+03
 
 TH 6
+        3.06E+00 -7.34E+00  3.74E+00 -8.11E+00 -2.28E+00  2.12E+02
 
 TH 7
+        3.98E-02  1.02E+01 -1.10E-01 -1.27E+00 -4.19E-02  8.50E-02  1.88E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.48E+00 -2.95E+01  1.27E+01 -1.19E+00 -8.92E+00  2.38E+00  2.10E+00  0.00E+00  2.13E+02
 
 TH10
+       -4.58E+00  1.78E+01 -9.04E+00 -4.73E+00 -6.82E+01 -8.11E-01  4.52E-01  0.00E+00  5.09E+00  1.12E+02
 
 TH11
+       -1.18E+01 -3.05E+00 -1.88E+01 -9.83E+00  9.74E+00  2.49E+00  1.75E-01  0.00E+00  1.15E+01  2.34E+01  7.30E+01
 
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
 #CPUT: Total CPU Time in Seconds,       15.923
Stop Time:
Sat Sep 25 08:03:04 CDT 2021
