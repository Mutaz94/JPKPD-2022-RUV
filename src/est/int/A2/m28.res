Fri Sep 24 21:18:18 CDT 2021
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
$DATA ../../../../data/int/A2/dat28.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2591.43110191606        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.8336E+01 -1.3175E+01  5.4438E+01 -1.5277E+02  8.7267E+01  4.4063E+00 -8.8216E+01 -6.3550E+00 -2.1165E+01 -4.9307E+01
            -2.3841E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3208.22491981399        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.7870E-01  9.4251E-01  9.2128E-01  1.1339E+00  8.9307E-01  9.4963E-01  1.1879E+00  5.0318E-01  9.7212E-01  7.8998E-01
             1.8260E+00
 PARAMETER:  7.8470E-02  4.0790E-02  1.8009E-02  2.2569E-01 -1.3095E-02  4.8320E-02  2.7217E-01 -5.8680E-01  7.1724E-02 -1.3575E-01
             7.0215E-01
 GRADIENT:   8.4620E+00  2.9237E+00  1.2299E+01  1.0760E+01  3.1001E+01 -1.4776E+01  1.5217E+00  3.3564E+00 -7.0578E+00 -1.8128E+01
             9.0154E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3213.17900422295        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9046E-01  8.3521E-01  7.8263E-01  1.1941E+00  7.5344E-01  9.8077E-01  1.1897E+00  3.1849E-01  9.7860E-01  8.1760E-01
             1.8121E+00
 PARAMETER:  9.0417E-02 -8.0074E-02 -1.4510E-01  2.7736E-01 -1.8310E-01  8.0585E-02  2.7369E-01 -1.0442E+00  7.8369E-02 -1.0138E-01
             6.9451E-01
 GRADIENT:   3.7510E+01  2.9596E+01 -1.0777E+01  2.1294E+01 -8.8431E+00 -1.8012E+00 -2.1823E-01  3.4872E+00 -2.6526E+00 -8.3617E-01
             4.5279E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3214.54178285916        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.7133E-01  7.2096E-01  7.0721E-01  1.2239E+00  6.6821E-01  9.8886E-01  1.2564E+00  1.7978E-01  9.8189E-01  7.6199E-01
             1.8010E+00
 PARAMETER:  7.0913E-02 -2.2718E-01 -2.4643E-01  3.0201E-01 -3.0315E-01  8.8794E-02  3.2824E-01 -1.6160E+00  8.1722E-02 -1.7183E-01
             6.8834E-01
 GRADIENT:  -5.0644E+00  8.4943E-01 -3.7553E+00  5.4349E+00  2.7004E+00  1.7185E+00 -1.7450E-01  1.4088E+00  4.3721E-01  6.3840E-01
            -1.8799E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3214.56153467479        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.7236E-01  7.0187E-01  6.9008E-01  1.2280E+00  6.5063E-01  9.8693E-01  1.2717E+00  1.3938E-01  9.8128E-01  7.4863E-01
             1.8004E+00
 PARAMETER:  7.1969E-02 -2.5401E-01 -2.7095E-01  3.0535E-01 -3.2982E-01  8.6844E-02  3.4033E-01 -1.8705E+00  8.1106E-02 -1.8952E-01
             6.8802E-01
 GRADIENT:  -2.5386E+00 -1.1976E-01 -2.4562E+00  1.5941E+00  1.2756E+00  1.0404E+00 -6.6088E-02  8.8514E-01  2.5691E-01  5.0825E-01
             6.1251E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3214.59703373364        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.7334E-01  6.8732E-01  6.7678E-01  1.2308E+00  6.3716E-01  9.8492E-01  1.2832E+00  8.4196E-02  9.8087E-01  7.3887E-01
             1.7999E+00
 PARAMETER:  7.2980E-02 -2.7495E-01 -2.9041E-01  3.0768E-01 -3.5074E-01  8.4810E-02  3.4937E-01 -2.3746E+00  8.0683E-02 -2.0263E-01
             6.8774E-01
 GRADIENT:  -1.3922E-01 -7.1765E-01 -1.1558E+00 -1.4623E+00  7.9582E-02  3.1311E-01  1.9260E-02  3.3281E-01  7.0286E-02  3.2988E-01
             9.1412E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3215.35720403167        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  9.7200E-01  7.3707E-01  7.2025E-01  1.2185E+00  6.8308E-01  9.8762E-01  1.2480E+00  1.0000E-02  9.8052E-01  7.7077E-01
             1.8024E+00
 PARAMETER:  7.1605E-02 -2.0507E-01 -2.2816E-01  2.9761E-01 -2.8114E-01  8.7547E-02  3.2154E-01 -6.8352E+00  8.0326E-02 -1.6037E-01
             6.8911E-01
 GRADIENT:  -1.4757E+01 -1.2798E+00 -4.0430E+00 -3.1802E+00 -5.7274E-01  2.2978E-01 -6.1049E-01  0.0000E+00 -4.7078E-01 -8.7734E-01
            -1.9969E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3215.54971660001        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  9.7862E-01  7.7274E-01  7.5402E-01  1.2074E+00  7.1696E-01  9.8674E-01  1.2218E+00  1.0000E-02  9.8219E-01  8.0454E-01
             1.8066E+00
 PARAMETER:  7.8391E-02 -1.5781E-01 -1.8233E-01  2.8844E-01 -2.3273E-01  8.6653E-02  3.0029E-01 -8.2515E+00  8.2026E-02 -1.1749E-01
             6.9146E-01
 GRADIENT:  -2.2175E-04 -2.1253E-03 -1.5674E-02 -4.2649E-03  1.4918E-02  1.1212E-03  4.5675E-04  0.0000E+00 -1.4488E-03 -2.5630E-03
            -3.8394E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -3215.54971660001        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  9.7862E-01  7.7274E-01  7.5402E-01  1.2074E+00  7.1696E-01  9.8674E-01  1.2218E+00  1.0000E-02  9.8219E-01  8.0454E-01
             1.8066E+00
 PARAMETER:  7.8391E-02 -1.5781E-01 -1.8233E-01  2.8844E-01 -2.3273E-01  8.6653E-02  3.0029E-01 -8.2515E+00  8.2026E-02 -1.1749E-01
             6.9146E-01
 GRADIENT:  -2.2175E-04 -2.1253E-03 -1.5674E-02 -4.2649E-03  1.4918E-02  1.1212E-03  4.5675E-04  0.0000E+00 -1.4488E-03 -2.5630E-03
            -3.8394E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      669
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6842E-04 -1.0284E-02 -1.4957E-04  1.8758E-03 -1.1483E-02
 SE:             2.9674E-02  2.2956E-02  2.1288E-04  2.8022E-02  2.3537E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7397E-01  6.5415E-01  4.8232E-01  9.4663E-01  6.2565E-01

 ETASHRINKSD(%)  5.8876E-01  2.3095E+01  9.9287E+01  6.1223E+00  2.1146E+01
 ETASHRINKVR(%)  1.1740E+00  4.0856E+01  9.9995E+01  1.1870E+01  3.7821E+01
 EBVSHRINKSD(%)  8.5462E-01  2.2422E+01  9.9191E+01  6.3250E+00  2.2570E+01
 EBVSHRINKVR(%)  1.7019E+00  3.9816E+01  9.9993E+01  1.2250E+01  4.0046E+01
 RELATIVEINF(%)  9.8281E+01  1.7593E+01  1.7146E-03  6.0087E+01  9.0283E+00
 EPSSHRINKSD(%)  1.8359E+01
 EPSSHRINKVR(%)  3.3347E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3215.5497166000087     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1561.4603568315979     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3215.550       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  7.73E-01  7.54E-01  1.21E+00  7.17E-01  9.87E-01  1.22E+00  1.00E-02  9.82E-01  8.05E-01  1.81E+00
 


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
+        1.17E+03
 
 TH 2
+       -4.23E+00  7.56E+02
 
 TH 3
+       -3.77E+00  3.82E+01  1.25E+03
 
 TH 4
+       -7.47E+00  1.51E+02 -7.98E+01  6.93E+02
 
 TH 5
+       -3.83E+00 -7.23E+02 -1.14E+03  2.07E+02  1.97E+03
 
 TH 6
+        1.97E+00 -1.76E+00  2.88E+00 -3.52E+00 -3.07E+00  1.98E+02
 
 TH 7
+        1.32E-01  2.02E+01  1.90E-01  5.49E+00 -2.79E+00 -2.28E-01  4.76E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.16E-01 -1.67E+01  3.44E+01  1.31E+01 -1.44E+01  4.92E-01  4.16E+00  0.00E+00  1.63E+02
 
 TH10
+        9.91E-01 -1.11E+01 -7.40E+01 -4.70E+00 -4.55E+00  5.01E-01  2.42E+01  0.00E+00  4.91E+00  1.19E+02
 
 TH11
+       -1.20E+01 -1.85E+01 -1.82E+01 -1.25E+01  4.66E-01  1.71E+00  7.35E+00  0.00E+00  5.51E+00  1.23E+01  3.49E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.923
Stop Time:
Fri Sep 24 21:18:45 CDT 2021
