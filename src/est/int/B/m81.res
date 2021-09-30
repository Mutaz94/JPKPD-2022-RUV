Tue Sep 28 20:59:52 CDT 2021
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
$DATA ../../../../data/int/B/dat81.csv ignore=@
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
Current Date:       28 SEP 2021
Days until program expires : 201
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3382.02968683753        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0765E+02  5.4648E+01  8.9119E+01  4.7486E+01  1.3798E+02  2.9831E+01 -3.7356E+01 -2.2247E+02 -3.6424E+00 -2.9868E+01
            -7.1370E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3857.41972188782        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  8.6588E-01  8.8574E-01  8.6450E-01  1.0236E+00  8.4029E-01  9.7607E-01  9.5053E-01  9.8250E-01  9.2517E-01  9.9630E-01
             1.2188E+00
 PARAMETER: -4.4004E-02 -2.1333E-02 -4.5607E-02  1.2335E-01 -7.4011E-02  7.5783E-02  4.9261E-02  8.2348E-02  2.2218E-02  9.6294E-02
             2.9785E-01
 GRADIENT:  -6.8198E+01 -2.8535E+01 -5.4356E+00 -1.7991E+00 -2.5937E+00 -2.5758E+01  1.6549E+00  2.4985E+01  8.8458E+00  7.8634E+00
             3.9835E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3862.58760145297        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  8.6543E-01  8.4657E-01  8.4080E-01  1.0363E+00  8.2236E-01  9.8107E-01  1.1415E+00  5.8599E-01  8.6503E-01  9.5276E-01
             1.2064E+00
 PARAMETER: -4.4530E-02 -6.6557E-02 -7.3402E-02  1.3567E-01 -9.5581E-02  8.0886E-02  2.3236E-01 -4.3445E-01 -4.4991E-02  5.1609E-02
             2.8767E-01
 GRADIENT:  -6.0162E+01 -3.0264E+01  5.6429E+00  1.0559E+01  1.7693E+01 -2.1687E+01  2.7044E+01  2.7195E-02 -1.0270E+01  1.1315E+01
             3.6791E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3920.24235519457        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  9.2772E-01  9.8943E-01  9.1768E-01  1.0097E+00  9.2105E-01  9.7872E-01  7.4068E-01  6.3241E-01  8.2138E-01  1.1256E+00
             1.0737E+00
 PARAMETER:  2.4974E-02  8.9370E-02  1.4091E-02  1.0962E-01  1.7754E-02  7.8486E-02 -2.0019E-01 -3.5821E-01 -9.6773E-02  2.1834E-01
             1.7115E-01
 GRADIENT:  -2.2937E+02 -1.8477E+01  1.0929E+01 -5.2275E+01 -1.8151E+01 -7.2805E+01 -2.0768E+01 -8.8254E+00 -4.5790E+01  2.4638E+00
             1.6062E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3942.43312757024        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      460
 NPARAMETR:  1.0083E+00  9.3681E-01  8.6943E-01  1.0340E+00  8.7124E-01  1.0105E+00  8.0421E-01  7.6586E-01  8.9333E-01  1.0401E+00
             1.0209E+00
 PARAMETER:  1.0827E-01  3.4722E-02 -3.9915E-02  1.3346E-01 -3.7841E-02  1.1041E-01 -1.1790E-01 -1.6675E-01 -1.2802E-02  1.3934E-01
             1.2072E-01
 GRADIENT:  -3.3165E+01 -2.9261E+01 -1.2651E+01 -3.9780E+01 -2.3489E+01 -3.5254E+01 -1.6399E+01 -6.4248E-02 -1.3164E+01 -5.9104E+00
             8.5270E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3944.01087651251        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  1.0091E+00  9.3754E-01  8.7003E-01  1.0336E+00  8.9040E-01  1.0802E+00  8.0590E-01  7.6522E-01  8.9341E-01  1.0398E+00
             1.0216E+00
 PARAMETER:  1.0906E-01  3.5501E-02 -3.9229E-02  1.3302E-01 -1.6085E-02  1.7718E-01 -1.1580E-01 -1.6759E-01 -1.2711E-02  1.3902E-01
             1.2133E-01
 GRADIENT:   4.7668E+02  2.5858E+01 -1.6454E+01  1.1368E+02  6.3349E+01  1.4911E+02 -1.1340E+01  1.5627E-01 -5.1582E+00  9.2518E-01
             8.8673E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3945.18342240211        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      796            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0110E+00  9.4035E-01  8.7207E-01  1.0354E+00  8.8452E-01  1.0928E+00  8.2312E-01  7.6590E-01  8.9768E-01  1.0432E+00
             1.0161E+00
 PARAMETER:  1.1091E-01  3.8497E-02 -3.6881E-02  1.3480E-01 -2.2707E-02  1.8870E-01 -9.4653E-02 -1.6671E-01 -7.9437E-03  1.4227E-01
             1.1602E-01
 GRADIENT:   4.9491E+02  3.8978E+01 -1.0876E+01  1.2284E+02  4.9939E+01  1.6505E+02 -9.3850E+00  4.7035E-02 -2.9701E+00  3.4767E+00
             7.8701E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3945.69536184895        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  1.0116E+00  9.4153E-01  8.7428E-01  1.0369E+00  8.8444E-01  1.0912E+00  8.3132E-01  7.6619E-01  9.0195E-01  1.0446E+00
             1.0136E+00
 PARAMETER:  1.1154E-01  3.9754E-02 -3.4358E-02  1.3625E-01 -2.2802E-02  1.8732E-01 -8.4743E-02 -1.6632E-01 -3.1997E-03  1.4366E-01
             1.1348E-01
 GRADIENT:  -2.1761E+01 -3.0439E+01 -1.6720E+01 -2.4258E+01 -7.2320E+00 -9.1186E-01 -1.2939E+01 -8.5373E-01 -9.0994E+00 -3.2251E+00
             7.1883E+01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -3945.69536184895        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:      995
 NPARAMETR:  1.0116E+00  9.4161E-01  8.7428E-01  1.0369E+00  8.8436E-01  1.0933E+00  8.3132E-01  7.6619E-01  9.0195E-01  1.0445E+00
             1.0137E+00
 PARAMETER:  1.1154E-01  3.9754E-02 -3.4358E-02  1.3625E-01 -2.2802E-02  1.8732E-01 -8.4743E-02 -1.6632E-01 -3.1997E-03  1.4366E-01
             1.1348E-01
 GRADIENT:  -2.0905E+01 -1.5320E+05  8.4790E+01 -7.2245E+01  3.0636E+05 -2.1608E+00 -7.1292E+01 -2.7813E+02 -7.4744E+01  2.1323E+05
            -2.7030E+05
 NUMSIGDIG:         6.7         2.3         6.2         6.1         2.3         1.2         6.2         5.4         6.2         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      995
 NO. OF SIG. DIGITS IN FINAL EST.:  1.2
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2094E-02 -1.6124E-02 -1.2565E-02  2.0874E-02 -1.7590E-02
 SE:             3.0067E-02  2.2340E-02  1.6514E-02  2.8712E-02  2.6209E-02
 N:                     100         100         100         100         100

 P VAL.:         6.8749E-01  4.7045E-01  4.4674E-01  4.6721E-01  5.0211E-01

 ETASHRINKSD(%)  1.0000E-10  2.5159E+01  4.4676E+01  3.8123E+00  1.2198E+01
 ETASHRINKVR(%)  1.0000E-10  4.3988E+01  6.9393E+01  7.4792E+00  2.2908E+01
 EBVSHRINKSD(%)  2.2022E-01  2.9815E+01  4.5970E+01  6.6315E+00  1.3110E+01
 EBVSHRINKVR(%)  4.3995E-01  5.0741E+01  7.0808E+01  1.2823E+01  2.4502E+01
 RELATIVEINF(%)  9.9557E+01  2.1259E+01  1.7061E+01  6.0138E+01  2.9374E+01
 EPSSHRINKSD(%)  2.3643E+01
 EPSSHRINKVR(%)  4.1696E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3945.6953618489524     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2291.6060020805417     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3945.695       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  9.42E-01  8.74E-01  1.04E+00  8.84E-01  1.09E+00  8.31E-01  7.66E-01  9.02E-01  1.04E+00  1.01E+00
 


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
+        1.20E+08
 
 TH 2
+       -8.35E+02  1.73E+08
 
 TH 3
+       -7.77E+07 -3.09E+04  2.00E+08
 
 TH 4
+        4.80E+07  1.27E+04  6.20E+07  3.84E+07
 
 TH 5
+        8.85E+02  9.20E+07  9.91E+07  6.13E+07  1.96E+08
 
 TH 6
+        1.94E+01 -1.36E+03  1.67E+00  4.61E-01  1.44E+03  1.68E+02
 
 TH 7
+        8.17E+07  1.87E+04 -1.05E+08 -6.52E+07 -1.04E+08 -6.19E-01  2.22E+08
 
 TH 8
+       -5.33E+07  6.39E+07  6.88E+07 -6.40E+04 -6.80E+07 -7.43E+00 -6.57E+00  9.43E+07
 
 TH 9
+        4.57E+00  1.93E+04 -9.71E+07 -6.01E+07 -9.60E+07  5.45E+00 -1.94E+04  1.84E+02  9.41E+07
 
 TH10
+        1.77E-01  4.99E+03  5.84E+07  7.34E+00 -5.38E+03  8.52E+02 -6.14E+07 -3.79E-01 -5.66E+07  3.40E+07
 
 TH11
+        5.89E+07 -7.06E+07 -7.61E+07  4.71E+07  7.52E+07 -1.11E+03  1.53E+04  1.23E+02  1.58E+04  4.12E+03  5.81E+07
 
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
 #CPUT: Total CPU Time in Seconds,       38.612
Stop Time:
Tue Sep 28 21:00:32 CDT 2021
