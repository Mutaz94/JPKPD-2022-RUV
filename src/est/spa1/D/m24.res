Thu Sep 30 02:49:45 CDT 2021
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
$DATA ../../../../data/spa1/D/dat24.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   5570.37113985138        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3051E+02  1.3458E+01 -4.8678E+00 -1.5656E+02  2.4572E+02 -1.0760E+03 -3.8541E+02 -9.9302E+01 -7.4438E+02 -4.4066E+02
            -1.2492E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -771.820154703468        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0148E+00  1.0281E+00  1.2133E+00  3.1938E+00  1.1808E+00  3.4145E+00  2.2394E+00  1.1003E+00  7.1850E+00  2.0928E+00
             9.0476E+00
 PARAMETER:  1.1473E-01  1.2767E-01  2.9331E-01  1.2612E+00  2.6622E-01  1.3280E+00  9.0619E-01  1.9555E-01  2.0720E+00  8.3851E-01
             2.3025E+00
 GRADIENT:  -3.2661E+01 -3.5672E+01 -5.1578E+01  5.1938E+01 -2.9213E+01  1.5629E+02  1.8079E+01  4.5147E+00  1.6079E+02  2.8256E+01
             3.3112E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -871.059128440899        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.4249E-01  5.6611E-01  2.1300E+01  4.1731E+00  3.4541E+02  4.2388E+00  1.1733E+01  3.0079E+00  4.1718E+00  1.0495E+02
             6.4933E+00
 PARAMETER:  4.0766E-02 -4.6897E-01  3.1587E+00  1.5287E+00  5.9447E+00  1.5443E+00  2.5624E+00  1.2012E+00  1.5283E+00  4.7535E+00
             1.9708E+00
 GRADIENT:  -1.6187E+01  1.0513E+01  1.9772E+00  1.6105E+02  3.2850E-01  2.7599E+02  1.1275E+01 -1.8545E+00  1.3794E+02  3.3794E+01
             1.8832E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -983.118892363349        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0638E+00  6.4390E-01  8.1020E+00  2.1618E+00  8.0955E+02  2.7353E+00  6.6342E+00  6.9620E+00  2.1081E+00  1.0582E+02
             7.1370E+00
 PARAMETER:  1.6182E-01 -3.4022E-01  2.1921E+00  8.7092E-01  6.7965E+00  1.1062E+00  1.9922E+00  2.0405E+00  8.4580E-01  4.7618E+00
             2.0653E+00
 GRADIENT:  -5.5449E+00  1.8564E+01  1.1941E+00  7.5663E+01  5.9245E-01  1.3716E+02  7.8454E+01  7.1648E-01  2.2279E+01  5.8944E+01
             2.3754E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -997.336978547920        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      329
 NPARAMETR:  1.0601E+00  6.8891E-01  7.1471E+00  1.9564E+00  7.9507E+02  2.5401E+00  5.6674E+00  6.1549E+00  2.0096E+00  1.0310E+02
             6.9246E+00
 PARAMETER:  1.5833E-01 -2.7265E-01  2.0667E+00  7.7112E-01  6.7784E+00  1.0322E+00  1.8347E+00  1.9173E+00  7.9795E-01  4.7357E+00
             2.0351E+00
 GRADIENT:  -2.1500E+01  1.2212E+01  1.3749E+00  2.9827E+01 -4.6069E-01  5.8623E+01  1.8417E-02  4.4336E-01  1.0601E+01  5.4541E+01
             2.0630E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1070.85470312193        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      467
 NPARAMETR:  1.1514E+00  3.8950E-01  1.3622E+00  1.7419E+00  8.4648E+02  1.8942E+00  5.7423E+00  5.2159E+00  1.8026E+00  7.9874E+01
             4.6610E+00
 PARAMETER:  2.4094E-01 -8.4290E-01  4.0909E-01  6.5497E-01  6.8411E+00  7.3879E-01  1.8479E+00  1.7517E+00  6.8924E-01  4.4804E+00
             1.6392E+00
 GRADIENT:   1.1793E+02  3.2530E+01 -1.1557E+01  9.7997E+01 -7.1918E-03  2.9648E+01  1.0335E+01  1.1055E+02  2.8031E+01  1.2623E-02
            -1.4206E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1157.46259924132        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  7.8584E-01  1.7333E-02  2.0152E-01  1.0541E+00  4.2495E+06  1.7538E+00  6.7143E+00  3.0813E+00  7.9540E-01  5.0954E+01
             5.4474E+00
 PARAMETER: -1.4100E-01 -3.9551E+00 -1.5019E+00  1.5265E-01  1.5362E+01  6.6180E-01  2.0042E+00  1.2254E+00 -1.2891E-01  4.0309E+00
             1.7951E+00
 GRADIENT:  -4.3182E-01 -2.1952E+00 -1.8198E+01  4.3087E+01  4.6338E-07  2.7633E+01 -2.4248E+00  3.8272E+01  2.2078E+01  1.6922E-11
             5.5565E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1183.12751418389        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  6.6055E-01  1.0000E-02  1.1150E-01  7.5983E-01  6.1321E+02  1.6063E+00  7.5676E+00  1.9751E+00  4.1974E-01  4.7019E+01
             5.2725E+00
 PARAMETER: -3.1469E-01 -5.0271E+00 -2.0937E+00 -1.7467E-01  6.5187E+00  5.7394E-01  2.1239E+00  7.8061E-01 -7.6811E-01  3.9506E+00
             1.7625E+00
 GRADIENT:  -1.8112E+01  0.0000E+00  5.8586E+00  6.6688E-01 -3.2484E-03  6.0670E+00 -1.0076E+00  6.4507E+00  4.0365E+00  2.9629E-03
             2.4287E+01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1183.26422974228        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      929
 NPARAMETR:  6.5939E-01  1.0000E-02  1.1024E-01  7.5611E-01  6.4050E+02  1.6052E+00  7.5682E+00  1.9641E+00  4.1652E-01  4.5523E+01
             5.2718E+00
 PARAMETER: -3.1627E-01 -5.0427E+00 -2.1040E+00 -1.7966E-01  6.4973E+00  5.7295E-01  2.1251E+00  7.7547E-01 -7.7622E-01  3.9487E+00
             1.7617E+00
 GRADIENT:   5.6343E+02  0.0000E+00  4.8624E+01 -1.0226E+03 -2.9641E-03 -3.1553E+02  8.6063E+01  1.2410E+02 -2.3300E+02  2.8358E-03
            -8.0880E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      929
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5659E-02  9.5618E-04 -1.4732E-02 -2.1383E-02 -4.4766E-04
 SE:             2.8825E-02  3.3644E-03  2.4108E-02  1.0530E-02  1.6291E-04
 N:                     100         100         100         100         100

 P VAL.:         5.8697E-01  7.7625E-01  5.4115E-01  4.2284E-02  5.9965E-03

 ETASHRINKSD(%)  3.4339E+00  8.8729E+01  1.9234E+01  6.4724E+01  9.9454E+01
 ETASHRINKVR(%)  6.7498E+00  9.8730E+01  3.4768E+01  8.7556E+01  9.9997E+01
 EBVSHRINKSD(%)  2.4933E+00  9.0128E+01  1.5373E+01  6.3369E+01  9.9565E+01
 EBVSHRINKVR(%)  4.9245E+00  9.9025E+01  2.8383E+01  8.6582E+01  9.9998E+01
 RELATIVEINF(%)  9.7083E+00  4.9459E-01  4.0730E+00  3.9224E-01  2.1677E-04
 EPSSHRINKSD(%)  1.7424E+01
 EPSSHRINKVR(%)  3.1813E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1183.2642297422842     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -264.32569653761152     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.28
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1183.264       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.59E-01  1.00E-02  1.10E-01  7.56E-01  6.00E+02  1.60E+00  7.58E+00  1.96E+00  4.16E-01  4.69E+01  5.27E+00
 


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
+        1.07E+05
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.37E+01  0.00E+00  1.02E+05
 
 TH 4
+       -3.29E+02  0.00E+00 -4.87E+03  2.51E+05
 
 TH 5
+        6.26E-06  0.00E+00 -4.50E-04  2.02E-04  1.35E-08
 
 TH 6
+        2.00E+01  0.00E+00  2.81E+01 -2.65E+01  1.28E+00  5.51E+03
 
 TH 7
+        1.50E+01  0.00E+00  3.22E+00 -2.36E+01 -9.05E-07  4.99E-02  1.80E+01
 
 TH 8
+       -1.52E+01  0.00E+00 -1.45E+02  3.23E+00 -2.06E-05  1.65E+00 -2.68E-01  2.00E+03
 
 TH 9
+       -1.53E+01  0.00E+00  5.88E+01 -3.80E+01  2.97E-05  1.86E-01 -2.13E-01  1.32E+01  4.41E+04
 
 TH10
+        9.37E-04  0.00E+00  1.24E-03 -1.29E-03 -2.34E-07  2.70E+01 -7.75E-07  1.52E-04 -6.08E-04  4.58E-06
 
 TH11
+       -8.10E+00  0.00E+00  4.22E+01 -1.48E+01 -1.48E-06  5.04E-01 -3.23E-01  2.49E+00  2.77E+00 -2.87E-05  7.29E+01
 
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
 #CPUT: Total CPU Time in Seconds,       25.616
Stop Time:
Thu Sep 30 02:50:12 CDT 2021
