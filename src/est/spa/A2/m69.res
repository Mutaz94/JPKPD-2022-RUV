Sat Sep 18 10:00:42 CDT 2021
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
$DATA ../../../../data/spa/A2/dat69.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -741.927381443108        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3463E-01  8.7700E+01  2.9539E+01  8.6991E+01  1.8345E+02  2.6632E-01 -5.2112E+01 -1.8880E+01 -9.8317E+01 -1.5632E+02
            -1.5499E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1342.85865590355        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0296E+00  8.4001E-01  9.4378E-01  1.0575E+00  8.2816E-01  9.2272E-01  1.0594E+00  9.2584E-01  1.1809E+00  1.1549E+00
             2.9196E+00
 PARAMETER:  1.2916E-01 -7.4340E-02  4.2139E-02  1.5592E-01 -8.8546E-02  1.9573E-02  1.5767E-01  2.2948E-02  2.6631E-01  2.4401E-01
             1.1714E+00
 GRADIENT:  -1.8159E+01  3.6013E+00 -1.4727E+00 -3.6362E+00  3.9516E+00 -1.6470E+01  1.0387E+00  7.1491E+00  9.4387E+00  5.5040E+00
             2.1739E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1350.78266542144        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0199E+00  4.9278E-01  6.5383E-01  1.3292E+00  5.6214E-01  9.6011E-01  1.1672E+00  3.2982E-01  1.1454E+00  8.5789E-01
             2.8146E+00
 PARAMETER:  1.1970E-01 -6.0770E-01 -3.2491E-01  3.8459E-01 -4.7600E-01  5.9296E-02  2.5462E-01 -1.0092E+00  2.3576E-01 -5.3276E-02
             1.1348E+00
 GRADIENT:  -4.9043E+01  3.0392E+01 -2.9748E+00  1.2176E+02  4.0890E+00 -7.8377E+00 -1.1511E+00  1.7660E+00  1.7264E+01  1.8598E-01
             9.4100E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1366.62802436786        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0435E+00  4.0089E-01  2.9346E-01  1.1078E+00  3.2258E-01  1.0128E+00  1.3900E+00  2.2259E-02  1.0209E+00  5.6561E-01
             2.6354E+00
 PARAMETER:  1.4255E-01 -8.1406E-01 -1.1260E+00  2.0236E-01 -1.0314E+00  1.1275E-01  4.2932E-01 -3.7050E+00  1.2065E-01 -4.6985E-01
             1.0690E+00
 GRADIENT:  -3.4300E+00  2.9150E+00 -3.1370E+00  1.2673E+01  1.2098E+01  7.1799E+00 -2.9830E+00  1.6298E-03 -6.1259E+00 -1.8067E+00
            -1.0479E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1366.74557059743        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0457E+00  3.6494E-01  2.6436E-01  1.0865E+00  2.9608E-01  9.9939E-01  1.5061E+00  1.3049E-02  1.0432E+00  5.4989E-01
             2.6459E+00
 PARAMETER:  1.4467E-01 -9.0801E-01 -1.2304E+00  1.8294E-01 -1.1171E+00  9.9390E-02  5.0950E-01 -4.2390E+00  1.4230E-01 -4.9804E-01
             1.0730E+00
 GRADIENT:  -8.2977E-01  1.5675E+00 -1.7375E+00  5.0693E+00  3.0637E+00  1.7197E+00 -6.7043E-01  6.5710E-05 -2.4827E+00 -1.2819E+00
            -4.5910E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1367.80653677919        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0409E+00  2.5478E-01  3.1100E-01  1.1624E+00  3.0778E-01  9.9532E-01  2.0799E+00  1.0000E-02  1.0100E+00  6.0798E-01
             2.6568E+00
 PARAMETER:  1.4009E-01 -1.2674E+00 -1.0680E+00  2.5052E-01 -1.0784E+00  9.5305E-02  8.3231E-01 -4.8401E+00  1.0992E-01 -3.9762E-01
             1.0771E+00
 GRADIENT:  -2.8875E+00 -1.5408E+00 -7.1748E+00 -4.0370E+00  1.3877E+01  3.4171E+00 -2.9739E-01  0.0000E+00 -7.4244E-03 -9.1143E-01
            -5.8652E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1368.52847873038        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      509
 NPARAMETR:  1.0427E+00  2.1643E-01  3.6630E-01  1.2302E+00  3.3659E-01  9.7968E-01  2.5953E+00  1.0000E-02  9.7627E-01  6.3388E-01
             2.7073E+00
 PARAMETER:  1.4177E-01 -1.4305E+00 -9.0431E-01  3.0719E-01 -9.8889E-01  7.9474E-02  1.0537E+00 -4.9213E+00  7.5982E-02 -3.5590E-01
             1.0959E+00
 GRADIENT:   2.6525E+00  2.1906E+01 -2.3622E+01 -9.1227E+00  2.1187E+01  6.2353E-01  2.7395E+01  0.0000E+00 -9.9096E+00 -7.3266E+00
            -9.8069E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1368.63711034708        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      684
 NPARAMETR:  1.0414E+00  2.0237E-01  3.6582E-01  1.2327E+00  3.3538E-01  9.8004E-01  2.7140E+00  1.0000E-02  9.7494E-01  6.3281E-01
             2.7064E+00
 PARAMETER:  1.4056E-01 -1.4977E+00 -9.0562E-01  3.0924E-01 -9.9250E-01  7.9833E-02  1.0984E+00 -5.1144E+00  7.4616E-02 -3.5758E-01
             1.0956E+00
 GRADIENT:  -3.0240E-01 -1.1699E+00  1.7640E+00  9.4408E-01 -1.9678E+00 -1.3712E-01 -1.5725E+00  0.0000E+00  7.0062E-01  4.7425E-01
             4.6396E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1368.65629494548        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  1.0407E+00  1.9029E-01  3.6645E-01  1.2371E+00  3.3490E-01  9.7961E-01  2.8997E+00  1.0000E-02  9.6888E-01  6.2570E-01
             2.7098E+00
 PARAMETER:  1.3985E-01 -1.5592E+00 -9.0389E-01  3.1273E-01 -9.9393E-01  7.9399E-02  1.1646E+00 -5.2779E+00  6.8387E-02 -3.6889E-01
             1.0969E+00
 GRADIENT:  -4.0645E-03  1.1094E-02 -5.6945E-03  1.3768E-03  5.2875E-04 -9.0258E-04  3.4086E-03  0.0000E+00 -6.1282E-03 -6.6012E-03
            -1.5522E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1368.65629719619        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      919
 NPARAMETR:  1.0406E+00  1.9008E-01  3.6649E-01  1.2372E+00  3.3491E-01  9.7960E-01  2.9030E+00  1.0000E-02  9.6878E-01  6.2564E-01
             2.7099E+00
 PARAMETER:  1.3984E-01 -1.5603E+00 -9.0377E-01  3.1281E-01 -9.9391E-01  7.9386E-02  1.1657E+00 -5.2808E+00  6.8286E-02 -3.6898E-01
             1.0969E+00
 GRADIENT:   2.9436E-03  1.5262E-03  8.2409E-03 -1.9462E-03 -8.8928E-03 -1.6353E-03  3.1714E-03  0.0000E+00 -4.7123E-04 -4.7131E-03
             1.4845E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      919
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.5565E-04  3.1490E-02 -1.3656E-04 -1.8618E-02  8.1652E-03
 SE:             2.9019E-02  1.4378E-02  2.0846E-04  2.6191E-02  1.8551E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8472E-01  2.8513E-02  5.1242E-01  4.7716E-01  6.5983E-01

 ETASHRINKSD(%)  2.7817E+00  5.1832E+01  9.9302E+01  1.2258E+01  3.7851E+01
 ETASHRINKVR(%)  5.4861E+00  7.6799E+01  9.9995E+01  2.3014E+01  6.1376E+01
 EBVSHRINKSD(%)  2.7844E+00  6.2849E+01  9.9236E+01  1.0814E+01  3.3774E+01
 EBVSHRINKVR(%)  5.4913E+00  8.6198E+01  9.9994E+01  2.0458E+01  5.6141E+01
 RELATIVEINF(%)  9.2474E+01  5.0599E+00  2.3287E-04  2.5606E+01  1.6923E+00
 EPSSHRINKSD(%)  3.3569E+01
 EPSSHRINKVR(%)  5.5869E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1368.6562971961860     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -633.50547063244778     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.55
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
 





 #OBJV:********************************************    -1368.656       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.90E-01  3.66E-01  1.24E+00  3.35E-01  9.80E-01  2.90E+00  1.00E-02  9.69E-01  6.26E-01  2.71E+00
 


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
+        1.01E+03
 
 TH 2
+        2.54E+01  5.10E+03
 
 TH 3
+       -4.14E+01 -1.32E+03  5.85E+03
 
 TH 4
+       -4.88E+01 -1.60E+02 -2.63E+02  6.56E+02
 
 TH 5
+        1.65E+02  5.32E+02 -7.85E+03 -2.30E+02  1.22E+04
 
 TH 6
+       -9.74E-01  4.10E+01 -9.63E+00 -1.38E+01  2.71E+01  1.89E+02
 
 TH 7
+        7.60E+00  3.22E+02 -1.24E+02 -2.53E+01  1.08E+02  3.30E+00  2.18E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.14E+00 -3.82E+02  1.79E+02  2.17E+01 -1.22E+02 -6.49E+00 -2.11E+01  0.00E+00  1.68E+02
 
 TH10
+       -1.36E+01 -2.77E+02 -3.10E+01  1.72E+01  2.95E+01 -2.99E+00 -1.55E+01  0.00E+00  1.71E+01  1.23E+02
 
 TH11
+       -1.57E+01 -1.03E+02  2.40E+01 -9.10E-01 -1.50E+01  2.32E+00 -6.42E+00  0.00E+00  1.55E+01  2.64E+01  4.15E+01
 
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
 #CPUT: Total CPU Time in Seconds,       17.745
Stop Time:
Sat Sep 18 10:01:02 CDT 2021
