Fri Sep 24 21:25:25 CDT 2021
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
$DATA ../../../../data/int/A2/dat39.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2591.53749738032        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.9147E+01  1.3407E+02  1.1397E+02 -4.1611E+01  2.1878E+01  2.7132E+01 -8.6382E+01  1.0344E+01  1.9311E+01 -6.5323E+01
            -2.3848E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3196.07305343496        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8305E-01  8.6521E-01  8.7418E-01  1.1194E+00  8.5465E-01  9.4077E-01  9.9734E-01  5.1705E-01  7.4433E-01  7.8345E-01
             1.8288E+00
 PARAMETER:  8.2900E-02 -4.4787E-02 -3.4473E-02  2.1277E-01 -5.7065E-02  3.8947E-02  9.7338E-02 -5.5961E-01 -1.9526E-01 -1.4404E-01
             7.0367E-01
 GRADIENT:  -9.5310E+00  7.2948E+01  3.8341E+01  3.5114E+01 -5.4049E+01  4.2068E+00 -1.7345E+01  8.1557E+00 -5.0430E+01 -2.0206E+01
            -1.9291E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3214.18022598622        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.8965E-01  7.0278E-01  7.2313E-01  1.2205E+00  6.9727E-01  9.2321E-01  9.5065E-01  1.6403E-01  8.9157E-01  9.2735E-01
             1.7962E+00
 PARAMETER:  8.9593E-02 -2.5271E-01 -2.2416E-01  2.9927E-01 -2.6059E-01  2.0100E-02  4.9395E-02 -1.7077E+00 -1.4772E-02  2.4574E-02
             6.8569E-01
 GRADIENT:   1.0368E+01  6.6416E+01  4.6260E+01  1.2828E+02 -2.2702E+01 -2.1296E+00 -2.0697E+01  1.2981E+00 -1.3893E+01  7.0606E+00
             2.0748E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3236.91892921636        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      237
 NPARAMETR:  9.8963E-01  3.6876E-01  3.5531E-01  1.2652E+00  3.6293E-01  9.4680E-01  1.1315E+00  3.2285E-02  9.6449E-01  6.3401E-01
             1.7213E+00
 PARAMETER:  8.9578E-02 -8.9762E-01 -9.3478E-01  3.3520E-01 -9.1354E-01  4.5334E-02  2.2356E-01 -3.3332E+00  6.3843E-02 -3.5570E-01
             6.4305E-01
 GRADIENT:   1.4236E+01  1.9154E+01  2.6380E+01  1.4342E+02 -1.5997E+01  5.9722E+00 -1.4927E+01  2.4704E-02 -1.3783E+01 -3.1313E-03
             7.5990E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3237.15791927444        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  9.8745E-01  3.2801E-01  3.0669E-01  1.2404E+00  3.2151E-01  9.3817E-01  1.2049E+00  2.1200E-02  9.7614E-01  5.9186E-01
             1.7024E+00
 PARAMETER:  8.7370E-02 -1.0147E+00 -1.0819E+00  3.1547E-01 -1.0347E+00  3.6174E-02  2.8640E-01 -3.7538E+00  7.5850E-02 -4.2449E-01
             6.3202E-01
 GRADIENT:   9.0388E+00  1.6936E+01  3.1497E+01  1.3664E+02 -2.8324E+01  2.1394E+00 -9.5494E+00 -4.8112E-03 -1.8644E+01 -8.3831E+00
            -1.3021E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3241.73102113546        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      465
 NPARAMETR:  9.8739E-01  3.4714E-01  3.3008E-01  1.1859E+00  3.4599E-01  9.2800E-01  1.2667E+00  1.9435E-02  9.8533E-01  5.9755E-01
             1.6997E+00
 PARAMETER:  8.7315E-02 -9.5804E-01 -1.0084E+00  2.7046E-01 -9.6135E-01  2.5277E-02  3.3641E-01 -3.8407E+00  8.5222E-02 -4.1492E-01
             6.3045E-01
 GRADIENT:  -2.8894E+00 -1.0569E+00  2.4733E+00  1.3245E+00 -1.8239E+00 -2.4527E+00  1.3587E+00 -1.0407E-03 -4.1480E+00 -6.9292E+00
            -6.0626E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3241.90111755199        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  9.8853E-01  3.4971E-01  3.3225E-01  1.1860E+00  3.4805E-01  9.3376E-01  1.2492E+00  1.8581E-02  9.9534E-01  6.1892E-01
             1.7041E+00
 PARAMETER:  8.8462E-02 -9.5065E-01 -1.0019E+00  2.7055E-01 -9.5542E-01  3.1462E-02  3.2252E-01 -3.8856E+00  9.5328E-02 -3.7977E-01
             6.3304E-01
 GRADIENT:  -2.0856E-03 -6.1218E-03  7.8065E-03 -2.1920E-02 -1.0569E-02 -4.4925E-03  2.9202E-03  2.2747E-03 -8.1891E-03  3.5434E-03
            -2.8025E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3241.90191882122        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      738
 NPARAMETR:  9.8854E-01  3.4969E-01  3.3222E-01  1.1860E+00  3.4803E-01  9.3378E-01  1.2492E+00  1.0000E-02  9.9540E-01  6.1895E-01
             1.7042E+00
 PARAMETER:  8.8472E-02 -9.5070E-01 -1.0020E+00  2.7055E-01 -9.5546E-01  3.1483E-02  3.2248E-01 -5.6669E+00  9.5391E-02 -3.7973E-01
             6.3307E-01
 GRADIENT:   1.3062E+01  9.7035E+00  8.1543E+00  1.1776E+01  3.5589E+01  9.9844E-01  4.7287E-01  0.0000E+00  5.4327E-01  5.1777E-01
             1.1483E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3241.90192265316        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  9.8853E-01  3.4970E-01  3.3223E-01  1.1860E+00  3.4803E-01  9.3377E-01  1.2492E+00  1.0000E-02  9.9537E-01  6.1898E-01
             1.7042E+00
 PARAMETER:  8.8460E-02 -9.5069E-01 -1.0019E+00  2.7055E-01 -9.5546E-01  3.1473E-02  3.2250E-01 -5.5779E+00  9.5357E-02 -3.7968E-01
             6.3307E-01
 GRADIENT:  -6.1531E-03 -1.3731E-03 -1.8704E-03 -3.5881E-03 -3.3840E-02  1.2446E-04 -3.0963E-04  0.0000E+00  1.4075E-03 -1.2950E-03
             3.7869E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -3241.90192265316        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:      952
 NPARAMETR:  9.8863E-01  3.4980E-01  3.3225E-01  1.1859E+00  3.4803E-01  9.3368E-01  1.2493E+00  1.0000E-02  9.9527E-01  6.1912E-01
             1.7041E+00
 PARAMETER:  8.8460E-02 -9.5069E-01 -1.0019E+00  2.7055E-01 -9.5546E-01  3.1473E-02  3.2250E-01 -5.5779E+00  9.5357E-02 -3.7968E-01
             6.3307E-01
 GRADIENT:  -9.0146E-03 -1.0951E-02 -4.7043E-03  3.2163E-03 -1.4698E-03  1.3967E-03 -8.2870E-04  0.0000E+00  1.4767E-03 -1.6583E-03
             4.2377E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      952
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4549E-04  2.0841E-03 -1.2419E-04 -2.7421E-03  1.3903E-04
 SE:             2.9668E-02  2.4693E-02  3.3456E-04  2.9049E-02  2.5277E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7726E-01  9.3274E-01  7.1048E-01  9.2480E-01  9.9561E-01

 ETASHRINKSD(%)  6.0814E-01  1.7274E+01  9.8879E+01  2.6809E+00  1.5319E+01
 ETASHRINKVR(%)  1.2126E+00  3.1565E+01  9.9987E+01  5.2900E+00  2.8291E+01
 EBVSHRINKSD(%)  8.4296E-01  1.5881E+01  9.8883E+01  2.6117E+00  1.5988E+01
 EBVSHRINKVR(%)  1.6788E+00  2.9240E+01  9.9988E+01  5.1551E+00  2.9421E+01
 RELATIVEINF(%)  9.8308E+01  1.7387E+01  9.6375E-04  8.2221E+01  4.9679E+00
 EPSSHRINKSD(%)  2.0100E+01
 EPSSHRINKVR(%)  3.6160E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3241.9019226531559     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1587.8125628847451     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3241.902       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.89E-01  3.50E-01  3.32E-01  1.19E+00  3.48E-01  9.34E-01  1.25E+00  1.00E-02  9.95E-01  6.19E-01  1.70E+00
 


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
+        1.29E+03
 
 TH 2
+       -2.31E+00  4.16E+03
 
 TH 3
+        1.48E+01 -9.75E+02  1.16E+04
 
 TH 4
+       -5.74E+00  2.67E+01 -3.27E+02  7.57E+02
 
 TH 5
+       -6.92E+00 -2.81E+03 -1.21E+04 -2.96E+01  1.70E+04
 
 TH 6
+       -1.82E+00 -4.03E+00  1.48E+01 -3.07E+00 -8.69E+00  2.25E+02
 
 TH 7
+       -1.03E+00  3.98E+01  2.07E+01  1.21E+00 -3.58E+01 -4.59E-02  6.28E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.68E+00 -1.65E+01  1.18E+02 -3.84E+00 -1.52E+01  1.32E+00  8.42E-01  0.00E+00  1.90E+02
 
 TH10
+       -2.08E+00 -1.58E+01 -8.09E+01  8.97E+00  9.87E+00 -6.40E-01  1.93E+01  0.00E+00  3.10E+00  2.79E+02
 
 TH11
+       -1.35E+01 -1.59E+01 -1.35E+02 -9.47E+00  7.88E+01  1.76E+00  1.05E+01  0.00E+00  4.78E+00  1.80E+01  3.73E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.042
Stop Time:
Fri Sep 24 21:25:59 CDT 2021
