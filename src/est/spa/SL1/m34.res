Wed Sep 29 15:00:39 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat34.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1582.46872738134        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9430E+02 -8.4853E+01 -3.4384E+01 -2.2377E+01  1.1923E+02 -2.4555E+00 -5.8738E+00 -1.1058E+01  3.7412E+01 -3.2759E+01
            -3.6046E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1595.32546797283        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      170
 NPARAMETR:  9.8271E-01  1.1284E+00  9.4287E-01  9.9700E-01  9.3494E-01  1.2501E+00  1.0227E+00  1.1112E+00  7.6281E-01  1.1606E+00
             1.2091E+00
 PARAMETER:  8.2560E-02  2.2080E-01  4.1172E-02  9.6995E-02  3.2731E-02  3.2325E-01  1.2243E-01  2.0542E-01 -1.7075E-01  2.4894E-01
             2.8988E-01
 GRADIENT:   4.9418E+01  3.8006E+00 -1.0400E+00  1.3394E+01 -2.1689E+01  5.1586E+01 -8.7491E+00 -8.4789E-01  2.7920E+00  1.4659E+01
             4.4245E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1598.68409828796        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8563E-01  1.2706E+00  9.7580E-01  9.1443E-01  9.8688E-01  1.1504E+00  1.1895E+00  1.3352E+00  5.8544E-01  1.1661E+00
             1.1838E+00
 PARAMETER:  8.5530E-02  3.3946E-01  7.5502E-02  1.0550E-02  8.6789E-02  2.4007E-01  2.7355E-01  3.8908E-01 -4.3540E-01  2.5367E-01
             2.6874E-01
 GRADIENT:   6.3618E+01  3.7540E+01  1.3187E+01  4.1474E-01 -4.5936E+01  2.4926E+01  1.2763E+01 -3.4256E+00 -6.0791E+00  7.0855E+00
             3.3571E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1603.46109183435        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.7515E-01  1.2364E+00  1.5305E+00  9.4482E-01  1.1893E+00  1.1221E+00  1.2035E+00  2.0334E+00  5.6023E-01  1.3201E+00
             1.1666E+00
 PARAMETER:  7.4833E-02  3.1224E-01  5.2560E-01  4.3237E-02  2.7334E-01  2.1516E-01  2.8527E-01  8.0973E-01 -4.7941E-01  3.7769E-01
             2.5407E-01
 GRADIENT:   5.1704E+01  2.2864E+01  3.0969E-01  1.0590E+01 -1.0356E+01  1.7531E+01  8.9697E+00 -1.1231E+00 -3.8462E+00  2.5741E+00
             2.7339E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1604.06839316210        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  9.7443E-01  1.2354E+00  1.5284E+00  9.4439E-01  1.2046E+00  1.0719E+00  1.2009E+00  2.0292E+00  5.6180E-01  1.3207E+00
             1.1643E+00
 PARAMETER:  7.4096E-02  3.1142E-01  5.2424E-01  4.2781E-02  2.8612E-01  1.6945E-01  2.8307E-01  8.0764E-01 -4.7662E-01  3.7813E-01
             2.5215E-01
 GRADIENT:   5.5009E+01  1.9847E+01 -2.3699E+00  1.1439E+01  1.7642E-01  4.8469E-01  8.7286E+00 -1.2994E+00 -3.7470E+00  9.5015E-01
             2.6203E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1605.21315794240        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  9.7155E-01  1.2309E+00  1.5266E+00  9.4320E-01  1.2020E+00  1.0588E+00  1.0946E+00  2.0203E+00  6.7910E-01  1.3050E+00
             1.1585E+00
 PARAMETER:  7.1141E-02  3.0777E-01  5.2307E-01  4.1525E-02  2.8399E-01  1.5714E-01  1.9038E-01  8.0324E-01 -2.8699E-01  3.6622E-01
             2.4710E-01
 GRADIENT:   3.3203E+02  9.5225E+01 -4.7202E-01  4.4864E+01  1.4192E+01  3.8774E+01  6.3374E+00  2.8614E+00  4.4217E+00  3.8802E+00
             2.7006E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1606.37154627948        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1044            RESET HESSIAN, TYPE II
 NPARAMETR:  9.5936E-01  1.2285E+00  1.5417E+00  9.4173E-01  1.1982E+00  1.0702E+00  1.0854E+00  2.0133E+00  6.8815E-01  1.3010E+00
             1.1198E+00
 PARAMETER:  5.8508E-02  3.0582E-01  5.3291E-01  3.9964E-02  2.8084E-01  1.6785E-01  1.8193E-01  7.9979E-01 -2.7375E-01  3.6315E-01
             2.1319E-01
 GRADIENT:   3.2746E+02  9.8444E+01  2.4881E+00  4.1657E+01  1.1306E+01  5.1661E+01  5.0386E+00  2.0291E+00  4.3226E+00  2.9182E+00
             1.4024E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1606.53018101200        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1206
 NPARAMETR:  9.4356E-01  1.2181E+00  1.5371E+00  9.3956E-01  1.1959E+00  1.0674E+00  1.0826E+00  2.0060E+00  6.8447E-01  1.3008E+00
             1.1132E+00
 PARAMETER:  4.1902E-02  2.9730E-01  5.2990E-01  3.7652E-02  2.7893E-01  1.6523E-01  1.7936E-01  7.9614E-01 -2.7910E-01  3.6301E-01
             2.0728E-01
 GRADIENT:  -5.4414E+00 -1.0227E+01 -1.7369E-01 -1.4049E+01 -1.2950E+00  8.7652E-02 -2.1458E+00 -5.3524E-01 -8.8626E-01 -9.8274E-01
             1.0436E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1606.61305464067        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1386
 NPARAMETR:  9.4419E-01  1.2178E+00  1.5375E+00  9.4121E-01  1.1962E+00  1.0675E+00  1.0853E+00  2.0046E+00  6.8663E-01  1.3024E+00
             1.1073E+00
 PARAMETER:  4.2572E-02  2.9704E-01  5.3014E-01  3.9409E-02  2.7917E-01  1.6531E-01  1.8189E-01  7.9547E-01 -2.7596E-01  3.6419E-01
             2.0195E-01
 GRADIENT:  -2.0495E+04  1.3783E+04  7.7442E+03  4.0963E+04 -9.5809E-01 -1.9932E-01  2.2526E+04  5.0826E+03 -7.4269E+03 -5.6274E+03
            -2.0290E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1386
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0513E-03 -1.9276E-04 -4.9314E-02 -2.8449E-03 -4.6130E-02
 SE:             2.9810E-02  2.3469E-02  1.5912E-02  1.8083E-02  2.1425E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1847E-01  9.9345E-01  1.9411E-03  8.7499E-01  3.1314E-02

 ETASHRINKSD(%)  1.3163E-01  2.1375E+01  4.6693E+01  3.9420E+01  2.8222E+01
 ETASHRINKVR(%)  2.6309E-01  3.8181E+01  7.1583E+01  6.3300E+01  4.8479E+01
 EBVSHRINKSD(%)  4.9095E-01  2.1879E+01  5.4971E+01  4.1068E+01  2.2714E+01
 EBVSHRINKVR(%)  9.7949E-01  3.8972E+01  7.9724E+01  6.5271E+01  4.0269E+01
 RELATIVEINF(%)  9.8679E+01  1.7666E+00  3.0295E+00  9.5771E-01  2.1626E+01
 EPSSHRINKSD(%)  4.5381E+01
 EPSSHRINKVR(%)  7.0168E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1606.6130546406741     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -871.46222807693596     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.14
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1606.613       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.44E-01  1.22E+00  1.54E+00  9.41E-01  1.20E+00  1.07E+00  1.09E+00  2.00E+00  6.87E-01  1.30E+00  1.11E+00
 


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
+        1.15E+07
 
 TH 2
+        8.44E+01  7.83E+05
 
 TH 3
+       -4.91E+03  1.29E+03  1.55E+05
 
 TH 4
+       -1.15E+07 -4.52E+02  1.34E+06  1.16E+07
 
 TH 5
+       -7.05E+02  1.24E+02  2.55E+01  7.61E+02  2.86E+02
 
 TH 6
+       -4.16E+02  1.07E+02  4.98E+01  4.17E+02 -1.39E+00  1.71E+02
 
 TH 7
+        1.27E+03 -2.98E+02  6.38E+05 -1.29E+03  3.35E+02  2.01E+02  2.63E+06
 
 TH 8
+        2.12E+01 -1.62E+02  2.75E+02 -2.06E+02  3.91E+01  2.47E+01 -7.54E+01  3.92E+04
 
 TH 9
+        1.42E+03 -3.81E+02 -6.65E+05 -5.74E+06 -3.46E+02 -2.10E+02  6.71E+02 -7.75E+01  2.85E+06
 
 TH10
+        1.11E+02 -3.17E+01 -2.66E+05 -1.10E+02 -6.47E+05 -8.38E+01 -1.09E+06 -3.53E+00  5.65E+01  4.55E+05
 
 TH11
+       -1.61E+02 -6.67E+02 -2.08E+03  1.59E+03 -2.98E+02 -1.74E+02  5.44E+02 -1.46E+02  6.09E+02  5.57E+01  2.05E+06
 
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
 #CPUT: Total CPU Time in Seconds,       25.615
Stop Time:
Wed Sep 29 15:01:07 CDT 2021
