Wed Sep 29 03:51:47 CDT 2021
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
$DATA ../../../../data/int/SL3/dat5.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      982
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      882
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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -932.514423625569        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3417E+02  1.4141E+02  1.3336E+02  7.7211E+01  9.8448E+01  4.1542E+01 -5.4148E+01 -1.8345E+02 -3.9964E+01  3.6157E+00
            -5.2986E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2696.29361289831        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.4562E-01  1.1338E+00  9.7491E-01  9.6480E-01  1.0648E+00  8.0485E-01  1.0052E+00  9.6206E-01  8.3968E-01  7.7782E-01
             3.0569E+00
 PARAMETER:  4.4082E-02  2.2558E-01  7.4593E-02  6.4170E-02  1.6283E-01 -1.1710E-01  1.0520E-01  6.1321E-02 -7.4730E-02 -1.5126E-01
             1.2174E+00
 GRADIENT:  -1.5264E+01  2.3863E+01 -1.1281E+01  1.6839E+01  6.5416E-01 -4.5909E+01 -5.2523E+00  8.5504E+00 -3.6269E+00  1.9774E+00
             3.7413E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2703.42370979292        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.2676E-01  8.8058E-01  7.0948E-01  1.1051E+00  7.9842E-01  8.6896E-01  1.3703E+00  7.0504E-01  8.5014E-01  3.8302E-01
             2.9820E+00
 PARAMETER:  2.3938E-02 -2.7175E-02 -2.4322E-01  1.9992E-01 -1.2512E-01 -4.0457E-02  4.1503E-01 -2.4950E-01 -6.2360E-02 -8.5968E-01
             1.1926E+00
 GRADIENT:  -5.8106E+01  1.7702E+01 -5.4552E+01  8.1365E+01  6.0728E+01 -1.9079E+01  1.1313E+01  1.0637E+01  3.5253E+00  2.1419E+00
             3.3755E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2738.26429226038        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      273
 NPARAMETR:  9.5556E-01  9.7057E-01  8.7872E-01  1.0172E+00  9.4177E-01  9.0605E-01  1.2772E+00  3.2434E-01  8.7228E-01  6.7460E-01
             2.5740E+00
 PARAMETER:  5.4540E-02  7.0133E-02 -2.9291E-02  1.1709E-01  4.0011E-02  1.3357E-03  3.4469E-01 -1.0260E+00 -3.6646E-02 -2.9363E-01
             1.0455E+00
 GRADIENT:   3.3382E+00 -1.7365E+01 -1.3512E+01 -2.2075E+01  2.1394E+01 -4.2286E+00  9.5115E+00  9.1532E-01 -1.0931E+00  3.4011E+00
             3.9828E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2739.82073196685        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  9.5393E-01  8.2235E-01  7.5388E-01  1.1034E+00  7.8229E-01  9.1713E-01  1.4053E+00  2.1608E-01  8.5906E-01  4.7492E-01
             2.5315E+00
 PARAMETER:  5.2839E-02 -9.5586E-02 -1.8252E-01  1.9836E-01 -1.4554E-01  1.3492E-02  4.4022E-01 -1.4321E+00 -5.1911E-02 -6.4460E-01
             1.0288E+00
 GRADIENT:   5.9606E-01  2.2732E+00 -1.0490E+00  8.5352E+00  2.3066E+00 -1.7794E-02  3.6766E-01  3.0129E-01  7.8374E-01  3.5785E-01
             8.2082E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2739.97206374199        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  9.5342E-01  7.8687E-01  7.1761E-01  1.1128E+00  7.4638E-01  9.1685E-01  1.4540E+00  9.8842E-02  8.5257E-01  4.1129E-01
             2.5317E+00
 PARAMETER:  5.2299E-02 -1.3969E-01 -2.3182E-01  2.0686E-01 -1.9253E-01  1.3183E-02  4.7434E-01 -2.2142E+00 -5.9497E-02 -7.8845E-01
             1.0289E+00
 GRADIENT:   6.0748E+01  1.3062E+01  6.6711E+00  3.5193E+01  1.3507E+01  4.4537E+00  6.5143E+00  6.9189E-02  1.0061E+00  6.6735E-01
             1.6264E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2739.99174586943        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  9.5327E-01  7.8570E-01  7.1538E-01  1.1114E+00  7.4635E-01  9.1653E-01  1.4504E+00  5.0671E-02  8.5197E-01  4.1965E-01
             2.5338E+00
 PARAMETER:  5.2145E-02 -1.4118E-01 -2.3495E-01  2.0561E-01 -1.9256E-01  1.2837E-02  4.7186E-01 -2.8824E+00 -6.0204E-02 -7.6833E-01
             1.0297E+00
 GRADIENT:  -7.9051E-01 -1.3232E+00 -1.1329E+00 -1.6547E+00  1.0971E+00 -1.7845E-01 -3.4245E-01  1.1657E-02 -8.3604E-02  4.0606E-02
             1.7923E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2739.99964519096        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  9.5348E-01  7.8638E-01  7.1690E-01  1.1123E+00  7.4686E-01  9.1700E-01  1.4533E+00  1.0000E-02  8.5206E-01  4.2038E-01
             2.5338E+00
 PARAMETER:  5.2360E-02 -1.4032E-01 -2.3282E-01  2.0639E-01 -1.9187E-01  1.3355E-02  4.7385E-01 -4.6896E+00 -6.0095E-02 -7.6660E-01
             1.0297E+00
 GRADIENT:  -2.6263E-01 -2.4206E-01  9.7963E-02  5.6983E-03 -3.2763E-01  4.1850E-03  6.3524E-02  0.0000E+00 -3.4402E-03 -1.7302E-02
            -6.4720E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2740.00236243507        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1098
 NPARAMETR:  9.5360E-01  7.9414E-01  7.2385E-01  1.1093E+00  7.5425E-01  9.1701E-01  1.4428E+00  1.0000E-02  8.5319E-01  4.3191E-01
             2.5337E+00
 PARAMETER:  5.2484E-02 -1.3049E-01 -2.2317E-01  2.0375E-01 -1.8203E-01  1.3358E-02  4.6656E-01 -1.0650E+01 -5.8772E-02 -7.3954E-01
             1.0297E+00
 GRADIENT:   1.1241E-03  1.0083E-02 -3.0352E-02 -4.8706E-02  2.5043E-02 -2.3075E-04  9.7021E-03  0.0000E+00  6.2037E-03 -4.1316E-03
             6.7246E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1098
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8383E-03  3.9688E-03 -2.3788E-04 -5.9416E-03  1.6164E-03
 SE:             2.9338E-02  2.5717E-02  2.0616E-04  2.6030E-02  1.4409E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5004E-01  8.7735E-01  2.4856E-01  8.1945E-01  9.1068E-01

 ETASHRINKSD(%)  1.7137E+00  1.3843E+01  9.9309E+01  1.2796E+01  5.1728E+01
 ETASHRINKVR(%)  3.3981E+00  2.5770E+01  9.9995E+01  2.3954E+01  7.6698E+01
 EBVSHRINKSD(%)  1.8250E+00  1.2750E+01  9.9316E+01  1.2770E+01  5.3277E+01
 EBVSHRINKVR(%)  3.6167E+00  2.3874E+01  9.9995E+01  2.3909E+01  7.8169E+01
 RELATIVEINF(%)  9.6336E+01  1.3336E+01  7.6203E-04  2.8192E+01  1.3914E+00
 EPSSHRINKSD(%)  1.5917E+01
 EPSSHRINKVR(%)  2.9300E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          882
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1621.0075725730426     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2740.0023624350729     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1118.9947898620303     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2740.002       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.54E-01  7.94E-01  7.24E-01  1.11E+00  7.54E-01  9.17E-01  1.44E+00  1.00E-02  8.53E-01  4.32E-01  2.53E+00
 


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
+        1.40E+03
 
 TH 2
+       -5.95E+00  6.50E+02
 
 TH 3
+        7.58E+00  2.97E+02  1.02E+03
 
 TH 4
+       -2.25E+01  1.43E+02 -2.68E+02  9.46E+02
 
 TH 5
+       -1.37E+01 -8.17E+02 -1.30E+03  5.12E+02  2.31E+03
 
 TH 6
+        7.53E+00 -1.60E+00  5.28E+00 -9.03E+00 -9.31E+00  2.19E+02
 
 TH 7
+        4.06E-01  2.60E+01 -3.89E+01  9.08E-01  3.72E+01 -2.26E-01  5.41E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.03E+00 -1.72E+01  1.76E+01  3.05E+01 -2.89E+01 -8.97E-02 -1.86E+00  0.00E+00  1.56E+02
 
 TH10
+        1.01E-01 -1.72E+01 -9.29E+01 -2.59E+01  3.99E+01  1.05E+00  3.05E+01  0.00E+00  1.09E+01  6.37E+01
 
 TH11
+       -1.83E+01 -1.35E+01 -6.62E+00 -1.92E+01  3.41E+00  2.95E+00  4.62E+00  0.00E+00  1.47E+01  1.04E+01  1.79E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.701
Stop Time:
Wed Sep 29 03:52:25 CDT 2021
