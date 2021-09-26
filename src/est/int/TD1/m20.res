Sat Sep 25 03:51:30 CDT 2021
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
$DATA ../../../../data/int/TD1/dat20.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3256.13533424710        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1067E+02 -2.2059E+01  1.0887E+02  9.1726E+00  1.0167E+02  1.7498E+00 -1.6054E+01 -2.5267E+02 -4.1147E+01 -1.2046E+01
            -8.7440E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3682.09951735662        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  7.8194E-01  9.6190E-01  9.3042E-01  9.9937E-01  9.1458E-01  1.3747E+00  1.1242E+00  9.0403E-01  9.5613E-01  8.8383E-01
             1.4714E+00
 PARAMETER: -1.4598E-01  6.1155E-02  2.7885E-02  9.9371E-02  1.0705E-02  4.1825E-01  2.1704E-01 -8.9181E-04  5.5138E-02 -2.3491E-02
             4.8621E-01
 GRADIENT:  -2.2524E+02 -2.4455E+01  5.9300E-01 -1.2291E+00  1.0937E+01  6.3674E+01  1.2684E+01  7.2608E+00 -6.3737E+00  5.4142E+00
             5.6687E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3767.26018605042        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  9.4540E-01  1.0013E+00  9.4433E-01  9.5976E-01  9.6168E-01  7.6243E-01  7.1657E-01  3.7280E-02  1.0491E+00  1.0954E+00
             1.2344E+00
 PARAMETER:  4.3856E-02  1.0134E-01  4.2717E-02  5.8928E-02  6.0923E-02 -1.7124E-01 -2.3328E-01 -3.1893E+00  1.4796E-01  1.9114E-01
             3.1058E-01
 GRADIENT:  -7.7359E+01 -6.1394E+01  5.4410E+01 -2.9310E+01  3.3552E+01 -1.4073E+02 -2.5754E+01 -9.7920E-02 -6.9237E-01  1.3291E+01
             2.4114E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3819.93038293014        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      239
 NPARAMETR:  9.2908E-01  8.9042E-01  8.1459E-01  1.0132E+00  8.3984E-01  9.5282E-01  9.9015E-01  3.4415E-01  1.0479E+00  9.7269E-01
             1.0359E+00
 PARAMETER:  2.6440E-02 -1.6065E-02 -1.0506E-01  1.1313E-01 -7.4543E-02  5.1672E-02  9.0104E-02 -9.6667E-01  1.4675E-01  7.2313E-02
             1.3531E-01
 GRADIENT:  -6.3496E+01 -3.7327E+01 -2.1154E+01 -8.3996E+00  2.9618E+01 -2.2244E+01 -1.8443E+01 -1.2111E+01  6.9834E+00  1.7682E+01
            -2.5256E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3827.32081244783        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  9.6270E-01  9.1032E-01  8.2093E-01  1.0216E+00  8.2932E-01  1.0084E+00  1.1245E+00  3.8183E-01  1.0280E+00  8.7756E-01
             1.0433E+00
 PARAMETER:  6.1982E-02  6.0388E-03 -9.7318E-02  1.2142E-01 -8.7151E-02  1.0839E-01  2.1736E-01 -8.6278E-01  1.2762E-01 -3.0612E-02
             1.4239E-01
 GRADIENT:   2.4406E+01  1.1292E+01  3.2510E+00  1.2433E+01 -1.0424E+01  6.0081E+00  1.7188E+00 -1.3901E+01  2.3217E+00  1.3664E+00
            -4.3992E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3827.52478698497        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      557
 NPARAMETR:  9.7102E-01  9.1953E-01  8.2954E-01  1.0175E+00  8.4460E-01  1.0154E+00  1.1170E+00  3.8183E-01  1.0275E+00  8.8215E-01
             1.0441E+00
 PARAMETER:  7.0591E-02  1.6109E-02 -8.6888E-02  1.1732E-01 -6.8894E-02  1.1529E-01  2.1064E-01 -8.6278E-01  1.2710E-01 -2.5396E-02
             1.4317E-01
 GRADIENT:  -1.7074E-01  2.6790E-01 -1.3411E+00  2.0501E+00  6.0030E-01  3.5069E-01 -5.0237E-01 -1.4108E+01  3.8261E-01 -1.6672E-01
            -9.9508E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3835.65785703156        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      717
 NPARAMETR:  9.7091E-01  9.2077E-01  8.3128E-01  1.0157E+00  8.4594E-01  1.0142E+00  1.1195E+00  6.1380E-01  1.0261E+00  8.8376E-01
             1.0446E+00
 PARAMETER:  7.0477E-02  1.7456E-02 -8.4786E-02  1.1562E-01 -6.7309E-02  1.1413E-01  2.1287E-01 -3.8808E-01  1.2580E-01 -2.3573E-02
             1.4366E-01
 GRADIENT:   4.2488E+01  5.8028E+00 -9.7706E+00  9.2872E+00  1.8284E+00  8.6664E+00  3.8608E+00 -1.7828E+01  4.0462E+00  6.5081E+00
             4.0510E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3837.46597045008        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  9.7050E-01  9.2077E-01  8.3128E-01  1.0157E+00  8.4606E-01  1.0138E+00  1.1187E+00  6.7831E-01  1.0261E+00  8.8185E-01
             1.0431E+00
 PARAMETER:  7.0054E-02  1.7457E-02 -8.4786E-02  1.1559E-01 -6.7170E-02  1.1369E-01  2.1219E-01 -2.8814E-01  1.2580E-01 -2.5731E-02
             1.4217E-01
 GRADIENT:  -1.5937E+00 -1.9979E+00 -1.4168E+01 -2.9883E+00 -4.7445E+00 -3.3446E-01  2.4898E+00 -1.6439E+01  2.9682E+00  6.8747E+00
             4.6585E+01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -3837.47134305528        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      965
 NPARAMETR:  9.7050E-01  9.2077E-01  8.3129E-01  1.0157E+00  8.4606E-01  1.0138E+00  1.1187E+00  6.7851E-01  1.0261E+00  8.8184E-01
             1.0431E+00
 PARAMETER:  7.0054E-02  1.7457E-02 -8.4781E-02  1.1559E-01 -6.7169E-02  1.1369E-01  2.1219E-01 -2.8786E-01  1.2579E-01 -2.5746E-02
             1.4216E-01
 GRADIENT:  -1.5583E+06 -1.5583E+06 -1.5582E+06 -1.3481E+06 -1.5583E+06 -1.3707E+06 -1.4688E+06  1.0825E+06  1.2387E+06  1.5583E+06
            -2.1926E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      965
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0771E-03 -2.1137E-02 -1.1129E-02  1.2596E-02 -2.3004E-02
 SE:             2.9911E-02  2.3782E-02  1.7668E-02  2.7844E-02  2.3230E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7127E-01  3.7411E-01  5.2877E-01  6.5099E-01  3.2205E-01

 ETASHRINKSD(%)  1.0000E-10  2.0327E+01  4.0809E+01  6.7194E+00  2.2177E+01
 ETASHRINKVR(%)  1.0000E-10  3.6522E+01  6.4964E+01  1.2987E+01  3.9436E+01
 EBVSHRINKSD(%)  2.7808E-01  1.9080E+01  4.9126E+01  6.4795E+00  2.0179E+01
 EBVSHRINKVR(%)  5.5538E-01  3.4520E+01  7.4119E+01  1.2539E+01  3.6287E+01
 RELATIVEINF(%)  9.9442E+01  3.3056E+01  1.4169E+01  6.5514E+01  2.3446E+01
 EPSSHRINKSD(%)  2.2790E+01
 EPSSHRINKVR(%)  4.0386E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3837.4713430552792     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2183.3819832868685     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.54
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3837.471       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.70E-01  9.21E-01  8.31E-01  1.02E+00  8.46E-01  1.01E+00  1.12E+00  6.79E-01  1.03E+00  8.82E-01  1.04E+00
 


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
+        8.27E+09
 
 TH 2
+       -2.08E+04  9.19E+09
 
 TH 3
+       -3.69E+05 -3.89E+05  1.13E+10
 
 TH 4
+       -1.63E+04  1.92E+05 -3.05E+05  5.65E+09
 
 TH 5
+       -4.95E+04 -5.26E+04 -5.84E+04 -4.07E+04  1.09E+10
 
 TH 6
+       -3.10E+04 -3.26E+04 -5.91E+08 -2.56E+04 -4.17E+04  5.87E+09
 
 TH 7
+       -6.50E+01 -5.03E+01 -1.51E+05 -4.76E+01 -2.02E+04 -1.27E+04  1.38E+09
 
 TH 8
+        9.79E+03 -4.33E+09  1.83E+05 -3.40E+09  2.46E+04  1.54E+04  3.23E+01  2.04E+09
 
 TH 9
+       -6.22E+09  3.87E+05  2.77E+05  3.04E+05  7.13E+09  2.33E+04 -2.54E+09  3.09E+09  1.43E+02
 
 TH10
+        3.40E+04  3.58E+04  3.97E+04  2.81E+04  5.44E+04  2.86E+04  1.39E+04 -1.69E+04 -2.56E+04  1.00E+10
 
 TH11
+        5.42E+09  2.11E+06  6.26E+09  1.65E+06 -6.21E+09 -4.56E+09  2.21E+09 -2.69E+09  4.07E+09  5.96E+09  8.98E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.371
Stop Time:
Sat Sep 25 03:52:09 CDT 2021
