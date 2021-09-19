Sat Sep 18 07:08:36 CDT 2021
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
$DATA ../../../../data/int/D/dat54.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   37330.6454068572        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.4379E+02  6.0186E+02 -4.6252E+01  5.0983E+02  6.3177E+01 -2.7242E+03 -1.2930E+03 -2.5614E+01 -1.9077E+03 -7.1898E+02
            -7.5327E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -810.785176074940        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.5801E+00  1.7855E+00  8.3048E-01  1.8744E+00  1.0456E+00  5.0670E+00  4.8848E+00  9.6696E-01  2.5506E+00  1.3065E+00
             1.2743E+01
 PARAMETER:  5.5752E-01  6.7973E-01 -8.5756E-02  7.2831E-01  1.4455E-01  1.7227E+00  1.6861E+00  6.6401E-02  1.0363E+00  3.6737E-01
             2.6449E+00
 GRADIENT:   6.8801E+00  1.4542E+01 -6.4529E+01  1.1650E+02  1.8496E+01  1.4899E+02  6.7372E+01  3.3825E+00  3.5083E+01  2.6325E+01
             3.0667E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -867.322009141230        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3870E+00  3.7978E+00  7.9909E+01  3.1335E+00  2.6797E+00  2.4030E+00  3.8920E+01  5.6063E-01  2.3373E+00  1.0086E+00
             1.3099E+01
 PARAMETER:  4.2716E-01  1.4344E+00  4.4809E+00  1.2421E+00  1.0857E+00  9.7670E-01  3.7615E+00 -4.7869E-01  9.4898E-01  1.0857E-01
             2.6725E+00
 GRADIENT:   6.8896E+00  1.0981E+01 -4.2185E+00  4.8480E+00  2.7286E+01  5.7671E+01  4.6170E+01  2.1435E-04  3.5427E+01  1.5842E+01
             3.6587E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1038.24835286758        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1048E+00  1.1237E+00  6.3662E+00  1.0448E+00  1.9698E+00  2.1311E+00  5.8626E+00  7.8371E-01  7.0047E-01  9.1957E-02
             1.0925E+01
 PARAMETER:  1.9966E-01  2.1667E-01  1.9510E+00  1.4387E-01  7.7793E-01  8.5666E-01  1.8686E+00 -1.4372E-01 -2.5601E-01 -2.2864E+00
             2.4911E+00
 GRADIENT:  -4.8364E+01 -7.1800E+00 -3.0452E+00 -8.3624E+01 -2.5339E+00  2.0162E+01  5.2220E+01  5.1051E-01  1.2113E+01  1.8859E-01
             2.4399E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1061.56949987131        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.1950E+00  8.1604E-01  1.4207E+01  1.3360E+00  2.3132E+00  1.8507E+00  5.9126E+00  7.4043E-01  8.5688E-01  9.9052E-02
             1.0128E+01
 PARAMETER:  2.7812E-01 -1.0329E-01  2.7538E+00  3.8965E-01  9.3861E-01  7.1559E-01  1.8771E+00 -2.0053E-01 -5.4458E-02 -2.2121E+00
             2.4153E+00
 GRADIENT:  -2.2671E+00 -1.3921E+00 -3.7415E+00  8.1355E+00  1.8432E+01 -1.3868E+01  4.2316E+00  1.5784E-01  5.3437E+00  2.0167E-01
             5.6373E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1061.65493100400        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      395
 NPARAMETR:  1.1946E+00  8.1564E-01  1.4190E+01  1.3338E+00  2.3111E+00  1.8524E+00  5.9093E+00  7.4067E-01  8.5381E-01  9.9061E-02
             1.0121E+01
 PARAMETER:  2.7777E-01 -1.0379E-01  2.7525E+00  3.8800E-01  9.3774E-01  7.1646E-01  1.8765E+00 -2.0020E-01 -5.8048E-02 -2.2120E+00
             2.4146E+00
 GRADIENT:  -2.2385E+00 -1.5015E+00 -3.6957E+00  7.7912E+00  1.8132E+01 -1.3520E+01  4.2100E+00  1.5805E-01  5.3005E+00  2.0180E-01
             5.5396E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1062.79736488038        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  1.1944E+00  8.2343E-01  2.4951E+01  1.3340E+00  2.3124E+00  1.8517E+00  5.9157E+00  7.4005E-01  7.8138E-01  1.0000E-02
             1.0124E+01
 PARAMETER:  2.7762E-01 -9.4272E-02  3.3169E+00  3.8821E-01  9.3828E-01  7.1608E-01  1.8776E+00 -2.0103E-01 -1.4669E-01 -7.0189E+00
             2.4149E+00
 GRADIENT:  -1.8308E+00  3.2110E-01 -1.9361E-02  1.9244E+01 -5.0768E+00 -1.2859E+01  8.9305E-01  4.8934E-02  9.6765E-02  0.0000E+00
             4.9866E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1063.02704830820        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      599
 NPARAMETR:  1.1944E+00  8.1109E-01  2.1255E+01  1.3339E+00  2.3126E+00  1.8520E+00  5.9147E+00  7.3659E-01  7.5192E-01  1.0000E-02
             1.0045E+01
 PARAMETER:  2.7763E-01 -1.0937E-01  3.1566E+00  3.8813E-01  9.3839E-01  7.1627E-01  1.8774E+00 -2.0572E-01 -1.8512E-01 -6.8336E+00
             2.4071E+00
 GRADIENT:  -3.9919E-01  4.4437E-02 -4.7986E-01  2.5493E+01 -1.0816E+00 -1.2806E+01 -6.1877E-01  6.9024E-02 -1.5809E+00  0.0000E+00
             3.3056E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1063.04806108505        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  1.1943E+00  7.9980E-01  2.4856E+01  1.3341E+00  2.3132E+00  1.8517E+00  5.9179E+00  7.3381E-01  7.5266E-01  1.0000E-02
             1.0049E+01
 PARAMETER:  2.7755E-01 -1.2340E-01  3.3131E+00  3.8823E-01  9.3864E-01  7.1609E-01  1.8780E+00 -2.0951E-01 -1.8414E-01 -6.8336E+00
             2.4075E+00
 GRADIENT:  -2.8195E-01 -4.5166E-01  1.5989E-02  2.4298E+01 -5.5917E+00 -1.2757E+01 -1.2233E+00  4.9271E-02 -1.7704E+00  0.0000E+00
             3.4548E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1063.07533149258        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      914
 NPARAMETR:  1.1943E+00  8.1438E-01  2.4938E+01  1.3341E+00  2.3135E+00  1.8516E+00  5.9190E+00  7.2394E-01  7.5265E-01  1.0000E-02
             1.0042E+01
 PARAMETER:  2.7754E-01 -1.0531E-01  3.3197E+00  3.8823E-01  9.3872E-01  7.1609E-01  1.8781E+00 -2.2289E-01 -1.8414E-01 -6.8336E+00
             2.4067E+00
 GRADIENT:   2.5829E+03  3.7370E-01  4.5237E-02 -1.8262E+03 -3.8849E+02  9.8189E+02 -2.0187E+02  5.6251E-02  3.9008E+03  0.0000E+00
            -2.6620E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      914
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0327E-02  3.3903E-02 -2.3782E-04 -7.2818E-02  6.8853E-05
 SE:             2.9649E-02  2.4858E-02  3.8042E-04  1.2590E-02  1.0305E-04
 N:                     100         100         100         100         100

 P VAL.:         4.9297E-01  1.7262E-01  5.3186E-01  7.3221E-09  5.0403E-01

 ETASHRINKSD(%)  6.7368E-01  1.6722E+01  9.8726E+01  5.7822E+01  9.9655E+01
 ETASHRINKVR(%)  1.3428E+00  3.0648E+01  9.9984E+01  8.2210E+01  9.9999E+01
 EBVSHRINKSD(%)  6.1430E+00  1.4199E+01  9.7935E+01  6.2200E+01  9.9536E+01
 EBVSHRINKVR(%)  1.1909E+01  2.6381E+01  9.9957E+01  8.5712E+01  9.9998E+01
 RELATIVEINF(%)  8.7697E+01  3.8602E+01  8.2102E-03  7.5231E+00  4.0513E-04
 EPSSHRINKSD(%)  5.7877E+00
 EPSSHRINKVR(%)  1.1240E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1063.0753314925846     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       591.01402827582615     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1063.075       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  8.14E-01  2.50E+01  1.33E+00  2.31E+00  1.85E+00  5.92E+00  7.24E-01  7.53E-01  1.00E-02  1.00E+01
 


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
+        1.63E+06
 
 TH 2
+        1.22E+02 -1.78E+00
 
 TH 3
+       -7.31E-01  4.20E-02  2.01E-03
 
 TH 4
+       -3.17E+01 -4.96E+01  5.30E-01  6.70E+05
 
 TH 5
+        1.00E+03 -2.75E+01 -1.44E-01 -1.68E+01  3.80E+04
 
 TH 6
+       -1.36E+03  3.75E+01 -1.83E-01 -1.66E+00  2.49E+02  1.01E+05
 
 TH 7
+        1.64E+02  7.76E-02  1.58E-02 -1.06E+01 -2.91E+01  3.01E+01  1.45E+03
 
 TH 8
+       -2.29E+00  3.24E+00 -2.79E-02 -5.34E+00  7.47E-01 -1.51E+00 -1.60E-01 -1.40E+01
 
 TH 9
+        5.41E+00  3.27E+02  1.56E+04 -7.22E+01  7.56E+00 -3.59E+00  2.95E+00  1.30E+01  9.35E+06
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        6.68E+01 -5.55E+00  6.70E-03 -1.36E+01 -1.24E+01  8.51E+01 -1.20E+00  9.58E-02  4.53E+00  0.00E+00  3.08E+02
 
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
 #CPUT: Total CPU Time in Seconds,       45.654
Stop Time:
Sat Sep 18 07:09:23 CDT 2021
