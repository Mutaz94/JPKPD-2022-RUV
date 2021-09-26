Sat Sep 25 12:57:40 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat53.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1679.73980869997        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.9389E+01 -7.4066E+01 -3.2533E+01 -4.9328E+01  2.9267E+01 -2.4806E+01  1.2315E+01  7.4077E+00  5.4236E+01  1.8313E+01
            -8.7283E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1694.86539335832        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7713E-01  1.1055E+00  1.0567E+00  9.9295E-01  1.0271E+00  1.0538E+00  9.1585E-01  9.8078E-01  7.3306E-01  9.0795E-01
             1.0222E+00
 PARAMETER:  7.6867E-02  2.0034E-01  1.5515E-01  9.2923E-02  1.2677E-01  1.5239E-01  1.2101E-02  8.0594E-02 -2.1052E-01  3.4393E-03
             1.2197E-01
 GRADIENT:   3.5431E+01  2.2138E+01  8.2496E+00  1.8551E+01 -5.5338E+00  2.6701E+00 -1.1429E+00 -2.9847E+00  5.2783E-01 -3.7801E-01
            -6.4032E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1695.18517932976        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.7428E-01  1.0644E+00  9.3697E-01  1.0116E+00  9.5449E-01  1.0631E+00  9.8178E-01  9.0512E-01  6.9702E-01  8.3216E-01
             1.0296E+00
 PARAMETER:  7.3948E-02  1.6238E-01  3.4895E-02  1.1154E-01  5.3427E-02  1.6117E-01  8.1612E-02  3.1172E-04 -2.6094E-01 -8.3725E-02
             1.2920E-01
 GRADIENT:   2.7235E+01  1.4795E+01 -1.1427E+00  2.3336E+01  1.2048E+00  6.2719E+00  2.2407E-01  3.8131E-01  6.3456E-01 -4.5684E-02
             4.6890E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1695.18783529697        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.7234E-01  1.0659E+00  9.1001E-01  1.0076E+00  9.4219E-01  1.0615E+00  9.8308E-01  8.6875E-01  6.9570E-01  8.2214E-01
             1.0275E+00
 PARAMETER:  7.1949E-02  1.6386E-01  5.7056E-03  1.0756E-01  4.0451E-02  1.5965E-01  8.2937E-02 -4.0704E-02 -2.6284E-01 -9.5849E-02
             1.2710E-01
 GRADIENT:   2.2999E+01  1.2822E+01 -1.2169E+00  2.0385E+01  1.5032E+00  5.4825E+00  1.1999E-01  3.6920E-01  6.2280E-01 -1.3361E-01
             3.9694E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1695.45037474003        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      392
 NPARAMETR:  9.8401E-01  1.1250E+00  7.2732E-01  9.5154E-01  8.7045E-01  1.0841E+00  9.6034E-01  6.2727E-01  7.0370E-01  7.6460E-01
             1.0109E+00
 PARAMETER:  8.3881E-02  2.1779E-01 -2.1839E-01  5.0330E-02 -3.8750E-02  1.8078E-01  5.9527E-02 -3.6638E-01 -2.5140E-01 -1.6841E-01
             1.1082E-01
 GRADIENT:   4.7319E+00 -4.7907E-01  1.4936E+00 -2.0598E+00 -3.9483E+00  2.4399E+00  4.9359E-01  4.0104E-01  3.4252E-01  6.3946E-01
            -2.0182E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1695.71370487969        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      568
 NPARAMETR:  9.8231E-01  1.2812E+00  5.5470E-01  8.4106E-01  8.4679E-01  1.0780E+00  8.6402E-01  3.2202E-01  7.4570E-01  7.2250E-01
             1.0148E+00
 PARAMETER:  8.2154E-02  3.4776E-01 -4.8933E-01 -7.3091E-02 -6.6303E-02  1.7511E-01 -4.6165E-02 -1.0331E+00 -1.9343E-01 -2.2503E-01
             1.1471E-01
 GRADIENT:  -8.7970E-01 -2.0769E+00 -2.0453E+00  1.8545E+00  4.1181E+00 -3.7870E-01 -4.2941E-02  1.1302E-01 -4.7144E-01 -2.5273E-01
            -6.2764E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1695.75788451515        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      743
 NPARAMETR:  9.8279E-01  1.3810E+00  5.0788E-01  7.7728E-01  8.7030E-01  1.0790E+00  8.1182E-01  2.1466E-01  7.9010E-01  7.3390E-01
             1.0165E+00
 PARAMETER:  8.2641E-02  4.2282E-01 -5.7750E-01 -1.5195E-01 -3.8916E-02  1.7606E-01 -1.0848E-01 -1.4387E+00 -1.3559E-01 -2.0938E-01
             1.1634E-01
 GRADIENT:  -7.8979E-03  1.2796E+00  4.3745E-02  6.7295E-01 -5.7030E-01 -2.9300E-02  1.1990E-02  4.7650E-02  4.2507E-02  9.3353E-02
            -2.0672E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1695.76685968037        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  9.8278E-01  1.3891E+00  5.0036E-01  7.7019E-01  8.7103E-01  1.0797E+00  8.0754E-01  1.3677E-01  7.9665E-01  7.3661E-01
             1.0163E+00
 PARAMETER:  8.2630E-02  4.2863E-01 -5.9242E-01 -1.6111E-01 -3.8073E-02  1.7673E-01 -1.1377E-01 -1.8894E+00 -1.2734E-01 -2.0569E-01
             1.1612E-01
 GRADIENT:   1.5748E-02 -2.0594E+00 -2.3285E-01 -1.5032E+00  8.6718E-02  2.2869E-01  1.7181E-01  2.0132E-02  4.3206E-01  4.7553E-01
             4.7987E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1695.78001876813        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1097
 NPARAMETR:  9.8282E-01  1.3595E+00  5.0451E-01  7.8884E-01  8.5575E-01  1.0792E+00  8.2213E-01  2.1405E-02  7.8050E-01  7.2573E-01
             1.0163E+00
 PARAMETER:  8.2670E-02  4.0708E-01 -5.8417E-01 -1.3719E-01 -5.5773E-02  1.7625E-01 -9.5855E-02 -3.7441E+00 -1.4782E-01 -2.2057E-01
             1.1613E-01
 GRADIENT:   7.9360E-03  3.5322E-01  9.0803E-02  1.7328E-01 -1.9180E-01  9.5413E-03  1.9025E-02  3.4900E-04  1.5489E-02  1.6206E-02
            -8.3493E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1695.78023352572        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  9.8281E-01  1.3615E+00  5.0294E-01  7.8733E-01  8.5591E-01  1.0792E+00  8.2098E-01  1.0000E-02  7.8116E-01  7.2514E-01
             1.0165E+00
 PARAMETER:  8.2662E-02  4.0855E-01 -5.8729E-01 -1.3911E-01 -5.5592E-02  1.7623E-01 -9.7257E-02 -4.6486E+00 -1.4698E-01 -2.2139E-01
             1.1633E-01
 GRADIENT:  -5.2158E-03  1.3653E-02  4.2980E-03  1.7002E-03 -1.3789E-02  7.9934E-04 -3.6228E-03  0.0000E+00 -2.5512E-03 -1.8397E-03
             3.4205E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1224
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.6934E-05 -1.0461E-02 -3.3165E-04  6.2355E-03 -1.7381E-02
 SE:             2.9870E-02  2.4171E-02  1.5974E-04  2.2563E-02  2.1616E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9875E-01  6.6516E-01  3.7876E-02  7.8227E-01  4.2134E-01

 ETASHRINKSD(%)  1.0000E-10  1.9025E+01  9.9465E+01  2.4412E+01  2.7584E+01
 ETASHRINKVR(%)  1.0000E-10  3.4431E+01  9.9997E+01  4.2865E+01  4.7559E+01
 EBVSHRINKSD(%)  3.5975E-01  1.9080E+01  9.9510E+01  2.5282E+01  2.7151E+01
 EBVSHRINKVR(%)  7.1821E-01  3.4519E+01  9.9998E+01  4.4172E+01  4.6931E+01
 RELATIVEINF(%)  9.9255E+01  2.6242E+00  1.2297E-04  1.9722E+00  3.8749E+00
 EPSSHRINKSD(%)  4.3300E+01
 EPSSHRINKVR(%)  6.7852E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1695.7802335257168     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -960.62940696197859     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.84
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1695.780       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.36E+00  5.03E-01  7.87E-01  8.56E-01  1.08E+00  8.21E-01  1.00E-02  7.81E-01  7.25E-01  1.02E+00
 


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
+        9.80E+02
 
 TH 2
+       -7.29E+00  5.94E+02
 
 TH 3
+        7.80E+00  3.01E+02  1.12E+03
 
 TH 4
+       -1.99E+01  4.88E+02 -7.89E+02  1.65E+03
 
 TH 5
+       -5.68E+00 -4.65E+02 -1.17E+03  7.11E+02  1.53E+03
 
 TH 6
+       -1.38E-01 -1.34E+00  2.37E+00 -4.84E+00 -1.71E+00  1.69E+02
 
 TH 7
+       -7.32E-02  2.44E+01 -3.97E+01 -1.39E+01  6.43E+00 -3.68E-01  1.29E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.06E+00 -2.13E+01 -5.41E+01  5.85E+01  1.89E-01 -8.75E-01  2.76E+01  0.00E+00  1.18E+02
 
 TH10
+       -1.04E+00 -1.75E+01 -6.27E+01 -1.74E+01 -8.06E+01 -2.53E-01  2.60E+01  0.00E+00  1.69E+01  1.19E+02
 
 TH11
+       -5.34E+00 -1.65E+01 -3.04E+01 -1.71E+00 -3.17E+00  2.74E+00  1.07E+01  0.00E+00  1.69E+01  2.36E+01  2.07E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.479
Stop Time:
Sat Sep 25 12:58:01 CDT 2021
