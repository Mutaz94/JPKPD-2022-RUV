Fri Sep 24 22:17:28 CDT 2021
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
$DATA ../../../../data/int/A3/dat37.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1990.52728471018        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.1050E+01  8.1864E+01  2.1570E+02 -7.6354E+01  7.3081E+01  4.4578E+00 -7.0994E+01 -1.3041E+02  6.8937E+00 -6.6455E+01
            -3.4334E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3006.24867643096        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9459E-01  9.6754E-01  8.8767E-01  1.0586E+00  9.3857E-01  9.8016E-01  1.2733E+00  8.4977E-01  8.0050E-01  1.0564E+00
             1.9613E+00
 PARAMETER:  9.4573E-02  6.7003E-02 -1.9157E-02  1.5697E-01  3.6607E-02  7.9958E-02  3.4161E-01 -6.2785E-02 -1.2252E-01  1.5486E-01
             7.7359E-01
 GRADIENT:   1.0723E+01  1.0350E+01  2.0179E+01 -4.5969E+00 -8.8242E+00  1.1476E+00 -2.9720E+00  8.1683E+00 -6.6040E+00  5.6011E+00
            -1.5350E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3011.04781597299        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0012E+00  8.1042E-01  6.8066E-01  1.1555E+00  7.4749E-01  1.0250E+00  1.3859E+00  4.4567E-01  8.7077E-01  9.2538E-01
             1.9463E+00
 PARAMETER:  1.0123E-01 -1.1020E-01 -2.8469E-01  2.4452E-01 -1.9104E-01  1.2468E-01  4.2638E-01 -7.0817E-01 -3.8382E-02  2.2448E-02
             7.6592E-01
 GRADIENT:   2.5403E+01  3.7154E+01 -2.6363E+01  5.0826E+01  1.9732E+01  1.6527E+01  5.4848E+00  2.2066E+00  5.1976E+00 -1.2198E+00
            -1.4824E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3020.59454380265        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.9013E-01  5.4921E-01  5.0210E-01  1.2206E+00  5.2424E-01  9.6150E-01  1.2941E+00  1.0207E-01  8.7435E-01  8.7895E-01
             2.1210E+00
 PARAMETER:  9.0077E-02 -4.9928E-01 -5.8896E-01  2.9932E-01 -5.4581E-01  6.0740E-02  3.5782E-01 -2.1821E+00 -3.4278E-02 -2.9024E-02
             8.5191E-01
 GRADIENT:  -6.9016E+00 -1.5696E+00  2.1507E+00  1.2214E+01  3.1257E+01 -6.0963E+00  2.6922E+00  1.5049E-01  5.9496E-01 -1.1061E+01
             9.1296E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3020.96844431883        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  9.8800E-01  3.7751E-01  3.3830E-01  1.2505E+00  3.5996E-01  9.4953E-01  1.1961E+00  2.4519E-02  9.1553E-01  8.6155E-01
             2.1697E+00
 PARAMETER:  8.7929E-02 -8.7417E-01 -9.8383E-01  3.2358E-01 -9.2176E-01  4.8215E-02  2.7908E-01 -3.6083E+00  1.1750E-02 -4.9018E-02
             8.7461E-01
 GRADIENT:  -1.7897E+01 -2.6986E+01  4.4957E+01  8.7180E+01  1.4029E+01 -1.1922E+01 -1.5749E+01  1.4100E-02  6.4405E+00 -8.4350E+00
             1.8001E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3021.02339690491        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      405
 NPARAMETR:  9.8765E-01  3.3924E-01  2.9754E-01  1.2473E+00  3.2187E-01  9.4820E-01  1.1683E+00  1.4960E-02  9.3110E-01  8.7167E-01
             2.1702E+00
 PARAMETER:  8.7568E-02 -9.8106E-01 -1.1122E+00  3.2102E-01 -1.0336E+00  4.6806E-02  2.5556E-01 -4.1024E+00  2.8615E-02 -3.7341E-02
             8.7483E-01
 GRADIENT:  -2.7957E+01 -4.0850E+01  5.0369E+01  1.0783E+02 -2.2176E+01 -1.3677E+01 -2.2714E+01  4.8992E-03  4.1542E+00 -3.2544E+00
             1.9217E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3029.03992957664        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      581
 NPARAMETR:  9.9363E-01  3.1024E-01  2.6620E-01  1.2572E+00  2.8916E-01  9.7151E-01  1.2307E+00  1.1650E-02  9.4061E-01  8.8277E-01
             2.0403E+00
 PARAMETER:  9.3605E-02 -1.0704E+00 -1.2235E+00  3.2890E-01 -1.1408E+00  7.1096E-02  3.0761E-01 -4.3524E+00  3.8768E-02 -2.4686E-02
             8.1309E-01
 GRADIENT:  -1.0156E+01 -3.8696E+01  8.6517E+01  1.6931E+02 -8.0062E+01 -5.2850E+00 -2.2829E+01  1.4180E-03 -7.8443E+00  2.0960E+00
             9.4638E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3073.39096194699        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  9.9803E-01  1.7778E-01  9.9083E-02  8.4699E-01  1.6084E-01  9.8896E-01  1.6289E+00  1.0000E-02  1.2776E+00  8.8591E-01
             1.7880E+00
 PARAMETER:  9.8030E-02 -1.6272E+00 -2.2118E+00 -6.6070E-02 -1.7273E+00  8.8902E-02  5.8790E-01 -9.9235E+00  3.4495E-01 -2.1136E-02
             6.8112E-01
 GRADIENT:   3.4611E+00  1.3035E+01 -8.4002E+00 -1.7379E+01  1.5795E+01 -1.5455E+00  5.8637E+00  0.0000E+00 -1.7675E+01 -2.9813E+00
             1.1464E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3075.33452695462        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      933
 NPARAMETR:  9.9407E-01  1.5852E-01  8.2429E-02  7.8895E-01  1.4544E-01  9.9723E-01  1.6229E+00  1.0000E-02  1.5564E+00  9.0617E-01
             1.7586E+00
 PARAMETER:  9.4056E-02 -1.7419E+00 -2.3958E+00 -1.3706E-01 -1.8280E+00  9.7225E-02  5.8420E-01 -1.1074E+01  5.4236E-01  1.4699E-03
             6.6451E-01
 GRADIENT:  -1.5572E-01 -1.2080E+00  5.1135E-01  9.3015E-01  4.2219E-01  2.0086E-01  4.8897E-02  0.0000E+00  2.5153E-01 -1.1385E-01
            -6.6292E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -3075.33799265804        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1025
 NPARAMETR:  9.9410E-01  1.5811E-01  8.1726E-02  7.8474E-01  1.4482E-01  9.9664E-01  1.6234E+00  1.0000E-02  1.5648E+00  9.0717E-01
             1.7586E+00
 PARAMETER:  9.4083E-02 -1.7445E+00 -2.4044E+00 -1.4240E-01 -1.8323E+00  9.6632E-02  5.8451E-01 -1.1137E+01  5.4777E-01  2.5781E-03
             6.6453E-01
 GRADIENT:  -5.4342E-04  1.9197E-02 -6.4939E-02  4.7909E-02  7.1056E-02 -1.3655E-03 -4.0468E-03  0.0000E+00  1.2500E-02 -1.5762E-02
             1.8652E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1025
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8861E-03  1.0847E-02  3.8331E-04 -2.8912E-03  7.4548E-03
 SE:             2.9533E-02  2.7316E-02  2.4458E-04  2.5880E-02  2.8899E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4908E-01  6.9129E-01  1.1706E-01  9.1105E-01  7.9643E-01

 ETASHRINKSD(%)  1.0598E+00  8.4891E+00  9.9181E+01  1.3299E+01  3.1860E+00
 ETASHRINKVR(%)  2.1084E+00  1.6258E+01  9.9993E+01  2.4829E+01  6.2705E+00
 EBVSHRINKSD(%)  1.0192E+00  7.5923E+00  9.9317E+01  8.0927E+00  3.9859E+00
 EBVSHRINKVR(%)  2.0280E+00  1.4608E+01  9.9995E+01  1.5531E+01  7.8129E+00
 RELATIVEINF(%)  9.7884E+01  3.9336E+01  7.1366E-04  3.6138E+01  1.8944E+01
 EPSSHRINKSD(%)  2.1317E+01
 EPSSHRINKVR(%)  3.8090E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3075.3379926580355     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1421.2486328896248     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3075.338       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  1.58E-01  8.17E-02  7.85E-01  1.45E-01  9.97E-01  1.62E+00  1.00E-02  1.56E+00  9.07E-01  1.76E+00
 


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
+        1.10E+03
 
 TH 2
+        5.11E+00  8.43E+07
 
 TH 3
+       -7.62E+00 -1.18E+08  1.65E+08
 
 TH 4
+        8.18E+00 -2.06E+02 -2.48E+03  6.15E+02
 
 TH 5
+       -1.20E+01 -2.63E+03 -4.48E+04 -4.29E+02  5.91E+04
 
 TH 6
+        1.14E+01 -2.07E+00  6.00E+00 -6.78E+00 -4.40E-01  1.77E+02
 
 TH 7
+       -1.04E+00  5.36E+01  1.25E+02 -4.76E+00 -3.74E+01 -5.65E-01  5.12E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.24E+00  4.09E+01  5.99E+02 -3.21E+00  1.13E+02 -2.89E+00  2.77E+00  0.00E+00  4.79E+01
 
 TH10
+        2.40E+00 -1.46E+01  1.98E+02  5.68E+00 -1.78E+01  3.69E+00 -2.71E-01  0.00E+00 -7.72E-01  1.98E+02
 
 TH11
+       -1.66E+01 -3.37E+01 -1.77E+02  3.08E+00  5.23E+01  4.87E+00  8.25E+00  0.00E+00  9.69E+00  9.33E+00  3.32E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.185
Stop Time:
Fri Sep 24 22:18:06 CDT 2021
