Sat Sep 25 09:41:17 CDT 2021
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
$DATA ../../../../data/spa/S1/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1159.66392529750        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.6786E+01 -1.2541E+01  4.3492E+01  5.2748E+00  7.7867E+01 -3.0088E+01  1.5397E+01 -1.9846E+02 -4.1927E+00 -1.2055E+01
            -7.6944E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1653.45461036225        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:       97
 NPARAMETR:  9.9295E-01  9.9868E-01  9.9562E-01  9.7114E-01  9.6916E-01  1.0474E+00  9.5700E-01  1.0677E+00  9.6028E-01  9.8783E-01
             1.1548E+00
 PARAMETER:  9.2924E-02  9.8679E-02  9.5611E-02  7.0713E-02  6.8673E-02  1.4627E-01  5.6044E-02  1.6554E-01  5.9469E-02  8.7751E-02
             2.4393E-01
 GRADIENT:   3.5567E+01 -2.4345E+01 -9.6272E+00 -2.0349E+01  1.7879E+01 -1.1892E+01  1.0092E+01  3.9782E+00  9.8559E+00  3.3075E+00
             3.8365E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1655.17116277911        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:      296             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9295E-01  9.9864E-01  9.9566E-01  9.7133E-01  9.6900E-01  1.1284E+00  6.9338E-01  1.0677E+00  9.6032E-01  9.8783E-01
             1.1548E+00
 PARAMETER:  9.2923E-02  9.8638E-02  9.5650E-02  7.0913E-02  6.8506E-02  2.2083E-01 -2.6618E-01  1.6554E-01  5.9509E-02  8.7751E-02
             2.4393E-01
 GRADIENT:   3.5173E+01 -2.8459E+01 -1.3075E+01 -2.2600E+01  2.3175E+01  2.2017E+01  8.7654E-01  2.5441E+00 -4.5891E+00 -1.4408E+00
             3.7126E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1655.30379465016        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      453
 NPARAMETR:  9.9295E-01  9.9864E-01  9.9566E-01  9.7133E-01  9.6900E-01  1.1006E+00  6.6732E-01  1.0677E+00  9.6032E-01  9.8783E-01
             1.1548E+00
 PARAMETER:  9.2923E-02  9.8638E-02  9.5651E-02  7.0913E-02  6.8506E-02  1.9581E-01 -3.0449E-01  1.6554E-01  5.9509E-02  8.7751E-02
             2.4393E-01
 GRADIENT:   9.7631E-02 -3.1178E+01 -1.3673E+01 -2.8347E+01  2.3041E+01  8.0718E-02  1.7279E-02  2.3570E+00 -6.8011E+00 -2.1126E+00
             3.6665E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1656.02389996455        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      631             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9292E-01  1.0068E+00  1.0044E+00  9.7524E-01  9.6475E-01  1.1003E+00  6.6594E-01  1.0511E+00  9.6752E-01  9.9139E-01
             1.1472E+00
 PARAMETER:  9.2894E-02  1.0674E-01  1.0439E-01  7.4926E-02  6.4110E-02  1.9563E-01 -3.0655E-01  1.4981E-01  6.6983E-02  9.1355E-02
             2.3731E-01
 GRADIENT:   3.5498E+01 -9.4203E+00 -4.7531E+00 -6.4624E+00  7.5442E+00  1.1057E+01  6.1703E-01  1.0646E+00 -4.7652E+00 -1.6169E+00
             3.3948E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1656.13632836827        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.8373E-01  1.0113E+00  1.0130E+00  9.7554E-01  9.6483E-01  1.0866E+00  6.0179E-01  1.0409E+00  9.9207E-01  1.0066E+00
             1.1472E+00
 PARAMETER:  8.3596E-02  1.1120E-01  1.1296E-01  7.5234E-02  6.4200E-02  1.8306E-01 -4.0784E-01  1.4008E-01  9.2036E-02  1.0653E-01
             2.3731E-01
 GRADIENT:   1.8183E+01 -5.5565E-01 -9.0211E-01  1.1303E+00 -1.6785E+00  5.4705E+00  6.7534E-01  1.0874E-01 -2.7820E+00 -1.1464E+00
             3.3885E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1656.15328486837        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      836             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8247E-01  1.0115E+00  1.0217E+00  9.7500E-01  9.6999E-01  1.1174E+00  5.4723E-01  1.0450E+00  9.9810E-01  1.0108E+00
             1.1473E+00
 PARAMETER:  8.2319E-02  1.1148E-01  1.2147E-01  7.4684E-02  6.9527E-02  2.1102E-01 -5.0289E-01  1.4397E-01  9.8098E-02  1.1076E-01
             2.3743E-01
 GRADIENT:   1.6928E+01  2.5573E-01 -1.2097E+00  1.1977E+00  7.1171E-02  1.8072E+01  2.8751E-01 -5.1228E-01 -4.6772E+00 -2.3683E+00
             3.3535E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1656.19394821385        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  9.8247E-01  1.0115E+00  1.0237E+00  9.7500E-01  9.7029E-01  1.0978E+00  5.4305E-01  1.0450E+00  9.9810E-01  1.0108E+00
             1.1473E+00
 PARAMETER:  8.2319E-02  1.1148E-01  1.2338E-01  7.4684E-02  6.9836E-02  1.9331E-01 -5.1055E-01  1.4397E-01  9.8098E-02  1.1076E-01
             2.3743E-01
 GRADIENT:  -1.9336E+01 -2.9922E+00 -9.0523E-01 -4.4650E+00 -1.1826E+00 -1.1904E+00 -2.3090E-01 -6.5399E-01 -5.5066E+00 -2.6348E+00
             3.3191E+01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1656.19689584020        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1173
 NPARAMETR:  9.8251E-01  1.0115E+00  1.0236E+00  9.7499E-01  9.7032E-01  1.1011E+00  5.6172E-01  1.0450E+00  9.9817E-01  1.0109E+00
             1.1474E+00
 PARAMETER:  8.2294E-02  1.1145E-01  1.2335E-01  7.4659E-02  6.9861E-02  1.9629E-01 -4.7676E-01  1.4393E-01  9.8073E-02  1.1073E-01
             2.3749E-01
 GRADIENT:  -1.9205E+01 -4.0408E+00 -7.7711E-01 -5.0229E+00 -1.1346E+00 -1.0289E-02 -1.5309E-03 -5.1283E-01 -4.3287E+00 -2.2279E+00
            -8.2995E+05
 NUMSIGDIG:         0.7         1.2         1.5         1.3         1.8         3.9         3.0         0.5         0.5         0.6
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1173
 NO. OF SIG. DIGITS IN FINAL EST.:  0.5
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0022E-02 -2.4732E-02 -2.5733E-02  7.2701E-03 -3.1288E-02
 SE:             2.9842E-02  1.1754E-02  1.3672E-02  2.7306E-02  2.2948E-02
 N:                     100         100         100         100         100

 P VAL.:         7.3700E-01  3.5365E-02  5.9804E-02  7.9005E-01  1.7275E-01

 ETASHRINKSD(%)  2.5033E-02  6.0623E+01  5.4198E+01  8.5201E+00  2.3122E+01
 ETASHRINKVR(%)  5.0059E-02  8.4494E+01  7.9022E+01  1.6314E+01  4.0897E+01
 EBVSHRINKSD(%)  4.5227E-01  6.1647E+01  5.7763E+01  9.7774E+00  2.1786E+01
 EBVSHRINKVR(%)  9.0250E-01  8.5291E+01  8.2160E+01  1.8599E+01  3.8825E+01
 RELATIVEINF(%)  9.8690E+01  5.8979E-01  2.9327E+00  4.1406E+00  1.0613E+01
 EPSSHRINKSD(%)  4.7553E+01
 EPSSHRINKVR(%)  7.2493E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1656.1968958402006     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -921.04606927646239     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1656.197       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.01E+00  1.02E+00  9.75E-01  9.70E-01  1.10E+00  5.62E-01  1.04E+00  9.98E-01  1.01E+00  1.15E+00
 


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
+        9.40E+02
 
 TH 2
+       -1.44E+01  5.27E+02
 
 TH 3
+        7.36E+00  1.31E+02  1.84E+02
 
 TH 4
+       -1.48E+01  5.80E+02 -4.48E+01  9.51E+02
 
 TH 5
+        2.13E+00 -2.96E+02 -2.90E+02  1.40E+00  7.12E+02
 
 TH 6
+        1.68E+01 -4.04E+00  1.43E+00 -1.97E+00 -4.44E-01  1.61E+02
 
 TH 7
+        1.24E+00 -2.40E+01  3.64E+00 -1.41E+01 -2.43E+00  1.68E-01  1.10E+01
 
 TH 8
+        3.33E+09  2.91E+09 -2.72E+01 -7.07E-01  6.31E+00 -6.68E-01  3.81E+00  2.18E+09
 
 TH 9
+        3.56E+00 -4.88E+01  1.34E+00  2.97E+01 -1.81E+00  8.31E-01  3.04E+01  7.52E-01  1.36E+02
 
 TH10
+       -1.56E+00  2.48E+00 -1.40E+01  2.67E+00 -6.29E+01  5.96E-01  1.15E+01  2.93E+09 -2.54E+00  7.74E+01
 
 TH11
+       -3.43E+03  1.21E+04 -5.59E+04  3.01E+04 -1.61E+04 -4.00E+03 -1.98E+03  1.64E+05  1.54E+04  5.00E+03  6.64E+08
 
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
 #CPUT: Total CPU Time in Seconds,       21.303
Stop Time:
Sat Sep 25 09:41:40 CDT 2021
