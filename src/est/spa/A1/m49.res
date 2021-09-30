Wed Sep 29 12:10:17 CDT 2021
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
$DATA ../../../../data/spa/A1/dat49.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -888.960682679994        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0004E+02 -1.8146E+01  5.7038E+01 -1.8451E+01  7.2891E+01  2.2162E+01 -4.4495E+01 -2.2218E+02 -8.3575E+01 -1.3230E+01
            -1.0874E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1480.32468156211        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      129
 NPARAMETR:  9.9077E-01  1.0996E+00  1.0822E+00  9.6890E-01  1.1159E+00  8.7633E-01  1.5966E+00  8.7550E-01  1.2973E+00  9.4579E-01
             1.4946E+00
 PARAMETER:  9.0727E-02  1.9491E-01  1.7897E-01  6.8408E-02  2.0965E-01 -3.2017E-02  5.6785E-01 -3.2965E-02  3.6031E-01  4.4270E-02
             5.0184E-01
 GRADIENT:   1.1576E+02 -1.3669E+01 -1.4738E+01  1.4904E+01  6.5005E+01 -6.6669E+01  1.7557E+01 -2.3573E+00  2.8224E+01 -8.5853E+00
            -2.5792E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1485.77991152785        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.6478E-01  1.0377E+00  1.0291E+00  1.0443E+00  1.0202E+00  9.8442E-01  1.8336E+00  4.4053E-01  1.1366E+00  7.6188E-01
             1.7192E+00
 PARAMETER:  6.4143E-02  1.3702E-01  1.2869E-01  1.4331E-01  1.1999E-01  8.4296E-02  7.0626E-01 -7.1978E-01  2.2805E-01 -1.7197E-01
             6.4186E-01
 GRADIENT:   2.7892E+01  2.0117E+01 -6.6368E+00  3.7032E+01  4.3820E+01 -9.4910E+00  2.5969E+01 -1.0657E+00  1.6967E+01 -8.9527E+00
             2.3933E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1494.46368449048        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      482
 NPARAMETR:  9.5257E-01  1.1655E+00  6.9314E-01  9.0004E-01  8.8516E-01  1.0085E+00  1.4446E+00  5.4115E-01  1.0471E+00  6.6975E-01
             1.5855E+00
 PARAMETER:  5.1404E-02  2.5314E-01 -2.6653E-01 -5.3198E-03 -2.1983E-02  1.0850E-01  4.6782E-01 -5.1406E-01  1.4598E-01 -3.0086E-01
             5.6090E-01
 GRADIENT:  -4.0661E+00  2.2000E+00 -2.2134E+00  1.6816E+00  3.7042E+00 -9.4206E-01  2.6353E+00  8.8354E-01  5.2755E-01 -9.7178E-01
             6.8191E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1494.66024071830        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  9.5322E-01  1.0429E+00  7.4909E-01  9.7534E-01  8.5966E-01  1.0100E+00  1.5506E+00  5.3609E-01  1.0013E+00  7.0061E-01
             1.5660E+00
 PARAMETER:  5.2095E-02  1.4205E-01 -1.8890E-01  7.5033E-02 -5.1224E-02  1.0992E-01  5.3867E-01 -5.2345E-01  1.0127E-01 -2.5581E-01
             5.4853E-01
 GRADIENT:  -8.2034E-01 -3.5269E-01 -4.6240E-01  5.8791E-01  3.2680E-01 -1.1065E-01  1.1017E-01  6.0539E-01  1.8242E-01  1.3801E-02
             3.2214E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1494.84911986117        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  9.5429E-01  1.0937E+00  6.5692E-01  9.3428E-01  8.3053E-01  1.0110E+00  1.4973E+00  2.7801E-01  1.0088E+00  6.7056E-01
             1.5533E+00
 PARAMETER:  5.3217E-02  1.8961E-01 -3.2019E-01  3.2025E-02 -8.5687E-02  1.1093E-01  5.0366E-01 -1.1801E+00  1.0875E-01 -2.9965E-01
             5.4036E-01
 GRADIENT:  -5.7275E-02 -3.9868E-01 -9.6033E-01 -2.4095E-01 -6.3519E-02 -1.1612E-01  3.1167E-01  1.5761E-01  5.1252E-02  4.9038E-01
             3.0591E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1494.89212068375        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1010
 NPARAMETR:  9.5451E-01  1.1045E+00  6.3914E-01  9.2579E-01  8.2532E-01  1.0105E+00  1.4849E+00  7.6430E-02  1.0115E+00  6.6393E-01
             1.5545E+00
 PARAMETER:  5.3439E-02  1.9940E-01 -3.4763E-01  2.2891E-02 -9.1981E-02  1.1047E-01  4.9532E-01 -2.4714E+00  1.1145E-01 -3.0957E-01
             5.4115E-01
 GRADIENT:   6.4705E-02 -1.8854E-01 -1.4573E-01 -1.1508E-01 -7.9522E-02 -3.4610E-01  3.8994E-02  5.8450E-03  9.6634E-03  1.2360E-01
             1.6420E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1494.89480928705        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  9.5446E-01  1.1108E+00  6.3693E-01  9.2205E-01  8.2698E-01  1.0114E+00  1.4784E+00  1.2790E-02  1.0143E+00  6.6326E-01
             1.5545E+00
 PARAMETER:  5.3395E-02  2.0506E-01 -3.5109E-01  1.8841E-02 -8.9979E-02  1.1134E-01  4.9098E-01 -4.2591E+00  1.1418E-01 -3.1059E-01
             5.4114E-01
 GRADIENT:  -7.6519E-02 -2.6157E-02 -2.4199E-02 -9.4935E-03  3.4379E-02 -2.4356E-02  3.1984E-03  1.6207E-04  4.6105E-03  2.6990E-03
             2.0855E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1494.89491400492        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  9.5451E-01  1.1101E+00  6.3712E-01  9.2251E-01  8.2671E-01  1.0115E+00  1.4792E+00  1.0000E-02  1.0139E+00  6.6325E-01
             1.5545E+00
 PARAMETER:  5.3447E-02  2.0445E-01 -3.5081E-01  1.9342E-02 -9.0302E-02  1.1143E-01  4.9148E-01 -4.7976E+00  1.1383E-01 -3.1060E-01
             5.4116E-01
 GRADIENT:   3.7149E-02  3.1322E-02  1.4280E-02  1.8413E-02 -4.8079E-02  1.2919E-02  6.3995E-04  0.0000E+00 -1.4488E-03  5.4988E-03
             1.1534E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1242
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2443E-04  3.5245E-03 -3.6134E-04 -8.0752E-03 -1.0200E-02
 SE:             2.9652E-02  2.4569E-02  1.5181E-04  2.2800E-02  1.7375E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8320E-01  8.8593E-01  1.7300E-02  7.2321E-01  5.5717E-01

 ETASHRINKSD(%)  6.6197E-01  1.7690E+01  9.9491E+01  2.3617E+01  4.1790E+01
 ETASHRINKVR(%)  1.3196E+00  3.2250E+01  9.9997E+01  4.1656E+01  6.6116E+01
 EBVSHRINKSD(%)  9.8523E-01  1.6694E+01  9.9528E+01  2.3932E+01  4.2364E+01
 EBVSHRINKVR(%)  1.9608E+00  3.0601E+01  9.9998E+01  4.2137E+01  6.6781E+01
 RELATIVEINF(%)  9.7849E+01  7.8864E+00  2.2327E-04  6.6442E+00  3.1090E+00
 EPSSHRINKSD(%)  3.9292E+01
 EPSSHRINKVR(%)  6.3145E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1494.8949140049249     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -759.74408744118671     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1494.895       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.55E-01  1.11E+00  6.37E-01  9.23E-01  8.27E-01  1.01E+00  1.48E+00  1.00E-02  1.01E+00  6.63E-01  1.55E+00
 


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
+        1.17E+03
 
 TH 2
+       -8.18E+00  2.93E+02
 
 TH 3
+        1.79E+01  1.91E+02  6.66E+02
 
 TH 4
+       -2.16E+01  1.85E+02 -3.79E+02  7.57E+02
 
 TH 5
+       -8.43E+00 -3.11E+02 -8.04E+02  4.39E+02  1.22E+03
 
 TH 6
+        1.01E+00 -1.68E+00  4.76E+00 -5.87E+00 -4.19E+00  1.87E+02
 
 TH 7
+        1.41E+00  2.30E+01 -2.14E+01 -1.58E+01  1.18E+01 -1.27E-01  4.50E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.64E+00 -1.28E+01 -3.22E+01  2.75E+01 -1.76E+00 -6.17E-01  1.05E+01  0.00E+00  7.45E+01
 
 TH10
+       -1.64E+00 -8.22E+00 -3.92E+01 -2.22E+01 -4.43E+01  9.15E-01  8.99E+00  0.00E+00  1.11E+01  5.96E+01
 
 TH11
+       -9.83E+00 -1.04E+01 -2.78E+01 -6.90E+00 -1.12E+01  2.58E+00  4.87E+00  0.00E+00  1.03E+01  2.34E+01  9.84E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.343
Stop Time:
Wed Sep 29 12:10:40 CDT 2021
