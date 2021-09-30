Wed Sep 29 20:05:24 CDT 2021
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
$DATA ../../../../data/spa/D/dat48.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14075.9601763257        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3530E+02  2.7944E+02 -2.0966E+01  2.6172E+02  1.7351E+02 -1.6488E+03 -6.5809E+02 -5.6708E+01 -1.0076E+03 -4.5225E+02
            -2.7359E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -621.079466303319        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2877E+00  1.1573E+00  9.6014E-01  1.4832E+00  1.1911E+00  1.3901E+00  1.1777E+00  9.7337E-01  1.1373E+00  1.1331E+00
             1.5212E+01
 PARAMETER:  3.5285E-01  2.4611E-01  5.9321E-02  4.9421E-01  2.7490E-01  4.2937E-01  2.6358E-01  7.3014E-02  2.2864E-01  2.2498E-01
             2.8221E+00
 GRADIENT:  -5.5102E+01  2.3408E+01 -5.3516E+00  4.6634E+01 -9.8990E+00 -8.0811E-01  3.1595E+00  3.2307E+00  1.3795E+01  4.0052E+00
             1.6711E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -631.721405334074        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.2800E+00  1.0501E+00  2.7021E+00  1.6737E+00  2.9945E+00  1.4006E+00  2.1375E+00  4.4859E-01  1.1426E+00  6.3867E+00
             1.3568E+01
 PARAMETER:  3.4683E-01  1.4890E-01  1.0940E+00  6.1502E-01  1.1968E+00  4.3693E-01  8.5962E-01 -7.0164E-01  2.3334E-01  1.9542E+00
             2.7077E+00
 GRADIENT:  -1.6488E+01  3.0462E+01  4.6447E+00  6.3886E+01 -1.8298E+01 -1.1206E+01  5.9834E+00  1.6572E-02  8.7404E+00  1.5523E+01
             9.5651E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -652.928743745598        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1377E+00  1.4381E+00  2.2270E+00  8.4650E-01  7.0360E+00  1.2316E+00  7.4928E-01  5.2200E-01  6.7392E-01  1.0225E+01
             1.1790E+01
 PARAMETER:  2.2898E-01  4.6331E-01  9.0066E-01 -6.6647E-02  2.0510E+00  3.0834E-01 -1.8864E-01 -5.5009E-01 -2.9464E-01  2.4249E+00
             2.5673E+00
 GRADIENT:   1.1951E+01 -4.5216E+00  2.9120E+00 -1.2738E+01 -3.6693E+00 -2.8653E+00  1.9871E+00  4.5202E-04  1.0702E+00  5.6131E+00
             2.5309E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -656.565746851190        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.1283E+00  1.1277E+00  7.8035E-01  1.1030E+00  5.3263E+00  1.3226E+00  7.8694E-01  4.7936E-02  6.1176E-01  8.2613E+00
             1.1247E+01
 PARAMETER:  2.2070E-01  2.2015E-01 -1.4801E-01  1.9805E-01  1.7727E+00  3.7962E-01 -1.3961E-01 -2.9379E+00 -3.9142E-01  2.2116E+00
             2.5201E+00
 GRADIENT:  -9.5367E-01  2.3442E+01 -5.5770E+00  2.6567E+01 -6.2399E+00  1.1647E+00 -6.2350E-01  4.1730E-03  1.9997E-01  4.7286E+00
            -2.4376E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -668.195449903673        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0318E+00  4.5989E-01  6.7045E-01  1.2915E+00  9.8000E+00  1.1546E+00  2.9930E-01  1.0000E-02  1.5388E-01  8.3490E+00
             1.1464E+01
 PARAMETER:  1.3134E-01 -6.7676E-01 -2.9981E-01  3.5579E-01  2.3824E+00  2.4377E-01 -1.1063E+00 -5.5869E+00 -1.7716E+00  2.2221E+00
             2.5392E+00
 GRADIENT:  -3.8799E+01  1.1393E+01 -2.3726E+00  3.6205E+01  6.4190E-01 -2.7773E+01  6.8204E-02  0.0000E+00  2.0770E-01  4.5269E+00
             6.9023E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -686.649857638610        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.9196E-01  1.7026E-01  3.8751E-01  1.2844E+00  1.6429E+01  1.4175E+00  1.2620E-01  1.0000E-02  2.9820E-02  6.8080E+00
             1.0845E+01
 PARAMETER:  9.1923E-02 -1.6704E+00 -8.4801E-01  3.5032E-01  2.8991E+00  4.4889E-01 -1.9699E+00 -8.3612E+00 -3.4126E+00  2.0181E+00
             2.4837E+00
 GRADIENT:   1.9454E+01  1.9086E+00  6.4926E+00  4.3033E+01  1.0316E+01  1.6556E+01  2.0366E-03  0.0000E+00 -6.9998E-03 -2.0250E+01
            -4.4491E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -688.090913687682        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      541
 NPARAMETR:  9.8605E-01  1.7513E-01  3.7598E-01  1.2702E+00  1.6335E+01  1.4030E+00  1.2418E-01  1.0000E-02  3.0505E-02  7.1086E+00
             1.0878E+01
 PARAMETER:  8.5949E-02 -1.6422E+00 -8.7821E-01  3.3918E-01  2.8933E+00  4.3860E-01 -1.9861E+00 -8.3975E+00 -3.3899E+00  2.0613E+00
             2.4868E+00
 GRADIENT:   8.8996E+00  2.3259E+00  5.1345E+00  2.0185E+01  9.7439E+00  8.7632E+00  2.2908E-03  0.0000E+00 -3.4040E-03 -1.9894E+01
            -5.4367E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -690.085223666722        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      641
 NPARAMETR:  9.6537E-01  1.7226E-01  3.7638E-01  1.2608E+00  1.6086E+01  1.3163E+00  1.1725E-02  1.0000E-02  1.0000E-02  7.3656E+00
             1.1229E+01
 PARAMETER:  6.4758E-02 -1.6587E+00 -8.7716E-01  3.3173E-01  2.8780E+00  3.7481E-01 -4.3460E+00 -8.3975E+00 -5.6061E+00  2.0968E+00
             2.5185E+00
 GRADIENT:  -1.5421E+01  3.4557E+00  9.5861E+00  2.0017E+01  9.1469E+00 -3.8608E+00  2.4035E-05  0.0000E+00  0.0000E+00 -1.8062E+01
             6.2046E+00

0ITERATION NO.:   43    OBJECTIVE VALUE:  -690.112822860600        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:      748
 NPARAMETR:  9.6450E-01  1.6972E-01  3.7852E-01  1.2571E+00  1.6344E+01  1.3204E+00  1.0859E-02  1.0000E-02  1.0000E-02  7.2973E+00
             1.1450E+01
 PARAMETER:  6.4584E-02 -1.6590E+00 -8.7725E-01  3.3159E-01  2.8772E+00  3.7482E-01 -4.3790E+00 -8.3975E+00 -5.6622E+00  2.0980E+00
             2.5185E+00
 GRADIENT:   5.6965E+02  2.0719E+01 -5.9113E+01  9.8324E+01 -3.0189E+01 -8.6083E+01  6.3490E-06  0.0000E+00  0.0000E+00  3.6442E+01
            -2.8889E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      748
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8698E-02 -2.6715E-05  2.0697E-05 -5.2963E-04 -2.5009E-02
 SE:             2.8573E-02  1.0088E-05  1.0286E-04  3.1424E-04  8.3635E-03
 N:                     100         100         100         100         100

 P VAL.:         5.1285E-01  8.0948E-03  8.4053E-01  9.1907E-02  2.7875E-03

 ETASHRINKSD(%)  4.2763E+00  9.9966E+01  9.9655E+01  9.8947E+01  7.1981E+01
 ETASHRINKVR(%)  8.3696E+00  1.0000E+02  9.9999E+01  9.9989E+01  9.2149E+01
 EBVSHRINKSD(%)  6.2048E+00  9.9954E+01  9.9649E+01  9.8973E+01  7.7644E+01
 EBVSHRINKVR(%)  1.2025E+01  1.0000E+02  9.9999E+01  9.9989E+01  9.5002E+01
 RELATIVEINF(%)  2.6038E+01  4.7887E-06  1.1259E-04  7.4323E-04  2.9298E+00
 EPSSHRINKSD(%)  3.8275E+00
 EPSSHRINKVR(%)  7.5085E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -690.11282286059998     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       45.038003703138202     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.08
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -690.113       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.65E-01  1.72E-01  3.76E-01  1.26E+00  1.61E+01  1.32E+00  1.13E-02  1.00E-02  1.00E-02  7.37E+00  1.12E+01
 


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
+        1.59E+05
 
 TH 2
+       -1.65E+02  3.53E+04
 
 TH 3
+        1.27E+01 -1.50E+04  2.68E+04
 
 TH 4
+       -2.01E+02  1.14E+04 -2.03E+02  8.80E+03
 
 TH 5
+        8.84E-01 -2.21E+02 -5.27E-01  4.23E+00  3.50E+00
 
 TH 6
+        6.16E+01  6.93E+00  3.78E+01 -4.11E+01 -3.30E-02  6.14E+03
 
 TH 7
+       -2.46E-04  1.69E-02 -1.05E-03  1.32E-02 -1.79E-05 -2.63E-03  2.90E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -3.28E+00  3.19E+02 -2.78E+02  2.08E+02 -6.42E+00  2.60E-01 -9.35E-04  0.00E+00  0.00E+00  1.22E+01
 
 TH11
+       -1.11E+01 -1.63E+02  1.53E+02 -6.36E+00  4.40E+00  1.68E+00 -1.00E-04  0.00E+00  0.00E+00 -2.85E+00  9.13E+00
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       17.043
Stop Time:
Wed Sep 29 20:05:42 CDT 2021
