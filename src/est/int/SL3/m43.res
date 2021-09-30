Wed Sep 29 04:18:26 CDT 2021
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
$DATA ../../../../data/int/SL3/dat43.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      986
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

 TOT. NO. OF OBS RECS:      886
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -922.006607015781        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4369E+02  5.1105E+01  8.9480E+01  2.0058E+02  1.7735E+02  4.6390E+01 -1.1669E+02 -1.8379E+02 -8.2876E+01 -4.2809E+01
            -5.3134E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2765.79481293460        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9813E-01  1.2424E+00  9.6986E-01  8.8146E-01  1.0709E+00  9.3828E-01  1.0067E+00  9.1712E-01  7.5690E-01  1.0973E+00
             2.3948E+00
 PARAMETER:  9.8125E-02  3.1706E-01  6.9397E-02 -2.6180E-02  1.6852E-01  3.6291E-02  1.0670E-01  1.3484E-02 -1.7852E-01  1.9281E-01
             9.7332E-01
 GRADIENT:   6.6890E+01  5.3048E+01 -8.8245E+00  1.1039E+01 -3.8521E+00 -8.6906E+00  2.2379E+01  5.5507E+00 -1.4205E+01 -2.3952E+01
            -1.3916E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2781.35793399326        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      183
 NPARAMETR:  9.9837E-01  1.4018E+00  1.1759E+00  8.1184E-01  1.2541E+00  9.8463E-01  8.6793E-01  4.2591E-01  8.3878E-01  1.3519E+00
             2.4433E+00
 PARAMETER:  9.8367E-02  4.3774E-01  2.6204E-01 -1.0846E-01  3.2642E-01  8.4506E-02 -4.1644E-02 -7.5352E-01 -7.5808E-02  4.0148E-01
             9.9335E-01
 GRADIENT:  -7.5482E+00  1.9017E+01 -2.1976E+00  4.0173E+01 -7.7857E+00  2.4601E+00  1.7775E+01  5.4339E-01 -5.7523E+00 -6.4005E+00
            -9.5372E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2794.55985725900        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      362
 NPARAMETR:  1.0028E+00  1.7854E+00  1.3176E+00  5.8534E-01  1.6311E+00  9.7546E-01  4.9388E-01  3.0320E-01  1.3414E+00  1.6040E+00
             2.5573E+00
 PARAMETER:  1.0283E-01  6.7962E-01  3.7580E-01 -4.3557E-01  5.8927E-01  7.5155E-02 -6.0547E-01 -1.0933E+00  3.9373E-01  5.7250E-01
             1.0390E+00
 GRADIENT:  -1.8128E+00  5.8707E+01 -4.4031E+00  5.6116E+01  1.9269E+01 -3.5774E-01 -2.6846E+00 -1.6178E-02  4.6258E+00 -8.5685E+00
             3.0280E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2800.42215218004        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      537
 NPARAMETR:  1.0033E+00  2.0692E+00  1.4690E+00  3.6177E-01  1.8127E+00  9.7645E-01  4.7582E-01  3.2115E-01  1.7794E+00  1.8061E+00
             2.5084E+00
 PARAMETER:  1.0333E-01  8.2716E-01  4.8456E-01 -9.1674E-01  6.9482E-01  7.6169E-02 -6.4272E-01 -1.0359E+00  6.7630E-01  6.9118E-01
             1.0196E+00
 GRADIENT:   2.7092E+00  1.7135E+01  3.5464E+00  5.4027E+00 -5.8452E+00  2.9545E-01 -1.2074E+00 -1.7763E-01 -1.1651E+00 -6.4315E-02
             6.2649E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2801.63175782920        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  1.0027E+00  2.2675E+00  9.1312E-01  2.2203E-01  1.8991E+00  9.7645E-01  4.8126E-01  3.2766E-01  2.3930E+00  1.8374E+00
             2.5004E+00
 PARAMETER:  1.0265E-01  9.1869E-01  9.1113E-03 -1.4049E+00  7.4135E-01  7.6165E-02 -6.3134E-01 -1.0158E+00  9.7256E-01  7.0832E-01
             1.0165E+00
 GRADIENT:   7.7188E-01  3.7915E+00 -3.3009E-01  1.3529E+00  8.6157E-01 -1.8733E-01 -2.8242E-01  1.4388E-02  6.7121E-02 -1.2225E-01
            -1.6138E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2801.64493612722        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  1.0022E+00  2.2679E+00  9.1956E-01  2.1978E-01  1.8978E+00  9.7687E-01  4.8175E-01  1.5420E-01  2.4020E+00  1.8393E+00
             2.5018E+00
 PARAMETER:  1.0219E-01  9.1883E-01  1.6142E-02 -1.4151E+00  7.4068E-01  7.6600E-02 -6.3034E-01 -1.7695E+00  9.7632E-01  7.0939E-01
             1.0170E+00
 GRADIENT:  -1.8634E-01 -2.1715E+00  4.7477E-03  1.0471E-01 -7.8655E-02  1.6609E-02 -4.8204E-02  1.7398E-03 -2.3805E-01 -2.0717E-02
             7.9429E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2801.64619197187        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  1.0024E+00  2.2711E+00  9.1221E-01  2.1735E-01  1.8988E+00  9.7689E-01  4.8221E-01  5.9590E-02  2.4250E+00  1.8413E+00
             2.5016E+00
 PARAMETER:  1.0239E-01  9.2028E-01  8.1132E-03 -1.4262E+00  7.4121E-01  7.6620E-02 -6.2938E-01 -2.7203E+00  9.8581E-01  7.1050E-01
             1.0169E+00
 GRADIENT:   2.8008E-01 -3.0985E+00 -4.0398E-02  2.1095E-01 -2.2447E-01  1.0119E-02  1.3411E-01  3.3544E-04  2.0092E-01  1.9382E-01
             6.6413E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2801.64646672884        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1236
 NPARAMETR:  1.0024E+00  2.2710E+00  9.1368E-01  2.1727E-01  1.8987E+00  9.7690E-01  4.8186E-01  1.0000E-02  2.4243E+00  1.8412E+00
             2.5014E+00
 PARAMETER:  1.0237E-01  9.2024E-01  9.7281E-03 -1.4266E+00  7.4117E-01  7.6630E-02 -6.3010E-01 -6.6488E+00  9.8555E-01  7.1040E-01
             1.0169E+00
 GRADIENT:   7.1344E+01  3.7853E+02  3.3634E-02  1.4734E+01  3.3913E+01  7.6049E+00  5.5807E+00  0.0000E+00  7.8603E+00  1.1110E+01
             1.7954E+01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -2801.64646906070        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:     1356
 NPARAMETR:  1.0025E+00  2.2707E+00  9.1459E-01  2.1744E-01  1.8986E+00  9.7701E-01  4.8196E-01  1.0000E-02  2.4220E+00  1.8412E+00
             2.5016E+00
 PARAMETER:  1.0235E-01  9.2020E-01  9.7258E-03 -1.4265E+00  7.4112E-01  7.6613E-02 -6.3016E-01 -7.0760E+00  9.8550E-01  7.1036E-01
             1.0169E+00
 GRADIENT:  -5.5288E-02  8.7241E-02 -2.9278E-03 -1.6673E-02  5.4963E-04 -8.5496E-03 -6.0715E-03  0.0000E+00  2.0274E-02 -1.5401E-03
            -3.4302E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1356
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6499E-03 -4.1206E-02 -3.8362E-05  3.7328E-02 -1.9085E-02
 SE:             2.9512E-02  2.1013E-02  3.0846E-05  1.9590E-02  2.6954E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5542E-01  4.9876E-02  2.1363E-01  5.6714E-02  4.7893E-01

 ETASHRINKSD(%)  1.1327E+00  2.9605E+01  9.9897E+01  3.4372E+01  9.6996E+00
 ETASHRINKVR(%)  2.2525E+00  5.0445E+01  1.0000E+02  5.6930E+01  1.8458E+01
 EBVSHRINKSD(%)  1.4243E+00  2.6380E+01  9.9915E+01  4.0294E+01  7.1081E+00
 EBVSHRINKVR(%)  2.8284E+00  4.5801E+01  1.0000E+02  6.4352E+01  1.3711E+01
 RELATIVEINF(%)  9.7088E+01  9.1986E+00  5.0861E-05  6.1150E+00  6.2155E+01
 EPSSHRINKSD(%)  1.6101E+01
 EPSSHRINKVR(%)  2.9610E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          886
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1628.3590808386800     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2801.6464690606977     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1173.2873882220176     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2801.646       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.27E+00  9.14E-01  2.17E-01  1.90E+00  9.77E-01  4.82E-01  1.00E-02  2.42E+00  1.84E+00  2.50E+00
 


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
+        1.13E+03
 
 TH 2
+       -1.97E+01  5.31E+02
 
 TH 3
+        5.32E-01  9.45E+00  9.55E+00
 
 TH 4
+       -2.89E+01  6.78E+02 -2.24E+01  1.50E+03
 
 TH 5
+       -1.77E+00 -3.41E+01 -9.18E+00  3.05E+01  7.85E+01
 
 TH 6
+        4.07E+00 -4.42E+00  4.09E-01 -1.02E+01 -5.43E-01  1.97E+02
 
 TH 7
+        3.39E+00 -5.83E+01  1.30E+00 -1.37E+01 -5.65E+00 -2.61E+00  2.84E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.87E-02 -7.24E+00 -1.92E+00  6.29E+01  7.11E-01 -6.20E-01  1.02E+01  0.00E+00  1.08E+01
 
 TH10
+        4.85E-01 -7.34E+00 -3.83E+00  1.69E+01 -5.90E+00  3.95E-01  4.67E+00  0.00E+00  8.40E-01  4.04E+01
 
 TH11
+       -1.40E+01 -2.38E+01  2.18E-01 -2.53E+01  1.19E+00  3.17E+00  1.83E+01  0.00E+00  1.68E+00  3.50E+00  1.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       49.051
Stop Time:
Wed Sep 29 04:19:17 CDT 2021
