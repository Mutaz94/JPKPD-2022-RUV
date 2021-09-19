Sat Sep 18 09:37:39 CDT 2021
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
$DATA ../../../../data/spa/A2/dat6.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -834.326360370397        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3800E+01 -4.1150E+01 -3.5015E+01 -1.3759E+01  2.3169E+02  1.3212E+01 -3.8701E+01 -4.0826E+00 -7.2769E+01 -8.6048E+01
            -1.4605E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1359.24162326744        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0018E+00  9.4701E-01  1.0693E+00  1.0785E+00  8.8929E-01  8.8251E-01  9.9040E-01  9.1943E-01  1.1601E+00  8.4467E-01
             2.1867E+00
 PARAMETER:  1.0177E-01  4.5552E-02  1.6701E-01  1.7559E-01 -1.7330E-02 -2.4989E-02  9.0353E-02  1.5993E-02  2.4853E-01 -6.8810E-02
             8.8238E-01
 GRADIENT:  -2.1367E+00  2.6395E+00 -4.5680E+00  1.2233E+01  5.5270E+01 -2.8766E+01 -9.5408E-01  2.3409E+00 -1.5429E+00 -1.8282E+01
            -1.5601E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1366.02712126712        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0012E+00  7.3569E-01  5.5905E-01  1.1722E+00  5.5186E-01  8.8530E-01  6.4217E-01  4.2222E-01  1.1827E+00  5.1630E-01
             2.2402E+00
 PARAMETER:  1.0117E-01 -2.0695E-01 -4.8152E-01  2.5891E-01 -4.9446E-01 -2.1834E-02 -3.4291E-01 -7.6224E-01  2.6777E-01 -5.6107E-01
             9.0658E-01
 GRADIENT:  -1.9912E+01  1.4790E+01 -2.0681E+01  5.1560E+01  4.3927E+01 -3.0336E+01 -1.1272E+01 -1.3735E+00  6.3746E+00 -2.2031E+01
            -1.4171E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1388.11316908272        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0132E+00  5.9132E-01  5.3474E-01  1.2351E+00  4.8907E-01  9.5845E-01  1.1755E+00  1.6356E-01  9.7972E-01  5.0909E-01
             2.7220E+00
 PARAMETER:  1.1310E-01 -4.2540E-01 -5.2598E-01  3.1113E-01 -6.1525E-01  5.7561E-02  2.6166E-01 -1.7106E+00  7.9510E-02 -5.7513E-01
             1.1014E+00
 GRADIENT:  -7.5644E+00  1.6067E+01 -2.9435E+00  3.3020E+01  2.5434E+00  3.6577E+00 -2.9347E+00  2.0896E-01 -5.2149E+00 -2.2066E+00
             3.2924E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1391.45149476262        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.0157E+00  3.5140E-01  4.3858E-01  1.3018E+00  3.7356E-01  9.5132E-01  1.7113E+00  1.4436E-02  9.6147E-01  5.7492E-01
             2.6152E+00
 PARAMETER:  1.1560E-01 -9.4584E-01 -7.2421E-01  3.6372E-01 -8.8469E-01  5.0100E-02  6.3724E-01 -4.1380E+00  6.0707E-02 -4.5352E-01
             1.0613E+00
 GRADIENT:  -9.2532E-01  1.5096E+01  1.7251E+01  2.2882E+01 -3.5069E+01 -1.2348E+00  2.3344E+00  1.9111E-03 -1.6436E-02  1.8889E+00
             6.4248E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1393.40406091135        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0097E+00  1.5694E-01  4.0632E-01  1.3448E+00  3.2950E-01  9.5382E-01  1.9193E+00  1.0000E-02  9.4256E-01  5.7919E-01
             2.6032E+00
 PARAMETER:  1.0970E-01 -1.7519E+00 -8.0062E-01  3.9625E-01 -1.0102E+00  5.2722E-02  7.5198E-01 -8.9653E+00  4.0848E-02 -4.4613E-01
             1.0567E+00
 GRADIENT:  -4.5257E-01  3.0416E+00  1.3686E+01  2.8162E+00 -2.7495E+01  5.8393E-01 -1.6873E+00  0.0000E+00 -2.3362E+00 -4.6348E-01
            -3.1131E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1396.03512901975        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  1.0016E+00  4.1218E-02  4.6214E-01  1.4279E+00  3.5081E-01  9.4825E-01  3.9233E+00  1.0000E-02  9.1047E-01  5.7746E-01
             2.6287E+00
 PARAMETER:  1.0158E-01 -3.0889E+00 -6.7190E-01  4.5618E-01 -9.4751E-01  4.6865E-02  1.4669E+00 -1.8249E+01  6.2096E-03 -4.4911E-01
             1.0665E+00
 GRADIENT:  -1.9066E+00  3.5852E-01  3.5054E+00  4.9215E+00 -5.8242E+00  6.8203E-01 -3.3381E-01  0.0000E+00 -6.0840E-01 -3.9547E-01
            -5.6597E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1396.25724876142        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      510
 NPARAMETR:  1.0008E+00  1.0000E-02  4.5558E-01  1.4330E+00  3.4475E-01  9.4547E-01  1.1747E+01  1.0000E-02  9.0804E-01  5.7916E-01
             2.6273E+00
 PARAMETER:  1.0083E-01 -4.6481E+00 -6.8618E-01  4.5977E-01 -9.6493E-01  4.3924E-02  2.5636E+00 -2.9760E+01  3.5301E-03 -4.4618E-01
             1.0660E+00
 GRADIENT:   1.5919E-01  0.0000E+00 -2.9895E-01 -4.3226E-01  4.7346E-01 -8.5597E-02  6.4860E-04  0.0000E+00  3.7628E-02  4.4124E-02
            -3.1768E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1396.70055250363        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  1.0025E+00  1.0000E-02  5.2092E-01  1.4775E+00  3.7933E-01  9.4423E-01  1.1394E+01  1.0000E-02  8.9219E-01  5.8896E-01
             2.6472E+00
 PARAMETER:  1.0248E-01 -4.6721E+00 -5.5216E-01  4.9035E-01 -8.6936E-01  4.2610E-02  2.5331E+00 -2.9857E+01 -1.4074E-02 -4.2939E-01
             1.0735E+00
 GRADIENT:   3.7760E-01  0.0000E+00 -2.2842E-01 -2.3593E-01  4.1750E-01 -5.7771E-02 -3.0039E-02  0.0000E+00 -7.8974E-03 -7.8678E-03
             9.6336E-03

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1396.70114735586        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:      799
 NPARAMETR:  1.0023E+00  1.0000E-02  5.2063E-01  1.4775E+00  3.7912E-01  9.4439E-01  1.1738E+01  1.0000E-02  8.9228E-01  5.8900E-01
             2.6470E+00
 PARAMETER:  1.0233E-01 -4.7134E+00 -5.5272E-01  4.9034E-01 -8.6991E-01  4.2789E-02  2.5628E+00 -3.0163E+01 -1.3980E-02 -4.2932E-01
             1.0734E+00
 GRADIENT:   1.5749E-03  0.0000E+00 -4.4889E-04  6.2629E-04  4.4297E-04 -1.3557E-03 -7.8663E-04  0.0000E+00  6.8353E-04 -3.9336E-04
             2.5949E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      799
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.9578E-04  1.0249E-03  7.3633E-05 -1.0154E-02 -6.3097E-03
 SE:             2.8997E-02  1.9343E-03  2.2476E-04  2.7078E-02  1.8903E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9186E-01  5.9622E-01  7.4321E-01  7.0766E-01  7.3853E-01

 ETASHRINKSD(%)  2.8561E+00  9.3520E+01  9.9247E+01  9.2851E+00  3.6674E+01
 ETASHRINKVR(%)  5.6307E+00  9.9580E+01  9.9994E+01  1.7708E+01  5.9898E+01
 EBVSHRINKSD(%)  2.7809E+00  9.4133E+01  9.9225E+01  8.7348E+00  3.6323E+01
 EBVSHRINKVR(%)  5.4845E+00  9.9656E+01  9.9994E+01  1.6707E+01  5.9452E+01
 RELATIVEINF(%)  8.1295E+01  1.8661E-02  2.4285E-04  8.4793E+00  1.3543E+00
 EPSSHRINKSD(%)  3.0683E+01
 EPSSHRINKVR(%)  5.1951E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1396.7011473558593     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -661.55032079212117     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.54
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1396.701       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  5.21E-01  1.48E+00  3.79E-01  9.44E-01  1.17E+01  1.00E-02  8.92E-01  5.89E-01  2.65E+00
 


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
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.11E+01  0.00E+00  2.42E+03
 
 TH 4
+       -4.19E+01  0.00E+00 -1.56E+02  5.35E+02
 
 TH 5
+        1.20E+02  0.00E+00 -4.41E+03 -2.29E+02  8.89E+03
 
 TH 6
+        1.39E+00  0.00E+00  1.37E+01 -1.07E+01  2.40E+00  1.91E+02
 
 TH 7
+       -1.70E-02  0.00E+00  7.81E-03 -1.98E-02  3.63E-02  9.37E-03  5.08E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.86E+00  0.00E+00  3.67E+01 -1.24E+01  9.83E+00  7.45E+00  1.85E-02  0.00E+00  1.83E+02
 
 TH10
+       -9.92E+00  0.00E+00 -5.23E+01 -3.31E+00  5.53E+01  2.79E+00  1.81E-02  0.00E+00  5.10E+00  1.00E+02
 
 TH11
+       -1.47E+01  0.00E+00 -6.01E+00 -6.35E+00 -7.65E+00  3.80E+00 -1.54E-03  0.00E+00  7.82E+00  2.73E+01  4.40E+01
 
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
 #CPUT: Total CPU Time in Seconds,       13.831
Stop Time:
Sat Sep 18 09:37:54 CDT 2021
