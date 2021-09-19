Sat Sep 18 08:20:58 CDT 2021
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
$DATA ../../../../data/spa/B/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1739.44513721827        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.1197E+01 -4.9543E+01 -8.3569E+00 -4.1434E+01  2.4694E+00  1.8993E+01  6.5501E+00  7.0860E+00  4.4181E+01 -3.3663E+00
             2.7412E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1749.09514556975        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0189E+00  1.0902E+00  1.0032E+00  9.7195E-01  1.0168E+00  9.6455E-01  9.8651E-01  9.5033E-01  7.5872E-01  1.0873E+00
             9.7391E-01
 PARAMETER:  1.1871E-01  1.8634E-01  1.0319E-01  7.1549E-02  1.1670E-01  6.3902E-02  8.6419E-02  4.9056E-02 -1.7613E-01  1.8369E-01
             7.3561E-02
 GRADIENT:   2.8283E+01  6.0974E-01  1.6900E+01 -2.3910E+01 -3.8178E+01  7.0063E+00 -4.2498E+00  2.9389E+00  1.4318E+00  3.5931E+00
             1.4042E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1750.22985556339        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0111E+00  8.8974E-01  9.8270E-01  1.1046E+00  9.4224E-01  9.8977E-01  1.2977E+00  5.8372E-01  6.4379E-01  1.0892E+00
             9.5563E-01
 PARAMETER:  1.1107E-01 -1.6830E-02  8.2548E-02  1.9947E-01  4.0503E-02  8.9713E-02  3.6056E-01 -4.3834E-01 -3.4039E-01  1.8544E-01
             5.4619E-02
 GRADIENT:   1.2572E+01  1.1720E+01  8.8764E+00  2.5925E+01 -1.4490E+01  1.7350E+01  2.0439E+00 -5.1054E-01 -1.2366E+00  4.7670E+00
             6.4031E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1750.74782728263        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0099E+00  1.0107E+00  8.3951E-01  1.0174E+00  9.2564E-01  9.5691E-01  1.1360E+00  4.3182E-01  6.9989E-01  1.0195E+00
             9.4455E-01
 PARAMETER:  1.0989E-01  1.1061E-01 -7.4934E-02  1.1722E-01  2.2728E-02  5.5955E-02  2.2748E-01 -7.3976E-01 -2.5683E-01  1.1936E-01
             4.2958E-02
 GRADIENT:   2.3621E+00  9.5031E-01 -1.4098E+00  7.3265E+00  5.4991E-01  3.2581E+00  2.7057E-01  3.5122E-01  1.8567E-01  1.1377E+00
             2.5879E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1751.42942340166        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  1.0349E+00  1.1226E+00  7.7266E-01  9.5267E-01  9.3710E-01  9.6118E-01  1.0326E+00  3.0872E-01  7.5010E-01  1.0121E+00
             9.4164E-01
 PARAMETER:  1.3430E-01  2.1563E-01 -1.5791E-01  5.1509E-02  3.5034E-02  6.0409E-02  1.3212E-01 -1.0753E+00 -1.8755E-01  1.1207E-01
             3.9870E-02
 GRADIENT:   3.7691E+00  5.8756E+00  5.6798E-01  7.6110E+00 -3.0177E+00  1.7905E-01 -4.4173E-01  1.3797E-01  1.4505E-01  1.9480E-01
             7.5631E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1751.52348408309        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      587
 NPARAMETR:  1.0341E+00  1.2570E+00  7.0878E-01  8.6314E-01  9.7601E-01  9.6163E-01  9.3483E-01  1.6041E-01  8.0665E-01  1.0276E+00
             9.3941E-01
 PARAMETER:  1.3350E-01  3.2873E-01 -2.4421E-01 -4.7183E-02  7.5714E-02  6.0875E-02  3.2609E-02 -1.7300E+00 -1.1486E-01  1.2720E-01
             3.7501E-02
 GRADIENT:   2.2938E-01 -1.0851E+00  5.1074E-02 -1.1027E+00  6.6968E-01 -9.9972E-02 -1.9834E-01  4.1832E-02 -3.8646E-02 -2.3967E-01
            -1.6871E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1751.54435428781        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      766
 NPARAMETR:  1.0339E+00  1.2264E+00  7.1334E-01  8.8274E-01  9.6026E-01  9.6167E-01  9.5710E-01  6.2838E-02  7.9258E-01  1.0193E+00
             9.3958E-01
 PARAMETER:  1.3337E-01  3.0409E-01 -2.3779E-01 -2.4728E-02  5.9453E-02  6.0917E-02  5.6151E-02 -2.6672E+00 -1.3247E-01  1.1911E-01
             3.7679E-02
 GRADIENT:   1.1522E-02  8.9815E-01  3.4709E-01  6.9422E-01 -2.8574E-01 -2.1825E-02 -3.4284E-02  5.0567E-03 -9.9182E-02 -1.2378E-01
            -8.8470E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1751.54813669989        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      941
 NPARAMETR:  1.0340E+00  1.2378E+00  7.0695E-01  8.7493E-01  9.6304E-01  9.6179E-01  9.4897E-01  1.0000E-02  7.9814E-01  1.0200E+00
             9.3958E-01
 PARAMETER:  1.3344E-01  3.1336E-01 -2.4680E-01 -3.3615E-02  6.2336E-02  6.1038E-02  4.7619E-02 -4.6008E+00 -1.2547E-01  1.1985E-01
             3.7677E-02
 GRADIENT:   3.1251E-02 -7.9110E-02 -3.4846E-03 -7.2017E-02  3.3230E-02 -1.9530E-02  1.6024E-03  0.0000E+00 -4.0415E-03 -3.7461E-03
            -9.0201E-04

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1751.54813950410        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  1.0340E+00  1.2381E+00  7.0679E-01  8.7478E-01  9.6306E-01  9.6184E-01  9.4875E-01  1.0000E-02  7.9829E-01  1.0200E+00
             9.3958E-01
 PARAMETER:  1.3343E-01  3.1359E-01 -2.4702E-01 -3.3780E-02  6.2358E-02  6.1088E-02  4.7392E-02 -4.5988E+00 -1.2529E-01  1.1983E-01
             3.7675E-02
 GRADIENT:   1.4786E-03 -9.0553E-03 -2.6164E-03 -5.2612E-03  2.6261E-03 -1.5437E-03 -2.6496E-04  0.0000E+00 -9.3973E-05  1.1149E-03
             1.9480E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      998
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.0909E-04 -1.1033E-02 -4.1865E-04  4.7358E-03 -2.0629E-02
 SE:             2.9841E-02  2.2644E-02  1.6460E-04  2.2352E-02  2.3828E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9441E-01  6.2609E-01  1.0978E-02  8.3221E-01  3.8665E-01

 ETASHRINKSD(%)  2.9839E-02  2.4140E+01  9.9449E+01  2.5119E+01  2.0172E+01
 ETASHRINKVR(%)  5.9669E-02  4.2453E+01  9.9997E+01  4.3928E+01  3.6274E+01
 EBVSHRINKSD(%)  4.0913E-01  2.3858E+01  9.9526E+01  2.6423E+01  1.7899E+01
 EBVSHRINKVR(%)  8.1659E-01  4.2024E+01  9.9998E+01  4.5864E+01  3.2594E+01
 RELATIVEINF(%)  9.8853E+01  2.3814E+00  2.1897E-04  2.1678E+00  8.1215E+00
 EPSSHRINKSD(%)  4.3811E+01
 EPSSHRINKVR(%)  6.8428E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1751.5481395041027     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1016.3973129403645     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.25
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1751.548       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.24E+00  7.07E-01  8.75E-01  9.63E-01  9.62E-01  9.49E-01  1.00E-02  7.98E-01  1.02E+00  9.40E-01
 


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
+        1.11E+03
 
 TH 2
+       -6.95E+00  4.84E+02
 
 TH 3
+        1.77E+01  1.70E+02  5.73E+02
 
 TH 4
+       -1.41E+01  4.96E+02 -3.77E+02  1.23E+03
 
 TH 5
+       -2.34E+00 -2.45E+02 -5.39E+02  3.42E+02  7.73E+02
 
 TH 6
+       -5.56E-01 -1.61E+00  6.00E+00 -5.74E+00  9.05E-01  2.11E+02
 
 TH 7
+        3.25E+00  2.74E+01 -1.38E+01 -1.32E+01 -4.69E+00 -2.40E+00  7.96E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.73E-01 -2.00E+01 -3.39E+01  3.34E+01  1.27E+01  7.17E-01  3.35E+01  0.00E+00  1.00E+02
 
 TH10
+       -5.40E-01 -1.59E+01 -5.50E+01 -1.53E+00 -5.54E+01  3.05E-01  8.26E+00  0.00E+00  1.26E+01  8.51E+01
 
 TH11
+       -7.99E+00 -1.74E+01 -4.28E+01  3.41E+00  7.33E+00  1.03E+00  6.30E+00  0.00E+00  1.37E+01  1.79E+01  2.42E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.171
Stop Time:
Sat Sep 18 08:21:16 CDT 2021
