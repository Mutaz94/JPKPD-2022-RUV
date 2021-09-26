Sat Sep 25 13:25:44 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat26.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1177.55754361742        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7393E+01 -4.6303E+01  3.9659E+01 -4.5024E+01  3.1295E+01  3.0572E+01 -6.4134E+00 -1.9013E+02 -1.2897E+01 -5.3998E+00
            -7.5688E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.75745123402        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.2629E-01  1.0764E+00  1.0374E+00  1.0086E+00  1.0846E+00  1.0266E+00  1.0984E+00  1.0985E+00  9.7204E-01  1.0822E+00
             1.1812E+00
 PARAMETER:  2.3427E-02  1.7361E-01  1.3672E-01  1.0856E-01  1.8123E-01  1.2621E-01  1.9385E-01  1.9398E-01  7.1645E-02  1.7897E-01
             2.6654E-01
 GRADIENT:  -1.2987E+02  6.6552E+00 -4.6721E+01  7.1015E+01  4.0891E+01  2.3900E+01  1.2573E+00  8.8279E+00  1.5289E+01  5.3560E+00
             7.4980E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1680.60284072126        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.4618E-01  8.4970E-01  1.4647E+00  1.1839E+00  1.1183E+00  9.6659E-01  1.4421E+00  9.9669E-01  8.5402E-01  1.1951E+00
             1.0906E+00
 PARAMETER:  4.4677E-02 -6.2867E-02  4.8164E-01  2.6877E-01  2.1185E-01  6.6016E-02  4.6610E-01  9.6688E-02 -5.7801E-02  2.7820E-01
             1.8673E-01
 GRADIENT:  -8.4873E+01  4.2042E+01 -1.9224E+00  1.0400E+02  1.2849E+00  8.7399E+00  9.8869E+00 -3.2882E+00  7.9458E+00  4.1158E+00
             4.1934E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1687.15676741161        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.6407E-01  9.4495E-01  1.4085E+00  1.0839E+00  1.1510E+00  9.3784E-01  1.2211E+00  1.0843E+00  9.0374E-01  1.2097E+00
             9.9674E-01
 PARAMETER:  6.3410E-02  4.3378E-02  4.4254E-01  1.8061E-01  2.4065E-01  3.5824E-02  2.9973E-01  1.8092E-01 -1.2136E-03  2.9034E-01
             9.6732E-02
 GRADIENT:  -3.8435E+01  1.3803E+01 -2.0564E+00  4.1030E+01  8.3738E+00  7.5507E-01  3.6465E+00 -1.9867E+00  9.2197E+00  3.3154E+00
             8.9250E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1689.31091942681        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  9.9734E-01  1.0006E+00  1.3273E+00  1.0313E+00  1.1332E+00  9.4371E-01  1.1906E+00  1.1688E+00  8.6798E-01  1.1467E+00
             9.7233E-01
 PARAMETER:  9.7334E-02  1.0064E-01  3.8315E-01  1.3078E-01  2.2509E-01  4.2069E-02  2.7449E-01  2.5595E-01 -4.1591E-02  2.3686E-01
             7.1936E-02
 GRADIENT:   1.7855E+00 -6.4909E-02  3.6834E-01  1.1430E+00 -4.3202E-01  1.0020E-01 -4.9944E-02 -7.2308E-02 -3.3916E-01 -3.6354E-01
            -1.4191E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1689.51431475427        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  9.9698E-01  1.2123E+00  1.1333E+00  8.9424E-01  1.1559E+00  9.4453E-01  1.0613E+00  1.0935E+00  9.3447E-01  1.1472E+00
             9.6907E-01
 PARAMETER:  9.6975E-02  2.9255E-01  2.2513E-01 -1.1784E-02  2.4489E-01  4.2929E-02  1.5953E-01  1.8940E-01  3.2223E-02  2.3732E-01
             6.8577E-02
 GRADIENT:  -3.4016E+00  4.6586E+00  1.8220E+00  3.7818E+00 -3.2326E+00 -2.2602E-01 -2.3926E-01 -2.0395E-01 -7.3009E-01  1.6484E-01
            -1.3301E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1689.53369532643        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      751
 NPARAMETR:  9.9796E-01  1.2221E+00  1.1145E+00  8.8738E-01  1.1595E+00  9.4479E-01  1.0577E+00  1.0912E+00  9.4134E-01  1.1463E+00
             9.7183E-01
 PARAMETER:  9.7953E-02  3.0060E-01  2.0840E-01 -1.9484E-02  2.4797E-01  4.3206E-02  1.5606E-01  1.8731E-01  3.9546E-02  2.3658E-01
             7.1426E-02
 GRADIENT:   4.4393E+01  1.8231E+01 -1.1477E+00  1.1491E+01  3.6853E+00  5.2756E+00  9.2678E-01  3.4686E-01  5.6567E-01  3.9616E-01
             3.0485E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1689.53820010451        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      931
 NPARAMETR:  9.9841E-01  1.2221E+00  1.1167E+00  8.8738E-01  1.1582E+00  9.4509E-01  1.0567E+00  1.0831E+00  9.4189E-01  1.1463E+00
             9.7177E-01
 PARAMETER:  9.8409E-02  3.0060E-01  2.1035E-01 -1.9484E-02  2.4688E-01  4.3529E-02  1.5513E-01  1.7980E-01  4.0138E-02  2.3658E-01
             7.1361E-02
 GRADIENT:  -1.7709E-01  3.0103E+00 -5.1047E-02  4.6387E+00  2.3101E-01 -1.9406E-02  1.6846E-02  4.5870E-02 -2.5787E-02  1.0778E-01
             6.8322E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1689.53821327730        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:      992
 NPARAMETR:  9.9841E-01  1.2221E+00  1.1166E+00  8.8738E-01  1.1582E+00  9.4510E-01  1.0567E+00  1.0829E+00  9.4192E-01  1.1463E+00
             9.7177E-01
 PARAMETER:  9.8413E-02  3.0060E-01  2.1031E-01 -1.9484E-02  2.4686E-01  4.3532E-02  1.5512E-01  1.7962E-01  4.0165E-02  2.3658E-01
             7.1364E-02
 GRADIENT:  -1.3138E+06 -8.7415E+05 -6.2471E+05 -2.6277E+06  5.3223E+05 -1.1642E-02  8.4697E+05 -1.4630E+06  2.6277E+06 -1.1107E+06
             2.6273E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      992
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2591E-04 -1.3873E-02 -3.2514E-02  2.9877E-03 -3.7094E-02
 SE:             2.9834E-02  2.1462E-02  1.2184E-02  2.2079E-02  2.2288E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9663E-01  5.1801E-01  7.6178E-03  8.9236E-01  9.6059E-02

 ETASHRINKSD(%)  5.0741E-02  2.8100E+01  5.9182E+01  2.6033E+01  2.5331E+01
 ETASHRINKVR(%)  1.0146E-01  4.8303E+01  8.3339E+01  4.5288E+01  4.4245E+01
 EBVSHRINKSD(%)  4.4889E-01  2.7797E+01  6.4466E+01  2.7727E+01  2.1471E+01
 EBVSHRINKVR(%)  8.9577E-01  4.7868E+01  8.7373E+01  4.7767E+01  3.8332E+01
 RELATIVEINF(%)  9.8322E+01  1.3739E+00  1.2850E+00  1.3943E+00  1.3367E+01
 EPSSHRINKSD(%)  4.4644E+01
 EPSSHRINKVR(%)  6.9357E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1689.5382132773050     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -954.38738671356680     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.00
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1689.538       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  1.22E+00  1.12E+00  8.87E-01  1.16E+00  9.45E-01  1.06E+00  1.08E+00  9.42E-01  1.15E+00  9.72E-01
 


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
+        6.59E+09
 
 TH 2
+        1.79E+09  4.87E+08
 
 TH 3
+       -3.86E+03  7.61E+08  1.19E+09
 
 TH 4
+       -1.02E+04  2.02E+09 -1.57E+05  8.34E+09
 
 TH 5
+        3.17E+03 -6.25E+08  4.86E+04 -2.59E+09  8.04E+08
 
 TH 6
+       -1.56E+04 -4.23E+03 -6.62E+03 -1.75E+04  5.44E+03  2.14E+02
 
 TH 7
+       -7.77E+03 -2.09E+03 -3.30E+03 -8.75E+03  2.71E+03  9.48E+03  2.45E+09
 
 TH 8
+       -4.67E+03  9.19E+08  1.44E+09  5.37E+05 -1.67E+05 -8.00E+03 -3.99E+03  1.74E+09
 
 TH 9
+        4.81E+04 -1.90E+09  2.04E+04  5.41E+04 -1.68E+04  1.65E+04  4.25E+09  2.47E+04  7.40E+09
 
 TH10
+       -3.35E+03  6.59E+08  1.03E+09  7.13E+03 -2.27E+03 -5.73E+03 -2.86E+03  1.76E+05  1.77E+04  8.93E+08
 
 TH11
+        2.10E+06 -1.84E+09  8.91E+05  2.36E+06 -7.32E+05  1.60E+04  7.98E+03  1.08E+06 -4.94E+04  7.71E+05  6.95E+09
 
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
 #CPUT: Total CPU Time in Seconds,       17.826
Stop Time:
Sat Sep 25 13:26:04 CDT 2021
