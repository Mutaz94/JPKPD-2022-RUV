Fri Sep 24 22:25:22 CDT 2021
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
$DATA ../../../../data/int/A3/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -620.941417672164        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2500E+02  1.1112E+02  1.8920E+02 -1.0014E+02  1.3521E+02  9.9908E+00 -1.1860E+02 -1.4574E+02 -4.4565E+01 -1.0434E+02
            -6.0036E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2766.77755427483        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7786E-01  1.0274E+00  1.0155E+00  1.0854E+00  9.9523E-01  8.9481E-01  7.9793E-01  9.2280E-01  9.7933E-01  7.6096E-01
             2.6581E+00
 PARAMETER:  7.7613E-02  1.2701E-01  1.1539E-01  1.8198E-01  9.5217E-02 -1.1147E-02 -1.2574E-01  1.9653E-02  7.9109E-02 -1.7318E-01
             1.0776E+00
 GRADIENT:   1.0761E+00  1.6395E+01  1.0158E+01  1.3579E+01 -1.5709E+00 -2.1978E+01  4.8764E+00  4.1207E+00 -7.9173E+00 -3.1174E+00
            -3.0293E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2771.06838902890        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8422E-01  9.1399E-01  8.1359E-01  1.1447E+00  8.5198E-01  9.3078E-01  6.5107E-01  4.2516E-01  1.0211E+00  8.5001E-01
             2.6568E+00
 PARAMETER:  8.4095E-02  1.0070E-02 -1.0629E-01  2.3513E-01 -6.0192E-02  2.8262E-02 -3.2913E-01 -7.5528E-01  1.2090E-01 -6.2509E-02
             1.0771E+00
 GRADIENT:   1.8445E+01  1.5491E+01 -1.8459E+01  2.7228E+01  1.3265E+01 -6.3319E+00 -3.0824E-01  1.9146E+00  2.7926E+00  3.9018E+00
             3.8449E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2772.97789743798        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.7484E-01  7.0460E-01  6.7580E-01  1.2268E+00  6.5916E-01  9.4992E-01  7.6469E-01  1.7567E-01  9.4236E-01  7.5927E-01
             2.6427E+00
 PARAMETER:  7.4518E-02 -2.5012E-01 -2.9186E-01  3.0442E-01 -3.1679E-01  4.8618E-02 -1.6828E-01 -1.6391E+00  4.0630E-02 -1.7540E-01
             1.0718E+00
 GRADIENT:  -3.6294E+00  6.7983E+00  3.6584E+00  1.0650E+01 -2.4760E+00  1.2135E+00  3.4233E-02  3.5510E-01 -5.0577E+00  4.0388E-01
             1.1533E-03

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2773.96939287827        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.6998E-01  4.5179E-01  4.3211E-01  1.3300E+00  4.1172E-01  9.6620E-01  8.8452E-01  2.8297E-02  8.7897E-01  6.5825E-01
             2.6126E+00
 PARAMETER:  6.9516E-02 -6.9454E-01 -7.3908E-01  3.8522E-01 -7.8740E-01  6.5613E-02 -2.2715E-02 -3.4650E+00 -2.9003E-02 -3.1817E-01
             1.0603E+00
 GRADIENT:  -1.6875E+01  3.8419E+01  6.0503E+01  1.3740E+02 -5.8380E+01  5.8504E+00 -8.1206E+00  3.4913E-03 -4.6900E+01 -1.0643E+01
             2.0911E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2774.08660668935        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.7143E-01  4.1481E-01  3.9496E-01  1.3348E+00  3.7801E-01  9.6392E-01  8.9039E-01  1.7861E-02  8.6965E-01  6.3080E-01
             2.6260E+00
 PARAMETER:  7.1012E-02 -7.7992E-01 -8.2897E-01  3.8882E-01 -8.7284E-01  6.3250E-02 -1.6100E-02 -3.9251E+00 -3.9662E-02 -3.6076E-01
             1.0655E+00
 GRADIENT:  -1.4937E+01  3.3015E+01  7.2782E+01  1.6108E+02 -7.6050E+01  4.8883E+00 -1.0625E+01 -5.9494E-04 -5.9037E+01 -1.6869E+01
             4.1664E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2775.32536740842        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  9.7590E-01  3.5816E-01  3.3373E-01  1.3227E+00  3.2741E-01  9.5524E-01  8.9422E-01  1.0000E-02  8.7426E-01  5.9146E-01
             2.6550E+00
 PARAMETER:  7.5609E-02 -9.2679E-01 -9.9742E-01  3.7969E-01 -1.0166E+00  5.4209E-02 -1.1806E-02 -5.0794E+00 -3.4381E-02 -4.2516E-01
             1.0764E+00
 GRADIENT:  -7.0637E+00  1.4062E+01  7.5820E+01  1.7737E+02 -8.7808E+01  1.4341E+00 -1.4189E+01  0.0000E+00 -7.2114E+01 -2.5831E+01
             8.1936E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2789.01333173939        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  9.8819E-01  2.4558E-01  2.0361E-01  1.1779E+00  2.2756E-01  9.2838E-01  9.5743E-01  1.0000E-02  1.0666E+00  5.7066E-01
             2.6492E+00
 PARAMETER:  8.8115E-02 -1.3041E+00 -1.4916E+00  2.6377E-01 -1.3803E+00  2.5682E-02  5.6498E-02 -1.0405E+01  1.6448E-01 -4.6097E-01
             1.0743E+00
 GRADIENT:   2.2724E+01 -4.6146E+01  4.5086E+01  9.1689E+01 -5.3161E+01 -1.0071E+01 -1.6033E+01  0.0000E+00 -4.6660E+01 -3.2799E+01
             1.3854E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2800.81450911447        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      611
 NPARAMETR:  9.7999E-01  2.1130E-01  1.4939E-01  1.0378E+00  1.8972E-01  9.4723E-01  1.1200E+00  1.0000E-02  1.3242E+00  6.7315E-01
             2.4684E+00
 PARAMETER:  7.9785E-02 -1.4545E+00 -1.8012E+00  1.3706E-01 -1.5622E+00  4.5785E-02  2.1335E-01 -1.4739E+01  3.8079E-01 -2.9579E-01
             1.0036E+00
 GRADIENT:   4.9634E+00 -1.9763E+01 -2.8919E+00  1.4490E+01 -4.4203E+01 -5.0997E+00 -3.9374E+00  0.0000E+00 -7.3357E+00 -7.4218E+00
             4.6289E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2807.68278270836        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  9.7614E-01  2.6410E-01  2.0542E-01  1.1275E+00  2.3727E-01  9.5972E-01  1.1179E+00  1.0000E-02  1.2234E+00  6.7113E-01
             2.4238E+00
 PARAMETER:  7.5854E-02 -1.2314E+00 -1.4827E+00  2.2001E-01 -1.3386E+00  5.8883E-02  2.1147E-01 -1.1450E+01  3.0167E-01 -2.9879E-01
             9.8532E-01
 GRADIENT:  -1.5444E-01  2.7231E-01 -3.2806E-01 -1.1456E-01  2.8128E-01  1.1767E-02 -2.6561E-02  0.0000E+00 -1.4233E-01  1.8759E-01
             3.0987E-01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -2807.68301905483        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.7620E-01  2.6392E-01  2.0531E-01  1.1274E+00  2.3715E-01  9.5970E-01  1.1183E+00  1.0000E-02  1.2242E+00  6.7059E-01
             2.4234E+00
 PARAMETER:  7.5910E-02 -1.2321E+00 -1.4832E+00  2.1996E-01 -1.3391E+00  5.8870E-02  2.1182E-01 -1.1456E+01  3.0231E-01 -2.9960E-01
             9.8517E-01
 GRADIENT:   1.4319E-04  6.3762E-03 -1.0912E-02 -2.8270E-03  3.8103E-03  2.0920E-04 -1.3554E-03  0.0000E+00  1.0956E-03 -1.1107E-04
            -1.1859E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      879
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6941E-03  5.6468E-03  8.6664E-05 -6.0909E-03  3.4062E-03
 SE:             2.9304E-02  2.2740E-02  3.0479E-04  2.7892E-02  2.5618E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5390E-01  8.0388E-01  7.7615E-01  8.2714E-01  8.9423E-01

 ETASHRINKSD(%)  1.8282E+00  2.3820E+01  9.8979E+01  6.5580E+00  1.4175E+01
 ETASHRINKVR(%)  3.6229E+00  4.1966E+01  9.9990E+01  1.2686E+01  2.6341E+01
 EBVSHRINKSD(%)  1.7659E+00  2.2888E+01  9.9081E+01  5.1522E+00  1.4652E+01
 EBVSHRINKVR(%)  3.5005E+00  4.0538E+01  9.9992E+01  1.0039E+01  2.7157E+01
 RELATIVEINF(%)  9.6479E+01  1.0225E+01  5.4352E-04  4.8791E+01  4.2139E+00
 EPSSHRINKSD(%)  1.8614E+01
 EPSSHRINKVR(%)  3.3763E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2807.6830190548349     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1153.5936592864241     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2807.683       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  2.64E-01  2.05E-01  1.13E+00  2.37E-01  9.60E-01  1.12E+00  1.00E-02  1.22E+00  6.71E-01  2.42E+00
 


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
+        1.23E+03
 
 TH 2
+        7.74E-01  7.69E+03
 
 TH 3
+        5.23E+01 -1.70E+03  1.96E+04
 
 TH 4
+       -3.72E+00 -7.32E+01 -4.44E+02  5.27E+02
 
 TH 5
+       -1.32E+01 -5.73E+03 -2.03E+04 -2.58E+02  3.22E+04
 
 TH 6
+        8.59E+00 -1.08E+01  4.00E+01 -8.62E+00 -1.88E+01  1.97E+02
 
 TH 7
+       -2.96E+00  3.92E+01  6.00E+01 -2.09E+00 -7.54E+01 -1.55E+00  5.58E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.80E+00 -9.40E+00  1.75E+02 -7.12E+00  1.82E+01 -2.31E+00  6.68E-01  0.00E+00  1.02E+02
 
 TH10
+       -2.33E+00 -2.20E+01  8.00E+01  9.50E+00 -2.04E+01 -1.67E-01  8.13E-01  0.00E+00  2.84E+00  2.46E+02
 
 TH11
+       -1.80E+01 -5.39E+00 -1.48E+02 -4.88E+00  1.25E+02  3.73E+00  1.51E+01  0.00E+00  5.32E+00  1.79E+01  1.83E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.785
Stop Time:
Fri Sep 24 22:25:54 CDT 2021
