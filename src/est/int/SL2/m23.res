Sat Sep 25 01:02:30 CDT 2021
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
$DATA ../../../../data/int/SL2/dat23.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      997
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

 TOT. NO. OF OBS RECS:      897
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1645.67404281510        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0649E+01  1.1623E+02  5.3983E+01  9.0206E+01  5.4283E+01  1.7392E+01 -1.0879E+02 -1.5066E+02 -4.5242E+01 -7.8671E+00
            -4.1227E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2924.40261016820        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0184E+00  1.0825E+00  9.8441E-01  9.2878E-01  1.0563E+00  9.2533E-01  1.0645E+00  1.0230E+00  8.6503E-01  9.6774E-01
             2.0315E+00
 PARAMETER:  1.1819E-01  1.7928E-01  8.4292E-02  2.6112E-02  1.5475E-01  2.2397E-02  1.6248E-01  1.2273E-01 -4.4985E-02  6.7205E-02
             8.0878E-01
 GRADIENT:  -5.3500E+00  9.5725E+00 -9.3056E+00  2.6864E+00  9.4465E+00 -1.2145E+01 -5.6204E+00 -7.6990E+00 -8.2664E+00 -6.5615E+00
            -2.7735E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2932.01388930738        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0199E+00  1.2751E+00  1.3071E+00  8.4897E-01  1.3043E+00  9.4873E-01  9.8688E-01  2.0949E+00  8.6438E-01  1.2001E+00
             2.0780E+00
 PARAMETER:  1.1971E-01  3.4306E-01  3.6784E-01 -6.3737E-02  3.6570E-01  4.7371E-02  8.6793E-02  8.3952E-01 -4.5741E-02  2.8240E-01
             8.3139E-01
 GRADIENT:  -3.3566E+00  3.2983E+01 -4.9699E+00  3.9259E+01  1.4352E+01 -2.2122E+00 -4.5962E-01  2.8244E+00 -2.1711E+00  8.2651E+00
            -1.9260E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2941.46235866194        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0237E+00  1.1741E+00  1.2149E+00  8.8594E-01  1.1927E+00  9.5100E-01  1.0282E+00  1.5507E+00  8.6530E-01  1.0589E+00
             2.2615E+00
 PARAMETER:  1.2341E-01  2.6050E-01  2.9466E-01 -2.1109E-02  2.7624E-01  4.9758E-02  1.2778E-01  5.3870E-01 -4.4683E-02  1.5724E-01
             9.1601E-01
 GRADIENT:   6.1013E-02  1.2299E-01 -4.3399E-01  1.6157E-01  4.6786E-02 -1.0814E-01 -1.1129E-01  2.2023E-02  2.1330E-01  3.9690E-02
            -2.0424E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2945.04573264904        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.0236E+00  1.2146E+00  1.8686E+00  8.7705E-01  1.3876E+00  9.5005E-01  1.0495E+00  2.6829E+00  7.7915E-01  1.0929E+00
             2.2329E+00
 PARAMETER:  1.2337E-01  2.9444E-01  7.2520E-01 -3.1190E-02  4.2761E-01  4.8762E-02  1.4831E-01  1.0869E+00 -1.4955E-01  1.8882E-01
             9.0330E-01
 GRADIENT:   1.6746E+00 -1.0213E-01 -6.8750E+00 -6.9063E-01  1.2674E+01  3.4704E-02 -2.9901E+00 -9.9529E-01  3.6479E+00 -9.7219E-01
            -4.2190E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2945.08187957274        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  1.0236E+00  1.2128E+00  1.8734E+00  8.7826E-01  1.3869E+00  9.5004E-01  1.0518E+00  2.6885E+00  7.7716E-01  1.0918E+00
             2.2328E+00
 PARAMETER:  1.2337E-01  2.9291E-01  7.2777E-01 -2.9816E-02  4.2710E-01  4.8748E-02  1.5048E-01  1.0890E+00 -1.5211E-01  1.8785E-01
             9.0327E-01
 GRADIENT:  -1.0610E+01 -5.8308E+00 -7.3716E+00 -2.6052E+00  9.2190E+00 -1.0823E+00 -3.2305E+00 -1.4290E+00  3.4783E+00 -1.1400E+00
            -5.6521E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2945.09653985324        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      628
 NPARAMETR:  1.0236E+00  1.2215E+00  1.8734E+00  8.7821E-01  1.3866E+00  9.5420E-01  1.0524E+00  2.6885E+00  7.7716E-01  1.0918E+00
             2.2328E+00
 PARAMETER:  1.2337E-01  3.0005E-01  7.2777E-01 -2.9867E-02  4.2689E-01  5.3115E-02  1.5110E-01  1.0890E+00 -1.5211E-01  1.8785E-01
             9.0327E-01
 GRADIENT:   1.6084E+00  8.9833E+00 -6.2990E+00  6.3712E+00  9.7624E+00  1.6363E+00 -2.3746E+00 -1.0640E+00  3.4512E+00 -1.1302E+00
            -4.6414E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2945.14309144608        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      771
 NPARAMETR:  1.0322E+00  1.2214E+00  2.0562E+00  8.7826E-01  1.3863E+00  9.5142E-01  1.0837E+00  2.6871E+00  7.7722E-01  1.0917E+00
             2.2318E+00
 PARAMETER:  1.3173E-01  3.0002E-01  8.2086E-01 -2.9817E-02  4.2667E-01  5.0197E-02  1.8037E-01  1.0884E+00 -1.5203E-01  1.8776E-01
             9.0282E-01
 GRADIENT:   9.7108E+00  1.0089E+01  5.9412E+00 -9.4799E+00 -1.1974E+01 -4.7358E-01  1.8086E+00 -6.0611E+00  4.7571E+00 -1.1005E+00
            -7.1642E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2945.16071581324        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      933
 NPARAMETR:  1.0323E+00  1.2215E+00  2.0566E+00  8.7828E-01  1.3862E+00  9.5248E-01  1.0713E+00  2.6863E+00  7.7725E-01  1.0917E+00
             2.2394E+00
 PARAMETER:  1.3176E-01  3.0009E-01  8.2106E-01 -2.9792E-02  4.2657E-01  5.1310E-02  1.6884E-01  1.0882E+00 -1.5199E-01  1.8771E-01
             9.0622E-01
 GRADIENT:   9.5290E+00  9.3317E+00  6.0307E+00 -8.8662E+00 -1.2414E+01  6.2582E-03 -7.2051E-03 -5.8141E+00  4.3080E+00 -1.0686E+00
            -5.4392E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2945.37161528437        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  1.0313E+00  1.2195E+00  2.0404E+00  8.7927E-01  1.3888E+00  9.5237E-01  1.0719E+00  2.7216E+00  7.6685E-01  1.0932E+00
             2.2393E+00
 PARAMETER:  1.3080E-01  2.9842E-01  8.1315E-01 -2.8664E-02  4.2843E-01  5.1199E-02  1.6939E-01  1.1012E+00 -1.6546E-01  1.8913E-01
             9.0614E-01
 GRADIENT:   7.1522E+00  7.4214E+00  3.7158E+00 -7.3639E+00 -8.4098E+00 -4.1973E-02 -6.6680E-01 -4.4916E+00  3.5387E+00 -8.9462E-01
             7.3406E-01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -2945.37161528437        NO. OF FUNC. EVALS.:  39
 CUMULATIVE NO. OF FUNC. EVALS.:     1098
 NPARAMETR:  1.0313E+00  1.2195E+00  2.0402E+00  8.7927E-01  1.3888E+00  9.5237E-01  1.0719E+00  2.7216E+00  7.6672E-01  1.0932E+00
             2.2393E+00
 PARAMETER:  1.3080E-01  2.9842E-01  8.1315E-01 -2.8664E-02  4.2843E-01  5.1199E-02  1.6939E-01  1.1012E+00 -1.6546E-01  1.8913E-01
             9.0614E-01
 GRADIENT:   7.0664E+00 -1.0297E+05  3.6889E+00 -1.5365E+05  7.1714E+04 -3.4464E-02 -6.5967E-01  2.7847E+04  3.5341E+00 -8.8785E-01
             7.3442E-01
 NUMSIGDIG:         1.4         3.3         1.2         3.3         3.3         2.8         1.3         3.3         0.3         1.2
                    3.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1098
 NO. OF SIG. DIGITS IN FINAL EST.:  0.3
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.8087E-03 -1.6254E-02 -3.4873E-02  1.1009E-02 -2.6370E-02
 SE:             2.9587E-02  2.3536E-02  1.9726E-02  1.9781E-02  2.2050E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5126E-01  4.8981E-01  7.7080E-02  5.7783E-01  2.3173E-01

 ETASHRINKSD(%)  8.7966E-01  2.1151E+01  3.3916E+01  3.3731E+01  2.6129E+01
 ETASHRINKVR(%)  1.7516E+00  3.7829E+01  5.6329E+01  5.6085E+01  4.5430E+01
 EBVSHRINKSD(%)  1.1956E+00  2.1809E+01  3.7917E+01  3.5040E+01  2.4296E+01
 EBVSHRINKVR(%)  2.3770E+00  3.8861E+01  6.1457E+01  5.7803E+01  4.2690E+01
 RELATIVEINF(%)  9.7573E+01  1.1555E+01  1.9216E+01  7.9372E+00  2.3830E+01
 EPSSHRINKSD(%)  1.7975E+01
 EPSSHRINKVR(%)  3.2720E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2945.3716152843717     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1296.7958867151890     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.60
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.69
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2945.372       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.22E+00  2.04E+00  8.79E-01  1.39E+00  9.52E-01  1.07E+00  2.72E+00  7.67E-01  1.09E+00  2.24E+00
 


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
+        4.22E+08
 
 TH 2
+       -3.61E+02  5.80E+07
 
 TH 3
+        3.43E+07 -5.59E+03  2.79E+06
 
 TH 4
+        6.48E+08  1.07E+04  5.26E+07  9.94E+08
 
 TH 5
+        2.08E+02 -3.55E+07  3.39E+03 -6.13E+03  2.17E+07
 
 TH 6
+       -1.88E+00 -5.16E+02  3.54E-01 -2.14E+03  3.12E+02  2.16E+02
 
 TH 7
+        9.12E-01  8.38E+02 -8.76E-01  3.34E+03 -4.95E+02 -4.64E-02  7.55E+01
 
 TH 8
+       -1.90E+07  2.64E+04 -1.54E+06 -2.91E+07 -1.61E+04  6.23E+01 -9.89E+01  8.52E+05
 
 TH 9
+        4.49E+08  1.85E+03  3.65E+07  6.89E+08 -1.13E+03  1.64E+00  2.91E+01 -2.02E+07  4.77E+08
 
 TH10
+       -2.75E+08  6.24E+02 -2.24E+07 -4.23E+08 -4.11E+02 -1.95E-01  3.29E+00  1.24E+07 -2.93E+08  1.80E+08
 
 TH11
+       -1.34E+01  2.88E+03 -6.83E-01  1.20E+04 -1.78E+03  2.99E+00  5.28E+00  1.26E+06 -2.98E+07  1.83E+07  2.18E+02
 
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
 #CPUT: Total CPU Time in Seconds,       45.398
Stop Time:
Sat Sep 25 01:03:17 CDT 2021
