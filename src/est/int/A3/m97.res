Wed Sep 29 00:40:21 CDT 2021
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
$DATA ../../../../data/int/A3/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1536.96759783356        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5455E+02  1.6189E+02  2.3734E+02 -6.8651E+01  2.1917E+02  2.6327E+01 -1.3908E+02 -2.8509E+02  5.9580E-01 -1.1494E+02
            -3.9079E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2869.55650744658        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  8.9898E-01  1.0277E+00  9.6136E-01  1.1248E+00  9.3378E-01  9.0525E-01  9.6517E-01  9.3496E-01  8.0874E-01  8.1473E-01
             2.3768E+00
 PARAMETER: -6.4903E-03  1.2737E-01  6.0589E-02  2.1764E-01  3.1488E-02  4.5788E-04  6.4547E-02  3.2753E-02 -1.1228E-01 -1.0490E-01
             9.6575E-01
 GRADIENT:   7.6102E+00  8.4454E+01  6.3083E+00  1.1851E+02 -2.6190E+01 -3.1640E+00 -7.3001E+00  1.3624E+01  1.1692E+00  1.1440E+01
             1.8548E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2874.22220264375        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.0708E-01  9.0049E-01  8.1982E-01  1.1889E+00  8.0238E-01  9.1212E-01  1.1191E+00  6.5838E-01  8.1544E-01  7.0327E-01
             2.3861E+00
 PARAMETER:  2.4708E-03 -4.8171E-03 -9.8672E-02  2.7303E-01 -1.2017E-01  8.0159E-03  2.1254E-01 -3.1798E-01 -1.0403E-01 -2.5202E-01
             9.6965E-01
 GRADIENT:   3.0346E+01  7.2897E+01 -5.4335E-01  1.3199E+02 -1.7061E+01 -2.1054E-01  5.4608E+00  8.1965E+00  9.1756E+00  6.4037E+00
             3.7617E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2878.02854983795        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  8.8231E-01  5.5628E-01  4.8558E-01  1.2887E+00  4.8658E-01  8.9949E-01  1.0573E+00  2.8852E-01  7.5468E-01  5.9730E-01
             2.2633E+00
 PARAMETER: -2.5208E-02 -4.8648E-01 -6.2242E-01  3.5366E-01 -6.2035E-01 -5.9225E-03  1.5572E-01 -1.1430E+00 -1.8147E-01 -4.1534E-01
             9.1683E-01
 GRADIENT:  -4.0068E+01  9.1399E+01 -3.5866E+01  1.5781E+02  1.5190E+02 -1.0303E+01 -1.2587E+01 -9.2998E-01 -3.2636E+01 -3.1242E+01
            -1.2453E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2907.12285557002        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  9.2970E-01  3.5777E-01  3.1008E-01  1.2960E+00  3.1487E-01  9.4250E-01  9.4393E-01  1.3024E-01  9.4435E-01  8.1673E-01
             2.1647E+00
 PARAMETER:  2.7112E-02 -9.2786E-01 -1.0709E+00  3.5927E-01 -1.0556E+00  4.0777E-02  4.2300E-02 -1.9384E+00  4.2739E-02 -1.0245E-01
             8.7229E-01
 GRADIENT:   3.5116E+01  2.3228E+01  1.3546E+01  1.2127E+02 -1.2293E+01  4.2996E+00 -1.5392E+01 -6.0576E-01  5.9077E+00  1.0842E+01
            -2.1500E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2933.33077659480        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  9.1409E-01  1.8555E-01  1.2046E-01  9.6102E-01  1.6322E-01  9.3420E-01  1.2608E+00  2.9856E-02  1.1765E+00  8.6435E-01
             2.0170E+00
 PARAMETER:  1.0176E-02 -1.5844E+00 -2.0165E+00  6.0236E-02 -1.7127E+00  3.1936E-02  3.3172E-01 -3.4114E+00  2.6258E-01 -4.5777E-02
             8.0162E-01
 GRADIENT:  -1.3898E+01 -9.9613E+00  7.8920E+00  1.9967E+01 -1.5934E+01 -1.5633E+00 -5.0695E+00 -1.4294E-01 -1.5899E+01  9.4811E+00
             2.5462E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2936.66797180385        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      760
 NPARAMETR:  9.1677E-01  1.6426E-01  9.5008E-02  8.5873E-01  1.4384E-01  9.3786E-01  1.3282E+00  2.1231E-02  1.4735E+00  8.6414E-01
             1.9767E+00
 PARAMETER:  1.3097E-02 -1.7063E+00 -2.2538E+00 -5.2302E-02 -1.8391E+00  3.5844E-02  3.8385E-01 -3.7523E+00  4.8765E-01 -4.6023E-02
             7.8143E-01
 GRADIENT:  -1.1040E+00 -1.9822E-01 -9.0738E-04  4.4092E-01 -8.8968E-02 -2.3319E-02  5.7397E-01 -6.1854E-02  4.0373E-02 -1.0439E-01
             9.6846E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2936.69518111586        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      939
 NPARAMETR:  9.1714E-01  1.6437E-01  9.4809E-02  8.5802E-01  1.4352E-01  9.3789E-01  1.3247E+00  2.5220E-02  1.4743E+00  8.6426E-01
             1.9761E+00
 PARAMETER:  1.3503E-02 -1.7056E+00 -2.2559E+00 -5.3123E-02 -1.8413E+00  3.5873E-02  3.8118E-01 -3.5801E+00  4.8815E-01 -4.5882E-02
             7.8110E-01
 GRADIENT:  -6.6514E-02  1.6624E+00  1.8653E-01  6.1269E-01 -3.3562E+00 -1.1835E-02 -7.8137E-02 -8.7162E-02 -2.0491E-01 -6.2512E-02
             1.8092E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2980.89213521882        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1122
 NPARAMETR:  9.1027E-01  1.5901E-01  9.5034E-02  8.9275E-01  1.4432E-01  9.3446E-01  1.2810E+00  1.1709E+00  1.6351E+00  7.4548E-01
             1.8333E+00
 PARAMETER:  5.9880E-03 -1.7388E+00 -2.2535E+00 -1.3451E-02 -1.8357E+00  3.2216E-02  3.4761E-01  2.5781E-01  5.9172E-01 -1.9373E-01
             7.0610E-01
 GRADIENT:  -2.1047E+00 -3.1394E+01 -1.1196E+01 -1.2040E+01 -3.3799E+00  1.6021E-02  1.1342E+00 -9.2360E-01 -1.1576E+00  8.6646E+00
             2.3569E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2981.80540436979        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1308             RESET HESSIAN, TYPE I
 NPARAMETR:  9.0887E-01  1.6045E-01  9.6989E-02  9.1705E-01  1.4505E-01  9.3470E-01  1.2699E+00  1.1953E+00  1.6581E+00  7.1993E-01
             1.8306E+00
 PARAMETER:  4.4434E-03 -1.7298E+00 -2.2332E+00  1.3404E-02 -1.8307E+00  3.2475E-02  3.3897E-01  2.7836E-01  6.0568E-01 -2.2860E-01
             7.0465E-01
 GRADIENT:   7.2540E+01  9.9927E+01  1.6224E+02  9.9442E+00  1.0325E+03  6.0732E+00  4.5297E+00  3.0153E+00  1.7334E+01  1.0907E+01
             1.4091E+01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -2981.82211992997        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1402
 NPARAMETR:  9.0887E-01  1.6049E-01  9.7215E-02  9.1637E-01  1.4501E-01  9.3456E-01  1.2692E+00  1.1916E+00  1.6581E+00  7.1993E-01
             1.8306E+00
 PARAMETER:  4.4460E-03 -1.7295E+00 -2.2308E+00  1.2663E-02 -1.8309E+00  3.2316E-02  3.3835E-01  2.7531E-01  6.0566E-01 -2.2860E-01
             7.0463E-01
 GRADIENT:   4.4021E+03 -3.3183E+01  1.2092E+01 -2.8221E-01  2.2890E+02 -5.8315E-02 -1.3012E+03 -4.3360E-02 -7.3111E+02  9.6994E+02
            -6.1823E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1402
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8534E-03  3.9149E-02  2.8474E-02  3.0267E-03  3.3256E-02
 SE:             2.9438E-02  2.5436E-02  2.1017E-02  2.6199E-02  2.3069E-02
 N:                     100         100         100         100         100

 P VAL.:         8.9586E-01  1.2377E-01  1.7548E-01  9.0803E-01  1.4943E-01

 ETASHRINKSD(%)  1.3781E+00  1.4786E+01  2.9591E+01  1.2229E+01  2.2716E+01
 ETASHRINKVR(%)  2.7372E+00  2.7386E+01  5.0426E+01  2.2963E+01  4.0271E+01
 EBVSHRINKSD(%)  1.5080E+00  1.4304E+01  2.8709E+01  9.4291E+00  2.2060E+01
 EBVSHRINKVR(%)  2.9932E+00  2.6562E+01  4.9176E+01  1.7969E+01  3.9254E+01
 RELATIVEINF(%)  9.6914E+01  3.1817E+01  1.3160E+01  5.0161E+01  1.4503E+01
 EPSSHRINKSD(%)  2.2222E+01
 EPSSHRINKVR(%)  3.9506E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2981.8221199299674     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1327.7327601615566     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.80
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2981.822       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.09E-01  1.60E-01  9.72E-02  9.16E-01  1.45E-01  9.35E-01  1.27E+00  1.19E+00  1.66E+00  7.20E-01  1.83E+00
 


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
+        2.67E+06
 
 TH 2
+        5.54E+00  1.63E+05
 
 TH 3
+        5.82E+05 -1.96E+05  4.09E+04
 
 TH 4
+       -4.78E+03 -2.96E+02 -7.44E+02  3.97E+02
 
 TH 5
+       -4.63E+05  1.40E+05  1.55E+05 -9.50E+02  3.75E+05
 
 TH 6
+        6.89E+02 -1.70E+01  2.88E+01 -3.90E+00  2.13E+02  2.15E+02
 
 TH 7
+        4.04E+02  1.96E+02 -9.55E+03  1.01E+03 -6.98E+01 -1.44E+02  5.97E+04
 
 TH 8
+       -6.79E+03  2.96E+01 -1.60E+05 -4.81E+00 -2.39E+03 -3.57E-01  1.44E+03  4.06E+01
 
 TH 9
+       -1.92E+03  6.99E+01 -3.94E+03  4.28E+02  2.81E+02 -5.93E+01 -3.28E+01  6.22E+02  1.11E+04
 
 TH10
+        1.95E+03  4.41E+01  2.52E+04 -2.64E+03  3.71E+02  3.74E+02  2.22E+02 -3.73E+03 -6.73E+04  4.08E+05
 
 TH11
+        1.73E+03 -1.88E+01 -3.18E+03  3.50E+02  1.39E+02 -4.54E+01 -1.81E+01  4.86E+02  8.50E+03 -5.15E+04  6.79E+03
 
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
 #CPUT: Total CPU Time in Seconds,       52.161
Stop Time:
Wed Sep 29 00:41:15 CDT 2021
