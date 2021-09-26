Fri Sep 24 23:02:45 CDT 2021
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
$DATA ../../../../data/int/S1/dat11.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m11.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3328.15237783033        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9298E+02 -6.2014E+01  8.7657E+01 -2.0344E+01  9.7996E+01  2.1194E+01 -1.7233E+01 -4.6445E+02 -1.6598E+02  2.4964E+01
            -3.7868E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3499.13702382434        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.2021E-01  1.0449E+00  9.6025E-01  1.0109E+00  9.4800E-01  9.7533E-01  1.0124E+00  1.4921E+00  1.1028E+00  9.7020E-01
             1.2760E+00
 PARAMETER:  1.6842E-02  1.4392E-01  5.9439E-02  1.1083E-01  4.6601E-02  7.5021E-02  1.1232E-01  5.0018E-01  1.9783E-01  6.9745E-02
             3.4372E-01
 GRADIENT:  -1.1659E+01  1.0049E+01  2.3935E+01  1.6407E+01 -1.3175E+01  2.1303E+01  4.9408E+00 -2.1805E+02 -4.1190E+01  2.5275E+01
             1.9946E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3642.04857210816        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      167
 NPARAMETR:  9.1259E-01  1.2614E+00  1.1373E+00  8.8208E-01  1.1995E+00  1.0696E+00  1.1961E+00  1.7072E+00  1.1634E+00  1.1165E+00
             1.2772E+00
 PARAMETER:  8.5268E-03  3.3222E-01  2.2863E-01 -2.5469E-02  2.8193E-01  1.6725E-01  2.7904E-01  6.3487E-01  2.5132E-01  2.1016E-01
             3.4463E-01
 GRADIENT:  -2.1975E+01  4.1222E+01 -3.5064E+01  2.6753E+01  6.5817E+01  5.2111E+01  3.5751E+01  1.9598E+01 -2.1986E+00  2.8794E+01
             2.0891E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3650.24161471802        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      271
 NPARAMETR:  9.1265E-01  1.2308E+00  1.1772E+00  8.9800E-01  1.1700E+00  1.0700E+00  1.1699E+00  1.7184E+00  1.2082E+00  1.0088E+00
             1.2601E+00
 PARAMETER:  8.6001E-03  3.0765E-01  2.6315E-01 -7.5837E-03  2.5699E-01  1.6762E-01  2.5693E-01  6.4140E-01  2.8916E-01  1.0874E-01
             3.3117E-01
 GRADIENT:  -4.5839E+01  2.1247E+01 -2.6961E+01  2.0914E+01  5.0599E+01  4.6446E+01  2.9918E+01  2.0781E+01  2.2491E+00  1.8214E+01
             1.8835E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3653.44199931067        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  9.1267E-01  1.2308E+00  1.1772E+00  8.8693E-01  1.1700E+00  9.6837E-01  1.1699E+00  1.7184E+00  1.2083E+00  1.0088E+00
             1.2601E+00
 PARAMETER:  8.6157E-03  3.0770E-01  2.6315E-01 -1.9992E-02  2.5699E-01  6.7863E-02  2.5693E-01  6.4140E-01  2.8920E-01  1.0876E-01
             3.3117E-01
 GRADIENT:  -5.5985E+01  1.6003E+01 -2.5808E+01  7.3071E+00  4.6453E+01  1.3151E+01  3.0106E+01  2.0857E+01  1.6836E+00  1.8155E+01
             1.8810E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3653.59429599494        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  9.1267E-01  1.2308E+00  1.1772E+00  8.8691E-01  1.1700E+00  9.6736E-01  1.1645E+00  1.7184E+00  1.2083E+00  1.0088E+00
             1.2601E+00
 PARAMETER:  8.6156E-03  3.0770E-01  2.6315E-01 -2.0012E-02  2.5699E-01  6.6812E-02  2.5228E-01  6.4140E-01  2.8920E-01  1.0876E-01
             3.3117E-01
 GRADIENT:  -3.1066E+01  3.5094E+01 -2.4929E+01  1.1645E+01  5.3933E+01  1.6301E+01  3.1017E+01  2.1324E+01  3.2518E+00  1.8394E+01
             1.8847E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3653.69795615177        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      732             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1267E-01  1.2309E+00  1.1772E+00  8.8691E-01  1.1700E+00  9.6272E-01  1.1629E+00  1.7184E+00  1.2071E+00  1.0087E+00
             1.2601E+00
 PARAMETER:  8.6188E-03  3.0771E-01  2.6316E-01 -2.0015E-02  2.5698E-01  6.2008E-02  2.5091E-01  6.4138E-01  2.8824E-01  1.0865E-01
             3.3116E-01
 GRADIENT:  -3.1609E+01  3.5115E+01 -2.4931E+01  1.1613E+01  5.3911E+01  1.4526E+01  3.0801E+01  2.1316E+01  3.0164E+00  1.8344E+01
             1.8834E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3653.77880431813        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      862             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1262E-01  1.2310E+00  1.1772E+00  8.8695E-01  1.1698E+00  9.1431E-01  1.1630E+00  1.7185E+00  1.2064E+00  1.0087E+00
             1.2598E+00
 PARAMETER:  8.5688E-03  3.0786E-01  2.6313E-01 -1.9965E-02  2.5686E-01  1.0412E-02  2.5103E-01  6.4145E-01  2.8767E-01  1.0871E-01
             3.3099E-01
 GRADIENT:  -3.7977E+01  3.5513E+01 -2.4947E+01  1.1764E+01  5.3764E+01 -5.6342E+00  3.0809E+01  2.1330E+01  2.9002E+00  1.8307E+01
             1.8793E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3653.88603198360        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  9.1262E-01  1.2310E+00  1.1772E+00  8.8695E-01  1.1698E+00  9.3902E-01  1.1630E+00  1.7185E+00  1.2064E+00  1.0087E+00
             1.2598E+00
 PARAMETER:  8.5692E-03  3.0786E-01  2.6313E-01 -1.9966E-02  2.5685E-01  3.7079E-02  2.5103E-01  6.4145E-01  2.8767E-01  1.0871E-01
             3.3096E-01
 GRADIENT:  -5.9708E+01  1.6254E+01 -2.5742E+01  7.4056E+00  4.6164E+01  1.4992E+00  2.9282E+01  2.0864E+01  1.2302E+00  1.8048E+01
             1.8740E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3654.90100807780        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1208
 NPARAMETR:  9.1273E-01  1.2308E+00  1.1774E+00  8.8694E-01  1.1685E+00  9.8415E-01  1.1626E+00  1.7168E+00  1.2064E+00  1.0087E+00
             1.2501E+00
 PARAMETER:  8.6825E-03  3.0769E-01  2.6333E-01 -1.9979E-02  2.5570E-01  8.4027E-02  2.5068E-01  6.4047E-01  2.8765E-01  1.0867E-01
             3.2319E-01
 GRADIENT:  -5.3877E+01  1.7495E+01 -2.5185E+01  6.6897E+00  4.5195E+01  1.8996E+01  2.8929E+01  2.0716E+01  8.1245E-01  1.8013E+01
             1.7407E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3655.70138948728        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1364
 NPARAMETR:  9.1291E-01  1.2307E+00  1.1776E+00  8.8693E-01  1.1677E+00  9.8609E-01  1.1624E+00  1.7158E+00  1.2064E+00  1.0078E+00
             1.2445E+00
 PARAMETER:  8.8842E-03  3.0759E-01  2.6345E-01 -1.9987E-02  2.5504E-01  8.5992E-02  2.5048E-01  6.3990E-01  2.8761E-01  1.0779E-01
             3.1873E-01
 GRADIENT:  -2.7448E+01  3.7931E+01 -2.4067E+01  1.0738E+01  5.2486E+01  2.3337E+01  3.0272E+01  2.1104E+01  2.2646E+00  1.8129E+01
             1.6678E+02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -3655.70530786604        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:     1449
 NPARAMETR:  9.1292E-01  1.2307E+00  1.1776E+00  8.8693E-01  1.1677E+00  9.8602E-01  1.1624E+00  1.7158E+00  1.2064E+00  1.0077E+00
             1.2445E+00
 PARAMETER:  8.8957E-03  3.0759E-01  2.6345E-01 -1.9987E-02  2.5504E-01  8.5922E-02  2.5048E-01  6.3990E-01  2.8761E-01  1.0768E-01
             3.1873E-01
 GRADIENT:  -2.7847E+05  1.1017E+01 -1.0569E+05 -5.5683E+05  4.1550E+01 -5.5681E+05  2.3533E+01 -8.7086E+04  9.6799E+04 -5.1709E+05
             1.4796E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1449
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4155E-02 -4.0807E-02 -7.8514E-03  1.8434E-02 -6.1859E-02
 SE:             2.8236E-02  2.0906E-02  1.6346E-02  2.6226E-02  2.0903E-02
 N:                     100         100         100         100         100

 P VAL.:         3.9230E-01  5.0948E-02  6.3099E-01  4.8213E-01  3.0827E-03

 ETASHRINKSD(%)  5.4043E+00  2.9962E+01  4.5240E+01  1.2140E+01  2.9974E+01
 ETASHRINKVR(%)  1.0517E+01  5.0947E+01  7.0013E+01  2.2806E+01  5.0963E+01
 EBVSHRINKSD(%)  4.1031E-01  2.0481E+01  3.7823E+01  1.3567E+01  2.2038E+01
 EBVSHRINKVR(%)  8.1895E-01  3.6767E+01  6.1341E+01  2.5294E+01  3.9219E+01
 RELATIVEINF(%)  9.9179E+01  3.2973E+01  2.9797E+01  4.5248E+01  2.9086E+01
 EPSSHRINKSD(%)  2.7197E+01
 EPSSHRINKVR(%)  4.6997E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3655.7053078660429     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2001.6159480976321     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.12
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3655.705       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.13E-01  1.23E+00  1.18E+00  8.87E-01  1.17E+00  9.86E-01  1.16E+00  1.72E+00  1.21E+00  1.01E+00  1.24E+00
 


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
+        1.67E+09
 
 TH 2
+       -2.22E+00  1.94E+08
 
 TH 3
+        4.92E+08  2.67E+01  2.89E+08
 
 TH 4
+       -2.53E+03  2.03E+02  5.06E+08  3.54E+09
 
 TH 5
+       -3.29E+00 -1.83E+02 -9.72E+01  1.67E+02  3.14E+08
 
 TH 6
+       -3.81E+03 -8.01E-03 -1.14E+03 -3.98E+03  1.08E+00  1.43E+09
 
 TH 7
+        1.55E-01  1.16E+01  1.54E+08  5.39E+08  6.44E+00 -9.42E-01  3.28E+08
 
 TH 8
+       -2.04E+02 -3.14E+00  4.09E+07  1.43E+08 -2.03E+01 -3.22E+02  4.36E+07  2.31E+07
 
 TH 9
+        6.48E+02 -7.12E+00 -1.29E+08  4.54E+04  3.91E+00  1.02E+03  8.58E+00 -3.65E+07  1.16E+08
 
 TH10
+       -1.60E+03 -1.17E+01  4.14E+08 -1.64E+03 -5.51E+01 -3.25E+03  3.96E+00 -1.32E+02  4.26E+02  1.18E+09
 
 TH11
+       -8.77E+00 -3.90E+01 -6.03E+00  9.70E+00 -2.52E+01 -1.70E-02 -1.21E+08  1.29E+01  1.75E+01  1.14E+01  1.77E+08
 
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
 #CPUT: Total CPU Time in Seconds,       54.913
Stop Time:
Fri Sep 24 23:03:41 CDT 2021
