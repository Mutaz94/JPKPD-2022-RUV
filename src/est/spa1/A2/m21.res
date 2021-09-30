Wed Sep 29 23:09:06 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat21.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -925.012869097185        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5636E+02  6.0187E+01  1.0856E+02  3.3074E+01  1.4073E+02  4.9451E+01 -8.1821E+01 -8.2783E+01 -5.2957E+01 -7.8785E+01
            -2.1079E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1686.23621099473        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8049E-01  1.1039E+00  1.1757E+00  1.0178E+00  1.1467E+00  8.2039E-01  8.5776E-01  9.2315E-01  8.4394E-01  6.6998E-01
             2.7687E+00
 PARAMETER:  8.0293E-02  1.9887E-01  2.6187E-01  1.1761E-01  2.3690E-01 -9.7975E-02 -5.3432E-02  2.0042E-02 -6.9675E-02 -3.0050E-01
             1.1184E+00
 GRADIENT:  -1.8845E+02  1.3825E+01 -1.1103E+01  8.2813E+00  7.9021E+00 -6.3613E+01 -1.2284E+01  5.0361E+00 -3.2347E+01  1.0010E+01
            -2.0502E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1696.86532621172        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.7987E-01  1.0128E+00  1.0248E+00  1.1015E+00  9.5129E-01  8.9716E-01  9.2978E-01  6.8694E-01  1.1324E+00  3.0779E-01
             2.8234E+00
 PARAMETER:  7.9663E-02  1.1269E-01  1.2446E-01  1.9671E-01  5.0066E-02 -8.5175E-03  2.7190E-02 -2.7551E-01  2.2434E-01 -1.0783E+00
             1.1379E+00
 GRADIENT:  -1.5694E+02  5.6256E+01  2.5686E+01  6.2580E+01 -5.6151E+01 -2.7206E+01  3.9910E+00  3.6127E+00  1.8623E+01  2.0294E+00
             2.0690E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1710.32461625638        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0283E+00  9.1919E-01  6.6073E-01  1.0800E+00  7.5014E-01  9.3832E-01  1.1342E+00  1.5211E-01  9.8477E-01  1.9288E-01
             2.7535E+00
 PARAMETER:  1.2787E-01  1.5739E-02 -3.1440E-01  1.7695E-01 -1.8749E-01  3.6335E-02  2.2591E-01 -1.7831E+00  8.4650E-02 -1.5457E+00
             1.1129E+00
 GRADIENT:  -7.4944E+00  6.0498E-01 -6.1340E+00  1.2728E+01  2.0459E+01  1.2571E+00  3.4152E+00  2.7157E-01  4.8071E+00  1.1963E+00
             1.1158E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1712.40314122501        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0477E+00  7.7849E-01  5.3678E-01  1.1807E+00  6.0543E-01  9.8277E-01  1.4069E+00  5.8301E-02  9.6958E-01  1.0894E-01
             2.7085E+00
 PARAMETER:  1.4658E-01 -1.5040E-01 -5.2217E-01  2.6613E-01 -4.0181E-01  8.2621E-02  4.4141E-01 -2.7421E+00  6.9104E-02 -2.1170E+00
             1.0964E+00
 GRADIENT:  -3.2551E+01  2.2091E+01 -2.0870E+01  6.5540E+01  3.5667E+01  1.2549E+01  1.1071E+01  3.6590E-02  8.3959E+00  2.7574E-01
             2.7707E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1725.83036847331        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0469E+00  3.6258E-01  3.3962E-01  1.2428E+00  3.5152E-01  9.4429E-01  2.4742E+00  1.0000E-02  9.1727E-01  1.4841E-02
             2.7197E+00
 PARAMETER:  1.4586E-01 -9.1452E-01 -9.7991E-01  3.1741E-01 -9.4548E-01  4.2677E-02  1.0059E+00 -8.2808E+00  1.3642E-02 -4.1104E+00
             1.1005E+00
 GRADIENT:  -2.5564E+01  1.6014E+01  1.2104E+01  4.1068E+01 -1.3022E+01 -6.4382E-01  2.1322E+01  0.0000E+00 -1.4436E+01 -3.3990E-02
             6.6939E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1736.29555392724        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.0330E+00  4.0911E-01  2.6734E-01  1.1564E+00  3.1515E-01  1.0091E+00  1.9457E+00  1.0000E-02  1.1287E+00  1.7305E-02
             2.4122E+00
 PARAMETER:  1.3247E-01 -7.9377E-01 -1.2193E+00  2.4528E-01 -1.0547E+00  1.0908E-01  7.6563E-01 -9.7693E+00  2.2110E-01 -3.9568E+00
             9.8055E-01
 GRADIENT:  -4.6535E+01  2.3093E+00  8.3117E+00  2.3664E+01 -8.1456E+00  1.6119E+01  1.3160E+01  0.0000E+00  1.6547E+01 -6.9352E-02
            -1.8919E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1739.29948036638        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  1.0571E+00  4.5512E-01  2.5650E-01  1.1168E+00  3.1677E-01  9.6531E-01  1.6315E+00  1.0000E-02  1.0852E+00  1.6313E-02
             2.4117E+00
 PARAMETER:  1.5552E-01 -6.8719E-01 -1.2606E+00  2.1046E-01 -1.0496E+00  6.4691E-02  5.8950E-01 -1.0956E+01  1.8176E-01 -4.0158E+00
             9.8034E-01
 GRADIENT:  -1.1696E+00  3.8717E-01  8.0794E-01  1.3770E+00 -1.2276E+00  6.5631E-01 -1.1329E-01  0.0000E+00  3.0581E-01 -5.4333E-02
            -4.9242E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1745.03462441917        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  1.0418E+00  4.4311E-01  2.5538E-01  1.1015E+00  3.1004E-01  9.5611E-01  1.5561E+00  1.0000E-02  1.0846E+00  3.7297E-01
             2.3520E+00
 PARAMETER:  1.4099E-01 -7.1394E-01 -1.2650E+00  1.9669E-01 -1.0711E+00  5.5113E-02  5.4216E-01 -1.0984E+01  1.8118E-01 -8.8627E-01
             9.5527E-01
 GRADIENT:   5.5439E+01  1.9461E+01  1.8067E+01  5.6009E+00  1.0243E+02  1.3437E+00  7.6517E+00  0.0000E+00 -6.0093E-01  1.9278E-01
             2.8530E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1746.76581195982        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1153
 NPARAMETR:  1.0529E+00  4.1575E-01  2.5511E-01  1.1284E+00  3.0186E-01  9.6719E-01  1.4167E+00  1.0000E-02  1.1233E+00  4.5764E-01
             2.2647E+00
 PARAMETER:  1.5155E-01 -7.7768E-01 -1.2660E+00  2.2079E-01 -1.0978E+00  6.6641E-02  4.4833E-01 -1.0984E+01  2.1626E-01 -6.8168E-01
             9.1742E-01
 GRADIENT:  -8.6192E+00  5.7396E+00  2.1755E+00  3.0780E-01 -1.9061E+00  3.4643E-01 -1.0453E-01  0.0000E+00  4.1908E-01 -3.3021E-01
            -9.3833E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1747.03047371168        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1328
 NPARAMETR:  1.0566E+00  3.8655E-01  2.5373E-01  1.1339E+00  2.9532E-01  9.6437E-01  1.4537E+00  1.0000E-02  1.1195E+00  4.7425E-01
             2.2861E+00
 PARAMETER:  1.5505E-01 -8.5048E-01 -1.2715E+00  2.2564E-01 -1.1197E+00  6.3719E-02  4.7410E-01 -1.0984E+01  2.1285E-01 -6.4602E-01
             9.2685E-01
 GRADIENT:  -3.0144E-01 -2.1563E-01  4.9031E-01 -1.4271E-01  1.1720E+00 -1.0167E-01  3.2442E-02  0.0000E+00  1.9985E-02  1.4868E-02
            -4.2659E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1747.03047371168        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1350
 NPARAMETR:  1.0566E+00  3.8655E-01  2.5373E-01  1.1339E+00  2.9532E-01  9.6437E-01  1.4537E+00  1.0000E-02  1.1195E+00  4.7425E-01
             2.2861E+00
 PARAMETER:  1.5505E-01 -8.5048E-01 -1.2715E+00  2.2564E-01 -1.1197E+00  6.3719E-02  4.7410E-01 -1.0984E+01  2.1285E-01 -6.4602E-01
             9.2685E-01
 GRADIENT:  -3.0144E-01 -2.1563E-01  4.9031E-01 -1.4271E-01  1.1720E+00 -1.0167E-01  3.2442E-02  0.0000E+00  1.9985E-02  1.4868E-02
            -4.2659E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1350
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1944E-03  2.5781E-02 -2.1982E-04 -1.0941E-02  1.1664E-02
 SE:             2.9364E-02  1.7786E-02  2.4610E-04  2.7730E-02  1.7889E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6756E-01  1.4719E-01  3.7174E-01  6.9319E-01  5.1438E-01

 ETASHRINKSD(%)  1.6276E+00  4.0415E+01  9.9176E+01  7.0993E+00  4.0069E+01
 ETASHRINKVR(%)  3.2287E+00  6.4497E+01  9.9993E+01  1.3695E+01  6.4083E+01
 EBVSHRINKSD(%)  1.7458E+00  4.2796E+01  9.9175E+01  6.2247E+00  3.8617E+01
 EBVSHRINKVR(%)  3.4612E+00  6.7277E+01  9.9993E+01  1.2062E+01  6.2322E+01
 RELATIVEINF(%)  9.6374E+01  7.0512E+00  4.6364E-04  5.4264E+01  1.6325E+00
 EPSSHRINKSD(%)  2.8722E+01
 EPSSHRINKVR(%)  4.9195E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1747.0304737116783     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -828.09194050700557     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1747.030       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  3.87E-01  2.54E-01  1.13E+00  2.95E-01  9.64E-01  1.45E+00  1.00E-02  1.12E+00  4.74E-01  2.29E+00
 


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
+        1.03E+03
 
 TH 2
+       -4.03E+01  1.15E+03
 
 TH 3
+       -1.08E+01  2.43E+03  1.16E+04
 
 TH 4
+       -1.46E+01  1.42E+02 -6.23E+02  6.09E+02
 
 TH 5
+        8.96E+01 -4.12E+03 -1.46E+04 -1.94E+01  2.18E+04
 
 TH 6
+        2.90E+00 -8.66E+00  2.88E+01 -6.60E+00 -1.58E+00  1.98E+02
 
 TH 7
+        1.41E+00  4.86E+01 -1.14E+01 -1.15E+00 -6.56E+01  5.79E-01  1.43E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.12E+00 -1.30E+01  8.26E+01 -5.45E+00  4.65E+01 -1.70E+00  7.01E+00  0.00E+00  1.19E+02
 
 TH10
+       -3.78E+00  4.81E+01 -2.59E+02 -3.65E+00  3.64E+02  2.26E+00  2.46E+01  0.00E+00 -4.19E+00  1.34E+02
 
 TH11
+       -1.48E+01 -3.00E-02 -1.07E+02 -7.70E+00  8.12E+01  3.06E+00  3.32E+00  0.00E+00  5.35E+00  2.30E+01  8.51E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.018
Stop Time:
Wed Sep 29 23:09:40 CDT 2021
