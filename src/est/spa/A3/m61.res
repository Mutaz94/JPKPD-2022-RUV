Wed Sep 29 13:40:29 CDT 2021
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
$DATA ../../../../data/spa/A3/dat61.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   308.248716908946        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1939E+02  7.2110E+01  1.2700E+02 -7.2469E+01  1.9494E+02  5.9176E+01 -4.0686E+01 -8.9310E+01 -2.7118E+02 -1.1490E+02
            -3.4036E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1116.87443100887        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1716E+00  9.5453E-01  8.2697E-01  1.2889E+00  8.0495E-01  7.6631E-01  9.6485E-01  1.0412E+00  1.3987E+00  9.7151E-01
             7.3334E+00
 PARAMETER:  2.5838E-01  5.3460E-02 -8.9988E-02  3.5377E-01 -1.1697E-01 -1.6617E-01  6.4215E-02  1.4036E-01  4.3554E-01  7.1094E-02
             2.0924E+00
 GRADIENT:   1.5496E+02 -1.2609E+01 -1.1710E+01  1.4084E+01 -5.5175E+00 -1.1881E+01  9.1038E+00  8.8580E+00  4.0162E+01  2.0159E+01
             2.8946E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1145.84602497576        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.1059E+00  4.2237E-01  1.9079E-01  1.5400E+00  2.1718E-01  1.0060E+00  4.0599E-01  7.4089E-01  2.3611E+00  2.2238E-01
             5.5989E+00
 PARAMETER:  2.0069E-01 -7.6188E-01 -1.5566E+00  5.3176E-01 -1.4270E+00  1.0602E-01 -8.0141E-01 -1.9990E-01  9.5913E-01 -1.4034E+00
             1.8226E+00
 GRADIENT:   4.5304E+01  1.1328E+02  1.3655E+02  1.1092E+02 -2.5334E+02  2.9414E+01 -3.1266E-01  2.9294E+00  4.9113E+01 -3.0369E-01
             1.9352E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1170.98384441892        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.0041E+00  3.2475E-01  8.6031E-02  1.5352E+00  1.5961E-01  1.0913E+00  1.8910E-01  2.5959E-01  3.1012E+00  1.1123E-01
             4.3819E+00
 PARAMETER:  1.0413E-01 -1.0247E+00 -2.3530E+00  5.2868E-01 -1.7350E+00  1.8733E-01 -1.5655E+00 -1.2487E+00  1.2318E+00 -2.0962E+00
             1.5775E+00
 GRADIENT:   3.3629E+01  8.6549E+01  7.7350E+01  7.5950E+01 -1.5267E+02  4.6829E+01 -4.1141E-01 -7.9231E-01 -3.0740E+01 -2.1859E+00
             4.3376E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1191.24131710069        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      416
 NPARAMETR:  9.5349E-01  3.2551E-01  8.5133E-02  1.4085E+00  1.6205E-01  9.3359E-01  6.6236E-02  4.7912E-02  3.7572E+00  3.4857E-01
             3.9450E+00
 PARAMETER:  5.2378E-02 -1.0223E+00 -2.3635E+00  4.4252E-01 -1.7198E+00  3.1285E-02 -2.6145E+00 -2.9384E+00  1.4237E+00 -9.5390E-01
             1.4724E+00
 GRADIENT:  -3.3998E+01  2.3252E+01  1.7165E+01  3.2669E+01 -3.4378E+01  1.2545E+01  3.2580E-02 -1.4811E-02  1.3556E+01  9.8667E-01
             2.0662E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1192.95388950814        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      559
 NPARAMETR:  9.5461E-01  3.2205E-01  8.2200E-02  1.3770E+00  1.6072E-01  9.1748E-01  1.0000E-02  6.2257E-02  3.6718E+00  3.4706E-01
             3.8608E+00
 PARAMETER:  5.3544E-02 -1.0331E+00 -2.3986E+00  4.1988E-01 -1.7281E+00  1.3873E-02 -4.8050E+00 -2.6765E+00  1.4007E+00 -9.5826E-01
             1.4509E+00
 GRADIENT:  -1.7606E+01  2.2848E+01  2.4098E+01  4.1459E+01  1.1094E+02  6.2181E+00  0.0000E+00 -3.1142E-02  3.6600E+01 -5.4931E-01
             1.2756E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1194.57843713943        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  9.5489E-01  3.1848E-01  8.1975E-02  1.3701E+00  1.6135E-01  9.1751E-01  1.0000E-02  6.8908E-01  3.6021E+00  3.4705E-01
             3.8336E+00
 PARAMETER:  5.3836E-02 -1.0442E+00 -2.4013E+00  4.1485E-01 -1.7241E+00  1.3908E-02 -5.3408E+00 -2.7240E-01  1.3815E+00 -9.5827E-01
             1.4438E+00
 GRADIENT:  -4.5758E+01 -2.1006E+00 -6.5670E+00  3.5573E+01  9.8699E+00  7.2746E+00  0.0000E+00 -2.3860E-01  1.2409E+01  2.1063E+00
             1.7263E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1195.19579256886        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  9.5653E-01  3.1336E-01  8.4770E-02  1.3584E+00  1.6358E-01  9.1895E-01  1.0000E-02  7.0691E-01  3.4991E+00  3.5200E-01
             3.7180E+00
 PARAMETER:  5.5557E-02 -1.0604E+00 -2.3678E+00  4.0633E-01 -1.7104E+00  1.5480E-02 -5.4512E+00 -2.4685E-01  1.3525E+00 -9.4413E-01
             1.4132E+00
 GRADIENT:  -4.0341E+01 -2.0294E+01 -9.9047E+00  3.4388E+01  4.7966E+01  1.0383E+01  0.0000E+00 -1.2929E+00  9.7106E+00  2.9676E+00
             6.8424E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1240.45802720062        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1018
 NPARAMETR:  1.0039E+00  2.9239E-01  7.4906E-02  9.2798E-01  1.5357E-01  9.8061E-01  1.0000E-02  1.1302E+00  2.0709E+00  5.4374E-01
             3.6122E+00
 PARAMETER:  1.0390E-01 -1.1297E+00 -2.4915E+00  2.5256E-02 -1.7736E+00  8.0418E-02 -5.4512E+00  2.2241E-01  8.2797E-01 -5.0928E-01
             1.3843E+00
 GRADIENT:  -5.4221E-01 -9.6337E-01 -2.1500E+01  4.3042E+01  1.1322E+02  1.4260E+01  0.0000E+00  1.4030E+01 -1.1297E+01 -3.2273E+00
             4.5756E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1253.52266885588        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  1.0126E+00  2.9443E-01  9.0797E-02  9.0313E-01  1.6111E-01  9.0956E-01  1.0000E-02  7.1182E-01  2.4322E+00  5.2740E-01
             3.2701E+00
 PARAMETER:  1.1254E-01 -1.1227E+00 -2.2991E+00 -1.8922E-03 -1.7257E+00  5.2092E-03 -5.4512E+00 -2.3993E-01  9.8879E-01 -5.3980E-01
             1.2848E+00
 GRADIENT:   1.0280E+01  4.6307E+00 -1.0282E+01  4.5846E+00 -5.8067E+01 -1.0701E+01  0.0000E+00  4.2664E+00  3.1111E+01 -9.9543E+00
            -1.2791E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1271.29609460034        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1333
 NPARAMETR:  1.0254E+00  3.3426E-01  1.5717E-01  1.0792E+00  2.1530E-01  9.2175E-01  1.0000E-02  1.2194E-01  1.6226E+00  4.6290E-01
             3.4116E+00
 PARAMETER:  1.2511E-01 -9.9584E-01 -1.7504E+00  1.7622E-01 -1.4357E+00  1.8519E-02 -5.4512E+00 -2.0042E+00  5.8406E-01 -6.7025E-01
             1.3272E+00
 GRADIENT:  -6.3971E+00 -1.0463E+01  1.4568E+00  8.8925E+00 -4.5382E+00 -4.3050E+00  0.0000E+00  1.2161E-01  3.5891E+00 -5.1872E-01
             6.4061E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1273.81571547626        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1509
 NPARAMETR:  1.0166E+00  4.6245E-01  1.5079E-01  1.0336E+00  2.4157E-01  9.3320E-01  1.0000E-02  7.7400E-02  1.6191E+00  4.1799E-01
             3.4201E+00
 PARAMETER:  1.1649E-01 -6.7122E-01 -1.7918E+00  1.3307E-01 -1.3206E+00  3.0864E-02 -5.4512E+00 -2.4588E+00  5.8186E-01 -7.7229E-01
             1.3297E+00
 GRADIENT:  -3.1120E-02  6.8613E-01  1.0741E+00 -2.9179E-01 -3.1374E+00 -1.6245E+00  0.0000E+00  2.8093E-02  6.1776E-01 -1.6367E+00
             1.9974E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1273.89342662497        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1685
 NPARAMETR:  1.0159E+00  4.6267E-01  1.4916E-01  1.0311E+00  2.4039E-01  9.3728E-01  1.0000E-02  4.3295E-02  1.6231E+00  4.4708E-01
             3.3851E+00
 PARAMETER:  1.1577E-01 -6.7074E-01 -1.8027E+00  1.3061E-01 -1.3255E+00  3.5224E-02 -5.4512E+00 -3.0397E+00  5.8431E-01 -7.0502E-01
             1.3194E+00
 GRADIENT:   8.0034E-02 -5.6626E-01 -5.5176E-01 -5.5613E-02  1.3162E+00  3.7215E-03  0.0000E+00  6.9850E-03  4.6791E-02  2.1794E-02
             1.1756E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1273.89676892884        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1847
 NPARAMETR:  1.0160E+00  4.6243E-01  1.4935E-01  1.0315E+00  2.4037E-01  9.3724E-01  1.0000E-02  1.0000E-02  1.6221E+00  4.4742E-01
             3.3847E+00
 PARAMETER:  1.1585E-01 -6.7126E-01 -1.8015E+00  1.3106E-01 -1.3256E+00  3.5179E-02 -5.4512E+00 -4.8172E+00  5.8372E-01 -7.0426E-01
             1.3193E+00
 GRADIENT:  -1.0738E-02  2.4085E-04  2.0776E-03  3.2092E-02  1.1749E-01 -1.1734E-02  0.0000E+00  0.0000E+00 -1.5363E-02 -4.9104E-03
            -1.4211E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1847
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0867E-03 -1.1060E-04  1.8740E-04 -1.4916E-02  9.3541E-03
 SE:             2.8499E-02  1.9629E-04  1.6574E-04  2.6302E-02  1.7496E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6958E-01  5.7313E-01  2.5819E-01  5.7064E-01  5.9289E-01

 ETASHRINKSD(%)  4.5236E+00  9.9342E+01  9.9445E+01  1.1886E+01  4.1387E+01
 ETASHRINKVR(%)  8.8426E+00  9.9996E+01  9.9997E+01  2.2359E+01  6.5645E+01
 EBVSHRINKSD(%)  4.3310E+00  9.9349E+01  9.9426E+01  9.1320E+00  4.1858E+01
 EBVSHRINKVR(%)  8.4744E+00  9.9996E+01  9.9997E+01  1.7430E+01  6.6195E+01
 RELATIVEINF(%)  8.4394E+01  1.7883E-04  4.8637E-04  5.8058E+01  9.3290E-01
 EPSSHRINKSD(%)  2.7616E+01
 EPSSHRINKVR(%)  4.7605E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1273.8967689288368     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -538.74594236509859     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1273.897       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  4.62E-01  1.49E-01  1.03E+00  2.40E-01  9.37E-01  1.00E-02  1.00E-02  1.62E+00  4.47E-01  3.38E+00
 


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
+       -6.28E+01  2.18E+03
 
 TH 3
+       -6.38E+02  4.74E+03  1.58E+04
 
 TH 4
+       -1.81E+01  9.71E+01 -5.08E+02  3.26E+02
 
 TH 5
+        5.84E+02 -8.08E+03 -1.96E+04 -1.27E+02  3.21E+04
 
 TH 6
+       -1.17E+00 -1.91E+00  5.27E+01 -8.71E+00  1.49E+01  1.89E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.38E+01 -3.02E+01  1.10E+02 -9.56E+00  1.12E+02  2.40E+00  0.00E+00  0.00E+00  4.37E+01
 
 TH10
+       -1.04E+01 -1.27E+02 -2.33E+02 -1.95E+00  5.97E+02  7.96E+00  0.00E+00  0.00E+00  3.18E+00  1.14E+02
 
 TH11
+       -1.99E+01 -1.27E+01 -1.54E+01 -2.38E+00  2.69E+01  2.07E+00  0.00E+00  0.00E+00  4.75E+00  2.51E+01  2.97E+01
 
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
 #CPUT: Total CPU Time in Seconds,       33.935
Stop Time:
Wed Sep 29 13:41:05 CDT 2021
