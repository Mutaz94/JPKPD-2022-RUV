Thu Sep 30 02:32:11 CDT 2021
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
$DATA ../../../../data/spa1/D/dat1.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2082.17408139387        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4456E+02  2.3637E+01  2.1966E+01  5.8910E+01 -4.4013E+01  4.1703E+01 -1.0163E+01 -2.9886E+00  2.0579E+01 -1.5447E+01
            -4.0922E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2087.41876899492        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0170E+00  1.2113E+00  1.0458E+00  9.3651E-01  1.1250E+00  1.0410E+00  1.1621E+00  1.0059E+00  8.7774E-01  1.1856E+00
             1.0512E+00
 PARAMETER:  1.1690E-01  2.9166E-01  1.4478E-01  3.4403E-02  2.1777E-01  1.4018E-01  2.5020E-01  1.0592E-01 -3.0400E-02  2.7028E-01
             1.4996E-01
 GRADIENT:   3.7248E+01  6.7488E+01  2.2165E+01  5.1795E+01 -4.4366E+01  8.0809E+00  9.8444E-03 -4.7296E+00 -6.4887E+00 -4.1940E+00
             2.1938E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2089.01755063452        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  1.0123E+00  1.2893E+00  1.1878E+00  8.7614E-01  1.3144E+00  1.0116E+00  1.1090E+00  1.1743E+00  9.1457E-01  1.4375E+00
             1.0076E+00
 PARAMETER:  1.1224E-01  3.5408E-01  2.7214E-01 -3.2234E-02  3.7334E-01  1.1149E-01  2.0342E-01  2.6067E-01  1.0702E-02  4.6290E-01
             1.0757E-01
 GRADIENT:   3.1788E+01  4.4676E+01  1.2278E+01  3.3788E+01 -2.7770E+00 -2.5594E+00  2.3472E+00 -8.6471E+00 -5.3219E+00  5.4996E+00
            -3.1064E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2091.65788676630        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.9870E-01  1.2788E+00  1.1408E+00  8.5549E-01  1.2854E+00  1.0208E+00  1.0506E+00  1.3823E+00  9.5984E-01  1.3267E+00
             1.0449E+00
 PARAMETER:  9.8698E-02  3.4594E-01  2.3173E-01 -5.6081E-02  3.5105E-01  1.2058E-01  1.4935E-01  4.2377E-01  5.9016E-02  3.8270E-01
             1.4389E-01
 GRADIENT:   3.8853E-01  2.9421E+00 -3.2485E-01  6.0024E+00 -1.6237E+00  1.2888E+00 -3.1429E-01  4.1718E-01  2.1726E-03  3.5495E-01
             2.8116E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2092.04946397623        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  1.0012E+00  1.6085E+00  7.6396E-01  6.4342E-01  1.2980E+00  1.0192E+00  9.0316E-01  1.0961E+00  1.1450E+00  1.3064E+00
             1.0371E+00
 PARAMETER:  1.0117E-01  5.7530E-01 -1.6924E-01 -3.4096E-01  3.6081E-01  1.1900E-01 -1.8563E-03  1.9174E-01  2.3542E-01  3.6728E-01
             1.3644E-01
 GRADIENT:   8.2591E-02  2.1034E+01  4.3629E+00  8.9303E+00 -7.7001E+00 -5.8110E-01 -9.1888E-01 -1.1588E-01 -1.3410E+00 -6.4286E-02
            -2.8151E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2092.60689413072        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  1.0038E+00  1.8366E+00  4.4884E-01  4.8789E-01  1.2714E+00  1.0214E+00  8.3730E-01  6.3055E-01  1.3422E+00  1.2407E+00
             1.0369E+00
 PARAMETER:  1.0384E-01  7.0793E-01 -7.0109E-01 -6.1767E-01  3.4012E-01  1.2116E-01 -7.7576E-02 -3.6116E-01  3.9432E-01  3.1565E-01
             1.3624E-01
 GRADIENT:   1.3968E+00  1.7411E+01 -6.3218E-01  1.0314E+01 -3.3545E+00 -7.9266E-01  1.5664E+00  9.1982E-01  1.7602E+00  1.8449E-01
             2.2932E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2092.79476114410        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1073             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0042E+00  1.8854E+00  3.8073E-01  4.4062E-01  1.2794E+00  1.0244E+00  8.1514E-01  5.2081E-01  1.3922E+00  1.2384E+00
             1.0294E+00
 PARAMETER:  1.0422E-01  7.3416E-01 -8.6566E-01 -7.1957E-01  3.4641E-01  1.2408E-01 -1.0440E-01 -5.5237E-01  4.3090E-01  3.1384E-01
             1.2894E-01
 GRADIENT:   4.4152E+02  9.0363E+02  2.9702E+00  1.2144E+02  2.0434E+01  6.0265E+01  7.9232E+00  3.1382E-01  1.1850E+01  3.0012E+00
             1.7284E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2092.80360952270        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1248
 NPARAMETR:  1.0011E+00  1.8900E+00  3.8577E-01  4.3769E-01  1.2864E+00  1.0221E+00  8.1414E-01  5.0178E-01  1.3934E+00  1.2471E+00
             1.0299E+00
 PARAMETER:  1.0109E-01  7.3658E-01 -8.5251E-01 -7.2624E-01  3.5182E-01  1.2190E-01 -1.0562E-01 -5.8960E-01  4.3174E-01  3.2081E-01
             1.2950E-01
 GRADIENT:  -4.8962E+00 -9.5753E+00 -1.5758E-01 -9.3406E-01 -8.8608E-01 -6.5337E-01 -1.0452E-01  3.1081E-01 -3.2211E-01 -2.7022E-01
            -3.7046E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2092.83229925843        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  9.9628E-01  1.8970E+00  4.0774E-01  4.3076E-01  1.3283E+00  1.0184E+00  8.1198E-01  3.7307E-01  1.4301E+00  1.2956E+00
             1.0355E+00
 PARAMETER:  9.6269E-02  7.4025E-01 -7.9712E-01 -7.4222E-01  3.8387E-01  1.1827E-01 -1.0828E-01 -8.8599E-01  4.5775E-01  3.5896E-01
             1.3488E-01
 GRADIENT:  -1.4559E+01 -1.6563E+01 -7.0127E-01 -2.6982E+00  2.2040E+00 -1.9558E+00 -1.7489E-01  2.9875E-01  2.8446E-01  4.6515E-01
             5.1026E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2093.06759679136        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  1.0027E+00  1.8847E+00  4.1684E-01  4.4459E-01  1.3141E+00  1.0229E+00  8.1868E-01  1.9910E-01  1.4024E+00  1.2872E+00
             1.0363E+00
 PARAMETER:  1.0272E-01  7.3379E-01 -7.7506E-01 -7.1060E-01  3.7317E-01  1.2268E-01 -1.0006E-01 -1.5139E+00  4.3819E-01  3.5247E-01
             1.3561E-01
 GRADIENT:  -7.8349E-01 -2.9619E+00  3.7236E-01  1.9796E-01 -7.1311E-02 -7.2296E-02 -6.9714E-02  8.2869E-02 -1.6638E-01  1.6743E-01
             4.3791E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2093.07753180737        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1781
 NPARAMETR:  1.0014E+00  1.8907E+00  4.0625E-01  4.3965E-01  1.3118E+00  1.0221E+00  8.1631E-01  1.0516E-01  1.4072E+00  1.2843E+00
             1.0356E+00
 PARAMETER:  1.0138E-01  7.3693E-01 -8.0079E-01 -7.2178E-01  3.7141E-01  1.2183E-01 -1.0296E-01 -2.1522E+00  4.4161E-01  3.5024E-01
             1.3502E-01
 GRADIENT:  -3.8017E+00 -5.1439E+00 -1.7601E-01 -8.0579E-03  1.8742E-01 -4.6383E-01 -1.9012E-01  2.4718E-02 -7.9958E-02  2.5790E-01
             3.5225E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2093.09398148339        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1960
 NPARAMETR:  1.0060E+00  1.8883E+00  4.0447E-01  4.3838E-01  1.3113E+00  1.0250E+00  8.1650E-01  2.9463E-02  1.4108E+00  1.2820E+00
             1.0351E+00
 PARAMETER:  1.0596E-01  7.3568E-01 -8.0519E-01 -7.2468E-01  3.7099E-01  1.2469E-01 -1.0272E-01 -3.4246E+00  4.4415E-01  3.4843E-01
             1.3452E-01
 GRADIENT:   6.0136E+00 -1.0698E+01 -3.5888E-01 -1.5824E+00  7.3655E-01  6.7453E-01 -2.4359E-02  2.0143E-03  1.9810E-01  1.4922E-01
             5.5873E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2093.11438307247        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2141
 NPARAMETR:  1.0053E+00  1.8845E+00  4.0820E-01  4.4214E-01  1.3074E+00  1.0243E+00  8.1778E-01  1.0000E-02  1.4034E+00  1.2789E+00
             1.0349E+00
 PARAMETER:  1.0527E-01  7.3365E-01 -7.9601E-01 -7.1613E-01  3.6808E-01  1.2396E-01 -1.0117E-01 -5.4365E+00  4.3891E-01  3.4599E-01
             1.3432E-01
 GRADIENT:   4.5548E+00 -7.4069E+00  1.6419E-01 -1.2893E+00 -3.6584E-01  3.9269E-01 -1.3415E-01  0.0000E+00 -1.4155E-02 -1.3549E-02
            -2.5545E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -2093.11455155042        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     2204
 NPARAMETR:  1.0036E+00  1.8768E+00  4.0746E-01  4.4324E-01  1.3079E+00  1.0240E+00  8.1872E-01  1.0000E-02  1.4042E+00  1.2798E+00
             1.0355E+00
 PARAMETER:  1.0422E-01  7.3414E-01 -7.9650E-01 -7.1562E-01  3.6820E-01  1.2340E-01 -1.0102E-01 -5.4365E+00  4.3883E-01  3.4599E-01
             1.3443E-01
 GRADIENT:   5.2253E-01  5.0324E+00  7.4541E-02 -4.3318E-01 -8.7206E-02 -5.2625E-02 -1.0148E-01  0.0000E+00 -3.5450E-02 -5.3766E-02
            -1.6429E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2204
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7668E-04 -3.0087E-02 -2.7177E-04  3.5366E-02 -4.0524E-02
 SE:             2.9893E-02  2.5762E-02  8.9552E-05  2.1294E-02  2.2554E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9262E-01  2.4285E-01  2.4076E-03  9.6750E-02  7.2378E-02

 ETASHRINKSD(%)  1.0000E-10  1.3695E+01  9.9700E+01  2.8662E+01  2.4440E+01
 ETASHRINKVR(%)  1.0000E-10  2.5514E+01  9.9999E+01  4.9109E+01  4.2907E+01
 EBVSHRINKSD(%)  3.3850E-01  1.3520E+01  9.9745E+01  3.2933E+01  2.0301E+01
 EBVSHRINKVR(%)  6.7586E-01  2.5212E+01  9.9999E+01  5.5020E+01  3.6481E+01
 RELATIVEINF(%)  9.9260E+01  8.0426E+00  1.1122E-04  3.7660E+00  1.9610E+01
 EPSSHRINKSD(%)  3.4004E+01
 EPSSHRINKVR(%)  5.6445E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2093.1145515504236     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1174.1760183457509     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.58
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2093.115       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.89E+00  4.08E-01  4.42E-01  1.31E+00  1.02E+00  8.18E-01  1.00E-02  1.40E+00  1.28E+00  1.04E+00
 


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
+        1.04E+03
 
 TH 2
+       -3.74E+00  3.65E+02
 
 TH 3
+        8.56E+00  1.12E+02  4.05E+02
 
 TH 4
+       -8.54E+00  3.57E+02 -4.41E+02  1.32E+03
 
 TH 5
+        1.69E+00 -7.83E+01 -2.24E+02  2.62E+02  2.58E+02
 
 TH 6
+        1.67E-02 -7.19E-01  2.81E+00 -3.08E+00  1.85E-01  1.88E+02
 
 TH 7
+        6.42E-01  7.84E+00 -1.86E+01 -1.94E+01 -1.64E+00 -4.92E-01  1.77E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.79E-01 -1.46E+01 -4.14E+01  7.29E+01  2.05E+00 -3.54E-01  1.67E+01  0.00E+00  3.25E+01
 
 TH10
+        1.40E+00 -1.06E+01 -2.63E+01  1.80E+01 -4.04E+01  5.23E-01  1.18E+00  0.00E+00  7.85E+00  5.35E+01
 
 TH11
+       -7.43E+00 -1.89E+01 -4.30E+01  1.38E+01 -3.08E+00  1.74E+00  1.26E+01  0.00E+00  4.64E+00  9.93E+00  3.75E+02
 
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
 #CPUT: Total CPU Time in Seconds,       43.150
Stop Time:
Thu Sep 30 02:32:55 CDT 2021
