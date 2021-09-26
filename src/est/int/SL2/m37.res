Sat Sep 25 01:11:44 CDT 2021
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
$DATA ../../../../data/int/SL2/dat37.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      999
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      899
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1457.39035384288        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.1230E+01  3.6495E+01  5.8902E+01  6.4949E+01  4.0714E+01  7.1065E+00 -1.0992E+02 -1.4001E+02 -1.1889E+02 -5.0841E+01
            -4.4347E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2852.85768574236        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0078E+00  1.2247E+00  1.0047E+00  8.5087E-01  1.1912E+00  9.5379E-01  1.3956E+00  1.0152E+00  1.1422E+00  1.3144E+00
             2.0039E+00
 PARAMETER:  1.0779E-01  3.0269E-01  1.0467E-01 -6.1498E-02  2.7496E-01  5.2685E-02  4.3332E-01  1.1508E-01  2.3297E-01  3.7339E-01
             7.9510E-01
 GRADIENT:  -5.8032E+00 -2.2949E+01 -1.5634E+01 -4.1154E+00  6.4135E+00 -1.2288E+01  2.8539E+01 -7.7807E+00  1.1042E+01  1.2738E+01
            -3.6918E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2871.65748314662        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0192E+00  1.7601E+00  1.0069E+00  6.4478E-01  1.5195E+00  1.0722E+00  1.1072E+00  1.8141E+00  1.0470E+00  1.5780E+00
             2.0691E+00
 PARAMETER:  1.1901E-01  6.6537E-01  1.0685E-01 -3.3884E-01  5.1835E-01  1.6969E-01  2.0187E-01  6.9562E-01  1.4596E-01  5.5616E-01
             8.2712E-01
 GRADIENT:   1.6386E+01  1.3655E+02 -3.8731E+00  6.8780E+01 -7.9646E+00  3.1447E+01  1.0238E+01 -4.3138E+00 -7.1902E+00  4.4694E+00
            -2.6433E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2891.18813922001        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.0159E+00  1.6072E+00  1.1798E+00  6.6416E-01  1.5201E+00  9.8058E-01  1.0198E+00  1.8352E+00  1.1492E+00  1.5097E+00
             2.2842E+00
 PARAMETER:  1.1578E-01  5.7447E-01  2.6537E-01 -3.0923E-01  5.1879E-01  8.0392E-02  1.1960E-01  7.0713E-01  2.3909E-01  5.1193E-01
             9.2604E-01
 GRADIENT:   1.8912E+00  1.8707E+01  1.8844E+00  1.4212E+01 -1.0389E+00 -2.3036E-02 -8.4160E-01 -3.2164E+00  2.1461E+00  1.8078E+00
            -1.0044E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2891.79343697066        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      314
 NPARAMETR:  1.0151E+00  1.6152E+00  1.2401E+00  6.5206E-01  1.5491E+00  9.8009E-01  1.0226E+00  2.1689E+00  1.1260E+00  1.5208E+00
             2.2884E+00
 PARAMETER:  1.1494E-01  5.7944E-01  3.1518E-01 -3.2762E-01  5.3764E-01  7.9893E-02  1.2236E-01  8.7421E-01  2.1868E-01  5.1923E-01
             9.2786E-01
 GRADIENT:   2.1618E-02  7.4697E+00 -1.2651E+00  8.3021E+00 -1.7730E+00 -1.3633E-01 -2.7156E-01 -1.4621E-01  2.0449E+00  1.2290E+00
            -8.8964E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2892.05947878464        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  1.0193E+00  1.6306E+00  1.2689E+00  6.3955E-01  1.5722E+00  9.8328E-01  1.0294E+00  2.1679E+00  1.0812E+00  1.5210E+00
             2.2895E+00
 PARAMETER:  1.1914E-01  5.8892E-01  3.3814E-01 -3.4699E-01  5.5247E-01  8.3134E-02  1.2901E-01  8.7377E-01  1.7810E-01  5.1939E-01
             9.2832E-01
 GRADIENT:  -1.0781E+00 -1.0724E+01  1.0721E+00 -3.6316E-01 -4.4715E+00  5.3562E-02 -5.1093E-01 -2.8066E+00 -1.2342E+00 -1.5019E+00
            -2.8704E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2892.51577539001        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      645
 NPARAMETR:  1.0197E+00  1.6831E+00  1.2057E+00  6.1289E-01  1.6051E+00  9.8295E-01  1.0140E+00  2.1679E+00  1.1386E+00  1.5578E+00
             2.2895E+00
 PARAMETER:  1.1954E-01  6.2062E-01  2.8708E-01 -3.8956E-01  5.7318E-01  8.2800E-02  1.1393E-01  8.7377E-01  2.2980E-01  5.4327E-01
             9.2833E-01
 GRADIENT:  -4.5261E-01 -1.7462E+00  7.1727E-02  6.7357E+00 -6.5565E-01 -1.8983E-01  1.6656E+00 -2.1622E+00  5.5948E-01  6.1463E-01
            -2.2447E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2892.87804927207        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  1.0199E+00  1.7326E+00  1.1408E+00  5.8286E-01  1.6349E+00  9.8300E-01  9.9337E-01  2.1679E+00  1.1833E+00  1.5850E+00
             2.2895E+00
 PARAMETER:  1.1975E-01  6.4964E-01  2.3169E-01 -4.3981E-01  5.9161E-01  8.2852E-02  9.3347E-02  8.7377E-01  2.6830E-01  5.6058E-01
             9.2833E-01
 GRADIENT:  -1.1630E-01  1.7199E+00 -2.5448E-01  8.2858E+00  9.0633E-01 -2.1982E-01  1.7644E+00 -1.8225E+00  9.5292E-01  1.1791E+00
            -2.0803E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2893.51044975947        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  1.0203E+00  1.8510E+00  9.8416E-01  5.0891E-01  1.7069E+00  9.8346E-01  9.4686E-01  2.1679E+00  1.3006E+00  1.6430E+00
             2.2895E+00
 PARAMETER:  1.2014E-01  7.1573E-01  8.4033E-02 -5.7549E-01  6.3470E-01  8.3319E-02  4.5394E-02  8.7377E-01  3.6282E-01  5.9654E-01
             9.2833E-01
 GRADIENT:   5.0760E-01  8.6894E+00 -7.4547E-01  9.6584E+00  3.3141E+00 -1.3786E-01  1.2653E+00 -1.3582E+00  1.2714E+00  1.8012E+00
            -2.0577E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2893.54725468988        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1190
 NPARAMETR:  1.0204E+00  1.8597E+00  9.7228E-01  5.0319E-01  1.7125E+00  9.8341E-01  9.4377E-01  2.1679E+00  1.3108E+00  1.6468E+00
             2.2895E+00
 PARAMETER:  1.2016E-01  7.2044E-01  7.1884E-02 -5.8679E-01  6.3796E-01  8.3269E-02  4.2123E-02  8.7376E-01  3.7062E-01  5.9886E-01
             9.2833E-01
 GRADIENT:   5.3309E-01  8.9052E+00 -7.8027E-01  9.5870E+00  3.4855E+00 -1.6235E-01  1.2656E+00 -1.3254E+00  1.3033E+00  1.7792E+00
            -2.0568E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2893.54991714133        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1335
 NPARAMETR:  1.0204E+00  1.8597E+00  9.7232E-01  5.0334E-01  1.7120E+00  9.8336E-01  9.3779E-01  2.1670E+00  1.3105E+00  1.6463E+00
             2.2906E+00
 PARAMETER:  1.2022E-01  7.2044E-01  7.1934E-02 -5.8649E-01  6.3764E-01  8.3220E-02  3.5772E-02  8.7333E-01  3.7043E-01  5.9856E-01
             9.2880E-01
 GRADIENT:   6.0061E-01  9.0108E+00 -7.3858E-01  9.7634E+00  3.3054E+00 -1.6665E-01 -1.5646E-01 -1.3503E+00  1.0587E+00  1.6991E+00
            -1.1497E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2893.56767448834        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:     1458
 NPARAMETR:  1.0198E+00  1.8597E+00  9.7627E-01  5.0333E-01  1.7120E+00  9.8340E-01  9.3836E-01  2.1832E+00  1.3067E+00  1.6441E+00
             2.2905E+00
 PARAMETER:  1.1960E-01  7.2043E-01  7.5989E-02 -5.8650E-01  6.3765E-01  8.3260E-02  3.6376E-02  8.8079E-01  3.6753E-01  5.9717E-01
             9.2878E-01
 GRADIENT:   1.0298E+01  3.6825E+01 -6.7458E-01  1.2781E+01  7.0517E+00  1.0208E+00  1.9360E-01 -1.1979E+00  1.2292E+00  2.2810E+00
             5.6343E-01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -2893.56767448834        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1532
 NPARAMETR:  1.0198E+00  1.8597E+00  9.7627E-01  5.0333E-01  1.7120E+00  9.8346E-01  9.3845E-01  2.1831E+00  1.3067E+00  1.6441E+00
             2.2906E+00
 PARAMETER:  1.1960E-01  7.2043E-01  7.5989E-02 -5.8650E-01  6.3765E-01  8.3260E-02  3.6376E-02  8.8079E-01  3.6753E-01  5.9717E-01
             9.2878E-01
 GRADIENT:  -4.5418E-01  1.4050E+00  4.4228E+05 -8.5888E+00  2.3328E+00 -1.4413E-01 -1.2160E-01  5.0136E+04 -1.2034E+05  9.5554E-01
            -2.3814E+04
 NUMSIGDIG:         9.5         8.2         3.3         7.5         8.1         2.1         1.9         3.3         3.3         8.5
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1532
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5148E-03 -2.2975E-02 -2.0307E-02  1.6528E-02 -2.3304E-02
 SE:             2.9546E-02  2.5540E-02  1.1468E-02  1.8183E-02  2.5673E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5911E-01  3.6836E-01  7.6589E-02  3.6336E-01  3.6401E-01

 ETASHRINKSD(%)  1.0163E+00  1.4438E+01  6.1582E+01  3.9086E+01  1.3994E+01
 ETASHRINKVR(%)  2.0223E+00  2.6791E+01  8.5240E+01  6.2894E+01  2.6029E+01
 EBVSHRINKSD(%)  1.2687E+00  1.4975E+01  6.5908E+01  4.3334E+01  1.0330E+01
 EBVSHRINKVR(%)  2.5213E+00  2.7708E+01  8.8377E+01  6.7890E+01  1.9593E+01
 RELATIVEINF(%)  9.7447E+01  1.0784E+01  5.8420E+00  4.1700E+00  3.3843E+01
 EPSSHRINKSD(%)  1.7546E+01
 EPSSHRINKVR(%)  3.2013E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          899
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1652.2514827020016     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2893.5676744883449     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1241.3161917863433     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2893.568       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.86E+00  9.76E-01  5.03E-01  1.71E+00  9.83E-01  9.38E-01  2.18E+00  1.31E+00  1.64E+00  2.29E+00
 


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
+        1.49E+09
 
 TH 2
+       -7.62E+00  1.23E+07
 
 TH 3
+        2.98E+00  1.36E+01  1.16E+09
 
 TH 4
+       -1.43E+01  3.25E+02 -7.11E+01  2.54E+08
 
 TH 5
+       -1.98E+00 -3.18E+01 -2.46E+01 -3.43E+07  1.86E+07
 
 TH 6
+        4.83E+00 -2.07E+00  1.14E+03 -4.88E+00 -1.37E+00  1.95E+02
 
 TH 7
+        6.15E+00  8.48E+00  1.49E+04 -1.88E+01  1.28E+00 -1.91E+00  1.26E+09
 
 TH 8
+       -4.72E+07 -3.34E+00 -5.90E+07  1.92E+01 -1.61E+00  5.79E+01 -6.14E+07  5.98E+06
 
 TH 9
+       -1.89E+08 -7.48E+00 -2.36E+08  2.73E+01  6.41E-01 -2.32E+02 -3.01E+03  1.20E+07  9.59E+07
 
 TH10
+       -9.23E+07 -4.23E+00 -1.15E+08  1.56E+01 -8.59E+00 -1.86E-02  1.75E+00  5.87E+06  2.35E+07  2.29E+07
 
 TH11
+        4.26E+07 -1.03E+01 -3.24E+03 -6.31E+00 -1.94E-02 -4.97E+01 -6.79E+02  8.18E+03  1.01E+03 -5.29E+06  2.44E+06
 
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
 #CPUT: Total CPU Time in Seconds,       55.959
Stop Time:
Sat Sep 25 01:12:41 CDT 2021
