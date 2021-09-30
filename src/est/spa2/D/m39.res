Thu Sep 30 09:01:18 CDT 2021
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
$DATA ../../../../data/spa2/D/dat39.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20003.4232611949        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.1850E+02  3.4290E+02  4.3803E+01  1.9503E+02  2.0581E+01 -1.5706E+03 -7.9132E+02 -3.3542E+01 -1.2319E+03 -4.3081E+02
            -4.0422E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -632.013499054046        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3212E+00  1.3040E+00  9.8389E-01  1.6011E+00  9.9424E-01  1.9590E+00  1.3128E+00  9.8098E-01  1.2797E+00  1.0454E+00
             1.4473E+01
 PARAMETER:  3.7856E-01  3.6545E-01  8.3757E-02  5.7071E-01  9.4219E-02  7.7246E-01  3.7215E-01  8.0797E-02  3.4666E-01  1.4435E-01
             2.7723E+00
 GRADIENT:  -1.6370E+01  3.1750E+01 -1.5073E+00  6.1437E+01 -1.0147E+01  3.4418E+01 -1.7016E+01  3.9025E+00 -8.5929E+00  1.1504E+01
             2.6096E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -703.566508696313        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.3104E+00  1.5875E+00  2.4017E+00  1.3905E+00  1.8526E+00  1.8619E+00  3.7304E+00  1.2159E-01  1.8613E+00  9.6042E-01
             1.3098E+01
 PARAMETER:  3.7033E-01  5.6216E-01  9.7619E-01  4.2966E-01  7.1657E-01  7.2157E-01  1.4165E+00 -2.0071E+00  7.2130E-01  5.9615E-02
             2.6725E+00
 GRADIENT:   2.8602E+00 -4.6281E+00 -1.1829E+01  8.7499E+00  4.5359E+00 -9.0060E+00 -5.9310E+00  1.6414E-02  3.2251E+01  7.2360E+00
             3.3516E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -785.278616253301        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0435E+00  8.1668E-01  1.1419E+01  1.3850E+00  2.9517E+00  2.1226E+00  4.5594E+00  8.1672E-01  9.2953E-01  5.6965E-01
             9.4335E+00
 PARAMETER:  1.4255E-01 -1.0250E-01  2.5353E+00  4.2569E-01  1.1824E+00  8.5264E-01  1.6172E+00 -1.0245E-01  2.6920E-02 -4.6274E-01
             2.3443E+00
 GRADIENT:  -4.3868E+01 -1.8380E+01 -8.0765E-01 -1.5153E+01  1.7717E+01  2.8693E+01  4.0511E+00  2.8023E-02  1.9238E+00  1.5256E+00
             1.1799E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -792.626184016952        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.1456E+00  1.3930E+00  2.8396E+00  1.0175E+00  1.8626E+00  1.9973E+00  3.6482E+00  6.3760E-01  6.7634E-01  3.3507E-01
             8.5622E+00
 PARAMETER:  2.3596E-01  4.3147E-01  1.1437E+00  1.1733E-01  7.2200E-01  7.9182E-01  1.3942E+00 -3.5004E-01 -2.9107E-01 -9.9341E-01
             2.2474E+00
 GRADIENT:   1.3121E+01  1.0578E+00  2.2168E+00 -2.8427E-01 -5.1021E+00 -5.9290E-01 -2.0459E+00  2.9478E-01  7.0031E-01  1.2862E+00
            -9.1710E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -793.029213167481        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.1167E+00  1.4325E+00  1.8513E+00  9.8586E-01  1.6601E+00  1.9970E+00  3.6035E+00  1.4660E-01  5.9718E-01  2.2964E-01
             8.6126E+00
 PARAMETER:  2.1034E-01  4.5945E-01  7.1587E-01  8.5759E-02  6.0689E-01  7.9165E-01  1.3819E+00 -1.8200E+00 -4.1554E-01 -1.3712E+00
             2.2532E+00
 GRADIENT:  -3.0651E+00  2.1091E+00  6.6645E-01  3.8374E-01 -1.9889E+00 -7.9398E-01 -3.0185E-01  3.5286E-02  1.1810E+00  7.2341E-01
             1.8456E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -793.048953315517        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.1182E+00  1.4240E+00  1.6525E+00  9.7121E-01  1.6101E+00  1.9988E+00  3.5905E+00  8.9874E-02  5.4170E-01  1.7953E-01
             8.5952E+00
 PARAMETER:  2.1170E-01  4.5350E-01  6.0229E-01  7.0792E-02  5.7628E-01  7.9253E-01  1.3783E+00 -2.3093E+00 -5.1304E-01 -1.6174E+00
             2.2512E+00
 GRADIENT:  -1.9122E+00 -2.8912E-01 -6.7283E-01 -6.7134E-01  4.3409E-01 -5.1954E-01  5.5771E-01  1.5708E-02  1.0415E+00  4.6832E-01
             3.1078E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -793.061890793525        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      511
 NPARAMETR:  1.1218E+00  1.4366E+00  1.5437E+00  9.5052E-01  1.5901E+00  2.0009E+00  3.5610E+00  5.8091E-02  4.6623E-01  1.2402E-01
             8.5843E+00
 PARAMETER:  2.1498E-01  4.6225E-01  5.3418E-01  4.9250E-02  5.6382E-01  7.9362E-01  1.3701E+00 -2.7457E+00 -6.6307E-01 -1.9873E+00
             2.2499E+00
 GRADIENT:   2.8704E-01 -1.7375E+00 -1.1877E+00 -1.5109E+00  1.6904E+00 -9.6237E-02  8.1052E-01  6.9679E-03  6.9465E-01  2.3084E-01
            -8.6010E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -793.077350107220        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      581
 NPARAMETR:  1.1231E+00  1.4683E+00  1.4886E+00  9.2974E-01  1.5856E+00  2.0019E+00  3.5244E+00  3.8775E-02  3.8285E-01  7.6977E-02
             8.5881E+00
 PARAMETER:  2.1609E-01  4.8410E-01  4.9782E-01  2.7151E-02  5.6096E-01  7.9409E-01  1.3597E+00 -3.1500E+00 -8.6010E-01 -2.4643E+00
             2.2504E+00
 GRADIENT:   9.0613E-01 -1.3570E+00 -7.8166E-01 -1.3505E+00  1.2835E+00  1.1351E-01  5.2274E-01  3.0913E-03  3.5082E-01  9.0334E-02
            -5.4161E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -796.676169847613        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      717
 NPARAMETR:  1.1637E+00  1.2614E+00  2.4482E+00  1.1165E+00  1.7894E+00  2.1354E+00  4.3350E+00  1.2028E-01  6.7898E-01  2.3980E-01
             8.8562E+00
 PARAMETER:  2.5163E-01  3.3223E-01  9.9534E-01  2.1016E-01  6.8187E-01  8.5867E-01  1.5667E+00 -2.0180E+00 -2.8717E-01 -1.3280E+00
             2.2811E+00
 GRADIENT:  -9.9695E-01  2.3048E-01 -1.6345E+00 -9.1156E+00  2.6677E+00  2.4473E+00  2.2934E+00  1.9182E-02  1.1700E+00  6.5538E-01
             7.7956E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -798.092121076579        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  1.1607E+00  8.7993E-01  5.3998E+00  1.3726E+00  2.0709E+00  2.1048E+00  5.1638E+00  2.9415E-01  8.9879E-01  2.3764E-01
             8.8437E+00
 PARAMETER:  2.4900E-01 -2.7918E-02  1.7864E+00  4.1668E-01  8.2800E-01  8.4423E-01  1.7417E+00 -1.1237E+00 -6.7049E-03 -1.3370E+00
             2.2797E+00
 GRADIENT:  -1.2075E+00  8.8701E-01  2.4252E-01  3.2298E-01  2.4448E-02 -1.6028E-01  5.5918E-01  3.1620E-02 -1.3269E+00  4.9707E-01
            -1.5625E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -798.444597648668        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1070
 NPARAMETR:  1.1722E+00  8.1967E-01  5.3079E+00  1.4093E+00  2.0354E+00  2.1093E+00  5.2654E+00  1.0463E-01  9.6423E-01  5.4635E-02
             8.8439E+00
 PARAMETER:  2.5886E-01 -9.8858E-02  1.7692E+00  4.4308E-01  8.1071E-01  8.4636E-01  1.7611E+00 -2.1574E+00  6.3573E-02 -2.8071E+00
             2.2797E+00
 GRADIENT:   3.1955E+00  1.2270E-01  1.2269E-01  1.5342E-02 -3.7054E-01  3.9017E-01 -3.7666E-01  4.4992E-03  1.4644E-03  2.6698E-02
            -1.7448E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -798.503633499832        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  1.1651E+00  7.8530E-01  5.3709E+00  1.4280E+00  2.0373E+00  2.1103E+00  5.3821E+00  2.8584E-02  9.6528E-01  1.0602E-02
             8.8584E+00
 PARAMETER:  2.5282E-01 -1.4169E-01  1.7810E+00  4.5629E-01  8.1163E-01  8.4685E-01  1.7831E+00 -3.4549E+00  6.4658E-02 -4.4467E+00
             2.2814E+00
 GRADIENT:   4.5711E-01 -5.5345E-02  1.6903E-02  1.8822E-01 -4.5981E-02  8.8833E-01  5.6660E-01  3.3997E-04 -4.0587E-01  9.9786E-04
             4.8990E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -798.510407127333        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1431             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1643E+00  7.8345E-01  5.3963E+00  1.4305E+00  2.0383E+00  2.1078E+00  5.4141E+00  1.0000E-02  9.7708E-01  1.0000E-02
             8.8533E+00
 PARAMETER:  2.5213E-01 -1.4405E-01  1.7857E+00  4.5800E-01  8.1212E-01  8.4566E-01  1.7890E+00 -4.9110E+00  7.6814E-02 -4.6038E+00
             2.2808E+00
 GRADIENT:   1.7720E+01  1.1128E+00  8.6391E-02  1.1407E+01  1.0486E+00  2.0181E+01  4.3626E+01  0.0000E+00  3.5801E-01  0.0000E+00
             2.7572E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -798.513471531385        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1609
 NPARAMETR:  1.1645E+00  7.7766E-01  5.4151E+00  1.4339E+00  2.0383E+00  2.1082E+00  5.4257E+00  1.0000E-02  9.7673E-01  1.0000E-02
             8.8545E+00
 PARAMETER:  2.5228E-01 -1.5147E-01  1.7892E+00  4.6039E-01  8.1210E-01  8.4582E-01  1.7911E+00 -4.9055E+00  7.6455E-02 -4.6038E+00
             2.2809E+00
 GRADIENT:   3.1792E-01  1.2747E-01 -2.0317E-02 -5.0368E-01  1.2399E-02  5.4953E-01  1.5245E+00  0.0000E+00 -1.8060E-03  0.0000E+00
             1.6097E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -798.517493172155        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1798             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1638E+00  7.6618E-01  5.4534E+00  1.4391E+00  2.0393E+00  2.1078E+00  5.4524E+00  1.0000E-02  9.8018E-01  1.0000E-02
             8.8521E+00
 PARAMETER:  2.5167E-01 -1.6633E-01  1.7962E+00  4.6402E-01  8.1260E-01  8.4563E-01  1.7961E+00 -4.9055E+00  7.9977E-02 -4.6038E+00
             2.2807E+00
 GRADIENT:   1.7553E+01  9.5074E-01  7.9957E-02  1.1890E+01  1.0327E+00  2.0251E+01  4.4061E+01  0.0000E+00  1.7706E-01  0.0000E+00
             2.7305E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -798.518349352901        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1978
 NPARAMETR:  1.1644E+00  7.6675E-01  5.4753E+00  1.4413E+00  2.0399E+00  2.1080E+00  5.4586E+00  1.0000E-02  9.8202E-01  1.0000E-02
             8.8553E+00
 PARAMETER:  2.5220E-01 -1.6559E-01  1.8002E+00  4.6556E-01  8.1290E-01  8.4576E-01  1.7972E+00 -4.9055E+00  8.1860E-02 -4.6038E+00
             2.2810E+00
 GRADIENT:   2.8834E-01  1.0872E-01 -1.6404E-02 -2.2330E-01 -7.0771E-02  5.8438E-01  1.6083E+00  0.0000E+00 -4.1218E-02  0.0000E+00
             7.7258E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -798.520184440438        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2167
 NPARAMETR:  1.1638E+00  7.6207E-01  5.5187E+00  1.4444E+00  2.0419E+00  2.1075E+00  5.4815E+00  1.0000E-02  9.8538E-01  1.0000E-02
             8.8523E+00
 PARAMETER:  2.5169E-01 -1.7172E-01  1.8081E+00  4.6769E-01  8.1390E-01  8.4549E-01  1.8014E+00 -4.9055E+00  8.5267E-02 -4.6038E+00
             2.2807E+00
 GRADIENT:   1.3064E-01  1.6755E-01 -1.3985E-02 -4.0281E-01 -7.1045E-02  5.0104E-01  1.9867E+00  0.0000E+00  2.5268E-02  0.0000E+00
            -2.8536E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -798.521035957138        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2355
 NPARAMETR:  1.1642E+00  7.5897E-01  5.5395E+00  1.4466E+00  2.0431E+00  2.1077E+00  5.4847E+00  1.0000E-02  9.8637E-01  1.0000E-02
             8.8548E+00
 PARAMETER:  2.5203E-01 -1.7579E-01  1.8119E+00  4.6921E-01  8.1447E-01  8.4560E-01  1.8020E+00 -4.9055E+00  8.6277E-02 -4.6038E+00
             2.2810E+00
 GRADIENT:   2.3768E-01  1.0085E-01 -1.2381E-02 -1.8387E-01 -8.5408E-02  5.7124E-01  1.7712E+00  0.0000E+00 -3.2429E-02  0.0000E+00
            -6.4403E-02

0ITERATION NO.:   93    OBJECTIVE VALUE:  -798.521339259623        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     2460
 NPARAMETR:  1.1637E+00  7.5785E-01  5.5946E+00  1.4500E+00  2.0428E+00  2.1072E+00  5.4987E+00  1.0000E-02  9.8612E-01  1.0000E-02
             8.8504E+00
 PARAMETER:  2.5184E-01 -1.7896E-01  1.8131E+00  4.6938E-01  8.1499E-01  8.4536E-01  1.8024E+00 -4.9055E+00  8.7021E-02 -4.6038E+00
             2.2808E+00
 GRADIENT:   3.1263E-02 -2.8396E-02 -3.0112E-02 -4.2652E-01  3.3955E-02  4.3898E-05 -1.3576E-01  0.0000E+00  1.7683E-02  0.0000E+00
             1.7027E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2460
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3123E-02  4.0643E-02  5.8087E-06 -7.1343E-02  9.7771E-06
 SE:             2.8902E-02  2.3004E-02  9.2117E-06  1.3567E-02  5.3522E-05
 N:                     100         100         100         100         100

 P VAL.:         6.4979E-01  7.7266E-02  5.2832E-01  1.4532E-07  8.5505E-01

 ETASHRINKSD(%)  3.1743E+00  2.2933E+01  9.9969E+01  5.4550E+01  9.9821E+01
 ETASHRINKVR(%)  6.2479E+00  4.0606E+01  1.0000E+02  7.9343E+01  1.0000E+02
 EBVSHRINKSD(%)  4.1571E+00  1.8987E+01  9.9944E+01  5.5711E+01  9.9727E+01
 EBVSHRINKVR(%)  8.1414E+00  3.4370E+01  1.0000E+02  8.0385E+01  9.9999E+01
 RELATIVEINF(%)  9.1181E+01  2.9035E+01  3.1505E-06  8.2738E+00  7.4053E-05
 EPSSHRINKSD(%)  9.5196E+00
 EPSSHRINKVR(%)  1.8133E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -798.52133925962278     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       304.20490058598432     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    53.15
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    11.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -798.521       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.16E+00  7.57E-01  5.55E+00  1.45E+00  2.04E+00  2.11E+00  5.49E+00  1.00E-02  9.87E-01  1.00E-02  8.85E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.69E+02
 
 TH 2
+       -1.08E+02  5.59E+01
 
 TH 3
+       -3.51E-02  1.12E-01  4.48E-04
 
 TH 4
+       -2.52E+02  1.21E+02  2.05E-01  3.00E+02
 
 TH 5
+        1.75E+01 -9.92E+00 -2.05E-02 -2.39E+01  1.97E+00
 
 TH 6
+       -7.01E+01  3.31E+01  6.82E-02  5.51E+01 -4.51E+00  2.86E+01
 
 TH 7
+        6.14E+00 -1.95E+00  1.12E-03 -8.98E+00  6.54E-01  1.04E+00  6.90E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.25E+01 -2.93E+01 -4.48E-02 -7.74E+01  6.11E+00 -1.08E+01  2.83E+00  0.00E+00  2.06E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        4.60E+00 -4.15E+00 -8.76E-03 -1.18E+01  9.82E-01 -7.62E-01  4.99E-01  0.00E+00  3.28E+00  0.00E+00  9.81E-01
 
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
+        1.63E+02
 
 TH 2
+       -2.84E+00  4.38E+01
 
 TH 3
+        7.24E-02  2.84E-01  1.67E-01
 
 TH 4
+       -7.13E+00  4.14E+01 -1.45E-01  1.40E+02
 
 TH 5
+       -7.04E-01 -6.70E+00 -1.47E+00 -1.01E+01  1.85E+01
 
 TH 6
+       -2.68E+00 -6.91E-01  5.24E-02  1.47E+00 -1.29E-02  3.37E+01
 
 TH 7
+        7.20E-01  5.12E+00 -3.52E-02 -8.06E+00  5.91E-01 -1.42E-01  3.25E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.05E-01 -3.71E+00 -4.08E-01 -3.58E+01  4.06E+00 -1.01E+00  2.32E+00  0.00E+00  2.74E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.38E+00 -3.42E+00  1.67E-02 -1.01E+01  5.56E-01  1.86E+00  2.54E-01  0.00E+00  3.05E+00  0.00E+00  8.33E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.70E+02
 
 TH 2
+        6.51E+01  4.42E+01
 
 TH 3
+        1.00E+00  4.00E-01  4.14E-02
 
 TH 4
+        9.00E+01  4.24E+01  7.37E-01  1.45E+02
 
 TH 5
+       -1.92E+01 -7.13E+00 -5.43E-01 -1.84E+01  8.26E+00
 
 TH 6
+        7.96E+00  7.99E+00 -1.22E-01 -2.30E+01  1.66E+00  3.18E+01
 
 TH 7
+        1.37E+00  5.27E+00 -5.34E-02 -8.84E+00  1.34E+00  4.46E+00  3.50E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.40E+01 -4.37E+00 -1.65E-01 -3.86E+01  5.67E+00  8.85E+00  2.46E+00  0.00E+00  2.25E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.44E+01 -8.26E+00 -2.51E-01 -3.43E+01  5.14E+00  6.01E+00  3.95E+00  0.00E+00  9.59E+00  0.00E+00  1.78E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       64.832
Stop Time:
Thu Sep 30 09:02:25 CDT 2021
