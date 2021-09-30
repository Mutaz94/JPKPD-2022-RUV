Thu Sep 30 03:36:26 CDT 2021
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
$DATA ../../../../data/spa1/D/dat75.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17897.1661474681        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1719E+02  2.8051E+02 -1.2527E+02 -1.1996E+01  3.1617E+02 -2.1318E+03 -9.4894E+02 -7.6556E+01 -1.9150E+03 -5.1548E+02
            -3.4118E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -625.599718818745        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.7148E+00  1.0829E+00  9.6592E-01  2.7696E+00  1.2153E+00  3.4434E+00  1.4026E+00  9.5497E-01  3.4295E+00  1.0816E+00
             1.1783E+01
 PARAMETER:  6.3927E-01  1.7962E-01  6.5325E-02  1.1187E+00  2.9501E-01  1.3365E+00  4.3835E-01  5.3924E-02  1.3324E+00  1.7842E-01
             2.5667E+00
 GRADIENT:   6.3558E+01  2.3589E+01 -5.3252E+01  9.3148E+01  3.3235E+00  1.0747E+02  4.1769E+00  7.2072E+00  1.9856E+01  4.5706E+00
             9.5362E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -683.377858043322        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4482E+00  1.8001E+00  3.5209E+00  1.9302E+00  6.4546E+00  2.7612E+00  1.3995E+00  3.5830E-01  3.9676E+00  8.5191E+00
             1.1075E+01
 PARAMETER:  4.7034E-01  6.8784E-01  1.3587E+00  7.5762E-01  1.9648E+00  1.1157E+00  4.3612E-01 -9.2639E-01  1.4782E+00  2.2423E+00
             2.5047E+00
 GRADIENT:   3.4575E+01  5.9997E+01  5.9600E+00  5.7625E+01 -3.9424E+00  4.6752E+01  7.6236E-01 -4.4229E-02 -4.4676E+00  1.3997E+01
             9.4770E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -709.275684629415        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.1686E+00  1.0863E+00  2.1601E+00  1.4377E+00  4.1512E+00  2.2784E+00  1.0975E+00  1.4129E+00  3.2825E+00  4.9998E+00
             1.0450E+01
 PARAMETER:  2.5584E-01  1.8279E-01  8.7016E-01  4.6302E-01  1.5234E+00  9.2346E-01  1.9302E-01  4.4561E-01  1.2886E+00  1.7094E+00
             2.4466E+00
 GRADIENT:  -3.1254E+01 -2.3045E+00  1.1523E+01 -5.7362E+00 -1.0257E+01 -1.1517E+01  5.9735E+00 -1.1575E+00  1.9220E+01  1.6584E+00
             1.2963E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -742.279943444864        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  1.1497E+00  4.4613E-01  5.1163E-01  1.4051E+00  1.8300E+01  2.2114E+00  5.3445E-01  5.0046E-01  2.2169E+00  1.3266E+01
             8.5035E+00
 PARAMETER:  2.3954E-01 -7.0714E-01 -5.7015E-01  4.4012E-01  3.0069E+00  8.9364E-01 -5.2651E-01 -5.9223E-01  8.9610E-01  2.6852E+00
             2.2405E+00
 GRADIENT:   4.4038E+01  3.2320E+01 -2.3434E+01 -5.1564E+00 -4.0709E+00  5.1955E+00  1.8313E+00  1.1359E+00  1.4416E+01  1.2914E+01
            -7.9770E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -831.888111250919        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  7.2923E-01  4.6532E-02  1.4567E-01  8.3203E-01  6.6473E+01  1.6863E+00  1.0000E-02  2.5192E+00  9.8808E-01  6.7436E+00
             8.2649E+00
 PARAMETER: -2.1576E-01 -2.9676E+00 -1.8264E+00 -8.3886E-02  4.2968E+00  6.2256E-01 -6.0209E+00  1.0239E+00  8.8004E-02  2.0086E+00
             2.2120E+00
 GRADIENT:  -6.5865E+00 -2.6054E-01  5.5466E+01 -2.3472E+01 -2.2429E-02  7.1357E+00  0.0000E+00  1.1138E+01  3.2846E+01  3.8436E-03
            -6.1259E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -846.617080183333        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  5.5208E-01  1.4767E-02  5.0950E-02  5.1668E-01  6.9580E+01  1.4847E+00  1.0000E-02  1.8925E+00  6.2417E-01  3.0099E+00
             8.5855E+00
 PARAMETER: -4.9407E-01 -4.1154E+00 -2.8769E+00 -5.6034E-01  4.3425E+00  4.9523E-01 -9.3011E+00  7.3791E-01 -3.7134E-01  1.2019E+00
             2.2501E+00
 GRADIENT:  -5.8000E+01  6.8724E-02 -2.8845E+01  6.4751E+01 -2.1864E-02 -4.2575E+01  0.0000E+00  1.7559E+01  1.0803E+01 -3.4291E-05
            -4.0691E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -857.747025310069        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  5.5167E-01  1.1791E-02  4.2548E-02  4.3449E-01  5.2728E+01  1.6898E+00  1.0000E-02  1.4994E+00  4.6871E-01  2.3262E+00
             8.9911E+00
 PARAMETER: -4.9480E-01 -4.3404E+00 -3.0571E+00 -7.3357E-01  4.0651E+00  6.2461E-01 -9.9234E+00  5.0504E-01 -6.5777E-01  9.4422E-01
             2.2962E+00
 GRADIENT:  -1.8934E+00 -1.3556E-01 -2.2993E+00  1.1096E+00  2.1796E-02  7.8983E-01  0.0000E+00  1.0271E+00  3.6762E+00 -5.0784E-05
             4.1142E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -858.686459380774        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      826
 NPARAMETR:  5.8756E-01  1.4410E-02  5.0641E-02  4.8660E-01  6.0580E+01  1.6882E+00  1.0000E-02  1.6608E+00  2.6432E-01  4.1303E+00
             8.9556E+00
 PARAMETER: -4.3178E-01 -4.1398E+00 -2.8830E+00 -6.2032E-01  4.2040E+00  6.2368E-01 -7.6214E+00  6.0728E-01 -1.2306E+00  1.5183E+00
             2.2923E+00
 GRADIENT:   3.0077E+00 -1.9847E-01 -1.4776E+00  1.0697E+00  3.1663E-02 -1.3606E-01  0.0000E+00 -1.3553E+00  8.7320E-01 -3.7637E-04
            -2.8005E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -859.021966471112        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  5.9366E-01  1.5174E-02  5.3275E-02  5.0084E-01  7.9637E+01  1.6877E+00  3.9002E-02  1.7470E+00  4.9761E-02  1.1889E+01
             8.9847E+00
 PARAMETER: -4.2145E-01 -4.0882E+00 -2.8323E+00 -5.9147E-01  4.4775E+00  6.2339E-01 -3.1441E+00  6.5789E-01 -2.9005E+00  2.5756E+00
             2.2955E+00
 GRADIENT:  -3.2100E-01 -2.1602E-01  1.0186E+00 -1.4264E+00  2.4572E-02  2.7714E-01  3.1524E-07  3.3988E-01  2.9015E-02 -2.7300E-03
             1.2015E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -859.089942257638        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1171
 NPARAMETR:  5.9344E-01  1.7534E-02  5.3168E-02  5.0017E-01  6.8551E+01  1.6853E+00  2.2185E-02  1.7451E+00  1.5655E-02  1.3000E+01
             8.9839E+00
 PARAMETER: -4.2183E-01 -3.9436E+00 -2.8343E+00 -5.9281E-01  4.3276E+00  6.2194E-01 -3.7083E+00  6.5679E-01 -4.0570E+00  2.6649E+00
             2.2954E+00
 GRADIENT:   3.6604E+01 -5.1591E-01  5.7017E+01  1.6966E+01  3.6328E-02  1.2786E+01  1.7797E-06  3.6026E+00  6.0553E-03 -2.2478E-03
             2.3983E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -859.308824279692        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1310
 NPARAMETR:  5.8888E-01  2.2404E-02  5.2984E-02  4.9959E-01  6.9051E+01  1.6727E+00  2.1123E-02  1.7452E+00  1.0000E-02  1.2941E+01
             8.9881E+00
 PARAMETER: -4.2954E-01 -3.6985E+00 -2.8378E+00 -5.9396E-01  4.3348E+00  6.1444E-01 -3.7574E+00  6.5688E-01 -6.3177E+00  2.6604E+00
             2.2959E+00
 GRADIENT:  -7.7560E+00 -1.7225E+00  5.1420E-01  5.3345E+00  4.3408E-02 -2.4398E+00  4.9487E-05 -1.3335E+00  0.0000E+00 -1.7165E-03
             2.3457E+00

0ITERATION NO.:   58    OBJECTIVE VALUE:  -859.357120629035        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:     1412
 NPARAMETR:  5.8803E-01  2.3083E-02  5.3056E-02  4.9932E-01  6.9317E+01  1.6719E+00  1.9766E-02  1.7447E+00  1.0000E-02  1.2910E+01
             8.9756E+00
 PARAMETER: -4.3070E-01 -3.6663E+00 -2.8382E+00 -5.9412E-01  4.3358E+00  6.1353E-01 -3.7859E+00  6.5699E-01 -6.5094E+00  2.6598E+00
             2.2960E+00
 GRADIENT:   6.7768E+02  1.5407E+02 -2.3105E+02  7.1965E+02 -1.3079E+02 -6.9189E+02  6.5744E-05  1.2878E+03  0.0000E+00  2.1115E+02
             2.4270E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1412
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5707E-02 -3.5501E-05  5.3322E-03 -4.4943E-04 -4.7550E-04
 SE:             2.9122E-02  7.9672E-06  2.5084E-02  2.3791E-04  3.3064E-04
 N:                     100         100         100         100         100

 P VAL.:         5.8963E-01  8.3625E-06  8.3166E-01  5.8885E-02  1.5040E-01

 ETASHRINKSD(%)  2.4380E+00  9.9973E+01  1.5964E+01  9.9203E+01  9.8892E+01
 ETASHRINKVR(%)  4.8166E+00  1.0000E+02  2.9380E+01  9.9994E+01  9.9988E+01
 EBVSHRINKSD(%)  3.4607E+00  9.9942E+01  1.4656E+01  9.9207E+01  9.9139E+01
 EBVSHRINKVR(%)  6.8016E+00  1.0000E+02  2.7164E+01  9.9994E+01  9.9993E+01
 RELATIVEINF(%)  6.3398E+00  1.0666E-05  2.0570E+00  1.1312E-04  7.3864E-04
 EPSSHRINKSD(%)  1.0940E+01
 EPSSHRINKVR(%)  2.0684E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -859.35712062903451     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       59.581412575638183     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -859.357       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.88E-01  2.31E-02  5.30E-02  5.00E-01  6.91E+01  1.67E+00  2.05E-02  1.75E+00  1.00E-02  1.29E+01  8.99E+00
 


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
+        2.23E+05
 
 TH 2
+        3.17E+05  2.95E+06
 
 TH 3
+        6.42E+02 -3.04E+05  1.48E+06
 
 TH 4
+       -2.18E+03  3.25E+04 -2.52E+05  2.79E+05
 
 TH 5
+       -1.72E+02 -5.40E+02  2.16E+02 -1.53E+02  1.57E-01
 
 TH 6
+        5.18E+02 -7.93E+04 -3.60E+02 -2.28E+04  4.42E+01  2.66E+04
 
 TH 7
+        2.17E-02  1.01E+00  8.22E-01 -2.50E-02 -2.22E-05 -1.50E-03  2.02E-01
 
 TH 8
+       -4.38E+02 -5.29E+03  1.77E+04  2.03E+04 -1.99E+01 -5.84E+03  1.47E-03  1.61E+04
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        7.68E+02  4.65E+03 -1.88E+03  2.30E+03 -6.02E-01 -3.75E+02  3.02E-04  1.71E+02  0.00E+00  2.00E+01
 
 TH11
+        1.29E+03  3.95E+03  1.17E+03  6.30E+02 -2.21E+00 -3.22E+02 -1.05E-03  2.93E+02  0.00E+00  2.42E+01  6.28E+01
 
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
 #CPUT: Total CPU Time in Seconds,       35.040
Stop Time:
Thu Sep 30 03:37:03 CDT 2021
