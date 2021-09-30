Wed Sep 29 19:37:42 CDT 2021
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
$DATA ../../../../data/spa/D/dat8.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m8.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1286.90065461848        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0859E+02 -1.2281E+01 -3.7043E+00  1.9636E+00  6.2062E+01 -2.1867E+02 -1.6377E+02 -2.1052E+01 -1.9772E+02 -3.7645E+01
            -4.5366E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1416.42458263399        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      170
 NPARAMETR:  1.0463E+00  1.1255E+00  9.5813E-01  9.2780E-01  1.0731E+00  2.8562E+00  2.7334E+00  1.1381E+00  1.9008E+00  1.0903E+00
             8.9910E-01
 PARAMETER:  1.4529E-01  2.1818E-01  5.7230E-02  2.5064E-02  1.7056E-01  1.1495E+00  1.1055E+00  2.2932E-01  7.4228E-01  1.8644E-01
            -6.3586E-03
 GRADIENT:   1.4660E+01 -6.2761E+00 -2.5205E+01  1.7284E+01  2.4770E+00  1.9889E+02  3.6623E+01  4.3897E+00  6.2026E+01  1.8573E+01
            -8.6906E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1425.19903371486        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  1.0745E+00  8.2804E-01  1.6513E+00  1.1066E+00  1.1356E+00  2.4455E+00  3.3692E+00  2.2875E+00  1.9198E+00  1.1939E+00
             9.0145E-01
 PARAMETER:  1.7183E-01 -8.8695E-02  6.0157E-01  2.0125E-01  2.2713E-01  9.9424E-01  1.3147E+00  9.2747E-01  7.5220E-01  2.7721E-01
            -3.7553E-03
 GRADIENT:   3.1402E+01 -4.8048E+00 -1.2196E+01 -1.5022E+01 -3.2187E+01  1.5220E+02  3.6059E+01  2.1931E+01  6.8919E+01  1.9738E+01
            -7.9035E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1442.65993680962        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  1.0775E+00  5.1747E-01  2.2488E+00  1.4858E+00  1.0743E+00  2.2048E+00  4.3134E+00  2.8125E+00  1.3521E+00  1.1815E+00
             9.5191E-01
 PARAMETER:  1.7468E-01 -5.5881E-01  9.1039E-01  4.9597E-01  1.7170E-01  8.9066E-01  1.5617E+00  1.1341E+00  4.0166E-01  2.6679E-01
             5.0716E-02
 GRADIENT:   3.9590E+01  2.1995E+01 -1.5489E+00  4.5312E+01 -7.6895E+01  1.1757E+02  2.0493E+01  3.0807E+01  2.5109E+01  1.8330E+01
             1.2173E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1475.65186907524        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  1.0117E+00  5.9127E-01  1.4942E+00  1.2819E+00  1.0204E+00  1.6434E+00  3.4479E+00  1.8020E+00  1.1421E+00  9.5106E-01
             9.2865E-01
 PARAMETER:  1.1160E-01 -4.2549E-01  5.0162E-01  3.4832E-01  1.2023E-01  5.9676E-01  1.3378E+00  6.8889E-01  2.3289E-01  4.9820E-02
             2.5980E-02
 GRADIENT:   1.0181E+01 -2.6648E-01 -2.2772E+00  1.1191E+01 -1.0443E+01  1.9139E+01 -2.1524E-02  8.5591E+00  6.8369E+00  2.1496E+00
             6.8837E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1475.92726889286        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      893             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0023E+00  5.9567E-01  1.5216E+00  1.2600E+00  1.0281E+00  1.6381E+00  3.4232E+00  1.8083E+00  1.1408E+00  9.4036E-01
             9.2814E-01
 PARAMETER:  1.0232E-01 -4.1806E-01  5.1976E-01  3.3111E-01  1.2769E-01  5.9356E-01  1.3306E+00  6.9237E-01  2.3171E-01  3.8512E-02
             2.5425E-02
 GRADIENT:   5.0858E+02  1.0279E+02  1.3044E+01  4.7165E+02 -3.3618E+00  7.8306E+02  5.1055E+02  1.2232E+01  4.1447E+01  1.1126E+00
             9.3267E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1476.05157145872        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  9.9933E-01  6.3330E-01  1.5273E+00  1.2477E+00  1.0421E+00  1.5969E+00  3.4011E+00  1.8122E+00  1.1399E+00  9.4997E-01
             9.2827E-01
 PARAMETER:  9.9328E-02 -3.5680E-01  5.2352E-01  3.2130E-01  1.4121E-01  5.6808E-01  1.3241E+00  6.9455E-01  2.3097E-01  4.8673E-02
             2.5572E-02
 GRADIENT:  -8.0728E-01  1.9519E-01  2.4692E+00  1.5189E+00 -8.6316E+00  6.6683E+00  3.4442E+00  5.6836E+00  7.1148E+00 -1.0783E-01
            -5.9874E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1477.54152032536        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.0023E+00  9.2824E-01  1.3976E+00  1.0653E+00  1.1169E+00  1.5880E+00  2.6828E+00  1.8397E+00  1.1341E+00  9.9039E-01
             9.2783E-01
 PARAMETER:  1.0230E-01  2.5532E-02  4.3473E-01  1.6328E-01  2.1052E-01  5.6247E-01  1.0869E+00  7.0963E-01  2.2586E-01  9.0348E-02
             2.5092E-02
 GRADIENT:  -7.6332E-01  2.0246E-01  2.2181E+00  1.6623E+00  3.6141E+00  4.0925E+00 -1.5242E+00  5.4585E+00  2.3074E+00 -1.9256E+00
            -8.3781E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1477.83013854259        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1394
 NPARAMETR:  1.0115E+00  9.2344E-01  1.3723E+00  1.0595E+00  1.1098E+00  1.6509E+00  2.6457E+00  1.6278E+00  1.0777E+00  1.0182E+00
             9.3047E-01
 PARAMETER:  1.1146E-01  2.0352E-02  4.1647E-01  1.5779E-01  2.0415E-01  6.0132E-01  1.0729E+00  5.8721E-01  1.7484E-01  1.1801E-01
             2.7930E-02
 GRADIENT:   7.3018E+00 -4.7529E-01  1.4384E+01 -4.6520E+00 -2.2679E+00  2.1722E+01 -7.7269E+00 -4.6006E+00 -5.1655E+00 -2.9753E-01
            -2.6656E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1478.72115487919        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1571
 NPARAMETR:  1.0029E+00  9.2145E-01  1.1662E+00  1.0581E+00  1.0496E+00  1.5715E+00  2.7844E+00  1.4648E+00  1.0979E+00  9.5163E-01
             9.2560E-01
 PARAMETER:  1.0288E-01  1.8197E-02  2.5377E-01  1.5652E-01  1.4839E-01  5.5200E-01  1.1240E+00  4.8174E-01  1.9344E-01  5.0423E-02
             2.2690E-02
 GRADIENT:  -6.2115E-01  1.4983E+00 -1.5249E+00  4.1967E+00  1.0681E+01 -1.4361E+00  4.7958E+00  2.3015E+00  1.4119E+00 -1.4460E-01
            -1.0168E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1479.49112176204        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1748
 NPARAMETR:  9.9908E-01  9.2706E-01  9.2695E-01  1.0309E+00  9.4443E-01  1.5456E+00  2.6955E+00  1.1207E+00  1.0791E+00  8.4105E-01
             9.2721E-01
 PARAMETER:  9.9077E-02  2.4263E-02  2.4144E-02  1.3045E-01  4.2823E-02  5.3542E-01  1.0916E+00  2.1392E-01  1.7610E-01 -7.3110E-02
             2.4424E-02
 GRADIENT:  -5.5009E+00 -1.8691E+00 -6.6728E-01  1.9373E+00  1.8632E+00 -1.0028E+01  1.4978E+00  1.3287E+00  3.4307E-01  3.2663E-01
            -1.4196E-02

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1479.49112176204        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1776
 NPARAMETR:  9.9908E-01  9.2706E-01  9.2695E-01  1.0309E+00  9.4443E-01  1.5456E+00  2.6955E+00  1.1207E+00  1.0791E+00  8.4105E-01
             9.2721E-01
 PARAMETER:  9.9077E-02  2.4263E-02  2.4144E-02  1.3045E-01  4.2823E-02  5.3542E-01  1.0916E+00  2.1392E-01  1.7610E-01 -7.3110E-02
             2.4424E-02
 GRADIENT:   8.6200E+04  8.6205E+04  8.6240E+04 -6.6098E+04 -8.6205E+04  1.6077E+04 -5.2529E-01 -4.0422E+04 -4.8980E+04  4.3104E+04
             8.6187E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1776
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.7376E-03  1.8868E-02 -4.7057E-02 -2.2672E-02 -2.3503E-02
 SE:             3.1776E-02  2.4516E-02  1.4999E-02  2.1759E-02  1.8190E-02
 N:                     100         100         100         100         100

 P VAL.:         8.5671E-01  4.4153E-01  1.7045E-03  2.9743E-01  1.9635E-01

 ETASHRINKSD(%)  1.0000E-10  1.7869E+01  4.9753E+01  2.7104E+01  3.9060E+01
 ETASHRINKVR(%)  1.0000E-10  3.2545E+01  7.4753E+01  4.6862E+01  6.2863E+01
 EBVSHRINKSD(%)  1.6755E-01  1.4780E+01  5.3707E+01  2.8332E+01  3.6306E+01
 EBVSHRINKVR(%)  3.3483E-01  2.7375E+01  7.8569E+01  4.8636E+01  5.9431E+01
 RELATIVEINF(%)  9.9526E+01  1.8967E+01  4.3666E+00  1.0762E+01  9.0437E+00
 EPSSHRINKSD(%)  4.6494E+01
 EPSSHRINKVR(%)  7.1371E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1479.4911217620427     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -744.34029519830449     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.40
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.07
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1479.491       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  9.27E-01  9.27E-01  1.03E+00  9.44E-01  1.55E+00  2.70E+00  1.12E+00  1.08E+00  8.41E-01  9.27E-01
 


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
+        2.16E+07
 
 TH 2
+        2.33E+07  5.02E+07
 
 TH 3
+        2.33E+07 -3.08E+02  2.51E+07
 
 TH 4
+       -1.60E+07 -1.73E+07  1.73E+07  2.38E+07
 
 TH 5
+        2.17E+01  2.55E+02 -2.46E+07 -1.70E+07  2.42E+07
 
 TH 6
+        2.61E+06 -4.00E+01  2.81E+06  1.94E+06  1.99E+00  3.15E+05
 
 TH 7
+        4.98E-01  1.41E+01 -7.13E+00 -1.57E+01  8.06E+00 -1.42E-01  2.49E+04
 
 TH 8
+       -6.88E+01  1.32E+02 -7.73E+03  3.10E+03 -9.78E+00 -1.10E+02  1.29E+00  3.77E+06
 
 TH 9
+       -1.13E+07 -1.22E+07 -9.70E+03  8.43E+06 -1.20E+07 -1.38E+02  6.95E+00 -4.71E+06  1.19E+07
 
 TH10
+        2.56E+07 -3.94E+02  2.77E+07  1.91E+07 -6.40E+01  3.10E+06 -7.94E-01 -4.76E+01  1.35E+07  3.05E+07
 
 TH11
+        1.74E+02 -3.59E+02  1.98E+04 -7.98E+03  1.04E+01  2.85E+02  9.46E-01 -9.72E+06  1.22E+07  1.69E+02  2.51E+07
 
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
 #CPUT: Total CPU Time in Seconds,       34.518
Stop Time:
Wed Sep 29 19:38:18 CDT 2021
