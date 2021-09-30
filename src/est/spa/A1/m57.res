Wed Sep 29 12:13:28 CDT 2021
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
$DATA ../../../../data/spa/A1/dat57.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1076.47895276326        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6398E+02 -2.2707E+01 -2.9660E+01  1.9364E+01  1.0497E+02  6.9548E+01 -9.9524E+00  1.1236E+01  1.1980E+00 -2.0840E+01
            -1.2082E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1478.89239959542        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0901E+00  1.0469E+00  1.2476E+00  1.0220E+00  1.0547E+00  8.7466E-01  8.9749E-01  8.3583E-01  8.9021E-01  7.0003E-01
             2.2355E+00
 PARAMETER:  1.8627E-01  1.4583E-01  3.2120E-01  1.2176E-01  1.5328E-01 -3.3918E-02 -8.1534E-03 -7.9330E-02 -1.6301E-02 -2.5664E-01
             9.0448E-01
 GRADIENT:   2.9452E+02 -1.3772E+01 -4.2341E+00 -1.4397E+01  3.6589E+01 -2.4386E+00  5.5562E-01  1.5762E+00 -6.7496E+00 -4.1014E+00
            -9.8901E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1487.18042222177        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0490E+00  1.0649E+00  7.2325E-01  9.8230E-01  8.0472E-01  8.9248E-01  8.2778E-01  2.6251E-01  9.8570E-01  4.8000E-01
             2.2744E+00
 PARAMETER:  1.4779E-01  1.6285E-01 -2.2401E-01  8.2140E-02 -1.1726E-01 -1.3750E-02 -8.9004E-02 -1.2375E+00  8.5598E-02 -6.3397E-01
             9.2171E-01
 GRADIENT:   1.3455E+02 -5.0246E-02 -2.2092E-01 -1.2124E+00  9.0397E+00  9.1632E+00 -6.1271E+00  4.8961E-01  3.3260E+00 -8.7665E-01
            -7.2782E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1490.47634982020        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0101E+00  9.2621E-01  5.0982E-01  1.0315E+00  6.2049E-01  8.8352E-01  1.0827E+00  1.6208E-01  8.4857E-01  2.3757E-01
             2.4579E+00
 PARAMETER:  1.1010E-01  2.3341E-02 -5.7369E-01  1.3104E-01 -3.7725E-01 -2.3846E-02  1.7943E-01 -1.7197E+00 -6.4205E-02 -1.3373E+00
             9.9931E-01
 GRADIENT:  -2.7134E+01  2.1231E+00 -1.0103E+01  1.7514E+01  2.1634E+01  1.2655E+00  3.6709E+00  4.6351E-01  4.7211E+00  1.7357E+00
            -1.5342E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1494.56326379787        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      407
 NPARAMETR:  1.0536E+00  8.3993E-01  7.1897E-01  1.1237E+00  7.0995E-01  8.6193E-01  1.1854E+00  1.3416E-01  7.7465E-01  2.2013E-01
             2.6222E+00
 PARAMETER:  1.5222E-01 -7.4435E-02 -2.2994E-01  2.1661E-01 -2.4256E-01 -4.8576E-02  2.7007E-01 -1.9087E+00 -1.5534E-01 -1.4135E+00
             1.0640E+00
 GRADIENT:   1.6447E+01  2.5605E+00  6.4100E-01 -2.1337E+00 -6.0288E+00 -6.1714E+00  3.2017E-01  1.7382E-01 -5.4998E-01  8.8956E-01
             5.2579E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1495.17522683103        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  1.0428E+00  7.4740E-01  9.2046E-01  1.2016E+00  7.8729E-01  8.7258E-01  1.2475E+00  9.0317E-02  7.5381E-01  1.3088E-01
             2.6771E+00
 PARAMETER:  1.4194E-01 -1.9116E-01  1.7118E-02  2.8368E-01 -1.3916E-01 -3.6303E-02  3.2113E-01 -2.3044E+00 -1.8261E-01 -1.9335E+00
             1.0847E+00
 GRADIENT:  -3.5504E+00  4.2675E-01 -4.4857E-01  2.3305E+00  4.8430E-03  5.8776E-01 -5.3482E-02  4.5785E-02 -3.6680E-01  1.8495E-01
             1.2531E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1495.36888667692        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  1.0398E+00  6.6333E-01  1.1971E+00  1.2722E+00  8.7994E-01  8.6851E-01  1.0049E+00  3.5705E-02  7.9373E-01  6.1427E-02
             2.7248E+00
 PARAMETER:  1.3899E-01 -3.1048E-01  2.7987E-01  3.4073E-01 -2.7899E-02 -4.0976E-02  1.0484E-01 -3.2325E+00 -1.3101E-01 -2.6899E+00
             1.1024E+00
 GRADIENT:  -2.8610E+00  5.7739E-01  9.4681E-01  5.4443E-01 -1.7592E+00  3.5482E-01  9.7882E-01  6.2856E-03  1.2298E+00  2.7917E-02
            -6.4614E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1495.74689680363        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  1.0395E+00  5.8122E-01  1.3656E+00  1.3236E+00  9.1999E-01  8.6905E-01  4.0790E-01  1.0000E-02  8.0755E-01  1.3926E-02
             2.7345E+00
 PARAMETER:  1.3878E-01 -4.4262E-01  4.1161E-01  3.8034E-01  1.6605E-02 -4.0351E-02 -7.9673E-01 -5.5057E+00 -1.1375E-01 -4.1740E+00
             1.1060E+00
 GRADIENT:   1.8130E+00 -3.3334E+00 -1.3289E+00 -8.1574E+00  3.4773E+00  4.7829E-01  1.8131E-01  0.0000E+00 -3.7449E-01  1.3338E-03
            -2.6981E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1495.80129885459        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1114
 NPARAMETR:  1.0393E+00  5.5968E-01  1.3616E+00  1.3392E+00  9.1351E-01  8.6836E-01  1.5520E-01  1.0000E-02  8.0741E-01  1.0000E-02
             2.7415E+00
 PARAMETER:  1.3858E-01 -4.8040E-01  4.0869E-01  3.9209E-01  9.5379E-03 -4.1151E-02 -1.7630E+00 -7.5529E+00 -1.1392E-01 -5.4320E+00
             1.1085E+00
 GRADIENT:   1.0165E+00 -1.7412E+00 -2.0035E+00 -2.2933E+00  3.3390E+00  3.9983E-03  2.2218E-02  0.0000E+00 -5.9070E-01  0.0000E+00
            -9.8985E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1495.99189048301        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1292
 NPARAMETR:  1.0409E+00  7.4174E-01  1.4030E+00  1.2290E+00  9.8317E-01  8.7021E-01  1.0000E-02  1.0000E-02  8.7952E-01  1.0000E-02
             2.7496E+00
 PARAMETER:  1.4009E-01 -1.9876E-01  4.3858E-01  3.0620E-01  8.3028E-02 -3.9022E-02 -6.4452E+00 -1.5235E+01 -2.8382E-02 -9.7892E+00
             1.1115E+00
 GRADIENT:  -2.0907E-01 -9.8702E-01 -5.1384E-01 -3.5678E-01  1.6906E+00  9.4596E-02  0.0000E+00  0.0000E+00  5.7017E-01  0.0000E+00
            -4.0494E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1496.01020681872        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1467
 NPARAMETR:  1.0417E+00  8.0616E-01  1.3940E+00  1.1883E+00  1.0007E+00  8.7069E-01  1.0000E-02  1.0000E-02  9.0498E-01  1.0000E-02
             2.7547E+00
 PARAMETER:  1.4084E-01 -1.1547E-01  4.3216E-01  2.7256E-01  1.0072E-01 -3.8473E-02 -7.6930E+00 -1.7277E+01  1.5618E-04 -1.0950E+01
             1.1133E+00
 GRADIENT:  -7.7714E-02  1.6877E-01 -2.6045E-02  2.9600E-01 -6.1267E-02 -2.1993E-02  0.0000E+00  0.0000E+00 -4.4721E-02  0.0000E+00
             8.8950E-02

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1496.01020681872        NO. OF FUNC. EVALS.:  32
 CUMULATIVE NO. OF FUNC. EVALS.:     1499
 NPARAMETR:  1.0424E+00  8.0523E-01  1.3946E+00  1.1873E+00  1.0009E+00  8.7099E-01  1.0000E-02  1.0000E-02  9.0583E-01  1.0000E-02
             2.7531E+00
 PARAMETER:  1.4084E-01 -1.1547E-01  4.3216E-01  2.7256E-01  1.0072E-01 -3.8473E-02 -7.6930E+00 -1.7277E+01  1.5618E-04 -1.0950E+01
             1.1133E+00
 GRADIENT:  -4.8711E-01  1.7072E-01 -2.5819E-02  4.4252E-01 -6.0932E-02 -3.0956E-02  0.0000E+00  0.0000E+00 -4.5726E-02  0.0000E+00
             1.2477E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1499
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4708E-03 -3.1674E-04  8.5298E-05 -5.0112E-03  1.0560E-05
 SE:             2.8967E-02  1.5921E-04  9.1326E-05  2.5363E-02  2.1379E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5950E-01  4.6659E-02  3.5030E-01  8.4338E-01  9.6061E-01

 ETASHRINKSD(%)  2.9583E+00  9.9467E+01  9.9694E+01  1.5030E+01  9.9284E+01
 ETASHRINKVR(%)  5.8292E+00  9.9997E+01  9.9999E+01  2.7801E+01  9.9995E+01
 EBVSHRINKSD(%)  3.1929E+00  9.9493E+01  9.9665E+01  1.4535E+01  9.9255E+01
 EBVSHRINKVR(%)  6.2839E+00  9.9997E+01  9.9999E+01  2.6957E+01  9.9994E+01
 RELATIVEINF(%)  9.1810E+01  9.1909E-05  1.3553E-04  3.8756E+00  3.7853E-04
 EPSSHRINKSD(%)  2.2170E+01
 EPSSHRINKVR(%)  3.9425E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1496.0102068187234     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -760.85938025498524     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1496.010       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  8.06E-01  1.39E+00  1.19E+00  1.00E+00  8.71E-01  1.00E-02  1.00E-02  9.05E-01  1.00E-02  2.75E+00
 


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
+        1.27E+03
 
 TH 2
+       -7.36E+01  4.51E+02
 
 TH 3
+        9.45E+00  7.03E+01  5.60E+01
 
 TH 4
+       -9.46E+01  4.97E+02  2.33E+01  7.00E+02
 
 TH 5
+       -6.20E+00 -3.20E+02 -1.71E+02 -1.68E+02  6.04E+02
 
 TH 6
+       -3.95E-02 -1.56E+01  5.14E+00 -1.91E+01 -7.77E+00  2.29E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.02E+01 -5.66E+01  5.41E+00  4.21E+00  2.01E+01  8.42E+00  0.00E+00  0.00E+00  1.17E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.72E+01 -2.00E+01 -9.47E-01 -1.42E+01  6.76E+00  2.81E+00  0.00E+00  0.00E+00  1.48E+01  0.00E+00  5.70E+01
 
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
 #CPUT: Total CPU Time in Seconds,       23.107
Stop Time:
Wed Sep 29 12:13:52 CDT 2021
