Wed Sep 29 14:00:44 CDT 2021
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
$DATA ../../../../data/spa/S1/dat6.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1696.46209520445        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3089E+02 -5.6200E+01 -8.8576E+01  4.4751E+01  1.4308E+02  5.1799E+01  3.9734E+00  1.3778E+01  7.5942E+00  5.9885E+00
             3.7253E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1707.44738797750        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.9747E-01  1.0815E+00  1.1940E+00  9.8360E-01  9.9370E-01  1.0006E+00  9.5517E-01  9.1367E-01  1.0126E+00  8.9428E-01
             8.3166E-01
 PARAMETER:  9.7465E-02  1.7834E-01  2.7734E-01  8.3461E-02  9.3679E-02  1.0057E-01  5.4138E-02  9.7193E-03  1.1248E-01 -1.1734E-02
            -8.4336E-02
 GRADIENT:  -3.0773E+00  3.3491E+00  1.2932E+01 -1.1889E+01  3.5660E+00  6.2282E-01  4.2276E+00  1.8775E+00 -9.9133E+00 -2.1171E+01
            -4.7882E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1711.60728403646        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0048E+00  1.1944E+00  1.1743E+00  9.1628E-01  1.0286E+00  9.9415E-01  5.5308E-01  5.6992E-01  1.1821E+00  1.0408E+00
             8.5567E-01
 PARAMETER:  1.0479E-01  2.7761E-01  2.6067E-01  1.2562E-02  1.2825E-01  9.4138E-02 -4.9226E-01 -4.6226E-01  2.6730E-01  1.3997E-01
            -5.5865E-02
 GRADIENT:   1.1837E+01  2.9930E+01  2.4364E+01  1.0419E+01 -2.2856E+01 -2.1244E+00 -1.6042E+00 -3.3507E+00 -1.1669E+01 -1.4148E+01
            -3.4561E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1714.33212765370        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.9852E-01  1.3372E+00  1.0136E+00  8.1660E-01  1.0528E+00  9.9839E-01  5.6045E-01  4.9428E-01  1.3348E+00  1.0698E+00
             9.2946E-01
 PARAMETER:  9.8523E-02  3.9056E-01  1.1349E-01 -1.0260E-01  1.5144E-01  9.8385E-02 -4.7902E-01 -6.0466E-01  3.8877E-01  1.6745E-01
             2.6846E-02
 GRADIENT:  -5.7994E+00  1.5858E-01 -1.2110E+00  8.5330E+00  3.4376E+00 -4.9583E-01 -2.2966E-01 -8.8174E-02  7.7142E-01 -1.5426E+00
             7.4755E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1715.06171774096        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  1.0044E+00  1.6116E+00  8.3203E-01  6.3703E-01  1.1096E+00  1.0011E+00  5.9384E-01  4.1997E-01  1.5616E+00  1.0953E+00
             9.0142E-01
 PARAMETER:  1.0435E-01  5.7725E-01 -8.3887E-02 -3.5094E-01  2.0401E-01  1.0112E-01 -4.2114E-01 -7.6757E-01  5.4571E-01  1.9105E-01
            -3.7854E-03
 GRADIENT:   6.0010E+00  1.0643E+01  5.5262E-01  8.4325E+00 -6.4507E+00  2.9954E-01 -2.3763E-01  5.1507E-01 -1.4562E+00  3.6126E+00
            -4.6167E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1715.09052096154        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0040E+00  1.6938E+00  7.7741E-01  5.8359E-01  1.1297E+00  1.0009E+00  5.9506E-01  3.9293E-01  1.6553E+00  1.0934E+00
             9.0419E-01
 PARAMETER:  1.0397E-01  6.2695E-01 -1.5179E-01 -4.3855E-01  2.2196E-01  1.0092E-01 -4.1909E-01 -8.3412E-01  6.0396E-01  1.8929E-01
            -7.1513E-04
 GRADIENT:   4.7137E+00  1.4463E+01  6.8370E-01  8.9026E+00 -6.6225E+00  2.0833E-01 -3.0737E-01  5.3776E-01 -1.5112E+00  2.7772E+00
            -3.5213E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1715.09429644782        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  1.0038E+00  1.7203E+00  7.5884E-01  5.6528E-01  1.1367E+00  1.0008E+00  5.9531E-01  3.8009E-01  1.6904E+00  1.0930E+00
             9.0531E-01
 PARAMETER:  1.0380E-01  6.4251E-01 -1.7597E-01 -4.7043E-01  2.2814E-01  1.0085E-01 -4.1867E-01 -8.6734E-01  6.2499E-01  1.8895E-01
             5.1715E-04
 GRADIENT:   4.2167E+00  1.3607E+01  5.4337E-01  8.2369E+00 -6.1609E+00  1.7828E-01 -2.2449E-01  5.2940E-01 -1.3695E+00  2.5243E+00
            -3.0534E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1715.34799441115        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1240             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0033E+00  1.7183E+00  7.5216E-01  5.5331E-01  1.1424E+00  1.0009E+00  5.9656E-01  2.7665E-01  1.7101E+00  1.0834E+00
             9.1021E-01
 PARAMETER:  1.0334E-01  6.4133E-01 -1.8481E-01 -4.9184E-01  2.3314E-01  1.0088E-01 -4.1658E-01 -1.1850E+00  6.3653E-01  1.8007E-01
             5.9243E-03
 GRADIENT:   5.5950E+02  8.6424E+02  9.3092E-01  1.1914E+02  2.1604E+01  6.4815E+01  1.7902E+01  3.2390E-01  3.8486E+01  2.5691E+00
             2.4283E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1715.44264325690        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  1.0022E+00  1.7425E+00  7.2565E-01  5.4020E-01  1.1426E+00  9.9997E-01  5.9948E-01  9.2092E-02  1.7358E+00  1.0774E+00
             9.1170E-01
 PARAMETER:  1.0215E-01  6.5532E-01 -2.2069E-01 -5.1582E-01  2.3328E-01  9.9970E-02 -4.1169E-01 -2.2850E+00  6.5144E-01  1.7452E-01
             7.5511E-03
 GRADIENT:   4.0701E-01 -6.7975E+00  5.8704E-02  6.8533E-01  2.1513E-01 -1.4467E-01  2.1941E-02  3.1783E-02 -2.7392E-01  2.6778E-01
             5.2942E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1715.46567788438        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1607
 NPARAMETR:  1.0030E+00  1.7393E+00  7.2496E-01  5.3986E-01  1.1417E+00  1.0009E+00  5.9912E-01  1.2318E-02  1.7400E+00  1.0756E+00
             9.1157E-01
 PARAMETER:  1.0296E-01  6.5350E-01 -2.2164E-01 -5.1645E-01  2.3253E-01  1.0088E-01 -4.1229E-01 -4.2967E+00  6.5389E-01  1.7291E-01
             7.4171E-03
 GRADIENT:   2.2888E+00 -1.2083E+01 -2.6644E-01 -8.7928E-01  1.0715E+00  2.1939E-01  1.6609E-01  5.8884E-04  3.9960E-01  1.2475E-01
             9.1415E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1715.46810242221        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1794
 NPARAMETR:  1.0027E+00  1.7386E+00  7.2602E-01  5.4154E-01  1.1405E+00  1.0007E+00  5.9869E-01  1.0000E-02  1.7368E+00  1.0749E+00
             9.1158E-01
 PARAMETER:  1.0270E-01  6.5308E-01 -2.2017E-01 -5.1334E-01  2.3149E-01  1.0074E-01 -4.1301E-01 -5.2801E+00  6.5204E-01  1.7219E-01
             7.4211E-03
 GRADIENT:   1.6787E+00 -9.0565E+00  1.0867E-01 -4.6105E-02 -2.9876E-03  1.6185E-01  2.4690E-02  0.0000E+00  2.5217E-01  7.9527E-03
             1.5012E-02

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1715.46810242221        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1822
 NPARAMETR:  1.0037E+00  1.7309E+00  7.2454E-01  5.4074E-01  1.1404E+00  1.0013E+00  6.0073E-01  1.0000E-02  1.7355E+00  1.0748E+00
             9.1145E-01
 PARAMETER:  1.0270E-01  6.5308E-01 -2.2017E-01 -5.1334E-01  2.3149E-01  1.0074E-01 -4.1301E-01 -5.2801E+00  6.5204E-01  1.7219E-01
             7.4211E-03
 GRADIENT:  -5.8844E-01  3.2506E+00  1.1050E-01  1.9195E-01  5.6148E-02 -5.6889E-02 -6.6701E-02  0.0000E+00  3.9456E-02  3.3309E-03
             1.5242E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1822
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8673E-04 -4.6302E-02 -2.2983E-04  2.6237E-02 -3.3813E-02
 SE:             2.9872E-02  1.8946E-02  9.7146E-05  2.4574E-02  2.4668E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8967E-01  1.4528E-02  1.7989E-02  2.8568E-01  1.7045E-01

 ETASHRINKSD(%)  1.0000E-10  3.6529E+01  9.9675E+01  1.7673E+01  1.7360E+01
 ETASHRINKVR(%)  1.0000E-10  5.9715E+01  9.9999E+01  3.2223E+01  3.1706E+01
 EBVSHRINKSD(%)  3.6880E-01  3.5488E+01  9.9709E+01  1.8321E+01  1.5238E+01
 EBVSHRINKVR(%)  7.3623E-01  5.8383E+01  9.9999E+01  3.3285E+01  2.8154E+01
 RELATIVEINF(%)  9.9219E+01  3.7058E+00  2.1151E-04  7.3141E+00  2.2971E+01
 EPSSHRINKSD(%)  4.4456E+01
 EPSSHRINKVR(%)  6.9149E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1715.4681024222093     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -980.31727585847113     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.71
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1715.468       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.74E+00  7.26E-01  5.42E-01  1.14E+00  1.00E+00  5.99E-01  1.00E-02  1.74E+00  1.07E+00  9.12E-01
 


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
+        1.10E+03
 
 TH 2
+       -6.25E+00  4.66E+02
 
 TH 3
+        4.64E+00  1.34E+02  1.97E+02
 
 TH 4
+       -1.03E+01  4.11E+02 -9.59E+01  8.49E+02
 
 TH 5
+       -2.80E+00 -2.05E+02 -2.02E+02  1.09E+02  5.30E+02
 
 TH 6
+        8.17E-01 -7.80E-01  8.02E-01 -2.84E+00 -6.88E-01  1.96E+02
 
 TH 7
+        1.72E+00 -3.12E+01  1.60E+01 -2.68E+01 -2.77E+01  7.22E-02  1.05E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.12E+00 -2.60E+01 -1.50E+01  4.57E+01  2.60E+00 -3.07E-01  2.59E+01  0.00E+00  3.31E+01
 
 TH10
+       -2.23E-01 -1.64E+01 -3.36E+01 -2.23E+00 -4.89E+01  1.57E-02  2.30E+01  0.00E+00  1.96E+00  8.91E+01
 
 TH11
+       -7.89E+00 -2.33E+01 -2.16E+01 -7.41E-01  6.43E+00  2.50E+00  9.52E+00  0.00E+00  3.87E+00  1.69E+01  2.60E+02
 
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
 #CPUT: Total CPU Time in Seconds,       29.859
Stop Time:
Wed Sep 29 14:01:16 CDT 2021
