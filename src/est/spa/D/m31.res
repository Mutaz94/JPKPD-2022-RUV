Wed Sep 29 19:57:01 CDT 2021
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
$DATA ../../../../data/spa/D/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12309.7508158601        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5085E+02  1.7221E+02 -4.7855E+01  6.5654E+01  2.3313E+02 -1.3580E+03 -6.6805E+02 -2.8270E+01 -1.1042E+03 -5.1545E+02
            -2.3978E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -646.587147331626        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2809E+00  1.1526E+00  1.0307E+00  1.5489E+00  1.1735E+00  1.5385E+00  1.3161E+00  9.6287E-01  1.3493E+00  1.1173E+00
             1.4604E+01
 PARAMETER:  3.4754E-01  2.4203E-01  1.3026E-01  5.3754E-01  2.5995E-01  5.3080E-01  3.7469E-01  6.2160E-02  3.9961E-01  2.1094E-01
             2.7813E+00
 GRADIENT:  -4.8160E+01  8.3603E+00 -7.7934E+00  1.9701E+01 -6.9406E+00  2.5671E+01  3.2367E+00  3.4702E+00  1.8349E+01  4.3632E+00
             1.9756E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -659.848874704758        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2752E+00  1.2966E+00  4.2567E+00  1.7409E+00  4.9946E+00  1.6163E+00  2.2318E+00  3.7210E-01  1.7927E+00  8.3513E+00
             1.2469E+01
 PARAMETER:  3.4309E-01  3.5972E-01  1.5485E+00  6.5439E-01  1.7084E+00  5.8012E-01  9.0279E-01 -8.8860E-01  6.8372E-01  2.2224E+00
             2.6233E+00
 GRADIENT:  -2.0469E+01  2.1013E+01  2.0201E+00  4.7269E+01 -2.8781E+00  1.2340E+01  8.6239E+00 -3.5187E-04  2.4662E+01  3.4616E+00
             1.2003E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -682.076880392949        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1732E+00  1.1384E+00  1.4757E+00  1.1590E+00  4.6754E+00  1.3071E+00  1.3548E+00  5.0537E-01  9.9202E-01  9.1118E+00
             1.1393E+01
 PARAMETER:  2.5976E-01  2.2962E-01  4.8912E-01  2.4757E-01  1.6423E+00  3.6779E-01  4.0364E-01 -5.8247E-01  9.1993E-02  2.3096E+00
             2.5330E+00
 GRADIENT:   1.9437E+01 -1.5881E+01  4.4999E+00 -3.8760E+01 -9.3722E+00  7.5723E-01  5.2402E+00  2.5678E-02  7.5745E+00  1.2542E+01
             6.0272E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -693.621777983930        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.1109E+00  7.3717E-01  1.0079E+00  1.3547E+00  6.6087E+00  1.2781E+00  1.6959E+00  2.5000E-01  7.2156E-01  8.0202E+00
             1.0442E+01
 PARAMETER:  2.0515E-01 -2.0494E-01  1.0789E-01  4.0356E-01  1.9884E+00  3.4535E-01  6.2819E-01 -1.2863E+00 -2.2634E-01  2.1820E+00
             2.4458E+00
 GRADIENT:   8.0460E-01  9.4272E+00  7.3725E+00 -3.9112E+00 -1.9850E+00 -7.1902E+00  9.1145E-01  3.9559E-02  1.7744E+00 -3.0404E+00
            -1.5411E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -739.307266510476        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  8.5271E-01  9.7598E-02  1.0218E-01  8.6464E-01  2.6661E+01  1.5935E+00  2.0643E-01  1.0000E-02  4.3615E-02  5.3641E+00
             1.1795E+01
 PARAMETER: -5.9333E-02 -2.2269E+00 -2.1811E+00 -4.5447E-02  3.3832E+00  5.6594E-01 -1.4778E+00 -1.5136E+01 -3.0324E+00  1.7797E+00
             2.5677E+00
 GRADIENT:   6.2778E+01  4.8130E+01 -4.5062E+01  2.0022E+01 -2.3105E+00  3.8842E+01  9.5034E-01  0.0000E+00  4.3921E-02  1.6548E+00
             6.9407E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -763.383625063035        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  5.8556E-01  3.2596E-02  4.6526E-02  4.5459E-01  7.9298E+01  1.3279E+00  3.7726E-02  1.0000E-02  1.0000E-02  8.0213E+00
             1.0393E+01
 PARAMETER: -4.3518E-01 -3.3236E+00 -2.9678E+00 -6.8836E-01  4.4732E+00  3.8360E-01 -3.1774E+00 -1.9725E+01 -4.7952E+00  2.1821E+00
             2.4412E+00
 GRADIENT:   3.4291E+01  3.0849E+00  4.7577E+01 -8.1609E+01 -4.5757E-02 -9.0853E-01  1.2445E-03  0.0000E+00  0.0000E+00  5.1465E-03
            -4.9930E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -767.910390491480        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  4.6020E-01  1.4502E-02  2.6531E-02  3.1188E-01  1.4814E+02  1.3166E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0165E+01
             1.0344E+01
 PARAMETER: -6.7609E-01 -4.1335E+00 -3.5294E+00 -1.0651E+00  5.0981E+00  3.7508E-01 -4.6751E+00 -2.4032E+01 -5.9293E+00  2.4190E+00
             2.4364E+00
 GRADIENT:   4.1250E+01  2.3927E-01  5.9103E+01  1.6622E+01 -2.5227E-04  3.3485E+00  0.0000E+00  0.0000E+00  0.0000E+00 -1.0401E-04
             2.0017E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -767.941634031351        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      808
 NPARAMETR:  4.5916E-01  1.3249E-02  2.6358E-02  3.1139E-01  1.4797E+02  1.3162E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0168E+01
             1.0329E+01
 PARAMETER: -6.7835E-01 -4.2238E+00 -3.5360E+00 -1.0667E+00  5.0970E+00  3.7475E-01 -4.7523E+00 -2.4032E+01 -5.9293E+00  2.4192E+00
             2.4349E+00
 GRADIENT:  -3.2740E+00  9.0790E-02 -4.4184E-01  2.4830E+00  1.1449E-03 -8.6133E-01  0.0000E+00  0.0000E+00  0.0000E+00 -2.6721E-04
            -1.4821E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -768.064261171062        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      970             RESET HESSIAN, TYPE I
 NPARAMETR:  4.5371E-01  1.0000E-02  2.5354E-02  3.0204E-01  2.8362E+01  1.3212E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6944E+01
             1.0349E+01
 PARAMETER: -6.9030E-01 -7.0417E+00 -3.5748E+00 -1.0972E+00  3.4451E+00  3.7851E-01 -4.7523E+00 -2.4032E+01 -5.9293E+00  2.9299E+00
             2.4369E+00
 GRADIENT:   4.7069E+01  0.0000E+00  5.7337E+01  1.9659E+01  3.4569E-02  4.6232E+00  0.0000E+00  0.0000E+00  0.0000E+00  6.2495E-03
             1.9627E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -768.069161941943        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1138
 NPARAMETR:  4.5208E-01  1.0000E-02  2.5254E-02  3.0184E-01  2.8031E+01  1.3171E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7148E+01
             1.0329E+01
 PARAMETER: -6.9390E-01 -7.0417E+00 -3.5788E+00 -1.0979E+00  3.4333E+00  3.7545E-01 -4.7523E+00 -2.4032E+01 -5.9293E+00  2.9419E+00
             2.4349E+00
 GRADIENT:  -1.8557E+00  0.0000E+00 -1.9000E+00  3.0389E+00  1.8721E-02 -3.5407E-01  0.0000E+00  0.0000E+00  0.0000E+00 -7.1946E-03
            -1.4201E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -768.086151168575        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1286
 NPARAMETR:  4.5090E-01  1.0000E-02  2.5020E-02  2.9998E-01  2.9361E+01  1.3178E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.8165E+01
             1.0326E+01
 PARAMETER: -6.9650E-01 -7.0417E+00 -3.5881E+00 -1.1040E+00  3.4797E+00  3.7597E-01 -4.7523E+00 -2.4032E+01 -5.9293E+00  2.9995E+00
             2.4346E+00
 GRADIENT:   4.6321E+01  0.0000E+00  5.2697E+01  2.7295E+01  3.6116E-02  3.6534E+00  0.0000E+00  0.0000E+00  0.0000E+00  1.1502E-02
             1.7562E+01

0ITERATION NO.:   58    OBJECTIVE VALUE:  -768.086479493107        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:     1390
 NPARAMETR:  4.5266E-01  1.0000E-02  2.5478E-02  2.9996E-01  2.9769E+01  1.3166E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7112E+01
             1.0333E+01
 PARAMETER: -6.9731E-01 -7.0417E+00 -3.5868E+00 -1.1067E+00  3.4589E+00  3.7587E-01 -4.7523E+00 -2.4032E+01 -5.9293E+00  2.9695E+00
             2.4361E+00
 GRADIENT:  -2.7085E+00  0.0000E+00 -2.6214E+01 -2.9102E+00 -3.4909E+01  1.3060E-01  0.0000E+00  0.0000E+00  0.0000E+00  4.1587E+01
             3.4183E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1390
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.6753E-03  2.0256E-06  8.1334E-05 -1.9581E-04 -7.0790E-04
 SE:             2.8574E-02  1.2711E-06  3.1360E-04  3.9105E-04  7.8209E-04
 N:                     100         100         100         100         100

 P VAL.:         8.9765E-01  1.1104E-01  7.9536E-01  6.1655E-01  3.6539E-01

 ETASHRINKSD(%)  4.2729E+00  9.9996E+01  9.8949E+01  9.8690E+01  9.7380E+01
 ETASHRINKVR(%)  8.3633E+00  1.0000E+02  9.9989E+01  9.9983E+01  9.9931E+01
 EBVSHRINKSD(%)  4.4954E+00  9.9995E+01  9.8982E+01  9.8727E+01  9.7612E+01
 EBVSHRINKVR(%)  8.7886E+00  1.0000E+02  9.9990E+01  9.9984E+01  9.9943E+01
 RELATIVEINF(%)  3.7021E-01  9.9180E-09  1.6302E-05  1.8047E-05  5.1952E-04
 EPSSHRINKSD(%)  7.2344E+00
 EPSSHRINKVR(%)  1.3945E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -768.08647949310716     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -32.935652929368985     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -768.086       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.51E-01  1.00E-02  2.51E-02  2.99E-01  2.88E+01  1.32E+00  1.00E-02  1.00E-02  1.00E-02  1.76E+01  1.03E+01
 


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
+        2.84E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.13E+05  0.00E+00  1.40E+06
 
 TH 4
+       -5.66E+02  0.00E+00 -2.00E+05  1.25E+04
 
 TH 5
+        1.25E+00  0.00E+00  3.14E+02  5.63E+00  3.05E-01
 
 TH 6
+        5.16E+00  0.00E+00  6.38E+04 -9.37E+01 -2.81E-01  9.41E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.25E+00  0.00E+00 -6.13E+02 -1.06E+01 -5.92E-01  5.36E-01  0.00E+00  0.00E+00  0.00E+00  1.15E+00
 
 TH11
+       -3.95E+02  0.00E+00  1.57E+03 -2.20E+01  1.17E+00  1.30E+00  0.00E+00  0.00E+00  0.00E+00 -2.27E+00  8.59E+00
 
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
 #CPUT: Total CPU Time in Seconds,       24.132
Stop Time:
Wed Sep 29 19:57:26 CDT 2021
