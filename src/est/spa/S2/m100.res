Sat Sep 18 13:47:06 CDT 2021
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
$DATA ../../../../data/spa/S2/dat100.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1665.84482911333        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2728E+02 -5.1401E+01 -3.2004E+01 -8.1587E+00  4.1633E+01 -1.1191E+01 -2.1004E+00 -6.0457E-01  3.1867E+01  3.8692E+00
             1.8121E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1673.41741079403        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6123E-01  1.1938E+00  1.0075E+00  8.7109E-01  1.1117E+00  1.0176E+00  1.1136E+00  1.1042E+00  6.7257E-01  9.9085E-01
             9.7140E-01
 PARAMETER:  6.0461E-02  2.7715E-01  1.0750E-01 -3.8008E-02  2.0588E-01  1.1740E-01  2.0758E-01  1.9911E-01 -2.9666E-01  9.0811E-02
             7.0987E-02
 GRADIENT:   4.1882E+01 -2.9982E+01 -5.7483E+00 -2.2937E+01  4.5600E+01  9.4546E-01 -5.6516E+00 -6.5039E+00 -5.3218E+00 -1.4193E+01
             8.5006E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1674.65994631729        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.7325E-01  1.3640E+00  9.4635E-01  8.1087E-01  1.1336E+00  1.0296E+00  1.1029E+00  1.3467E+00  6.6903E-01  1.0597E+00
             9.5659E-01
 PARAMETER:  7.2887E-02  4.1042E-01  4.4858E-02 -1.0965E-01  2.2537E-01  1.2915E-01  1.9793E-01  3.9767E-01 -3.0193E-01  1.5794E-01
             5.5618E-02
 GRADIENT:   2.0846E+01  3.5575E+01 -8.4498E+00  3.8660E+01  1.5170E+01 -2.2006E+00  6.4164E+00 -2.3339E-01 -4.1635E+00 -2.5900E+00
            -3.2117E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1675.38998226716        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  9.7328E-01  1.3649E+00  9.4809E-01  7.8675E-01  1.1373E+00  1.0465E+00  1.0558E+00  1.3502E+00  6.6916E-01  1.0781E+00
             9.6503E-01
 PARAMETER:  7.2914E-02  4.1111E-01  4.6692E-02 -1.3984E-01  2.2865E-01  1.4541E-01  1.5432E-01  4.0028E-01 -3.0173E-01  1.7517E-01
             6.4402E-02
 GRADIENT:   2.0526E+01  1.0342E+00 -1.2382E+00  9.5665E-01  2.5039E+00  4.3004E+00 -9.5573E-02 -1.4787E+00 -5.1533E+00 -9.9408E-02
             2.8905E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1675.41758804676        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      577
 NPARAMETR:  9.7328E-01  1.3649E+00  9.5273E-01  7.8714E-01  1.1366E+00  1.0349E+00  1.0568E+00  1.3502E+00  6.6916E-01  1.0791E+00
             9.6461E-01
 PARAMETER:  7.2914E-02  4.1111E-01  5.1579E-02 -1.3935E-01  2.2806E-01  1.3435E-01  1.5522E-01  4.0028E-01 -3.0173E-01  1.7608E-01
             6.3967E-02
 GRADIENT:   2.1045E+01  2.2801E+00 -2.2638E-02  8.8464E-02  1.1815E-04  1.6257E-02  3.6605E-02 -1.7186E+00 -5.1906E+00  6.1939E-03
             2.1754E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1676.03278608020        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  9.6349E-01  1.3657E+00  9.6372E-01  7.8714E-01  1.1366E+00  1.0349E+00  1.0268E+00  1.3505E+00  7.7120E-01  1.0789E+00
             9.6027E-01
 PARAMETER:  6.2804E-02  4.1168E-01  6.3044E-02 -1.3935E-01  2.2802E-01  1.3432E-01  1.2648E-01  4.0048E-01 -1.5981E-01  1.7596E-01
             5.9458E-02
 GRADIENT:   4.5714E-01 -1.1054E+00  1.9861E-01 -2.5432E+00 -2.1239E+00  1.7240E-01  2.2760E+00 -8.8473E-01 -1.1220E+00  1.7245E+00
            -2.1087E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1676.12923471993        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      911            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6331E-01  1.3696E+00  9.6625E-01  7.8657E-01  1.1394E+00  1.0346E+00  9.9544E-01  1.3715E+00  8.1470E-01  1.0621E+00
             9.6016E-01
 PARAMETER:  6.2618E-02  4.1449E-01  6.5669E-02 -1.4007E-01  2.3054E-01  1.3400E-01  9.5435E-02  4.1588E-01 -1.0493E-01  1.6026E-01
             5.9342E-02
 GRADIENT:   4.6315E+01  2.5730E+01 -1.2053E+00  8.6752E+00  3.7615E+00  8.2297E+00  7.5769E-01 -2.7350E-01  3.5266E-01  6.8308E-02
            -2.2008E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1676.14365356923        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  9.6261E-01  1.3725E+00  9.7414E-01  7.8491E-01  1.1399E+00  1.0339E+00  9.9213E-01  1.4067E+00  8.1576E-01  1.0602E+00
             9.6131E-01
 PARAMETER:  6.1890E-02  4.1661E-01  7.3798E-02 -1.4218E-01  2.3094E-01  1.3331E-01  9.2096E-02  4.4122E-01 -1.0364E-01  1.5849E-01
             6.0546E-02
 GRADIENT:  -1.4901E+00 -1.1959E+00 -3.4848E-01  3.6674E-01 -1.3741E+00 -2.6865E-01  1.3883E-02  9.7105E-02 -3.6026E-02 -7.1176E-02
             1.4301E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1676.17368580717        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1268
 NPARAMETR:  9.6143E-01  1.3825E+00  1.0162E+00  7.7917E-01  1.1636E+00  1.0328E+00  9.8932E-01  1.4780E+00  8.1516E-01  1.0841E+00
             9.6127E-01
 PARAMETER:  6.0665E-02  4.2387E-01  1.1604E-01 -1.4953E-01  2.5148E-01  1.3224E-01  8.9264E-02  4.9070E-01 -1.0437E-01  1.8079E-01
             6.0496E-02
 GRADIENT:  -3.6516E+00 -1.7842E+00  3.4854E-01 -1.2402E+00 -8.8687E-01 -6.5625E-01  2.0396E-01 -2.0356E-01 -1.2234E-01  5.6995E-02
            -9.6342E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1676.18222839279        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1443
 NPARAMETR:  9.6249E-01  1.3865E+00  1.0462E+00  7.7830E-01  1.1777E+00  1.0339E+00  9.8511E-01  1.5358E+00  8.1807E-01  1.0968E+00
             9.6189E-01
 PARAMETER:  6.1769E-02  4.2677E-01  1.4516E-01 -1.5064E-01  2.6353E-01  1.3333E-01  8.5001E-02  5.2907E-01 -1.0081E-01  1.9241E-01
             6.1148E-02
 GRADIENT:  -1.1762E+00 -9.5698E-01  6.1314E-02  9.8932E-02 -3.2119E-02 -2.0841E-01 -1.5382E-02  5.8722E-02 -2.6503E-02  1.2733E-01
             5.2548E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1676.18311303969        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1629             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6280E-01  1.3868E+00  1.0442E+00  7.7807E-01  1.1770E+00  1.0342E+00  9.8415E-01  1.5325E+00  8.1990E-01  1.0958E+00
             9.6182E-01
 PARAMETER:  6.2088E-02  4.2700E-01  1.4321E-01 -1.5094E-01  2.6300E-01  1.3362E-01  8.4023E-02  5.2687E-01 -9.8575E-02  1.9152E-01
             6.1068E-02
 GRADIENT:   4.5808E+01  2.8087E+01  1.6772E-01  7.8164E+00  1.6872E+00  8.1401E+00  5.6277E-01  1.5738E-01  3.3557E-01  2.5694E-01
             1.0633E-01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1676.18311303969        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1695
 NPARAMETR:  9.6278E-01  1.3867E+00  1.0440E+00  7.7805E-01  1.1770E+00  1.0343E+00  9.8423E-01  1.5326E+00  8.1989E-01  1.0959E+00
             9.6179E-01
 PARAMETER:  6.2088E-02  4.2700E-01  1.4321E-01 -1.5094E-01  2.6300E-01  1.3362E-01  8.4023E-02  5.2687E-01 -9.8575E-02  1.9152E-01
             6.1068E-02
 GRADIENT:   7.7875E+05  1.8237E+05  4.5682E-02  2.5796E+05  3.5764E-02 -8.6116E-02 -2.6441E-02 -1.4794E+05  3.8937E+05 -4.0662E+05
             1.5852E-02
 NUMSIGDIG:         3.3         3.3         2.5         3.3         3.9         2.5         2.6         3.3         3.3         3.3
                    3.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1695
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0281E-04 -8.9226E-03 -4.0564E-02  3.9883E-03 -4.0028E-02
 SE:             2.9864E-02  2.4089E-02  1.4148E-02  1.8872E-02  2.1509E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8390E-01  7.1108E-01  4.1411E-03  8.3263E-01  6.2747E-02

 ETASHRINKSD(%)  1.0000E-10  1.9300E+01  5.2604E+01  3.6776E+01  2.7942E+01
 ETASHRINKVR(%)  1.0000E-10  3.4875E+01  7.7536E+01  6.0028E+01  4.8076E+01
 EBVSHRINKSD(%)  3.7749E-01  1.9164E+01  5.7784E+01  3.9094E+01  2.3569E+01
 EBVSHRINKVR(%)  7.5355E-01  3.4655E+01  8.2178E+01  6.2904E+01  4.1583E+01
 RELATIVEINF(%)  9.8774E+01  1.7979E+00  1.3636E+00  9.2052E-01  1.5897E+01
 EPSSHRINKSD(%)  4.5018E+01
 EPSSHRINKVR(%)  6.9769E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1676.1831130396870     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -941.03228647594881     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.11
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1676.183       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.63E-01  1.39E+00  1.04E+00  7.78E-01  1.18E+00  1.03E+00  9.84E-01  1.53E+00  8.20E-01  1.10E+00  9.62E-01
 


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
+        2.10E+09
 
 TH 2
+        4.90E+02  5.55E+07
 
 TH 3
+        1.30E+05  2.12E+04  9.23E+01
 
 TH 4
+        2.50E+03 -1.92E+04  1.06E+05  1.41E+09
 
 TH 5
+        7.60E+03  1.15E+03 -1.35E+02  6.39E+03  2.03E+08
 
 TH 6
+        2.94E+03  4.75E+02  8.06E-01  2.40E+03 -6.65E-01  1.85E+02
 
 TH 7
+        1.62E+04  2.66E+03  1.16E+00  1.33E+04 -3.17E+00 -7.54E-01  9.79E+01
 
 TH 8
+       -3.66E+02 -7.54E+04 -1.55E+04 -3.80E+05 -9.06E+02 -3.51E+02 -1.93E+03  2.99E+07
 
 TH 9
+        2.47E+09 -7.30E+03  1.53E+05 -3.68E+04  8.93E+03  3.45E+03  1.91E+04 -2.94E+08  2.90E+09
 
 TH10
+       -1.40E+03 -3.46E+02 -5.97E+04 -1.71E+03 -3.55E+03 -1.35E+03 -7.43E+03  1.15E+08  2.06E+04  4.42E+08
 
 TH11
+       -6.53E+04 -1.06E+04 -8.52E+00 -5.35E+04 -4.41E+00  1.58E+00  6.96E+00  7.80E+03 -7.67E+04  3.00E+04  2.32E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.570
Stop Time:
Sat Sep 18 13:47:36 CDT 2021
