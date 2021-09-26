Sat Sep 25 10:07:11 CDT 2021
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
$DATA ../../../../data/spa/S1/dat81.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1685.71264958426        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.2902E+00 -6.6757E+01 -3.1054E+01 -2.9526E+01  5.6541E+01 -3.7021E+01 -8.7083E+00 -7.8038E-01  2.9754E+01 -3.0840E+00
             5.5745E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1691.64048265220        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9168E-01  1.0775E+00  9.8542E-01  9.6492E-01  1.0468E+00  1.0932E+00  1.0717E+00  1.0250E+00  8.5327E-01  9.9851E-01
             9.8786E-01
 PARAMETER:  9.1643E-02  1.7466E-01  8.5312E-02  6.4285E-02  1.4571E-01  1.8915E-01  1.6925E-01  1.2472E-01 -5.8681E-02  9.8508E-02
             8.7785E-02
 GRADIENT:  -1.2108E+01 -4.8748E+01 -3.1196E+01 -5.4889E+00  7.6298E+01  1.0711E+01 -1.0814E+01 -1.1957E+00  7.3315E+00 -6.5698E+00
            -9.4926E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1692.18779086941        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  9.8588E-01  1.1527E+00  7.4625E-01  9.0307E-01  9.2504E-01  1.1149E+00  1.1962E+00  1.0744E+00  8.0315E-01  8.0396E-01
             9.4210E-01
 PARAMETER:  8.5776E-02  2.4213E-01 -1.9270E-01 -1.9525E-03  2.2081E-02  2.0879E-01  2.7918E-01  1.7178E-01 -1.1921E-01 -1.1820E-01
             4.0361E-02
 GRADIENT:  -1.9285E+01 -2.5102E+01 -3.0206E+01  2.5359E+00  5.2605E+01  2.0581E+01  9.3292E+00  5.8839E+00  6.8261E+00 -4.6857E+00
            -2.0439E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1693.79807001092        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:      246
 NPARAMETR:  9.9232E-01  1.3030E+00  6.2406E-01  8.2761E-01  9.0128E-01  1.1088E+00  1.1310E+00  8.1018E-01  7.9606E-01  8.2374E-01
             9.5082E-01
 PARAMETER:  9.2294E-02  3.6463E-01 -3.7151E-01 -8.9215E-02 -3.9398E-03  2.0331E-01  2.2312E-01 -1.1050E-01 -1.2808E-01 -9.3906E-02
             4.9567E-02
 GRADIENT:  -1.2081E+01  3.3778E+01 -1.0767E+01  1.6514E+01  1.2244E+01  1.7130E+01  1.4508E+01  1.7225E+00 -1.2679E+00 -5.5018E-01
            -1.7209E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1694.87868544757        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  1.0011E+00  1.2904E+00  6.5328E-01  8.3691E-01  9.1255E-01  1.1090E+00  1.1215E+00  8.1547E-01  8.0124E-01  8.4125E-01
             9.6087E-01
 PARAMETER:  1.0112E-01  3.5494E-01 -3.2576E-01 -7.8037E-02  8.4869E-03  2.0345E-01  2.1467E-01 -1.0400E-01 -1.2159E-01 -7.2871E-02
             6.0079E-02
 GRADIENT:  -4.7364E+01  6.9019E+00 -8.7986E+00  7.1927E+00  8.6015E+00 -2.8404E+00  9.7617E+00  1.2627E+00 -1.2047E+00 -4.8894E-01
            -1.2758E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1694.89741075044        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      599
 NPARAMETR:  1.0012E+00  1.2904E+00  6.5373E-01  8.3711E-01  9.1259E-01  1.1123E+00  1.1214E+00  8.1557E-01  8.0131E-01  8.4137E-01
             9.6091E-01
 PARAMETER:  1.0125E-01  3.5494E-01 -3.2506E-01 -7.7797E-02  8.5340E-03  2.0640E-01  2.1454E-01 -1.0387E-01 -1.2151E-01 -7.2727E-02
             6.0125E-02
 GRADIENT:  -4.6848E+01  7.1823E+00 -8.6370E+00  7.2678E+00  8.3267E+00 -1.6540E+00  9.7125E+00  1.2505E+00 -1.2095E+00 -4.9379E-01
            -1.2738E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1696.01270831042        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      755
 NPARAMETR:  1.0193E+00  1.2901E+00  6.6204E-01  8.3432E-01  9.1071E-01  1.1026E+00  1.0645E+00  7.2215E-01  8.1818E-01  8.4963E-01
             9.9086E-01
 PARAMETER:  1.1910E-01  3.5476E-01 -3.1243E-01 -8.1139E-02  6.4660E-03  1.9769E-01  1.6252E-01 -2.2552E-01 -1.0067E-01 -6.2951E-02
             9.0815E-02
 GRADIENT:  -1.5474E+01  5.6698E+00  3.9938E+00 -5.0605E+00 -9.7937E+00 -3.8678E+00  1.7145E-01 -3.6445E-01 -2.0308E+00 -9.5239E-01
            -6.9203E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1696.21860873454        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      932
 NPARAMETR:  1.0253E+00  1.2901E+00  7.0872E-01  8.3587E-01  9.4701E-01  1.1102E+00  1.0542E+00  8.0436E-01  8.3751E-01  9.0052E-01
             9.9087E-01
 PARAMETER:  1.2495E-01  3.5476E-01 -2.4429E-01 -7.9288E-02  4.5552E-02  2.0450E-01  1.5275E-01 -1.1771E-01 -7.7320E-02 -4.7834E-03
             9.0827E-02
 GRADIENT:  -4.0021E+00 -3.8722E+00 -1.0184E+00 -4.2122E+00  3.9219E-01 -8.0554E-01  2.8544E-02 -1.5044E-01  2.8148E-01  1.1910E+00
            -7.4723E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1696.29134141050        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1107
 NPARAMETR:  1.0270E+00  1.2901E+00  7.7770E-01  8.4473E-01  9.8109E-01  1.1115E+00  1.0535E+00  9.8156E-01  8.3498E-01  9.1874E-01
             9.9203E-01
 PARAMETER:  1.2666E-01  3.5476E-01 -1.5141E-01 -6.8738E-02  8.0908E-02  2.0569E-01  1.5209E-01  8.1390E-02 -8.0348E-02  1.5253E-02
             9.1995E-02
 GRADIENT:   8.9747E-03  1.7910E+00  5.9537E-03 -8.5694E-04 -1.3282E-02  1.8458E-03  4.8502E-03  2.2406E-04  3.1769E-03  2.5194E-03
            -2.4539E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1696.29366195432        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1285
 NPARAMETR:  1.0270E+00  1.2884E+00  7.7844E-01  8.4575E-01  9.8040E-01  1.1114E+00  1.0542E+00  9.7958E-01  8.3441E-01  9.1845E-01
             9.9191E-01
 PARAMETER:  1.2663E-01  3.5341E-01 -1.5046E-01 -6.7530E-02  8.0210E-02  2.0565E-01  1.5282E-01  7.9371E-02 -8.1035E-02  1.4930E-02
             9.1878E-02
 GRADIENT:  -3.1243E-02  1.7903E+00  1.8858E-01 -2.0795E-01 -3.0534E-01 -1.4891E-02 -6.5842E-02 -3.1431E-02 -1.5066E-02 -4.3135E-03
            -5.5753E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1696.36303504781        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1469
 NPARAMETR:  1.0264E+00  1.2194E+00  8.0347E-01  8.8942E-01  9.5982E-01  1.1104E+00  1.1020E+00  9.3922E-01  8.1050E-01  9.0875E-01
             9.9154E-01
 PARAMETER:  1.2608E-01  2.9835E-01 -1.1882E-01 -1.7181E-02  5.8995E-02  2.0472E-01  1.9712E-01  3.7298E-02 -1.1010E-01  4.3153E-03
             9.1502E-02
 GRADIENT:  -5.8464E-01  8.9347E-01 -1.4392E+00  2.5895E+00  2.2479E+00 -3.1413E-01 -2.8808E-01 -2.9025E-02  5.1906E-01  1.4897E-01
             1.3672E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1696.36628811293        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1653
 NPARAMETR:  1.0265E+00  1.2181E+00  8.0467E-01  8.9006E-01  9.5911E-01  1.1107E+00  1.1034E+00  9.3891E-01  8.1002E-01  9.0828E-01
             9.9144E-01
 PARAMETER:  1.2617E-01  2.9727E-01 -1.1733E-01 -1.6470E-02  5.8251E-02  2.0496E-01  1.9839E-01  3.6969E-02 -1.1069E-01  3.7997E-03
             9.1407E-02
 GRADIENT:  -1.1298E+05  4.2497E+05  2.1535E+06  1.2632E+06  1.2634E+06 -6.1637E+05  5.4032E+03 -1.2633E+06  1.1412E+06 -2.5267E+06
             1.2631E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1653
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1970E-04 -6.1723E-03 -2.9297E-02 -2.8958E-04 -2.7923E-02
 SE:             2.9891E-02  2.3642E-02  1.2456E-02  2.1017E-02  2.1289E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8880E-01  7.9404E-01  1.8674E-02  9.8901E-01  1.8965E-01

 ETASHRINKSD(%)  1.0000E-10  2.0795E+01  5.8270E+01  2.9590E+01  2.8680E+01
 ETASHRINKVR(%)  1.0000E-10  3.7265E+01  8.2586E+01  5.0424E+01  4.9134E+01
 EBVSHRINKSD(%)  3.5098E-01  2.1010E+01  6.1663E+01  3.0452E+01  2.6175E+01
 EBVSHRINKVR(%)  7.0072E-01  3.7605E+01  8.5303E+01  5.1630E+01  4.5498E+01
 RELATIVEINF(%)  9.8945E+01  2.9715E+00  1.1021E+00  2.0854E+00  8.2927E+00
 EPSSHRINKSD(%)  4.4876E+01
 EPSSHRINKVR(%)  6.9613E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1696.3662881129269     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -961.21546154918872     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1696.366       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.22E+00  8.05E-01  8.90E-01  9.59E-01  1.11E+00  1.10E+00  9.39E-01  8.10E-01  9.08E-01  9.91E-01
 


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
+        3.77E+09
 
 TH 2
+       -2.24E+08  4.82E+08
 
 TH 3
+        3.58E+05  1.85E+09  7.09E+09
 
 TH 4
+        1.92E+04 -2.78E+04  2.08E+05  7.97E+09
 
 TH 5
+       -4.95E+09  1.76E+09  6.15E+09  1.82E+04  6.87E+09
 
 TH 6
+       -8.81E+04 -7.66E+08 -2.94E+09  3.29E+04 -7.49E+07  1.22E+09
 
 TH 7
+       -1.13E+09 -1.06E+04 -4.59E+04  1.34E+03  5.19E+08  1.14E+04  1.41E+08
 
 TH 8
+       -2.37E+04  3.19E+04 -1.78E+05 -1.40E+05 -3.52E+03 -4.12E+04 -3.42E+02  7.17E+09
 
 TH 9
+        2.86E+04 -3.71E+04  1.73E+05 -4.10E+05 -5.72E+03  5.02E+04 -2.66E+02 -2.74E+09  7.86E+09
 
 TH10
+        5.22E+09 -9.16E+08 -7.22E+09  1.97E+04 -7.25E+09  1.31E+09  1.14E+05 -3.21E+04  4.31E+04  7.66E+09
 
 TH11
+        1.02E+05 -1.25E+05 -1.75E+05 -2.28E+06 -2.04E+05 -4.76E+08 -1.30E+04  2.20E+06 -2.33E+06  2.36E+05  6.42E+09
 
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
 #CPUT: Total CPU Time in Seconds,       30.502
Stop Time:
Sat Sep 25 10:07:43 CDT 2021
