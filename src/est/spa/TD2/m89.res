Sat Sep 18 14:58:24 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat89.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1641.71342223120        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9739E+01 -7.1313E+01 -1.9146E+01 -6.0563E+01  6.7665E+01  3.7238E+01 -5.5560E+00 -3.0214E+00  9.8219E+00 -1.5162E+01
            -2.1834E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1648.32074595060        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9291E-01  1.0439E+00  9.6696E-01  1.0290E+00  9.5823E-01  8.9239E-01  1.0194E+00  1.0243E+00  9.3482E-01  1.0331E+00
             1.0510E+00
 PARAMETER:  9.2883E-02  1.4294E-01  6.6397E-02  1.2854E-01  5.7336E-02 -1.3849E-02  1.1923E-01  1.2396E-01  3.2601E-02  1.3256E-01
             1.4979E-01
 GRADIENT:   2.5559E+01  3.8994E+00 -5.9092E+00  1.6503E+01  1.1897E+01 -4.0569E+00 -3.0819E+00  1.0322E-01  8.0947E-01  9.5465E-01
            -7.8722E-02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1648.77080616042        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9153E-01  9.2766E-01  8.7539E-01  1.1015E+00  8.4581E-01  9.0095E-01  1.2702E+00  9.0501E-01  8.0822E-01  8.6827E-01
             1.0618E+00
 PARAMETER:  9.1497E-02  2.4909E-02 -3.3090E-02  1.9668E-01 -6.7463E-02 -4.3009E-03  3.3914E-01  1.8760E-04 -1.1292E-01 -4.1250E-02
             1.5995E-01
 GRADIENT:   1.9468E+01  2.0566E+01 -5.0003E-01  3.5071E+01  2.5759E+00 -5.7284E-01  3.1980E-01  1.5540E+00 -7.3880E+00 -3.3562E+00
             2.3349E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1649.35519397189        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.8640E-01  9.1997E-01  6.3657E-01  1.0663E+00  7.1295E-01  9.0584E-01  1.2727E+00  5.5384E-01  8.1968E-01  7.5735E-01
             1.0502E+00
 PARAMETER:  8.6310E-02  1.6588E-02 -3.5165E-01  1.6419E-01 -2.3835E-01  1.1101E-03  3.4117E-01 -4.9087E-01 -9.8844E-02 -1.7792E-01
             1.4896E-01
 GRADIENT:  -3.7090E-01  6.4176E+00 -1.0765E+01  1.8400E+01  1.0249E+01  3.8958E-01  1.6281E+00  2.8123E+00 -2.7931E+00  2.9194E+00
            -1.0358E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1649.45681400654        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  9.8564E-01  8.7222E-01  5.4876E-01  1.0701E+00  6.3762E-01  9.0662E-01  1.3145E+00  4.0950E-01  8.0523E-01  6.6207E-01
             1.0527E+00
 PARAMETER:  8.5536E-02 -3.6717E-02 -5.0009E-01  1.6774E-01 -3.5001E-01  1.9700E-03  3.7348E-01 -7.9283E-01 -1.1663E-01 -3.1238E-01
             1.5137E-01
 GRADIENT:  -4.7569E+00  2.7828E+00 -9.4568E+00  1.0185E+01  8.7590E+00  3.5474E-01  8.6864E-01  2.2533E+00 -1.3841E+00  2.9532E+00
            -1.2402E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1650.90935709124        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      466
 NPARAMETR:  1.0028E+00  6.9604E-01  6.1471E-01  1.1799E+00  6.1749E-01  9.1280E-01  1.6106E+00  3.1416E-01  7.6905E-01  6.9206E-01
             1.0598E+00
 PARAMETER:  1.0280E-01 -2.6235E-01 -3.8661E-01  2.6541E-01 -3.8209E-01  8.7657E-03  5.7663E-01 -1.0579E+00 -1.6260E-01 -2.6809E-01
             1.5804E-01
 GRADIENT:   1.1223E+01  1.6124E+00 -2.1604E-01  3.3843E+00  2.6128E+00  1.1064E+00  1.1177E+00 -2.2617E-01  1.0048E+00 -8.9981E-01
            -8.1322E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1650.96639930086        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      641
 NPARAMETR:  9.9877E-01  6.7833E-01  6.0154E-01  1.1852E+00  6.0303E-01  9.1025E-01  1.6257E+00  2.9614E-01  7.5985E-01  6.8598E-01
             1.0605E+00
 PARAMETER:  9.8769E-02 -2.8813E-01 -4.0827E-01  2.6987E-01 -4.0579E-01  5.9648E-03  5.8593E-01 -1.1169E+00 -1.7464E-01 -2.7691E-01
             1.5875E-01
 GRADIENT:   3.1765E-01  6.4914E-01  3.3380E-01  3.6890E-01 -5.9253E-01  4.1506E-02  8.4455E-03 -9.4366E-02 -2.3481E-01  2.5275E-02
             4.3300E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1650.99850923708        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  9.9931E-01  6.4601E-01  6.5227E-01  1.2155E+00  6.2116E-01  9.1052E-01  1.7069E+00  4.4345E-01  7.4841E-01  7.0861E-01
             1.0611E+00
 PARAMETER:  9.9310E-02 -3.3694E-01 -3.2730E-01  2.9512E-01 -3.7617E-01  6.2616E-03  6.3468E-01 -7.1317E-01 -1.8980E-01 -2.4445E-01
             1.5928E-01
 GRADIENT:   4.3077E+00  3.5575E+00  8.7585E-01  2.0815E+00 -4.0023E+00  4.0796E-01  1.2340E+00  2.9096E-01 -4.6257E-01  6.1528E-01
             1.0353E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1651.30347394434        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  9.9294E-01  4.6754E-01  8.3760E-01  1.3469E+00  6.7068E-01  9.0739E-01  2.1238E+00  6.8391E-01  7.0883E-01  7.9974E-01
             1.0550E+00
 PARAMETER:  9.2919E-02 -6.6028E-01 -7.7221E-02  3.9780E-01 -2.9947E-01  2.8146E-03  8.5322E-01 -2.7992E-01 -2.4414E-01 -1.2347E-01
             1.5350E-01
 GRADIENT:   2.9714E-01  3.8748E+00  7.5126E+00  8.7202E+00 -9.9014E+00  2.1854E-01 -3.0087E-01 -1.5497E+00 -3.2223E-01 -7.8691E-01
            -3.4713E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1651.58103260933        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  9.8875E-01  3.4679E-01  1.0860E+00  1.4397E+00  7.5990E-01  9.0417E-01  2.5601E+00  1.0026E+00  6.8433E-01  8.8574E-01
             1.0521E+00
 PARAMETER:  8.8683E-02 -9.5904E-01  1.8247E-01  4.6441E-01 -1.7456E-01 -7.3353E-04  1.0400E+00  1.0264E-01 -2.7931E-01 -2.1328E-02
             1.5082E-01
 GRADIENT:   3.8932E-01  8.9134E-01  3.0345E-01  8.8899E-01 -6.3060E-01  3.2984E-02  2.0501E-01  4.5635E-01 -1.7900E-01 -1.6560E-01
             8.2122E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1651.61527963022        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1353
 NPARAMETR:  9.8739E-01  2.8497E-01  1.1492E+00  1.4791E+00  7.7119E-01  9.0345E-01  2.8538E+00  1.0660E+00  6.7611E-01  9.0409E-01
             1.0518E+00
 PARAMETER:  8.7306E-02 -1.1554E+00  2.3908E-01  4.9145E-01 -1.5983E-01 -1.5310E-03  1.1486E+00  1.6393E-01 -2.9140E-01 -8.2688E-04
             1.5047E-01
 GRADIENT:   4.4066E-01 -2.0406E-01  2.6279E-01 -1.7259E+00  1.2900E-01  3.6654E-02  3.8980E-02 -2.9527E-02  9.6383E-02 -1.4578E-01
             7.1279E-02

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1651.61900707761        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1481
 NPARAMETR:  9.8764E-01  2.9793E-01  1.1069E+00  1.4684E+00  7.5707E-01  9.0366E-01  2.7867E+00  1.0211E+00  6.7798E-01  8.9235E-01
             1.0517E+00
 PARAMETER:  8.7561E-02 -1.1109E+00  2.0159E-01  4.8416E-01 -1.7829E-01 -1.2970E-03  1.1249E+00  1.2084E-01 -2.8863E-01 -1.3898E-02
             1.5044E-01
 GRADIENT:  -7.7515E-03 -2.9215E-03 -5.1178E-03  1.5206E-02  4.9112E-03  2.9677E-05 -1.0195E-02  2.1119E-03 -3.1386E-03 -7.3532E-04
            -5.5562E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1481
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0325E-03  3.1342E-02 -3.1643E-02 -2.6616E-02 -2.0105E-02
 SE:             2.9813E-02  1.7327E-02  1.6930E-02  2.4554E-02  2.0275E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7237E-01  7.0480E-02  6.1614E-02  2.7837E-01  3.2138E-01

 ETASHRINKSD(%)  1.2139E-01  4.1951E+01  4.3282E+01  1.7742E+01  3.2077E+01
 ETASHRINKVR(%)  2.4263E-01  6.6303E+01  6.7831E+01  3.2336E+01  5.3865E+01
 EBVSHRINKSD(%)  5.6926E-01  4.8246E+01  4.5159E+01  1.4822E+01  2.8317E+01
 EBVSHRINKVR(%)  1.1353E+00  7.3215E+01  6.9924E+01  2.7447E+01  4.8616E+01
 RELATIVEINF(%)  9.8018E+01  4.2477E+00  4.4976E+00  1.3124E+01  7.0697E+00
 EPSSHRINKSD(%)  4.4811E+01
 EPSSHRINKVR(%)  6.9542E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1651.6190070776079     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -916.46818051386970     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1651.619       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  2.98E-01  1.11E+00  1.47E+00  7.57E-01  9.04E-01  2.79E+00  1.02E+00  6.78E-01  8.92E-01  1.05E+00
 


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
+        1.37E+03
 
 TH 2
+       -2.76E+01  5.11E+02
 
 TH 3
+        2.66E+00  7.68E+01  2.78E+02
 
 TH 4
+       -8.83E+00  4.27E+02 -8.41E+01  8.40E+02
 
 TH 5
+       -1.58E+00 -2.84E+02 -5.08E+02  2.85E+01  1.27E+03
 
 TH 6
+        1.10E+00 -3.85E+00  2.00E+00 -1.74E+00  2.86E+00  2.37E+02
 
 TH 7
+        4.81E-01  3.44E+01 -8.15E-01 -5.94E+00  1.82E+00 -1.77E-01  7.35E+00
 
 TH 8
+        3.23E+00  3.77E-01 -4.87E+01 -7.02E+00  1.32E+00  2.53E+00  5.39E-01  3.67E+01
 
 TH 9
+        5.10E+00 -2.99E+01  2.65E+00 -1.35E+01  7.06E+00 -2.88E+00  2.50E+00 -9.94E-01  2.70E+02
 
 TH10
+        5.30E+00  1.35E+01 -7.13E+00 -2.37E+01 -7.94E+01 -5.84E+00  2.23E+00  1.96E+01 -1.84E+00  7.43E+01
 
 TH11
+       -9.29E+00 -3.93E+00 -6.35E+00 -8.63E+00 -4.56E+00  1.19E-01  1.92E-01  2.28E+00  1.13E+01  1.13E+01  1.92E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.227
Stop Time:
Sat Sep 18 14:58:50 CDT 2021
