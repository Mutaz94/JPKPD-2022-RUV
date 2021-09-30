Thu Sep 30 00:45:55 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat99.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   321.963449988074        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2029E+02  7.6264E+01  1.1811E+02  6.7135E+01  2.5546E+02  6.4293E+01 -5.0670E+01 -1.0748E+02 -4.5326E+00 -1.7459E+02
            -4.4942E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1467.19135062252        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0708E+00  1.0652E+00  9.9437E-01  1.1136E+00  8.9313E-01  7.1060E-01  9.5922E-01  9.5382E-01  8.5975E-01  9.7412E-01
             5.1151E+00
 PARAMETER:  1.6837E-01  1.6318E-01  9.4358E-02  2.0757E-01 -1.3021E-02 -2.4165E-01  5.8364E-02  5.2724E-02 -5.1112E-02  7.3778E-02
             1.7322E+00
 GRADIENT:   6.5617E+01  3.8802E+01 -3.8093E+00  6.1401E+01 -3.1218E+01 -4.3030E+01  6.0701E+00  5.3641E+00  1.1945E+01  2.1076E+01
             2.9258E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1515.05318299351        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0108E+00  7.3385E-01  5.0607E-01  1.2173E+00  5.8299E-01  8.1468E-01  1.5509E+00  2.0817E-02  9.0442E-01  3.8515E-01
             4.1394E+00
 PARAMETER:  1.1073E-01 -2.0944E-01 -5.8108E-01  2.9660E-01 -4.3959E-01 -1.0497E-01  5.3883E-01 -3.7720E+00 -4.5971E-04 -8.5413E-01
             1.5206E+00
 GRADIENT:  -5.1593E+01 -2.4913E+00 -8.8358E+01  1.1129E+02  1.0643E+02 -1.3987E+01  2.1934E+01  5.6493E-03  2.1369E+01  4.2087E+00
             1.7540E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1543.04749411500        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.0019E+00  6.4969E-01  5.3136E-01  1.1838E+00  5.2267E-01  8.7171E-01  1.4289E+00  1.2894E-02  7.6722E-01  4.5284E-01
             3.3844E+00
 PARAMETER:  1.0189E-01 -3.3126E-01 -5.3232E-01  2.6871E-01 -5.4880E-01 -3.7293E-02  4.5690E-01 -4.2510E+00 -1.6498E-01 -6.9221E-01
             1.3192E+00
 GRADIENT:  -2.0623E+01  1.6877E+01  7.0564E+00  3.4870E+01  7.8555E-01  5.3857E+00 -2.0695E-01  9.0749E-04 -8.5813E+00 -8.4630E-01
            -2.6170E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1544.83859647879        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      416
 NPARAMETR:  1.0228E+00  4.9765E-01  6.0172E-01  1.2794E+00  5.1385E-01  8.4459E-01  1.6790E+00  1.0000E-02  7.9417E-01  5.2143E-01
             3.4011E+00
 PARAMETER:  1.2255E-01 -5.9786E-01 -4.0797E-01  3.4638E-01 -5.6582E-01 -6.8906E-02  6.1822E-01 -4.5694E+00 -1.3046E-01 -5.5119E-01
             1.3241E+00
 GRADIENT:   1.0639E+01  1.3077E+01  2.0772E+01  1.1082E+01 -2.7981E+01 -5.2995E+00  4.9400E-01  0.0000E+00 -7.1990E-01 -2.0666E+00
            -1.1011E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1548.96907256423        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      592
 NPARAMETR:  1.0176E+00  2.4733E-01  3.1696E-01  1.2414E+00  2.9910E-01  8.4871E-01  1.6159E+00  1.0000E-02  8.8995E-01  6.4125E-01
             3.2490E+00
 PARAMETER:  1.1744E-01 -1.2970E+00 -1.0490E+00  3.1621E-01 -1.1070E+00 -6.4038E-02  5.7992E-01 -6.2144E+00 -1.6591E-02 -3.4434E-01
             1.2784E+00
 GRADIENT:  -1.3959E+00  1.7925E+00 -1.6060E+01  4.1087E+01  7.7455E+00 -7.2202E+00 -2.3595E+00  0.0000E+00 -3.5389E+00 -4.5210E+00
            -3.8332E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1550.08916141583        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      767
 NPARAMETR:  1.0144E+00  1.9470E-01  3.1299E-01  1.2264E+00  2.8703E-01  8.5977E-01  2.0516E+00  1.0000E-02  8.9051E-01  6.9431E-01
             3.2278E+00
 PARAMETER:  1.1432E-01 -1.5363E+00 -1.0616E+00  3.0409E-01 -1.1482E+00 -5.1095E-02  8.1863E-01 -6.9442E+00 -1.5966E-02 -2.6484E-01
             1.2718E+00
 GRADIENT:  -5.2010E-01  2.1936E+00  5.0693E+00 -8.0593E-01 -9.9885E+00 -9.8720E-01 -5.8268E-02  0.0000E+00 -5.5376E-01  1.9148E-01
            -1.2918E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1551.31107649461        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  1.0040E+00  8.6584E-02  3.9998E-01  1.3547E+00  3.2963E-01  8.5658E-01  4.1128E+00  1.0000E-02  8.3806E-01  6.8474E-01
             3.2889E+00
 PARAMETER:  1.0396E-01 -2.3466E+00 -8.1633E-01  4.0359E-01 -1.0098E+00 -5.4807E-02  1.5141E+00 -9.1594E+00 -7.6666E-02 -2.7872E-01
             1.2905E+00
 GRADIENT:  -6.9837E+00 -4.1609E-01  7.2014E+00  2.1099E+01 -1.2165E+01  2.3707E-01 -3.7738E+00  0.0000E+00  9.7049E-01  2.5208E+00
             8.1204E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1552.34538365968        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1120
 NPARAMETR:  1.0031E+00  4.6994E-02  3.4810E-01  1.3178E+00  2.9741E-01  8.5549E-01  5.9692E+00  1.0000E-02  8.5399E-01  6.5924E-01
             3.2329E+00
 PARAMETER:  1.0308E-01 -2.9577E+00 -9.5525E-01  3.7597E-01 -1.1127E+00 -5.6075E-02  1.8866E+00 -1.1090E+01 -5.7838E-02 -3.1667E-01
             1.2734E+00
 GRADIENT:  -1.3053E+00 -1.9392E+00  7.3164E+00  1.9371E+01 -1.5412E+01 -2.9095E-01 -5.3893E+00  0.0000E+00  6.4053E-01  5.7157E-01
            -4.7250E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1552.59905947641        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1297
 NPARAMETR:  1.0025E+00  4.0303E-02  3.2742E-01  1.2913E+00  2.8641E-01  8.5808E-01  6.4505E+00  1.0000E-02  8.7037E-01  6.5728E-01
             3.2162E+00
 PARAMETER:  1.0254E-01 -3.1113E+00 -1.0165E+00  3.5565E-01 -1.1503E+00 -5.3059E-02  1.9642E+00 -1.1653E+01 -3.8835E-02 -3.1965E-01
             1.2682E+00
 GRADIENT:   3.0959E+00  8.2029E+00 -5.4885E+00 -1.0539E+01  8.0265E+00 -1.7344E+00  1.5742E+01  0.0000E+00 -8.5286E+00 -3.1714E+00
            -5.2249E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1552.64399272920        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1475
 NPARAMETR:  1.0021E+00  3.7367E-02  3.2549E-01  1.2894E+00  2.8493E-01  8.5774E-01  6.7157E+00  1.0000E-02  8.7279E-01  6.6165E-01
             3.2218E+00
 PARAMETER:  1.0212E-01 -3.1870E+00 -1.0224E+00  3.5417E-01 -1.1555E+00 -5.3450E-02  2.0045E+00 -1.1888E+01 -3.6059E-02 -3.1301E-01
             1.2700E+00
 GRADIENT:   1.2234E+00  5.1023E+00 -3.0629E+00 -7.0036E+00  4.5059E+00 -1.0471E+00  9.8414E+00  0.0000E+00 -4.7303E+00 -1.5500E+00
            -1.4016E+00

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1552.66179656412        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1603
 NPARAMETR:  1.0022E+00  3.5425E-02  3.3221E-01  1.2976E+00  2.8883E-01  8.5749E-01  6.8927E+00  1.0000E-02  8.6557E-01  6.6050E-01
             3.2245E+00
 PARAMETER:  1.0220E-01 -3.2403E+00 -1.0020E+00  3.6055E-01 -1.1419E+00 -5.3750E-02  2.0305E+00 -1.2038E+01 -4.4364E-02 -3.1476E-01
             1.2708E+00
 GRADIENT:   6.3820E+00  1.9188E+01 -1.0454E+01 -2.7787E+01  1.7387E+01 -4.4768E+00  3.8462E+01  0.0000E+00 -2.0707E+01 -5.9105E+00
            -6.5600E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1603
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.7501E-04  1.3181E-02  1.1757E-04 -1.5536E-02 -2.3775E-03
 SE:             2.8630E-02  7.7352E-03  1.9723E-04  2.5967E-02  1.9845E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7562E-01  8.8368E-02  5.5109E-01  5.4964E-01  9.0464E-01

 ETASHRINKSD(%)  4.0859E+00  7.4086E+01  9.9339E+01  1.3006E+01  3.3518E+01
 ETASHRINKVR(%)  8.0049E+00  9.3285E+01  9.9996E+01  2.4321E+01  5.5801E+01
 EBVSHRINKSD(%)  3.8001E+00  8.1986E+01  9.9283E+01  1.1551E+01  3.3161E+01
 EBVSHRINKVR(%)  7.4558E+00  9.6755E+01  9.9995E+01  2.1767E+01  5.5326E+01
 RELATIVEINF(%)  9.0447E+01  1.9818E+00  2.2917E-04  1.6482E+01  2.0422E+00
 EPSSHRINKSD(%)  2.3219E+01
 EPSSHRINKVR(%)  4.1046E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1552.6617965641210     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -633.72326335944831     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1552.662       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  3.54E-02  3.32E-01  1.30E+00  2.89E-01  8.57E-01  6.89E+00  1.00E-02  8.66E-01  6.60E-01  3.22E+00
 


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
+        1.94E+05
 
 TH 2
+       -4.06E+02  1.50E+05
 
 TH 3
+        5.38E+00  6.56E+03  2.25E+04
 
 TH 4
+       -1.61E+01  1.87E+03 -8.87E+02  9.75E+03
 
 TH 5
+        1.48E+02 -1.25E+04 -6.69E+03 -3.22E+02  3.07E+04
 
 TH 6
+        3.17E+02  1.24E+02 -1.71E+01 -3.81E+01  3.05E+01  2.74E+05
 
 TH 7
+       -2.05E+00 -3.32E+01  5.11E+01  1.21E+01 -9.31E+01  1.04E+00  1.01E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.67E+02  1.09E+03 -1.03E+02 -1.58E+02  2.15E+02 -2.21E+02  8.51E+00  0.00E+00  2.68E+05
 
 TH10
+        1.01E+02  9.03E+02 -1.43E+02 -4.70E+01  1.76E+02 -9.35E+01  7.27E+00  0.00E+00 -4.80E+02  4.66E+04
 
 TH11
+       -1.79E+01 -1.97E+01 -2.60E+01 -1.90E+01  2.60E+01  8.12E-01  2.31E+00  0.00E+00 -8.32E+00  1.15E+01  1.68E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.853
Stop Time:
Thu Sep 30 00:46:33 CDT 2021
