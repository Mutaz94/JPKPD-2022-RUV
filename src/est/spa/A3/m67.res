Sat Sep 18 10:37:48 CDT 2021
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
$DATA ../../../../data/spa/A3/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -223.976766522495        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.8376E+01 -2.8034E+01  1.0475E+02 -2.0728E+02  1.2233E+02 -2.5593E+00 -5.8291E+01 -4.6442E+01 -1.8202E+02 -1.0605E+02
            -2.4367E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1205.16648663419        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0388E+00  9.8584E-01  8.7102E-01  1.2788E+00  8.8343E-01  8.8325E-01  1.0534E+00  1.0151E+00  1.2657E+00  9.9935E-01
             5.2809E+00
 PARAMETER:  1.3804E-01  8.5737E-02 -3.8090E-02  3.4594E-01 -2.3939E-02 -2.4149E-02  1.5206E-01  1.1503E-01  3.3563E-01  9.9351E-02
             1.7641E+00
 GRADIENT:  -5.9129E+00  3.1014E+00 -1.3633E+01  2.7816E+01 -1.3971E+01 -4.4238E+00  1.0420E+01  9.0457E+00  3.4519E+01  2.2953E+01
             2.1137E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1242.34793161511        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0107E+00  6.1917E-01  3.4848E-01  1.3846E+00  3.8604E-01  9.6224E-01  1.5350E+00  2.1226E-01  1.3286E+00  2.7506E-01
             4.2404E+00
 PARAMETER:  1.1069E-01 -3.7937E-01 -9.5416E-01  4.2540E-01 -8.5180E-01  6.1512E-02  5.2852E-01 -1.4499E+00  3.8416E-01 -1.1907E+00
             1.5447E+00
 GRADIENT:  -3.1879E+01  4.8199E+01 -7.2583E+00  8.4455E+01 -2.2740E+01  1.8089E+00  2.3269E+01  3.7789E-01  3.4892E+01  2.4943E+00
             1.3972E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1279.66730095288        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.7734E-01  5.0371E-01  2.0174E-01  1.1633E+00  2.6920E-01  1.0064E+00  1.2997E+00  8.3572E-01  1.1708E+00  1.6828E-01
             2.9719E+00
 PARAMETER:  7.7078E-02 -5.8575E-01 -1.5008E+00  2.5129E-01 -1.2123E+00  1.0638E-01  3.6217E-01 -7.9458E-02  2.5769E-01 -1.6821E+00
             1.1892E+00
 GRADIENT:  -4.0083E+01  3.9777E+01  2.0332E+01  2.9678E+01 -6.4170E+01  4.1151E+00  9.3094E+00 -3.6139E+00 -2.0843E+01 -4.1546E-01
             3.5069E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1286.28829741147        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.8170E-01  5.0639E-01  1.6870E-01  1.0964E+00  2.5490E-01  9.8240E-01  9.6149E-01  1.3427E+00  1.4412E+00  1.6080E-01
             2.6716E+00
 PARAMETER:  8.1534E-02 -5.8045E-01 -1.6796E+00  1.9200E-01 -1.2669E+00  8.2246E-02  6.0725E-02  3.9465E-01  4.6550E-01 -1.7276E+00
             1.0827E+00
 GRADIENT:   4.3535E+00  4.9154E+00  7.4211E+00 -1.9345E+00 -1.3589E+01 -9.9547E-01 -1.4920E+00 -5.0676E-01 -3.6027E-01 -8.5903E-01
            -4.6743E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1286.46720270123        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  9.7797E-01  5.0622E-01  1.6581E-01  1.0963E+00  2.5484E-01  9.8367E-01  9.7662E-01  1.3378E+00  1.4508E+00  1.9890E-01
             2.6816E+00
 PARAMETER:  7.7723E-02 -5.8078E-01 -1.6969E+00  1.9195E-01 -1.2671E+00  8.3539E-02  7.6343E-02  3.9103E-01  4.7214E-01 -1.5149E+00
             1.0864E+00
 GRADIENT:  -6.2514E+00 -4.3252E+00 -4.2622E+00 -1.5384E+00  5.4084E-01 -8.1892E-01  1.0763E+00  1.6211E-01  1.0211E+00 -6.3329E-01
            -1.1612E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1286.69291068174        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      576
 NPARAMETR:  9.8058E-01  5.3566E-01  1.7160E-01  1.0983E+00  2.6508E-01  9.8400E-01  9.3877E-01  1.3538E+00  1.4225E+00  1.9904E-01
             2.6961E+00
 PARAMETER:  8.0389E-02 -5.2426E-01 -1.6626E+00  1.9374E-01 -1.2277E+00  8.3872E-02  3.6816E-02  4.0294E-01  4.5238E-01 -1.5143E+00
             1.0918E+00
 GRADIENT:  -2.7431E-01  7.9085E-01  3.1878E-01  9.3076E-01 -1.7430E+00  1.1890E-02 -2.5603E-01 -9.4482E-02 -5.4394E-01 -3.3236E-03
            -4.3719E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1286.71097022373        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      755
 NPARAMETR:  9.8019E-01  5.5635E-01  1.7466E-01  1.0945E+00  2.7276E-01  9.8242E-01  9.5616E-01  1.3782E+00  1.4176E+00  1.3092E-01
             2.7022E+00
 PARAMETER:  7.9994E-02 -4.8637E-01 -1.6449E+00  1.9034E-01 -1.1992E+00  8.2259E-02  5.5170E-02  4.2078E-01  4.4897E-01 -1.9332E+00
             1.0941E+00
 GRADIENT:   1.3830E-01 -1.8311E-01 -3.4461E-02 -2.6791E-01  1.4389E+00  8.7280E-03  2.2452E-01  1.6369E-01  3.0839E-01  8.0931E-02
             2.7519E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1286.72976449543        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      930
 NPARAMETR:  9.7990E-01  5.6446E-01  1.7582E-01  1.0937E+00  2.7572E-01  9.8198E-01  9.6505E-01  1.3871E+00  1.4145E+00  4.4143E-02
             2.7044E+00
 PARAMETER:  7.9695E-02 -4.7188E-01 -1.6383E+00  1.8954E-01 -1.1884E+00  8.1812E-02  6.4425E-02  4.2725E-01  4.4677E-01 -3.0203E+00
             1.0949E+00
 GRADIENT:  -1.2016E-01  1.1154E-01  2.7266E-02  9.7706E-02  2.1317E-01  1.3328E-02  3.3739E-02  4.9058E-03  9.0937E-02  5.4839E-03
             5.7573E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1286.73198118640        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1105
 NPARAMETR:  9.7996E-01  5.6462E-01  1.7590E-01  1.0936E+00  2.7585E-01  9.8191E-01  9.6689E-01  1.3880E+00  1.4139E+00  1.0000E-02
             2.7044E+00
 PARAMETER:  7.9760E-02 -4.7161E-01 -1.6378E+00  1.8947E-01 -1.1879E+00  8.1745E-02  6.6325E-02  4.2788E-01  4.4638E-01 -4.7041E+00
             1.0949E+00
 GRADIENT:   7.7365E-03 -8.9415E-03 -1.1700E-02  1.7640E-02  6.6644E-02 -2.7950E-03  1.4099E-02  5.0128E-03 -2.7456E-03  0.0000E+00
            -1.3087E-03

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1286.73198697047        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1197
 NPARAMETR:  9.7996E-01  5.6451E-01  1.7588E-01  1.0936E+00  2.7580E-01  9.8193E-01  9.6676E-01  1.3879E+00  1.4141E+00  1.0000E-02
             2.7044E+00
 PARAMETER:  7.9759E-02 -4.7180E-01 -1.6379E+00  1.8946E-01 -1.1881E+00  8.1762E-02  6.6190E-02  4.2780E-01  4.4646E-01 -4.6805E+00
             1.0949E+00
 GRADIENT:   4.5789E-05  4.9813E-04  1.4066E-04  4.6229E-04 -7.2214E-04  2.8535E-04 -1.3703E-04 -1.5854E-05 -3.8738E-04  0.0000E+00
            -3.8610E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1197
 NO. OF SIG. DIGITS IN FINAL EST.:  4.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8407E-03 -6.5363E-03  2.1109E-03 -7.8893E-03  3.2302E-04
 SE:             2.8889E-02  1.8208E-02  1.9589E-02  2.6420E-02  3.5026E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4920E-01  7.1961E-01  9.1419E-01  7.6524E-01  3.5641E-01

 ETASHRINKSD(%)  3.2169E+00  3.9001E+01  3.4373E+01  1.1488E+01  9.8827E+01
 ETASHRINKVR(%)  6.3304E+00  6.2791E+01  5.6931E+01  2.1657E+01  9.9986E+01
 EBVSHRINKSD(%)  3.1400E+00  3.8539E+01  3.4885E+01  1.0337E+01  9.8898E+01
 EBVSHRINKVR(%)  6.1813E+00  6.2225E+01  5.7600E+01  1.9605E+01  9.9988E+01
 RELATIVEINF(%)  8.9772E+01  1.8495E+00  1.1006E+01  5.3079E+01  4.5585E-04
 EPSSHRINKSD(%)  3.6918E+01
 EPSSHRINKVR(%)  6.0207E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1286.7319869704711     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -551.58116040673292     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.28
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1286.732       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  5.65E-01  1.76E-01  1.09E+00  2.76E-01  9.82E-01  9.67E-01  1.39E+00  1.41E+00  1.00E-02  2.70E+00
 


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
+        1.13E+03
 
 TH 2
+       -2.20E+00  1.47E+03
 
 TH 3
+       -2.79E+02  2.17E+03  8.09E+03
 
 TH 4
+       -2.02E+01  1.44E+02 -3.59E+02  3.71E+02
 
 TH 5
+        1.99E+02 -4.95E+03 -9.15E+03 -1.33E+02  1.95E+04
 
 TH 6
+        5.03E+00 -1.04E+01  1.87E+00 -8.01E+00  5.44E+01  1.80E+02
 
 TH 7
+       -3.05E-01  9.17E+00 -5.16E+01 -3.32E+00  7.10E+01  3.79E-01  3.77E+01
 
 TH 8
+        3.40E+00 -2.97E-01 -4.64E+01 -1.47E+00 -3.49E+01  3.93E+00  3.64E+00  3.00E+01
 
 TH 9
+        1.22E+01 -2.92E+01  5.97E+01 -8.14E+00  1.72E+02  8.61E-01  9.43E+00 -3.78E+00  5.93E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.99E+01 -1.24E+01 -8.27E+00 -3.55E-01  3.61E+01  2.18E+00  6.96E+00  6.18E+00  7.32E+00  0.00E+00  3.28E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.979
Stop Time:
Sat Sep 18 10:38:11 CDT 2021
