Thu Sep 30 00:42:40 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat94.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1528.26510099977        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2759E+02  8.0598E+01  7.2258E+01  4.5673E+01  3.2096E+02  5.8081E+01 -1.0857E+02 -4.7623E+01 -1.3952E+02 -1.9258E+02
            -6.7580E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -893.890864002365        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.4448E+00  1.2466E+00  9.6956E-01  1.5696E+00  6.7487E-01  5.1498E-01  1.1009E+00  9.7508E-01  1.2044E+00  1.1080E+00
             1.5802E+01
 PARAMETER:  4.6794E-01  3.2039E-01  6.9090E-02  5.5083E-01 -2.9323E-01 -5.6362E-01  1.9614E-01  7.4768E-02  2.8595E-01  2.0256E-01
             2.8601E+00
 GRADIENT:   1.4124E+02 -1.0665E+01  6.1875E+00 -1.7671E+01 -2.2566E+01 -1.6656E+01  1.1550E+01  3.9116E+00  3.7729E+01  1.5750E+01
             4.6685E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1031.35924537234        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.1858E+00  4.5477E-01  5.1384E-02  1.5308E+00  1.0041E-01  5.4893E-01  7.5461E-01  1.3808E-01  2.3147E-01  8.2172E-01
             1.0411E+01
 PARAMETER:  2.7040E-01 -6.8795E-01 -2.8684E+00  5.2580E-01 -2.1985E+00 -4.9979E-01 -1.8155E-01 -1.8799E+00 -1.3633E+00 -9.6361E-02
             2.4429E+00
 GRADIENT:   2.6176E+02  2.3040E+02  4.9152E+01  3.1282E+01 -3.2860E+02 -3.4340E+01 -4.4980E+00  6.8781E-02 -2.8211E-01 -6.6093E+00
             2.2744E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1279.11094101962        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.3284E-01  2.5550E-01  4.8118E-02  1.3523E+00  1.1713E-01  8.6256E-01  4.2882E-01  1.0000E-02  9.3300E-02  8.1051E-02
             5.9322E+00
 PARAMETER:  3.0475E-02 -1.2645E+00 -2.9341E+00  4.0181E-01 -2.0445E+00 -4.7852E-02 -7.4672E-01 -7.5119E+00 -2.2719E+00 -2.4127E+00
             1.8804E+00
 GRADIENT:  -6.6037E+01  4.4262E+01  3.4616E+01 -1.6956E-01 -6.6559E+01  3.1691E+00  1.4922E+00  0.0000E+00 -2.5221E-02 -1.3386E+00
             1.5967E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1307.36461702606        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  9.7606E-01  4.3831E-01  1.8580E-01  1.1010E+00  2.2784E-01  8.2759E-01  2.7359E-01  1.0000E-02  1.8462E-01  1.3093E-01
             5.4753E+00
 PARAMETER:  7.5765E-02 -7.2484E-01 -1.5831E+00  1.9622E-01 -1.3791E+00 -8.9243E-02 -1.1961E+00 -7.9798E+00 -1.5895E+00 -1.9331E+00
             1.8003E+00
 GRADIENT:  -2.0025E+02  1.0175E+02  1.3606E+02  2.5336E+01 -3.3674E+02 -2.8318E+01 -2.2690E+00  0.0000E+00 -2.7172E+00 -4.4242E+00
            -4.0491E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1440.84123987113        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.0468E+00  3.4450E-01  2.4750E-01  1.0757E+00  2.7084E-01  8.6354E-01  1.0563E-01  1.0000E-02  1.5723E+00  1.7309E-01
             3.8496E+00
 PARAMETER:  1.4570E-01 -9.6566E-01 -1.2963E+00  1.7294E-01 -1.2062E+00 -4.6715E-02 -2.1478E+00 -1.0498E+01  5.5254E-01 -1.6540E+00
             1.4480E+00
 GRADIENT:   6.5749E+01 -1.1391E+01 -3.9477E+01 -3.4093E+01  5.4373E+01 -1.6041E+01 -7.3029E-02  0.0000E+00  2.8008E+00 -3.7985E+00
            -3.0781E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1446.63929734993        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      764
 NPARAMETR:  1.0156E+00  2.7824E-01  2.7451E-01  1.2042E+00  2.6725E-01  8.8507E-01  7.9886E-02  1.0000E-02  1.4809E+00  4.8784E-01
             3.8394E+00
 PARAMETER:  1.1548E-01 -1.1793E+00 -1.1928E+00  2.8583E-01 -1.2196E+00 -2.2089E-02 -2.4272E+00 -1.1720E+01  4.9268E-01 -6.1777E-01
             1.4453E+00
 GRADIENT:  -2.8063E+01  6.8821E+00  5.3747E+00  2.8864E+01 -2.0552E+01 -5.1463E+00  7.9241E-05  0.0000E+00  3.7091E+00  5.2931E+00
             1.3586E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1447.97233291627        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      942
 NPARAMETR:  1.0276E+00  2.3765E-01  2.7549E-01  1.1784E+00  2.6314E-01  8.9486E-01  6.2992E-02  1.0000E-02  1.4315E+00  4.4257E-01
             3.8025E+00
 PARAMETER:  1.2725E-01 -1.3370E+00 -1.1892E+00  2.6412E-01 -1.2351E+00 -1.1093E-02 -2.6647E+00 -1.2334E+01  4.5870E-01 -7.1517E-01
             1.4357E+00
 GRADIENT:   1.1164E+01  2.9365E+00  4.9407E+00  6.7108E-02 -1.5992E+01 -6.8534E-01 -1.2013E-03  0.0000E+00 -6.0047E-02  1.5668E+00
            -3.9697E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1449.54643186977        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1118
 NPARAMETR:  1.0089E+00  1.0620E-01  3.1107E-01  1.2612E+00  2.7308E-01  8.8676E-01  1.1870E-02  1.0000E-02  1.3106E+00  3.1076E-01
             3.8830E+00
 PARAMETER:  1.0891E-01 -2.1424E+00 -1.0677E+00  3.3209E-01 -1.1980E+00 -2.0184E-02 -4.3338E+00 -1.6845E+01  3.7045E-01 -1.0687E+00
             1.4566E+00
 GRADIENT:  -1.1259E+01  1.4666E+00  1.0211E+01 -3.9204E-01 -1.6124E+01 -4.9572E-01 -9.7427E-06  0.0000E+00 -7.6730E-01 -1.4425E+00
             1.5133E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1450.03781831453        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1293
 NPARAMETR:  1.0081E+00  3.9954E-02  3.1508E-01  1.2905E+00  2.7130E-01  8.8440E-01  1.0000E-02  1.0000E-02  1.2708E+00  3.6995E-01
             3.8635E+00
 PARAMETER:  1.0811E-01 -3.1200E+00 -1.0549E+00  3.5502E-01 -1.2045E+00 -2.2851E-02 -6.4759E+00 -2.2403E+01  3.3963E-01 -8.9440E-01
             1.4516E+00
 GRADIENT:  -1.6075E+00  1.0370E-01 -3.3101E-01 -3.6780E-01  1.2268E+00 -4.9750E-01  0.0000E+00  0.0000E+00 -8.3362E-01  4.5205E-01
             1.7936E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1450.08125471361        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1469
 NPARAMETR:  1.0070E+00  1.4449E-02  3.1659E-01  1.3010E+00  2.7043E-01  8.8482E-01  1.0000E-02  1.0000E-02  1.2611E+00  3.5803E-01
             3.8642E+00
 PARAMETER:  1.0693E-01 -4.1371E+00 -1.0501E+00  3.6314E-01 -1.2077E+00 -2.2375E-02 -8.7165E+00 -2.8205E+01  3.3197E-01 -9.2713E-01
             1.4518E+00
 GRADIENT:   6.3594E-01  2.6940E-03 -6.5324E-01 -4.2491E-01  1.6528E+00  4.3865E-02  0.0000E+00  0.0000E+00  4.7724E-02  4.3909E-02
             1.0522E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1450.08301299886        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1649
 NPARAMETR:  1.0065E+00  1.0000E-02  3.1650E-01  1.3023E+00  2.7005E-01  8.8459E-01  1.0000E-02  1.0000E-02  1.2594E+00  3.5770E-01
             3.8644E+00
 PARAMETER:  1.0645E-01 -4.5401E+00 -1.0504E+00  3.6411E-01 -1.2091E+00 -2.2637E-02 -9.5063E+00 -3.0250E+01  3.3067E-01 -9.2806E-01
             1.4518E+00
 GRADIENT:   2.5129E-01  1.5413E-04 -6.7999E-01 -5.0281E-01  1.4145E+00  2.3102E-02  0.0000E+00  0.0000E+00  9.0354E-02  3.0819E-02
             8.1990E-02

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1450.08301299886        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1671
 NPARAMETR:  1.0065E+00  1.0000E-02  3.1650E-01  1.3023E+00  2.7005E-01  8.8459E-01  1.0000E-02  1.0000E-02  1.2594E+00  3.5770E-01
             3.8644E+00
 PARAMETER:  1.0645E-01 -4.5401E+00 -1.0504E+00  3.6411E-01 -1.2091E+00 -2.2637E-02 -9.5063E+00 -3.0250E+01  3.3067E-01 -9.2806E-01
             1.4518E+00
 GRADIENT:   2.5129E-01  1.5413E-04 -6.7999E-01 -5.0281E-01  1.4145E+00  2.3102E-02  0.0000E+00  0.0000E+00  9.0354E-02  3.0819E-02
             8.1990E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1671
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5095E-04 -2.7132E-06  1.9954E-04 -1.1717E-02  2.9475E-03
 SE:             2.8284E-02  1.4054E-06  2.0969E-04  2.7072E-02  1.1753E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9292E-01  5.3532E-02  3.4131E-01  6.6516E-01  8.0199E-01

 ETASHRINKSD(%)  5.2443E+00  9.9995E+01  9.9297E+01  9.3043E+00  6.0625E+01
 ETASHRINKVR(%)  1.0214E+01  1.0000E+02  9.9995E+01  1.7743E+01  8.4496E+01
 EBVSHRINKSD(%)  5.0494E+00  9.9996E+01  9.9255E+01  7.5525E+00  6.0930E+01
 EBVSHRINKVR(%)  9.8438E+00  1.0000E+02  9.9994E+01  1.4535E+01  8.4736E+01
 RELATIVEINF(%)  7.2509E+01  2.8602E-08  3.1270E-04  2.9341E+01  5.7404E-01
 EPSSHRINKSD(%)  1.9224E+01
 EPSSHRINKVR(%)  3.4753E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1450.0830129988581     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -531.14447979418537     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1450.083       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.00E-02  3.16E-01  1.30E+00  2.70E-01  8.85E-01  1.00E-02  1.00E-02  1.26E+00  3.58E-01  3.86E+00
 


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
+       -1.56E+01 -6.81E+00
 
 TH 3
+       -9.92E+01  7.78E+01  5.71E+03
 
 TH 4
+       -3.25E+01  1.90E+01 -1.17E+02  3.48E+02
 
 TH 5
+        2.77E+02 -1.90E+02 -9.26E+03 -5.22E+02  1.91E+04
 
 TH 6
+       -1.06E+00 -1.38E+00  3.12E+01 -1.01E+01 -1.73E+00  2.03E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.44E+01 -5.58E+00  5.05E+01 -2.24E+00  1.08E+02  3.63E+00  0.00E+00  0.00E+00  7.97E+01
 
 TH10
+       -1.39E+01 -1.51E+00 -8.92E+01 -3.42E+00  2.07E+02  5.09E-01  0.00E+00  0.00E+00  1.69E-01  4.98E+01
 
 TH11
+       -2.40E+01 -2.61E-01 -1.20E+01 -6.44E+00  2.88E+01  3.68E+00  0.00E+00  0.00E+00  4.45E+00  1.57E+01  3.83E+01
 
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
 #CPUT: Total CPU Time in Seconds,       37.708
Stop Time:
Thu Sep 30 00:43:19 CDT 2021
