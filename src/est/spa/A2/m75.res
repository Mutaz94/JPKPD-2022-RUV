Sat Sep 25 08:51:04 CDT 2021
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
$DATA ../../../../data/spa/A2/dat75.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1240.45352290260        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0449E+01 -5.9082E+01  1.0547E+01 -1.3445E+02  6.0615E+01 -4.0714E+00 -2.7610E+01  5.7313E-01 -8.8163E+01 -2.1381E+01
            -7.0882E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1450.64874067486        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0584E+00  8.5016E-01  1.0725E+00  1.1408E+00  9.8112E-01  9.6317E-01  9.8606E-01  8.3077E-01  1.2202E+00  8.4506E-01
             1.9326E+00
 PARAMETER:  1.5676E-01 -6.2330E-02  1.6996E-01  2.3177E-01  8.0938E-02  6.2472E-02  8.5958E-02 -8.5397E-02  2.9901E-01 -6.8346E-02
             7.5885E-01
 GRADIENT:   1.3531E+02 -4.6062E+01 -2.7622E+01 -4.9517E+01  7.5500E+01 -2.2427E+01  1.4036E+00  2.1771E+00  8.0006E+00 -1.3082E+01
            -4.3445E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1456.94645225325        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0333E+00  6.6698E-01  8.0811E-01  1.2657E+00  7.3131E-01  9.5285E-01  6.6358E-01  1.5271E-01  1.2330E+00  9.5884E-01
             1.9528E+00
 PARAMETER:  1.3279E-01 -3.0500E-01 -1.1306E-01  3.3561E-01 -2.1292E-01  5.1701E-02 -3.1010E-01 -1.7792E+00  3.0947E-01  5.7971E-02
             7.6925E-01
 GRADIENT:   7.2611E+01 -1.0832E+01 -3.1566E+01 -3.3266E+00  2.5920E+01 -2.3116E+01  3.5383E-01  2.2653E-01  2.6126E+01  1.3126E+01
            -1.7279E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1464.60580675526        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.9706E-01  5.2048E-01  8.3186E-01  1.3674E+00  6.8718E-01  9.9019E-01  1.8531E+00  1.3247E-01  9.7586E-01  7.4704E-01
             2.0737E+00
 PARAMETER:  9.7051E-02 -5.5301E-01 -8.4095E-02  4.1293E-01 -2.7516E-01  9.0141E-02  7.1684E-01 -1.9214E+00  7.5566E-02 -1.9164E-01
             8.2931E-01
 GRADIENT:  -1.3691E+01  6.6105E+00  2.5339E+00  9.4567E+00  3.7510E+00 -4.5999E+00  3.7834E+00  1.2510E-01  3.4169E-01  1.2368E+00
             8.4116E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1465.67088440401        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0058E+00  4.8955E-01  5.8242E-01  1.3299E+00  5.2579E-01  1.0079E+00  1.9217E+00  6.3375E-02  9.5206E-01  6.1341E-01
             2.0043E+00
 PARAMETER:  1.0578E-01 -6.1428E-01 -4.4057E-01  3.8510E-01 -5.4285E-01  1.0788E-01  7.5322E-01 -2.6587E+00  5.0872E-02 -3.8873E-01
             7.9527E-01
 GRADIENT:  -3.5758E-01  6.8594E+00  4.1707E+00  7.8233E+00 -8.3663E+00  2.0101E-01  6.7194E-01  7.1654E-02  2.1542E+00  3.3553E+00
             4.2322E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1465.67270580156        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0059E+00  4.8217E-01  5.6117E-01  1.3265E+00  5.1097E-01  1.0085E+00  1.9398E+00  5.7958E-02  9.4605E-01  5.9485E-01
             1.9984E+00
 PARAMETER:  1.0587E-01 -6.2946E-01 -4.7774E-01  3.8253E-01 -5.7144E-01  1.0843E-01  7.6258E-01 -2.7480E+00  4.4535E-02 -4.1945E-01
             7.9233E-01
 GRADIENT:  -6.5647E-01  6.3900E+00  3.9056E+00  7.4865E+00 -7.8097E+00  2.2473E-01  5.6964E-01  6.1002E-02  1.7248E+00  2.6930E+00
             3.4062E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1466.18987561812        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  1.0073E+00  4.2651E-01  6.3085E-01  1.3781E+00  5.4492E-01  1.0068E+00  2.1547E+00  5.4470E-02  9.2423E-01  6.0652E-01
             2.0080E+00
 PARAMETER:  1.0726E-01 -7.5213E-01 -3.6068E-01  4.2068E-01 -5.0712E-01  1.0681E-01  8.6767E-01 -2.8101E+00  2.1200E-02 -4.0001E-01
             7.9715E-01
 GRADIENT:  -2.7011E+00  1.7369E+00 -3.2587E+00  5.6828E+00  5.4349E+00 -4.9045E-01  1.1189E-01  2.8404E-02 -2.2965E+00 -2.1744E+00
            -1.6589E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1467.73713863663        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  1.0060E+00  2.4088E-01  7.1286E-01  1.4882E+00  5.4171E-01  1.0200E+00  2.5130E+00  1.0000E-02  9.7136E-01  7.1591E-01
             1.9809E+00
 PARAMETER:  1.0601E-01 -1.3234E+00 -2.3847E-01  4.9757E-01 -5.1302E-01  1.1980E-01  1.0215E+00 -4.6975E+00  7.0937E-02 -2.3421E-01
             7.8357E-01
 GRADIENT:   8.2317E+00  3.2181E+00  3.4703E+01 -4.6526E+00 -4.1713E+01  6.3390E+00 -3.4830E+00  0.0000E+00  1.2281E+01 -1.9805E+00
            -9.4480E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1469.90266475634        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  9.9965E-01  1.3248E-01  6.8226E-01  1.5372E+00  5.1537E-01  9.9031E-01  3.6203E+00  1.0000E-02  8.9822E-01  7.0612E-01
             2.0196E+00
 PARAMETER:  9.9652E-02 -1.9213E+00 -2.8235E-01  5.2995E-01 -5.6287E-01  9.0267E-02  1.3866E+00 -7.2547E+00 -7.3401E-03 -2.4797E-01
             8.0289E-01
 GRADIENT:   8.9645E-02 -1.0119E+00  8.3477E+00  5.4593E+00 -8.6633E+00 -3.7215E+00 -5.7720E+00  0.0000E+00 -3.5088E+00  2.7574E+00
             3.6956E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1470.85202729189        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1085
 NPARAMETR:  9.9874E-01  7.0252E-02  5.6059E-01  1.5411E+00  4.3618E-01  9.9791E-01  4.8653E+00  1.0000E-02  9.1360E-01  6.2429E-01
             1.9820E+00
 PARAMETER:  9.8740E-02 -2.5557E+00 -4.7877E-01  5.3247E-01 -7.2971E-01  9.7905E-02  1.6821E+00 -1.0115E+01  9.6375E-03 -3.7114E-01
             7.8408E-01
 GRADIENT:   5.4762E-01 -2.8045E+00  1.4198E+01  4.1719E+01 -2.6545E+01 -1.3723E+00 -1.1406E+01  0.0000E+00  1.5357E+00 -1.8622E+00
            -2.2276E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1471.82973812678        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1265
 NPARAMETR:  9.9568E-01  3.6617E-02  4.9117E-01  1.4947E+00  3.9374E-01  1.0097E+00  6.6043E+00  1.0000E-02  9.1826E-01  6.2242E-01
             1.9841E+00
 PARAMETER:  9.5671E-02 -3.2072E+00 -6.1097E-01  5.0192E-01 -8.3206E-01  1.0963E-01  1.9877E+00 -1.3490E+01  1.4724E-02 -3.7415E-01
             7.8517E-01
 GRADIENT:  -3.4008E+00 -3.2098E+00  4.0330E+00  1.1890E+00 -1.2863E+01  3.1008E+00 -9.6274E+00  0.0000E+00  3.0051E+00  2.4911E+00
             6.5272E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1473.20721502910        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1446
 NPARAMETR:  9.9551E-01  1.1454E-02  5.4334E-01  1.5292E+00  4.2386E-01  9.9715E-01  1.1792E+01  1.0000E-02  8.9581E-01  6.1626E-01
             1.9791E+00
 PARAMETER:  9.5503E-02 -4.3694E+00 -5.1001E-01  5.2474E-01 -7.5836E-01  9.7143E-02  2.5674E+00 -1.8263E+01 -1.0028E-02 -3.8409E-01
             7.8266E-01
 GRADIENT:   4.9775E-01 -1.0718E+00 -1.1094E+00  9.9901E-02  6.0952E-01 -8.7711E-01 -2.3130E+00  0.0000E+00 -1.9439E+00 -7.9925E-01
            -6.0803E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1473.24076535928        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1621
 NPARAMETR:  9.9509E-01  1.0823E-02  5.4893E-01  1.5326E+00  4.2703E-01  9.9983E-01  1.2234E+01  1.0000E-02  9.0407E-01  6.2165E-01
             1.9811E+00
 PARAMETER:  9.5074E-02 -4.4261E+00 -4.9979E-01  5.2699E-01 -7.5091E-01  9.9834E-02  2.6043E+00 -1.8501E+01 -8.5229E-04 -3.7538E-01
             7.8367E-01
 GRADIENT:   1.5016E-01  1.3602E+00 -2.3933E+00 -2.0625E+00  2.6594E+00 -1.0764E-01  2.8226E+00  0.0000E+00 -1.1072E+00 -9.1022E-01
            -5.3805E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1473.24082360561        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1684
 NPARAMETR:  9.9500E-01  1.0791E-02  5.4929E-01  1.5335E+00  4.2675E-01  9.9993E-01  1.2194E+01  1.0000E-02  9.0402E-01  6.2192E-01
             1.9828E+00
 PARAMETER:  9.5087E-02 -4.4246E+00 -4.9963E-01  5.2701E-01 -7.5080E-01  9.9833E-02  2.6035E+00 -1.8495E+01 -1.0042E-03 -3.7531E-01
             7.8372E-01
 GRADIENT:   3.9248E+03  1.7563E+02 -7.8442E+02 -1.4856E+03  1.0408E+03 -3.9248E+03  2.9808E+02  0.0000E+00 -3.9241E+03 -1.0458E+03
            -1.0006E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1684
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3913E-04  1.2730E-02 -7.3144E-06 -9.0202E-03 -8.2627E-03
 SE:             2.9425E-02  6.8715E-03  2.4974E-04  2.7908E-02  2.1327E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9080E-01  6.3941E-02  9.7663E-01  7.4653E-01  6.9844E-01

 ETASHRINKSD(%)  1.4222E+00  7.6980E+01  9.9163E+01  6.5052E+00  2.8550E+01
 ETASHRINKVR(%)  2.8242E+00  9.4701E+01  9.9993E+01  1.2587E+01  4.8949E+01
 EBVSHRINKSD(%)  1.4996E+00  8.1509E+01  9.9115E+01  5.8858E+00  2.7944E+01
 EBVSHRINKVR(%)  2.9766E+00  9.6581E+01  9.9992E+01  1.1425E+01  4.8079E+01
 RELATIVEINF(%)  9.6769E+01  2.9777E+00  3.2471E-04  4.6885E+01  2.1723E+00
 EPSSHRINKSD(%)  3.5598E+01
 EPSSHRINKVR(%)  5.8524E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1473.2408236056112     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -738.08999704187306     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.95
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1473.241       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  1.08E-02  5.49E-01  1.53E+00  4.27E-01  1.00E+00  1.22E+01  1.00E-02  9.04E-01  6.22E-01  1.98E+00
 


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
+        1.98E+07
 
 TH 2
+       -2.23E+03  8.49E+07
 
 TH 3
+        4.67E+02  7.17E+04  2.60E+06
 
 TH 4
+        1.33E+02  1.96E+04 -4.61E+03  3.00E+05
 
 TH 5
+       -3.18E+02 -1.22E+05  6.53E+03  3.18E+03  1.90E+06
 
 TH 6
+       -1.97E+07  1.47E+03 -3.36E+02 -1.11E+02  2.69E+02  1.96E+07
 
 TH 7
+       -3.21E+00 -5.34E+02  1.07E+02  2.92E+01 -1.83E+02  2.23E+00  1.93E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.93E+03  1.79E+04 -3.82E+03 -1.23E+03  3.04E+03 -1.32E+03  2.68E+01  0.00E+00  2.40E+07
 
 TH10
+        6.79E+02  1.20E+04 -2.56E+03 -8.08E+02  2.02E+03 -4.76E+02  1.80E+01  0.00E+00 -5.29E+03  3.60E+06
 
 TH11
+        7.61E+01  2.72E+03 -6.90E+02 -2.04E+02  4.87E+02 -5.82E+01  4.06E+00  0.00E+00 -6.73E+02 -4.12E+02  8.13E+04
 
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
 #CPUT: Total CPU Time in Seconds,       27.608
Stop Time:
Sat Sep 25 08:51:33 CDT 2021
