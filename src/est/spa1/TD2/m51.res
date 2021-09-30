Thu Sep 30 02:10:50 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat51.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1603.93177336429        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1444E+02  2.0065E+01  3.6561E+01  7.5725E+01  9.4662E+01  6.1114E+01  8.1404E+00 -2.1119E+02 -3.2150E+01  1.2215E+01
            -7.1203E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2019.98689708727        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:      122
 NPARAMETR:  8.2432E-01  9.9814E-01  9.7568E-01  9.6135E-01  9.3191E-01  1.0987E+00  9.8882E-01  1.1879E+00  1.0257E+00  9.7253E-01
             1.5411E+00
 PARAMETER: -9.3196E-02  9.8138E-02  7.5383E-02  6.0581E-02  2.9480E-02  1.9417E-01  8.8757E-02  2.7215E-01  1.2542E-01  7.2146E-02
             5.3250E-01
 GRADIENT:  -2.7811E+02 -2.6208E+01 -1.2792E+01 -3.6033E+01 -1.1795E+01  1.0354E+01  7.6347E+00  9.9729E+00  1.0049E+00  3.0446E+01
             2.7659E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2047.46816948539        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.7325E-01  9.9774E-01  1.0382E+00  9.6096E-01  9.3229E-01  9.5623E-01  5.5763E-01  1.1865E+00  1.0917E+00  9.7214E-01
             1.5443E+00
 PARAMETER:  7.2890E-02  9.7734E-02  1.3746E-01  6.0177E-02  2.9885E-02  5.5248E-02 -4.8406E-01  2.7104E-01  1.8774E-01  7.1741E-02
             5.3459E-01
 GRADIENT:   2.1660E+02  6.1363E+00  3.4686E+00 -5.8776E+00 -3.1766E+01  3.0397E+01  4.3003E+00  5.2245E+00 -2.6100E+00  2.5730E+01
             2.7489E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2048.36277430802        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      362
 NPARAMETR:  9.5344E-01  9.9772E-01  1.0337E+00  9.6095E-01  9.3230E-01  9.3357E-01  4.2716E-01  1.1865E+00  1.1438E+00  9.7212E-01
             1.5443E+00
 PARAMETER:  5.2317E-02  9.7721E-02  1.3313E-01  6.0164E-02  2.9898E-02  3.1259E-02 -7.5060E-01  2.7101E-01  2.3433E-01  7.1729E-02
             5.3460E-01
 GRADIENT:  -1.9564E+01 -1.2558E+01  1.9278E-01 -2.8037E+01 -3.3050E+01 -3.3894E+00  1.5498E+00  4.5974E+00 -3.4427E+00  2.4537E+01
             2.7188E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2049.13606544916        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      540
 NPARAMETR:  9.6540E-01  9.9773E-01  1.0231E+00  9.6095E-01  9.3230E-01  9.4796E-01  2.2369E-01  1.1865E+00  1.2138E+00  9.7213E-01
             1.5438E+00
 PARAMETER:  6.4788E-02  9.7725E-02  1.2287E-01  6.0169E-02  2.9897E-02  4.6557E-02 -1.3975E+00  2.7102E-01  2.9374E-01  7.1730E-02
             5.3425E-01
 GRADIENT:   1.1430E+01 -1.6621E+01 -3.4880E+00 -1.8835E+01 -2.7492E+01  2.4705E+00  7.2278E-01  4.3827E+00  7.0196E+00  2.3761E+01
             2.7206E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2049.20804961709        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  9.6526E-01  9.9770E-01  1.0235E+00  9.6098E-01  9.3233E-01  9.3179E-01  2.2068E-01  1.1494E+00  1.2135E+00  9.7209E-01
             1.5440E+00
 PARAMETER:  6.4644E-02  9.7693E-02  1.2320E-01  6.0201E-02  2.9930E-02  2.9352E-02 -1.4110E+00  2.3927E-01  2.9354E-01  7.1698E-02
             5.3441E-01
 GRADIENT:   1.9765E+02  6.6794E+00 -7.3164E-01  1.0919E+01 -2.4413E+01  2.1162E+01  2.0033E+00  3.1613E+00  1.9508E+01  2.3274E+01
             2.7550E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2049.46394081416        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  9.6523E-01  9.9767E-01  1.0234E+00  9.6101E-01  9.3236E-01  9.4822E-01  2.2076E-01  9.9261E-01  1.2136E+00  9.7206E-01
             1.5426E+00
 PARAMETER:  6.4616E-02  9.7670E-02  1.2317E-01  6.0231E-02  2.9963E-02  4.6833E-02 -1.4107E+00  9.2586E-02  2.9360E-01  7.1667E-02
             5.3350E-01
 GRADIENT:   1.0991E+01 -9.1757E+00  6.4889E+00 -1.8393E+01 -2.8935E+01  2.5637E+00  5.8774E-01 -2.7400E+00  7.0167E+00  1.9409E+01
             2.6909E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2084.13508369081        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1022
 NPARAMETR:  9.5375E-01  9.8781E-01  1.0099E+00  9.7328E-01  9.4531E-01  9.4367E-01  2.5457E-01  2.4846E-02  1.2477E+00  9.5950E-01
             1.0592E+00
 PARAMETER:  5.2643E-02  8.7735E-02  1.0990E-01  7.2914E-02  4.3761E-02  4.2020E-02 -1.2682E+00 -3.5950E+00  3.2133E-01  5.8657E-02
             1.5755E-01
 GRADIENT:  -6.8001E+00  2.3841E+01  1.0336E+02  4.8467E+00 -1.4021E+01 -1.2778E-01 -1.0657E+00 -1.6016E-01  1.6523E+01 -3.8288E+01
            -1.2375E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2098.52995498124        NO. OF FUNC. EVALS.: 113
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  9.0540E-01  9.8176E-01  7.3934E-01  9.3443E-01  9.1760E-01  8.9303E-01  1.0370E+00  3.7174E-02  1.0707E+00  1.1190E+00
             1.0025E+00
 PARAMETER:  6.2620E-04  8.1587E-02 -2.0200E-01  3.2187E-02  1.4008E-02 -1.3135E-02  1.3635E-01 -3.1922E+00  1.6832E-01  2.1242E-01
             1.0252E-01
 GRADIENT:  -1.5361E+02 -1.0120E+02 -3.9000E+01 -1.3267E+01  1.2250E+02 -3.2318E+01  1.3303E-01 -2.0581E-01  7.5840E+00  3.4865E+01
             5.0241E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2128.28422528644        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1314
 NPARAMETR:  9.6724E-01  1.1175E+00  6.0928E-01  8.9411E-01  7.8518E-01  9.5529E-01  8.6921E-01  4.3853E-01  1.0304E+00  8.0586E-01
             8.9407E-01
 PARAMETER:  6.6691E-02  2.1114E-01 -3.9548E-01 -1.1925E-02 -1.4184E-01  5.4256E-02 -4.0165E-02 -7.2433E-01  1.2997E-01 -1.1585E-01
            -1.1971E-02
 GRADIENT:   2.2167E+01  2.8505E+01  6.2002E+00  2.6174E+01 -3.6386E+01  2.6199E+00 -8.5711E+00 -7.8539E-01 -2.1283E+01  1.2294E+01
            -2.0812E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2131.22099322778        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1497             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6173E-01  1.1108E+00  6.2389E-01  8.8068E-01  8.0436E-01  9.5192E-01  9.3689E-01  4.8717E-01  1.1019E+00  7.3451E-01
             9.1771E-01
 PARAMETER:  6.0978E-02  2.0511E-01 -3.7177E-01 -2.7056E-02 -1.1771E-01  5.0731E-02  3.4815E-02 -6.1915E-01  1.9707E-01 -2.0855E-01
             1.4127E-02
 GRADIENT:   5.1981E+02  1.2704E+02  1.1292E+01  8.4676E+01  1.0266E+01  7.0172E+01  3.8879E+00 -2.5879E-01  1.6970E+01  4.1072E+00
             1.6562E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2131.40905230146        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1644
 NPARAMETR:  9.5876E-01  1.1234E+00  6.2390E-01  8.7094E-01  8.1004E-01  9.4893E-01  9.5235E-01  4.8716E-01  1.1154E+00  7.1320E-01
             9.1768E-01
 PARAMETER:  5.7888E-02  2.1637E-01 -3.7177E-01 -3.8181E-02 -1.1067E-01  4.7584E-02  5.1179E-02 -6.1916E-01  2.0924E-01 -2.3799E-01
             1.4097E-02
 GRADIENT:   9.3005E-01 -4.2156E+00  3.3970E+00  1.9801E+00 -5.2227E+00  3.3975E-01  8.5836E-01 -9.7789E-01 -2.8511E-01 -6.4940E-01
             2.3613E-01

0ITERATION NO.:   58    OBJECTIVE VALUE:  -2131.46553204577        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1744
 NPARAMETR:  9.5894E-01  1.1329E+00  6.2455E-01  8.6698E-01  8.1352E-01  9.4868E-01  9.4926E-01  4.8631E-01  1.1225E+00  7.1224E-01
             9.1734E-01
 PARAMETER:  5.7795E-02  2.2259E-01 -3.7177E-01 -4.2736E-02 -1.0639E-01  4.7600E-02  4.7927E-02 -6.1915E-01  2.1493E-01 -2.3934E-01
             1.3724E-02
 GRADIENT:  -1.5012E+05 -2.6700E+00 -8.0695E+04 -5.5863E+01  1.4597E+01  3.0026E+05  2.1585E-01  4.8414E+04 -1.3978E+05 -6.8369E+00
            -8.8117E+02
 NUMSIGDIG:         2.3         1.8         2.3         6.3         6.9         2.3         8.7         2.3         2.3         6.9
                    5.1

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1744
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.6668E-04 -1.1958E-02 -1.5886E-02  6.6785E-03 -1.8850E-02
 SE:             2.9890E-02  2.1575E-02  8.7005E-03  2.6069E-02  2.1610E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7687E-01  5.7942E-01  6.7879E-02  7.9781E-01  3.8307E-01

 ETASHRINKSD(%)  1.0000E-10  2.7722E+01  7.0852E+01  1.2666E+01  2.7603E+01
 ETASHRINKVR(%)  1.0000E-10  4.7759E+01  9.1504E+01  2.3727E+01  4.7586E+01
 EBVSHRINKSD(%)  3.0885E-01  2.7417E+01  7.3784E+01  1.2368E+01  2.8220E+01
 EBVSHRINKVR(%)  6.1675E-01  4.7316E+01  9.3127E+01  2.3206E+01  4.8476E+01
 RELATIVEINF(%)  9.9265E+01  4.8802E+00  1.3143E+00  1.0324E+01  7.3705E+00
 EPSSHRINKSD(%)  3.4519E+01
 EPSSHRINKVR(%)  5.7123E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2131.4655320457673     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1212.5269988410946     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2131.466       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.59E-01  1.13E+00  6.24E-01  8.67E-01  8.14E-01  9.49E-01  9.49E-01  4.87E-01  1.12E+00  7.12E-01  9.17E-01
 


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
+        8.17E+07
 
 TH 2
+        3.11E+07  2.37E+07
 
 TH 3
+       -2.26E+03  1.02E+03  1.39E+07
 
 TH 4
+        1.13E+00  3.02E+02  6.90E+03  9.98E+07
 
 TH 5
+        9.05E+07 -5.78E+02 -3.51E+03  3.70E+02  2.00E+08
 
 TH 6
+       -5.01E+03  3.14E+07 -2.07E+03 -8.88E-01 -2.75E-01  1.67E+08
 
 TH 7
+        8.25E+07  3.38E+01  6.64E+01 -5.42E+00 -9.14E+07  8.33E+07  1.67E+08
 
 TH 8
+        1.72E+03 -5.24E+02  1.08E+07 -5.54E+03  1.96E+03  1.59E+03 -5.77E+01  8.22E+06
 
 TH 9
+        2.72E+07 -2.94E+01  1.58E+04  7.04E+01 -2.23E+01 -1.99E+03  1.31E+01 -1.21E+04  2.58E+07
 
 TH10
+        4.59E+07 -1.75E+07  7.91E+02 -1.68E+01 -5.09E+07  7.32E-02 -4.64E+07 -6.33E+02 -1.83E+07  2.58E+07
 
 TH11
+       -8.01E+00  3.25E+07  1.03E+05 -9.42E+07 -9.47E+07  1.22E+00  7.35E+00 -7.96E+04 -9.96E+04  4.81E+07  1.78E+08
 
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
 #CPUT: Total CPU Time in Seconds,       31.768
Stop Time:
Thu Sep 30 02:11:24 CDT 2021
