Sat Sep 18 11:45:42 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat56.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1653.53066929372        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5044E+01 -6.1630E+01 -1.0224E+01 -8.8744E+01 -2.7028E+01 -4.1310E+00 -1.8515E+01  1.2479E+01 -1.1029E+01  8.1571E+00
            -1.4429E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1660.63605329438        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0243E+00  1.1481E+00  1.2000E+00  9.8312E-01  1.2260E+00  1.0110E+00  1.2801E+00  8.5554E-01  1.0372E+00  9.9220E-01
             1.0753E+00
 PARAMETER:  1.2400E-01  2.3811E-01  2.8233E-01  8.2978E-02  3.0379E-01  1.1091E-01  3.4697E-01 -5.6018E-02  1.3652E-01  9.2173E-02
             1.7258E-01
 GRADIENT:   8.0493E+01  6.7118E+00  1.3821E+00  1.4299E+01  4.9089E+01 -2.3129E-01  6.7100E+00 -4.2162E+00  9.0324E+00 -1.9634E+01
             8.4594E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1663.60065430847        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0149E+00  9.5375E-01  9.3741E-01  1.0852E+00  9.8869E-01  1.0231E+00  1.5786E+00  4.1563E-01  8.7438E-01  8.5824E-01
             1.0564E+00
 PARAMETER:  1.1476E-01  5.2650E-02  3.5370E-02  1.8176E-01  8.8623E-02  1.2287E-01  5.5651E-01 -7.7796E-01 -3.4236E-02 -5.2871E-02
             1.5488E-01
 GRADIENT:   5.6141E+01  1.5811E+00 -1.6523E+01  2.8206E+01  4.7589E+01  5.4576E+00  1.3263E+01 -5.1957E-01  2.7148E+00 -1.1297E+01
             5.3296E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.60193646676        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.9504E-01  1.0188E+00  8.5358E-01  1.0322E+00  9.5161E-01  1.0110E+00  1.3979E+00  3.6123E-01  9.0464E-01  8.6847E-01
             1.0329E+00
 PARAMETER:  9.5024E-02  1.1858E-01 -5.8312E-02  1.3166E-01  5.0398E-02  1.1091E-01  4.3500E-01 -9.1824E-01 -2.2315E-04 -4.1023E-02
             1.3241E-01
 GRADIENT:   8.0847E+00 -1.5365E+00 -7.6741E+00  5.6255E+00  1.0254E+01 -1.7155E-01  2.9752E+00  7.8977E-01  2.2314E+00 -1.1143E+00
            -7.1575E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.60300993156        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.9393E-01  1.0194E+00  8.5094E-01  1.0313E+00  9.4893E-01  1.0111E+00  1.3915E+00  3.4771E-01  9.0242E-01  8.6963E-01
             1.0335E+00
 PARAMETER:  9.3914E-02  1.1917E-01 -6.1409E-02  1.3078E-01  4.7583E-02  1.1106E-01  4.3039E-01 -9.5640E-01 -2.6793E-03 -3.9687E-02
             1.3296E-01
 GRADIENT:   5.5943E+00 -1.2157E+00 -5.8901E+00  4.1127E+00  7.3832E+00 -1.5768E-01  2.1908E+00  6.9625E-01  1.6201E+00 -7.3218E-01
            -4.8988E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.60384945397        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.9315E-01  1.0204E+00  8.4676E-01  1.0300E+00  9.4592E-01  1.0113E+00  1.3859E+00  3.2819E-01  9.0099E-01  8.6970E-01
             1.0339E+00
 PARAMETER:  9.3128E-02  1.2021E-01 -6.6339E-02  1.2959E-01  4.4402E-02  1.1121E-01  4.2636E-01 -1.0142E+00 -4.2646E-03 -3.9609E-02
             1.3338E-01
 GRADIENT:   3.7784E+00 -9.4182E-01 -4.4318E+00  2.9402E+00  5.2018E+00 -1.3556E-01  1.5801E+00  5.8716E-01  1.1535E+00 -4.5823E-01
            -3.2298E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1664.60442291096        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  9.9259E-01  1.0215E+00  8.4191E-01  1.0288E+00  9.4278E-01  1.0114E+00  1.3813E+00  3.0608E-01  8.9998E-01  8.6904E-01
             1.0343E+00
 PARAMETER:  9.2558E-02  1.2131E-01 -7.2087E-02  1.2836E-01  4.1082E-02  1.1136E-01  4.2301E-01 -1.0839E+00 -5.3777E-03 -4.0366E-02
             1.3371E-01
 GRADIENT:   2.4207E+00 -7.1781E-01 -3.2683E+00  2.0340E+00  3.5328E+00 -1.1261E-01  1.1034E+00  4.8504E-01  7.9427E-01 -2.5593E-01
            -1.9678E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1664.60487466165        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      510
 NPARAMETR:  9.9207E-01  1.0228E+00  8.3550E-01  1.0272E+00  9.3884E-01  1.0116E+00  1.3765E+00  2.7667E-01  8.9905E-01  8.6760E-01
             1.0346E+00
 PARAMETER:  9.2036E-02  1.2259E-01 -7.9725E-02  1.2688E-01  3.6891E-02  1.1153E-01  4.1957E-01 -1.1849E+00 -6.4126E-03 -4.2020E-02
             1.3402E-01
 GRADIENT:   1.1313E+00 -4.9551E-01 -2.1254E+00  1.1534E+00  1.9251E+00 -8.6373E-02  6.3735E-01  3.7481E-01  4.4712E-01 -6.2051E-02
            -7.3829E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1664.60574557392        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      580
 NPARAMETR:  9.9160E-01  1.0243E+00  8.2743E-01  1.0255E+00  9.3400E-01  1.0118E+00  1.3718E+00  2.3744E-01  8.9817E-01  8.6527E-01
             1.0349E+00
 PARAMETER:  9.1569E-02  1.2400E-01 -8.9431E-02  1.2517E-01  3.1718E-02  1.1172E-01  4.1612E-01 -1.3379E+00 -7.4014E-03 -4.4717E-02
             1.3431E-01
 GRADIENT:  -7.5705E-02 -2.8461E-01 -1.0319E+00  3.1563E-01  4.0870E-01 -5.7941E-02  1.8856E-01  2.6074E-01  1.1688E-01  1.2219E-01
             4.4710E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1664.60903183705        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  9.9120E-01  1.0260E+00  8.1670E-01  1.0233E+00  9.2770E-01  1.0120E+00  1.3668E+00  1.7763E-01  8.9732E-01  8.6155E-01
             1.0352E+00
 PARAMETER:  9.1166E-02  1.2569E-01 -1.0248E-01  1.2305E-01  2.4956E-02  1.1193E-01  4.1248E-01 -1.6280E+00 -8.3456E-03 -4.9025E-02
             1.3458E-01
 GRADIENT:  -1.2054E+00 -9.2645E-02  5.9100E-03 -4.9518E-01 -1.0240E+00 -2.4793E-02 -2.4997E-01  1.3867E-01 -1.9851E-01  3.0772E-01
             1.6628E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1664.62735167913        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  9.9118E-01  1.0280E+00  8.0467E-01  1.0209E+00  9.2103E-01  1.0122E+00  1.3636E+00  8.2455E-02  8.9717E-01  8.5605E-01
             1.0352E+00
 PARAMETER:  9.1138E-02  1.2758E-01 -1.1732E-01  1.2072E-01  1.7737E-02  1.1210E-01  4.1012E-01 -2.3955E+00 -8.5136E-03 -5.5432E-02
             1.3462E-01
 GRADIENT:  -1.5476E+00  1.7545E-02  4.7842E-01 -8.9067E-01 -1.6010E+00 -9.9677E-04 -4.0957E-01  3.0204E-02 -3.1588E-01  3.7423E-01
             2.2214E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1664.99552412818        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      840
 NPARAMETR:  1.0015E+00  1.0166E+00  8.0768E-01  1.0315E+00  9.1384E-01  1.0228E+00  1.3810E+00  1.0000E-02  8.8517E-01  8.5685E-01
             1.0375E+00
 PARAMETER:  1.0154E-01  1.1643E-01 -1.1359E-01  1.3098E-01  9.9024E-03  1.2255E-01  4.2278E-01 -4.6860E+00 -2.1975E-02 -5.4493E-02
             1.3686E-01
 GRADIENT:  -1.6856E+01  6.7191E-01  3.8480E+00 -6.5037E+00 -7.2217E+00 -1.1628E+00 -2.2864E+00  0.0000E+00 -1.8862E+00  5.4128E-01
             6.3544E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1665.22013714434        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1015
 NPARAMETR:  1.0088E+00  8.8467E-01  8.2386E-01  1.1116E+00  8.6687E-01  1.0238E+00  1.5671E+00  1.0000E-02  8.4149E-01  8.3221E-01
             1.0371E+00
 PARAMETER:  1.0873E-01 -2.2535E-02 -9.3760E-02  2.0578E-01 -4.2865E-02  1.2351E-01  5.4924E-01 -9.0570E+00 -7.2582E-02 -8.3673E-02
             1.3641E-01
 GRADIENT:   5.7924E-02  3.2336E-01 -3.8154E-01  1.4012E+00  1.7273E+00 -3.0865E-01  4.1601E-01  0.0000E+00  3.4475E-01 -6.6502E-01
            -3.8451E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1665.22436789091        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1177
 NPARAMETR:  1.0087E+00  8.7517E-01  8.2304E-01  1.1163E+00  8.6169E-01  1.0245E+00  1.5759E+00  1.0000E-02  8.3654E-01  8.3333E-01
             1.0377E+00
 PARAMETER:  1.0866E-01 -3.3341E-02 -9.4745E-02  2.1003E-01 -4.8858E-02  1.2421E-01  5.5482E-01 -8.5702E+00 -7.8486E-02 -8.2324E-02
             1.3702E-01
 GRADIENT:  -3.4405E-03 -1.7559E-02  3.1319E-02 -6.8517E-02 -3.0720E-02 -3.1114E-04  1.7686E-04  0.0000E+00  2.5740E-03  6.3519E-04
             1.1072E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1177
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6859E-05  9.7548E-03 -4.8865E-04 -1.2974E-02 -8.7363E-03
 SE:             2.9852E-02  2.2446E-02  2.0486E-04  2.3598E-02  2.1896E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9928E-01  6.6387E-01  1.7064E-02  5.8246E-01  6.8990E-01

 ETASHRINKSD(%)  1.0000E-10  2.4802E+01  9.9314E+01  2.0943E+01  2.6646E+01
 ETASHRINKVR(%)  1.0000E-10  4.3452E+01  9.9995E+01  3.7500E+01  4.6192E+01
 EBVSHRINKSD(%)  4.3293E-01  2.4302E+01  9.9381E+01  2.1122E+01  2.5182E+01
 EBVSHRINKVR(%)  8.6398E-01  4.2699E+01  9.9996E+01  3.7782E+01  4.4023E+01
 RELATIVEINF(%)  9.8768E+01  5.9121E+00  5.3457E-04  7.1513E+00  5.5229E+00
 EPSSHRINKSD(%)  4.2716E+01
 EPSSHRINKVR(%)  6.7186E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.2243678909106     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.07354132717239     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.02
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.224       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  8.75E-01  8.23E-01  1.12E+00  8.62E-01  1.02E+00  1.58E+00  1.00E-02  8.37E-01  8.33E-01  1.04E+00
 


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
+        1.03E+03
 
 TH 2
+       -7.79E+00  3.47E+02
 
 TH 3
+        1.81E+01  1.84E+02  7.13E+02
 
 TH 4
+       -7.73E+00  2.80E+02 -3.05E+02  8.07E+02
 
 TH 5
+       -5.77E+00 -3.11E+02 -7.90E+02  3.54E+02  1.19E+03
 
 TH 6
+        7.02E-01 -1.48E+00  2.87E+00 -1.74E+00  8.51E-01  1.85E+02
 
 TH 7
+        1.03E+00  3.03E+01 -8.90E+00 -1.40E+01  1.98E+00  4.22E-02  3.24E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.19E-01 -2.31E+01 -2.98E+01  1.90E+01  1.10E+01 -6.45E-01  1.41E+01  0.00E+00  1.23E+02
 
 TH10
+       -1.48E+00 -5.42E+00 -5.98E+01 -2.62E+01 -5.67E+01  2.19E-02  5.79E+00  0.00E+00  1.53E+01  9.66E+01
 
 TH11
+       -7.44E+00 -1.03E+01 -4.62E+01 -7.70E+00  1.08E+01  3.47E+00  3.38E+00  0.00E+00  9.35E+00  2.30E+01  2.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.228
Stop Time:
Sat Sep 18 11:46:00 CDT 2021
