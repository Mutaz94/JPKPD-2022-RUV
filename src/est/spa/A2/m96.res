Wed Sep 29 13:11:04 CDT 2021
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
$DATA ../../../../data/spa/A2/dat96.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -976.843786631395        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7512E+02  2.6186E+01  1.3039E+01  3.9823E+01  6.1293E+01  2.3217E+01  1.5953E+00 -6.1156E+00 -8.6328E+00 -1.6052E+01
            -1.2430E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1369.84506675022        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1601E+00  1.0193E+00  1.0775E+00  1.0135E+00  9.9798E-01  1.2907E+00  8.9451E-01  9.4677E-01  9.3689E-01  8.1896E-01
             2.4977E+00
 PARAMETER:  2.4851E-01  1.1915E-01  1.7468E-01  1.1343E-01  9.7975E-02  3.5515E-01 -1.1475E-02  4.5301E-02  3.4813E-02 -9.9725E-02
             1.0154E+00
 GRADIENT:   3.8843E+02 -8.7690E+00  7.9580E+00 -2.6185E+01 -1.4653E+01  5.9042E+01  8.2747E+00  2.9600E+00  3.0092E+00  1.2405E+01
            -4.3678E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1399.47464279237        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0151E+00  9.4657E-01  5.3457E-01  1.0312E+00  6.8760E-01  9.6717E-01  7.5637E-01  1.8515E-01  1.0038E+00  4.7306E-01
             2.3659E+00
 PARAMETER:  1.1494E-01  4.5094E-02 -5.2628E-01  1.3074E-01 -2.7455E-01  6.6624E-02 -1.7922E-01 -1.5866E+00  1.0381E-01 -6.4854E-01
             9.6116E-01
 GRADIENT:   1.6599E+02 -2.9898E+01 -3.8286E+01  4.5493E+01  6.8236E+01 -4.9534E+00 -6.2955E+00  1.1969E-01  1.9174E+01  8.2563E-02
            -4.3443E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1403.30849776516        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      251
 NPARAMETR:  9.7629E-01  9.4907E-01  4.4095E-01  1.0057E+00  6.0416E-01  9.1798E-01  9.9960E-01  2.3003E-01  8.9461E-01  3.2889E-01
             2.3487E+00
 PARAMETER:  7.6008E-02  4.7730E-02 -7.1881E-01  1.0566E-01 -4.0391E-01  1.4420E-02  9.9599E-02 -1.3695E+00 -1.1362E-02 -1.0120E+00
             9.5386E-01
 GRADIENT:  -1.3183E+00  1.7054E+01 -3.1136E+00  3.2189E+01 -3.6158E+00 -2.7674E+01  5.1932E+00  8.7874E-02  1.2375E+01  1.1986E+00
            -4.1694E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1407.39471662578        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  9.7942E-01  9.8349E-01  4.3260E-01  9.6774E-01  6.1703E-01  9.7617E-01  9.3694E-01  2.4590E-01  7.9938E-01  2.2353E-01
             2.5790E+00
 PARAMETER:  7.9208E-02  8.3349E-02 -7.3794E-01  6.7206E-02 -3.8284E-01  7.5879E-02  3.4868E-02 -1.3028E+00 -1.2392E-01 -1.3982E+00
             1.0474E+00
 GRADIENT:   3.2627E+00 -3.4520E-01  1.0950E+00 -2.3681E+00 -3.8536E+00  7.2039E-02  7.4334E-01  1.3825E-01  8.7879E-01  8.2398E-01
             4.6568E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1407.87431720494        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.7481E-01  1.2190E+00  4.5178E-01  8.4352E-01  7.5157E-01  9.7755E-01  8.2686E-01  2.9834E-01  8.7783E-01  1.0849E-01
             2.6092E+00
 PARAMETER:  7.4485E-02  2.9807E-01 -6.9456E-01 -7.0173E-02 -1.8559E-01  7.7290E-02 -9.0124E-02 -1.1095E+00 -3.0301E-02 -2.1211E+00
             1.0590E+00
 GRADIENT:  -5.5890E+00 -5.6852E+00 -3.1787E+00  3.9402E-01  7.5433E+00  1.8729E+00  7.9393E-01  9.9091E-02 -7.6541E-01  1.4071E-01
             2.9791E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1407.92729840048        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      780
 NPARAMETR:  9.7670E-01  1.2820E+00  4.7962E-01  8.1412E-01  8.0488E-01  9.7008E-01  7.9442E-01  2.0271E-01  9.2561E-01  9.8285E-02
             2.6353E+00
 PARAMETER:  7.6429E-02  3.4842E-01 -6.3477E-01 -1.0565E-01 -1.1706E-01  6.9623E-02 -1.3014E-01 -1.4960E+00  2.2696E-02 -2.2199E+00
             1.0690E+00
 GRADIENT:  -1.3529E+00 -3.1810E-01  2.0859E-02 -2.3423E-01  3.0064E-01 -2.5164E-01  3.1544E-01  1.5482E-02 -2.3870E-03  6.5190E-02
             5.3159E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1407.94562586978        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      961
 NPARAMETR:  9.7915E-01  1.4018E+00  6.2089E-01  7.6704E-01  9.8062E-01  9.6546E-01  7.2533E-01  1.5489E-02  1.0417E+00  4.0106E-02
             2.7110E+00
 PARAMETER:  7.8926E-02  4.3774E-01 -3.7661E-01 -1.6522E-01  8.0429E-02  6.4848E-02 -2.2113E-01 -4.0676E+00  1.4086E-01 -3.1162E+00
             1.0973E+00
 GRADIENT:   5.0445E+00  2.7442E+00 -4.0643E+00  6.5079E+00  6.8421E+00 -1.8052E-01  1.7756E+00  3.1164E-04  1.6036E-01  1.0623E-03
             1.8432E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1409.05973403041        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1143
 NPARAMETR:  9.7416E-01  1.2133E+00  1.0576E+00  8.9841E-01  1.0916E+00  9.5756E-01  6.2288E-01  1.0000E-02  1.0225E+00  1.0000E-02
             2.7694E+00
 PARAMETER:  7.3825E-02  2.9332E-01  1.5601E-01 -7.1308E-03  1.8764E-01  5.6632E-02 -3.7340E-01 -1.0960E+01  1.2227E-01 -5.4717E+00
             1.1186E+00
 GRADIENT:   2.2955E+00  1.0228E+01  4.5899E+00  1.4460E+00 -1.1248E+01  1.4397E-01  2.0352E+00  0.0000E+00 -3.4277E-01  0.0000E+00
            -2.5800E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1409.71722456950        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1321
 NPARAMETR:  9.6969E-01  1.0488E+00  1.1864E+00  1.0072E+00  1.0901E+00  9.5380E-01  4.6902E-01  1.0000E-02  9.8974E-01  6.0435E-02
             2.7850E+00
 PARAMETER:  6.9224E-02  1.4766E-01  2.7092E-01  1.0721E-01  1.8631E-01  5.2697E-02 -6.5710E-01 -7.8530E+00  8.9685E-02 -2.7062E+00
             1.1242E+00
 GRADIENT:  -6.3188E+00  1.8000E+00 -3.0779E+00  4.9852E+00  5.9548E+00 -1.1271E+00  6.2447E-01  0.0000E+00  3.8453E-01 -1.8485E-02
             1.1467E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1410.12844343878        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1497
 NPARAMETR:  9.6944E-01  9.3557E-01  1.2259E+00  1.0751E+00  1.0432E+00  9.5710E-01  4.3840E-01  1.0000E-02  9.4071E-01  3.6316E-01
             2.7415E+00
 PARAMETER:  6.8964E-02  3.3405E-02  3.0370E-01  1.7237E-01  1.4229E-01  5.6151E-02 -7.2463E-01 -5.9644E+00  3.8875E-02 -9.1291E-01
             1.1085E+00
 GRADIENT:  -4.4606E+00  2.2062E+00  1.2834E+00  1.8964E-02 -2.0822E+00 -1.4637E-01  5.9214E-01  0.0000E+00 -6.7707E-02 -2.3051E-02
            -3.0376E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1410.16730874033        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1675
 NPARAMETR:  9.6824E-01  8.4284E-01  1.2940E+00  1.1416E+00  1.0336E+00  9.5289E-01  4.5514E-01  1.0000E-02  8.8040E-01  3.6504E-01
             2.7524E+00
 PARAMETER:  6.7722E-02 -7.0975E-02  3.5771E-01  2.3242E-01  1.3304E-01  5.1743E-02 -6.8714E-01 -6.7699E+00 -2.7380E-02 -9.0774E-01
             1.1125E+00
 GRADIENT:  -6.4561E+00  6.5299E+00  1.5609E+00  8.4918E+00 -1.8544E+00 -1.3103E+00  4.5669E-01  0.0000E+00 -1.5387E+00  4.7562E-02
             9.2485E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1410.31480225248        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1851
 NPARAMETR:  9.6855E-01  6.6423E-01  1.1767E+00  1.2479E+00  9.1919E-01  9.5730E-01  6.5459E-01  1.0000E-02  8.0786E-01  3.7833E-01
             2.7327E+00
 PARAMETER:  6.8044E-02 -3.0913E-01  2.6268E-01  3.2145E-01  1.5738E-02  5.6358E-02 -3.2375E-01 -8.4487E+00 -1.1336E-01 -8.7198E-01
             1.1053E+00
 GRADIENT:  -3.2396E+00  4.7839E+00  1.1761E+00  7.7362E+00 -3.5019E+00 -1.2188E-01  5.2500E-01  0.0000E+00 -2.8574E-01  2.0574E-02
             4.5097E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1410.38210218967        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2026
 NPARAMETR:  9.6626E-01  5.4029E-01  1.2026E+00  1.3237E+00  8.8985E-01  9.5461E-01  8.6650E-01  1.0000E-02  7.5001E-01  3.6780E-01
             2.7366E+00
 PARAMETER:  6.5678E-02 -5.1566E-01  2.8449E-01  3.8040E-01 -1.6699E-02  5.3549E-02 -4.3289E-02 -1.1016E+01 -1.8767E-01 -9.0022E-01
             1.1067E+00
 GRADIENT:  -5.0698E+00  2.7701E+00  5.8121E-01  4.7391E+00 -4.0610E-01 -5.3310E-01  7.2882E-01  0.0000E+00 -5.6784E-01  7.2666E-02
             9.1713E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1410.52459212838        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2207             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6944E-01  4.8613E-01  1.1799E+00  1.3474E+00  8.6627E-01  9.5646E-01  7.2767E-01  1.0000E-02  7.3361E-01  3.4982E-01
             2.7310E+00
 PARAMETER:  6.8961E-02 -6.2127E-01  2.6543E-01  3.9819E-01 -4.3557E-02  5.5479E-02 -2.1790E-01 -1.2148E+01 -2.0978E-01 -9.5034E-01
             1.1047E+00
 GRADIENT:   5.7021E+01  5.3085E+00 -1.3004E+00  4.7819E+01  4.3237E+00  4.6964E+00  2.0994E-01  0.0000E+00 -2.6232E+00 -1.8816E-01
             5.5261E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1410.64193817100        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     2302
 NPARAMETR:  9.6450E-01  4.7723E-01  1.1819E+00  1.3447E+00  8.6378E-01  9.5376E-01  3.5409E-01  1.0000E-02  7.5513E-01  3.8505E-01
             2.7229E+00
 PARAMETER:  6.3849E-02 -6.3975E-01  2.6714E-01  3.9619E-01 -4.6437E-02  5.2657E-02 -9.3821E-01 -1.2148E+01 -1.8086E-01 -8.5438E-01
             1.1017E+00
 GRADIENT:  -5.1632E+00 -3.7897E+00 -1.9632E+00 -2.0952E+01  3.5741E+00 -6.0652E-01  5.5027E-02  0.0000E+00 -2.4595E+00 -2.6205E-01
            -2.1483E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1410.77579979228        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2480
 NPARAMETR:  9.6791E-01  4.7626E-01  1.1912E+00  1.3537E+00  8.6493E-01  9.5702E-01  2.5033E-01  1.0000E-02  7.6952E-01  4.2984E-01
             2.7080E+00
 PARAMETER:  6.7379E-02 -6.4179E-01  2.7497E-01  4.0287E-01 -4.5109E-02  5.6069E-02 -1.2850E+00 -1.2148E+01 -1.6198E-01 -7.4434E-01
             1.0962E+00
 GRADIENT:   1.8282E+00 -3.6345E-01  3.3499E-01 -4.9755E+00 -8.9037E-01  6.2910E-02  4.9599E-02  0.0000E+00  5.7043E-01 -4.0720E-02
            -2.2740E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1410.86531068029        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2655
 NPARAMETR:  9.6536E-01  4.0107E-01  1.2688E+00  1.4072E+00  8.7499E-01  9.5495E-01  1.0000E-02  1.0000E-02  7.3594E-01  4.0834E-01
             2.7310E+00
 PARAMETER:  6.4749E-02 -8.1361E-01  3.3804E-01  4.4162E-01 -3.3538E-02  5.3904E-02 -5.4135E+00 -1.2148E+01 -2.0661E-01 -7.9565E-01
             1.1047E+00
 GRADIENT:  -1.9903E+00  6.6696E-01  6.8033E-01  1.5421E+00  6.4728E-02  5.3819E-02  0.0000E+00  0.0000E+00  1.3891E-01 -1.4959E-02
             7.7752E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1410.87366707696        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2833
 NPARAMETR:  9.6548E-01  3.3215E-01  1.3055E+00  1.4506E+00  8.6879E-01  9.5371E-01  1.0000E-02  1.0000E-02  7.0880E-01  4.1232E-01
             2.7249E+00
 PARAMETER:  6.4869E-02 -1.0022E+00  3.6655E-01  4.7197E-01 -4.0655E-02  5.2603E-02 -8.8551E+00 -1.2148E+01 -2.4418E-01 -7.8595E-01
             1.1024E+00
 GRADIENT:   1.0054E+00  3.5728E-01  1.0557E+00  3.0471E-01  7.7399E-01 -7.1782E-02  0.0000E+00  0.0000E+00 -5.2104E-01 -9.3842E-03
            -1.3191E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1410.99862673323        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3013
 NPARAMETR:  9.6634E-01  1.8655E-01  1.2064E+00  1.5305E+00  7.9026E-01  9.5389E-01  1.0000E-02  1.0000E-02  6.6743E-01  4.5594E-01
             2.6907E+00
 PARAMETER:  6.5762E-02 -1.5790E+00  2.8762E-01  5.2561E-01 -1.3539E-01  5.2791E-02 -1.6228E+01 -1.2148E+01 -3.0432E-01 -6.8538E-01
             1.0898E+00
 GRADIENT:   8.1985E+00  4.7736E-01  8.8808E-02  2.2124E+00  1.0916E+00 -3.1303E-01  0.0000E+00  0.0000E+00 -2.2517E+00  9.3004E-02
            -5.4281E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1411.11696999559        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3190
 NPARAMETR:  9.6175E-01  1.2833E-01  1.0771E+00  1.5504E+00  7.2113E-01  9.5526E-01  1.0000E-02  1.0000E-02  6.6522E-01  4.4513E-01
             2.6997E+00
 PARAMETER:  6.1004E-02 -1.9531E+00  1.7425E-01  5.3849E-01 -2.2693E-01  5.4231E-02 -1.9961E+01 -1.2148E+01 -3.0764E-01 -7.0938E-01
             1.0932E+00
 GRADIENT:  -1.0496E+00  3.3655E-01  2.3505E-02  1.4004E+00 -1.2260E+00 -2.5239E-02  0.0000E+00  0.0000E+00  1.0998E-01  2.5698E-07
             7.2634E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1411.15340867271        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3370
 NPARAMETR:  9.5976E-01  6.9593E-02  1.1028E+00  1.5886E+00  7.1816E-01  9.5418E-01  1.0000E-02  1.0000E-02  6.4794E-01  4.3348E-01
             2.7116E+00
 PARAMETER:  5.8923E-02 -2.5651E+00  1.9790E-01  5.6283E-01 -2.3106E-01  5.3099E-02 -2.9412E+01 -1.2148E+01 -3.3396E-01 -7.3591E-01
             1.0976E+00
 GRADIENT:  -3.2790E+00  3.1393E-01  7.5247E-01  6.6224E+00 -2.1621E+00  3.8110E-02  0.0000E+00  0.0000E+00  2.6376E-01 -3.1783E-02
             2.3965E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1411.17544711632        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3547
 NPARAMETR:  9.6087E-01  4.1441E-02  1.1159E+00  1.6050E+00  7.1679E-01  9.5373E-01  1.0000E-02  1.0000E-02  6.3834E-01  4.4924E-01
             2.6926E+00
 PARAMETER:  6.0084E-02 -3.0835E+00  2.0962E-01  5.7315E-01 -2.3297E-01  5.2628E-02 -3.7064E+01 -1.2148E+01 -3.4889E-01 -7.0019E-01
             1.0905E+00
 GRADIENT:   1.3282E+00  1.5342E-01  8.9571E-01  4.3825E+00 -1.4504E+00 -9.0870E-02  0.0000E+00  0.0000E+00 -7.5216E-01 -3.9489E-02
            -1.8688E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1411.21016875069        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3725
 NPARAMETR:  9.5931E-01  1.3796E-02  1.0651E+00  1.6132E+00  6.8952E-01  9.5406E-01  1.0000E-02  1.0000E-02  6.3651E-01  4.4983E-01
             2.6951E+00
 PARAMETER:  5.8455E-02 -4.1834E+00  1.6308E-01  5.7821E-01 -2.7176E-01  5.2974E-02 -5.2089E+01 -1.2148E+01 -3.5176E-01 -6.9888E-01
             1.0914E+00
 GRADIENT:  -1.0088E+00  3.1338E-02  2.5434E-01  1.3875E+00 -8.9968E-01 -1.5015E-02  0.0000E+00  0.0000E+00  1.0473E-01 -4.7124E-03
             7.4531E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1411.21388107011        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     3897
 NPARAMETR:  9.5971E-01  1.0000E-02  1.0754E+00  1.6152E+00  6.9362E-01  9.5391E-01  1.0000E-02  1.0000E-02  6.3464E-01  4.4840E-01
             2.6940E+00
 PARAMETER:  5.8871E-02 -4.5408E+00  1.7269E-01  5.7948E-01 -2.6584E-01  5.2816E-02 -5.6719E+01 -1.2148E+01 -3.5470E-01 -7.0208E-01
             1.0910E+00
 GRADIENT:   4.8459E-01  7.0497E-04 -3.5763E-01 -2.5106E+00  7.9756E-01  2.8711E-02  0.0000E+00  0.0000E+00 -3.3434E-04 -2.1876E-02
            -8.9516E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3897
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.2220E-04 -2.8239E-06  1.2195E-04 -1.0876E-02 -7.0036E-03
 SE:             2.9069E-02  1.8068E-06  1.5384E-04  2.5591E-02  1.1223E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9116E-01  1.1806E-01  4.2794E-01  6.7084E-01  5.3261E-01

 ETASHRINKSD(%)  2.6152E+00  9.9994E+01  9.9485E+01  1.4267E+01  6.2401E+01
 ETASHRINKVR(%)  5.1620E+00  1.0000E+02  9.9997E+01  2.6499E+01  8.5863E+01
 EBVSHRINKSD(%)  2.5888E+00  9.9994E+01  9.9462E+01  1.4177E+01  6.2727E+01
 EBVSHRINKVR(%)  5.1106E+00  1.0000E+02  9.9997E+01  2.6344E+01  8.6107E+01
 RELATIVEINF(%)  8.7623E+01  8.4571E-09  9.1903E-05  2.8812E+00  2.8200E-01
 EPSSHRINKSD(%)  2.4890E+01
 EPSSHRINKVR(%)  4.3585E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1411.2138810701085     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -676.06305450637035     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    46.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1411.214       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.60E-01  1.00E-02  1.08E+00  1.62E+00  6.94E-01  9.54E-01  1.00E-02  1.00E-02  6.35E-01  4.48E-01  2.69E+00
 


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
+        1.26E+03
 
 TH 2
+       -8.60E+00  1.71E+02
 
 TH 3
+       -8.15E-01  1.38E+01  2.59E+02
 
 TH 4
+       -7.86E+01  5.46E+01 -7.60E+00  7.78E+02
 
 TH 5
+        2.96E+01 -5.14E+01 -5.99E+02 -2.21E+02  1.50E+03
 
 TH 6
+        1.06E-01 -1.05E+00  6.86E+00 -1.72E+01 -6.04E+00  1.95E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.15E+00 -8.98E+00  1.13E+01 -2.31E+01  1.42E+01 -3.64E+00  0.00E+00  0.00E+00  2.62E+02
 
 TH10
+       -7.04E+00 -1.38E+00 -7.37E+00 -1.13E+01  2.36E+01 -8.20E-01  0.00E+00  0.00E+00  3.29E+00  2.22E+01
 
 TH11
+       -1.42E+01 -1.10E+00  1.78E+00 -1.08E+01 -9.69E+00  4.93E+00  0.00E+00  0.00E+00  2.30E+01  1.81E+01  5.10E+01
 
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
 #CPUT: Total CPU Time in Seconds,       51.623
Stop Time:
Wed Sep 29 13:12:08 CDT 2021
