Wed Sep 29 09:05:10 CDT 2021
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
$DATA ../../../../data/int/D/dat51.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25263.3715041249        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.1370E+02  5.0303E+02 -7.7068E+01  3.4492E+02  2.5744E+02 -1.7976E+03 -8.9036E+02 -8.4717E+01 -1.4827E+03 -6.7455E+02
            -5.2611E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -816.500483279937        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  8.4932E-01  1.8948E+00  9.7882E-01  2.8538E+00  7.6959E-01  3.9376E+00  2.9585E+00  1.0014E+00  2.6712E+00  1.5362E+00
             1.2937E+01
 PARAMETER: -6.3322E-02  7.3909E-01  7.8592E-02  1.1486E+00 -1.6190E-01  1.4706E+00  1.1847E+00  1.0141E-01  1.0825E+00  5.2933E-01
             2.6601E+00
 GRADIENT:  -6.5470E+01  8.6233E+01 -1.1111E+01  1.8394E+02 -6.8549E+01  1.5731E+02 -6.4798E+01  4.7995E+00 -3.1799E+01  3.1320E+01
             5.0068E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -884.232102099040        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.5199E-01  7.6534E-01  1.9948E+02  1.1042E+01  3.3137E+00  2.3523E+00  7.2805E+00  9.5046E-01  6.2447E+00  1.5455E+00
             1.1959E+01
 PARAMETER:  5.0796E-02 -1.6744E-01  5.3957E+00  2.5017E+00  1.2981E+00  9.5539E-01  2.0852E+00  4.9189E-02  1.9317E+00  5.3534E-01
             2.5815E+00
 GRADIENT:  -9.0256E+01  2.9198E+00 -1.6410E+00  1.6235E+02  4.9408E+01  5.6176E+01  1.4295E+01  1.0710E-01  1.7202E+02  2.5632E+01
             5.2660E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -952.346330485521        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  5.7154E-01  9.5459E-01  8.4499E+02  3.5796E+00  2.2717E+00  3.5112E+00  1.1767E+01  3.0378E+01  2.4887E+00  1.3708E+00
             1.0209E+01
 PARAMETER: -4.5941E-01  5.3524E-02  6.8393E+00  1.3753E+00  9.2054E-01  1.3560E+00  2.5653E+00  3.5137E+00  1.0118E+00  4.1543E-01
             2.4233E+00
 GRADIENT:  -9.4056E+01  2.4147E+01  4.0200E+00  2.1224E+02 -2.4928E+01  2.8263E+02  1.9532E+02  8.5586E+01 -1.4622E+01  3.4420E+01
             3.9898E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1120.55191026246        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      371             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0439E+00  3.3489E-01  7.1616E+02  1.6844E+00  1.7250E+00  3.7230E+00  9.4235E+00  2.7269E+01  1.1217E+00  1.9105E+00
             9.1377E+00
 PARAMETER:  1.4293E-01 -9.9395E-01  6.6739E+00  6.2144E-01  6.4520E-01  1.4145E+00  2.3432E+00  3.4057E+00  2.1481E-01  7.4737E-01
             2.3124E+00
 GRADIENT:  -6.9215E+00  1.3749E+01  2.7731E+00  3.1283E+01 -7.0178E+01  2.1154E+02  2.5293E+02  4.1206E+01 -3.9670E+00  5.7402E+01
             3.2706E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1141.82440772375        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      486
 NPARAMETR:  1.3986E+00  1.0228E-01  7.2563E+02  1.6620E+00  1.7311E+00  2.6947E+00  8.8846E+00  2.6887E+01  9.8016E-01  1.9124E+00
             8.7971E+00
 PARAMETER:  4.3546E-01 -2.1800E+00  6.6870E+00  6.0802E-01  6.4877E-01  1.0913E+00  2.2843E+00  3.3916E+00  7.9959E-02  7.4838E-01
             2.2744E+00
 GRADIENT:   7.3525E+01 -1.7504E+01  1.6354E+00 -3.9147E+00 -6.9857E+01  7.8641E+01 -1.3276E+01  1.1692E+01 -1.9922E+01  5.9805E+01
             1.8241E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1192.07552911454        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  1.1005E+00  3.3607E-01  7.4234E+02  1.5486E+00  1.7633E+00  2.2379E+00  8.0174E+00  2.5925E+01  1.1627E+00  1.9102E+00
             7.6250E+00
 PARAMETER:  1.9579E-01 -9.9042E-01  6.7098E+00  5.3732E-01  6.6717E-01  9.0554E-01  2.1816E+00  3.3552E+00  2.5071E-01  7.4719E-01
             2.1314E+00
 GRADIENT:   1.7904E-01 -1.0787E-01  7.6378E-01  4.3069E+00 -8.2728E+01 -9.9588E-01  2.5891E+01  4.8899E+00 -1.1698E+00  5.2221E+01
            -4.0840E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1208.78874592271        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      838
 NPARAMETR:  1.1177E+00  1.3197E+00  7.3526E+02  9.7104E-01  2.0595E+00  2.3163E+00  4.3795E+00  2.2841E+01  5.5750E-01  1.8445E+00
             7.7555E+00
 PARAMETER:  2.1123E-01  3.7743E-01  6.7002E+00  7.0616E-02  8.2245E-01  9.3998E-01  1.5769E+00  3.2286E+00 -4.8428E-01  7.1222E-01
             2.1484E+00
 GRADIENT:   3.3783E+00  9.4026E-01  6.7148E-02 -3.4551E+00 -4.1526E+01  9.9693E+00 -9.5973E+00  1.5074E-02  4.8579E-01  4.5767E+01
            -3.6491E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1213.29854196376        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1015
 NPARAMETR:  1.0972E+00  1.1814E+00  7.2741E+02  1.0837E+00  2.2929E+00  2.1189E+00  4.9791E+00  2.2487E+01  6.2075E-01  1.7694E+00
             7.8693E+00
 PARAMETER:  1.9280E-01  2.6674E-01  6.6895E+00  1.8036E-01  9.2983E-01  8.5091E-01  1.7052E+00  3.2129E+00 -3.7683E-01  6.7064E-01
             2.1630E+00
 GRADIENT:  -7.0095E+00  4.3489E+00  2.5751E-02 -1.8764E+00 -5.2508E+00 -2.2934E+01  6.6870E+00  1.9002E-02 -1.2533E+00  4.1436E+01
             1.9326E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1217.43876743916        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1190
 NPARAMETR:  1.1114E+00  6.8264E-01  7.2148E+02  1.3877E+00  2.3923E+00  2.2509E+00  5.8714E+00  2.2872E+01  1.0635E+00  1.7209E+00
             7.7751E+00
 PARAMETER:  2.0564E-01 -2.8179E-01  6.6813E+00  4.2765E-01  9.7224E-01  9.1132E-01  1.8701E+00  3.2299E+00  1.6157E-01  6.4286E-01
             2.1509E+00
 GRADIENT:   1.9013E-01 -1.8800E-01  1.3992E-02  5.6938E+00  9.6965E+00 -4.1306E-02 -8.1842E-01  1.8987E-02  6.1950E-02  3.7876E+01
            -1.7809E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1217.46691075169        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1365
 NPARAMETR:  1.1110E+00  7.3707E-01  7.2089E+02  1.3420E+00  2.4144E+00  2.2564E+00  5.7418E+00  2.2732E+01  1.0117E+00  1.7134E+00
             7.7769E+00
 PARAMETER:  2.0525E-01 -2.0507E-01  6.6805E+00  3.9416E-01  9.8145E-01  9.1375E-01  1.8478E+00  3.2238E+00  1.1162E-01  6.3846E-01
             2.1512E+00
 GRADIENT:   7.8319E-02 -2.8287E-01  3.5003E-03  3.8981E-01  1.3316E+01  6.9712E-01  7.0242E-01  2.0220E-02  1.2185E-02  3.7061E+01
             1.5616E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1217.58266949818        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1543
 NPARAMETR:  1.1101E+00  9.1642E-01  7.2058E+02  1.2305E+00  2.4397E+00  2.2644E+00  5.3437E+00  2.2351E+01  9.0111E-01  1.6726E+00
             7.7855E+00
 PARAMETER:  2.0448E-01  1.2715E-02  6.6801E+00  3.0740E-01  9.9189E-01  9.1731E-01  1.7759E+00  3.2069E+00 -4.1225E-03  6.1441E-01
             2.1523E+00
 GRADIENT:  -3.0517E-01  6.7995E-01 -6.3960E-03 -5.8260E+00  1.7852E+01  1.6872E+00  2.5414E+00  1.8662E-02  5.8214E-01  3.4463E+01
             2.9198E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1222.78824999535        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1725
 NPARAMETR:  1.1145E+00  1.3193E+00  7.2108E+02  8.9204E-01  2.3165E+00  2.1951E+00  4.5204E+00  2.0423E+01  4.0512E-01  1.0421E+00
             7.7249E+00
 PARAMETER:  2.0839E-01  3.7707E-01  6.6808E+00 -1.4250E-02  9.4007E-01  8.8623E-01  1.6086E+00  3.1167E+00 -8.0356E-01  1.4128E-01
             2.1444E+00
 GRADIENT:   3.9341E+00 -8.4544E+00  1.6252E-02 -3.8072E+01  8.2470E+00 -1.1162E+01  1.3553E+01  3.6470E-03  2.0605E+00  7.4564E+00
            -1.5000E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1225.94694900117        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1901
 NPARAMETR:  1.1098E+00  1.3588E+00  7.1313E+02  9.4011E-01  2.3364E+00  2.2587E+00  4.4857E+00  1.9317E+01  2.7249E-01  7.1265E-01
             7.9006E+00
 PARAMETER:  2.0422E-01  4.0659E-01  6.6697E+00  3.8241E-02  9.4861E-01  9.1477E-01  1.6009E+00  3.0610E+00 -1.2002E+00 -2.3877E-01
             2.1669E+00
 GRADIENT:   1.1161E-02  6.2470E-02  3.5869E-02  1.0577E-01  1.3946E-01  9.0192E-01 -2.8759E-01 -5.8660E-03 -2.1264E-01  4.6490E-01
             9.5904E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1226.07123478025        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2087             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1091E+00  1.3365E+00  6.9104E+02  9.6486E-01  2.3374E+00  2.2892E+00  4.5477E+00  1.9621E+01  3.4113E-01  6.7620E-01
             7.8910E+00
 PARAMETER:  2.0356E-01  3.9003E-01  6.6382E+00  6.4227E-02  9.4905E-01  9.2822E-01  1.6146E+00  3.0766E+00 -9.7548E-01 -2.9127E-01
             2.1657E+00
 GRADIENT:   1.7939E+01  1.0483E+01  4.7303E-02  6.8646E+00  3.2381E+00  6.2030E+01  8.5742E+01 -1.2907E-02 -4.4038E-01 -6.6860E-02
             3.4627E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1226.51432490013        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2263
 NPARAMETR:  1.1153E+00  1.1879E+00  7.1137E+02  1.0524E+00  2.3410E+00  2.1947E+00  4.7275E+00  1.9260E+01  6.7074E-01  7.0553E-01
             7.9163E+00
 PARAMETER:  2.0916E-01  2.7221E-01  6.6672E+00  1.5105E-01  9.5058E-01  8.8603E-01  1.6534E+00  3.0580E+00 -2.9937E-01 -2.4881E-01
             2.1689E+00
 GRADIENT:   1.8392E+00  6.7309E-02  5.4534E-02 -1.0180E+01  2.3803E+00 -9.2899E+00  7.7016E-01 -2.3410E-02  1.0240E+00  2.4592E-01
             8.6238E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1227.23244613058        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2441
 NPARAMETR:  1.0971E+00  9.9493E-01  7.1296E+02  1.1849E+00  2.3086E+00  2.2183E+00  4.9087E+00  1.8589E+01  9.0987E-01  7.7837E-01
             7.8251E+00
 PARAMETER:  1.9263E-01  9.4913E-02  6.6694E+00  2.6965E-01  9.3664E-01  8.9675E-01  1.6910E+00  3.0226E+00  5.5444E-03 -1.5055E-01
             2.1573E+00
 GRADIENT:  -4.0023E+00  3.5381E-01  8.5971E-02  1.1522E+00 -3.8614E+00 -5.2493E+00 -9.4216E+00 -3.9649E-02  7.4121E-01  6.1814E-01
            -1.1289E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1229.71754604652        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2624
 NPARAMETR:  1.1157E+00  8.1201E-01  5.9707E+02  1.2907E+00  2.2994E+00  2.2668E+00  5.5116E+00  1.7537E+01  9.8765E-01  6.3414E-01
             7.8598E+00
 PARAMETER:  2.0947E-01 -1.0824E-01  6.4920E+00  3.5519E-01  9.3263E-01  9.1837E-01  1.8069E+00  2.9643E+00  8.7571E-02 -3.5548E-01
             2.1618E+00
 GRADIENT:   2.7608E+00  3.5571E-02  5.7144E-01 -2.8047E-01 -1.0282E+00  3.1077E+00  1.4921E+00 -1.0700E+00 -5.4365E-03  3.5986E-01
            -5.3764E-01

0ITERATION NO.:   88    OBJECTIVE VALUE:  -1229.72511282798        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     2722
 NPARAMETR:  1.1152E+00  8.0596E-01  5.9336E+02  1.2940E+00  2.3035E+00  2.2480E+00  5.5207E+00  1.7570E+01  9.9102E-01  6.3444E-01
             7.8730E+00
 PARAMETER:  2.0935E-01 -1.1640E-01  6.4787E+00  3.5816E-01  9.3330E-01  9.1100E-01  1.8099E+00  2.9696E+00  9.1984E-02 -3.5548E-01
             2.1629E+00
 GRADIENT:   5.3392E+01 -6.9075E-02 -2.6740E+00  5.8495E+01 -1.8654E+01  1.9788E+01  8.4978E+00  6.4853E+00  1.0125E-01 -5.9126E+01
            -6.1310E+00
 NUMSIGDIG:         2.3         1.7         2.4         2.3         2.3         2.4         2.5         2.4         1.4         2.3
                    3.0

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2722
 NO. OF SIG. DIGITS IN FINAL EST.:  1.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.9551E-03  3.8063E-02 -5.0395E-03 -7.2257E-02 -2.8658E-03
 SE:             2.9256E-02  2.4470E-02  3.8033E-03  1.3909E-02  1.0080E-02
 N:                     100         100         100         100         100

 P VAL.:         8.1209E-01  1.1983E-01  1.8516E-01  2.0531E-07  7.7618E-01

 ETASHRINKSD(%)  1.9897E+00  1.8021E+01  8.7259E+01  5.3402E+01  6.6230E+01
 ETASHRINKVR(%)  3.9398E+00  3.2794E+01  9.8377E+01  7.8286E+01  8.8596E+01
 EBVSHRINKSD(%)  3.0248E+00  1.4802E+01  8.9882E+01  5.5455E+01  6.6505E+01
 EBVSHRINKVR(%)  5.9581E+00  2.7412E+01  9.8976E+01  8.0157E+01  8.8781E+01
 RELATIVEINF(%)  9.3903E+01  3.4239E+01  9.8052E-01  9.3228E+00  1.0435E+01
 EPSSHRINKSD(%)  7.6097E+00
 EPSSHRINKVR(%)  1.4640E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1229.7251128279827     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       424.36424694042807     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    85.76
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1229.725       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.12E+00  8.05E-01  5.89E+02  1.29E+00  2.30E+00  2.25E+00  5.53E+00  1.76E+01  9.92E-01  6.34E-01  7.87E+00
 


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
+        1.98E+04
 
 TH 2
+       -8.67E+01  4.08E+01
 
 TH 3
+        1.90E-03  6.71E-03  3.18E-05
 
 TH 4
+       -1.09E+01  3.19E+00  2.33E-02  2.61E+03
 
 TH 5
+       -1.74E-01  2.97E+00 -5.41E-02 -2.92E+00  1.19E+02
 
 TH 6
+        1.12E+03 -2.80E+03 -1.23E-03  5.70E+00 -1.67E+00 -9.39E+01
 
 TH 7
+        2.49E+02  2.74E+00 -3.89E-03  8.84E+00 -3.64E+00  2.83E+01  1.00E+01
 
 TH 8
+       -1.37E-01 -5.07E-01  1.27E-04 -1.63E+00  2.65E+00  9.04E-02  2.87E-01  2.01E-01
 
 TH 9
+        1.46E+01 -4.06E+00 -8.25E-03 -3.67E+01  2.46E+00 -2.65E+03  2.94E+00  6.05E-01  2.76E+01
 
 TH10
+        1.41E+01  8.68E+01 -3.46E-02  1.29E+02 -2.35E+01 -9.41E+00 -2.54E+02  2.52E+00 -1.45E+01  1.03E+04
 
 TH11
+       -5.76E+00 -2.39E+00 -4.63E-03 -1.01E+01 -6.29E+00  3.64E+01 -6.10E+00  3.20E-01  3.24E+00 -1.31E+00  1.36E+01
 
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
 #CPUT: Total CPU Time in Seconds,      102.251
Stop Time:
Wed Sep 29 09:06:54 CDT 2021
