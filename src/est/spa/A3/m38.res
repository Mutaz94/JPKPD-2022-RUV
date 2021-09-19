Sat Sep 18 10:25:47 CDT 2021
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
$DATA ../../../../data/spa/A3/dat38.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1160.40316178538        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.9375E+01  1.1912E+02  4.8576E+01  1.0687E+02  2.0563E+02  1.4397E+00 -9.2533E+01 -3.9223E+01 -1.7702E+02 -2.1251E+02
            -5.0965E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -803.851880317744        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3405E+00  1.0209E+00  9.0292E-01  1.2923E+00  7.1666E-01  7.2171E-01  1.1029E+00  1.0324E+00  1.1254E+00  1.2980E+00
             1.5692E+01
 PARAMETER:  3.9306E-01  1.2071E-01 -2.1209E-03  3.5640E-01 -2.3316E-01 -2.2614E-01  1.9797E-01  1.3186E-01  2.1818E-01  3.6085E-01
             2.8531E+00
 GRADIENT:   2.2939E+01 -3.8911E+01 -7.4754E+00 -6.4008E+01 -7.0110E-02  9.5533E-01  8.2832E+00  3.6911E+00  2.3469E+01  2.0952E+01
             4.3028E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1018.87874997390        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.7340E-01  3.9307E-01  2.3041E-01  1.0520E+00  2.6101E-01  6.6603E-01  2.7972E+00  3.0958E-01  1.8016E+00  8.0878E-01
             7.0589E+00
 PARAMETER:  7.3043E-02 -8.3377E-01 -1.3679E+00  1.5071E-01 -1.2432E+00 -3.0641E-01  1.1286E+00 -1.0725E+00  6.8865E-01 -1.1222E-01
             2.0543E+00
 GRADIENT:  -3.1396E+02  3.6632E+01 -3.3801E+01  1.0566E+01  1.3677E+01 -5.6952E+01  4.1060E+01  1.0803E+00  5.6471E+01  2.8614E+01
             2.6880E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1172.73696739605        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9373E-01  1.3866E+00  1.1841E+00  7.1509E-01  1.1773E+00  8.8207E-01  2.0566E-01  8.1792E-02  7.4792E-01  5.9083E-01
             5.0279E+00
 PARAMETER:  9.3709E-02  4.2684E-01  2.6898E-01 -2.3534E-01  2.6322E-01 -2.5480E-02 -1.4815E+00 -2.4036E+00 -1.9046E-01 -4.2623E-01
             1.7150E+00
 GRADIENT:  -4.7474E+01 -5.1315E+01 -1.3789E-01 -4.8193E+01 -4.0907E+00  7.4417E-01 -1.3698E+00  6.6828E-03 -7.7831E+00  4.2406E+00
             2.1927E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1180.66268480854        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0230E+00  1.3583E+00  1.4747E+00  8.0550E-01  1.2723E+00  8.7148E-01  2.7031E-01  3.9881E-01  1.1236E+00  2.9512E-01
             4.7963E+00
 PARAMETER:  1.2273E-01  4.0623E-01  4.8848E-01 -1.1629E-01  3.4086E-01 -3.7560E-02 -1.2082E+00 -8.1928E-01  2.1651E-01 -1.1204E+00
             1.6678E+00
 GRADIENT:  -3.6234E+00  2.2663E+00  6.3292E-01  8.1851E-01 -3.0432E+00 -4.6050E+00 -2.7088E-01  1.5733E-01  1.9978E+00  7.6807E-01
             2.1489E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1181.33741788309        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0243E+00  1.1309E+00  1.2556E+00  9.5387E-01  1.1462E+00  8.8783E-01  5.1692E-01  8.4897E-02  9.0876E-01  1.5389E-01
             4.7706E+00
 PARAMETER:  1.2402E-01  2.2300E-01  3.2765E-01  5.2775E-02  2.3644E-01 -1.8973E-02 -5.5987E-01 -2.3663E+00  4.3215E-03 -1.7715E+00
             1.6625E+00
 GRADIENT:   1.2299E+00 -3.4117E+00 -2.0877E+00 -1.8969E+00  3.0106E+00 -4.3268E-01  1.9936E-01  1.7853E-02  1.7464E+00  2.5347E-01
            -9.7547E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1181.57479403148        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.0231E+00  1.0466E+00  1.6041E+00  1.0209E+00  1.2019E+00  8.8762E-01  6.8577E-01  1.0000E-02  7.5394E-01  1.4700E-01
             4.7848E+00
 PARAMETER:  1.2281E-01  1.4554E-01  5.7253E-01  1.2067E-01  2.8391E-01 -1.9214E-02 -2.7721E-01 -4.5785E+00 -1.8244E-01 -1.8173E+00
             1.6654E+00
 GRADIENT:  -6.2197E-01  8.9553E-01 -4.1576E-01  1.2050E+00  1.2452E+00 -3.8363E-01  3.4807E-01  0.0000E+00  2.3828E-01  2.1688E-01
             1.6770E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1181.71296472472        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  1.0222E+00  9.2291E-01  1.5950E+00  1.0990E+00  1.1484E+00  8.8809E-01  6.4525E-01  3.7967E-02  7.3714E-01  4.5390E-02
             4.7864E+00
 PARAMETER:  1.2195E-01  1.9776E-02  5.6687E-01  1.9436E-01  2.3833E-01 -1.8683E-02 -3.3811E-01 -3.1710E+00 -2.0498E-01 -2.9925E+00
             1.6658E+00
 GRADIENT:  -1.0207E+00  2.9951E+00  1.6415E-01  4.3817E+00 -1.5553E+00 -3.0180E-01 -3.0632E-01  2.6610E-03 -6.8538E-01  2.0392E-02
            -5.3758E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1181.75397239151        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  1.0216E+00  8.5963E-01  1.6671E+00  1.1380E+00  1.1439E+00  8.8847E-01  7.1841E-01  6.3658E-02  7.1023E-01  2.4230E-02
             4.7819E+00
 PARAMETER:  1.2134E-01 -5.1252E-02  6.1107E-01  2.2930E-01  2.3442E-01 -1.8256E-02 -2.3072E-01 -2.6542E+00 -2.4217E-01 -3.6202E+00
             1.6648E+00
 GRADIENT:  -1.7638E-02  1.6623E-01  9.1612E-02  9.7460E-02 -2.0620E-01  1.8636E-02  6.3246E-03  7.2869E-03 -4.7855E-02  5.9626E-03
            -4.5240E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1181.76038821019        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  1.0214E+00  8.4025E-01  1.6555E+00  1.1505E+00  1.1354E+00  8.8848E-01  7.3462E-01  6.6529E-02  7.0372E-01  2.0291E-02
             4.7808E+00
 PARAMETER:  1.2114E-01 -7.4060E-02  6.0412E-01  2.4023E-01  2.2696E-01 -1.8245E-02 -2.0840E-01 -2.6101E+00 -2.5137E-01 -3.7976E+00
             1.6646E+00
 GRADIENT:  -3.3632E-01  3.6848E-01  1.0457E-02  6.3199E-01 -1.1644E-01 -4.4584E-02  1.1941E-02  8.3817E-03 -6.5108E-02  4.2181E-03
            -2.0107E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1181.77395873043        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  1.0210E+00  8.0544E-01  1.6432E+00  1.1730E+00  1.1219E+00  8.8845E-01  7.6001E-01  1.6025E-02  6.9159E-01  1.9037E-02
             4.7794E+00
 PARAMETER:  1.2078E-01 -1.1636E-01  5.9662E-01  2.5955E-01  2.1498E-01 -1.8276E-02 -1.7443E-01 -4.0336E+00 -2.6876E-01 -3.8614E+00
             1.6643E+00
 GRADIENT:  -8.4716E-01  7.1722E-01 -1.0730E-01  1.5487E+00 -1.0565E-02 -1.5712E-01 -2.2693E-02  5.2636E-04 -2.1385E-01  3.7505E-03
            -5.7431E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1181.77931306358        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  1.0211E+00  7.9698E-01  1.6664E+00  1.1787E+00  1.1228E+00  8.8857E-01  7.5226E-01  1.5098E-02  6.9285E-01  1.7571E-02
             4.7802E+00
 PARAMETER:  1.2085E-01 -1.2693E-01  6.1069E-01  2.6441E-01  2.1587E-01 -1.8137E-02 -1.8468E-01 -4.0932E+00 -2.6694E-01 -3.9415E+00
             1.6645E+00
 GRADIENT:  -2.7144E+00  7.0719E-01  9.2805E-02  6.3873E-01 -5.3037E-01 -1.9939E-01 -3.4660E-02  4.5287E-04 -1.7123E-01  3.1453E-03
            -1.6457E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1181.83890278626        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  1.0199E+00  6.0794E-01  1.7455E+00  1.2983E+00  1.0840E+00  8.8920E-01  8.8687E-01  1.0000E-02  6.4644E-01  1.0000E-02
             4.7792E+00
 PARAMETER:  1.1969E-01 -3.9768E-01  6.5704E-01  3.6108E-01  1.8062E-01 -1.7438E-02 -2.0058E-02 -6.2382E+00 -3.3627E-01 -5.8422E+00
             1.6643E+00
 GRADIENT:  -9.8138E-01 -8.3942E-01 -1.8476E-01 -2.9815E+00  9.0296E-01  1.0198E-01  4.6365E-02  0.0000E+00  1.0120E-01  0.0000E+00
            -8.8606E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1181.87398584440        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1187
 NPARAMETR:  1.0199E+00  4.8192E-01  1.7458E+00  1.3811E+00  1.0452E+00  8.8887E-01  1.0049E+00  1.0000E-02  6.1967E-01  1.0000E-02
             4.7890E+00
 PARAMETER:  1.1970E-01 -6.2998E-01  6.5720E-01  4.2285E-01  1.4424E-01 -1.7806E-02  1.0490E-01 -7.6197E+00 -3.7857E-01 -7.6685E+00
             1.6663E+00
 GRADIENT:   9.9574E-02  8.0674E-01 -2.5225E-02  2.8897E+00 -4.5879E-01 -2.1552E-01 -4.1490E-02  0.0000E+00 -1.8735E-01  0.0000E+00
             3.8947E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1181.89747828040        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  1.0186E+00  3.8182E-01  1.8260E+00  1.4458E+00  1.0336E+00  8.8876E-01  1.1126E+00  1.0000E-02  6.0247E-01  1.0000E-02
             4.7855E+00
 PARAMETER:  1.1843E-01 -8.6281E-01  7.0212E-01  4.6869E-01  1.3304E-01 -1.7928E-02  2.0673E-01 -8.3954E+00 -4.0672E-01 -9.7357E+00
             1.6656E+00
 GRADIENT:  -3.9885E-01  7.6815E-01  2.4419E-01  3.5486E+00 -1.0288E+00 -1.7711E-01 -4.1535E-02  0.0000E+00 -1.3963E-01  0.0000E+00
             2.8443E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1181.92236904260        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1537
 NPARAMETR:  1.0170E+00  2.7354E-01  1.8566E+00  1.5131E+00  1.0110E+00  8.8896E-01  1.3970E+00  1.0000E-02  5.8249E-01  1.0000E-02
             4.7800E+00
 PARAMETER:  1.1682E-01 -1.1963E+00  7.1873E-01  5.1413E-01  1.1092E-01 -1.7706E-02  4.3431E-01 -9.0080E+00 -4.4045E-01 -1.2874E+01
             1.6644E+00
 GRADIENT:  -8.9480E-01  3.8293E-01  1.7455E-01  2.4277E+00 -7.1035E-01 -7.2001E-02 -7.6568E-03  0.0000E+00 -1.0885E-01  0.0000E+00
            -4.8648E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1181.93640646797        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1713
 NPARAMETR:  1.0161E+00  1.9998E-01  1.8666E+00  1.5579E+00  9.9358E-01  8.8892E-01  1.7116E+00  1.0000E-02  5.7255E-01  1.0000E-02
             4.7787E+00
 PARAMETER:  1.1592E-01 -1.5095E+00  7.2414E-01  5.4335E-01  9.3555E-02 -1.7748E-02  6.3740E-01 -9.1970E+00 -4.5766E-01 -1.5975E+01
             1.6642E+00
 GRADIENT:  -6.4162E-01  1.4959E-01  1.5176E-01  1.0315E+00 -4.7347E-01 -2.4276E-02  2.1953E-02  0.0000E+00  1.1296E-01  0.0000E+00
            -1.8383E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1181.94540737416        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1888
 NPARAMETR:  1.0158E+00  1.4939E-01  1.9420E+00  1.5921E+00  9.9904E-01  8.8869E-01  1.9219E+00  1.0000E-02  5.6352E-01  1.0000E-02
             4.7803E+00
 PARAMETER:  1.1564E-01 -1.8012E+00  7.6371E-01  5.6505E-01  9.9037E-02 -1.8010E-02  7.5333E-01 -9.0834E+00 -4.7355E-01 -1.8993E+01
             1.6645E+00
 GRADIENT:  -5.7457E-02  3.3821E-02 -1.2828E-02  1.8649E-01  1.0893E-01 -7.1624E-03  4.1690E-03  0.0000E+00  2.0688E-02  0.0000E+00
            -1.8164E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1181.95224807945        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2065
 NPARAMETR:  1.0149E+00  9.2151E-02  1.9511E+00  1.6285E+00  9.8542E-01  8.8868E-01  2.1788E+00  1.0000E-02  5.5828E-01  1.0000E-02
             4.7789E+00
 PARAMETER:  1.1484E-01 -2.2843E+00  7.6841E-01  5.8765E-01  8.5310E-02 -1.8018E-02  8.7877E-01 -8.8479E+00 -4.8290E-01 -2.4057E+01
             1.6642E+00
 GRADIENT:  -6.1933E-01  1.1498E-01  1.3208E-01  2.6357E+00 -6.1882E-01 -4.9333E-02 -1.4383E-02  0.0000E+00 -6.7991E-02  0.0000E+00
            -2.4483E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1181.95731273815        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2244
 NPARAMETR:  1.0144E+00  5.3694E-02  1.9572E+00  1.6517E+00  9.7704E-01  8.8865E-01  2.5388E+00  1.0000E-02  5.5363E-01  1.0000E-02
             4.7784E+00
 PARAMETER:  1.1435E-01 -2.8244E+00  7.7151E-01  6.0179E-01  7.6769E-02 -1.8047E-02  1.0317E+00 -8.2332E+00 -4.9126E-01 -2.9851E+01
             1.6641E+00
 GRADIENT:  -3.8050E-01  4.5263E-02  1.2803E-01  1.6887E+00 -5.3798E-01 -1.6877E-02 -1.0930E-02  0.0000E+00 -8.2244E-02  0.0000E+00
            -1.7850E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1181.96121777621        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2419
 NPARAMETR:  1.0142E+00  2.6401E-02  1.9831E+00  1.6697E+00  9.7638E-01  8.8852E-01  3.4105E+00  1.0000E-02  5.4912E-01  1.0000E-02
             4.7806E+00
 PARAMETER:  1.1411E-01 -3.5344E+00  7.8465E-01  6.1263E-01  7.6102E-02 -1.8195E-02  1.3269E+00 -7.0341E+00 -4.9945E-01 -3.7575E+01
             1.6646E+00
 GRADIENT:  -3.3138E-01  1.9918E-02  7.7610E-02  1.5971E+00 -3.6987E-01 -1.7333E-02 -5.0572E-03  0.0000E+00 -2.4020E-02  0.0000E+00
             2.9850E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1181.96382969465        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2595
 NPARAMETR:  1.0140E+00  1.3397E-02  1.9849E+00  1.6767E+00  9.7422E-01  8.8845E-01  5.0418E+00  1.0000E-02  5.4692E-01  1.0000E-02
             4.7784E+00
 PARAMETER:  1.1391E-01 -4.2127E+00  7.8555E-01  6.1681E-01  7.3878E-02 -1.8274E-02  1.7178E+00 -5.6578E+00 -5.0345E-01 -4.5006E+01
             1.6641E+00
 GRADIENT:   1.0923E-01 -2.4626E-03 -2.7126E-02 -5.8476E-01  1.3347E-01  1.8953E-04 -2.3902E-03  0.0000E+00 -1.8635E-02  0.0000E+00
             2.1881E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1181.96444087252        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2772
 NPARAMETR:  1.0140E+00  1.0000E-02  1.9840E+00  1.6790E+00  9.7300E-01  8.8849E-01  6.0321E+00  1.0000E-02  5.4667E-01  1.0000E-02
             4.7782E+00
 PARAMETER:  1.1386E-01 -4.5113E+00  7.8510E-01  6.1817E-01  7.2634E-02 -1.8229E-02  1.8971E+00 -5.0318E+00 -5.0391E-01 -4.8281E+01
             1.6641E+00
 GRADIENT:   6.4183E-03  0.0000E+00 -1.0442E-03 -2.5274E-02  3.7579E-03  8.8069E-04 -1.7203E-03  0.0000E+00 -3.3807E-03  0.0000E+00
            -5.4974E-03

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1181.96467897981        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2950
 NPARAMETR:  1.0140E+00  1.0000E-02  1.9829E+00  1.6789E+00  9.7275E-01  8.8849E-01  7.3273E+00  1.0161E-02  5.4643E-01  1.0000E-02
             4.7782E+00
 PARAMETER:  1.1386E-01 -4.7447E+00  7.8454E-01  6.1815E-01  7.2374E-02 -1.8231E-02  2.0916E+00 -4.4892E+00 -5.0435E-01 -5.0846E+01
             1.6641E+00
 GRADIENT:  -3.6988E-03  0.0000E+00 -1.1009E-03  3.2842E-03 -3.1990E-04 -1.7081E-03 -4.7002E-04  2.8627E-04 -5.0695E-03  0.0000E+00
            -1.2140E-02

0ITERATION NO.:  119    OBJECTIVE VALUE:  -1181.96469203686        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     3090
 NPARAMETR:  1.0140E+00  1.0000E-02  1.9828E+00  1.6789E+00  9.7288E-01  8.8856E-01  7.4957E+00  1.0000E-02  5.4641E-01  1.0000E-02
             4.7781E+00
 PARAMETER:  1.1386E-01 -4.7729E+00  7.8464E-01  6.1816E-01  7.2407E-02 -1.8226E-02  2.1153E+00 -4.5223E+00 -5.0436E-01 -5.1139E+01
             1.6641E+00
 GRADIENT:  -1.8335E-04  0.0000E+00  1.5170E-03  6.6225E-03 -1.2341E-02 -8.4952E-03  1.3903E-04  0.0000E+00  6.0770E-04  0.0000E+00
             5.5869E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3090
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8164E-03 -1.0056E-03  5.8458E-05 -1.6277E-02 -9.4433E-06
 SE:             2.7943E-02  1.1290E-03  3.8179E-05  2.0366E-02  1.0066E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4817E-01  3.7311E-01  1.2573E-01  4.2416E-01  9.2526E-01

 ETASHRINKSD(%)  6.3887E+00  9.6218E+01  9.9872E+01  3.1772E+01  9.9663E+01
 ETASHRINKVR(%)  1.2369E+01  9.9857E+01  1.0000E+02  5.3450E+01  9.9999E+01
 EBVSHRINKSD(%)  6.3540E+00  9.6228E+01  9.9827E+01  3.2023E+01  9.9592E+01
 EBVSHRINKVR(%)  1.2304E+01  9.9858E+01  1.0000E+02  5.3792E+01  9.9998E+01
 RELATIVEINF(%)  5.1823E+01  1.6606E-04  8.1621E-06  6.0258E-02  4.5175E-05
 EPSSHRINKSD(%)  1.5826E+01
 EPSSHRINKVR(%)  2.9148E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1181.9646920368632     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -446.81386547312502     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    38.91
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1181.965       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.00E-02  1.98E+00  1.68E+00  9.73E-01  8.88E-01  7.50E+00  1.00E-02  5.46E-01  1.00E-02  4.78E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.26E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.79E+00  0.00E+00  3.45E+00
 
 TH 4
+       -1.54E+02  0.00E+00  5.42E+00  5.48E+02
 
 TH 5
+       -4.28E+00  0.00E+00 -2.12E+01 -1.51E+02  1.63E+02
 
 TH 6
+        9.29E+01  0.00E+00  8.29E+00 -3.14E+01  7.40E-03  2.75E+02
 
 TH 7
+        2.49E-02  0.00E+00 -5.19E-03 -1.16E-02  1.92E-02 -8.54E-02  3.12E-05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.21E+01  0.00E+00 -1.76E+01 -2.56E+01  1.21E+02  4.35E+01  2.33E-03  0.00E+00  1.19E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.28E+01  0.00E+00 -1.95E+00 -1.19E+01  1.66E+01  9.97E+00 -1.81E-03  0.00E+00  1.43E+01  0.00E+00  2.28E+00
 
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
+        1.22E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.30E+00  0.00E+00  7.52E+00
 
 TH 4
+       -1.80E+02  0.00E+00  2.84E+00  5.80E+02
 
 TH 5
+       -6.28E+00  0.00E+00 -3.42E+01 -1.46E+02  2.17E+02
 
 TH 6
+       -1.13E+01  0.00E+00  3.20E+00 -2.86E+01  1.91E+01  2.34E+02
 
 TH 7
+        5.59E-02  0.00E+00  6.64E-04 -9.01E-03  3.82E-02 -7.09E-02  2.04E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.99E+00  0.00E+00  9.01E-01 -4.03E+01  2.92E+01  1.18E+01 -1.47E-02  0.00E+00  1.37E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.87E+01  0.00E+00  1.55E-01 -1.33E+01  3.34E+00  3.80E+00  3.73E-03  0.00E+00  1.61E+01  0.00E+00  2.11E+01
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.20E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.00E-01  0.00E+00  4.62E+00
 
 TH 4
+       -1.92E+02  0.00E+00  6.08E-01  6.14E+02
 
 TH 5
+        2.04E+01  0.00E+00 -2.13E+01 -1.37E+02  1.34E+02
 
 TH 6
+       -1.06E+02  0.00E+00  3.83E+00 -1.87E+01 -8.19E+00  1.89E+02
 
 TH 7
+       -1.18E-02  0.00E+00 -4.51E-04 -1.50E-02  1.04E-02  3.03E-03  7.03E-06
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.27E+01  0.00E+00  2.52E-02 -5.46E+01  3.05E+01  1.47E+01  2.80E-02  0.00E+00  1.16E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        8.64E+01  0.00E+00  6.71E-01 -5.87E+01  1.04E+01  1.26E+01  3.20E-03  0.00E+00  1.20E+01  0.00E+00  6.85E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       46.092
Stop Time:
Sat Sep 18 10:26:35 CDT 2021
