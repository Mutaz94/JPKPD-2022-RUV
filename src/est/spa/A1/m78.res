Sat Sep 18 09:26:41 CDT 2021
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
$DATA ../../../../data/spa/A1/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1465.68828302923        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0953E+02 -7.8430E+01 -5.2180E+01 -5.7471E+01  8.7202E+01  5.9122E+00 -3.5113E+01  1.5283E+01 -2.8221E+01 -1.6377E+01
            -3.1417E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1534.41529903242        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.3981E-01  1.0957E+00  1.2293E+00  1.0346E+00  1.0158E+00  9.4480E-01  1.1700E+00  8.7824E-01  1.1058E+00  9.6273E-01
             1.7331E+00
 PARAMETER:  3.7921E-02  1.9137E-01  3.0645E-01  1.3399E-01  1.1564E-01  4.3222E-02  2.5698E-01 -2.9837E-02  2.0061E-01  6.2014E-02
             6.4990E-01
 GRADIENT:  -7.6130E+01  3.4919E+01  3.2383E+00  2.9926E+01 -3.6894E+01 -1.5907E+01  6.4773E-01  6.1086E+00  5.8812E+00 -2.6763E-01
             4.6953E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1540.39666458995        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.5929E-01  9.6278E-01  1.4920E+00  1.1579E+00  1.1085E+00  9.3829E-01  9.6633E-01  2.6181E-01  1.1543E+00  1.1545E+00
             1.5916E+00
 PARAMETER:  5.8443E-02  6.2072E-02  5.0015E-01  2.4662E-01  2.0302E-01  3.6308E-02  6.5751E-02 -1.2401E+00  2.4352E-01  2.4365E-01
             5.6472E-01
 GRADIENT:  -1.9297E+01  3.6912E+01 -9.5154E+00  7.5924E+01  7.1182E+00 -1.6765E+01 -5.7090E+00  3.7930E-01  5.3094E+00  2.0215E+00
             1.3955E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1544.95277886051        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.6636E-01  7.7864E-01  1.5112E+00  1.2178E+00  1.0325E+00  9.8345E-01  1.5684E+00  2.2853E-01  9.8061E-01  1.0873E+00
             1.5275E+00
 PARAMETER:  6.5784E-02 -1.5021E-01  5.1288E-01  2.9701E-01  1.3200E-01  8.3316E-02  5.5007E-01 -1.3761E+00  8.0417E-02  1.8370E-01
             5.2361E-01
 GRADIENT:   8.6237E+00  1.0776E+01  7.4251E+00  9.4311E+00 -9.7276E+00  3.1755E+00  6.5566E+00  2.7938E-01  3.3433E+00 -2.1725E+00
            -2.0540E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1545.89256377571        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.6099E-01  6.5520E-01  1.3428E+00  1.2785E+00  9.4641E-01  9.7022E-01  1.5886E+00  8.3193E-02  9.3201E-01  1.0466E+00
             1.5252E+00
 PARAMETER:  6.0212E-02 -3.2281E-01  3.9479E-01  3.4570E-01  4.4918E-02  6.9766E-02  5.6283E-01 -2.3866E+00  2.9591E-02  1.4552E-01
             5.2209E-01
 GRADIENT:  -3.7622E+00  1.1529E+00 -4.1466E+00  4.9699E+00  3.8903E+00 -1.8793E+00 -1.7372E+00  6.3807E-02 -3.4858E+00  1.5074E+00
             7.4936E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1547.34084954103        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      417
 NPARAMETR:  9.6528E-01  4.1246E-01  1.3074E+00  1.4452E+00  8.4540E-01  9.7398E-01  2.2118E+00  1.3782E-02  8.8160E-01  9.8662E-01
             1.5149E+00
 PARAMETER:  6.4667E-02 -7.8561E-01  3.6806E-01  4.6823E-01 -6.7948E-02  7.3636E-02  8.9382E-01 -4.1844E+00 -2.6013E-02  8.6527E-02
             5.1532E-01
 GRADIENT:  -5.1034E+00  8.2286E+00  4.6496E+00  2.3263E+01 -1.2945E+01 -1.0160E+00  1.5002E-01  1.9760E-03 -1.9965E+00 -1.1685E+00
            -1.4401E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1548.34769077049        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      593
 NPARAMETR:  9.6611E-01  2.4826E-01  1.4648E+00  1.5451E+00  8.6720E-01  9.7343E-01  2.8513E+00  1.0000E-02  8.5966E-01  1.0479E+00
             1.5134E+00
 PARAMETER:  6.5520E-02 -1.2933E+00  4.8172E-01  5.3512E-01 -4.2483E-02  7.3067E-02  1.1478E+00 -6.7079E+00 -5.1218E-02  1.4678E-01
             5.1435E-01
 GRADIENT:   2.5604E+00  1.8575E+00 -8.6198E-02  1.1335E+01  2.1005E-01 -4.1220E-01 -1.3648E-01  0.0000E+00  2.5349E-01 -3.2278E-01
             1.5423E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1548.68601535954        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      768
 NPARAMETR:  9.6295E-01  1.1952E-01  1.4229E+00  1.6139E+00  8.1892E-01  9.7195E-01  4.1384E+00  1.0000E-02  8.3917E-01  1.0357E+00
             1.5056E+00
 PARAMETER:  6.2243E-02 -2.0242E+00  4.5272E-01  5.7863E-01 -9.9775E-02  7.1548E-02  1.5203E+00 -1.0913E+01 -7.5341E-02  1.3512E-01
             5.0918E-01
 GRADIENT:  -3.1268E-01  8.0054E-01  8.0440E-01  6.3811E+00 -2.9218E+00 -3.1400E-01  4.1208E-01  0.0000E+00 -5.5319E-01  2.3432E-01
            -2.0227E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1548.95609412609        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      951
 NPARAMETR:  9.5970E-01  3.8043E-02  1.5740E+00  1.6700E+00  8.4551E-01  9.7468E-01  7.1623E+00  1.0000E-02  8.1376E-01  1.0901E+00
             1.5011E+00
 PARAMETER:  5.8867E-02 -3.1690E+00  5.5364E-01  6.1281E-01 -6.7819E-02  7.4358E-02  2.0688E+00 -1.7630E+01 -1.0609E-01  1.8626E-01
             5.0619E-01
 GRADIENT:  -6.0589E+00  3.0500E+00  3.3475E+00  3.8029E-01 -7.5034E+00  1.9086E+00  5.9893E+00  0.0000E+00 -1.2900E+00 -3.3957E+00
            -2.4946E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1550.17206266657        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1127
 NPARAMETR:  9.6035E-01  1.4575E-02  1.4712E+00  1.6678E+00  8.1360E-01  9.7142E-01  9.9004E+00  1.0000E-02  8.1321E-01  1.0563E+00
             1.5015E+00
 PARAMETER:  5.9538E-02 -4.1285E+00  4.8610E-01  6.1153E-01 -1.0629E-01  7.0999E-02  2.3926E+00 -2.3805E+01 -1.0677E-01  1.5476E-01
             5.0650E-01
 GRADIENT:  -2.5283E+00 -9.4392E-01  1.2265E+00 -5.7525E+00 -8.2828E-01  2.1253E-02 -2.9402E+00  0.0000E+00 -1.8691E+00  2.1980E+00
            -7.0768E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1550.19055041758        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1322             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6247E-01  1.3922E-02  1.4707E+00  1.6681E+00  8.1451E-01  9.7135E-01  1.0105E+01  1.0000E-02  8.1346E-01  1.0563E+00
             1.5017E+00
 PARAMETER:  6.1751E-02 -4.1743E+00  4.8574E-01  6.1168E-01 -1.0516E-01  7.0927E-02  2.4130E+00 -2.4102E+01 -1.0645E-01  1.5473E-01
             5.0660E-01
 GRADIENT:   2.0299E+01 -5.9815E-01  5.3196E-01  4.3704E+01  1.4956E+00  1.4028E+00  2.9754E-01  0.0000E+00 -6.9581E-01  2.1499E+00
             3.7600E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1550.19299471660        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1458
 NPARAMETR:  9.6140E-01  1.3924E-02  1.4707E+00  1.6681E+00  8.1350E-01  9.7161E-01  1.0107E+01  1.0000E-02  8.1346E-01  1.0563E+00
             1.5018E+00
 PARAMETER:  6.0638E-02 -4.1741E+00  4.8573E-01  6.1167E-01 -1.0641E-01  7.1196E-02  2.4132E+00 -2.4102E+01 -1.0645E-01  1.5473E-01
             5.0666E-01
 GRADIENT:   3.6949E-02 -9.3297E-01  1.0393E+00 -5.9802E+00 -5.0096E-01  1.0547E-01 -2.8612E+00  0.0000E+00 -1.6970E+00  2.2112E+00
             1.8726E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1550.23598719091        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1556
 NPARAMETR:  9.6032E-01  1.3747E-02  1.4687E+00  1.6671E+00  8.1291E-01  9.7104E-01  1.0268E+01  1.0000E-02  8.1866E-01  1.0486E+00
             1.5020E+00
 PARAMETER:  5.9513E-02 -4.1869E+00  4.8436E-01  6.1111E-01 -1.0714E-01  7.0611E-02  2.4290E+00 -2.4102E+01 -1.0009E-01  1.4745E-01
             5.0681E-01
 GRADIENT:   1.5254E+01  3.2357E-02  1.1347E+00  4.1447E+01  1.3450E+00  1.2801E+00  1.7159E+00  0.0000E+00  8.9906E-01  4.9648E-01
             2.1003E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1550.24050909799        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1649
 NPARAMETR:  9.5915E-01  1.3173E-02  1.4629E+00  1.6652E+00  8.1076E-01  9.7032E-01  1.0464E+01  1.0000E-02  8.2105E-01  1.0433E+00
             1.5024E+00
 PARAMETER:  5.8288E-02 -4.2296E+00  4.8040E-01  6.0995E-01 -1.0978E-01  6.9868E-02  2.4479E+00 -2.4102E+01 -9.7172E-02  1.4242E-01
             5.0708E-01
 GRADIENT:  -5.1364E+00 -3.0936E-01  8.9489E-01 -1.1522E+01  1.3241E+00 -4.1380E-01 -1.4242E+00  0.0000E+00  9.2987E-01 -1.8309E-01
            -7.3404E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1550.32608894647        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1828
 NPARAMETR:  9.6040E-01  1.0048E-02  1.3883E+00  1.6600E+00  7.8418E-01  9.7070E-01  1.1801E+01  1.0000E-02  8.1400E-01  1.0323E+00
             1.4981E+00
 PARAMETER:  5.9590E-02 -4.5004E+00  4.2805E-01  6.0683E-01 -1.4312E-01  7.0265E-02  2.5682E+00 -2.4102E+01 -1.0579E-01  1.3177E-01
             5.0419E-01
 GRADIENT:  -1.7615E+00 -6.0372E-01  3.8063E-01 -1.2280E+01 -8.7387E-01 -1.4451E-01 -1.7686E+00  0.0000E+00 -2.1847E+00  1.7147E+00
            -2.5703E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1550.35557344748        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2009
 NPARAMETR:  9.6086E-01  1.0000E-02  1.3679E+00  1.6630E+00  7.7698E-01  9.7120E-01  1.1837E+01  1.0000E-02  8.1926E-01  1.0224E+00
             1.4989E+00
 PARAMETER:  6.0075E-02 -4.5087E+00  4.1328E-01  6.0860E-01 -1.5234E-01  7.0778E-02  2.5712E+00 -2.4102E+01 -9.9353E-02  1.2217E-01
             5.0475E-01
 GRADIENT:  -6.0732E-01  1.9396E-01 -1.5409E+00 -7.2990E+00  9.8964E-01 -1.0102E-02  3.0532E+00  0.0000E+00 -2.1119E+00 -1.4199E+00
            -3.8487E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1550.38577912743        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2190
 NPARAMETR:  9.6140E-01  1.0000E-02  1.4240E+00  1.6712E+00  7.9670E-01  9.7140E-01  1.1976E+01  1.0000E-02  8.2011E-01  1.0363E+00
             1.5020E+00
 PARAMETER:  6.0637E-02 -4.5140E+00  4.5346E-01  6.1354E-01 -1.2727E-01  7.0980E-02  2.5829E+00 -2.4102E+01 -9.8315E-02  1.3562E-01
             5.0683E-01
 GRADIENT:   3.3743E-01  0.0000E+00 -2.2099E-01  1.3396E+00  1.0069E-01  3.4439E-03 -1.5973E-01  0.0000E+00  1.6668E-01  7.1946E-02
             1.5437E-01

0ITERATION NO.:   82    OBJECTIVE VALUE:  -1550.38624490065        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2247
 NPARAMETR:  9.6128E-01  1.0000E-02  1.4283E+00  1.6709E+00  7.9803E-01  9.7138E-01  1.1981E+01  1.0000E-02  8.1976E-01  1.0374E+00
             1.5018E+00
 PARAMETER:  6.0510E-02 -4.5142E+00  4.5649E-01  6.1338E-01 -1.2560E-01  7.0963E-02  2.5833E+00 -2.4102E+01 -9.8738E-02  1.3674E-01
             5.0665E-01
 GRADIENT:   6.3238E-02  0.0000E+00 -2.7340E-01 -7.0368E-01  3.9201E-01 -3.5050E-03  1.0428E+00  0.0000E+00 -4.3649E-01 -5.0173E-01
            -6.9973E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2247
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.6514E-05  4.6421E-03 -4.4235E-05 -1.0046E-02 -2.7192E-02
 SE:             2.9560E-02  4.7448E-03  1.1941E-04  2.8492E-02  2.3423E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9820E-01  3.2791E-01  7.1105E-01  7.2441E-01  2.4569E-01

 ETASHRINKSD(%)  9.7163E-01  8.4104E+01  9.9600E+01  4.5466E+00  2.1529E+01
 ETASHRINKVR(%)  1.9338E+00  9.7473E+01  9.9998E+01  8.8865E+00  3.8424E+01
 EBVSHRINKSD(%)  9.9779E-01  8.8377E+01  9.9534E+01  4.1597E+00  1.9465E+01
 EBVSHRINKVR(%)  1.9856E+00  9.8649E+01  9.9998E+01  8.1464E+00  3.5142E+01
 RELATIVEINF(%)  9.7919E+01  9.9609E-01  1.9183E-04  5.7458E+01  5.7605E+00
 EPSSHRINKSD(%)  3.8500E+01
 EPSSHRINKVR(%)  6.2178E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1550.3862449006508     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -815.23541833691263     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.08
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1550.386       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  1.00E-02  1.43E+00  1.67E+00  7.98E-01  9.71E-01  1.20E+01  1.00E-02  8.20E-01  1.04E+00  1.50E+00
 


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
+        1.75E-02
 
 TH 2
+        4.66E-18  5.83E-33
 
 TH 3
+        1.40E+00  1.75E-15  1.12E+02
 
 TH 4
+        2.96E-01  3.70E-16  2.37E+01  5.00E+00
 
 TH 5
+       -6.43E+02 -8.05E-13 -5.15E+04 -1.09E+04  2.36E+07
 
 TH 6
+       -5.14E-03 -6.43E-18 -4.12E-01 -8.69E-02  1.89E+02  1.51E-03
 
 TH 7
+        3.51E-03  4.39E-18  2.81E-01  5.93E-02 -1.29E+02 -1.03E-03  7.04E-04
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.86E+02  9.84E-13  6.30E+04  1.33E+04 -2.89E+07 -2.31E+02  1.58E+02  0.00E+00  3.54E+07
 
 TH10
+        4.54E+02  5.68E-13  3.64E+04  7.68E+03 -1.67E+07 -1.33E+02  9.11E+01  0.00E+00  2.04E+07  1.18E+07
 
 TH11
+        8.47E+01  1.06E-13  6.79E+03  1.43E+03 -3.11E+06 -2.49E+01  1.70E+01  0.00E+00  3.81E+06  2.20E+06  4.10E+05
 
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
+        1.29E+03
 
 TH 2
+        0.00E+00  1.78E+03
 
 TH 3
+        4.29E+01  0.00E+00  2.38E+05
 
 TH 4
+        2.89E+00  0.00E+00 -5.52E+02  9.56E+04
 
 TH 5
+       -2.71E+02  0.00E+00  2.31E+03  1.55E+03  1.01E+07
 
 TH 6
+       -1.95E+01  0.00E+00 -1.24E+01 -1.12E+01  6.85E+01  2.53E+02
 
 TH 7
+       -6.50E-01  0.00E+00  1.66E+01  4.71E+01 -5.06E+01  1.44E-01  1.05E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.62E+02  0.00E+00 -8.93E+02 -4.89E+02 -1.23E+07  1.60E+02  1.40E+01  0.00E+00  1.51E+07
 
 TH10
+        1.92E+02  0.00E+00 -1.91E+03 -1.13E+03 -7.12E+06 -6.99E+01  3.54E+01  0.00E+00  8.71E+06  5.03E+06
 
 TH11
+       -5.27E+01  0.00E+00 -1.30E+02 -7.35E+01 -1.33E+06  1.90E+01  1.93E+00  0.00E+00  1.62E+06  9.38E+05  1.75E+05
 
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
+        1.25E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -7.01E+02  0.00E+00  9.80E+04
 
 TH 4
+       -5.88E+02  0.00E+00  6.24E+04  4.03E+04
 
 TH 5
+        4.46E+03  0.00E+00 -6.36E+05 -4.05E+05  4.12E+06
 
 TH 6
+       -1.80E+02  0.00E+00 -1.73E+02 -1.12E+02  1.17E+03  2.02E+02
 
 TH 7
+        1.41E+01  0.00E+00 -2.05E+03 -1.30E+03  1.33E+04  4.01E+00  4.27E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.47E+01  0.00E+00  4.65E+02  1.39E+02 -3.02E+03 -1.17E+01 -9.22E+00  0.00E+00  3.44E+02
 
 TH10
+       -3.19E+03  0.00E+00  4.48E+05  2.85E+05 -2.91E+06 -8.25E+02 -9.36E+03  0.00E+00  2.05E+03  2.05E+06
 
 TH11
+        1.62E+02  0.00E+00 -9.64E+01 -2.80E+02  5.42E+02  2.48E+01  2.17E+00  0.00E+00  2.40E+02 -4.14E+02  3.57E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       34.849
Stop Time:
Sat Sep 18 09:27:17 CDT 2021
