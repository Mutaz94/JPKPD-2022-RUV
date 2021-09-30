Wed Sep 29 08:47:17 CDT 2021
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
$DATA ../../../../data/int/D/dat41.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   28922.5809020732        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.6887E+02  4.7076E+02 -5.1052E+01  2.6920E+02  1.4648E+02 -1.5970E+03 -1.0342E+03 -1.0663E+02 -1.6041E+03 -4.5052E+02
            -6.0102E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -875.054849496831        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1503E+00  1.8189E+00  9.1474E-01  2.0297E+00  9.0916E-01  2.8492E+00  4.0833E+00  1.0034E+00  2.4450E+00  1.5367E+00
             1.3451E+01
 PARAMETER:  2.3999E-01  6.9824E-01  1.0882E-02  8.0787E-01  4.7671E-03  1.1471E+00  1.5069E+00  1.0343E-01  9.9404E-01  5.2961E-01
             2.6990E+00
 GRADIENT:  -4.9271E+01  3.4139E+01 -3.6168E+01  1.1714E+02 -9.6645E+00  9.1961E+01  5.3291E+01  4.3098E+00  3.6520E+01  3.5342E+01
             5.2964E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -951.787842940255        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.1626E+00  1.7715E+00  2.2641E+01  2.6492E+00  2.3329E+00  1.8743E+00  1.1363E+01  7.4677E-01  2.0495E+00  1.0318E+00
             1.3420E+01
 PARAMETER:  2.5068E-01  6.7183E-01  3.2198E+00  1.0742E+00  9.4713E-01  7.2823E-01  2.5304E+00 -1.9199E-01  8.1758E-01  1.3135E-01
             2.6967E+00
 GRADIENT:  -3.7937E+01  2.9400E+01 -7.5099E+00  7.9760E+01  1.1461E+01  5.3931E+01  8.2861E+01  6.5340E-02  1.8787E+01  1.3274E+01
             5.2624E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1116.78300606606        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1104E+00  1.9779E+00  5.0325E+00  6.5003E-01  1.9479E+00  2.1668E+00  2.3186E+00  4.0611E+00  2.0030E+00  5.8213E-01
             9.7104E+00
 PARAMETER:  2.0468E-01  7.8201E-01  1.7159E+00 -3.3074E-01  7.6673E-01  8.7327E-01  9.4095E-01  1.5015E+00  7.9464E-01 -4.4105E-01
             2.3732E+00
 GRADIENT:  -3.5709E+01 -2.5772E+01  3.5508E+00 -3.6922E+00 -2.5115E+01  3.1103E+01 -6.7059E+01  1.1607E+00  6.0533E+00  5.1036E+00
             3.2274E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1150.30812818400        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.1208E+00  2.0944E+00  5.6050E+00  4.5264E-01  2.1226E+00  1.9825E+00  2.6095E+00  6.4969E+00  1.7415E+00  4.1616E-01
             8.1720E+00
 PARAMETER:  2.1404E-01  8.3927E-01  1.8237E+00 -6.9267E-01  8.5266E-01  7.8434E-01  1.0591E+00  1.9713E+00  6.5472E-01 -7.7668E-01
             2.2007E+00
 GRADIENT:  -5.1164E-01 -4.6953E+00  1.4248E+00 -7.2672E+00 -1.0888E+01  1.9548E+00  6.6182E-01  9.9849E+00  4.6548E+00  2.4144E+00
             2.1196E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1158.22026116270        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.1146E+00  1.9801E+00  2.8096E+00  5.9242E-01  1.9498E+00  1.9990E+00  2.7965E+00  2.8002E+00  1.6413E+00  1.6949E-01
             8.1290E+00
 PARAMETER:  2.0851E-01  7.8313E-01  1.1330E+00 -4.2355E-01  7.6772E-01  7.9263E-01  1.1284E+00  1.1297E+00  5.9547E-01 -1.6750E+00
             2.1954E+00
 GRADIENT:  -4.0110E+00  8.6340E+00  2.3724E+00 -3.3993E+00  9.9793E-01  5.7176E+00  3.8839E+00  1.1864E+00  4.9869E+00  5.2189E-01
             1.1487E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1159.28113986710        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.1238E+00  1.7926E+00  1.6100E+00  6.7297E-01  1.6941E+00  1.9664E+00  2.9067E+00  1.4232E+00  1.3759E+00  1.0457E-01
             8.0443E+00
 PARAMETER:  2.1674E-01  6.8369E-01  5.7624E-01 -2.9606E-01  6.2713E-01  7.7619E-01  1.1670E+00  4.5293E-01  4.1910E-01 -2.1579E+00
             2.1850E+00
 GRADIENT:   1.4073E+00  1.8755E+00 -1.7888E+00  3.3978E+00  4.1072E+00  1.0815E+00  4.0382E+00  1.0987E+00  5.0056E+00  2.4197E-01
            -1.3116E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1159.56235727196        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  1.1265E+00  1.7538E+00  1.4331E+00  6.7155E-01  1.6412E+00  1.9486E+00  2.9022E+00  1.0730E+00  1.1091E+00  8.3388E-02
             8.1139E+00
 PARAMETER:  2.1916E-01  6.6179E-01  4.5984E-01 -2.9816E-01  5.9543E-01  7.6710E-01  1.1655E+00  1.7047E-01  2.0357E-01 -2.3843E+00
             2.1936E+00
 GRADIENT:   1.4622E+00 -2.3934E+00 -9.6687E-01 -5.7488E-01 -2.7767E-01 -1.4955E+00  1.8882E-02  4.6270E-01  9.0358E-01  1.5184E-01
             6.7986E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1163.76461851630        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.1579E+00  1.8672E+00  2.0564E+00  7.0723E-01  1.7808E+00  2.0536E+00  3.2615E+00  9.6384E-01  9.2460E-01  5.0098E-02
             8.3499E+00
 PARAMETER:  2.4663E-01  7.2446E-01  8.2097E-01 -2.4640E-01  6.7707E-01  8.1957E-01  1.2822E+00  6.3171E-02  2.1610E-02 -2.8938E+00
             2.2223E+00
 GRADIENT:  -4.9957E+00 -9.0155E-01  6.6782E-01 -3.4697E+00 -6.9497E-01 -1.1909E+00  4.5545E-02  4.8965E-02  2.6736E-01  4.2385E-02
             6.2469E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1163.88494641186        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  1.1692E+00  1.8018E+00  2.1210E+00  7.4880E-01  1.7547E+00  2.0587E+00  3.3520E+00  8.5754E-01  9.3356E-01  4.5951E-02
             8.3363E+00
 PARAMETER:  2.5629E-01  6.8881E-01  8.5187E-01 -1.8928E-01  6.6228E-01  8.2209E-01  1.3096E+00 -5.3682E-02  3.1248E-02 -2.9802E+00
             2.2206E+00
 GRADIENT:  -4.4201E-02  1.2065E-01 -1.1006E-01  3.5451E-01  2.6952E-01  1.6942E-01  2.2320E-01  2.9648E-02 -8.6005E-02  3.6015E-02
            -3.5904E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1163.90237105677        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      980
 NPARAMETR:  1.1686E+00  1.7976E+00  2.1206E+00  7.4933E-01  1.7537E+00  2.0576E+00  3.3473E+00  7.7276E-01  9.4469E-01  1.0000E-02
             8.3341E+00
 PARAMETER:  2.5582E-01  6.8646E-01  8.5168E-01 -1.8857E-01  6.6174E-01  8.2155E-01  1.3082E+00 -1.5778E-01  4.3100E-02 -4.9960E+00
             2.2204E+00
 GRADIENT:  -2.1178E-01 -1.9766E-01  5.3476E-02  1.0314E-01  2.2756E-02 -6.6352E-02 -2.7519E-01 -4.5059E-03 -1.1138E-02  0.0000E+00
            -6.9079E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1163.90393439295        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1158
 NPARAMETR:  1.1697E+00  1.8014E+00  2.1143E+00  7.4824E-01  1.7538E+00  2.0602E+00  3.3545E+00  7.8060E-01  9.4284E-01  1.0000E-02
             8.3365E+00
 PARAMETER:  2.5673E-01  6.8856E-01  8.4874E-01 -1.9003E-01  6.6180E-01  8.2279E-01  1.3103E+00 -1.4769E-01  4.1141E-02 -4.6500E+00
             2.2206E+00
 GRADIENT:   2.2680E-01  1.6364E-02 -2.6788E-02 -1.2283E-01  1.0275E-01  4.1181E-01  5.8560E-01  2.0884E-03  2.9221E-02  0.0000E+00
            -1.2345E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1163.90405673711        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1345
 NPARAMETR:  1.1695E+00  1.8006E+00  2.1148E+00  7.4855E-01  1.7536E+00  2.0601E+00  3.3552E+00  7.7511E-01  9.4065E-01  1.0000E-02
             8.3367E+00
 PARAMETER:  2.5661E-01  6.8811E-01  8.4898E-01 -1.8962E-01  6.6165E-01  8.2278E-01  1.3105E+00 -1.5475E-01  3.8818E-02 -4.6500E+00
             2.2207E+00
 GRADIENT:   1.6300E-01  1.0292E-02 -1.4922E-02 -9.1139E-02  5.8019E-02  4.0874E-01  5.3899E-01  8.4438E-04  7.1884E-03  0.0000E+00
            -1.4517E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1163.90414788571        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1533
 NPARAMETR:  1.1696E+00  1.8001E+00  2.1159E+00  7.4904E-01  1.7531E+00  2.0604E+00  3.3563E+00  7.7649E-01  9.3920E-01  1.0000E-02
             8.3363E+00
 PARAMETER:  2.5665E-01  6.8782E-01  8.4950E-01 -1.8897E-01  6.6141E-01  8.2289E-01  1.3108E+00 -1.5297E-01  3.7275E-02 -4.6500E+00
             2.2206E+00
 GRADIENT:   1.9702E-01  6.5912E-02 -3.5536E-03 -2.3530E-02 -1.9276E-02  4.6236E-01  5.1369E-01  9.9765E-04 -1.5082E-02  0.0000E+00
            -3.2384E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1163.90421700070        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1721             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1696E+00  1.7989E+00  2.1168E+00  7.4939E-01  1.7530E+00  2.0603E+00  3.3575E+00  7.7223E-01  9.4104E-01  1.0000E-02
             8.3366E+00
 PARAMETER:  2.5662E-01  6.8720E-01  8.4988E-01 -1.8849E-01  6.6131E-01  8.2284E-01  1.3112E+00 -1.5848E-01  3.9230E-02 -4.6500E+00
             2.2207E+00
 GRADIENT:   2.0116E+01  1.6886E+01  3.3636E-01  2.3706E+00  3.9063E+00  2.0089E+01  3.1169E+01  2.5853E-03  9.3217E-02  0.0000E+00
             3.9789E+01

0ITERATION NO.:   72    OBJECTIVE VALUE:  -1163.90421700070        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1778
 NPARAMETR:  1.1696E+00  1.7989E+00  2.1168E+00  7.4939E-01  1.7530E+00  2.0603E+00  3.3575E+00  7.7223E-01  9.4104E-01  1.0000E-02
             8.3366E+00
 PARAMETER:  2.5662E-01  6.8720E-01  8.4988E-01 -1.8849E-01  6.6131E-01  8.2284E-01  1.3112E+00 -1.5848E-01  3.9230E-02 -4.6500E+00
             2.2207E+00
 GRADIENT:   1.7388E-01  9.4291E-03 -1.5437E-02 -6.9160E-02  2.4215E-02  4.3443E-01  5.7001E-01  1.3209E-04  4.3074E-03  0.0000E+00
            -1.8052E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1778
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9969E-03  5.2558E-03 -4.4196E-03 -3.2733E-02 -1.3126E-05
 SE:             2.9209E-02  2.7495E-02  3.1116E-03  1.0162E-02  1.4699E-04
 N:                     100         100         100         100         100

 P VAL.:         7.3216E-01  8.4840E-01  1.5550E-01  1.2779E-03  9.2885E-01

 ETASHRINKSD(%)  2.1473E+00  7.8896E+00  8.9576E+01  6.5954E+01  9.9508E+01
 ETASHRINKVR(%)  4.2486E+00  1.5157E+01  9.8913E+01  8.8409E+01  9.9998E+01
 EBVSHRINKSD(%)  2.8532E+00  5.2345E+00  8.9882E+01  6.9702E+01  9.9445E+01
 EBVSHRINKVR(%)  5.6249E+00  1.0195E+01  9.8976E+01  9.0821E+01  9.9997E+01
 RELATIVEINF(%)  9.4124E+01  3.5725E+01  1.8464E-01  2.5509E+00  6.6430E-04
 EPSSHRINKSD(%)  7.5710E+00
 EPSSHRINKVR(%)  1.4569E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1163.9042170007042     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       490.18514276770657     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    53.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.92
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1163.904       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.17E+00  1.80E+00  2.12E+00  7.49E-01  1.75E+00  2.06E+00  3.36E+00  7.72E-01  9.41E-01  1.00E-02  8.34E+00
 


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
+        1.82E+02
 
 TH 2
+       -1.26E+00  2.73E+01
 
 TH 3
+        1.18E+00  2.23E+00  5.70E+00
 
 TH 4
+       -9.22E+00  3.28E+01 -2.05E+01  2.44E+02
 
 TH 5
+       -5.50E+00 -1.25E+01 -2.06E+01  7.66E+01  1.01E+02
 
 TH 6
+        6.98E+00 -2.71E-01  2.58E-02  1.56E+00 -2.73E+00  4.66E+01
 
 TH 7
+       -1.05E-01  2.06E+00 -1.55E+00 -1.99E+01  3.51E+00 -1.22E+00  1.28E+01
 
 TH 8
+       -1.56E-01 -2.10E-01 -4.41E-01  9.62E-01  1.10E+00  3.89E-02  1.94E-01  2.04E-01
 
 TH 9
+       -2.46E-01 -1.79E+00 -7.36E-01 -1.05E+01  3.32E+00 -5.72E-01  2.60E+00 -4.39E-02  5.24E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.24E+00 -3.07E+00 -9.38E-02 -9.88E+00 -6.87E-02 -4.11E+01  1.11E+00  1.61E-01  1.95E+00  0.00E+00  1.88E+01
 
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
 #CPUT: Total CPU Time in Seconds,       69.409
Stop Time:
Wed Sep 29 08:48:30 CDT 2021
