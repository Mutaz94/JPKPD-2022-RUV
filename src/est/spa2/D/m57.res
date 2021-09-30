Thu Sep 30 09:20:31 CDT 2021
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
$DATA ../../../../data/spa2/D/dat57.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   21209.5797403979        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0756E+02  3.8461E+02  3.6776E+01  1.9581E+02  2.0434E+02 -2.2732E+03 -1.0842E+03 -9.1534E+01 -1.6657E+03 -8.8103E+02
            -4.0954E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -637.618354945458        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.6701E+00  1.2620E+00  9.4520E-01  1.5435E+00  1.0345E+00  2.1699E+00  1.4907E+00  9.7218E-01  1.3093E+00  1.1083E+00
             1.4154E+01
 PARAMETER:  6.1287E-01  3.3267E-01  4.3637E-02  5.3404E-01  1.3392E-01  8.7469E-01  4.9922E-01  7.1786E-02  3.6947E-01  2.0285E-01
             2.7500E+00
 GRADIENT:   8.6147E+01 -1.4130E+01 -1.2607E+01  1.2737E+01  2.3412E+01  5.2250E+01 -1.6258E+01  4.7629E+00 -1.7100E+01  1.5105E+01
             1.6400E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -719.054350251320        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.6075E+00  1.7427E+00  1.8511E+00  1.3083E+00  1.6458E+00  2.8443E+00  5.1873E+00  5.0511E-01  4.1860E+00  9.5542E-01
             1.2195E+01
 PARAMETER:  5.7469E-01  6.5544E-01  7.1577E-01  3.6872E-01  5.9820E-01  1.1453E+00  1.7462E+00 -5.8297E-01  1.5317E+00  5.4400E-02
             2.6010E+00
 GRADIENT:   5.0084E+01 -8.2680E+00 -2.9974E+01  2.2806E+01  2.8978E+00  5.6508E+01  7.8716E+01  1.9415E-01  7.5596E+01  8.7508E+00
             2.5526E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -745.029581232281        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      283
 NPARAMETR:  1.5571E+00  7.4107E-01  2.8553E+00  2.0788E+00  9.9455E-01  2.6529E+00  7.9540E+00  5.3098E-01  2.2551E+00  1.4048E+00
             1.2039E+01
 PARAMETER:  5.4284E-01 -1.9966E-01  1.1492E+00  8.3178E-01  9.4540E-02  1.0756E+00  2.1737E+00 -5.3302E-01  9.1321E-01  4.3993E-01
             2.5881E+00
 GRADIENT:   2.8389E+01  1.2599E+01  3.6016E+00  1.7121E+01 -6.6922E+01  2.1144E+01  1.3925E+01  4.1224E-01  5.2708E+01  1.6899E+01
             1.7839E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -823.901948797748        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  1.1920E+00  7.7368E-01  1.9126E+01  1.7354E+00  3.9152E+00  2.1011E+00  6.9069E+00  1.3588E+00  1.0252E+00  1.4836E+00
             1.0382E+01
 PARAMETER:  2.7566E-01 -1.5660E-01  3.0511E+00  6.5124E-01  1.4649E+00  8.4246E-01  2.0325E+00  4.0662E-01  1.2493E-01  4.9449E-01
             2.4401E+00
 GRADIENT:  -2.6707E+00  1.5725E+00 -9.0192E-01  1.1974E+01  8.0524E+00  3.0263E+01 -7.4348E+00  1.8107E-02 -2.9182E+00  4.3955E+00
             4.4515E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -833.173136587108        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  1.1327E+00  5.9591E-01  1.1742E+01  1.7183E+00  2.5127E+00  1.6424E+00  7.6524E+00  1.7145E+00  1.0931E+00  5.5408E-01
             1.0273E+01
 PARAMETER:  2.2457E-01 -4.1766E-01  2.5632E+00  6.4136E-01  1.0214E+00  5.9614E-01  2.1350E+00  6.3911E-01  1.8900E-01 -4.9044E-01
             2.4295E+00
 GRADIENT:   5.0782E-01  2.7128E-01 -2.3214E+00 -9.9924E+00  4.6657E+00  3.7872E+00  9.5642E+00  2.0231E-01  2.2881E+00  1.6684E+00
             2.4381E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -835.380503266152        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.1405E+00  4.0974E-01  2.1690E+01  1.8236E+00  2.5187E+00  1.6297E+00  7.8784E+00  5.5692E-01  1.1734E+00  1.0769E-01
             1.0040E+01
 PARAMETER:  2.3143E-01 -7.9224E-01  3.1768E+00  7.0081E-01  1.0238E+00  5.8843E-01  2.1641E+00 -4.8533E-01  2.5990E-01 -2.1285E+00
             2.4066E+00
 GRADIENT:   1.0363E+01 -2.8286E+00 -9.9601E-01 -1.9088E+00  2.6645E+00  3.3579E+00  6.7243E-02  7.1872E-03 -2.0507E+00  6.6212E-02
            -8.2153E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -837.168406704623        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1000             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1269E+00  3.9559E-01  8.4969E+01  1.8745E+00  2.5508E+00  1.6144E+00  8.7587E+00  5.3262E-01  1.2706E+00  1.9525E-02
             1.0104E+01
 PARAMETER:  2.1944E-01 -8.2738E-01  4.5423E+00  7.2836E-01  1.0364E+00  5.7897E-01  2.2700E+00 -5.2995E-01  3.3948E-01 -3.8360E+00
             2.4129E+00
 GRADIENT:   1.2424E+01  5.4580E+00 -1.7987E-01  1.0494E+01  2.0265E+00  7.0698E+00  1.9607E+02  4.0706E-04  3.9513E+00  2.4190E-03
             3.5846E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -837.711139686430        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1181             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1265E+00  3.3220E-01  1.1323E+03  1.9299E+00  2.5681E+00  1.6081E+00  9.3850E+00  4.8622E-01  1.3074E+00  1.0000E-02
             1.0090E+01
 PARAMETER:  2.1916E-01 -1.0020E+00  7.1320E+00  7.5749E-01  1.0432E+00  5.7508E-01  2.3391E+00 -6.2109E-01  3.6801E-01 -4.7758E+00
             2.4115E+00
 GRADIENT:   1.2930E+01  6.8659E+00 -1.3094E-02  1.1815E+01  1.3603E+00  6.6647E+00  2.2338E+02  6.4109E-07  3.5395E+00  0.0000E+00
             3.3576E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -837.846516964209        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1344
 NPARAMETR:  1.1239E+00  2.6015E-01  7.4696E+03  1.9741E+00  2.5630E+00  1.6107E+00  9.0214E+00  4.8925E-01  1.2796E+00  1.0000E-02
             1.0064E+01
 PARAMETER:  2.1681E-01 -1.2465E+00  9.0186E+00  7.8011E-01  1.0412E+00  5.7669E-01  2.2996E+00 -6.1488E-01  3.4657E-01 -4.7769E+00
             2.4089E+00
 GRADIENT:  -1.0300E+00 -1.2258E+00 -2.0510E-03  3.6291E+00 -1.1351E+00  8.4947E-01  1.4192E+01  1.0479E-06 -4.3400E+00  0.0000E+00
            -5.3406E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -838.507447483449        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1531             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1274E+00  2.2962E-01  1.8313E+04  2.0028E+00  2.5923E+00  1.6095E+00  9.7814E+00  5.0193E-01  1.3659E+00  1.0000E-02
             1.0106E+01
 PARAMETER:  2.1989E-01 -1.3713E+00  9.9154E+00  7.9452E-01  1.0526E+00  5.7589E-01  2.3805E+00 -5.8929E-01  4.1182E-01 -4.7769E+00
             2.4132E+00
 GRADIENT:   1.1747E+01  6.0996E+00 -1.0153E-03  1.5147E+01  2.0122E+00  7.6586E+00  2.5046E+02  4.1479E-06  3.3273E+00  0.0000E+00
             3.3643E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -838.739652389325        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1712
 NPARAMETR:  1.1231E+00  1.8410E-01  3.9393E+05  2.0273E+00  2.5956E+00  1.6081E+00  1.0156E+01  5.0785E-01  1.3848E+00  1.0000E-02
             1.0098E+01
 PARAMETER:  2.1612E-01 -1.5923E+00  1.2984E+01  8.0671E-01  1.0538E+00  5.7504E-01  2.4181E+00 -5.7757E-01  4.2556E-01 -4.7769E+00
             2.4124E+00
 GRADIENT:  -4.0772E+01 -9.6576E-02 -5.3766E-05  4.0064E-01  1.1612E+00  8.5653E+00  3.3295E+01  1.6721E-05 -3.5808E+00  0.0000E+00
             1.6664E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -838.904246216559        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1903             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1236E+00  1.7371E-01  4.2365E+05  2.0562E+00  2.5949E+00  1.6117E+00  1.0354E+01  5.1291E-01  1.3994E+00  1.0000E-02
             1.0116E+01
 PARAMETER:  2.1658E-01 -1.6504E+00  1.3057E+01  8.2085E-01  1.0536E+00  5.7730E-01  2.4373E+00 -5.6765E-01  4.3602E-01 -4.7769E+00
             2.4141E+00
 GRADIENT:  -4.2523E+00  5.9465E+00 -4.3382E-05  7.1490E+00  2.3225E+00  2.6249E+01  2.8200E+02  6.5067E-05 -3.4720E+00  0.0000E+00
             4.6445E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -838.975079408167        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2080
 NPARAMETR:  1.1304E+00  1.6075E-01  7.7927E+06  2.0743E+00  2.5764E+00  1.6211E+00  1.0462E+01  5.1193E-01  1.3897E+00  1.0000E-02
             1.0093E+01
 PARAMETER:  2.2261E-01 -1.7279E+00  1.5969E+01  8.2963E-01  1.0464E+00  5.8309E-01  2.4477E+00 -5.6957E-01  4.2909E-01 -4.7769E+00
             2.4119E+00
 GRADIENT:  -3.7590E+01 -1.7785E+01 -3.8550E-05  1.4668E+01 -5.2305E-01 -5.4732E+00  3.2370E+01  2.8983E-04  2.6707E+01  0.0000E+00
             1.7584E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -839.049949531793        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     2251             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1278E+00  1.5208E-01  7.3240E+06  2.0793E+00  2.5822E+00  1.6118E+00  1.0787E+01  5.2368E-01  1.3794E+00  1.0000E-02
             1.0101E+01
 PARAMETER:  2.2025E-01 -1.7833E+00  1.5907E+01  8.3201E-01  1.0486E+00  5.7736E-01  2.4783E+00 -5.4688E-01  4.2167E-01 -4.7769E+00
             2.4127E+00
 GRADIENT:   1.1184E+01  9.6762E+00 -3.0628E-05 -5.3371E+01  2.9514E+00  7.2213E+01  3.1548E+02 -4.0135E-03 -3.1139E+01  0.0000E+00
             3.2021E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -839.105818463189        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     2408
 NPARAMETR:  1.1296E+00  1.4878E-01  1.1631E+08  2.0797E+00  2.5841E+00  1.6139E+00  1.0689E+01  5.2311E-01  1.3873E+00  1.0000E-02
             1.0103E+01
 PARAMETER:  2.2189E-01 -1.8053E+00  1.8672E+01  8.3222E-01  1.0494E+00  5.7864E-01  2.4692E+00 -5.4797E-01  4.2732E-01 -4.7769E+00
             2.4129E+00
 GRADIENT:  -4.6897E+01  1.0826E+01 -5.3103E-04 -1.0549E+02  6.9401E+00  3.4991E+00  7.2192E+01  9.9922E-03 -3.5264E+01  0.0000E+00
             3.5755E+01

0ITERATION NO.:   76    OBJECTIVE VALUE:  -839.105818463189        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     2437
 NPARAMETR:  1.1290E+00  1.5149E-01  1.1098E+08  2.0786E+00  2.5771E+00  1.6137E+00  1.0752E+01  5.2217E-01  1.3876E+00  1.0000E-02
             1.0103E+01
 PARAMETER:  2.2189E-01 -1.8053E+00  1.8672E+01  8.3222E-01  1.0494E+00  5.7864E-01  2.4692E+00 -5.4797E-01  4.2732E-01 -4.7769E+00
             2.4129E+00
 GRADIENT:   1.6541E+01 -1.6722E+00  2.9397E-01  3.0404E+00  6.9453E+00  4.9939E-01 -2.8131E+00  1.2539E+01 -3.2868E+00  0.0000E+00
             1.1044E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2437
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5772E-02  8.1187E-02  4.4736E-09 -8.8849E-02 -5.2295E-06
 SE:             2.7474E-02  1.9868E-02  5.7138E-09  1.7185E-02  4.1488E-05
 N:                     100         100         100         100         100

 P VAL.:         5.6591E-01  4.3860E-05  4.3366E-01  2.3445E-07  8.9970E-01

 ETASHRINKSD(%)  7.9599E+00  3.3439E+01  1.0000E+02  4.2427E+01  9.9861E+01
 ETASHRINKVR(%)  1.5286E+01  5.5696E+01  1.0000E+02  6.6854E+01  1.0000E+02
 EBVSHRINKSD(%)  1.0785E+01  4.1879E+01  1.0000E+02  2.9166E+01  9.9775E+01
 EBVSHRINKVR(%)  2.0406E+01  6.6220E+01  1.0000E+02  4.9826E+01  9.9999E+01
 RELATIVEINF(%)  7.5328E+01  2.4877E+01  0.0000E+00  3.4972E+01  4.6735E-04
 EPSSHRINKSD(%)  8.4191E+00
 EPSSHRINKVR(%)  1.6129E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -839.10581846318928     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       263.62042138241782     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    65.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    22.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -839.106       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.13E+00  1.49E-01  1.16E+08  2.08E+00  2.58E+00  1.61E+00  1.07E+01  5.23E-01  1.39E+00  1.00E-02  1.01E+01
 


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
+        2.95E+03
 
 TH 2
+        2.56E+02  8.43E+03
 
 TH 3
+       -1.83E-07  7.54E-07  4.49E-17
 
 TH 4
+       -2.92E+01  2.08E+02  3.23E-08  1.32E+02
 
 TH 5
+        2.25E+02 -1.03E+01 -1.15E-08  2.16E+01  2.76E+01
 
 TH 6
+        2.27E+02  6.68E+01 -3.59E-08 -7.01E+00 -5.53E+01  3.89E+01
 
 TH 7
+        1.68E+00 -7.47E+00  1.84E-10 -9.25E+00 -2.68E+00  7.00E+00  1.25E+00
 
 TH 8
+       -1.41E+03 -1.23E+03  3.54E-08  1.98E+02 -3.33E+02 -5.52E+02  9.72E+00  3.71E+03
 
 TH 9
+        2.20E+05 -9.21E+02  1.84E-08  4.81E+01 -3.88E+01  3.02E+02  1.61E+00  8.84E+02  4.60E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.74E+01 -6.23E+00 -7.01E-09 -1.12E+01 -1.56E+00  1.66E+01  9.24E-02  4.70E+03 -4.49E+00  0.00E+00  6.72E+00
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       88.235
Stop Time:
Thu Sep 30 09:22:01 CDT 2021
