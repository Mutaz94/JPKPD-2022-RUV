Sat Sep 25 14:01:03 CDT 2021
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
$DATA ../../../../data/spa/D/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1424.31365708769        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.7683E+01  1.6289E+01 -5.3426E+01  4.6304E+01  9.0943E+01 -1.1726E+02 -9.2576E+01 -8.1508E+00 -1.8448E+02 -1.9127E+01
            -5.3815E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1497.33271025085        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0281E+00  1.0425E+00  1.5613E+00  9.8473E-01  1.1030E+00  1.1801E+00  1.8625E+00  1.2030E+00  1.7628E+00  1.0315E+00
             1.0845E+00
 PARAMETER:  1.2772E-01  1.4160E-01  5.4552E-01  8.4609E-02  1.9801E-01  2.6561E-01  7.2194E-01  2.8484E-01  6.6693E-01  1.3105E-01
             1.8109E-01
 GRADIENT:   1.2268E+02  2.0601E+01  3.1621E+01  1.7705E+01 -3.3722E+01 -2.2175E+01  1.0546E+01 -1.3846E+01  4.4521E+01 -3.1686E+00
            -1.4438E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1503.39021627220        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9186E-01  8.9719E-01  2.2203E+00  1.1343E+00  1.2084E+00  1.2737E+00  1.8678E+00  2.1829E+00  1.5746E+00  1.1122E+00
             1.1499E+00
 PARAMETER:  9.1823E-02 -8.4875E-03  8.9762E-01  2.2600E-01  2.8932E-01  3.4194E-01  7.2475E-01  8.8064E-01  5.5397E-01  2.0634E-01
             2.3971E-01
 GRADIENT:   4.8189E+01  2.3258E+01  1.9581E+00  5.2171E+01  1.1850E+01  1.8280E+01  1.3520E+01  4.0567E+00  2.2458E+01  4.9820E-01
             1.5676E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1506.67765737703        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      274
 NPARAMETR:  9.8381E-01  8.9780E-01  1.7348E+00  1.0551E+00  1.1129E+00  1.2771E+00  1.7210E+00  1.7653E+00  1.5187E+00  1.0053E+00
             1.1012E+00
 PARAMETER:  8.3681E-02 -7.8062E-03  6.5089E-01  1.5366E-01  2.0692E-01  3.4459E-01  6.4288E-01  6.6832E-01  5.1785E-01  1.0529E-01
             1.9637E-01
 GRADIENT:   9.9120E-01  2.1249E+00  4.0395E+00 -5.4252E+00 -3.6543E+00  6.2803E-01  9.5511E-01 -9.3882E-01  6.2289E-01 -2.8022E+00
             1.0009E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1507.95770136519        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.7539E-01  5.1943E-01  2.0929E+00  1.3256E+00  1.0770E+00  1.2628E+00  1.8196E+00  1.8777E+00  1.3637E+00  1.0562E+00
             1.1007E+00
 PARAMETER:  7.5078E-02 -5.5502E-01  8.3856E-01  3.8187E-01  1.7415E-01  3.3333E-01  6.9860E-01  7.3005E-01  4.1024E-01  1.5465E-01
             1.9599E-01
 GRADIENT:  -6.3898E+00  4.7512E+00  2.0654E+00  9.9757E+00 -2.9116E+00 -2.4919E+00  7.0113E-01 -6.8302E-01  5.8946E-01  7.7108E-02
             2.0025E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1508.18248644441        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      638
 NPARAMETR:  9.7741E-01  4.7132E-01  2.1303E+00  1.3491E+00  1.0707E+00  1.2649E+00  1.8016E+00  1.9029E+00  1.3434E+00  1.0510E+00
             1.0981E+00
 PARAMETER:  7.7153E-02 -6.5221E-01  8.5626E-01  3.9941E-01  1.6829E-01  3.3501E-01  6.8866E-01  7.4338E-01  3.9519E-01  1.4972E-01
             1.9362E-01
 GRADIENT:  -2.8108E+00  2.9353E+00  2.0445E+00  4.0293E+00 -2.6064E+00 -1.6327E+00  6.1107E-01 -4.0368E-01 -3.0288E-01 -4.4658E-01
             1.2076E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1508.18789244084        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      830            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7745E-01  4.7126E-01  2.1308E+00  1.3489E+00  1.0707E+00  1.2702E+00  1.8000E+00  1.9025E+00  1.3447E+00  1.0546E+00
             1.0982E+00
 PARAMETER:  7.7193E-02 -6.5235E-01  8.5648E-01  3.9929E-01  1.6827E-01  3.3916E-01  6.8777E-01  7.4319E-01  3.9616E-01  1.5314E-01
             1.9364E-01
 GRADIENT:   3.7022E+01  1.0169E+01  2.8147E+00  5.5599E+01 -2.3430E+00  2.0151E+01  2.0728E+00 -2.3711E-02  9.6844E+00  6.2426E-02
             1.4964E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1508.22261051277        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:      981
 NPARAMETR:  9.7842E-01  4.7131E-01  2.1317E+00  1.3485E+00  1.0765E+00  1.2670E+00  1.3223E+00  1.9018E+00  1.3442E+00  1.0569E+00
             1.0912E+00
 PARAMETER:  7.8181E-02 -6.5224E-01  8.5690E-01  3.9901E-01  1.7370E-01  3.3669E-01  3.7937E-01  7.4282E-01  3.9576E-01  1.5537E-01
             1.8727E-01
 GRADIENT:  -1.5602E+00  4.2853E+00  4.4546E-01  6.9510E+00  2.3135E+00 -9.2199E-01 -2.9129E-01 -1.3315E+00 -7.0706E+00 -1.9589E+00
            -1.8926E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1508.32523615370        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1157
 NPARAMETR:  9.7936E-01  4.7131E-01  2.1317E+00  1.3485E+00  1.0773E+00  1.2697E+00  1.0654E+00  1.9018E+00  1.3864E+00  1.0770E+00
             1.0938E+00
 PARAMETER:  7.9141E-02 -6.5224E-01  8.5690E-01  3.9901E-01  1.7450E-01  3.3877E-01  1.6338E-01  7.4282E-01  4.2673E-01  1.7419E-01
             1.8964E-01
 GRADIENT:  -3.2291E-01  3.9114E+00  2.1370E-01  9.2434E+00 -7.8374E-03 -1.0419E-01  6.0288E-03 -1.0866E+00  1.4212E-01  1.1715E-02
            -3.6444E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1508.32550152229        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1333
 NPARAMETR:  9.7960E-01  4.7131E-01  2.1317E+00  1.3485E+00  1.0773E+00  1.2700E+00  1.0634E+00  1.9018E+00  1.3859E+00  1.0767E+00
             1.0948E+00
 PARAMETER:  7.9390E-02 -6.5224E-01  8.5690E-01  3.9901E-01  1.7449E-01  3.3902E-01  1.6146E-01  7.4282E-01  4.2638E-01  1.7391E-01
             1.9059E-01
 GRADIENT:   6.2293E-03  3.9344E+00  1.8819E-01  9.2429E+00 -5.3721E-03  1.8538E-03 -4.7678E-04 -1.0730E+00 -3.5105E-04  3.6030E-04
            -9.4465E-05

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1508.44790824426        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:     1482
 NPARAMETR:  9.7982E-01  4.4872E-01  2.1643E+00  1.3496E+00  1.0730E+00  1.2709E+00  1.0399E+00  1.9364E+00  1.3784E+00  1.0758E+00
             1.0931E+00
 PARAMETER:  7.9614E-02 -7.0136E-01  8.7208E-01  3.9982E-01  1.7045E-01  3.3969E-01  1.3914E-01  7.6085E-01  4.2093E-01  1.7302E-01
             1.8904E-01
 GRADIENT:   7.1451E-01  1.3225E+00  1.2567E+00 -2.1429E+00 -3.3970E+00  3.3425E-01  1.8172E-01 -1.3025E-01  8.5249E-01  6.7244E-01
            -1.3014E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1508.75420740394        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1659
 NPARAMETR:  9.7868E-01  3.4851E-01  2.2721E+00  1.4244E+00  1.0606E+00  1.2704E+00  8.4428E-01  2.0017E+00  1.3236E+00  1.0554E+00
             1.0888E+00
 PARAMETER:  7.8454E-02 -9.5408E-01  9.2069E-01  4.5374E-01  1.5882E-01  3.3932E-01 -6.9267E-02  7.9398E-01  3.8034E-01  1.5391E-01
             1.8509E-01
 GRADIENT:   5.4230E-01  2.3159E+00  2.7421E+00  4.5860E+00 -5.9744E+00  4.3496E-01  1.6077E-01  1.7494E-01 -5.5010E-01 -4.7997E-01
            -1.9259E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1508.75989311448        NO. OF FUNC. EVALS.: 204
 CUMULATIVE NO. OF FUNC. EVALS.:     1863
 NPARAMETR:  9.7871E-01  3.4738E-01  2.2735E+00  1.4148E+00  1.0604E+00  1.2679E+00  8.3650E-01  2.0008E+00  1.3267E+00  1.0553E+00
             1.0890E+00
 PARAMETER:  7.8484E-02 -9.5733E-01  9.2130E-01  4.4698E-01  1.5860E-01  3.3736E-01 -7.8531E-02  7.9354E-01  3.8269E-01  1.5383E-01
             1.8523E-01
 GRADIENT:   6.6667E-01  7.1433E-01  3.0386E+00 -4.2649E+00 -5.7313E+00 -3.2441E-01  2.0156E-01  1.2649E-01  3.9254E-01 -4.8538E-01
            -1.7462E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1508.87333228511        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2044
 NPARAMETR:  9.7871E-01  3.4738E-01  2.2734E+00  1.4171E+00  1.0604E+00  1.2691E+00  2.1866E-01  2.0008E+00  1.3381E+00  1.0553E+00
             1.0890E+00
 PARAMETER:  7.8484E-02 -9.5734E-01  9.2130E-01  4.4865E-01  1.5860E-01  3.3827E-01 -1.4202E+00  7.9354E-01  3.9128E-01  1.5383E-01
             1.8523E-01
 GRADIENT:   5.7203E-01  1.3382E+00  3.0269E+00 -2.4649E-01 -6.1078E+00  3.0148E-02  2.1430E-02 -9.7459E-02  7.0584E-01 -7.0441E-01
            -1.8511E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1508.88438118507        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2219
 NPARAMETR:  9.7871E-01  3.4738E-01  2.2734E+00  1.4171E+00  1.0604E+00  1.2689E+00  2.2185E-02  2.0008E+00  1.3368E+00  1.0553E+00
             1.0890E+00
 PARAMETER:  7.8484E-02 -9.5734E-01  9.2129E-01  4.4863E-01  1.5860E-01  3.3813E-01 -3.7084E+00  7.9354E-01  3.9031E-01  1.5383E-01
             1.8523E-01
 GRADIENT:   5.5911E-01  1.4589E+00  3.0196E+00 -9.9558E-02 -6.1287E+00 -2.4681E-02  2.2757E-04 -1.0973E-01  9.8996E-02 -7.2279E-01
            -1.8758E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1508.90457577805        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2402             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7837E-01  3.3851E-01  2.2370E+00  1.4172E+00  1.0662E+00  1.2690E+00  1.0000E-02  2.0025E+00  1.3366E+00  1.0604E+00
             1.0938E+00
 PARAMETER:  7.8137E-02 -9.8320E-01  9.0514E-01  4.4866E-01  1.6410E-01  3.3822E-01 -5.0601E+00  7.9442E-01  3.9009E-01  1.5861E-01
             1.8966E-01
 GRADIENT:   4.0287E+01  5.8451E+00 -1.6737E+00  5.7447E+01  6.7777E+00  2.0169E+01  0.0000E+00  1.7407E+00  1.2141E+01 -7.0815E-01
             5.6581E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1508.94174945814        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2587
 NPARAMETR:  9.7810E-01  3.3880E-01  2.2484E+00  1.4201E+00  1.0629E+00  1.2688E+00  1.0000E-02  1.9824E+00  1.3338E+00  1.0653E+00
             1.0938E+00
 PARAMETER:  7.7859E-02 -9.8234E-01  9.1020E-01  4.5070E-01  1.6103E-01  3.3807E-01 -5.0601E+00  7.8429E-01  3.8801E-01  1.6323E-01
             1.8968E-01
 GRADIENT:  -1.5532E-01  2.3607E-01  1.3754E-01 -2.6353E+00  6.1572E-01  4.3552E-03  0.0000E+00  2.7260E-01  8.1865E-01 -1.0517E-01
             2.4087E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1508.94369804646        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2778
 NPARAMETR:  9.7822E-01  3.3855E-01  2.2477E+00  1.4203E+00  1.0622E+00  1.2688E+00  1.0000E-02  1.9820E+00  1.3310E+00  1.0653E+00
             1.0931E+00
 PARAMETER:  7.7974E-02 -9.8309E-01  9.0992E-01  4.5086E-01  1.6033E-01  3.3806E-01 -5.0601E+00  7.8409E-01  3.8591E-01  1.6321E-01
             1.8904E-01
 GRADIENT:   1.4646E-02  4.1345E-01  2.8920E-01 -2.5550E+00  3.9048E-02  5.8169E-03  0.0000E+00  2.8119E-01  5.4395E-02 -4.1050E-02
            -1.0997E-02

0ITERATION NO.:   88    OBJECTIVE VALUE:  -1508.94542042053        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:     2887
 NPARAMETR:  9.7823E-01  3.3892E-01  2.2454E+00  1.4210E+00  1.0622E+00  1.2689E+00  1.0000E-02  1.9837E+00  1.3307E+00  1.0651E+00
             1.0932E+00
 PARAMETER:  7.7974E-02 -9.8297E-01  9.0981E-01  4.5181E-01  1.6032E-01  3.3807E-01 -5.0601E+00  7.8419E-01  3.8588E-01  1.6319E-01
             1.8905E-01
 GRADIENT:  -1.1591E-02 -5.8620E+04  3.1679E+04  1.2754E+05 -4.1947E-03 -1.3135E-02  0.0000E+00 -7.3582E+04  2.8478E-02  3.5312E+05
            -2.1029E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2887
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2352E-04 -1.5468E-04 -3.9467E-02 -3.4321E-03 -4.7761E-02
 SE:             2.9892E-02  4.9827E-05  1.8197E-02  2.9642E-02  2.0370E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9670E-01  1.9069E-03  3.0096E-02  9.0782E-01  1.9046E-02

 ETASHRINKSD(%)  1.0000E-10  9.9833E+01  3.9037E+01  6.9710E-01  3.1757E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.2835E+01  1.3893E+00  5.3429E+01
 EBVSHRINKSD(%)  3.1158E-01  9.9857E+01  4.4877E+01  1.1403E+00  2.7904E+01
 EBVSHRINKVR(%)  6.2218E-01  1.0000E+02  6.9614E+01  2.2677E+00  4.8022E+01
 RELATIVEINF(%)  9.8898E+01  2.7888E-05  1.1918E+01  1.5613E+01  1.6633E+01
 EPSSHRINKSD(%)  4.5746E+01
 EPSSHRINKVR(%)  7.0565E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1508.9454204205294     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -773.79459385679127     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1508.945       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  3.39E-01  2.25E+00  1.42E+00  1.06E+00  1.27E+00  1.00E-02  1.98E+00  1.33E+00  1.07E+00  1.09E+00
 


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
+        1.51E+09
 
 TH 2
+       -9.69E+02  1.30E+08
 
 TH 3
+        1.55E+02  4.57E+03  3.45E+06
 
 TH 4
+        4.93E+02 -6.74E+07 -2.10E+03  3.49E+07
 
 TH 5
+        8.65E+08 -2.55E+04  4.07E+03  1.32E+04  4.97E+08
 
 TH 6
+        3.57E-01 -6.70E+02  1.09E+02  3.47E+02 -8.64E-01  1.22E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -2.05E+02 -5.99E+03 -1.16E+04  2.77E+03 -5.45E+03 -1.43E+02  0.00E+00  5.98E+06
 
 TH 9
+       -5.98E-02  9.39E+03 -1.54E+03 -4.89E+03 -3.17E-01 -6.04E-01  0.00E+00  2.03E+03  1.08E+02
 
 TH10
+        1.83E+03 -2.49E+08 -3.62E+03 -2.47E+04  4.86E+04  1.28E+03  0.00E+00  4.77E+03 -1.81E+04  4.77E+08
 
 TH11
+       -6.23E+00 -6.56E+03  1.06E+03 -1.09E+08 -1.25E+01  1.55E+00  0.00E+00 -1.39E+03  4.11E+00  1.25E+04  1.66E+02
 
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
 #CPUT: Total CPU Time in Seconds,       45.390
Stop Time:
Sat Sep 25 14:01:53 CDT 2021
