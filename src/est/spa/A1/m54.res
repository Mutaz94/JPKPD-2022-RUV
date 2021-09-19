Sat Sep 18 09:18:21 CDT 2021
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
$DATA ../../../../data/spa/A1/dat54.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1246.29270258990        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.2781E+01  4.4278E+01  4.2669E+01  1.7736E+01  2.1358E+01 -3.7700E+01 -8.7865E+00 -2.7772E+01 -3.5540E+01 -1.4730E+01
            -7.8802E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1480.07369969074        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0715E+00  9.7295E-01  9.4213E-01  1.0071E+00  9.1059E-01  1.1338E+00  9.3479E-01  1.2861E+00  1.1406E+00  6.5563E-01
             1.8325E+00
 PARAMETER:  1.6911E-01  7.2576E-02  4.0393E-02  1.0703E-01  6.3397E-03  2.2562E-01  3.2565E-02  3.5163E-01  2.3153E-01 -3.2216E-01
             7.0569E-01
 GRADIENT:   3.6852E+01  1.4078E+01  2.3545E+01  1.5706E+00 -3.4421E+01  2.1033E+01  2.2916E-02 -5.8489E+00  1.0428E+01 -3.2960E+00
            -8.8668E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1483.67744874685        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0548E+00  1.0911E+00  7.9604E-01  9.3529E-01  9.0757E-01  1.0459E+00  1.2294E+00  1.1741E+00  1.0055E+00  3.8473E-01
             1.9214E+00
 PARAMETER:  1.5332E-01  1.8716E-01 -1.2811E-01  3.3098E-02  3.0125E-03  1.4491E-01  3.0649E-01  2.6051E-01  1.0544E-01 -8.5520E-01
             7.5305E-01
 GRADIENT:  -1.6573E+00  2.6755E+01  1.6831E+01  4.6040E+00  2.2636E+00 -1.0119E+01  1.1415E+01 -7.2273E+00 -2.1750E+00 -1.9110E+00
            -7.8102E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1492.03788305880        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0711E+00  1.1297E+00  4.3838E-01  8.6100E-01  7.1262E-01  1.0714E+00  1.0894E+00  7.6235E-01  9.4353E-01  2.3599E-01
             2.2450E+00
 PARAMETER:  1.6867E-01  2.2195E-01 -7.2468E-01 -4.9662E-02 -2.3880E-01  1.6899E-01  1.8567E-01 -1.7135E-01  4.1870E-02 -1.3440E+00
             9.0872E-01
 GRADIENT:   1.3117E+01  1.1136E+01  2.0305E+00  7.8172E+00 -2.0861E+00  2.0227E+00  6.0047E+00 -4.8101E-01  1.8116E+00  1.6814E+00
             2.4739E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1493.92021201273        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0484E+00  1.0483E+00  2.7179E-01  8.2630E-01  5.5084E-01  1.0556E+00  1.0526E+00  6.2923E-01  8.5912E-01  1.2353E-01
             2.0343E+00
 PARAMETER:  1.4726E-01  1.4713E-01 -1.2027E+00 -9.0797E-02 -4.9631E-01  1.5411E-01  1.5127E-01 -3.6326E-01 -5.1845E-02 -1.9912E+00
             8.1013E-01
 GRADIENT:  -1.3293E+01  9.0652E+00  3.7207E+00  2.3979E+00 -3.7057E+00 -4.6859E+00  2.7736E+00 -2.0406E+00 -4.0657E+00  7.8298E-01
            -9.5832E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1507.97776245532        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0587E+00  1.1001E+00  2.0403E-01  7.4742E-01  5.1880E-01  1.0798E+00  8.9156E-01  1.6057E+00  8.7851E-01  3.1781E-02
             2.0416E+00
 PARAMETER:  1.5703E-01  1.9541E-01 -1.4895E+00 -1.9113E-01 -5.5623E-01  1.7680E-01 -1.4788E-02  5.7357E-01 -2.9530E-02 -3.3489E+00
             8.1375E-01
 GRADIENT:   2.6401E+01 -8.3106E-01  8.2344E+00 -6.3968E+00 -6.0820E+01  1.0366E+01 -9.1228E+00  5.1030E+00 -2.2820E-01  6.4932E-02
             7.3179E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1533.72216650895        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.0324E+00  1.5762E+00  2.1492E-01  5.4921E-01  7.8389E-01  1.0364E+00  8.4249E-01  2.3017E+00  1.1531E+00  4.1248E-02
             1.6680E+00
 PARAMETER:  1.3184E-01  5.5499E-01 -1.4375E+00 -4.9928E-01 -1.4349E-01  1.3574E-01 -7.1390E-02  9.3366E-01  2.4245E-01 -3.0881E+00
             6.1163E-01
 GRADIENT:  -1.8582E+01  4.9272E+01  4.5571E+00  3.2889E+01 -7.9770E+00 -5.3430E+00  7.4898E+00 -5.1093E-01 -3.7188E+00  7.8157E-03
            -3.6965E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1537.31400859108        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      513
 NPARAMETR:  1.0395E+00  1.6765E+00  1.6206E-01  4.5454E-01  8.2830E-01  1.0456E+00  7.6966E-01  2.2014E+00  1.3363E+00  1.7774E-02
             1.6983E+00
 PARAMETER:  1.3876E-01  6.1672E-01 -1.7198E+00 -6.8848E-01 -8.8379E-02  1.4455E-01 -1.6181E-01  8.8911E-01  3.8988E-01 -3.9300E+00
             6.2965E-01
 GRADIENT:   7.7120E-01 -1.0066E+00 -5.2226E+00  9.4180E+00  6.7119E+00 -1.0592E+00 -1.1985E+00 -4.5373E-01 -2.4410E+00  1.9735E-04
            -4.6372E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1537.48736195091        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  1.0515E+00  1.6779E+00  1.6207E-01  4.5350E-01  8.2952E-01  1.0568E+00  7.7665E-01  2.2030E+00  1.3397E+00  1.7612E-02
             1.6997E+00
 PARAMETER:  1.5023E-01  6.1757E-01 -1.7197E+00 -6.9077E-01 -8.6907E-02  1.5524E-01 -1.5277E-01  8.8981E-01  3.9248E-01 -3.9392E+00
             6.3046E-01
 GRADIENT:   2.5414E+01 -2.0047E+00 -4.7465E+00  8.7867E+00  7.3415E+00  3.3675E+00  1.2124E+00 -8.2734E-01 -2.0053E+00 -5.2312E-05
            -3.9107E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1537.51035014842        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      804
 NPARAMETR:  1.0527E+00  1.6791E+00  1.6237E-01  4.5345E-01  8.2952E-01  1.0554E+00  7.6807E-01  2.2032E+00  1.3398E+00  1.7620E-02
             1.7000E+00
 PARAMETER:  1.5139E-01  6.1823E-01 -1.7179E+00 -6.9087E-01 -8.6912E-02  1.5392E-01 -1.6388E-01  8.8989E-01  3.9251E-01 -3.9387E+00
             6.3065E-01
 GRADIENT:   2.9047E+00 -2.6793E+01 -5.0922E+00  3.6970E+00  6.0272E+00  4.1179E-01 -1.9221E+00 -2.0344E+00 -2.8087E+00 -9.1893E-04
            -4.3914E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1537.67673864774        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      981
 NPARAMETR:  1.0515E+00  1.6917E+00  1.6578E-01  4.5288E-01  8.2947E-01  1.0544E+00  7.7340E-01  2.2053E+00  1.3402E+00  1.7717E-02
             1.7037E+00
 PARAMETER:  1.5025E-01  6.2571E-01 -1.6971E+00 -6.9212E-01 -8.6964E-02  1.5298E-01 -1.5696E-01  8.9087E-01  3.9281E-01 -3.9332E+00
             6.3283E-01
 GRADIENT:  -4.2692E+00 -7.5238E-01 -1.0883E-02  4.1867E+00 -6.8963E+00  9.5315E-01 -1.5110E-01 -3.5184E+00 -2.8758E+00  3.3775E-02
            -2.8282E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1537.81669550024        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0509E+00  1.6918E+00  1.6560E-01  4.5118E-01  8.3278E-01  1.0541E+00  7.7303E-01  2.2586E+00  1.3620E+00  1.7771E-02
             1.7093E+00
 PARAMETER:  1.4966E-01  6.2582E-01 -1.6982E+00 -6.9589E-01 -8.2985E-02  1.5270E-01 -1.5744E-01  9.1476E-01  4.0895E-01 -3.9302E+00
             6.3607E-01
 GRADIENT:   2.3238E+01  1.7690E+01 -7.2434E-03  7.4991E+00 -8.1970E+00  2.3855E+00  1.0305E+00 -3.9084E-01 -6.0378E-01 -3.5297E-04
             4.1272E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1537.83442194457        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1283
 NPARAMETR:  1.0515E+00  1.6967E+00  1.6569E-01  4.5128E-01  8.3281E-01  1.0542E+00  7.7065E-01  2.2580E+00  1.3649E+00  3.1007E-02
             1.7096E+00
 PARAMETER:  1.5023E-01  6.2868E-01 -1.6976E+00 -6.9568E-01 -8.2953E-02  1.5277E-01 -1.6053E-01  9.1448E-01  4.1110E-01 -3.3735E+00
             6.3627E-01
 GRADIENT:  -2.0824E-01 -4.9043E-01 -4.9001E-01  3.7396E+00 -1.3101E+01 -3.6722E-02  7.1442E-02 -1.3261E+00 -9.7580E-01 -2.7953E-03
            -1.0555E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1537.83471994481        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1356
 NPARAMETR:  1.0516E+00  1.6968E+00  1.6574E-01  4.5132E-01  8.3282E-01  1.0542E+00  7.7051E-01  2.2577E+00  1.3648E+00  3.4119E-02
             1.7098E+00
 PARAMETER:  1.5023E-01  6.2869E-01 -1.6976E+00 -6.9568E-01 -8.2953E-02  1.5277E-01 -1.6054E-01  9.1448E-01  4.1110E-01 -3.2774E+00
             6.3627E-01
 GRADIENT:  -1.4753E-01 -4.5643E-01 -6.6829E+01 -3.0937E+02 -1.0953E+03 -2.1850E-02  6.3922E-02  2.4659E+02  5.1629E+02  6.1113E+01
            -3.5291E+02
 NUMSIGDIG:         3.0         3.6         3.3         3.3         3.3         3.1         2.5         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1356
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.4513E-04 -1.3212E-02 -1.3660E-02  1.3626E-02 -2.1603E-03
 SE:             2.9634E-02  2.6604E-02  1.6973E-02  2.2095E-02  9.7787E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8263E-01  6.1947E-01  4.2093E-01  5.3741E-01  2.7159E-02

 ETASHRINKSD(%)  7.2377E-01  1.0873E+01  4.3137E+01  2.5980E+01  9.6724E+01
 ETASHRINKVR(%)  1.4423E+00  2.0564E+01  6.7666E+01  4.5210E+01  9.9893E+01
 EBVSHRINKSD(%)  1.0558E+00  1.1082E+01  4.3858E+01  2.5910E+01  9.6813E+01
 EBVSHRINKVR(%)  2.1004E+00  2.0937E+01  6.8480E+01  4.5107E+01  9.9898E+01
 RELATIVEINF(%)  9.6772E+01  1.3067E+01  1.1566E+01  9.8220E+00  1.6691E-02
 EPSSHRINKSD(%)  4.0431E+01
 EPSSHRINKVR(%)  6.4515E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1537.8347199448081     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -802.68389338106988     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.92
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1537.835       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.70E+00  1.66E-01  4.51E-01  8.33E-01  1.05E+00  7.71E-01  2.26E+00  1.36E+00  3.41E-02  1.71E+00
 


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
+        8.87E+02
 
 TH 2
+       -3.01E+05  4.60E+04
 
 TH 3
+        1.28E+02  6.69E+02  7.27E+05
 
 TH 4
+        3.93E+02  4.59E+02  5.06E+03  5.62E+05
 
 TH 5
+        2.04E+03  1.33E+02  1.52E+04  5.98E+04  7.90E+06
 
 TH 6
+       -1.47E+00 -2.94E+05  1.44E+02  3.72E+02  1.82E+03  1.73E+02
 
 TH 7
+        3.77E-04  6.28E+00  5.92E+01  3.99E+01  1.37E+02 -1.47E+00  2.14E+02
 
 TH 8
+       -2.37E+01 -4.47E+01 -7.07E+03 -1.01E+03 -2.63E+03 -2.13E+01 -1.13E+01  1.36E+04
 
 TH 9
+       -8.11E+02  5.05E+01 -9.92E+02 -2.36E+03 -1.15E+04 -7.11E+02  6.47E+01  1.52E+02  1.70E+05
 
 TH10
+       -2.93E+06  4.37E+05 -7.45E+02 -1.91E+03 -9.09E+03 -2.87E+06  3.51E+06  1.14E+02  3.54E+03  4.18E+06
 
 TH11
+        4.97E+01  6.72E+01  8.85E+03  2.59E+03  6.60E+03  5.53E+01  3.21E+01 -1.39E+03 -3.38E+02 -2.69E+02  4.86E+04
 
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
 #CPUT: Total CPU Time in Seconds,       26.114
Stop Time:
Sat Sep 18 09:18:49 CDT 2021
