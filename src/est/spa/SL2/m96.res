Sat Sep 25 11:25:57 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat96.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1631.45853687810        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2414E+02 -9.3328E+00 -3.6478E+01  3.6334E+01  3.4942E+01 -1.0040E+01  9.2620E+00  9.8062E+00  2.4640E+01  1.3627E+01
            -5.4798E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1640.47906554770        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.5986E-01  1.0333E+00  1.0834E+00  9.7611E-01  1.0212E+00  1.0182E+00  9.4895E-01  9.4867E-01  8.9816E-01  9.2535E-01
             1.1391E+00
 PARAMETER:  5.9030E-02  1.3272E-01  1.8007E-01  7.5825E-02  1.2093E-01  1.1806E-01  4.7597E-02  4.7301E-02 -7.4123E-03  2.2416E-02
             2.3025E-01
 GRADIENT:   2.3490E+01  8.5175E+00 -4.0807E+00  1.3299E+01  6.2076E+00  1.2585E+00  1.9316E+00  3.3978E+00  8.6674E-01 -4.7204E-01
             2.1229E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1640.93019837447        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.5819E-01  9.4513E-01  1.1077E+00  1.0312E+00  1.0004E+00  1.0352E+00  9.4720E-01  6.3041E-01  9.1340E-01  9.8757E-01
             1.1395E+00
 PARAMETER:  5.7287E-02  4.3571E-02  2.0231E-01  1.3074E-01  1.0041E-01  1.3464E-01  4.5760E-02 -3.6138E-01  9.4142E-03  8.7487E-02
             2.3060E-01
 GRADIENT:   2.1364E+01  5.2518E+00  5.2415E-03  1.5873E+01  4.3019E+00  8.4588E+00  1.8175E+00  1.0181E-01  7.5666E+00  3.3694E+00
             1.4268E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1641.33910265142        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.5253E-01  9.3286E-01  9.6650E-01  1.0273E+00  9.2967E-01  1.0203E+00  1.0417E+00  4.7524E-01  8.6133E-01  9.0492E-01
             1.1327E+00
 PARAMETER:  5.1362E-02  3.0495E-02  6.5923E-02  1.2689E-01  2.7078E-02  1.2007E-01  1.4090E-01 -6.4394E-01 -4.9282E-02  8.5761E-05
             2.2462E-01
 GRADIENT:   6.8713E+00  2.0919E+00 -2.3143E+00  8.0017E+00  3.9273E+00  1.7259E+00  1.8941E+00  4.8779E-01  2.3469E+00  1.3411E+00
             1.4099E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1641.34252960513        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.5148E-01  9.4193E-01  9.3849E-01  1.0184E+00  9.1941E-01  1.0187E+00  1.0338E+00  4.2739E-01  8.6026E-01  8.9321E-01
             1.1327E+00
 PARAMETER:  5.0267E-02  4.0174E-02  3.6513E-02  1.1820E-01  1.5974E-02  1.1853E-01  1.3323E-01 -7.5006E-01 -5.0524E-02 -1.2937E-02
             2.2464E-01
 GRADIENT:   3.9081E+00  7.6779E-01 -1.8303E+00  4.3761E+00  2.6343E+00  9.1434E-01  1.1969E+00  4.2160E-01  1.4225E+00  1.0013E+00
             1.5265E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1641.34321006734        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.5103E-01  9.4795E-01  9.2171E-01  1.0129E+00  9.1346E-01  1.0181E+00  1.0292E+00  3.9067E-01  8.6010E-01  8.8667E-01
             1.1328E+00
 PARAMETER:  4.9793E-02  4.6542E-02  1.8480E-02  1.1285E-01  9.4881E-03  1.1790E-01  1.2882E-01 -8.3989E-01 -5.0712E-02 -2.0280E-02
             2.2471E-01
 GRADIENT:   2.5062E+00  2.9297E-01 -1.4455E+00  2.7295E+00  1.8784E+00  5.5373E-01  8.2643E-01  3.5226E-01  9.5313E-01  7.6879E-01
             1.3832E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1641.34395584672        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.5069E-01  9.5358E-01  9.0617E-01  1.0079E+00  9.0792E-01  1.0176E+00  1.0252E+00  3.4967E-01  8.6008E-01  8.8076E-01
             1.1329E+00
 PARAMETER:  4.9435E-02  5.2464E-02  1.4752E-03  1.0792E-01  3.3977E-03  1.1744E-01  1.2489E-01 -9.5078E-01 -5.0728E-02 -2.6972E-02
             2.2480E-01
 GRADIENT:   1.3771E+00 -1.2112E-02 -9.9067E-01  1.3810E+00  1.1357E+00  2.7245E-01  4.9400E-01  2.7447E-01  5.4902E-01  5.2017E-01
             1.0086E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1641.34439667015        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      517
 NPARAMETR:  9.5051E-01  9.5715E-01  8.9544E-01  1.0048E+00  9.0400E-01  1.0174E+00  1.0232E+00  3.1776E-01  8.6003E-01  8.7690E-01
             1.1330E+00
 PARAMETER:  4.9244E-02  5.6207E-02 -1.0437E-02  1.0474E-01 -9.2775E-04  1.1720E-01  1.2293E-01 -1.0465E+00 -5.0789E-02 -3.1357E-02
             2.2486E-01
 GRADIENT:   7.1685E-01 -1.9575E-01 -7.8379E-01  6.4032E-01  7.6355E-01  1.1229E-01  3.1115E-01  2.2408E-01  3.2215E-01  3.9622E-01
             9.6742E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1641.34499529380        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      587
 NPARAMETR:  9.5032E-01  9.6138E-01  8.8229E-01  1.0010E+00  8.9911E-01  1.0171E+00  1.0211E+00  2.7214E-01  8.5996E-01  8.7226E-01
             1.1331E+00
 PARAMETER:  4.9049E-02  6.0614E-02 -2.5233E-02  1.0096E-01 -6.3519E-03  1.1696E-01  1.2085E-01 -1.2014E+00 -5.0870E-02 -3.6663E-02
             2.2494E-01
 GRADIENT:  -1.2040E-02 -3.8843E-01 -5.3444E-01 -1.7803E-01  3.4135E-01 -6.1586E-02  1.0396E-01  1.6111E-01  6.8343E-02  2.4828E-01
             8.9011E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1641.34675513801        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  9.5024E-01  9.6881E-01  8.6743E-01  9.9511E-01  8.9474E-01  1.0170E+00  1.0171E+00  2.1192E-01  8.6127E-01  8.6777E-01
             1.1332E+00
 PARAMETER:  4.8954E-02  6.8309E-02 -4.2226E-02  9.5103E-02 -1.1224E-02  1.1684E-01  1.1697E-01 -1.4516E+00 -4.9353E-02 -4.1824E-02
             2.2504E-01
 GRADIENT:  -5.7196E-01 -5.1777E-01 -2.9090E-01 -7.8948E-01 -1.0454E-01 -1.8663E-01 -5.5138E-02  9.6860E-02 -1.2429E-01  1.4019E-01
             1.0050E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1641.34822814199        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  9.5025E-01  9.8165E-01  8.5568E-01  9.8637E-01  8.9451E-01  1.0170E+00  1.0085E+00  1.5409E-01  8.6545E-01  8.6566E-01
             1.1333E+00
 PARAMETER:  4.8969E-02  8.1478E-02 -5.5861E-02  8.6279E-02 -1.1483E-02  1.1685E-01  1.0846E-01 -1.7702E+00 -4.4500E-02 -4.4259E-02
             2.2510E-01
 GRADIENT:  -8.4361E-01 -4.9486E-01 -1.1206E-01 -1.0586E+00 -2.8974E-01 -2.4105E-01 -1.5993E-01  5.0664E-02 -2.4092E-01  2.9373E-02
             7.0451E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1641.35503098131        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  9.5041E-01  9.9821E-01  8.4378E-01  9.7543E-01  8.9640E-01  1.0171E+00  9.9844E-01  5.6797E-02  8.7182E-01  8.6544E-01
             1.1333E+00
 PARAMETER:  4.9137E-02  9.8206E-02 -6.9859E-02  7.5121E-02 -9.3705E-03  1.1698E-01  9.8440E-02 -2.7683E+00 -3.7167E-02 -4.4519E-02
             2.2514E-01
 GRADIENT:  -7.9300E-01 -6.9904E-01 -4.7107E-01 -1.0546E+00  3.3882E-01 -2.3329E-01 -9.7527E-02  7.4276E-03 -1.7718E-01  1.0530E-01
             1.1816E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1641.38447073943        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      873
 NPARAMETR:  9.5063E-01  9.8062E-01  8.5186E-01  9.8719E-01  8.9243E-01  1.0175E+00  1.0116E+00  1.1896E-02  8.6562E-01  8.6629E-01
             1.1334E+00
 PARAMETER:  4.9373E-02  8.0425E-02 -6.0328E-02  8.7108E-02 -1.3812E-02  1.1737E-01  1.1153E-01 -4.3316E+00 -4.4306E-02 -4.3532E-02
             2.2521E-01
 GRADIENT:  -5.0626E-02 -2.9671E-03 -3.4262E-02 -2.8509E-02  1.1738E-02 -1.6036E-02  4.6757E-03  3.0183E-04 -1.2775E-02  1.6169E-02
             1.1861E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1641.67060270119        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1036
 NPARAMETR:  9.6521E-01  9.1929E-01  8.7378E-01  1.0284E+00  8.7585E-01  1.0271E+00  1.0620E+00  1.0000E-02  8.4159E-01  8.6576E-01
             1.1350E+00
 PARAMETER:  6.4595E-02  1.5842E-02 -3.4921E-02  1.2798E-01 -3.2557E-02  1.2678E-01  1.6011E-01 -1.0500E+01 -7.2457E-02 -4.4143E-02
             2.2667E-01
 GRADIENT:   5.6467E-01 -1.0971E-01  6.8708E-01 -1.3079E+00 -1.5080E+00 -8.9131E-02 -1.6707E-01  0.0000E+00 -4.1983E-01 -1.1701E-01
             3.5377E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1641.67411789495        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1211
 NPARAMETR:  9.6501E-01  9.3612E-01  8.7515E-01  1.0189E+00  8.8477E-01  1.0274E+00  1.0469E+00  1.0000E-02  8.5035E-01  8.7165E-01
             1.1350E+00
 PARAMETER:  6.4382E-02  3.3987E-02 -3.3363E-02  1.1872E-01 -2.2433E-02  1.2705E-01  1.4582E-01 -1.0465E+01 -6.2102E-02 -3.7364E-02
             2.2663E-01
 GRADIENT:   5.1371E-03  9.4648E-03  1.0228E-03  1.3931E-02 -3.6526E-03 -4.8618E-04 -2.7168E-03  0.0000E+00 -4.0532E-03 -1.8693E-03
            -1.5452E-03

0ITERATION NO.:   71    OBJECTIVE VALUE:  -1641.67411789495        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.6501E-01  9.3612E-01  8.7515E-01  1.0189E+00  8.8477E-01  1.0274E+00  1.0469E+00  1.0000E-02  8.5035E-01  8.7165E-01
             1.1350E+00
 PARAMETER:  6.4382E-02  3.3987E-02 -3.3363E-02  1.1872E-01 -2.2433E-02  1.2705E-01  1.4582E-01 -1.0465E+01 -6.2102E-02 -3.7364E-02
             2.2663E-01
 GRADIENT:   5.1371E-03  9.4648E-03  1.0228E-03  1.3931E-02 -3.6526E-03 -4.8618E-04 -2.7168E-03  0.0000E+00 -4.0532E-03 -1.8693E-03
            -1.5452E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1233
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.3365E-04 -5.0217E-03 -3.0392E-04 -2.9860E-03 -1.8562E-02
 SE:             2.9796E-02  1.9237E-02  1.6514E-04  2.4277E-02  2.3276E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9374E-01  7.9406E-01  6.5722E-02  9.0211E-01  4.2516E-01

 ETASHRINKSD(%)  1.7901E-01  3.5553E+01  9.9447E+01  1.8668E+01  2.2024E+01
 ETASHRINKVR(%)  3.5770E-01  5.8466E+01  9.9997E+01  3.3852E+01  3.9197E+01
 EBVSHRINKSD(%)  5.0676E-01  3.5567E+01  9.9474E+01  1.8822E+01  2.0732E+01
 EBVSHRINKVR(%)  1.0109E+00  5.8484E+01  9.9997E+01  3.4101E+01  3.7166E+01
 RELATIVEINF(%)  9.8441E+01  1.2810E+00  2.5467E-04  2.5242E+00  4.9022E+00
 EPSSHRINKSD(%)  4.1669E+01
 EPSSHRINKVR(%)  6.5975E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1641.6741178949492     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -906.52329133121100     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1641.674       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.65E-01  9.36E-01  8.75E-01  1.02E+00  8.85E-01  1.03E+00  1.05E+00  1.00E-02  8.50E-01  8.72E-01  1.13E+00
 


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
+        1.12E+03
 
 TH 2
+       -1.21E+01  5.00E+02
 
 TH 3
+        1.43E+01  1.98E+02  4.29E+02
 
 TH 4
+       -1.48E+01  4.79E+02 -1.90E+02  9.81E+02
 
 TH 5
+       -5.80E+00 -3.93E+02 -6.23E+02  2.19E+02  1.17E+03
 
 TH 6
+       -5.78E-01 -3.13E+00  3.11E+00 -3.54E+00 -2.79E+00  1.87E+02
 
 TH 7
+        6.49E-01  2.34E+01  7.76E+00 -7.05E+00 -1.33E+01  3.67E-01  3.82E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.34E+00 -2.33E+01 -1.83E+01  2.42E+01  1.04E+01  5.93E-01  2.73E+01  0.00E+00  1.38E+02
 
 TH10
+       -1.85E+00  1.75E-01 -3.80E+01 -1.64E+01 -5.83E+01 -1.39E+00  1.26E+01  0.00E+00  4.00E+00  1.09E+02
 
 TH11
+       -7.05E+00 -1.61E+01 -3.45E+01 -6.87E+00  9.81E+00  2.20E+00  5.26E+00  0.00E+00  1.09E+01  2.67E+01  1.76E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.149
Stop Time:
Sat Sep 25 11:26:14 CDT 2021
