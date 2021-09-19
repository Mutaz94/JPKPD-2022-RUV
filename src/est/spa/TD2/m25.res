Sat Sep 18 14:33:53 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1699.40111734719        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6959E+01 -7.6379E+01 -1.4972E+01 -1.1150E+02 -4.2694E+01  3.0445E+01 -3.9349E+00  1.7748E+01  3.0835E+00  2.8862E+01
             7.1941E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1709.46101516531        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0060E+00  1.1451E+00  1.2491E+00  1.0103E+00  1.1935E+00  8.8958E-01  1.0446E+00  8.4880E-01  9.5406E-01  8.4228E-01
             1.0293E+00
 PARAMETER:  1.0597E-01  2.3550E-01  3.2241E-01  1.1024E-01  2.7692E-01 -1.7010E-02  1.4365E-01 -6.3936E-02  5.2976E-02 -7.1638E-02
             1.2890E-01
 GRADIENT:   7.8102E+01  2.9863E+01  2.0922E+01  1.8438E+01  1.4282E+01 -1.4486E+01 -7.0323E-01 -5.4466E+00 -7.9697E+00 -2.1705E+01
             1.6298E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1712.57427819046        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9521E-01  1.0332E+00  9.1303E-01  1.0491E+00  1.0047E+00  9.0905E-01  1.1299E+00  2.7379E-01  9.0921E-01  9.1973E-01
             9.8682E-01
 PARAMETER:  9.5203E-02  1.3267E-01  9.0138E-03  1.4795E-01  1.0473E-01  4.6428E-03  2.2213E-01 -1.1954E+00  4.8256E-03  1.6322E-02
             8.6736E-02
 GRADIENT:   4.4939E+01 -1.2509E+01 -2.2932E+01  1.5467E+01  2.4473E+01 -5.4133E+00 -1.3759E+00  4.6704E-01 -3.8245E+00  1.0955E+01
            -3.6230E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1713.45805529510        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8209E-01  1.0866E+00  9.1274E-01  1.0124E+00  1.0138E+00  9.2241E-01  1.1013E+00  3.2962E-01  9.4957E-01  8.6622E-01
             9.8799E-01
 PARAMETER:  8.1929E-02  1.8308E-01  8.6990E-03  1.1237E-01  1.1368E-01  1.9236E-02  1.9649E-01 -1.0098E+00  4.8258E-02 -4.3619E-02
             8.7921E-02
 GRADIENT:   1.0374E+01 -2.2340E+00 -2.8149E+00  1.1899E+00  1.4907E+00  2.1281E-01 -1.1623E-01  3.8602E-01 -3.6832E-01  3.5218E+00
            -1.7448E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1713.45885070545        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.8068E-01  1.0907E+00  9.0567E-01  1.0093E+00  1.0108E+00  9.2260E-01  1.1006E+00  3.2285E-01  9.5093E-01  8.5465E-01
             9.9020E-01
 PARAMETER:  8.0487E-02  1.8685E-01  9.1634E-04  1.0924E-01  1.1070E-01  1.9441E-02  1.9588E-01 -1.0306E+00  4.9690E-02 -5.7061E-02
             9.0149E-02
 GRADIENT:   6.2559E+00 -1.7058E+00 -2.2344E+00  7.3306E-01  9.5215E-01  1.7689E-01 -1.2495E-01  3.9131E-01 -2.9913E-01  2.5313E+00
            -9.3400E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1713.45951232624        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.7980E-01  1.0940E+00  8.9841E-01  1.0066E+00  1.0076E+00  9.2273E-01  1.1004E+00  3.0982E-01  9.5194E-01  8.4561E-01
             9.9159E-01
 PARAMETER:  7.9596E-02  1.8988E-01 -7.1300E-03  1.0661E-01  1.0757E-01  1.9581E-02  1.9568E-01 -1.0718E+00  5.0745E-02 -6.7696E-02
             9.1553E-02
 GRADIENT:   3.6447E+00 -1.3142E+00 -1.8046E+00  4.4954E-01  6.2154E-01  1.4234E-01 -1.2304E-01  3.6989E-01 -2.4084E-01  1.8266E+00
            -4.3024E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1713.46001464618        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  9.7896E-01  1.0979E+00  8.8880E-01  1.0035E+00  1.0033E+00  9.2287E-01  1.1004E+00  2.8804E-01  9.5304E-01  8.3516E-01
             9.9295E-01
 PARAMETER:  7.8732E-02  1.9339E-01 -1.7888E-02  1.0347E-01  1.0334E-01  1.9729E-02  1.9570E-01 -1.1447E+00  5.1905E-02 -8.0130E-02
             9.2925E-02
 GRADIENT:   1.0531E+00 -8.7536E-01 -1.3179E+00  1.7007E-01  2.9643E-01  9.8863E-02 -1.1131E-01  3.2598E-01 -1.6852E-01  1.0649E+00
             5.5425E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1713.46115255051        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  9.7833E-01  1.1013E+00  8.7929E-01  1.0006E+00  9.9907E-01  9.2298E-01  1.1007E+00  2.6224E-01  9.5398E-01  8.2582E-01
             9.9397E-01
 PARAMETER:  7.8089E-02  1.9651E-01 -2.8640E-02  1.0058E-01  9.9072E-02  1.9849E-02  1.9591E-01 -1.2385E+00  5.2890E-02 -9.1373E-02
             9.3954E-02
 GRADIENT:  -9.3000E-01 -4.9999E-01 -8.9627E-01 -4.2970E-02  6.1351E-02  5.7597E-02 -9.4463E-02  2.7307E-01 -1.0056E-01  4.2701E-01
             4.1558E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1713.46942096233        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  9.7742E-01  1.1101E+00  8.5932E-01  9.9340E-01  9.9111E-01  9.2316E-01  1.0996E+00  1.9824E-01  9.5668E-01  8.0859E-01
             9.9549E-01
 PARAMETER:  7.7161E-02  2.0446E-01 -5.1617E-02  9.3375E-02  9.1073E-02  2.0051E-02  1.9490E-01 -1.5183E+00  5.5712E-02 -1.1246E-01
             9.5478E-02
 GRADIENT:  -3.9564E+00  1.6899E-01 -2.3802E-01 -2.8710E-01 -2.9232E-01 -2.5375E-02 -5.0101E-02  1.6473E-01  2.4493E-02 -6.3714E-01
             9.5884E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1713.57609337741        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  9.7865E-01  1.1104E+00  8.4479E-01  9.9214E-01  9.8329E-01  9.2341E-01  1.1021E+00  3.1541E-02  9.5535E-01  8.1107E-01
             9.9290E-01
 PARAMETER:  7.8423E-02  2.0474E-01 -6.8662E-02  9.2113E-02  8.3150E-02  2.0320E-02  1.9719E-01 -3.3565E+00  5.4320E-02 -1.0940E-01
             9.2875E-02
 GRADIENT:  -7.6626E-01  1.0419E-01 -2.7522E-01 -5.3125E-02 -7.1421E-01  1.3250E-01  1.2730E-01  4.0480E-03  8.6165E-03  4.7556E-01
             2.3453E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1714.35012612445        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  1.0002E+00  1.2706E+00  8.2057E-01  8.9921E-01  1.0584E+00  9.3238E-01  9.9844E-01  1.0000E-02  1.0393E+00  8.4657E-01
             9.9469E-01
 PARAMETER:  1.0022E-01  3.3951E-01 -9.7755E-02 -6.2346E-03  1.5673E-01  2.9982E-02  9.8434E-02 -6.4486E+00  1.3856E-01 -6.6568E-02
             9.4675E-02
 GRADIENT:   1.3487E+01 -1.5127E+00 -8.3856E-01  2.4316E+00  1.4379E+00  7.1229E-01  2.0274E-01  0.0000E+00  1.0430E-01 -2.1468E-01
             3.1733E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1714.67843816491        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      994
 NPARAMETR:  9.9586E-01  1.5425E+00  7.5368E-01  7.2782E-01  1.1902E+00  9.3150E-01  8.5910E-01  1.0000E-02  1.2239E+00  9.2461E-01
             9.9577E-01
 PARAMETER:  9.5851E-02  5.3342E-01 -1.8279E-01 -2.1769E-01  2.7414E-01  2.9042E-02 -5.1870E-02 -8.7699E+00  3.0203E-01  2.1622E-02
             9.5763E-02
 GRADIENT:   9.6879E-01  5.1193E+00  1.6025E+00  2.1859E+00 -1.7789E+00  2.3001E-01 -3.1172E-01  0.0000E+00 -1.3626E-01 -2.9022E-01
            -4.3026E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1714.68902852013        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  9.9550E-01  1.5546E+00  7.4422E-01  7.1703E-01  1.1956E+00  9.3097E-01  8.5525E-01  1.0000E-02  1.2340E+00  9.2689E-01
             9.9623E-01
 PARAMETER:  9.5490E-02  5.4120E-01 -1.9542E-01 -2.3264E-01  2.7861E-01  2.8474E-02 -5.6360E-02 -9.0936E+00  3.1028E-01  2.4078E-02
             9.6227E-02
 GRADIENT:  -5.2094E-02  1.7235E-02 -5.8117E-03  2.3986E-02  6.6373E-03 -2.3411E-03 -1.0324E-02  0.0000E+00 -3.4423E-03 -3.6469E-04
             4.9626E-03

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1714.68903025468        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1271
 NPARAMETR:  9.9552E-01  1.5545E+00  7.4423E-01  7.1704E-01  1.1955E+00  9.3098E-01  8.5526E-01  1.0000E-02  1.2340E+00  9.2686E-01
             9.9621E-01
 PARAMETER:  9.5517E-02  5.4117E-01 -1.9541E-01 -2.3267E-01  2.7860E-01  2.8483E-02 -5.6246E-02 -9.0936E+00  3.1030E-01  2.4070E-02
             9.6213E-02
 GRADIENT:   2.7985E-03 -7.7891E-03 -3.0252E-04 -2.3388E-02 -2.3507E-03 -5.0876E-04  4.1663E-03  0.0000E+00  2.8040E-03  1.3609E-03
             1.5567E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1271
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.5163E-05 -2.2626E-02 -3.2569E-04  1.6147E-02 -3.3947E-02
 SE:             2.9816E-02  2.2664E-02  1.2803E-04  2.3249E-02  2.1766E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9772E-01  3.1814E-01  1.0964E-02  4.8735E-01  1.1886E-01

 ETASHRINKSD(%)  1.1351E-01  2.4072E+01  9.9571E+01  2.2112E+01  2.7080E+01
 ETASHRINKVR(%)  2.2690E-01  4.2350E+01  9.9998E+01  3.9335E+01  4.6826E+01
 EBVSHRINKSD(%)  4.8683E-01  2.3179E+01  9.9606E+01  2.3244E+01  2.6014E+01
 EBVSHRINKVR(%)  9.7128E-01  4.0986E+01  9.9998E+01  4.1085E+01  4.5261E+01
 RELATIVEINF(%)  9.8874E+01  3.3716E+00  2.1555E-04  3.7552E+00  9.2677E+00
 EPSSHRINKSD(%)  4.2797E+01
 EPSSHRINKVR(%)  6.7278E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.6890302546781     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -979.53820369093989     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.13
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.689       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.55E+00  7.44E-01  7.17E-01  1.20E+00  9.31E-01  8.55E-01  1.00E-02  1.23E+00  9.27E-01  9.96E-01
 


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
+        1.28E+03
 
 TH 2
+       -6.27E+00  3.89E+02
 
 TH 3
+        1.32E+01  1.58E+02  3.53E+02
 
 TH 4
+       -1.39E+01  3.22E+02 -2.19E+02  8.45E+02
 
 TH 5
+       -4.09E+00 -1.94E+02 -3.01E+02  2.17E+02  4.97E+02
 
 TH 6
+        1.51E+00 -7.20E-01  3.60E+00 -4.83E+00 -3.11E-01  2.24E+02
 
 TH 7
+        3.63E+00  1.69E+01 -7.34E+00 -1.29E+01 -8.35E+00  4.13E+00  1.07E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.93E-01 -2.15E+01 -2.98E+01  4.71E+01 -8.25E-01 -9.25E-01  1.64E+01  0.00E+00  5.84E+01
 
 TH10
+       -2.67E+00 -1.07E+01 -3.30E+01 -9.22E+00 -6.24E+01 -2.62E+00  1.19E+01  0.00E+00  6.17E+00  8.31E+01
 
 TH11
+       -8.75E+00 -1.97E+01 -3.02E+01  3.54E+00  5.33E-01 -6.72E-02  1.10E+01  0.00E+00  5.45E+00  2.14E+01  2.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.014
Stop Time:
Sat Sep 18 14:34:12 CDT 2021
