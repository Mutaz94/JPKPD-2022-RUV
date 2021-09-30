Wed Sep 29 13:28:45 CDT 2021
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
$DATA ../../../../data/spa/A3/dat36.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   182.038453910194        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3170E+02  6.0426E+01  9.7528E+01 -9.5673E+00  2.7871E+02  2.5352E+01 -6.6602E+01 -6.1315E+01 -1.6744E+02 -1.9799E+02
            -3.1250E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1198.06379173943        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0505E+00  9.6137E-01  8.8229E-01  1.1337E+00  8.7322E-01  8.8045E-01  1.0100E+00  1.0279E+00  1.1460E+00  1.0054E+00
             4.3155E+00
 PARAMETER:  1.4931E-01  6.0601E-02 -2.5238E-02  2.2546E-01 -3.5565E-02 -2.7325E-02  1.0997E-01  1.2752E-01  2.3629E-01  1.0542E-01
             1.5622E+00
 GRADIENT:   9.1392E+01 -1.1284E+01 -1.7065E+01  7.7695E+00  1.8106E+01 -1.2182E+01  6.9015E+00  5.2444E+00  1.2593E+01  1.9558E+01
             3.5098E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1208.63065457757        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0541E+00  7.3559E-01  4.1688E-01  1.2442E+00  4.5289E-01  9.7142E-01  7.6913E-01  6.9725E-01  1.2328E+00  3.9910E-01
             4.1573E+00
 PARAMETER:  1.5270E-01 -2.0709E-01 -7.7495E-01  3.1853E-01 -6.9211E-01  7.1004E-02 -1.6250E-01 -2.6062E-01  3.0930E-01 -8.1854E-01
             1.5249E+00
 GRADIENT:   6.5668E+01  6.2424E+01  3.4576E+01  8.7924E+01 -4.9573E+01  4.6927E+00 -2.2493E+00  7.3014E-01  6.9509E+00  1.5968E+00
             4.5928E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1211.71720344948        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0225E+00  5.9882E-01  3.0544E-01  1.2913E+00  3.4557E-01  1.0193E+00  1.0941E+00  6.1149E-01  1.1468E+00  2.8945E-01
             3.8031E+00
 PARAMETER:  1.2223E-01 -4.1280E-01 -1.0860E+00  3.5564E-01 -9.6257E-01  1.1914E-01  1.8991E-01 -3.9185E-01  2.3695E-01 -1.1398E+00
             1.4358E+00
 GRADIENT:   1.5441E+00  7.9121E+01  4.1970E+01  1.5807E+02 -6.8794E+01  1.3290E+01 -7.1620E+00 -3.9361E+00 -1.8623E+01 -1.1989E+00
             1.7624E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1215.90229742064        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  9.9080E-01  6.1101E-01  3.0396E-01  1.2311E+00  3.5604E-01  1.0133E+00  1.5824E+00  4.3245E-01  1.0944E+00  2.2667E-01
             3.6569E+00
 PARAMETER:  9.0758E-02 -3.9265E-01 -1.0908E+00  3.0790E-01 -9.3272E-01  1.1320E-01  5.5893E-01 -7.3828E-01  1.9024E-01 -1.3843E+00
             1.3966E+00
 GRADIENT:  -7.5839E+01  4.8590E+01  9.8855E+00  1.0143E+02 -3.0637E+01  8.9741E+00  4.2238E+00 -2.4286E+00 -1.9730E+01 -4.2328E-01
             3.4771E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1226.98045455147        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      543
 NPARAMETR:  1.0130E+00  5.7588E-01  1.8687E-01  1.0297E+00  2.8908E-01  9.9731E-01  1.3484E+00  3.1388E-02  1.3412E+00  9.4251E-02
             3.4324E+00
 PARAMETER:  1.1291E-01 -4.5185E-01 -1.5773E+00  1.2922E-01 -1.1411E+00  9.7308E-02  3.9888E-01 -3.3613E+00  3.9353E-01 -2.2618E+00
             1.3333E+00
 GRADIENT:   5.4100E+00 -5.1845E+00 -4.2643E+00  6.4121E+00  5.3938E+00 -3.5539E-02 -1.9178E+00 -2.2806E-02  5.6564E-01 -1.9392E-01
             5.3886E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1227.21109617340        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      680
 NPARAMETR:  1.0107E+00  5.8625E-01  1.8782E-01  1.0172E+00  2.9274E-01  9.9251E-01  1.3463E+00  3.8258E-02  1.3403E+00  2.5271E-01
             3.3922E+00
 PARAMETER:  1.1066E-01 -4.3401E-01 -1.5723E+00  1.1704E-01 -1.1285E+00  9.2478E-02  3.9737E-01 -3.1634E+00  3.9289E-01 -1.2755E+00
             1.3215E+00
 GRADIENT:   3.3812E+01 -5.9107E-01  2.1878E+00  3.0448E+00  7.1748E+01  1.1587E+00  4.1214E+00 -2.2357E-02  5.1523E+00  5.2056E-01
             1.1774E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1227.48838717277        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.0097E+00  5.8113E-01  1.8797E-01  1.0260E+00  2.8723E-01  9.9963E-01  1.2663E+00  4.8558E-02  1.3452E+00  3.2934E-01
             3.3824E+00
 PARAMETER:  1.0967E-01 -4.4279E-01 -1.5715E+00  1.2565E-01 -1.1475E+00  9.9626E-02  3.3610E-01 -2.9250E+00  3.9654E-01 -1.0107E+00
             1.3186E+00
 GRADIENT:  -1.9587E+00  5.8385E+00  1.5311E+00  3.2590E-01 -3.3220E+00  7.7351E-01 -4.7286E-01 -4.0868E-02  2.2179E-02 -1.9038E-01
            -1.9152E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1230.55561809806        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  1.0091E+00  4.6456E-01  1.8016E-01  1.0463E+00  2.5128E-01  9.8145E-01  1.0210E+00  1.0693E+00  1.3482E+00  6.4350E-01
             3.3412E+00
 PARAMETER:  1.0902E-01 -6.6667E-01 -1.6139E+00  1.4530E-01 -1.2812E+00  8.1279E-02  1.2077E-01  1.6700E-01  3.9876E-01 -3.4084E-01
             1.3063E+00
 GRADIENT:  -1.3617E+01 -6.9759E+00 -4.3421E+00 -4.7449E+00  2.7455E+01 -6.8017E-01  4.3386E+00  4.5486E-01 -2.7843E+00  1.2588E+01
             2.8031E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1233.75251327204        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1172
 NPARAMETR:  1.0104E+00  4.7767E-01  1.7417E-01  1.0439E+00  2.4930E-01  9.8622E-01  1.0284E+00  1.2437E+00  1.4499E+00  4.2489E-01
             3.1720E+00
 PARAMETER:  1.1035E-01 -6.3883E-01 -1.6478E+00  1.4299E-01 -1.2891E+00  8.6120E-02  1.2797E-01  3.1809E-01  4.7152E-01 -7.5593E-01
             1.2544E+00
 GRADIENT:   1.8374E+00 -8.1384E-01 -3.4103E-01  1.2665E+00 -2.3486E+00  7.1632E-01 -5.5106E-01  1.0908E-01  2.3172E+00  2.7949E-01
             1.8072E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1234.66715742832        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1348
 NPARAMETR:  1.0063E+00  5.6390E-01  1.8771E-01  1.0340E+00  2.8161E-01  9.7934E-01  1.2344E+00  1.3265E+00  1.3727E+00  1.9659E-01
             3.1605E+00
 PARAMETER:  1.0626E-01 -4.7287E-01 -1.5728E+00  1.3346E-01 -1.1672E+00  7.9125E-02  3.1058E-01  3.8256E-01  4.1676E-01 -1.5266E+00
             1.2507E+00
 GRADIENT:   6.1225E-02 -2.8559E+00 -4.3212E+00 -6.8714E-01  1.0705E+01 -7.9426E-01  2.2657E+00  9.4512E-01  2.9442E-02  9.1892E-01
            -6.6259E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1235.01880771238        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1524
 NPARAMETR:  1.0059E+00  5.7506E-01  1.8980E-01  1.0332E+00  2.8523E-01  9.7951E-01  1.2003E+00  1.3392E+00  1.3725E+00  7.4411E-02
             3.1858E+00
 PARAMETER:  1.0591E-01 -4.5329E-01 -1.5618E+00  1.3269E-01 -1.1545E+00  7.9299E-02  2.8256E-01  3.9211E-01  4.1662E-01 -2.4982E+00
             1.2587E+00
 GRADIENT:  -1.0322E+00 -9.3146E-01 -1.1609E+00 -5.9962E-02  2.5393E+00 -5.7403E-01 -2.2527E-01  2.1720E-01 -1.5151E-01  1.0306E-01
             5.0937E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1235.07377549700        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1686
 NPARAMETR:  1.0067E+00  5.7233E-01  1.9058E-01  1.0348E+00  2.8474E-01  9.8132E-01  1.2193E+00  1.3338E+00  1.3699E+00  1.0000E-02
             3.1823E+00
 PARAMETER:  1.0664E-01 -4.5803E-01 -1.5577E+00  1.3419E-01 -1.1562E+00  8.1141E-02  2.9825E-01  3.8807E-01  4.1473E-01 -4.6875E+00
             1.2576E+00
 GRADIENT:  -1.3643E-01  1.1105E-02 -1.3043E-01 -1.1980E-01 -5.5347E-02 -8.5653E-02 -6.5584E-03  4.0090E-02 -1.3019E-01  0.0000E+00
             7.4805E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1686
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.8550E-04  1.2006E-03 -5.6757E-03 -1.1731E-02  2.9255E-04
 SE:             2.8737E-02  1.9741E-02  1.6726E-02  2.5469E-02  2.9820E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7542E-01  9.5150E-01  7.3436E-01  6.4509E-01  3.2657E-01

 ETASHRINKSD(%)  3.7285E+00  3.3865E+01  4.3966E+01  1.4677E+01  9.9001E+01
 ETASHRINKVR(%)  7.3180E+00  5.6261E+01  6.8602E+01  2.7199E+01  9.9990E+01
 EBVSHRINKSD(%)  3.7160E+00  3.3708E+01  4.3447E+01  1.3272E+01  9.9031E+01
 EBVSHRINKVR(%)  7.2940E+00  5.6054E+01  6.8018E+01  2.4783E+01  9.9991E+01
 RELATIVEINF(%)  8.7715E+01  2.9209E+00  7.2505E+00  3.9183E+01  4.9048E-04
 EPSSHRINKSD(%)  3.4141E+01
 EPSSHRINKVR(%)  5.6626E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1235.0737754970010     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -499.92294893326277     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.85
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.00
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1235.074       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  5.72E-01  1.91E-01  1.03E+00  2.85E-01  9.81E-01  1.22E+00  1.33E+00  1.37E+00  1.00E-02  3.18E+00
 


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
+        1.06E+03
 
 TH 2
+       -2.37E+01  1.02E+03
 
 TH 3
+       -2.36E+02  1.54E+03  5.79E+03
 
 TH 4
+       -3.44E+01  1.45E+02 -4.49E+02  4.14E+02
 
 TH 5
+        2.51E+02 -3.43E+03 -6.70E+03 -6.83E+01  1.44E+04
 
 TH 6
+       -2.50E-01 -8.24E+00  1.80E+00 -9.45E+00  5.14E+01  1.77E+02
 
 TH 7
+       -2.71E+00  1.70E+01 -3.98E+01 -8.38E-01  5.45E+01 -3.16E-02  2.22E+01
 
 TH 8
+        3.71E+00 -6.65E+00 -5.87E+01 -2.09E+00  1.45E+01  4.43E+00  5.83E+00  1.70E+01
 
 TH 9
+        9.62E+00 -2.50E+01  5.74E+01 -9.68E+00  1.24E+02  1.48E+00  8.06E+00 -1.42E+00  5.29E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.82E+01 -5.13E+00 -1.87E+01 -3.66E+00  2.18E+01  1.64E+00  7.19E+00  4.22E+00  6.62E+00  0.00E+00  2.61E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.927
Stop Time:
Wed Sep 29 13:29:15 CDT 2021
