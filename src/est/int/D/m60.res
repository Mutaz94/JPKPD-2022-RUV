Wed Sep 29 09:20:10 CDT 2021
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
$DATA ../../../../data/int/D/dat60.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26342.9441611638        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9372E+02  4.7060E+02  1.3120E+01  7.2020E+01  4.3943E+02 -2.2186E+03 -1.0509E+03 -3.6459E+02 -1.9717E+03 -9.2822E+02
            -5.3444E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1019.23504855742        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.8353E+00  1.6920E+00  1.0338E+00  2.1358E+00  7.1821E-01  4.2545E+00  5.5298E+00  9.7611E-01  4.1213E+00  2.5639E+00
             1.1161E+01
 PARAMETER:  7.0721E-01  6.2591E-01  1.3321E-01  8.5883E-01 -2.3100E-01  1.5480E+00  1.8102E+00  7.5824E-02  1.5162E+00  1.0415E+00
             2.5124E+00
 GRADIENT:   6.8649E+01  2.3837E+01 -2.7545E+01  6.3394E+01 -6.8299E+01  1.4430E+02  8.6700E+01  4.7485E+00  1.2035E+02  5.9964E+01
             4.5736E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1124.33339573811        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.5419E+00  1.5486E+00  7.2263E+01  3.9372E+00  3.3876E+00  3.9561E+00  1.2109E+01  4.4745E-01  3.8646E+00  2.0612E+00
             1.1000E+01
 PARAMETER:  5.3305E-01  5.3737E-01  4.3803E+00  1.4705E+00  1.3201E+00  1.4753E+00  2.5939E+00 -7.0418E-01  1.4519E+00  8.2328E-01
             2.4979E+00
 GRADIENT:   4.0914E+01  2.4584E+01 -6.0105E+00  1.2404E+02  5.5456E+01  1.3491E+02  1.0270E+02 -6.0173E-05  8.3847E+01  5.1037E+01
             4.6128E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1315.49255839379        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  8.5431E-01  1.3968E+00  1.3613E+01  9.1321E-01  2.3519E+00  3.1389E+00  4.4236E+00  1.5050E+00  2.3632E+00  5.4826E-01
             8.2236E+00
 PARAMETER: -5.7456E-02  4.3416E-01  2.7111E+00  9.2107E-03  9.5522E-01  1.2439E+00  1.5869E+00  5.0879E-01  9.6003E-01 -5.0100E-01
             2.2070E+00
 GRADIENT:  -8.1142E+01 -1.7211E+01 -1.6462E+00 -1.3708E+01  4.2336E+01  8.2274E+01  6.9669E+01 -1.1215E-02  2.4305E+01  6.3644E+00
             2.3424E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1330.87386936589        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.6054E-01  1.6361E+00  7.4842E+00  7.5450E-01  2.2727E+00  2.5849E+00  3.9101E+00  4.6350E+00  2.1697E+00  4.5777E-01
             7.7984E+00
 PARAMETER:  5.9745E-02  5.9233E-01  2.1128E+00 -1.8169E-01  9.2095E-01  1.0497E+00  1.4636E+00  1.6336E+00  8.7458E-01 -6.8138E-01
             2.1539E+00
 GRADIENT:  -6.8945E+01 -1.5586E+01 -5.6038E+00  9.8341E-02  4.2808E+01  4.7356E+00  4.1670E+01  5.5569E+00  1.8738E+01  5.2224E+00
             1.2181E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1345.22731301455        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  1.2190E+00  1.6416E+00  7.4547E+00  7.6635E-01  2.0980E+00  2.5775E+00  3.9399E+00  3.5203E+00  8.5471E-01  1.7268E-01
             7.7537E+00
 PARAMETER:  2.9800E-01  5.9567E-01  2.1088E+00 -1.6612E-01  8.4098E-01  1.0468E+00  1.4712E+00  1.3585E+00 -5.6991E-02 -1.6563E+00
             2.1482E+00
 GRADIENT:   1.6515E+01 -2.0781E+00  2.6682E+00  5.3981E+00 -1.8878E+00 -1.4564E+00  2.7589E+01 -4.4922E-01 -9.9548E-01  6.9617E-01
             7.9707E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1350.06016301194        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      630
 NPARAMETR:  1.3568E+00  1.6478E+00  7.4489E+00  7.9396E-01  2.1198E+00  2.6307E+00  4.2314E+00  3.7669E+00  1.2439E+00  3.2521E-02
             7.5836E+00
 PARAMETER:  4.0513E-01  5.9942E-01  2.1081E+00 -1.3072E-01  8.5132E-01  1.0672E+00  1.5425E+00  1.4263E+00  3.1826E-01 -3.3259E+00
             2.1260E+00
 GRADIENT:   2.2824E+01 -1.3333E+01 -2.1470E-01 -4.0738E+00  6.6534E+00 -5.8819E+01 -3.2448E+01  7.8804E-01  3.2481E+00  2.5451E-02
             2.9851E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1357.35321323434        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      806
 NPARAMETR:  1.2471E+00  1.6632E+00  7.4307E+00  8.9218E-01  2.0876E+00  2.8113E+00  5.0459E+00  2.8069E+00  9.3999E-01  1.0000E-02
             7.5773E+00
 PARAMETER:  3.2080E-01  6.0876E-01  2.1056E+00 -1.4085E-02  8.3603E-01  1.1336E+00  1.7186E+00  1.1321E+00  3.8109E-02 -7.2620E+00
             2.1252E+00
 GRADIENT:  -2.5416E+00  5.2179E+00  1.1151E+00  5.1592E+00  4.8120E-01 -2.6654E+01 -3.2901E+00 -1.2203E+00 -7.8639E-01  0.0000E+00
            -7.7348E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1358.18149505975        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      984
 NPARAMETR:  1.2494E+00  1.6592E+00  7.4257E+00  8.6884E-01  2.0791E+00  2.9037E+00  5.0609E+00  2.8658E+00  1.0133E+00  1.0000E-02
             7.5762E+00
 PARAMETER:  3.2264E-01  6.0633E-01  2.1050E+00 -4.0595E-02  8.3192E-01  1.1660E+00  1.7215E+00  1.1529E+00  1.1324E-01 -7.8404E+00
             2.1250E+00
 GRADIENT:  -1.9066E+00  2.9927E+00  1.5655E+00 -5.6740E+00 -1.9960E+00 -1.4285E+01  3.3310E+00 -1.3448E+00  1.2796E+00  0.0000E+00
             7.7816E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1358.93823095422        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1166
 NPARAMETR:  1.2633E+00  1.5540E+00  7.3474E+00  8.8438E-01  2.0799E+00  3.0622E+00  5.0649E+00  3.0868E+00  9.5824E-01  1.0000E-02
             7.5785E+00
 PARAMETER:  3.3371E-01  5.4085E-01  2.0943E+00 -2.2871E-02  8.3232E-01  1.2191E+00  1.7223E+00  1.2271E+00  5.7340E-02 -8.3957E+00
             2.1253E+00
 GRADIENT:   7.3563E-01 -6.5444E-01 -1.9761E-01 -2.6527E+00  2.3443E+00  5.6547E+00 -1.2829E+00  2.8620E-01  8.8707E-02  0.0000E+00
             1.3345E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1358.96103044651        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1341
 NPARAMETR:  1.2589E+00  1.5267E+00  7.4276E+00  9.0497E-01  2.0694E+00  3.0105E+00  5.1489E+00  3.0225E+00  9.7197E-01  1.0000E-02
             7.5782E+00
 PARAMETER:  3.3025E-01  5.2312E-01  2.1052E+00  1.4825E-04  8.2725E-01  1.2021E+00  1.7388E+00  1.2061E+00  7.1574E-02 -8.3957E+00
             2.1253E+00
 GRADIENT:  -7.6508E-02  2.2697E-01  1.4036E-02 -5.7668E-01  5.1855E-01 -5.2506E-01 -8.2879E-01  2.9391E-02 -9.9521E-02  0.0000E+00
             1.2570E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1358.96343292644        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1518
 NPARAMETR:  1.2571E+00  1.5038E+00  7.4328E+00  9.1413E-01  2.0634E+00  3.0055E+00  5.1707E+00  3.0043E+00  9.9027E-01  1.0000E-02
             7.5741E+00
 PARAMETER:  3.2880E-01  5.0801E-01  2.1059E+00  1.0216E-02  8.2434E-01  1.2004E+00  1.7430E+00  1.2001E+00  9.0222E-02 -8.3957E+00
             2.1247E+00
 GRADIENT:  -3.9207E-01 -5.3066E-02  4.3201E-03  1.9066E-02 -5.4028E-02 -1.1172E+00 -1.3881E+00 -1.0260E-02 -3.4225E-02  0.0000E+00
            -4.8148E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1359.00743560919        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1697
 NPARAMETR:  1.2567E+00  1.4735E+00  7.7089E+00  9.2865E-01  2.0657E+00  3.0185E+00  5.2629E+00  3.0460E+00  1.0038E+00  1.0000E-02
             7.5764E+00
 PARAMETER:  3.2853E-01  4.8764E-01  2.1424E+00  2.5972E-02  8.2548E-01  1.2048E+00  1.7607E+00  1.2138E+00  1.0382E-01 -8.3957E+00
             2.1250E+00
 GRADIENT:  -4.5175E-01  5.4961E-02 -2.7390E-02 -2.1567E-01  3.2662E-02  5.3157E-01  4.8467E-01  4.8088E-04  4.2198E-02  0.0000E+00
             5.9541E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1359.04123073314        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1864             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2615E+00  1.4746E+00  7.7007E+00  9.2732E-01  2.0660E+00  3.0604E+00  5.2659E+00  3.0511E+00  1.0097E+00  1.0000E-02
             7.5729E+00
 PARAMETER:  3.3226E-01  4.8841E-01  2.1413E+00  2.4540E-02  8.2563E-01  1.2185E+00  1.7613E+00  1.2155E+00  1.0963E-01 -8.3957E+00
             2.1246E+00
 GRADIENT:   3.3482E+01  1.2116E+01  3.9158E-01  2.4626E+00  4.8515E+00  8.2974E+01  1.1900E+02  2.6508E-01  4.3156E-01  0.0000E+00
             3.8573E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1359.04282864599        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2048             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2615E+00  1.4714E+00  7.7176E+00  9.2880E-01  2.0658E+00  3.0595E+00  5.2666E+00  3.0468E+00  1.0054E+00  1.0000E-02
             7.5731E+00
 PARAMETER:  3.3229E-01  4.8619E-01  2.1435E+00  2.6137E-02  8.2553E-01  1.2183E+00  1.7614E+00  1.2141E+00  1.0534E-01 -8.3957E+00
             2.1246E+00
 GRADIENT:   3.3484E+01  1.2034E+01  4.3496E-01  2.9165E+00  4.6794E+00  8.2846E+01  1.1873E+02  2.3108E-01  2.8517E-01  0.0000E+00
             3.8276E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1359.04305679676        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2231
 NPARAMETR:  1.2614E+00  1.4692E+00  7.7276E+00  9.2948E-01  2.0658E+00  3.0603E+00  5.2702E+00  3.0466E+00  1.0058E+00  1.0000E-02
             7.5731E+00
 PARAMETER:  3.3226E-01  4.8473E-01  2.1448E+00  2.6875E-02  8.2551E-01  1.2185E+00  1.7621E+00  1.2140E+00  1.0579E-01 -8.3957E+00
             2.1246E+00
 GRADIENT:   4.6121E-01 -2.9313E-02 -1.6592E-02 -3.8741E-01 -4.4990E-03  5.5330E+00  6.1589E-01 -1.1159E-02  6.2356E-02  0.0000E+00
             7.4304E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1359.04334282408        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     2404
 NPARAMETR:  1.2614E+00  1.4685E+00  7.7333E+00  9.3026E-01  2.0658E+00  3.0603E+00  5.2713E+00  3.0476E+00  1.0047E+00  1.0000E-02
             7.5728E+00
 PARAMETER:  3.3226E-01  4.8427E-01  2.1455E+00  2.7709E-02  8.2549E-01  1.2185E+00  1.7623E+00  1.2143E+00  1.0473E-01 -8.3957E+00
             2.1246E+00
 GRADIENT:   4.6118E-01  1.2197E-02 -2.1794E-02 -1.4636E-01 -2.6779E-02  5.5336E+00  5.0280E-01 -4.4868E-03  9.2174E-03  0.0000E+00
            -1.4762E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2404
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1783E-03  1.8348E-02 -1.7747E-02 -5.0889E-02 -4.1846E-05
 SE:             2.9529E-02  2.6264E-02  8.3051E-03  1.1698E-02  1.3478E-04
 N:                     100         100         100         100         100

 P VAL.:         8.0793E-01  4.8481E-01  3.2606E-02  1.3611E-05  7.5621E-01

 ETASHRINKSD(%)  1.0750E+00  1.2012E+01  7.2177E+01  6.0810E+01  9.9548E+01
 ETASHRINKVR(%)  2.1385E+00  2.2582E+01  9.2259E+01  8.4642E+01  9.9998E+01
 EBVSHRINKSD(%)  1.1773E+00  6.8071E+00  7.7302E+01  6.8176E+01  9.9437E+01
 EBVSHRINKVR(%)  2.3407E+00  1.3151E+01  9.4848E+01  8.9873E+01  9.9997E+01
 RELATIVEINF(%)  9.7592E+01  4.4557E+01  1.8678E+00  4.9990E+00  1.2235E-03
 EPSSHRINKSD(%)  8.0542E+00
 EPSSHRINKVR(%)  1.5460E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1359.0433428240804     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       295.04601694433040     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    79.96
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    17.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1359.043       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.26E+00  1.47E+00  7.73E+00  9.30E-01  2.07E+00  3.06E+00  5.27E+00  3.05E+00  1.00E+00  1.00E-02  7.57E+00
 


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
+        4.11E+01
 
 TH 2
+        1.62E+00  3.84E+00
 
 TH 3
+        7.20E-02 -1.85E-01  9.49E-03
 
 TH 4
+       -1.20E+01  2.73E+01 -1.40E+00  2.08E+02
 
 TH 5
+        3.64E-01  1.09E+00 -5.33E-02  7.82E+00  3.16E-01
 
 TH 6
+        5.26E+00 -2.42E-01  3.12E-02 -4.83E+00 -7.77E-02  7.29E-01
 
 TH 7
+        2.36E+00 -2.26E+00  1.21E-01 -1.80E+01 -6.50E-01  5.82E-01  1.60E+00
 
 TH 8
+       -1.71E-01  2.47E-01 -1.29E-02  1.91E+00  7.14E-02 -5.15E-02 -1.67E-01  1.78E-02
 
 TH 9
+        2.35E+00 -3.62E+00  1.88E-01 -2.80E+01 -1.03E+00  7.45E-01  2.45E+00 -2.58E-01  3.79E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.03E+00 -1.03E+00  4.31E-02 -6.91E+00 -2.39E-01  1.22E-02  5.50E-01 -5.59E-02  9.63E-01  0.00E+00  6.79E-01
 
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
+        7.18E+01
 
 TH 2
+       -2.27E-01  1.62E+01
 
 TH 3
+        5.32E-02  1.83E-01  1.97E-01
 
 TH 4
+       -1.86E+00  2.48E+01 -1.44E+00  1.85E+02
 
 TH 5
+       -9.09E-01 -5.14E+00 -2.47E+00  1.18E+01  7.34E+01
 
 TH 6
+       -8.78E-03  6.11E-02  1.78E-03  4.02E-01 -3.78E-01  1.88E+01
 
 TH 7
+        1.78E-01  1.81E+00 -1.11E-01 -1.53E+01  1.35E+00  6.30E-03  4.73E+00
 
 TH 8
+       -6.15E-02 -2.11E-01 -3.16E-01  1.73E+00  1.10E+00  2.04E-02  1.85E-01  1.01E+00
 
 TH 9
+       -2.77E-02 -2.11E+00 -1.08E-01 -2.33E+01  3.62E+00 -2.39E-02  2.15E+00 -2.65E-01  9.79E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.20E+00 -1.91E+00 -1.18E-02 -9.72E+00 -7.96E-01  6.83E-01  6.02E-01  1.32E-01  3.39E+00  0.00E+00  1.81E+01
 
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
+        7.29E+01
 
 TH 2
+        3.01E+01  1.63E+01
 
 TH 3
+        3.34E-02  1.49E-01  1.16E-01
 
 TH 4
+        4.02E+01  2.55E+01 -1.08E-01  1.78E+02
 
 TH 5
+       -1.21E+00 -4.33E+00 -1.90E+00 -7.61E+00  5.40E+01
 
 TH 6
+        1.86E+01  7.79E+00 -3.23E-02 -2.00E+01  2.73E+00  2.50E+01
 
 TH 7
+        5.74E+00  1.37E+00 -1.25E-01 -1.55E+01  2.79E+00  8.03E+00  4.06E+00
 
 TH 8
+       -2.10E-01 -1.40E-01 -1.54E-01  5.15E-01  7.63E-01  5.76E-02  9.26E-02  4.39E-01
 
 TH 9
+       -8.17E+00 -2.80E+00 -1.35E-01 -4.98E+01  8.11E+00  8.62E+00  3.79E+00 -1.84E-01  2.85E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.92E+01 -2.17E+01 -3.13E-01 -1.61E+02  2.85E+01  2.64E+01  1.22E+01 -7.59E-01  8.04E+01  0.00E+00  7.82E+02
 
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
 #CPUT: Total CPU Time in Seconds,       97.195
Stop Time:
Wed Sep 29 09:21:49 CDT 2021
