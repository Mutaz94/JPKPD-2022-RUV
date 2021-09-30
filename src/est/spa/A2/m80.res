Wed Sep 29 13:03:57 CDT 2021
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
$DATA ../../../../data/spa/A2/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -947.253469068550        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7179E+02  3.6383E+01  7.9207E+01  5.3687E+00  2.5172E+01  7.2781E+01  1.1194E+01 -4.0120E+01  1.8804E+01 -2.2386E+01
            -1.4226E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1429.02726645142        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0406E+00  1.0596E+00  9.1951E-01  1.0765E+00  9.6538E-01  8.1937E-01  8.5325E-01  1.0408E+00  7.8942E-01  8.6665E-01
             2.0342E+00
 PARAMETER:  1.3982E-01  1.5792E-01  1.6090E-02  1.7368E-01  6.4769E-02 -9.9217E-02 -5.8700E-02  1.3994E-01 -1.3646E-01 -4.3124E-02
             8.1013E-01
 GRADIENT:   1.6166E+02  1.0331E+02  1.9237E+01  1.3876E+02 -1.4393E+01 -2.0158E+01  3.8530E-01 -8.6125E+00 -1.0497E+01  4.0328E+00
            -1.9603E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1448.18757622359        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0143E+00  1.0392E+00  1.1331E+00  1.0490E+00  1.0311E+00  9.0852E-01  8.5013E-01  1.5445E+00  9.3175E-01  1.1616E-01
             2.2688E+00
 PARAMETER:  1.1424E-01  1.3841E-01  2.2495E-01  1.4780E-01  1.3065E-01  4.0585E-03 -6.2364E-02  5.3471E-01  2.9306E-02 -2.0528E+00
             9.1925E-01
 GRADIENT:   5.2098E+01  3.9759E+01  2.9781E+01  4.1440E+01 -1.9004E+01  1.9844E+01  7.8823E+00 -1.2910E+01  2.2920E+01 -1.3062E-01
            -1.2672E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1464.83756949275        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.0086E+00  1.0435E+00  3.0150E-01  9.3696E-01  5.5950E-01  8.5646E-01  7.8982E-01  6.1724E-01  7.8747E-01  1.0914E-01
             2.5916E+00
 PARAMETER:  1.0857E-01  1.4261E-01 -1.0990E+00  3.4888E-02 -4.8070E-01 -5.4950E-02 -1.3595E-01 -3.8249E-01 -1.3894E-01 -2.1151E+00
             1.0523E+00
 GRADIENT:  -1.7548E+01  2.6838E+01 -6.6499E+00  9.0579E+01  3.2681E+01 -6.7495E+00  6.6859E-01 -1.0770E+00  7.9422E+00  2.8528E-01
             9.0925E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1465.06242620722        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.0087E+00  1.0382E+00  2.2950E-01  9.0816E-01  4.9425E-01  8.5526E-01  7.8516E-01  5.0169E-01  7.7460E-01  1.1374E-01
             2.6115E+00
 PARAMETER:  1.0869E-01  1.3744E-01 -1.3718E+00  3.6706E-03 -6.0471E-01 -5.6355E-02 -1.4187E-01 -5.8977E-01 -1.5540E-01 -2.0738E+00
             1.0599E+00
 GRADIENT:  -7.5174E+00  6.9772E+01 -1.2023E+00  1.2048E+02  4.2004E+00 -8.9112E+00  2.3855E+00 -1.9759E+00  1.9786E+00  4.2576E-01
             2.5028E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1465.16152724615        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0094E+00  1.0211E+00  1.9828E-01  8.8402E-01  4.5911E-01  8.5892E-01  7.9266E-01  4.5878E-01  7.4689E-01  1.1089E-01
             2.6259E+00
 PARAMETER:  1.0937E-01  1.2092E-01 -1.5181E+00 -2.3272E-02 -6.7847E-01 -5.2077E-02 -1.3235E-01 -6.7919E-01 -1.9183E-01 -2.0993E+00
             1.0654E+00
 GRADIENT:   5.3294E+00  8.9570E+01  3.6952E+00  1.0788E+02 -2.6113E+01 -7.9480E+00  3.5678E+00 -2.5614E+00 -4.6620E+00  4.7194E-01
             3.4130E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1469.03138496786        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  1.0381E+00  1.1395E+00  2.2978E-01  8.2645E-01  5.3441E-01  8.7282E-01  7.3831E-01  5.8933E-01  7.1387E-01  1.2918E-01
             2.6269E+00
 PARAMETER:  1.3743E-01  2.3061E-01 -1.3706E+00 -9.0621E-02 -5.2659E-01 -3.6025E-02 -2.0339E-01 -4.2878E-01 -2.3705E-01 -1.9465E+00
             1.0658E+00
 GRADIENT:   2.2651E+01  5.9631E+01  2.0919E+01  1.6281E+01 -5.6511E+01 -1.6341E+00 -4.9616E-01 -3.3262E+00 -7.1951E+00  5.3752E-01
             1.2720E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1492.93733589296        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      717
 NPARAMETR:  1.0225E+00  1.5934E+00  1.3469E-01  5.6228E-01  7.0380E-01  8.5902E-01  6.2065E-01  2.1026E+00  9.1382E-01  1.0000E-02
             2.5617E+00
 PARAMETER:  1.2221E-01  5.6585E-01 -1.9048E+00 -4.7576E-01 -2.5126E-01 -5.1965E-02 -3.7700E-01  8.4317E-01  9.8766E-03 -4.6424E+00
             1.0407E+00
 GRADIENT:   9.1119E+00  5.3852E+01 -5.9200E+00  3.9098E+01 -2.5937E+01 -1.1117E+00  3.2016E+00  5.3970E+00 -4.3863E+00  0.0000E+00
             7.4572E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1501.02722434373        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      892
 NPARAMETR:  1.0134E+00  1.5193E+00  1.4644E-01  5.8250E-01  6.8013E-01  8.6009E-01  6.4156E-01  2.1037E+00  1.0238E+00  1.1705E-02
             2.1243E+00
 PARAMETER:  1.1327E-01  5.1822E-01 -1.8212E+00 -4.4043E-01 -2.8547E-01 -5.0719E-02 -3.4386E-01  8.4369E-01  1.2354E-01 -4.3478E+00
             8.5343E-01
 GRADIENT:   5.4322E+00  1.0282E+01 -8.0335E-03  1.5267E+01 -1.0738E+01 -4.5446E-01  8.8199E-01  3.7968E+00 -9.1128E-01  5.7553E-03
             4.4464E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1506.56491973982        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1070
 NPARAMETR:  1.0099E+00  2.1131E+00  4.7292E-02  2.5687E-01  9.7504E-01  8.4839E-01  5.8752E-01  1.3807E+00  2.2898E+00  1.0000E-02
             2.2964E+00
 PARAMETER:  1.0981E-01  8.4816E-01 -2.9514E+00 -1.2592E+00  7.4726E-02 -6.4421E-02 -4.3185E-01  4.2261E-01  9.2846E-01 -1.1348E+01
             9.3135E-01
 GRADIENT:  -1.4091E+01  9.4734E+01 -2.2316E-01  4.0703E+00 -5.9186E+01 -6.3466E+00  1.0014E+01  5.5711E+00 -5.0050E+00  0.0000E+00
             2.0195E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1508.49259618995        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  1.0098E+00  2.1817E+00  3.8891E-02  2.2585E-01  1.0251E+00  8.4757E-01  5.7470E-01  1.2141E+00  2.7151E+00  1.0000E-02
             2.2778E+00
 PARAMETER:  1.0978E-01  8.8010E-01 -3.1470E+00 -1.3879E+00  1.2480E-01 -6.5383E-02 -4.5390E-01  2.9398E-01  1.0988E+00 -1.2606E+01
             9.2320E-01
 GRADIENT:  -1.4329E+01  1.0739E+02 -3.0944E+00  6.7518E+00 -5.7496E+01 -6.8906E+00  7.3206E+00  5.4462E+00 -4.1899E+00  0.0000E+00
             1.3273E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1512.11254045239        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:     1374
 NPARAMETR:  1.0110E+00  2.1518E+00  4.0040E-02  2.2864E-01  1.0271E+00  8.5680E-01  5.7166E-01  1.8266E-01  2.6909E+00  1.0000E-02
             2.2981E+00
 PARAMETER:  1.1092E-01  8.6630E-01 -3.1179E+00 -1.3756E+00  1.2677E-01 -5.4548E-02 -4.5922E-01 -1.6001E+00  1.0899E+00 -1.2622E+01
             9.3209E-01
 GRADIENT:   7.5419E+01  2.8392E+02  9.9244E+00  1.6063E+01 -2.0558E+01  2.3659E+00  1.1637E+01  1.4639E-01 -1.0776E-01  0.0000E+00
             2.7697E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1512.23303845020        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:     1499
 NPARAMETR:  1.0111E+00  2.1487E+00  4.0165E-02  2.2900E-01  1.0270E+00  8.4779E-01  5.7198E-01  4.1494E-02  2.6948E+00  1.0000E-02
             2.2940E+00
 PARAMETER:  1.1106E-01  8.6486E-01 -3.1148E+00 -1.3740E+00  1.2661E-01 -6.5119E-02 -4.5864E-01 -3.0822E+00  1.0913E+00 -1.2622E+01
             9.3030E-01
 GRADIENT:   7.6526E+01  2.7889E+02  9.9552E+00  1.6171E+01 -1.7991E+01 -1.5348E+00  1.1662E+01  8.6414E-03  2.8093E-01  0.0000E+00
             2.6778E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1512.44321077371        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:     1607
 NPARAMETR:  1.0112E+00  2.1369E+00  4.0199E-02  2.2918E-01  1.0268E+00  8.4375E-01  5.7228E-01  1.0000E-02  2.6986E+00  1.0000E-02
             2.2878E+00
 PARAMETER:  1.1114E-01  8.5936E-01 -3.1139E+00 -1.3732E+00  1.2646E-01 -6.9893E-02 -4.5813E-01 -6.1075E+00  1.0927E+00 -1.2622E+01
             9.2761E-01
 GRADIENT:  -6.8487E+00  1.8363E+01 -3.2870E+00  3.5605E+00 -1.1058E+01 -8.7126E+00  8.0122E+00  0.0000E+00 -4.0176E+00  0.0000E+00
             1.8750E+01

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1512.44321077371        NO. OF FUNC. EVALS.:  32
 CUMULATIVE NO. OF FUNC. EVALS.:     1639
 NPARAMETR:  1.0123E+00  2.1545E+00  4.1444E-02  2.3227E-01  1.0255E+00  8.4458E-01  5.7484E-01  1.0000E-02  2.7278E+00  1.0000E-02
             2.2667E+00
 PARAMETER:  1.1114E-01  8.5936E-01 -3.1139E+00 -1.3732E+00  1.2646E-01 -6.9893E-02 -4.5813E-01 -6.1075E+00  1.0927E+00 -1.2622E+01
             9.2761E-01
 GRADIENT:  -6.4684E+03 -1.6519E+03 -4.5200E+02 -1.0496E+03  1.1352E+04 -7.1918E+03 -3.1274E+03  0.0000E+00 -6.7420E+02  0.0000E+00
             7.8857E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1639
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.5559E-03 -2.4818E-02  1.9900E-04  9.3416E-03 -2.3321E-04
 SE:             2.9832E-02  2.4706E-02  7.2151E-05  2.0205E-02  2.3818E-04
 N:                     100         100         100         100         100

 P VAL.:         9.0512E-01  3.1512E-01  5.8138E-03  6.4383E-01  3.2750E-01

 ETASHRINKSD(%)  5.7881E-02  1.7232E+01  9.9758E+01  3.2311E+01  9.9202E+01
 ETASHRINKVR(%)  1.1573E-01  3.1494E+01  9.9999E+01  5.4182E+01  9.9994E+01
 EBVSHRINKSD(%)  2.4585E+00  1.5914E+01  9.9733E+01  2.6849E+01  9.9161E+01
 EBVSHRINKVR(%)  4.8566E+00  2.9295E+01  9.9999E+01  4.6489E+01  9.9993E+01
 RELATIVEINF(%)  9.1600E+01  2.2123E+01  5.6171E-04  2.3534E+01  2.0234E-03
 EPSSHRINKSD(%)  3.1873E+01
 EPSSHRINKVR(%)  5.3587E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1512.4432107737114     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -777.29238420997319     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.07
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1512.443       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  2.14E+00  4.02E-02  2.29E-01  1.03E+00  8.44E-01  5.72E-01  1.00E-02  2.70E+00  1.00E-02  2.29E+00
 


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
+        2.84E+06
 
 TH 2
+       -1.73E+02  5.32E+02
 
 TH 3
+       -2.21E+03  1.18E+03  2.24E+06
 
 TH 4
+        1.33E+04  9.67E+02  9.37E+03  3.69E+05
 
 TH 5
+       -3.71E+02 -3.96E+02 -3.30E+02 -1.11E+04  2.13E+06
 
 TH 6
+       -9.91E+02 -7.05E+01 -1.01E+03 -3.76E+02  8.58E+02  5.04E+06
 
 TH 7
+       -6.41E+02 -6.02E+01 -5.70E+02  5.72E+03  6.03E+02 -4.20E+02  5.23E+05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.02E+01 -6.52E+03  3.86E+03  5.39E+02 -2.27E+01 -3.80E+01 -2.30E+01  0.00E+00  4.30E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        7.03E+01  9.25E+03  1.43E+03 -7.17E+02  6.75E+00  5.62E+01  5.72E+01  0.00E+00 -2.20E+02  0.00E+00  7.97E+03
 
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
 #CPUT: Total CPU Time in Seconds,       27.313
Stop Time:
Wed Sep 29 13:04:26 CDT 2021
