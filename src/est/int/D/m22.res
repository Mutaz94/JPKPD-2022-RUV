Wed Sep 29 08:15:06 CDT 2021
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
$DATA ../../../../data/int/D/dat22.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m22.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1356.32714772719        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6028E+02  4.1580E+01 -9.0579E+01 -3.4247E+00  4.1828E+02 -1.2544E+03 -3.7330E+02 -7.4987E+01 -6.5842E+02 -5.6769E+02
            -7.5775E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1802.68425677492        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      109
 NPARAMETR:  1.5577E+00  1.2017E+00  1.0448E+00  2.3909E+00  4.3424E-01  6.5772E+00  2.5042E+00  1.2664E+00  7.4281E+00  5.2883E+00
             5.1322E+00
 PARAMETER:  5.4320E-01  2.8370E-01  1.4386E-01  9.7166E-01 -7.3416E-01  1.9836E+00  1.0180E+00  3.3617E-01  2.1053E+00  1.7655E+00
             1.7355E+00
 GRADIENT:   1.5283E+01 -2.8385E+01 -1.4154E+01  3.6029E+01 -6.2850E+01  1.6576E+02  7.3338E+01  9.0964E+00  1.3120E+02  9.2458E+01
             7.6089E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2045.53328597939        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      284
 NPARAMETR:  2.5721E+00  3.7784E+00  8.4232E+02  9.9854E-02  4.1521E+00  1.7731E+01  2.6450E+00  1.7667E+01  2.0385E+01  3.7102E+00
             3.2115E+00
 PARAMETER:  1.0447E+00  1.4293E+00  6.8362E+00 -2.2040E+00  1.5236E+00  2.9753E+00  1.0727E+00  2.9717E+00  3.1148E+00  1.4111E+00
             1.2667E+00
 GRADIENT:   2.3175E+01  4.3186E+00 -1.7598E+00  1.3107E+01  8.8276E+01  1.3123E+02  2.5160E+01  3.9560E+00  2.2705E+02  9.3544E+01
             2.1250E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2072.61040131665        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      466
 NPARAMETR:  3.8509E-01  3.6213E+00  1.0477E+02  6.6439E-01  3.3908E+00  1.5905E+01  2.5792E+00  5.6984E+01  1.5502E+01  2.2188E+00
             3.2212E+00
 PARAMETER: -8.5428E-01  1.3868E+00  4.7517E+00 -3.0889E-01  1.3211E+00  2.8666E+00  1.0475E+00  4.1428E+00  2.8410E+00  8.9697E-01
             1.2697E+00
 GRADIENT:   2.0013E+00  3.5687E+01 -2.0071E+00  2.0752E+01  1.1077E+02  8.9439E+02  3.1785E+01  1.4391E+01  1.1634E+02  6.7914E+01
             1.5776E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2220.87013411132        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.1032E+00  3.6207E+00  1.0359E+02  7.2059E-01  3.2933E+00  1.8589E+00  2.5787E+00  5.5561E+01  2.8235E+00  2.1850E+00
             3.2127E+00
 PARAMETER:  1.9819E-01  1.3867E+00  4.7405E+00 -2.2769E-01  1.2919E+00  7.2000E-01  1.0473E+00  4.1175E+00  1.1380E+00  8.8159E-01
             1.2671E+00
 GRADIENT:   4.8716E+01  4.1513E+02 -1.5554E+00  9.3957E+01  1.4659E+02 -1.0916E+02  2.5174E+01  1.0266E+02 -2.4119E+00  8.9988E+01
             1.0229E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2339.45463189716        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  1.0487E+00  3.1387E+00  1.0384E+02  7.2027E-01  2.4509E+00  2.0758E+00  2.5685E+00  3.1516E+01  3.6445E+00  2.0219E+00
             3.1709E+00
 PARAMETER:  1.4752E-01  1.2438E+00  4.7428E+00 -2.2812E-01  9.9644E-01  8.3033E-01  1.0433E+00  3.5505E+00  1.3932E+00  8.0403E-01
             1.2540E+00
 GRADIENT:   1.2592E+01  3.6800E+02 -3.1353E-01  6.8251E+01  7.8673E+01  2.8469E+01  5.5518E+01  2.1942E+01  2.5739E+01  7.5222E+01
             1.8333E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2358.49880872205        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      785
 NPARAMETR:  1.0478E+00  2.3317E+00  1.0497E+02  7.1953E-01  1.9030E+00  1.9477E+00  2.5358E+00  2.0721E+01  3.2213E+00  1.6976E+00
             3.0673E+00
 PARAMETER:  1.4671E-01  9.4658E-01  4.7537E+00 -2.2915E-01  7.4344E-01  7.6664E-01  1.0305E+00  3.1311E+00  1.2698E+00  6.2924E-01
             1.2208E+00
 GRADIENT:   1.2520E+01  2.2600E+02  2.0446E+00  4.8966E+01 -3.8151E+01 -2.2771E+01  5.9777E+01  2.5431E+01  4.4894E+01  5.4140E+01
             1.3451E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2420.52955756818        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  1.2243E+00  2.0165E+00  1.0475E+02  7.1917E-01  1.8861E+00  2.5299E+00  2.5575E+00  1.4943E+01  3.3044E+00  1.5421E+00
             2.9966E+00
 PARAMETER:  3.0236E-01  8.0135E-01  4.7516E+00 -2.2966E-01  7.3452E-01  1.0282E+00  1.0390E+00  2.8043E+00  1.2953E+00  5.3316E-01
             1.1975E+00
 GRADIENT:   2.2669E+01 -1.0536E+00  1.4834E+00  1.3136E+01 -7.2991E+01 -4.5576E+01 -1.2297E+02  8.1699E-01  8.0247E+00  2.2063E+01
             7.7568E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2442.05545604063        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1097
 NPARAMETR:  1.0797E+00  1.8968E+00  1.0349E+02  7.1880E-01  2.1073E+00  2.7269E+00  2.7526E+00  1.3824E+01  2.7256E+00  1.3936E+00
             2.8970E+00
 PARAMETER:  1.7668E-01  7.4018E-01  4.7394E+00 -2.3017E-01  8.4543E-01  1.1032E+00  1.1125E+00  2.7264E+00  1.1027E+00  4.3187E-01
             1.1637E+00
 GRADIENT:  -1.8131E+01 -7.2032E+00  9.9155E-01  7.0758E+00  1.2494E+01 -3.3761E+00 -1.0299E+02  2.6707E-01  9.4486E+00 -5.3747E+00
             2.4935E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2445.73205783431        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  1.1327E+00  1.9070E+00  1.0514E+02  7.1776E-01  2.1086E+00  2.6877E+00  2.8508E+00  1.3629E+01  2.7293E+00  1.3975E+00
             2.9005E+00
 PARAMETER:  2.2464E-01  7.4553E-01  4.7553E+00 -2.3162E-01  8.4601E-01  1.0887E+00  1.1476E+00  2.7122E+00  1.1041E+00  4.3472E-01
             1.1649E+00
 GRADIENT:  -3.6962E+00 -3.9892E+00  1.1237E+00  5.2169E+00  1.2620E+01 -1.1387E+01 -8.6011E+01 -1.8122E-02  1.0860E+01 -4.5550E+00
             4.9938E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2459.33414486765        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  1.1390E+00  1.8788E+00  8.0559E+01  7.2800E-01  2.0705E+00  2.7409E+00  3.5520E+00  1.2979E+01  2.2648E+00  1.4075E+00
             2.8700E+00
 PARAMETER:  2.3017E-01  7.3065E-01  4.4890E+00 -2.1745E-01  8.2777E-01  1.1083E+00  1.3675E+00  2.6633E+00  9.1747E-01  4.4183E-01
             1.1543E+00
 GRADIENT:   1.2238E+02  1.3859E+02  1.9272E+00  1.1187E+01  2.6948E+01  3.2519E+02  3.2214E+02  3.3289E+00  3.3391E+01 -7.9044E-01
             1.7476E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2461.13468273976        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1589             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1384E+00  1.8658E+00  8.1220E+01  7.3102E-01  2.0711E+00  2.7497E+00  3.5585E+00  1.3365E+01  2.0077E+00  1.4086E+00
             2.8716E+00
 PARAMETER:  2.2966E-01  7.2371E-01  4.4972E+00 -2.1331E-01  8.2809E-01  1.1115E+00  1.3693E+00  2.6927E+00  7.9698E-01  4.4262E-01
             1.1549E+00
 GRADIENT:   1.2182E+02  1.3472E+02  1.6472E+00  1.0591E+01  2.6979E+01  3.2747E+02  3.1913E+02  3.8263E+00  2.6390E+01 -7.9897E-01
             1.6806E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2464.48346803385        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1683
 NPARAMETR:  1.1404E+00  1.8808E+00  8.4818E+01  7.2906E-01  2.0917E+00  2.7136E+00  3.5202E+00  1.2982E+01  1.2338E+00  1.4167E+00
             2.9133E+00
 PARAMETER:  2.3136E-01  7.3171E-01  4.5405E+00 -2.1600E-01  8.3796E-01  1.0983E+00  1.3585E+00  2.6635E+00  3.1008E-01  4.4830E-01
             1.1693E+00
 GRADIENT:  -1.5977E+00  1.3300E+01  1.0514E+00 -3.2098E+00  4.0522E+00 -7.2964E+00 -2.5419E+01  6.7762E-01  4.3441E+00 -2.6076E+00
             6.4052E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2465.21109086445        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1862
 NPARAMETR:  1.1288E+00  1.8793E+00  8.5062E+01  7.2884E-01  2.0926E+00  2.6415E+00  3.6457E+00  1.2921E+01  1.0474E+00  1.4177E+00
             2.9085E+00
 PARAMETER:  2.2118E-01  7.3093E-01  4.5434E+00 -2.1630E-01  8.3841E-01  1.0714E+00  1.3936E+00  2.6588E+00  1.4628E-01  4.4906E-01
             1.1676E+00
 GRADIENT:  -5.0894E+00  1.4341E+01 -1.9710E-01 -3.2182E+00  2.8239E+00 -2.0182E+01 -1.5710E+01  3.0013E+00  2.7075E+00 -4.5684E+00
             3.2455E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2468.83645119230        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2043
 NPARAMETR:  1.1644E+00  1.7738E+00  8.7856E+01  7.3369E-01  2.0358E+00  2.7571E+00  3.9402E+00  1.1901E+01  6.8944E-01  1.4038E+00
             2.8910E+00
 PARAMETER:  2.5223E-01  6.7314E-01  4.5757E+00 -2.0967E-01  8.1090E-01  1.1142E+00  1.4712E+00  2.5766E+00 -2.7188E-01  4.3920E-01
             1.1616E+00
 GRADIENT:   5.0009E+00  9.5503E+00  1.1911E+00 -1.1252E+01 -1.9441E+01 -4.7585E-01  2.3889E+00  3.9220E-01  1.9141E-01 -5.0501E+00
            -1.8364E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2471.00591470883        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2205
 NPARAMETR:  1.1445E+00  1.6117E+00  8.7191E+01  7.8869E-01  2.0844E+00  2.8858E+00  4.0056E+00  1.1641E+01  6.3840E-01  1.4352E+00
             2.9189E+00
 PARAMETER:  2.3499E-01  5.7727E-01  4.5681E+00 -1.3739E-01  8.3447E-01  1.1598E+00  1.4877E+00  2.5545E+00 -3.4880E-01  4.6130E-01
             1.1712E+00
 GRADIENT:   1.2214E+02  8.8505E+01  8.6701E-01  1.2823E+01  2.9509E+01  3.5026E+02  3.5714E+02  4.0116E+00 -1.4066E+00  2.4675E+00
             2.3845E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2472.19643706081        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     2343
 NPARAMETR:  1.1345E+00  1.4250E+00  9.0052E+01  8.8436E-01  2.0877E+00  2.7108E+00  4.0680E+00  1.1400E+01  1.0284E+00  1.4446E+00
             2.9056E+00
 PARAMETER:  2.2616E-01  4.5415E-01  4.6004E+00 -2.2889E-02  8.3605E-01  1.0973E+00  1.5032E+00  2.5336E+00  1.2800E-01  4.6782E-01
             1.1666E+00
 GRADIENT:  -3.1194E+00  8.1963E-01  8.6411E-01 -1.0225E+01  4.1672E+00 -7.6011E+00 -1.0510E+01  1.5444E+00  8.8478E-01  1.5456E+00
            -1.0038E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2472.36359052780        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2520
 NPARAMETR:  1.1390E+00  1.4067E+00  9.0232E+01  9.0076E-01  2.0827E+00  2.6974E+00  4.0921E+00  1.1385E+01  1.0596E+00  1.4524E+00
             2.9073E+00
 PARAMETER:  2.3019E-01  4.4124E-01  4.6024E+00 -4.5155E-03  8.3369E-01  1.0923E+00  1.5091E+00  2.5323E+00  1.5790E-01  4.7322E-01
             1.1672E+00
 GRADIENT:  -1.9293E+00  1.1496E+00 -7.6422E-01 -7.1327E+00 -1.3363E-02 -1.0115E+01 -1.0373E+01  4.8164E+00  6.6879E-01  4.8519E-01
             1.4011E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2472.77846897019        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     2690
 NPARAMETR:  1.1503E+00  1.3565E+00  9.1237E+01  9.3156E-01  2.0735E+00  2.6677E+00  4.2077E+00  1.1326E+01  1.0272E+00  1.4932E+00
             2.9173E+00
 PARAMETER:  2.4001E-01  4.0489E-01  4.6135E+00  2.9103E-02  8.2923E-01  1.0812E+00  1.5369E+00  2.5271E+00  1.2682E-01  5.0091E-01
             1.1707E+00
 GRADIENT:  -1.1451E-01 -3.7142E+02 -6.1089E+01  1.5079E+03 -3.4586E+02  1.0665E+02  7.8009E+01  5.7928E+01 -5.9715E+02 -5.6525E+02
            -1.2194E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2690
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2879E-03  1.4267E-02 -1.4534E-02 -5.4734E-02 -1.1501E-02
 SE:             3.2224E-02  2.8436E-02  3.9370E-03  1.2044E-02  2.5205E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6812E-01  6.1586E-01  2.2297E-04  5.5116E-06  6.4818E-01

 ETASHRINKSD(%)  1.0000E-10  4.7358E+00  8.6811E+01  5.9652E+01  1.5560E+01
 ETASHRINKVR(%)  1.0000E-10  9.2474E+00  9.8260E+01  8.3720E+01  2.8700E+01
 EBVSHRINKSD(%)  2.8008E-01  6.7200E+00  8.6797E+01  6.5943E+01  1.1212E+01
 EBVSHRINKVR(%)  5.5937E-01  1.2988E+01  9.8257E+01  8.8401E+01  2.1167E+01
 RELATIVEINF(%)  9.9434E+01  3.5196E+01  1.3292E+00  4.6741E+00  6.0481E+01
 EPSSHRINKSD(%)  1.6099E+01
 EPSSHRINKVR(%)  2.9606E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2472.7784689701857     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -818.68910920177495     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    89.87
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2472.778       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.15E+00  1.36E+00  9.12E+01  9.32E-01  2.07E+00  2.67E+00  4.21E+00  1.13E+01  1.03E+00  1.49E+00  2.92E+00
 


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
+        1.18E+02
 
 TH 2
+        4.64E-03  1.26E+04
 
 TH 3
+       -8.70E-04  4.13E-02  6.47E-02
 
 TH 4
+       -1.50E+00  3.83E+01 -3.66E-01  4.36E+05
 
 TH 5
+       -7.36E-01  7.36E+00 -3.39E-01 -7.42E+01  2.45E+03
 
 TH 6
+       -4.75E+03 -1.03E-01  1.02E-03  9.84E-01  1.31E-01  4.86E+02
 
 TH 7
+        9.94E-03  1.97E+00 -2.79E-02  5.79E+03 -6.11E+00 -4.25E-02  1.84E+02
 
 TH 8
+        1.06E-02 -6.40E-01 -6.43E-01  5.73E+00 -3.64E+01 -1.44E-02  4.24E-01  1.42E+01
 
 TH 9
+        4.82E-01 -2.48E+01  1.99E-02  7.06E+01  3.37E+00 -3.15E-01  4.40E+00 -2.90E-01  2.23E+05
 
 TH10
+       -5.89E-01  2.02E+01 -2.68E-01 -1.61E+02  2.56E+03  6.57E-01 -7.19E+02  3.93E+00  5.57E+00  1.28E+04
 
 TH11
+       -2.33E+00 -5.92E-01 -1.54E-02 -2.45E+01 -5.45E+00  1.21E+00 -1.56E+02  1.93E-01  3.34E+00  1.36E+03  7.52E+02
 
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
 #CPUT: Total CPU Time in Seconds,      106.075
Stop Time:
Wed Sep 29 08:16:54 CDT 2021
