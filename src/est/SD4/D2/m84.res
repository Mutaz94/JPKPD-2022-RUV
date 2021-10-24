Sun Oct 24 04:43:28 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat84.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1619.65077589098        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9132E+02  5.8850E+00 -3.0090E+01  5.7977E+01  2.0360E+01  3.5458E+01 -3.8216E+01  8.6261E+00 -3.3800E+01 -1.6127E+00
            -2.6406E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1633.58304844715        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0010E+00  1.1828E+00  1.2871E+00  9.0141E-01  1.2042E+00  9.8504E-01  1.5773E+00  8.9014E-01  1.2367E+00  1.0275E+00
             1.0609E+00
 PARAMETER:  1.0100E-01  2.6790E-01  3.5241E-01 -3.7933E-03  2.8583E-01  8.4932E-02  5.5569E-01 -1.6372E-02  3.1247E-01  1.2708E-01
             1.5914E-01
 GRADIENT:  -3.5057E+01  6.7936E+00  8.3446E+00 -6.9891E+00  9.4221E+00 -1.3780E+00  2.1432E+01  5.8159E-01  2.0589E+01 -1.8150E+01
            -6.8580E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1635.20747887486        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.9533E-01  1.1325E+00  1.5018E+00  9.5168E-01  1.2865E+00  9.7516E-01  1.5591E+00  7.4452E-01  1.1595E+00  1.2211E+00
             1.0887E+00
 PARAMETER:  9.5317E-02  2.2446E-01  5.0667E-01  5.0476E-02  3.5195E-01  7.4850E-02  5.4408E-01 -1.9502E-01  2.4800E-01  2.9974E-01
             1.8500E-01
 GRADIENT:  -4.7772E+01  8.9074E+00  4.8232E+00  3.9847E+00  1.9971E+01 -5.2517E+00  1.5149E+01 -7.0021E-01  1.3775E+01 -6.1515E+00
             2.9047E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1640.45386136635        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0177E+00  8.9657E-01  1.0987E+00  1.0811E+00  1.0156E+00  9.8800E-01  1.8003E+00  2.6261E-01  8.7668E-01  1.0482E+00
             1.0606E+00
 PARAMETER:  1.1751E-01 -9.1822E-03  1.9410E-01  1.7797E-01  1.1547E-01  8.7930E-02  6.8798E-01 -1.2371E+00 -3.1610E-02  1.4706E-01
             1.5880E-01
 GRADIENT:   3.7041E+00  8.9330E+00 -4.0725E-01  1.4990E+01  1.9527E+00  4.4016E-01 -1.0942E+00  7.7370E-02 -3.8764E-01 -1.3901E+00
            -1.3120E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1641.19392704482        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0133E+00  5.6832E-01  1.2955E+00  1.2803E+00  9.7703E-01  9.8404E-01  2.3945E+00  1.3862E-01  8.2810E-01  1.1139E+00
             1.0687E+00
 PARAMETER:  1.1325E-01 -4.6507E-01  3.5890E-01  3.4706E-01  7.6766E-02  8.3910E-02  9.7318E-01 -1.8761E+00 -8.8618E-02  2.0790E-01
             1.6647E-01
 GRADIENT:   2.9601E+00  3.3058E+00  2.3454E+00  5.4074E+00 -2.7042E+00  2.4795E-01  9.2489E-02 -7.0149E-03  1.2863E+00  3.8672E-01
            -1.4795E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1641.26695252590        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  1.0124E+00  5.4087E-01  1.2923E+00  1.2879E+00  9.6807E-01  9.8264E-01  2.4680E+00  1.3974E-01  8.2128E-01  1.1100E+00
             1.0681E+00
 PARAMETER:  1.1231E-01 -5.1457E-01  3.5645E-01  3.5305E-01  6.7546E-02  8.2483E-02  1.0034E+00 -1.8680E+00 -9.6887E-02  2.0435E-01
             1.6587E-01
 GRADIENT:   4.2108E+02  5.4100E+01  6.7912E+00  3.8231E+02  3.4378E+00  2.5114E+01  6.0356E+01  3.4388E-02  8.5085E+00  1.9988E+00
             1.2737E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1641.28057361792        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1051
 NPARAMETR:  1.0122E+00  5.4226E-01  1.2855E+00  1.2895E+00  9.6814E-01  9.8333E-01  2.4803E+00  1.7602E-01  8.1099E-01  1.1082E+00
             1.0682E+00
 PARAMETER:  1.1214E-01 -5.1202E-01  3.5114E-01  3.5427E-01  6.7624E-02  8.3186E-02  1.0084E+00 -1.6372E+00 -1.0950E-01  2.0272E-01
             1.6595E-01
 GRADIENT:   9.4920E-01  1.1475E+00 -3.9063E-01 -6.6813E-01  8.2166E-01 -1.0211E-03  6.6853E-01  7.7679E-03 -7.2657E-02  3.3582E-01
             2.7541E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1641.28775288298        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1234            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0125E+00  5.3747E-01  1.2839E+00  1.2894E+00  9.6655E-01  9.8337E-01  2.4883E+00  1.4459E-01  8.1230E-01  1.1070E+00
             1.0681E+00
 PARAMETER:  1.1243E-01 -5.2088E-01  3.4987E-01  3.5419E-01  6.5979E-02  8.3230E-02  1.0116E+00 -1.8339E+00 -1.0789E-01  2.0162E-01
             1.6587E-01
 GRADIENT:   4.2143E+02  5.3994E+01  4.3258E+00  3.8611E+02  7.0301E+00  2.5355E+01  6.1819E+01  4.3576E-02  7.8932E+00  1.7513E+00
             1.5425E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1641.28993568666        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  1.0123E+00  5.3678E-01  1.2834E+00  1.2909E+00  9.6543E-01  9.8334E-01  2.4890E+00  1.2350E-01  8.1111E-01  1.1065E+00
             1.0678E+00
 PARAMETER:  1.1223E-01 -5.2217E-01  3.4954E-01  3.5537E-01  6.4817E-02  8.3201E-02  1.0119E+00 -1.9915E+00 -1.0935E-01  2.0120E-01
             1.6564E-01
 GRADIENT:   1.3132E+00  4.5251E-01  3.1888E-01 -3.2604E+00  4.0738E-01  3.1619E-02  4.8480E-01  3.4462E-04 -3.8667E-02 -3.8007E-03
            -4.6421E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1641.29290084359        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1592
 NPARAMETR:  1.0126E+00  5.3436E-01  1.2821E+00  1.2914E+00  9.6353E-01  9.8337E-01  2.5034E+00  1.0520E-01  8.1085E-01  1.1063E+00
             1.0678E+00
 PARAMETER:  1.1256E-01 -5.2668E-01  3.4847E-01  3.5576E-01  6.2852E-02  8.3234E-02  1.0176E+00 -2.1519E+00 -1.0968E-01  2.0099E-01
             1.6557E-01
 GRADIENT:   2.1530E+00  4.6858E-01  9.0642E-01 -4.9300E+00 -5.0976E-01  4.8781E-02  1.0023E+00 -1.6296E-04  1.3767E-01  7.2590E-02
            -5.2703E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1641.29565178569        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1776
 NPARAMETR:  1.0125E+00  5.3232E-01  1.2805E+00  1.2925E+00  9.6241E-01  9.8335E-01  2.5061E+00  8.4381E-02  8.1034E-01  1.1056E+00
             1.0677E+00
 PARAMETER:  1.1238E-01 -5.3052E-01  3.4726E-01  3.5657E-01  6.1683E-02  8.3215E-02  1.0187E+00 -2.3724E+00 -1.1031E-01  2.0042E-01
             1.6549E-01
 GRADIENT:   1.7588E+00  3.3511E-01  7.7052E-01 -4.8535E+00 -2.8529E-01  4.7624E-02  8.3383E-01 -1.0262E-04  4.5929E-02 -7.6724E-03
            -1.0282E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1641.29625078322        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1954
 NPARAMETR:  1.0122E+00  5.3171E-01  1.2798E+00  1.2935E+00  9.6232E-01  9.8331E-01  2.5024E+00  5.4489E-02  8.1039E-01  1.1057E+00
             1.0678E+00
 PARAMETER:  1.1214E-01 -5.3166E-01  3.4671E-01  3.5738E-01  6.1588E-02  8.3171E-02  1.0172E+00 -2.8098E+00 -1.1024E-01  2.0052E-01
             1.6559E-01
 GRADIENT:   1.1970E+00  3.7027E-01  3.6266E-01 -3.4390E+00  3.7035E-01  3.5263E-02  5.0015E-01 -3.6410E-05 -6.7454E-02 -9.0616E-02
            -1.0303E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1641.29630660435        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2131
 NPARAMETR:  1.0117E+00  5.2921E-01  1.2795E+00  1.2953E+00  9.6136E-01  9.8314E-01  2.5013E+00  3.3678E-02  8.1103E-01  1.1060E+00
             1.0679E+00
 PARAMETER:  1.1159E-01 -5.3637E-01  3.4644E-01  3.5874E-01  6.0597E-02  8.2995E-02  1.0168E+00 -3.2909E+00 -1.0945E-01  2.0074E-01
             1.6569E-01
 GRADIENT:  -9.5804E-03  2.3157E-01  1.2451E-01 -2.8424E+00  5.1964E-01 -1.9177E-02  8.8945E-02  7.5313E-05 -4.7656E-02 -3.8635E-02
            -4.0839E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1641.29990683730        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     2295
 NPARAMETR:  1.0123E+00  5.2797E-01  1.2781E+00  1.2951E+00  9.6082E-01  9.8327E-01  2.5140E+00  7.0200E-02  8.0870E-01  1.1057E+00
             1.0678E+00
 PARAMETER:  1.1223E-01 -5.3872E-01  3.4539E-01  3.5861E-01  6.0037E-02  8.3124E-02  1.0219E+00 -2.5564E+00 -1.1232E-01  2.0052E-01
             1.6559E-01
 GRADIENT:   1.4932E+00  1.3348E-01 -1.3568E-01 -3.9739E+00  7.7296E-01  2.2282E-02  5.3433E-01  6.9887E-04 -1.6193E-01  6.6440E-02
             4.9451E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1641.30159175315        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2478
 NPARAMETR:  1.0124E+00  5.2651E-01  1.2776E+00  1.2958E+00  9.5990E-01  9.8329E-01  2.5198E+00  4.6417E-02  8.0994E-01  1.1055E+00
             1.0678E+00
 PARAMETER:  1.1232E-01 -5.4149E-01  3.4500E-01  3.5910E-01  5.9075E-02  8.3154E-02  1.0242E+00 -2.9701E+00 -1.1079E-01  2.0025E-01
             1.6563E-01
 GRADIENT:   1.7334E+00  1.3194E-01  4.5566E-02 -4.5203E+00  5.0212E-01  3.9993E-02  7.4082E-01  3.1343E-04  1.2739E-01  6.0400E-02
             7.4414E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1641.30263818573        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2655
 NPARAMETR:  1.0124E+00  5.2553E-01  1.2777E+00  1.2966E+00  9.5927E-01  9.8328E-01  2.5234E+00  2.1406E-02  8.0999E-01  1.1053E+00
             1.0678E+00
 PARAMETER:  1.1230E-01 -5.4335E-01  3.4505E-01  3.5973E-01  5.8412E-02  8.3134E-02  1.0256E+00 -3.7441E+00 -1.1073E-01  2.0015E-01
             1.6559E-01
 GRADIENT:   1.7340E+00  2.8478E-01  2.6756E-01 -4.2910E+00  1.4290E-01  3.6983E-02  8.0212E-01  6.9671E-05  1.3910E-01  3.7886E-02
             1.5094E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1641.30273095137        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2830
 NPARAMETR:  1.0120E+00  5.2453E-01  1.2778E+00  1.2977E+00  9.5898E-01  9.8320E-01  2.5202E+00  1.5694E-02  8.1014E-01  1.1053E+00
             1.0678E+00
 PARAMETER:  1.1196E-01 -5.4525E-01  3.4514E-01  3.6056E-01  5.8115E-02  8.3057E-02  1.0243E+00 -4.0544E+00 -1.1055E-01  2.0013E-01
             1.6561E-01
 GRADIENT:   9.6284E-01  2.7332E-01  8.5449E-02 -3.3936E+00  3.0246E-01  1.5083E-02  4.4881E-01  4.8200E-05  5.1776E-02  1.7372E-02
             2.5446E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1641.30524632766        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3012
 NPARAMETR:  1.0126E+00  5.2173E-01  1.2764E+00  1.2982E+00  9.5771E-01  9.8330E-01  2.5388E+00  1.0067E-02  8.0749E-01  1.1047E+00
             1.0676E+00
 PARAMETER:  1.1249E-01 -5.5060E-01  3.4405E-01  3.6098E-01  5.6786E-02  8.3159E-02  1.0317E+00 -4.4985E+00 -1.1383E-01  1.9957E-01
             1.6539E-01
 GRADIENT:   2.2647E+00  1.9465E-01  1.5512E-01 -4.9398E+00  2.5036E-01  5.2663E-02  9.9427E-01  2.5229E-04 -8.6519E-02  1.3897E-02
            -3.3408E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1641.30624129703        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     3181
 NPARAMETR:  1.0123E+00  5.2048E-01  1.2763E+00  1.2995E+00  9.5715E-01  9.8324E-01  2.5499E+00  1.0000E-02  8.0825E-01  1.1044E+00
             1.0676E+00
 PARAMETER:  1.1225E-01 -5.5226E-01  3.4421E-01  3.6165E-01  5.6287E-02  8.3109E-02  1.0312E+00 -4.5602E+00 -1.1176E-01  1.9958E-01
             1.6550E-01
 GRADIENT:   3.0279E-02  1.3211E-01  1.4656E-01 -7.9458E-01  1.0886E-01  2.9307E-03 -1.8498E+03  0.0000E+00  8.9852E-02  5.1858E-02
             3.0602E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3181
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3649E-04  2.7993E-02 -2.9058E-04 -3.2389E-02 -1.1579E-02
 SE:             2.9833E-02  2.0977E-02  1.4003E-04  2.2283E-02  2.3151E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9635E-01  1.8206E-01  3.7976E-02  1.4607E-01  6.1697E-01

 ETASHRINKSD(%)  5.6149E-02  2.9724E+01  9.9531E+01  2.5349E+01  2.2441E+01
 ETASHRINKVR(%)  1.1227E-01  5.0612E+01  9.9998E+01  4.4272E+01  3.9847E+01
 EBVSHRINKSD(%)  4.9444E-01  3.2541E+01  9.9554E+01  2.2375E+01  1.8201E+01
 EBVSHRINKVR(%)  9.8643E-01  5.4493E+01  9.9998E+01  3.9744E+01  3.3090E+01
 RELATIVEINF(%)  9.8514E+01  6.9403E+00  2.8862E-04  9.3636E+00  1.0168E+01
 EPSSHRINKSD(%)  4.1467E+01
 EPSSHRINKVR(%)  6.5739E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1641.3062412970301     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -906.15541473329188     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1641.306       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  5.21E-01  1.28E+00  1.30E+00  9.57E-01  9.83E-01  2.54E+00  1.00E-02  8.09E-01  1.10E+00  1.07E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,      111.173
Stop Time:
Sun Oct 24 04:43:48 CDT 2021
