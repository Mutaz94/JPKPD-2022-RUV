Thu Sep 30 00:36:05 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat83.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   874.606065574264        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5785E+02  1.2957E+02  2.0860E+02  9.1416E+00  9.1526E+01  3.0104E+01 -4.6150E+01 -1.4706E+02 -2.5405E+01 -1.4374E+02
            -5.5092E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1408.18257362368        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8636E-01  9.7839E-01  8.9804E-01  1.1318E+00  1.0914E+00  8.3029E-01  9.1088E-01  9.8818E-01  8.1979E-01  1.0047E+00
             5.3347E+00
 PARAMETER:  8.6262E-02  7.8151E-02 -7.5423E-03  2.2379E-01  1.8747E-01 -8.5978E-02  6.6553E-03  8.8106E-02 -9.8705E-02  1.0466E-01
             1.7742E+00
 GRADIENT:  -9.9233E+01 -2.1296E+01 -2.9872E+01 -5.0060E+00  7.2937E+00 -3.0820E+01  1.2201E+01  7.9379E+00  2.3457E+01  1.8604E+01
             3.0816E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1468.04517908950        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0054E+00  1.2838E+00  1.8260E+00  9.6754E-01  2.2825E+00  8.8791E-01  3.3968E-01  9.1539E-02  9.9981E-01  1.2498E+00
             4.2129E+00
 PARAMETER:  1.0537E-01  3.4986E-01  7.0214E-01  6.7005E-02  9.2526E-01 -1.8884E-02 -9.7976E-01 -2.2910E+00  9.9810E-02  3.2295E-01
             1.5381E+00
 GRADIENT:   3.3899E+01  6.5383E+00 -3.3964E+00  1.8785E+01  1.0601E+01 -1.1967E+01  1.8830E+00  5.3414E-03  1.5633E+01  2.5940E+00
             8.3553E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1479.30761364947        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.8343E-01  7.1135E-01  6.9474E-01  1.2160E+00  7.3051E-01  9.4079E-01  5.5464E-01  2.8563E-02  7.9050E-01  3.5589E-01
             3.7503E+00
 PARAMETER:  8.3296E-02 -2.4059E-01 -2.6421E-01  2.9554E-01 -2.1401E-01  3.8967E-02 -4.8944E-01 -3.4556E+00 -1.3509E-01 -9.3313E-01
             1.4218E+00
 GRADIENT:  -1.1934E+01  1.6984E+00  5.5262E+00 -3.1573E+00 -4.0899E+00  1.0520E+00  1.0788E+00  1.3028E-02  2.0083E-01  3.1687E+00
            -1.7815E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1481.23721627199        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0042E+00  6.8416E-01  3.7237E-01  1.1849E+00  4.7958E-01  9.5188E-01  5.2026E-01  1.0000E-02  8.6957E-01  1.8322E-01
             3.8591E+00
 PARAMETER:  1.0415E-01 -2.7957E-01 -8.8786E-01  2.6965E-01 -6.3485E-01  5.0686E-02 -5.5343E-01 -4.7737E+00 -3.9761E-02 -1.5970E+00
             1.4504E+00
 GRADIENT:  -1.2929E+01  1.5950E+01 -4.0204E-01  3.4844E+01 -5.9586E+00 -2.7754E+00  2.4100E-01  0.0000E+00  9.5739E-01  8.5788E-02
             9.7436E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1482.41244271486        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0085E+00  5.1451E-01  3.2639E-01  1.2165E+00  3.9231E-01  9.5936E-01  4.8106E-01  1.0000E-02  8.6957E-01  1.0071E-01
             3.8198E+00
 PARAMETER:  1.0848E-01 -5.6455E-01 -1.0197E+00  2.9596E-01 -8.3571E-01  5.8506E-02 -6.3176E-01 -6.2730E+00 -3.9761E-02 -2.1955E+00
             1.4402E+00
 GRADIENT:   2.3696E+00  9.0849E+00  5.5753E+00  1.8153E+01 -1.2421E+01 -1.2122E+00 -7.3308E-01  0.0000E+00 -2.2056E+00 -5.2109E-01
             2.1876E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1482.87494393827        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0055E+00  4.3384E-01  2.7113E-01  1.1895E+00  3.3087E-01  9.6739E-01  6.2806E-01  1.0000E-02  9.3836E-01  5.9016E-02
             3.7691E+00
 PARAMETER:  1.0548E-01 -7.3508E-01 -1.2052E+00  2.7356E-01 -1.0060E+00  6.6850E-02 -3.6512E-01 -8.1536E+00  3.6377E-02 -2.7300E+00
             1.4268E+00
 GRADIENT:  -2.2106E-01  8.1484E-01  4.7449E-01  1.5687E+00 -4.0579E-01 -1.7069E-01 -5.1775E-02  0.0000E+00  2.5916E-01 -3.3519E-01
            -9.9768E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1484.43748565600        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      904
 NPARAMETR:  1.0073E+00  5.3223E-01  2.8775E-01  1.1812E+00  3.6731E-01  9.6428E-01  3.9527E-01  1.0793E-02  9.4646E-01  4.1717E-01
             3.7376E+00
 PARAMETER:  1.0728E-01 -5.3067E-01 -1.1457E+00  2.6656E-01 -9.0154E-01  6.3628E-02 -8.2818E-01 -4.4289E+00  4.4970E-02 -7.7425E-01
             1.4185E+00
 GRADIENT:  -3.4101E+00 -3.1400E-01 -1.1779E+01  1.5261E+01  2.8281E+01 -2.4878E+00  9.1498E-01  7.8367E-04 -6.2282E-02  2.0249E+00
             2.3424E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1490.33126338101        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0022E+00  3.8594E-01  1.8987E-01  1.1045E+00  2.5667E-01  9.8445E-01  2.4079E-01  1.0000E-02  1.0808E+00  5.2168E-01
             3.4445E+00
 PARAMETER:  1.0221E-01 -8.5207E-01 -1.5614E+00  1.9940E-01 -1.2600E+00  8.4327E-02 -1.3238E+00 -5.5920E+00  1.7771E-01 -5.5070E-01
             1.3368E+00
 GRADIENT:  -7.6075E+00  1.9145E+00 -7.2392E+00  5.0202E+00  3.9887E+00 -2.2946E+00  5.3598E-01  0.0000E+00 -7.3787E+00  2.1838E+00
            -4.9438E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1490.78615774982        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1257
 NPARAMETR:  1.0058E+00  3.7215E-01  1.8725E-01  1.0969E+00  2.5162E-01  9.9028E-01  9.2511E-02  1.0000E-02  1.1229E+00  5.0641E-01
             3.4597E+00
 PARAMETER:  1.0576E-01 -8.8846E-01 -1.5753E+00  1.9250E-01 -1.2798E+00  9.0235E-02 -2.2804E+00 -4.9312E+00  2.1596E-01 -5.8040E-01
             1.3412E+00
 GRADIENT:   6.3896E-01 -5.9267E-01 -7.9398E-01 -2.5172E+00  8.3369E-01  2.0370E-02  7.1935E-02  0.0000E+00 -4.1978E-01 -4.0672E-02
             6.9586E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1491.23218497329        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1436
 NPARAMETR:  1.0058E+00  3.7643E-01  1.8892E-01  1.0999E+00  2.5342E-01  9.9093E-01  1.0000E-02  5.3035E-01  1.1200E+00  5.0763E-01
             3.4559E+00
 PARAMETER:  1.0580E-01 -8.7702E-01 -1.5664E+00  1.9519E-01 -1.2727E+00  9.0890E-02 -6.6888E+00 -5.3422E-01  2.1336E-01 -5.7801E-01
             1.3401E+00
 GRADIENT:  -2.3318E+00 -3.1673E+00  3.3028E+00 -4.4097E-01  5.4401E+00 -3.6232E-01  0.0000E+00  4.3862E-01  2.1512E-01  6.9933E+00
             1.6707E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1501.72177941706        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1619
 NPARAMETR:  1.0016E+00  3.4852E-01  1.4078E-01  1.0391E+00  2.1455E-01  1.0147E+00  1.0000E-02  1.4423E+00  1.3798E+00  2.8453E-01
             3.0676E+00
 PARAMETER:  1.0162E-01 -9.5405E-01 -1.8605E+00  1.3838E-01 -1.4392E+00  1.1455E-01 -9.6322E+00  4.6623E-01  4.2197E-01 -1.1569E+00
             1.2209E+00
 GRADIENT:   3.6918E+00  3.5714E+00  1.7251E+01  1.4875E+00 -9.1076E+00  1.7389E+00  0.0000E+00 -6.4455E+00  3.5148E-01  2.4566E+00
             7.1813E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1503.42878668868        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1804
 NPARAMETR:  9.9684E-01  3.3631E-01  1.2998E-01  1.0286E+00  2.0647E-01  1.0051E+00  1.0000E-02  1.5628E+00  1.4553E+00  2.2956E-01
             2.9819E+00
 PARAMETER:  9.6837E-02 -9.8971E-01 -1.9404E+00  1.2818E-01 -1.4776E+00  1.0504E-01 -1.0209E+01  5.4648E-01  4.7524E-01 -1.3716E+00
             1.1926E+00
 GRADIENT:   8.3491E-01  6.7933E+00  3.2568E+00  4.9345E-01 -1.7394E+01 -1.0515E+00  0.0000E+00  1.4786E+00 -1.5808E+00 -1.6708E-01
            -1.1094E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1503.53056110020        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1979
 NPARAMETR:  9.9686E-01  3.3446E-01  1.3029E-01  1.0310E+00  2.0709E-01  1.0069E+00  1.0000E-02  1.5587E+00  1.4572E+00  2.2113E-01
             3.0153E+00
 PARAMETER:  9.6854E-02 -9.9523E-01 -1.9380E+00  1.3055E-01 -1.4746E+00  1.0685E-01 -1.0209E+01  5.4382E-01  4.7649E-01 -1.4090E+00
             1.2037E+00
 GRADIENT:  -3.8600E-02 -4.9671E-01 -1.0706E-01 -8.7903E-05 -9.7700E-01 -4.5943E-02  0.0000E+00  8.3536E-01 -6.5433E-02  1.3631E-02
             1.1019E-01

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1503.53056110020        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2001
 NPARAMETR:  9.9686E-01  3.3446E-01  1.3029E-01  1.0310E+00  2.0709E-01  1.0069E+00  1.0000E-02  1.5587E+00  1.4572E+00  2.2113E-01
             3.0153E+00
 PARAMETER:  9.6854E-02 -9.9523E-01 -1.9380E+00  1.3055E-01 -1.4746E+00  1.0685E-01 -1.0209E+01  5.4382E-01  4.7649E-01 -1.4090E+00
             1.2037E+00
 GRADIENT:  -3.8600E-02 -4.9671E-01 -1.0706E-01 -8.7903E-05 -9.7700E-01 -4.5943E-02  0.0000E+00  8.3536E-01 -6.5433E-02  1.3631E-02
             1.1019E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2001
 NO. OF SIG. DIGITS IN FINAL EST.:  2.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1291E-03 -3.5175E-04  1.5685E-02 -8.5225E-03  1.9567E-02
 SE:             2.8994E-02  1.5095E-04  2.4116E-02  2.6213E-02  8.6623E-03
 N:                     100         100         100         100         100

 P VAL.:         9.1406E-01  1.9796E-02  5.1545E-01  7.4509E-01  2.3891E-02

 ETASHRINKSD(%)  2.8668E+00  9.9494E+01  1.9208E+01  1.2182E+01  7.0980E+01
 ETASHRINKVR(%)  5.6514E+00  9.9997E+01  3.4727E+01  2.2879E+01  9.1579E+01
 EBVSHRINKSD(%)  2.8732E+00  9.9485E+01  1.8122E+01  1.0323E+01  7.2592E+01
 EBVSHRINKVR(%)  5.6638E+00  9.9997E+01  3.2961E+01  1.9580E+01  9.2488E+01
 RELATIVEINF(%)  9.4183E+01  6.4982E-04  1.8299E+01  6.2743E+01  8.4280E-01
 EPSSHRINKSD(%)  2.6775E+01
 EPSSHRINKVR(%)  4.6381E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1503.5305611002038     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -584.59202789553115     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1503.531       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.97E-01  3.34E-01  1.30E-01  1.03E+00  2.07E-01  1.01E+00  1.00E-02  1.56E+00  1.46E+00  2.21E-01  3.02E+00
 


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
+        1.05E+03
 
 TH 2
+       -3.78E+01  8.29E+03
 
 TH 3
+       -1.57E+02  1.59E+03  3.03E+04
 
 TH 4
+       -1.34E+01  1.79E+02 -4.67E+02  4.33E+02
 
 TH 5
+        1.81E+02 -1.01E+04 -1.63E+04 -1.02E+03  5.32E+04
 
 TH 6
+        6.45E+00 -3.32E+01  4.30E-01 -1.74E+01  6.33E+01  1.76E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -5.48E+01  2.25E+03 -3.44E+03  4.79E+03 -5.06E+02 -7.00E+01  0.00E+00  9.65E+02
 
 TH 9
+        1.02E+01 -7.03E+01  2.62E+02 -8.44E+00  5.02E+02 -8.75E-02  0.00E+00 -1.15E+02  5.55E+01
 
 TH10
+       -4.48E-01 -1.37E+02  1.97E+02 -1.27E+00  7.13E+02  1.61E+00  0.00E+00 -2.06E+03  1.38E+01  6.56E+01
 
 TH11
+       -1.93E+01 -1.77E+02  8.04E+02 -3.41E+01  3.65E+02  2.15E+00  0.00E+00 -2.23E+02  1.90E+01  2.67E+01  1.06E+02
 
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
 #CPUT: Total CPU Time in Seconds,       43.852
Stop Time:
Thu Sep 30 00:36:50 CDT 2021
