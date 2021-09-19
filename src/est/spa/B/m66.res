Sat Sep 18 08:39:11 CDT 2021
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
$DATA ../../../../data/spa/B/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1640.17380987063        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4094E+02 -7.7780E+01 -4.6661E+00 -9.5374E+01 -2.7417E+00  1.1901E+01 -8.0841E+00  8.3924E+00  2.0462E+01 -3.0177E+00
            -7.7278E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1650.47611957212        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.5805E-01  1.1667E+00  1.0411E+00  9.5404E-01  1.1168E+00  9.4543E-01  1.1448E+00  8.9521E-01  7.5082E-01  1.0914E+00
             1.0469E+00
 PARAMETER:  5.7140E-02  2.5420E-01  1.4023E-01  5.2953E-02  2.1043E-01  4.3882E-02  2.3526E-01 -1.0692E-02 -1.8659E-01  1.8743E-01
             1.4583E-01
 GRADIENT:   4.1898E+01 -1.9902E+00  1.0980E+01 -2.6199E+01  2.5723E+00 -4.9220E+00 -2.6831E+00 -1.8970E-01 -8.0357E+00 -6.0282E+00
             4.4674E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1652.22952444093        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.6484E-01  8.7474E-01  1.0229E+00  1.1624E+00  9.5238E-01  9.2261E-01  1.5598E+00  5.0357E-01  6.6788E-01  1.0291E+00
             1.0431E+00
 PARAMETER:  6.4206E-02 -3.3832E-02  1.2261E-01  2.5050E-01  5.1205E-02  1.9454E-02  5.4457E-01 -5.8604E-01 -3.0365E-01  1.2864E-01
             1.4219E-01
 GRADIENT:   6.3525E+01  2.8659E+01  1.6465E+01  4.2496E+01 -1.0048E+01 -1.5125E+01  6.5784E+00 -8.3897E-01 -5.6071E+00 -3.3389E+00
             3.4870E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1653.91784722474        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.4631E-01  9.2852E-01  8.3233E-01  1.0978E+00  8.7855E-01  9.4753E-01  1.3912E+00  3.1086E-01  7.3368E-01  9.3768E-01
             1.0235E+00
 PARAMETER:  4.4817E-02  2.5840E-02 -8.3525E-02  1.9332E-01 -2.9479E-02  4.6100E-02  4.3015E-01 -1.0684E+00 -2.0969E-01  3.5650E-02
             1.2321E-01
 GRADIENT:   1.0510E+01  3.5811E+00 -4.8010E+00  1.1659E+01  5.8978E+00 -4.4226E+00  2.1093E+00  6.2368E-01  6.8926E-01  6.4608E-01
             1.4305E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1653.92352639592        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.4447E-01  9.2842E-01  8.2006E-01  1.0946E+00  8.7030E-01  9.5189E-01  1.3821E+00  2.8212E-01  7.3588E-01  9.2861E-01
             1.0221E+00
 PARAMETER:  4.2865E-02  2.5730E-02 -9.8382E-02  1.9037E-01 -3.8913E-02  5.0689E-02  4.2358E-01 -1.1654E+00 -2.0669E-01  2.5931E-02
             1.2190E-01
 GRADIENT:   5.7928E+00  1.6635E+00 -3.8034E+00  7.1400E+00  4.1702E+00 -2.6778E+00  1.3311E+00  5.3506E-01  7.2359E-01  6.7580E-01
             9.7769E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1653.92816108305        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.4299E-01  9.2676E-01  8.0821E-01  1.0925E+00  8.6160E-01  9.5565E-01  1.3766E+00  2.4210E-01  7.3673E-01  9.1963E-01
             1.0211E+00
 PARAMETER:  4.1301E-02  2.3940E-02 -1.1293E-01  1.8851E-01 -4.8958E-02  5.4637E-02  4.1965E-01 -1.3184E+00 -2.0553E-01  1.6212E-02
             1.2089E-01
 GRADIENT:   2.0523E+00  1.6722E-01 -2.5692E+00  3.1209E+00  2.4197E+00 -1.2007E+00  6.6646E-01  4.0232E-01  6.2042E-01  5.9384E-01
             5.6331E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1653.93313121648        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  9.4185E-01  9.2444E-01  7.9737E-01  1.0913E+00  8.5329E-01  9.5878E-01  1.3735E+00  1.9133E-01  7.3686E-01  9.1152E-01
             1.0203E+00
 PARAMETER:  4.0087E-02  2.1434E-02 -1.2643E-01  1.8741E-01 -5.8650E-02  5.7901E-02  4.1733E-01 -1.5537E+00 -2.0535E-01  7.3546E-03
             1.2012E-01
 GRADIENT:  -8.3404E-01 -9.2026E-01 -1.3538E+00 -1.3828E-01  8.3389E-01  3.6997E-03  1.0494E-01  2.5418E-01  4.5542E-01  4.6197E-01
             1.9728E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1653.94896907730        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      514
 NPARAMETR:  9.4100E-01  9.2157E-01  7.8747E-01  1.0909E+00  8.4540E-01  9.6132E-01  1.3721E+00  1.1379E-01  7.3644E-01  9.0446E-01
             1.0198E+00
 PARAMETER:  3.9184E-02  1.8319E-02 -1.3893E-01  1.8697E-01 -6.7948E-02  6.0551E-02  4.1635E-01 -2.0734E+00 -2.0593E-01 -4.1263E-04
             1.1957E-01
 GRADIENT:  -2.9814E+00 -1.6559E+00 -1.3745E-01 -2.7525E+00 -6.1653E-01  9.7171E-01 -3.6859E-01  9.0632E-02  2.3583E-01  2.8731E-01
            -1.2659E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1654.49364405559        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  9.5334E-01  8.9798E-01  8.1469E-01  1.1157E+00  8.5291E-01  9.6878E-01  1.4116E+00  1.0000E-02  7.3398E-01  9.2442E-01
             1.0247E+00
 PARAMETER:  5.2213E-02 -7.6102E-03 -1.0495E-01  2.0948E-01 -5.9101E-02  6.8278E-02  4.4471E-01 -7.5122E+00 -2.0927E-01  2.1411E-02
             1.2441E-01
 GRADIENT:  -6.0589E+00  2.5973E+00 -2.1186E+00  2.1511E+00  2.7455E+00  8.6281E-01 -1.2103E+00  0.0000E+00  5.7570E-01 -2.5619E-01
             1.2937E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1655.28410733292        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  9.5363E-01  6.6040E-01  8.4421E-01  1.2516E+00  7.7700E-01  9.6198E-01  1.8353E+00  4.7947E-02  6.6785E-01  9.0337E-01
             1.0216E+00
 PARAMETER:  5.2526E-02 -3.1490E-01 -6.9350E-02  3.2445E-01 -1.5232E-01  6.1235E-02  7.0720E-01 -2.9377E+00 -3.0370E-01 -1.6237E-03
             1.2138E-01
 GRADIENT:  -1.3122E-01  2.0454E+00  5.4527E-01  2.1648E+00 -7.7654E-01 -1.0729E+00  2.8823E-01  1.0707E-02 -3.4097E-01 -4.8735E-02
             4.2809E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1655.32841707097        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  9.5294E-01  6.0596E-01  8.4586E-01  1.2799E+00  7.6051E-01  9.6376E-01  1.9529E+00  2.2349E-02  6.5984E-01  8.9914E-01
             1.0213E+00
 PARAMETER:  5.1800E-02 -4.0094E-01 -6.7406E-02  3.4680E-01 -1.7377E-01  6.3088E-02  7.6930E-01 -3.7010E+00 -3.1576E-01 -6.3210E-03
             1.2110E-01
 GRADIENT:  -3.4628E-01 -1.1103E-01 -3.5407E-01  1.1238E-01  4.1928E-01 -1.0431E-01  1.7476E-02  2.3973E-03  9.4691E-02 -6.1101E-04
             6.5667E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1655.32954216766        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1153
 NPARAMETR:  9.5310E-01  6.0770E-01  8.4611E-01  1.2790E+00  7.6108E-01  9.6404E-01  1.9487E+00  1.0000E-02  6.5980E-01  8.9969E-01
             1.0212E+00
 PARAMETER:  5.1959E-02 -3.9807E-01 -6.7102E-02  3.4611E-01 -1.7302E-01  6.3377E-02  7.6718E-01 -4.6569E+00 -3.1582E-01 -5.7103E-03
             1.2102E-01
 GRADIENT:  -7.2674E-03  9.6965E-03  8.0153E-03  1.9320E-02 -1.5570E-02  3.0865E-03 -1.3039E-03  0.0000E+00 -3.4455E-03 -4.7853E-05
            -1.1188E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1153
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3762E-04  2.2647E-02 -4.8134E-04 -2.5838E-02 -1.1848E-03
 SE:             2.9826E-02  2.1393E-02  2.1184E-04  2.2833E-02  2.3493E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8829E-01  2.8977E-01  2.3075E-02  2.5779E-01  9.5978E-01

 ETASHRINKSD(%)  7.9339E-02  2.8331E+01  9.9290E+01  2.3507E+01  2.1296E+01
 ETASHRINKVR(%)  1.5862E-01  4.8635E+01  9.9995E+01  4.1488E+01  3.8057E+01
 EBVSHRINKSD(%)  4.8159E-01  2.9322E+01  9.9327E+01  2.2247E+01  1.8835E+01
 EBVSHRINKVR(%)  9.6086E-01  5.0046E+01  9.9995E+01  3.9544E+01  3.4123E+01
 RELATIVEINF(%)  9.8415E+01  6.6550E+00  4.4302E-04  8.6973E+00  5.7646E+00
 EPSSHRINKSD(%)  4.2650E+01
 EPSSHRINKVR(%)  6.7110E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1655.3295421676598     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -920.17871560392166     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1655.330       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  6.08E-01  8.46E-01  1.28E+00  7.61E-01  9.64E-01  1.95E+00  1.00E-02  6.60E-01  9.00E-01  1.02E+00
 


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
+        1.31E+03
 
 TH 2
+       -1.31E+01  4.19E+02
 
 TH 3
+        2.14E+01  1.69E+02  7.70E+02
 
 TH 4
+       -7.87E+00  3.84E+02 -3.11E+02  9.70E+02
 
 TH 5
+       -2.94E+00 -3.36E+02 -9.94E+02  3.34E+02  1.56E+03
 
 TH 6
+       -4.41E+00 -3.10E+00  6.60E+00 -3.53E+00 -1.01E+00  2.15E+02
 
 TH 7
+        1.75E+00  3.56E+01 -1.55E+00 -1.32E+01  7.50E-02  9.09E-02  1.91E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.18E+00 -2.00E+01 -4.12E+01 -2.13E+01  4.16E+01  1.69E+00  1.62E+01  0.00E+00  1.81E+02
 
 TH10
+       -2.52E-02 -4.74E+00 -6.41E+01 -2.76E+01 -5.07E+01 -3.39E+00  3.32E+00  0.00E+00  1.48E+01  1.08E+02
 
 TH11
+       -8.08E+00 -1.09E+01 -4.15E+01 -6.51E+00  1.88E+01  1.17E+00  1.76E+00  0.00E+00  1.60E+01  2.75E+01  2.14E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.065
Stop Time:
Sat Sep 18 08:39:31 CDT 2021
