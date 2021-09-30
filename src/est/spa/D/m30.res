Wed Sep 29 19:56:16 CDT 2021
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
$DATA ../../../../data/spa/D/dat30.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   10408.8121057926        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.8571E+02  7.7190E+01 -4.1883E+01 -3.5889E+01  1.7740E+02 -1.0630E+03 -5.4348E+02 -1.1320E+02 -9.9625E+02 -2.8984E+02
            -2.0929E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -720.919594067320        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4971E+00  1.2921E+00  1.0434E+00  1.7901E+00  1.1424E+00  2.1788E+00  1.3535E+00  9.7407E-01  1.5821E+00  1.0482E+00
             1.3837E+01
 PARAMETER:  5.0350E-01  3.5628E-01  1.4249E-01  6.8224E-01  2.3316E-01  8.7877E-01  4.0271E-01  7.3724E-02  5.5875E-01  1.4712E-01
             2.7273E+00
 GRADIENT:   1.3263E+01  3.3733E+01 -3.5201E+00  5.8777E+01 -1.6440E+01  6.3481E+01 -3.9822E+00  4.3517E+00  4.5634E+00  3.3614E+00
             1.6584E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -740.688645347236        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4247E+00  9.3490E-01  1.4695E+00  1.9455E+00  1.7069E+00  1.8203E+00  3.7668E+00  4.6793E-01  1.4044E+00  3.7425E+00
             1.2523E+01
 PARAMETER:  4.5394E-01  3.2684E-02  4.8495E-01  7.6553E-01  6.3469E-01  6.9902E-01  1.4262E+00 -6.5944E-01  4.3962E-01  1.4198E+00
             2.6276E+00
 GRADIENT:   1.9528E+01  2.3428E+01 -9.0171E+00  5.0834E+01 -3.7353E+01  1.0336E+01  1.2248E+01  4.9516E-01  1.8027E+01  3.3251E+01
             1.3740E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -777.483377576390        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1469E+00  9.5827E-01  1.7007E+00  1.2817E+00  5.4763E+00  1.6271E+00  1.9287E+00  3.4494E-01  1.1753E+00  3.1909E+00
             9.2065E+00
 PARAMETER:  2.3710E-01  5.7379E-02  6.3106E-01  3.4817E-01  1.8004E+00  5.8681E-01  7.5687E-01 -9.6440E-01  2.6154E-01  1.2603E+00
             2.3199E+00
 GRADIENT:  -1.6642E+01  1.5201E+00  9.7233E+00 -1.3735E+01 -3.2593E+00  1.1962E+01  1.7080E+00 -7.7305E-03 -2.4397E+00 -4.6154E-01
            -1.2514E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -804.176958350889        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0428E+00  2.0328E-01  7.0350E-01  1.6322E+00  1.2260E+01  1.7168E+00  2.2070E+00  1.0000E-02  1.1780E+00  1.5830E+00
             8.6480E+00
 PARAMETER:  1.4187E-01 -1.4932E+00 -2.5168E-01  5.8993E-01  2.6063E+00  6.4045E-01  8.9164E-01 -5.3593E+00  2.6381E-01  5.5931E-01
             2.2573E+00
 GRADIENT:  -4.7103E+01  7.1418E+00  2.6623E+01  6.2140E+01 -2.1281E-01 -3.2825E+00 -1.9637E-02  0.0000E+00 -1.3512E+01 -1.4431E-02
            -1.4168E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -845.150660378525        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  8.4804E-01  7.5237E-02  9.6139E-02  8.4075E-01  1.4027E+02  1.7960E+00  4.2194E-01  1.0000E-02  1.2999E+00  1.2348E+00
             7.7868E+00
 PARAMETER: -6.4822E-02 -2.4871E+00 -2.2420E+00 -7.3456E-02  5.0436E+00  6.8558E-01 -7.6290E-01 -1.3604E+01  3.6232E-01  3.1093E-01
             2.1524E+00
 GRADIENT:   6.8098E+01  2.2308E+01 -6.6214E+01  1.1516E+02 -6.6854E-02 -8.6151E+00  1.6300E+00  0.0000E+00  1.7930E+01  6.2341E-04
            -4.1942E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -863.526058312912        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      507
 NPARAMETR:  6.8722E-01  4.2238E-02  5.8878E-02  5.4807E-01  4.1272E+02  1.8258E+00  1.8736E-01  1.0000E-02  1.0731E+00  1.6699E+00
             7.8421E+00
 PARAMETER: -2.7511E-01 -3.0644E+00 -2.7323E+00 -5.0135E-01  6.1228E+00  7.0200E-01 -1.5747E+00 -1.6482E+01  1.7057E-01  6.1278E-01
             2.1595E+00
 GRADIENT:   3.0721E+01  1.1130E+01 -3.2570E+01  2.3134E+01 -1.4256E-02  2.5017E-01  1.1541E-01  0.0000E+00 -2.2941E+00  4.3056E-05
            -3.1582E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -867.035022653283        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      683            RESET HESSIAN, TYPE II
 NPARAMETR:  5.9032E-01  2.1409E-02  4.4637E-02  4.3898E-01  8.0682E+02  1.7844E+00  3.1899E-02  1.0000E-02  1.0538E+00  1.9064E+00
             8.0400E+00
 PARAMETER: -4.2710E-01 -3.7439E+00 -3.0092E+00 -7.2330E-01  6.7931E+00  6.7911E-01 -3.3452E+00 -1.8346E+01  1.5244E-01  7.4520E-01
             2.1844E+00
 GRADIENT:   4.1149E+01  3.8131E-01  5.5763E+01  2.2346E+01  1.1784E-03  1.6648E+01  1.9418E-05  0.0000E+00  1.4479E+00  5.3365E-08
             1.7396E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -867.049782062718        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  5.8408E-01  1.8896E-02  4.4022E-02  4.3500E-01  8.0560E+02  1.7763E+00  1.0000E-02  1.0000E-02  1.0499E+00  1.9064E+00
             8.0335E+00
 PARAMETER: -4.3772E-01 -3.8688E+00 -3.0231E+00 -7.3242E-01  6.7916E+00  6.7454E-01 -5.2184E+00 -1.8346E+01  1.4867E-01  7.4524E-01
             2.1836E+00
 GRADIENT:  -2.6813E+00  4.7198E-02  1.9531E+00 -2.2679E+00  1.2882E-03 -1.0258E+00  0.0000E+00  0.0000E+00 -5.0884E-02 -1.3925E-07
            -3.0180E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -867.054477447396        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      966
 NPARAMETR:  5.8762E-01  1.9619E-02  4.4325E-02  4.3779E-01  8.0435E+02  1.7813E+00  1.2957E-02  1.0000E-02  1.0506E+00  1.9064E+00
             8.0396E+00
 PARAMETER: -4.3167E-01 -3.8313E+00 -3.0162E+00 -7.2601E-01  6.7900E+00  6.7733E-01 -4.2461E+00 -1.8346E+01  1.4938E-01  7.4522E-01
             2.1844E+00
 GRADIENT:  -6.8194E-01  7.2384E-02 -3.7004E-02 -4.6790E-01  1.4565E-03 -2.8992E-01  8.1130E-07  0.0000E+00 -1.1808E-02 -1.2796E-07
            -2.2869E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -867.058016195863        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1153             RESET HESSIAN, TYPE I
 NPARAMETR:  5.8833E-01  1.8630E-02  4.4279E-02  4.3785E-01  5.6662E+02  1.7848E+00  1.0931E-02  1.0000E-02  1.0504E+00  1.9397E+00
             8.0416E+00
 PARAMETER: -4.3046E-01 -3.8830E+00 -3.0172E+00 -7.2588E-01  6.4397E+00  6.7930E-01 -4.4162E+00 -1.8346E+01  1.4921E-01  7.6256E-01
             2.1846E+00
 GRADIENT:   4.1829E+01  6.5738E-02  5.3279E+01  2.6680E+01  2.2657E-03  1.6619E+01  1.0253E-06  0.0000E+00  7.6347E-01 -8.5337E-08
             1.7182E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -867.059543796374        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1334
 NPARAMETR:  5.8812E-01  1.8626E-02  4.4273E-02  4.3759E-01  3.8726E+02  1.7847E+00  1.0620E-02  1.0000E-02  1.0506E+00  1.9002E+00
             8.0426E+00
 PARAMETER: -4.3083E-01 -3.8832E+00 -3.0174E+00 -7.2648E-01  6.0591E+00  6.7922E-01 -4.4450E+00 -1.8346E+01  1.4933E-01  7.4198E-01
             2.1847E+00
 GRADIENT:   3.6261E-01  2.2376E-03 -5.5952E-01 -2.4283E-01  3.3300E-03  2.6410E-01  3.5873E-07  0.0000E+00  1.0266E-02 -6.1981E-07
            -1.1605E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -867.123383998720        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1511
 NPARAMETR:  5.8850E-01  1.5244E-02  4.4602E-02  4.3962E-01  1.2909E+01  1.7748E+00  1.0366E-02  1.0000E-02  1.0483E+00  1.6226E+00
             8.0467E+00
 PARAMETER: -4.3018E-01 -4.0836E+00 -3.0100E+00 -7.2185E-01  2.6579E+00  6.7366E-01 -4.4693E+00 -1.8346E+01  1.4716E-01  5.8403E-01
             2.1853E+00
 GRADIENT:  -2.3508E+00 -2.9665E-02  2.9159E+00 -3.5151E+00  5.5915E-02 -2.0709E+00  2.2653E-07  0.0000E+00  2.2324E-01 -3.2896E-04
             8.4774E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -867.149432439187        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1696             RESET HESSIAN, TYPE I
 NPARAMETR:  5.9113E-01  1.5916E-02  4.4536E-02  4.4057E-01  1.1354E+01  1.7886E+00  1.0000E-02  1.0000E-02  1.0483E+00  1.2940E+00
             8.0411E+00
 PARAMETER: -4.2571E-01 -4.0405E+00 -3.0114E+00 -7.1969E-01  2.5296E+00  6.8146E-01 -4.5619E+00 -1.8346E+01  1.4714E-01  3.5776E-01
             2.1846E+00
 GRADIENT:   4.0806E+01  6.2120E-02  5.3837E+01  2.6091E+01  1.0038E-02  1.6765E+01  0.0000E+00  0.0000E+00  7.7320E-01  3.9537E-04
             1.7312E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -867.150193615256        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1885             RESET HESSIAN, TYPE I
 NPARAMETR:  5.9097E-01  1.5920E-02  4.4494E-02  4.4018E-01  1.1411E+01  1.7884E+00  1.0000E-02  1.0000E-02  1.0483E+00  1.1585E+00
             8.0420E+00
 PARAMETER: -4.2600E-01 -4.0402E+00 -3.0124E+00 -7.2057E-01  2.5346E+00  6.8131E-01 -4.6495E+00 -1.8346E+01  1.4717E-01  2.4713E-01
             2.1847E+00
 GRADIENT:   4.0876E+01  5.9886E-02  5.4029E+01  2.5861E+01  1.4453E-02  1.6735E+01  0.0000E+00  0.0000E+00  7.9653E-01  2.5520E-04
             1.7408E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -867.150827082146        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2077
 NPARAMETR:  5.9071E-01  1.5903E-02  4.4451E-02  4.3981E-01  1.1397E+01  1.7883E+00  1.0000E-02  1.0000E-02  1.0483E+00  1.1277E+00
             8.0421E+00
 PARAMETER: -4.2643E-01 -4.0412E+00 -3.0134E+00 -7.2142E-01  2.5334E+00  6.8125E-01 -4.6495E+00 -1.8346E+01  1.4721E-01  2.2020E-01
             2.1847E+00
 GRADIENT:   1.9730E-01  1.0918E-02 -5.2627E-01 -3.4358E-01 -1.1575E-02  3.2763E-01  0.0000E+00  0.0000E+00  1.4514E-02  1.8033E-04
            -2.9882E-02

0ITERATION NO.:   76    OBJECTIVE VALUE:  -867.150827082146        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2099
 NPARAMETR:  5.9071E-01  1.5903E-02  4.4451E-02  4.3981E-01  1.1397E+01  1.7883E+00  1.0000E-02  1.0000E-02  1.0483E+00  1.1277E+00
             8.0421E+00
 PARAMETER: -4.2643E-01 -4.0412E+00 -3.0134E+00 -7.2142E-01  2.5334E+00  6.8125E-01 -4.6495E+00 -1.8346E+01  1.4721E-01  2.2020E-01
             2.1847E+00
 GRADIENT:   1.9730E-01  1.0918E-02 -5.2627E-01 -3.4358E-01 -1.1575E-02  3.2763E-01  0.0000E+00  0.0000E+00  1.4514E-02  1.8033E-04
            -2.9882E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2099
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.5823E-04  1.8456E-07  6.2196E-05 -2.1111E-02 -9.5884E-06
 SE:             2.9080E-02  1.0881E-06  2.0040E-04  2.5244E-02  1.3190E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7920E-01  8.6531E-01  7.5629E-01  4.0301E-01  9.4205E-01

 ETASHRINKSD(%)  2.5770E+00  9.9996E+01  9.9329E+01  1.5429E+01  9.9558E+01
 ETASHRINKVR(%)  5.0877E+00  1.0000E+02  9.9995E+01  2.8478E+01  9.9998E+01
 EBVSHRINKSD(%)  2.2358E+00  9.9995E+01  9.9392E+01  1.5548E+01  9.9560E+01
 EBVSHRINKVR(%)  4.4215E+00  1.0000E+02  9.9996E+01  2.8679E+01  9.9998E+01
 RELATIVEINF(%)  6.6188E+00  3.5057E-08  5.0600E-05  8.6270E-01  1.9220E-04
 EPSSHRINKSD(%)  1.6270E+01
 EPSSHRINKVR(%)  2.9893E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -867.15082708214561     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -132.00000051840743     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.15
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     8.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -867.151       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.91E-01  1.59E-02  4.45E-02  4.40E-01  1.14E+01  1.79E+00  1.00E-02  1.00E-02  1.05E+00  1.13E+00  8.04E+00
 


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
+        5.54E+01
 
 TH 2
+        1.77E+00  5.64E-02
 
 TH 3
+       -2.84E+03 -9.06E+01  1.46E+05
 
 TH 4
+        3.59E+02  1.14E+01 -1.84E+04  2.32E+03
 
 TH 5
+        4.54E-02  1.45E-03 -2.33E+00  2.94E-01  3.72E-05
 
 TH 6
+       -4.95E-01 -1.58E-02  2.54E+01 -3.21E+00 -4.06E-04  4.43E-03
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.95E+00 -1.90E-01  3.05E+02 -3.85E+01 -4.88E-03  5.32E-02  0.00E+00  0.00E+00  6.39E-01
 
 TH10
+        5.39E-04  1.72E-05 -2.76E-02  3.49E-03  4.42E-07 -4.82E-06  0.00E+00  0.00E+00 -5.79E-05  5.25E-09
 
 TH11
+       -1.62E+00 -5.15E-02  8.29E+01 -1.05E+01 -1.32E-03  1.45E-02  0.00E+00  0.00E+00  1.74E-01 -1.57E-05  4.71E-02
 
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
+        9.42E+02
 
 TH 2
+       -2.89E+02  1.52E+03
 
 TH 3
+       -3.81E+03 -1.15E+02  1.94E+05
 
 TH 4
+       -7.76E+01  1.17E+02 -2.44E+04  3.51E+03
 
 TH 5
+        2.36E-01 -1.06E+00 -3.09E+00  3.56E-01  5.06E-03
 
 TH 6
+       -5.07E+00  1.65E+01  3.25E+01 -2.03E+01  5.12E-03  5.06E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.12E+00 -1.79E+01  4.04E+02 -6.05E+01 -1.88E-02 -2.61E+00  0.00E+00  0.00E+00  9.08E+01
 
 TH10
+        1.58E-03  3.48E-02 -3.67E-02  3.49E-03 -2.10E-04  2.70E-03  0.00E+00  0.00E+00 -2.25E-02 -7.90E-03
 
 TH11
+       -1.26E+01  3.23E+00  1.10E+02 -9.67E+00 -4.27E-03  6.96E-01  0.00E+00  0.00E+00  3.36E+00 -6.47E-06  5.77E+00
 
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
+        9.52E+02
 
 TH 2
+       -2.10E+02  7.76E+01
 
 TH 3
+       -6.74E+03 -2.24E+02  2.58E+05
 
 TH 4
+        1.56E+02  1.71E+02 -2.88E+04  3.76E+03
 
 TH 5
+        2.54E-01 -3.62E-02 -4.34E+00  4.11E-01  1.25E-04
 
 TH 6
+        2.76E+02 -1.14E+01 -3.27E+03  1.54E+02  7.92E-02  2.14E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.12E+01 -3.87E+01  7.34E+02 -1.76E+02 -2.60E-02  1.90E+01  0.00E+00  0.00E+00  1.04E+02
 
 TH10
+        5.60E-03 -8.18E-04 -3.60E-02 -2.16E-03 -3.33E-08  4.27E-03  0.00E+00  0.00E+00  2.69E-03  2.01E-07
 
 TH11
+        2.64E+01 -5.21E+00  3.78E+02 -7.13E+01  2.87E-03  1.09E+01  0.00E+00  0.00E+00 -1.03E+01  2.43E-04  8.46E+01
 
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
 #CPUT: Total CPU Time in Seconds,       42.754
Stop Time:
Wed Sep 29 19:57:00 CDT 2021
