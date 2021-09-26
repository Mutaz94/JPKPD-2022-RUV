Sat Sep 25 14:22:33 CDT 2021
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
$DATA ../../../../data/spa/D/dat47.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13964.8305697274        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9203E+02  2.5869E+02 -5.7394E+01  3.5025E+01  1.3026E+02 -1.5034E+03 -5.8816E+02 -4.1100E+01 -1.3007E+03 -3.9335E+02
            -2.7238E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -568.965222233935        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1948E+00  1.0447E+00  9.5796E-01  1.5343E+00  1.2868E+00  1.8203E+00  1.1437E+00  9.6567E-01  1.2350E+00  1.0506E+00
             1.4566E+01
 PARAMETER:  2.7802E-01  1.4375E-01  5.7047E-02  5.2811E-01  3.5217E-01  6.9900E-01  2.3425E-01  6.5063E-02  3.1108E-01  1.4936E-01
             2.7787E+00
 GRADIENT:  -3.8841E+01  8.9940E+00 -5.5637E+00  1.1507E+01 -5.8355E+00  3.9046E+01  1.2224E+00  4.0916E+00  1.1273E+01  2.4137E+00
             1.5346E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -577.952196510349        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.2394E+00  9.3744E-01  1.0696E+00  1.7879E+00  2.4290E+00  1.6449E+00  2.8926E+00  5.1802E-01  1.2550E+00  2.4101E+00
             1.3728E+01
 PARAMETER:  3.1465E-01  3.5392E-02  1.6725E-01  6.8102E-01  9.8746E-01  5.9770E-01  1.1622E+00 -5.5773E-01  3.2713E-01  9.7968E-01
             2.7194E+00
 GRADIENT:  -1.2417E+01  1.9192E+01 -8.9204E+00  5.0063E+01 -1.0030E+00 -1.5912E+00  5.3352E+00  8.9661E-01  1.4209E+01  3.2211E+00
             1.0430E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -604.992016970480        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1150E+00  4.8984E-01  1.5174E+00  1.6568E+00  4.4518E+00  1.4575E+00  1.7161E+00  5.8960E-01  8.4571E-01  6.5960E+00
             1.1965E+01
 PARAMETER:  2.0889E-01 -6.1367E-01  5.1701E-01  6.0489E-01  1.5933E+00  4.7673E-01  6.4004E-01 -4.2832E-01 -6.7579E-02  1.9865E+00
             2.5820E+00
 GRADIENT:  -8.2334E+00  1.0269E+01  9.9672E+00  3.3248E+01  2.8406E-01 -4.4280E+00  1.3264E-01 -9.9824E-03  1.8889E+00 -3.6218E+00
             3.2631E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -622.377143821941        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  9.6744E-01  2.2484E-01  5.6165E-01  1.3126E+00  7.5788E+00  1.2875E+00  8.5119E-01  1.0000E-02  3.4498E-01  6.3145E+00
             1.1277E+01
 PARAMETER:  6.6901E-02 -1.3924E+00 -4.7687E-01  3.7202E-01  2.1254E+00  3.5269E-01 -6.1119E-02 -5.1135E+00 -9.6426E-01  1.9428E+00
             2.5228E+00
 GRADIENT:   1.1387E+01  1.9683E+00  1.1664E+01 -5.4534E+01 -1.0129E+01 -1.4791E+01  1.1849E-01  0.0000E+00  1.6836E+00  1.2526E+01
            -5.3227E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -680.879374480325        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  7.3153E-01  7.1965E-02  8.4590E-02  7.7713E-01  1.8516E+01  1.5144E+00  1.0000E-02  1.0000E-02  2.1244E-01  6.8467E-01
             1.2243E+01
 PARAMETER: -2.1261E-01 -2.5316E+00 -2.3699E+00 -1.5215E-01  3.0187E+00  5.1500E-01 -5.3142E+00 -1.6964E+01 -1.4491E+00 -2.7882E-01
             2.6049E+00
 GRADIENT:   2.1607E+01  8.2034E+00 -3.4253E+00  1.7461E+01 -4.2972E-01  6.5205E+00  0.0000E+00  0.0000E+00  2.3861E+00  2.6026E-02
             8.3571E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -700.076653312601        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  5.1177E-01  3.4347E-02  3.6331E-02  4.3438E-01  4.5122E+01  1.2581E+00  1.0000E-02  1.0000E-02  5.3288E-02  7.3981E-01
             1.0140E+01
 PARAMETER: -5.6988E-01 -3.2712E+00 -3.2151E+00 -7.3384E-01  3.9094E+00  3.2960E-01 -5.9336E+00 -2.1839E+01 -2.8320E+00 -2.0136E-01
             2.4165E+00
 GRADIENT:   1.5146E+01  1.0349E+01 -3.2489E+01  4.8449E+01 -1.9990E-01 -3.5245E+01  0.0000E+00  0.0000E+00  1.8792E-01  3.7234E-03
            -2.6569E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -708.655110301986        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  3.0378E-01  1.0000E-02  1.0000E-02  1.4996E-01  2.9925E+02  1.2409E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0857E+00
             1.0060E+01
 PARAMETER: -1.0914E+00 -4.5992E+00 -4.5064E+00 -1.7974E+00  5.8013E+00  3.1587E-01 -9.5122E+00 -3.2493E+01 -6.1676E+00  1.8225E-01
             2.4086E+00
 GRADIENT:   2.9614E+01  0.0000E+00 -1.0082E+01  1.1439E+01 -1.9123E-02 -1.8819E+01  0.0000E+00  0.0000E+00  0.0000E+00  1.5986E-05
            -1.5308E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -709.402654625938        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      611
 NPARAMETR:  2.9761E-01  1.0000E-02  1.0000E-02  1.4968E-01  3.0430E+02  1.3103E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0326E+00
             1.0151E+01
 PARAMETER: -1.1120E+00 -4.6079E+00 -4.5244E+00 -1.7993E+00  5.8180E+00  3.7027E-01 -9.7564E+00 -3.2861E+01 -6.1403E+00  1.3210E-01
             2.4176E+00
 GRADIENT:  -8.0268E+00  0.0000E+00  0.0000E+00  1.7577E+00 -2.5946E-02 -6.0143E-01  0.0000E+00  0.0000E+00  0.0000E+00  1.4350E-05
            -2.4811E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -709.443616466956        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      792
 NPARAMETR:  3.0008E-01  1.0000E-02  1.0000E-02  1.4962E-01  3.3828E+02  1.3152E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0987E+00
             1.0194E+01
 PARAMETER: -1.1037E+00 -4.6027E+00 -4.5228E+00 -1.7997E+00  5.9239E+00  3.7401E-01 -9.8707E+00 -3.2908E+01 -6.1430E+00  1.9409E-01
             2.4218E+00
 GRADIENT:   5.2856E-01  0.0000E+00  0.0000E+00  3.1874E-02 -2.2157E-02  5.5285E-01  0.0000E+00  0.0000E+00  0.0000E+00  1.2904E-05
            -7.9196E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -709.466277226834        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      961
 NPARAMETR:  2.9988E-01  1.0000E-02  1.0000E-02  1.4959E-01  1.5565E+05  1.3130E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.0404E+00
             1.0201E+01
 PARAMETER: -1.1044E+00 -4.5124E+00 -4.5235E+00 -1.7999E+00  1.2055E+01  3.7232E-01 -1.1712E+01 -3.3368E+01 -5.9295E+00  1.4963E+00
             2.4225E+00
 GRADIENT:   8.1628E+00  0.0000E+00  0.0000E+00  3.0877E+00 -4.6454E-05  4.3079E-01  0.0000E+00  0.0000E+00  0.0000E+00  7.5977E-10
             1.8039E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -709.466352591827        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1140
 NPARAMETR:  2.9994E-01  1.0000E-02  1.0000E-02  1.4959E-01  7.5215E+05  1.3129E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.0404E+00
             1.0203E+01
 PARAMETER: -1.1042E+00 -4.5129E+00 -4.5235E+00 -1.7998E+00  1.3631E+01  3.7225E-01 -1.1712E+01 -3.3368E+01 -5.9295E+00  1.4963E+00
             2.4226E+00
 GRADIENT:   7.3518E-03  0.0000E+00  0.0000E+00 -5.3816E-04 -9.6522E-06  6.5680E-03  0.0000E+00  0.0000E+00  0.0000E+00  3.7988E-11
            -2.3988E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -709.466362428693        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  2.9994E-01  1.0000E-02  1.0000E-02  1.4959E-01  6.6769E+07  1.3129E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.0404E+00
             1.0203E+01
 PARAMETER: -1.1042E+00 -4.5153E+00 -4.5235E+00 -1.7998E+00  1.8117E+01  3.7223E-01 -1.1712E+01 -3.3368E+01 -5.9295E+00  1.4963E+00
             2.4227E+00
 GRADIENT:  -3.3822E-05  0.0000E+00  0.0000E+00  1.3006E-03 -1.0873E-07  1.4303E-04  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
             4.0569E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -709.466362613590        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1494
 NPARAMETR:  2.9995E-01  1.0000E-02  1.0000E-02  1.4959E-01  5.2838E+08  1.3129E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.0420E+00
             1.0203E+01
 PARAMETER: -1.1041E+00 -4.5165E+00 -4.5235E+00 -1.7998E+00  2.0185E+01  3.7223E-01 -1.1712E+01 -3.3368E+01 -5.9295E+00  1.4967E+00
             2.4227E+00
 GRADIENT:   2.6235E-02  0.0000E+00  0.0000E+00 -9.4688E-03 -1.3740E-08 -6.5877E-05  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -7.6545E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1494
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0767E-03  7.8832E-06  5.9116E-05 -1.3149E-04  1.2715E-12
 SE:             2.8770E-02  2.0385E-05  2.9473E-04  3.3578E-04  1.5742E-11
 N:                     100         100         100         100         100

 P VAL.:         9.7015E-01  6.9896E-01  8.4103E-01  6.9535E-01  9.3563E-01

 ETASHRINKSD(%)  3.6185E+00  9.9932E+01  9.9013E+01  9.8875E+01  1.0000E+02
 ETASHRINKVR(%)  7.1061E+00  1.0000E+02  9.9990E+01  9.9987E+01  1.0000E+02
 EBVSHRINKSD(%)  4.0617E+00  9.9890E+01  9.8951E+01  9.8783E+01  1.0000E+02
 EBVSHRINKVR(%)  7.9585E+00  1.0000E+02  9.9989E+01  9.9985E+01  1.0000E+02
 RELATIVEINF(%)  1.9805E+00  1.8280E-05  2.2102E-05  2.9253E-05  0.0000E+00
 EPSSHRINKSD(%)  7.5961E+00
 EPSSHRINKVR(%)  1.4615E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -709.46636261358969     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       25.684463950148483     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -709.466       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.00E-01  1.00E-02  1.00E-02  1.50E-01  5.28E+08  1.31E+00  1.00E-02  1.00E-02  1.00E-02  4.04E+00  1.02E+01
 


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
+        6.47E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        0.00E+00  0.00E+00  0.00E+00
 
 TH 4
+       -6.05E+02  0.00E+00  0.00E+00  5.22E+04
 
 TH 5
+       -1.05E-12  0.00E+00  0.00E+00  1.09E-12  3.85E-24
 
 TH 6
+       -2.31E+01  0.00E+00  0.00E+00 -1.23E+02 -2.06E-13  9.04E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        1.00E-03  0.00E+00  0.00E+00  1.06E-03 -8.26E-15  6.84E-04  0.00E+00  0.00E+00  0.00E+00 -1.90E-05
 
 TH11
+       -4.17E+01  0.00E+00  0.00E+00 -1.78E+01  1.40E-14  6.09E-01  0.00E+00  0.00E+00  0.00E+00  1.07E-05  4.17E+00
 
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
 #CPUT: Total CPU Time in Seconds,       24.087
Stop Time:
Sat Sep 25 14:23:08 CDT 2021
