Sat Sep 25 05:48:08 CDT 2021
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
$DATA ../../../../data/int/D/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   29060.5497563292        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6973E+02  4.3568E+02  1.5580E+01  3.6185E+02 -3.8917E+01 -1.8839E+03 -1.1844E+03 -4.8503E+01 -1.6064E+03 -2.8406E+02
            -6.0163E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -776.551972951087        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.1618E+00  1.6093E+00  9.5125E-01  1.6153E+00  9.8796E-01  2.8078E+00  2.7455E+00  9.8566E-01  1.8949E+00  1.1341E+00
             1.3704E+01
 PARAMETER:  2.4998E-01  5.7578E-01  5.0026E-02  5.7952E-01  8.7892E-02  1.1324E+00  1.1100E+00  8.5551E-02  7.3916E-01  2.2581E-01
             2.7177E+00
 GRADIENT:  -3.5688E+01  2.3127E+01 -1.9395E+01  1.0125E+02 -1.0955E+01  7.1189E+01 -2.2450E+01  3.8923E+00  1.4567E+01  2.2132E+01
             5.0877E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -831.535734749635        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0147E+00  1.7227E+00  1.7078E+01  2.0332E+00  2.0133E+00  2.4656E+00  7.2076E+00  7.0263E-01  2.2501E+00  7.0623E-01
             1.3130E+01
 PARAMETER:  1.1464E-01  6.4391E-01  2.9378E+00  8.0960E-01  7.9976E-01  1.0024E+00  2.0751E+00 -2.5293E-01  9.1098E-01 -2.4781E-01
             2.6749E+00
 GRADIENT:  -6.6095E+01  3.6736E+01  7.4244E-01  2.9028E+01 -2.5617E+01  5.5537E+01  3.7858E+01  7.3964E-02  2.5997E+01  7.5445E+00
             5.2076E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1005.14742340306        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.1263E+00  2.6946E-01  8.8725E+01  1.8650E+00  2.3836E+00  2.1785E+00  5.8242E+00  3.7145E+00  2.1319E+00  4.9585E-01
             8.3461E+00
 PARAMETER:  2.1890E-01 -1.2113E+00  4.5855E+00  7.2325E-01  9.6860E-01  8.7865E-01  1.8620E+00  1.4123E+00  8.5702E-01 -6.0148E-01
             2.2218E+00
 GRADIENT:   1.4725E+01 -3.7277E-02 -2.0800E+00 -1.0062E+00  1.6981E+01 -7.1899E+00  1.4774E+00  4.9328E-01  2.2915E+00  2.7632E+00
             5.7555E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1006.76622660059        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0808E+00  2.4435E-01  9.5999E+01  1.8601E+00  2.3215E+00  2.2068E+00  5.7024E+00  8.9355E-01  2.1352E+00  3.2755E-01
             8.3050E+00
 PARAMETER:  1.7768E-01 -1.3092E+00  4.6643E+00  7.2063E-01  9.4222E-01  8.9153E-01  1.8409E+00 -1.2555E-02  8.5856E-01 -1.0161E+00
             2.2169E+00
 GRADIENT:  -2.0520E+00  6.2991E-02 -1.0490E+00 -2.2751E+00  1.3918E+00 -2.8284E+00  1.5979E+00  4.2646E-02  1.9655E+00  9.3702E-01
             1.7230E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1008.28913922536        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0862E+00  8.8487E-02  4.3207E+02  1.9727E+00  2.3901E+00  2.2244E+00  4.0782E+00  1.0000E-02  2.0646E+00  1.0000E-02
             8.3106E+00
 PARAMETER:  1.8268E-01 -2.3249E+00  6.1686E+00  7.7942E-01  9.7136E-01  8.9947E-01  1.5057E+00 -8.9811E+00  8.2492E-01 -5.0134E+00
             2.2175E+00
 GRADIENT:   2.1932E-01 -4.5210E-01 -1.6214E-01 -8.3850E-01  6.8726E-01  5.6992E-02  5.5475E-01  0.0000E+00  4.6150E-01  0.0000E+00
             1.3020E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1009.46789014901        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      455
 NPARAMETR:  1.0381E+00  3.0318E-01  1.4777E+02  1.8763E+00  2.3117E+00  2.2030E+00  2.1716E-02  1.0000E-02  2.1664E+00  1.0000E-02
             8.3341E+00
 PARAMETER:  1.3739E-01 -1.0934E+00  5.0957E+00  7.2928E-01  9.3799E-01  8.8984E-01 -3.7297E+00 -2.8851E+01  8.7305E-01 -7.3665E+00
             2.2204E+00
 GRADIENT:  -1.9980E+01  1.0953E+00 -4.1532E-01  2.0601E+01 -1.4744E+01 -3.8787E+00  3.9634E-04  0.0000E+00 -1.3896E+01  0.0000E+00
             6.6401E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1015.47655465369        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.0883E+00  8.4551E-01  1.7022E+02  1.3510E+00  2.3816E+00  2.2545E+00  1.0000E-02  1.0000E-02  2.9123E+00  1.0000E-02
             8.3756E+00
 PARAMETER:  1.8459E-01 -6.7815E-02  5.2371E+00  4.0082E-01  9.6779E-01  9.1292E-01 -1.7500E+01 -9.5317E+01  1.1689E+00 -2.1680E+01
             2.2253E+00
 GRADIENT:  -2.7699E-01  6.1144E+00 -6.6938E-01  4.0300E+00 -3.7374E+00  3.9796E+00  0.0000E+00  0.0000E+00 -7.7756E+00  0.0000E+00
             3.9885E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1016.15682575455        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  1.0899E+00  8.6053E-01  7.9308E+02  1.3151E+00  2.4149E+00  2.2408E+00  1.0000E-02  1.0000E-02  2.9943E+00  1.0000E-02
             8.3438E+00
 PARAMETER:  1.8604E-01 -5.0203E-02  6.7759E+00  3.7392E-01  9.8168E-01  9.0683E-01 -2.6919E+01 -1.4600E+02  1.1967E+00 -3.3620E+01
             2.2215E+00
 GRADIENT:   7.8940E-01  1.4673E+00 -1.2748E-01  3.6872E-01  2.0077E-01  1.5034E+00  0.0000E+00  0.0000E+00 -1.8580E+00  0.0000E+00
            -5.3881E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1016.28088245375        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  1.0887E+00  8.5708E-01  7.1955E+04  1.3143E+00  2.4182E+00  2.2303E+00  1.0000E-02  1.0000E-02  3.0052E+00  1.0000E-02
             8.3387E+00
 PARAMETER:  1.8499E-01 -5.4223E-02  1.1284E+01  3.7330E-01  9.8304E-01  9.0212E-01 -5.3395E+01 -2.8925E+02  1.2003E+00 -6.7353E+01
             2.2209E+00
 GRADIENT:  -1.2142E+00  1.1455E-01 -1.3040E-03  2.6504E-01  2.2853E-02 -2.9576E-01  0.0000E+00  0.0000E+00  8.5171E-02  0.0000E+00
             1.6204E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1016.31073009093        NO. OF FUNC. EVALS.: 152
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  1.0925E+00  8.6618E-01  1.5360E+05  1.3179E+00  2.4205E+00  2.2519E+00  1.0000E-02  1.0000E-02  3.0451E+00  1.0000E-02
             8.3569E+00
 PARAMETER:  1.8851E-01 -4.3661E-02  1.2042E+01  3.7605E-01  9.8399E-01  9.1177E-01 -5.7998E+01 -3.1404E+02  1.2135E+00 -7.3192E+01
             2.2231E+00
 GRADIENT:  -5.4770E+00  1.0841E+00 -6.0327E-04  3.0111E+00 -8.0932E-02  4.3949E+00  0.0000E+00  0.0000E+00  2.1875E+00  0.0000E+00
             5.3263E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1016.31225778591        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  1.0927E+00  8.6593E-01  1.5545E+05  1.3149E+00  2.4227E+00  2.2473E+00  1.0000E-02  1.0000E-02  3.0410E+00  1.0000E-02
             8.3578E+00
 PARAMETER:  1.8861E-01 -4.3955E-02  1.2054E+01  3.7375E-01  9.8487E-01  9.0971E-01 -5.7998E+01 -3.1404E+02  1.2122E+00 -7.3192E+01
             2.2232E+00
 GRADIENT:  -5.6586E+00  4.6903E-01 -6.1288E-04  1.5038E+00  1.0935E-01 -9.3361E-01  0.0000E+00  0.0000E+00  1.8417E-01  0.0000E+00
             1.1773E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1016.31265133902        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1095
 NPARAMETR:  1.0924E+00  8.6794E-01  1.5329E+05  1.3122E+00  2.4205E+00  2.2478E+00  1.0000E-02  1.0000E-02  3.0413E+00  1.0000E-02
             8.3565E+00
 PARAMETER:  1.8837E-01 -4.1632E-02  1.2040E+01  3.7174E-01  9.8399E-01  9.0997E-01 -5.7998E+01 -3.1404E+02  1.2123E+00 -7.3192E+01
             2.2230E+00
 GRADIENT:  -9.2399E+00  1.6524E+00 -6.1522E-04  1.9209E+00 -4.6867E-01  1.3458E+00  0.0000E+00  0.0000E+00 -8.7508E-01  0.0000E+00
             1.3516E+00

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1016.31293636597        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:     1206
 NPARAMETR:  1.0922E+00  8.6794E-01  1.5256E+05  1.3124E+00  2.4218E+00  2.2482E+00  1.0000E-02  1.0000E-02  3.0455E+00  1.0000E-02
             8.3563E+00
 PARAMETER:  1.8821E-01 -4.1728E-02  1.2040E+01  3.7182E-01  9.8445E-01  9.0980E-01 -5.7998E+01 -3.1404E+02  1.2140E+00 -7.3192E+01
             2.2232E+00
 GRADIENT:   3.1371E-02 -2.6214E-01  1.0906E-03 -2.1572E-03 -2.2297E-02 -1.0297E-01  0.0000E+00  0.0000E+00  8.4522E-02  0.0000E+00
             2.0143E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1206
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4519E-03 -4.8880E-04  6.4522E-09 -3.3149E-04 -1.6398E-05
 SE:             2.9017E-02  1.1445E-04  3.0824E-09  2.8232E-02  1.5180E-04
 N:                     100         100         100         100         100

 P VAL.:         7.7084E-01  1.9496E-05  3.6330E-02  9.9063E-01  9.1397E-01

 ETASHRINKSD(%)  2.7908E+00  9.9617E+01  1.0000E+02  5.4184E+00  9.9491E+01
 ETASHRINKVR(%)  5.5038E+00  9.9999E+01  1.0000E+02  1.0543E+01  9.9997E+01
 EBVSHRINKSD(%)  2.2438E+00  9.9735E+01  1.0000E+02  3.8023E+00  9.9479E+01
 EBVSHRINKVR(%)  4.4372E+00  9.9999E+01  1.0000E+02  7.4600E+00  9.9997E+01
 RELATIVEINF(%)  9.5369E+01  3.0530E-04  0.0000E+00  4.0492E+01  1.6389E-03
 EPSSHRINKSD(%)  6.8310E+00
 EPSSHRINKVR(%)  1.3195E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1016.3129363659739     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       637.77642340243688     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1016.313       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.09E+00  8.68E-01  1.53E+05  1.31E+00  2.42E+00  2.25E+00  1.00E-02  1.00E-02  3.05E+00  1.00E-02  8.36E+00
 


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
+        4.81E+02
 
 TH 2
+       -6.46E+02  1.82E+03
 
 TH 3
+        1.40E-05  2.35E-05  3.01E-12
 
 TH 4
+       -5.12E+01 -1.23E+02 -4.16E-06  1.27E+02
 
 TH 5
+       -5.37E+01 -6.62E+01  1.78E-06 -8.01E+00  4.75E+01
 
 TH 6
+       -4.62E+00 -4.98E+00  1.03E-07  4.07E+00  6.95E-01  3.39E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.55E+00 -5.01E+01 -1.48E-06  5.27E+00  1.36E+00 -4.69E-01  0.00E+00  0.00E+00  1.66E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.45E+00 -5.81E+00 -9.80E-08 -5.01E+00  9.49E-01  1.56E+00  0.00E+00  0.00E+00  2.06E+00  0.00E+00  1.52E+01
 
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
 #CPUT: Total CPU Time in Seconds,       42.849
Stop Time:
Sat Sep 25 05:48:52 CDT 2021
