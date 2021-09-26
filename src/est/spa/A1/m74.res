Sat Sep 25 08:13:42 CDT 2021
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
$DATA ../../../../data/spa/A1/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1543.92500999939        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.1936E+01 -4.6826E+01  2.1586E+01 -9.8541E+01 -7.7790E+00  2.2529E+01 -1.6519E+01 -4.2416E+00 -2.1615E+01  9.4259E-01
            -2.0253E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1576.02490600951        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6930E-01  1.0624E+00  1.0024E+00  1.0846E+00  1.0141E+00  8.7518E-01  1.0736E+00  9.8434E-01  9.8556E-01  9.2656E-01
             1.4484E+00
 PARAMETER:  6.8817E-02  1.6049E-01  1.0237E-01  1.8122E-01  1.1401E-01 -3.3329E-02  1.7098E-01  8.4212E-02  8.5457E-02  2.3719E-02
             4.7043E-01
 GRADIENT:  -5.5922E+01  5.0194E+01  2.0298E+00  7.4497E+01 -4.0884E+00 -3.1236E+01 -8.5989E+00  4.0475E-01 -8.4856E+00  1.7885E+00
             8.4580E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1578.44138619334        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.7153E-01  9.1550E-01  7.3704E-01  1.1473E+00  7.7918E-01  8.9913E-01  1.4039E+00  7.5918E-01  8.8726E-01  5.9939E-01
             1.4516E+00
 PARAMETER:  7.1121E-02  1.1714E-02 -2.0511E-01  2.3740E-01 -1.4951E-01 -6.3255E-03  4.3927E-01 -1.7552E-01 -1.9617E-02 -4.1184E-01
             4.7263E-01
 GRADIENT:  -5.0724E+01  4.7448E+01  4.3879E-01  7.8985E+01 -9.2860E+00 -2.0116E+01  5.3255E+00  1.3770E+00 -4.9473E+00 -2.7877E+00
             1.1963E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1582.16099508864        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8837E-01  8.9640E-01  6.7290E-01  1.0953E+00  7.5028E-01  9.4219E-01  1.3130E+00  5.6070E-01  9.1589E-01  6.3973E-01
             1.3962E+00
 PARAMETER:  8.8300E-02 -9.3679E-03 -2.9616E-01  1.9106E-01 -1.8730E-01  4.0450E-02  3.7232E-01 -4.7857E-01  1.2137E-02 -3.4671E-01
             4.3374E-01
 GRADIENT:   6.2397E-01  8.4125E-01 -4.5041E-01  3.2316E+00  2.0794E+00  1.3271E-01 -4.9768E-01 -4.2049E-01  4.1351E-01 -8.0132E-01
            -6.1736E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1582.89212984570        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.9242E-01  7.1009E-01  8.0408E-01  1.2310E+00  7.4956E-01  9.4146E-01  1.5516E+00  6.6210E-01  8.7140E-01  7.0384E-01
             1.3966E+00
 PARAMETER:  9.2395E-02 -2.4236E-01 -1.1806E-01  3.0785E-01 -1.8827E-01  3.9676E-02  5.3931E-01 -3.1234E-01 -3.7649E-02 -2.5120E-01
             4.3401E-01
 GRADIENT:  -2.8316E+00  1.0446E+01  1.0003E+01  9.4485E+00 -1.7238E+01 -5.4742E-01 -1.4453E+00 -3.6543E-01 -6.6193E-01  7.5414E-01
            -6.0146E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1583.73318592368        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8547E-01  4.3211E-01  9.0425E-01  1.4028E+00  7.1926E-01  9.3549E-01  2.1898E+00  7.5311E-01  8.1671E-01  7.1400E-01
             1.3946E+00
 PARAMETER:  8.5365E-02 -7.3908E-01 -6.5339E-04  4.3845E-01 -2.2953E-01  3.3313E-02  8.8383E-01 -1.8355E-01 -1.0247E-01 -2.3687E-01
             4.3263E-01
 GRADIENT:  -7.3752E+00  7.9970E+00  1.0258E+01  1.9269E+01 -1.5729E+01 -8.3646E-01  9.1528E-01 -2.7985E-01 -1.9076E+00 -7.6309E-01
            -1.1534E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1584.81145452029        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.8593E-01  2.1602E-01  9.0673E-01  1.5143E+00  6.7187E-01  9.4099E-01  2.9919E+00  7.7437E-01  8.0385E-01  7.3584E-01
             1.3942E+00
 PARAMETER:  8.5827E-02 -1.4324E+00  2.0899E-03  5.1498E-01 -2.9769E-01  3.9182E-02  1.1959E+00 -1.5571E-01 -1.1834E-01 -2.0674E-01
             4.3233E-01
 GRADIENT:   6.1344E+00  8.2845E-01  3.1919E+00  1.1374E+01 -1.0176E+01  3.3680E+00 -4.1360E+00  1.3096E+00  1.7651E+00  1.4463E+00
             1.6653E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1585.73902064628        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.7913E-01  9.7529E-02  9.9855E-01  1.5860E+00  6.9488E-01  9.2772E-01  4.4059E+00  8.7282E-01  7.8965E-01  7.5262E-01
             1.3936E+00
 PARAMETER:  7.8905E-02 -2.2276E+00  9.8551E-02  5.6123E-01 -2.6401E-01  2.4974E-02  1.5830E+00 -3.6027E-02 -1.3616E-01 -1.8420E-01
             4.3187E-01
 GRADIENT:  -2.1447E+00 -5.0354E-01  1.4976E+00  3.6387E+00 -2.4857E+00 -8.1100E-01 -2.9538E+00  4.8924E-01  3.7444E+00  6.4590E-01
             7.1730E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1586.84661457118        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.7936E-01  1.6841E-02  1.0610E+00  1.6444E+00  7.0516E-01  9.2799E-01  9.9603E+00  9.1740E-01  7.5069E-01  7.0077E-01
             1.3921E+00
 PARAMETER:  7.9148E-02 -3.9839E+00  1.5922E-01  5.9738E-01 -2.4934E-01  2.5262E-02  2.3986E+00  1.3784E-02 -1.8676E-01 -2.5558E-01
             4.3078E-01
 GRADIENT:   3.2043E+00 -1.6588E+00  1.2264E+01  4.2789E+01 -9.8288E+00  3.5948E-01 -5.0495E+00 -2.4978E+00 -2.8709E+00 -2.9425E+00
            -4.7148E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1587.46671888409        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1243
 NPARAMETR:  9.7916E-01  1.0273E-02  9.9346E-01  1.6259E+00  6.7812E-01  9.2698E-01  1.2333E+01  8.6058E-01  7.4655E-01  6.8125E-01
             1.3908E+00
 PARAMETER:  7.8940E-02 -4.4782E+00  9.3434E-02  5.8608E-01 -2.8843E-01  2.4175E-02  2.6122E+00 -5.0151E-02 -1.9229E-01 -2.8383E-01
             4.2991E-01
 GRADIENT:   3.4728E+00 -1.4662E+00  3.2966E+00  1.6767E+01  3.9206E+00 -2.4412E-02 -3.2731E+00 -1.9852E+00 -5.4388E+00 -3.4420E+00
            -4.0541E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1587.93030822697        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1429
 NPARAMETR:  9.7842E-01  1.0032E-02  8.9796E-01  1.6035E+00  6.3374E-01  9.2829E-01  1.2344E+01  7.5761E-01  7.6404E-01  6.9724E-01
             1.3891E+00
 PARAMETER:  7.8181E-02 -4.5020E+00 -7.6339E-03  5.7217E-01 -3.5612E-01  2.5594E-02  2.6131E+00 -1.7759E-01 -1.6913E-01 -2.6062E-01
             4.2868E-01
 GRADIENT:   1.6623E+00 -9.1387E-01 -1.3516E+00  2.4941E-01  5.2581E+00  3.4368E-01 -3.0540E+00 -2.5224E-01  8.2256E-02 -5.5532E-01
            -5.6426E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1587.97909584391        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1605
 NPARAMETR:  9.7777E-01  1.0037E-02  8.6234E-01  1.5975E+00  6.1407E-01  9.2816E-01  1.2320E+01  7.2287E-01  7.6979E-01  6.9959E-01
             1.3870E+00
 PARAMETER:  7.7522E-02 -4.5015E+00 -4.8105E-02  5.6846E-01 -3.8764E-01  2.5447E-02  2.6112E+00 -2.2453E-01 -1.6163E-01 -2.5726E-01
             4.2713E-01
 GRADIENT:   1.4380E-01  3.2190E+00 -2.3262E+00 -6.6989E+00  2.9668E+00 -2.7103E-01  6.4295E+00 -5.1472E-01 -4.5137E+00 -6.9233E-01
            -3.7121E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1587.97966540441        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1778
 NPARAMETR:  9.7779E-01  1.0000E-02  8.6294E-01  1.5977E+00  6.1439E-01  9.2806E-01  1.2336E+01  7.2146E-01  7.6942E-01  7.0059E-01
             1.3871E+00
 PARAMETER:  7.7539E-02 -4.5053E+00 -4.7360E-02  5.6855E-01 -3.8709E-01  2.5403E-02  2.6128E+00 -2.2635E-01 -1.6213E-01 -2.5558E-01
             4.2722E-01
 GRADIENT:  -8.8647E-02  9.6252E+01  5.6080E-01 -1.5284E+03  2.2503E+03  5.8643E-02  3.2968E+02  8.0022E-02 -5.3792E+03  3.0160E-01
             7.9679E-02
 NUMSIGDIG:         3.5         3.3         2.5         3.3         3.3         2.4         3.3         2.5         3.3         2.2
                    3.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1778
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.7037E-04  1.4015E-02 -1.2284E-02 -8.1164E-03 -1.7423E-02
 SE:             2.9657E-02  7.1878E-03  1.5466E-02  2.8336E-02  1.9829E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9004E-01  5.1200E-02  4.2702E-01  7.7455E-01  3.7960E-01

 ETASHRINKSD(%)  6.4618E-01  7.5920E+01  4.8188E+01  5.0703E+00  3.3569E+01
 ETASHRINKVR(%)  1.2882E+00  9.4202E+01  7.3155E+01  9.8834E+00  5.5870E+01
 EBVSHRINKSD(%)  8.9097E-01  8.1448E+01  4.8654E+01  4.8329E+00  3.2514E+01
 EBVSHRINKVR(%)  1.7740E+00  9.6558E+01  7.3636E+01  9.4323E+00  5.4456E+01
 RELATIVEINF(%)  9.8124E+01  2.9684E+00  2.1064E+00  6.1245E+01  3.6310E+00
 EPSSHRINKSD(%)  4.1007E+01
 EPSSHRINKVR(%)  6.5199E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1587.9796654044080     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -852.82883884066985     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.93
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.69
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1587.980       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  1.00E-02  8.63E-01  1.60E+00  6.14E-01  9.28E-01  1.23E+01  7.22E-01  7.69E-01  7.01E-01  1.39E+00
 


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
+        1.32E+03
 
 TH 2
+       -1.66E+03  1.06E+08
 
 TH 3
+       -1.82E+02  2.90E+04  2.79E+03
 
 TH 4
+        4.08E+01  5.35E+06 -1.25E+03  2.64E+05
 
 TH 5
+       -2.31E+02 -2.04E+07  3.44E+03  2.20E+03  3.85E+06
 
 TH 6
+       -5.87E+01  2.21E+03  2.19E+02 -8.53E+01  3.43E+02  2.64E+02
 
 TH 7
+       -1.46E+00 -3.17E+02  3.08E+01  3.85E+01 -6.07E+01  1.95E+00  2.09E+02
 
 TH 8
+       -1.59E+01  4.08E+03  4.22E+02 -1.36E+02  6.01E+02  6.01E+01  3.16E+00  1.44E+02
 
 TH 9
+       -1.79E+07  3.88E+07  2.03E+07 -1.41E+03 -7.34E+06  1.88E+07  3.50E+01  1.07E+07  1.40E+07
 
 TH10
+       -6.08E+01  5.67E+03  8.38E+02 -1.79E+02  7.03E+02  8.21E+01  3.84E+00  2.39E+02  9.76E+06  4.27E+02
 
 TH11
+       -2.19E+01  1.69E+03  1.71E+02 -6.91E+01  2.47E+02  2.78E+01  1.43E+00  6.20E+01  2.95E+06  1.02E+02  1.32E+02
 
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
 #CPUT: Total CPU Time in Seconds,       29.686
Stop Time:
Sat Sep 25 08:14:13 CDT 2021
