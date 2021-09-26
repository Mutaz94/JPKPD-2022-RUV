Fri Sep 24 21:17:43 CDT 2021
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
$DATA ../../../../data/int/A2/dat27.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m27.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2735.40945203082        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.7520E+01  6.3585E+01  8.2058E+01  2.7476E+01  1.6625E+02  2.4033E+01 -9.8747E+01 -4.9661E+01 -2.3876E+01 -6.7322E+01
            -1.9996E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3235.15183058307        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8425E-01  8.4575E-01  7.9652E-01  1.0856E+00  7.8579E-01  8.8685E-01  1.0156E+00  7.7208E-01  1.0803E+00  8.7829E-01
             1.7134E+00
 PARAMETER:  8.4124E-02 -6.7534E-02 -1.2751E-01  1.8212E-01 -1.4106E-01 -2.0079E-02  1.1547E-01 -1.5867E-01  1.7726E-01 -2.9779E-02
             6.3846E-01
 GRADIENT:   2.2631E+01 -2.4565E+00 -1.5461E+01  2.9831E+01  3.4622E+01 -2.0501E+01 -3.6172E+00  1.6271E+01  1.0980E+01 -1.3195E+00
            -2.4469E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3240.28887743552        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0011E+00  7.4354E-01  6.7923E-01  1.1100E+00  6.4570E-01  9.3708E-01  1.2217E+00  3.0348E-01  1.0821E+00  8.5295E-01
             1.7130E+00
 PARAMETER:  1.0113E-01 -1.9633E-01 -2.8679E-01  2.0436E-01 -3.3741E-01  3.5010E-02  3.0021E-01 -1.0924E+00  1.7889E-01 -5.9050E-02
             6.3824E-01
 GRADIENT:   6.6480E+01  3.5921E+01  2.8488E+01 -4.0784E+00 -6.0855E+01  7.1983E-01  1.5717E+01  2.4426E+00  1.4655E+01  8.7297E+00
            -1.1432E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3243.99446678241        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8002E-01  7.3378E-01  6.7340E-01  1.1137E+00  6.5826E-01  9.2468E-01  1.0720E+00  1.8278E-01  1.0436E+00  8.2353E-01
             1.7295E+00
 PARAMETER:  7.9820E-02 -2.0954E-01 -2.9542E-01  2.0771E-01 -3.1815E-01  2.1690E-02  1.6949E-01 -1.5995E+00  1.4269E-01 -9.4158E-02
             6.4784E-01
 GRADIENT:   1.2129E+01  1.0031E+01  3.8586E-01  1.7155E+00 -4.2883E+00 -3.4330E+00 -2.2261E+00  7.4924E-01  1.9280E+00 -1.8837E+00
             2.1664E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3243.99803439960        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.7858E-01  7.2374E-01  6.6566E-01  1.1166E+00  6.5051E-01  9.2683E-01  1.0820E+00  1.6970E-01  1.0413E+00  8.1914E-01
             1.7284E+00
 PARAMETER:  7.8345E-02 -2.2332E-01 -3.0697E-01  2.1025E-01 -3.2999E-01  2.4012E-02  1.7877E-01 -1.6737E+00  1.4047E-01 -9.9504E-02
             6.4719E-01
 GRADIENT:   8.5398E+00  7.2274E+00 -2.0686E-01  1.2315E+00 -2.6689E+00 -2.5467E+00 -1.7063E+00  6.6298E-01  1.4336E+00 -1.3652E+00
             1.7806E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3243.99863042023        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.7810E-01  7.1993E-01  6.6268E-01  1.1176E+00  6.4750E-01  9.2761E-01  1.0858E+00  1.6140E-01  1.0405E+00  8.1751E-01
             1.7280E+00
 PARAMETER:  7.7853E-02 -2.2860E-01 -3.1146E-01  2.1122E-01 -3.3463E-01  2.4859E-02  1.8232E-01 -1.7238E+00  1.3968E-01 -1.0149E-01
             6.4695E-01
 GRADIENT:   7.3592E+00  6.2648E+00 -2.9777E-01  1.0660E+00 -2.2187E+00 -2.2263E+00 -1.5013E+00  6.0277E-01  1.2545E+00 -1.1851E+00
             1.5939E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3243.99902578389        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.7781E-01  7.1757E-01  6.6084E-01  1.1183E+00  6.4563E-01  9.2810E-01  1.0882E+00  1.5520E-01  1.0400E+00  8.1653E-01
             1.7277E+00
 PARAMETER:  7.7562E-02 -2.3188E-01 -3.1425E-01  2.1182E-01 -3.3753E-01  2.5383E-02  1.8452E-01 -1.7631E+00  1.3920E-01 -1.0270E-01
             6.4680E-01
 GRADIENT:   6.6661E+00  5.6884E+00 -3.1595E-01  9.6735E-01 -1.9791E+00 -2.0286E+00 -1.3717E+00  5.5827E-01  1.1437E+00 -1.0767E+00
             1.4672E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3243.99918106841        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.7765E-01  7.1618E-01  6.5974E-01  1.1187E+00  6.4452E-01  9.2839E-01  1.0896E+00  1.5112E-01  1.0397E+00  8.1595E-01
             1.7276E+00
 PARAMETER:  7.7394E-02 -2.3382E-01 -3.1591E-01  2.1217E-01 -3.3925E-01  2.5693E-02  1.8582E-01 -1.7897E+00  1.3891E-01 -1.0340E-01
             6.4672E-01
 GRADIENT:   6.2673E+00  5.3530E+00 -3.1511E-01  9.0990E-01 -1.8487E+00 -1.9118E+00 -1.2942E+00  5.2965E-01  1.0781E+00 -1.0134E+00
             1.3885E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3243.99924783796        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  9.7755E-01  7.1531E-01  6.5906E-01  1.1190E+00  6.4383E-01  9.2857E-01  1.0905E+00  1.4841E-01  1.0395E+00  8.1560E-01
             1.7275E+00
 PARAMETER:  7.7291E-02 -2.3503E-01 -3.1694E-01  2.1239E-01 -3.4032E-01  2.5887E-02  1.8663E-01 -1.8077E+00  1.3874E-01 -1.0383E-01
             6.4667E-01
 GRADIENT:   6.0223E+00  5.1459E+00 -3.1037E-01  8.7454E-01 -1.7715E+00 -1.8390E+00 -1.2455E+00  5.1096E-01  1.0371E+00 -9.7423E-01
             1.3381E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3243.99927793872        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      675
 NPARAMETR:  9.7748E-01  7.1473E-01  6.5860E-01  1.1191E+00  6.4336E-01  9.2869E-01  1.0911E+00  1.4651E-01  1.0394E+00  8.1537E-01
             1.7274E+00
 PARAMETER:  7.7222E-02 -2.3586E-01 -3.1764E-01  2.1254E-01 -3.4105E-01  2.6018E-02  1.8718E-01 -1.8207E+00  1.3862E-01 -1.0412E-01
             6.4663E-01
 GRADIENT:   5.8581E+00  5.0066E+00 -3.0539E-01  8.5081E-01 -1.7210E+00 -1.7897E+00 -1.2124E+00  4.9796E-01  1.0094E+00 -9.4787E-01
             1.3034E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3244.24319124176        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  9.7759E-01  7.3018E-01  6.7849E-01  1.1165E+00  6.6345E-01  9.3665E-01  1.0892E+00  1.4225E-01  1.0376E+00  8.3487E-01
             1.7314E+00
 PARAMETER:  7.7336E-02 -2.1447E-01 -2.8788E-01  2.1021E-01 -3.1030E-01  3.4552E-02  1.8541E-01 -1.8502E+00  1.3689E-01 -8.0479E-02
             6.4891E-01
 GRADIENT:  -7.9786E+00 -4.9311E+00 -6.9787E-01 -1.7590E+00  5.0763E-01  4.2307E-01 -1.3252E-01  4.1069E-01 -3.9174E-05  2.5979E-01
             3.7649E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3244.47877110640        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      971
 NPARAMETR:  9.8100E-01  7.4403E-01  6.8765E-01  1.1131E+00  6.7405E-01  9.3572E-01  1.0838E+00  3.0450E-02  1.0388E+00  8.4407E-01
             1.7305E+00
 PARAMETER:  8.0821E-02 -1.9567E-01 -2.7447E-01  2.0711E-01 -2.9445E-01  3.3564E-02  1.8043E-01 -3.3917E+00  1.3802E-01 -6.9518E-02
             6.4841E-01
 GRADIENT:   6.0488E-01  6.0054E-02 -7.1030E-01  2.5883E-01  6.4945E-01  6.9238E-02  2.6223E-01  1.5276E-02 -2.8578E-02 -1.4690E-01
             1.4300E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3244.48644105288        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1133
 NPARAMETR:  9.8078E-01  7.4464E-01  6.8854E-01  1.1127E+00  6.7459E-01  9.3554E-01  1.0803E+00  1.0000E-02  1.0389E+00  8.4606E-01
             1.7305E+00
 PARAMETER:  8.0594E-02 -1.9485E-01 -2.7319E-01  2.0683E-01 -2.9365E-01  3.3365E-02  1.7723E-01 -4.6785E+00  1.3815E-01 -6.7160E-02
             6.4843E-01
 GRADIENT:   3.6904E-02  8.4741E-03 -1.3999E-02  1.0740E-02  8.3307E-03 -1.5945E-03  1.4908E-03  0.0000E+00  1.6536E-03  7.1640E-04
             6.0016E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1133
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1803E-04 -1.4043E-02 -1.3181E-04  2.8827E-03 -1.2601E-02
 SE:             2.9671E-02  2.1975E-02  2.2595E-04  2.8319E-02  2.4654E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7801E-01  5.2280E-01  5.5964E-01  9.1892E-01  6.0926E-01

 ETASHRINKSD(%)  5.9775E-01  2.6382E+01  9.9243E+01  5.1280E+00  1.7407E+01
 ETASHRINKVR(%)  1.1919E+00  4.5803E+01  9.9994E+01  9.9930E+00  3.1783E+01
 EBVSHRINKSD(%)  8.5224E-01  2.5923E+01  9.9202E+01  5.2124E+00  1.8416E+01
 EBVSHRINKVR(%)  1.6972E+00  4.5126E+01  9.9994E+01  1.0153E+01  3.3441E+01
 RELATIVEINF(%)  9.8277E+01  1.7050E+01  1.8121E-03  6.8307E+01  1.1610E+01
 EPSSHRINKSD(%)  1.8595E+01
 EPSSHRINKVR(%)  3.3732E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3244.4864410528780     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1590.3970812844673     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.94
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3244.486       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  7.45E-01  6.89E-01  1.11E+00  6.75E-01  9.36E-01  1.08E+00  1.00E-02  1.04E+00  8.46E-01  1.73E+00
 


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
+        1.30E+03
 
 TH 2
+       -3.93E+00  9.49E+02
 
 TH 3
+        3.80E-01 -1.55E+01  1.38E+03
 
 TH 4
+       -8.62E+00  1.76E+02 -7.92E+01  7.47E+02
 
 TH 5
+       -4.35E+00 -8.38E+02 -1.24E+03  1.75E+02  2.25E+03
 
 TH 6
+        4.73E+00 -2.99E+00  3.38E+00 -2.07E+00 -2.83E+00  2.20E+02
 
 TH 7
+        4.66E-01  1.95E+01 -2.11E+00  3.97E+00 -3.52E+00 -4.66E-02  5.08E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.10E+00 -2.19E+01  3.60E+01  1.24E+01 -1.60E+01 -1.66E+00  4.42E+00  0.00E+00  1.49E+02
 
 TH10
+       -1.37E+00 -1.37E+01 -5.82E+01 -8.91E-01 -1.70E+01  7.14E-02  2.31E+01  0.00E+00  2.41E-01  1.28E+02
 
 TH11
+       -1.33E+01 -2.05E+01 -2.46E+01 -1.24E+01  3.75E+00  2.33E+00  1.11E+01  0.00E+00  5.79E+00  1.28E+01  3.75E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       32.727
Stop Time:
Fri Sep 24 21:18:17 CDT 2021
