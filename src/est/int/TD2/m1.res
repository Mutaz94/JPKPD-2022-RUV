Wed Sep 29 06:53:52 CDT 2021
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
$DATA ../../../../data/int/TD2/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3464.24150620592        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4199E+02 -1.4517E+01  1.6547E+02  1.6063E+02  5.2236E+01  5.0026E+01 -1.4609E+01 -8.6875E+01 -2.5571E+01 -1.0520E+01
            -6.7246E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3562.90387045254        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0448E+00  1.0490E+00  7.6319E-01  9.4165E-01  9.3467E-01  1.0353E+00  1.0159E+00  1.1594E+00  1.0515E+00  1.0525E+00
             1.2566E+00
 PARAMETER:  1.4384E-01  1.4788E-01 -1.7025E-01  3.9878E-02  3.2437E-02  1.3470E-01  1.1582E-01  2.4786E-01  1.5017E-01  1.5114E-01
             3.2842E-01
 GRADIENT:   4.8808E+02  5.4236E+01 -1.8011E+01  2.0434E+01 -4.9958E+00  4.9091E+01  6.5882E+00 -1.1263E+01  3.5648E+00 -3.2725E+00
            -1.6669E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3571.44781532996        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      263
 NPARAMETR:  1.0539E+00  1.3719E+00  9.3319E-01  8.2875E-01  1.3273E+00  9.1987E-01  8.8105E-01  1.6991E+00  1.2201E+00  1.4449E+00
             1.3316E+00
 PARAMETER:  1.5249E-01  4.1618E-01  3.0849E-02 -8.7832E-02  3.8316E-01  1.6474E-02 -2.6642E-02  6.3009E-01  2.9894E-01  4.6802E-01
             3.8635E-01
 GRADIENT:   1.1870E+02 -1.8537E+01 -5.0381E+00  7.1522E+01  4.9518E+01 -3.9591E+01  1.3928E+01  2.1257E+00  1.7013E+01  2.1464E+01
             9.6139E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3579.01640143078        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.0375E+00  1.3874E+00  9.4324E-01  8.0799E-01  1.3099E+00  9.3469E-01  8.3306E-01  1.8261E+00  1.1646E+00  1.3300E+00
             1.3044E+00
 PARAMETER:  1.3682E-01  4.2742E-01  4.1565E-02 -1.1320E-01  3.6993E-01  3.2462E-02 -8.2649E-02  7.0220E-01  2.5235E-01  3.8515E-01
             3.6574E-01
 GRADIENT:   7.5837E+01 -8.9295E+00 -3.4503E+00  5.1055E+01  3.2529E+01 -2.9563E+01  4.8997E+00  7.0327E+00  5.1604E+00  5.4810E-01
             5.6173E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3582.91642065883        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      631
 NPARAMETR:  1.0070E+00  1.4109E+00  9.3647E-01  7.8686E-01  1.3098E+00  9.9862E-01  7.9469E-01  1.7813E+00  1.1658E+00  1.3351E+00
             1.2920E+00
 PARAMETER:  1.0698E-01  4.4425E-01  3.4363E-02 -1.3970E-01  3.6989E-01  9.8622E-02 -1.2981E-01  6.7732E-01  2.5342E-01  3.8898E-01
             3.5620E-01
 GRADIENT:   5.7708E-01 -2.3889E+00 -5.3111E-03  3.6409E+01  1.9302E+01  9.9149E-02  8.9842E-03  3.8650E+00  1.3384E+00 -4.7586E-01
             3.3803E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3584.09398149681        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.0046E+00  1.4395E+00  9.1330E-01  7.5006E-01  1.3122E+00  9.9690E-01  7.9449E-01  1.7509E+00  1.1701E+00  1.3422E+00
             1.2897E+00
 PARAMETER:  1.0455E-01  4.6429E-01  9.3133E-03 -1.8760E-01  3.7171E-01  9.6895E-02 -1.3005E-01  6.6015E-01  2.5711E-01  3.9433E-01
             3.5444E-01
 GRADIENT:  -4.6448E+00 -1.2023E+01  2.2409E-01  3.9519E+00  2.6997E+00 -5.2567E-01  5.8965E-01  2.3000E+00 -1.8936E-01 -1.5701E+00
             2.6957E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3584.14702286024        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  1.0032E+00  1.4413E+00  9.1012E-01  7.4801E-01  1.3146E+00  9.9902E-01  7.8461E-01  1.7567E+00  1.1723E+00  1.3488E+00
             1.2875E+00
 PARAMETER:  1.0322E-01  4.6553E-01  5.8176E-03 -1.9034E-01  3.7354E-01  9.9017E-02 -1.4257E-01  6.6341E-01  2.5897E-01  3.9922E-01
             3.5267E-01
 GRADIENT:  -7.5494E+00 -1.3859E+01 -4.2511E-01  2.9849E+00  3.4903E+00  2.7917E-01 -5.7730E-01  2.4983E+00 -6.5611E-01 -8.5841E-01
             2.3388E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3584.25053547804        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1152
 NPARAMETR:  1.0061E+00  1.4540E+00  9.0820E-01  7.4073E-01  1.3134E+00  9.9710E-01  7.8760E-01  1.7538E+00  1.1730E+00  1.3475E+00
             1.2886E+00
 PARAMETER:  1.0610E-01  4.7433E-01  3.7109E-03 -2.0011E-01  3.7263E-01  9.7098E-02 -1.3877E-01  6.6179E-01  2.5960E-01  3.9824E-01
             3.5352E-01
 GRADIENT:  -1.1925E+00 -6.5000E+00  7.5792E-01  4.6031E-01 -4.1198E+00 -4.3268E-01 -9.9985E-02  2.1469E+00 -1.0338E+00 -1.7264E+00
             2.3742E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3584.28578371970        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1313
 NPARAMETR:  1.0056E+00  1.4533E+00  9.0588E-01  7.4015E-01  1.3159E+00  9.9818E-01  7.8910E-01  1.7596E+00  1.1715E+00  1.3502E+00
             1.2863E+00
 PARAMETER:  1.0561E-01  4.7385E-01  1.1569E-03 -2.0090E-01  3.7450E-01  9.8177E-02 -1.3687E-01  6.6511E-01  2.5830E-01  4.0025E-01
             3.5173E-01
 GRADIENT:  -2.1951E+00 -9.2728E+00 -1.1050E-01 -1.9258E-01 -2.0592E+00 -1.9375E-02  3.7512E-02  2.4443E+00 -1.1616E+00 -1.3964E+00
             2.0393E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3584.34921417925        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1495             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0082E+00  1.4606E+00  9.0022E-01  7.3482E-01  1.3172E+00  9.9897E-01  7.8319E-01  1.7626E+00  1.1827E+00  1.3642E+00
             1.2851E+00
 PARAMETER:  1.0812E-01  4.7887E-01 -5.1188E-03 -2.0813E-01  3.7547E-01  9.8965E-02 -1.4438E-01  6.6682E-01  2.6782E-01  4.1059E-01
             3.5080E-01
 GRADIENT:   3.0336E+02  3.9232E+02  1.4695E+00  6.4559E+01  9.1507E+01  3.0715E+01  6.3402E+00  6.3323E+00  1.0301E+01  1.9851E+01
             2.4495E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3584.43262166669        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1654            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0069E+00  1.4630E+00  9.0288E-01  7.3628E-01  1.3202E+00  9.9846E-01  7.8332E-01  1.7597E+00  1.1831E+00  1.3626E+00
             1.2794E+00
 PARAMETER:  1.0691E-01  4.8050E-01 -2.1626E-03 -2.0615E-01  3.7781E-01  9.8457E-02 -1.4421E-01  6.6516E-01  2.6817E-01  4.0936E-01
             3.4638E-01
 GRADIENT:   2.9995E+02  4.0190E+02  2.0109E+00  6.8233E+01  9.4525E+01  3.0850E+01  6.2838E+00  5.6757E+00  1.0351E+01  1.9538E+01
             1.4873E+01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -3584.43262166669        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     1715
 NPARAMETR:  1.0070E+00  1.4648E+00  9.0280E-01  7.3498E-01  1.3187E+00  9.9852E-01  7.8445E-01  1.7561E+00  1.1841E+00  1.3624E+00
             1.2808E+00
 PARAMETER:  1.0691E-01  4.8050E-01 -2.1626E-03 -2.0615E-01  3.7781E-01  9.8457E-02 -1.4421E-01  6.6516E-01  2.6817E-01  4.0936E-01
             3.4638E-01
 GRADIENT:  -3.9028E-01 -3.5607E+00  3.1315E-02  1.4651E+00  1.3255E+04 -3.8880E-02 -1.2402E-01  3.7297E+03 -2.3024E-01  6.0544E-02
            -1.4465E+04
 NUMSIGDIG:         2.8         2.4         2.9         1.9         2.3         3.0         1.8         2.3         2.3         3.2
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1715
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1480E-04 -4.3399E-02 -2.6023E-02  3.0851E-02 -3.1284E-02
 SE:             2.9855E-02  2.0972E-02  1.8901E-02  2.5316E-02  2.6088E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7823E-01  3.8508E-02  1.6857E-01  2.2297E-01  2.3047E-01

 ETASHRINKSD(%)  1.0000E-10  2.9741E+01  3.6679E+01  1.5190E+01  1.2601E+01
 ETASHRINKVR(%)  1.0000E-10  5.0637E+01  5.9905E+01  2.8072E+01  2.3614E+01
 EBVSHRINKSD(%)  4.1446E-01  2.9797E+01  3.8907E+01  1.7599E+01  1.1239E+01
 EBVSHRINKVR(%)  8.2720E-01  5.0715E+01  6.2677E+01  3.2101E+01  2.1214E+01
 RELATIVEINF(%)  9.9169E+01  1.3365E+01  3.3554E+01  2.0861E+01  3.4138E+01
 EPSSHRINKSD(%)  2.1211E+01
 EPSSHRINKVR(%)  3.7922E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3584.4326216666900     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1930.3432618982793     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    54.51
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3584.433       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.46E+00  9.03E-01  7.36E-01  1.32E+00  9.98E-01  7.83E-01  1.76E+00  1.18E+00  1.36E+00  1.28E+00
 


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
+        1.09E+03
 
 TH 2
+       -3.36E+00  2.54E+05
 
 TH 3
+        2.49E-01  3.67E+01  1.37E+02
 
 TH 4
+       -5.17E+00 -1.17E+06 -5.28E+01  9.57E+02
 
 TH 5
+        3.96E+01 -8.22E+02  1.64E+03 -1.02E+03  5.03E+05
 
 TH 6
+        2.16E+00 -1.19E+00  2.31E-01 -2.17E+00  5.78E+01  1.97E+02
 
 TH 7
+        5.18E-01 -3.10E+00 -3.58E+00 -2.43E+00 -2.65E+02 -3.26E-01  8.72E+01
 
 TH 8
+        1.74E+01 -3.00E+02  6.99E+02 -4.78E+02 -1.17E+02  2.45E+01 -1.10E+02  1.83E+05
 
 TH 9
+        3.67E+06  5.62E+05  4.37E+06 -2.60E+06 -7.92E+05 -5.74E-01  3.21E+01 -2.89E+02  1.25E+06
 
 TH10
+       -2.09E+06 -3.20E+05  1.29E+01  1.48E+06  4.50E+05  1.62E-01  6.82E+00 -1.38E+02 -7.08E+05  4.02E+05
 
 TH11
+       -5.44E+01  4.03E+05 -1.92E+03  1.30E+03 -5.67E+05 -6.30E+01  3.03E+02 -2.88E+02  8.92E+05 -5.07E+05  6.39E+05
 
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
 #CPUT: Total CPU Time in Seconds,       70.009
Stop Time:
Wed Sep 29 06:55:05 CDT 2021
