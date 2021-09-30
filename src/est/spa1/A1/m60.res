Wed Sep 29 22:40:13 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat60.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1330.13563605658        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8623E+02 -3.1245E+01  1.2044E+02  4.1580E+00  1.2439E+02  2.9541E+01 -2.1134E+01 -3.1127E+02  1.2259E+01 -2.1558E+01
            -1.1854E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1760.26317662464        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0804E+00  1.1239E+00  1.0303E+00  9.8940E-01  1.0443E+00  1.0100E+00  1.0825E+00  1.0328E+00  8.7608E-01  7.7405E-01
             2.3468E+00
 PARAMETER:  1.7737E-01  2.1684E-01  1.2990E-01  8.9348E-02  1.4334E-01  1.0997E-01  1.7924E-01  1.3224E-01 -3.2297E-02 -1.5612E-01
             9.5306E-01
 GRADIENT:   2.6588E+02 -1.4250E+01 -6.3792E+00 -1.1184E+01  1.7750E+01  1.1753E+01 -8.0133E+00  3.7476E+00 -3.2333E-02  8.6252E+00
             3.7843E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1762.90178825125        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0645E+00  6.9143E-01  6.1785E-01  1.2999E+00  6.0416E-01  9.8411E-01  1.8237E+00  6.8734E-01  8.2655E-01  2.5685E-01
             2.2430E+00
 PARAMETER:  1.6251E-01 -2.6899E-01 -3.8150E-01  3.6232E-01 -4.0392E-01  8.3983E-02  7.0087E-01 -2.7492E-01 -9.0493E-02 -1.2592E+00
             9.0782E-01
 GRADIENT:   2.2662E+02  5.0746E+01 -1.8656E+01  2.5230E+02  2.0271E+01  6.1551E-01  1.0082E+01  4.9445E-01  1.5651E+01 -3.2238E-01
             1.6808E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1774.52962798804        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      335
 NPARAMETR:  1.0083E+00  7.8209E-01  5.7081E-01  1.1447E+00  6.1738E-01  9.8283E-01  1.5793E+00  5.7957E-01  7.1475E-01  3.1462E-01
             2.1844E+00
 PARAMETER:  1.0822E-01 -1.4579E-01 -4.6069E-01  2.3518E-01 -3.8228E-01  8.2683E-02  5.5700E-01 -4.4547E-01 -2.3582E-01 -1.0564E+00
             8.8132E-01
 GRADIENT:  -1.6036E+01  6.0736E+00 -3.1514E+00  2.0159E+01  6.7996E+00 -1.2662E+00 -8.4900E-01 -9.3497E-01 -6.9094E+00  5.5301E-02
            -8.9936E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1775.31044991376        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  1.0169E+00  8.1367E-01  5.1671E-01  1.1025E+00  5.9207E-01  9.8484E-01  1.5077E+00  7.6587E-01  7.4463E-01  1.6447E-01
             2.1745E+00
 PARAMETER:  1.1677E-01 -1.0620E-01 -5.6028E-01  1.9762E-01 -4.2412E-01  8.4728E-02  5.1056E-01 -1.6674E-01 -1.9487E-01 -1.7050E+00
             8.7682E-01
 GRADIENT:   7.1268E-01 -2.6133E+00  7.3469E-01 -6.2958E+00  2.7897E-01 -2.0880E+00  1.1518E+00  8.4162E-01 -6.2421E-01  4.2167E-01
            -1.4446E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1775.66965053050        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  1.0185E+00  9.7033E-01  4.7642E-01  1.0184E+00  6.2432E-01  9.9373E-01  1.3043E+00  8.4245E-01  7.8041E-01  2.3271E-02
             2.1718E+00
 PARAMETER:  1.1836E-01  6.9880E-02 -6.4146E-01  1.1822E-01 -3.7109E-01  9.3707E-02  3.6565E-01 -7.1436E-02 -1.4793E-01 -3.6605E+00
             8.7556E-01
 GRADIENT:   6.7391E-01  8.9963E-01  3.2629E-01  4.9694E-01 -1.0643E+00 -3.0153E-01  4.8727E-02 -1.0960E-02 -8.5521E-02  1.4187E-02
             2.7320E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1809.65906285794        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  1.0160E+00  1.6124E+00  4.2243E-01  7.3297E-01  8.6684E-01  1.0210E+00  8.5999E-01  2.1210E+00  9.2821E-01  1.0000E-02
             2.0894E+00
 PARAMETER:  1.1584E-01  5.7775E-01 -7.6172E-01 -2.1065E-01 -4.2903E-02  1.2081E-01 -5.0838E-02  8.5188E-01  2.5499E-02 -1.2074E+01
             8.3689E-01
 GRADIENT:  -1.1053E+01  1.2653E+02  1.4822E+01  8.2380E+01 -2.4877E+01  5.0183E+00 -2.5408E+01 -4.1380E+01 -1.9182E+01  0.0000E+00
             1.0961E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1830.60605238551        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1045
 NPARAMETR:  1.0132E+00  1.5610E+00  4.2722E-01  6.9402E-01  9.0439E-01  1.0041E+00  9.0904E-01  2.6321E+00  1.1457E+00  1.0000E-02
             1.7574E+00
 PARAMETER:  1.1307E-01  5.4530E-01 -7.5046E-01 -2.6525E-01 -4.9779E-04  1.0408E-01  4.6360E-03  1.0678E+00  2.3599E-01 -1.3588E+01
             6.6382E-01
 GRADIENT:  -3.6252E+00 -1.4927E+01 -8.6296E+00  1.8678E+01  3.6384E+01 -6.7211E-01 -3.8746E-02 -7.0501E-01 -1.6250E+00  0.0000E+00
            -1.2241E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1832.08907997181        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  1.0160E+00  1.7917E+00  4.1520E-01  5.5696E-01  1.0164E+00  1.0050E+00  8.3414E-01  2.9397E+00  1.3332E+00  1.0000E-02
             1.8300E+00
 PARAMETER:  1.1585E-01  6.8315E-01 -7.7900E-01 -4.8526E-01  1.1626E-01  1.0495E-01 -8.1358E-02  1.1783E+00  3.8761E-01 -1.6488E+01
             7.0432E-01
 GRADIENT:   7.3411E-01  1.6519E+01  5.1258E+00  5.6418E+00 -5.1895E+00 -2.5485E-01 -2.4394E+00 -1.3930E+00 -1.0144E+00  0.0000E+00
            -1.6903E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1832.35789657025        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  1.0157E+00  1.8132E+00  3.9607E-01  5.3552E-01  1.0215E+00  1.0060E+00  8.2662E-01  2.9754E+00  1.3770E+00  1.0000E-02
             1.8354E+00
 PARAMETER:  1.1562E-01  6.9512E-01 -8.2617E-01 -5.2451E-01  1.2129E-01  1.0601E-01 -9.0409E-02  1.1904E+00  4.1993E-01 -1.7121E+01
             7.0726E-01
 GRADIENT:   2.9232E-01  7.5477E+00  4.8410E+00  8.2942E-03 -1.0334E+01  1.8329E-01 -2.2315E+00  4.4157E-01 -3.9628E-02  0.0000E+00
            -1.2179E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1832.42597397935        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1575            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0159E+00  1.8103E+00  3.8548E-01  5.3226E-01  1.0299E+00  1.0057E+00  8.3504E-01  2.9492E+00  1.3733E+00  1.0000E-02
             1.8384E+00
 PARAMETER:  1.1576E-01  6.9349E-01 -8.5326E-01 -5.3063E-01  1.2946E-01  1.0564E-01 -8.0279E-02  1.1815E+00  4.1724E-01 -1.7121E+01
             7.0891E-01
 GRADIENT:   1.3530E+02  1.9887E+02  6.9297E+00  3.1388E+01  1.1607E+01  1.1289E+01  2.0336E+00  5.0322E+01  3.7537E+00  0.0000E+00
             6.8926E+00

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1832.47442019242        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1714
 NPARAMETR:  1.0148E+00  1.8153E+00  3.8607E-01  5.3050E-01  1.0235E+00  1.0061E+00  8.3407E-01  2.9165E+00  1.3799E+00  1.0000E-02
             1.8385E+00
 PARAMETER:  1.1469E-01  6.9621E-01 -8.6030E-01 -5.3395E-01  1.2443E-01  1.0515E-01 -8.0611E-02  1.1821E+00  4.1786E-01 -1.7121E+01
             7.0896E-01
 GRADIENT:  -1.0026E+00 -7.4219E+00 -5.5817E+02 -5.4388E+00  3.8823E+03 -1.7818E-01  9.1893E-02  3.9421E+02 -5.7943E+02  0.0000E+00
             9.0654E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1714
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1657E-03 -1.9176E-02 -2.1157E-02  2.7891E-02 -6.1216E-04
 SE:             2.9656E-02  2.5176E-02  1.9419E-02  2.0953E-02  2.7512E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4178E-01  4.4626E-01  2.7594E-01  1.8316E-01  2.6079E-02

 ETASHRINKSD(%)  6.4832E-01  1.5656E+01  3.4943E+01  2.9804E+01  9.9078E+01
 ETASHRINKVR(%)  1.2924E+00  2.8861E+01  5.7676E+01  5.0725E+01  9.9992E+01
 EBVSHRINKSD(%)  1.0442E+00  1.5480E+01  3.5002E+01  3.1764E+01  9.9169E+01
 EBVSHRINKVR(%)  2.0774E+00  2.8564E+01  5.7753E+01  5.3439E+01  9.9993E+01
 RELATIVEINF(%)  9.7870E+01  1.3611E+01  1.5886E+01  8.6546E+00  1.8787E-03
 EPSSHRINKSD(%)  2.9938E+01
 EPSSHRINKVR(%)  5.0913E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1832.4744201924245     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -913.53588698775184     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1832.474       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.82E+00  3.83E-01  5.30E-01  1.02E+00  1.01E+00  8.35E-01  2.95E+00  1.37E+00  1.00E-02  1.84E+00
 


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
+        1.79E+06
 
 TH 2
+       -1.27E+01  1.55E+04
 
 TH 3
+       -8.70E+01 -2.90E+04  1.11E+05
 
 TH 4
+        3.66E+05  3.40E+04 -1.30E+05  3.03E+05
 
 TH 5
+        4.15E+00 -2.19E+02  7.36E+02  3.22E+02  7.42E+05
 
 TH 6
+        1.84E+00 -3.32E+00 -5.55E+01 -3.31E+00  1.45E+02  1.88E+02
 
 TH 7
+        2.37E+00  1.20E+01 -3.84E+02 -5.39E+00 -1.13E+06  1.92E-01  1.74E+06
 
 TH 8
+        8.79E+00 -2.75E+01  7.01E+01 -7.26E+01 -1.19E+02  5.81E+00  3.71E+01  9.45E+02
 
 TH 9
+        8.28E-01 -9.74E+00  1.44E+01  3.48E+01 -1.65E+05 -3.24E+01 -2.03E+02 -5.72E-01  3.67E+04
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        7.98E+04 -1.31E+01 -1.24E+01  1.03E+01 -7.27E+04  2.09E+00  6.22E+00  2.16E+00 -1.62E+04  0.00E+00  1.44E+04
 
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
 #CPUT: Total CPU Time in Seconds,       36.435
Stop Time:
Wed Sep 29 22:40:51 CDT 2021
