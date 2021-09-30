Thu Sep 30 00:33:31 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17.6191053478040        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0434E+02  8.5289E+01  1.2689E+02  2.0149E+01  2.9963E+02  2.7442E+01 -1.0359E+02 -5.6865E+01 -7.9730E+01 -2.2637E+02
            -3.7083E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1321.74858334320        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0034E+00  1.0036E+00  9.3855E-01  1.1781E+00  7.8737E-01  8.3288E-01  1.0670E+00  9.4588E-01  1.0701E+00  1.1647E+00
             6.6141E+00
 PARAMETER:  1.0336E-01  1.0362E-01  3.6581E-02  2.6391E-01 -1.3906E-01 -8.2862E-02  1.6488E-01  4.4365E-02  1.6777E-01  2.5242E-01
             1.9892E+00
 GRADIENT:  -1.3177E+02  1.0876E+01  9.2199E-01  1.1747E+01 -3.9261E+01 -3.1567E+01  8.9907E+00  5.5101E+00  2.2961E+01  2.9553E+01
             4.5665E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1389.05982033072        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.8397E-01  4.3283E-01  2.1880E-01  1.5086E+00  2.3190E-01  8.4932E-01  1.6241E+00  5.4771E-02  2.1091E+00  6.7399E-01
             4.7364E+00
 PARAMETER:  8.3837E-02 -7.3741E-01 -1.4196E+00  5.1119E-01 -1.3615E+00 -6.3321E-02  5.8493E-01 -2.8046E+00  8.4628E-01 -2.9453E-01
             1.6553E+00
 GRADIENT:  -8.3456E+01  8.2202E+01  4.9309E+01  1.5286E+02 -1.0016E+02 -4.2747E+01  1.4975E+01  4.3784E-02  6.4746E+01  1.6993E+01
             3.1615E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1496.18442836533        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      273
 NPARAMETR:  9.5597E-01  5.1769E-01  6.7280E-01  1.4100E+00  6.4135E-01  9.3384E-01  9.3060E-01  1.0000E-02  1.1357E+00  4.9519E-01
             3.7752E+00
 PARAMETER:  5.4975E-02 -5.5838E-01 -2.9630E-01  4.4356E-01 -3.4418E-01  3.1550E-02  2.8073E-02 -8.0966E+00  2.2724E-01 -6.0281E-01
             1.4285E+00
 GRADIENT:  -7.3474E+01 -1.9301E+01 -1.2412E+02  6.6890E+01  1.4720E+02 -1.0487E+01  1.3185E+00  0.0000E+00  2.7740E+01  6.2127E+00
             9.5289E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1522.40451303819        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.8136E-01  3.0632E-01  5.6065E-01  1.4260E+00  4.3217E-01  9.6310E-01  1.3782E+00  1.0000E-02  9.4520E-01  4.5975E-01
             3.2610E+00
 PARAMETER:  8.1186E-02 -1.0831E+00 -4.7865E-01  4.5484E-01 -7.3894E-01  6.2404E-02  4.2077E-01 -9.1829E+00  4.3637E-02 -6.7708E-01
             1.2820E+00
 GRADIENT:   2.4376E+00  2.1776E+01  5.5475E+01  7.1002E+01 -7.2463E+01  5.6894E-01 -3.7923E+00  0.0000E+00 -1.8129E+01 -1.0072E+01
            -1.7649E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1548.11388261326        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  9.7246E-01  3.1617E-01  2.5497E-01  1.1716E+00  2.5738E-01  9.5899E-01  1.3487E+00  1.0000E-02  1.2746E+00  7.3006E-01
             2.8539E+00
 PARAMETER:  7.2076E-02 -1.0515E+00 -1.2666E+00  2.5837E-01 -1.2572E+00  5.8128E-02  3.9911E-01 -9.3917E+00  3.4261E-01 -2.1462E-01
             1.1487E+00
 GRADIENT:  -2.1081E+01  1.6497E+01  1.1326E+01  2.2727E+01 -2.5674E+01 -3.7668E+00  1.4946E+00  0.0000E+00  1.1474E+00  7.4858E+00
             5.1711E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1549.72425364927        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      802
 NPARAMETR:  9.7995E-01  2.6464E-01  2.3180E-01  1.1277E+00  2.3868E-01  9.5598E-01  1.1693E+00  1.0000E-02  1.2967E+00  7.2160E-01
             2.8324E+00
 PARAMETER:  7.9741E-02 -1.2294E+00 -1.3619E+00  2.2021E-01 -1.3327E+00  5.4984E-02  2.5643E-01 -1.0732E+01  3.5985E-01 -2.2629E-01
             1.1411E+00
 GRADIENT:  -9.7157E-01  3.9664E-01 -4.5926E+00 -1.5473E+00  6.8948E+00 -3.3228E+00 -3.1847E-01  0.0000E+00  8.9765E-01  5.6687E-01
             4.8763E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1550.30606321751        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  9.7543E-01  2.0001E-01  2.6566E-01  1.1931E+00  2.5099E-01  9.5894E-01  1.8487E+00  1.0000E-02  1.2214E+00  7.2041E-01
             2.8555E+00
 PARAMETER:  7.5125E-02 -1.5094E+00 -1.2255E+00  2.7658E-01 -1.2824E+00  5.8077E-02  7.1448E-01 -1.1351E+01  2.9998E-01 -2.2793E-01
             1.1492E+00
 GRADIENT:  -2.0713E+00  1.7813E+00  1.0180E+00  2.8929E+00 -3.6542E+00 -4.5111E-02  2.7968E-01  0.0000E+00 -7.3355E-01  9.1853E-01
            -1.0087E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1550.81085756844        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1162
 NPARAMETR:  9.7182E-01  1.2516E-01  2.8537E-01  1.2326E+00  2.5609E-01  9.5679E-01  2.7268E+00  1.0000E-02  1.1904E+00  7.1711E-01
             2.8908E+00
 PARAMETER:  7.1411E-02 -1.9782E+00 -1.1540E+00  3.0916E-01 -1.2622E+00  5.5833E-02  1.1031E+00 -1.3300E+01  2.7428E-01 -2.3252E-01
             1.1615E+00
 GRADIENT:   2.1625E+00  7.1349E+00  1.0169E+00 -6.6018E+00 -3.5604E+00  8.3706E-01  6.5409E+00  0.0000E+00 -6.8805E-01 -2.8131E+00
             3.9233E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1551.86737898246        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1346             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7078E-01  6.6426E-02  3.0828E-01  1.2844E+00  2.6490E-01  9.3071E-01  4.8667E+00  1.0000E-02  1.1459E+00  6.9303E-01
             2.8864E+00
 PARAMETER:  7.0341E-02 -2.6117E+00 -1.0767E+00  3.5033E-01 -1.2284E+00  2.8189E-02  1.6824E+00 -1.5998E+01  2.3622E-01 -2.6668E-01
             1.1600E+00
 GRADIENT:   5.0109E+01 -5.8283E-01  4.6330E+01  6.1228E+01  4.2454E+01 -5.4953E+00 -2.4619E+00  0.0000E+00  5.9056E+00  7.1992E+00
             1.5853E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1552.31531992259        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1482
 NPARAMETR:  9.6547E-01  6.6201E-02  3.0712E-01  1.2832E+00  2.6622E-01  9.5816E-01  4.8600E+00  1.0000E-02  1.1641E+00  6.7591E-01
             2.8839E+00
 PARAMETER:  6.4859E-02 -2.6151E+00 -1.0805E+00  3.4932E-01 -1.2234E+00  5.7258E-02  1.6810E+00 -1.5998E+01  2.5197E-01 -2.9170E-01
             1.1591E+00
 GRADIENT:  -3.1467E-01 -3.0293E+00  1.8957E+01  5.7359E+00 -2.1587E+01  1.7275E+00 -8.1829E+00  0.0000E+00  4.3453E+00  3.9669E+00
             3.6984E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1552.34660285918        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1547
 NPARAMETR:  9.6536E-01  6.6061E-02  3.0717E-01  1.2834E+00  2.6625E-01  9.5739E-01  4.8544E+00  1.0000E-02  1.1642E+00  6.7408E-01
             2.8856E+00
 PARAMETER:  6.4727E-02 -2.6155E+00 -1.0810E+00  3.4925E-01 -1.2227E+00  5.7460E-02  1.6809E+00 -1.5998E+01  2.5183E-01 -2.9462E-01
             1.1590E+00
 GRADIENT:  -2.9048E-01  4.3508E+01 -9.8758E+01 -1.7656E+02  8.2561E+01  1.4133E+00  6.4858E+01  0.0000E+00 -2.4958E+02 -2.1315E+02
            -5.1440E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1547
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5901E-03  1.6899E-02  3.4525E-05 -1.4103E-02  6.5257E-03
 SE:             2.8727E-02  1.0113E-02  2.1854E-04  2.6836E-02  2.0685E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5586E-01  9.4719E-02  8.7447E-01  5.9922E-01  7.5240E-01

 ETASHRINKSD(%)  3.7613E+00  6.6120E+01  9.9268E+01  1.0096E+01  3.0703E+01
 ETASHRINKVR(%)  7.3812E+00  8.8522E+01  9.9995E+01  1.9173E+01  5.1979E+01
 EBVSHRINKSD(%)  2.8939E+00  7.8278E+01  9.9241E+01  7.1677E+00  2.8394E+01
 EBVSHRINKVR(%)  5.7041E+00  9.5282E+01  9.9994E+01  1.3822E+01  4.8726E+01
 RELATIVEINF(%)  9.2589E+01  2.8840E+00  3.3756E-04  3.2190E+01  2.9359E+00
 EPSSHRINKSD(%)  2.5785E+01
 EPSSHRINKVR(%)  4.4921E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1552.3466028591786     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -633.40806965450588     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.60
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1552.347       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.65E-01  6.62E-02  3.07E-01  1.28E+00  2.66E-01  9.58E-01  4.86E+00  1.00E-02  1.16E+00  6.74E-01  2.88E+00
 


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
+        1.28E+03
 
 TH 2
+       -1.56E+02  1.02E+05
 
 TH 3
+        1.05E+01  2.45E+03  3.46E+04
 
 TH 4
+        5.08E+00  3.92E+02  1.02E+04  8.33E+03
 
 TH 5
+        9.56E+01  2.39E+04 -2.35E+04 -1.35E+02  4.63E+04
 
 TH 6
+       -3.37E+01  5.21E+01 -1.55E+01 -2.50E+01  1.81E+01  2.21E+02
 
 TH 7
+       -1.11E+00  1.06E+03  4.35E+01  5.45E+00  5.31E+02  1.33E+00  4.74E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.97E+01  1.74E+02  1.62E+04 -7.41E+01  1.44E+02 -2.98E+01  4.43E+00  0.00E+00  1.86E+04
 
 TH10
+        4.60E+01  2.99E+02  2.37E+04 -1.06E+02  3.66E+02 -4.05E+01  7.46E+00  0.00E+00  2.73E+04  4.04E+04
 
 TH11
+       -1.78E+01  1.76E+01  1.37E+03 -1.25E+01 -1.38E+03  1.18E+00  4.26E-01  0.00E+00 -3.28E+00  7.50E+00  2.01E+02
 
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
 #CPUT: Total CPU Time in Seconds,       34.344
Stop Time:
Thu Sep 30 00:34:10 CDT 2021
