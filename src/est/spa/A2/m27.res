Sat Sep 18 09:45:34 CDT 2021
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
$DATA ../../../../data/spa/A2/dat27.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -887.858676421124        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.5645E+01  1.7325E+01  1.3448E+01 -9.6802E+00  1.2444E+02  2.2923E+01 -2.2273E+01 -6.6402E+00 -7.7076E+01 -6.0000E+01
            -1.3602E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1368.50768478346        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0119E+00  9.5664E-01  1.0217E+00  1.0969E+00  9.0429E-01  8.3907E-01  9.5762E-01  9.4156E-01  1.1505E+00  8.3508E-01
             3.2571E+00
 PARAMETER:  1.1185E-01  5.5673E-02  1.2150E-01  1.9248E-01 -6.0367E-04 -7.5466E-02  5.6691E-02  3.9786E-02  2.4018E-01 -8.0222E-02
             1.2808E+00
 GRADIENT:   3.2616E+01  2.2568E+01  1.9146E+00  3.5521E+01 -2.1447E+01 -2.3592E+01  7.8582E+00  6.0676E+00  1.5165E+01  1.9699E+01
             9.7723E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1378.36284234519        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0151E+00  9.6255E-01  5.9016E-01  1.0591E+00  6.9496E-01  9.1117E-01  1.0055E+00  2.7670E-01  1.1481E+00  4.7647E-01
             3.1444E+00
 PARAMETER:  1.1495E-01  6.1827E-02 -4.2737E-01  1.5738E-01 -2.6390E-01  6.9709E-03  1.0545E-01 -1.1848E+00  2.3808E-01 -6.4136E-01
             1.2456E+00
 GRADIENT:   3.0774E+01  2.0279E+01 -2.1985E+01  4.3070E+01  8.1647E+00  1.1399E+00  1.1997E+01  9.4195E-01  1.0950E+01  9.0259E+00
             8.8346E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1389.58535555883        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.8830E-01  6.5073E-01  6.7768E-01  1.1994E+00  6.2641E-01  9.0882E-01  1.1364E+00  1.9783E-01  1.0012E+00  3.7690E-01
             2.7114E+00
 PARAMETER:  8.8234E-02 -3.2966E-01 -2.8908E-01  2.8179E-01 -3.6775E-01  4.3954E-03  2.2791E-01 -1.5204E+00  1.0117E-01 -8.7577E-01
             1.0975E+00
 GRADIENT:  -7.3658E+00  8.8739E+00  8.9036E+00  8.0134E+00 -1.3097E+01 -9.8861E-01 -3.9805E-02 -5.6237E-02 -1.7263E+00 -3.9008E-01
            -1.0628E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1390.04948911600        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.8954E-01  5.2698E-01  6.9394E-01  1.2731E+00  6.0301E-01  9.0745E-01  1.2505E+00  2.0098E-01  9.6301E-01  3.2178E-01
             2.7918E+00
 PARAMETER:  8.9486E-02 -5.4059E-01 -2.6537E-01  3.4146E-01 -4.0583E-01  2.8871E-03  3.2357E-01 -1.5046E+00  6.2313E-02 -1.0339E+00
             1.1267E+00
 GRADIENT:  -1.9156E+00  4.7676E+00  1.9045E+00  1.2362E+01 -4.5290E+00 -3.2972E-01  6.8046E-03 -8.4759E-02  1.6969E-01 -6.6724E-01
            -1.6135E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1390.72060435545        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  9.8692E-01  2.5610E-01  6.0612E-01  1.3989E+00  4.8451E-01  9.1229E-01  2.0637E+00  2.9992E-01  9.2064E-01  5.0302E-01
             2.6993E+00
 PARAMETER:  8.6834E-02 -1.2622E+00 -4.0067E-01  4.3572E-01 -6.2462E-01  8.2058E-03  8.2448E-01 -1.1042E+00  1.7319E-02 -5.8713E-01
             1.0930E+00
 GRADIENT:  -1.2217E-01  4.6423E+00  6.0970E+00  2.2692E+01 -1.1623E+01  2.6786E-01  1.1490E+00  6.6725E-01  6.2600E+00  1.2055E+00
             1.0538E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1391.02406344321        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  9.8549E-01  2.1632E-01  5.4954E-01  1.3860E+00  4.4571E-01  9.1337E-01  2.1696E+00  2.1286E-01  8.9535E-01  5.4653E-01
             2.6189E+00
 PARAMETER:  8.5385E-02 -1.4310E+00 -4.9867E-01  4.2641E-01 -7.0809E-01  9.3903E-03  8.7453E-01 -1.4471E+00 -1.0545E-02 -5.0416E-01
             1.0628E+00
 GRADIENT:  -8.7460E-01  2.1358E+00  1.4906E+00  5.0246E+00 -4.0035E+00 -2.4793E-01 -8.7008E-02  3.8090E-01 -1.3961E+00  8.9385E-01
            -1.3089E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1391.12936543012        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  9.8395E-01  1.5798E-01  5.1208E-01  1.3941E+00  4.1449E-01  9.1427E-01  2.4952E+00  1.1133E-01  8.8628E-01  5.6627E-01
             2.5940E+00
 PARAMETER:  8.3819E-02 -1.7453E+00 -5.6927E-01  4.3223E-01 -7.8071E-01  1.0369E-02  1.0143E+00 -2.0953E+00 -2.0724E-02 -4.6868E-01
             1.0532E+00
 GRADIENT:  -7.5485E-01  9.9610E-01  3.7063E+00  1.5235E+00 -6.9786E+00 -1.2031E-01 -8.5905E-01  1.0951E-01 -2.5830E+00  4.8257E-01
            -1.9677E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1391.61007840165        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      595
 NPARAMETR:  9.7981E-01  5.6661E-02  5.0848E-01  1.4378E+00  4.0018E-01  9.1249E-01  4.2320E+00  1.0000E-02  8.7593E-01  5.7229E-01
             2.6016E+00
 PARAMETER:  7.9605E-02 -2.7707E+00 -5.7633E-01  4.6312E-01 -8.1584E-01  8.4211E-03  1.5427E+00 -4.9862E+00 -3.2466E-02 -4.5811E-01
             1.0561E+00
 GRADIENT:  -2.5596E-02  9.5869E-02  6.2431E-01  4.6712E+00 -1.4702E+00  1.6987E-02 -3.1563E-01  0.0000E+00  1.2511E-01 -2.7215E-01
            -2.9987E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1391.92004341506        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  9.8024E-01  2.4803E-02  5.6591E-01  1.4855E+00  4.2861E-01  9.1088E-01  6.6299E+00  1.0000E-02  8.5778E-01  5.7736E-01
             2.6296E+00
 PARAMETER:  8.0042E-02 -3.5968E+00 -4.6932E-01  4.9577E-01 -7.4721E-01  6.6559E-03  1.9916E+00 -7.7289E+00 -5.3403E-02 -4.4929E-01
             1.0668E+00
 GRADIENT:   6.9851E-02  8.1381E-02  8.1199E-01  7.5830E+00 -2.9859E+00  4.5274E-02 -8.4665E-02  0.0000E+00 -4.2512E-01 -7.7491E-02
            -1.0652E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1391.95682033355        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  9.7936E-01  1.0098E-02  5.6644E-01  1.4883E+00  4.2732E-01  9.1024E-01  1.1566E+01  1.0000E-02  8.5585E-01  5.8122E-01
             2.6273E+00
 PARAMETER:  7.9142E-02 -4.4954E+00 -4.6838E-01  4.9764E-01 -7.5023E-01  5.9545E-03  2.5481E+00 -1.0818E+01 -5.5663E-02 -4.4262E-01
             1.0659E+00
 GRADIENT:  -4.7478E-02  4.5851E-02 -1.2738E-01  6.4497E-01  2.1588E-02 -7.6749E-03  5.6701E-02  0.0000E+00 -4.4546E-02  4.6858E-02
            -3.9905E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1391.95799473400        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1088            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7935E-01  1.0000E-02  5.6705E-01  1.4882E+00  4.2764E-01  9.1019E-01  1.1266E+01  1.0000E-02  8.5541E-01  5.8007E-01
             2.6284E+00
 PARAMETER:  7.9137E-02 -4.5074E+00 -4.6730E-01  4.9759E-01 -7.4948E-01  5.9015E-03  2.5218E+00 -1.0858E+01 -5.6175E-02 -4.4460E-01
             1.0664E+00
 GRADIENT:   5.2935E+00  1.5839E-02  5.4634E-01  1.1277E+01  2.1980E+00  3.8913E-01  2.1752E-02  0.0000E+00  1.1495E-01  2.5282E-02
             8.3994E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1391.95808388969        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  9.7935E-01  1.0000E-02  5.6704E-01  1.4882E+00  4.2763E-01  9.1020E-01  1.1141E+01  1.0000E-02  8.5583E-01  5.8032E-01
             2.6283E+00
 PARAMETER:  7.9133E-02 -4.5090E+00 -4.6733E-01  4.9756E-01 -7.4950E-01  5.9070E-03  2.5107E+00 -1.0858E+01 -5.5684E-02 -4.4418E-01
             1.0663E+00
 GRADIENT:  -8.3389E-03  1.2753E-03  4.3514E-03 -1.4323E-02  2.7748E-02 -1.5596E-03 -3.8642E-05  0.0000E+00  7.4782E-04 -1.8856E-03
             1.0446E-02

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1391.95808388969        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1246
 NPARAMETR:  9.7935E-01  1.0000E-02  5.6704E-01  1.4882E+00  4.2763E-01  9.1020E-01  1.1141E+01  1.0000E-02  8.5583E-01  5.8032E-01
             2.6283E+00
 PARAMETER:  7.9133E-02 -4.5090E+00 -4.6733E-01  4.9756E-01 -7.4950E-01  5.9070E-03  2.5107E+00 -1.0858E+01 -5.5684E-02 -4.4418E-01
             1.0663E+00
 GRADIENT:  -8.3389E-03  1.2753E-03  4.3514E-03 -1.4323E-02  2.7748E-02 -1.5596E-03 -3.8642E-05  0.0000E+00  7.4782E-04 -1.8856E-03
             1.0446E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1246
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.0326E-04  7.0094E-04  7.8824E-05 -1.0699E-02 -8.5105E-03
 SE:             2.8961E-02  1.8196E-03  2.2352E-04  2.7036E-02  1.8107E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9165E-01  7.0008E-01  7.2435E-01  6.9230E-01  6.3834E-01

 ETASHRINKSD(%)  2.9760E+00  9.3904E+01  9.9251E+01  9.4251E+00  3.9340E+01
 ETASHRINKVR(%)  5.8634E+00  9.9628E+01  9.9994E+01  1.7962E+01  6.3203E+01
 EBVSHRINKSD(%)  2.9099E+00  9.4310E+01  9.9235E+01  9.0566E+00  3.9261E+01
 EBVSHRINKVR(%)  5.7350E+00  9.9676E+01  9.9994E+01  1.7293E+01  6.3107E+01
 RELATIVEINF(%)  8.0776E+01  1.6203E-02  2.2368E-04  8.1773E+00  1.0188E+00
 EPSSHRINKSD(%)  3.0018E+01
 EPSSHRINKVR(%)  5.1025E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1391.9580838896868     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -656.80725732594863     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1391.958       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.00E-02  5.67E-01  1.49E+00  4.28E-01  9.10E-01  1.11E+01  1.00E-02  8.56E-01  5.80E-01  2.63E+00
 


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
+       -2.27E+01  4.67E+03
 
 TH 3
+       -1.96E+01  2.94E+01  1.98E+03
 
 TH 4
+       -4.94E+01  1.94E+01 -1.21E+02  5.66E+02
 
 TH 5
+        1.05E+02 -7.68E+01 -3.42E+03 -2.09E+02  6.43E+03
 
 TH 6
+       -2.68E+00 -2.46E+01  9.26E+00 -1.29E+01 -4.41E+00  2.08E+02
 
 TH 7
+       -4.53E-02  1.73E-01  6.71E-03 -1.47E-02  1.30E-02  3.12E-02  3.42E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.65E+01 -4.22E+01  3.22E+01 -1.30E+01  8.83E+00 -2.95E+00  5.14E-02  0.00E+00  2.03E+02
 
 TH10
+       -1.09E+01  5.14E+00 -3.99E+01 -5.84E+00  5.37E+01 -2.91E+00  2.47E-02  0.00E+00  1.79E+00  7.08E+01
 
 TH11
+       -1.59E+01 -3.00E-01 -7.58E+00 -6.11E+00 -1.25E+01  3.41E+00  2.00E-03  0.00E+00  9.50E+00  3.06E+01  4.46E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.410
Stop Time:
Sat Sep 18 09:45:55 CDT 2021
