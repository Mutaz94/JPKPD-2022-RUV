Thu Sep 30 03:45:39 CDT 2021
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
$DATA ../../../../data/spa1/D/dat89.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   29751.1262998774        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.6239E+02  6.2000E+02 -2.0432E+01  6.7508E+02  3.2908E+02 -3.2314E+03 -1.2737E+03 -7.1161E+01 -1.7323E+03 -8.7889E+02
            -5.5809E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -446.050177763986        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2827E+00  1.0834E+00  8.9179E-01  1.5472E+00  1.3236E+00  2.0119E+00  1.1228E+00  9.5702E-01  1.0596E+00  9.8744E-01
             1.4676E+01
 PARAMETER:  3.4897E-01  1.8006E-01 -1.4524E-02  5.3643E-01  3.8038E-01  7.9907E-01  2.1578E-01  5.6073E-02  1.5786E-01  8.7357E-02
             2.7862E+00
 GRADIENT:   1.4443E+01  3.4933E+01 -2.1615E+00  5.0137E+01 -5.0722E+00  4.6600E+01 -6.5699E+00  4.9021E+00 -1.6894E+01  1.8654E+00
            -4.0641E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -463.375303571345        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2530E+00  8.7587E-01  1.1322E+00  1.6581E+00  6.1109E+00  1.5598E+00  6.9830E-01  4.4052E-01  1.3982E+00  8.5640E-01
             1.5311E+01
 PARAMETER:  3.2557E-01 -3.2543E-02  2.2413E-01  6.0567E-01  1.9101E+00  5.4454E-01 -2.5910E-01 -7.1981E-01  4.3522E-01 -5.5021E-02
             2.8286E+00
 GRADIENT:  -2.1934E+00  2.8187E+01  1.0934E+01  5.5162E+01 -2.5880E+00 -3.8969E+01 -2.7620E-01  7.6439E-02  8.7781E+00  3.0747E-02
             1.3832E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -543.213724706703        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  8.2842E-01  9.0222E-02  1.5620E-01  1.2045E+00  3.8674E+01  1.8019E+00  3.5453E+00  1.0000E-02  1.4841E-01  4.2202E+00
             1.4106E+01
 PARAMETER: -8.8233E-02 -2.3055E+00 -1.7566E+00  2.8604E-01  3.7552E+00  6.8885E-01  1.3656E+00 -1.7386E+01 -1.8078E+00  1.5399E+00
             2.7466E+00
 GRADIENT:  -6.8022E+00  4.1953E+01 -1.2458E+02  3.2942E+02 -1.1806E+00 -3.9710E+01  5.1065E+01  0.0000E+00 -1.2148E+00  1.0860E-01
            -9.5696E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -612.972659946887        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  5.9840E-01  3.0389E-02  4.3636E-02  5.5862E-01  3.5847E+02  1.7985E+00  1.4005E+00  1.0000E-02  3.3345E-02  3.0595E+00
             1.4766E+01
 PARAMETER: -4.1349E-01 -3.3937E+00 -3.0319E+00 -4.8229E-01  5.9819E+00  6.8697E-01  4.3683E-01 -2.8474E+01 -3.3008E+00  1.2182E+00
             2.7923E+00
 GRADIENT:   2.5505E+01  4.8273E+00 -1.3840E+02  2.7284E+02  7.0303E-03 -1.7465E+01  8.5250E+00  0.0000E+00  7.0566E-02  2.8587E-04
             5.4342E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -642.572954128175        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      436
 NPARAMETR:  5.2524E-01  1.6841E-02  2.6286E-02  3.3212E-01  7.5139E+02  1.9870E+00  1.4051E+00  1.0000E-02  1.8750E-02  2.5155E+00
             1.4196E+01
 PARAMETER: -5.4390E-01 -3.9840E+00 -3.5387E+00 -1.0023E+00  6.7219E+00  7.8665E-01  4.4011E-01 -3.3893E+01 -3.8766E+00  1.0225E+00
             2.7530E+00
 GRADIENT:   3.8392E+01  2.1385E+00 -4.4729E+01  4.1200E+01  1.9796E-03  2.8346E+01  8.8244E-01  0.0000E+00  2.6568E-02  3.9970E-06
             5.2865E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -652.218144348440        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      614
 NPARAMETR:  4.0253E-01  1.0654E-02  1.6980E-02  2.2659E-01  8.8037E+02  1.7811E+00  1.0575E+00  1.0000E-02  2.2583E-02  2.3217E+00
             1.2720E+01
 PARAMETER: -8.0998E-01 -4.4419E+00 -3.9757E+00 -1.3846E+00  6.8803E+00  6.7724E-01  1.5592E-01 -3.5767E+01 -3.6905E+00  9.4229E-01
             2.6432E+00
 GRADIENT:   1.0044E+01  3.4543E-01 -2.4566E+00 -3.9993E+00  7.4255E-04  4.8042E-01  6.8062E-02  0.0000E+00  3.7569E-02  1.4799E-07
            -5.1088E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -652.327639011304        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  3.9691E-01  1.0000E-02  1.6980E-02  2.2686E-01  6.4991E+02  1.7791E+00  6.0040E-01  1.0000E-02  1.0000E-02  2.2912E+00
             1.2752E+01
 PARAMETER: -8.2404E-01 -4.5320E+00 -3.9757E+00 -1.3834E+00  6.5768E+00  6.7612E-01 -4.1016E-01 -3.5767E+01 -4.6769E+00  9.2906E-01
             2.6457E+00
 GRADIENT:   6.9638E-01  1.0205E-02 -8.2853E-01 -1.0653E+00  8.7559E-04  4.2542E-03  5.7256E-03  0.0000E+00  0.0000E+00 -1.8467E-07
            -3.5421E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -652.337905730300        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:      997
 NPARAMETR:  3.9683E-01  1.0000E-02  1.6916E-02  2.2656E-01  4.2009E+02  1.7807E+00  3.3307E-02  1.0000E-02  1.0000E-02  2.2081E+00
             1.2753E+01
 PARAMETER: -8.2426E-01 -4.5547E+00 -3.9795E+00 -1.3847E+00  6.1405E+00  6.7698E-01 -3.3020E+00 -3.5767E+01 -4.6769E+00  8.9213E-01
             2.6458E+00
 GRADIENT:   1.4951E+00  0.0000E+00 -4.6261E+00  3.1210E+00  1.6445E-03  2.2691E-01  1.5965E-05  0.0000E+00  0.0000E+00 -4.5346E-07
            -6.7756E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -652.353499054865        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     1198             RESET HESSIAN, TYPE I
 NPARAMETR:  3.9626E-01  1.0000E-02  1.6763E-02  2.2525E-01  1.3128E+02  1.7820E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.5297E+00
             1.2751E+01
 PARAMETER: -8.2568E-01 -4.5547E+00 -3.9886E+00 -1.3906E+00  4.9774E+00  6.7772E-01 -5.9511E+00 -3.5767E+01 -4.6769E+00  1.0281E+00
             2.6456E+00
 GRADIENT:   5.5859E+01  0.0000E+00  7.4618E+01  3.0579E+01  5.9179E-03  1.2523E+01  0.0000E+00  0.0000E+00  0.0000E+00 -7.2330E-07
             2.5168E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -652.366000815943        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1373
 NPARAMETR:  3.9482E-01  1.0000E-02  1.6836E-02  2.2521E-01  3.0528E+01  1.7771E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.2916E+00
             1.2749E+01
 PARAMETER: -8.2932E-01 -4.5547E+00 -3.9842E+00 -1.3907E+00  3.5187E+00  6.7497E-01 -5.9511E+00 -3.5767E+01 -4.6769E+00  1.2914E+00
             2.6454E+00
 GRADIENT:  -1.4124E+00  0.0000E+00  1.4670E+00 -3.0710E+00  1.3928E-02 -2.9492E-01  0.0000E+00  0.0000E+00  0.0000E+00 -2.3437E-06
             4.3763E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -652.383249933809        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1552
 NPARAMETR:  3.9691E-01  1.0000E-02  1.6886E-02  2.2620E-01  1.3295E+01  1.7795E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7571E-02
             1.2744E+01
 PARAMETER: -8.2405E-01 -4.5547E+00 -3.9813E+00 -1.3863E+00  2.6874E+00  6.7631E-01 -5.9511E+00 -3.5767E+01 -4.6769E+00 -3.9415E+00
             2.6451E+00
 GRADIENT:   8.7115E-02  0.0000E+00 -1.6221E+00 -1.5923E-01  9.2023E-03 -5.3515E-02  0.0000E+00  0.0000E+00  0.0000E+00  2.8794E-07
            -3.9429E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -652.394771119536        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1750
 NPARAMETR:  3.9713E-01  1.0000E-02  1.6796E-02  2.2516E-01  1.3179E+01  1.7820E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0174E-02
             1.2763E+01
 PARAMETER: -8.2349E-01 -4.5547E+00 -3.9866E+00 -1.3910E+00  2.6786E+00  6.7774E-01 -5.9511E+00 -3.5767E+01 -4.6769E+00 -4.4879E+00
             2.6465E+00
 GRADIENT:   1.5282E+00  0.0000E+00 -1.4456E+00 -1.2815E+00  7.1274E-03  6.3075E-01  0.0000E+00  0.0000E+00  0.0000E+00  9.4823E-08
             5.2698E-01

0ITERATION NO.:   61    OBJECTIVE VALUE:  -652.394771119536        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1772
 NPARAMETR:  3.9713E-01  1.0000E-02  1.6796E-02  2.2516E-01  1.3179E+01  1.7820E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0174E-02
             1.2763E+01
 PARAMETER: -8.2349E-01 -4.5547E+00 -3.9866E+00 -1.3910E+00  2.6786E+00  6.7774E-01 -5.9511E+00 -3.5767E+01 -4.6769E+00 -4.4879E+00
             2.6465E+00
 GRADIENT:   1.5282E+00  0.0000E+00 -1.4456E+00 -1.2815E+00  7.1274E-03  6.3075E-01  0.0000E+00  0.0000E+00  0.0000E+00  9.4823E-08
             5.2698E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1772
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1378E-03  7.9011E-06  5.2551E-05 -1.6681E-04 -5.7912E-07
 SE:             2.9076E-02  4.1557E-06  2.8108E-04  3.2097E-04  1.1520E-06
 N:                     100         100         100         100         100

 P VAL.:         9.6878E-01  5.7267E-02  8.5169E-01  6.0327E-01  6.1518E-01

 ETASHRINKSD(%)  2.5924E+00  9.9986E+01  9.9058E+01  9.8925E+01  9.9996E+01
 ETASHRINKVR(%)  5.1176E+00  1.0000E+02  9.9991E+01  9.9988E+01  1.0000E+02
 EBVSHRINKSD(%)  2.6894E+00  9.9980E+01  9.8969E+01  9.8780E+01  9.9996E+01
 EBVSHRINKVR(%)  5.3064E+00  1.0000E+02  9.9989E+01  9.9985E+01  1.0000E+02
 RELATIVEINF(%)  4.8008E+00  2.4127E-07  4.7075E-05  5.8814E-05  6.0085E-09
 EPSSHRINKSD(%)  4.8284E+00
 EPSSHRINKVR(%)  9.4236E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -652.39477111953647     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       266.54376208513622     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.67
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -652.395       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.97E-01  1.00E-02  1.68E-02  2.25E-01  1.32E+01  1.78E+00  1.00E-02  1.00E-02  1.00E-02  1.02E-02  1.28E+01
 


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
+        2.08E+03
 
 TH 2
+        0.00E+00  5.01E+03
 
 TH 3
+       -1.95E+04  0.00E+00  3.05E+06
 
 TH 4
+        8.54E+00  0.00E+00 -2.61E+05  2.40E+04
 
 TH 5
+        2.29E-01  0.00E+00 -9.88E+00  8.13E-01  8.12E-04
 
 TH 6
+        2.21E+00  0.00E+00  8.02E+02 -1.01E+02 -1.38E-03  5.49E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        3.06E-03  0.00E+00 -2.18E-03  1.46E-03 -1.99E-05 -1.50E-03  0.00E+00  0.00E+00  0.00E+00 -3.46E-02
 
 TH11
+       -2.07E+01  0.00E+00  4.45E+02 -2.71E+01 -4.11E-03  8.67E-01  0.00E+00  0.00E+00  0.00E+00 -4.48E-05  2.72E+00
 
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
 #CPUT: Total CPU Time in Seconds,       36.324
Stop Time:
Thu Sep 30 03:46:17 CDT 2021
