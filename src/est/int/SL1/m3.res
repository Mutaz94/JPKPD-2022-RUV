Fri Sep 24 23:54:10 CDT 2021
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
$DATA ../../../../data/int/SL1/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3152.82272695461        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0312E+00 -3.8639E+01  6.2038E+01 -2.0569E+01  5.9910E+01  3.7125E+00 -5.1476E+01 -1.0536E+02 -5.9627E+01 -1.1776E+01
            -1.3107E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3403.96584424584        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0258E+00  1.0580E+00  9.8330E-01  9.7954E-01  1.0049E+00  9.9670E-01  1.0904E+00  1.0816E+00  1.1031E+00  9.8442E-01
             1.4217E+00
 PARAMETER:  1.2548E-01  1.5642E-01  8.3161E-02  7.9331E-02  1.0488E-01  9.6692E-02  1.8657E-01  1.7843E-01  1.9810E-01  8.4298E-02
             4.5183E-01
 GRADIENT:   3.7523E+01 -1.9540E+01  5.5908E+00 -2.0900E+01  8.6766E+00  3.1803E-01 -9.7133E+00 -1.6739E+01 -7.3116E+00  1.9890E+00
            -1.6479E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3406.00767407121        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0119E+00  1.0886E+00  1.0895E+00  9.8485E-01  1.0593E+00  9.9433E-01  1.1333E+00  1.7016E+00  1.0935E+00  9.3236E-01
             1.3940E+00
 PARAMETER:  1.1187E-01  1.8491E-01  1.8572E-01  8.4732E-02  1.5758E-01  9.4313E-02  2.2514E-01  6.3157E-01  1.8942E-01  2.9960E-02
             4.3214E-01
 GRADIENT:   5.2992E+00 -4.5010E+00 -7.5730E+00  2.6376E+00  1.0580E+01 -9.4291E-01 -2.5615E+00 -3.4409E+00 -7.0344E+00 -1.7826E+01
            -1.5920E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3406.51176665430        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      324
 NPARAMETR:  1.0112E+00  1.0992E+00  1.1007E+00  9.8114E-01  1.0705E+00  1.0060E+00  1.1285E+00  1.7694E+00  1.0878E+00  9.5474E-01
             1.3923E+00
 PARAMETER:  1.1111E-01  1.9459E-01  1.9597E-01  8.0955E-02  1.6816E-01  1.0602E-01  2.2090E-01  6.7065E-01  1.8419E-01  5.3683E-02
             4.3096E-01
 GRADIENT:  -2.2906E+01 -1.1566E+01 -9.1345E+00 -2.8201E+00  6.0110E+00 -1.9793E-01 -1.8110E+00 -4.0997E+00 -7.9252E+00 -1.7556E+01
            -1.5819E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3406.51963764171        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      504
 NPARAMETR:  1.0112E+00  1.0992E+00  1.1007E+00  9.8113E-01  1.0705E+00  1.0008E+00  1.1285E+00  1.7694E+00  1.0878E+00  9.5474E-01
             1.3924E+00
 PARAMETER:  1.1112E-01  1.9458E-01  1.9597E-01  8.0954E-02  1.6816E-01  1.0081E-01  2.2090E-01  6.7065E-01  1.8419E-01  5.3686E-02
             4.3105E-01
 GRADIENT:  -2.3154E+01 -1.1564E+01 -9.1361E+00 -2.8079E+00  6.0082E+00 -2.2538E+00 -1.8081E+00 -4.0893E+00 -7.9169E+00 -1.7550E+01
            -1.5801E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3406.84747814297        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      686
 NPARAMETR:  1.0113E+00  1.0992E+00  1.1009E+00  9.8109E-01  1.0704E+00  9.8393E-01  1.1286E+00  1.7692E+00  1.0880E+00  9.5482E-01
             1.3962E+00
 PARAMETER:  1.1122E-01  1.9454E-01  1.9613E-01  8.0911E-02  1.6805E-01  8.3799E-02  2.2102E-01  6.7053E-01  1.8433E-01  5.3766E-02
             4.3375E-01
 GRADIENT:  -2.3824E+01 -1.1782E+01 -9.0846E+00 -2.8915E+00  5.7549E+00 -9.0929E+00 -1.6856E+00 -3.8152E+00 -7.7299E+00 -1.7376E+01
            -1.5189E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3411.74620600695        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      866
 NPARAMETR:  1.0145E+00  1.0977E+00  1.1061E+00  9.7985E-01  1.0669E+00  9.5728E-01  1.1326E+00  1.7627E+00  1.0924E+00  9.5710E-01
             1.5132E+00
 PARAMETER:  1.1439E-01  1.9319E-01  2.0084E-01  7.9640E-02  1.6477E-01  5.6345E-02  2.2453E-01  6.6686E-01  1.8840E-01  5.6152E-02
             5.1425E-01
 GRADIENT:  -1.9916E+01 -1.8687E+01 -6.7496E+00 -7.1797E+00 -1.6594E+00 -2.0057E+01  1.7193E+00  3.9011E+00 -2.5355E+00 -1.2414E+01
             1.9778E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3413.69262119900        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1021
 NPARAMETR:  1.0120E+00  1.1197E+00  1.1045E+00  9.7583E-01  1.0834E+00  9.9577E-01  1.0547E+00  1.7714E+00  1.0927E+00  1.0861E+00
             1.4981E+00
 PARAMETER:  1.1194E-01  2.1306E-01  1.9938E-01  7.5538E-02  1.8008E-01  9.5761E-02  1.5325E-01  6.7178E-01  1.8863E-01  1.8256E-01
             5.0423E-01
 GRADIENT:   1.8259E-02 -1.4939E+01 -6.2543E+00  7.2797E+00 -1.5155E+00 -6.2562E-01  5.4543E-01  1.1250E+00  5.9067E-01 -8.8873E-01
             2.3661E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3415.97982203092        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     1099
 NPARAMETR:  1.0102E+00  1.2093E+00  1.1677E+00  9.3158E-01  1.1505E+00  9.9760E-01  9.5938E-01  1.9477E+00  1.0907E+00  1.1913E+00
             1.4966E+00
 PARAMETER:  1.1015E-01  2.9005E-01  2.5503E-01  2.9126E-02  2.4016E-01  9.7594E-02  5.8537E-02  7.6667E-01  1.8678E-01  2.7507E-01
             5.0323E-01
 GRADIENT:  -4.3208E+00  3.5898E+00 -1.9907E+00  2.2277E+00 -1.2174E+01 -8.9878E-03 -1.3478E+00 -1.5388E+00 -1.1809E+00  1.4622E+00
             1.3059E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3415.99602240909        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1243
 NPARAMETR:  1.0103E+00  1.2091E+00  1.1678E+00  9.3153E-01  1.1506E+00  1.0061E+00  9.5934E-01  1.9470E+00  1.0906E+00  1.1912E+00
             1.4970E+00
 PARAMETER:  1.1021E-01  2.8991E-01  2.5516E-01  2.9076E-02  2.4028E-01  1.0608E-01  5.8487E-02  7.6629E-01  1.8669E-01  2.7493E-01
             5.0348E-01
 GRADIENT:  -2.7033E+01 -9.6780E+00 -2.5986E+00 -2.0040E+00 -1.6937E+01 -8.8265E-02 -1.7458E+00 -2.1496E+00 -2.0464E+00  7.4713E-01
             9.2050E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3416.53250737092        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1386
 NPARAMETR:  1.0223E+00  1.2186E+00  1.1827E+00  9.3228E-01  1.1695E+00  1.0056E+00  9.7416E-01  1.9927E+00  1.1004E+00  1.1861E+00
             1.4962E+00
 PARAMETER:  1.2206E-01  2.9772E-01  2.6782E-01  2.9876E-02  2.5658E-01  1.0556E-01  7.3819E-02  7.8950E-01  1.9570E-01  2.7071E-01
             5.0295E-01
 GRADIENT:  -1.3514E+00 -6.8822E+00 -4.5287E+00  7.3381E+00 -7.3220E+00  4.3710E-02  4.9860E-01 -1.8270E+00  7.0514E-01 -6.4749E-01
             7.5256E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -3416.53250737092        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1416
 NPARAMETR:  1.0223E+00  1.2185E+00  1.1828E+00  9.3225E-01  1.1694E+00  1.0055E+00  9.7419E-01  1.9927E+00  1.1005E+00  1.1861E+00
             1.4965E+00
 PARAMETER:  1.2206E-01  2.9772E-01  2.6782E-01  2.9876E-02  2.5658E-01  1.0556E-01  7.3819E-02  7.8950E-01  1.9570E-01  2.7071E-01
             5.0295E-01
 GRADIENT:  -3.8377E+05  1.4731E+05 -1.7490E+05  2.1930E+05  1.9418E+05  3.2157E-02 -4.9826E+05 -7.7178E+03 -4.7874E+05  2.2038E+04
            -8.7212E+04
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.8         3.3         4.5         3.3         4.5
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1416
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3919E-03 -3.2114E-02 -2.6466E-02  1.7240E-02 -3.4547E-02
 SE:             2.9812E-02  2.1251E-02  2.0270E-02  2.5946E-02  2.4053E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6276E-01  1.3074E-01  1.9167E-01  5.0640E-01  1.5093E-01

 ETASHRINKSD(%)  1.2648E-01  2.8808E+01  3.2092E+01  1.3078E+01  1.9418E+01
 ETASHRINKVR(%)  2.5280E-01  4.9316E+01  5.3885E+01  2.4445E+01  3.5065E+01
 EBVSHRINKSD(%)  5.4694E-01  2.9223E+01  3.4002E+01  1.4430E+01  1.8417E+01
 EBVSHRINKVR(%)  1.0909E+00  4.9906E+01  5.6443E+01  2.6778E+01  3.3442E+01
 RELATIVEINF(%)  9.8903E+01  1.8361E+01  3.5486E+01  3.3271E+01  2.8502E+01
 EPSSHRINKSD(%)  2.0618E+01
 EPSSHRINKVR(%)  3.6986E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3416.5325073709214     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1762.4431476025106     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.77
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3416.533       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.22E+00  1.18E+00  9.32E-01  1.17E+00  1.01E+00  9.74E-01  1.99E+00  1.10E+00  1.19E+00  1.50E+00
 


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
+        1.50E+09
 
 TH 2
+        9.48E+01  1.78E+08
 
 TH 3
+        2.77E+08  1.83E+04  2.33E+08
 
 TH 4
+        1.07E+09  3.24E+08  2.27E+04  1.26E+09
 
 TH 5
+        2.93E+08 -1.14E+08  1.15E+08  3.91E+03  2.60E+08
 
 TH 6
+       -2.35E+03  3.42E+02 -9.23E+02  1.33E+03  5.58E+02  1.91E+02
 
 TH 7
+       -1.67E+04  3.58E+03 -7.59E+08 -1.37E+09  2.63E+03 -1.73E+03  1.31E+09
 
 TH 8
+        2.60E+02 -2.18E+07  6.76E+04 -4.57E+03 -6.71E+04  7.94E+01  2.07E+05  4.42E+06
 
 TH 9
+       -4.64E+08 -2.12E+03 -3.43E+08 -6.20E+08 -4.26E+03 -1.36E+03  5.23E+08  9.28E+04  2.36E+08
 
 TH10
+        3.11E+08 -8.46E+01  1.47E+07  4.16E+08 -4.46E+03  1.35E+02  3.51E+08 -6.25E+04  1.59E+08  1.06E+08
 
 TH11
+       -3.05E+02  4.56E+07 -8.83E+03  9.56E+03 -5.52E+07 -1.64E+02  1.70E+08  9.27E+06  7.69E+07 -5.15E+07  1.94E+07
 
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
 #CPUT: Total CPU Time in Seconds,       56.369
Stop Time:
Fri Sep 24 23:55:15 CDT 2021
