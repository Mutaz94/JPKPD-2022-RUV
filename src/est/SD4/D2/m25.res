Sun Oct 24 04:32:29 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat25.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1520.82552092181        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9953E+02 -1.2433E+02 -3.4311E+01 -1.6410E+02  3.8470E+01 -3.4372E+01 -8.3491E+01  2.8750E+00 -1.2325E+02 -3.6862E+00
            -3.2944E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1542.73251733767        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      109
 NPARAMETR:  9.6539E-01  1.0748E+00  1.2188E+00  8.0911E-01  1.2622E+00  1.0454E+00  1.9601E+00  9.7103E-01  1.8462E+00  9.8034E-01
             1.1748E+00
 PARAMETER:  6.4775E-02  1.7210E-01  2.9785E-01 -1.1182E-01  3.3284E-01  1.4439E-01  7.7298E-01  7.0599E-02  7.1313E-01  8.0140E-02
             2.6108E-01
 GRADIENT:  -1.6437E+02 -1.1392E+02 -2.1032E+01 -9.6246E+01  4.7690E+01 -7.8685E+01  8.7644E+00  3.3038E+00  3.0490E+01 -9.9645E+00
             3.8187E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1568.87554149512        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      288
 NPARAMETR:  1.1590E+00  1.4004E+00  1.3303E+00  9.3733E-01  1.3106E+00  1.2296E+00  2.0685E+00  1.0426E+00  1.7399E+00  1.2280E+00
             1.0608E+00
 PARAMETER:  2.4757E-01  4.3679E-01  3.8543E-01  3.5284E-02  3.7050E-01  3.0671E-01  8.2682E-01  1.4169E-01  6.5380E-01  3.0541E-01
             1.5899E-01
 GRADIENT:   1.5513E+02  1.7902E+01  3.3666E+00  1.8921E+00 -3.1876E+00 -7.7328E+00  5.1676E+01 -5.1034E+00  3.5193E+01  4.9351E+00
            -1.4947E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1590.30456044194        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      465
 NPARAMETR:  1.0384E+00  1.2210E+00  1.6353E+00  1.0490E+00  1.3204E+00  1.2203E+00  1.5925E+00  1.5413E+00  1.3435E+00  1.2295E+00
             1.0764E+00
 PARAMETER:  1.3768E-01  2.9970E-01  5.9182E-01  1.4784E-01  3.7792E-01  2.9908E-01  5.6533E-01  5.3263E-01  3.9530E-01  3.0662E-01
             1.7359E-01
 GRADIENT:  -1.0578E+01  5.2955E+00  2.1738E+00  9.9300E-01 -5.1690E+00  5.6295E+00  1.0802E+00  1.4875E-01  7.9443E-01  8.7469E-01
             4.9376E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1590.77858373924        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.0507E+00  9.0897E-01  1.9294E+00  1.2716E+00  1.2793E+00  1.1837E+00  1.8362E+00  1.5989E+00  1.2042E+00  1.2319E+00
             1.0613E+00
 PARAMETER:  1.4943E-01  4.5574E-03  7.5723E-01  3.4025E-01  3.4633E-01  2.6861E-01  7.0767E-01  5.6930E-01  2.8585E-01  3.0856E-01
             1.5953E-01
 GRADIENT:   1.1513E+01  9.5613E+00  2.3038E+00  1.1641E+01 -2.7492E+00 -5.8808E+00 -7.9801E-01 -1.8734E+00 -3.9986E+00  1.3798E-01
            -1.3707E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1591.52483734463        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  1.0272E+00  5.8236E-01  2.3223E+00  1.4962E+00  1.2607E+00  1.1972E+00  2.3019E+00  1.8345E+00  1.1544E+00  1.2449E+00
             1.0670E+00
 PARAMETER:  1.2688E-01 -4.4066E-01  9.4256E-01  5.0295E-01  3.3169E-01  2.8002E-01  9.3375E-01  7.0677E-01  2.4354E-01  3.1902E-01
             1.6482E-01
 GRADIENT:  -1.9607E+01  8.5812E+00  1.7557E+00  1.4223E+01 -3.9363E+00 -1.2187E-01  2.7555E+00 -2.5010E-01  4.2158E+00  3.0617E-03
             1.2647E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1591.96950986991        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      996
 NPARAMETR:  1.0266E+00  3.6402E-01  2.6288E+00  1.6600E+00  1.2581E+00  1.2076E+00  2.8447E+00  1.9836E+00  1.0760E+00  1.2598E+00
             1.0621E+00
 PARAMETER:  1.2624E-01 -9.1054E-01  1.0665E+00  6.0685E-01  3.2964E-01  2.8860E-01  1.1454E+00  7.8491E-01  1.7323E-01  3.3092E-01
             1.6026E-01
 GRADIENT:  -1.6920E+01  8.8826E+00  4.0630E+00  3.0396E+01 -3.3309E+00  4.4923E+00  1.4213E+00 -3.4990E+00 -2.2038E+00 -1.5429E-01
            -1.2987E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1593.65782702416        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1175
 NPARAMETR:  1.0434E+00  1.2532E-01  3.0626E+00  1.7991E+00  1.2732E+00  1.1961E+00  4.0475E+00  2.3092E+00  1.0120E+00  1.2690E+00
             1.0611E+00
 PARAMETER:  1.4249E-01 -1.9769E+00  1.2193E+00  6.8731E-01  3.4156E-01  2.7904E-01  1.4981E+00  9.3689E-01  1.1191E-01  3.3820E-01
             1.5929E-01
 GRADIENT:   1.2196E+01  2.3205E+00  2.0796E+00  5.8547E+00 -5.2152E-01  2.2196E+00 -1.5982E-01 -1.4442E+00 -5.3947E+00 -7.9372E-01
            -2.2710E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1594.04650053905        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1353
 NPARAMETR:  1.0348E+00  3.8494E-02  3.1211E+00  1.8652E+00  1.2593E+00  1.1842E+00  5.5995E+00  2.3916E+00  1.0069E+00  1.2673E+00
             1.0651E+00
 PARAMETER:  1.3424E-01 -3.1572E+00  1.2382E+00  7.2335E-01  3.3059E-01  2.6909E-01  1.8227E+00  9.7197E-01  1.0684E-01  3.3691E-01
             1.6306E-01
 GRADIENT:   5.8904E-01  9.0333E-01  4.2094E-01  1.3412E+01 -3.8745E+00 -1.3548E+00  8.1487E-04  8.6354E-01 -2.7462E-01  4.4073E-01
            -2.4103E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1594.56887822971        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1536             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0359E+00  1.0000E-02  3.0947E+00  1.8391E+00  1.2535E+00  1.1921E+00  8.0162E+00  2.3633E+00  9.9269E-01  1.2645E+00
             1.0631E+00
 PARAMETER:  1.3525E-01 -4.6997E+00  1.2297E+00  7.0928E-01  3.2593E-01  2.7571E-01  2.1815E+00  9.6004E-01  9.2660E-02  3.3467E-01
             1.6120E-01
 GRADIENT:   4.9094E+02  0.0000E+00  7.8185E+00  1.1563E+03  1.2825E+01  1.3054E+02  3.2018E-01  4.7270E+00  1.6899E+01  2.8129E+00
             1.3527E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1594.60635731750        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1728             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0354E+00  1.0000E-02  3.0777E+00  1.8482E+00  1.2513E+00  1.1903E+00  7.2282E+00  2.3584E+00  9.9606E-01  1.2622E+00
             1.0630E+00
 PARAMETER:  1.3475E-01 -4.6997E+00  1.2242E+00  7.1422E-01  3.2416E-01  2.7423E-01  2.0780E+00  9.5799E-01  9.6049E-02  3.3289E-01
             1.6110E-01
 GRADIENT:   4.8813E+02  0.0000E+00  7.3525E+00  1.1751E+03  1.2494E+01  1.2932E+02  2.6051E-01  4.9069E+00  1.8374E+01  2.7295E+00
             1.2471E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1594.61444092129        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1922
 NPARAMETR:  1.0353E+00  1.0000E-02  3.0653E+00  1.8478E+00  1.2496E+00  1.1903E+00  1.2260E-01  2.3517E+00  9.9804E-01  1.2609E+00
             1.0630E+00
 PARAMETER:  1.3473E-01 -4.6997E+00  1.2201E+00  7.1402E-01  3.2279E-01  2.7422E-01 -1.9988E+00  9.5514E-01  9.8033E-02  3.3184E-01
             1.6107E-01
 GRADIENT:   2.3876E+00  0.0000E+00  2.4683E-01 -2.6243E+01  3.0027E-01  8.8388E-01  6.5386E-06  2.3808E-01 -2.2455E-01  1.0261E-01
            -2.0446E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1594.61826259071        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2116             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0353E+00  1.0000E-02  3.0425E+00  1.8474E+00  1.2466E+00  1.1903E+00  1.2286E-01  2.3396E+00  9.9906E-01  1.2583E+00
             1.0629E+00
 PARAMETER:  1.3470E-01 -4.6997E+00  1.2127E+00  7.1377E-01  3.2044E-01  2.7422E-01 -1.9967E+00  9.4996E-01  9.9059E-02  3.2976E-01
             1.6096E-01
 GRADIENT:   4.8789E+02  0.0000E+00  7.0556E+00  1.1734E+03  1.2743E+01  1.2929E+02  9.9488E-05  4.8678E+00  1.9261E+01  2.5691E+00
             1.3105E+00

0ITERATION NO.:   64    OBJECTIVE VALUE:  -1594.61927988762        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     2260
 NPARAMETR:  1.0353E+00  1.0000E-02  3.0153E+00  1.8470E+00  1.2454E+00  1.1903E+00  1.1068E-01  2.3318E+00  9.9922E-01  1.2543E+00
             1.0627E+00
 PARAMETER:  1.3468E-01 -4.6997E+00  1.2115E+00  7.1366E-01  3.1870E-01  2.7424E-01 -2.0906E+00  9.4699E-01  9.8774E-02  3.2986E-01
             1.6085E-01
 GRADIENT:  -6.7683E-03  0.0000E+00  4.7411E-01  9.2048E-02 -2.3180E-01  1.9343E-03  2.3330E-05  1.9839E-02 -7.0919E-02  1.8414E-01
            -1.2341E-05

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2260
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.6963E-04 -3.7851E-05 -4.2479E-02 -7.2138E-03 -6.1900E-02
 SE:             2.9873E-02  1.4462E-05  1.8260E-02  2.9549E-02  1.9309E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9280E-01  8.8613E-03  2.0003E-02  8.0713E-01  1.3475E-03

 ETASHRINKSD(%)  1.0000E-10  9.9952E+01  3.8825E+01  1.0082E+00  3.5311E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.2576E+01  2.0062E+00  5.8153E+01
 EBVSHRINKSD(%)  3.4802E-01  9.9956E+01  4.3457E+01  1.4739E+00  3.0726E+01
 EBVSHRINKVR(%)  6.9483E-01  1.0000E+02  6.8029E+01  2.9260E+00  5.2011E+01
 RELATIVEINF(%)  9.8227E+01  2.0262E-06  1.4636E+01  1.1352E+01  1.7207E+01
 EPSSHRINKSD(%)  4.5065E+01
 EPSSHRINKVR(%)  6.9822E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1594.6192798876159     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -859.46845332387772     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1594.619       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  3.04E+00  1.85E+00  1.24E+00  1.19E+00  1.12E-01  2.33E+00  9.99E-01  1.26E+00  1.06E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       73.814
Stop Time:
Sun Oct 24 04:32:43 CDT 2021
