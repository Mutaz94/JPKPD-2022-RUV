Sat Sep 18 07:43:04 CDT 2021
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
$DATA ../../../../data/int/D/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   62004.3899512590        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0007E+03  7.0133E+02 -5.6144E+01  8.4467E+02 -2.2498E+02 -3.6705E+03 -2.8506E+03 -3.8569E+01 -2.7557E+03 -4.2851E+02
            -1.2067E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -401.100913279948        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0756E+00  1.9411E+00  9.3548E-01  1.5112E+00  1.0006E+00  4.8841E+00  4.7490E+00  9.7395E-01  1.9434E+00  1.1415E+00
             1.3146E+01
 PARAMETER:  1.7286E-01  7.6324E-01  3.3308E-02  5.1290E-01  1.0058E-01  1.6860E+00  1.6579E+00  7.3600E-02  7.6446E-01  2.3236E-01
             2.6761E+00
 GRADIENT:  -1.1047E+01  2.8729E+01 -4.0923E+01  8.9038E+01 -1.3979E+01  1.3158E+02  3.4819E+01  3.8535E+00  2.4216E+01  2.1594E+01
            -7.7180E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -446.942507637034        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.5083E-01  1.8282E+00  3.9700E+01  3.4961E+00  2.8256E+00  3.2614E+00  8.9668E+00  6.6831E-01  4.1395E+00  1.0089E+00
             1.2958E+01
 PARAMETER:  4.9580E-02  7.0334E-01  3.7813E+00  1.3516E+00  1.1387E+00  1.2822E+00  2.2935E+00 -3.0301E-01  1.5206E+00  1.0884E-01
             2.6617E+00
 GRADIENT:  -2.6317E+01  2.5424E+01 -7.5839E+00  6.3844E+01  3.7388E+01  8.4599E+01  3.2185E+01  1.0233E-02  4.1230E+01  1.3026E+01
            -8.8972E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -548.966731401153        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0922E+00  1.1692E+00  1.3705E+01  1.4384E+00  2.1981E+00  2.1683E+00  5.3557E+00  4.8079E-01  1.2431E+00  3.3227E-01
             1.3994E+01
 PARAMETER:  1.8817E-01  2.5633E-01  2.7178E+00  4.6355E-01  8.8760E-01  8.7393E-01  1.7782E+00 -6.3233E-01  3.1759E-01 -1.0018E+00
             2.7386E+00
 GRADIENT:  -2.3684E+01 -1.0943E+00 -2.0636E+00  2.8062E+00  9.5043E+00 -1.3123E+01 -3.2694E+00  3.5053E-02  7.6896E-01  1.5674E+00
             6.2064E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -551.949373827126        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.1275E+00  9.0131E-01  2.5331E+01  1.5536E+00  2.2192E+00  2.2587E+00  5.9458E+00  3.7373E-01  1.3782E+00  2.1547E-01
             1.3547E+01
 PARAMETER:  2.2002E-01 -3.9049E-03  3.3320E+00  5.4059E-01  8.9713E-01  9.1481E-01  1.8827E+00 -8.8422E-01  4.2075E-01 -1.4349E+00
             2.7062E+00
 GRADIENT:  -3.6295E+00 -1.8656E+00 -1.2134E-01 -1.6890E+00 -1.8428E-01  3.6633E-01  2.0873E+00  6.8297E-03  1.8408E+00  6.4113E-01
             1.3689E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -551.985442090770        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.1363E+00  1.0068E+00  2.2342E+01  1.4797E+00  2.2114E+00  2.2460E+00  5.6924E+00  4.5265E-01  1.3198E+00  1.7785E-01
             1.3501E+01
 PARAMETER:  2.2782E-01  1.0680E-01  3.2065E+00  4.9185E-01  8.9364E-01  9.0914E-01  1.8391E+00 -6.9263E-01  3.7746E-01 -1.6268E+00
             2.7028E+00
 GRADIENT:  -1.1364E-01 -1.5704E+00 -1.4160E-01 -2.3363E+00 -1.3837E-01 -1.4969E+00  1.1519E+00  1.1100E-02  1.4214E+00  4.3819E-01
             7.7058E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -551.988217010857        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  1.1353E+00  1.0526E+00  2.1715E+01  1.4524E+00  2.2143E+00  2.2488E+00  5.6018E+00  4.7732E-01  1.2933E+00  1.5930E-01
             1.3473E+01
 PARAMETER:  2.2692E-01  1.5122E-01  3.1780E+00  4.7322E-01  8.9494E-01  9.1038E-01  1.8231E+00 -6.3957E-01  3.5716E-01 -1.7369E+00
             2.7007E+00
 GRADIENT:  -3.9189E-02 -1.1270E+00 -1.4111E-01 -1.6675E+00 -2.7224E-02 -1.2108E+00  8.0348E-01  1.2364E-02  1.0249E+00  3.5182E-01
             4.3333E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -551.988444654414        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  1.1349E+00  1.0680E+00  2.1616E+01  1.4432E+00  2.2162E+00  2.2499E+00  5.5725E+00  4.8203E-01  1.2832E+00  1.5137E-01
             1.3463E+01
 PARAMETER:  2.2653E-01  1.6579E-01  3.1734E+00  4.6686E-01  8.9581E-01  9.1089E-01  1.8178E+00 -6.2975E-01  3.4933E-01 -1.7880E+00
             2.7000E+00
 GRADIENT:  -4.6666E-02 -9.7302E-01 -1.3511E-01 -1.4046E+00  5.9684E-03 -1.0845E+00  6.7482E-01  1.2490E-02  8.5590E-01  3.1767E-01
             3.1948E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -551.988610757531        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.1345E+00  1.0841E+00  2.1590E+01  1.4338E+00  2.2188E+00  2.2512E+00  5.5427E+00  4.8285E-01  1.2725E+00  1.4188E-01
             1.3454E+01
 PARAMETER:  2.2616E-01  1.8079E-01  3.1723E+00  4.6030E-01  8.9698E-01  9.1145E-01  1.8125E+00 -6.2806E-01  3.4102E-01 -1.8528E+00
             2.6993E+00
 GRADIENT:  -5.3616E-02 -8.0943E-01 -1.2472E-01 -1.1309E+00  3.2739E-02 -9.4196E-01  5.4060E-01  1.2318E-02  6.7935E-01  2.7900E-01
             2.1068E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -551.988729260845        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  1.1342E+00  1.0955E+00  2.1629E+01  1.4273E+00  2.2210E+00  2.2521E+00  5.5223E+00  4.8025E-01  1.2651E+00  1.3438E-01
             1.3449E+01
 PARAMETER:  2.2593E-01  1.9117E-01  3.1740E+00  4.5575E-01  8.9798E-01  9.1186E-01  1.8088E+00 -6.3344E-01  3.3517E-01 -1.9071E+00
             2.6989E+00
 GRADIENT:  -5.6642E-02 -6.9258E-01 -1.1473E-01 -9.4222E-01  4.6468E-02 -8.3653E-01  4.4753E-01  1.1976E-02  5.5772E-01  2.5027E-01
             1.4117E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -551.988926104361        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  1.1340E+00  1.1084E+00  2.1742E+01  1.4199E+00  2.2241E+00  2.2532E+00  5.4994E+00  4.7334E-01  1.2567E+00  1.2484E-01
             1.3443E+01
 PARAMETER:  2.2571E-01  2.0291E-01  3.1792E+00  4.5062E-01  8.9934E-01  9.1234E-01  1.8046E+00 -6.4793E-01  3.2846E-01 -1.9807E+00
             2.6985E+00
 GRADIENT:  -5.4396E-02 -5.5697E-01 -1.0035E-01 -7.2996E-01  5.7191E-02 -7.0984E-01  3.4145E-01  1.1331E-02  4.2067E-01  2.1591E-01
             6.8337E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -551.989105376505        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  1.1338E+00  1.1193E+00  2.1907E+01  1.4139E+00  2.2271E+00  2.2541E+00  5.4805E+00  4.6354E-01  1.2495E+00  1.1583E-01
             1.3439E+01
 PARAMETER:  2.2556E-01  2.1274E-01  3.1868E+00  4.4634E-01  9.0071E-01  9.1276E-01  1.8012E+00 -6.6886E-01  3.2277E-01 -2.0557E+00
             2.6982E+00
 GRADIENT:  -4.7404E-02 -4.4061E-01 -8.5548E-02 -5.5385E-01  6.2062E-02 -5.9664E-01  2.5189E-01  1.0557E-02  3.0634E-01  1.8579E-01
             1.2160E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -551.989323225565        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      872
 NPARAMETR:  1.1336E+00  1.1303E+00  2.2156E+01  1.4079E+00  2.2307E+00  2.2551E+00  5.4621E+00  4.4914E-01  1.2425E+00  1.0576E-01
             1.3435E+01
 PARAMETER:  2.2543E-01  2.2244E-01  3.1981E+00  4.4213E-01  9.0232E-01  9.1319E-01  1.7978E+00 -7.0041E-01  3.1710E-01 -2.1466E+00
             2.6979E+00
 GRADIENT:  -3.7929E-02 -3.2226E-01 -6.7996E-02 -3.8065E-01  6.1757E-02 -4.8058E-01  1.6236E-01  9.5579E-03  1.9373E-01  1.5484E-01
            -3.8424E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -551.989513624916        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  1.1335E+00  1.1405E+00  2.2493E+01  1.4025E+00  2.2347E+00  2.2560E+00  5.4451E+00  4.3033E-01  1.2359E+00  9.5115E-02
             1.3433E+01
 PARAMETER:  2.2535E-01  2.3146E-01  3.2132E+00  4.3825E-01  9.0412E-01  9.1361E-01  1.7947E+00 -7.4321E-01  3.1180E-01 -2.2527E+00
             2.6977E+00
 GRADIENT:  -2.5408E-02 -2.0903E-01 -4.8710E-02 -2.2083E-01  5.7470E-02 -3.6668E-01  7.8149E-02  8.3995E-03  8.9215E-02  1.2517E-01
            -8.0961E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -551.989880970514        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  1.1335E+00  1.1517E+00  2.3048E+01  1.3967E+00  2.2403E+00  2.2571E+00  5.4270E+00  4.0087E-01  1.2288E+00  8.1415E-02
             1.3430E+01
 PARAMETER:  2.2529E-01  2.4127E-01  3.2376E+00  4.3411E-01  9.0659E-01  9.1409E-01  1.7914E+00 -8.1412E-01  3.0602E-01 -2.4082E+00
             2.6975E+00
 GRADIENT:  -4.5840E-03 -8.1871E-02 -2.3109E-02 -4.9512E-02  4.5389E-02 -2.3535E-01 -1.5801E-02  6.8380E-03 -2.4891E-02  9.1652E-02
            -1.2123E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -551.990823299903        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.1335E+00  1.1621E+00  2.3895E+01  1.3917E+00  2.2472E+00  2.2582E+00  5.4110E+00  3.5933E-01  1.2224E+00  6.5576E-02
             1.3429E+01
 PARAMETER:  2.2530E-01  2.5025E-01  3.2737E+00  4.3049E-01  9.0970E-01  9.1455E-01  1.7884E+00 -9.2351E-01  3.0085E-01 -2.6245E+00
             2.6974E+00
 GRADIENT:   2.0706E-02  4.2089E-02  7.0486E-03  1.0619E-01  2.3955E-02 -1.0673E-01 -1.0448E-01  5.0354E-03 -1.2945E-01  5.9414E-02
            -1.4861E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -551.997729909754        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1156
 NPARAMETR:  1.1338E+00  1.1736E+00  2.6857E+01  1.3877E+00  2.2651E+00  2.2596E+00  5.3961E+00  2.4426E-01  1.2169E+00  3.3160E-02
             1.3433E+01
 PARAMETER:  2.2557E-01  2.6010E-01  3.3905E+00  4.2762E-01  9.1763E-01  9.1521E-01  1.7857E+00 -1.3095E+00  2.9632E-01 -3.3064E+00
             2.6977E+00
 GRADIENT:   8.5317E-02  2.0970E-01  7.1163E-02  2.6808E-01 -4.5993E-02  8.1419E-02 -2.1396E-01  1.8006E-03 -2.4333E-01  1.5170E-02
            -1.3606E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -552.375172631228        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1318
 NPARAMETR:  1.1336E+00  8.9525E-01  5.9269E+01  1.5767E+00  2.3275E+00  2.2770E+00  6.1557E+00  1.4879E-02  1.3421E+00  1.0000E-02
             1.3496E+01
 PARAMETER:  2.2544E-01 -1.0652E-02  4.1821E+00  5.5533E-01  9.4479E-01  9.2284E-01  1.9174E+00 -4.1078E+00  3.9426E-01 -7.0855E+00
             2.7024E+00
 GRADIENT:  -1.5099E+00  1.6339E-01  1.1287E-01  7.5916E-01 -4.0468E-01  7.6342E-01  4.6027E-02  1.9368E-06 -5.2391E-01  0.0000E+00
            -1.1674E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -552.418628122747        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1493
 NPARAMETR:  1.1394E+00  8.9860E-01  3.5384E+01  1.5715E+00  2.2776E+00  2.2702E+00  6.1262E+00  8.1857E-02  1.3491E+00  1.2709E-02
             1.3511E+01
 PARAMETER:  2.3047E-01 -6.9198E-03  3.6663E+00  5.5204E-01  9.2311E-01  9.1985E-01  1.9126E+00 -2.4028E+00  3.9947E-01 -4.2654E+00
             2.7035E+00
 GRADIENT:   3.5690E-02 -1.2975E-02  1.6132E-02 -5.0355E-03 -5.0438E-02 -6.1553E-02 -1.4659E-02  1.6552E-04 -1.0046E-02  2.1723E-03
             1.4171E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -552.418814164070        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1673
 NPARAMETR:  1.1388E+00  9.0457E-01  3.4728E+01  1.5691E+00  2.2764E+00  2.2719E+00  6.1127E+00  7.3806E-02  1.3481E+00  1.1070E-02
             1.3511E+01
 PARAMETER:  2.3000E-01 -2.9378E-04  3.6475E+00  5.5048E-01  9.2258E-01  9.2062E-01  1.9104E+00 -2.5063E+00  3.9870E-01 -4.4035E+00
             2.7035E+00
 GRADIENT:  -1.3551E-01  7.3532E-02  6.0820E-03  1.4937E-01  9.8848E-03  1.6590E-01 -4.0841E-02  1.3897E-04 -8.1216E-03  1.6489E-03
             1.2304E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -552.419320433461        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1857             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1392E+00  9.0027E-01  3.4172E+01  1.5702E+00  2.2736E+00  2.2707E+00  6.1229E+00  6.2717E-02  1.3486E+00  1.0000E-02
             1.3509E+01
 PARAMETER:  2.3031E-01 -5.0656E-03  3.6314E+00  5.5121E-01  9.2138E-01  9.2010E-01  1.9120E+00 -2.6691E+00  3.9904E-01 -4.6118E+00
             2.7034E+00
 GRADIENT:   7.8884E-01  7.0415E-02  3.3029E-03  1.0747E+00  2.2327E-01  1.8884E+00  6.3392E+00  1.0459E-04  1.1186E-01  0.0000E+00
             4.8594E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -552.419328493103        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2032
 NPARAMETR:  1.1392E+00  8.9985E-01  3.4201E+01  1.5705E+00  2.2738E+00  2.2706E+00  6.1226E+00  5.8335E-02  1.3488E+00  1.0000E-02
             1.3510E+01
 PARAMETER:  2.3034E-01 -5.5275E-03  3.6323E+00  5.5139E-01  9.2144E-01  9.2005E-01  1.9120E+00 -2.7416E+00  3.9924E-01 -4.6118E+00
             2.7034E+00
 GRADIENT:   7.9111E-03 -2.8964E-03 -2.3159E-04  1.0746E-03  7.3872E-03 -3.7988E-03 -4.8407E-04  9.0172E-05  3.3891E-03  0.0000E+00
             2.8148E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -552.419363790465        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2209
 NPARAMETR:  1.1392E+00  8.9963E-01  3.4003E+01  1.5704E+00  2.2729E+00  2.2706E+00  6.1228E+00  1.7276E-02  1.3486E+00  1.0000E-02
             1.3509E+01
 PARAMETER:  2.3033E-01 -5.7683E-03  3.6265E+00  5.5131E-01  9.2104E-01  9.2003E-01  1.9120E+00 -3.9585E+00  3.9906E-01 -4.6118E+00
             2.7034E+00
 GRADIENT:   1.0143E-02 -1.0482E-02 -1.5979E-03 -8.5391E-03 -8.1943E-03 -1.0348E-02 -1.0267E-03  8.0331E-06 -2.6681E-03  0.0000E+00
             8.8634E-03

0ITERATION NO.:  115    OBJECTIVE VALUE:  -552.419373381912        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2385
 NPARAMETR:  1.1392E+00  9.0005E-01  3.4175E+01  1.5704E+00  2.2737E+00  2.2707E+00  6.1225E+00  1.0789E-02  1.3487E+00  1.0000E-02
             1.3510E+01
 PARAMETER:  2.3032E-01 -5.3048E-03  3.6315E+00  5.5130E-01  9.2139E-01  9.2007E-01  1.9120E+00 -4.4292E+00  3.9912E-01 -4.6118E+00
             2.7034E+00
 GRADIENT:   5.6533E-03  1.6889E-03 -1.7376E-04  2.7948E-03  2.6779E-03  2.8713E-03  8.2214E-03  3.0883E-06  6.7793E-04  0.0000E+00
             6.6739E-03

0ITERATION NO.:  117    OBJECTIVE VALUE:  -552.419373478027        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     2451
 NPARAMETR:  1.1392E+00  8.9995E-01  3.4184E+01  1.5703E+00  2.2736E+00  2.2709E+00  6.1254E+00  1.0000E-02  1.3487E+00  1.0000E-02
             1.3509E+01
 PARAMETER:  2.3031E-01 -5.3140E-03  3.6316E+00  5.5130E-01  9.2139E-01  9.2007E-01  1.9120E+00 -4.5886E+00  3.9911E-01 -4.6118E+00
             2.7034E+00
 GRADIENT:  -1.0954E-03  5.8897E-03 -6.8558E-05  2.7590E-03  3.7905E-03 -9.7946E-03 -2.9866E-02  0.0000E+00 -2.5369E-04  0.0000E+00
             1.3591E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2451
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6064E-02  4.1483E-02 -8.4670E-07 -7.8462E-02 -3.6134E-06
 SE:             2.7954E-02  2.3134E-02  3.9566E-06  1.4595E-02  8.8421E-05
 N:                     100         100         100         100         100

 P VAL.:         5.6551E-01  7.2952E-02  8.3055E-01  7.6293E-08  9.6740E-01

 ETASHRINKSD(%)  6.3506E+00  2.2497E+01  9.9987E+01  5.1106E+01  9.9704E+01
 ETASHRINKVR(%)  1.2298E+01  3.9933E+01  1.0000E+02  7.6094E+01  9.9999E+01
 EBVSHRINKSD(%)  6.7040E+00  1.8780E+01  9.9981E+01  5.1070E+01  9.9608E+01
 EBVSHRINKVR(%)  1.2959E+01  3.4033E+01  1.0000E+02  7.6059E+01  9.9998E+01
 RELATIVEINF(%)  8.6858E+01  3.0658E+01  6.9619E-07  1.1215E+01  2.7780E-04
 EPSSHRINKSD(%)  3.7748E+00
 EPSSHRINKVR(%)  7.4071E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -552.41937347802741     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       1101.6699862903833     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    63.06
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    16.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -552.419       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  9.00E-01  3.42E+01  1.57E+00  2.27E+00  2.27E+00  6.12E+00  1.00E-02  1.35E+00  1.00E-02  1.35E+01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.44E+02
 
 TH 2
+       -1.00E+02  6.82E+01
 
 TH 3
+       -4.41E-02 -3.77E-02  3.28E-04
 
 TH 4
+       -1.68E+02  7.84E+01  6.48E-02  1.52E+02
 
 TH 5
+        1.88E+01  6.07E+00 -8.91E-02 -2.40E+01  2.45E+01
 
 TH 6
+       -4.29E+01  1.82E+01  4.71E-02  2.53E+01 -1.37E+01  2.63E+01
 
 TH 7
+        2.14E+00  2.39E+00 -1.50E-02 -4.65E+00  4.13E+00 -1.70E-01  9.77E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.38E+01 -2.25E+01 -3.33E-02 -5.08E+01  1.11E+01 -8.64E+00  2.21E+00  0.00E+00  1.75E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        9.24E-01 -3.34E+00 -1.81E-04 -6.01E+00  3.41E-01  2.26E-02  1.49E-01  0.00E+00  1.96E+00  0.00E+00  6.48E-01
 
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
+        1.44E+02
 
 TH 2
+        3.30E+00  4.25E+01
 
 TH 3
+        3.41E-03  1.78E-02  6.26E-04
 
 TH 4
+       -3.76E+00  1.99E+01  1.67E-02  7.76E+01
 
 TH 5
+       -1.46E+00 -5.79E+00 -1.13E-01 -8.20E+00  3.04E+01
 
 TH 6
+        1.27E+00  6.15E-02  1.99E-03  3.52E-01 -7.93E-01  3.02E+01
 
 TH 7
+        3.33E-01  3.82E+00 -3.02E-03 -5.58E+00  4.64E-01 -1.72E-01  2.72E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.14E-01 -6.34E+00 -1.04E-02 -2.35E+01  3.53E+00 -5.90E-01  1.85E+00  0.00E+00  1.98E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.85E+00 -2.50E+00 -1.74E-03 -6.07E+00  7.06E-01  5.53E-01  2.00E-01  0.00E+00  1.75E+00  0.00E+00  4.45E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.45E+02
 
 TH 2
+        3.99E+01  2.52E+01
 
 TH 3
+        2.47E-02  1.17E-02  2.41E-04
 
 TH 4
+        6.34E+01  2.35E+01  3.06E-02  8.11E+01
 
 TH 5
+       -1.32E+01 -2.46E+00 -6.07E-02 -1.27E+01  1.86E+01
 
 TH 6
+        1.94E+01  6.29E+00  6.63E-03 -1.12E+01 -3.52E+00  3.51E+01
 
 TH 7
+       -1.66E+00  3.78E+00 -2.91E-03 -5.98E+00  2.64E+00  2.05E+00  3.06E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -9.77E+00 -2.27E+00 -9.20E-03 -2.21E+01  4.43E+00  6.78E+00  1.69E+00  0.00E+00  1.50E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.96E+01 -3.56E+00 -1.47E-02 -1.27E+01  5.66E+00 -1.31E+00  1.92E+00  0.00E+00  3.26E+00  0.00E+00  8.51E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       79.708
Stop Time:
Sat Sep 18 07:44:25 CDT 2021
