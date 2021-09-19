Sat Sep 18 01:01:25 CDT 2021
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
$DATA ../../../../data/int/A2/dat72.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2906.89636430109        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.5260E+01  8.0295E+01  1.4071E+02 -8.7669E+01 -1.1618E+01  2.8933E+01 -3.7734E+01 -1.0476E+02  2.7016E+01 -3.6701E+01
            -1.7903E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3320.20867440250        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0589E+00  9.7116E-01  9.0904E-01  1.1026E+00  9.6568E-01  9.0118E-01  9.6349E-01  9.9960E-01  7.9763E-01  9.7630E-01
             1.5920E+00
 PARAMETER:  1.5725E-01  7.0741E-02  4.6315E-03  1.9770E-01  6.5075E-02 -4.0545E-03  6.2810E-02  9.9598E-02 -1.2611E-01  7.6018E-02
             5.6501E-01
 GRADIENT:   8.5033E+01  6.0987E+01  8.8761E+00  8.4586E+01 -8.6839E+00 -9.3456E+00 -5.5799E+00  1.0504E+00 -1.1616E+01 -1.1406E+00
            -1.3468E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3322.52627936739        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0606E+00  8.5228E-01  7.9270E-01  1.1568E+00  8.4978E-01  9.3055E-01  1.0258E+00  7.1102E-01  7.7947E-01  8.8608E-01
             1.5943E+00
 PARAMETER:  1.5886E-01 -5.9837E-02 -1.3231E-01  2.4564E-01 -6.2775E-02  2.8025E-02  1.2546E-01 -2.4105E-01 -1.4914E-01 -2.0948E-02
             5.6643E-01
 GRADIENT:   8.5516E+01  4.0401E+01 -1.9027E+00  9.5679E+01  1.8382E+01  2.9948E+00 -4.7881E+00 -3.9996E+00 -1.8414E+01 -8.3182E+00
            -1.2934E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3330.30811169838        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0239E+00  6.6543E-01  6.2983E-01  1.2062E+00  6.5469E-01  9.1749E-01  1.1539E+00  5.0315E-01  8.2236E-01  8.0866E-01
             1.6820E+00
 PARAMETER:  1.2358E-01 -3.0732E-01 -3.6231E-01  2.8746E-01 -3.2360E-01  1.3882E-02  2.4311E-01 -5.8687E-01 -9.5576E-02 -1.1238E-01
             6.1997E-01
 GRADIENT:  -1.4653E+01  1.7871E+01  1.3872E+01  3.5776E+01 -5.2450E+00 -1.6765E+00  2.3503E+00 -2.2935E-01  1.1993E+00  1.0867E+00
             1.7392E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3330.37430848768        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.0211E+00  6.1300E-01  5.8153E-01  1.2272E+00  6.0061E-01  9.1857E-01  1.2061E+00  4.5050E-01  8.2268E-01  7.7083E-01
             1.6870E+00
 PARAMETER:  1.2084E-01 -3.8939E-01 -4.4209E-01  3.0473E-01 -4.0980E-01  1.5063E-02  2.8740E-01 -6.9739E-01 -9.5190E-02 -1.6029E-01
             6.2297E-01
 GRADIENT:  -2.2793E+01  1.8552E+01  2.4853E+01  5.0836E+01 -1.2198E+01 -1.6889E+00  4.2756E+00 -5.6335E-01  1.2797E+00  9.4703E-01
             3.0142E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3330.37688681998        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:      393
 NPARAMETR:  1.0203E+00  6.0000E-01  5.6908E-01  1.2316E+00  5.8734E-01  9.1895E-01  1.2183E+00  4.3781E-01  8.2311E-01  7.6225E-01
             1.6883E+00
 PARAMETER:  1.2013E-01 -4.1082E-01 -4.6374E-01  3.0828E-01 -4.3215E-01  1.5480E-02  2.9742E-01 -7.2597E-01 -9.4671E-02 -1.7148E-01
             6.2374E-01
 GRADIENT:  -2.4889E+01  1.8423E+01  2.7340E+01  5.3766E+01 -1.3721E+01 -1.6584E+00  4.7340E+00 -6.4035E-01  1.3228E+00  9.0628E-01
             3.3572E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3330.37741422008        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      478
 NPARAMETR:  1.0200E+00  5.9381E-01  5.6308E-01  1.2336E+00  5.8102E-01  9.1914E-01  1.2239E+00  4.3186E-01  8.2335E-01  7.5828E-01
             1.6889E+00
 PARAMETER:  1.1980E-01 -4.2120E-01 -4.7433E-01  3.0990E-01 -4.4297E-01  1.5686E-02  3.0206E-01 -7.3965E-01 -9.4376E-02 -1.7670E-01
             6.2407E-01
 GRADIENT:  -2.5866E+01  1.8348E+01  2.8508E+01  5.5134E+01 -1.4441E+01 -1.6428E+00  4.9486E+00 -6.7683E-01  1.3415E+00  8.8374E-01
             3.5176E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3330.37750894189        NO. OF FUNC. EVALS.:  87
 CUMULATIVE NO. OF FUNC. EVALS.:      565
 NPARAMETR:  1.0198E+00  5.8998E-01  5.5937E-01  1.2348E+00  5.7712E-01  9.1926E-01  1.2274E+00  4.2821E-01  8.2351E-01  7.5587E-01
             1.6892E+00
 PARAMETER:  1.1959E-01 -4.2766E-01 -4.8095E-01  3.1087E-01 -4.4970E-01  1.5817E-02  3.0488E-01 -7.4814E-01 -9.4178E-02 -1.7988E-01
             6.2427E-01
 GRADIENT:  -2.6464E+01  1.8297E+01  2.9224E+01  5.5967E+01 -1.4883E+01 -1.6327E+00  5.0803E+00 -6.9923E-01  1.3527E+00  8.6937E-01
             3.6161E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3330.37755482654        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:      656
 NPARAMETR:  1.0197E+00  5.8754E-01  5.5699E-01  1.2355E+00  5.7464E-01  9.1934E-01  1.2296E+00  4.2589E-01  8.2362E-01  7.5436E-01
             1.6894E+00
 PARAMETER:  1.1946E-01 -4.3181E-01 -4.8522E-01  3.1148E-01 -4.5402E-01  1.5901E-02  3.0666E-01 -7.5358E-01 -9.4044E-02 -1.8189E-01
             6.2439E-01
 GRADIENT:  -2.6844E+01  1.8261E+01  2.9679E+01  5.6493E+01 -1.5164E+01 -1.6259E+00  5.1640E+00 -7.1347E-01  1.3598E+00  8.6013E-01
             3.6788E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3330.37759580528        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      749
 NPARAMETR:  1.0195E+00  5.8549E-01  5.5498E-01  1.2361E+00  5.7254E-01  9.1941E-01  1.2314E+00  4.2394E-01  8.2372E-01  7.5309E-01
             1.6896E+00
 PARAMETER:  1.1936E-01 -4.3531E-01 -4.8883E-01  3.1198E-01 -4.5767E-01  1.5972E-02  3.0815E-01 -7.5816E-01 -9.3928E-02 -1.8357E-01
             6.2449E-01
 GRADIENT:  -2.7162E+01  1.8231E+01  3.0060E+01  5.6932E+01 -1.5400E+01 -1.6202E+00  5.2342E+00 -7.2542E-01  1.3656E+00  8.5230E-01
             3.7313E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3330.59344281576        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      916
 NPARAMETR:  1.0276E+00  5.7846E-01  5.5743E-01  1.2342E+00  5.7487E-01  9.2575E-01  1.2459E+00  4.3248E-01  8.2055E-01  7.4097E-01
             1.6941E+00
 PARAMETER:  1.2722E-01 -4.4739E-01 -4.8441E-01  3.1038E-01 -4.5361E-01  2.2848E-02  3.1987E-01 -7.3821E-01 -9.7776E-02 -1.9980E-01
             6.2715E-01
 GRADIENT:  -2.3874E+01 -2.8669E+00  2.9569E+01  3.5549E+01 -1.6208E+01  3.3900E-01  5.1434E+00 -9.3090E-01  5.5566E-01 -5.9383E-01
             4.1979E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3340.68274211893        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1098
 NPARAMETR:  1.0417E+00  3.6031E-01  3.0741E-01  1.2612E+00  3.4315E-01  9.2068E-01  1.6337E+00  1.1637E+00  8.4212E-01  3.8717E-01
             1.5514E+00
 PARAMETER:  1.4085E-01 -9.2079E-01 -1.0796E+00  3.3203E-01 -9.6958E-01  1.7352E-02  5.9084E-01  2.5159E-01 -7.1833E-02 -8.4890E-01
             5.3915E-01
 GRADIENT:   5.6561E+00  2.2810E+01  2.4945E+01  1.6262E+02 -1.1220E+01 -7.7240E+00  2.4082E+01  1.9925E+01 -2.7156E+01  4.6302E+00
             7.1697E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3376.66443819139        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1275
 NPARAMETR:  1.0136E+00  2.2629E-01  1.6159E-01  1.0510E+00  2.2109E-01  9.5294E-01  1.5477E+00  1.2785E+00  1.1915E+00  1.8365E-01
             1.4032E+00
 PARAMETER:  1.1355E-01 -1.3859E+00 -1.7227E+00  1.4973E-01 -1.4092E+00  5.1796E-02  5.3675E-01  3.4567E-01  2.7521E-01 -1.5947E+00
             4.3876E-01
 GRADIENT:  -5.5213E+01 -8.6845E+00  4.2009E+00  4.0520E+00 -1.2356E+01  4.9768E+00 -1.4016E+01 -8.5089E+00  8.4557E+00 -4.1366E+00
             2.1648E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3379.18396248862        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1451
 NPARAMETR:  1.0302E+00  2.3468E-01  1.6647E-01  1.0737E+00  2.2636E-01  9.4461E-01  1.5945E+00  1.2949E+00  1.1302E+00  2.6986E-01
             1.3884E+00
 PARAMETER:  1.2978E-01 -1.3495E+00 -1.6929E+00  1.7110E-01 -1.3856E+00  4.3017E-02  5.6656E-01  3.5844E-01  2.2242E-01 -1.2098E+00
             4.2819E-01
 GRADIENT:  -1.6634E+01  3.5713E+00 -5.0827E+00  2.0157E+01  1.1344E+01  2.4927E+00  1.4194E+00  4.3539E+00  4.4781E+00 -2.4874E+00
            -2.6105E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3382.23803781388        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1627
 NPARAMETR:  1.0402E+00  1.9624E-01  1.2880E-01  9.7345E-01  1.8977E-01  9.3569E-01  1.5090E+00  1.1071E+00  1.1488E+00  5.6565E-01
             1.3888E+00
 PARAMETER:  1.3944E-01 -1.5284E+00 -1.9495E+00  7.3092E-02 -1.5619E+00  3.3528E-02  5.1146E-01  2.0173E-01  2.3871E-01 -4.6977E-01
             4.2841E-01
 GRADIENT:   6.7174E+00 -3.4780E+00  9.2157E+00  7.0993E-01 -1.1478E+01 -3.4121E-01 -2.8254E+00  1.2072E+00 -1.0478E+00 -1.2072E+00
            -2.2332E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3382.40273787411        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1802
 NPARAMETR:  1.0375E+00  1.8917E-01  1.1954E-01  9.4769E-01  1.8243E-01  9.3597E-01  1.5242E+00  1.0660E+00  1.1712E+00  6.0918E-01
             1.3897E+00
 PARAMETER:  1.3677E-01 -1.5651E+00 -2.0241E+00  4.6275E-02 -1.6014E+00  3.3827E-02  5.2148E-01  1.6388E-01  2.5801E-01 -3.9563E-01
             4.2910E-01
 GRADIENT:   1.1291E-02 -4.9663E-03 -2.7758E-02  2.7276E-01 -7.1750E-02 -1.6625E-02  5.0641E-03 -3.7536E-02 -2.5016E-02  2.1737E-02
            -8.7739E-02

0ITERATION NO.:   77    OBJECTIVE VALUE:  -3382.40278126310        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1859
 NPARAMETR:  1.0374E+00  1.8914E-01  1.1951E-01  9.4745E-01  1.8241E-01  9.3600E-01  1.5242E+00  1.0661E+00  1.1714E+00  6.0910E-01
             1.3898E+00
 PARAMETER:  1.3677E-01 -1.5653E+00 -2.0243E+00  4.6016E-02 -1.6015E+00  3.3859E-02  5.2148E-01  1.6402E-01  2.5819E-01 -3.9577E-01
             4.2914E-01
 GRADIENT:   1.9035E-03 -1.0526E-02 -2.9503E-02  2.0926E-02  4.8769E-02 -2.9525E-03  5.5178E-04 -1.3371E-02 -1.0875E-02 -2.3325E-03
             7.3513E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1859
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4535E-04  1.2649E-02  1.6828E-02  3.8957E-03  2.1136E-02
 SE:             2.9789E-02  2.7106E-02  2.2106E-02  2.7729E-02  2.2283E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9075E-01  6.4074E-01  4.4652E-01  8.8827E-01  3.4285E-01

 ETASHRINKSD(%)  2.0263E-01  9.1908E+00  2.5942E+01  7.1049E+00  2.5350E+01
 ETASHRINKVR(%)  4.0486E-01  1.7537E+01  4.5153E+01  1.3705E+01  4.4274E+01
 EBVSHRINKSD(%)  5.6395E-01  7.9146E+00  2.6403E+01  6.6461E+00  2.6165E+01
 EBVSHRINKVR(%)  1.1247E+00  1.5203E+01  4.5835E+01  1.2850E+01  4.5484E+01
 RELATIVEINF(%)  9.8868E+01  3.7904E+01  1.0267E+01  4.3490E+01  9.6345E+00
 EPSSHRINKSD(%)  2.3467E+01
 EPSSHRINKVR(%)  4.1427E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3382.4027812631048     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1728.3134214946940     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.23
 Elapsed covariance  time in seconds:    11.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3382.403       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.89E-01  1.20E-01  9.47E-01  1.82E-01  9.36E-01  1.52E+00  1.07E+00  1.17E+00  6.09E-01  1.39E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.95E-02  6.12E-02  6.61E-02  2.01E-01  5.84E-02  5.60E-02  1.20E-01  2.97E-01  1.89E-01  3.08E-01  4.67E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.68E-04
 
 TH 2
+        3.56E-04  3.75E-03
 
 TH 3
+        2.74E-04  3.99E-03  4.37E-03
 
 TH 4
+        7.57E-04  1.19E-02  1.31E-02  4.05E-02
 
 TH 5
+        3.21E-04  3.53E-03  3.85E-03  1.15E-02  3.41E-03
 
 TH 6
+       -2.52E-06  1.07E-03  1.10E-03  3.42E-03  1.03E-03  3.13E-03
 
 TH 7
+       -2.50E-05 -1.56E-03 -1.51E-03 -4.05E-03 -1.25E-03  4.32E-04  1.44E-02
 
 TH 8
+        3.01E-03  1.39E-02  1.43E-02  4.55E-02  1.32E-02  6.36E-03 -2.22E-03  8.79E-02
 
 TH 9
+        3.05E-04 -8.89E-03 -9.93E-03 -2.91E-02 -8.56E-03 -2.06E-03  3.21E-03 -2.00E-02  3.57E-02
 
 TH10
+       -2.08E-03 -1.79E-02 -1.93E-02 -5.90E-02 -1.73E-02 -4.87E-03  4.77E-03 -7.64E-02  4.07E-02  9.48E-02
 
 TH11
+        1.26E-04  1.80E-04  1.26E-04  3.61E-04  1.42E-04  3.50E-04 -5.13E-05  1.85E-03 -5.00E-04 -1.17E-04  2.18E-03
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.95E-02
 
 TH 2
+        1.97E-01  6.12E-02
 
 TH 3
+        1.41E-01  9.87E-01  6.61E-02
 
 TH 4
+        1.28E-01  9.70E-01  9.87E-01  2.01E-01
 
 TH 5
+        1.87E-01  9.88E-01  9.96E-01  9.82E-01  5.84E-02
 
 TH 6
+       -1.53E-03  3.13E-01  2.99E-01  3.03E-01  3.14E-01  5.60E-02
 
 TH 7
+       -7.07E-03 -2.12E-01 -1.90E-01 -1.68E-01 -1.79E-01  6.44E-02  1.20E-01
 
 TH 8
+        3.44E-01  7.68E-01  7.32E-01  7.63E-01  7.63E-01  3.83E-01 -6.26E-02  2.97E-01
 
 TH 9
+        5.49E-02 -7.69E-01 -7.96E-01 -7.66E-01 -7.76E-01 -1.94E-01  1.42E-01 -3.57E-01  1.89E-01
 
 TH10
+       -2.29E-01 -9.49E-01 -9.50E-01 -9.53E-01 -9.63E-01 -2.83E-01  1.29E-01 -8.36E-01  6.99E-01  3.08E-01
 
 TH11
+        9.14E-02  6.28E-02  4.08E-02  3.85E-02  5.20E-02  1.34E-01 -9.17E-03  1.34E-01 -5.67E-02 -8.13E-03  4.67E-02
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.92E+03
 
 TH 2
+       -1.41E+03  1.65E+04
 
 TH 3
+        6.30E+03 -2.28E+04  1.29E+05
 
 TH 4
+       -8.10E+01  1.71E+03 -9.81E+03  1.55E+03
 
 TH 5
+       -6.12E+03  5.60E+03 -1.03E+05  5.39E+03  1.08E+05
 
 TH 6
+        2.06E+02 -1.77E+02  1.38E+03 -4.71E+01 -1.71E+03  4.44E+02
 
 TH 7
+       -1.56E+01  2.42E+02  1.08E+01 -7.98E+00 -8.28E+01 -2.40E+01  8.06E+01
 
 TH 8
+        7.31E+00 -5.01E+02  1.42E+03 -1.52E+02 -6.97E+02 -4.80E+01 -4.70E+00  9.65E+01
 
 TH 9
+       -5.77E+01  1.81E+02  2.05E+02  4.10E+00 -9.43E+01  1.68E+01  9.51E+00 -4.95E+01  1.19E+02
 
 TH10
+       -6.31E+01  9.71E+00 -1.60E+03  1.44E+02  2.49E+03 -1.10E+02  1.44E+01  7.13E+01 -3.05E+01  2.93E+02
 
 TH11
+       -4.07E+00 -1.03E+02  9.01E+02 -3.20E+01 -1.13E+03  8.64E+00 -2.71E-01 -5.31E+01  4.74E+01 -1.24E+02  5.43E+02
 
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
 #CPUT: Total CPU Time in Seconds,       52.280
Stop Time:
Sat Sep 18 01:02:19 CDT 2021
