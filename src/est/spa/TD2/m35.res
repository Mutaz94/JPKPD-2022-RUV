Wed Sep 29 18:56:55 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1699.98639806319        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4675E+02  7.1852E+00 -4.5225E+01  8.4729E+01  8.3311E+01  6.4463E+01 -8.5093E-01  1.1490E+01  5.4455E+00 -1.7208E+01
             1.1005E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1707.12233477545        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0557E+00  1.0167E+00  1.2165E+00  1.0362E+00  1.0652E+00  9.0853E-01  1.0455E+00  7.7663E-01  1.0283E+00  1.2813E+00
             9.2770E-01
 PARAMETER:  1.5422E-01  1.1652E-01  2.9599E-01  1.3555E-01  1.6320E-01  4.0701E-03  1.4450E-01 -1.5279E-01  1.2792E-01  3.4786E-01
             2.4955E-02
 GRADIENT:   3.5432E+00  2.9608E+01 -2.4906E+00  4.4321E+01  2.0819E+00 -2.7627E+01  3.6178E+00 -1.5427E-01  3.4733E+00  5.8179E+00
            -1.9599E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1708.85398403133        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0627E+00  9.1107E-01  1.0198E+00  1.0722E+00  9.3383E-01  9.4522E-01  1.0625E+00  4.2640E-01  1.0267E+00  1.0879E+00
             9.6199E-01
 PARAMETER:  1.6079E-01  6.8616E-03  1.1960E-01  1.6967E-01  3.1539E-02  4.3658E-02  1.6058E-01 -7.5238E-01  1.2632E-01  1.8422E-01
             6.1246E-02
 GRADIENT:   1.7545E+01  4.1149E+00 -3.3660E+00  2.1913E+01  8.5338E+00 -1.0870E+01 -6.8509E-01 -7.6003E-01  1.3890E+01 -4.9404E+00
            -6.3064E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1710.37073403445        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  1.0522E+00  7.4613E-01  9.7890E-01  1.1660E+00  8.3768E-01  9.7347E-01  1.5081E+00  3.9289E-01  8.1846E-01  1.0042E+00
             9.7996E-01
 PARAMETER:  1.5086E-01 -1.9285E-01  7.8674E-02  2.5360E-01 -7.7122E-02  7.3109E-02  5.1088E-01 -8.3421E-01 -1.0033E-01  1.0415E-01
             7.9754E-02
 GRADIENT:  -5.6890E+00  1.4815E+01  8.8347E+00  1.4899E+01 -1.2976E+01  1.5775E+00  1.2758E-01 -1.6877E-01 -4.1620E+00 -2.4537E+00
             5.8539E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1711.05463744151        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0528E+00  5.3899E-01  9.6473E-01  1.2776E+00  7.6905E-01  9.6621E-01  1.9142E+00  2.5437E-01  7.6926E-01  9.9296E-01
             9.8088E-01
 PARAMETER:  1.5149E-01 -5.1805E-01  6.4093E-02  3.4496E-01 -1.6260E-01  6.5625E-02  7.4928E-01 -1.2690E+00 -1.6232E-01  9.2939E-02
             8.0691E-02
 GRADIENT:   5.8106E-01  5.6613E+00  3.6996E+00  9.5601E+00 -6.6393E+00 -4.2714E-01  8.3200E-02 -3.4818E-01 -6.4671E-01 -3.1231E-01
            -1.2925E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1711.14712463568        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0525E+00  4.5070E-01  9.8333E-01  1.3186E+00  7.5975E-01  9.6614E-01  2.1401E+00  2.2620E-01  7.5108E-01  1.0086E+00
             9.8688E-01
 PARAMETER:  1.5119E-01 -6.9696E-01  8.3188E-02  3.7658E-01 -1.7477E-01  6.5552E-02  8.6083E-01 -1.3863E+00 -1.8624E-01  1.0853E-01
             8.6793E-02
 GRADIENT:   3.2397E+00 -1.6262E+00 -2.3551E+00 -9.6050E+00  4.4102E+00  1.8628E-01 -5.4094E-01 -3.2313E-01  3.9527E-01  2.9190E-01
             3.9008E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1711.19756124832        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  1.0470E+00  4.1570E-01  9.9706E-01  1.3523E+00  7.4923E-01  9.6389E-01  2.2552E+00  3.3493E-01  7.4209E-01  1.0016E+00
             9.8376E-01
 PARAMETER:  1.4588E-01 -7.7780E-01  9.7052E-02  4.0183E-01 -1.8871E-01  6.3224E-02  9.1324E-01 -9.9384E-01 -1.9828E-01  1.0158E-01
             8.3623E-02
 GRADIENT:  -8.0291E+00  4.4772E+00  1.1002E+00  1.6005E+01 -6.6654E+00 -5.2938E-01 -9.4486E-01 -5.1978E-01 -4.5468E-01 -7.4505E-02
            -2.1658E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1711.23116484919        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1260
 NPARAMETR:  1.0449E+00  3.7771E-01  1.0175E+00  1.3788E+00  7.4507E-01  9.6285E-01  2.3972E+00  4.7787E-01  7.3418E-01  9.9035E-01
             9.7977E-01
 PARAMETER:  1.4389E-01 -8.7363E-01  1.1738E-01  4.2124E-01 -1.9428E-01  6.2142E-02  9.7429E-01 -6.3843E-01 -2.0900E-01  9.0299E-02
             7.9561E-02
 GRADIENT:  -1.0875E+01  5.5088E+00 -9.7670E-01  2.2264E+01 -9.9973E+00 -6.8809E-01 -9.1714E-01  2.1572E-01 -6.5665E-01 -6.6853E-02
            -2.7267E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1711.65713497486        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1440
 NPARAMETR:  1.0450E+00  2.2447E-01  1.1783E+00  1.4777E+00  7.7243E-01  9.6148E-01  3.2201E+00  7.3440E-01  7.0784E-01  1.0256E+00
             9.7600E-01
 PARAMETER:  1.4405E-01 -1.3940E+00  2.6405E-01  4.9049E-01 -1.5822E-01  6.0723E-02  1.2694E+00 -2.0870E-01 -2.4554E-01  1.2528E-01
             7.5705E-02
 GRADIENT:  -1.3706E+00  3.4760E+00 -7.4054E-01  1.1050E+01 -1.1166E+01  8.9555E-02  1.2363E+00  1.5808E+00 -1.5112E+00 -1.7546E+00
            -8.4985E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1711.99453257500        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1616
 NPARAMETR:  1.0405E+00  1.3525E-01  1.3729E+00  1.5500E+00  8.2792E-01  9.5871E-01  4.0653E+00  8.3315E-01  7.0019E-01  1.1208E+00
             9.7598E-01
 PARAMETER:  1.3969E-01 -1.9007E+00  4.1695E-01  5.3824E-01 -8.8833E-02  5.7835E-02  1.5025E+00 -8.2542E-02 -2.5640E-01  2.1407E-01
             7.5682E-02
 GRADIENT:  -7.2846E+00  1.8581E+00  1.3202E+00  1.9924E+01 -6.9080E+00 -4.5867E-01  1.4539E-01 -1.2097E+00  8.7446E-01 -8.7655E-02
            -4.4569E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1712.10142966812        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1795
 NPARAMETR:  1.0406E+00  9.8552E-02  1.5403E+00  1.5805E+00  8.7932E-01  9.5811E-01  4.6820E+00  9.9470E-01  6.8926E-01  1.1742E+00
             9.7514E-01
 PARAMETER:  1.3984E-01 -2.2172E+00  5.3196E-01  5.5775E-01 -2.8609E-02  5.7210E-02  1.6437E+00  9.4689E-02 -2.7214E-01  2.6060E-01
             7.4824E-02
 GRADIENT:  -5.2621E+00  5.0608E-01 -9.9398E-01  1.4092E+01 -3.7715E-01 -4.6357E-01 -2.4703E-01 -8.6701E-01 -3.0963E-02  1.2003E+00
             3.0699E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1712.17256532811        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1972
 NPARAMETR:  1.0421E+00  8.9668E-02  1.6044E+00  1.5850E+00  8.9447E-01  9.5895E-01  4.8777E+00  1.0737E+00  6.8730E-01  1.1763E+00
             9.7502E-01
 PARAMETER:  1.4126E-01 -2.3116E+00  5.7276E-01  5.6058E-01 -1.1522E-02  5.8089E-02  1.6847E+00  1.7114E-01 -2.7498E-01  2.6238E-01
             7.4706E-02
 GRADIENT:  -1.4013E+00  3.9894E-02  2.8068E-01  2.5648E+00 -1.1744E+00 -5.8055E-02 -1.8095E-01 -1.1312E-01  4.1877E-01  3.1464E-01
             5.0554E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1712.23809013084        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2158
 NPARAMETR:  1.0440E+00  8.8351E-02  1.6130E+00  1.5795E+00  8.9669E-01  9.5940E-01  4.9107E+00  1.0805E+00  6.8907E-01  1.1764E+00
             9.7393E-01
 PARAMETER:  1.4310E-01 -2.3264E+00  5.7810E-01  5.5708E-01 -9.0418E-03  5.8553E-02  1.6914E+00  1.7746E-01 -2.7241E-01  2.6249E-01
             7.3582E-02
 GRADIENT:   3.3649E+00 -7.8675E-02  1.2341E+00 -1.6283E+01 -3.5388E-01  1.1976E-01  9.1316E-01 -1.0353E-01  7.9554E-01  1.3318E-01
             1.6676E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1712.24553641985        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2343
 NPARAMETR:  1.0441E+00  8.8650E-02  1.6087E+00  1.5783E+00  8.9551E-01  9.5948E-01  4.9091E+00  1.0775E+00  6.8954E-01  1.1753E+00
             9.7351E-01
 PARAMETER:  1.4315E-01 -2.3231E+00  5.7545E-01  5.5632E-01 -1.0364E-02  5.8633E-02  1.6911E+00  1.7465E-01 -2.7173E-01  2.6150E-01
             7.3154E-02
 GRADIENT:   3.5403E+00  1.0008E-01  1.1830E+00 -1.9042E+01 -1.5238E-01  1.4367E-01  1.4955E+00 -5.5442E-02  5.3108E-01  1.2304E-01
            -7.2598E-04

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1712.24553641985        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     2367
 NPARAMETR:  1.0441E+00  8.8650E-02  1.6087E+00  1.5783E+00  8.9551E-01  9.5948E-01  4.9091E+00  1.0775E+00  6.8954E-01  1.1753E+00
             9.7351E-01
 PARAMETER:  1.4315E-01 -2.3231E+00  5.7545E-01  5.5632E-01 -1.0364E-02  5.8633E-02  1.6911E+00  1.7465E-01 -2.7173E-01  2.6150E-01
             7.3154E-02
 GRADIENT:  -2.0576E-01 -3.2958E-01  6.1926E-01  2.9473E+00 -1.4184E-01  8.5464E-03 -3.1727E-01 -3.3107E-02  4.0311E-01  1.1209E-01
            -5.6352E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2367
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.6140E-04  2.8284E-02 -2.9283E-02 -2.3806E-02 -2.7846E-02
 SE:             2.9849E-02  1.3384E-02  1.4191E-02  2.6494E-02  2.2390E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9034E-01  3.4573E-02  3.9066E-02  3.6890E-01  2.1363E-01

 ETASHRINKSD(%)  2.7046E-03  5.5162E+01  5.2458E+01  1.1243E+01  2.4990E+01
 ETASHRINKVR(%)  5.4091E-03  7.9896E+01  7.7397E+01  2.1222E+01  4.3735E+01
 EBVSHRINKSD(%)  4.1518E-01  6.6990E+01  5.6298E+01  7.3968E+00  1.9357E+01
 EBVSHRINKVR(%)  8.2863E-01  8.9103E+01  8.0902E+01  1.4246E+01  3.4967E+01
 RELATIVEINF(%)  9.9000E+01  4.1197E+00  4.1522E+00  3.2150E+01  1.4338E+01
 EPSSHRINKSD(%)  4.4409E+01
 EPSSHRINKVR(%)  6.9096E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1712.2455364198474     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -977.09470985610926     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.88
 Elapsed covariance  time in seconds:     7.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1712.246       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  8.86E-02  1.61E+00  1.58E+00  8.96E-01  9.59E-01  4.91E+00  1.08E+00  6.90E-01  1.18E+00  9.74E-01
 


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
 
         3.03E-02  1.94E-02  6.08E-01  5.25E-02  1.93E-01  5.96E-02  1.97E-01  5.34E-01  5.48E-02  2.18E-01  6.24E-02
 


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
+        9.19E-04
 
 TH 2
+       -4.00E-05  3.77E-04
 
 TH 3
+        1.85E-03  3.38E-03  3.70E-01
 
 TH 4
+        3.01E-04  2.42E-04  2.50E-02  2.76E-03
 
 TH 5
+        5.22E-04  9.83E-04  1.14E-01  8.01E-03  3.71E-02
 
 TH 6
+        1.57E-04 -5.19E-05 -4.71E-04 -1.96E-04  2.77E-04  3.55E-03
 
 TH 7
+        3.10E-04 -2.73E-03 -6.49E-02 -3.95E-03 -1.88E-02  7.91E-04  3.87E-02
 
 TH 8
+        1.39E-03  2.63E-03  2.98E-01  2.08E-02  9.09E-02 -1.89E-03 -5.60E-02  2.85E-01
 
 TH 9
+       -3.69E-04  5.95E-04 -1.40E-02 -1.48E-03 -4.07E-03  1.66E-04 -1.85E-04 -1.17E-02  3.01E-03
 
 TH10
+        3.30E-04  1.71E-03  1.13E-01  7.48E-03  3.70E-02  1.16E-03 -2.38E-02  8.39E-02 -2.44E-03  4.76E-02
 
 TH11
+        2.81E-04 -5.25E-05 -1.41E-03 -2.63E-04 -1.18E-03 -2.76E-04  6.81E-04 -5.61E-04 -2.26E-04 -2.65E-03  3.89E-03
 
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
+        3.03E-02
 
 TH 2
+       -6.80E-02  1.94E-02
 
 TH 3
+        1.00E-01  2.86E-01  6.08E-01
 
 TH 4
+        1.89E-01  2.37E-01  7.83E-01  5.25E-02
 
 TH 5
+        8.94E-02  2.63E-01  9.77E-01  7.91E-01  1.93E-01
 
 TH 6
+        8.71E-02 -4.48E-02 -1.30E-02 -6.27E-02  2.41E-02  5.96E-02
 
 TH 7
+        5.20E-02 -7.16E-01 -5.42E-01 -3.82E-01 -4.95E-01  6.75E-02  1.97E-01
 
 TH 8
+        8.59E-02  2.54E-01  9.18E-01  7.41E-01  8.85E-01 -5.94E-02 -5.34E-01  5.34E-01
 
 TH 9
+       -2.22E-01  5.59E-01 -4.19E-01 -5.13E-01 -3.86E-01  5.07E-02 -1.71E-02 -4.00E-01  5.48E-02
 
 TH10
+        4.98E-02  4.05E-01  8.51E-01  6.52E-01  8.81E-01  8.95E-02 -5.55E-01  7.20E-01 -2.04E-01  2.18E-01
 
 TH11
+        1.49E-01 -4.34E-02 -3.71E-02 -8.03E-02 -9.79E-02 -7.42E-02  5.55E-02 -1.69E-02 -6.61E-02 -1.95E-01  6.24E-02
 
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
+        1.24E+03
 
 TH 2
+       -1.09E+03  5.29E+04
 
 TH 3
+        1.45E+01 -1.29E+03  1.32E+02
 
 TH 4
+        1.98E+01 -1.06E+04  2.87E+02  3.31E+03
 
 TH 5
+       -2.01E+01  4.34E+03 -3.41E+02 -1.15E+03  1.27E+03
 
 TH 6
+       -7.58E+01  1.75E+02  1.01E+01  3.76E+01 -3.64E+01  3.07E+02
 
 TH 7
+       -7.38E+01  2.65E+03 -4.55E+01 -5.06E+02  1.45E+02  1.23E+00  1.82E+02
 
 TH 8
+       -4.44E+00  3.50E+02 -3.12E+01 -8.54E+01  2.53E+01  4.69E+00  2.34E+01  2.78E+01
 
 TH 9
+        3.94E+02 -1.50E+04  4.28E+02  3.27E+03 -1.34E+03 -4.59E+01 -7.05E+02 -1.03E+02  4.82E+03
 
 TH10
+       -2.58E+00 -6.51E+02  6.85E+00  1.48E+02 -1.91E+02 -2.46E+01 -8.96E+00  1.24E+01  1.38E+02  1.38E+02
 
 TH11
+       -7.60E+01 -7.94E+02 -2.08E+01  2.06E+02  9.83E+00  5.82E+00 -4.11E+01 -1.66E+00  2.20E+02  4.96E+01  3.15E+02
 
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
 #CPUT: Total CPU Time in Seconds,       41.558
Stop Time:
Wed Sep 29 18:57:38 CDT 2021
