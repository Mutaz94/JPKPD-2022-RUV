Sat Sep 18 15:20:02 CDT 2021
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
$DATA ../../../../data/spa/D/dat41.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13987.1606812327        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1993E+02  1.7259E+02 -1.1619E+01  1.1689E+02  1.8542E+02 -1.4295E+03 -7.1179E+02 -6.7235E+01 -1.0716E+03 -4.3191E+02
            -2.7444E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -621.779512657361        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4324E+00  1.2732E+00  9.5038E-01  1.6670E+00  1.1505E+00  1.3236E+00  1.1971E+00  9.8183E-01  1.1302E+00  1.0975E+00
             1.5171E+01
 PARAMETER:  4.5936E-01  3.4156E-01  4.9109E-02  6.1103E-01  2.4022E-01  3.8037E-01  2.7994E-01  8.1659E-02  2.2241E-01  1.9305E-01
             2.8194E+00
 GRADIENT:  -1.7835E-01  3.6949E+01 -2.9846E+00  6.9322E+01 -1.2228E+01 -2.5439E+00 -5.8175E-01  3.3293E+00  7.0728E+00  3.4353E+00
             9.8327E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -639.580824822285        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.3253E+00  1.0807E+00  3.4747E+00  1.6477E+00  3.0775E+00  1.3045E+00  1.5249E+00  3.0446E-01  1.0674E+00  5.3799E+00
             1.3870E+01
 PARAMETER:  3.8165E-01  1.7763E-01  1.3455E+00  5.9941E-01  1.2241E+00  3.6579E-01  5.2191E-01 -1.0892E+00  1.6525E-01  1.7827E+00
             2.7297E+00
 GRADIENT:  -7.3441E+00  2.4761E+01  2.7244E+00  5.0994E+01 -7.7898E+00  1.1804E+00  1.5326E+00  2.3199E-03  5.3755E+00  4.1129E+00
             8.9523E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -662.236093430475        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1505E+00  4.4204E-01  1.9192E+00  1.6302E+00  3.1982E+00  1.2078E+00  2.4351E+00  2.2572E-01  5.4963E-01  5.0206E+00
             1.1254E+01
 PARAMETER:  2.4017E-01 -7.1636E-01  7.5189E-01  5.8869E-01  1.2626E+00  2.8879E-01  9.9001E-01 -1.3885E+00 -4.9852E-01  1.7136E+00
             2.5207E+00
 GRADIENT:  -6.9489E+00  1.0287E+01  7.4720E+00  3.8349E+01 -4.7665E+00 -2.1976E+00  7.3434E-01  8.6813E-04 -5.9375E+00 -1.4877E+00
            -3.5360E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -677.144031146315        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0459E+00  3.7132E-01  6.4462E-01  1.2700E+00  5.8697E+00  1.0614E+00  8.2481E-01  2.0721E-02  4.0110E-01  7.4212E+00
             1.1195E+01
 PARAMETER:  1.4492E-01 -8.9069E-01 -3.3910E-01  3.3904E-01  1.8698E+00  1.5957E-01 -9.2599E-02 -3.7766E+00 -8.1355E-01  2.1043E+00
             2.5154E+00
 GRADIENT:   1.9549E+01 -8.5721E-01  8.2153E+00 -9.1534E+01 -9.5351E+00 -8.6617E+00  4.1502E-01  7.6295E-05  3.8960E+00 -1.1818E+00
             2.1905E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -684.000334967196        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.0525E+00  2.5175E-01  6.3622E-01  1.3813E+00  6.9259E+00  1.0735E+00  7.3172E-01  1.1699E-02  3.7748E-01  7.6959E+00
             1.1667E+01
 PARAMETER:  1.5113E-01 -1.2793E+00 -3.5221E-01  4.2304E-01  2.0353E+00  1.7090E-01 -2.1236E-01 -4.3483E+00 -8.7423E-01  2.1407E+00
             2.5567E+00
 GRADIENT:   2.1237E+00  1.1888E+00  7.8798E+00 -4.2047E+01 -4.0025E+00 -1.0217E+01  1.3373E-01 -8.6487E-05  3.7727E+00  2.4565E+00
             2.6871E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -688.884387337655        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.0797E+00  8.6654E-02  6.1594E-01  1.5241E+00  7.4178E+00  1.1260E+00  4.6530E-01  1.4314E-02  2.4521E-01  7.2260E+00
             1.2003E+01
 PARAMETER:  1.7667E-01 -2.3458E+00 -3.8461E-01  5.2140E-01  2.1039E+00  2.1866E-01 -6.6508E-01 -4.1465E+00 -1.3056E+00  2.0777E+00
             2.5851E+00
 GRADIENT:   2.1308E+01  1.0426E+00  5.2830E+00 -4.6893E+00 -5.7816E-01 -4.2247E+00  3.9977E-03 -1.7436E-04  1.7399E+00 -2.6325E+00
             2.9963E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -689.530624835751        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  1.0741E+00  6.6681E-02  6.1772E-01  1.5374E+00  7.4713E+00  1.1385E+00  4.3600E-01  1.6268E-02  2.2962E-01  7.1939E+00
             1.1898E+01
 PARAMETER:  1.7151E-01 -2.6078E+00 -3.8172E-01  5.3012E-01  2.1111E+00  2.2969E-01 -7.3011E-01 -4.0186E+00 -1.3713E+00  2.0732E+00
             2.5764E+00
 GRADIENT:   1.3902E+01  4.2112E-01  1.0753E+01 -8.9180E+00  8.4723E+00 -3.0957E+00  2.7451E-03 -8.8436E-05  1.7552E+00 -8.4453E+00
             3.2041E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -689.695583864173        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  1.0722E+00  6.1686E-02  6.1827E-01  1.5410E+00  7.4856E+00  1.1422E+00  4.2805E-01  1.6931E-02  2.2535E-01  7.1863E+00
             1.1863E+01
 PARAMETER:  1.6974E-01 -2.6857E+00 -3.8083E-01  5.3240E-01  2.1130E+00  2.3297E-01 -7.4851E-01 -3.9786E+00 -1.3901E+00  2.0722E+00
             2.5734E+00
 GRADIENT:   3.7108E+01  3.7577E+00 -3.5967E+01  1.0041E+02 -4.5976E+01  1.2411E+01 -5.2500E-03 -2.1617E-03 -2.2937E+00  1.6841E+01
            -2.4211E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -690.318632024989        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  1.0721E+00  6.1613E-02  6.1836E-01  1.5407E+00  7.5270E+00  1.1421E+00  4.1765E-01  1.6978E-02  2.2548E-01  7.1837E+00
             1.1801E+01
 PARAMETER:  1.6966E-01 -2.6869E+00 -3.8069E-01  5.3224E-01  2.1185E+00  2.3287E-01 -7.7312E-01 -3.9758E+00 -1.3895E+00  2.0718E+00
             2.5681E+00
 GRADIENT:   1.8367E+01  1.3034E+00  5.0126E+00  1.1524E+01 -6.3241E+00 -5.7205E+00  2.4664E-03 -1.6676E-04  1.6766E+00 -2.9609E+00
             1.7340E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -690.360275017906        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      925
 NPARAMETR:  1.0716E+00  6.1607E-02  6.1836E-01  1.5407E+00  7.5275E+00  1.1422E+00  4.1764E-01  1.7118E-02  2.2019E-01  7.1833E+00
             1.1801E+01
 PARAMETER:  1.6915E-01 -2.6870E+00 -3.8068E-01  5.3222E-01  2.1186E+00  2.3296E-01 -7.7314E-01 -3.9676E+00 -1.4133E+00  2.0718E+00
             2.5682E+00
 GRADIENT:   3.8930E+01  3.8474E+00 -3.3537E+01  9.6721E+01 -4.9454E+01  1.0894E+01 -4.7359E-03 -2.1591E-03 -2.0710E+00  1.6333E+01
            -2.4505E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -690.664461371334        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1000
 NPARAMETR:  1.0463E+00  6.1294E-02  6.1842E-01  1.5382E+00  7.5566E+00  1.1218E+00  4.2014E-01  2.8139E-01  1.7533E-01  7.2029E+00
             1.1712E+01
 PARAMETER:  1.4524E-01 -2.6921E+00 -3.8058E-01  5.3059E-01  2.1224E+00  2.1498E-01 -7.6718E-01 -1.1680E+00 -1.6411E+00  2.0745E+00
             2.5606E+00
 GRADIENT:  -2.0379E+01  3.3837E-01  1.8716E+01 -1.1612E+01  1.4262E+01 -9.1151E+00  2.5700E-03  4.7192E-02  1.2646E+00 -1.5850E+01
             3.4454E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -691.132756016292        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1073
 NPARAMETR:  1.0444E+00  5.9807E-02  6.1880E-01  1.5265E+00  7.6879E+00  1.1227E+00  4.3060E-01  2.2468E-01  1.7649E-01  7.3023E+00
             1.1303E+01
 PARAMETER:  1.4344E-01 -2.7166E+00 -3.7998E-01  5.2299E-01  2.1396E+00  2.1571E-01 -7.4257E-01 -1.3931E+00 -1.6345E+00  2.0882E+00
             2.5251E+00
 GRADIENT:  -4.6746E+00  4.9350E-01  1.6696E+01 -1.1096E+01  8.9871E+00 -8.0512E+00  1.9932E-03 -1.2113E-02  1.0221E+00 -1.6290E+01
            -7.0769E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -691.188486181564        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1145
 NPARAMETR:  1.0409E+00  5.9707E-02  6.1875E-01  1.5257E+00  7.6957E+00  1.1465E+00  4.3118E-01  2.7300E-01  1.8728E-01  7.3219E+00
             1.1290E+01
 PARAMETER:  1.4007E-01 -2.7183E+00 -3.8005E-01  5.2244E-01  2.1407E+00  2.3673E-01 -7.4123E-01 -1.1983E+00 -1.5752E+00  2.0909E+00
             2.5239E+00
 GRADIENT:  -7.5934E+00  4.7875E-01  1.7010E+01 -1.2440E+01  8.6083E+00 -1.9386E+00  1.7595E-03 -3.7054E-02  1.0332E+00 -1.5327E+01
             1.1696E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -691.407546027409        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1219
 NPARAMETR:  1.0378E+00  5.9035E-02  6.1831E-01  1.5205E+00  7.7420E+00  1.1691E+00  4.3416E-01  3.6663E-01  1.8937E-01  7.4924E+00
             1.1268E+01
 PARAMETER:  1.3711E-01 -2.7296E+00 -3.8076E-01  5.1904E-01  2.1467E+00  2.5623E-01 -7.3435E-01 -9.0341E-01 -1.5641E+00  2.1139E+00
             2.5219E+00
 GRADIENT:  -4.6516E+00  8.7218E-01  1.4909E+01 -1.7894E+01 -5.2611E+00  3.3396E+00  1.5445E-03 -8.9344E-02  1.1077E+00 -7.7022E+00
            -8.1638E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -692.050091133685        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1291
 NPARAMETR:  1.0352E+00  5.8465E-02  6.1655E-01  1.5167E+00  7.7531E+00  1.1292E+00  4.3244E-01  3.3023E-01  1.3614E-01  7.9268E+00
             1.1616E+01
 PARAMETER:  1.3456E-01 -2.7393E+00 -3.8361E-01  5.1652E-01  2.1481E+00  2.2149E-01 -7.3831E-01 -1.0080E+00 -1.8941E+00  2.1702E+00
             2.5524E+00
 GRADIENT:  -2.4347E+01  3.9000E-01  1.8440E+01 -3.5858E+01  4.8881E+00 -1.8451E+00  2.5270E-03  4.4370E-02  9.1437E-01 -1.5736E+01
             3.4589E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -692.266760553113        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1367
 NPARAMETR:  1.0223E+00  5.6982E-02  6.1428E-01  1.5055E+00  7.8684E+00  1.1166E+00  4.3743E-01  7.6900E-01  9.5414E-02  8.5725E+00
             1.1807E+01
 PARAMETER:  1.2210E-01 -2.7650E+00 -3.8731E-01  5.0915E-01  2.1629E+00  2.1031E-01 -7.2683E-01 -1.6266E-01 -2.2495E+00  2.2486E+00
             2.5687E+00
 GRADIENT:  -4.3615E+01 -2.5438E-01  2.2537E+01 -6.0289E+01  9.9655E+00 -7.6536E+00  3.1339E-03  5.7104E-01  5.3971E-01 -9.8332E+00
             6.6855E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -693.301237583922        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:     1446
 NPARAMETR:  1.0294E+00  5.7500E-02  6.1509E-01  1.5094E+00  7.8568E+00  1.1087E+00  4.3766E-01  7.5874E-01  1.0537E-01  8.2856E+00
             1.1653E+01
 PARAMETER:  1.2900E-01 -2.7560E+00 -3.8599E-01  5.1172E-01  2.1614E+00  2.0321E-01 -7.2632E-01 -1.7609E-01 -2.1503E+00  2.2145E+00
             2.5556E+00
 GRADIENT:  -4.4049E+01  2.0635E-01  2.3134E+01 -4.7819E+01  1.1529E+01  8.9691E-01  2.4188E-03  7.1049E-01  5.9934E-01 -3.4581E+01
             5.7977E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -693.370287825046        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:     1555
 NPARAMETR:  1.0302E+00  5.7503E-02  6.1504E-01  1.5095E+00  7.8596E+00  1.1027E+00  4.3777E-01  7.5467E-01  1.0310E-01  8.2888E+00
             1.1656E+01
 PARAMETER:  1.2975E-01 -2.7559E+00 -3.8607E-01  5.1175E-01  2.1617E+00  1.9775E-01 -7.2607E-01 -1.8148E-01 -2.1720E+00  2.2149E+00
             2.5558E+00
 GRADIENT:  -6.3898E+01  2.0763E-01  2.7951E+01 -6.0162E+01  2.0617E+01  1.3088E+01  2.4132E-03  1.1507E+00  6.5045E-01 -6.7906E+01
             7.5295E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -693.426799953395        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1695
 NPARAMETR:  1.0302E+00  5.6558E-02  6.1505E-01  1.5094E+00  7.8602E+00  1.1026E+00  4.3777E-01  7.6923E-01  8.9067E-02  8.2885E+00
             1.1656E+01
 PARAMETER:  1.2977E-01 -2.7725E+00 -3.8606E-01  5.1174E-01  2.1618E+00  1.9770E-01 -7.2607E-01 -1.6236E-01 -2.3184E+00  2.2149E+00
             2.5558E+00
 GRADIENT:  -3.4684E+02 -1.2168E-02  8.8408E+01 -2.4620E+02  1.5313E+02  2.1304E+02  2.6035E-03  6.6298E+00  1.4303E+00 -5.3720E+02
             3.6682E+02

0ITERATION NO.:   96    OBJECTIVE VALUE:  -693.426799953395        NO. OF FUNC. EVALS.:  36
 CUMULATIVE NO. OF FUNC. EVALS.:     1731
 NPARAMETR:  1.0302E+00  5.6715E-02  6.1509E-01  1.5093E+00  7.8633E+00  1.1027E+00  4.3764E-01  7.6926E-01  8.9109E-02  8.2854E+00
             1.1661E+01
 PARAMETER:  1.2977E-01 -2.7725E+00 -3.8606E-01  5.1174E-01  2.1618E+00  1.9770E-01 -7.2607E-01 -1.6236E-01 -2.3184E+00  2.2149E+00
             2.5558E+00
 GRADIENT:   1.2299E+03 -1.1225E-02 -4.0562E+02  1.1535E+03 -2.7468E+02 -8.3880E+02  1.2351E-02 -1.0164E+03 -7.0686E+01  2.6392E+02
            -1.9170E+02
 NUMSIGDIG:         3.3         2.6         3.3         3.3         3.3         3.3         2.9         3.3         3.3         3.4
                    3.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1731
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3923E-02 -2.8022E-04 -1.3182E-02 -2.1697E-03 -4.0104E-02
 SE:             2.8521E-02  1.7510E-04  4.4027E-03  2.3471E-03  8.9930E-03
 N:                     100         100         100         100         100

 P VAL.:         4.0159E-01  1.0951E-01  2.7518E-03  3.5528E-01  8.2218E-06

 ETASHRINKSD(%)  4.4526E+00  9.9413E+01  8.5251E+01  9.2137E+01  6.9872E+01
 ETASHRINKVR(%)  8.7069E+00  9.9997E+01  9.7825E+01  9.9382E+01  9.0923E+01
 EBVSHRINKSD(%)  9.1189E+00  9.9319E+01  8.0861E+01  9.1362E+01  6.7218E+01
 EBVSHRINKVR(%)  1.7406E+01  9.9995E+01  9.6337E+01  9.9254E+01  8.9253E+01
 RELATIVEINF(%)  1.6408E+01  7.7464E-05  4.2694E-01  9.9760E-03  5.1926E+00
 EPSSHRINKSD(%)  7.4181E+00
 EPSSHRINKVR(%)  1.4286E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -693.42679995339472     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       41.724026610343458     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -693.427       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  5.66E-02  6.15E-01  1.51E+00  7.86E+00  1.10E+00  4.38E-01  7.69E-01  8.91E-02  8.29E+00  1.17E+01
 


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
+        2.42E+06
 
 TH 2
+       -1.36E+04  1.67E+06
 
 TH 3
+       -6.46E+04  6.49E+03  7.74E+05
 
 TH 4
+        1.47E+04  9.40E+05 -6.34E+05  2.64E+05
 
 TH 5
+       -5.56E+02  5.57E+01  4.88E+02 -1.09E+02  5.73E+02
 
 TH 6
+       -5.94E+04  8.60E+03  3.88E+04 -9.83E+03  3.56E+02  9.11E+05
 
 TH 7
+        2.23E+01  8.26E+05 -4.52E+00 -1.23E+00 -4.92E-02 -6.46E+00  2.99E+00
 
 TH 8
+       -4.51E+03  1.53E+04  2.66E+03 -7.75E+02  3.09E+01  2.72E+03  3.02E+00  2.66E+06
 
 TH 9
+       -9.11E+02  9.25E+03  6.85E+02 -1.67E+02  6.13E+00  5.13E+02  6.53E-02  4.30E+03  9.72E+05
 
 TH10
+        3.78E+02 -4.56E+01 -2.65E+04  8.12E+03 -2.60E+01 -2.44E+02  3.21E-02 -2.13E+01 -4.05E+00  5.47E+02
 
 TH11
+       -2.93E+02  2.35E+01  2.48E+02 -6.36E+01  5.21E+02  1.81E+02 -2.75E-02  1.59E+01  3.66E+00 -2.67E+01  3.58E+02
 
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
 #CPUT: Total CPU Time in Seconds,       38.664
Stop Time:
Sat Sep 18 15:20:43 CDT 2021
