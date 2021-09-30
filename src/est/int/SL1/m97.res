Wed Sep 29 02:34:44 CDT 2021
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
$DATA ../../../../data/int/SL1/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3029.10338991665        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6988E+02 -4.3135E+00  1.5904E+02  1.2492E+02  1.7218E+02  3.1641E+01 -3.4227E+01 -1.3134E+02 -5.7224E+01  9.1125E+00
            -1.3665E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3307.58142042044        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  8.7680E-01  1.0535E+00  6.9144E-01  9.6426E-01  8.2973E-01  9.7900E-01  1.0060E+00  1.2150E+00  1.0869E+00  8.0246E-01
             1.7012E+00
 PARAMETER: -3.1472E-02  1.5212E-01 -2.6897E-01  6.3610E-02 -8.6651E-02  7.8773E-02  1.0594E-01  2.9478E-01  1.8332E-01 -1.2008E-01
             6.3135E-01
 GRADIENT:  -3.7866E+00  9.1613E+01 -4.4163E+01 -3.0319E+01 -4.1619E+01 -2.2404E+00 -2.1110E+00  7.6186E-01  1.2927E+01  6.3149E+00
             2.0371E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3310.96633343914        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      168
 NPARAMETR:  8.7832E-01  9.3654E-01  6.3599E-01  1.0512E+00  7.4502E-01  9.7982E-01  1.0723E+00  1.3522E+00  1.0815E+00  7.0046E-01
             1.6750E+00
 PARAMETER: -2.9747E-02  3.4434E-02 -3.5258E-01  1.4994E-01 -1.9435E-01  7.9611E-02  1.6982E-01  4.0171E-01  1.7832E-01 -2.5601E-01
             6.1583E-01
 GRADIENT:   3.0686E+00  8.1064E+01 -4.0967E+01  4.2454E+01  4.8597E+00 -1.1780E+00 -1.0080E+01  1.3008E+01  2.1252E+01  5.4684E+00
             2.1366E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3311.72732065452        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      346
 NPARAMETR:  8.7856E-01  9.3537E-01  6.3546E-01  1.0515E+00  7.4462E-01  1.0522E+00  1.0730E+00  1.3513E+00  1.0806E+00  7.0021E-01
             1.6773E+00
 PARAMETER: -2.9473E-02  3.3183E-02 -3.5341E-01  1.5019E-01 -1.9488E-01  1.5084E-01  1.7049E-01  4.0107E-01  1.7756E-01 -2.5638E-01
             6.1717E-01
 GRADIENT:  -1.1844E+02  5.3305E+01 -5.3735E+01 -9.6284E+00 -3.8926E+01  1.1162E+00 -1.4618E+01  1.1011E+01  1.5178E+01  3.0233E+00
             2.0629E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3315.41105050349        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      492             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5598E-01  9.3540E-01  6.3553E-01  1.0514E+00  7.4457E-01  1.0536E+00  1.0991E+00  1.3511E+00  1.0597E+00  6.3323E-01
             1.6776E+00
 PARAMETER:  5.4986E-02  3.3216E-02 -3.5329E-01  1.5014E-01 -1.9495E-01  1.5221E-01  1.9450E-01  4.0094E-01  1.5797E-01 -3.5693E-01
             6.1737E-01
 GRADIENT:   1.8914E+02  8.2966E+01 -4.1126E+01  4.1033E+01  9.7156E+00  4.1533E+01 -1.2224E+01  1.4809E+01  1.5053E+01  2.3010E+00
             2.1064E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3316.10905897434        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  9.3314E-01  9.3540E-01  6.3553E-01  1.0514E+00  7.4457E-01  1.0295E+00  1.0991E+00  1.3511E+00  1.0597E+00  6.4359E-01
             1.6775E+00
 PARAMETER:  3.0797E-02  3.3216E-02 -3.5329E-01  1.5014E-01 -1.9495E-01  1.2912E-01  1.9450E-01  4.0094E-01  1.5797E-01 -3.4069E-01
             6.1733E-01
 GRADIENT:   6.5417E-01  5.6256E+01 -5.3816E+01 -1.0370E+01 -3.5123E+01  1.3580E-01 -1.6717E+01  1.2495E+01  9.7687E+00  4.0381E-02
             2.0194E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3317.59169174653        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      743
 NPARAMETR:  9.2505E-01  9.3492E-01  6.3667E-01  1.0539E+00  7.4384E-01  1.0205E+00  1.2071E+00  1.1812E+00  1.0201E+00  6.2555E-01
             1.6827E+00
 PARAMETER:  2.2096E-02  3.2711E-02 -3.5150E-01  1.5249E-01 -1.9593E-01  1.2031E-01  2.8819E-01  2.6652E-01  1.1986E-01 -3.6912E-01
             6.2040E-01
 GRADIENT:  -1.7777E+01  5.5371E+01 -6.2904E+01 -7.8702E+00 -1.4335E+01 -3.6486E+00 -1.0386E+00  1.5920E+00  5.9197E-01  6.1565E+00
             1.9397E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3318.61065038149        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      923
 NPARAMETR:  9.3297E-01  9.3491E-01  6.3683E-01  1.0584E+00  7.4383E-01  1.0182E+00  1.2468E+00  1.1855E+00  1.0114E+00  5.1480E-01
             1.6793E+00
 PARAMETER:  3.0618E-02  3.2692E-02 -3.5126E-01  1.5672E-01 -1.9594E-01  1.1799E-01  3.2062E-01  2.7020E-01  1.1131E-01 -5.6398E-01
             6.1838E-01
 GRADIENT:   5.0370E-02  5.8762E+01 -6.3148E+01  3.6418E+00  1.1196E+00 -4.5076E+00 -3.5326E+00  3.4102E+00 -3.1433E+00 -6.9318E-01
             1.8598E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3324.25331397712        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  9.2023E-01  8.7943E-01  6.3922E-01  1.0576E+00  7.3660E-01  1.0138E+00  1.2550E+00  1.3190E+00  1.0228E+00  5.5563E-01
             1.6405E+00
 PARAMETER:  1.6866E-02 -2.8481E-02 -3.4750E-01  1.5601E-01 -2.0571E-01  1.1374E-01  3.2716E-01  3.7684E-01  1.2251E-01 -4.8765E-01
             5.9499E-01
 GRADIENT:   1.1817E+02  3.4040E+01 -3.8669E+01  3.6866E+01  7.8373E+01  2.2242E+01  1.0260E+01  1.1812E+01  3.4197E+00  5.2233E+00
             1.8114E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3324.66169683572        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1136
 NPARAMETR:  9.0030E-01  8.2243E-01  6.3933E-01  1.0667E+00  6.9325E-01  1.0003E+00  1.2835E+00  1.3187E+00  1.0134E+00  4.5782E-01
             1.6409E+00
 PARAMETER: -5.0287E-03 -9.5486E-02 -3.4734E-01  1.6458E-01 -2.6636E-01  1.0034E-01  3.4960E-01  3.7666E-01  1.1329E-01 -6.8127E-01
             5.9522E-01
 GRADIENT:   6.8860E+01  2.7986E+01 -2.7121E+00  1.4858E+01  6.0111E+01  1.2592E+01  4.5829E+00  1.6764E+01  3.0364E+00  3.9126E+00
             1.8862E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3327.49677213944        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1317
 NPARAMETR:  9.3249E-01  7.8403E-01  6.3931E-01  1.1118E+00  6.7433E-01  1.0288E+00  1.3817E+00  1.3188E+00  1.0028E+00  3.7505E-01
             1.6404E+00
 PARAMETER:  3.0105E-02 -1.4330E-01 -3.4737E-01  2.0596E-01 -2.9404E-01  1.2835E-01  4.2331E-01  3.7673E-01  1.0276E-01 -8.8070E-01
             5.9494E-01
 GRADIENT:  -3.5461E-01 -4.8099E+00 -7.2611E+00  1.0056E+01  1.6849E+01 -7.6301E-02  2.1522E+00  1.4801E+01  2.4078E-02 -6.4784E-02
             1.8080E+02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3336.48054227019        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1498
 NPARAMETR:  9.3420E-01  7.8085E-01  6.3797E-01  1.1100E+00  6.6782E-01  1.0313E+00  1.3794E+00  1.1689E+00  1.0028E+00  3.7812E-01
             1.4734E+00
 PARAMETER:  3.1940E-02 -1.4738E-01 -3.4946E-01  2.0435E-01 -3.0373E-01  1.3080E-01  4.2164E-01  2.5606E-01  1.0276E-01 -8.7255E-01
             4.8756E-01
 GRADIENT:   6.8769E+00  5.0337E+00  9.3484E+00  5.6898E+00  1.6046E+01  3.2410E-01 -2.4984E+00 -1.2763E+01 -4.4368E+00 -2.1885E+00
            -2.7411E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3336.95715426085        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     1619
 NPARAMETR:  9.2500E-01  7.8030E-01  6.2651E-01  1.1096E+00  6.5874E-01  1.0220E+00  1.3790E+00  1.1821E+00  1.0152E+00  4.2256E-01
             1.4758E+00
 PARAMETER:  2.2039E-02 -1.4808E-01 -3.6759E-01  2.0400E-01 -3.1743E-01  1.2180E-01  4.2134E-01  2.6729E-01  1.1511E-01 -7.6144E-01
             4.8918E-01
 GRADIENT:   1.6533E+02  5.7347E+01  2.7619E+01  9.9351E+01  9.6493E+01  3.3641E+01  2.0004E+01 -8.7381E+00  5.6050E+00  4.4111E+00
            -6.5491E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3337.29153566497        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1812             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2894E-01  7.8037E-01  6.2535E-01  1.1094E+00  6.5890E-01  1.0275E+00  1.3793E+00  1.2142E+00  1.0149E+00  4.1581E-01
             1.4795E+00
 PARAMETER:  2.6291E-02 -1.4798E-01 -3.6944E-01  2.0385E-01 -3.1718E-01  1.2716E-01  4.2160E-01  2.9412E-01  1.1480E-01 -7.7752E-01
             4.9168E-01
 GRADIENT:   1.7337E+02  5.7055E+01  2.5509E+01  9.8768E+01  9.7444E+01  3.6943E+01  1.9568E+01 -5.1126E+00  5.0971E+00  4.1623E+00
             2.1487E+00

0ITERATION NO.:   69    OBJECTIVE VALUE:  -3337.30938187249        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1951
 NPARAMETR:  9.2871E-01  7.8037E-01  6.2535E-01  1.1094E+00  6.5890E-01  1.0305E+00  1.3793E+00  1.2151E+00  1.0152E+00  4.1115E-01
             1.4828E+00
 PARAMETER:  2.6041E-02 -1.4798E-01 -3.6944E-01  2.0385E-01 -3.1718E-01  1.3005E-01  4.2160E-01  2.9485E-01  1.1509E-01 -7.8881E-01
             4.9390E-01
 GRADIENT:   1.2565E+05  5.5505E+00  4.8790E+01 -1.9972E+01 -2.7114E+01 -1.2344E-01 -5.5503E+00  4.2564E+04  1.0918E+05 -1.5930E+04
            -2.5193E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1951
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8286E-03 -1.2798E-02 -1.7610E-02  3.5414E-03 -2.2278E-02
 SE:             2.9798E-02  2.6710E-02  2.3561E-02  2.8293E-02  1.2687E-02
 N:                     100         100         100         100         100

 P VAL.:         8.9777E-01  6.3184E-01  4.5480E-01  9.0039E-01  7.9092E-02

 ETASHRINKSD(%)  1.7192E-01  1.0517E+01  2.1068E+01  5.2146E+00  5.7497E+01
 ETASHRINKVR(%)  3.4355E-01  1.9927E+01  3.7698E+01  1.0157E+01  8.1935E+01
 EBVSHRINKSD(%)  5.8861E-01  9.5432E+00  2.3245E+01  5.0853E+00  5.9117E+01
 EBVSHRINKVR(%)  1.1738E+00  1.8176E+01  4.1087E+01  9.9120E+00  8.3286E+01
 RELATIVEINF(%)  9.8821E+01  2.5888E+01  2.1746E+01  6.2192E+01  2.9639E+00
 EPSSHRINKSD(%)  2.1017E+01
 EPSSHRINKVR(%)  3.7618E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3337.3093818724860     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1683.2200221040753     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    51.86
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3337.309       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.29E-01  7.80E-01  6.25E-01  1.11E+00  6.59E-01  1.03E+00  1.38E+00  1.22E+00  1.02E+00  4.11E-01  1.48E+00
 


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
+        3.64E+07
 
 TH 2
+       -7.16E-01  4.71E+07
 
 TH 3
+        4.17E-01  1.04E+02  1.18E+07
 
 TH 4
+       -2.70E+00 -1.20E+07 -8.53E+01  1.23E+07
 
 TH 5
+       -2.93E+00 -9.24E+02 -1.12E+03  3.52E+02  1.44E+07
 
 TH 6
+        1.47E+03 -1.02E+00  3.03E+00  1.05E-02 -4.53E+00  1.83E+02
 
 TH 7
+        4.41E-02  4.68E+06 -4.61E+00 -2.39E+06  1.53E+01 -1.46E-01  1.86E+06
 
 TH 8
+       -9.45E+06  5.30E+00 -4.48E+01  1.11E+01 -3.68E+01  3.78E+02 -3.09E+00  4.89E+06
 
 TH 9
+       -2.90E+07 -8.33E+00  1.17E+01  1.89E+01 -5.91E+01  1.16E+03 -3.68E+00  7.50E+06  4.60E+07
 
 TH10
+       -2.89E+02 -1.16E+01 -4.06E+00 -1.95E+01 -5.11E+01 -4.18E+02  2.15E+01  3.31E+01  1.06E+02  2.99E+06
 
 TH11
+       -9.54E+00  3.69E+06 -5.42E+01 -8.53E+00 -1.42E+01  1.60E+00  1.34E+01 -1.85E+04 -3.70E+06  1.97E+01  5.95E+05
 
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
 #CPUT: Total CPU Time in Seconds,       64.630
Stop Time:
Wed Sep 29 02:35:51 CDT 2021
