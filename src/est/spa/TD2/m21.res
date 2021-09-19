Sat Sep 18 14:32:11 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1742.03409203237        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.2564E+01 -4.0666E+01 -1.0531E+01 -7.0354E+01 -4.5818E+01  1.8811E+00  2.6232E-01  1.6568E+01  3.3538E+00  1.8301E+01
             3.6075E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1751.92618107663        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0519E+00  1.0531E+00  1.1912E+00  1.0033E+00  1.1531E+00  9.7539E-01  1.0092E+00  8.7452E-01  9.7167E-01  9.5357E-01
             8.9647E-01
 PARAMETER:  1.5059E-01  1.5172E-01  2.7499E-01  1.0326E-01  2.4247E-01  7.5087E-02  1.0914E-01 -3.4080E-02  7.1261E-02  5.2459E-02
            -9.2927E-03
 GRADIENT:   9.3594E+01 -1.0390E+01  8.7316E+00 -2.7074E+01  1.3147E+01 -2.6023E+00 -2.0758E+00 -3.9528E-01 -6.5787E+00 -1.8084E+01
            -1.8847E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1752.95340952766        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0418E+00  9.7658E-01  1.0499E+00  1.0510E+00  1.0421E+00  1.0065E+00  9.9900E-01  5.4653E-01  9.4906E-01  9.4329E-01
             8.9840E-01
 PARAMETER:  1.4092E-01  7.6302E-02  1.4865E-01  1.4974E-01  1.4126E-01  1.0645E-01  9.9003E-02 -5.0416E-01  4.7713E-02  4.1614E-02
            -7.1422E-03
 GRADIENT:   6.4407E+01 -8.5715E+00  7.1654E+00 -1.0591E+01  8.1127E-02  9.9444E+00 -7.4876E+00 -1.2132E+00 -5.6243E+00 -9.8215E+00
            -1.5718E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1753.42739257526        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0296E+00  9.9309E-01  9.7127E-01  1.0424E+00  1.0178E+00  9.9648E-01  1.1432E+00  4.8501E-01  9.0870E-01  9.4373E-01
             9.0640E-01
 PARAMETER:  1.2916E-01  9.3070E-02  7.0854E-02  1.4150E-01  1.1762E-01  9.6470E-02  2.3387E-01 -6.2358E-01  4.2609E-03  4.2082E-02
             1.7214E-03
 GRADIENT:   2.7668E+01 -6.3863E+00 -6.7495E+00  9.7076E-01  8.0957E+00  4.7348E+00 -1.7584E+00  9.1177E-01 -2.9408E+00 -1.6745E+00
            -7.5371E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1754.36777215441        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      409
 NPARAMETR:  1.0509E+00  8.4741E-01  1.0773E+00  1.1475E+00  1.0035E+00  9.9698E-01  1.3164E+00  5.2319E-01  8.6771E-01  9.7139E-01
             9.2221E-01
 PARAMETER:  1.4965E-01 -6.5574E-02  1.7447E-01  2.3760E-01  1.0348E-01  9.6975E-02  3.7493E-01 -5.4780E-01 -4.1894E-02  7.0977E-02
             1.9015E-02
 GRADIENT:   5.3449E+00  1.8062E+00  5.0656E-02  1.3856E+00  9.4954E-01  3.7537E-01  7.2645E-01 -2.2876E-01 -9.2092E-02 -4.6356E-01
            -1.0262E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1754.54320380925        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0460E+00  6.5552E-01  1.2062E+00  1.2717E+00  9.8650E-01  9.9285E-01  1.5238E+00  6.7681E-01  8.2296E-01  9.9174E-01
             9.2315E-01
 PARAMETER:  1.4497E-01 -3.2233E-01  2.8749E-01  3.4032E-01  8.6405E-02  9.2822E-02  5.2124E-01 -2.9036E-01 -9.4842E-02  9.1703E-02
             2.0040E-02
 GRADIENT:   1.2259E+00  2.0626E+00  3.2053E-01  2.9907E+00 -2.9551E+00 -2.2415E-02  3.0926E-01  4.2320E-01 -5.0959E-02  6.2999E-01
             2.3057E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1754.60597987052        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      760
 NPARAMETR:  1.0435E+00  5.4813E-01  1.3338E+00  1.3420E+00  1.0048E+00  9.9121E-01  1.5826E+00  7.9984E-01  8.1277E-01  1.0215E+00
             9.2648E-01
 PARAMETER:  1.4257E-01 -5.0124E-01  3.8803E-01  3.9412E-01  1.0482E-01  9.1167E-02  5.5909E-01 -1.2334E-01 -1.0730E-01  1.2131E-01
             2.3635E-02
 GRADIENT:   2.9278E-01  9.6432E-01  1.7934E+00 -1.1437E+00 -3.0685E+00  9.2724E-02  1.6572E-01  1.5272E-01 -9.7184E-01  5.9987E-01
             1.1298E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1754.78525416850        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  1.0405E+00  3.8521E-01  1.4639E+00  1.4484E+00  1.0050E+00  9.8882E-01  1.5358E+00  9.0930E-01  8.0762E-01  1.0444E+00
             9.2475E-01
 PARAMETER:  1.3967E-01 -8.5396E-01  4.8112E-01  4.7048E-01  1.0497E-01  8.8754E-02  5.2907E-01  4.9172E-03 -1.1366E-01  1.4341E-01
             2.1771E-02
 GRADIENT:   6.8175E-02  1.2043E+00  2.2259E+00  3.6271E+00 -1.8147E+00  4.1028E-02  3.3238E-01 -7.6712E-01  1.2705E+00  9.9568E-03
             3.6316E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1754.97811059517        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  1.0369E+00  2.1329E-01  1.5387E+00  1.5593E+00  9.8023E-01  9.8578E-01  1.1790E+00  9.9171E-01  7.8016E-01  1.0503E+00
             9.2178E-01
 PARAMETER:  1.3628E-01 -1.4451E+00  5.3092E-01  5.4425E-01  8.0035E-02  8.5674E-02  2.6470E-01  9.1679E-02 -1.4826E-01  1.4904E-01
             1.8555E-02
 GRADIENT:  -9.9604E-01  1.6719E+00  1.1580E+00  1.6816E+01 -4.2045E+00 -1.8425E-01  1.6322E-01 -3.1832E-01  1.2061E+00  6.3580E-01
            -5.2951E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1755.09407282304        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  1.0361E+00  1.3901E-01  1.6223E+00  1.6017E+00  9.9040E-01  9.8490E-01  8.2318E-01  1.0729E+00  7.5854E-01  1.0518E+00
             9.2265E-01
 PARAMETER:  1.3542E-01 -1.8732E+00  5.8385E-01  5.7107E-01  9.0355E-02  8.4783E-02 -9.4585E-02  1.7035E-01 -1.7636E-01  1.5053E-01
             1.9493E-02
 GRADIENT:   7.0682E-02  1.3411E-01  9.0803E-01  2.1415E+00 -3.2135E-01 -1.7008E-02  5.8685E-02 -3.3650E-01  1.4072E-01 -4.5437E-01
            -3.1795E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1755.10235452873        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1470             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0357E+00  1.1837E-01  1.6140E+00  1.6127E+00  9.8192E-01  9.8460E-01  6.3116E-01  1.0639E+00  7.5278E-01  1.0492E+00
             9.2232E-01
 PARAMETER:  1.3505E-01 -2.0339E+00  5.7872E-01  5.7790E-01  8.1757E-02  8.4479E-02 -3.6020E-01  1.6192E-01 -1.8398E-01  1.4799E-01
             1.9139E-02
 GRADIENT:   7.2029E+01  1.4433E+00  1.6801E+00  1.3328E+02  5.0756E-01  5.5164E+00  3.5645E-02 -3.2319E-01  2.0141E+00 -2.4872E-01
            -1.9580E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1755.11126012289        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1630
 NPARAMETR:  1.0358E+00  1.2357E-01  1.6126E+00  1.6103E+00  9.8319E-01  9.8480E-01  4.1762E-01  1.0727E+00  7.5571E-01  1.0534E+00
             9.2412E-01
 PARAMETER:  1.3514E-01 -1.9909E+00  5.7783E-01  5.7644E-01  8.3046E-02  8.4681E-02 -7.7318E-01  1.7016E-01 -1.8009E-01  1.5206E-01
             2.1089E-02
 GRADIENT:   8.2497E-02  4.8084E-02 -4.5262E-01  2.1368E+00 -5.3827E-01  5.1457E-02  1.4010E-02  3.3988E-01  3.6090E-01  4.6876E-01
             8.5991E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1755.14937797829        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1805
 NPARAMETR:  1.0375E+00  1.9035E-01  1.6046E+00  1.5677E+00  9.9996E-01  9.8592E-01  1.0000E-02  1.0653E+00  7.7673E-01  1.0621E+00
             9.2164E-01
 PARAMETER:  1.3680E-01 -1.5589E+00  5.7289E-01  5.4958E-01  9.9962E-02  8.5821E-02 -5.3062E+00  1.6326E-01 -1.5266E-01  1.6028E-01
             1.8403E-02
 GRADIENT:   1.1431E+00 -2.1778E-01 -3.0815E-02 -1.9120E+00 -4.1832E-01  4.8856E-02  0.0000E+00  2.5361E-01  4.8357E-01  5.6196E-01
            -6.0205E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1755.16784794166        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1984
 NPARAMETR:  1.0381E+00  2.5405E-01  1.5912E+00  1.5295E+00  1.0145E+00  9.8670E-01  1.0000E-02  1.0490E+00  7.9606E-01  1.0645E+00
             9.2450E-01
 PARAMETER:  1.3735E-01 -1.2702E+00  5.6451E-01  5.2496E-01  1.1439E-01  8.6607E-02 -8.8980E+00  1.4782E-01 -1.2807E-01  1.6251E-01
             2.1498E-02
 GRADIENT:  -1.9800E-01  1.1254E-01  1.0092E-01  1.0968E+00  9.9636E-02 -3.4157E-02  0.0000E+00 -5.3764E-02 -1.6864E-01 -1.0867E-01
             1.2054E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1755.17146648166        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2168
 NPARAMETR:  1.0387E+00  2.8401E-01  1.5814E+00  1.5102E+00  1.0201E+00  9.8723E-01  1.0000E-02  1.0412E+00  8.0673E-01  1.0676E+00
             9.2453E-01
 PARAMETER:  1.3793E-01 -1.1587E+00  5.5829E-01  5.1224E-01  1.1992E-01  8.7143E-02 -1.0569E+01  1.4040E-01 -1.1477E-01  1.6540E-01
             2.1533E-02
 GRADIENT:  -2.4303E-02  7.1276E-03  6.5156E-02 -4.3342E-02 -8.0156E-02 -9.0804E-03  0.0000E+00  1.9368E-02 -2.8218E-03  2.4639E-02
             5.8253E-02

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1755.17147990601        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     2260
 NPARAMETR:  1.0387E+00  2.8399E-01  1.5812E+00  1.5102E+00  1.0202E+00  9.8726E-01  1.0000E-02  1.0408E+00  8.0674E-01  1.0675E+00
             9.2436E-01
 PARAMETER:  1.3793E-01 -1.1588E+00  5.5818E-01  5.1224E-01  1.1995E-01  8.7181E-02 -1.0569E+01  1.3998E-01 -1.1475E-01  1.6535E-01
             2.1351E-02
 GRADIENT:  -2.1073E-02 -4.2936E-03  3.2513E-02 -5.1626E-02  4.6920E-02  4.8565E-03  0.0000E+00 -8.2599E-04  4.8355E-03 -4.0861E-03
            -2.3771E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2260
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0978E-04 -1.1079E-04 -2.8740E-02 -5.5679E-03 -3.7669E-02
 SE:             2.9857E-02  5.0946E-05  1.5179E-02  2.9370E-02  2.2108E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9439E-01  2.9654E-02  5.8301E-02  8.4964E-01  8.8404E-02

 ETASHRINKSD(%)  1.0000E-10  9.9829E+01  4.9150E+01  1.6067E+00  2.5937E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  7.4142E+01  3.1875E+00  4.5146E+01
 EBVSHRINKSD(%)  3.7205E-01  9.9840E+01  5.2603E+01  1.9799E+00  2.2800E+01
 EBVSHRINKVR(%)  7.4271E-01  1.0000E+02  7.7535E+01  3.9205E+00  4.0401E+01
 RELATIVEINF(%)  9.7571E+01  1.4721E-05  5.7542E+00  6.6616E+00  9.6026E+00
 EPSSHRINKSD(%)  4.4304E+01
 EPSSHRINKVR(%)  6.8980E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1755.1714799060064     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1020.0206533422682     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.03
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1755.171       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.84E-01  1.58E+00  1.51E+00  1.02E+00  9.87E-01  1.00E-02  1.04E+00  8.07E-01  1.07E+00  9.24E-01
 


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
+        1.05E+03
 
 TH 2
+       -2.60E+01  3.70E+02
 
 TH 3
+       -8.79E-02  4.47E+01  9.84E+01
 
 TH 4
+       -8.62E+00  4.53E+02 -2.12E+01  7.19E+02
 
 TH 5
+        3.39E+00 -1.92E+02 -1.85E+02 -3.80E+01  5.70E+02
 
 TH 6
+        2.64E+00 -4.23E+00  3.28E-01 -2.24E+00  6.35E-01  1.97E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -7.81E-01 -2.50E+00 -2.54E+01 -3.41E+00 -3.15E+00  1.13E+00  0.00E+00  2.50E+01
 
 TH 9
+        4.29E+00 -9.72E+01  6.55E+00 -4.07E-02 -4.49E+00 -3.19E+00  0.00E+00  3.77E+00  2.89E+02
 
 TH10
+       -7.49E-01  7.85E+00 -1.95E+00 -4.38E-01 -6.13E+01 -4.98E+00  0.00E+00  1.22E+01  1.44E-01  7.20E+01
 
 TH11
+       -5.91E+00 -1.31E+01 -1.69E+01 -9.23E+00 -5.05E+00 -2.80E+00  0.00E+00  1.72E+01  2.55E+00  1.56E+01  2.53E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.752
Stop Time:
Sat Sep 18 14:32:46 CDT 2021
