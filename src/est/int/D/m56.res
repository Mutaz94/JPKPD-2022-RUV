Wed Sep 29 09:12:11 CDT 2021
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
$DATA ../../../../data/int/D/dat56.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   37352.2078150354        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6174E+02  5.4472E+02 -3.2143E+01  3.4225E+02  1.9772E+02 -2.7071E+03 -1.2234E+03 -7.8673E+01 -1.9918E+03 -7.6075E+02
            -7.4903E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -838.553872823693        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2292E+00  2.0289E+00  8.3435E-01  2.3123E+00  9.4914E-01  4.8968E+00  5.7846E+00  9.9997E-01  2.9926E+00  1.6329E+00
             1.2273E+01
 PARAMETER:  3.0636E-01  8.0748E-01 -8.1106E-02  9.3824E-01  4.7797E-02  1.6886E+00  1.8552E+00  9.9968E-02  1.1961E+00  5.9034E-01
             2.6074E+00
 GRADIENT:  -5.5449E+00  2.7225E+01 -5.2775E+01  1.1734E+02 -1.5499E+01  2.1010E+02  1.1272E+02  4.3588E+00  4.8291E+01  3.0409E+01
             3.0849E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -906.724811615270        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.1970E+00  4.0138E+00  8.2188E+00  2.7924E+00  4.5904E+00  4.3319E+00  4.2162E+00  9.0664E-01  9.8476E+00  1.4524E+00
             1.1985E+01
 PARAMETER:  2.7982E-01  1.4897E+00  2.2064E+00  1.1269E+00  1.6240E+00  1.5660E+00  1.5389E+00  1.9882E-03  2.3872E+00  4.7319E-01
             2.5837E+00
 GRADIENT:   1.9205E-01  5.7642E+01 -9.5952E+00  5.1459E+01  5.3216E+01  1.7534E+02  5.4979E+01  1.8929E+00  1.9003E+01  1.7257E+01
             3.1942E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1019.47077420987        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.1893E+00  7.3735E-01  5.2633E+00  1.6507E+00  1.3508E+00  2.5938E+00  4.9272E+00  3.1487E-01  2.1444E+00  1.5404E+00
             1.1532E+01
 PARAMETER:  2.7334E-01 -2.0469E-01  1.7608E+00  6.0119E-01  4.0069E-01  1.0531E+00  1.6948E+00 -1.0556E+00  8.6285E-01  5.3202E-01
             2.5451E+00
 GRADIENT:  -1.1507E+01 -1.9238E+01  1.1287E+00  6.3619E+00 -7.8777E+01  5.6510E+01 -2.0200E+00  2.0445E-01  3.7340E+01  3.3833E+01
             3.1715E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1098.40527562793        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  1.1072E+00  1.5134E+00  1.1641E+01  1.1144E+00  2.3380E+00  2.0965E+00  4.0572E+00  2.8747E-02  1.8204E+00  3.1485E-01
             9.0554E+00
 PARAMETER:  2.0184E-01  5.1435E-01  2.5545E+00  2.0830E-01  9.4928E-01  8.4027E-01  1.5005E+00 -3.4492E+00  6.9905E-01 -1.0556E+00
             2.3034E+00
 GRADIENT:   2.6384E+00 -2.9860E+00 -1.9602E+00  1.5300E+01 -1.1657E+01  3.4244E+00  2.0941E+00  8.4239E-05  1.5574E+01  1.4895E+00
            -6.3740E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1110.37661232390        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  1.1162E+00  1.8948E+00  1.4496E+01  8.7375E-01  2.5401E+00  2.0595E+00  3.9197E+00  1.0000E-02  5.3843E-01  1.6357E-01
             9.6074E+00
 PARAMETER:  2.0992E-01  7.3909E-01  2.7739E+00 -3.4966E-02  1.0322E+00  8.2249E-01  1.4660E+00 -5.6076E+00 -5.1910E-01 -1.7105E+00
             2.3625E+00
 GRADIENT:  -4.8708E+00  7.5323E+00  1.5982E-01  9.8402E+00  3.7621E+00 -6.4854E+00 -1.2289E+01  0.0000E+00  1.4391E+00  4.0006E-01
             4.4945E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1120.11187726976        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      499
 NPARAMETR:  1.1074E+00  1.8249E+00  1.7741E+01  9.0179E-01  2.5588E+00  2.0947E+00  4.5280E+00  1.0000E-02  2.2821E-01  9.3063E-02
             9.5818E+00
 PARAMETER:  2.0204E-01  7.0152E-01  2.9759E+00 -3.3701E-03  1.0395E+00  8.3939E-01  1.6103E+00 -8.1823E+00 -1.3775E+00 -2.2745E+00
             2.3599E+00
 GRADIENT:  -1.8441E+01 -9.3834E-01  3.5694E-01  3.5346E+00 -1.1895E+00 -2.6167E+01 -3.8064E+01  0.0000E+00 -1.4056E-01  1.1308E-01
            -3.1962E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1125.55582636115        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      676
 NPARAMETR:  1.1605E+00  1.5477E+00  1.6995E+01  1.0007E+00  2.5517E+00  2.2697E+00  5.3255E+00  1.0000E-02  3.6660E-01  1.5180E-01
             9.6238E+00
 PARAMETER:  2.4888E-01  5.3676E-01  2.9329E+00  1.0067E-01  1.0368E+00  9.1964E-01  1.7725E+00 -7.0681E+00 -9.0349E-01 -1.7852E+00
             2.3642E+00
 GRADIENT:   5.3999E+00 -1.9536E+00 -9.9616E-02 -5.3678E+00 -5.4350E-01  1.1874E+01 -5.2270E+00  0.0000E+00 -6.2520E-01  2.7620E-01
            -1.4173E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1126.65959990358        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  1.1470E+00  1.3785E+00  1.8220E+01  1.1189E+00  2.5405E+00  2.2291E+00  5.7422E+00  1.0000E-02  6.7496E-01  2.3775E-01
             9.6029E+00
 PARAMETER:  2.3715E-01  4.2099E-01  3.0025E+00  2.1230E-01  1.0324E+00  9.0159E-01  1.8478E+00 -5.8842E+00 -2.9311E-01 -1.3365E+00
             2.3621E+00
 GRADIENT:  -2.0956E+00  7.1281E-01 -8.3238E-01 -2.6096E+00  1.5591E+00 -1.0485E+00  9.2767E-01  0.0000E+00  9.0354E-01  6.9205E-01
            -1.0534E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1127.56452046495        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1033
 NPARAMETR:  1.1557E+00  1.2056E+00  4.7058E+01  1.1831E+00  2.6640E+00  2.2310E+00  5.9678E+00  1.0000E-02  7.5005E-01  1.2444E-01
             9.6151E+00
 PARAMETER:  2.4467E-01  2.8699E-01  3.9514E+00  2.6814E-01  1.0798E+00  9.0244E-01  1.8864E+00 -8.5597E+00 -1.8761E-01 -1.9839E+00
             2.3633E+00
 GRADIENT:   4.4638E-01 -2.0705E+00 -3.2764E-01 -3.1551E+00  5.4951E+00 -1.2455E+00  3.6198E-01  0.0000E+00 -1.9061E-02  1.7943E-01
             1.7287E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1128.13265338028        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1218
 NPARAMETR:  1.1578E+00  1.1788E+00  2.4334E+02  1.2205E+00  2.6656E+00  2.2456E+00  6.3388E+00  1.0000E-02  8.0717E-01  2.2716E-02
             9.6057E+00
 PARAMETER:  2.4649E-01  2.6453E-01  5.5945E+00  2.9929E-01  1.0804E+00  9.0899E-01  1.9467E+00 -1.0754E+01 -1.1422E-01 -3.6847E+00
             2.3624E+00
 GRADIENT:   1.7422E+00  9.3219E-01 -1.2295E-02 -4.7151E+00 -3.0890E+00  1.4208E+00  1.1305E+01  0.0000E+00  1.2090E-02  5.6590E-03
            -3.5507E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1128.22246974996        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1393
 NPARAMETR:  1.1560E+00  1.1305E+00  5.6236E+02  1.2548E+00  2.7080E+00  2.2444E+00  6.2686E+00  1.0000E-02  8.6595E-01  1.4752E-02
             9.6324E+00
 PARAMETER:  2.4500E-01  2.2268E-01  6.4321E+00  3.2694E-01  1.0962E+00  9.0844E-01  1.9356E+00 -1.0754E+01 -4.3927E-02 -4.1163E+00
             2.3651E+00
 GRADIENT:   3.4391E-01 -6.8908E-02 -1.9608E-02 -4.1427E-01  1.8446E+00  9.3952E-01  4.9482E+00  0.0000E+00  1.4947E-02  2.4886E-03
             4.2112E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1128.32043956550        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1564
 NPARAMETR:  1.1541E+00  1.0954E+00  1.6973E+04  1.2653E+00  2.6995E+00  2.2462E+00  6.4056E+00  1.0000E-02  8.5548E-01  1.0000E-02
             9.6071E+00
 PARAMETER:  2.4333E-01  1.9115E-01  9.8394E+00  3.3531E-01  1.0931E+00  9.0926E-01  1.9572E+00 -1.0754E+01 -5.6089E-02 -1.6924E+01
             2.3625E+00
 GRADIENT:   5.1006E-02  1.9001E-02 -4.5574E-04  6.1847E-01 -6.6333E-01  1.3930E+00  7.7218E+00  0.0000E+00 -8.9400E-01  0.0000E+00
            -9.9390E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1128.33387835976        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1740
 NPARAMETR:  1.1569E+00  1.0936E+00  3.7169E+04  1.2750E+00  2.7035E+00  2.2458E+00  6.3993E+00  1.0000E-02  8.9179E-01  1.0000E-02
             9.6210E+00
 PARAMETER:  2.4577E-01  1.8946E-01  1.0623E+01  3.4298E-01  1.0945E+00  9.0908E-01  1.9562E+00 -1.0754E+01 -1.4522E-02 -1.6924E+01
             2.3640E+00
 GRADIENT:  -6.8510E-01  8.2577E-02 -2.6879E-04 -3.4004E-02  2.3936E-01  2.5600E+00  7.3628E+00  0.0000E+00 -1.0682E-01  0.0000E+00
             2.5521E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1128.35985207932        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1921
 NPARAMETR:  1.1601E+00  1.0821E+00  4.0255E+04  1.2814E+00  2.7026E+00  2.2466E+00  6.4823E+00  1.0000E-02  8.9228E-01  1.0000E-02
             9.6199E+00
 PARAMETER:  2.4847E-01  1.7888E-01  1.0703E+01  3.4792E-01  1.0942E+00  9.0940E-01  1.9691E+00 -1.0754E+01 -1.3980E-02 -1.6924E+01
             2.3638E+00
 GRADIENT:   8.1878E-01  4.2017E-01 -2.5300E-04 -2.5118E-01 -2.0659E-02  2.0293E+00  9.5104E+00  0.0000E+00 -2.0485E-01  0.0000E+00
             1.7155E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1128.38167071827        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2108
 NPARAMETR:  1.1539E+00  1.0550E+00  6.2125E+04  1.2874E+00  2.7005E+00  2.2478E+00  6.5350E+00  1.0000E-02  8.9585E-01  1.0000E-02
             9.6017E+00
 PARAMETER:  2.4310E-01  1.5352E-01  1.1137E+01  3.5266E-01  1.0934E+00  9.0994E-01  1.9772E+00 -1.0754E+01 -9.9839E-03 -1.6924E+01
             2.3619E+00
 GRADIENT:  -3.0151E+00 -1.5262E-02 -1.6779E-04 -2.5521E-01 -1.4874E-01  3.8298E+00  1.0189E+01  0.0000E+00 -4.2637E-01  0.0000E+00
            -2.3957E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1128.39030197814        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2302             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1538E+00  1.0500E+00  7.5696E+04  1.2937E+00  2.7022E+00  2.2462E+00  6.5645E+00  1.0000E-02  9.0779E-01  1.0000E-02
             9.6078E+00
 PARAMETER:  2.4309E-01  1.4880E-01  1.1334E+01  3.5749E-01  1.0941E+00  9.0923E-01  1.9817E+00 -1.0754E+01  3.2564E-03 -1.6924E+01
             2.3626E+00
 GRADIENT:   8.6263E+00  2.0080E+00 -1.4191E-04  9.7305E+00  3.4963E+00  3.3047E+01  1.5466E+02  0.0000E+00 -1.2164E-01  0.0000E+00
             4.4181E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1128.39573931360        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2496
 NPARAMETR:  1.1538E+00  1.0400E+00  4.7226E+04  1.2986E+00  2.7009E+00  2.2462E+00  6.5786E+00  1.0000E-02  9.2015E-01  1.0000E-02
             9.6054E+00
 PARAMETER:  2.4307E-01  1.3925E-01  1.0863E+01  3.6130E-01  1.0936E+00  9.0926E-01  1.9838E+00 -1.0754E+01  1.6777E-02 -1.6924E+01
             2.3623E+00
 GRADIENT:  -7.2190E-01  1.2420E-02 -2.5019E-04 -1.2054E+00  1.0837E-01  1.7623E+00  1.0775E+01  0.0000E+00  1.0467E-01  0.0000E+00
             3.8586E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1128.39904250360        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2693             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1539E+00  1.0308E+00  6.6212E+04  1.3028E+00  2.7020E+00  2.2460E+00  6.6033E+00  1.0000E-02  9.2436E-01  1.0000E-02
             9.6080E+00
 PARAMETER:  2.4319E-01  1.3037E-01  1.1201E+01  3.6450E-01  1.0940E+00  9.0914E-01  1.9876E+00 -1.0754E+01  2.1345E-02 -1.6924E+01
             2.3626E+00
 GRADIENT:   9.7928E+00  1.5084E+00 -1.7592E-04  9.6285E+00  3.6081E+00  3.2626E+01  1.5670E+02  0.0000E+00 -9.5300E-02  0.0000E+00
             4.4497E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1128.40067256048        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2881
 NPARAMETR:  1.1539E+00  1.0209E+00  6.7171E+04  1.3062E+00  2.7013E+00  2.2462E+00  6.6090E+00  1.0000E-02  9.2492E-01  1.0000E-02
             9.6063E+00
 PARAMETER:  2.4314E-01  1.2067E-01  1.1215E+01  3.6716E-01  1.0937E+00  9.0923E-01  1.9884E+00 -1.0754E+01  2.1947E-02 -1.6924E+01
             2.3624E+00
 GRADIENT:  -6.2924E+00 -1.7988E-01 -1.8293E-04  1.1089E+00  9.2633E-02  5.2496E+00  1.0371E+01  0.0000E+00 -6.2918E-01  0.0000E+00
             1.9090E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1128.40270137888        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3069             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1539E+00  1.0216E+00  7.4586E+04  1.3084E+00  2.7013E+00  2.2462E+00  6.6191E+00  1.0000E-02  9.2879E-01  1.0000E-02
             9.6062E+00
 PARAMETER:  2.4311E-01  1.2135E-01  1.1320E+01  3.6877E-01  1.0937E+00  9.0923E-01  1.9900E+00 -1.0754E+01  2.6125E-02 -1.6924E+01
             2.3624E+00
 GRADIENT:   5.7926E+00  1.5067E+00 -1.5697E-04  1.1415E+01  3.3316E+00  3.5435E+01  1.5730E+02  0.0000E+00 -6.5552E-01  0.0000E+00
             4.4617E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1128.40339358848        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3260
 NPARAMETR:  1.1539E+00  1.0208E+00  8.3074E+04  1.3097E+00  2.7009E+00  2.2462E+00  6.6248E+00  1.0000E-02  9.2995E-01  1.0000E-02
             9.6054E+00
 PARAMETER:  2.4316E-01  1.2058E-01  1.1427E+01  3.6982E-01  1.0936E+00  9.0924E-01  1.9908E+00 -1.0754E+01  2.7378E-02 -1.6924E+01
             2.3623E+00
 GRADIENT:  -4.3420E+00  3.2574E-02 -1.5087E-04  4.7421E-01 -6.9664E-02  3.2963E+00  1.0907E+01  0.0000E+00 -2.5944E-01  0.0000E+00
             9.6415E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1128.40433584569        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     3439
 NPARAMETR:  1.1538E+00  1.0175E+00  4.4956E+04  1.3140E+00  2.7003E+00  2.2463E+00  6.6343E+00  1.0000E-02  9.3510E-01  1.0000E-02
             9.6040E+00
 PARAMETER:  2.4303E-01  1.1731E-01  1.0813E+01  3.7307E-01  1.0934E+00  9.0929E-01  1.9922E+00 -1.0754E+01  3.2896E-02 -1.6924E+01
             2.3622E+00
 GRADIENT:  -7.0368E-01  6.7130E-02 -2.8019E-04  4.7626E-01 -2.3646E-01  1.7156E+00  1.0889E+01  0.0000E+00 -2.3473E-01  0.0000E+00
            -4.4242E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1128.40526071060        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3627
 NPARAMETR:  1.1539E+00  1.0175E+00  4.0038E+05  1.3150E+00  2.7008E+00  2.2462E+00  6.6405E+00  1.0000E-02  9.3947E-01  1.0000E-02
             9.6053E+00
 PARAMETER:  2.4303E-01  1.1616E-01  1.2972E+01  3.7334E-01  1.0935E+00  9.0924E-01  1.9930E+00 -1.0754E+01  3.7668E-02 -1.6924E+01
             2.3623E+00
 GRADIENT:  -3.1369E-02 -1.5575E-01 -3.9148E-04 -2.6653E-01 -7.4780E-03  2.9377E-03 -4.5054E-02  0.0000E+00  1.0169E-02  0.0000E+00
            -1.3426E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3627
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5553E-02  3.3587E-02  8.3139E-11 -7.0662E-02  4.2476E-07
 SE:             2.8783E-02  2.4454E-02  3.6784E-10  1.2909E-02  1.1173E-04
 N:                     100         100         100         100         100

 P VAL.:         5.8894E-01  1.6960E-01  8.2118E-01  4.4125E-08  9.9697E-01

 ETASHRINKSD(%)  3.5740E+00  1.8077E+01  1.0000E+02  5.6754E+01  9.9626E+01
 ETASHRINKVR(%)  7.0203E+00  3.2887E+01  1.0000E+02  8.1298E+01  9.9999E+01
 EBVSHRINKSD(%)  5.1066E+00  1.1957E+01  1.0000E+02  6.1763E+01  9.9558E+01
 EBVSHRINKVR(%)  9.9524E+00  2.2485E+01  1.0000E+02  8.5379E+01  9.9998E+01
 RELATIVEINF(%)  8.9903E+01  4.3048E+01  0.0000E+00  8.0280E+00  1.5258E-03
 EPSSHRINKSD(%)  5.5911E+00
 EPSSHRINKVR(%)  1.0870E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1128.4052607105973     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       525.68409905781346     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   128.78
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1128.405       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.15E+00  1.02E+00  3.89E+05  1.31E+00  2.70E+00  2.25E+00  6.64E+00  1.00E-02  9.40E-01  1.00E-02  9.61E+00
 


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
+        1.45E+02
 
 TH 2
+       -6.90E+00  2.22E+01
 
 TH 3
+       -1.13E-07 -4.56E-07  1.69E-14
 
 TH 4
+        1.09E+00  2.10E+01 -2.82E-08  1.51E+02
 
 TH 5
+       -2.33E+00 -2.84E+00  4.99E-09 -1.34E+01  2.83E+01
 
 TH 6
+        7.32E+00 -2.29E-01 -3.24E-08  1.56E+00 -1.91E+00  4.22E+01
 
 TH 7
+        1.00E-01  2.42E+00 -9.32E-11 -8.63E+00  4.32E-01 -2.13E-01  2.59E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.66E+00  7.32E+00 -1.07E-06 -4.74E+01  5.02E+00 -1.69E+00  2.23E+00  0.00E+00  2.20E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.16E+00 -2.04E+00 -4.94E-10 -9.62E+00  6.87E-01  7.68E+00  3.88E-01  0.00E+00  3.48E+00  0.00E+00  9.77E+00
 
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
 #CPUT: Total CPU Time in Seconds,      146.060
Stop Time:
Wed Sep 29 09:14:39 CDT 2021
