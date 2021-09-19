Sat Sep 18 01:45:43 CDT 2021
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
$DATA ../../../../data/int/A3/dat54.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   85.9705474253466        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4561E+01  3.1056E+02  3.6302E+02 -1.2438E+01  2.9690E+02  1.4918E+01 -3.2109E+02 -3.8410E+02 -1.4322E+02 -2.7400E+02
            -6.7582E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2369.28673443828        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1010E+00  9.3133E-01  7.4715E-01  1.1129E+00  8.3368E-01  8.6037E-01  9.5040E-01  1.4563E+00  9.8587E-01  8.3219E-01
             5.1028E+00
 PARAMETER:  1.9619E-01  2.8856E-02 -1.9149E-01  2.0695E-01 -8.1901E-02 -5.0396E-02  4.9127E-02  4.7591E-01  8.5772E-02 -8.3697E-02
             1.7298E+00
 GRADIENT:   5.4520E+01 -2.2186E+01 -3.1450E+01  1.1042E+01  1.5846E+01 -1.8763E+01  1.1452E+01  2.5914E+01  2.2538E+01  2.4406E+01
             7.7235E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2425.08154071591        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0490E+00  3.6767E-01  2.1610E-01  1.5614E+00  2.5357E-01  1.0225E+00  1.1674E+00  1.2028E+00  1.9354E+00  5.7642E-01
             4.3148E+00
 PARAMETER:  1.4786E-01 -9.0058E-01 -1.4320E+00  5.4556E-01 -1.2721E+00  1.2221E-01  2.5474E-01  2.8466E-01  7.6031E-01 -4.5092E-01
             1.5620E+00
 GRADIENT:  -3.2200E+01  1.8977E+02  5.0742E+01  1.7921E+02 -2.2017E+02  1.9444E+01 -6.6534E-01  3.0126E+01  3.1127E+01  8.4945E+00
             6.9685E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2539.59459989074        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.5560E-01  1.9951E-01  7.8272E-02  1.8204E+00  1.3927E-01  1.0279E+00  1.9922E+00  7.3058E-01  5.6236E+00  3.4743E-01
             2.8908E+00
 PARAMETER:  5.4586E-02 -1.5119E+00 -2.4476E+00  6.9906E-01 -1.8713E+00  1.2752E-01  7.8923E-01 -2.1392E-01  1.8270E+00 -9.5720E-01
             1.1615E+00
 GRADIENT:   4.3605E+01  1.1513E+02  2.1508E+00  1.0353E+01 -1.1414E+02  4.6159E+01  6.8970E+01 -2.9084E+01  2.9336E+01 -3.2468E+01
             1.5788E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2594.17485030364        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  9.2847E-01  1.5243E-01  6.8629E-02  1.9226E+00  1.2728E-01  8.4563E-01  1.0996E+00  6.8305E-01  5.4835E+00  7.4471E-01
             2.5504E+00
 PARAMETER:  2.5783E-02 -1.7810E+00 -2.5790E+00  7.5370E-01 -1.9614E+00 -6.7669E-02  1.9497E-01 -2.8119E-01  1.8017E+00 -1.9477E-01
             1.0362E+00
 GRADIENT:   4.7098E+00  2.0300E+01  3.3604E+01  1.1634E+01 -2.2316E+01  6.2628E+00 -1.6781E+01 -3.1273E+00 -1.4855E+01  6.7188E+00
            -6.6126E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2594.55029805628        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      385
 NPARAMETR:  9.3037E-01  1.5113E-01  6.7460E-02  1.9235E+00  1.2639E-01  8.4073E-01  1.1194E+00  6.7462E-01  5.5135E+00  7.4577E-01
             2.5589E+00
 PARAMETER:  2.7831E-02 -1.7896E+00 -2.5962E+00  7.5415E-01 -1.9684E+00 -7.3480E-02  2.1281E-01 -2.9361E-01  1.8072E+00 -1.9333E-01
             1.0396E+00
 GRADIENT:   5.5488E+00  1.8293E+01  2.8376E+01  1.2328E+01 -2.3140E+01  4.4901E+00 -1.3877E+01 -2.9010E+00 -1.3651E+01  5.3801E+00
            -5.9946E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2597.94764465130        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      566
 NPARAMETR:  9.2389E-01  1.5194E-01  6.9004E-02  1.8645E+00  1.3064E-01  8.2729E-01  1.1780E+00  6.8459E-01  5.9274E+00  7.1805E-01
             2.6236E+00
 PARAMETER:  2.0839E-02 -1.7843E+00 -2.5736E+00  7.2299E-01 -1.9353E+00 -8.9600E-02  2.6381E-01 -2.7893E-01  1.8796E+00 -2.3122E-01
             1.0646E+00
 GRADIENT:  -9.3429E+00  8.2967E-01  2.0188E+00  4.4558E+00 -3.6655E+00  1.3406E+00 -1.1533E+00 -2.1537E+00  6.4158E-01  3.3211E+00
             3.1378E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2598.16695301763        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      744
 NPARAMETR:  9.2528E-01  1.5265E-01  6.9568E-02  1.8500E+00  1.3132E-01  8.2669E-01  1.1794E+00  7.2775E-01  5.9038E+00  6.9892E-01
             2.6206E+00
 PARAMETER:  2.2342E-02 -1.7796E+00 -2.5655E+00  7.1520E-01 -1.9301E+00 -9.0331E-02  2.6502E-01 -2.1779E-01  1.8756E+00 -2.5822E-01
             1.0634E+00
 GRADIENT:  -5.9400E+00  1.5762E-01  3.2529E+00  3.6719E+00 -2.2677E+00  1.7258E+00 -1.4966E+00 -2.1120E+00  1.2720E+00  1.3718E-01
            -2.0528E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2598.18415003492        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      926
 NPARAMETR:  9.2546E-01  1.5268E-01  6.9563E-02  1.8486E+00  1.3135E-01  8.2645E-01  1.1800E+00  7.3147E-01  5.9000E+00  6.9836E-01
             2.6204E+00
 PARAMETER:  2.2538E-02 -1.7794E+00 -2.5655E+00  7.1445E-01 -1.9299E+00 -9.0610E-02  2.6553E-01 -2.1270E-01  1.8749E+00 -2.5902E-01
             1.0633E+00
 GRADIENT:  -5.7863E+00  1.9290E-01  3.1847E+00  3.6471E+00 -2.2167E+00  1.6849E+00 -1.4345E+00 -2.0714E+00  1.2417E+00  1.0393E-01
            -1.9731E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2598.28895441026        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  9.2542E-01  1.5254E-01  6.9436E-02  1.8493E+00  1.3119E-01  8.2641E-01  1.1802E+00  8.0897E-01  5.9049E+00  6.9827E-01
             2.6218E+00
 PARAMETER:  2.2488E-02 -1.7804E+00 -2.5674E+00  7.1481E-01 -1.9311E+00 -9.0661E-02  2.6567E-01 -1.1200E-01  1.8758E+00 -2.5915E-01
             1.0639E+00
 GRADIENT:  -6.7676E+00  6.9251E+00  1.5452E+01  5.9308E+00  5.8906E+01  3.2077E+00 -9.9316E-01  8.2296E-02  1.3674E+01  2.2886E+00
             6.5795E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2598.28897298596        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     1104
 NPARAMETR:  9.2542E-01  1.5253E-01  6.9403E-02  1.8493E+00  1.3115E-01  8.2641E-01  1.1802E+00  8.1127E-01  5.9043E+00  6.9827E-01
             2.6218E+00
 PARAMETER:  2.2487E-02 -1.7804E+00 -2.5678E+00  7.1480E-01 -1.9314E+00 -9.0662E-02  2.6567E-01 -1.0915E-01  1.8757E+00 -2.5915E-01
             1.0639E+00
 GRADIENT:  -7.0289E+00  7.0496E+00  1.5350E+01  5.9506E+00  5.8503E+01  3.2289E+00 -9.9150E-01  1.5390E-01  1.3670E+01  2.3052E+00
             6.6019E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2598.29414656741        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  9.2542E-01  1.5249E-01  6.9179E-02  1.8492E+00  1.3106E-01  8.2641E-01  1.1802E+00  7.9924E-01  5.9004E+00  6.9826E-01
             2.6213E+00
 PARAMETER:  2.2489E-02 -1.7806E+00 -2.5711E+00  7.1474E-01 -1.9321E+00 -9.0666E-02  2.6568E-01 -1.2410E-01  1.8750E+00 -2.5917E-01
             1.0637E+00
 GRADIENT:  -1.1120E+01  3.5816E-01  1.7992E+00  3.9073E+00 -6.6086E+00  2.6744E+00 -1.2462E+00 -2.5927E-01  2.2247E+00  1.4337E+00
             3.3743E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2598.30689661475        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1430
 NPARAMETR:  9.2543E-01  1.5253E-01  6.9041E-02  1.8488E+00  1.3132E-01  8.2641E-01  1.1802E+00  8.1109E-01  5.8976E+00  6.9825E-01
             2.6202E+00
 PARAMETER:  2.2501E-02 -1.7804E+00 -2.5731E+00  7.1456E-01 -1.9301E+00 -9.0669E-02  2.6569E-01 -1.0938E-01  1.8745E+00 -2.5918E-01
             1.0633E+00
 GRADIENT:  -1.3055E+01 -4.9222E-01 -1.0227E+00  3.9860E+00  2.3961E-01  2.6133E+00 -1.3681E+00 -2.5510E-02  2.6357E+00  1.7543E+00
             3.2862E+00

0ITERATION NO.:   62    OBJECTIVE VALUE:  -2598.30692650899        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1493
 NPARAMETR:  9.2538E-01  1.5239E-01  6.9132E-02  1.8482E+00  1.3120E-01  8.2636E-01  1.1800E+00  8.1129E-01  5.9031E+00  6.9816E-01
             2.6216E+00
 PARAMETER:  2.2501E-02 -1.7804E+00 -2.5730E+00  7.1456E-01 -1.9301E+00 -9.0669E-02  2.6569E-01 -1.0924E-01  1.8746E+00 -2.5918E-01
             1.0632E+00
 GRADIENT:   6.4697E+03  3.8225E+02 -1.2902E+02  4.6074E+02  3.6144E+02  3.2340E+03  2.4669E+03 -2.0874E-02 -3.5305E+02  1.2579E+03
            -6.2258E+02
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         3.3         3.3         2.9         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1493
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1142E-03  2.1382E-02  3.4098E-02 -9.4157E-03  2.5590E-02
 SE:             2.6227E-02  2.4764E-02  1.5858E-02  2.7564E-02  2.4897E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6611E-01  3.8790E-01  3.1538E-02  7.3266E-01  3.0402E-01

 ETASHRINKSD(%)  1.2135E+01  1.7036E+01  4.6874E+01  7.6579E+00  1.6592E+01
 ETASHRINKVR(%)  2.2798E+01  3.1170E+01  7.1776E+01  1.4729E+01  3.0432E+01
 EBVSHRINKSD(%)  1.0070E+01  1.8648E+01  4.9856E+01  4.7450E+00  1.6073E+01
 EBVSHRINKVR(%)  1.9126E+01  3.3818E+01  7.4855E+01  9.2649E+00  2.9562E+01
 RELATIVEINF(%)  7.3194E+01  5.1216E+01  1.4134E+01  8.4808E+01  4.0925E+01
 EPSSHRINKSD(%)  2.0002E+01
 EPSSHRINKVR(%)  3.6003E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2598.3069265089866     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -944.21756674057588     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    52.10
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    22.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2598.307       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.25E-01  1.53E-01  6.90E-02  1.85E+00  1.31E-01  8.26E-01  1.18E+00  8.11E-01  5.90E+00  6.98E-01  2.62E+00
 


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
+        1.90E+07
 
 TH 2
+        3.32E+03  2.38E+06
 
 TH 3
+       -6.46E+03 -1.89E+04  5.34E+06
 
 TH 4
+        3.09E+03  7.31E+02 -3.05E+03  9.39E+04
 
 TH 5
+        1.70E+03  5.96E+04 -3.93E+04  3.95E+02  5.04E+06
 
 TH 6
+        3.85E+04  1.15E+03 -2.80E+03  1.04E+03  6.94E+02  2.37E+07
 
 TH 7
+        9.59E+03  7.88E+03 -2.59E+04  3.37E+03  5.23E+03  3.25E+03  1.68E+06
 
 TH 8
+        2.21E+04  3.56E+02 -1.39E+03  5.60E+02  5.20E+01  2.21E+07  1.67E+03  2.06E+07
 
 TH 9
+       -1.17E+02 -7.40E+02  1.04E+03 -3.51E+01 -2.42E+02 -4.22E+01 -3.17E+02 -1.67E+01  1.38E+03
 
 TH10
+        2.99E+04  5.45E+03 -1.49E+04  5.22E+03  4.32E+03  1.01E+04  1.56E+04  5.61E+03 -1.99E+02  4.99E+06
 
 TH11
+       -6.91E+02 -4.15E+03  1.22E+04 -1.98E+02 -3.20E+03 -2.20E+02 -1.77E+03 -9.40E+01  1.47E+02 -1.04E+03  2.19E+04
 
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
 #CPUT: Total CPU Time in Seconds,       74.280
Stop Time:
Sat Sep 18 01:46:59 CDT 2021
