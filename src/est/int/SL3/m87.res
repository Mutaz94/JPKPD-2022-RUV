Sat Sep 25 02:47:30 CDT 2021
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
$DATA ../../../../data/int/SL3/dat87.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF DATA RECS IN DATA SET:      986
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

 TOT. NO. OF OBS RECS:      886
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1203.44108986238        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7717E+01 -5.8039E+01  3.0040E+02  5.0221E+01  2.3298E+02 -3.6450E-01 -1.5497E+02 -7.6430E+02 -2.5256E+02 -7.3813E+01
            -8.8938E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2365.49420381862        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0878E+00  1.4058E+00  9.6269E-01  8.9622E-01  1.0811E+00  9.4139E-01  1.0322E+00  1.0132E+00  9.7627E-01  8.4843E-01
             5.2023E+00
 PARAMETER:  1.8416E-01  4.4062E-01  6.1973E-02 -9.5727E-03  1.7801E-01  3.9598E-02  1.3171E-01  1.1307E-01  7.5985E-02 -6.4371E-02
             1.7491E+00
 GRADIENT:   1.4435E+01  1.5180E+01 -7.3173E+00  6.9385E+00 -1.8022E+01 -8.2347E+00  1.8657E+01  5.5606E+00  8.6788E+00  7.9066E+00
             7.4695E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2510.73900627610        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0264E+00  1.4172E+00  1.6515E+00  8.1690E-01  1.3502E+00  9.5458E-01  8.4017E-01  7.1990E-01  1.6077E+00  5.8298E-01
             3.9014E+00
 PARAMETER:  1.2604E-01  4.4869E-01  6.0168E-01 -1.0224E-01  4.0026E-01  5.3511E-02 -7.4148E-02 -2.2865E-01  5.7481E-01 -4.3961E-01
             1.4613E+00
 GRADIENT:  -4.7264E+01 -7.0838E+01  1.5392E+01 -2.8949E+00  4.6358E+01 -1.0484E+01 -2.9806E+01 -5.2162E+00  5.3362E+01 -4.1568E+01
             3.8440E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2603.30531048052        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0410E+00  1.7582E+00  1.0381E+00  6.1264E-01  1.4236E+00  1.0043E+00  7.6971E-01  7.2916E-01  1.3783E+00  1.2675E+00
             3.1122E+00
 PARAMETER:  1.4021E-01  6.6428E-01  1.3743E-01 -3.8997E-01  4.5322E-01  1.0432E-01 -1.6174E-01 -2.1587E-01  4.2082E-01  3.3704E-01
             1.2353E+00
 GRADIENT:   1.3120E+01  7.3379E+00 -7.3774E+00  1.9952E+01  6.3835E+00  4.7148E+00  3.4826E+00  1.2001E+00  5.0741E-01 -5.8668E-01
             3.5670E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2608.24224457059        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0339E+00  1.9797E+00  1.3686E+00  4.5220E-01  1.6405E+00  9.8460E-01  6.6786E-01  1.8961E-01  1.7796E+00  1.4423E+00
             3.0124E+00
 PARAMETER:  1.3329E-01  7.8296E-01  4.1379E-01 -6.9362E-01  5.9498E-01  8.4476E-02 -3.0368E-01 -1.5628E+00  6.7640E-01  4.6627E-01
             1.2027E+00
 GRADIENT:   2.8152E+00  1.4209E+00  1.4332E+00  5.5968E+00 -8.6067E-01 -2.7799E+00 -4.8224E+00 -9.2358E-03  3.7120E+00 -3.4935E+00
            -2.2982E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2610.99215021387        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0329E+00  2.2004E+00  8.2396E-01  2.9253E-01  1.7047E+00  9.9385E-01  6.8399E-01  1.4853E-02  2.1421E+00  1.4903E+00
             3.0300E+00
 PARAMETER:  1.3233E-01  8.8862E-01 -9.3633E-02 -1.1292E+00  6.3339E-01  9.3827E-02 -2.7981E-01 -4.1095E+00  8.6178E-01  4.9899E-01
             1.2086E+00
 GRADIENT:   4.9283E-01 -1.0089E+01 -2.5212E-01 -3.2437E+00 -5.4176E+00  7.3748E-01  1.7201E+00  1.5778E-04  1.1753E-01 -2.0183E+00
            -1.4418E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2612.24126176984        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  1.0374E+00  2.3632E+00  6.5160E-01  2.1014E-01  1.8132E+00  9.9455E-01  6.6931E-01  1.0000E-02  2.5432E+00  1.5705E+00
             3.0301E+00
 PARAMETER:  1.3667E-01  9.6000E-01 -3.2833E-01 -1.4600E+00  6.9510E-01  9.4536E-02 -3.0150E-01 -6.1178E+00  1.0334E+00  5.5137E-01
             1.2086E+00
 GRADIENT:   1.0814E+00  1.1924E+01 -1.0559E+00  4.8561E+00  4.4049E+00 -3.2463E-01 -7.9498E-01  0.0000E+00  4.6660E-01 -1.5226E-01
            -1.3852E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2612.82651873170        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  1.0363E+00  2.4978E+00  4.5909E-01  1.1334E-01  1.8791E+00  9.9481E-01  6.6383E-01  1.0000E-02  3.3427E+00  1.6375E+00
             3.0257E+00
 PARAMETER:  1.3561E-01  1.0154E+00 -6.7851E-01 -2.0774E+00  7.3079E-01  9.4793E-02 -3.0973E-01 -1.0680E+01  1.3068E+00  5.9316E-01
             1.2072E+00
 GRADIENT:  -6.0334E-01  1.8369E+00 -1.5887E-01  2.2238E-01  3.1496E-01 -1.6115E-01  3.6544E-01  0.0000E+00 -1.4147E-01 -4.8309E-02
            -8.8322E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2612.84461528790        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  1.0363E+00  2.5231E+00  4.5579E-01  9.6036E-02  1.8943E+00  9.9487E-01  6.6110E-01  1.0000E-02  3.5616E+00  1.6530E+00
             3.0257E+00
 PARAMETER:  1.3565E-01  1.0255E+00 -6.8571E-01 -2.2430E+00  7.3884E-01  9.4853E-02 -3.1385E-01 -1.1930E+01  1.3702E+00  6.0257E-01
             1.2071E+00
 GRADIENT:  -5.2858E-01  1.5922E+00 -7.9992E-02  1.4027E-01 -8.9912E-03 -1.6243E-01 -2.0144E-01  0.0000E+00 -6.2148E-01  6.3009E-02
            -1.5838E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2613.28794083247        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1070
 NPARAMETR:  1.0369E+00  2.6035E+00  3.5385E-01  4.3762E-02  1.9209E+00  9.9427E-01  6.6583E-01  1.0000E-02  5.3891E+00  1.6661E+00
             3.0162E+00
 PARAMETER:  1.3621E-01  1.0568E+00 -9.3889E-01 -3.0290E+00  7.5281E-01  9.4256E-02 -3.0672E-01 -1.8202E+01  1.7844E+00  6.1052E-01
             1.2040E+00
 GRADIENT:   4.3439E+00 -4.2147E+01 -4.6170E-01  8.8108E+00 -1.3300E+01  3.1913E+00 -8.1173E+00  0.0000E+00  2.9132E+01 -1.2862E+01
            -1.1744E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2614.23586741282        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1247
 NPARAMETR:  1.0362E+00  2.5965E+00  3.9296E-01  4.5360E-02  1.9262E+00  9.9430E-01  6.5843E-01  1.0000E-02  5.1794E+00  1.6802E+00
             3.0209E+00
 PARAMETER:  1.3558E-01  1.0541E+00 -8.3403E-01 -2.9931E+00  7.5553E-01  9.4283E-02 -3.1789E-01 -1.7900E+01  1.7447E+00  6.1890E-01
             1.2056E+00
 GRADIENT:  -1.8535E-01  4.8810E+00 -2.0136E-01 -2.2888E-01 -8.5137E-02 -3.6509E-01  5.5206E-01  0.0000E+00 -1.5933E+00  6.8061E-01
            -7.1974E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2617.36361101218        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1423
 NPARAMETR:  1.0359E+00  2.6574E+00  5.7048E-01  1.0000E-02  1.9758E+00  1.0008E+00  6.3922E-01  1.0000E-02  1.0565E+01  1.7418E+00
             3.0303E+00
 PARAMETER:  1.3531E-01  1.0773E+00 -4.6127E-01 -4.7165E+00  7.8099E-01  1.0079E-01 -3.4751E-01 -3.1831E+01  2.4575E+00  6.5492E-01
             1.2087E+00
 GRADIENT:  -2.0875E+00  2.4332E+01  2.7581E-01  0.0000E+00  1.0185E+01  1.1499E+00 -1.5943E+00  0.0000E+00 -2.8334E+00  6.2137E+00
             6.5604E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2617.64259493223        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1600
 NPARAMETR:  1.0368E+00  2.6473E+00  6.0783E-01  1.0000E-02  1.9456E+00  9.9788E-01  6.4436E-01  1.0000E-02  1.0778E+01  1.7114E+00
             3.0263E+00
 PARAMETER:  1.3616E-01  1.0735E+00 -3.9786E-01 -4.8071E+00  7.6555E-01  9.7880E-02 -3.3950E-01 -3.2566E+01  2.4775E+00  6.3730E-01
             1.2073E+00
 GRADIENT:   4.7269E-01  9.2002E+00  5.4029E-01  0.0000E+00  2.2679E+00  3.3481E-01  5.9468E-01  0.0000E+00 -2.2208E+00  2.5967E+00
             3.2277E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2617.69475152229        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1779
 NPARAMETR:  1.0365E+00  2.6431E+00  5.7625E-01  1.0000E-02  1.9344E+00  9.9679E-01  6.4707E-01  1.0000E-02  1.0816E+01  1.6930E+00
             3.0226E+00
 PARAMETER:  1.3585E-01  1.0720E+00 -4.5121E-01 -4.8334E+00  7.5980E-01  9.6783E-02 -3.3530E-01 -3.2779E+01  2.4810E+00  6.2650E-01
             1.2061E+00
 GRADIENT:   2.6433E-01 -5.8014E-01  1.8129E-01  0.0000E+00 -7.1309E-01 -1.3905E-02 -2.4747E-01  0.0000E+00  3.1162E-01  4.4245E-02
            -4.2466E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2617.69829007998        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1954
 NPARAMETR:  1.0364E+00  2.6435E+00  5.5724E-01  1.0000E-02  1.9368E+00  9.9690E-01  6.4693E-01  1.0000E-02  1.0811E+01  1.6932E+00
             3.0226E+00
 PARAMETER:  1.3573E-01  1.0721E+00 -4.8475E-01 -4.8300E+00  7.6104E-01  9.6894E-02 -3.3552E-01 -3.2755E+01  2.4806E+00  6.2662E-01
             1.2061E+00
 GRADIENT:  -6.3910E-03 -2.9343E-01 -3.8125E-03  0.0000E+00  2.0460E-03  4.5025E-04 -2.1547E-01  0.0000E+00  2.4060E-01 -1.2927E-04
            -2.4668E-02

0ITERATION NO.:   71    OBJECTIVE VALUE:  -2617.69829007998        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1976
 NPARAMETR:  1.0364E+00  2.6435E+00  5.5724E-01  1.0000E-02  1.9368E+00  9.9690E-01  6.4693E-01  1.0000E-02  1.0811E+01  1.6932E+00
             3.0226E+00
 PARAMETER:  1.3573E-01  1.0721E+00 -4.8475E-01 -4.8300E+00  7.6104E-01  9.6894E-02 -3.3552E-01 -3.2755E+01  2.4806E+00  6.2662E-01
             1.2061E+00
 GRADIENT:  -6.3910E-03 -2.9343E-01 -3.8125E-03  0.0000E+00  2.0460E-03  4.5025E-04 -2.1547E-01  0.0000E+00  2.4060E-01 -1.2927E-04
            -2.4668E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1976
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0944E-03 -1.1469E-02  1.7961E-06  6.4729E-03 -1.5954E-02
 SE:             2.9166E-02  2.7443E-02  1.7750E-06  5.3681E-03  2.6559E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7007E-01  6.7601E-01  3.1158E-01  2.2790E-01  5.4805E-01

 ETASHRINKSD(%)  2.2886E+00  8.0610E+00  9.9994E+01  8.2016E+01  1.1022E+01
 ETASHRINKVR(%)  4.5247E+00  1.5472E+01  1.0000E+02  9.6766E+01  2.0830E+01
 EBVSHRINKSD(%)  2.0779E+00  7.4220E+00  9.9949E+01  8.6745E+01  9.7319E+00
 EBVSHRINKVR(%)  4.1125E+00  1.4293E+01  1.0000E+02  9.8243E+01  1.8517E+01
 RELATIVEINF(%)  9.5778E+01  4.8437E+01  2.5747E-05  1.0038E+00  8.0692E+01
 EPSSHRINKSD(%)  1.5145E+01
 EPSSHRINKVR(%)  2.7997E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          886
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1628.3590808386800     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2617.6982900799776     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -989.33920924129757     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    50.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2617.698       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.64E+00  5.57E-01  1.00E-02  1.94E+00  9.97E-01  6.47E-01  1.00E-02  1.08E+01  1.69E+00  3.02E+00
 


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
+        9.98E+02
 
 TH 2
+       -1.40E+01  2.69E+04
 
 TH 3
+       -1.77E+00  5.07E+01  1.22E+01
 
 TH 4
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 5
+       -3.42E+00 -1.89E+01  1.12E+00  0.00E+00  8.34E+01
 
 TH 6
+        3.38E+00 -5.17E+00 -1.30E+00  0.00E+00 -1.79E+00  1.84E+02
 
 TH 7
+        4.84E+01 -3.01E+02  5.88E+02  0.00E+00 -3.95E+01  5.15E+00  4.54E+06
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.45E-01  1.15E+01 -5.04E+00  0.00E+00  3.15E-01 -7.93E-02  2.82E+01  0.00E+00  2.97E+02
 
 TH10
+       -7.63E-01  4.87E-01  2.65E-01  0.00E+00 -5.54E+00  6.69E-01  4.47E+00  0.00E+00 -4.63E-02  4.28E+01
 
 TH11
+       -1.94E+01 -1.57E+01  7.30E+00  0.00E+00  2.01E+00  2.48E+00 -4.07E+01  0.00E+00  1.75E-01  5.90E+00  1.62E+04
 
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
 #CPUT: Total CPU Time in Seconds,       63.695
Stop Time:
Sat Sep 25 02:48:35 CDT 2021
