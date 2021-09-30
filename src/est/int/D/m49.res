Wed Sep 29 09:01:01 CDT 2021
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
$DATA ../../../../data/int/D/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   24901.5295192380        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1283E+02  3.6666E+02  1.8679E+01  3.3982E+02  2.8344E+02 -1.9502E+03 -8.4268E+02 -8.4042E+01 -1.1157E+03 -6.7103E+02
            -5.2288E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -978.152509309190        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4256E+00  1.7081E+00  9.2475E-01  1.6156E+00  9.0370E-01  4.2689E+00  6.0405E+00  1.0081E+00  2.9507E+00  2.0816E+00
             1.2137E+01
 PARAMETER:  4.5456E-01  6.3540E-01  2.1764E-02  5.7972E-01 -1.2616E-03  1.5514E+00  1.8985E+00  1.0802E-01  1.1820E+00  8.3316E-01
             2.5962E+00
 GRADIENT:   2.1830E+01  1.5559E+01 -3.9916E+01  5.3556E+01 -2.7401E+01  2.4179E+02  1.4099E+02  3.8595E+00  8.3236E+01  5.0989E+01
             5.6421E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1031.75254162761        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      264
 NPARAMETR:  1.2164E+00  1.3050E+00  4.5012E+01  2.6546E+00  3.1479E+00  5.5913E+00  1.2930E+01  8.8763E-01  1.7725E+00  1.2901E+00
             1.2090E+01
 PARAMETER:  2.9593E-01  3.6623E-01  3.9069E+00  1.0763E+00  1.2467E+00  1.8212E+00  2.6596E+00 -1.9199E-02  6.7240E-01  3.5471E-01
             2.5924E+00
 GRADIENT:  -4.2552E+00  2.0261E+01 -7.6171E+00  5.9278E+01  5.6224E+01  2.0383E+02  3.6269E+01  1.9036E-02 -2.9806E+01  1.9610E+01
             5.0821E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1232.31097219603        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.0539E+00  6.3945E-01  5.1792E+01  2.1650E+00  2.1085E+00  1.6804E+00  7.7039E+00  2.0772E+00  3.5985E+00  1.4881E+00
             9.7340E+00
 PARAMETER:  1.5250E-01 -3.4714E-01  4.0472E+00  8.7242E-01  8.4596E-01  6.1905E-01  2.1417E+00  8.3101E-01  1.3805E+00  4.9751E-01
             2.3756E+00
 GRADIENT:  -2.6005E+00  4.4509E+00 -1.1653E+00  3.1821E+00 -2.5650E+01 -1.4514E+01  2.0384E+01  2.3116E+00  4.6633E+01  3.8242E+01
             3.3695E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1297.90315995384        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      626
 NPARAMETR:  9.8358E-01  4.5997E-01  4.6208E+01  1.5523E+00  2.2925E+00  1.8999E+00  6.9714E+00  7.2828E+00  2.1179E+00  6.3448E-01
             8.0134E+00
 PARAMETER:  8.3445E-02 -6.7658E-01  3.9332E+00  5.3971E-01  9.2966E-01  7.4179E-01  2.0418E+00  2.0855E+00  8.5041E-01 -3.5496E-01
             2.1811E+00
 GRADIENT:  -8.0073E+00 -1.1679E+01 -8.7922E+00 -3.7875E+01 -1.1523E+01  2.0184E+01  1.0759E+01  2.1765E+01  1.8391E+01  3.1161E+00
             3.4316E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1306.01745564918        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  9.8194E-01  4.7100E-01  5.3075E+01  1.7865E+00  2.3583E+00  1.7762E+00  7.0036E+00  5.0972E+00  1.9335E+00  5.0419E-01
             7.9323E+00
 PARAMETER:  8.1778E-02 -6.5290E-01  4.0717E+00  6.8027E-01  9.5796E-01  6.7446E-01  2.0464E+00  1.7287E+00  7.5932E-01 -5.8479E-01
             2.1709E+00
 GRADIENT:  -1.0616E+01 -3.5376E+00 -3.1715E+00  6.8571E+00  1.0182E+01 -9.0597E-01  3.4761E+00  4.1187E-01 -2.0021E+00  2.1318E+00
            -2.0885E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1306.49072070983        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  1.0037E+00  6.7310E-01  7.0861E+01  1.5814E+00  2.3369E+00  1.7741E+00  6.2858E+00  4.9214E+00  1.8667E+00  3.7118E-01
             7.9796E+00
 PARAMETER:  1.0369E-01 -2.9587E-01  4.3607E+00  5.5831E-01  9.4882E-01  6.7328E-01  1.9383E+00  1.6936E+00  7.2416E-01 -8.9106E-01
             2.1769E+00
 GRADIENT:   2.5107E+00 -4.1636E+00 -9.9454E-01 -3.4888E+00 -5.3705E+00 -8.8376E-01 -3.7452E+00  6.4466E-02  1.3378E+00  8.7057E-01
             1.2508E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1306.72394347141        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  1.0013E+00  8.8190E-01  1.0234E+02  1.4684E+00  2.3932E+00  1.7804E+00  5.9951E+00  5.0266E+00  1.7840E+00  2.7195E-01
             7.9477E+00
 PARAMETER:  1.0128E-01 -2.5671E-02  4.7283E+00  4.8419E-01  9.7263E-01  6.7685E-01  1.8909E+00  1.7147E+00  6.7885E-01 -1.2021E+00
             2.1729E+00
 GRADIENT:   1.8732E+00 -2.5664E-01 -4.8568E-01 -1.3246E+00 -1.0573E+00 -1.1331E-02 -9.9413E-01  1.0042E-02  5.5133E-01  4.3656E-01
             2.0174E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1308.58563258064        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1335
 NPARAMETR:  9.9571E-01  5.9487E-01  5.9305E+02  1.7522E+00  2.4719E+00  1.7901E+00  6.6671E+00  6.1771E+00  1.9413E+00  1.1678E-01
             7.9799E+00
 PARAMETER:  9.5705E-02 -4.1941E-01  6.4853E+00  6.6087E-01  1.0050E+00  6.8229E-01  1.9972E+00  1.9208E+00  7.6334E-01 -2.0474E+00
             2.1769E+00
 GRADIENT:  -2.8208E+00 -9.2687E-03 -1.8607E-01  9.1776E+00  5.6537E+00  1.5761E+00  4.1519E-01  3.0144E-03 -3.1057E+00  7.5297E-02
             4.6613E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1309.87260083367        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1519             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9932E-01  5.0534E-01  2.1194E+04  1.7519E+00  2.4530E+00  1.7830E+00  7.6478E+00  3.7203E+00  1.9853E+00  2.2638E-02
             7.9585E+00
 PARAMETER:  9.9322E-02 -5.8253E-01  1.0061E+01  6.6068E-01  9.9731E-01  6.7832E-01  2.1344E+00  1.4138E+00  7.8575E-01 -3.6881E+00
             2.1742E+00
 GRADIENT:   6.8626E+00  6.5022E+00 -4.9185E-03  2.4429E+01  5.9608E+00  1.3528E+01  2.4006E+02 -1.5238E-05  1.4364E+01  3.1023E-03
             4.2284E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1309.99728597743        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1681
 NPARAMETR:  9.9865E-01  4.7517E-01  4.4483E+05  1.8209E+00  2.4484E+00  1.7834E+00  7.4735E+00  3.9652E+00  1.9228E+00  1.0000E-02
             7.9486E+00
 PARAMETER:  9.8648E-02 -6.4408E-01  1.3105E+01  6.9936E-01  9.9544E-01  6.7853E-01  2.1114E+00  1.4775E+00  7.5380E-01 -4.7242E+00
             2.1730E+00
 GRADIENT:  -5.5152E+01  2.2024E+00 -2.5924E-04  1.6818E+01 -2.9914E+00  1.1039E+01  1.6321E+01  1.0138E-05 -8.6322E+00  0.0000E+00
             5.9071E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1310.28378019997        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1863
 NPARAMETR:  9.9331E-01  4.3196E-01  5.1345E+05  1.8160E+00  2.4536E+00  1.7840E+00  7.8118E+00  4.1628E+00  1.9958E+00  1.0000E-02
             7.9616E+00
 PARAMETER:  9.3283E-02 -7.3941E-01  1.3249E+01  6.9664E-01  9.9754E-01  6.7886E-01  2.1556E+00  1.5262E+00  7.9103E-01 -4.7100E+00
             2.1746E+00
 GRADIENT:  -2.4057E+01  1.0680E+00 -2.3915E-04 -7.5144E-01  8.2857E-01  4.4179E+00  2.5291E+01  3.5045E-06 -2.4499E+00  0.0000E+00
             6.9639E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1310.37595752543        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2048             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8899E-01  4.0759E-01  6.2323E+05  1.8522E+00  2.4527E+00  1.7811E+00  7.8330E+00  4.2292E+00  2.0045E+00  1.0000E-02
             7.9649E+00
 PARAMETER:  8.8928E-02 -7.9748E-01  1.3443E+01  7.1637E-01  9.9720E-01  6.7726E-01  2.1584E+00  1.5420E+00  7.9539E-01 -4.7100E+00
             2.1750E+00
 GRADIENT:  -2.3978E+01  7.7503E+00 -1.7842E-04  4.7485E+01  4.1268E+00  8.3404E+00  2.4839E+02  3.6820E-05  8.7265E+00  0.0000E+00
             4.6066E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1310.52132760322        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     2215
 NPARAMETR:  9.9651E-01  3.7891E-01  8.1434E+05  1.8933E+00  2.4497E+00  1.7833E+00  8.1092E+00  1.0594E+00  1.9964E+00  1.0000E-02
             7.9470E+00
 PARAMETER:  9.6501E-02 -8.7047E-01  1.3710E+01  7.3830E-01  9.9596E-01  6.7847E-01  2.1930E+00  1.5771E-01  7.9134E-01 -4.7100E+00
             2.1728E+00
 GRADIENT:  -4.9599E+01  2.0221E+00 -1.2487E-04  9.7300E+00 -1.2487E-01 -5.2239E+00  2.7276E+01  1.3798E-04 -6.8939E+00  0.0000E+00
             4.7519E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1310.61431995053        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     2416
 NPARAMETR:  1.0020E+00  3.5140E-01  1.0185E+06  1.9022E+00  2.4516E+00  1.7843E+00  8.2504E+00  1.0586E+00  2.0214E+00  1.0000E-02
             7.9562E+00
 PARAMETER:  1.0198E-01 -9.4582E-01  1.3934E+01  7.4303E-01  9.9672E-01  6.7903E-01  2.2103E+00  1.5692E-01  8.0377E-01 -4.7100E+00
             2.1739E+00
 GRADIENT:  -1.0230E+02  6.4768E+00 -1.4648E-04  4.7669E+01 -3.0075E+00 -2.6131E+01  2.7921E+01  2.2263E-04 -1.0198E+01  0.0000E+00
             8.8979E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1310.63761598167        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2612
 NPARAMETR:  1.0008E+00  3.2173E-01  1.3901E+06  1.9071E+00  2.4525E+00  1.7845E+00  8.3042E+00  1.0530E+00  2.0217E+00  1.0000E-02
             7.9584E+00
 PARAMETER:  1.0078E-01 -1.0340E+00  1.4245E+01  7.4558E-01  9.9712E-01  6.7912E-01  2.2168E+00  1.5162E-01  8.0394E-01 -4.7100E+00
             2.1742E+00
 GRADIENT:  -7.1234E+01  6.8358E+00 -1.2156E-04 -5.4970E+01  1.1383E+01 -3.8565E+00  5.1696E+01 -2.0116E-04 -4.0511E+01  0.0000E+00
             2.3929E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1310.66593228942        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     2815             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0014E+00  3.2196E-01  1.6070E+06  1.9199E+00  2.4512E+00  1.7850E+00  8.3797E+00  1.0522E+00  2.0261E+00  1.0000E-02
             7.9573E+00
 PARAMETER:  1.0143E-01 -1.0333E+00  1.4390E+01  7.5229E-01  9.9657E-01  6.7941E-01  2.2258E+00  1.5088E-01  8.0610E-01 -4.7100E+00
             2.1741E+00
 GRADIENT:  -1.2033E+02 -3.8239E+00 -7.0494E-05  5.1661E+01  8.9126E-01 -1.0562E+01  2.3923E+02  1.4261E-03  7.9095E+00  0.0000E+00
             5.8566E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1310.68252948102        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3012
 NPARAMETR:  9.9945E-01  3.0911E-01  1.1232E+06  1.9243E+00  2.4527E+00  1.7842E+00  8.4599E+00  1.0482E+00  2.0245E+00  1.0000E-02
             7.9594E+00
 PARAMETER:  9.9454E-02 -1.0740E+00  1.4032E+01  7.5455E-01  9.9719E-01  6.7897E-01  2.2353E+00  1.4708E-01  8.0535E-01 -4.7100E+00
             2.1744E+00
 GRADIENT:  -6.0778E+01  1.3971E+00 -1.1707E-04  2.4063E+01 -8.3375E-01  3.7567E+00  2.8217E+01  1.0723E-04 -8.7627E+00  0.0000E+00
             8.3034E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1310.69367765731        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3207
 NPARAMETR:  9.9992E-01  3.0388E-01  1.9683E+06  1.9277E+00  2.4527E+00  1.7843E+00  8.5264E+00  1.0489E+00  2.0242E+00  1.0000E-02
             7.9591E+00
 PARAMETER:  9.9923E-02 -1.0911E+00  1.4593E+01  7.5634E-01  9.9719E-01  6.7905E-01  2.2432E+00  1.4776E-01  8.0517E-01 -4.7100E+00
             2.1743E+00
 GRADIENT:  -8.4705E+01  6.4119E+00 -1.6325E-04 -9.3978E+00  2.5609E-01 -2.7297E+01  4.8185E+01 -5.8879E-04 -2.0016E+01  0.0000E+00
             1.6485E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1310.69865163428        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3399             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9795E-01  2.9861E-01  3.7930E+06  1.9333E+00  2.4529E+00  1.7842E+00  8.5719E+00  1.0458E+00  2.0256E+00  1.0000E-02
             7.9603E+00
 PARAMETER:  9.7949E-02 -1.1086E+00  1.5249E+01  7.5923E-01  9.9728E-01  6.7899E-01  2.2485E+00  1.4482E-01  8.0589E-01 -4.7100E+00
             2.1745E+00
 GRADIENT:  -4.8911E+01  7.3811E+00 -1.5000E-05  5.1297E+01  4.2365E+00  1.7204E+01  2.8825E+02  4.1104E-03  9.6084E+00  0.0000E+00
             4.7983E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1310.70323330264        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3592             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0001E+00  2.9621E-01  2.8251E+06  1.9364E+00  2.4528E+00  1.7845E+00  8.6129E+00  1.0480E+00  2.0249E+00  1.0000E-02
             7.9595E+00
 PARAMETER:  1.0011E-01 -1.1167E+00  1.4954E+01  7.6082E-01  9.9724E-01  6.7913E-01  2.2533E+00  1.4685E-01  8.0553E-01 -4.7100E+00
             2.1744E+00
 GRADIENT:  -1.0520E+02  9.9481E+00 -3.0974E-05  7.3071E+01  1.7708E+00  1.5769E+01  2.9128E+02  2.7128E-03  7.7974E+00  0.0000E+00
             5.5763E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1310.70655304641        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     3792
 NPARAMETR:  9.9953E-01  2.9380E-01  4.7170E+06  1.9366E+00  2.4530E+00  1.7846E+00  8.6465E+00  1.0404E+00  2.0235E+00  1.0000E-02
             7.9602E+00
 PARAMETER:  9.9525E-02 -1.1248E+00  1.5467E+01  7.6094E-01  9.9733E-01  6.7917E-01  2.2572E+00  1.3961E-01  8.0485E-01 -4.7100E+00
             2.1745E+00
 GRADIENT:  -6.9911E+01 -2.4908E+00 -3.9471E-05  3.5695E+00  6.8715E-02  2.0929E+01  2.6512E+01  1.5597E-03 -1.1515E+01  0.0000E+00
             1.5757E+01

0ITERATION NO.:  106    OBJECTIVE VALUE:  -1310.70655304641        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     3822
 NPARAMETR:  9.9975E-01  2.9310E-01  5.5061E+06  1.9441E+00  2.4526E+00  1.7812E+00  8.7059E+00  1.0400E+00  2.0260E+00  1.0000E-02
             7.9590E+00
 PARAMETER:  9.9525E-02 -1.1248E+00  1.5467E+01  7.6094E-01  9.9733E-01  6.7917E-01  2.2572E+00  1.3961E-01  8.0485E-01 -4.7100E+00
             2.1745E+00
 GRADIENT:  -6.9681E-02  2.0550E-02 -2.5876E-03 -5.3240E-01  2.6417E-02  1.5869E-01 -4.1157E-01  9.2264E-02 -6.7901E-02  0.0000E+00
             8.3617E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3822
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0715E-02  6.6546E-02  1.0642E-08 -8.5666E-02 -4.5570E-05
 SE:             2.8149E-02  2.0386E-02  6.2459E-09  1.8614E-02  1.5268E-04
 N:                     100         100         100         100         100

 P VAL.:         7.0345E-01  1.0974E-03  8.8405E-02  4.1812E-06  7.6535E-01

 ETASHRINKSD(%)  5.6979E+00  3.1705E+01  1.0000E+02  3.7642E+01  9.9488E+01
 ETASHRINKVR(%)  1.1071E+01  5.3357E+01  1.0000E+02  6.1115E+01  9.9997E+01
 EBVSHRINKSD(%)  6.7974E+00  3.9588E+01  1.0000E+02  2.3051E+01  9.9467E+01
 EBVSHRINKVR(%)  1.3133E+01  6.3503E+01  1.0000E+02  4.0788E+01  9.9997E+01
 RELATIVEINF(%)  8.6652E+01  2.1516E+01  0.0000E+00  3.4469E+01  6.5199E-04
 EPSSHRINKSD(%)  7.7361E+00
 EPSSHRINKVR(%)  1.4874E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1310.7065530464051     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       343.38280672200563     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   163.68
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    23.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1310.707       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.94E-01  4.72E+06  1.94E+00  2.45E+00  1.78E+00  8.65E+00  1.04E+00  2.02E+00  1.00E-02  7.96E+00
 


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
+        5.87E+02
 
 TH 2
+        2.16E+01  1.14E+02
 
 TH 3
+        2.71E-07  1.62E-08  2.51E-16
 
 TH 4
+       -2.71E+00  3.35E+00 -1.37E-08  6.64E+01
 
 TH 5
+       -3.47E+00 -2.03E+00 -1.53E-08 -7.74E+00  4.84E+01
 
 TH 6
+        7.48E-01  5.02E+00  5.43E-09  2.13E-01  6.52E-02  5.15E+01
 
 TH 7
+        1.81E+00  9.23E+00 -1.05E-09 -4.04E+00  1.99E-01  4.74E-02  1.52E+00
 
 TH 8
+        2.41E+00  2.94E+01  1.11E-07 -5.92E-01 -6.39E+00 -6.23E+01  1.48E-01  4.77E+01
 
 TH 9
+       -7.18E+00 -4.43E+00  2.21E-09 -4.48E+00  1.62E+00 -6.32E-01 -5.93E-01 -3.28E+00  2.36E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -9.26E+00 -2.00E+00 -1.91E-09 -4.41E+00 -2.08E-01  1.24E+00  7.10E-02  3.25E-01  1.22E+00  0.00E+00  1.67E+01
 
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
 #CPUT: Total CPU Time in Seconds,      187.586
Stop Time:
Wed Sep 29 09:04:10 CDT 2021
