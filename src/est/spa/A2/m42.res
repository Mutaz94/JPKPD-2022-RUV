Sat Sep 18 09:51:19 CDT 2021
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
$DATA ../../../../data/spa/A2/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -961.737782434081        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.1896E+01  3.1911E+01  5.1413E+01 -3.0516E+01  4.3285E+01  8.1278E+00 -1.0934E+01 -1.2677E+01 -5.1030E+01 -3.3193E+01
            -1.3507E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1440.22505292434        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0382E+00  9.5204E-01  1.0965E+00  1.0713E+00  1.0169E+00  8.8254E-01  8.7841E-01  8.8622E-01  1.0013E+00  7.4989E-01
             2.5551E+00
 PARAMETER:  1.3744E-01  5.0850E-02  1.9216E-01  1.6890E-01  1.1671E-01 -2.4954E-02 -2.9642E-02 -2.0795E-02  1.0135E-01 -1.8782E-01
             1.0381E+00
 GRADIENT:  -3.1824E+01  9.9699E+00  1.2260E+00  6.0465E+00  8.1644E+00 -2.8042E+01  4.3563E+00  3.5023E+00 -5.3581E+00  3.8230E+00
            -3.9489E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1445.71198496317        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0415E+00  8.4586E-01  8.5662E-01  1.1400E+00  8.4735E-01  9.7639E-01  6.1250E-01  6.5346E-01  1.0378E+00  5.2112E-01
             2.5693E+00
 PARAMETER:  1.4066E-01 -6.7404E-02 -5.4767E-02  2.3099E-01 -6.5647E-02  7.6112E-02 -3.9020E-01 -3.2547E-01  1.3711E-01 -5.5177E-01
             1.0436E+00
 GRADIENT:  -2.3744E+01  1.8234E+01 -6.4895E+00  3.8506E+01  9.8751E+00  7.2485E+00 -2.3844E-01  2.0721E+00 -3.6898E+00  7.2369E-01
            -3.7071E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1449.89437092052        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.0516E+00  6.2975E-01  8.2150E-01  1.2476E+00  7.4451E-01  9.4736E-01  6.7498E-01  1.4142E-01  9.4256E-01  4.2690E-01
             2.7656E+00
 PARAMETER:  1.5028E-01 -3.6243E-01 -9.6620E-02  3.2118E-01 -1.9503E-01  4.5919E-02 -2.9308E-01 -1.8560E+00  4.0844E-02 -7.5121E-01
             1.1172E+00
 GRADIENT:  -1.0503E+00  8.5843E+00 -1.0122E+00  2.0777E+01 -2.2684E+00 -5.0350E-01  1.0215E-01  1.2162E-01 -4.0857E-01  1.0065E-01
            -1.8524E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1450.79882737723        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  1.0471E+00  3.7161E-01  8.4578E-01  1.3889E+00  6.8016E-01  9.4493E-01  2.2549E-01  1.2473E-02  8.7967E-01  4.6872E-01
             2.7519E+00
 PARAMETER:  1.4602E-01 -8.8991E-01 -6.7500E-02  4.2850E-01 -2.8543E-01  4.3356E-02 -1.3895E+00 -4.2842E+00 -2.8206E-02 -6.5775E-01
             1.1123E+00
 GRADIENT:   9.5995E-01  3.9747E+00  6.0260E+00  2.0857E+01 -8.3212E+00 -2.1028E-01  2.3214E-02  9.4320E-04  4.3003E+00 -3.9329E-01
             8.1168E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1451.06614251281        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  1.0447E+00  2.6331E-01  6.9787E-01  1.4091E+00  5.6759E-01  9.4900E-01  6.1487E-02  1.0000E-02  8.4364E-01  5.4925E-01
             2.6583E+00
 PARAMETER:  1.4370E-01 -1.2344E+00 -2.5972E-01  4.4298E-01 -4.6636E-01  4.7651E-02 -2.6889E+00 -6.3885E+00 -7.0024E-02 -4.9920E-01
             1.0777E+00
 GRADIENT:  -1.7354E+00  2.5957E+00  4.2123E+00  7.7123E+00 -7.1277E+00 -1.0779E-01 -1.0443E-04  0.0000E+00 -4.2729E-02  3.0219E-01
            -3.7921E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1451.07141235224        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  1.0446E+00  2.1841E-01  6.3675E-01  1.4144E+00  5.2167E-01  9.5133E-01  2.8807E-02  1.0000E-02  8.3571E-01  5.5831E-01
             2.6293E+00
 PARAMETER:  1.4359E-01 -1.4214E+00 -3.5138E-01  4.4668E-01 -5.5071E-01  5.0102E-02 -3.4471E+00 -7.5533E+00 -7.9478E-02 -4.8284E-01
             1.0667E+00
 GRADIENT:  -8.4664E-01  2.3899E+00  4.2769E+00  6.5442E+00 -7.8465E+00  1.7302E-01 -7.7235E-05  0.0000E+00 -2.2443E-02  1.7844E-01
            -1.6972E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1451.07704694027        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  1.0441E+00  1.8043E-01  5.9427E-01  1.4180E+00  4.8897E-01  9.5266E-01  1.2709E-02  1.0000E-02  8.2877E-01  5.6456E-01
             2.6074E+00
 PARAMETER:  1.4312E-01 -1.6124E+00 -4.2042E-01  4.4927E-01 -6.1545E-01  5.1507E-02 -4.2654E+00 -8.7638E+00 -8.7815E-02 -4.7172E-01
             1.0584E+00
 GRADIENT:  -3.3733E-01  2.0277E+00  4.3106E+00  5.2582E+00 -8.0243E+00  3.1239E-01 -1.5164E-05  0.0000E+00 -1.0605E-01  9.4878E-02
            -3.4379E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1451.09187573476        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.0431E+00  1.3644E-01  5.5335E-01  1.4214E+00  4.5706E-01  9.5356E-01  1.0000E-02  1.0000E-02  8.2141E-01  5.7048E-01
             2.5851E+00
 PARAMETER:  1.4217E-01 -1.8919E+00 -4.9177E-01  4.5163E-01 -6.8294E-01  5.2451E-02 -5.5233E+00 -1.0567E+01 -9.6733E-02 -4.6127E-01
             1.0498E+00
 GRADIENT:   2.0285E-01  1.3948E+00  3.5008E+00  2.8592E+00 -6.7159E+00  3.9400E-01  0.0000E+00  0.0000E+00 -1.2771E-01  6.8103E-03
             1.8920E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1451.12086952368        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.0414E+00  8.8844E-02  5.1891E-01  1.4263E+00  4.2919E-01  9.5374E-01  1.0000E-02  1.0000E-02  8.1395E-01  5.7820E-01
             2.5632E+00
 PARAMETER:  1.4053E-01 -2.3209E+00 -5.5603E-01  4.5510E-01 -7.4585E-01  5.2632E-02 -7.5250E+00 -1.3370E+01 -1.0586E-01 -4.4783E-01
             1.0413E+00
 GRADIENT:   2.7899E-01  8.1335E-01  2.8976E+00  1.8108E+00 -5.8367E+00  3.7671E-01  0.0000E+00  0.0000E+00 -1.6101E-01  3.2812E-02
             3.3880E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1451.50480017362        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      733
 NPARAMETR:  1.0363E+00  1.0428E-02  5.0609E-01  1.4571E+00  4.1097E-01  9.4994E-01  1.0000E-02  1.0000E-02  7.9578E-01  5.8813E-01
             2.5448E+00
 PARAMETER:  1.3565E-01 -4.4632E+00 -5.8103E-01  4.7643E-01 -7.8923E-01  4.8647E-02 -1.7898E+01 -2.7589E+01 -1.2843E-01 -4.3080E-01
             1.0340E+00
 GRADIENT:  -2.6649E+00  1.0530E-01  3.6061E+00  1.1916E+01 -8.0751E+00 -4.2662E-01  0.0000E+00  0.0000E+00 -9.2380E-01 -2.1616E-02
            -2.8831E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1451.56482402901        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  1.0369E+00  1.0000E-02  5.1126E-01  1.4554E+00  4.1482E-01  9.5031E-01  1.0000E-02  1.0000E-02  7.9628E-01  5.8890E-01
             2.5472E+00
 PARAMETER:  1.3621E-01 -4.5127E+00 -5.7088E-01  4.7527E-01 -7.7991E-01  4.9028E-02 -1.8178E+01 -2.7940E+01 -1.2781E-01 -4.2949E-01
             1.0350E+00
 GRADIENT:  -8.7481E+00  0.0000E+00  4.8449E-01 -1.0665E+01 -4.8259E+00 -5.7409E-01  0.0000E+00  0.0000E+00 -5.9424E-01  1.4190E-01
            -7.4893E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1452.08531788761        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  1.0402E+00  1.0000E-02  6.2274E-01  1.5102E+00  4.8092E-01  9.4743E-01  1.0000E-02  1.0000E-02  7.8068E-01  5.7856E-01
             2.6036E+00
 PARAMETER:  1.3938E-01 -4.6051E+00 -3.7362E-01  5.1226E-01 -6.3206E-01  4.6002E-02 -1.8350E+01 -2.8355E+01 -1.4759E-01 -4.4722E-01
             1.0569E+00
 GRADIENT:  -3.5933E-03  0.0000E+00 -8.3467E-03 -1.0482E-02  1.2816E-02 -5.7327E-04  0.0000E+00  0.0000E+00 -1.9818E-03 -9.7323E-04
             5.3195E-03

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1452.08531788761        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  1.0402E+00  1.0000E-02  6.2274E-01  1.5102E+00  4.8092E-01  9.4743E-01  1.0000E-02  1.0000E-02  7.8068E-01  5.7856E-01
             2.6036E+00
 PARAMETER:  1.3938E-01 -4.6051E+00 -3.7362E-01  5.1226E-01 -6.3206E-01  4.6002E-02 -1.8350E+01 -2.8355E+01 -1.4759E-01 -4.4722E-01
             1.0569E+00
 GRADIENT:  -3.5933E-03  0.0000E+00 -8.3467E-03 -1.0482E-02  1.2816E-02 -5.7327E-04  0.0000E+00  0.0000E+00 -1.9818E-03 -9.7323E-04
             5.3195E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1027
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.5282E-04 -2.5925E-06  7.5018E-05 -1.0114E-02 -9.6330E-03
 SE:             2.9101E-02  1.6493E-06  2.1427E-04  2.6820E-02  1.7453E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8484E-01  1.1598E-01  7.2625E-01  7.0610E-01  5.8100E-01

 ETASHRINKSD(%)  2.5088E+00  9.9994E+01  9.9282E+01  1.0149E+01  4.1529E+01
 ETASHRINKVR(%)  4.9547E+00  1.0000E+02  9.9995E+01  1.9267E+01  6.5811E+01
 EBVSHRINKSD(%)  2.5688E+00  9.9995E+01  9.9255E+01  9.8919E+00  4.1683E+01
 EBVSHRINKVR(%)  5.0716E+00  1.0000E+02  9.9994E+01  1.8805E+01  6.5992E+01
 RELATIVEINF(%)  8.1048E+01  1.2184E-08  1.7793E-04  7.9360E+00  6.9008E-01
 EPSSHRINKSD(%)  2.9481E+01
 EPSSHRINKVR(%)  5.0270E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1452.0853178876130     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -716.93449132387479     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1452.085       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  6.23E-01  1.51E+00  4.81E-01  9.47E-01  1.00E-02  1.00E-02  7.81E-01  5.79E-01  2.60E+00
 


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
+        1.10E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.03E+01  0.00E+00  1.55E+03
 
 TH 4
+       -5.45E+01  0.00E+00 -1.30E+02  6.51E+02
 
 TH 5
+        7.90E+01  0.00E+00 -2.58E+03 -1.92E+02  4.69E+03
 
 TH 6
+        4.33E+00  0.00E+00  9.07E+00 -1.14E+01 -7.75E+00  2.03E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.05E+01  0.00E+00  2.58E+01 -1.58E+01  9.51E+00  1.19E+01  0.00E+00  0.00E+00  2.22E+02
 
 TH10
+       -7.14E+00  0.00E+00 -3.41E+01 -5.66E+00  4.29E+01  2.46E+00  0.00E+00  0.00E+00  1.91E-01  6.54E+01
 
 TH11
+       -1.30E+01  0.00E+00 -4.14E+00 -7.92E+00 -1.63E+01  2.71E+00  0.00E+00  0.00E+00  1.26E+01  2.92E+01  4.66E+01
 
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
 #CPUT: Total CPU Time in Seconds,       14.658
Stop Time:
Sat Sep 18 09:51:36 CDT 2021
