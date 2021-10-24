Sun Oct 24 00:31:07 CDT 2021
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
$DATA ../../../../data/SD3/TD1/dat63.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2157.76670072063        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3619E+02 -5.4303E+01 -2.4370E+01 -4.1504E+01  4.3799E+01  2.8453E+01 -4.0316E+00  5.9371E+00 -9.0178E-01  1.0868E+01
             3.2621E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2167.52085912965        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0454E+00  1.0385E+00  1.0534E+00  1.0630E+00  1.0006E+00  1.0486E+00  1.0257E+00  9.8245E-01  1.0027E+00  9.4189E-01
             9.6401E-01
 PARAMETER:  1.4439E-01  1.3779E-01  1.5202E-01  1.6114E-01  1.0060E-01  1.4749E-01  1.2541E-01  8.2292E-02  1.0270E-01  4.0135E-02
             6.3343E-02
 GRADIENT:   1.2567E+00 -2.6303E-02 -1.9762E+00 -8.7547E-01  7.1662E+00  4.8993E+00 -3.4737E+00 -9.8406E-01  1.0365E+00  1.0393E-01
            -1.6988E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2168.22287427830        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0441E+00  8.7960E-01  1.1041E+00  1.1767E+00  9.4310E-01  1.0229E+00  1.2814E+00  9.7455E-01  9.1776E-01  8.7158E-01
             9.6020E-01
 PARAMETER:  1.4318E-01 -2.8282E-02  1.9902E-01  2.6269E-01  4.1418E-02  1.2268E-01  3.4792E-01  7.4219E-02  1.4182E-02 -3.7444E-02
             5.9384E-02
 GRADIENT:   1.8489E+00  1.5247E+01  5.7464E+00  1.8923E+01 -4.6221E+00 -4.5138E+00  1.0235E+00 -1.7317E+00  2.2738E+00 -3.6779E+00
            -5.6579E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2169.17109161938        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  1.0392E+00  6.2613E-01  1.2658E+00  1.3303E+00  9.2000E-01  1.0328E+00  1.3638E+00  1.0674E+00  8.7221E-01  9.4784E-01
             9.6424E-01
 PARAMETER:  1.3843E-01 -3.6819E-01  3.3571E-01  3.8541E-01  1.6622E-02  1.3225E-01  4.1027E-01  1.6522E-01 -3.6720E-02  4.6427E-02
             6.3585E-02
 GRADIENT:  -8.7537E-01  7.4006E+00  5.6475E+00  5.6739E+00 -1.4255E+01  1.3382E+00  2.8492E-01  2.2156E-02  1.6538E+00  2.0864E+00
             1.7566E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2169.86426386894        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.0350E+00  3.5580E-01  1.5156E+00  1.5155E+00  9.3709E-01  1.0211E+00  1.3286E+00  1.2633E+00  8.1894E-01  9.8617E-01
             9.5984E-01
 PARAMETER:  1.3439E-01 -9.3340E-01  5.1583E-01  5.1578E-01  3.5019E-02  1.2087E-01  3.8410E-01  3.3371E-01 -9.9750E-02  8.6072E-02
             5.9014E-02
 GRADIENT:  -1.5342E-01  6.7932E+00  3.4237E+00  2.4422E+01 -5.1774E+00 -1.1553E+00  2.0879E-01 -8.5898E-01 -3.8230E-01 -5.5779E-02
            -1.6130E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2170.19885231367        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  1.0323E+00  2.0714E-01  1.6040E+00  1.6055E+00  9.2630E-01  1.0210E+00  1.1588E+00  1.3518E+00  7.8646E-01  9.8785E-01
             9.6094E-01
 PARAMETER:  1.3181E-01 -1.4743E+00  5.7249E-01  5.7346E-01  2.3444E-02  1.2079E-01  2.4738E-01  4.0143E-01 -1.4021E-01  8.7772E-02
             6.0157E-02
 GRADIENT:  -2.9844E-01  3.3079E+00  2.0361E+00  1.5804E+01 -3.9728E+00 -1.3581E-01  9.1618E-02 -2.5145E-01 -4.2776E-02  2.3453E-01
            -4.5015E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2170.41784440093        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  1.0323E+00  1.7168E-01  1.6114E+00  1.6142E+00  9.2347E-01  1.0209E+00  7.0470E-01  1.3612E+00  7.8064E-01  9.8871E-01
             9.6122E-01
 PARAMETER:  1.3179E-01 -1.6621E+00  5.7712E-01  5.7887E-01  2.0388E-02  1.2071E-01 -2.4999E-01  4.0835E-01 -1.4764E-01  8.8647E-02
             6.0451E-02
 GRADIENT:   1.1773E+00  6.1584E-01 -1.3268E-03 -1.2176E+01  1.8170E+00  1.4697E-01  3.4812E-02 -2.1657E-02  5.8543E-01  7.5452E-05
             1.0689E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2170.44524408298        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  1.0322E+00  1.5675E-01  1.6117E+00  1.6205E+00  9.2005E-01  1.0206E+00  1.7061E-01  1.3632E+00  7.7652E-01  9.8729E-01
             9.6107E-01
 PARAMETER:  1.3174E-01 -1.7531E+00  5.7731E-01  5.8273E-01  1.6669E-02  1.2037E-01 -1.6684E+00  4.0982E-01 -1.5293E-01  8.7208E-02
             6.0289E-02
 GRADIENT:   1.6904E+00  2.0584E-01 -3.2093E-01 -1.8158E+01  2.3739E+00  1.3675E-01  2.4652E-03  5.1566E-02  8.5071E-02 -4.0838E-02
             3.7024E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2170.45064802350        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1458             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0327E+00  1.5346E-01  1.6101E+00  1.6221E+00  9.1753E-01  1.0205E+00  5.4010E-02  1.3629E+00  7.7582E-01  9.8494E-01
             9.6101E-01
 PARAMETER:  1.3218E-01 -1.7743E+00  5.7628E-01  5.8371E-01  1.3930E-02  1.2032E-01 -2.8186E+00  4.0961E-01 -1.5384E-01  8.4830E-02
             6.0226E-02
 GRADIENT:   6.1876E+02  1.8699E+01  1.2659E+01  1.1912E+03  8.1103E+00  5.0935E+01  8.0971E-03  2.2154E+00  2.1092E+01  7.5034E-01
             9.8721E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2170.45271922500        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1642
 NPARAMETR:  1.0327E+00  1.5234E-01  1.6074E+00  1.6231E+00  9.1556E-01  1.0205E+00  2.4371E-02  1.3610E+00  7.7529E-01  9.8527E-01
             9.6103E-01
 PARAMETER:  1.3222E-01 -1.7816E+00  5.7462E-01  5.8437E-01  1.1779E-02  1.2028E-01 -3.6144E+00  4.0821E-01 -1.5452E-01  8.5155E-02
             6.0253E-02
 GRADIENT:   2.9222E+00  3.5177E-01  6.6785E-01 -1.7382E+01 -2.8957E-01  1.3863E-01  1.0289E-04  1.4524E-01  7.0601E-02  2.0165E-01
             5.2241E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2170.45376366938        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1830
 NPARAMETR:  1.0327E+00  1.5096E-01  1.6053E+00  1.6233E+00  9.1388E-01  1.0205E+00  1.7011E-02  1.3586E+00  7.7495E-01  9.8520E-01
             9.6103E-01
 PARAMETER:  1.3221E-01 -1.7908E+00  5.7330E-01  5.8449E-01  9.9492E-03  1.2027E-01 -3.9739E+00  4.0647E-01 -1.5496E-01  8.5090E-02
             6.0247E-02
 GRADIENT:   2.9581E+00  3.2466E-01  1.1122E+00 -1.8383E+01 -1.1959E+00  1.4858E-01  5.6458E-05  1.2301E-01  9.9372E-02  3.3665E-01
             8.5831E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2170.45585873531        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2015
 NPARAMETR:  1.0327E+00  1.4997E-01  1.6027E+00  1.6239E+00  9.1363E-01  1.0204E+00  1.1065E-02  1.3570E+00  7.7459E-01  9.8340E-01
             9.6098E-01
 PARAMETER:  1.3218E-01 -1.7973E+00  5.7169E-01  5.8485E-01  9.6679E-03  1.2023E-01 -4.4039E+00  4.0531E-01 -1.5542E-01  8.3260E-02
             6.0203E-02
 GRADIENT:   2.9546E+00  2.8865E-01  4.7937E-01 -1.8221E+01 -1.4358E-02  1.4052E-01  2.7427E-05  8.9953E-02  4.4014E-02  8.0929E-02
             1.4737E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2170.45684800802        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2204             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0327E+00  1.4886E-01  1.6008E+00  1.6244E+00  9.1256E-01  1.0204E+00  1.0000E-02  1.3550E+00  7.7430E-01  9.8301E-01
             9.6097E-01
 PARAMETER:  1.3216E-01 -1.8047E+00  5.7051E-01  5.8511E-01  8.4934E-03  1.2021E-01 -4.5984E+00  4.0378E-01 -1.5580E-01  8.2863E-02
             6.0187E-02
 GRADIENT:   6.1864E+02  1.7963E+01  1.3018E+01  1.1986E+03  6.8027E+00  5.0834E+01  0.0000E+00  2.1442E+00  2.1264E+01  9.0084E-01
             1.0112E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2170.45746406373        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2391
 NPARAMETR:  1.0327E+00  1.4869E-01  1.5990E+00  1.6245E+00  9.1185E-01  1.0204E+00  1.0000E-02  1.3537E+00  7.7424E-01  9.8258E-01
             9.6096E-01
 PARAMETER:  1.3214E-01 -1.8059E+00  5.6935E-01  5.8518E-01  7.7207E-03  1.2020E-01 -4.5984E+00  4.0286E-01 -1.5588E-01  8.2429E-02
             6.0177E-02
 GRADIENT:   2.9431E+00  3.0165E-01  5.7861E-01 -1.8186E+01 -4.2510E-01  1.3971E-01  0.0000E+00  6.3291E-02  3.8002E-02  1.1040E-01
             1.1486E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2170.45804910302        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2581             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0327E+00  1.4772E-01  1.5952E+00  1.6243E+00  9.1159E-01  1.0204E+00  1.0000E-02  1.3517E+00  7.7424E-01  9.8094E-01
             9.6093E-01
 PARAMETER:  1.3214E-01 -1.8124E+00  5.6703E-01  5.8508E-01  7.4390E-03  1.2019E-01 -4.5984E+00  4.0139E-01 -1.5588E-01  8.0753E-02
             6.0150E-02
 GRADIENT:   6.1849E+02  1.7634E+01  1.2104E+01  1.1984E+03  8.3443E+00  5.0810E+01  0.0000E+00  2.1625E+00  2.1355E+01  6.5655E-01
             9.8529E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2170.45861314626        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2768
 NPARAMETR:  1.0327E+00  1.4758E-01  1.5940E+00  1.6244E+00  9.1104E-01  1.0204E+00  1.0000E-02  1.3505E+00  7.7418E-01  9.8094E-01
             9.6093E-01
 PARAMETER:  1.3213E-01 -1.8134E+00  5.6626E-01  5.8513E-01  6.8329E-03  1.2018E-01 -4.5984E+00  4.0049E-01 -1.5596E-01  8.0753E-02
             6.0142E-02
 GRADIENT:   2.9950E+00  1.6888E-01 -2.3359E-01 -1.9354E+01  1.0151E+00  1.4595E-01  0.0000E+00  7.6026E-02  1.2865E-01 -1.0012E-01
            -1.1033E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2170.45912110710        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2946
 NPARAMETR:  1.0326E+00  1.4804E-01  1.5931E+00  1.6248E+00  9.0958E-01  1.0204E+00  1.0000E-02  1.3483E+00  7.7395E-01  9.8179E-01
             9.6094E-01
 PARAMETER:  1.3212E-01 -1.8135E+00  5.6580E-01  5.8518E-01  6.2287E-03  1.2017E-01 -4.5984E+00  3.9983E-01 -1.5603E-01  8.0729E-02
             6.0144E-02
 GRADIENT:   2.4060E-02 -2.2292E-02  2.0310E-02 -3.6490E-01  5.6699E-01  3.1336E-03  0.0000E+00  5.0196E-02  3.8485E-02 -5.1037E-02
            -3.9858E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2946
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8213E-04 -4.7211E-05 -3.3748E-02 -5.3455E-03 -3.8977E-02
 SE:             2.9888E-02  2.5563E-05  1.8988E-02  2.9499E-02  2.0861E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9247E-01  6.4768E-02  7.5509E-02  8.5620E-01  6.1707E-02

 ETASHRINKSD(%)  1.0000E-10  9.9914E+01  3.6389E+01  1.1737E+00  3.0112E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  5.9537E+01  2.3337E+00  5.1157E+01
 EBVSHRINKSD(%)  2.9629E-01  9.9920E+01  3.9100E+01  1.5913E+00  2.7727E+01
 EBVSHRINKVR(%)  5.9169E-01  1.0000E+02  6.2912E+01  3.1573E+00  4.7767E+01
 RELATIVEINF(%)  9.7365E+01  4.5910E-06  1.0062E+01  8.2723E+00  1.0969E+01
 EPSSHRINKSD(%)  3.4555E+01
 EPSSHRINKVR(%)  5.7170E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2170.4591211070997     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1251.5205879024270     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2170.459       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.48E-01  1.59E+00  1.62E+00  9.10E-01  1.02E+00  1.00E-02  1.35E+00  7.74E-01  9.81E-01  9.61E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,      101.793
Stop Time:
Sun Oct 24 00:31:25 CDT 2021
