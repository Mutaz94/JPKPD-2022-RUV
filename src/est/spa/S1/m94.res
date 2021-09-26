Sat Sep 25 10:12:06 CDT 2021
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
$DATA ../../../../data/spa/S1/dat94.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1660.87134369135        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8536E+01 -6.1118E+01 -3.2439E+01 -7.5781E+01  4.6679E+01 -6.7758E+00 -7.8264E+00  9.5660E+00 -4.0500E+01  2.6270E+01
             1.4643E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1670.76347011526        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0024E+00  9.4623E-01  1.0143E+00  1.0811E+00  9.1464E-01  9.8065E-01  9.8078E-01  9.1599E-01  1.2149E+00  7.0357E-01
             1.0247E+00
 PARAMETER:  1.0241E-01  4.4728E-02  1.1424E-01  1.7797E-01  1.0776E-02  8.0461E-02  8.0596E-02  1.2248E-02  2.9464E-01 -2.5159E-01
             1.2440E-01
 GRADIENT:   5.3685E+01  7.3848E+00  1.6824E+01  8.0529E+00 -9.2731E+00 -1.5686E+01  1.5202E+00 -1.5136E+00  1.8060E+01 -3.7782E+00
             1.8814E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1672.73460211237        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0036E+00  9.6398E-01  7.0170E-01  1.0604E+00  7.7423E-01  1.0093E+00  1.2023E+00  6.4149E-01  1.1256E+00  5.5288E-01
             9.7878E-01
 PARAMETER:  1.0359E-01  6.3312E-02 -2.5425E-01  1.5861E-01 -1.5589E-01  1.0928E-01  2.8425E-01 -3.4396E-01  2.1829E-01 -4.9262E-01
             7.8554E-02
 GRADIENT:   5.3577E+01  1.4084E+00 -2.3674E+01  3.6861E+01  2.4699E+01 -4.5045E+00  5.7912E+00  4.0552E+00  1.3732E+01  1.3321E+00
             6.6744E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1673.49335313264        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.8351E-01  9.9770E-01  7.1125E-01  1.0215E+00  7.9330E-01  1.0150E+00  1.1472E+00  5.6095E-01  1.1014E+00  5.9771E-01
             9.6668E-01
 PARAMETER:  8.3370E-02  9.7697E-02 -2.4073E-01  1.2125E-01 -1.3156E-01  1.1490E-01  2.3730E-01 -4.7812E-01  1.9661E-01 -4.1465E-01
             6.6108E-02
 GRADIENT:   1.0249E+01 -1.1723E+00 -5.9534E+00  2.9076E+00  5.7034E+00 -1.8565E+00  1.4609E+00  1.9830E+00  2.5419E+00  1.4331E+00
             1.2418E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1673.51173452738        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.7980E-01  9.8126E-01  7.0192E-01  1.0281E+00  7.8132E-01  1.0182E+00  1.1648E+00  4.9284E-01  1.0828E+00  5.9725E-01
             9.6531E-01
 PARAMETER:  7.9595E-02  8.1084E-02 -2.5394E-01  1.2770E-01 -1.4678E-01  1.1805E-01  2.5255E-01 -6.0758E-01  1.7951E-01 -4.1542E-01
             6.4697E-02
 GRADIENT:   2.3587E+00 -1.2693E+00 -3.2596E+00  5.4298E-01  3.2688E+00 -6.8665E-01  7.9781E-01  1.2188E+00  6.1878E-01  1.0823E+00
             6.4527E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1673.51388173447        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.7893E-01  9.7944E-01  6.8692E-01  1.0267E+00  7.7227E-01  1.0194E+00  1.1707E+00  4.3603E-01  1.0765E+00  5.9232E-01
             9.6468E-01
 PARAMETER:  7.8707E-02  7.9224E-02 -2.7553E-01  1.2632E-01 -1.5843E-01  1.1918E-01  2.5762E-01 -7.3004E-01  1.7367E-01 -4.2371E-01
             6.4039E-02
 GRADIENT:   2.9617E-01 -1.2093E+00 -2.1788E+00 -8.9986E-02  2.1971E+00 -3.1269E-01  4.9019E-01  8.2920E-01  8.9019E-02  8.0355E-01
             3.9028E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1673.51450286821        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  9.7876E-01  9.8383E-01  6.6940E-01  1.0217E+00  7.6430E-01  1.0200E+00  1.1716E+00  3.7203E-01  1.0742E+00  5.8677E-01
             9.6424E-01
 PARAMETER:  7.8534E-02  8.3702E-02 -3.0138E-01  1.2145E-01 -1.6879E-01  1.1982E-01  2.5838E-01 -8.8877E-01  1.7159E-01 -4.3313E-01
             6.3588E-02
 GRADIENT:  -3.5333E-01 -8.8266E-01 -1.4013E+00 -2.4071E-01  1.2966E+00 -1.2317E-01  2.8912E-01  5.3271E-01 -8.8460E-02  5.6038E-01
             2.3714E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1673.51568709818        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      509
 NPARAMETR:  9.7893E-01  9.8993E-01  6.5490E-01  1.0162E+00  7.5871E-01  1.0203E+00  1.1700E+00  3.1301E-01  1.0743E+00  5.8246E-01
             9.6398E-01
 PARAMETER:  7.8704E-02  8.9876E-02 -3.2327E-01  1.1612E-01 -1.7613E-01  1.2009E-01  2.5702E-01 -1.0615E+00  1.7166E-01 -4.4049E-01
             6.3320E-02
 GRADIENT:  -2.4511E-01 -5.3981E-01 -9.3170E-01 -1.4786E-01  7.3538E-01 -6.2622E-02  1.7796E-01  3.4280E-01 -7.0863E-02  3.8000E-01
             1.5947E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1673.51601104965        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  9.7910E-01  9.9436E-01  6.4594E-01  1.0125E+00  7.5554E-01  1.0204E+00  1.1684E+00  2.7017E-01  1.0749E+00  5.8001E-01
             9.6384E-01
 PARAMETER:  7.8878E-02  9.4340E-02 -3.3705E-01  1.1242E-01 -1.8033E-01  1.2018E-01  2.5567E-01 -1.2087E+00  1.7219E-01 -4.4471E-01
             6.3174E-02
 GRADIENT:  -3.1274E-02 -3.3643E-01 -6.8953E-01 -4.6895E-02  4.7674E-01 -5.1653E-02  1.2438E-01  2.4037E-01 -1.9417E-02  2.7037E-01
             1.2168E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1674.17264773201        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      738
 NPARAMETR:  9.9845E-01  9.3398E-01  7.1419E-01  1.0651E+00  7.6945E-01  1.0362E+00  1.2289E+00  3.3256E-01  1.0542E+00  6.2943E-01
             9.6531E-01
 PARAMETER:  9.8452E-02  3.1701E-02 -2.3661E-01  1.6308E-01 -1.6208E-01  1.3552E-01  3.0611E-01 -1.0009E+00  1.5274E-01 -3.6294E-01
             6.4693E-02
 GRADIENT:  -1.3488E+00  1.9322E+00  1.4389E+00 -7.5879E-02 -1.5932E+00  9.7475E-01  4.2227E-01 -1.2658E-01 -1.1156E+00  1.9024E-01
            -1.5072E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1674.24617717190        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  9.9934E-01  8.2994E-01  8.0276E-01  1.1393E+00  7.7347E-01  1.0337E+00  1.3023E+00  5.7269E-01  1.0306E+00  6.3129E-01
             9.6643E-01
 PARAMETER:  9.9341E-02 -8.6408E-02 -1.1971E-01  2.3039E-01 -1.5687E-01  1.3318E-01  3.6414E-01 -4.5741E-01  1.3012E-01 -3.5999E-01
             6.5856E-02
 GRADIENT:   3.0584E+00  2.3166E+00  9.1728E-01  1.5811E+00 -1.6004E+00  5.2687E-01  6.6157E-01  1.7390E-01  6.7838E-01 -8.5674E-02
             2.1837E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1674.30762709373        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  9.9630E-01  7.1159E-01  8.6549E-01  1.2169E+00  7.6198E-01  1.0308E+00  1.3918E+00  6.4042E-01  9.8989E-01  6.5054E-01
             9.6531E-01
 PARAMETER:  9.6290E-02 -2.4026E-01 -4.4462E-02  2.9627E-01 -1.7183E-01  1.3030E-01  4.3059E-01 -3.4564E-01  8.9838E-02 -3.2996E-01
             6.4691E-02
 GRADIENT:  -2.9635E-01  2.0643E+00  2.5054E+00  1.5971E+00 -4.1060E+00 -9.9010E-02 -5.2776E-01 -1.9923E-01 -8.2517E-01 -8.6945E-02
            -3.3916E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1674.38167710886        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1268
 NPARAMETR:  9.9037E-01  5.1345E-01  9.7901E-01  1.3495E+00  7.5267E-01  1.0256E+00  1.6279E+00  7.7088E-01  9.3886E-01  6.7919E-01
             9.6579E-01
 PARAMETER:  9.0321E-02 -5.6661E-01  7.8784E-02  3.9977E-01 -1.8413E-01  1.2532E-01  5.8729E-01 -1.6023E-01  3.6907E-02 -2.8686E-01
             6.5192E-02
 GRADIENT:  -5.7275E+00  5.0922E+00  5.5508E+00  7.0634E+00 -9.0306E+00 -9.3476E-01  1.6563E+00 -3.8140E-02 -1.1202E+00  5.3906E-01
             6.6206E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1674.70835617376        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1445
 NPARAMETR:  9.9168E-01  3.5989E-01  1.0310E+00  1.4447E+00  7.3536E-01  1.0271E+00  1.6597E+00  8.2861E-01  9.1357E-01  6.9241E-01
             9.6557E-01
 PARAMETER:  9.1642E-02 -9.2197E-01  1.3051E-01  4.6792E-01 -2.0739E-01  1.2670E-01  6.0663E-01 -8.8008E-02  9.6004E-03 -2.6757E-01
             6.4961E-02
 GRADIENT:   3.1372E+00  3.4325E+00  6.1061E+00  1.1185E+01 -1.0284E+01  4.5216E-01  7.5337E-02 -2.1960E-01  2.7934E-01  4.0918E-01
            -4.8598E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1674.93766605299        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1620
 NPARAMETR:  9.8798E-01  1.9802E-01  1.0575E+00  1.5332E+00  7.1182E-01  1.0251E+00  1.4922E+00  8.8719E-01  8.7487E-01  6.8784E-01
             9.6560E-01
 PARAMETER:  8.7903E-02 -1.5194E+00  1.5595E-01  5.2735E-01 -2.3993E-01  1.2475E-01  5.0028E-01 -1.9693E-02 -3.3683E-02 -2.7420E-01
             6.4999E-02
 GRADIENT:   2.3997E+00  4.3230E-01 -1.3849E+00  2.4147E+00  2.6878E-01  4.9174E-01 -1.7035E-01  5.1856E-01 -1.2072E+00 -1.8496E-01
             1.5326E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1675.02277654470        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1795
 NPARAMETR:  9.8440E-01  9.3778E-02  1.0695E+00  1.5920E+00  6.9424E-01  1.0224E+00  1.2333E+00  9.0368E-01  8.4801E-01  6.8810E-01
             9.6522E-01
 PARAMETER:  8.4282E-02 -2.2668E+00  1.6720E-01  5.6496E-01 -2.6494E-01  1.2219E-01  3.0969E-01 -1.2753E-03 -6.4868E-02 -2.7383E-01
             6.4603E-02
 GRADIENT:  -1.6504E-01  9.7888E-02 -1.3492E+00  8.3640E-01  9.2140E-01  1.3955E-02 -2.7309E-02  3.4978E-01 -6.0149E-01  1.0201E-01
             1.7859E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1675.04873973800        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1972
 NPARAMETR:  9.8337E-01  4.8606E-02  1.0804E+00  1.6185E+00  6.8826E-01  1.0216E+00  1.0221E+00  9.1313E-01  8.3562E-01  6.8393E-01
             9.6472E-01
 PARAMETER:  8.3227E-02 -2.9240E+00  1.7736E-01  5.8149E-01 -2.7359E-01  1.2134E-01  1.2184E-01  9.1283E-03 -7.9581E-02 -2.7989E-01
             6.4079E-02
 GRADIENT:  -4.2572E-02  5.8741E-02  1.2025E+00  9.1280E-01 -9.2382E-01 -8.0224E-02 -3.9534E-03 -3.0967E-01 -1.9101E-02 -2.5256E-01
            -1.9311E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1675.06203712681        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2147
 NPARAMETR:  9.8263E-01  1.9345E-02  1.0705E+00  1.6332E+00  6.7846E-01  1.0213E+00  7.6116E-01  9.0883E-01  8.2719E-01  6.8232E-01
             9.6471E-01
 PARAMETER:  8.2475E-02 -3.8453E+00  1.6813E-01  5.9051E-01 -2.8793E-01  1.2111E-01 -1.7291E-01  4.3990E-03 -8.9716E-02 -2.8225E-01
             6.4072E-02
 GRADIENT:  -8.5868E-02  8.6627E-03 -1.5402E-01 -4.4358E-02 -1.8992E-01 -1.1146E-02 -3.1855E-04  5.5203E-02 -2.8129E-02  1.3935E-01
             4.3819E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1675.06484540485        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2324
 NPARAMETR:  9.8245E-01  1.0007E-02  1.0759E+00  1.6393E+00  6.7856E-01  1.0212E+00  6.1619E-01  9.1700E-01  8.2418E-01  6.7920E-01
             9.6468E-01
 PARAMETER:  8.2298E-02 -4.5044E+00  1.7319E-01  5.9425E-01 -2.8778E-01  1.2098E-01 -3.8421E-01  1.3349E-02 -9.3364E-02 -2.8684E-01
             6.4046E-02
 GRADIENT:   5.1633E-02  1.2626E-02  1.4741E-01  3.3597E-01 -8.0338E-02 -5.0565E-04 -5.0191E-05 -3.0358E-02 -3.3664E-02 -6.0953E-02
            -3.1128E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1675.06493416988        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2506
 NPARAMETR:  9.8243E-01  1.0000E-02  1.0757E+00  1.6391E+00  6.7852E-01  1.0212E+00  6.1624E-01  9.1656E-01  8.2425E-01  6.7977E-01
             9.6471E-01
 PARAMETER:  8.2275E-02 -4.5098E+00  1.7297E-01  5.9415E-01 -2.8784E-01  1.2102E-01 -3.8413E-01  1.2877E-02 -9.3287E-02 -2.8600E-01
             6.4069E-02
 GRADIENT:  -1.3805E-02  0.0000E+00  6.1709E-02  7.0231E-02  6.1147E-03  6.0805E-03 -5.0088E-05  2.4635E-03 -2.9596E-03 -2.9679E-03
            -4.0091E-04

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1675.06501152417        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2687
 NPARAMETR:  9.8244E-01  1.0000E-02  1.0747E+00  1.6388E+00  6.7814E-01  1.0212E+00  6.1665E-01  9.1514E-01  8.2428E-01  6.8000E-01
             9.6470E-01
 PARAMETER:  8.2280E-02 -4.5192E+00  1.7204E-01  5.9398E-01 -2.8841E-01  1.2101E-01 -3.8346E-01  1.1322E-02 -9.3243E-02 -2.8567E-01
             6.4058E-02
 GRADIENT:   1.8695E-02  0.0000E+00  4.5560E-02 -2.3316E-01  1.1196E-02  7.4015E-03 -4.9555E-05 -7.1473E-03 -3.5684E-03  6.5139E-03
             7.1238E-05

0ITERATION NO.:  103    OBJECTIVE VALUE:  -1675.06501368784        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     2786
 NPARAMETR:  9.8243E-01  1.0000E-02  1.0747E+00  1.6388E+00  6.7813E-01  1.0212E+00  6.1688E-01  9.1522E-01  8.2429E-01  6.7996E-01
             9.6469E-01
 PARAMETER:  8.2278E-02 -4.5192E+00  1.7200E-01  5.9398E-01 -2.8841E-01  1.2100E-01 -3.8307E-01  1.1406E-02 -9.3235E-02 -2.8573E-01
             6.4056E-02
 GRADIENT:   3.4428E-03  0.0000E+00 -8.1892E-04 -2.2461E-01  6.5295E-02 -3.7403E-03 -5.0461E-05  1.6697E-04  2.1011E-03  6.4507E-03
             1.3192E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2786
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6874E-04 -1.6087E-04 -1.8340E-02 -3.6795E-03 -2.3759E-02
 SE:             2.9848E-02  1.0809E-04  1.8776E-02  2.9444E-02  2.0010E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9549E-01  1.3667E-01  3.2866E-01  9.0055E-01  2.3507E-01

 ETASHRINKSD(%)  5.8210E-03  9.9638E+01  3.7099E+01  1.3605E+00  3.2965E+01
 ETASHRINKVR(%)  1.1642E-02  9.9999E+01  6.0435E+01  2.7025E+00  5.5063E+01
 EBVSHRINKSD(%)  3.8636E-01  9.9666E+01  3.8363E+01  1.7784E+00  3.1950E+01
 EBVSHRINKVR(%)  7.7122E-01  9.9999E+01  6.2009E+01  3.5253E+00  5.3692E+01
 RELATIVEINF(%)  9.5344E+01  6.2160E-05  4.4798E+00  8.9848E+00  2.7989E+00
 EPSSHRINKSD(%)  4.4902E+01
 EPSSHRINKVR(%)  6.9642E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1675.0650136878367     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.91418712409848     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.27
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1675.065       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.00E-02  1.07E+00  1.64E+00  6.78E-01  1.02E+00  6.17E-01  9.15E-01  8.24E-01  6.80E-01  9.65E-01
 


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
+        1.11E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.29E+00  0.00E+00  4.47E+02
 
 TH 4
+       -7.32E+00  0.00E+00 -3.51E+01  5.88E+02
 
 TH 5
+        9.36E+00  0.00E+00 -9.39E+02 -9.75E+01  2.45E+03
 
 TH 6
+       -1.12E+01  0.00E+00 -2.43E+00 -3.01E+00  1.18E+00  2.00E+02
 
 TH 7
+        1.81E+00  0.00E+00 -2.89E-01 -4.79E-01 -1.83E+00  1.10E+00  1.68E+00
 
 TH 8
+        1.81E+01  0.00E+00 -4.66E+01 -1.79E+00 -1.16E+01 -6.37E-02 -5.81E+00  5.41E+01
 
 TH 9
+        1.95E+00  0.00E+00  1.01E+01  2.45E-01  4.99E+00 -1.46E+01  2.07E+00 -1.15E+01  2.88E+02
 
 TH10
+       -8.24E+00  0.00E+00  2.01E+00 -2.24E-01 -1.32E+02  6.23E+00  6.89E-02  3.77E+01  1.11E+00  1.01E+02
 
 TH11
+       -3.86E+00  0.00E+00 -1.56E+01 -5.53E+00 -8.00E+00  5.77E+00  9.36E-01  1.43E+01 -3.27E-01  2.40E+01  2.38E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.505
Stop Time:
Sat Sep 25 10:12:45 CDT 2021
