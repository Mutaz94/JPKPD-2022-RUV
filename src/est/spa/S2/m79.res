Sat Sep 25 12:29:39 CDT 2021
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
$DATA ../../../../data/spa/S2/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1670.25281044064        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.8483E+00 -4.3692E+01 -2.9888E-01 -4.3990E+01 -2.3743E+01  2.7422E+01 -1.3079E+01 -6.6359E-01  1.6334E+01  1.0643E+01
            -5.7540E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1678.55646234704        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0051E+00  1.1963E+00  1.0563E+00  9.4166E-01  1.1483E+00  8.8483E-01  1.1818E+00  1.0236E+00  8.4472E-01  8.9820E-01
             1.2356E+00
 PARAMETER:  1.0506E-01  2.7926E-01  1.5480E-01  3.9892E-02  2.3830E-01 -2.2364E-02  2.6701E-01  1.2337E-01 -6.8747E-02 -7.3674E-03
             3.1158E-01
 GRADIENT:  -8.9481E+00  3.1374E+01 -3.6737E+00  4.1414E+01  3.3934E+01 -2.0084E+01  2.3113E+00 -6.8290E+00 -3.3836E+00 -6.4023E+00
             2.7457E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1679.85340863989        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0106E+00  1.0751E+00  1.1512E+00  1.0089E+00  1.0813E+00  9.0834E-01  1.3662E+00  1.4308E+00  7.3310E-01  7.6327E-01
             1.2220E+00
 PARAMETER:  1.1054E-01  1.7238E-01  2.4079E-01  1.0885E-01  1.7821E-01  3.8640E-03  4.1200E-01  4.5824E-01 -2.1047E-01 -1.7014E-01
             3.0050E-01
 GRADIENT:   1.3130E+01  2.7052E+01  4.2974E+00  2.1355E+01 -5.8822E+00 -8.0224E+00  9.3433E+00  5.2382E-01 -4.7183E+00 -6.5445E+00
             2.3943E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1682.14107396742        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0037E+00  1.0737E+00  1.1808E+00  9.9575E-01  1.1168E+00  9.2560E-01  1.2541E+00  1.3611E+00  8.1189E-01  9.0788E-01
             1.1408E+00
 PARAMETER:  1.0367E-01  1.7114E-01  2.6621E-01  9.5740E-02  2.1046E-01  2.2688E-02  3.2641E-01  4.0827E-01 -1.0839E-01  3.3561E-03
             2.3169E-01
 GRADIENT:   1.0805E+00  3.8764E+00  1.1926E+00  2.7530E+00 -8.7303E-01 -8.3266E-01  1.0909E+00 -5.0875E-02 -3.6583E-01 -1.6204E-01
             1.0731E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1682.43274409677        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      406
 NPARAMETR:  1.0204E+00  1.1887E+00  1.0275E+00  9.2515E-01  1.1027E+00  9.3700E-01  1.1560E+00  1.2612E+00  8.6883E-01  8.8533E-01
             1.1386E+00
 PARAMETER:  1.2019E-01  2.7283E-01  1.2717E-01  2.2203E-02  1.9778E-01  3.4925E-02  2.4495E-01  3.3203E-01 -4.0602E-02 -2.1793E-02
             2.2976E-01
 GRADIENT:   4.6811E+00  4.8972E+00  8.8561E-01  6.1778E+00 -3.4172E+00  1.0038E+00 -3.1673E-02 -6.1987E-02  6.0768E-01  4.5517E-02
             4.2341E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1682.72024603315        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0205E+00  1.4433E+00  8.0113E-01  7.6044E-01  1.1298E+00  9.3381E-01  9.8315E-01  1.1758E+00  9.9830E-01  8.8350E-01
             1.1373E+00
 PARAMETER:  1.2032E-01  4.6692E-01 -1.2173E-01 -1.7386E-01  2.2206E-01  3.1515E-02  8.3007E-02  2.6197E-01  9.8301E-02 -2.3864E-02
             2.2864E-01
 GRADIENT:   2.3010E-01  1.1150E+01  2.1081E+00  8.5263E+00 -4.5085E+00 -1.1339E+00 -1.1707E+00 -5.1866E-01  3.2706E-01 -2.9294E-01
            -2.9455E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1683.24928225812        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  1.0184E+00  1.7425E+00  5.4386E-01  5.5839E-01  1.1797E+00  9.3685E-01  8.5216E-01  1.1518E+00  1.1827E+00  9.0177E-01
             1.1338E+00
 PARAMETER:  1.1827E-01  6.5530E-01 -5.0906E-01 -4.8270E-01  2.6522E-01  3.4768E-02 -5.9983E-02  2.4132E-01  2.6780E-01 -3.4012E-03
             2.2561E-01
 GRADIENT:  -7.9114E+00  1.2922E+01  2.4506E+00  2.5708E+00 -8.2385E+00 -2.2522E-01 -1.1965E+00  1.6354E-01 -1.8865E+00  5.6828E-01
            -1.9398E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1683.64508522756        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  1.0223E+00  1.9588E+00  3.7012E-01  4.1090E-01  1.2472E+00  9.3825E-01  7.8305E-01  9.8862E-01  1.4426E+00  9.3868E-01
             1.1424E+00
 PARAMETER:  1.2204E-01  7.7233E-01 -8.9392E-01 -7.8941E-01  3.2093E-01  3.6261E-02 -1.4456E-01  8.8556E-02  4.6644E-01  3.6724E-02
             2.3316E-01
 GRADIENT:   6.5844E-01  6.5833E+00 -9.6771E-01  3.5224E+00  6.7888E-01  4.0305E-01 -4.1457E-01  5.9022E-01 -3.9039E-02  2.6440E-02
             5.2528E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1683.75501432508        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1110
 NPARAMETR:  1.0229E+00  2.0498E+00  2.9973E-01  3.4170E-01  1.2839E+00  9.3690E-01  7.6078E-01  6.7450E-01  1.6189E+00  9.6416E-01
             1.1464E+00
 PARAMETER:  1.2260E-01  8.1773E-01 -1.1049E+00 -9.7384E-01  3.4991E-01  3.4824E-02 -1.7342E-01 -2.9379E-01  5.8176E-01  6.3505E-02
             2.3666E-01
 GRADIENT:   2.3729E+00 -6.2287E+00 -7.3275E-02 -2.0944E+00 -5.4629E-01 -1.7481E-02  7.1634E-01  3.1025E-01  5.3160E-01 -2.1564E-01
             5.8898E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1683.89001222195        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1286
 NPARAMETR:  1.0220E+00  2.0252E+00  2.9973E-01  3.5822E-01  1.2590E+00  9.3728E-01  7.6778E-01  2.7636E-01  1.5672E+00  9.4867E-01
             1.1435E+00
 PARAMETER:  1.2180E-01  8.0566E-01 -1.1049E+00 -9.2660E-01  3.3035E-01  3.5225E-02 -1.6425E-01 -1.1860E+00  5.4927E-01  4.7305E-02
             2.3413E-01
 GRADIENT:   1.8552E-01 -2.5698E+00 -7.9275E-02 -8.3362E-01  1.1629E-01 -3.1363E-02  4.0004E-01  3.9958E-02 -1.4192E-01  7.1436E-02
             1.7550E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1683.91217773872        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1461
 NPARAMETR:  1.0220E+00  2.0212E+00  3.0041E-01  3.6204E-01  1.2545E+00  9.3750E-01  7.6839E-01  4.8119E-02  1.5622E+00  9.4581E-01
             1.1429E+00
 PARAMETER:  1.2176E-01  8.0370E-01 -1.1026E+00 -9.1600E-01  3.2671E-01  3.5456E-02 -1.6345E-01 -2.9341E+00  5.4607E-01  4.4284E-02
             2.3360E-01
 GRADIENT:  -4.8278E-02  1.7715E-01 -1.1312E-01  2.2758E-01  1.0108E-01 -6.2210E-03  8.5447E-02  1.2383E-03  7.4959E-02  6.0142E-02
             1.3460E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1683.91282768172        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1641             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0220E+00  2.0216E+00  3.0040E-01  3.6153E-01  1.2550E+00  9.3750E-01  7.6803E-01  1.1123E-02  1.5632E+00  9.4581E-01
             1.1430E+00
 PARAMETER:  1.2179E-01  8.0390E-01 -1.1026E+00 -9.1742E-01  3.2711E-01  3.5458E-02 -1.6392E-01 -4.3987E+00  5.4672E-01  4.4283E-02
             2.3365E-01
 GRADIENT:   4.1659E+01  8.5197E+01  3.4013E-01  7.7129E+00  1.3673E+00  2.7167E+00  1.1534E+00  8.1815E-05  1.1482E+00  2.8217E-02
             1.7900E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1683.91283533576        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1821
 NPARAMETR:  1.0220E+00  2.0216E+00  3.0041E-01  3.6156E-01  1.2550E+00  9.3752E-01  7.6802E-01  1.0320E-02  1.5632E+00  9.4590E-01
             1.1430E+00
 PARAMETER:  1.2178E-01  8.0387E-01 -1.1026E+00 -9.1733E-01  3.2712E-01  3.5487E-02 -1.6394E-01 -4.4736E+00  5.4671E-01  4.4386E-02
             2.3366E-01
 GRADIENT:  -5.6589E-02 -2.2166E-01 -5.2662E-03 -4.3790E-02  2.0694E-02  1.5686E-02 -3.9259E-03  6.0983E-05 -8.3235E-04  5.0868E-04
             3.8563E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1683.91284065172        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1987
 NPARAMETR:  1.0220E+00  2.0216E+00  3.0042E-01  3.6159E-01  1.2550E+00  9.3750E-01  7.6804E-01  1.0000E-02  1.5632E+00  9.4591E-01
             1.1430E+00
 PARAMETER:  1.2179E-01  8.0388E-01 -1.1026E+00 -9.1726E-01  3.2711E-01  3.5464E-02 -1.6392E-01 -4.8089E+00  5.4671E-01  4.4393E-02
             2.3366E-01
 GRADIENT:  -4.1705E-02 -1.3238E-01 -8.3580E-03 -1.3884E-02  1.8984E-02  7.8417E-03 -3.4980E-04  0.0000E+00  4.4246E-03  1.5373E-03
             2.3574E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1987
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2318E-04 -2.6217E-02 -2.0914E-04  3.3834E-02 -4.1184E-02
 SE:             2.9777E-02  2.6485E-02  7.3114E-05  2.0383E-02  2.0470E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9670E-01  3.2223E-01  4.2312E-03  9.6930E-02  4.4228E-02

 ETASHRINKSD(%)  2.4311E-01  1.1273E+01  9.9755E+01  3.1714E+01  3.1423E+01
 ETASHRINKVR(%)  4.8562E-01  2.1276E+01  9.9999E+01  5.3370E+01  5.2971E+01
 EBVSHRINKSD(%)  6.3875E-01  1.1711E+01  9.9788E+01  3.4659E+01  2.9524E+01
 EBVSHRINKVR(%)  1.2734E+00  2.2050E+01  1.0000E+02  5.7306E+01  5.0332E+01
 RELATIVEINF(%)  9.8656E+01  6.9924E+00  3.4399E-05  2.4931E+00  1.2467E+01
 EPSSHRINKSD(%)  4.3228E+01
 EPSSHRINKVR(%)  6.7770E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1683.9128406517239     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -948.76201408798568     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1683.913       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.02E+00  3.00E-01  3.62E-01  1.25E+00  9.38E-01  7.68E-01  1.00E-02  1.56E+00  9.46E-01  1.14E+00
 


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
+        1.20E+03
 
 TH 2
+       -9.06E+00  3.71E+02
 
 TH 3
+        8.34E+00  1.40E+02  6.39E+02
 
 TH 4
+       -2.55E+01  3.18E+02 -6.85E+02  1.56E+03
 
 TH 5
+       -7.20E+00 -1.20E+02 -3.56E+02  3.77E+02  4.02E+02
 
 TH 6
+        6.34E-01 -1.99E+00  1.03E+00 -8.27E+00 -2.53E+00  2.29E+02
 
 TH 7
+       -4.59E-01  9.05E+00 -3.08E+01 -2.26E+01  1.55E+00  1.92E+00  2.20E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.87E-01 -1.17E+01 -4.75E+01  8.22E+01 -3.83E+00  5.96E-02  9.27E+00  0.00E+00  2.76E+01
 
 TH10
+       -1.10E+00 -1.21E+01 -3.01E+01  5.22E+00 -6.16E+01  1.63E+00  1.47E+01  0.00E+00  6.21E+00  7.54E+01
 
 TH11
+       -7.97E+00 -1.61E+01 -2.72E+01  1.49E+01 -8.43E+00  3.63E+00  1.33E+01  0.00E+00  4.34E+00  1.79E+01  1.66E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.365
Stop Time:
Sat Sep 25 12:30:10 CDT 2021
