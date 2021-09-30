Wed Sep 29 19:47:25 CDT 2021
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
$DATA ../../../../data/spa/D/dat16.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1037.56252255264        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3684E+02 -1.7181E+02 -1.0543E+02 -1.5673E+02  3.0296E+02 -3.1179E+02 -3.0047E+02 -1.8865E+01 -4.3350E+02 -1.7190E+02
            -1.4173E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1429.37917179060        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      131
 NPARAMETR:  1.2665E+00  1.0922E+00  1.0484E+00  1.1971E+00  1.0074E+00  1.5219E+00  2.0260E+00  1.0307E+00  1.8380E+00  1.3406E+00
             9.3555E-01
 PARAMETER:  3.3629E-01  1.8821E-01  1.4726E-01  2.7987E-01  1.0736E-01  5.1995E-01  8.0606E-01  1.3027E-01  7.0867E-01  3.9314E-01
             3.3382E-02
 GRADIENT:   1.2084E+02 -5.9520E+01 -6.7086E+01  8.4551E+01  5.4926E+01 -5.7869E+01 -1.2673E+02 -7.2134E-01  1.4265E+01  2.6766E+01
             8.2157E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1466.93757638819        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.2107E+00  6.9696E-01  3.2814E+00  1.6477E+00  1.3007E+00  1.6045E+00  3.4839E+00  1.3578E+00  1.7001E+00  1.6202E+00
             1.0038E+00
 PARAMETER:  2.9119E-01 -2.6103E-01  1.2883E+00  5.9940E-01  3.6287E-01  5.7284E-01  1.3482E+00  4.0585E-01  6.3068E-01  5.8254E-01
             1.0379E-01
 GRADIENT:   6.9314E+01  1.3392E+01 -1.1135E+01  8.5041E+01  3.4070E+01 -1.8120E+01 -5.9288E+01 -2.1310E-01 -4.4351E+00  1.6681E+01
            -8.5907E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1474.95925680821        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      486
 NPARAMETR:  1.1646E+00  5.6043E-01  3.4497E+00  1.6580E+00  1.1959E+00  1.6320E+00  4.4087E+00  7.7009E-01  1.6985E+00  1.4059E+00
             1.0007E+00
 PARAMETER:  2.5237E-01 -4.7905E-01  1.3383E+00  6.0562E-01  2.7892E-01  5.8978E-01  1.5836E+00 -1.6125E-01  6.2975E-01  4.4068E-01
             1.0066E-01
 GRADIENT:   3.2560E+01  1.3643E+01  1.0515E+01  4.3024E+01 -2.6091E+00 -4.2271E+00 -5.6835E+01 -4.1659E+00  7.6371E+00 -5.2321E+00
            -2.5829E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1479.78779665043        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:      659
 NPARAMETR:  1.1309E+00  4.8488E-01  3.1665E+00  1.5631E+00  1.1905E+00  1.6462E+00  4.5584E+00  7.7101E-01  1.6589E+00  1.4226E+00
             1.0300E+00
 PARAMETER:  2.2298E-01 -6.2385E-01  1.2526E+00  5.4669E-01  2.7439E-01  5.9848E-01  1.6170E+00 -1.6005E-01  6.0618E-01  4.5250E-01
             1.2960E-01
 GRADIENT:   5.6630E+00  2.7327E+00  4.2559E+00 -1.2125E+00  1.3727E+01  7.4988E-01 -6.4140E+01 -3.2060E+00  6.7376E+00 -4.8672E-02
            -5.1629E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1481.06083887723        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  1.1212E+00  4.6818E-01  2.8794E+00  1.5184E+00  1.1477E+00  1.7439E+00  4.5435E+00  7.7394E-01  1.6631E+00  1.4195E+00
             1.0487E+00
 PARAMETER:  2.1438E-01 -6.5889E-01  1.1576E+00  5.1763E-01  2.3774E-01  6.5611E-01  1.6137E+00 -1.5626E-01  6.0868E-01  4.5027E-01
             1.4755E-01
 GRADIENT:  -1.1329E+00 -2.4317E+00  3.3118E+00 -1.5253E+01  2.5802E+00  2.8680E+01 -6.8191E+01 -2.2485E+00  1.2032E+01  4.7206E+00
             7.6176E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1485.07645913833        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  1.1274E+00  4.6517E-01  2.9252E+00  1.5394E+00  1.1425E+00  1.6750E+00  4.6585E+00  1.8749E+00  1.6623E+00  1.3319E+00
             9.8786E-01
 PARAMETER:  2.1991E-01 -6.6535E-01  1.1734E+00  5.3142E-01  2.3324E-01  6.1583E-01  1.6387E+00  7.2853E-01  6.0823E-01  3.8664E-01
             8.7784E-02
 GRADIENT:   4.1282E+00 -2.6576E+00  1.1955E+00 -5.1824E+00  4.7259E+00  1.3464E+01 -8.8161E+01 -6.1978E+00  2.0762E+01  6.7389E+00
             1.5876E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1486.64947585957        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1192
 NPARAMETR:  1.1449E+00  4.9777E-01  2.9085E+00  1.5227E+00  1.1308E+00  1.8019E+00  4.6700E+00  2.1511E+00  1.4127E+00  1.2218E+00
             9.1086E-01
 PARAMETER:  2.3529E-01 -5.9762E-01  1.1676E+00  5.2051E-01  2.2295E-01  6.8884E-01  1.6412E+00  8.6597E-01  4.4549E-01  3.0033E-01
             6.6380E-03
 GRADIENT:   1.8086E+01  1.2108E+00  4.5297E+00  4.3138E+00  4.6142E+00  5.5273E+01 -1.0781E+02 -5.9574E+00 -1.4509E+01 -1.4863E+00
            -7.2417E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1487.69415620966        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1384
 NPARAMETR:  1.1330E+00  4.7870E-01  2.8192E+00  1.4918E+00  1.1257E+00  1.7400E+00  4.7509E+00  2.3338E+00  1.5570E+00  1.2263E+00
             9.2746E-01
 PARAMETER:  2.2484E-01 -6.3667E-01  1.1365E+00  5.0001E-01  2.1842E-01  6.5388E-01  1.6583E+00  9.4749E-01  5.4273E-01  3.0398E-01
             2.4691E-02
 GRADIENT:   9.8283E+00 -2.4740E+00 -6.6037E+00 -2.0037E+01  7.8429E+00  3.6991E+01 -1.0153E+02  4.2753E+00  5.0620E+00  2.6849E+00
             3.9494E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1488.15151698183        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1581             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1324E+00  5.0799E-01  2.9103E+00  1.5027E+00  1.1189E+00  1.7438E+00  4.5533E+00  2.3364E+00  1.5543E+00  1.1977E+00
             9.2027E-01
 PARAMETER:  2.2437E-01 -5.7729E-01  1.1683E+00  5.0727E-01  2.1236E-01  6.5608E-01  1.6158E+00  9.4860E-01  5.4104E-01  2.8041E-01
             1.6912E-02
 GRADIENT:   1.1755E+03  1.3718E+02  9.9227E+00  8.7251E+02  1.2551E+01  1.1089E+03  1.2901E+03  6.7687E+00  2.2445E+02  3.4984E+00
             1.0332E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1488.17844979869        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1768
 NPARAMETR:  1.1326E+00  5.2042E-01  2.9234E+00  1.4988E+00  1.1228E+00  1.7444E+00  4.4962E+00  2.3531E+00  1.5546E+00  1.1961E+00
             9.2014E-01
 PARAMETER:  2.2453E-01 -5.5311E-01  1.1727E+00  5.0467E-01  2.1581E-01  6.5640E-01  1.6032E+00  9.5572E-01  5.4121E-01  2.7904E-01
             1.6772E-02
 GRADIENT:   9.4837E+00  1.3931E-02 -3.6814E-01 -7.3829E+00 -1.6405E+00  3.8888E+01 -9.0755E+01  9.3145E-01  2.1305E+00  2.8546E-02
             5.5649E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1488.20484833811        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1939
 NPARAMETR:  1.1310E+00  5.3332E-01  2.9204E+00  1.4877E+00  1.1253E+00  1.7381E+00  4.4103E+00  2.3652E+00  1.5505E+00  1.1952E+00
             9.1943E-01
 PARAMETER:  2.2313E-01 -5.2863E-01  1.1717E+00  4.9720E-01  2.1808E-01  6.5280E-01  1.5839E+00  9.6086E-01  5.3860E-01  2.7832E-01
             1.5999E-02
 GRADIENT:   8.2383E+00 -6.4969E-01 -6.8286E-01 -7.3314E+00 -1.9109E+00  3.7219E+01 -8.8901E+01  1.8978E+00  1.4198E+00 -4.1559E-01
            -5.6195E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1488.24430447930        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2130             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1344E+00  5.4596E-01  2.9297E+00  1.4710E+00  1.1297E+00  1.7447E+00  4.3801E+00  2.3584E+00  1.5584E+00  1.2004E+00
             9.2101E-01
 PARAMETER:  2.2609E-01 -5.0521E-01  1.1749E+00  4.8594E-01  2.2191E-01  6.5660E-01  1.5771E+00  9.5798E-01  5.4368E-01  2.8266E-01
             1.7719E-02
 GRADIENT:   1.1846E+03  1.2757E+02  1.1177E+01  8.0804E+02  1.4462E+01  1.1091E+03  1.2356E+03  6.2315E+00  2.1685E+02  3.6563E+00
             1.6744E+00

0ITERATION NO.:   64    OBJECTIVE VALUE:  -1488.25032851413        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     2267
 NPARAMETR:  1.1337E+00  5.5389E-01  2.9146E+00  1.4718E+00  1.1303E+00  1.7447E+00  4.3974E+00  2.3640E+00  1.5563E+00  1.1967E+00
             9.1927E-01
 PARAMETER:  2.2552E-01 -4.9079E-01  1.1697E+00  4.8651E-01  2.2246E-01  6.5660E-01  1.5810E+00  9.6037E-01  5.4232E-01  2.7960E-01
             1.5820E-02
 GRADIENT:  -1.3586E-02  3.6702E-02 -5.5079E-02  8.8047E-01 -2.0580E-01 -3.6874E-02  1.0399E+00 -2.5653E+03  1.4828E-02 -2.1074E-01
            -1.7619E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2267
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1214E-04  4.0270E-02 -5.8702E-02 -4.2820E-02 -3.7807E-02
 SE:             3.0036E-02  2.0985E-02  1.7155E-02  2.1581E-02  2.0356E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7843E-01  5.4988E-02  6.2206E-04  4.7245E-02  6.3268E-02

 ETASHRINKSD(%)  1.0000E-10  2.9697E+01  4.2528E+01  2.7699E+01  3.1806E+01
 ETASHRINKVR(%)  1.0000E-10  5.0575E+01  6.6970E+01  4.7726E+01  5.3496E+01
 EBVSHRINKSD(%)  1.2924E-01  3.3413E+01  5.1056E+01  2.1652E+01  2.4589E+01
 EBVSHRINKVR(%)  2.5831E-01  5.5662E+01  7.6045E+01  3.8617E+01  4.3131E+01
 RELATIVEINF(%)  9.9641E+01  1.3048E+01  1.1250E+01  1.8315E+01  2.7573E+01
 EPSSHRINKSD(%)  4.7218E+01
 EPSSHRINKVR(%)  7.2140E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1488.2503285141252     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -753.09950195038698     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.82
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1488.250       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.13E+00  5.54E-01  2.91E+00  1.47E+00  1.13E+00  1.74E+00  4.40E+00  2.36E+00  1.56E+00  1.20E+00  9.19E-01
 


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
+        9.35E+05
 
 TH 2
+       -2.12E+00  8.70E+01
 
 TH 3
+        1.24E-01  6.59E+00  4.83E+01
 
 TH 4
+        4.30E-01  7.08E+01 -7.16E+00  1.37E+02
 
 TH 5
+       -8.19E-01 -3.77E+01 -1.52E+01  3.03E+01  3.78E+02
 
 TH 6
+       -4.16E-02 -7.87E-01  7.17E-02  4.72E-01 -1.77E-01  6.54E+01
 
 TH 7
+        1.22E-01  9.82E+00 -5.29E-01 -3.99E+00  9.61E-01 -1.15E-01  1.27E+03
 
 TH 8
+       -1.07E+05  1.01E+05 -8.13E+03  3.79E+04 -1.10E+05 -8.05E+00  3.93E+03  1.24E+04
 
 TH 9
+        2.83E+05 -6.30E-01  2.13E+04 -1.01E+05  5.67E+00  4.35E-01 -1.04E+04 -3.21E+04  8.57E+04
 
 TH10
+        3.90E-01  2.91E+00  1.06E+01 -2.93E+00 -3.64E+01 -1.53E-01 -1.16E+00  8.10E+04 -2.16E+05  6.71E+01
 
 TH11
+       -1.50E+00  1.39E+00  1.27E+01 -2.14E+00  1.35E+01  6.32E-01 -7.19E-01 -2.26E+03  5.83E-01  3.24E+01  2.64E+02
 
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
 #CPUT: Total CPU Time in Seconds,       49.213
Stop Time:
Wed Sep 29 19:48:16 CDT 2021
