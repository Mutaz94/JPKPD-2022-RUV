Sun Oct 24 03:53:22 CDT 2021
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
$DATA ../../../../data/SD4/TD1/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1664.43103152295        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5398E+02 -6.6764E+01 -3.6343E+01 -4.7835E+01  1.1655E+01  3.3272E+01 -1.6257E+01  1.5616E+01 -7.7957E+00  1.2458E+01
             1.3697E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1676.72808641389        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.7532E-01  1.1490E+00  1.2392E+00  9.9549E-01  1.1570E+00  9.8880E-01  1.1245E+00  8.9630E-01  1.0531E+00  9.4036E-01
             9.4181E-01
 PARAMETER:  7.5015E-02  2.3887E-01  3.1443E-01  9.5482E-02  2.4579E-01  8.8733E-02  2.1735E-01 -9.4840E-03  1.5170E-01  3.8505E-02
             4.0052E-02
 GRADIENT:  -7.0825E-01 -4.2432E-01  8.9327E+00 -1.6991E+00  7.9872E+00 -2.2110E+00 -3.9071E+00 -1.1557E+00  9.2194E-01 -2.4259E+01
            -2.2870E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1678.19247472001        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.7449E-01  1.1703E+00  1.2854E+00  9.8924E-01  1.2244E+00  1.0051E+00  1.1287E+00  7.6212E-01  1.1015E+00  1.1070E+00
             9.6868E-01
 PARAMETER:  7.4163E-02  2.5725E-01  3.5108E-01  8.9185E-02  3.0242E-01  1.0504E-01  2.2105E-01 -1.7165E-01  1.9664E-01  2.0167E-01
             6.8175E-02
 GRADIENT:  -2.7115E+00 -4.0105E+00  2.6909E-01  8.2449E+00  2.4435E+01  4.2066E+00  1.1859E+00 -2.2204E+00  7.7601E+00 -7.1195E+00
            -8.4307E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1679.22069405081        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.7816E-01  1.2932E+00  1.0566E+00  9.1024E-01  1.1598E+00  9.9468E-01  1.1008E+00  5.8654E-01  1.0758E+00  1.0738E+00
             9.8175E-01
 PARAMETER:  7.7916E-02  3.5715E-01  1.5502E-01  5.9479E-03  2.4825E-01  9.4664E-02  1.9603E-01 -4.3352E-01  1.7303E-01  1.7122E-01
             8.1577E-02
 GRADIENT:   1.9975E+00  1.0877E+01  1.9773E+00  9.9921E+00 -5.0642E+00 -7.2796E-02  2.4180E-01  2.6531E-02 -3.4803E-01  1.5882E-01
             3.1129E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1679.37309713816        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.7805E-01  1.4661E+00  9.3721E-01  7.9245E-01  1.2025E+00  9.9535E-01  1.0020E+00  4.8473E-01  1.1770E+00  1.0876E+00
             9.7956E-01
 PARAMETER:  7.7804E-02  4.8263E-01  3.5149E-02 -1.3263E-01  2.8437E-01  9.5344E-02  1.0195E-01 -6.2416E-01  2.6298E-01  1.8401E-01
             7.9346E-02
 GRADIENT:  -1.5071E-01  7.7256E+00  2.0899E+00  4.9649E+00 -5.2645E+00 -1.3198E-01 -7.8476E-01  4.0217E-02 -4.9540E-01  9.7471E-02
            -9.2732E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1679.39413530510        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  9.7838E-01  1.5693E+00  8.6159E-01  7.2529E-01  1.2317E+00  9.9610E-01  9.6421E-01  3.9921E-01  1.2406E+00  1.0960E+00
             9.8136E-01
 PARAMETER:  7.8143E-02  5.5066E-01 -4.8980E-02 -2.2118E-01  3.0842E-01  9.6089E-02  6.3549E-02 -8.1826E-01  3.1556E-01  1.9164E-01
             8.1184E-02
 GRADIENT:  -4.5552E-01  8.0997E+00  6.0407E-01  6.0476E+00 -2.0495E+00 -2.6398E-02 -2.0146E-01  6.9829E-02 -5.6978E-01 -2.6279E-01
            -1.7722E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1679.43781535667        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.7914E-01  1.5851E+00  8.4311E-01  7.0360E-01  1.2402E+00  9.9633E-01  9.5550E-01  3.4108E-01  1.2626E+00  1.1003E+00
             9.8169E-01
 PARAMETER:  7.8921E-02  5.6063E-01 -7.0661E-02 -2.5154E-01  3.1527E-01  9.6321E-02  5.4478E-02 -9.7564E-01  3.3317E-01  1.9560E-01
             8.1517E-02
 GRADIENT:   1.2501E+00 -6.3458E+00  3.7353E-01 -2.9714E+00 -2.2939E-01  8.6675E-02 -2.9427E-02  4.5399E-02 -2.4803E-01 -5.7104E-02
             6.1705E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1679.44489967141        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  9.7848E-01  1.6004E+00  8.2858E-01  6.9844E-01  1.2415E+00  9.9618E-01  9.5253E-01  2.4376E-01  1.2719E+00  1.1025E+00
             9.8227E-01
 PARAMETER:  7.8249E-02  5.7027E-01 -8.8040E-02 -2.5890E-01  3.1635E-01  9.6169E-02  5.1369E-02 -1.3116E+00  3.4048E-01  1.9759E-01
             8.2109E-02
 GRADIENT:  -4.8551E-01 -3.6156E-01  2.3674E-02  1.8131E+00  3.2355E-01 -3.9591E-02  1.1204E-01  8.3151E-03  9.4613E-02  1.6416E-01
             1.1331E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1679.45079149485        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1423
 NPARAMETR:  9.7939E-01  1.5984E+00  8.2576E-01  6.9712E-01  1.2407E+00  9.9645E-01  9.5217E-01  2.1593E-01  1.2728E+00  1.1014E+00
             9.8227E-01
 PARAMETER:  7.9180E-02  5.6901E-01 -9.1448E-02 -2.6080E-01  3.1570E-01  9.6447E-02  5.0990E-02 -1.4328E+00  3.4119E-01  1.9660E-01
             8.2112E-02
 GRADIENT:   1.5909E+00 -3.8983E+00 -2.1278E-01 -2.4320E-02  9.6972E-01  7.8036E-02  1.1602E-02  6.7372E-03  2.0110E-01  1.7099E-01
             1.4580E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1679.45183056879        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1605             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7932E-01  1.5988E+00  8.2420E-01  6.9694E-01  1.2398E+00  9.9644E-01  9.5227E-01  2.0064E-01  1.2722E+00  1.1004E+00
             9.8225E-01
 PARAMETER:  7.9106E-02  5.6928E-01 -9.3343E-02 -2.6105E-01  3.1495E-01  9.6434E-02  5.1089E-02 -1.5063E+00  3.4077E-01  1.9570E-01
             8.2086E-02
 GRADIENT:   4.2467E+02  4.5797E+02  6.0903E-01  8.1144E+01  1.7850E+01  3.6007E+01  6.4826E+00  3.9107E-02  1.3995E+01  1.7464E+00
             9.0140E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1679.45259491202        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1788             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7941E-01  1.5987E+00  8.2299E-01  6.9685E-01  1.2390E+00  9.9647E-01  9.5242E-01  1.8407E-01  1.2720E+00  1.0997E+00
             9.8221E-01
 PARAMETER:  7.9197E-02  5.6919E-01 -9.4806E-02 -2.6119E-01  3.1431E-01  9.6463E-02  5.1252E-02 -1.5924E+00  3.4057E-01  1.9503E-01
             8.2048E-02
 GRADIENT:   4.2488E+02  4.5778E+02  6.7659E-01  8.1047E+01  1.7718E+01  3.6007E+01  6.4758E+00  3.3597E-02  1.3972E+01  1.7169E+00
             8.5985E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1679.45302760692        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1967
 NPARAMETR:  9.7935E-01  1.5993E+00  8.2207E-01  6.9683E-01  1.2379E+00  9.9647E-01  9.5269E-01  1.5808E-01  1.2712E+00  1.0986E+00
             9.8215E-01
 PARAMETER:  7.9135E-02  5.6959E-01 -9.5929E-02 -2.6122E-01  3.1339E-01  9.6465E-02  5.1530E-02 -1.7447E+00  3.3996E-01  1.9407E-01
             8.1994E-02
 GRADIENT:   1.4581E+00 -2.5530E+00  4.2727E-01 -5.7678E-02 -1.7256E-01  7.5112E-02  2.6160E-03 -6.3399E-04  5.3558E-03  1.3606E-02
            -7.9283E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1679.45366055240        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2151
 NPARAMETR:  9.7943E-01  1.5987E+00  8.2083E-01  6.9653E-01  1.2374E+00  9.9649E-01  9.5273E-01  1.4501E-01  1.2710E+00  1.0981E+00
             9.8213E-01
 PARAMETER:  7.9216E-02  5.6922E-01 -9.7437E-02 -2.6165E-01  3.1300E-01  9.6485E-02  5.1574E-02 -1.8310E+00  3.3978E-01  1.9361E-01
             8.1964E-02
 GRADIENT:   1.6328E+00 -3.4794E+00  3.1253E-01 -5.0150E-01  3.1774E-02  8.3564E-02 -2.4755E-02 -1.7790E-04  1.4309E-02  2.5003E-02
            -6.5737E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1679.45410615991        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2333
 NPARAMETR:  9.7957E-01  1.5986E+00  8.1959E-01  6.9656E-01  1.2364E+00  9.9655E-01  9.5302E-01  1.2596E-01  1.2706E+00  1.0973E+00
             9.8211E-01
 PARAMETER:  7.9356E-02  5.6910E-01 -9.8955E-02 -2.6161E-01  3.1222E-01  9.6542E-02  5.1878E-02 -1.9718E+00  3.3947E-01  1.9285E-01
             8.1949E-02
 GRADIENT:   1.9320E+00 -3.4625E+00  3.5804E-01 -5.7388E-01 -1.4453E-01  1.0376E-01 -1.3131E-02 -2.7301E-04  1.8129E-02  2.0268E-02
            -8.0380E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1679.45449889403        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2518
 NPARAMETR:  9.7944E-01  1.5986E+00  8.1855E-01  6.9661E-01  1.2360E+00  9.9650E-01  9.5313E-01  1.1926E-01  1.2703E+00  1.0969E+00
             9.8215E-01
 PARAMETER:  7.9225E-02  5.6910E-01 -1.0022E-01 -2.6153E-01  3.1192E-01  9.6496E-02  5.1992E-02 -2.0264E+00  3.3924E-01  1.9249E-01
             8.1989E-02
 GRADIENT:   1.6250E+00 -3.4872E+00  1.9533E-01 -3.9372E-01  5.3399E-02  8.2597E-02 -1.5523E-02  2.1982E-04  4.1233E-02  3.2581E-02
            -2.8613E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1679.45473931317        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2698
 NPARAMETR:  9.7944E-01  1.5985E+00  8.1777E-01  6.9664E-01  1.2355E+00  9.9651E-01  9.5331E-01  1.0530E-01  1.2699E+00  1.0964E+00
             9.8213E-01
 PARAMETER:  7.9230E-02  5.6904E-01 -1.0117E-01 -2.6149E-01  3.1145E-01  9.6502E-02  5.2181E-02 -2.1509E+00  3.3894E-01  1.9202E-01
             8.1970E-02
 GRADIENT:   1.6285E+00 -3.4462E+00  2.2001E-01 -4.2048E-01 -4.1430E-02  8.3272E-02 -1.3189E-02  1.2045E-04  2.8176E-02  2.6877E-02
            -4.1573E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1679.45491979247        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2881
 NPARAMETR:  9.7945E-01  1.5984E+00  8.1715E-01  6.9674E-01  1.2351E+00  9.9651E-01  9.5344E-01  8.9666E-02  1.2697E+00  1.0960E+00
             9.8214E-01
 PARAMETER:  7.9231E-02  5.6901E-01 -1.0193E-01 -2.6134E-01  3.1118E-01  9.6502E-02  5.2318E-02 -2.3117E+00  3.3882E-01  1.9170E-01
             8.1977E-02
 GRADIENT:   1.6221E+00 -3.3526E+00  1.8648E-01 -3.0316E-01  2.0380E-02  8.0749E-02 -1.3001E-02  5.1199E-05  4.1548E-02  1.4021E-02
            -4.5429E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1679.45506079052        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3067
 NPARAMETR:  9.7949E-01  1.5981E+00  8.1650E-01  6.9669E-01  1.2346E+00  9.9652E-01  9.5357E-01  7.8175E-02  1.2693E+00  1.0955E+00
             9.8211E-01
 PARAMETER:  7.9275E-02  5.6880E-01 -1.0273E-01 -2.6141E-01  3.1072E-01  9.6519E-02  5.2455E-02 -2.4488E+00  3.3846E-01  1.9120E-01
             8.1944E-02
 GRADIENT:   1.7147E+00 -3.6056E+00  2.0958E-01 -5.2797E-01 -7.2483E-02  8.7499E-02 -2.1648E-02  7.9933E-05  1.8981E-02  7.6429E-03
            -5.3988E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1679.45520627248        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3251
 NPARAMETR:  9.7940E-01  1.5982E+00  8.1586E-01  6.9686E-01  1.2345E+00  9.9649E-01  9.5371E-01  7.7668E-02  1.2694E+00  1.0955E+00
             9.8220E-01
 PARAMETER:  7.9180E-02  5.6889E-01 -1.0351E-01 -2.6118E-01  3.1063E-01  9.6487E-02  5.2605E-02 -2.4553E+00  3.3852E-01  1.9117E-01
             8.2037E-02
 GRADIENT:   1.4877E+00 -3.4518E+00  8.7461E-03 -1.7059E-01  1.6274E-01  7.1521E-02  3.1555E-03  3.2680E-04  8.6275E-02  4.6403E-02
             3.0212E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1679.45529459340        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3435
 NPARAMETR:  9.7944E-01  1.5982E+00  8.1553E-01  6.9692E-01  1.2342E+00  9.9651E-01  9.5380E-01  7.0386E-02  1.2692E+00  1.0952E+00
             9.8218E-01
 PARAMETER:  7.9221E-02  5.6885E-01 -1.0392E-01 -2.6108E-01  3.1040E-01  9.6504E-02  5.2700E-02 -2.5538E+00  3.3835E-01  1.9091E-01
             8.2016E-02
 GRADIENT:   1.5751E+00 -3.3903E+00  1.3771E-02 -1.4821E-01  1.2464E-01  7.7696E-02  3.2808E-03  2.8015E-04  8.1937E-02  3.8998E-02
             2.0252E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1679.45535606764        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3623
 NPARAMETR:  9.7941E-01  1.5981E+00  8.1474E-01  6.9704E-01  1.2339E+00  9.9651E-01  9.5399E-01  1.8254E-02  1.2691E+00  1.0950E+00
             9.8227E-01
 PARAMETER:  7.9199E-02  5.6879E-01 -1.0489E-01 -2.6091E-01  3.1017E-01  9.6501E-02  5.2899E-02 -3.9034E+00  3.3830E-01  1.9077E-01
             8.2110E-02
 GRADIENT:   1.5150E+00 -3.3830E+00 -7.8763E-02  1.0458E-02  3.0106E-01  7.3400E-02  1.0435E-02  3.1922E-05  1.1887E-01  5.4385E-02
             5.6695E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1679.45550886945        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3815             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7946E-01  1.5979E+00  8.1456E-01  6.9701E-01  1.2333E+00  9.9652E-01  9.5409E-01  1.0000E-02  1.2685E+00  1.0944E+00
             9.8215E-01
 PARAMETER:  7.9244E-02  5.6866E-01 -1.0510E-01 -2.6096E-01  3.0970E-01  9.6517E-02  5.3000E-02 -4.6853E+00  3.3784E-01  1.9025E-01
             8.1988E-02
 GRADIENT:   4.2481E+02  4.5691E+02  8.3016E-01  8.0971E+01  1.6962E+01  3.6002E+01  6.4980E+00  0.0000E+00  1.3734E+01  1.5913E+00
             7.9100E-01

0ITERATION NO.:  107    OBJECTIVE VALUE:  -1679.45550886945        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:     3873
 NPARAMETR:  9.7946E-01  1.5979E+00  8.1456E-01  6.9701E-01  1.2333E+00  9.9652E-01  9.5409E-01  1.0000E-02  1.2685E+00  1.0944E+00
             9.8215E-01
 PARAMETER:  7.9244E-02  5.6866E-01 -1.0510E-01 -2.6096E-01  3.0970E-01  9.6517E-02  5.3000E-02 -4.6853E+00  3.3784E-01  1.9025E-01
             8.1988E-02
 GRADIENT:   1.6191E+00 -3.3651E+00  8.3415E-02 -2.4799E-01 -4.4039E-02  8.0987E-02  2.1050E-03  0.0000E+00  5.3742E-02  2.6318E-02
            -9.0136E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3873
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.4785E-04 -2.2290E-02 -3.0751E-04  1.6849E-02 -3.5517E-02
 SE:             2.9793E-02  2.3074E-02  1.1401E-04  2.2293E-02  2.2651E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9604E-01  3.3403E-01  6.9895E-03  4.4977E-01  1.1687E-01

 ETASHRINKSD(%)  1.8993E-01  2.2699E+01  9.9618E+01  2.5317E+01  2.4117E+01
 ETASHRINKVR(%)  3.7950E-01  4.0246E+01  9.9999E+01  4.4224E+01  4.2418E+01
 EBVSHRINKSD(%)  4.4582E-01  2.1340E+01  9.9677E+01  2.7658E+01  2.1828E+01
 EBVSHRINKVR(%)  8.8966E-01  3.8127E+01  9.9999E+01  4.7666E+01  3.8891E+01
 RELATIVEINF(%)  9.8911E+01  3.2516E+00  1.2735E-04  2.7637E+00  1.2930E+01
 EPSSHRINKSD(%)  4.3158E+01
 EPSSHRINKVR(%)  6.7690E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1679.4555088694492     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -944.30468230571103     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1679.456       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.60E+00  8.15E-01  6.97E-01  1.23E+00  9.97E-01  9.54E-01  1.00E-02  1.27E+00  1.09E+00  9.82E-01
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,      124.807
Stop Time:
Sun Oct 24 03:53:44 CDT 2021
