Sun Oct 24 04:30:04 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat9.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1490.65685851548        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4413E+02 -6.5375E+01 -2.3528E+01 -9.3177E+01  5.6119E+01 -5.8784E+01 -5.1704E+01 -1.2577E+00 -1.1059E+02 -1.4715E+01
             7.1161E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1528.17989217453        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  9.3660E-01  9.8348E-01  1.1205E+00  1.0694E+00  9.9739E-01  1.3596E+00  1.2909E+00  1.0197E+00  1.5456E+00  1.0786E+00
             9.4578E-01
 PARAMETER:  3.4503E-02  8.3346E-02  2.1373E-01  1.6705E-01  9.7389E-02  4.0718E-01  3.5531E-01  1.1952E-01  5.3543E-01  1.7563E-01
             4.4253E-02
 GRADIENT:  -3.0305E+01 -1.7640E+01 -4.4194E+00 -2.4223E+01 -2.2516E+01  3.6843E+01 -3.3826E-01 -3.6652E+00  3.2218E+01 -1.9729E+00
            -4.1473E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1531.91420841768        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.3626E-01  9.3388E-01  1.2999E+00  1.1344E+00  1.0894E+00  1.3231E+00  1.6319E+00  1.2761E+00  1.4474E+00  1.2007E+00
             9.4260E-01
 PARAMETER:  3.4136E-02  3.1591E-02  3.6227E-01  2.2613E-01  1.8561E-01  3.7998E-01  5.8973E-01  3.4382E-01  4.6975E-01  2.8288E-01
             4.0885E-02
 GRADIENT:  -3.1114E+01 -1.1245E+01 -1.7075E+01 -9.0333E+00  4.8330E+00  2.7574E+01  1.0441E+01  3.6821E+00  2.7685E+01  7.3106E+00
            -1.3619E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1538.23459622209        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.6428E-01  6.9141E-01  1.7435E+00  1.3300E+00  1.1247E+00  1.2150E+00  1.3778E+00  1.5835E+00  1.2384E+00  1.2192E+00
             9.5881E-01
 PARAMETER:  6.3631E-02 -2.6902E-01  6.5589E-01  3.8518E-01  2.1751E-01  2.9473E-01  4.2047E-01  5.5963E-01  3.1379E-01  2.9817E-01
             5.7940E-02
 GRADIENT:   1.0910E+01  8.2079E+00  2.5200E+00  7.6965E+00 -5.8317E+00 -2.6236E+00  1.9052E+00  4.2833E-01 -2.5139E-01 -3.7631E-01
             1.7949E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1540.05306254380        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.5257E-01  3.2567E-01  1.8092E+00  1.5601E+00  1.0387E+00  1.2039E+00  4.1936E-01  1.5506E+00  1.1202E+00  1.2170E+00
             9.5558E-01
 PARAMETER:  5.1409E-02 -1.0219E+00  6.9291E-01  5.4472E-01  1.3801E-01  2.8556E-01 -7.6902E-01  5.3867E-01  2.1347E-01  2.9642E-01
             5.4560E-02
 GRADIENT:  -1.2198E+00  3.6354E+00  1.9842E+00  7.7050E+00 -6.3187E+00 -5.1128E+00  3.7573E-02 -9.1832E-01  1.2037E+00  1.5399E+00
             6.8653E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1540.28802879266        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  9.5097E-01  2.1435E-01  1.8403E+00  1.6403E+00  1.0180E+00  1.2178E+00  1.7729E-01  1.5823E+00  1.0642E+00  1.1966E+00
             9.5380E-01
 PARAMETER:  4.9729E-02 -1.4401E+00  7.0991E-01  5.9485E-01  1.1785E-01  2.9701E-01 -1.6300E+00  5.5888E-01  1.6220E-01  2.7948E-01
             5.2698E-02
 GRADIENT:  -1.3448E+00  3.1386E+00  9.5759E-01  1.3804E+01 -3.9772E+00 -2.3363E-01  4.1182E-03 -2.5844E-01 -1.2146E+00  8.9862E-01
            -1.5696E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1540.49001737408        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  9.5004E-01  1.1063E-01  1.8672E+00  1.7058E+00  9.9898E-01  1.2233E+00  3.4776E-02  1.6109E+00  1.0185E+00  1.1733E+00
             9.5355E-01
 PARAMETER:  4.8750E-02 -2.1015E+00  7.2446E-01  6.3403E-01  9.8977E-02  3.0156E-01 -3.2588E+00  5.7681E-01  1.1834E-01  2.5983E-01
             5.2432E-02
 GRADIENT:  -4.1597E-01  1.4781E+00  1.9742E-01  8.7465E+00 -9.9527E-01  1.7749E+00  8.5050E-05  1.7002E-01 -2.1163E+00 -4.8513E-02
            -4.0697E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1540.70651871302        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.5015E-01  6.0874E-02  1.8723E+00  1.7163E+00  9.9005E-01  1.2209E+00  1.1762E-02  1.6137E+00  1.0081E+00  1.1680E+00
             9.5405E-01
 PARAMETER:  4.8870E-02 -2.6990E+00  7.2715E-01  6.4015E-01  9.0004E-02  2.9957E-01 -4.3429E+00  5.7851E-01  1.0805E-01  2.5530E-01
             5.2965E-02
 GRADIENT:   1.0583E+00  5.5991E-02 -1.4671E-01 -1.9212E+01  2.5261E+00  1.0729E+00  5.0270E-06  2.6880E-01  1.5966E+00  1.7875E-02
             1.0836E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1540.73194069488        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1421
 NPARAMETR:  9.5000E-01  4.6939E-02  1.8708E+00  1.7262E+00  9.8399E-01  1.2200E+00  1.1587E-02  1.6092E+00  9.9859E-01  1.1649E+00
             9.5393E-01
 PARAMETER:  4.8712E-02 -2.9589E+00  7.2636E-01  6.4591E-01  8.3862E-02  2.9887E-01 -4.3579E+00  5.7576E-01  9.8589E-02  2.5262E-01
             5.2840E-02
 GRADIENT:   1.1536E+00  1.1515E-01  7.1801E-01 -1.8235E+01  8.2908E-01  8.2695E-01  2.8218E-06  1.4221E-01  1.9733E-02  2.0259E-01
            -4.6268E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1540.74195803206        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  9.4976E-01  3.8152E-02  1.8676E+00  1.7332E+00  9.8137E-01  1.2195E+00  1.2304E-02  1.6071E+00  9.9708E-01  1.1628E+00
             9.5402E-01
 PARAMETER:  4.8450E-02 -3.1662E+00  7.2467E-01  6.4995E-01  8.1192E-02  2.9842E-01 -4.2979E+00  5.7442E-01  9.7074E-02  2.5079E-01
             5.2932E-02
 GRADIENT:   9.9375E-01  1.0312E-01  3.9961E-01 -1.6749E+01  1.2534E+00  6.5184E-01  2.0853E-06  1.5220E-01  6.9758E-01  1.5944E-01
            -8.1357E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1540.76292875504        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1785
 NPARAMETR:  9.4959E-01  2.3904E-02  1.8641E+00  1.7414E+00  9.7601E-01  1.2188E+00  1.2326E-02  1.6036E+00  9.9020E-01  1.1594E+00
             9.5404E-01
 PARAMETER:  4.8279E-02 -3.6337E+00  7.2279E-01  6.5467E-01  7.5714E-02  2.9791E-01 -4.2960E+00  5.7225E-01  9.0154E-02  2.4789E-01
             5.2954E-02
 GRADIENT:   1.0703E+00  6.1940E-02  6.5999E-01 -1.8059E+01  7.2654E-01  4.7079E-01  8.3380E-07  1.6555E-01  1.0009E-01  2.3796E-01
            -2.7480E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1540.77783481739        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1973
 NPARAMETR:  9.4976E-01  1.3195E-02  1.8544E+00  1.7423E+00  9.6955E-01  1.2219E+00  1.0265E-02  1.5957E+00  9.8463E-01  1.1552E+00
             9.5405E-01
 PARAMETER:  4.8454E-02 -4.2279E+00  7.1755E-01  6.5522E-01  6.9081E-02  3.0043E-01 -4.4790E+00  5.6732E-01  8.4506E-02  2.4424E-01
             5.2963E-02
 GRADIENT:   1.6358E+00  3.9184E-03  1.0810E+00 -2.5095E+01  1.2279E-01  1.5201E+00  9.6190E-07  2.5460E-01 -5.5836E-01  3.1514E-01
             2.6113E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1540.79085137642        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2157
 NPARAMETR:  9.4950E-01  1.0000E-02  1.8477E+00  1.7480E+00  9.6540E-01  1.2206E+00  1.0017E-02  1.5891E+00  9.8441E-01  1.1531E+00
             9.5405E-01
 PARAMETER:  4.8176E-02 -4.8333E+00  7.1394E-01  6.5845E-01  6.4789E-02  2.9931E-01 -4.5034E+00  5.6314E-01  8.4288E-02  2.4242E-01
             5.2959E-02
 GRADIENT:   1.3287E+00  0.0000E+00  1.5555E+00 -2.0185E+01 -1.4181E+00  1.0849E+00  7.4688E-07  1.9452E-01 -1.9071E-01  4.5536E-01
            -1.2208E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1540.79646945870        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     2356
 NPARAMETR:  9.4947E-01  1.0000E-02  1.8396E+00  1.7475E+00  9.6348E-01  1.2205E+00  1.0000E-02  1.5826E+00  9.8504E-01  1.1510E+00
             9.5413E-01
 PARAMETER:  4.8151E-02 -4.8333E+00  7.0953E-01  6.5821E-01  6.2800E-02  2.9930E-01 -4.5937E+00  5.5905E-01  8.4930E-02  2.4060E-01
             5.3046E-02
 GRADIENT:   1.3887E+00  0.0000E+00  1.3046E+00 -2.0273E+01 -1.0671E+00  1.0700E+00  0.0000E+00  1.7817E-01  3.9801E-02  3.3275E-01
             2.4245E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1540.80177675812        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2554             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4948E-01  1.0000E-02  1.8289E+00  1.7471E+00  9.6295E-01  1.2205E+00  1.0264E-02  1.5759E+00  9.8505E-01  1.1482E+00
             9.5408E-01
 PARAMETER:  4.8161E-02 -4.8333E+00  7.0369E-01  6.5795E-01  6.2241E-02  2.9928E-01 -4.4791E+00  5.5484E-01  8.4939E-02  2.3817E-01
             5.2995E-02
 GRADIENT:   3.9377E+02  0.0000E+00  1.0252E+01  1.1135E+03  7.4918E+00  1.4967E+02  1.5382E-06  2.3991E+00  1.4063E+01  2.0694E+00
             7.9995E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1540.80331614042        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2747
 NPARAMETR:  9.4947E-01  1.0000E-02  1.8225E+00  1.7468E+00  9.6187E-01  1.2205E+00  1.0000E-02  1.5713E+00  9.8507E-01  1.1461E+00
             9.5408E-01
 PARAMETER:  4.8152E-02 -4.8333E+00  7.0021E-01  6.5778E-01  6.1127E-02  2.9927E-01 -4.5830E+00  5.5190E-01  8.4958E-02  2.3641E-01
             5.2995E-02
 GRADIENT:   1.3493E+00  0.0000E+00 -4.3215E-01 -2.0412E+01  2.0207E+00  1.0685E+00  0.0000E+00  2.0365E-01 -6.0202E-03 -2.7460E-01
            -1.0458E-02

0ITERATION NO.:   78    OBJECTIVE VALUE:  -1540.80529594220        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:     2856
 NPARAMETR:  9.4946E-01  1.0000E-02  1.8156E+00  1.7464E+00  9.6036E-01  1.2205E+00  1.0000E-02  1.5649E+00  9.8511E-01  1.1451E+00
             9.5401E-01
 PARAMETER:  4.8148E-02 -4.8333E+00  7.0078E-01  6.5776E-01  5.9172E-02  2.9927E-01 -4.5372E+00  5.5100E-01  8.4997E-02  2.3788E-01
             5.3008E-02
 GRADIENT:   5.8252E-03  0.0000E+00  6.9587E-01  1.9856E-01 -1.6892E-01  3.0575E-03  0.0000E+00  1.5758E-01 -8.9592E-05  1.5393E-01
             1.5987E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2856
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.8255E-05 -3.1213E-06 -3.7281E-02 -5.9470E-03 -4.7739E-02
 SE:             2.9875E-02  1.3570E-06  1.8584E-02  2.9584E-02  2.0336E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9871E-01  2.1445E-02  4.4846E-02  8.4068E-01  1.8901E-02

 ETASHRINKSD(%)  1.0000E-10  9.9995E+01  3.7742E+01  8.8922E-01  3.1870E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.1239E+01  1.7705E+00  5.3583E+01
 EBVSHRINKSD(%)  2.6966E-01  9.9996E+01  4.0867E+01  1.2869E+00  2.7994E+01
 EBVSHRINKVR(%)  5.3860E-01  1.0000E+02  6.5032E+01  2.5573E+00  4.8152E+01
 RELATIVEINF(%)  9.7979E+01  1.5158E-08  1.0531E+01  1.0215E+01  1.0429E+01
 EPSSHRINKSD(%)  4.6052E+01
 EPSSHRINKVR(%)  7.0896E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1540.8052959421957     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -805.65446937845752     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1540.805       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  1.00E-02  1.82E+00  1.75E+00  9.60E-01  1.22E+00  1.00E-02  1.57E+00  9.85E-01  1.15E+00  9.54E-01
 


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
 #CPUT: Total CPU Time in Seconds,       88.029
Stop Time:
Sun Oct 24 04:30:20 CDT 2021
