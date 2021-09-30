Wed Sep 29 13:50:57 CDT 2021
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
$DATA ../../../../data/spa/A3/dat83.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   728.503554108023        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4889E+02  1.4169E+02  1.7700E+02 -2.4500E+01  5.3012E+01  1.9601E+01 -5.5025E+01 -6.9388E+01 -1.3906E+02 -1.7272E+02
            -4.2582E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1178.47477409175        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0359E+00  9.8293E-01  8.6791E-01  1.2253E+00  1.1891E+00  8.2908E-01  8.6151E-01  9.9263E-01  8.0447E-01  9.4686E-01
             5.4464E+00
 PARAMETER:  1.3523E-01  8.2785E-02 -4.1666E-02  3.0316E-01  2.7317E-01 -8.7433E-02 -4.9069E-02  9.2607E-02 -1.1757E-01  4.5397E-02
             1.7949E+00
 GRADIENT:   5.0080E-02 -5.0119E+00 -4.2218E+01  4.1915E+01  1.7515E+01 -3.2856E+01  1.0374E+01  7.7035E+00  1.9076E+01  1.1073E+01
             1.5653E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1201.96451700993        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0157E+00  1.3465E+00  1.2981E+00  1.0259E+00  2.0840E+00  9.1809E-01  1.8751E-01  8.0024E-01  8.2137E-01  1.6378E+00
             5.0812E+00
 PARAMETER:  1.1557E-01  3.9751E-01  3.6093E-01  1.2555E-01  8.3428E-01  1.4542E-02 -1.5739E+00 -1.2284E-01 -9.6782E-02  5.9336E-01
             1.7255E+00
 GRADIENT:  -2.2325E+01  4.5945E+01 -4.1926E+00  6.2014E+01 -6.2930E+00 -3.9863E+00  6.3649E-01  1.0005E+00  1.0354E+01  5.2409E+00
             9.7186E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1211.26516370764        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.9865E-01  1.1599E+00  2.1206E+00  1.0686E+00  2.7090E+00  9.2019E-01  2.3667E-01  1.6418E+00  4.9206E-01  1.5768E+00
             4.5972E+00
 PARAMETER:  9.8644E-02  2.4836E-01  8.5170E-01  1.6631E-01  1.0966E+00  1.6821E-02 -1.3411E+00  5.9582E-01 -6.0916E-01  5.5541E-01
             1.6254E+00
 GRADIENT:  -4.1382E-02 -1.9124E+00 -2.1687E+00 -1.8922E+00  1.5170E+00 -1.2365E+00 -6.0142E-03  1.9659E-01 -1.7880E-01  5.3663E-01
             2.0932E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1214.29136560250        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0038E+00  1.1880E+00  2.1370E+00  1.0395E+00  1.8690E+00  9.3063E-01  2.9826E-01  3.0517E+00  4.0608E-01  5.5611E-01
             4.5335E+00
 PARAMETER:  1.0375E-01  2.7227E-01  8.5940E-01  1.3877E-01  7.2541E-01  2.8106E-02 -1.1098E+00  1.2157E+00 -8.0120E-01 -4.8679E-01
             1.6115E+00
 GRADIENT:   3.7536E+00  7.6156E+00 -2.8345E+00  8.3603E+00 -5.9573E-01  6.8724E-01 -4.1951E-01 -9.0628E-01 -1.3707E+00  1.1943E+00
            -7.7544E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1216.06725042854        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0030E+00  1.1369E+00  3.9460E+00  1.0764E+00  1.8651E+00  9.3266E-01  6.9176E-01  4.5986E+00  2.6096E-01  1.7700E-01
             4.5068E+00
 PARAMETER:  1.0295E-01  2.2828E-01  1.4727E+00  1.7365E-01  7.2330E-01  3.0286E-02 -2.6851E-01  1.6258E+00 -1.2434E+00 -1.6316E+00
             1.6056E+00
 GRADIENT:  -7.7897E-01  1.9715E+00  3.4406E+00 -2.1072E+00  1.1477E+00  3.2022E+00  2.3250E+00 -5.8230E+00  4.2066E-01  1.7105E-01
             6.8529E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1216.66407758428        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.0038E+00  1.1572E+00  5.1278E+00  1.0620E+00  1.9224E+00  9.2469E-01  5.7109E-01  5.6260E+00  2.2989E-01  1.1570E-01
             4.5366E+00
 PARAMETER:  1.0382E-01  2.4597E-01  1.7347E+00  1.6018E-01  7.5355E-01  2.1703E-02 -4.6022E-01  1.8274E+00 -1.3702E+00 -2.0567E+00
             1.6122E+00
 GRADIENT:   1.9894E+00 -3.9101E+00  3.1561E+00 -7.6462E+00  2.9492E+00  4.9525E-01  7.7664E-02 -4.7780E+00 -1.1396E-01  7.1354E-02
             5.8693E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1216.68920371494        NO. OF FUNC. EVALS.:  90
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  1.0032E+00  1.1641E+00  5.4953E+00  1.0567E+00  1.8989E+00  9.2416E-01  5.7676E-01  5.8749E+00  2.2416E-01  9.7982E-02
             4.5262E+00
 PARAMETER:  1.0315E-01  2.5195E-01  1.8039E+00  1.5517E-01  7.4126E-01  2.1129E-02 -4.5033E-01  1.8707E+00 -1.3954E+00 -2.2230E+00
             1.6099E+00
 GRADIENT:  -1.9256E+01 -6.0901E+00  2.9610E+00 -9.7228E+00  1.2098E+00 -1.1870E+00 -1.9888E-01 -5.6609E+00 -3.1970E-01  5.3207E-02
            -9.1106E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1216.83897036407        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:      729             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0052E+00  1.1852E+00  5.7872E+00  1.0529E+00  1.9003E+00  9.2982E-01  5.9297E-01  6.2230E+00  2.1954E-01  8.6066E-02
             4.5033E+00
 PARAMETER:  1.0514E-01  2.6989E-01  1.8556E+00  1.5154E-01  7.4203E-01  2.7233E-02 -4.2261E-01  1.9283E+00 -1.4162E+00 -2.3526E+00
             1.6048E+00
 GRADIENT:   1.3242E+00  1.0739E+01  2.8245E+00  1.1019E+01 -7.7277E-01  1.6100E+00 -1.1552E-01 -3.9801E+00 -1.9035E-01  4.2568E-02
            -1.9603E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1216.88448615424        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  1.0052E+00  1.1768E+00  5.7802E+00  1.0529E+00  1.9002E+00  9.2840E-01  6.1340E-01  6.2357E+00  2.1952E-01  8.6045E-02
             4.5078E+00
 PARAMETER:  1.0516E-01  2.6282E-01  1.8544E+00  1.5155E-01  7.4195E-01  2.5702E-02 -3.8874E-01  1.9303E+00 -1.4163E+00 -2.3529E+00
             1.6058E+00
 GRADIENT:  -1.6509E+01 -1.6536E+00  2.6703E+00 -3.8829E+00  2.6354E-01  1.3940E-01  2.3424E-01 -5.1053E+00 -2.2801E-01  4.1333E-02
            -1.1299E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1217.22968719077        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1021
 NPARAMETR:  1.0057E+00  1.1771E+00  5.6195E+00  1.0531E+00  1.8996E+00  9.1771E-01  5.2270E-01  6.5792E+00  2.1957E-01  8.5915E-02
             4.5999E+00
 PARAMETER:  1.0565E-01  2.6304E-01  1.8262E+00  1.5171E-01  7.4162E-01  1.4129E-02 -5.4875E-01  1.9839E+00 -1.4161E+00 -2.3544E+00
             1.6260E+00
 GRADIENT:  -1.7774E+01 -1.4885E+00  1.7247E+00 -1.6477E+00 -7.1329E-01 -3.4419E+00 -4.9586E-01 -3.3905E+00 -3.7955E-01  3.9515E-02
             2.4265E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1217.40051261271        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     1142
 NPARAMETR:  1.0061E+00  1.1750E+00  5.6344E+00  1.0548E+00  1.8912E+00  9.2481E-01  4.7214E-01  6.5634E+00  3.4893E-01  2.5895E-02
             4.5873E+00
 PARAMETER:  1.0607E-01  2.6128E-01  1.8289E+00  1.5337E-01  7.3722E-01  2.1836E-02 -6.5048E-01  1.9815E+00 -9.5288E-01 -3.5537E+00
             1.6233E+00
 GRADIENT:  -1.6883E+01 -3.2597E+00  1.8243E+00 -4.7345E+00 -1.9188E-01 -7.2926E-01  1.9014E-01 -3.4592E+00 -3.0127E-02  3.6830E-03
             5.4384E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1217.52307303499        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1318
 NPARAMETR:  1.0149E+00  1.1803E+00  5.6282E+00  1.0590E+00  1.9166E+00  9.2687E-01  4.3933E-01  6.5969E+00  3.7556E-01  1.5841E-02
             4.5836E+00
 PARAMETER:  1.1476E-01  2.6579E-01  1.8278E+00  1.5732E-01  7.5053E-01  2.4053E-02 -7.2250E-01  1.9866E+00 -8.7935E-01 -4.0452E+00
             1.6225E+00
 GRADIENT:   8.5269E-01  8.1372E-01  1.7768E+00  8.1244E-01 -1.2558E-01 -2.6095E-01  2.0610E-02 -3.4738E+00 -9.1556E-02  1.2961E-03
             2.5316E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1217.58889032220        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1494
 NPARAMETR:  1.0145E+00  1.1771E+00  5.6504E+00  1.0625E+00  1.9504E+00  9.2815E-01  2.8914E-01  6.8209E+00  4.4821E-01  1.0000E-02
             4.5797E+00
 PARAMETER:  1.1444E-01  2.6308E-01  1.8317E+00  1.6061E-01  7.6803E-01  2.5441E-02 -1.1409E+00  2.0200E+00 -7.0249E-01 -7.8900E+00
             1.6216E+00
 GRADIENT:   3.8340E-01 -1.4763E-01  3.9509E-01  1.4274E+00 -1.5372E-01  5.0991E-02 -1.2051E-01 -1.5095E+00 -9.5079E-02  0.0000E+00
             5.6421E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1217.60044026671        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1680
 NPARAMETR:  1.0144E+00  1.1779E+00  5.6510E+00  1.0621E+00  1.9540E+00  9.2793E-01  3.1200E-01  6.8413E+00  4.5107E-01  1.0000E-02
             4.5789E+00
 PARAMETER:  1.1435E-01  2.6370E-01  1.8318E+00  1.6028E-01  7.6987E-01  2.5206E-02 -1.0647E+00  2.0230E+00 -6.9613E-01 -8.1572E+00
             1.6215E+00
 GRADIENT:   7.1111E-02 -8.8891E-01 -6.5147E-02  1.1523E+00 -4.5577E-01  3.0758E-02 -4.5045E-02 -9.1898E-01  1.6513E-01  0.0000E+00
             7.0092E-01

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1217.60236558764        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     1802
 NPARAMETR:  1.0146E+00  1.1787E+00  5.6431E+00  1.0619E+00  1.9529E+00  9.2788E-01  3.2770E-01  6.8528E+00  4.4796E-01  1.0000E-02
             4.5731E+00
 PARAMETER:  1.1441E-01  2.6477E-01  1.8319E+00  1.6028E-01  7.6990E-01  2.5132E-02 -1.0162E+00  2.0231E+00 -6.9609E-01 -8.1572E+00
             1.6214E+00
 GRADIENT:  -2.7415E-01  2.2764E+02  8.1906E+01  3.7814E+02  1.9181E+02 -2.0445E-02 -2.6123E-03 -7.4630E+01  2.2704E-01  0.0000E+00
             8.8216E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1802
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.8006E-03 -7.4228E-03 -1.1013E-02 -7.7750E-03 -1.9009E-05
 SE:             2.8474E-02  6.8939E-03  6.6370E-03  1.2245E-02  3.9503E-05
 N:                     100         100         100         100         100

 P VAL.:         8.6612E-01  2.8160E-01  9.7056E-02  5.2546E-01  6.3036E-01

 ETASHRINKSD(%)  4.6085E+00  7.6905E+01  7.7765E+01  5.8978E+01  9.9868E+01
 ETASHRINKVR(%)  9.0046E+00  9.4666E+01  9.5056E+01  8.3172E+01  1.0000E+02
 EBVSHRINKSD(%)  4.7012E+00  7.7259E+01  8.2723E+01  5.9170E+01  9.9796E+01
 EBVSHRINKVR(%)  9.1815E+00  9.4829E+01  9.7015E+01  8.3330E+01  1.0000E+02
 RELATIVEINF(%)  8.4112E+01  6.9489E-03  8.6986E-01  2.1906E-02  2.0323E-04
 EPSSHRINKSD(%)  1.3738E+01
 EPSSHRINKVR(%)  2.5589E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1217.6023655876397     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -482.45153902390155     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1217.602       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.18E+00  5.65E+00  1.06E+00  1.95E+00  9.28E-01  3.28E-01  6.84E+00  4.51E-01  1.00E-02  4.58E+00
 


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
+        1.15E+03
 
 TH 2
+       -1.72E+02  1.57E+04
 
 TH 3
+        5.05E+00  2.21E+01  7.61E+01
 
 TH 4
+       -2.50E+02  1.66E+02  2.98E+01  5.26E+04
 
 TH 5
+        4.50E+01  1.14E+02  3.28E+02  1.60E+02  2.93E+03
 
 TH 6
+        1.01E+01 -2.20E+01  2.39E+00 -3.08E+01  1.63E+01  1.87E+02
 
 TH 7
+       -3.65E+00 -1.87E+01 -4.90E+00 -2.85E+01 -3.23E+01  1.60E+00  6.19E+00
 
 TH 8
+       -3.82E+00 -1.74E+01 -4.63E+01 -2.25E+01 -3.17E+02 -1.68E+00  3.65E+00  4.27E+01
 
 TH 9
+       -4.16E+04 -1.54E+04 -1.03E+01 -2.79E+01 -6.76E+01  6.30E+00  1.40E+01  7.73E+00  1.54E+04
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.14E+00  1.45E+01  1.13E+01  2.17E+01 -5.14E+01  7.40E+00 -1.62E+00 -8.61E+00 -3.03E+00  0.00E+00  1.63E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.325
Stop Time:
Wed Sep 29 13:51:29 CDT 2021
