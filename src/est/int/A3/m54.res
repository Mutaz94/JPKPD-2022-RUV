Wed Sep 29 00:10:05 CDT 2021
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
$DATA ../../../../data/int/A3/dat54.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   85.9705474253466        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0130E+02  3.7253E+02  3.7147E+02  7.7699E+01  3.6011E+02  4.5694E+01 -3.1167E+02 -3.7708E+02 -1.3080E+02 -2.3935E+02
            -6.7505E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2365.45430210647        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0533E+00  8.5383E-01  6.9514E-01  1.0955E+00  7.9291E-01  7.9100E-01  9.2680E-01  1.5102E+00  9.6482E-01  8.0628E-01
             5.0750E+00
 PARAMETER:  1.5192E-01 -5.8018E-02 -2.6364E-01  1.9118E-01 -1.3205E-01 -1.3445E-01  2.3986E-02  5.1225E-01  6.4190E-02 -1.1533E-01
             1.7243E+00
 GRADIENT:  -4.2717E+01 -7.2385E+01 -4.4974E+01 -3.7205E+01  6.2988E+01 -4.7138E+01  8.0909E+00  2.9372E+01  1.9870E+01  2.1602E+01
             8.0262E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2395.27306870357        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0050E+00  4.4096E-01  2.9862E-01  1.3761E+00  3.3372E-01  1.0167E+00  7.6845E-01  2.2829E+00  1.0972E+00  4.6303E-01
             4.5180E+00
 PARAMETER:  1.0499E-01 -7.1880E-01 -1.1086E+00  4.1923E-01 -9.9746E-01  1.1655E-01 -1.6338E-01  9.2545E-01  1.9272E-01 -6.6995E-01
             1.6081E+00
 GRADIENT:  -1.1992E+02  2.2437E+02  1.1476E+01  2.4824E+02 -1.5477E+02  1.2740E+01 -9.7857E+00  8.9198E+01 -2.1601E+01  9.1093E-01
             6.9186E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2592.88030561654        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      322
 NPARAMETR:  9.7094E-01  3.4946E-01  2.8747E-01  1.4367E+00  3.0861E-01  1.1970E+00  8.8675E-01  1.9629E+00  1.9572E+00  2.5440E-01
             3.1294E+00
 PARAMETER:  7.0512E-02 -9.5135E-01 -1.1466E+00  4.6234E-01 -1.0757E+00  2.7982E-01 -2.0188E-02  7.7440E-01  7.7154E-01 -1.2689E+00
             1.2408E+00
 GRADIENT:  -1.0486E+02  5.6605E+01  2.6064E+01  1.1762E+02  6.8935E-01  5.7228E+01 -1.3857E+01  6.0268E+01  7.8461E+01 -6.0936E-01
             3.2628E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2651.54648530301        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      499
 NPARAMETR:  1.0269E+00  4.3980E-01  3.3265E-01  1.2292E+00  3.4717E-01  1.0686E+00  1.0210E+00  1.7326E+00  1.4008E+00  5.5704E-01
             2.5894E+00
 PARAMETER:  1.2658E-01 -7.2143E-01 -1.0007E+00  3.0633E-01 -9.5794E-01  1.6639E-01  1.2083E-01  6.4964E-01  4.3704E-01 -4.8511E-01
             1.0514E+00
 GRADIENT:  -1.3423E+01  1.7944E+02  5.9934E+01  6.2683E+01 -9.5679E+01  3.2018E+01 -2.6117E+01  6.3339E+01  4.0824E+01  1.0794E+00
            -3.3527E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2659.66382731874        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      685             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0228E+00  3.7552E-01  2.5413E-01  1.2249E+00  2.8000E-01  1.0562E+00  1.1048E+00  1.9658E+00  1.4736E+00  5.6151E-01
             2.4931E+00
 PARAMETER:  1.2251E-01 -8.7945E-01 -1.2699E+00  3.0288E-01 -1.1730E+00  1.5472E-01  1.9964E-01  7.7592E-01  4.8774E-01 -4.7712E-01
             1.0135E+00
 GRADIENT:   4.9198E+01  2.7474E+02  1.2266E+02  1.1670E+02  1.4934E+02  3.3745E+01 -2.9655E+01  4.7625E+01  4.2462E+01  7.0527E+00
             2.2496E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2671.73037469245        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      842
 NPARAMETR:  1.0499E+00  3.7275E-01  2.5351E-01  1.2247E+00  2.8266E-01  9.7866E-01  1.1048E+00  1.9624E+00  1.3184E+00  4.5889E-01
             2.4934E+00
 PARAMETER:  1.4868E-01 -8.8684E-01 -1.2723E+00  3.0270E-01 -1.1635E+00  7.8428E-02  1.9966E-01  7.7419E-01  3.7640E-01 -6.7893E-01
             1.0137E+00
 GRADIENT:   3.0701E+01  2.1100E+02  5.5467E+01  1.0732E+02 -1.3597E+02 -8.4646E-01 -3.5398E+01  7.4560E+01 -1.4353E+00 -2.0525E+00
            -3.2458E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2688.43903713078        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1021
 NPARAMETR:  1.0322E+00  3.4273E-01  2.4872E-01  1.2202E+00  2.9645E-01  1.0858E+00  1.1052E+00  1.8248E+00  1.2113E+00  5.2501E-01
             2.5242E+00
 PARAMETER:  1.3165E-01 -9.7082E-01 -1.2914E+00  2.9899E-01 -1.1159E+00  1.8228E-01  2.0005E-01  7.0146E-01  2.9165E-01 -5.4434E-01
             1.0259E+00
 GRADIENT:  -7.2196E+00  5.6669E+01  1.3063E+01  1.0936E+02  5.6377E+01  3.5127E+01 -1.1217E+01  7.0385E+01 -1.7103E+01  4.9376E+00
            -9.2293E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2713.07098959825        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1199
 NPARAMETR:  1.0482E+00  2.8489E-01  2.3662E-01  1.2029E+00  2.6294E-01  9.3219E-01  1.1066E+00  1.3404E+00  1.3271E+00  4.4834E-01
             2.6199E+00
 PARAMETER:  1.4706E-01 -1.1556E+00 -1.3413E+00  2.8477E-01 -1.2358E+00  2.9785E-02  2.0127E-01  3.9300E-01  3.8298E-01 -7.0220E-01
             1.0631E+00
 GRADIENT:   2.5006E+01 -3.0437E+01  9.3609E+01  1.0283E+02 -2.6019E+01 -1.8472E+01 -1.6898E+01  1.9799E+01  9.1748E-02 -1.8258E+00
             7.9874E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2732.65408609563        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1376
 NPARAMETR:  1.0315E+00  2.6463E-01  1.8779E-01  1.1742E+00  2.3810E-01  9.8418E-01  1.1086E+00  9.8015E-01  1.4707E+00  5.6230E-01
             2.3809E+00
 PARAMETER:  1.3101E-01 -1.2294E+00 -1.5724E+00  2.6062E-01 -1.3351E+00  8.4058E-02  2.0309E-01  7.9949E-02  4.8573E-01 -4.7571E-01
             9.6749E-01
 GRADIENT:  -3.4930E+00 -1.7852E+01  5.5245E+01  9.7646E+01  2.2283E+00  6.4375E-01 -1.4366E+01 -1.0854E+01  1.0400E+00  1.4722E+00
            -8.2854E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2738.19105342088        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1553
 NPARAMETR:  1.0345E+00  2.5605E-01  1.6274E-01  1.1546E+00  2.2246E-01  9.8257E-01  1.1099E+00  8.5686E-01  1.5473E+00  5.9573E-01
             2.4722E+00
 PARAMETER:  1.3389E-01 -1.2624E+00 -1.7156E+00  2.4376E-01 -1.4030E+00  8.2418E-02  2.0431E-01 -5.4484E-02  5.3649E-01 -4.1797E-01
             1.0051E+00
 GRADIENT:   5.1592E-01  1.5165E+01  1.0817E+00  9.8220E+01  8.1317E+00 -4.0151E-02 -9.4312E+00 -1.6370E+01  2.2870E-01  9.5240E-01
             7.6560E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2751.62753077658        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1747             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0231E+00  2.0464E-01  1.2098E-01  9.4412E-01  1.8578E-01  9.8188E-01  1.1182E+00  1.0298E+00  1.7007E+00  6.4346E-01
             2.3577E+00
 PARAMETER:  1.2285E-01 -1.4865E+00 -2.0121E+00  4.2500E-02 -1.5832E+00  8.1713E-02  2.1174E-01  1.2938E-01  6.3106E-01 -3.4089E-01
             9.5769E-01
 GRADIENT:   5.0632E+01  3.6406E+01  9.2407E+01  1.0077E+01  5.6723E+02  4.9641E+00 -1.4416E+01 -1.7608E+01  1.6881E+01  3.9948E+00
            -1.8429E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2751.78562606193        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1884
 NPARAMETR:  1.0232E+00  2.0492E-01  1.2085E-01  9.4724E-01  1.8501E-01  9.8509E-01  1.1185E+00  1.0297E+00  1.7015E+00  6.4726E-01
             2.3589E+00
 PARAMETER:  1.2295E-01 -1.4851E+00 -2.0132E+00  4.5794E-02 -1.5874E+00  8.4981E-02  2.1196E-01  1.2930E-01  6.3149E-01 -3.3501E-01
             9.5818E-01
 GRADIENT:  -1.7844E+01 -3.6313E+01  7.7232E+00  3.0068E+00  1.2352E+01  1.7267E+00 -1.5423E+01 -1.8016E+01 -1.9538E-01  4.1917E-01
            -3.5160E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2752.65804714387        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2061
 NPARAMETR:  1.0240E+00  2.0948E-01  1.1840E-01  9.3794E-01  1.8389E-01  9.7816E-01  1.1203E+00  1.0294E+00  1.7059E+00  6.4455E-01
             2.3967E+00
 PARAMETER:  1.2376E-01 -1.4631E+00 -2.0337E+00  3.5933E-02 -1.5934E+00  7.7915E-02  2.1364E-01  1.2898E-01  6.3410E-01 -3.3920E-01
             9.7408E-01
 GRADIENT:  -1.7622E+01 -3.6283E+00 -5.9350E+00 -1.7733E+00  1.9534E+00 -7.7035E-01 -1.2603E+01 -1.6853E+01 -5.0478E-01 -2.9717E-01
            -2.5414E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2754.59489404038        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     2199
 NPARAMETR:  1.0291E+00  2.0998E-01  1.1883E-01  9.4052E-01  1.8362E-01  9.8018E-01  1.1751E+00  1.1644E+00  1.7087E+00  6.4511E-01
             2.3979E+00
 PARAMETER:  1.2871E-01 -1.4607E+00 -2.0301E+00  3.8672E-02 -1.5949E+00  7.9985E-02  2.6136E-01  2.5218E-01  6.3576E-01 -3.3833E-01
             9.7461E-01
 GRADIENT:  -6.4626E+00  2.6560E+00  1.8987E+00 -2.4258E+00 -1.0281E+01  4.8029E-02 -4.1394E+00 -7.3516E+00  1.1566E+00  3.2414E+00
             1.2065E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2754.64973611683        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     2360
 NPARAMETR:  1.0298E+00  2.1000E-01  1.1999E-01  9.4005E-01  1.8507E-01  9.8061E-01  1.2013E+00  1.1629E+00  1.7033E+00  6.2428E-01
             2.4081E+00
 PARAMETER:  1.2936E-01 -1.4606E+00 -2.0203E+00  3.8175E-02 -1.5870E+00  8.0424E-02  2.8338E-01  2.5094E-01  6.3256E-01 -3.7115E-01
             9.7885E-01
 GRADIENT:  -5.5949E+00 -3.8999E+00  4.2610E-01 -5.0978E+00  5.6285E-01  2.3423E-01 -6.6638E-01 -7.1846E+00  1.6234E+00  3.1597E-01
             2.0304E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2754.94345697630        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     2521
 NPARAMETR:  1.0302E+00  2.1027E-01  1.2019E-01  9.4263E-01  1.8492E-01  9.8038E-01  1.2033E+00  1.2010E+00  1.7008E+00  6.2279E-01
             2.3989E+00
 PARAMETER:  1.2979E-01 -1.4593E+00 -2.0187E+00  4.0913E-02 -1.5878E+00  8.0190E-02  2.8509E-01  2.8316E-01  6.3107E-01 -3.7355E-01
             9.7501E-01
 GRADIENT:  -4.2149E+00 -1.6943E+00  3.6508E+00 -4.3768E+00 -5.5674E+00  1.3708E-01 -6.1455E-01 -5.4988E+00  9.4940E-01  4.9903E-01
             1.6256E+01

0ITERATION NO.:   83    OBJECTIVE VALUE:  -2754.99551221353        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     2622
 NPARAMETR:  1.0302E+00  2.1028E-01  1.2018E-01  9.4563E-01  1.8492E-01  9.8030E-01  1.2050E+00  1.2006E+00  1.6994E+00  6.2185E-01
             2.3916E+00
 PARAMETER:  1.2979E-01 -1.4593E+00 -2.0188E+00  4.4094E-02 -1.5878E+00  8.0100E-02  2.8647E-01  2.8281E-01  6.3029E-01 -3.7506E-01
             9.7196E-01
 GRADIENT:  -2.8230E+00 -4.1102E+00  2.2629E+01  5.6581E+03  1.0112E+00  5.3676E-02  1.9780E+03 -4.3498E+01  8.9224E+02  1.3889E-01
            -5.6106E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2622
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3515E-03  1.7301E-02  1.9286E-02 -1.5304E-03  2.7567E-02
 SE:             2.9424E-02  2.4158E-02  2.1262E-02  2.7077E-02  2.1253E-02
 N:                     100         100         100         100         100

 P VAL.:         9.0931E-01  4.7390E-01  3.6438E-01  9.5493E-01  1.9461E-01

 ETASHRINKSD(%)  1.4275E+00  1.9068E+01  2.8768E+01  9.2885E+00  2.8799E+01
 ETASHRINKVR(%)  2.8346E+00  3.4500E+01  4.9260E+01  1.7714E+01  4.9305E+01
 EBVSHRINKSD(%)  1.6351E+00  1.8500E+01  2.8225E+01  6.8825E+00  2.9902E+01
 EBVSHRINKVR(%)  3.2434E+00  3.3577E+01  4.8484E+01  1.3291E+01  5.0862E+01
 RELATIVEINF(%)  9.6677E+01  2.2507E+01  1.1284E+01  5.6491E+01  8.3357E+00
 EPSSHRINKSD(%)  2.1009E+01
 EPSSHRINKVR(%)  3.7605E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2754.9955122135339     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1100.9061524451231     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    88.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2754.996       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.10E-01  1.20E-01  9.46E-01  1.85E-01  9.80E-01  1.20E+00  1.20E+00  1.70E+00  6.22E-01  2.39E+00
 


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
+        1.59E+06
 
 TH 2
+        4.63E+00  3.13E+05
 
 TH 3
+       -7.88E+00 -1.83E+03  2.49E+05
 
 TH 4
+        4.59E+00 -1.21E+02 -2.52E+02  1.58E+06
 
 TH 5
+       -2.94E+01 -9.55E+03 -2.20E+04 -5.64E+02  2.05E+05
 
 TH 6
+        6.53E+00 -6.76E+00  9.22E+00  3.23E+02 -9.99E+00  1.93E+02
 
 TH 7
+        3.07E+05  8.17E+01  2.00E+01  4.33E+05 -2.32E+01  8.87E+01  1.19E+05
 
 TH 8
+       -3.15E+05  2.91E+01  7.54E+01  4.38E+05  2.35E+01 -1.23E-01 -2.26E+03  2.46E+05
 
 TH 9
+        9.81E+04  5.62E+00  1.36E+02  1.39E+05  1.70E+02  3.00E+01  3.80E+04 -7.22E+02  1.22E+04
 
 TH10
+       -3.39E-01  6.69E+01  1.27E+02  1.46E+03  3.50E+02 -7.84E-02  4.14E+02 -1.77E+05  1.36E+02  1.58E+02
 
 TH11
+       -1.97E+01 -3.34E+01 -2.66E+04  2.07E+02 -2.13E+04 -1.15E+01  1.69E+01  3.58E+02  1.02E+02 -5.21E+01  2.69E+03
 
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
 #CPUT: Total CPU Time in Seconds,      107.302
Stop Time:
Wed Sep 29 00:12:11 CDT 2021
