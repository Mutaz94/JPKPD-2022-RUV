Wed Sep 29 20:01:03 CDT 2021
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
$DATA ../../../../data/spa/D/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12249.4004084698        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9562E+02  2.4047E+02 -2.1571E+01  2.3585E+02  9.2039E+01 -1.3531E+03 -7.1864E+02 -8.6652E+01 -9.8759E+02 -4.1929E+02
            -2.3953E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -646.842092664352        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4274E+00  1.1404E+00  9.9229E-01  1.4216E+00  1.1935E+00  1.7005E+00  1.2697E+00  9.8731E-01  1.1781E+00  1.0156E+00
             1.4810E+01
 PARAMETER:  4.5583E-01  2.3141E-01  9.2265E-02  4.5178E-01  2.7690E-01  6.3093E-01  3.3878E-01  8.7228E-02  2.6387E-01  1.1550E-01
             2.7953E+00
 GRADIENT:   4.0634E+01  3.1081E+00 -2.2842E+00  1.7580E+00 -6.2712E+00  5.1908E+01  2.0485E+00  3.2572E+00  1.1484E+01  3.2192E+00
             1.7081E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -667.573515947735        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2872E+00  7.6043E-01  1.3853E+00  1.6531E+00  2.0417E+00  1.4921E+00  2.8586E+00  4.0000E-01  1.0200E+00  2.0697E+00
             1.3492E+01
 PARAMETER:  3.5247E-01 -1.7388E-01  4.2591E-01  6.0265E-01  8.1377E-01  5.0017E-01  1.1503E+00 -8.1629E-01  1.1983E-01  8.2739E-01
             2.7021E+00
 GRADIENT:  -3.5740E+00  1.5515E+01  3.2210E+00  2.4874E+01 -4.2242E+00  1.2995E+01  1.0241E+01  2.4933E-01  1.2188E+01  2.7666E-01
             1.4263E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -687.900370869580        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0975E+00  5.9306E-01  1.1007E+00  1.4056E+00  3.5319E+00  1.3447E+00  2.3577E+00  6.9043E-01  7.7133E-01  8.4253E+00
             1.1805E+01
 PARAMETER:  1.9301E-01 -4.2246E-01  1.9592E-01  4.4049E-01  1.3618E+00  3.9620E-01  9.5770E-01 -2.7043E-01 -1.5964E-01  2.2312E+00
             2.5685E+00
 GRADIENT:  -5.2355E+01  9.2364E+00  6.1512E+00 -6.7700E+00 -1.1105E+01  1.0485E+00  5.8240E+00  3.0334E-01  1.1399E+01  3.9583E+00
             9.8220E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -707.196944864548        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.8503E-01  2.0237E-01  4.2169E-01  1.1819E+00  1.1045E+01  1.1452E+00  1.7620E+00  2.2263E-01  2.0666E-01  4.8211E+00
             9.9323E+00
 PARAMETER:  8.4920E-02 -1.4977E+00 -7.6349E-01  2.6714E-01  2.5020E+00  2.3558E-01  6.6645E-01 -1.4022E+00 -1.4767E+00  1.6730E+00
             2.3958E+00
 GRADIENT:   4.0061E+01  7.9221E-01  3.4762E+01 -5.6733E+01  1.4859E+01 -5.1266E+01  7.1995E-01  1.4841E-01  8.9235E-01 -2.0412E+01
            -6.9781E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -733.245006930304        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  7.7772E-01  8.6023E-02  8.6015E-02  7.2862E-01  1.3412E+01  1.4204E+00  6.6660E-01  1.0000E-02  2.8871E-02  2.9148E+00
             1.0943E+01
 PARAMETER: -1.5138E-01 -2.3531E+00 -2.3532E+00 -2.1661E-01  2.6962E+00  4.5096E-01 -3.0557E-01 -7.6009E+00 -3.4449E+00  1.1698E+00
             2.4927E+00
 GRADIENT:   4.6418E+01  9.4597E+01 -9.6295E+01  4.6753E+01 -1.2053E+01 -1.7278E+01  1.0246E+01  0.0000E+00  3.9080E-02  4.0353E+00
             9.7962E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -742.305000812750        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  5.9207E-01  4.0592E-02  3.9272E-02  4.5325E-01  2.5591E+01  1.3982E+00  2.2748E-01  1.0000E-02  1.0000E-02  2.3183E+00
             1.0583E+01
 PARAMETER: -4.2414E-01 -3.1042E+00 -3.1372E+00 -6.9130E-01  3.3422E+00  4.3517E-01 -1.3807E+00 -1.1011E+01 -4.9778E+00  9.4082E-01
             2.4593E+00
 GRADIENT:   4.9706E+01  1.6294E+01 -1.3564E+02  1.4212E+02 -7.9871E-01 -3.0409E+01  1.1151E+00  0.0000E+00  0.0000E+00  2.1671E-01
            -2.5848E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -757.028116077777        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      643
 NPARAMETR:  5.0746E-01  2.5040E-02  3.2877E-02  3.4896E-01  3.5268E+01  1.4741E+00  1.4452E-01  1.0000E-02  1.0000E-02  1.8552E+00
             1.0413E+01
 PARAMETER: -5.7834E-01 -3.5873E+00 -3.3150E+00 -9.5279E-01  3.6630E+00  4.8804E-01 -1.8343E+00 -1.3037E+01 -5.9890E+00  7.1800E-01
             2.4431E+00
 GRADIENT:  -1.5189E+00  2.1249E+00 -1.9555E+00  1.1453E+00 -5.0146E-02 -1.4499E+00  2.5822E-02  0.0000E+00  0.0000E+00  2.9890E-03
            -2.5928E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -757.359650069117        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  4.8670E-01  1.8648E-02  2.9679E-02  3.2104E-01  4.7104E+01  1.4831E+00  7.0870E-02  1.0000E-02  1.0000E-02  1.6501E+00
             1.0447E+01
 PARAMETER: -6.2011E-01 -3.8820E+00 -3.4173E+00 -1.0362E+00  3.9524E+00  4.9412E-01 -2.5469E+00 -1.3728E+01 -6.3510E+00  6.0083E-01
             2.4463E+00
 GRADIENT:   1.6304E-01  3.4628E-01 -2.1699E-03 -8.1238E-01  2.8970E-03  4.3180E-03  4.6702E-04  0.0000E+00  0.0000E+00 -2.5009E-05
             2.6270E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -757.430963159540        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1001
 NPARAMETR:  4.8597E-01  1.3392E-02  2.9832E-02  3.2156E-01  9.0850E+00  1.4868E+00  4.0295E-02  1.0000E-02  1.0000E-02  1.6716E+00
             1.0443E+01
 PARAMETER: -6.2161E-01 -4.2131E+00 -3.4122E+00 -1.0346E+00  2.3066E+00  4.9662E-01 -3.1115E+00 -1.3728E+01 -6.3510E+00  6.1378E-01
             2.4459E+00
 GRADIENT:  -3.2466E+00  1.3066E-01  7.3675E+00 -8.6498E+00 -4.8644E-02  2.2880E-01  1.7985E-05  0.0000E+00  0.0000E+00  2.0584E-02
             1.2585E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -757.686493240665        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1177
 NPARAMETR:  4.8705E-01  1.2383E-02  2.9745E-02  3.2228E-01  8.9371E+00  1.4898E+00  3.6809E-02  1.0000E-02  1.0000E-02  2.0705E+00
             1.0443E+01
 PARAMETER: -6.1940E-01 -4.2914E+00 -3.4151E+00 -1.0323E+00  2.2902E+00  4.9861E-01 -3.2020E+00 -1.3728E+01 -6.3510E+00  8.2780E-01
             2.4460E+00
 GRADIENT:  -2.7485E+00 -8.9182E-02  3.4935E+00 -3.2533E+00  3.4187E+00 -1.8753E-01  1.0682E-05  0.0000E+00  0.0000E+00 -2.1833E+00
             3.4801E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -757.946316926394        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     1349
 NPARAMETR:  4.8833E-01  1.1015E-02  2.9744E-02  3.2252E-01  9.6189E+00  1.4937E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.7535E+00
             1.0420E+01
 PARAMETER: -6.1676E-01 -4.4085E+00 -3.4151E+00 -1.0316E+00  2.3637E+00  5.0125E-01 -4.5112E+00 -1.3728E+01 -6.3510E+00  1.1129E+00
             2.4437E+00
 GRADIENT:   1.0738E+00 -1.3725E-01  1.9358E-01 -8.3797E-01  2.7571E+00  2.3753E-01  1.7797E-07  0.0000E+00  0.0000E+00 -1.2568E+00
             8.4091E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -757.958646317659        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1525
 NPARAMETR:  4.8778E-01  1.1017E-02  2.9824E-02  3.2288E-01  9.5754E+00  1.4928E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.7688E+00
             1.0409E+01
 PARAMETER: -6.1788E-01 -4.4083E+00 -3.4124E+00 -1.0305E+00  2.3592E+00  5.0067E-01 -8.0695E+00 -1.3728E+01 -6.3510E+00  1.1184E+00
             2.4427E+00
 GRADIENT:  -3.7681E-01  3.7843E-01  7.0342E-01 -1.2961E+00  5.4181E-01 -3.3882E-02  0.0000E+00  0.0000E+00  0.0000E+00 -4.9406E-01
            -1.8289E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -757.989797029107        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1706
 NPARAMETR:  4.8678E-01  1.0401E-02  2.9749E-02  3.2286E-01  9.5238E+00  1.4952E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.8220E+00
             1.0402E+01
 PARAMETER: -6.1994E-01 -4.4659E+00 -3.4150E+00 -1.0306E+00  2.3538E+00  5.0228E-01 -1.5503E+01 -1.3728E+01 -6.3510E+00  1.1375E+00
             2.4420E+00
 GRADIENT:  -4.7877E-01  2.8882E+00 -8.8159E+00  8.8629E+00 -1.2561E+01 -8.6448E-01  0.0000E+00  0.0000E+00  0.0000E+00  3.8597E+00
            -3.3432E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -758.014314846115        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1885
 NPARAMETR:  4.8416E-01  1.0000E-02  3.0277E-02  3.1913E-01  9.6283E+00  1.5015E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.8255E+00
             1.0605E+01
 PARAMETER: -6.1939E-01 -4.5063E+00 -3.4168E+00 -1.0332E+00  2.3498E+00  5.0147E-01 -1.9455E+01 -1.3728E+01 -6.3510E+00  1.1501E+00
             2.4427E+00
 GRADIENT:   8.1432E+01  1.5027E+01 -3.3054E+01  8.3579E+01 -4.2040E+01 -1.0091E+02  0.0000E+00  0.0000E+00  0.0000E+00  7.8655E+01
            -3.5523E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1885
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4468E-03  3.7284E-06  1.0520E-04 -2.2378E-04 -5.4665E-03
 SE:             2.8844E-02  1.3514E-06  3.1462E-04  4.0467E-04  3.8359E-03
 N:                     100         100         100         100         100

 P VAL.:         9.3240E-01  5.7982E-03  7.3811E-01  5.8027E-01  1.5413E-01

 ETASHRINKSD(%)  3.3678E+00  9.9995E+01  9.8946E+01  9.8644E+01  8.7149E+01
 ETASHRINKVR(%)  6.6222E+00  1.0000E+02  9.9989E+01  9.9982E+01  9.8349E+01
 EBVSHRINKSD(%)  3.5520E+00  9.9995E+01  9.8970E+01  9.8678E+01  8.9970E+01
 EBVSHRINKVR(%)  6.9777E+00  1.0000E+02  9.9989E+01  9.9983E+01  9.8994E+01
 RELATIVEINF(%)  5.5382E+00  1.5494E-08  9.1525E-05  1.4945E-04  2.2901E-01
 EPSSHRINKSD(%)  7.4049E+00
 EPSSHRINKVR(%)  1.4262E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -758.01431484611544     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -22.863488282377261     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.17
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -758.014       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.87E-01  1.00E-02  2.97E-02  3.22E-01  9.49E+00  1.49E+00  1.00E-02  1.00E-02  1.00E-02  2.86E+00  1.04E+01
 


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
+        2.68E+04
 
 TH 2
+       -1.28E+04  1.65E+06
 
 TH 3
+       -5.73E+03  9.39E+03  1.03E+06
 
 TH 4
+       -7.91E+02 -8.71E+04 -7.75E+04  3.00E+04
 
 TH 5
+        1.23E+01 -1.29E+03  3.96E+02 -1.16E+02  1.69E+01
 
 TH 6
+       -1.00E+04  4.57E+03  6.84E+02 -6.50E+01 -7.27E-01  2.32E+03
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.38E+01 -7.86E+03  3.50E+01  1.87E+03 -1.80E+01  5.08E+00  0.00E+00  0.00E+00  0.00E+00  3.17E+02
 
 TH11
+       -1.73E+01  1.13E+03  1.63E+02 -1.20E+02  1.34E+00  2.86E-01  0.00E+00  0.00E+00  0.00E+00 -1.30E+01  8.89E+00
 
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
 #CPUT: Total CPU Time in Seconds,       32.531
Stop Time:
Wed Sep 29 20:01:37 CDT 2021
