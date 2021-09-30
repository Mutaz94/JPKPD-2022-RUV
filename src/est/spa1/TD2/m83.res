Thu Sep 30 02:23:55 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat83.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2101.73527141614        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4203E+02  5.5763E+01  2.9973E+01  5.8839E+01 -4.2887E+01  2.5996E+01  1.5340E+01 -3.8770E+00  2.7218E+01  2.8799E+01
            -5.6065E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2108.90932018924        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0197E+00  9.8214E-01  9.9435E-01  1.0266E+00  1.0012E+00  1.0696E+00  9.3129E-01  1.0154E+00  8.9433E-01  8.4786E-01
             1.0759E+00
 PARAMETER:  1.1949E-01  8.1975E-02  9.4333E-02  1.2624E-01  1.0124E-01  1.6732E-01  2.8813E-02  1.1532E-01 -1.1680E-02 -6.5043E-02
             1.7313E-01
 GRADIENT:   3.4154E+01  2.3404E+01  1.7674E+01  3.7087E+00 -1.5219E+01  8.8812E+00  1.8916E+00 -4.7687E+00 -5.0010E+00  8.7958E+00
             2.8275E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2111.65847325053        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  1.0032E+00  7.7489E-01  9.3121E-01  1.1355E+00  8.6973E-01  1.0550E+00  8.7980E-01  1.0019E+00  8.7324E-01  6.6370E-01
             1.0574E+00
 PARAMETER:  1.0318E-01 -1.5503E-01  2.8733E-02  2.2711E-01 -3.9577E-02  1.5358E-01 -2.8059E-02  1.0187E-01 -3.5539E-02 -3.0992E-01
             1.5580E-01
 GRADIENT:   4.9640E+00  4.6807E+00  9.4080E+00 -8.6587E+00 -1.2607E+01  4.7295E+00 -3.4818E+00  1.9981E+00  1.3935E+00 -5.0247E+00
            -1.0093E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2113.17350827987        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.0003E+00  6.9772E-01  8.3090E-01  1.1798E+00  7.8660E-01  1.0418E+00  1.2038E+00  7.8625E-01  8.0457E-01  6.3978E-01
             1.0698E+00
 PARAMETER:  1.0029E-01 -2.5994E-01 -8.5247E-02  2.6531E-01 -1.4004E-01  1.4099E-01  2.8552E-01 -1.4049E-01 -1.1745E-01 -3.4662E-01
             1.6746E-01
 GRADIENT:  -4.9265E-01  1.0986E+01  8.0971E+00  1.0849E+01 -1.2153E+01 -1.7942E-01 -8.9121E-01 -1.3057E+00 -8.8826E-01 -1.3094E-01
             2.1079E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2113.76133865208        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.9684E-01  5.2043E-01  9.2723E-01  1.2893E+00  7.7875E-01  1.0390E+00  1.4874E+00  8.8450E-01  7.5426E-01  6.4549E-01
             1.0672E+00
 PARAMETER:  9.6833E-02 -5.5311E-01  2.4441E-02  3.5408E-01 -1.5006E-01  1.3822E-01  4.9705E-01 -2.2735E-02 -1.8202E-01 -3.3775E-01
             1.6506E-01
 GRADIENT:  -2.0827E-01  6.6188E+00  6.2290E+00  7.5214E+00 -1.0614E+01 -1.9705E-01  6.2506E-01  6.2719E-02 -8.8686E-01 -2.7513E-01
             6.1300E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2114.04724604225        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.9197E-01  3.4407E-01  1.0756E+00  1.4064E+00  8.0309E-01  1.0366E+00  1.6822E+00  1.0289E+00  7.2204E-01  6.8763E-01
             1.0668E+00
 PARAMETER:  9.1937E-02 -9.6692E-01  1.7287E-01  4.4100E-01 -1.1928E-01  1.3590E-01  6.2010E-01  1.2848E-01 -2.2568E-01 -2.7451E-01
             1.6467E-01
 GRADIENT:  -1.5292E+00  3.9993E+00  3.8374E+00  1.0739E+01 -5.0032E+00 -3.0672E-02  4.4515E-02 -3.8592E-01 -1.4752E+00 -9.0077E-02
             2.7042E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2114.10256525097        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  9.9090E-01  2.8538E-01  1.0995E+00  1.4406E+00  7.9914E-01  1.0358E+00  1.7735E+00  1.0605E+00  7.1350E-01  6.9046E-01
             1.0663E+00
 PARAMETER:  9.0861E-02 -1.1539E+00  1.9489E-01  4.6508E-01 -1.2421E-01  1.3513E-01  6.7298E-01  1.5870E-01 -2.3758E-01 -2.7040E-01
             1.6419E-01
 GRADIENT:  -8.5811E-01  2.4124E+00  2.1038E+00  7.1723E+00 -2.8689E+00  1.0644E-03 -1.5513E-01 -2.6996E-01 -8.8013E-01 -7.5164E-02
             2.0850E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2114.14130359254        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.9137E-01  2.7377E-01  1.1011E+00  1.4418E+00  7.9882E-01  1.0358E+00  1.8269E+00  1.0653E+00  7.1375E-01  6.9134E-01
             1.0662E+00
 PARAMETER:  9.1332E-02 -1.1955E+00  1.9635E-01  4.6592E-01 -1.2461E-01  1.3518E-01  7.0263E-01  1.6321E-01 -2.3722E-01 -2.6913E-01
             1.6410E-01
 GRADIENT:   7.7638E-01  2.0671E-01 -1.2930E-01 -7.9169E+00  2.1120E+00  1.0708E-01  1.7187E-02  2.2149E-02  4.6182E-01 -5.0016E-02
             1.9029E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2114.14462719160        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1426             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9165E-01  2.7388E-01  1.0995E+00  1.4416E+00  7.9684E-01  1.0359E+00  1.8381E+00  1.0633E+00  7.1294E-01  6.9072E-01
             1.0659E+00
 PARAMETER:  9.1618E-02 -1.1951E+00  1.9483E-01  4.6573E-01 -1.2711E-01  1.3529E-01  7.0871E-01  1.6134E-01 -2.3836E-01 -2.7003E-01
             1.6383E-01
 GRADIENT:   3.8120E+02  3.7171E+01  7.3763E+00  7.0179E+02  8.3070E+00  4.9451E+01  5.1550E+00  6.3193E-01  2.3245E+01  1.1330E+00
             1.6457E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2114.14646371082        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  9.9165E-01  2.7373E-01  1.0978E+00  1.4415E+00  7.9607E-01  1.0359E+00  1.8420E+00  1.0618E+00  7.1284E-01  6.9019E-01
             1.0659E+00
 PARAMETER:  9.1616E-02 -1.1956E+00  1.9332E-01  4.6566E-01 -1.2807E-01  1.3529E-01  7.1086E-01  1.5994E-01 -2.3850E-01 -2.7080E-01
             1.6382E-01
 GRADIENT:   1.3500E+00  4.7571E-01  1.4499E+00 -7.8546E+00 -7.5398E-01  1.5513E-01  3.1334E-02  3.8962E-03  1.7388E-01  1.2575E-01
             5.6534E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2114.14838237423        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1787             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9166E-01  2.7286E-01  1.0941E+00  1.4411E+00  7.9557E-01  1.0359E+00  1.8458E+00  1.0603E+00  7.1275E-01  6.8889E-01
             1.0659E+00
 PARAMETER:  9.1624E-02 -1.1988E+00  1.8997E-01  4.6544E-01 -1.2870E-01  1.3529E-01  7.1294E-01  1.5859E-01 -2.3862E-01 -2.7268E-01
             1.6381E-01
 GRADIENT:   3.8115E+02  3.6617E+01  4.9388E+00  7.0016E+02  1.1814E+01  4.9424E+01  5.1965E+00  7.8869E-01  2.3207E+01  9.7897E-01
             1.6312E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2114.14969125263        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1968
 NPARAMETR:  9.9166E-01  2.7274E-01  1.0927E+00  1.4410E+00  7.9485E-01  1.0359E+00  1.8486E+00  1.0590E+00  7.1273E-01  6.8834E-01
             1.0659E+00
 PARAMETER:  9.1623E-02 -1.1992E+00  1.8862E-01  4.6536E-01 -1.2960E-01  1.3528E-01  7.1443E-01  1.5732E-01 -2.3866E-01 -2.7347E-01
             1.6381E-01
 GRADIENT:   1.4262E+00  1.7880E-03 -8.1898E-01 -9.2617E+00  2.5359E+00  1.6146E-01  1.7616E-02  1.6782E-01  1.3755E-01 -4.3843E-02
            -7.7730E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2114.15164426199        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2151             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9164E-01  2.7337E-01  1.0920E+00  1.4411E+00  7.9312E-01  1.0359E+00  1.8537E+00  1.0563E+00  7.1272E-01  6.8807E-01
             1.0659E+00
 PARAMETER:  9.1606E-02 -1.1969E+00  1.8798E-01  4.6540E-01 -1.3178E-01  1.3527E-01  7.1720E-01  1.5480E-01 -2.3866E-01 -2.7386E-01
             1.6381E-01
 GRADIENT:   3.8098E+02  3.7161E+01  7.1616E+00  7.0135E+02  8.3214E+00  4.9400E+01  5.3052E+00  5.8573E-01  2.3246E+01  1.1503E+00
             1.6446E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2114.15249491678        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2329
 NPARAMETR:  9.9164E-01  2.7335E-01  1.0907E+00  1.4410E+00  7.9245E-01  1.0359E+00  1.8561E+00  1.0550E+00  7.1272E-01  6.8761E-01
             1.0659E+00
 PARAMETER:  9.1604E-02 -1.1970E+00  1.8681E-01  4.6533E-01 -1.3262E-01  1.3526E-01  7.1849E-01  1.5358E-01 -2.3867E-01 -2.7453E-01
             1.6381E-01
 GRADIENT:   1.3339E+00  5.2039E-01  1.5470E+00 -7.6169E+00 -1.3870E+00  1.5311E-01  4.0428E-02 -2.9954E-02  1.6297E-01  1.3285E-01
             5.9442E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2114.15371232275        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2512             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9165E-01  2.7257E-01  1.0875E+00  1.4406E+00  7.9234E-01  1.0359E+00  1.8576E+00  1.0541E+00  7.1269E-01  6.8650E-01
             1.0659E+00
 PARAMETER:  9.1618E-02 -1.1999E+00  1.8386E-01  4.6507E-01 -1.3277E-01  1.3527E-01  7.1929E-01  1.5269E-01 -2.3871E-01 -2.7614E-01
             1.6380E-01
 GRADIENT:   3.8095E+02  3.6609E+01  4.6737E+00  6.9953E+02  1.1948E+01  4.9380E+01  5.3113E+00  7.3744E-01  2.3211E+01  9.7906E-01
             1.6297E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2114.15438650527        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2693
 NPARAMETR:  9.9165E-01  2.7258E-01  1.0864E+00  1.4405E+00  7.9183E-01  1.0359E+00  1.8593E+00  1.0531E+00  7.1269E-01  6.8612E-01
             1.0659E+00
 PARAMETER:  9.1618E-02 -1.1998E+00  1.8286E-01  4.6499E-01 -1.3340E-01  1.3526E-01  7.2022E-01  1.5171E-01 -2.3871E-01 -2.7671E-01
             1.6380E-01
 GRADIENT:   1.4148E+00  5.3146E-03 -9.3873E-01 -9.0912E+00  2.3105E+00  1.5996E-01  2.1663E-02  1.4039E-01  1.2914E-01 -5.3384E-02
            -9.0181E-03

0ITERATION NO.:   76    OBJECTIVE VALUE:  -2114.15438650527        NO. OF FUNC. EVALS.:  26
 CUMULATIVE NO. OF FUNC. EVALS.:     2719
 NPARAMETR:  9.9164E-01  2.7338E-01  1.0874E+00  1.4406E+00  7.9078E-01  1.0359E+00  1.8618E+00  1.0518E+00  7.1272E-01  6.8643E-01
             1.0659E+00
 PARAMETER:  9.1618E-02 -1.1998E+00  1.8286E-01  4.6499E-01 -1.3340E-01  1.3526E-01  7.2022E-01  1.5171E-01 -2.3871E-01 -2.7671E-01
             1.6380E-01
 GRADIENT:   4.0079E-02 -2.5223E-01 -5.3544E-01 -4.8319E-01  1.6963E+00  4.2932E-03 -1.5730E-02  9.4538E-02 -1.6918E-02 -5.6709E-02
            -7.8957E-03
 NUMSIGDIG:         3.7         2.5         2.2         3.6         1.9         4.1         2.6         1.9         3.6         2.7
                    4.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2719
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3730E-04  3.2574E-03 -2.2951E-02 -4.4822E-03 -2.7317E-02
 SE:             2.9870E-02  9.6245E-03  1.9659E-02  2.8170E-02  1.8676E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9099E-01  7.3502E-01  2.4304E-01  8.7358E-01  1.4355E-01

 ETASHRINKSD(%)  1.0000E-10  6.7757E+01  3.4139E+01  5.6281E+00  3.7433E+01
 ETASHRINKVR(%)  1.0000E-10  8.9604E+01  5.6623E+01  1.0939E+01  6.0854E+01
 EBVSHRINKSD(%)  3.4754E-01  6.8832E+01  3.4939E+01  5.9110E+00  3.7021E+01
 EBVSHRINKVR(%)  6.9387E-01  9.0286E+01  5.7670E+01  1.1473E+01  6.0336E+01
 RELATIVEINF(%)  9.5337E+01  3.8388E-01  7.3816E+00  4.5897E+00  4.4319E+00
 EPSSHRINKSD(%)  3.3519E+01
 EPSSHRINKVR(%)  5.5802E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2114.1543865052677     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1195.2158533005950     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.94
 Elapsed covariance  time in seconds:     7.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2114.154       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.92E-01  2.73E-01  1.09E+00  1.44E+00  7.92E-01  1.04E+00  1.86E+00  1.05E+00  7.13E-01  6.86E-01  1.07E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.20E-02  1.46E-01  2.71E-01  8.45E-02  1.41E-01  8.02E-02  7.91E-01  2.58E-01  7.27E-02  1.84E-01  5.87E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.03E-03
 
 TH 2
+        1.63E-03  2.12E-02
 
 TH 3
+        1.33E-04  5.11E-03  7.34E-02
 
 TH 4
+       -8.67E-04 -1.06E-02  4.75E-03  7.14E-03
 
 TH 5
+        3.80E-04  8.52E-03  3.58E-02 -8.37E-04  1.98E-02
 
 TH 6
+       -8.87E-04 -3.19E-03  3.15E-03  1.80E-03  7.32E-04  6.43E-03
 
 TH 7
+       -5.63E-03 -8.92E-02 -7.25E-02  3.84E-02 -5.99E-02  9.93E-03  6.25E-01
 
 TH 8
+        1.03E-04  5.68E-03  5.73E-02  3.21E-03  2.86E-02  1.96E-03 -5.85E-02  6.66E-02
 
 TH 9
+        2.66E-04  6.56E-03 -1.37E-04 -3.87E-03  2.00E-03 -1.31E-03 -2.20E-02 -1.75E-03  5.29E-03
 
 TH10
+        3.03E-04  4.12E-03  3.78E-02  1.38E-03  1.91E-02  5.57E-04 -6.31E-02  2.85E-02  1.99E-03  3.40E-02
 
 TH11
+        7.81E-05 -4.33E-04 -4.40E-04  1.90E-04 -1.22E-04 -7.86E-05 -1.81E-03 -1.66E-03 -3.77E-04  7.98E-04  3.45E-03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.20E-02
 
 TH 2
+        3.50E-01  1.46E-01
 
 TH 3
+        1.53E-02  1.29E-01  2.71E-01
 
 TH 4
+       -3.20E-01 -8.65E-01  2.07E-01  8.45E-02
 
 TH 5
+        8.44E-02  4.15E-01  9.39E-01 -7.04E-02  1.41E-01
 
 TH 6
+       -3.45E-01 -2.73E-01  1.45E-01  2.66E-01  6.49E-02  8.02E-02
 
 TH 7
+       -2.22E-01 -7.74E-01 -3.38E-01  5.75E-01 -5.38E-01  1.57E-01  7.91E-01
 
 TH 8
+        1.25E-02  1.51E-01  8.20E-01  1.47E-01  7.87E-01  9.50E-02 -2.87E-01  2.58E-01
 
 TH 9
+        1.14E-01  6.19E-01 -6.94E-03 -6.30E-01  1.96E-01 -2.25E-01 -3.83E-01 -9.33E-02  7.27E-02
 
 TH10
+        5.13E-02  1.53E-01  7.58E-01  8.89E-02  7.35E-01  3.77E-02 -4.33E-01  5.99E-01  1.49E-01  1.84E-01
 
 TH11
+        4.15E-02 -5.05E-02 -2.76E-02  3.83E-02 -1.48E-02 -1.67E-02 -3.90E-02 -1.09E-01 -8.82E-02  7.37E-02  5.87E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.32E+03
 
 TH 2
+       -2.12E+02  6.61E+02
 
 TH 3
+       -1.34E+02  1.94E+02  4.35E+02
 
 TH 4
+        7.68E+01  5.29E+02 -7.19E+01  1.01E+03
 
 TH 5
+        2.56E+02 -5.40E+02 -8.03E+02 -6.72E+01  1.79E+03
 
 TH 6
+        1.39E+02  4.93E+01 -4.85E+00  5.22E+01 -4.00E+01  1.96E+02
 
 TH 7
+       -1.26E+01  3.13E+01  3.31E-02  3.90E+00  5.18E-01  2.21E-01  6.02E+00
 
 TH 8
+        2.01E+01 -2.09E+01 -2.44E+01 -3.39E+00 -3.05E+01  6.43E+00 -2.32E+00  5.28E+01
 
 TH 9
+        1.51E+02 -9.83E+01  3.72E+01  1.29E+02 -1.17E+02  3.41E+01 -1.47E+01  4.64E+01  4.34E+02
 
 TH10
+       -3.43E+01  6.58E+01 -3.64E+01  1.31E+01 -7.26E+00  8.11E+00  9.76E+00 -5.20E+00 -6.40E+01  9.25E+01
 
 TH11
+       -3.81E+01  4.57E+01  5.90E+01  9.85E+00 -1.35E+02  7.63E+00  2.20E+00  2.35E+01  5.53E+01 -2.22E+01  3.22E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       52.148
Stop Time:
Thu Sep 30 02:24:48 CDT 2021
