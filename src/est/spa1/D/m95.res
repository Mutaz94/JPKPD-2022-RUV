Thu Sep 30 03:48:14 CDT 2021
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
$DATA ../../../../data/spa1/D/dat95.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   30180.0762664309        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.8889E+02  7.8888E+02  6.0640E+01  4.8522E+02 -6.2217E+01 -3.0506E+03 -1.5360E+03 -1.3677E+02 -2.3656E+03 -6.9978E+02
            -5.6028E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -373.316276823039        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1873E+00  9.5690E-01  8.2232E-01  1.7558E+00  1.7309E+00  2.9550E+00  1.2352E+00  9.1005E-01  1.4553E+00  8.6776E-01
             1.3660E+01
 PARAMETER:  2.7166E-01  5.5945E-02 -9.5622E-02  6.6295E-01  6.4864E-01  1.1835E+00  3.1121E-01  5.7471E-03  4.7525E-01 -4.1841E-02
             2.7144E+00
 GRADIENT:  -6.3904E+00  5.7572E+01  4.4833E+00  6.1408E+01 -8.6630E+00  8.2983E+01 -8.3352E+00 -1.0296E+00 -9.1888E+01  3.6770E-01
            -8.4208E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -444.229651711374        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.1759E+00  5.3256E-01  5.4704E-01  2.1959E+00  8.1013E+00  2.5019E+00  4.7363E-01  6.0466E-01  2.7086E+00  4.4364E-02
             1.4762E+01
 PARAMETER:  2.6201E-01 -5.3005E-01 -5.0324E-01  8.8658E-01  2.1920E+00  1.0171E+00 -6.4732E-01 -4.0309E-01  1.0964E+00 -3.0153E+00
             2.7921E+00
 GRADIENT:  -1.1084E+01  2.7087E+01 -1.8810E+01  7.9579E+01 -7.2216E+00  8.6855E+00  1.3053E+00  2.0935E+00  4.6873E+01  7.7036E-04
             6.4742E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -527.592432360291        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.1628E-01  8.6308E-02  3.7327E-01  1.3794E+00  9.8436E+00  1.4509E+00  9.1452E-02  3.5183E+00  6.8710E-01  5.3634E+00
             1.4310E+01
 PARAMETER:  1.2568E-02 -2.3498E+00 -8.8546E-01  4.2163E-01  2.3868E+00  4.7221E-01 -2.2919E+00  1.3580E+00 -2.7527E-01  1.7796E+00
             2.7610E+00
 GRADIENT:  -6.7378E+01  1.8922E+01  2.0800E+01  9.0781E+01 -6.6269E+00 -3.6702E+01  1.3486E-02 -8.2431E+00  3.4846E+00  1.8314E+00
             4.9042E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -547.856918804621        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  7.9008E-01  2.5794E-02  1.1635E-01  9.0901E-01  4.0399E+01  1.0937E+00  1.3601E-02  2.8413E+00  2.9407E-01  2.9253E+00
             1.5447E+01
 PARAMETER: -1.3563E-01 -3.5576E+00 -2.0512E+00  4.6032E-03  3.7988E+00  1.8955E-01 -4.1976E+00  1.1443E+00 -1.1239E+00  1.1734E+00
             2.8374E+00
 GRADIENT:   1.2373E+01 -3.8861E-01  1.2099E+01  2.9248E+01  7.2755E-02 -3.7907E+01  2.0272E-06 -2.2552E+00  3.1790E+00 -5.7821E-04
             9.7647E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -557.087928177453        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      396
 NPARAMETR:  5.2537E-01  1.3934E-02  3.5361E-02  4.9102E-01  2.7263E+02  1.6992E+00  1.3400E-02  1.6254E+00  1.8024E-01  2.1330E+00
             1.2878E+01
 PARAMETER: -5.4365E-01 -4.1735E+00 -3.2421E+00 -6.1126E-01  5.7081E+00  6.3014E-01 -4.2125E+00  5.8577E-01 -1.6135E+00  8.5754E-01
             2.6555E+00
 GRADIENT:  -4.6481E+01  3.9810E-01 -1.5202E+02  2.3818E+02  1.0607E-02 -2.3897E+01  2.2538E-07 -3.8047E+00  4.1565E-01 -3.1179E-06
            -2.0339E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -587.250613175481        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      572
 NPARAMETR:  6.3923E-01  3.2029E-02  6.0595E-02  5.6070E-01  1.4948E+02  1.8905E+00  9.0043E-02  1.2870E+00  2.5029E-01  1.4191E+00
             1.2357E+01
 PARAMETER: -3.4749E-01 -3.3411E+00 -2.7035E+00 -4.7857E-01  5.1071E+00  7.3687E-01 -2.3075E+00  3.5228E-01 -1.2851E+00  4.5000E-01
             2.6142E+00
 GRADIENT:   5.8677E+00 -1.6887E+00  3.0795E+00  2.5347E+01  2.9916E-02 -5.5157E+00  1.6112E-03 -8.2154E+00  1.5810E+00  7.8727E-06
            -3.6974E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -599.983491682071        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      748
 NPARAMETR:  3.9154E-01  1.1688E-02  1.6478E-02  2.1774E-01  5.8294E+02  1.8500E+00  2.0777E-02  6.2040E-01  6.7079E-02  3.8414E-01
             1.2874E+01
 PARAMETER: -8.3766E-01 -4.3492E+00 -4.0057E+00 -1.4245E+00  6.4681E+00  7.1516E-01 -3.7739E+00 -3.7738E-01 -2.6019E+00 -8.5674E-01
             2.6552E+00
 GRADIENT:  -3.1634E+00  1.1768E+00  5.8370E-01  1.3356E+00 -5.5802E-05  7.7045E-01  9.0887E-05 -1.6008E+00  1.2974E-01  6.5718E-08
             2.0033E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -600.864845028166        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      893
 NPARAMETR:  3.5536E-01  1.0000E-02  1.2932E-02  1.8033E-01  1.0154E+04  1.8622E+00  1.0000E-02  5.7427E-01  1.3995E-02  3.0526E-01
             1.2724E+01
 PARAMETER: -9.3462E-01 -4.6837E+00 -4.2481E+00 -1.6130E+00  9.3256E+00  7.2175E-01 -6.3963E+00 -4.5465E-01 -4.1691E+00 -1.0866E+00
             2.6435E+00
 GRADIENT:   1.0036E+00  0.0000E+00 -1.8670E+01  1.7992E+01 -2.4680E-05  5.8334E+00  0.0000E+00  7.1581E-02  6.0770E-03  2.1449E-10
            -8.6750E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -601.015804525002        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  3.5820E-01  1.0000E-02  1.3117E-02  1.8069E-01  4.7145E+04  1.8277E+00  1.0000E-02  5.8927E-01  1.0000E-02  3.0439E-01
             1.2834E+01
 PARAMETER: -9.2666E-01 -4.6837E+00 -4.2339E+00 -1.6110E+00  1.0861E+01  7.0303E-01 -5.9494E+00 -4.2887E-01 -5.0049E+00 -1.0894E+00
             2.6521E+00
 GRADIENT:   6.9067E-01  0.0000E+00 -3.6246E+00  4.2398E-01 -9.3887E-06  1.9166E-01  0.0000E+00  1.9248E-02  0.0000E+00  1.0435E-11
            -1.7987E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -601.023154017196        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     1280             RESET HESSIAN, TYPE I
 NPARAMETR:  3.5796E-01  1.0000E-02  1.3068E-02  1.8025E-01  1.1465E+06  1.8276E+00  1.0000E-02  5.8799E-01  1.0000E-02  3.1255E-01
             1.2838E+01
 PARAMETER: -9.2733E-01 -4.6837E+00 -4.2376E+00 -1.6134E+00  1.4052E+01  7.0301E-01 -5.9494E+00 -4.3104E-01 -5.0049E+00 -1.0630E+00
             2.6524E+00
 GRADIENT:   6.2667E+01  0.0000E+00  8.8585E+01  2.7473E+01 -3.5811E-07  1.2869E+01  0.0000E+00  3.4583E-01  0.0000E+00  0.0000E+00
             2.6363E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -601.027050856667        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1479
 NPARAMETR:  3.5751E-01  1.0000E-02  1.3025E-02  1.7992E-01  2.3279E+06  1.8275E+00  1.0000E-02  5.9071E-01  1.0000E-02  3.1344E-01
             1.2840E+01
 PARAMETER: -9.2860E-01 -4.6837E+00 -4.2409E+00 -1.6153E+00  1.4760E+01  7.0294E-01 -5.9494E+00 -4.2643E-01 -5.0049E+00 -1.0601E+00
             2.6526E+00
 GRADIENT:   1.0406E+00  0.0000E+00 -6.2717E+00  3.1858E+00 -1.9469E-07  2.6487E-01  0.0000E+00  2.2758E-01  0.0000E+00  0.0000E+00
             4.8076E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -601.032273064529        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1664
 NPARAMETR:  3.5705E-01  1.0000E-02  1.2994E-02  1.7952E-01  2.7159E+06  1.8274E+00  1.0000E-02  5.9137E-01  1.0000E-02  3.1277E-01
             1.2843E+01
 PARAMETER: -9.2988E-01 -4.6837E+00 -4.2433E+00 -1.6175E+00  1.4915E+01  7.0289E-01 -5.9494E+00 -4.2531E-01 -5.0049E+00 -1.0623E+00
             2.6528E+00
 GRADIENT:   8.0699E-01  0.0000E+00 -5.9476E+00  2.8155E+00 -1.8342E-07  3.1418E-01  0.0000E+00  2.7513E-01  0.0000E+00  5.3511E-12
             4.0461E-01

0ITERATION NO.:   64    OBJECTIVE VALUE:  -601.035578993761        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1804
 NPARAMETR:  3.5718E-01  1.0000E-02  1.2986E-02  1.7942E-01  3.3048E+06  1.8269E+00  1.0000E-02  5.8432E-01  1.0000E-02  3.0732E-01
             1.2833E+01
 PARAMETER: -9.2953E-01 -4.6837E+00 -4.2439E+00 -1.6180E+00  1.5111E+01  7.0263E-01 -5.9494E+00 -4.3730E-01 -5.0049E+00 -1.0799E+00
             2.6520E+00
 GRADIENT:   1.5151E+00  0.0000E+00 -5.5966E+00  2.2034E+00 -1.4368E-07  1.5190E-01  0.0000E+00  3.1559E-02  0.0000E+00  0.0000E+00
            -8.3969E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1804
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3019E-03 -3.6166E-06  9.1730E-03 -3.0967E-04 -3.7051E-12
 SE:             2.9039E-02  5.2348E-06  1.5763E-02  3.0060E-04  9.0068E-11
 N:                     100         100         100         100         100

 P VAL.:         9.0947E-01  4.8964E-01  5.6061E-01  3.0293E-01  9.6719E-01

 ETASHRINKSD(%)  2.7145E+00  9.9982E+01  4.7192E+01  9.8993E+01  1.0000E+02
 ETASHRINKVR(%)  5.3553E+00  1.0000E+02  7.2113E+01  9.9990E+01  1.0000E+02
 EBVSHRINKSD(%)  2.9748E+00  9.9958E+01  4.8361E+01  9.8953E+01  1.0000E+02
 EBVSHRINKVR(%)  5.8612E+00  1.0000E+02  7.3334E+01  9.9989E+01  1.0000E+02
 RELATIVEINF(%)  7.1952E+00  9.8116E-06  1.4261E-01  5.8769E-05  0.0000E+00
 EPSSHRINKSD(%)  6.0876E+00
 EPSSHRINKVR(%)  1.1805E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -601.03557899376108     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       317.90295421091162     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -601.036       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.57E-01  1.00E-02  1.30E-02  1.79E-01  3.30E+06  1.83E+00  1.00E-02  5.84E-01  1.00E-02  3.07E-01  1.28E+01
 


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
+        2.42E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.30E+04  0.00E+00  4.08E+06
 
 TH 4
+       -3.39E+02  0.00E+00 -3.30E+05  2.93E+04
 
 TH 5
+       -4.79E-12  0.00E+00 -9.64E-11  2.10E-11  7.03E-20
 
 TH 6
+        2.78E+00  0.00E+00  3.21E+02 -6.41E+01 -3.34E-13  5.06E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -3.71E+01  0.00E+00 -3.38E+03  2.67E+02  8.31E-12  2.18E+00  0.00E+00  2.47E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        9.56E-04  0.00E+00 -2.04E-03 -5.97E-04 -9.17E-13  1.95E-05  0.00E+00  1.50E-04  0.00E+00  1.78E-03
 
 TH11
+       -1.99E+01  0.00E+00  5.88E+02 -3.62E+01  6.61E-14  8.62E-01  0.00E+00  2.55E+00  0.00E+00 -4.01E-06  2.54E+00
 
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
 #CPUT: Total CPU Time in Seconds,       46.475
Stop Time:
Thu Sep 30 03:49:02 CDT 2021
