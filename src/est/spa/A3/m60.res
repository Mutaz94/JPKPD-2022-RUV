Sat Sep 18 10:34:55 CDT 2021
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
$DATA ../../../../data/spa/A3/dat60.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1189.78674624056        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.0801E+01  6.2351E+01  1.0312E+02 -3.8779E+01  2.5059E+02  2.1525E+01 -7.0374E+01 -7.0042E+01 -1.4812E+02 -1.9561E+02
            -5.1731E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -818.774767253344        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3191E+00  1.0846E+00  8.3334E-01  1.4610E+00  6.8067E-01  6.6041E-01  1.0382E+00  1.1187E+00  1.1735E+00  1.2045E+00
             1.5476E+01
 PARAMETER:  3.7694E-01  1.8117E-01 -8.2314E-02  4.7912E-01 -2.8468E-01 -3.1490E-01  1.3748E-01  2.1216E-01  2.5996E-01  2.8609E-01
             2.8393E+00
 GRADIENT:  -3.6223E+01 -2.2516E+01 -9.5222E+00 -2.4567E+01 -3.8892E+00  6.5023E+00  8.8680E+00  5.6988E+00  3.2988E+01  1.9190E+01
             4.3278E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1009.65467455721        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0751E+00  6.1484E-01  1.3442E-01  1.3805E+00  2.3922E-01  6.6360E-01  5.9657E-01  1.5342E+00  2.4667E-01  8.5344E-01
             9.3916E+00
 PARAMETER:  1.7240E-01 -3.8639E-01 -1.9068E+00  4.2245E-01 -1.3304E+00 -3.1008E-01 -4.1656E-01  5.2800E-01 -1.2997E+00 -5.8485E-02
             2.3398E+00
 GRADIENT:  -7.0420E+01  8.3958E+00 -9.0115E+00  7.5932E+01  3.8375E+00 -7.1507E+00  6.3299E+00  1.6894E+01  1.8322E+00  3.2852E+01
             2.8341E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1144.54991010896        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.8046E-01  8.8947E-01  1.2841E-01  8.9416E-01  3.2974E-01  9.5000E-01  3.9677E-01  9.2138E-01  4.4852E-01  4.3714E-01
             5.3680E+00
 PARAMETER:  8.0264E-02 -1.7124E-02 -1.9525E+00 -1.1872E-02 -1.0095E+00  4.8709E-02 -8.2440E-01  1.8116E-02 -7.0181E-01 -7.2749E-01
             1.7805E+00
 GRADIENT:  -6.3243E+01  8.8366E+01 -2.2893E+01  3.7455E+01 -9.6058E+01  2.4660E+01  1.2692E+00  4.9493E+00 -1.3266E+00  8.7373E+00
             1.8167E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1180.24854973868        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0056E+00  1.2679E+00  3.8072E-01  8.5013E-01  7.0143E-01  7.9887E-01  5.1079E-01  3.8522E-01  8.0442E-01  1.8320E-01
             4.9193E+00
 PARAMETER:  1.0556E-01  3.3736E-01 -8.6570E-01 -6.2361E-02 -2.5464E-01 -1.2456E-01 -5.7180E-01 -8.5393E-01 -1.1764E-01 -1.5972E+00
             1.6932E+00
 GRADIENT:  -5.0036E+01 -7.4599E+00 -1.0173E+01  7.7277E+00  8.8782E+00 -1.0469E+01  3.9404E-01  9.5119E-01 -2.2152E+00  1.1999E+00
             1.3415E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1185.81860945359        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0193E+00  1.4757E+00  8.5805E-01  7.8773E-01  1.1303E+00  8.0758E-01  3.7320E-01  1.9327E-01  1.0932E+00  1.0231E-01
             4.8525E+00
 PARAMETER:  1.1916E-01  4.8916E-01 -5.3092E-02 -1.3859E-01  2.2246E-01 -1.1371E-01 -8.8565E-01 -1.5437E+00  1.8911E-01 -2.1797E+00
             1.6795E+00
 GRADIENT:  -3.2206E+00  7.1648E-01 -2.7848E+00  2.2487E+00  5.0792E+00 -1.1301E+00  6.4619E-02  1.1409E-01 -2.2434E-02  1.5599E-01
             2.4219E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1185.94027821139        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  1.0184E+00  1.4318E+00  9.6800E-01  8.0976E-01  1.1294E+00  8.0718E-01  3.7306E-01  1.7972E-01  1.0670E+00  1.0014E-01
             4.8482E+00
 PARAMETER:  1.1827E-01  4.5893E-01  6.7472E-02 -1.1101E-01  2.2170E-01 -1.1420E-01 -8.8601E-01 -1.6164E+00  1.6483E-01 -2.2012E+00
             1.6786E+00
 GRADIENT:  -2.9201E-01  1.4839E+00  7.7646E-01 -1.5739E+00 -1.6094E+00 -2.0642E-01  7.1202E-02  8.6752E-02  1.5862E-01  1.4642E-01
             1.5494E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1186.02390854357        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  1.0155E+00  1.2750E+00  9.9817E-01  9.1868E-01  1.0674E+00  8.0650E-01  3.2607E-01  1.6932E-01  1.0232E+00  9.4852E-02
             4.8215E+00
 PARAMETER:  1.1542E-01  3.4297E-01  9.8173E-02  1.5184E-02  1.6519E-01 -1.1505E-01 -1.0206E+00 -1.6760E+00  1.2296E-01 -2.2554E+00
             1.6731E+00
 GRADIENT:  -1.1071E+01  9.9590E+00  8.6724E-01  7.8590E+00 -2.8856E+00 -1.4211E+00 -8.2602E-02  1.1097E-01  3.9457E-01  1.3637E-01
            -2.6920E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1186.26649871477        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      587
 NPARAMETR:  1.0174E+00  1.0120E+00  1.0631E+00  1.0788E+00  9.7888E-01  8.0970E-01  2.1911E-01  1.4221E-01  9.1176E-01  8.0884E-02
             4.8195E+00
 PARAMETER:  1.1721E-01  1.1189E-01  1.6121E-01  1.7583E-01  7.8656E-02 -1.1109E-01 -1.4182E+00 -1.8504E+00  7.6236E-03 -2.4147E+00
             1.6727E+00
 GRADIENT:   1.4873E-01  1.2142E+00  1.9168E-01  1.1865E+00 -9.9322E-03 -5.1418E-02 -5.6267E-02  1.1519E-01 -9.1067E-02  1.0653E-01
            -1.6324E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1186.36261259391        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  1.0160E+00  8.0052E-01  1.0400E+00  1.2045E+00  8.7906E-01  8.1046E-01  1.4454E-01  1.1791E-01  8.4534E-01  6.7859E-02
             4.7990E+00
 PARAMETER:  1.1585E-01 -1.2249E-01  1.3926E-01  2.8604E-01 -2.8906E-02 -1.1015E-01 -1.8342E+00 -2.0378E+00 -6.8021E-02 -2.5903E+00
             1.6684E+00
 GRADIENT:   2.0595E-01  1.3091E+00  9.5896E-01  1.6009E+00 -1.8542E+00 -4.8360E-02 -1.8444E-02  1.1026E-01 -2.0988E-01  8.3046E-02
            -4.3940E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1186.38976503074        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  1.0151E+00  7.3965E-01  1.0417E+00  1.2398E+00  8.6006E-01  8.1050E-01  1.2203E-01  1.0483E-01  8.2590E-01  6.0974E-02
             4.7946E+00
 PARAMETER:  1.1504E-01 -2.0158E-01  1.4089E-01  3.1494E-01 -5.0750E-02 -1.1010E-01 -2.0035E+00 -2.1554E+00 -9.1284E-02 -2.6973E+00
             1.6675E+00
 GRADIENT:   1.1448E-01 -4.6324E-01 -2.3188E-02 -5.3223E-01  4.7019E-01 -1.0544E-03 -7.9625E-03  9.3109E-02  1.5672E-02  6.8852E-02
             6.5793E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1186.42224912808        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      808
 NPARAMETR:  1.0139E+00  6.0732E-01  9.9607E-01  1.3130E+00  7.9072E-01  8.1174E-01  8.3984E-02  7.5107E-02  7.9704E-01  4.6573E-02
             4.7792E+00
 PARAMETER:  1.1384E-01 -3.9871E-01  9.6060E-02  3.7233E-01 -1.3481E-01 -1.0858E-01 -2.3771E+00 -2.4888E+00 -1.2685E-01 -2.9667E+00
             1.6643E+00
 GRADIENT:  -3.6043E-02 -2.8441E-01  6.2003E-01 -5.9480E-01 -6.5093E-01  2.0977E-01 -2.1162E-03  5.8045E-02  3.4408E-01  4.3268E-02
             2.7784E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1186.43541914092        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  1.0152E+00  5.7365E-01  9.9213E-01  1.3354E+00  7.7974E-01  8.1208E-01  7.4162E-02  6.4794E-02  7.8057E-01  4.1391E-02
             4.7917E+00
 PARAMETER:  1.1510E-01 -4.5574E-01  9.2101E-02  3.8923E-01 -1.4879E-01 -1.0816E-01 -2.5015E+00 -2.6365E+00 -1.4773E-01 -3.0847E+00
             1.6669E+00
 GRADIENT:   1.0168E+00 -1.2497E-01  2.4636E-03 -6.2927E-02 -6.8177E-03  4.9048E-02 -2.3619E-03  4.4321E-02  3.3483E-02  3.4174E-02
             6.8309E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1186.44126346460        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1119
 NPARAMETR:  1.0142E+00  5.1761E-01  9.7966E-01  1.3661E+00  7.5656E-01  8.1222E-01  6.0097E-02  4.7613E-02  7.6617E-01  3.2768E-02
             4.7836E+00
 PARAMETER:  1.1411E-01 -5.5854E-01  7.9451E-02  4.1198E-01 -1.7897E-01 -1.0798E-01 -2.7118E+00 -2.9447E+00 -1.6635E-01 -3.3183E+00
             1.6652E+00
 GRADIENT:   3.5177E-01 -2.2658E-01 -6.8464E-02 -2.2319E-01  2.1766E-01  2.2889E-02 -1.2384E-03  2.5195E-02  5.2413E-03  2.1873E-02
             1.5416E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1186.44901303266        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1298
 NPARAMETR:  1.0139E+00  4.8267E-01  9.6287E-01  1.3835E+00  7.3849E-01  8.1249E-01  5.5033E-02  2.1531E-02  7.5989E-01  1.9524E-02
             4.7811E+00
 PARAMETER:  1.1381E-01 -6.2842E-01  6.2166E-02  4.2461E-01 -2.0315E-01 -1.0765E-01 -2.7998E+00 -3.7383E+00 -1.7458E-01 -3.8361E+00
             1.6647E+00
 GRADIENT:   8.8566E-01 -6.4045E-01 -4.0575E-01 -1.4985E+00  9.2075E-01  7.3986E-02 -8.3703E-04  5.3953E-03  1.1975E-01  7.9286E-03
             6.5353E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1186.47856099684        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1475
 NPARAMETR:  1.0153E+00  6.2103E-01  9.9546E-01  1.3095E+00  7.9626E-01  8.1185E-01  1.1693E-01  1.0000E-02  7.9437E-01  1.0000E-02
             4.7905E+00
 PARAMETER:  1.1516E-01 -3.7638E-01  9.5446E-02  3.6962E-01 -1.2783E-01 -1.0844E-01 -2.0461E+00 -5.6323E+00 -1.3020E-01 -4.9044E+00
             1.6666E+00
 GRADIENT:  -6.9415E-01  8.1444E-01  3.2303E-01  2.1345E+00 -8.6404E-01 -1.0296E-01 -7.8685E-03  0.0000E+00 -9.4146E-02  0.0000E+00
            -3.7567E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1186.49468654345        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1652
 NPARAMETR:  1.0163E+00  7.2592E-01  1.0066E+00  1.2469E+00  8.3803E-01  8.1140E-01  2.9625E-01  1.0000E-02  8.1912E-01  1.0000E-02
             4.7987E+00
 PARAMETER:  1.1616E-01 -2.2032E-01  1.0660E-01  3.2069E-01 -7.6705E-02 -1.0899E-01 -1.1165E+00 -1.2002E+01 -9.9519E-02 -8.7866E+00
             1.6683E+00
 GRADIENT:   2.0702E-01 -4.7738E-01 -1.5927E-01 -1.0100E+00  4.3716E-01  2.4434E-02 -1.8698E-02  0.0000E+00  1.8334E-01  0.0000E+00
             2.3075E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1186.49794313297        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1830
 NPARAMETR:  1.0163E+00  7.2572E-01  1.0062E+00  1.2481E+00  8.3720E-01  8.1146E-01  3.3266E-01  1.0000E-02  8.1310E-01  1.0000E-02
             4.7987E+00
 PARAMETER:  1.1620E-01 -2.2059E-01  1.0616E-01  3.2160E-01 -7.7693E-02 -1.0892E-01 -1.0006E+00 -1.3109E+01 -1.0690E-01 -9.4719E+00
             1.6684E+00
 GRADIENT:  -9.9483E-02  1.6133E-01  3.9151E-02  1.2851E-01 -1.2034E-01 -3.0932E-02 -2.1535E-02  0.0000E+00 -6.3980E-02  0.0000E+00
            -1.8887E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1186.50204046202        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2005
 NPARAMETR:  1.0159E+00  6.7086E-01  9.9426E-01  1.2796E+00  8.1191E-01  8.1179E-01  3.9899E-01  1.0000E-02  7.9471E-01  1.0000E-02
             4.7945E+00
 PARAMETER:  1.1581E-01 -2.9919E-01  9.4242E-02  3.4653E-01 -1.0836E-01 -1.0852E-01 -8.1882E-01 -1.6371E+01 -1.2977E-01 -1.1519E+01
             1.6675E+00
 GRADIENT:   7.5155E-03 -1.7885E-02 -1.6565E-02 -3.1569E-02  3.6198E-02 -1.7645E-04  1.6946E-03  0.0000E+00  1.3028E-02  0.0000E+00
             2.9656E-02

0ITERATION NO.:   93    OBJECTIVE VALUE:  -1186.50205885721        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     2105
 NPARAMETR:  1.0159E+00  6.6673E-01  9.9350E-01  1.2819E+00  8.1016E-01  8.1184E-01  3.9870E-01  1.0000E-02  7.9357E-01  1.0000E-02
             4.7940E+00
 PARAMETER:  1.1577E-01 -3.0529E-01  9.3512E-02  3.4836E-01 -1.1056E-01 -1.0849E-01 -8.1971E-01 -1.6480E+01 -1.3109E-01 -1.1589E+01
             1.6674E+00
 GRADIENT:  -1.7010E-03  2.4914E-03  9.3993E-04  6.0908E-03 -2.7364E-03 -1.1688E-03 -1.1830E-05  0.0000E+00  1.3131E-03  0.0000E+00
             1.7708E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2105
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2633E-03 -6.8485E-03  7.9385E-05 -1.4408E-02 -1.4974E-05
 SE:             2.7453E-02  4.2087E-03  7.2444E-05  2.0229E-02  1.3382E-04
 N:                     100         100         100         100         100

 P VAL.:         9.3430E-01  1.0369E-01  2.7316E-01  4.7631E-01  9.1091E-01

 ETASHRINKSD(%)  8.0273E+00  8.5900E+01  9.9757E+01  3.2229E+01  9.9552E+01
 ETASHRINKVR(%)  1.5410E+01  9.8012E+01  9.9999E+01  5.4071E+01  9.9998E+01
 EBVSHRINKSD(%)  7.9226E+00  8.6370E+01  9.9671E+01  3.2278E+01  9.9452E+01
 EBVSHRINKVR(%)  1.5217E+01  9.8142E+01  9.9999E+01  5.4137E+01  9.9997E+01
 RELATIVEINF(%)  7.5558E+01  2.6304E-02  4.6983E-05  1.0797E+00  6.8007E-05
 EPSSHRINKSD(%)  1.5586E+01
 EPSSHRINKVR(%)  2.8742E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1186.5020588572133     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -451.35123229347516     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.67
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1186.502       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  6.67E-01  9.94E-01  1.28E+00  8.10E-01  8.12E-01  3.99E-01  1.00E-02  7.94E-01  1.00E-02  4.79E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.59E+03
 
 TH 2
+       -8.38E+01  3.30E+02
 
 TH 3
+        6.29E+01  1.52E+02  1.56E+02
 
 TH 4
+       -1.50E+02  3.24E+02  5.66E+01  4.21E+02
 
 TH 5
+       -5.09E+01 -3.86E+02 -3.37E+02 -2.06E+02  7.48E+02
 
 TH 6
+        2.17E+02 -5.44E+01  2.35E+00 -6.81E+01  3.04E+01  2.51E+02
 
 TH 7
+       -5.47E-01 -6.02E+00 -3.75E+00 -4.88E+00  8.84E+00  7.35E-01  1.22E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.34E+00 -2.41E+01 -1.49E+01 -1.96E+01  3.52E+01  3.86E+00  4.86E-01  0.00E+00  1.95E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.74E+01 -1.60E+01 -6.83E+00 -1.65E+01  1.76E+01  8.34E+00  3.05E-01  0.00E+00  1.26E+00  0.00E+00  1.42E+00
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.39E+03
 
 TH 2
+       -1.42E+02  2.95E+02
 
 TH 3
+        1.36E+01  9.44E+01  1.10E+02
 
 TH 4
+       -1.97E+02  3.29E+02  3.23E+01  4.67E+02
 
 TH 5
+        4.23E+01 -2.68E+02 -2.18E+02 -1.55E+02  4.97E+02
 
 TH 6
+        2.35E-01 -2.15E+01  1.42E+01 -3.66E+01 -5.70E+00  2.20E+02
 
 TH 7
+        6.74E-01 -5.36E+00 -9.60E-01 -3.47E+00  6.49E+00  4.56E-01  1.19E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.10E+00 -2.72E+01 -1.78E+00 -3.60E+00  2.55E+01  4.10E+00  4.39E+00  0.00E+00  5.52E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.18E+01 -1.42E+01 -2.91E+00 -1.25E+01  7.96E+00  5.02E+00  1.04E+00  0.00E+00  1.12E+01  0.00E+00  2.12E+01
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.39E+03
 
 TH 2
+       -2.25E+02  3.01E+02
 
 TH 3
+       -3.22E+01  5.81E+01  6.75E+01
 
 TH 4
+       -2.62E+02  3.59E+02  1.94E+01  4.95E+02
 
 TH 5
+        1.12E+02 -1.95E+02 -1.44E+02 -1.31E+02  3.39E+02
 
 TH 6
+       -2.31E+02  2.10E+01  2.24E+01 -5.39E+00 -5.14E+01  2.10E+02
 
 TH 7
+        2.24E+00 -7.49E+00 -7.98E-01 -6.44E+00  5.01E+00 -1.38E+00  6.97E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.77E+01 -6.27E+01  5.06E+00 -6.13E+01  1.23E+01 -1.06E+01  5.60E+00  0.00E+00  6.02E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.04E+02 -3.78E+01 -8.52E+00 -3.97E+01  2.64E+01  1.47E+01  1.24E+00  0.00E+00  1.23E+01  0.00E+00  6.60E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       31.587
Stop Time:
Sat Sep 18 10:35:28 CDT 2021
