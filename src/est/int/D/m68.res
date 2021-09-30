Wed Sep 29 09:30:59 CDT 2021
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
$DATA ../../../../data/int/D/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26421.3210914178        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5860E+02  5.0206E+02  1.7011E+01  2.2551E+02  4.4524E+02 -2.8353E+03 -1.1781E+03 -1.0556E+02 -2.1676E+03 -1.0739E+03
            -5.2809E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1000.35188868175        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  2.2538E+00  1.5408E+00  7.9700E-01  2.6936E+00  5.9958E-01  5.2235E+00  4.6165E+00  1.0087E+00  4.9109E+00  3.3660E+00
             1.0008E+01
 PARAMETER:  9.1263E-01  5.3230E-01 -1.2690E-01  1.0909E+00 -4.1153E-01  1.7532E+00  1.6296E+00  1.0867E-01  1.6915E+00  1.3137E+00
             2.4034E+00
 GRADIENT:   8.8745E+01  4.2822E+01 -2.2786E+01  7.5874E+01 -2.1420E+01  1.9760E+02  7.9135E+01  5.3591E+00  1.0652E+02  5.6012E+01
             3.3115E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1071.79839113609        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      266
 NPARAMETR:  2.6850E+00  5.2075E-01  1.0194E+01  5.3040E+00  1.8025E+00  3.9360E+00  4.1748E+00  2.7264E-01  4.3948E+00  2.8851E+00
             1.0340E+01
 PARAMETER:  1.0877E+00 -5.5248E-01  2.4218E+00  1.7685E+00  6.8916E-01  1.4702E+00  1.5291E+00 -1.1996E+00  1.5804E+00  1.1596E+00
             2.4360E+00
 GRADIENT:   8.1678E+01  7.1369E+00 -4.1499E+01  9.5063E+01  1.7175E+01  5.7428E+00  5.4845E+00 -3.5080E-03 -1.1300E+01  5.6798E+01
             3.3965E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1249.03452937116        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.1961E+00  3.8013E-01  1.8831E+02  2.4055E+00  2.4770E+00  3.3637E+00  7.3942E-01  1.0761E-02  4.3312E+00  1.4736E+00
             9.4865E+00
 PARAMETER:  2.7910E-01 -8.6725E-01  5.3381E+00  9.7777E-01  1.0071E+00  1.3130E+00 -2.0189E-01 -4.4319E+00  1.5658E+00  4.8769E-01
             2.3499E+00
 GRADIENT:  -3.7853E+01 -1.6155E+01 -1.7982E+00  6.1704E+00  1.9258E+01 -1.8087E+00  6.5286E-01  1.6968E-05  4.9161E+01  3.2169E+01
             3.4660E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1302.10738028706        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      618
 NPARAMETR:  1.3459E+00  1.2697E-01  5.9328E+02  2.3507E+00  2.3890E+00  3.3556E+00  9.3362E-01  1.0000E-02  3.5138E+00  6.2698E-01
             7.7176E+00
 PARAMETER:  3.9707E-01 -1.9638E+00  6.4857E+00  9.5473E-01  9.7090E-01  1.3106E+00  3.1311E-02 -5.5983E+00  1.3567E+00 -3.6684E-01
             2.1435E+00
 GRADIENT:  -2.2072E+00 -3.4957E+00 -2.0857E-01  8.2582E-01  7.1365E+00 -1.0814E+00  1.2665E-01  0.0000E+00  5.9545E-02  3.2021E+00
             1.5257E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1320.67377692309        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  1.2285E+00  1.2429E+00  2.3088E+02  1.4153E+00  2.3281E+00  3.3799E+00  1.4978E+00  1.0000E-02  5.0198E+00  4.3898E-01
             7.6329E+00
 PARAMETER:  3.0581E-01  3.1748E-01  5.5419E+00  4.4734E-01  9.4507E-01  1.3178E+00  5.0397E-01 -6.5504E+00  1.7134E+00 -7.2329E-01
             2.1325E+00
 GRADIENT:  -1.8818E+01  9.4197E+00 -1.3682E+00  1.5692E+01 -1.1834E+01  1.5837E+00  3.4492E+00  0.0000E+00 -1.8921E+01  5.4832E-02
            -3.5025E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1326.91124210897        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  1.3898E+00  1.4576E+00  2.9436E+02  9.6635E-01  2.3945E+00  3.3240E+00  1.0787E+00  1.0000E-02  6.0930E+00  4.3371E-01
             7.7893E+00
 PARAMETER:  4.2918E-01  4.7678E-01  5.7848E+00  6.5770E-02  9.7319E-01  1.3012E+00  1.7575E-01 -6.1924E+00  1.9071E+00 -7.3538E-01
             2.1527E+00
 GRADIENT:   3.8627E+00  6.6409E+00 -9.0782E-01  5.1730E+00  5.9738E-02 -4.3813E+00 -2.8590E+00  0.0000E+00 -6.9254E+00  2.2207E-01
            -1.4399E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1327.78187190290        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1152
 NPARAMETR:  1.3706E+00  1.7045E+00  5.4538E+02  7.1597E-01  2.4246E+00  3.3560E+00  1.3548E+00  1.0000E-02  6.9831E+00  2.6526E-01
             7.8030E+00
 PARAMETER:  4.1522E-01  6.3325E-01  6.4015E+00 -2.3411E-01  9.8567E-01  1.3108E+00  4.0363E-01 -8.7201E+00  2.0435E+00 -1.2271E+00
             2.1545E+00
 GRADIENT:   1.2794E+00  5.4499E+00 -1.6796E-01  2.9429E+00 -4.1783E-01  1.9757E-01 -6.5772E-01  0.0000E+00 -2.4005E+00 -1.7274E-02
            -3.6613E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1328.13453388133        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1341
 NPARAMETR:  1.3654E+00  1.6811E+00  8.0442E+02  6.8320E-01  2.4243E+00  3.4318E+00  1.3480E+00  1.0000E-02  6.8478E+00  2.6252E-01
             7.8032E+00
 PARAMETER:  4.1148E-01  6.1947E-01  6.7901E+00 -2.8097E-01  9.8554E-01  1.3331E+00  3.9859E-01 -8.7201E+00  2.0239E+00 -1.2374E+00
             2.1545E+00
 GRADIENT:   6.0304E-01  1.4847E+00 -1.4195E-02 -1.2950E-02 -9.7941E-01  9.1834E+00 -3.0461E-02  0.0000E+00 -6.8006E+00 -2.1389E-04
            -3.5729E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1328.15748433706        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1528
 NPARAMETR:  1.3648E+00  1.6682E+00  8.8607E+02  7.0135E-01  2.4283E+00  3.4389E+00  1.3336E+00  1.0000E-02  6.7999E+00  2.4631E-01
             7.8010E+00
 PARAMETER:  4.1103E-01  6.1173E-01  6.8868E+00 -2.5475E-01  9.8719E-01  1.3352E+00  3.8789E-01 -8.7201E+00  2.0169E+00 -1.3012E+00
             2.1542E+00
 GRADIENT:   5.1086E-01  2.0628E+00 -7.6203E-04  4.0482E-01 -7.9197E-01  1.0011E+01 -2.8162E-01  0.0000E+00 -6.5224E+00 -1.0186E-02
            -1.8040E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1328.17159791576        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1715
 NPARAMETR:  1.3649E+00  1.6590E+00  9.0208E+02  7.0637E-01  2.4291E+00  3.4391E+00  1.3253E+00  1.0000E-02  6.7655E+00  2.4196E-01
             7.8026E+00
 PARAMETER:  4.1106E-01  6.0620E-01  6.9047E+00 -2.4762E-01  9.8753E-01  1.3352E+00  3.8167E-01 -8.7201E+00  2.0118E+00 -1.3190E+00
             2.1545E+00
 GRADIENT:   5.3100E-01  1.7108E+00  5.7109E-05  2.4592E-01 -6.4570E-01  1.0031E+01 -2.8878E-01  0.0000E+00 -6.7759E+00 -1.0801E-02
            -1.5287E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1328.18343136971        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1902
 NPARAMETR:  1.3652E+00  1.6516E+00  9.0995E+02  7.0902E-01  2.4295E+00  3.4448E+00  1.3176E+00  1.0000E-02  6.7288E+00  2.4547E-01
             7.8037E+00
 PARAMETER:  4.1127E-01  6.0172E-01  6.9134E+00 -2.4387E-01  9.8769E-01  1.3369E+00  3.7585E-01 -8.7201E+00  2.0064E+00 -1.3046E+00
             2.1546E+00
 GRADIENT:   5.0813E-01  1.5683E+00  2.1629E-05  8.9862E-03 -4.6982E-01  1.0685E+01 -3.5126E-01  0.0000E+00 -7.3853E+00 -7.3213E-03
            -1.2292E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1328.18552765033        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     2097
 NPARAMETR:  1.3654E+00  1.6412E+00  9.0722E+02  7.0908E-01  2.4304E+00  3.4500E+00  1.3315E+00  1.0000E-02  6.7134E+00  2.4961E-01
             7.8066E+00
 PARAMETER:  4.1148E-01  5.9545E-01  6.9104E+00 -2.4378E-01  9.8805E-01  1.3384E+00  3.8631E-01 -8.7201E+00  2.0041E+00 -1.2878E+00
             2.1550E+00
 GRADIENT:   5.6602E-01 -1.2695E+00 -2.0156E-03 -4.7457E-01  4.8012E-01  1.1224E+01  6.8818E-01  0.0000E+00 -6.7264E+00  1.3313E-02
             1.2817E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1328.19680532277        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2289             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3654E+00  1.6388E+00  9.1659E+02  7.2020E-01  2.4299E+00  3.4501E+00  1.3145E+00  1.0000E-02  6.6810E+00  2.3368E-01
             7.8048E+00
 PARAMETER:  4.1148E-01  5.9394E-01  6.9207E+00 -2.2823E-01  9.8785E-01  1.3384E+00  3.7347E-01 -8.7201E+00  1.9993E+00 -1.3538E+00
             2.1547E+00
 GRADIENT:   3.8571E+01  2.4368E+01  6.9001E-03  4.2545E+00  4.2235E+00  1.2830E+02  1.2853E+00  0.0000E+00  1.9901E+02  2.0960E-02
             3.7934E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1328.20057407652        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2477
 NPARAMETR:  1.3655E+00  1.6340E+00  9.2607E+02  7.2526E-01  2.4304E+00  3.4501E+00  1.3106E+00  1.0000E-02  6.6623E+00  2.3792E-01
             7.8055E+00
 PARAMETER:  4.1150E-01  5.9106E-01  6.9309E+00 -2.2122E-01  9.8807E-01  1.3384E+00  3.7048E-01 -8.7201E+00  1.9965E+00 -1.3358E+00
             2.1548E+00
 GRADIENT:   5.4880E-01  8.9320E-01 -2.3632E-04  2.9685E-02 -2.6334E-01  1.1263E+01 -1.0899E-01  0.0000E+00 -7.3724E+00 -6.7292E-03
            -7.2259E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1328.20339315530        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2657
 NPARAMETR:  1.3655E+00  1.6219E+00  9.2361E+02  7.2876E-01  2.4310E+00  3.4501E+00  1.3119E+00  1.0000E-02  6.6376E+00  2.4552E-01
             7.8081E+00
 PARAMETER:  4.1155E-01  5.8843E-01  6.9344E+00 -2.1715E-01  9.8795E-01  1.3384E+00  3.6777E-01 -8.7201E+00  1.9941E+00 -1.2932E+00
             2.1548E+00
 GRADIENT:   1.0496E-03  8.5169E-01  5.5276E-04 -1.5256E-02 -9.9996E-02 -1.1238E-04 -4.6228E-02  0.0000E+00  1.6430E-01  1.2936E-03
            -3.4503E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2657
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0705E-03 -7.6834E-02  1.3667E-06  3.5679E-02 -6.1169E-04
 SE:             2.9599E-02  1.4171E-02  1.5026E-06  2.4377E-02  4.0385E-03
 N:                     100         100         100         100         100

 P VAL.:         8.3750E-01  5.9081E-08  3.6306E-01  1.4329E-01  8.7961E-01

 ETASHRINKSD(%)  8.3910E-01  5.2526E+01  9.9995E+01  1.8335E+01  8.6470E+01
 ETASHRINKVR(%)  1.6712E+00  7.7462E+01  1.0000E+02  3.3308E+01  9.8170E+01
 EBVSHRINKSD(%)  8.6193E-01  5.7378E+01  9.9991E+01  1.1910E+01  8.6623E+01
 EBVSHRINKVR(%)  1.7164E+00  8.1834E+01  1.0000E+02  2.2401E+01  9.8211E+01
 RELATIVEINF(%)  9.8206E+01  9.5065E+00  7.0299E-07  4.0896E+01  1.5351E+00
 EPSSHRINKSD(%)  7.7091E+00
 EPSSHRINKVR(%)  1.4824E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1328.2033931552999     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       325.88596661311090     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    99.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1328.203       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.37E+00  1.63E+00  9.29E+02  7.28E-01  2.43E+00  3.45E+00  1.31E+00  1.00E-02  6.65E+00  2.48E-01  7.81E+00
 


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
+        4.86E+01
 
 TH 2
+       -5.52E-01  6.43E+01
 
 TH 3
+        1.87E-05  1.99E-05  1.02E-07
 
 TH 4
+       -7.63E-02  2.26E+01 -9.51E-06  3.78E+01
 
 TH 5
+       -9.18E-01 -9.64E+00 -7.04E-04 -4.29E+00  4.93E+01
 
 TH 6
+        9.55E-02  3.07E-01  1.72E-05  6.25E-02 -1.02E-01  1.54E+01
 
 TH 7
+        1.66E-01 -1.99E+01  8.46E-06 -2.68E+00  3.43E+00 -1.45E-01  1.41E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.05E-01 -6.27E+00  1.70E-05  3.43E+00  5.91E-01 -4.24E-03  1.71E+00  0.00E+00  2.58E+00
 
 TH10
+       -3.32E-01 -1.67E+00  1.72E-05 -6.39E-01  3.56E+00  1.74E-02  8.46E-01  0.00E+00  1.29E-01  1.80E+00
 
 TH11
+       -2.39E+00 -5.72E+00  5.71E-06 -2.57E+00 -6.33E-02  6.77E-01  3.81E+00  0.00E+00  3.65E-01  1.05E+00  1.67E+01
 
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
 #CPUT: Total CPU Time in Seconds,      117.510
Stop Time:
Wed Sep 29 09:32:58 CDT 2021
