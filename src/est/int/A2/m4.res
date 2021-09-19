Sat Sep 18 00:25:49 CDT 2021
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
$DATA ../../../../data/int/A2/dat4.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m4.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2622.85280187314        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2088E+02  8.0645E+01  1.7310E+02  9.2103E+00  2.1373E+02 -1.8951E+01 -1.4088E+02 -2.7322E+02 -1.2619E+01 -1.2759E+02
            -1.8746E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3155.43747309380        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.7552E-01  7.9362E-01  7.1086E-01  1.1483E+00  6.5729E-01  1.0014E+00  1.2191E+00  9.8639E-01  9.9654E-01  9.8583E-01
             1.7357E+00
 PARAMETER:  7.5217E-02 -1.3115E-01 -2.4128E-01  2.3827E-01 -3.1963E-01  1.0137E-01  2.9811E-01  8.6298E-02  9.6536E-02  8.5724E-02
             6.5141E-01
 GRADIENT:   2.8163E+01  7.6339E+01  4.2270E+01  6.6703E+01 -7.1763E+01 -1.3510E+01 -4.3054E+00  8.5531E+00 -2.5120E+00  9.5463E+00
            -1.9368E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3157.31708278544        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      166
 NPARAMETR:  9.6024E-01  7.2613E-01  5.4316E-01  1.1786E+00  5.6037E-01  1.0705E+00  1.3775E+00  6.2816E-01  1.0082E+00  7.6962E-01
             1.7924E+00
 PARAMETER:  5.9429E-02 -2.2003E-01 -5.1035E-01  2.6429E-01 -4.7916E-01  1.6817E-01  4.2026E-01 -3.6496E-01  1.0813E-01 -1.6186E-01
             6.8355E-01
 GRADIENT:  -5.2587E+00  1.2078E+02 -3.1333E+01  7.7683E+01 -3.7294E+01  1.2712E+01 -3.5605E+00  1.9158E+00 -6.7458E+00 -2.7271E+00
             4.8513E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3167.42706870108        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      263
 NPARAMETR:  9.5844E-01  4.9802E-01  4.2333E-01  1.1889E+00  4.3188E-01  1.0384E+00  1.4174E+00  5.3048E-01  1.0241E+00  6.4028E-01
             1.7154E+00
 PARAMETER:  5.7550E-02 -5.9711E-01 -7.5961E-01  2.7304E-01 -7.3961E-01  1.3769E-01  4.4886E-01 -5.3398E-01  1.2383E-01 -3.4584E-01
             6.3962E-01
 GRADIENT:  -2.0107E+01  3.7709E+00  2.1812E+01  2.8604E+01 -3.3600E+01 -2.4556E+00 -6.4985E+00 -6.5534E+00 -3.0696E+00 -9.0250E+00
            -1.6470E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3174.94763104567        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  9.7205E-01  4.8088E-01  4.1067E-01  1.1586E+00  4.2944E-01  1.0549E+00  1.5986E+00  1.2278E+00  1.0465E+00  5.3718E-01
             1.6671E+00
 PARAMETER:  7.1652E-02 -6.3214E-01 -7.8997E-01  2.4723E-01 -7.4527E-01  1.5348E-01  5.6914E-01  3.0522E-01  1.4541E-01 -5.2142E-01
             6.1111E-01
 GRADIENT:   6.3270E+00 -7.6933E+00  3.3897E+00 -7.7271E+00  4.7241E+01  1.1352E+00  7.4278E+00 -1.4583E+01 -2.5902E+00  3.6749E+00
             7.8454E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3193.80478123090        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      619
 NPARAMETR:  9.7170E-01  3.7140E-01  3.0613E-01  1.1463E+00  3.3506E-01  1.0608E+00  1.6983E+00  1.4948E+00  1.1084E+00  4.9432E-01
             1.5488E+00
 PARAMETER:  7.1294E-02 -8.9048E-01 -1.0838E+00  2.3656E-01 -9.9345E-01  1.5905E-01  6.2963E-01  5.0200E-01  2.0292E-01 -6.0457E-01
             5.3749E-01
 GRADIENT:   8.0064E+00 -2.0422E+01  1.7113E+01  4.1137E-01  2.8501E+01  2.2081E+00  9.8682E+00  3.9522E+00 -2.3704E+00  9.8451E+00
             4.4951E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3195.81298968654        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      811
 NPARAMETR:  9.6870E-01  3.7098E-01  3.0581E-01  1.1396E+00  3.2668E-01  1.0568E+00  1.6979E+00  1.4653E+00  1.1096E+00  3.3515E-01
             1.5490E+00
 PARAMETER:  6.8203E-02 -8.9161E-01 -1.0848E+00  2.3071E-01 -1.0188E+00  1.5521E-01  6.2939E-01  4.8206E-01  2.0396E-01 -9.9316E-01
             5.3763E-01
 GRADIENT:   1.9554E+00 -4.7165E+00  4.4162E+01 -1.4181E+00 -3.0998E+00  5.6711E-01  5.8872E+00 -4.5055E-02 -1.3733E+00 -1.3744E-01
             3.3829E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3197.99973454581        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      946
 NPARAMETR:  9.6483E-01  3.7144E-01  2.7010E-01  1.1356E+00  3.1766E-01  1.0530E+00  1.6280E+00  1.4621E+00  1.1180E+00  3.3734E-01
             1.4923E+00
 PARAMETER:  6.4194E-02 -8.9036E-01 -1.2090E+00  2.2716E-01 -1.0468E+00  1.5165E-01  5.8738E-01  4.7984E-01  2.1157E-01 -9.8665E-01
             5.0032E-01
 GRADIENT:   1.0365E+01  2.8031E+01 -1.2053E+01  1.3925E+01  1.2525E+02  1.4751E+00 -6.2580E+00 -1.7824E-03 -2.3217E+00  6.3056E-01
            -3.0223E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3203.82750946495        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  9.6192E-01  3.4584E-01  2.2679E-01  1.1194E+00  2.8680E-01  1.0513E+00  1.6112E+00  1.4665E+00  1.1404E+00  2.7081E-01
             1.4531E+00
 PARAMETER:  6.1171E-02 -9.6178E-01 -1.3837E+00  2.1278E-01 -1.1490E+00  1.5005E-01  5.7700E-01  4.8290E-01  2.3141E-01 -1.2063E+00
             4.7371E-01
 GRADIENT:   6.6604E+00  4.9122E+01 -4.1993E+01  3.9024E+00  1.6318E+02  1.4410E+00 -1.6257E+01  1.1753E-01 -2.6127E+00 -1.0135E+00
            -6.7973E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3204.06346892527        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1193
 NPARAMETR:  9.6209E-01  3.4491E-01  2.2619E-01  1.1197E+00  2.8595E-01  1.0514E+00  1.6132E+00  1.4653E+00  1.1414E+00  2.7353E-01
             1.4530E+00
 PARAMETER:  6.1353E-02 -9.6448E-01 -1.3864E+00  2.1307E-01 -1.1519E+00  1.5008E-01  5.7820E-01  4.8207E-01  2.3227E-01 -1.1963E+00
             4.7363E-01
 GRADIENT:  -8.6309E+00  3.3097E+01 -5.6625E+01 -4.2980E+00  8.6331E+01 -1.1174E+00 -1.8024E+01 -1.4388E+00 -4.3248E+00 -1.2187E+00
            -6.8521E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3205.60215258682        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1359
 NPARAMETR:  9.6207E-01  3.4499E-01  2.2611E-01  1.1230E+00  2.8603E-01  1.0541E+00  1.7020E+00  1.4695E+00  1.1530E+00  2.7345E-01
             1.5002E+00
 PARAMETER:  6.1328E-02 -9.6424E-01 -1.3867E+00  2.1603E-01 -1.1517E+00  1.5271E-01  6.3180E-01  4.8494E-01  2.4240E-01 -1.1966E+00
             5.0562E-01
 GRADIENT:  -9.4893E+00  3.0252E+01 -5.7138E+01  1.9339E-02  8.8403E+01  1.2371E-02 -1.5459E-02 -5.3205E-03  1.4215E-03 -8.5064E-02
             7.6722E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3205.72918199974        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1541
 NPARAMETR:  9.6204E-01  3.4495E-01  2.2638E-01  1.1288E+00  2.8563E-01  1.0574E+00  1.7003E+00  1.4574E+00  1.1393E+00  2.7335E-01
             1.4985E+00
 PARAMETER:  6.1299E-02 -9.6436E-01 -1.3855E+00  2.2111E-01 -1.1530E+00  1.5577E-01  6.3080E-01  4.7663E-01  2.3045E-01 -1.1970E+00
             5.0445E-01
 GRADIENT:  -9.5388E+00  3.1226E+01 -5.5609E+01  8.6280E+00  8.1558E+01  1.1795E+00 -4.2062E-01 -2.2436E+00 -4.5124E+00 -9.7372E-02
            -2.7402E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3205.76371299096        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1724
 NPARAMETR:  9.6201E-01  3.4503E-01  2.2634E-01  1.1228E+00  2.8566E-01  1.0541E+00  1.7026E+00  1.4679E+00  1.1537E+00  2.7325E-01
             1.4986E+00
 PARAMETER:  6.1270E-02 -9.6413E-01 -1.3857E+00  2.1580E-01 -1.1529E+00  1.5272E-01  6.3214E-01  4.8385E-01  2.4294E-01 -1.1974E+00
             5.0454E-01
 GRADIENT:  -9.5637E+00  3.1233E+01 -5.4902E+01 -2.0844E-02  8.3951E+01  1.2992E-02 -1.6009E-02 -2.7361E-02 -6.8132E-03 -3.7275E-02
            -2.0557E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3205.80102100724        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1901
 NPARAMETR:  9.6200E-01  3.4501E-01  2.2641E-01  1.1219E+00  2.8556E-01  1.0590E+00  1.6963E+00  1.4574E+00  1.1528E+00  2.7322E-01
             1.4987E+00
 PARAMETER:  6.1262E-02 -9.6418E-01 -1.3854E+00  2.1499E-01 -1.1533E+00  1.5745E-01  6.2831E-01  4.7664E-01  2.4217E-01 -1.1975E+00
             5.0458E-01
 GRADIENT:   1.0729E+06 -1.1125E+05  7.7401E+04  4.9905E+05 -9.2912E+04  1.8332E+00 -1.1510E+00  2.2508E+05  4.4305E+05  8.9597E+04
            -1.0642E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         1.2         1.9         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1901
 NO. OF SIG. DIGITS IN FINAL EST.:  1.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.2294E-03 -1.9463E-02  3.1753E-02 -2.1551E-04 -7.8472E-03
 SE:             2.9635E-02  2.7937E-02  2.7320E-02  2.9081E-02  9.8697E-03
 N:                     100         100         100         100         100

 P VAL.:         8.5993E-01  4.8600E-01  2.4512E-01  9.9409E-01  4.2656E-01

 ETASHRINKSD(%)  7.1824E-01  6.4076E+00  8.4757E+00  2.5762E+00  6.6935E+01
 ETASHRINKVR(%)  1.4313E+00  1.2405E+01  1.6233E+01  5.0860E+00  8.9067E+01
 EBVSHRINKSD(%)  5.6714E-01  4.9121E+00  8.0818E+00  3.0668E+00  6.8636E+01
 EBVSHRINKVR(%)  1.1311E+00  9.5829E+00  1.5510E+01  6.0396E+00  9.0163E+01
 RELATIVEINF(%)  9.8865E+01  5.0939E+01  3.3174E+01  8.6186E+01  2.8596E+00
 EPSSHRINKSD(%)  2.3079E+01
 EPSSHRINKVR(%)  4.0831E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3205.8010210072357     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1551.7116612388249     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    49.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3205.801       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.62E-01  3.45E-01  2.26E-01  1.12E+00  2.86E-01  1.06E+00  1.70E+00  1.46E+00  1.15E+00  2.73E-01  1.50E+00
 


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
+        2.90E+09
 
 TH 2
+       -1.94E+04  2.42E+08
 
 TH 3
+        1.72E+03 -5.78E+03  2.73E+08
 
 TH 4
+       -2.88E+04  8.34E+03 -9.00E+03  4.61E+08
 
 TH 5
+       -1.67E+03  2.33E+03  2.60E+08  7.51E+03  2.47E+08
 
 TH 6
+        7.83E+03 -2.26E+03  2.39E+03  3.12E+03 -2.32E+03  1.73E+02
 
 TH 7
+        3.13E+03 -9.31E+02  9.62E+02  1.25E+03 -8.84E+02 -9.23E-03  5.27E+01
 
 TH 8
+       -5.65E+04  1.63E+04 -1.74E+04  1.60E+08  1.50E+04  1.09E+03  4.36E+02  5.56E+07
 
 TH 9
+        2.16E+04 -6.26E+03  6.59E+03  8.62E+03 -6.47E+03 -5.76E+08 -9.02E+07  3.04E+03  3.44E+08
 
 TH10
+       -6.28E+04  2.46E+08 -1.99E+04 -8.31E+03  1.40E+04  2.32E+03  9.38E+02 -1.61E+04  6.49E+03  2.51E+08
 
 TH11
+        3.68E+08 -1.07E+08 -2.08E+05  3.56E+03  2.01E+05 -1.00E+03 -3.95E+02  6.89E+03 -2.82E+03 -1.09E+08  4.70E+07
 
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
 #CPUT: Total CPU Time in Seconds,       63.774
Stop Time:
Sat Sep 18 00:26:54 CDT 2021
