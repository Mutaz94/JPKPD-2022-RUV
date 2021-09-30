Thu Sep 30 01:17:45 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2043.68442499146        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7325E+02 -2.3718E+01 -4.6767E+01  3.5502E+01  6.6856E+01  2.3317E+01 -1.8265E+01  7.5825E+00 -1.4470E+01  1.4504E+01
            -9.1211E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2056.24172708366        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9071E-01  1.0531E+00  1.0828E+00  1.0259E+00  1.0168E+00  1.0732E+00  1.0801E+00  9.7249E-01  1.0785E+00  9.2868E-01
             1.1046E+00
 PARAMETER:  9.0667E-02  1.5178E-01  1.7953E-01  1.2558E-01  1.1662E-01  1.7069E-01  1.7702E-01  7.2100E-02  1.7557E-01  2.6005E-02
             1.9951E-01
 GRADIENT:   7.3168E+00  1.1279E+00 -1.8006E+01  2.4854E+01  2.0695E+01  2.0853E+00 -9.8810E+00  2.6158E+00 -1.8102E+00  3.5637E+00
            -5.5725E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2057.44125502907        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.8868E-01  9.7435E-01  1.1049E+00  1.0663E+00  9.7653E-01  1.0555E+00  1.3460E+00  9.2269E-01  1.0103E+00  8.7792E-01
             1.0938E+00
 PARAMETER:  8.8615E-02  7.4011E-02  1.9976E-01  1.6421E-01  7.6249E-02  1.5401E-01  3.9716E-01  1.9541E-02  1.1021E-01 -3.0196E-02
             1.8962E-01
 GRADIENT:   5.1228E+00  7.2620E+00 -1.4286E+00  7.4885E+00 -2.1732E+00 -4.4550E+00  2.9696E+00  9.7487E-01  4.4857E-01  1.1947E+00
            -1.4030E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2057.87407621403        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8362E-01  8.1258E-01  1.2854E+00  1.1741E+00  9.8121E-01  1.0692E+00  1.3850E+00  9.8355E-01  9.8579E-01  9.0141E-01
             1.1232E+00
 PARAMETER:  8.3485E-02 -1.0755E-01  3.5111E-01  2.6053E-01  8.1028E-02  1.6690E-01  4.2567E-01  8.3417E-02  8.5691E-02 -3.7916E-03
             2.1622E-01
 GRADIENT:  -1.6436E+00  5.7181E+00  3.5915E+00  5.3529E+00 -5.3200E+00  1.5263E+00 -4.4652E-01 -7.8891E-01  4.7246E-01 -8.6463E-01
             6.1020E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2058.46025963287        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.7995E-01  4.1892E-01  1.6236E+00  1.4421E+00  9.6494E-01  1.0561E+00  1.7724E+00  1.1982E+00  8.9080E-01  9.4913E-01
             1.1058E+00
 PARAMETER:  7.9744E-02 -7.7008E-01  5.8466E-01  4.6613E-01  6.4314E-02  1.5462E-01  6.7235E-01  2.8085E-01 -1.5634E-02  4.7787E-02
             2.0059E-01
 GRADIENT:   9.6217E-01  6.5503E+00  2.8760E+00  2.0215E+01 -6.1290E+00 -1.7344E+00 -5.2251E-01 -1.0014E+00 -9.0019E-01  3.0732E-01
            -6.4230E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2058.76077736039        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.7737E-01  2.8133E-01  1.7266E+00  1.5296E+00  9.5522E-01  1.0569E+00  2.0696E+00  1.2822E+00  8.5933E-01  9.4958E-01
             1.1110E+00
 PARAMETER:  7.7108E-02 -1.1682E+00  6.4617E-01  5.2503E-01  5.4189E-02  1.5538E-01  8.2736E-01  3.4855E-01 -5.1599E-02  4.8263E-02
             2.0522E-01
 GRADIENT:  -7.2476E-01  4.1346E+00  2.0569E+00  1.6238E+01 -6.0212E+00 -9.1393E-01 -5.3246E-01 -5.0726E-01 -1.0747E+00  5.9025E-01
            -2.3910E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2058.83814905315        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.7609E-01  1.8691E-01  1.8394E+00  1.5953E+00  9.5771E-01  1.0590E+00  2.4586E+00  1.3747E+00  8.3762E-01  9.5207E-01
             1.1153E+00
 PARAMETER:  7.5803E-02 -1.5771E+00  7.0942E-01  5.6707E-01  5.6792E-02  1.5737E-01  9.9960E-01  4.1823E-01 -7.7188E-02  5.0884E-02
             2.0909E-01
 GRADIENT:  -1.0717E+00  3.3138E+00  2.5117E+00  1.9048E+01 -7.4097E+00  1.2042E-01 -3.8091E-01  3.0346E-02 -1.1685E+00  7.7482E-01
             8.2054E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2058.98882997216        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.7525E-01  1.0585E-01  1.9319E+00  1.6470E+00  9.6032E-01  1.0600E+00  3.1954E+00  1.4450E+00  8.2007E-01  9.5175E-01
             1.1174E+00
 PARAMETER:  7.4934E-02 -2.1458E+00  7.5851E-01  5.9893E-01  5.9511E-02  1.5828E-01  1.2617E+00  4.6814E-01 -9.8361E-02  5.0546E-02
             2.1098E-01
 GRADIENT:  -8.1220E-01  1.7455E+00  1.9542E+00  1.2506E+01 -5.2328E+00  6.8639E-01 -1.9203E-01  2.7774E-01 -7.2221E-01  4.6916E-01
             2.3724E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2059.24190954427        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1414             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7531E-01  3.6080E-02  1.9769E+00  1.6710E+00  9.6034E-01  1.0589E+00  5.8691E+00  1.4748E+00  8.0859E-01  9.4943E-01
             1.1147E+00
 PARAMETER:  7.4998E-02 -3.2220E+00  7.8151E-01  6.1340E-01  5.9530E-02  1.5718E-01  1.8697E+00  4.8851E-01 -1.1247E-01  4.8108E-02
             2.0857E-01
 GRADIENT:   3.5575E+02  5.3858E+00  8.2254E+00  1.0968E+03  1.0584E+01  6.6393E+01  1.9739E+00  1.5127E+00  1.8937E+01  7.7929E-02
             2.8177E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2059.27506578783        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1596             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7520E-01  2.2879E-02  1.9772E+00  1.6791E+00  9.5560E-01  1.0583E+00  6.9605E+00  1.4696E+00  8.0335E-01  9.5022E-01
             1.1138E+00
 PARAMETER:  7.4884E-02 -3.6775E+00  7.8168E-01  6.1825E-01  5.4583E-02  1.5663E-01  2.0403E+00  4.8500E-01 -1.1897E-01  4.8934E-02
             2.0781E-01
 GRADIENT:   3.5635E+02  3.3197E+00  9.2059E+00  1.1186E+03  8.2110E+00  6.6040E+01  1.2024E+00  1.2887E+00  1.8772E+01  4.4710E-01
             2.1571E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2059.29021255932        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1775
 NPARAMETR:  9.7470E-01  1.8385E-02  1.9765E+00  1.6893E+00  9.5292E-01  1.0577E+00  7.4566E+00  1.4689E+00  8.0283E-01  9.5061E-01
             1.1138E+00
 PARAMETER:  7.4372E-02 -3.8962E+00  7.8134E-01  6.2430E-01  5.1771E-02  1.5605E-01  2.1091E+00  4.8454E-01 -1.1961E-01  4.9348E-02
             2.0780E-01
 GRADIENT:   7.5880E-01  1.4883E-01  4.7351E-01 -1.6520E+01  2.6221E-01  1.4716E-01 -1.9936E-02  4.6145E-02  1.7969E-01  2.3288E-01
            -2.4405E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2059.30514376732        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1940
 NPARAMETR:  9.7463E-01  1.3430E-02  1.9752E+00  1.6890E+00  9.5249E-01  1.0577E+00  1.0532E+01  1.4685E+00  8.0239E-01  9.4885E-01
             1.1138E+00
 PARAMETER:  7.4301E-02 -4.2103E+00  7.8065E-01  6.2416E-01  5.1322E-02  1.5605E-01  2.4545E+00  4.8427E-01 -1.2016E-01  4.7495E-02
             2.0781E-01
 GRADIENT:   8.4745E-01  1.3199E-01 -5.2557E-03 -2.3180E+01  2.0579E+00  1.8601E-01  5.8682E-02  5.5881E-02  6.4222E-01  4.5905E-03
             2.7218E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2059.31424278568        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2129
 NPARAMETR:  9.7478E-01  1.0000E-02  1.9728E+00  1.6915E+00  9.5002E-01  1.0579E+00  1.1227E+01  1.4673E+00  8.0068E-01  9.4768E-01
             1.1138E+00
 PARAMETER:  7.4451E-02 -4.5179E+00  7.7943E-01  6.2560E-01  4.8731E-02  1.5624E-01  2.5183E+00  4.8343E-01 -1.2230E-01  4.6259E-02
             2.0777E-01
 GRADIENT:   1.2435E+00  4.6469E-01  4.1946E-01 -2.2277E+01  6.3578E-01  2.9452E-01  8.7222E-03  5.6155E-02  5.3209E-02  1.1950E-01
            -1.1234E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2059.31750742720        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2320             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7482E-01  1.0000E-02  1.9618E+00  1.6907E+00  9.4799E-01  1.0579E+00  1.2233E+01  1.4623E+00  8.0084E-01  9.4292E-01
             1.1138E+00
 PARAMETER:  7.4497E-02 -4.5859E+00  7.7388E-01  6.2514E-01  4.6589E-02  1.5629E-01  2.6042E+00  4.8000E-01 -1.2209E-01  4.1221E-02
             2.0774E-01
 GRADIENT:   3.5574E+02  0.0000E+00  8.9771E+00  1.1511E+03  6.9777E+00  6.5749E+01  9.5323E-01  1.2662E+00  1.9590E+01  3.0009E-01
             1.9898E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2059.31975586100        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2508             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7481E-01  1.0000E-02  1.9592E+00  1.6905E+00  9.4693E-01  1.0579E+00  1.1559E+01  1.4602E+00  8.0082E-01  9.4326E-01
             1.1138E+00
 PARAMETER:  7.4484E-02 -4.5859E+00  7.7254E-01  6.2504E-01  4.5467E-02  1.5629E-01  2.5475E+00  4.7859E-01 -1.2212E-01  4.1587E-02
             2.0777E-01
 GRADIENT:   3.5570E+02  0.0000E+00  9.2095E+00  1.1506E+03  6.2562E+00  6.5741E+01  8.0049E-01  1.2774E+00  1.9521E+01  4.4415E-01
             2.0415E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2059.32068475231        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2696
 NPARAMETR:  9.7479E-01  1.0000E-02  1.9550E+00  1.6903E+00  9.4618E-01  1.0579E+00  1.1591E+01  1.4572E+00  8.0089E-01  9.4202E-01
             1.1138E+00
 PARAMETER:  7.4470E-02 -4.5859E+00  7.7040E-01  6.2492E-01  4.4682E-02  1.5626E-01  2.5503E+00  4.7650E-01 -1.2203E-01  4.0268E-02
             2.0774E-01
 GRADIENT:   1.4053E+00  0.0000E+00 -5.3914E-02 -2.2864E+01  1.3302E+00  3.1308E-01  1.7422E-02  8.2619E-03  3.3373E-02 -2.4149E-01
            -9.2698E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2059.32220348421        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2888             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7477E-01  1.0000E-02  1.9519E+00  1.6901E+00  9.4444E-01  1.0579E+00  1.1615E+01  1.4545E+00  8.0101E-01  9.4360E-01
             1.1138E+00
 PARAMETER:  7.4448E-02 -4.5859E+00  7.6882E-01  6.2478E-01  4.2836E-02  1.5624E-01  2.5523E+00  4.7465E-01 -1.2188E-01  4.1945E-02
             2.0781E-01
 GRADIENT:   3.5549E+02  0.0000E+00  9.6114E+00  1.1498E+03  4.9486E+00  6.5696E+01  8.1105E-01  1.2826E+00  1.9520E+01  7.1569E-01
             2.1171E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2059.32296941309        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3076
 NPARAMETR:  9.7477E-01  1.0000E-02  1.9491E+00  1.6900E+00  9.4410E-01  1.0578E+00  1.1614E+01  1.4530E+00  8.0103E-01  9.4272E-01
             1.1138E+00
 PARAMETER:  7.4444E-02 -4.5859E+00  7.6738E-01  6.2471E-01  4.2473E-02  1.5623E-01  2.5522E+00  4.7362E-01 -1.2185E-01  4.1015E-02
             2.0779E-01
 GRADIENT:   1.3728E+00  0.0000E+00  2.6982E-01 -2.2843E+01  1.0016E-01  3.2345E-01  1.8174E-02  6.7149E-02  6.7035E-02  5.5492E-02
             9.2914E-03

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2059.32339749479        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3268             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7475E-01  1.0000E-02  1.9426E+00  1.6897E+00  9.4348E-01  1.0578E+00  1.1615E+01  1.4494E+00  8.0109E-01  9.4038E-01
             1.1138E+00
 PARAMETER:  7.4430E-02 -4.5859E+00  7.6405E-01  6.2455E-01  4.1821E-02  1.5619E-01  2.5523E+00  4.7112E-01 -1.2178E-01  3.8532E-02
             2.0777E-01
 GRADIENT:   3.5550E+02  0.0000E+00  8.7335E+00  1.1490E+03  6.7755E+00  6.5634E+01  8.1203E-01  1.2298E+00  1.9471E+01  3.3742E-01
             2.0256E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2059.32426423818        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3456
 NPARAMETR:  9.7475E-01  1.0000E-02  1.9417E+00  1.6896E+00  9.4275E-01  1.0578E+00  1.1624E+01  1.4480E+00  8.0114E-01  9.4092E-01
             1.1138E+00
 PARAMETER:  7.4421E-02 -4.5859E+00  7.6358E-01  6.2449E-01  4.1050E-02  1.5619E-01  2.5530E+00  4.7021E-01 -1.2172E-01  3.9105E-02
             2.0777E-01
 GRADIENT:   1.4070E+00  0.0000E+00 -5.7617E-02 -2.2903E+01  6.6599E-01  3.1541E-01  1.8015E-02  2.6870E-02  5.3473E-02 -8.8100E-02
            -3.5882E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2059.32463644046        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3648             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7473E-01  1.0000E-02  1.9402E+00  1.6895E+00  9.4185E-01  1.0578E+00  1.1635E+01  1.4465E+00  8.0121E-01  9.4168E-01
             1.1138E+00
 PARAMETER:  7.4408E-02 -4.5859E+00  7.6280E-01  6.2442E-01  4.0094E-02  1.5618E-01  2.5540E+00  4.6918E-01 -1.2163E-01  3.9909E-02
             2.0781E-01
 GRADIENT:   3.5550E+02  0.0000E+00  9.3742E+00  1.1485E+03  5.1141E+00  6.5628E+01  8.1558E-01  1.2413E+00  1.9507E+01  6.7000E-01
             2.1081E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2059.32483419306        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3836
 NPARAMETR:  9.7473E-01  1.0000E-02  1.9388E+00  1.6894E+00  9.4154E-01  1.0578E+00  1.1637E+01  1.4456E+00  8.0123E-01  9.4144E-01
             1.1138E+00
 PARAMETER:  7.4404E-02 -4.5859E+00  7.6208E-01  6.2438E-01  3.9759E-02  1.5618E-01  2.5542E+00  4.6850E-01 -1.2161E-01  3.9656E-02
             2.0781E-01
 GRADIENT:   1.3706E+00  0.0000E+00  2.2903E-01 -2.2898E+01 -2.1067E-01  3.3740E-01  1.8452E-02  4.4020E-02  5.9823E-02  9.6225E-02
             2.8444E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -2059.32502259028        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     4028             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7473E-01  1.0000E-02  1.9360E+00  1.6893E+00  9.4142E-01  1.0578E+00  1.1632E+01  1.4440E+00  8.0124E-01  9.3998E-01
             1.1138E+00
 PARAMETER:  7.4401E-02 -4.5859E+00  7.6062E-01  6.2432E-01  3.9629E-02  1.5616E-01  2.5538E+00  4.6744E-01 -1.2160E-01  3.8104E-02
             2.0778E-01
 GRADIENT:   3.5550E+02  0.0000E+00  8.9779E+00  1.1482E+03  5.9907E+00  6.5603E+01  8.1553E-01  1.1969E+00  1.9475E+01  4.6309E-01
             2.0411E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -2059.32510371741        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4216
 NPARAMETR:  9.7472E-01  1.0000E-02  1.9349E+00  1.6892E+00  9.4117E-01  1.0578E+00  1.1633E+01  1.4432E+00  8.0126E-01  9.3963E-01
             1.1138E+00
 PARAMETER:  7.4397E-02 -4.5859E+00  7.6008E-01  6.2428E-01  3.9368E-02  1.5615E-01  2.5538E+00  4.6686E-01 -1.2157E-01  3.7736E-02
             2.0777E-01
 GRADIENT:   1.3656E+00  0.0000E+00 -1.5447E-01 -2.2898E+01  6.7469E-01  3.0884E-01  1.7934E-02 -5.3914E-03  4.2926E-02 -1.2981E-01
            -4.9864E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -2059.32540576522        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     4408             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7470E-01  1.0000E-02  1.9312E+00  1.6890E+00  9.3985E-01  1.0578E+00  1.1654E+01  1.4414E+00  8.0135E-01  9.4096E-01
             1.1139E+00
 PARAMETER:  7.4379E-02 -4.5859E+00  7.5813E-01  6.2415E-01  3.7968E-02  1.5615E-01  2.5556E+00  4.6560E-01 -1.2145E-01  3.9140E-02
             2.0786E-01
 GRADIENT:   3.5539E+02  0.0000E+00  9.1273E+00  1.1475E+03  5.1254E+00  6.5599E+01  8.1959E-01  1.2988E+00  1.9482E+01  7.5470E-01
             2.1789E+00

0ITERATION NO.:  124    OBJECTIVE VALUE:  -2059.32557557441        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     4548
 NPARAMETR:  9.7470E-01  1.0000E-02  1.9292E+00  1.6889E+00  9.3899E-01  1.0577E+00  1.1659E+01  1.4398E+00  8.0141E-01  9.4060E-01
             1.1139E+00
 PARAMETER:  7.4379E-02 -4.5859E+00  7.5807E-01  6.2414E-01  3.8053E-02  1.5614E-01  2.5550E+00  4.6490E-01 -1.2146E-01  3.7772E-02
             2.0780E-01
 GRADIENT:   1.0619E-03  0.0000E+00  3.5675E-02  1.4813E-02  1.2770E-01  2.0648E-03 -3.5104E-05  2.9767E-03 -2.9042E-03 -1.1597E-02
            -8.7793E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4548
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3897E-04 -2.7085E-04 -2.9163E-02 -5.1197E-03 -3.7725E-02
 SE:             2.9827E-02  1.7426E-03  1.8097E-02  2.9419E-02  2.0209E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9361E-01  8.7648E-01  1.0708E-01  8.6185E-01  6.1933E-02

 ETASHRINKSD(%)  7.5309E-02  9.4162E+01  3.9372E+01  1.4411E+00  3.2299E+01
 ETASHRINKVR(%)  1.5056E-01  9.9659E+01  6.3243E+01  2.8615E+00  5.4166E+01
 EBVSHRINKSD(%)  3.9160E-01  9.4259E+01  4.2658E+01  1.7897E+00  2.9816E+01
 EBVSHRINKVR(%)  7.8168E-01  9.9670E+01  6.7119E+01  3.5474E+00  5.0741E+01
 RELATIVEINF(%)  9.5525E+01  8.2644E-03  6.8338E+00  2.8084E+00  8.2952E+00
 EPSSHRINKSD(%)  3.3331E+01
 EPSSHRINKVR(%)  5.5552E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2059.3255755744131     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1140.3870423697404     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    79.88
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     7.69
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2059.326       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.00E-02  1.93E+00  1.69E+00  9.40E-01  1.06E+00  1.16E+01  1.44E+00  8.01E-01  9.40E-01  1.11E+00
 


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
 
         3.10E-02  0.00E+00  3.57E-01  4.37E-02  8.99E-02  8.76E-02  1.54E+00  3.45E-01  5.29E-02  1.38E-01  5.53E-02
 


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
+        9.60E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.07E-03  0.00E+00  1.27E-01
 
 TH 4
+       -1.07E-04  0.00E+00  5.30E-03  1.91E-03
 
 TH 5
+       -3.02E-04  0.00E+00  2.84E-02  1.32E-03  8.08E-03
 
 TH 6
+       -3.64E-04  0.00E+00 -3.97E-05  2.85E-04 -5.29E-04  7.68E-03
 
 TH 7
+       -9.57E-03  0.00E+00  1.07E-02 -9.13E-03  4.84E-03 -1.29E-02  2.38E+00
 
 TH 8
+       -2.34E-03  0.00E+00  8.19E-02  1.79E-03  1.63E-02  2.62E-03  4.44E-02  1.19E-01
 
 TH 9
+        2.90E-04  0.00E+00 -2.27E-03 -5.51E-04 -3.70E-04 -2.37E-04  1.63E-02 -2.88E-03  2.80E-03
 
 TH10
+       -2.65E-04  0.00E+00  1.28E-02  9.03E-04  5.32E-03  1.24E-03 -3.02E-02 -8.00E-03 -5.21E-04  1.92E-02
 
 TH11
+        2.80E-04  0.00E+00  1.40E-04 -1.68E-04  3.28E-06 -3.80E-04  1.30E-02 -1.20E-03  3.80E-04 -6.32E-04  3.06E-03
 
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
+        3.10E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.72E-02  0.00E+00  3.57E-01
 
 TH 4
+       -7.92E-02  0.00E+00  3.40E-01  4.37E-02
 
 TH 5
+       -1.09E-01  0.00E+00  8.86E-01  3.35E-01  8.99E-02
 
 TH 6
+       -1.34E-01  0.00E+00 -1.27E-03  7.44E-02 -6.72E-02  8.76E-02
 
 TH 7
+       -2.00E-01  0.00E+00  1.94E-02 -1.35E-01  3.49E-02 -9.50E-02  1.54E+00
 
 TH 8
+       -2.19E-01  0.00E+00  6.65E-01  1.19E-01  5.27E-01  8.67E-02  8.33E-02  3.45E-01
 
 TH 9
+        1.77E-01  0.00E+00 -1.20E-01 -2.38E-01 -7.78E-02 -5.11E-02  2.00E-01 -1.58E-01  5.29E-02
 
 TH10
+       -6.18E-02  0.00E+00  2.60E-01  1.49E-01  4.28E-01  1.02E-01 -1.41E-01 -1.67E-01 -7.11E-02  1.38E-01
 
 TH11
+        1.63E-01  0.00E+00  7.08E-03 -6.94E-02  6.60E-04 -7.84E-02  1.52E-01 -6.27E-02  1.30E-01 -8.24E-02  5.53E-02
 
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
+        1.10E+03
 
 TH 2
+        1.11E-14  1.33E-30
 
 TH 3
+       -1.93E+01  1.46E-15  3.59E+01
 
 TH 4
+       -2.15E+01 -5.81E-15 -1.24E+01  6.19E+02
 
 TH 5
+        7.02E+01 -6.06E-15 -1.60E+02 -6.49E+01  7.49E+02
 
 TH 6
+        5.39E+01 -5.34E-15 -1.64E+01 -2.09E+01  6.27E+01  1.34E+02
 
 TH 7
+       -1.41E-02 -3.94E-18  3.21E-03  1.09E-02 -1.02E-02  8.03E-04  8.61E-06
 
 TH 8
+       -9.95E-01 -3.50E-16 -3.68E-01 -3.06E+00 -1.01E-02  2.72E+00 -5.50E-04  4.90E-01
 
 TH 9
+       -1.05E+02 -2.66E-14  1.75E+01  9.70E+01 -5.38E+01  2.61E+00  5.79E-02 -3.96E+00  3.91E+02
 
 TH10
+       -6.66E+00  5.93E-17  1.64E+01  1.41E+00 -7.96E+01  8.59E-01  3.31E-04  7.51E-01  1.01E-01  9.71E+00
 
 TH11
+       -6.89E+01 -1.41E-14 -9.15E+00  3.36E+01 -1.41E+01  4.19E+00 -3.77E-03  1.17E+01 -2.70E+01  1.75E+01  3.47E+02
 
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
+        1.04E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.98E+01  0.00E+00  4.87E+01
 
 TH 4
+       -3.54E+01  0.00E+00 -8.11E+00  5.85E+02
 
 TH 5
+       -2.23E+01  0.00E+00 -1.53E+02 -5.01E+01  7.28E+02
 
 TH 6
+       -5.95E+01  0.00E+00  5.73E+00  2.97E+01 -8.89E+01  2.33E+02
 
 TH 7
+       -4.13E-03  0.00E+00  3.45E-04 -2.45E-02  1.18E-02 -3.14E-03  1.16E-05
 
 TH 8
+       -2.97E+01  0.00E+00 -1.30E+01 -2.04E+01 -2.36E+00  1.08E+01  3.28E-04  1.86E+01
 
 TH 9
+        8.56E+01  0.00E+00  4.45E+00 -8.03E+01  4.07E+01 -1.17E+01  4.10E-02 -7.09E+00  2.35E+02
 
 TH10
+       -2.12E+01  0.00E+00 -2.24E+00 -3.76E+00 -6.31E+01  2.96E+01  5.39E-04  1.18E+01 -1.09E+01  4.95E+01
 
 TH11
+        6.23E+01  0.00E+00 -3.14E+00 -5.85E+01  7.04E-01 -7.97E+00  1.12E-02  8.45E+00  3.91E+01  5.00E+00  3.21E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,       87.637
Stop Time:
Thu Sep 30 01:19:14 CDT 2021
