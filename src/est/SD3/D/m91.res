Sun Oct 24 01:07:18 CDT 2021
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
$DATA ../../../../data/SD3/D/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2021.15928008679        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9516E+02 -5.8896E+01 -5.9376E+01 -2.3994E+01  9.9313E+01 -7.9911E+00 -2.5966E+01  1.4189E+01 -5.1218E+01  1.4187E+01
            -1.5972E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2039.28339757137        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  9.5547E-01  1.0054E+00  1.1427E+00  1.0776E+00  9.8082E-01  1.1984E+00  1.1534E+00  9.1707E-01  1.3387E+00  8.8358E-01
             1.0451E+00
 PARAMETER:  5.4452E-02  1.0536E-01  2.3340E-01  1.7475E-01  8.0630E-02  2.8095E-01  2.4269E-01  1.3426E-02  3.9168E-01 -2.3774E-02
             1.4415E-01
 GRADIENT:   8.1580E+00  2.2500E+00 -4.1359E+00  9.5974E+00  4.0644E+00  3.2561E+01  2.8643E+00  4.6622E+00  2.8864E+01 -8.9418E-01
             1.5441E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2040.72827580794        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  9.5778E-01  7.6919E-01  1.1817E+00  1.2735E+00  8.8780E-01  1.1495E+00  2.0391E+00  5.6432E-01  1.0190E+00  8.1571E-01
             1.0191E+00
 PARAMETER:  5.6867E-02 -1.6242E-01  2.6697E-01  3.4177E-01 -1.9004E-02  2.3930E-01  8.1253E-01 -4.7214E-01  1.1880E-01 -1.0369E-01
             1.1894E-01
 GRADIENT:   1.7723E+01  4.3681E+01  3.8002E+01  4.4534E+01 -2.4605E+01  1.8641E+01  1.8760E+01 -5.3341E+00  7.0516E+00 -1.2709E+01
            -9.9308E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2046.93568225421        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  9.5030E-01  7.8866E-01  9.1908E-01  1.1883E+00  8.0105E-01  1.0878E+00  1.6651E+00  3.3837E-01  9.9720E-01  7.9041E-01
             1.0196E+00
 PARAMETER:  4.9024E-02 -1.3742E-01  1.5617E-02  2.7253E-01 -1.2183E-01  1.8417E-01  6.0989E-01 -9.8362E-01  9.7194E-02 -1.3520E-01
             1.1945E-01
 GRADIENT:   1.3490E+00  8.7984E+00  3.4014E+00  5.0606E+00 -7.6843E+00 -2.7292E+00  5.9838E-01 -3.7002E-01 -5.6284E-01  2.7313E-01
            -8.6978E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2048.04685628649        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.4583E-01  4.7888E-01  1.1291E+00  1.3930E+00  7.9185E-01  1.0901E+00  2.1455E+00  7.1209E-01  9.5032E-01  8.1554E-01
             1.0194E+00
 PARAMETER:  4.4311E-02 -6.3631E-01  2.2140E-01  4.3149E-01 -1.3339E-01  1.8631E-01  8.6338E-01 -2.3956E-01  4.9041E-02 -1.0391E-01
             1.1925E-01
 GRADIENT:   1.5351E+00  4.9859E+00  8.3729E-01  8.2747E+00 -7.4884E+00  1.9990E-01  2.5751E-01  1.3298E+00  1.0037E+00 -2.1060E-01
             1.9955E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2048.30512687862        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  9.4138E-01  3.0272E-01  1.2870E+00  1.5217E+00  8.0737E-01  1.0883E+00  2.3922E+00  8.0968E-01  9.2694E-01  9.1088E-01
             1.0149E+00
 PARAMETER:  3.9591E-02 -1.0949E+00  3.5231E-01  5.1984E-01 -1.1398E-01  1.8466E-01  9.7220E-01 -1.1111E-01  2.4137E-02  6.6510E-03
             1.1482E-01
 GRADIENT:  -1.4689E+00  5.3134E+00  2.7288E-01  1.8493E+01 -6.3905E+00  1.0114E+00  3.7897E-01  4.6616E-01 -2.6051E+00  3.7088E+00
            -4.9571E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2048.82193024991        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  9.3608E-01  1.6939E-01  1.4081E+00  1.6043E+00  8.1421E-01  1.0742E+00  2.8485E+00  9.2306E-01  9.0352E-01  9.2896E-01
             1.0156E+00
 PARAMETER:  3.3943E-02 -1.6755E+00  4.4223E-01  5.7269E-01 -1.0553E-01  1.7156E-01  1.1468E+00  1.9935E-02 -1.4616E-03  2.6313E-02
             1.1552E-01
 GRADIENT:  -8.1524E+00  3.0176E+00  5.6009E+00  7.9754E+00 -9.9762E+00 -3.4217E+00  2.4277E-01  3.7191E-01 -2.3845E+00  3.9805E+00
             1.3054E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2049.21096085576        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  9.3749E-01  6.5000E-02  1.4659E+00  1.6706E+00  8.0724E-01  1.0785E+00  4.0412E+00  9.8884E-01  8.8428E-01  8.9643E-01
             1.0158E+00
 PARAMETER:  3.5453E-02 -2.6334E+00  4.8249E-01  6.1320E-01 -1.1413E-01  1.7554E-01  1.4965E+00  8.8776E-02 -2.2978E-02 -9.3388E-03
             1.1568E-01
 GRADIENT:  -1.8787E+00  1.1630E+00  6.4305E+00  6.0931E+00 -7.1563E+00 -1.2136E+00  1.4673E-01 -6.6194E-01 -9.3878E-01 -7.8142E-01
            -5.8517E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2049.44130740351        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1416             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3820E-01  1.0000E-02  1.4546E+00  1.6769E+00  7.9776E-01  1.0817E+00  2.6404E-02  1.0063E+00  8.7624E-01  9.0512E-01
             1.0162E+00
 PARAMETER:  3.6207E-02 -4.6399E+00  4.7474E-01  6.1695E-01 -1.2594E-01  1.7850E-01 -3.5342E+00  1.0626E-01 -3.2117E-02  3.0723E-04
             1.1606E-01
 GRADIENT:   3.6690E+02  0.0000E+00  8.7083E+00  1.1019E+03  1.5105E+01  6.8267E+01  9.6001E-06  1.4763E+00  1.6235E+01  2.0460E+00
             2.3801E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2049.52695472273        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1600             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3787E-01  1.0000E-02  1.4570E+00  1.6896E+00  7.9581E-01  1.0813E+00  2.4217E-02  9.9051E-01  8.7397E-01  8.9707E-01
             1.0156E+00
 PARAMETER:  3.5859E-02 -4.6399E+00  4.7635E-01  6.2451E-01 -1.2840E-01  1.7813E-01 -3.6207E+00  9.0461E-02 -3.4711E-02 -8.6161E-03
             1.1548E-01
 GRADIENT:   3.6648E+02  0.0000E+00  1.2117E+01  1.1382E+03  1.0233E+01  6.8093E+01  8.2087E-06 -2.2754E-02  1.5714E+01  4.2486E-01
             1.0075E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2049.53097303458        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1789             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3793E-01  1.0000E-02  1.4541E+00  1.6886E+00  7.9502E-01  1.0813E+00  2.1380E-02  9.9105E-01  8.7375E-01  8.9660E-01
             1.0156E+00
 PARAMETER:  3.5920E-02 -4.6399E+00  4.7439E-01  6.2392E-01 -1.2939E-01  1.7817E-01 -3.7453E+00  9.1007E-02 -3.4963E-02 -9.1402E-03
             1.1545E-01
 GRADIENT:   3.6664E+02  0.0000E+00  1.1731E+01  1.1357E+03  1.0582E+01  6.8123E+01  6.6374E-06  1.3023E-01  1.5576E+01  5.1294E-01
             1.0648E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2049.53339573879        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1981
 NPARAMETR:  9.3795E-01  1.0000E-02  1.4517E+00  1.6879E+00  7.9431E-01  1.0813E+00  2.5101E-02  9.8976E-01  8.7380E-01  8.9644E-01
             1.0156E+00
 PARAMETER:  3.5944E-02 -4.6399E+00  4.7274E-01  6.2347E-01 -1.3029E-01  1.7819E-01 -3.5849E+00  8.9709E-02 -3.4906E-02 -9.3245E-03
             1.1543E-01
 GRADIENT:   1.2815E+00  0.0000E+00  7.7175E-01 -1.8778E+01  9.1424E-01  2.8398E-01  6.8617E-07 -6.8362E-02  3.1307E-02 -1.2155E-01
            -6.6748E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2049.53621702896        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2173             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3794E-01  1.0000E-02  1.4480E+00  1.6876E+00  7.9317E-01  1.0813E+00  1.0000E-02  9.8963E-01  8.7386E-01  8.9707E-01
             1.0156E+00
 PARAMETER:  3.5931E-02 -4.6399E+00  4.7019E-01  6.2329E-01 -1.3172E-01  1.7818E-01 -2.4977E+01  8.9575E-02 -3.4833E-02 -8.6214E-03
             1.1548E-01
 GRADIENT:   3.6657E+02  0.0000E+00  1.1163E+01  1.1330E+03  1.0751E+01  6.8096E+01  0.0000E+00  3.7236E-01  1.5554E+01  8.4995E-01
             1.2397E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2049.53711293069        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2365
 NPARAMETR:  9.3793E-01  1.0000E-02  1.4454E+00  1.6874E+00  7.9238E-01  1.0813E+00  1.0000E-02  9.9008E-01  8.7390E-01  8.9768E-01
             1.0157E+00
 PARAMETER:  3.5924E-02 -4.6399E+00  4.6837E-01  6.2316E-01 -1.3271E-01  1.7817E-01 -2.4977E+01  9.0035E-02 -3.4784E-02 -7.9429E-03
             1.1554E-01
 GRADIENT:   1.2985E+00  0.0000E+00  2.2720E-02 -1.8763E+01  7.2369E-01  2.8712E-01  0.0000E+00  3.2373E-01  2.1286E-02  4.2681E-01
             2.3827E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2049.53936369088        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2563             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3792E-01  1.0000E-02  1.4437E+00  1.6872E+00  7.9167E-01  1.0813E+00  1.0000E-02  9.8435E-01  8.7396E-01  8.9517E-01
             1.0155E+00
 PARAMETER:  3.5914E-02 -4.6399E+00  4.6722E-01  6.2304E-01 -1.3361E-01  1.7816E-01 -2.4977E+01  8.4222E-02 -3.4724E-02 -1.0742E-02
             1.1541E-01
 GRADIENT:   3.6657E+02  0.0000E+00  1.1293E+01  1.1323E+03  1.0768E+01  6.8079E+01  0.0000E+00  1.7770E-01  1.5563E+01  5.8133E-01
             1.0981E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2049.54021433573        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2755
 NPARAMETR:  9.3792E-01  1.0000E-02  1.4420E+00  1.6870E+00  7.9103E-01  1.0813E+00  1.0000E-02  9.8191E-01  8.7400E-01  8.9446E-01
             1.0155E+00
 PARAMETER:  3.5907E-02 -4.6399E+00  4.6601E-01  6.2294E-01 -1.3442E-01  1.7815E-01 -2.4977E+01  8.1747E-02 -3.4674E-02 -1.1539E-02
             1.1536E-01
 GRADIENT:   1.2576E+00  0.0000E+00  6.2985E-01 -1.8836E+01  5.0970E-01  2.7942E-01  0.0000E+00 -1.4846E-01  3.7382E-02 -2.0465E-01
            -1.3391E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2049.54166556506        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2953             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3791E-01  1.0000E-02  1.4392E+00  1.6868E+00  7.9024E-01  1.0813E+00  1.0000E-02  9.8281E-01  8.7404E-01  8.9533E-01
             1.0156E+00
 PARAMETER:  3.5899E-02 -4.6399E+00  4.6411E-01  6.2281E-01 -1.3542E-01  1.7814E-01 -2.4977E+01  8.2657E-02 -3.4627E-02 -1.0566E-02
             1.1543E-01
 GRADIENT:   3.6644E+02  0.0000E+00  1.0951E+01  1.1315E+03  1.0763E+01  6.8052E+01  0.0000E+00  3.1565E-01  1.5559E+01  7.9152E-01
             1.2030E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2049.54212611002        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3145
 NPARAMETR:  9.3791E-01  1.0000E-02  1.4378E+00  1.6866E+00  7.8978E-01  1.0813E+00  1.0000E-02  9.8276E-01  8.7407E-01  8.9524E-01
             1.0156E+00
 PARAMETER:  3.5895E-02 -4.6399E+00  4.6308E-01  6.2274E-01 -1.3600E-01  1.7814E-01 -2.4977E+01  8.2606E-02 -3.4597E-02 -1.0666E-02
             1.1547E-01
 GRADIENT:   1.2570E+00  0.0000E+00  3.3332E-02 -1.8811E+01  4.7026E-01  2.7837E-01  0.0000E+00  1.4847E-01  3.0370E-02  1.7916E-01
             1.1155E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2049.54277660961        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     3343             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3790E-01  1.0000E-02  1.4363E+00  1.6864E+00  7.8910E-01  1.0813E+00  1.0000E-02  9.7808E-01  8.7412E-01  8.9349E-01
             1.0155E+00
 PARAMETER:  3.5886E-02 -4.6399E+00  4.6204E-01  6.2262E-01 -1.3687E-01  1.7813E-01 -2.4977E+01  7.7831E-02 -3.4542E-02 -1.2620E-02
             1.1534E-01
 GRADIENT:   3.6647E+02  0.0000E+00  1.1229E+01  1.1310E+03  1.0667E+01  6.8038E+01  0.0000E+00  9.2373E-02  1.5573E+01  4.8882E-01
             1.0272E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2049.54323485971        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3535
 NPARAMETR:  9.3790E-01  1.0000E-02  1.4354E+00  1.6864E+00  7.8886E-01  1.0813E+00  1.0000E-02  9.7858E-01  8.7413E-01  8.9394E-01
             1.0155E+00
 PARAMETER:  3.5883E-02 -4.6399E+00  4.6141E-01  6.2258E-01 -1.3717E-01  1.7813E-01 -2.4977E+01  7.8352E-02 -3.4525E-02 -1.2112E-02
             1.1537E-01
 GRADIENT:   1.2710E+00  0.0000E+00  2.9165E-01 -1.8861E+01  2.6094E-01  2.7834E-01  0.0000E+00 -4.3051E-02  3.8132E-02 -5.0807E-02
            -4.9890E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2049.54358063275        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     3733             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3789E-01  1.0000E-02  1.4337E+00  1.6863E+00  7.8845E-01  1.0813E+00  1.0000E-02  9.7926E-01  8.7415E-01  8.9415E-01
             1.0156E+00
 PARAMETER:  3.5880E-02 -4.6399E+00  4.6027E-01  6.2251E-01 -1.3769E-01  1.7812E-01 -2.4977E+01  7.9038E-02 -3.4504E-02 -1.1885E-02
             1.1544E-01
 GRADIENT:   3.6636E+02  0.0000E+00  1.0680E+01  1.1304E+03  1.0913E+01  6.8016E+01  0.0000E+00  3.2131E-01  1.5559E+01  7.6531E-01
             1.2232E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2049.54378722317        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3925
 NPARAMETR:  9.3789E-01  1.0000E-02  1.4333E+00  1.6862E+00  7.8822E-01  1.0813E+00  1.0000E-02  9.7814E-01  8.7417E-01  8.9393E-01
             1.0155E+00
 PARAMETER:  3.5877E-02 -4.6399E+00  4.5997E-01  6.2248E-01 -1.3798E-01  1.7812E-01 -2.4977E+01  7.7899E-02 -3.4484E-02 -1.2131E-02
             1.1541E-01
 GRADIENT:   1.2708E+00  0.0000E+00  9.1754E-02 -1.8861E+01  2.4858E-01  2.8443E-01  0.0000E+00  3.2618E-02  3.9120E-02  4.1075E-02
             2.0311E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -2049.54418341206        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     4123             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3788E-01  1.0000E-02  1.4311E+00  1.6860E+00  7.8737E-01  1.0812E+00  1.0000E-02  9.7540E-01  8.7423E-01  8.9298E-01
             1.0155E+00
 PARAMETER:  3.5869E-02 -4.6399E+00  4.5847E-01  6.2235E-01 -1.3905E-01  1.7811E-01 -2.4977E+01  7.5091E-02 -3.4414E-02 -1.3189E-02
             1.1536E-01
 GRADIENT:   3.6640E+02  0.0000E+00  1.0997E+01  1.1300E+03  1.0607E+01  6.8014E+01  0.0000E+00  1.6197E-01  1.5574E+01  5.9246E-01
             1.0953E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -2049.54433144039        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     4315
 NPARAMETR:  9.3788E-01  1.0000E-02  1.4304E+00  1.6859E+00  7.8715E-01  1.0812E+00  1.0000E-02  9.7540E-01  8.7423E-01  8.9308E-01
             1.0155E+00
 PARAMETER:  3.5864E-02 -4.6399E+00  4.5793E-01  6.2231E-01 -1.3933E-01  1.7810E-01 -2.4977E+01  7.5091E-02 -3.4412E-02 -1.3082E-02
             1.1535E-01
 GRADIENT:   1.2612E+00  0.0000E+00  1.6503E-01 -1.8890E+01  3.0537E-03  2.7026E-01  0.0000E+00 -2.7464E-02  3.9665E-02 -3.0809E-02
            -4.6164E-02

0ITERATION NO.:  116    OBJECTIVE VALUE:  -2049.54433144039        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     4345
 NPARAMETR:  9.3788E-01  1.0000E-02  1.4295E+00  1.6859E+00  7.8713E-01  1.0812E+00  1.0000E-02  9.7637E-01  8.7423E-01  8.9350E-01
             1.0156E+00
 PARAMETER:  3.5864E-02 -4.6399E+00  4.5793E-01  6.2231E-01 -1.3933E-01  1.7810E-01 -2.4977E+01  7.5091E-02 -3.4412E-02 -1.3082E-02
             1.1535E-01
 GRADIENT:  -4.0958E-03  0.0000E+00  1.9627E-01  1.2576E-02  2.3479E-02 -1.8036E-03  0.0000E+00 -2.8094E-02 -6.0606E-04 -3.1349E-02
            -4.5653E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4345
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8894E-04 -3.0763E-06 -2.1535E-02 -3.7222E-03 -2.6362E-02
 SE:             2.9861E-02  1.5407E-06  1.6448E-02  2.9608E-02  2.2038E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8961E-01  4.5855E-02  1.9045E-01  8.9996E-01  2.3161E-01

 ETASHRINKSD(%)  1.0000E-10  9.9995E+01  4.4896E+01  8.0823E-01  2.6170E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.9636E+01  1.6099E+00  4.5492E+01
 EBVSHRINKSD(%)  3.2358E-01  9.9995E+01  4.7124E+01  1.2111E+00  2.4426E+01
 EBVSHRINKVR(%)  6.4612E-01  1.0000E+02  7.2042E+01  2.4075E+00  4.2885E+01
 RELATIVEINF(%)  9.7558E+01  2.0734E-08  5.5163E+00  1.0680E+01  8.9085E+00
 EPSSHRINKSD(%)  3.3692E+01
 EPSSHRINKVR(%)  5.6032E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2049.5443314403906     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1130.6057982357179     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2049.544       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.38E-01  1.00E-02  1.43E+00  1.69E+00  7.87E-01  1.08E+00  1.00E-02  9.75E-01  8.74E-01  8.93E-01  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,      156.608
Stop Time:
Sun Oct 24 01:07:44 CDT 2021
