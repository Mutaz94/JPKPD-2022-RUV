Sun Oct 24 04:15:35 CDT 2021
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
$DATA ../../../../data/SD4/D/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1653.62327913869        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3339E+02  5.7981E+01  1.8612E+00  1.1063E+02  5.4461E+00  5.5367E+01 -1.4478E+01 -7.6439E+00 -1.3418E+01 -4.4755E-01
             7.9125E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1655.43653244625        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0010E+00  1.0624E+00  1.0082E+00  9.5702E-01  1.0390E+00  1.0036E+00  1.2863E+00  1.1486E+00  1.1823E+00  9.9468E-01
             9.6055E-01
 PARAMETER:  1.0096E-01  1.6050E-01  1.0816E-01  5.6064E-02  1.3830E-01  1.0358E-01  3.5179E-01  2.3857E-01  2.6748E-01  9.4666E-02
             5.9747E-02
 GRADIENT:  -2.5404E+01  7.0320E+00 -6.2080E+00  2.2121E+01  2.8742E+00  1.1440E-01  1.1646E+01 -1.8349E+00  2.3193E+01  2.5992E-01
            -4.0272E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1656.54243118007        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0037E+00  1.0884E+00  1.4015E+00  9.5712E-01  1.2048E+00  1.0155E+00  1.1267E+00  1.9441E+00  1.2034E+00  1.1456E+00
             9.4973E-01
 PARAMETER:  1.0372E-01  1.8470E-01  4.3751E-01  5.6175E-02  2.8627E-01  1.1538E-01  2.1929E-01  7.6479E-01  2.8516E-01  2.3593E-01
             4.8424E-02
 GRADIENT:  -1.6344E+01  5.2993E+00 -1.0943E+01  2.3734E+01  7.8752E+00  5.5489E+00  7.3810E+00  7.9921E+00  1.6483E+01  3.7752E+00
            -7.9449E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1658.29553729970        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0103E+00  9.4594E-01  1.4366E+00  1.0307E+00  1.1418E+00  9.9898E-01  1.1760E+00  1.6865E+00  1.0483E+00  1.0841E+00
             9.7227E-01
 PARAMETER:  1.1029E-01  4.4424E-02  4.6227E-01  1.3020E-01  2.3260E-01  9.8983E-02  2.6210E-01  6.2267E-01  1.4716E-01  1.8072E-01
             7.1875E-02
 GRADIENT:  -1.7268E-01  3.9087E-01  1.1058E-01 -8.3795E-01 -2.3765E-01 -1.0846E-02  4.1527E-01  4.6300E-02 -4.6996E-01  9.6385E-04
             1.2086E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1658.30947100563        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0102E+00  8.7215E-01  1.5442E+00  1.0841E+00  1.1393E+00  9.9845E-01  1.1899E+00  1.7378E+00  1.0325E+00  1.0932E+00
             9.7072E-01
 PARAMETER:  1.1012E-01 -3.6789E-02  5.3454E-01  1.8071E-01  2.3040E-01  9.8451E-02  2.7390E-01  6.5263E-01  1.3197E-01  1.8915E-01
             7.0281E-02
 GRADIENT:   1.1598E+00  4.4104E+00  2.0556E+00  4.0574E+00 -4.2933E+00  1.6487E-01  6.6174E-01 -4.9255E-01  2.6732E-01  5.8010E-01
            -4.9111E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1658.35481372212        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      905
 NPARAMETR:  1.0093E+00  7.7467E-01  1.6935E+00  1.1495E+00  1.1422E+00  9.9706E-01  1.1090E+00  1.8316E+00  1.0264E+00  1.1042E+00
             9.6992E-01
 PARAMETER:  1.0926E-01 -1.5532E-01  6.2680E-01  2.3933E-01  2.3300E-01  9.7051E-02  2.0344E-01  7.0520E-01  1.2610E-01  1.9915E-01
             6.9457E-02
 GRADIENT:   1.4469E+00  4.9265E+00  2.5499E+00  6.2247E+00 -5.0012E+00  1.4730E-01  8.3835E-03 -7.1883E-01  9.3327E-01  5.1973E-01
            -7.3283E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1658.40301330667        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  1.0076E+00  6.5638E-01  1.8084E+00  1.2238E+00  1.1344E+00  9.9507E-01  1.0144E+00  1.8881E+00  1.0003E+00  1.1064E+00
             9.7052E-01
 PARAMETER:  1.0757E-01 -3.2102E-01  6.9244E-01  3.0193E-01  2.2614E-01  9.5059E-02  1.1433E-01  7.3558E-01  1.0026E-01  2.0114E-01
             7.0072E-02
 GRADIENT:   3.8324E-01  2.7482E+00  1.3014E+00  3.3567E+00 -2.2623E+00  3.5110E-02 -4.7690E-01 -4.7787E-01  4.3901E-01  1.2672E-01
            -3.4595E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1658.44905706343        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1262
 NPARAMETR:  1.0071E+00  6.4839E-01  1.8015E+00  1.2272E+00  1.1331E+00  9.9476E-01  1.1347E+00  1.8778E+00  9.8294E-01  1.1042E+00
             9.7124E-01
 PARAMETER:  1.0706E-01 -3.3326E-01  6.8864E-01  3.0472E-01  2.2495E-01  9.4751E-02  2.2635E-01  7.3008E-01  8.2788E-02  1.9916E-01
             7.0820E-02
 GRADIENT:  -4.8204E-01  9.5409E-01  4.8096E-01 -1.1113E+00  1.7556E-01 -1.4571E-02 -2.9521E-02 -8.4661E-02  4.4104E-02  1.5826E-01
            -3.6062E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1658.45402759816        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1437
 NPARAMETR:  1.0067E+00  6.2131E-01  1.7887E+00  1.2453E+00  1.1212E+00  9.9417E-01  1.2273E+00  1.8444E+00  9.5990E-01  1.0965E+00
             9.7159E-01
 PARAMETER:  1.0670E-01 -3.7593E-01  6.8149E-01  3.1935E-01  2.1442E-01  9.4150E-02  3.0479E-01  7.1213E-01  5.9071E-02  1.9208E-01
             7.1178E-02
 GRADIENT:  -5.7709E-01  1.5640E+00  5.2511E-01 -4.9543E-01 -2.2782E-02 -8.8762E-02  3.3844E-02 -3.2825E-01 -1.5102E+00  7.3122E-02
            -5.5837E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1658.45439097068        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1616
 NPARAMETR:  1.0065E+00  6.0324E-01  1.7843E+00  1.2569E+00  1.1143E+00  9.9388E-01  1.2581E+00  1.8295E+00  9.5098E-01  1.0924E+00
             9.7168E-01
 PARAMETER:  1.0645E-01 -4.0544E-01  6.7905E-01  3.2864E-01  2.0825E-01  9.3859E-02  3.2964E-01  7.0406E-01  4.9738E-02  1.8840E-01
             7.1272E-02
 GRADIENT:  -6.5433E-01  1.7197E+00  4.8556E-01 -2.8738E-01 -2.2436E-01 -9.9432E-02  3.3594E-02 -4.3484E-01 -1.8722E+00  2.0929E-02
            -6.5080E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1658.46856133757        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1802             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0074E+00  5.8749E-01  1.7818E+00  1.2629E+00  1.1107E+00  9.9418E-01  1.2729E+00  1.8282E+00  9.5279E-01  1.0904E+00
             9.7180E-01
 PARAMETER:  1.0734E-01 -4.3189E-01  6.7765E-01  3.3339E-01  2.0502E-01  9.4166E-02  3.4133E-01  7.0332E-01  5.1638E-02  1.8654E-01
             7.1398E-02
 GRADIENT:   5.3175E+02  7.0453E+01  8.8346E+00  4.5713E+02  1.3177E+01  5.6855E+01  6.8177E+00  4.6509E+00  1.1459E+01  1.4620E+00
             9.8562E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1658.47721187182        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1981
 NPARAMETR:  1.0061E+00  5.7580E-01  1.7910E+00  1.2737E+00  1.1084E+00  9.9371E-01  1.2022E+00  1.8341E+00  9.5552E-01  1.0908E+00
             9.7131E-01
 PARAMETER:  1.0611E-01 -4.5199E-01  6.8277E-01  3.4191E-01  2.0291E-01  9.3689E-02  2.8418E-01  7.0655E-01  5.4502E-02  1.8694E-01
             7.0887E-02
 GRADIENT:  -6.8209E-01  1.1403E+00 -4.6095E-01 -3.6976E-01  1.8275E-01 -2.5474E-02 -4.2424E-02 -3.2857E-02 -7.0085E-02 -5.2598E-03
            -4.1266E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1658.48359883542        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2164
 NPARAMETR:  1.0068E+00  5.6850E-01  1.7969E+00  1.2767E+00  1.1076E+00  9.9382E-01  1.2004E+00  1.8369E+00  9.5524E-01  1.0910E+00
             9.7129E-01
 PARAMETER:  1.0680E-01 -4.6476E-01  6.8609E-01  3.4430E-01  2.0217E-01  9.3800E-02  2.8266E-01  7.0806E-01  5.4207E-02  1.8705E-01
             7.0874E-02
 GRADIENT:   1.1034E+00  3.0636E-01 -3.6132E-01 -3.0868E+00  1.9259E-01  6.7238E-02  1.0781E-02 -4.2554E-02  3.7835E-01  2.4432E-02
            -3.4197E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1658.48487734505        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2348             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0073E+00  5.6742E-01  1.8005E+00  1.2764E+00  1.1077E+00  9.9410E-01  1.2025E+00  1.8404E+00  9.5263E-01  1.0902E+00
             9.7131E-01
 PARAMETER:  1.0730E-01 -4.6665E-01  6.8804E-01  3.4406E-01  2.0231E-01  9.4082E-02  2.8437E-01  7.0997E-01  5.1474E-02  1.8641E-01
             7.0890E-02
 GRADIENT:   5.3246E+02  7.2481E+01  9.1473E+00  4.8101E+02  1.1860E+01  5.6937E+01  4.5402E+00  4.7414E+00  1.1329E+01  1.4067E+00
             7.7341E-01

0ITERATION NO.:   68    OBJECTIVE VALUE:  -1658.48501308701        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     2442
 NPARAMETR:  1.0070E+00  5.6748E-01  1.8009E+00  1.2772E+00  1.1077E+00  9.9397E-01  1.2059E+00  1.8404E+00  9.5317E-01  1.0904E+00
             9.7133E-01
 PARAMETER:  1.0702E-01 -4.6656E-01  6.8829E-01  3.4471E-01  2.0229E-01  9.3948E-02  2.8722E-01  7.0998E-01  5.2039E-02  1.8652E-01
             7.0908E-02
 GRADIENT:   1.6334E+00  3.0720E-01 -1.8713E-01 -3.5573E+00 -9.1836E-02  1.3631E-01 -8.8233E-03 -2.4401E-02 -2.1016E-02 -3.5311E-02
            -1.8164E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2442
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.3237E-04 -1.5171E-02 -4.3486E-02 -6.6808E-04 -5.0624E-02
 SE:             2.9893E-02  1.0427E-02  1.8211E-02  2.7730E-02  2.0411E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8579E-01  1.4568E-01  1.6943E-02  9.8078E-01  1.3131E-02

 ETASHRINKSD(%)  1.0000E-10  6.5068E+01  3.8992E+01  7.0997E+00  3.1619E+01
 ETASHRINKVR(%)  1.0000E-10  8.7797E+01  6.2781E+01  1.3695E+01  5.3241E+01
 EBVSHRINKSD(%)  4.0722E-01  6.5306E+01  4.2938E+01  7.4795E+00  2.7868E+01
 EBVSHRINKVR(%)  8.1278E-01  8.7964E+01  6.7440E+01  1.4400E+01  4.7969E+01
 RELATIVEINF(%)  9.7037E+01  4.1325E-01  1.0547E+01  3.3015E+00  1.2007E+01
 EPSSHRINKSD(%)  4.6032E+01
 EPSSHRINKVR(%)  7.0874E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1658.4850130870120     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -923.33418652327384     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1658.485       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  5.67E-01  1.80E+00  1.28E+00  1.11E+00  9.94E-01  1.21E+00  1.84E+00  9.53E-01  1.09E+00  9.71E-01
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       74.819
Stop Time:
Sun Oct 24 04:15:49 CDT 2021
