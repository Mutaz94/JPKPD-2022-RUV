Wed Sep 29 23:07:47 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat20.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1218.12401441602        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7442E+02  5.4052E-01  7.2449E+01  1.5738E+01  1.2108E+02  4.4000E+01 -5.8939E+00 -1.6967E+02 -6.1633E+01 -2.0354E+01
            -1.4488E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1701.76991677911        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0866E+00  1.1576E+00  7.8690E-01  8.6375E-01  8.3533E-01  9.9721E-01  9.0457E-01  2.7296E+00  1.0891E+00  6.8765E-01
             2.0180E+00
 PARAMETER:  1.8305E-01  2.4632E-01 -1.3965E-01 -4.6470E-02 -7.9929E-02  9.7210E-02 -2.9817E-04  1.1042E+00  1.8536E-01 -2.7448E-01
             8.0213E-01
 GRADIENT:   4.1428E+02 -5.1707E+01  2.6364E+01 -6.8036E+01 -5.5862E+01  3.9193E+00  1.2887E+01  2.4926E+01 -6.5137E-01  1.7905E+01
            -4.9733E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1725.04945378230        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0822E+00  1.4421E+00  6.8001E-01  7.4869E-01  9.1436E-01  9.7156E-01  6.4339E-01  2.7709E+00  1.3701E+00  3.1987E-01
             2.0501E+00
 PARAMETER:  1.7900E-01  4.6611E-01 -2.8565E-01 -1.8943E-01  1.0468E-02  7.1147E-02 -3.4100E-01  1.1192E+00  4.1489E-01 -1.0398E+00
             8.1791E-01
 GRADIENT:   4.0011E+02  9.2336E+01  3.0000E+01  1.4092E+01 -5.2396E+01 -3.9326E+00  8.2236E+00  4.3337E+00  7.9402E+00  2.0183E+00
            -2.1869E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1744.07241713402        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      256
 NPARAMETR:  9.8821E-01  1.2872E+00  3.4400E-01  7.9200E-01  6.8372E-01  8.8285E-01  7.7103E-01  2.0152E+00  1.1963E+00  1.6098E-01
             1.9809E+00
 PARAMETER:  8.8135E-02  3.5249E-01 -9.6713E-01 -1.3319E-01 -2.8020E-01 -2.4597E-02 -1.6003E-01  8.0072E-01  2.7922E-01 -1.7265E+00
             7.8357E-01
 GRADIENT:   2.8125E+01 -1.8414E+01  3.8463E+00 -1.5588E+01 -2.1110E+01 -2.7612E+01 -2.3760E+00 -1.5741E+01 -2.1234E+00  2.5000E-01
             1.6341E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1750.30566598495        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      422
 NPARAMETR:  9.7359E-01  1.3917E+00  4.0103E-01  7.4949E-01  7.7141E-01  9.3164E-01  6.9060E-01  2.5353E+00  1.3179E+00  1.9985E-01
             1.9319E+00
 PARAMETER:  7.3232E-02  4.3051E-01 -8.1371E-01 -1.8837E-01 -1.5954E-01  2.9190E-02 -2.7019E-01  1.0303E+00  3.7607E-01 -1.5102E+00
             7.5849E-01
 GRADIENT:   9.1019E+01  6.7095E+01  1.7401E+01  1.0955E+01 -2.2607E+01  2.2231E+00  2.0852E+00  1.1812E+01  3.1519E+00 -2.4484E-02
             1.1140E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1750.70811983371        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      599
 NPARAMETR:  9.7706E-01  1.3912E+00  4.0121E-01  7.5432E-01  7.7154E-01  9.4510E-01  6.3014E-01  2.6271E+00  1.3835E+00  3.6285E-01
             1.8964E+00
 PARAMETER:  7.6791E-02  4.3017E-01 -8.1327E-01 -1.8193E-01 -1.5937E-01  4.3535E-02 -3.6181E-01  1.0659E+00  4.2464E-01 -9.1378E-01
             7.3996E-01
 GRADIENT:   8.2724E-01 -7.1209E+00  7.4921E+00  9.4813E-01 -3.7001E+01  4.6001E-01  3.3943E-01  4.1442E-01 -2.5067E-01 -5.9601E-03
             4.9830E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1750.74013064884        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      775
 NPARAMETR:  9.7707E-01  1.3922E+00  4.0002E-01  7.5114E-01  7.7180E-01  9.4311E-01  6.1856E-01  2.6148E+00  1.3956E+00  4.0141E-01
             1.8901E+00
 PARAMETER:  7.6800E-02  4.3091E-01 -8.1623E-01 -1.8617E-01 -1.5903E-01  4.1425E-02 -3.8036E-01  1.0612E+00  4.3332E-01 -8.1277E-01
             7.3663E-01
 GRADIENT:   1.0691E+00 -9.7061E+00  7.4428E+00 -1.4032E+00 -3.7782E+01 -3.8477E-01  3.4747E-01 -3.7388E-01 -5.3035E-02  3.3789E-01
             3.4373E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1751.38656119762        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      964
 NPARAMETR:  9.8138E-01  1.4202E+00  3.6869E-01  7.3170E-01  7.8024E-01  9.4671E-01  6.7146E-01  2.4953E+00  1.3648E+00  3.7506E-01
             1.9007E+00
 PARAMETER:  8.1208E-02  4.5079E-01 -8.9779E-01 -2.1239E-01 -1.4815E-01  4.5236E-02 -2.9830E-01  1.0144E+00  4.1098E-01 -8.8068E-01
             7.4221E-01
 GRADIENT:   1.1128E+01 -2.7802E+01 -1.1645E+00 -3.3157E-02 -1.1506E+01  6.5983E-01  5.5273E-01 -4.9568E+00 -1.4752E+00  1.6951E-01
            -1.1673E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1751.68398683667        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1139
 NPARAMETR:  9.7792E-01  1.4346E+00  3.6912E-01  7.2270E-01  7.8118E-01  9.4440E-01  6.3403E-01  2.5255E+00  1.3904E+00  4.0560E-01
             1.9016E+00
 PARAMETER:  7.7677E-02  4.6089E-01 -8.9663E-01 -2.2477E-01 -1.4695E-01  4.2799E-02 -3.5566E-01  1.0265E+00  4.2959E-01 -8.0238E-01
             7.4271E-01
 GRADIENT:   2.3659E+00 -1.3606E+01  2.0885E+00 -2.6543E+00 -2.9582E+01 -2.0078E-01 -1.3539E+00 -4.2272E+00 -2.8194E+00  2.8916E-01
             5.5674E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1751.78690808185        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1318
 NPARAMETR:  9.7763E-01  1.4388E+00  3.6938E-01  7.2611E-01  7.8237E-01  9.4589E-01  6.3595E-01  2.5341E+00  1.3986E+00  4.0363E-01
             1.9032E+00
 PARAMETER:  7.7375E-02  4.6379E-01 -8.9593E-01 -2.2006E-01 -1.4543E-01  4.4373E-02 -3.5263E-01  1.0298E+00  4.3545E-01 -8.0726E-01
             7.4356E-01
 GRADIENT:   1.4051E+00 -6.0830E+00  2.1019E+00  3.6051E+00 -3.0898E+01  3.7268E-01 -9.2806E-01 -3.8918E+00 -1.6908E+00  2.7626E-01
             1.3554E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1751.95149055423        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1497
 NPARAMETR:  9.7802E-01  1.4437E+00  3.6890E-01  7.2207E-01  7.8458E-01  9.4541E-01  6.3776E-01  2.5477E+00  1.4104E+00  4.0247E-01
             1.9029E+00
 PARAMETER:  7.7770E-02  4.6724E-01 -8.9722E-01 -2.2563E-01 -1.4261E-01  4.3862E-02 -3.4979E-01  1.0352E+00  4.4386E-01 -8.1014E-01
             7.4340E-01
 GRADIENT:   2.4793E+00 -7.1418E+00  2.3061E+00  2.4285E+00 -3.2198E+01  2.2131E-01 -4.0020E-01 -3.3568E+00 -5.9484E-01  2.3574E-01
             1.3874E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1752.11165021395        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1676
 NPARAMETR:  9.7876E-01  1.4491E+00  3.6850E-01  7.1793E-01  7.8708E-01  9.4552E-01  6.3529E-01  2.5604E+00  1.4161E+00  4.0115E-01
             1.9041E+00
 PARAMETER:  7.8528E-02  4.7092E-01 -8.9830E-01 -2.3139E-01 -1.3942E-01  4.3980E-02 -3.5367E-01  1.0402E+00  4.4789E-01 -8.1341E-01
             7.4399E-01
 GRADIENT:   4.3219E+00 -7.4830E+00  2.4854E+00  1.2286E+00 -3.3179E+01  2.8935E-01 -5.0434E-01 -3.0484E+00 -6.5427E-01  1.5542E-01
             1.5413E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1752.26871953726        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1856
 NPARAMETR:  9.7428E-01  1.4543E+00  3.6787E-01  7.1081E-01  7.8987E-01  9.4460E-01  6.3866E-01  2.5722E+00  1.4152E+00  3.9993E-01
             1.9048E+00
 PARAMETER:  7.3945E-02  4.7456E-01 -9.0002E-01 -2.4135E-01 -1.3589E-01  4.3002E-02 -3.4839E-01  1.0447E+00  4.4730E-01 -8.1647E-01
             7.4436E-01
 GRADIENT:  -6.7007E+00 -1.1506E+01  2.8195E+00 -3.0599E+00 -3.3757E+01 -8.6945E-02 -4.2716E-01 -2.5761E+00 -1.2175E+00  1.3274E-01
             1.6336E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1752.36451150335        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2036
 NPARAMETR:  9.7686E-01  1.4570E+00  3.6732E-01  7.1111E-01  7.9120E-01  9.4506E-01  6.3776E-01  2.5758E+00  1.4360E+00  4.0339E-01
             1.9046E+00
 PARAMETER:  7.6583E-02  4.7635E-01 -9.0152E-01 -2.4093E-01 -1.3420E-01  4.3489E-02 -3.4979E-01  1.0462E+00  4.6189E-01 -8.0786E-01
             7.4426E-01
 GRADIENT:  -2.1762E-01 -1.0906E+01  2.4893E+00 -6.5340E-01 -3.3503E+01  1.6599E-01  3.4220E-01 -2.6534E+00  1.1003E+00  1.4921E-01
             2.3523E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1752.84892920233        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2215
 NPARAMETR:  9.7448E-01  1.4836E+00  3.6026E-01  6.9184E-01  7.9983E-01  9.4949E-01  6.6658E-01  2.6454E+00  1.4108E+00  3.6723E-01
             1.8933E+00
 PARAMETER:  7.4145E-02  4.9446E-01 -9.2094E-01 -2.6840E-01 -1.2336E-01  4.8169E-02 -3.0560E-01  1.0728E+00  4.4417E-01 -9.0176E-01
             7.3833E-01
 GRADIENT:  -5.8859E+00 -8.3055E+00  2.4853E+00 -5.4226E+00 -3.9584E+01  1.7727E+00  5.7179E-01  1.8635E+00 -2.5997E+00 -3.3738E-01
            -6.9005E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1753.42119271620        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2397
 NPARAMETR:  9.7706E-01  1.4892E+00  3.5568E-01  6.9796E-01  8.1818E-01  9.4487E-01  6.6158E-01  2.6293E+00  1.4318E+00  3.9312E-01
             1.9105E+00
 PARAMETER:  7.6791E-02  4.9823E-01 -9.3372E-01 -2.5960E-01 -1.0067E-01  4.3288E-02 -3.1312E-01  1.0667E+00  4.5896E-01 -8.3365E-01
             7.4735E-01
 GRADIENT:   6.7285E-02 -2.8267E+01 -6.2647E+00  9.6992E+00 -2.1021E+00  1.3909E-02  1.4542E+00 -1.8131E+00 -7.4820E-01 -1.5626E-01
             1.9398E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1755.62743370323        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2581             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7657E-01  1.5559E+00  3.8830E-01  6.5583E-01  8.6467E-01  9.4592E-01  6.0705E-01  2.7486E+00  1.5751E+00  5.2306E-01
             1.8845E+00
 PARAMETER:  7.6294E-02  5.4202E-01 -8.4597E-01 -3.2185E-01 -4.5405E-02  4.4401E-02 -3.9914E-01  1.1111E+00  5.5430E-01 -5.4806E-01
             7.3368E-01
 GRADIENT:   1.0629E+02  1.1660E+02  1.2871E+01  2.3524E+01 -1.9309E+01  8.0304E+00  5.4110E+00  1.3035E+01  1.1000E+01  6.9786E-01
             8.3278E+00

0ITERATION NO.:   84    OBJECTIVE VALUE:  -1755.66052849935        NO. OF FUNC. EVALS.: 112
 CUMULATIVE NO. OF FUNC. EVALS.:     2693
 NPARAMETR:  9.7671E-01  1.5583E+00  3.8732E-01  6.5440E-01  8.6442E-01  9.4467E-01  5.9085E-01  2.7398E+00  1.5687E+00  5.2232E-01
             1.8838E+00
 PARAMETER:  7.5430E-02  5.4208E-01 -8.4608E-01 -3.2444E-01 -4.5416E-02  4.2979E-02 -4.2664E-01  1.1110E+00  5.4555E-01 -5.4800E-01
             7.3298E-01
 GRADIENT:  -2.1457E+00 -8.5607E+02  5.4806E+02 -4.6356E-01  9.1841E+03 -6.8266E-02 -4.2028E-02  8.2756E+02 -7.2190E-01  4.4259E-02
            -4.4082E-01
 NUMSIGDIG:         1.8         2.3         2.3         2.7         2.3         2.7         2.7         2.3         1.8         2.3
                    3.1

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2693
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1464E-03 -3.8846E-02 -1.9838E-02  2.6316E-02 -2.3051E-02
 SE:             2.9598E-02  1.8315E-02  2.1706E-02  2.5321E-02  1.2640E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1534E-01  3.3922E-02  3.6075E-01  2.9866E-01  6.8210E-02

 ETASHRINKSD(%)  8.4246E-01  3.8643E+01  2.7281E+01  1.5171E+01  5.7654E+01
 ETASHRINKVR(%)  1.6778E+00  6.2353E+01  4.7120E+01  2.8041E+01  8.2068E+01
 EBVSHRINKSD(%)  1.2747E+00  3.8083E+01  3.0012E+01  1.5589E+01  5.7736E+01
 EBVSHRINKVR(%)  2.5332E+00  6.1663E+01  5.1017E+01  2.8748E+01  8.2138E+01
 RELATIVEINF(%)  9.7395E+01  5.8015E+00  2.3047E+01  1.7739E+01  4.2323E+00
 EPSSHRINKSD(%)  3.1987E+01
 EPSSHRINKVR(%)  5.3742E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1755.6605284993516     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -836.72199529467889     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    47.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1755.661       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.56E+00  3.88E-01  6.54E-01  8.65E-01  9.45E-01  5.91E-01  2.75E+00  1.56E+00  5.23E-01  1.88E+00
 


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
+        1.28E+03
 
 TH 2
+       -6.59E+01  3.28E+04
 
 TH 3
+        1.11E+02 -8.25E+04  2.12E+05
 
 TH 4
+       -1.75E+01  2.48E+02  1.49E+02  5.13E+05
 
 TH 5
+        5.80E+02  4.19E+02 -2.64E+03  1.25E+06  3.07E+06
 
 TH 6
+        4.89E+00 -6.29E+01  1.49E+02 -3.40E+00  7.42E+02  2.10E+02
 
 TH 7
+        2.83E+00 -8.48E+01  7.62E+01 -4.31E+05 -1.06E+06  1.38E-01  3.63E+05
 
 TH 8
+        1.28E+01 -9.01E+03  3.21E+02  3.06E+01 -2.76E+02  1.66E+01  1.70E+01  2.52E+03
 
 TH 9
+        3.92E+00 -7.45E+01  1.51E+02  1.28E+05  3.13E+05  1.98E+00 -1.07E+05  1.40E+01  3.56E+01
 
 TH10
+       -7.94E-01 -4.45E+02  1.11E+03 -3.78E+05 -9.26E+05  8.49E+05  3.18E+05  1.26E+02 -9.43E+04  2.79E+05
 
 TH11
+       -1.40E+01  9.64E+02 -2.53E+03  7.68E+04  1.88E+05  1.95E+00 -6.46E+04 -2.69E+02  1.03E+01 -5.67E+04  1.16E+04
 
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
 
 Elapsed finaloutput time in seconds:     1.31
 #CPUT: Total CPU Time in Seconds,       56.975
Stop Time:
Wed Sep 29 23:09:05 CDT 2021
