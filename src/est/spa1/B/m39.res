Wed Sep 29 21:05:52 CDT 2021
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
$DATA ../../../../data/spa1/B/dat39.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2147.69586474109        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4265E+02  1.2106E+01  1.2024E+01  6.9734E+00 -7.3619E+01  5.4651E+01  1.8102E+01  1.8259E+01  3.1653E+01  2.0614E+01
            -1.5948E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2156.49763053587        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.8854E-01  1.0720E+00  1.1537E+00  1.0008E+00  1.1858E+00  9.1812E-01  9.1005E-01  8.8218E-01  8.4970E-01  9.0718E-01
             1.0632E+00
 PARAMETER:  8.8472E-02  1.6954E-01  2.4300E-01  1.0084E-01  2.7046E-01  1.4571E-02  5.7448E-03 -2.5363E-02 -6.2874E-02  2.5867E-03
             1.6129E-01
 GRADIENT:  -7.3004E+00 -5.0690E+00 -1.6940E+00 -1.1657E+01  2.7318E+01 -1.6449E+01  1.1079E+00  4.5732E+00 -1.4141E+01 -2.4026E+01
             2.0983E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2160.80169035434        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.8511E-01  9.4805E-01  1.4001E+00  1.0801E+00  1.2693E+00  9.1551E-01  6.6022E-01  5.4241E-01  9.6956E-01  1.2199E+00
             1.0504E+00
 PARAMETER:  8.4998E-02  4.6656E-02  4.3652E-01  1.7709E-01  3.3849E-01  1.1725E-02 -3.1518E-01 -5.1173E-01  6.9089E-02  2.9878E-01
             1.4919E-01
 GRADIENT:  -9.1701E+00 -6.0406E+00  1.4266E+01 -8.9581E+00  8.1024E+00 -1.6734E+01  3.3354E+00 -2.1042E+00  9.2578E+00 -4.3931E+00
             1.0248E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2162.81625583725        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.9101E-01  1.0716E+00  1.1917E+00  1.0019E+00  1.2164E+00  9.5508E-01  5.1654E-01  5.3969E-01  1.0240E+00  1.1831E+00
             1.0305E+00
 PARAMETER:  9.0967E-02  1.6914E-01  2.7537E-01  1.0192E-01  2.9593E-01  5.4043E-02 -5.6060E-01 -5.1677E-01  1.2375E-01  2.6810E-01
             1.3000E-01
 GRADIENT:   1.8214E+00  5.3503E+00  2.3530E-01  6.9453E+00 -1.0583E+00 -2.0086E-01 -1.6409E-01  1.0522E-01 -7.4533E-01  1.3477E+00
             1.2698E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2162.87433701096        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.9139E-01  1.1531E+00  1.0817E+00  9.4546E-01  1.2019E+00  9.5621E-01  5.9987E-01  3.9029E-01  1.0483E+00  1.1450E+00
             1.0289E+00
 PARAMETER:  9.1349E-02  2.4245E-01  1.7853E-01  4.3915E-02  2.8392E-01  5.5222E-02 -4.1105E-01 -8.4087E-01  1.4721E-01  2.3541E-01
             1.2848E-01
 GRADIENT:   3.1754E-01 -3.1140E-02 -4.4039E-01  6.2323E-01 -2.7481E-01 -8.5461E-02  2.0250E-01  1.2568E-01  3.7263E-01  3.9240E-01
             6.0492E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2162.87838864233        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.9199E-01  1.2183E+00  1.0261E+00  9.0295E-01  1.2083E+00  9.5635E-01  6.1524E-01  3.0660E-01  1.0763E+00  1.1375E+00
             1.0278E+00
 PARAMETER:  9.1953E-02  2.9743E-01  1.2575E-01 -2.0926E-03  2.8918E-01  5.5370E-02 -3.8574E-01 -1.0822E+00  1.7355E-01  2.2881E-01
             1.2746E-01
 GRADIENT:   2.8764E-01  6.5623E-01  1.8322E-01  2.6393E-01 -1.0884E+00 -2.6205E-01  1.0335E-01  1.0203E-01  2.6195E-01  4.0365E-01
             1.1113E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2162.89299532000        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.9244E-01  1.2298E+00  1.0119E+00  8.9557E-01  1.2097E+00  9.5731E-01  6.2031E-01  2.1172E-01  1.0783E+00  1.1377E+00
             1.0283E+00
 PARAMETER:  9.2406E-02  3.0688E-01  1.1187E-01 -1.0290E-02  2.9035E-01  5.6370E-02 -3.7754E-01 -1.4525E+00  1.7541E-01  2.2904E-01
             1.2794E-01
 GRADIENT:   1.1110E+00  5.9185E-01  2.0261E-01  7.4573E-01  7.5400E-01  9.4198E-02 -1.1759E-01  1.4300E-02 -2.4375E-01  2.2461E-01
             1.2670E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2162.91708166464        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1236
 NPARAMETR:  9.9187E-01  1.2339E+00  9.8094E-01  8.9218E-01  1.1923E+00  9.5708E-01  6.4962E-01  6.4893E-02  1.0692E+00  1.1147E+00
             1.0276E+00
 PARAMETER:  9.1839E-02  3.1015E-01  8.0761E-02 -1.4082E-02  2.7589E-01  5.6130E-02 -3.3137E-01 -2.6350E+00  1.6695E-01  2.0859E-01
             1.2725E-01
 GRADIENT:  -7.5283E-01 -5.2795E-01 -2.0596E-01 -4.1775E-02  3.2870E-02 -6.7667E-02 -3.1249E-02  4.6852E-03 -8.3143E-02 -2.0115E-02
            -3.1371E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2162.91723756360        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1412
 NPARAMETR:  9.9163E-01  1.2350E+00  9.7994E-01  8.9139E-01  1.1920E+00  9.5693E-01  6.5016E-01  4.7769E-02  1.0704E+00  1.1151E+00
             1.0278E+00
 PARAMETER:  9.1595E-02  3.1107E-01  7.9736E-02 -1.4972E-02  2.7563E-01  5.5975E-02 -3.3054E-01 -2.9414E+00  1.6807E-01  2.0891E-01
             1.2742E-01
 GRADIENT:  -1.3680E+00 -4.2339E-01  7.0214E-03 -1.5160E-01 -4.8036E-01 -1.3288E-01  3.4467E-02  2.5924E-03  1.0169E-01  1.1473E-01
             1.3427E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2162.91734776994        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1587
 NPARAMETR:  9.9149E-01  1.2358E+00  9.8007E-01  8.9076E-01  1.1929E+00  9.5686E-01  6.4894E-01  3.4602E-02  1.0710E+00  1.1157E+00
             1.0278E+00
 PARAMETER:  9.1449E-02  3.1175E-01  7.9865E-02 -1.5682E-02  2.7639E-01  5.5905E-02 -3.3242E-01 -3.2639E+00  1.6864E-01  2.0947E-01
             1.2737E-01
 GRADIENT:  -1.7207E+00 -5.8863E-01 -1.6671E-02 -2.7903E-01 -2.2206E-01 -1.6250E-01 -4.3349E-03  1.2961E-03  2.2919E-02  5.4054E-02
             6.1882E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2162.91739289402        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1762
 NPARAMETR:  9.9143E-01  1.2364E+00  9.8012E-01  8.9037E-01  1.1935E+00  9.5685E-01  6.4824E-01  2.3149E-02  1.0713E+00  1.1160E+00
             1.0277E+00
 PARAMETER:  9.1389E-02  3.1218E-01  7.9922E-02 -1.6124E-02  2.7688E-01  5.5893E-02 -3.3349E-01 -3.6658E+00  1.6890E-01  2.0975E-01
             1.2733E-01
 GRADIENT:  -1.8679E+00 -6.8273E-01 -4.2577E-02 -3.3812E-01 -2.4082E-02 -1.6843E-01 -3.3170E-02  5.6972E-04 -4.3902E-02  4.4931E-03
             3.0037E-04

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2162.91897364135        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1943
 NPARAMETR:  9.9201E-01  1.2401E+00  9.7863E-01  8.8824E-01  1.1949E+00  9.5794E-01  6.4889E-01  1.0000E-02  1.0730E+00  1.1162E+00
             1.0277E+00
 PARAMETER:  9.1975E-02  3.1522E-01  7.8400E-02 -1.8514E-02  2.7806E-01  5.7027E-02 -3.3249E-01 -1.2545E+01  1.7045E-01  2.0990E-01
             1.2731E-01
 GRADIENT:  -5.0983E-01 -2.2266E-01  8.4471E-02  4.8473E-02 -5.2468E-02  2.7012E-01 -9.1109E-03  0.0000E+00 -2.5794E-02 -2.6469E-02
            -1.8335E-02

0ITERATION NO.:   57    OBJECTIVE VALUE:  -2162.91935073922        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     2007
 NPARAMETR:  9.9343E-01  1.2418E+00  9.7628E-01  8.8622E-01  1.1953E+00  9.5712E-01  6.5057E-01  1.0000E-02  1.0741E+00  1.1157E+00
             1.0276E+00
 PARAMETER:  9.2404E-02  3.1739E-01  7.6278E-02 -2.0276E-02  2.7833E-01  5.6940E-02 -3.3062E-01 -1.8407E+01  1.7131E-01  2.0939E-01
             1.2729E-01
 GRADIENT:  -1.0716E+00  5.5004E-01  6.2549E-02  3.3730E-01 -5.9466E-02  1.3370E-01 -9.5106E-03  0.0000E+00 -1.7398E-02 -5.9133E-03
             1.4190E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2007
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.7936E-04 -2.7499E-02 -3.2786E-04  1.0677E-02 -3.1543E-02
 SE:             2.9848E-02  1.5683E-02  1.4956E-04  2.6285E-02  2.3945E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7382E-01  7.9535E-02  2.8366E-02  6.8459E-01  1.8774E-01

 ETASHRINKSD(%)  5.7050E-03  4.7459E+01  9.9499E+01  1.1942E+01  1.9780E+01
 ETASHRINKVR(%)  1.1410E-02  7.2394E+01  9.9997E+01  2.2457E+01  3.5647E+01
 EBVSHRINKSD(%)  3.7161E-01  4.6810E+01  9.9539E+01  1.2287E+01  1.7055E+01
 EBVSHRINKVR(%)  7.4184E-01  7.1708E+01  9.9998E+01  2.3065E+01  3.1202E+01
 RELATIVEINF(%)  9.8840E+01  1.1498E+00  3.2495E-04  3.9345E+00  1.2706E+01
 EPSSHRINKSD(%)  3.1759E+01
 EPSSHRINKVR(%)  5.3431E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2162.9193507392206     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1243.9808175345479     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2162.919       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.92E-01  1.24E+00  9.77E-01  8.87E-01  1.20E+00  9.58E-01  6.50E-01  1.00E-02  1.07E+00  1.12E+00  1.03E+00
 


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
+        1.22E+03
 
 TH 2
+       -1.09E+01  5.06E+02
 
 TH 3
+        6.20E+00  1.62E+02  2.61E+02
 
 TH 4
+       -5.94E+00  5.11E+02 -9.69E+01  9.44E+02
 
 TH 5
+        3.15E+00 -2.21E+02 -2.59E+02  7.79E+01  4.34E+02
 
 TH 6
+        1.79E+00 -1.94E+00  1.49E+00 -2.10E+00 -2.29E-01  2.14E+02
 
 TH 7
+        1.13E+00 -1.72E+01  1.32E+01 -1.68E+01 -1.52E+01  1.77E-01  3.55E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.75E-01 -3.41E+01 -1.01E+01  3.62E+01  1.89E+00 -1.31E-01  3.90E+01  0.00E+00  1.03E+02
 
 TH10
+        1.58E+00 -3.43E+00 -2.92E+01 -5.34E-01 -4.66E+01  1.62E-01  1.25E+01  0.00E+00  1.55E+00  7.47E+01
 
 TH11
+       -8.82E+00 -2.77E+01 -4.17E+01 -7.57E+00  9.54E+00  1.36E+00  7.43E+00  0.00E+00  5.96E+00  2.22E+01  4.03E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.812
Stop Time:
Wed Sep 29 21:06:31 CDT 2021
