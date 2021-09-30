Wed Sep 29 09:34:01 CDT 2021
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
$DATA ../../../../data/int/D/dat70.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m70.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26495.0879883422        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9700E+02  4.8143E+02 -5.0626E+01  1.5714E+02  3.5044E+02 -2.7769E+03 -1.2765E+03 -6.7155E+01 -2.2025E+03 -9.7554E+02
            -5.2964E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -962.677088629594        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  2.2206E+00  1.6534E+00  9.0175E-01  3.3883E+00  8.0777E-01  4.9431E+00  5.3737E+00  9.7785E-01  6.2196E+00  2.8988E+00
             1.0110E+01
 PARAMETER:  8.9778E-01  6.0285E-01 -3.4214E-03  1.3203E+00 -1.1348E-01  1.6980E+00  1.7815E+00  7.7603E-02  1.9277E+00  1.1643E+00
             2.4135E+00
 GRADIENT:   9.5132E+01  2.0493E+01 -3.9492E+01  7.2556E+01 -4.1264E+01  1.8330E+02  1.0754E+02  4.7691E+00  1.3987E+02  5.3405E+01
             3.5277E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1108.63541177033        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.5598E+00  9.7718E-01  4.0869E+01  7.7593E+00  2.1060E+00  4.0678E+00  7.7586E+00  5.1147E-01  6.8867E+00  1.6465E+00
             1.0193E+01
 PARAMETER:  5.4456E-01  7.6911E-02  3.8104E+00  2.1489E+00  8.4479E-01  1.5031E+00  2.1488E+00 -5.7046E-01  2.0296E+00  5.9862E-01
             2.4217E+00
 GRADIENT:   4.6023E+01  3.6537E+00 -1.1710E+01  1.4177E+02 -3.3698E+01  1.5813E+02  2.1860E+01 -4.7529E-03  2.2536E+02  4.3575E+01
             4.0840E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1274.24392043589        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      238
 NPARAMETR:  1.4323E+00  1.7195E-01  1.4750E+02  1.6925E+00  2.6928E+00  2.7345E+00  1.0161E+01  5.9466E-02  1.7648E+00  8.6839E-01
             9.7964E+00
 PARAMETER:  4.5930E-01 -1.6605E+00  5.0938E+00  6.2622E-01  1.0906E+00  1.1059E+00  2.4185E+00 -2.7223E+00  6.6801E-01 -4.1113E-02
             2.3820E+00
 GRADIENT:   4.0612E+01 -1.1835E+00  1.3383E-01 -4.3781E+00  2.7186E+01  2.1086E+01  3.8281E+02  1.0354E-05  3.0797E+01  1.3205E+01
             3.8581E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1311.08975350317        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      417
 NPARAMETR:  1.0845E+00  1.3595E-01  3.2909E+02  1.6432E+00  2.5329E+00  3.2958E+00  1.1288E+01  2.7526E-01  1.5919E+00  7.9521E-01
             8.0094E+00
 PARAMETER:  1.8112E-01 -1.8954E+00  5.8963E+00  5.9665E-01  1.0294E+00  1.2926E+00  2.5238E+00 -1.1900E+00  5.6494E-01 -1.2915E-01
             2.1806E+00
 GRADIENT:  -4.5204E+01 -1.2123E+01  2.8808E-01  7.5365E+00  7.1176E+00 -1.0437E+00 -1.5293E+01  1.6254E-05  6.9269E+00  1.1155E+01
             7.9818E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1330.62464596437        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      592
 NPARAMETR:  1.3164E+00  4.2128E-01  1.4031E+01  1.4977E+00  2.4170E+00  3.3022E+00  1.0259E+01  1.4598E-01  1.6592E+00  3.2784E-01
             7.8636E+00
 PARAMETER:  3.7492E-01 -7.6446E-01  2.7413E+00  5.0395E-01  9.8251E-01  1.2946E+00  2.4282E+00 -1.8243E+00  6.0631E-01 -1.0152E+00
             2.1622E+00
 GRADIENT:   7.8312E-01  1.1413E-01 -7.6087E+00 -1.5363E+01  5.2077E+01  2.3317E+01  3.7767E+01 -1.1843E-03  1.7126E+01  1.9127E+00
             2.9372E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1340.87703053335        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      768
 NPARAMETR:  1.3095E+00  6.5520E-01  6.7465E+00  1.2581E+00  1.9643E+00  3.1292E+00  6.8475E+00  1.1031E+00  1.4035E+00  2.4440E-01
             7.7527E+00
 PARAMETER:  3.6965E-01 -3.2282E-01  2.0090E+00  3.2961E-01  7.7516E-01  1.2408E+00  2.0239E+00  1.9814E-01  4.3895E-01 -1.3090E+00
             2.1480E+00
 GRADIENT:   1.0470E+00 -9.6476E+00 -2.0468E-01 -5.8560E+00 -1.4741E+01  5.0958E+00 -1.0927E+01 -1.1287E-01  7.3601E+00  1.0327E+00
            -1.6558E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1343.82438490302        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      945
 NPARAMETR:  1.2695E+00  9.8799E-01  7.5161E+00  1.1834E+00  2.1369E+00  3.1364E+00  6.5041E+00  2.1781E+00  1.3169E+00  1.7264E-01
             7.7847E+00
 PARAMETER:  3.3866E-01  8.7922E-02  2.1170E+00  2.6842E-01  8.5936E-01  1.2431E+00  1.9724E+00  8.7843E-01  3.7529E-01 -1.6566E+00
             2.1522E+00
 GRADIENT:  -6.4760E+00 -7.3892E-01 -3.8360E+00  1.5503E+00  9.9752E+00  5.5582E+00  4.7898E+00  8.7237E-01  6.2888E+00  6.1114E-01
             4.8315E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1344.18130061187        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:     1096
 NPARAMETR:  1.2685E+00  9.8822E-01  7.6762E+00  1.1861E+00  2.1189E+00  3.0953E+00  6.4937E+00  2.1818E+00  1.3125E+00  3.5231E-02
             7.7420E+00
 PARAMETER:  3.3780E-01  8.8154E-02  2.1381E+00  2.7064E-01  8.5089E-01  1.2299E+00  1.9708E+00  8.8013E-01  3.7192E-01 -3.2458E+00
             2.1467E+00
 GRADIENT:   2.5389E+01  1.8723E+00 -1.4819E+00  1.2834E+01  7.2964E+00  9.6792E+01  1.9987E+02  5.5827E-01  6.2922E+00  2.7267E-02
             3.0957E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1344.19513262865        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  1.2677E+00  9.8803E-01  7.7238E+00  1.1866E+00  2.1146E+00  3.0828E+00  6.4489E+00  2.1853E+00  1.3113E+00  2.4650E-02
             7.7611E+00
 PARAMETER:  3.3720E-01  8.7962E-02  2.1443E+00  2.7110E-01  8.4885E-01  1.2258E+00  1.9639E+00  8.8176E-01  3.7101E-01 -3.6030E+00
             2.1491E+00
 GRADIENT:  -6.8225E+00 -3.5888E-01 -1.6113E+00  3.4659E+00  9.9933E-01 -1.1238E+00  1.5163E+00  4.1245E-01  4.9227E+00  1.2200E-02
            -3.4284E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1344.20554818732        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1432
 NPARAMETR:  1.2677E+00  9.8786E-01  7.8246E+00  1.1869E+00  2.1112E+00  3.0772E+00  6.4006E+00  2.1882E+00  1.3090E+00  1.5670E-02
             7.7865E+00
 PARAMETER:  3.3718E-01  8.7788E-02  2.1573E+00  2.7137E-01  8.4725E-01  1.2240E+00  1.9564E+00  8.8310E-01  3.6929E-01 -4.0560E+00
             2.1524E+00
 GRADIENT:  -7.0456E+00 -6.0096E-01 -9.1840E-01  3.6836E+00 -1.6011E+00 -1.6758E+00 -3.3842E-01  2.7525E-01  4.7508E+00  4.9017E-03
             2.5849E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1344.40865627249        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1612
 NPARAMETR:  1.2743E+00  9.8799E-01  8.2946E+00  1.1848E+00  2.1262E+00  3.1273E+00  6.3829E+00  2.1832E+00  1.2990E+00  1.0000E-02
             7.7772E+00
 PARAMETER:  3.4243E-01  8.7918E-02  2.2156E+00  2.6955E-01  8.5435E-01  1.2402E+00  1.9536E+00  8.8081E-01  3.6156E-01 -5.0407E+00
             2.1512E+00
 GRADIENT:  -5.4873E+00 -6.1357E-01  4.1278E-01  2.9382E+00 -4.1348E+00  4.6004E+00 -1.4688E+00 -2.3880E-01  4.1269E+00  0.0000E+00
            -1.8873E-01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1344.40913911842        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1676
 NPARAMETR:  1.2740E+00  9.8793E-01  8.2889E+00  1.1850E+00  2.1257E+00  3.1228E+00  6.3970E+00  2.2024E+00  1.2989E+00  1.0000E-02
             7.7761E+00
 PARAMETER:  3.4237E-01  8.7929E-02  2.2134E+00  2.6955E-01  8.5466E-01  1.2395E+00  1.9545E+00  8.8076E-01  3.6177E-01 -4.9903E+00
             2.1510E+00
 GRADIENT:   2.2393E+02  7.9767E+02 -3.5767E+01 -2.9318E+02  9.2185E+01  5.1212E+01 -4.1873E+01 -2.0070E-01  1.1445E+02  0.0000E+00
            -2.2612E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1676
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4228E-02  3.6159E-02 -1.3423E-02 -7.0283E-02  4.4198E-05
 SE:             2.9543E-02  2.4531E-02  6.1302E-03  1.4461E-02  1.3583E-04
 N:                     100         100         100         100         100

 P VAL.:         6.3009E-01  1.4048E-01  2.8548E-02  1.1747E-06  7.4489E-01

 ETASHRINKSD(%)  1.0257E+00  1.7818E+01  7.9463E+01  5.1553E+01  9.9545E+01
 ETASHRINKVR(%)  2.0409E+00  3.2461E+01  9.5782E+01  7.6529E+01  9.9998E+01
 EBVSHRINKSD(%)  1.5120E+00  1.2437E+01  8.1076E+01  5.3811E+01  9.9459E+01
 EBVSHRINKVR(%)  3.0011E+00  2.3327E+01  9.6419E+01  7.8666E+01  9.9997E+01
 RELATIVEINF(%)  9.6916E+01  3.9802E+01  9.7287E-01  1.0945E+01  8.3424E-04
 EPSSHRINKSD(%)  7.8149E+00
 EPSSHRINKVR(%)  1.5019E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1344.4091391184184     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       309.68022064999241     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    55.00
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1344.409       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.27E+00  9.88E-01  8.28E+00  1.18E+00  2.13E+00  3.13E+00  6.39E+00  2.18E+00  1.30E+00  1.00E-02  7.78E+00
 


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
+        1.03E+04
 
 TH 2
+        4.56E+04  2.04E+05
 
 TH 3
+        7.13E+00  8.47E-02  6.18E+00
 
 TH 4
+        4.04E+02  1.12E+01 -8.57E-01  1.97E+04
 
 TH 5
+       -7.23E+01 -4.20E+00 -2.96E+00 -1.75E+02  7.00E+02
 
 TH 6
+       -1.38E+03  9.07E-01 -1.40E-02  5.73E-01  1.44E+01  1.11E+02
 
 TH 7
+        1.04E+01  2.61E+00 -9.03E-02 -8.58E+00 -3.23E+00 -3.42E+01  1.52E+01
 
 TH 8
+       -2.35E+03 -2.36E-01 -2.49E-01  9.94E-01  1.41E+00  1.54E-02  5.07E-02  6.82E-01
 
 TH 9
+       -2.76E+02  5.01E+00 -4.69E-01 -2.06E+01  1.28E+02 -1.56E-01  1.63E+00 -4.66E-02  9.05E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.17E+02 -1.91E+00 -4.11E-02 -9.44E+00  5.06E+01 -1.88E+01  5.22E-01  2.06E-01  4.07E+00  0.00E+00  2.11E+01
 
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
 #CPUT: Total CPU Time in Seconds,       72.741
Stop Time:
Wed Sep 29 09:35:15 CDT 2021
