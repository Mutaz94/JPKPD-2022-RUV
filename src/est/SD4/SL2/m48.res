Sun Oct 24 03:16:50 CDT 2021
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
$DATA ../../../../data/SD4/SL2/dat48.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1630.32940456513        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3444E+02 -5.0686E+00  4.7270E+00  1.8439E+01  1.5971E+01  1.1424E+01 -3.2605E+00 -4.8357E+00  1.2617E+01  6.8189E+00
            -5.5694E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1637.20021963895        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0091E+00  1.0859E+00  9.6903E-01  9.8654E-01  1.0130E+00  1.0950E+00  1.0874E+00  1.0563E+00  9.1236E-01  8.8328E-01
             1.1211E+00
 PARAMETER:  1.0904E-01  1.8243E-01  6.8543E-02  8.6451E-02  1.1288E-01  1.9071E-01  1.8380E-01  1.5475E-01  8.2836E-03 -2.4114E-02
             2.1434E-01
 GRADIENT:   8.6766E+00  2.9716E-01  3.1858E+00  9.6667E-01  1.4832E+01  5.3913E+00 -3.0711E+00 -4.3545E+00 -6.9226E+00 -3.9564E+00
            -5.6897E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1637.59482293637        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0109E+00  1.0932E+00  8.8799E-01  9.7664E-01  9.6039E-01  1.0827E+00  1.1159E+00  1.0963E+00  9.3576E-01  7.9044E-01
             1.1222E+00
 PARAMETER:  1.1082E-01  1.8908E-01 -1.8791E-02  7.6361E-02  5.9584E-02  1.7947E-01  2.0969E-01  1.9195E-01  3.3605E-02 -1.3516E-01
             2.1525E-01
 GRADIENT:   1.1065E+01  2.1663E+00  1.7517E+00  1.2919E+00  2.2285E+00  7.9065E-01  1.0245E+00  3.9545E-01 -4.9558E-01 -4.6764E+00
            -4.0433E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1637.90918696922        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0052E+00  1.2213E+00  7.4317E-01  8.9081E-01  9.4831E-01  1.0810E+00  1.0215E+00  8.7025E-01  9.9451E-01  8.2758E-01
             1.1306E+00
 PARAMETER:  1.0521E-01  2.9988E-01 -1.9683E-01 -1.5629E-02  4.6922E-02  1.7790E-01  1.2123E-01 -3.8975E-02  9.4494E-02 -8.9250E-02
             2.2275E-01
 GRADIENT:  -2.3587E+00  5.6788E+00  2.3409E+00  3.6343E+00 -5.6263E+00 -2.7676E-01 -3.0609E-01 -1.6425E-01  3.7281E-03  9.6702E-01
             6.5909E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1637.99717063327        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0080E+00  1.3803E+00  5.9310E-01  7.8148E-01  9.4688E-01  1.0834E+00  9.4492E-01  6.7450E-01  1.0695E+00  7.9549E-01
             1.1288E+00
 PARAMETER:  1.0792E-01  4.2232E-01 -4.2240E-01 -1.4657E-01  4.5418E-02  1.8009E-01  4.3347E-02 -2.9378E-01  1.6718E-01 -1.2879E-01
             2.2118E-01
 GRADIENT:   9.9311E-01  7.6487E+00  1.7581E+00  4.0323E+00 -4.3698E+00  1.6139E-01  1.9371E-02 -1.6143E-03 -3.7296E-01 -3.1371E-01
            -4.7090E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1638.02329911374        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0081E+00  1.4695E+00  5.2367E-01  7.1905E-01  9.5817E-01  1.0838E+00  9.0301E-01  5.6264E-01  1.1227E+00  7.9099E-01
             1.1290E+00
 PARAMETER:  1.0812E-01  4.8495E-01 -5.4690E-01 -2.2982E-01  5.7269E-02  1.8050E-01 -2.0196E-03 -4.7511E-01  2.1572E-01 -1.3447E-01
             2.2130E-01
 GRADIENT:   1.0069E+00  4.2796E+00  1.3839E-01  3.0677E+00 -1.2144E+00  2.2026E-01 -1.2668E-01  1.5863E-01 -4.3036E-01 -4.9912E-01
            -5.2353E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1638.05320874638        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.0079E+00  1.4787E+00  5.1497E-01  7.0939E-01  9.6135E-01  1.0835E+00  8.9938E-01  4.6399E-01  1.1356E+00  8.0085E-01
             1.1302E+00
 PARAMETER:  1.0787E-01  4.9115E-01 -5.6364E-01 -2.4334E-01  6.0581E-02  1.8017E-01 -6.0484E-03 -6.6789E-01  2.2715E-01 -1.2208E-01
             2.2240E-01
 GRADIENT:   5.6261E-01 -4.9314E-01  7.1031E-01 -3.3322E-01 -4.8710E-01  1.0064E-01 -9.9974E-02  2.2242E-02 -1.0147E-01 -8.7614E-02
             2.7709E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1638.05563952786        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.0070E+00  1.4849E+00  5.0724E-01  7.0598E-01  9.6102E-01  1.0827E+00  8.9967E-01  3.8744E-01  1.1403E+00  8.0564E-01
             1.1298E+00
 PARAMETER:  1.0702E-01  4.9536E-01 -5.7876E-01 -2.4817E-01  6.0245E-02  1.7949E-01 -5.7268E-03 -8.4819E-01  2.3130E-01 -1.1612E-01
             2.2207E-01
 GRADIENT:  -1.1201E+00  4.9923E-01  7.4143E-02  1.6416E+00  8.1308E-01 -2.0252E-01  3.1811E-01  6.2801E-03  1.2241E-01  2.6780E-01
            -7.1123E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1638.05956965048        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1429
 NPARAMETR:  1.0060E+00  1.5072E+00  4.8127E-01  6.9127E-01  9.5772E-01  1.0819E+00  8.9426E-01  2.1421E-01  1.1538E+00  8.0544E-01
             1.1291E+00
 PARAMETER:  1.0602E-01  5.1026E-01 -6.3133E-01 -2.6923E-01  5.6798E-02  1.7874E-01 -1.1757E-02 -1.4408E+00  2.4306E-01 -1.1636E-01
             2.2141E-01
 GRADIENT:  -3.1187E+00  1.6866E+00 -1.5285E+00  4.5820E+00  2.1977E+00 -5.6112E-01  1.0987E+00  2.6522E-02  5.4707E-01  9.5968E-01
             5.9070E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1638.06319130365        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  1.0058E+00  1.5255E+00  4.6343E-01  6.7845E-01  9.5601E-01  1.0818E+00  8.8670E-01  1.3612E-01  1.1642E+00  7.9920E-01
             1.1288E+00
 PARAMETER:  1.0583E-01  5.2234E-01 -6.6909E-01 -2.8795E-01  5.5018E-02  1.7863E-01 -2.0246E-02 -1.8942E+00  2.5200E-01 -1.2414E-01
             2.2112E-01
 GRADIENT:  -3.4626E+00  2.0054E+00 -2.1543E+00  5.0728E+00  1.8225E+00 -6.2036E-01  1.2877E+00  2.0336E-02  6.7470E-01  1.1002E+00
             7.8058E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1638.07164118379        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1779
 NPARAMETR:  1.0067E+00  1.5488E+00  4.4560E-01  6.6091E-01  9.5724E-01  1.0825E+00  8.7405E-01  8.3752E-02  1.1783E+00  7.8916E-01
             1.1291E+00
 PARAMETER:  1.0664E-01  5.3747E-01 -7.0834E-01 -3.1414E-01  5.6300E-02  1.7932E-01 -3.4617E-02 -2.3799E+00  2.6411E-01 -1.3678E-01
             2.2146E-01
 GRADIENT:  -1.8060E+00  9.8419E-01 -1.4521E+00  2.8218E+00  2.5549E-01 -3.2183E-01  7.0903E-01  9.3090E-03  3.4577E-01  5.4893E-01
             3.5916E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1638.08359038113        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1958
 NPARAMETR:  1.0083E+00  1.5515E+00  4.4467E-01  6.5694E-01  9.5830E-01  1.0839E+00  8.6973E-01  2.9681E-02  1.1814E+00  7.8502E-01
             1.1295E+00
 PARAMETER:  1.0824E-01  5.3919E-01 -7.1043E-01 -3.2017E-01  5.7406E-02  1.8054E-01 -3.9568E-02 -3.4172E+00  2.6667E-01 -1.4204E-01
             2.2181E-01
 GRADIENT:   1.3000E+00 -8.5202E-01  1.7117E-01 -5.5844E-01 -1.2936E+00  1.9023E-01 -1.5663E-01  9.3719E-04 -1.1407E-01 -1.4751E-01
            -1.2411E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1638.08564188556        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2138
 NPARAMETR:  1.0084E+00  1.5505E+00  4.4544E-01  6.5710E-01  9.5890E-01  1.0841E+00  8.7029E-01  1.0000E-02  1.1816E+00  7.8640E-01
             1.1296E+00
 PARAMETER:  1.0838E-01  5.3855E-01 -7.0870E-01 -3.1991E-01  5.8035E-02  1.8075E-01 -3.8925E-02 -5.1119E+00  2.6690E-01 -1.4029E-01
             2.2187E-01
 GRADIENT:   1.5731E+00 -2.0414E+00 -2.3047E-02 -9.2347E-01 -7.6482E-01  2.7704E-01 -6.8727E-02  0.0000E+00 -4.6509E-02 -5.5328E-02
            -5.2999E-02

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1638.08564188556        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2160
 NPARAMETR:  1.0084E+00  1.5505E+00  4.4544E-01  6.5710E-01  9.5890E-01  1.0841E+00  8.7029E-01  1.0000E-02  1.1816E+00  7.8640E-01
             1.1296E+00
 PARAMETER:  1.0838E-01  5.3855E-01 -7.0870E-01 -3.1991E-01  5.8035E-02  1.8075E-01 -3.8925E-02 -5.1119E+00  2.6690E-01 -1.4029E-01
             2.2187E-01
 GRADIENT:   1.5731E+00 -2.0414E+00 -2.3047E-02 -9.2347E-01 -7.6482E-01  2.7704E-01 -6.8727E-02  0.0000E+00 -4.6509E-02 -5.5328E-02
            -5.2999E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2160
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0873E-04 -1.9684E-02 -2.9956E-04  1.6952E-02 -2.8110E-02
 SE:             2.9825E-02  2.4653E-02  1.2553E-04  2.3536E-02  2.0283E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9442E-01  4.2461E-01  1.7015E-02  4.7137E-01  1.6578E-01

 ETASHRINKSD(%)  8.3993E-02  1.7409E+01  9.9579E+01  2.1153E+01  3.2049E+01
 ETASHRINKVR(%)  1.6791E-01  3.1788E+01  9.9998E+01  3.7832E+01  5.3826E+01
 EBVSHRINKSD(%)  4.6884E-01  1.7374E+01  9.9627E+01  2.1786E+01  3.1902E+01
 EBVSHRINKVR(%)  9.3549E-01  3.1729E+01  9.9999E+01  3.8826E+01  5.3627E+01
 RELATIVEINF(%)  9.9031E+01  4.4701E+00  1.0883E-04  3.7404E+00  5.8817E+00
 EPSSHRINKSD(%)  4.3678E+01
 EPSSHRINKVR(%)  6.8278E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1638.0856418855572     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -902.93481532181897     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1638.086       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.55E+00  4.45E-01  6.57E-01  9.59E-01  1.08E+00  8.70E-01  1.00E-02  1.18E+00  7.86E-01  1.13E+00
 


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
 #CPUT: Total CPU Time in Seconds,       63.907
Stop Time:
Sun Oct 24 03:17:03 CDT 2021
