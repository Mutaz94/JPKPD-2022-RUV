Sun Oct 24 02:02:15 CDT 2021
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
$DATA ../../../../data/SD4/A1/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1315.04758950943        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9612E+02  3.1965E+01 -3.0175E+01  9.4319E+01  1.0234E+02  6.5851E+01 -8.5371E+00 -2.2702E-01 -3.5414E+01 -1.8140E+01
            -6.5771E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1497.32027970952        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0594E+00  9.8599E-01  1.0695E+00  1.0207E+00  9.5679E-01  9.9654E-01  9.9760E-01  9.4680E-01  1.1236E+00  9.5353E-01
             2.0088E+00
 PARAMETER:  1.5770E-01  8.5887E-02  1.6718E-01  1.2048E-01  5.5831E-02  9.6533E-02  9.7601E-02  4.5330E-02  2.1656E-01  5.2419E-02
             7.9756E-01
 GRADIENT:   2.3846E+02  1.3747E+01 -4.4999E-01  2.6893E+01  3.8232E-01  3.4997E+01  4.6680E+00  2.7251E+00  1.1557E+01 -2.8308E+00
             1.5686E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1499.47109690394        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0428E+00  9.1127E-01  9.3496E-01  1.0572E+00  8.7283E-01  9.2690E-01  8.4641E-01  5.5960E-01  1.1104E+00  9.8767E-01
             1.9701E+00
 PARAMETER:  1.4187E-01  7.0802E-03  3.2745E-02  1.5559E-01 -3.6009E-02  2.4093E-02 -6.6748E-02 -4.8054E-01  2.0469E-01  8.7595E-02
             7.7810E-01
 GRADIENT:   1.9339E+02  1.1652E+01 -6.3661E+00  3.9496E+01  2.6043E+00  1.1221E+01 -1.1536E+00  6.7656E-01  9.0196E+00  5.9405E+00
             4.7090E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1500.44492070957        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      294
 NPARAMETR:  1.0315E+00  8.2322E-01  1.1432E+00  1.1198E+00  9.2279E-01  9.2306E-01  9.2598E-01  8.4885E-01  1.0359E+00  9.6485E-01
             1.9949E+00
 PARAMETER:  1.3099E-01 -9.4528E-02  2.3386E-01  2.1311E-01  1.9648E-02  1.9938E-02  2.3101E-02 -6.3873E-02  1.3528E-01  6.4219E-02
             7.9062E-01
 GRADIENT:   8.6408E+00  5.8546E+00  1.4517E+00  2.7151E+00 -3.6160E+00  1.2400E+00 -2.6288E-01  6.2799E-01 -7.2166E-01 -2.6728E+00
             2.5503E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1501.70245666095        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  1.0241E+00  4.8225E-01  1.0749E+00  1.3250E+00  7.8251E-01  9.2386E-01  1.0818E+00  5.9613E-01  9.1584E-01  9.3626E-01
             1.9799E+00
 PARAMETER:  1.2386E-01 -6.2929E-01  1.7220E-01  3.8145E-01 -1.4525E-01  2.0803E-02  1.7860E-01 -4.1729E-01  1.2091E-02  3.4137E-02
             7.8306E-01
 GRADIENT:  -1.8506E+00  6.1345E+00  3.2746E+00  1.3242E+01 -8.9992E+00  2.3698E+00 -2.1583E-01 -4.5663E-02  4.8016E-01  1.5584E+00
            -3.2951E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1502.46980291692        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  1.0186E+00  2.5515E-01  1.0856E+00  1.4533E+00  7.3278E-01  9.1677E-01  1.2648E+00  5.7081E-01  8.3894E-01  9.0567E-01
             2.0137E+00
 PARAMETER:  1.1847E-01 -1.2659E+00  1.8215E-01  4.7385E-01 -2.1091E-01  1.3099E-02  3.3489E-01 -4.6070E-01 -7.5616E-02  9.1612E-04
             7.9998E-01
 GRADIENT:  -6.2277E+00  1.6060E+00 -1.6201E+00  7.2381E+00  1.5503E+00  1.0457E+00 -8.7936E-02  6.5954E-02 -7.3064E-01  1.2402E-01
             1.5027E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1502.77236978228        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  1.0200E+00  1.0290E-01  1.0597E+00  1.5362E+00  6.8479E-01  9.1085E-01  1.6151E+00  4.7253E-01  7.9574E-01  8.8773E-01
             2.0215E+00
 PARAMETER:  1.1978E-01 -2.1740E+00  1.5795E-01  5.2931E-01 -2.7865E-01  6.6256E-03  5.7939E-01 -6.4965E-01 -1.2848E-01 -1.9083E-02
             8.0383E-01
 GRADIENT:   5.0378E+00  7.1959E-01  7.8081E-01  8.1117E+00 -2.4018E+00 -3.7147E-01 -1.1518E-03 -8.9511E-02 -4.9731E-01 -4.0549E-01
            -1.0403E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1502.84540236405        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  1.0160E+00  4.8711E-02  1.0839E+00  1.5656E+00  6.8325E-01  9.1072E-01  1.9279E+00  4.7034E-01  7.8023E-01  8.9299E-01
             2.0327E+00
 PARAMETER:  1.1586E-01 -2.9219E+00  1.8053E-01  5.4828E-01 -2.8090E-01  6.4778E-03  7.5643E-01 -6.5431E-01 -1.4817E-01 -1.3180E-02
             8.0936E-01
 GRADIENT:  -1.7331E+00  1.6884E-01  9.7652E-01  3.0040E-01 -1.2591E+00  1.2341E-01  4.7451E-03 -9.7023E-02  2.5744E-02  4.5728E-02
             1.0016E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1502.85412698247        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  1.0156E+00  2.9195E-02  1.0920E+00  1.5788E+00  6.8204E-01  9.1025E-01  2.1685E+00  4.8981E-01  7.7398E-01  8.9248E-01
             2.0305E+00
 PARAMETER:  1.1550E-01 -3.4338E+00  1.8798E-01  5.5666E-01 -2.8266E-01  5.9624E-03  8.7402E-01 -6.1374E-01 -1.5620E-01 -1.3748E-02
             8.0829E-01
 GRADIENT:  -1.5184E+00  1.2920E-01  1.1071E+00  2.8287E+00 -1.7846E+00  5.4192E-02  1.6397E-03 -1.2200E-01 -2.2555E-01  6.3387E-02
             3.9742E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1502.87240484484        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1356
 NPARAMETR:  1.0156E+00  1.0536E-02  1.1177E+00  1.5923E+00  6.8802E-01  9.0975E-01  2.6898E+00  6.1142E-01  7.6894E-01  8.8769E-01
             2.0191E+00
 PARAMETER:  1.1544E-01 -4.4530E+00  2.1124E-01  5.6515E-01 -2.7394E-01  5.4187E-03  1.0895E+00 -3.9197E-01 -1.6274E-01 -1.9138E-02
             8.0265E-01
 GRADIENT:  -1.6094E-01  3.9290E-02  2.4814E-01  1.9980E+00 -8.1439E-01 -3.8320E-02 -1.2597E-04  1.3066E-02 -1.9438E-01  1.1514E-01
             3.4622E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1502.88031976697        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1545             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0157E+00  1.0000E-02  1.1182E+00  1.5893E+00  6.8832E-01  9.0982E-01  2.8257E+00  6.1113E-01  7.6924E-01  8.8731E-01
             2.0185E+00
 PARAMETER:  1.1559E-01 -4.5763E+00  2.1171E-01  5.6327E-01 -2.7349E-01  5.4960E-03  1.1387E+00 -3.9245E-01 -1.6235E-01 -1.9557E-02
             8.0237E-01
 GRADIENT:   1.2118E+02  0.0000E+00  1.9554E+00  2.6234E+02  8.6148E+00  8.9739E+00  8.8880E-03  1.0509E-01  5.1107E+00  8.2199E-01
             6.3971E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1502.88070189325        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1733
 NPARAMETR:  1.0157E+00  1.0000E-02  1.1167E+00  1.5890E+00  6.8735E-01  9.0981E-01  2.8264E+00  6.0625E-01  7.6933E-01  8.8704E-01
             2.0184E+00
 PARAMETER:  1.1556E-01 -4.5763E+00  2.1040E-01  5.6312E-01 -2.7492E-01  5.4847E-03  1.1390E+00 -4.0047E-01 -1.6224E-01 -1.9864E-02
             8.0231E-01
 GRADIENT:   4.9123E-01  0.0000E+00  8.5091E-01 -4.9243E+00 -6.5908E-01  2.8255E-02  1.2623E-04 -1.7624E-02  4.9339E-02  6.2392E-03
            -1.7098E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1502.88142739178        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1931             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0157E+00  1.0000E-02  1.1144E+00  1.5888E+00  6.8718E-01  9.0981E-01  2.7637E+00  6.0843E-01  7.6930E-01  8.8701E-01
             2.0192E+00
 PARAMETER:  1.1556E-01 -4.5763E+00  2.0828E-01  5.6300E-01 -2.7517E-01  5.4838E-03  1.1166E+00 -3.9688E-01 -1.6227E-01 -1.9905E-02
             8.0268E-01
 GRADIENT:   1.2094E+02  0.0000E+00  1.4500E+00  2.6186E+02  9.2497E+00  8.9595E+00  8.3455E-03  1.3268E-01  5.0883E+00  8.9059E-01
             6.5926E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1502.88184351619        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     2102             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0157E+00  1.0000E-02  1.1136E+00  1.5887E+00  6.8667E-01  9.0982E-01  2.6271E+00  6.0686E-01  7.6935E-01  8.8705E-01
             2.0193E+00
 PARAMETER:  1.1558E-01 -4.5763E+00  2.0757E-01  5.6294E-01 -2.7590E-01  5.4891E-03  1.0659E+00 -3.9945E-01 -1.6221E-01 -1.9858E-02
             8.0276E-01
 GRADIENT:   1.2097E+02  0.0000E+00  1.6707E+00  2.6185E+02  8.8701E+00  8.9593E+00  7.2052E-03  1.3362E-01  5.0964E+00  9.3314E-01
             6.6176E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1502.88214973383        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2291
 NPARAMETR:  1.0157E+00  1.0000E-02  1.1128E+00  1.5886E+00  6.8636E-01  9.0981E-01  2.6660E+00  6.0528E-01  7.6936E-01  8.8681E-01
             2.0193E+00
 PARAMETER:  1.1556E-01 -4.5763E+00  2.0689E-01  5.6285E-01 -2.7635E-01  5.4835E-03  1.0806E+00 -4.0206E-01 -1.6220E-01 -2.0121E-02
             8.0276E-01
 GRADIENT:   4.6065E-01  0.0000E+00  8.7051E-02 -5.0008E+00  2.9806E-01  2.9596E-02  5.9552E-05  2.4724E-02  2.0856E-02  1.0633E-01
             1.5948E-01

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1502.88228545227        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     2389
 NPARAMETR:  1.0157E+00  1.0000E-02  1.1127E+00  1.5886E+00  6.8628E-01  9.0981E-01  2.6676E+00  6.0203E-01  7.6936E-01  8.8625E-01
             2.0190E+00
 PARAMETER:  1.1557E-01 -4.5763E+00  2.0681E-01  5.6284E-01 -2.7646E-01  5.4787E-03  1.0812E+00 -4.0745E-01 -1.6220E-01 -2.0755E-02
             8.0259E-01
 GRADIENT:   4.8826E-01  0.0000E+00  2.0613E-01 -4.9525E+00  2.5539E-01  2.4386E-02  6.7927E-05 -4.7087E-03  1.6918E-02 -2.6430E-02
            -8.0194E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2389
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.2661E-04 -7.3193E-04 -2.0295E-03 -9.0799E-03 -2.2880E-02
 SE:             2.9373E-02  4.4088E-04  8.8798E-03  2.7757E-02  2.0532E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8841E-01  9.6887E-02  8.1922E-01  7.4357E-01  2.6512E-01

 ETASHRINKSD(%)  1.5980E+00  9.8523E+01  7.0252E+01  7.0120E+00  3.1214E+01
 ETASHRINKVR(%)  3.1705E+00  9.9978E+01  9.1150E+01  1.3532E+01  5.2685E+01
 EBVSHRINKSD(%)  1.7855E+00  9.8594E+01  7.0157E+01  6.8066E+00  3.0657E+01
 EBVSHRINKVR(%)  3.5391E+00  9.9980E+01  9.1094E+01  1.3150E+01  5.1916E+01
 RELATIVEINF(%)  9.0347E+01  7.3429E-04  6.3154E-01  4.8324E+00  2.0643E+00
 EPSSHRINKSD(%)  3.5320E+01
 EPSSHRINKVR(%)  5.8165E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1502.8822854522741     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -767.73145888853594     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1502.882       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E-02  1.11E+00  1.59E+00  6.86E-01  9.10E-01  2.67E+00  6.02E-01  7.69E-01  8.86E-01  2.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,       70.047
Stop Time:
Sun Oct 24 02:02:29 CDT 2021
