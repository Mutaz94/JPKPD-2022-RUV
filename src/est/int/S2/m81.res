Sat Sep 18 04:42:24 CDT 2021
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
$DATA ../../../../data/int/S2/dat81.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3474.10429581060        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.0184E+00 -1.3156E+01  1.0265E+02 -5.1533E+01  9.6193E+01 -4.1234E+01 -5.8430E+01 -5.0343E+02 -1.0535E+02 -4.6440E+01
            -1.3513E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3712.20900297813        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0379E+00  9.6421E-01  9.3981E-01  1.0362E+00  9.5476E-01  1.0930E+00  1.2224E+00  2.6196E+00  9.8935E-01  1.1233E+00
             9.8896E-01
 PARAMETER:  1.3715E-01  6.3553E-02  3.7928E-02  1.3553E-01  5.3703E-02  1.8897E-01  3.0083E-01  1.0630E+00  8.9295E-02  2.1629E-01
             8.8899E-02
 GRADIENT:   1.0029E+02 -6.6642E+00 -1.4859E+01  2.1309E+01  8.9210E-01  2.9774E+00 -1.1036E+01 -9.3764E+00  1.7035E+01 -7.3630E+00
            -6.5770E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3714.60414493559        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0284E+00  1.0156E+00  1.0488E+00  1.0106E+00  1.0218E+00  1.0483E+00  1.1419E+00  2.7338E+00  1.0084E+00  1.2650E+00
             1.0017E+00
 PARAMETER:  1.2796E-01  1.1546E-01  1.4769E-01  1.1057E-01  1.2161E-01  1.4718E-01  2.3270E-01  1.1057E+00  1.0839E-01  3.3510E-01
             1.0168E-01
 GRADIENT:   7.8869E+01 -6.9906E+00  2.0885E+00  1.5481E+01  1.7266E+01 -1.7046E+01 -2.4895E+01 -6.1517E+00  2.3197E+01  1.2665E+01
            -4.3644E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3715.50740681602        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  1.0230E+00  1.0148E+00  1.0521E+00  1.0097E+00  1.0206E+00  1.0461E+00  1.1858E+00  2.7325E+00  1.0114E+00  1.2494E+00
             1.0046E+00
 PARAMETER:  1.2274E-01  1.1464E-01  1.5077E-01  1.0969E-01  1.2035E-01  1.4507E-01  2.7044E-01  1.1052E+00  1.1135E-01  3.2267E-01
             1.0459E-01
 GRADIENT:   9.1315E+00 -1.5115E+01  1.3195E+00  2.8604E+00  7.0613E+00 -2.7372E+01 -1.8853E+01 -1.3804E+01  2.3264E+01 -2.2067E+00
            -3.6722E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3716.39149225475        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  1.0230E+00  1.0147E+00  1.0520E+00  1.0098E+00  1.0206E+00  1.1169E+00  1.1857E+00  2.7332E+00  1.0114E+00  1.2607E+00
             1.0046E+00
 PARAMETER:  1.2271E-01  1.1462E-01  1.5073E-01  1.0972E-01  1.2038E-01  2.1059E-01  2.7037E-01  1.1055E+00  1.1132E-01  3.3166E-01
             1.0461E-01
 GRADIENT:   8.0834E+00 -1.5805E+01  1.3849E+00  3.0491E+00  6.4388E+00  2.5956E-01 -1.8758E+01 -1.3886E+01  2.3446E+01 -2.6230E-01
            -3.6082E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3716.66895670080        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      578
 NPARAMETR:  1.0214E+00  1.0147E+00  1.0515E+00  1.0093E+00  1.0200E+00  1.1142E+00  1.1911E+00  2.7333E+00  1.0050E+00  1.2589E+00
             1.0057E+00
 PARAMETER:  1.2116E-01  1.1459E-01  1.5022E-01  1.0927E-01  1.1985E-01  2.0815E-01  2.7487E-01  1.1055E+00  1.0497E-01  3.3021E-01
             1.0565E-01
 GRADIENT:   6.0986E+01 -6.2594E+00  2.2718E+00  1.3054E+01  1.3782E+01  1.2900E+01 -1.5443E+01 -6.7099E+00  2.3004E+01  1.2283E+01
            -3.3723E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3717.66124559552        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      655
 NPARAMETR:  1.0139E+00  1.0139E+00  1.0493E+00  1.0082E+00  1.0163E+00  1.1053E+00  1.2160E+00  2.7333E+00  9.7505E-01  1.2428E+00
             1.0103E+00
 PARAMETER:  1.1385E-01  1.1385E-01  1.4814E-01  1.0820E-01  1.1621E-01  2.0007E-01  2.9560E-01  1.1055E+00  7.4737E-02  3.1739E-01
             1.1026E-01
 GRADIENT:   4.3880E+01 -4.4657E+00  1.6483E+00  9.1495E+00  1.0069E+01  9.2479E+00 -1.1155E+01 -6.8610E+00  1.6667E+01  8.8582E+00
            -2.4367E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3717.83831027251        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      789
 NPARAMETR:  1.0128E+00  1.0139E+00  1.0493E+00  1.0082E+00  1.0163E+00  1.1042E+00  1.2160E+00  2.7334E+00  9.6441E-01  1.2412E+00
             1.0118E+00
 PARAMETER:  1.1272E-01  1.1385E-01  1.4814E-01  1.0820E-01  1.1621E-01  1.9913E-01  2.9560E-01  1.1055E+00  6.3765E-02  3.1606E-01
             1.1173E-01
 GRADIENT:   4.1118E+01 -4.3804E+00  1.5426E+00  8.7647E+00  9.9926E+00  8.7871E+00 -1.1303E+01 -6.9774E+00  1.4353E+01  8.2518E+00
            -2.1704E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3718.10379571380        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      939
 NPARAMETR:  1.0141E+00  1.0139E+00  1.0493E+00  1.0082E+00  1.0163E+00  1.1072E+00  1.2160E+00  2.7334E+00  9.5124E-01  1.2504E+00
             1.0145E+00
 PARAMETER:  1.1403E-01  1.1385E-01  1.4814E-01  1.0820E-01  1.1621E-01  2.0179E-01  2.9560E-01  1.1055E+00  5.0016E-02  3.2348E-01
             1.1442E-01
 GRADIENT:  -7.6215E+00 -1.3991E+01  4.3713E-01 -2.0991E+00  1.5648E+00 -3.2315E+00 -1.3974E+01 -1.4451E+01  1.0960E+01 -2.5052E+00
            -1.6527E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3718.12734786758        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1137             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0142E+00  1.0139E+00  1.0493E+00  1.0082E+00  1.0163E+00  1.1074E+00  1.2160E+00  2.7334E+00  9.4974E-01  1.2510E+00
             1.0148E+00
 PARAMETER:  1.1414E-01  1.1385E-01  1.4814E-01  1.0820E-01  1.1621E-01  2.0201E-01  2.9560E-01  1.1055E+00  4.8433E-02  3.2397E-01
             1.1465E-01
 GRADIENT:   4.4038E+01 -5.0392E+00  1.4374E+00  8.5739E+00  9.1328E+00  1.0053E+01 -1.1381E+01 -7.3080E+00  1.1406E+01  9.8613E+00
            -1.5881E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3718.14351265571        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1334            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0143E+00  1.0140E+00  1.0494E+00  1.0082E+00  1.0163E+00  1.1164E+00  1.2164E+00  2.7326E+00  9.4956E-01  1.2512E+00
             1.0148E+00
 PARAMETER:  1.1418E-01  1.1393E-01  1.4817E-01  1.0818E-01  1.1618E-01  2.1008E-01  2.9589E-01  1.1053E+00  4.8245E-02  3.2409E-01
             1.1465E-01
 GRADIENT:   4.4261E+01 -4.9590E+00  1.4584E+00  8.5568E+00  9.0501E+00  1.3800E+01 -1.1314E+01 -7.4081E+00  1.1368E+01  9.8918E+00
            -1.5863E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3718.15036793468        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1502
 NPARAMETR:  1.0143E+00  1.0140E+00  1.0493E+00  1.0082E+00  1.0163E+00  1.1163E+00  1.2165E+00  2.7326E+00  9.4914E-01  1.2512E+00
             1.0148E+00
 PARAMETER:  1.1416E-01  1.1395E-01  1.4817E-01  1.0818E-01  1.1617E-01  2.1004E-01  2.9599E-01  1.1053E+00  4.7805E-02  3.2408E-01
             1.1468E-01
 GRADIENT:  -7.2434E+00 -1.3921E+01  4.4294E-01 -2.1386E+00  1.4037E+00  6.1431E-02 -1.3909E+01 -1.4583E+01  1.0519E+01 -2.4245E+00
            -1.6013E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3718.33535098877        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1682
 NPARAMETR:  1.0163E+00  1.0141E+00  1.0494E+00  1.0082E+00  1.0163E+00  1.1171E+00  1.2165E+00  2.7325E+00  9.3563E-01  1.2590E+00
             1.0171E+00
 PARAMETER:  1.1618E-01  1.1395E-01  1.4817E-01  1.0817E-01  1.1617E-01  2.1097E-01  2.9599E-01  1.1053E+00  3.3461E-02  3.3033E-01
             1.1694E-01
 GRADIENT:  -1.0579E+05 -1.0787E+05 -1.6589E+05  2.2722E+05  1.0579E+05  4.3787E-01 -4.1538E+04  2.2186E+04 -2.4580E+05 -3.7208E+04
             2.1018E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.0         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1682
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1364E-03 -3.0616E-02 -3.9034E-02  2.4115E-02 -5.6410E-02
 SE:             2.9880E-02  2.4451E-02  2.6465E-02  2.5861E-02  2.4190E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4300E-01  2.1051E-01  1.4024E-01  3.5110E-01  1.9703E-02

 ETASHRINKSD(%)  1.0000E-10  1.8087E+01  1.1339E+01  1.3361E+01  1.8960E+01
 ETASHRINKVR(%)  1.0000E-10  3.2903E+01  2.1392E+01  2.4937E+01  3.4325E+01
 EBVSHRINKSD(%)  2.1919E-01  2.1124E+01  1.5204E+01  1.3015E+01  1.9469E+01
 EBVSHRINKVR(%)  4.3789E-01  3.7786E+01  2.8097E+01  2.4336E+01  3.5148E+01
 RELATIVEINF(%)  9.9560E+01  3.7409E+01  6.4397E+01  5.2494E+01  4.0348E+01
 EPSSHRINKSD(%)  2.3316E+01
 EPSSHRINKVR(%)  4.1195E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3718.3353509887661     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2064.2459912203553     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    52.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3718.335       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.01E+00  1.05E+00  1.01E+00  1.02E+00  1.12E+00  1.22E+00  2.73E+00  9.36E-01  1.26E+00  1.02E+00
 


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
+        4.41E+08
 
 TH 2
+       -1.07E+03  4.60E+08
 
 TH 3
+       -7.78E+02  1.75E+04  2.54E+08
 
 TH 4
+        1.12E+03 -2.47E+04  2.27E+04  5.17E+08
 
 TH 5
+        1.03E+03 -2.32E+04  2.67E+04 -2.98E+04  4.41E+08
 
 TH 6
+       -1.18E+03 -1.21E+03 -8.97E+02  1.28E+03  1.19E+03  1.57E+02
 
 TH 7
+       -3.40E+02  7.58E+03 -8.78E+03  9.84E+03 -1.65E+04 -3.88E+02  4.74E+07
 
 TH 8
+        4.04E+01 -9.01E+02 -1.31E+07  1.87E+07  1.73E+07  4.59E+01 -5.66E+06  1.35E+06
 
 TH 9
+        3.48E+03  3.54E+03  2.65E+03 -3.73E+03 -3.46E+03 -1.49E+03  1.15E+03 -1.32E+02  7.02E+08
 
 TH10
+        3.66E+03  3.71E+03  2.78E+03 -3.95E+03 -3.68E+03 -3.33E+02  1.20E+03 -1.46E+02  1.58E+08  3.55E+07
 
 TH11
+        1.02E+03  1.66E+04  1.23E+04 -1.76E+04 -1.63E+04  1.17E+03  5.36E+03 -6.33E+02 -3.43E+03 -3.62E+03  4.34E+08
 
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
 #CPUT: Total CPU Time in Seconds,       68.195
Stop Time:
Sat Sep 18 04:43:34 CDT 2021
