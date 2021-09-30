Wed Sep 29 16:50:07 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat76.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1560.89269650926        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9875E+02  4.1596E+01  8.1334E+00  7.6317E+01  6.1496E+01  6.0688E-03  6.9805E+00 -1.0954E+01  4.2641E+00 -2.7782E+01
            -1.0431E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1563.78036117860        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.8953E-01  9.3372E-01  8.4130E-01  1.0390E+00  8.8868E-01  9.7151E-01  9.2882E-01  1.0633E+00  9.7184E-01  1.0845E+00
             1.4759E+00
 PARAMETER:  8.9472E-02  3.1421E-02 -7.2807E-02  1.3827E-01 -1.8015E-02  7.1092E-02  2.6156E-02  1.6133E-01  7.1436E-02  1.8116E-01
             4.8924E-01
 GRADIENT:   2.2225E+02 -1.4546E+00 -3.7479E+01  6.3498E+01  3.3921E+01 -3.0766E+01  2.6803E+00  9.2878E+00  6.4979E+00  1.9880E+01
             7.4166E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1566.08951797940        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.8036E-01  9.8064E-01  6.0249E-01  1.0025E+00  7.1698E-01  1.0136E+00  7.5808E-01  9.4864E-01  1.0312E+00  8.4099E-01
             1.3986E+00
 PARAMETER:  8.0162E-02  8.0448E-02 -4.0669E-01  1.0255E-01 -2.3270E-01  1.1354E-01 -1.7696E-01  4.7279E-02  1.3071E-01 -7.3170E-02
             4.3547E-01
 GRADIENT:   2.1486E+02  4.2952E+01 -3.7159E+00  7.8288E+01 -1.7071E+01 -7.4291E+00 -4.6307E+00  9.3132E+00  1.1455E+01  9.5203E+00
             5.5195E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1573.25940730513        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      342
 NPARAMETR:  9.5902E-01  8.9731E-01  7.2969E-01  1.0412E+00  7.6819E-01  1.0528E+00  1.0622E+00  9.1920E-01  8.9481E-01  8.2530E-01
             1.3225E+00
 PARAMETER:  5.8153E-02 -8.3575E-03 -2.1513E-01  1.4038E-01 -1.6372E-01  1.5145E-01  1.6035E-01  1.5751E-02 -1.1141E-02 -9.2012E-02
             3.7951E-01
 GRADIENT:  -4.1694E+01  3.3964E+00 -4.2223E+00 -4.6894E+00 -8.2731E+00 -1.8703E+01  1.9189E+00  4.7366E+00 -9.1227E+00 -6.8389E-01
             3.1100E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1577.08164241292        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.7840E-01  6.9049E-01  9.5763E-01  1.1822E+00  8.1926E-01  1.0992E+00  1.1294E+00  9.6804E-01  8.9955E-01  9.6281E-01
             1.2140E+00
 PARAMETER:  7.8165E-02 -2.7035E-01  5.6706E-02  2.6736E-01 -9.9348E-02  1.9462E-01  2.2166E-01  6.7518E-02 -5.8654E-03  6.2103E-02
             2.9390E-01
 GRADIENT:   7.2830E+00 -2.4846E+00 -1.0296E+00 -7.2557E+00  6.0680E-01 -1.6231E-01  1.0642E+00  5.6994E-01  3.5036E+00 -8.1193E-01
            -1.2803E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1577.17951310275        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  9.7366E-01  6.6622E-01  1.0378E+00  1.2067E+00  8.4939E-01  1.0999E+00  1.0709E+00  1.0318E+00  8.8732E-01  1.0086E+00
             1.2124E+00
 PARAMETER:  7.3310E-02 -3.0614E-01  1.3708E-01  2.8791E-01 -6.3237E-02  1.9525E-01  1.6847E-01  1.3129E-01 -1.9549E-02  1.0859E-01
             2.9262E-01
 GRADIENT:   4.4663E-02  1.4641E+00  1.2013E+00  6.3433E-01 -1.1231E+00  4.0558E-01 -1.6886E-01 -2.7655E-01  3.2874E-01 -1.0296E-01
            -1.0485E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1577.18314786847        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  9.7400E-01  6.3602E-01  1.0575E+00  1.2260E+00  8.4893E-01  1.0988E+00  1.0547E+00  1.0451E+00  8.8122E-01  1.0125E+00
             1.2115E+00
 PARAMETER:  7.3657E-02 -3.5252E-01  1.5590E-01  3.0373E-01 -6.3781E-02  1.9424E-01  1.5326E-01  1.4408E-01 -2.6450E-02  1.1245E-01
             2.9183E-01
 GRADIENT:   1.5403E+00  9.6045E-01  5.9223E-01  7.9431E-01 -7.6547E-02  1.4368E-01 -5.6655E-01 -3.2769E-01  3.3771E-01 -3.8184E-01
            -1.4288E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1577.18364468731        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1044
 NPARAMETR:  9.7367E-01  5.9124E-01  1.0872E+00  1.2550E+00  8.4773E-01  1.0976E+00  1.0559E+00  1.0667E+00  8.6878E-01  1.0183E+00
             1.2113E+00
 PARAMETER:  7.3319E-02 -4.2553E-01  1.8360E-01  3.2717E-01 -6.5191E-02  1.9310E-01  1.5437E-01  1.6455E-01 -4.0660E-02  1.1816E-01
             2.9169E-01
 GRADIENT:   2.2430E+00  1.0386E+00  4.6852E-01  1.6576E+00  4.1544E-02 -8.1769E-02 -7.5086E-01 -3.5790E-01  2.5216E-01 -4.8034E-01
            -1.5203E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1577.18469145571        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  9.7320E-01  5.5802E-01  1.1082E+00  1.2764E+00  8.4638E-01  1.0968E+00  1.0663E+00  1.0830E+00  8.5869E-01  1.0219E+00
             1.2115E+00
 PARAMETER:  7.2838E-02 -4.8336E-01  2.0270E-01  3.4407E-01 -6.6789E-02  1.9242E-01  1.6416E-01  1.7972E-01 -5.2346E-02  1.2169E-01
             2.9188E-01
 GRADIENT:   2.3625E+00  1.0976E+00  4.5206E-01  2.0713E+00 -2.8066E-02 -1.7358E-01 -7.8466E-01 -3.5179E-01  1.6922E-01 -4.8486E-01
            -1.4451E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1577.18537061534        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1404
 NPARAMETR:  9.7275E-01  5.3053E-01  1.1248E+00  1.2941E+00  8.4499E-01  1.0962E+00  1.0799E+00  1.0967E+00  8.5006E-01  1.0244E+00
             1.2118E+00
 PARAMETER:  7.2377E-02 -5.3387E-01  2.1759E-01  3.5780E-01 -6.8429E-02  1.9188E-01  1.7689E-01  1.9228E-01 -6.2453E-02  1.2415E-01
             2.9215E-01
 GRADIENT:   2.3564E+00  1.1390E+00  4.4448E-01  2.3607E+00 -1.1579E-01 -2.3556E-01 -7.8254E-01 -3.3646E-01  9.3180E-02 -4.6873E-01
            -1.3364E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1577.18571183780        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1584
 NPARAMETR:  9.7239E-01  5.0978E-01  1.1370E+00  1.3073E+00  8.4385E-01  1.0958E+00  1.0935E+00  1.1073E+00  8.4345E-01  1.0261E+00
             1.2121E+00
 PARAMETER:  7.1998E-02 -5.7377E-01  2.2841E-01  3.6798E-01 -6.9776E-02  1.9151E-01  1.8938E-01  2.0188E-01 -7.0258E-02  1.2578E-01
             2.9236E-01
 GRADIENT:   2.3075E+00  1.1410E+00  4.2813E-01  2.4844E+00 -1.6225E-01 -2.6864E-01 -7.6725E-01 -3.2026E-01  4.4181E-02 -4.5011E-01
            -1.2450E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1577.18583751152        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1764
 NPARAMETR:  9.7211E-01  4.9475E-01  1.1458E+00  1.3169E+00  8.4301E-01  1.0956E+00  1.1053E+00  1.1151E+00  8.3863E-01  1.0272E+00
             1.2123E+00
 PARAMETER:  7.1709E-02 -6.0370E-01  2.3607E-01  3.7527E-01 -7.0780E-02  1.9126E-01  2.0015E-01  2.0895E-01 -7.5987E-02  1.2686E-01
             2.9251E-01
 GRADIENT:   2.2521E+00  1.1261E+00  4.1231E-01  2.5159E+00 -1.8261E-01 -2.8258E-01 -7.4976E-01 -3.0742E-01  1.7198E-02 -4.3449E-01
            -1.1808E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1577.23891354942        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1946
 NPARAMETR:  9.7042E-01  4.8928E-01  1.1466E+00  1.3189E+00  8.4253E-01  1.0960E+00  1.2282E+00  1.1182E+00  8.3039E-01  1.0260E+00
             1.2151E+00
 PARAMETER:  6.9973E-02 -6.1481E-01  2.3680E-01  3.7682E-01 -7.1349E-02  1.9166E-01  3.0553E-01  2.1174E-01 -8.5856E-02  1.2569E-01
             2.9485E-01
 GRADIENT:  -7.5827E-01  2.7931E-01 -2.2183E-01 -1.3114E+00  1.1699E+00 -4.2605E-02 -1.0353E-01 -2.5230E-02  4.9911E-01 -9.5043E-02
            -1.4100E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1577.24538271430        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2125             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7125E-01  4.7676E-01  1.1479E+00  1.3251E+00  8.3831E-01  1.0967E+00  1.2705E+00  1.1186E+00  8.2295E-01  1.0233E+00
             1.2152E+00
 PARAMETER:  7.0827E-02 -6.4074E-01  2.3791E-01  3.8151E-01 -7.6372E-02  1.9232E-01  3.3943E-01  2.1210E-01 -9.4862E-02  1.2303E-01
             2.9488E-01
 GRADIENT:   2.9364E+02  4.1870E+01  4.9188E+00  3.3697E+02  5.8619E+00  6.1413E+01  2.5914E+00  5.0899E-01  7.0246E+00  7.8892E-01
             2.3339E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1577.24603390153        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2306
 NPARAMETR:  9.7123E-01  4.7682E-01  1.1467E+00  1.3254E+00  8.3788E-01  1.0967E+00  1.2693E+00  1.1175E+00  8.2336E-01  1.0230E+00
             1.2150E+00
 PARAMETER:  7.0813E-02 -6.4061E-01  2.3686E-01  3.8168E-01 -7.6885E-02  1.9234E-01  3.3845E-01  2.1109E-01 -9.4365E-02  1.2270E-01
             2.9477E-01
 GRADIENT:   1.1185E+00  8.4285E-02  4.4537E-01 -3.6170E+00  1.0149E-01  2.9456E-01 -2.3274E-02 -5.2767E-02 -1.0970E-01 -1.8278E-02
            -6.1190E-02

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1577.24625778019        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     2401
 NPARAMETR:  9.7123E-01  4.7671E-01  1.1457E+00  1.3254E+00  8.3780E-01  1.0967E+00  1.2745E+00  1.1183E+00  8.2372E-01  1.0231E+00
             1.2152E+00
 PARAMETER:  7.0810E-02 -6.4085E-01  2.3598E-01  3.8173E-01 -7.6977E-02  1.9233E-01  3.4255E-01  2.1179E-01 -9.3923E-02  1.2284E-01
             2.9491E-01
 GRADIENT:  -3.0375E-03 -5.1001E-02 -2.2027E-02 -6.4965E-02  5.8040E-01 -2.9622E-03  1.4024E-02  5.5317E-02  1.1471E-01  6.9041E-02
             6.1525E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2401
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3406E-04 -8.6337E-03 -2.9857E-02 -3.0812E-03 -3.4575E-02
 SE:             2.9781E-02  1.0216E-02  1.6528E-02  2.7351E-02  2.1384E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8837E-01  3.9805E-01  7.0847E-02  9.1030E-01  1.0591E-01

 ETASHRINKSD(%)  2.3151E-01  6.5775E+01  4.4630E+01  8.3714E+00  2.8361E+01
 ETASHRINKVR(%)  4.6248E-01  8.8286E+01  6.9341E+01  1.6042E+01  4.8678E+01
 EBVSHRINKSD(%)  5.5482E-01  6.6034E+01  4.7647E+01  8.4831E+00  2.5930E+01
 EBVSHRINKVR(%)  1.1066E+00  8.8463E+01  7.2591E+01  1.6247E+01  4.5136E+01
 RELATIVEINF(%)  9.6152E+01  2.9728E-01  5.3174E+00  2.6707E+00  5.9784E+00
 EPSSHRINKSD(%)  4.4027E+01
 EPSSHRINKVR(%)  6.8670E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1577.2462577801900     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -842.09543121645186     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.20
 Elapsed covariance  time in seconds:     5.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1577.246       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.71E-01  4.77E-01  1.15E+00  1.33E+00  8.38E-01  1.10E+00  1.27E+00  1.12E+00  8.24E-01  1.02E+00  1.22E+00
 


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
 
         3.37E-02  5.49E-01  5.47E-01  3.59E-01  1.43E-01  8.56E-02  8.25E-01  5.54E-01  1.84E-01  2.09E-01  8.28E-02
 


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
+        1.13E-03
 
 TH 2
+        5.07E-03  3.02E-01
 
 TH 3
+       -6.22E-03 -2.53E-01  3.00E-01
 
 TH 4
+       -3.59E-03 -1.96E-01  1.70E-01  1.29E-01
 
 TH 5
+       -1.43E-03 -2.18E-02  5.76E-02  1.69E-02  2.03E-02
 
 TH 6
+       -4.87E-04  1.66E-02 -1.58E-02 -1.05E-02 -1.48E-03  7.33E-03
 
 TH 7
+       -6.79E-03 -4.32E-01  3.61E-01  2.79E-01  3.16E-02 -2.19E-02  6.81E-01
 
 TH 8
+       -5.42E-03 -2.45E-01  2.85E-01  1.65E-01  5.30E-02 -1.56E-02  3.54E-01  3.07E-01
 
 TH 9
+        1.57E-03  9.51E-02 -8.22E-02 -6.25E-02 -7.80E-03  5.58E-03 -1.33E-01 -8.13E-02  3.40E-02
 
 TH10
+       -2.19E-03 -6.15E-02  9.13E-02  4.25E-02  2.31E-02 -3.72E-03  9.06E-02  8.15E-02 -2.11E-02  4.37E-02
 
 TH11
+       -2.59E-04 -1.08E-02  1.27E-02  6.73E-03  2.24E-03 -6.47E-04  1.74E-02  1.15E-02 -2.22E-03  3.68E-03  6.85E-03
 
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
+        3.37E-02
 
 TH 2
+        2.74E-01  5.49E-01
 
 TH 3
+       -3.38E-01 -8.40E-01  5.47E-01
 
 TH 4
+       -2.97E-01 -9.94E-01  8.64E-01  3.59E-01
 
 TH 5
+       -2.98E-01 -2.79E-01  7.38E-01  3.29E-01  1.43E-01
 
 TH 6
+       -1.69E-01  3.53E-01 -3.38E-01 -3.41E-01 -1.22E-01  8.56E-02
 
 TH 7
+       -2.44E-01 -9.54E-01  8.00E-01  9.42E-01  2.69E-01 -3.11E-01  8.25E-01
 
 TH 8
+       -2.91E-01 -8.05E-01  9.39E-01  8.27E-01  6.70E-01 -3.30E-01  7.73E-01  5.54E-01
 
 TH 9
+        2.53E-01  9.39E-01 -8.14E-01 -9.43E-01 -2.97E-01  3.53E-01 -8.76E-01 -7.95E-01  1.84E-01
 
 TH10
+       -3.11E-01 -5.36E-01  7.98E-01  5.66E-01  7.74E-01 -2.08E-01  5.26E-01  7.03E-01 -5.48E-01  2.09E-01
 
 TH11
+       -9.29E-02 -2.37E-01  2.81E-01  2.26E-01  1.90E-01 -9.14E-02  2.56E-01  2.51E-01 -1.46E-01  2.13E-01  8.28E-02
 
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
+        1.15E+03
 
 TH 2
+        1.31E+00  5.09E+02
 
 TH 3
+       -1.44E+01  1.16E+02  2.26E+02
 
 TH 4
+        1.01E+02  5.62E+02 -2.61E+01  9.08E+02
 
 TH 5
+        9.05E+01 -3.07E+02 -3.90E+02 -8.03E+01  9.51E+02
 
 TH 6
+        1.32E+02 -2.58E+01  3.76E+01 -5.15E+01 -6.14E+01  1.85E+02
 
 TH 7
+       -1.25E+01  4.35E+01  5.49E+00  1.86E+01 -5.95E+00 -6.61E+00  1.82E+01
 
 TH 8
+       -9.77E+00 -7.67E+00 -3.00E+01 -5.63E-01 -4.46E+00  1.19E-01 -2.64E+00  3.00E+01
 
 TH 9
+        2.80E+01 -2.51E+01  1.32E+01  1.04E+02 -7.79E+01 -8.15E+00 -1.19E+01  1.59E+01  3.08E+02
 
 TH10
+        1.23E+01 -2.78E+00 -2.14E+01  1.29E+01 -6.14E+01 -7.79E-02 -4.86E+00  1.14E+01  2.65E+01  8.60E+01
 
 TH11
+        7.27E+00  2.88E+01 -2.82E+01  4.72E+01  3.67E+01 -3.23E+00 -2.35E+00  5.98E-01 -5.04E+01 -1.09E+00  1.75E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       36.201
Stop Time:
Wed Sep 29 16:50:45 CDT 2021
