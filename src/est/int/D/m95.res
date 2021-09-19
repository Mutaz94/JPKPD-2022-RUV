Sat Sep 18 07:47:38 CDT 2021
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
$DATA ../../../../data/int/D/dat95.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   51023.0518911110        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.8046E+02  8.1078E+02  3.3392E+01  6.6182E+02  8.7294E+01 -4.2103E+03 -2.1939E+03 -8.7291E+01 -2.8765E+03 -1.1549E+03
            -9.8430E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -509.219453580942        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.2190E+00  1.8567E+00  8.7025E-01  2.1692E+00  1.0044E+00  5.0033E+00  3.7568E+00  9.8266E-01  2.1886E+00  1.2639E+00
             1.2618E+01
 PARAMETER:  2.9805E-01  7.1883E-01 -3.8976E-02  8.7434E-01  1.0443E-01  1.7101E+00  1.4236E+00  8.2504E-02  8.8328E-01  3.3418E-01
             2.6351E+00
 GRADIENT:  -7.4027E+00  3.4191E+01 -3.7786E+01  1.8506E+02 -1.4307E+01  1.0962E+02 -9.4310E+01  4.8596E+00 -4.9369E+01  2.4426E+01
            -1.6616E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -556.053405977075        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.0575E-01  3.1013E+00  1.0426E+02  1.0084E+01  3.4280E+00  1.9934E+00  8.9434E+00  6.8836E-01  7.4434E+00  8.4501E-01
             1.2973E+01
 PARAMETER:  1.0105E-03  1.2318E+00  4.7469E+00  2.4110E+00  1.3320E+00  7.8983E-01  2.2909E+00 -2.7344E-01  2.1073E+00 -6.8411E-02
             2.6629E+00
 GRADIENT:  -5.8249E+01  2.1964E+01 -5.2791E+00  7.3494E+01  2.6615E+01  1.4082E+01  1.5839E+01  4.9799E-02  1.0125E+01  7.0585E+00
             6.2626E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -644.755140849356        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1706E+00  1.8227E+00  6.7529E+00  2.0736E+00  2.1222E+00  2.9289E+00  4.1753E+00  2.2248E+00  3.3039E+00  7.6850E-01
             1.2900E+01
 PARAMETER:  2.5754E-01  7.0031E-01  2.0100E+00  8.2928E-01  8.5243E-01  1.1746E+00  1.5292E+00  8.9968E-01  1.2951E+00 -1.6331E-01
             2.6572E+00
 GRADIENT:  -1.4553E+01  1.0558E+01 -1.5651E+01  8.2699E+01 -1.3052E+01  2.8535E+01 -8.3976E+01  3.0757E+00  1.5777E+01  7.4259E+00
             9.3000E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -715.906477860968        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.1903E+00  1.0082E+00  4.7161E+01  1.5584E+00  2.6505E+00  2.1685E+00  7.6416E+00  1.8444E+00  1.3698E+00  4.4056E-01
             1.2721E+01
 PARAMETER:  2.7419E-01  1.0817E-01  3.9536E+00  5.4368E-01  1.0748E+00  8.7405E-01  2.1336E+00  7.1218E-01  4.1466E-01 -7.1970E-01
             2.6432E+00
 GRADIENT:   3.5347E+01  7.4465E+00  1.0733E-01  1.3149E+01 -1.1030E+01  1.9397E+01  1.9469E+01  1.8441E-02  1.3780E+00  1.8521E+00
             3.7369E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -721.172445985852        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0643E+00  8.6766E-01  5.1322E+01  1.4391E+00  2.7273E+00  2.0076E+00  6.8597E+00  2.7131E+00  1.2965E+00  3.5086E-01
             1.2425E+01
 PARAMETER:  1.6230E-01 -4.1954E-02  4.0381E+00  4.6403E-01  1.1033E+00  7.9692E-01  2.0257E+00  1.0981E+00  3.5970E-01 -9.4738E-01
             2.6197E+00
 GRADIENT:  -1.3753E+00 -3.1077E+00 -1.6462E-01 -2.0513E+00 -4.0356E-01  3.9684E-01 -4.0706E-01  3.3845E-02  1.8027E+00  1.2534E+00
             3.2013E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -721.565263145392        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.0496E+00  9.5566E-01  4.8339E+01  1.3955E+00  2.7493E+00  2.0219E+00  6.7690E+00  2.5808E+00  1.2773E+00  2.5536E-01
             1.2086E+01
 PARAMETER:  1.4844E-01  5.4646E-02  3.9782E+00  4.3327E-01  1.1113E+00  8.0404E-01  2.0124E+00  1.0481E+00  3.4472E-01 -1.2651E+00
             2.5921E+00
 GRADIENT:  -1.7832E+00 -5.9267E-01 -2.2056E-01  3.2693E-01  1.2785E+00  1.0351E+00  2.1133E+00  2.8338E-02  6.5930E-01  6.5322E-01
            -1.0970E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -721.575529720554        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      514
 NPARAMETR:  1.0514E+00  9.9598E-01  4.6699E+01  1.3825E+00  2.7420E+00  2.0175E+00  6.6949E+00  2.5443E+00  1.2662E+00  2.1875E-01
             1.2120E+01
 PARAMETER:  1.5016E-01  9.5970E-02  3.9437E+00  4.2389E-01  1.1087E+00  8.0186E-01  2.0014E+00  1.0339E+00  3.3602E-01 -1.4198E+00
             2.5949E+00
 GRADIENT:  -1.6235E+00 -6.7669E-03 -1.9827E-01  3.6914E-01  5.8927E-01  4.4372E-01  1.7847E+00  2.8064E-02  6.0964E-01  4.7662E-01
            -7.1838E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -721.576360534115        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  1.0527E+00  1.0075E+00  4.6434E+01  1.3771E+00  2.7411E+00  2.0162E+00  6.6657E+00  2.5248E+00  1.2562E+00  2.0149E-01
             1.2135E+01
 PARAMETER:  1.5132E-01  1.0744E-01  3.9380E+00  4.1996E-01  1.1084E+00  8.0119E-01  1.9970E+00  1.0262E+00  3.2811E-01 -1.5020E+00
             2.5961E+00
 GRADIENT:  -1.3627E+00  5.8519E-02 -1.9009E-01  3.9474E-01  4.5794E-01  3.0179E-01  1.4153E+00  2.7753E-02  3.7749E-01  4.0349E-01
            -5.7049E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -721.576603807401        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  1.0535E+00  1.0150E+00  4.6398E+01  1.3733E+00  2.7410E+00  2.0153E+00  6.6456E+00  2.5064E+00  1.2486E+00  1.8761E-01
             1.2145E+01
 PARAMETER:  1.5216E-01  1.1488E-01  3.9373E+00  4.1722E-01  1.1083E+00  8.0077E-01  1.9940E+00  1.0189E+00  3.2200E-01 -1.5734E+00
             2.5969E+00
 GRADIENT:  -1.1693E+00  8.7754E-02 -1.8344E-01  4.1846E-01  3.7874E-01  2.2113E-01  1.1349E+00  2.7306E-02  1.8888E-01  3.4924E-01
            -4.6642E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -721.576775678488        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.0542E+00  1.0201E+00  4.6473E+01  1.3707E+00  2.7411E+00  2.0147E+00  6.6319E+00  2.4902E+00  1.2433E+00  1.7688E-01
             1.2152E+01
 PARAMETER:  1.5274E-01  1.1987E-01  3.9389E+00  4.1535E-01  1.1083E+00  8.0048E-01  1.9919E+00  1.0123E+00  3.1779E-01 -1.6323E+00
             2.5975E+00
 GRADIENT:  -1.0303E+00  1.0455E-01 -1.7807E-01  4.3176E-01  3.2540E-01  1.6762E-01  9.3990E-01  2.6822E-02  6.0309E-02  3.1006E-01
            -3.9379E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -721.576940882019        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      802
 NPARAMETR:  1.0549E+00  1.0255E+00  4.6678E+01  1.3680E+00  2.7414E+00  2.0141E+00  6.6171E+00  2.4679E+00  1.2377E+00  1.6384E-01
             1.2160E+01
 PARAMETER:  1.5340E-01  1.2523E-01  3.9433E+00  4.1334E-01  1.1085E+00  8.0017E-01  1.9897E+00  1.0034E+00  3.1325E-01 -1.7089E+00
             2.5982E+00
 GRADIENT:  -8.7428E-01  1.2035E-01 -1.7122E-01  4.4163E-01  2.6704E-01  1.0977E-01  7.2463E-01  2.6080E-02 -7.6458E-02  2.6563E-01
            -3.1287E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -721.577155870705        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      872
 NPARAMETR:  1.0556E+00  1.0315E+00  4.7107E+01  1.3650E+00  2.7422E+00  2.0134E+00  6.6010E+00  2.4363E+00  1.2317E+00  1.4771E-01
             1.2169E+01
 PARAMETER:  1.5414E-01  1.3097E-01  3.9524E+00  4.1118E-01  1.1088E+00  7.9983E-01  1.9872E+00  9.9050E-01  3.0837E-01 -1.8125E+00
             2.5989E+00
 GRADIENT:  -6.9555E-01  1.3515E-01 -1.6219E-01  4.4614E-01  2.0282E-01  4.6733E-02  4.8808E-01  2.4933E-02 -2.2063E-01  2.1552E-01
            -2.2334E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -721.577359242637        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  1.0563E+00  1.0365E+00  4.7732E+01  1.3625E+00  2.7434E+00  2.0128E+00  6.5872E+00  2.4004E+00  1.2266E+00  1.3172E-01
             1.2177E+01
 PARAMETER:  1.5479E-01  1.3588E-01  3.9656E+00  4.0935E-01  1.1092E+00  7.9954E-01  1.9851E+00  9.7563E-01  3.0423E-01 -1.9270E+00
             2.5995E+00
 GRADIENT:  -5.3733E-01  1.4638E-01 -1.5266E-01  4.4530E-01  1.4652E-01 -7.9695E-03  2.8130E-01  2.3581E-02 -3.4102E-01  1.7112E-01
            -1.4500E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -721.577867199309        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  1.0572E+00  1.0428E+00  4.9067E+01  1.3595E+00  2.7459E+00  2.0121E+00  6.5701E+00  2.3376E+00  1.2205E+00  1.0823E-01
             1.2186E+01
 PARAMETER:  1.5561E-01  1.4192E-01  3.9932E+00  4.0713E-01  1.1101E+00  7.9918E-01  1.9825E+00  9.4910E-01  2.9923E-01 -2.1235E+00
             2.6003E+00
 GRADIENT:  -3.3251E-01  1.5910E-01 -1.3745E-01  4.3634E-01  7.5000E-02 -7.5895E-02  2.1276E-02  2.1204E-02 -4.8403E-01  1.1527E-01
            -4.7134E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -721.579557501565        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0580E+00  1.0490E+00  5.1700E+01  1.3567E+00  2.7503E+00  2.0114E+00  6.5535E+00  2.2364E+00  1.2148E+00  7.8756E-02
             1.2196E+01
 PARAMETER:  1.5642E-01  1.4782E-01  4.0455E+00  4.0505E-01  1.1117E+00  7.9884E-01  1.9800E+00  9.0488E-01  2.9456E-01 -2.4414E+00
             2.6011E+00
 GRADIENT:  -1.2845E-01  1.7165E-01 -1.1600E-01  4.1719E-01  2.0003E-03 -1.3828E-01 -2.3999E-01  1.7568E-02 -6.1524E-01  6.0903E-02
             4.8581E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -721.591000870539        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  1.0588E+00  1.0548E+00  6.0567E+01  1.3546E+00  2.7628E+00  2.0109E+00  6.5384E+00  1.9816E+00  1.2106E+00  3.3523E-02
             1.2204E+01
 PARAMETER:  1.5718E-01  1.5333E-01  4.2037E+00  4.0352E-01  1.1162E+00  7.9859E-01  1.9777E+00  7.8390E-01  2.9112E-01 -3.2955E+00
             2.6018E+00
 GRADIENT:   7.2679E-02  1.8737E-01 -7.4439E-02  3.6654E-01 -6.1812E-02 -1.8552E-01 -5.0106E-01  1.0200E-02 -7.1052E-01  1.1011E-02
             1.3505E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -721.821481609942        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1231
 NPARAMETR:  1.0569E+00  1.0131E+00  1.0272E+03  1.3852E+00  2.8816E+00  2.0145E+00  6.6479E+00  2.4284E-01  1.2690E+00  1.0000E-02
             1.2179E+01
 PARAMETER:  1.5534E-01  1.1298E-01  7.0346E+00  4.2585E-01  1.1583E+00  8.0039E-01  1.9943E+00 -1.3153E+00  3.3821E-01 -1.7896E+01
             2.5997E+00
 GRADIENT:  -3.8994E-01 -1.1929E-01 -6.4782E-03 -6.5629E-01  4.6373E+00  7.2638E-02  7.9708E-01  5.4123E-07  8.8521E-01  0.0000E+00
            -2.2140E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -722.799273469348        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1412
 NPARAMETR:  1.0610E+00  7.7563E-01  5.4809E+04  1.5100E+00  2.8409E+00  2.0152E+00  7.7405E+00  1.3156E-02  1.3049E+00  1.0000E-02
             1.2253E+01
 PARAMETER:  1.5922E-01 -1.5408E-01  1.1012E+01  5.1213E-01  1.1441E+00  8.0073E-01  2.1465E+00 -4.2308E+00  3.6610E-01 -3.7176E+01
             2.6058E+00
 GRADIENT:  -1.9969E+00  5.4408E-01 -2.4960E-05  4.4359E-01 -4.0866E-01  7.0562E-01  1.5192E+00  5.0130E-07 -2.7317E-01  0.0000E+00
             2.5699E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -722.807294814809        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1604             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0614E+00  7.7011E-01  6.7455E+04  1.5101E+00  2.8427E+00  2.0153E+00  7.6888E+00  1.1344E-02  1.3083E+00  1.0000E-02
             1.2236E+01
 PARAMETER:  1.5963E-01 -1.6122E-01  1.1219E+01  5.1219E-01  1.1448E+00  8.0079E-01  2.1398E+00 -4.3791E+00  3.6872E-01 -3.8252E+01
             2.6044E+00
 GRADIENT:  -1.2425E+00 -1.6553E+00 -3.8788E-05  1.0732E+00  2.7431E-01  1.2242E+00  1.8952E+01 -6.0141E-06  1.5107E+00  0.0000E+00
             9.2754E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -722.807510115850        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1804
 NPARAMETR:  1.0615E+00  7.7016E-01  6.9229E+04  1.5091E+00  2.8445E+00  2.0155E+00  7.6890E+00  1.1195E-02  1.3083E+00  1.0000E-02
             1.2236E+01
 PARAMETER:  1.5965E-01 -1.6116E-01  1.1245E+01  5.1151E-01  1.1454E+00  8.0089E-01  2.1398E+00 -4.3923E+00  3.6871E-01 -3.8252E+01
             2.6044E+00
 GRADIENT:  -4.5732E+00 -1.6756E+00 -3.1241E-05  1.2781E+00  2.3503E-01  1.8147E+00  2.6766E+00 -6.5067E-06  3.4014E-01  0.0000E+00
             4.8854E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -722.807541216435        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1997
 NPARAMETR:  1.0615E+00  7.7023E-01  7.0107E+04  1.5089E+00  2.8446E+00  2.0152E+00  7.6909E+00  1.1145E-02  1.3077E+00  1.0000E-02
             1.2237E+01
 PARAMETER:  1.5970E-01 -1.6107E-01  1.1258E+01  5.1141E-01  1.1454E+00  8.0073E-01  2.1400E+00 -4.3967E+00  3.6829E-01 -3.8252E+01
             2.6044E+00
 GRADIENT:  -1.0642E+01 -1.8733E+00 -3.0601E-05  2.2073E+00  2.6875E-01  8.5709E+00  2.9382E+00  1.5907E-06 -1.7805E-01  0.0000E+00
             7.4511E+00

0ITERATION NO.:  109    OBJECTIVE VALUE:  -722.807617004397        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     2119
 NPARAMETR:  1.0605E+00  7.7005E-01  7.1376E+04  1.5082E+00  2.8450E+00  2.0157E+00  7.6877E+00  1.1192E-02  1.3075E+00  1.0000E-02
             1.2228E+01
 PARAMETER:  1.5860E-01 -1.6116E-01  1.1264E+01  5.1137E-01  1.1454E+00  8.0076E-01  2.1399E+00 -4.3908E+00  3.6829E-01 -3.8252E+01
             2.6044E+00
 GRADIENT:  -1.7425E-01  1.6343E-01 -3.1106E-03  9.7081E-02 -1.2012E-02 -2.1019E-02  1.2181E-02  3.1064E-03  5.7564E-02  0.0000E+00
             1.9917E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2119
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8336E-02  4.6533E-02  4.5699E-10 -8.3499E-02 -2.9478E-06
 SE:             2.7679E-02  2.2844E-02  2.1518E-09  1.4945E-02  9.0202E-05
 N:                     100         100         100         100         100

 P VAL.:         5.0769E-01  4.1648E-02  8.3181E-01  2.3154E-08  9.7393E-01

 ETASHRINKSD(%)  7.2716E+00  2.3471E+01  1.0000E+02  4.9933E+01  9.9698E+01
 ETASHRINKVR(%)  1.4015E+01  4.1433E+01  1.0000E+02  7.4932E+01  9.9999E+01
 EBVSHRINKSD(%)  8.7811E+00  1.7896E+01  1.0000E+02  5.0961E+01  9.9632E+01
 EBVSHRINKVR(%)  1.6791E+01  3.2590E+01  1.0000E+02  7.5952E+01  9.9999E+01
 RELATIVEINF(%)  8.3021E+01  3.7208E+01  0.0000E+00  1.3065E+01  2.1189E-04
 EPSSHRINKSD(%)  4.3477E+00
 EPSSHRINKVR(%)  8.5064E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -722.80761700439712     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       931.28174276401364     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    58.13
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -722.808       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  7.70E-01  7.06E+04  1.51E+00  2.84E+00  2.02E+00  7.69E+00  1.12E-02  1.31E+00  1.00E-02  1.22E+01
 


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
+        1.20E+03
 
 TH 2
+        4.83E+02  1.01E+03
 
 TH 3
+       -9.13E-05  1.87E-04  4.30E-11
 
 TH 4
+       -2.25E+00  1.39E+02 -8.34E-06  1.21E+02
 
 TH 5
+       -1.25E+01 -1.05E+01 -3.43E-06 -1.18E+01  1.94E+01
 
 TH 6
+       -4.65E+01  1.76E+01 -4.97E-06  6.00E+00 -2.62E-01  4.84E+01
 
 TH 7
+       -1.04E+00  8.83E+00 -5.91E-07 -3.30E+00  3.85E-01 -3.06E-01  1.78E+00
 
 TH 8
+        1.54E+03  1.87E+03 -1.16E-04 -4.70E+02 -1.16E+02 -5.40E+01 -2.93E+01  1.14E+04
 
 TH 9
+       -4.23E+01  3.46E+02  2.88E-05 -2.47E+01  7.12E+00  1.26E+00 -5.96E-02  2.50E+02  1.49E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.45E+00 -1.58E+00 -3.21E-07 -6.78E+00  5.67E-01  8.15E-01  2.17E-01  2.59E+00  2.76E+00  0.00E+00  5.87E+00
 
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
 #CPUT: Total CPU Time in Seconds,       75.365
Stop Time:
Sat Sep 18 07:48:55 CDT 2021
