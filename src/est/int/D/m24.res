Wed Sep 29 08:17:46 CDT 2021
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
$DATA ../../../../data/int/D/dat24.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   8824.21178412119        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.6173E+02  1.2492E+02 -2.1098E+01 -1.8837E+02  2.4506E+02 -1.1034E+03 -4.5490E+02 -1.2587E+02 -8.7565E+02 -4.3582E+02
            -2.2297E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1211.17735166224        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  8.6585E-01  8.7100E-01  8.7009E-01  2.7509E+00  9.4721E-01  4.8311E+00  3.6115E+00  1.0685E+00  5.6848E+00  3.2513E+00
             9.8382E+00
 PARAMETER: -4.4042E-02 -3.8108E-02 -3.9157E-02  1.1119E+00  4.5771E-02  1.6751E+00  1.3841E+00  1.6622E-01  1.8378E+00  1.2791E+00
             2.3863E+00
 GRADIENT:  -3.4951E+01 -4.6021E+01 -5.1450E+01  5.5770E+01  3.9604E+01  2.1329E+02  6.5357E+01  6.1974E+00  1.5583E+02  7.9381E+01
             8.1681E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1473.12779200833        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  2.2916E+00  4.5249E-02  6.8993E+02  4.0420E+00  1.9880E+00  5.3588E+00  5.1949E+01  2.0021E+00  2.4497E+00  1.6097E+00
             8.1753E+00
 PARAMETER:  9.2926E-01 -2.9956E+00  6.6366E+00  1.4967E+00  7.8711E-01  1.7787E+00  4.0503E+00  7.9419E-01  9.9597E-01  5.7605E-01
             2.2011E+00
 GRADIENT:   1.2195E+02  4.2379E+00  8.4094E-01  2.4960E+02 -1.0606E+02  2.3088E+02  3.1000E+02 -1.3469E-03 -4.6767E+00  2.9250E+01
             7.1658E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1589.09512009101        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      252
 NPARAMETR:  1.5038E+00  1.9559E-01  3.2163E+02  1.1582E+00  2.3988E+00  5.0961E+00  1.9391E+01  1.8135E+01  8.5290E-01  1.9558E+00
             6.7804E+00
 PARAMETER:  5.0802E-01 -1.5317E+00  5.8734E+00  2.4690E-01  9.7495E-01  1.7285E+00  3.0648E+00  2.9978E+00 -5.9117E-02  7.7082E-01
             2.0140E+00
 GRADIENT:   1.9769E+01 -2.8315E+00  1.8016E-01 -2.7326E+02 -6.7088E+00  1.5017E+02  1.1017E+02  1.2064E-01  1.9171E+01  6.4445E+01
             6.2933E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1842.37623272821        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  1.0455E+00  2.0512E-01  2.3304E+02  1.5626E+00  2.3885E+00  2.1580E+00  2.2744E+01  6.6397E+00  1.1171E+00  1.1985E+00
             4.7140E+00
 PARAMETER:  1.4449E-01 -1.4842E+00  5.5512E+00  5.4638E-01  9.7068E-01  8.6916E-01  3.2243E+00  1.9931E+00  2.1072E-01  2.8110E-01
             1.6505E+00
 GRADIENT:  -2.2455E+01  5.0396E+00  6.6413E-01 -1.0275E+02 -8.4803E+00 -8.5671E+00  1.3782E+02 -2.4496E-02  3.2871E+00 -1.4880E+00
            -2.1237E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1912.42866258681        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      608
 NPARAMETR:  1.0267E+00  9.0498E-01  4.0844E+00  9.9147E-01  2.0762E+00  2.1538E+00  5.0631E+00  3.8196E-01  3.7159E-01  1.5723E-01
             5.1371E+00
 PARAMETER:  1.2635E-01  1.6275E-04  1.5072E+00  9.1431E-02  8.3056E-01  8.6723E-01  1.7220E+00 -8.6244E-01 -8.8997E-01 -1.7500E+00
             1.7365E+00
 GRADIENT:  -3.5130E+01 -2.5167E+01 -1.1308E+01 -5.0363E+00  4.8119E+01 -1.1281E+01  4.3665E+00 -1.9053E-01 -1.4740E+00 -7.8497E-01
             6.9401E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1921.14652756820        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      789
 NPARAMETR:  1.1079E+00  1.5022E+00  3.1455E+00  8.6672E-01  1.9826E+00  2.2059E+00  4.2508E+00  4.9387E-01  3.2073E-01  1.4865E-01
             5.0015E+00
 PARAMETER:  2.0245E-01  5.0693E-01  1.2460E+00 -4.3036E-02  7.8440E-01  8.9112E-01  1.5471E+00 -6.0548E-01 -1.0372E+00 -1.8061E+00
             1.7097E+00
 GRADIENT:   1.2244E+00  1.4319E-01  3.9240E-02 -4.2192E+00 -1.9112E+00  5.9415E-01  3.2657E+00 -5.4978E-01  5.1366E-01 -9.6123E-01
             5.8937E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1921.80704454018        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      969
 NPARAMETR:  1.0856E+00  1.5402E+00  3.5297E+00  8.3803E-01  2.0534E+00  2.2037E+00  4.0214E+00  1.7586E+00  2.8324E-01  1.7152E-01
             4.9673E+00
 PARAMETER:  1.8211E-01  5.3189E-01  1.3612E+00 -7.6700E-02  8.1952E-01  8.9013E-01  1.4916E+00  6.6452E-01 -1.1615E+00 -1.6631E+00
             1.7029E+00
 GRADIENT:  -7.7618E+00 -3.2148E+00 -1.8819E-01 -8.1608E+00  1.1541E+01  3.6476E-01 -5.6254E+00  6.3352E-01  6.3923E-01 -7.7873E-01
            -1.0159E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1922.65395918369        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1145
 NPARAMETR:  1.1048E+00  1.4713E+00  3.2333E+00  8.7636E-01  1.9592E+00  2.2024E+00  4.2475E+00  1.5972E+00  2.9967E-01  1.8155E-01
             4.9644E+00
 PARAMETER:  1.9963E-01  4.8613E-01  1.2735E+00 -3.1978E-02  7.7254E-01  8.8953E-01  1.5463E+00  5.6823E-01 -1.1051E+00 -1.6062E+00
             1.7023E+00
 GRADIENT:  -7.4653E-02 -6.3923E-01 -2.6044E-01 -1.0310E+00 -9.6340E-01  4.8283E-01  6.0536E-01  1.1812E-01  3.2071E-02 -9.8272E-01
            -1.6495E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1929.97451739744        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1325
 NPARAMETR:  1.0889E+00  1.5477E+00  2.6073E+00  8.2017E-01  1.8390E+00  2.1916E+00  4.2604E+00  8.8401E-01  2.5612E-01  1.0620E+00
             4.7531E+00
 PARAMETER:  1.8518E-01  5.3675E-01  1.0583E+00 -9.8244E-02  7.0920E-01  8.8462E-01  1.5494E+00 -2.3291E-02 -1.2621E+00  1.6019E-01
             1.6588E+00
 GRADIENT:  -5.4454E+00 -1.1252E+00 -1.3484E+00 -2.8473E+01 -1.3958E+01 -8.5211E-01  1.5590E+01  1.5246E+00  1.3437E+00  2.4185E+00
            -8.8428E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1932.10890086288        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1501
 NPARAMETR:  1.1056E+00  1.4564E+00  3.5091E+00  8.7964E-01  1.9648E+00  2.1967E+00  4.3070E+00  3.9488E-01  2.6878E-01  1.0676E+00
             4.7955E+00
 PARAMETER:  2.0037E-01  4.7599E-01  1.3554E+00 -2.8241E-02  7.7539E-01  8.8695E-01  1.5602E+00 -8.2917E-01 -1.2139E+00  1.6545E-01
             1.6677E+00
 GRADIENT:   6.2534E-01  2.4805E-01 -8.3425E-01  2.1965E+00  4.0123E+00 -8.3869E-02 -3.3293E-01  2.3268E-01  2.0347E-01  3.2383E-01
             2.4163E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1932.22181314015        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1680
 NPARAMETR:  1.1031E+00  1.4488E+00  3.4434E+00  8.7919E-01  1.9443E+00  2.1983E+00  4.3055E+00  1.3362E-01  2.7430E-01  1.0634E+00
             4.7861E+00
 PARAMETER:  1.9811E-01  4.7077E-01  1.3365E+00 -2.8760E-02  7.6492E-01  8.8766E-01  1.5599E+00 -1.9128E+00 -1.1935E+00  1.6148E-01
             1.6657E+00
 GRADIENT:  -2.0889E-01 -1.3622E-01  1.6558E-01 -3.0037E-01 -1.0971E+00  2.7139E-01 -4.9463E-01  2.6574E-02  2.6501E-01 -7.2949E-02
            -1.4894E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1932.33495767023        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1861             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1043E+00  1.4320E+00  3.4750E+00  8.8237E-01  1.9498E+00  2.2261E+00  4.3944E+00  1.2919E-02  2.3590E-01  1.0665E+00
             4.7890E+00
 PARAMETER:  1.9924E-01  4.5907E-01  1.3456E+00 -2.5149E-02  7.6773E-01  9.0024E-01  1.5803E+00 -4.2491E+00 -1.3444E+00  1.6437E-01
             1.6663E+00
 GRADIENT:   3.7901E+01  2.3675E+01  6.3849E-01  4.9938E+00  1.0499E+01  9.6775E+01  1.8615E+02  2.9980E-04  1.8819E-01  1.9513E-01
             2.9705E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1932.33909684070        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2046
 NPARAMETR:  1.1043E+00  1.4356E+00  3.4828E+00  8.8320E-01  1.9488E+00  2.2260E+00  4.4023E+00  1.0000E-02  2.4361E-01  1.0665E+00
             4.7888E+00
 PARAMETER:  1.9923E-01  4.6162E-01  1.3478E+00 -2.4208E-02  7.6719E-01  9.0019E-01  1.5821E+00 -5.7199E+00 -1.3122E+00  1.6440E-01
             1.6663E+00
 GRADIENT:   2.6800E-01  6.6272E-02 -1.1653E-01 -1.0603E-02  1.9506E-01  5.1768E+00  4.1948E+00  0.0000E+00 -3.4004E-02 -8.6412E-03
            -2.7759E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1932.34034427411        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2235             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1043E+00  1.4317E+00  3.4947E+00  8.8422E-01  1.9478E+00  2.2259E+00  4.4135E+00  1.0000E-02  2.5268E-01  1.0667E+00
             4.7892E+00
 PARAMETER:  1.9923E-01  4.5885E-01  1.3513E+00 -2.3046E-02  7.6670E-01  9.0017E-01  1.5847E+00 -5.7199E+00 -1.2756E+00  1.6455E-01
             1.6664E+00
 GRADIENT:   3.7941E+01  2.3871E+01  9.4459E-01  2.7229E+00  9.3703E+00  9.6761E+01  1.8787E+02  0.0000E+00  3.5403E-01  2.3736E-01
             3.0275E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1932.34183637944        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2419
 NPARAMETR:  1.1043E+00  1.4306E+00  3.4969E+00  8.8536E-01  1.9478E+00  2.2259E+00  4.4150E+00  1.0000E-02  2.5273E-01  1.0666E+00
             4.7890E+00
 PARAMETER:  1.9923E-01  4.5812E-01  1.3519E+00 -2.1755E-02  7.6671E-01  9.0018E-01  1.5850E+00 -5.7199E+00 -1.2754E+00  1.6449E-01
             1.6663E+00
 GRADIENT:   2.7489E-01  2.9261E-02 -3.4028E-02 -6.8796E-01 -1.2080E-01  5.1662E+00  4.6234E+00  0.0000E+00  1.0964E-02  6.4073E-03
            -2.1775E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1932.34277266667        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2605             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1043E+00  1.4281E+00  3.5031E+00  8.8702E-01  1.9478E+00  2.2259E+00  4.4192E+00  1.0000E-02  2.5422E-01  1.0666E+00
             4.7888E+00
 PARAMETER:  1.9921E-01  4.5635E-01  1.3536E+00 -1.9891E-02  7.6670E-01  9.0018E-01  1.5860E+00 -5.7199E+00 -1.2695E+00  1.6448E-01
             1.6663E+00
 GRADIENT:   3.7908E+01  2.3833E+01  8.1240E-01  4.6153E+00  9.7773E+00  9.6762E+01  1.8771E+02  0.0000E+00  2.7918E-01  2.0665E-01
             2.9712E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1932.34322235815        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2790
 NPARAMETR:  1.1043E+00  1.4265E+00  3.5073E+00  8.8770E-01  1.9477E+00  2.2259E+00  4.4221E+00  1.0000E-02  2.5602E-01  1.0666E+00
             4.7887E+00
 PARAMETER:  1.9921E-01  4.5526E-01  1.3549E+00 -1.9121E-02  7.6667E-01  9.0018E-01  1.5866E+00 -5.7199E+00 -1.2625E+00  1.6451E-01
             1.6663E+00
 GRADIENT:   2.6537E-01  9.0772E-02 -9.2844E-02  2.5222E-01  6.1036E-02  5.1727E+00  4.4315E+00  0.0000E+00 -2.6323E-02 -9.5378E-03
            -2.8718E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1932.34359355325        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2979             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1043E+00  1.4239E+00  3.5155E+00  8.8806E-01  1.9474E+00  2.2259E+00  4.4277E+00  1.0000E-02  2.6073E-01  1.0668E+00
             4.7890E+00
 PARAMETER:  1.9922E-01  4.5338E-01  1.3572E+00 -1.8716E-02  7.6650E-01  9.0017E-01  1.5879E+00 -5.7199E+00 -1.2443E+00  1.6466E-01
             1.6663E+00
 GRADIENT:   3.7930E+01  2.3488E+01  9.2941E-01  3.5333E+00  9.4259E+00  9.6757E+01  1.8874E+02  0.0000E+00  3.4688E-01  2.2744E-01
             3.0099E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1932.34383167340        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3161
 NPARAMETR:  1.1043E+00  1.4228E+00  3.5187E+00  8.8856E-01  1.9474E+00  2.2259E+00  4.4297E+00  1.0000E-02  2.6204E-01  1.0668E+00
             4.7890E+00
 PARAMETER:  1.9922E-01  4.5265E-01  1.3581E+00 -1.8154E-02  7.6650E-01  9.0017E-01  1.5883E+00 -5.7199E+00 -1.2393E+00  1.6469E-01
             1.6663E+00
 GRADIENT:   2.7832E-01 -5.6751E-02  2.0572E-02 -8.2835E-01 -2.7550E-01  5.1634E+00  4.8395E+00  0.0000E+00  3.1193E-02  1.1368E-02
             9.8071E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1932.34399522802        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3338
 NPARAMETR:  1.1043E+00  1.4223E+00  3.5203E+00  8.8986E-01  1.9478E+00  2.2259E+00  4.4298E+00  1.0000E-02  2.6074E-01  1.0667E+00
             4.7886E+00
 PARAMETER:  1.9922E-01  4.5201E-01  1.3589E+00 -1.7695E-02  7.6651E-01  9.0017E-01  1.5887E+00 -5.7199E+00 -1.2351E+00  1.6472E-01
             1.6663E+00
 GRADIENT:   7.6651E-03 -3.1084E-02  2.5970E-02 -5.6513E-01 -1.8538E-01 -4.6837E-03  1.3136E-01  0.0000E+00  1.8312E-02  1.0398E-02
             2.2180E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3338
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.2205E-03  7.4725E-03 -2.0615E-05 -2.4842E-02 -8.4967E-03
 SE:             2.9415E-02  2.7927E-02  4.0546E-05  6.1384E-03  1.9926E-02
 N:                     100         100         100         100         100

 P VAL.:         8.5914E-01  7.8903E-01  6.1115E-01  5.1928E-05  6.6980E-01

 ETASHRINKSD(%)  1.4554E+00  6.4397E+00  9.9864E+01  7.9436E+01  3.3247E+01
 ETASHRINKVR(%)  2.8896E+00  1.2465E+01  1.0000E+02  9.5771E+01  5.5440E+01
 EBVSHRINKSD(%)  1.6054E+00  3.3246E+00  9.9840E+01  8.5836E+01  3.3653E+01
 EBVSHRINKVR(%)  3.1851E+00  6.5387E+00  1.0000E+02  9.7994E+01  5.5981E+01
 RELATIVEINF(%)  9.6778E+01  5.3976E+01  6.3413E-05  9.3560E-01  1.2375E+01
 EPSSHRINKSD(%)  1.1507E+01
 EPSSHRINKVR(%)  2.1689E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1932.3439952280221     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -278.25463545961134     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   109.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1932.344       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.10E+00  1.42E+00  3.52E+00  8.89E-01  1.95E+00  2.23E+00  4.43E+00  1.00E-02  2.63E-01  1.07E+00  4.79E+00
 


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
+        1.76E+02
 
 TH 2
+       -5.36E-01  2.62E+01
 
 TH 3
+        2.60E-01  4.59E-01  2.47E+00
 
 TH 4
+       -4.68E+00  5.16E+01 -1.31E+01  5.92E+02
 
 TH 5
+       -1.49E+00 -5.23E+00 -1.39E+01  6.67E+01  1.12E+02
 
 TH 6
+        3.25E-02 -2.03E-02  3.51E-02  1.29E+00 -8.03E-01  3.82E+01
 
 TH 7
+        2.55E-01  1.47E+00 -2.74E-01 -3.57E+01  1.82E+00 -1.90E-01  7.95E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.06E-01 -2.69E+00 -2.36E-01 -7.91E+01  3.08E+00 -2.87E-01  3.94E+00  0.00E+00  2.41E+01
 
 TH10
+       -8.01E-01 -1.96E-01 -5.65E-01 -2.35E+00 -3.34E+00  3.11E-01 -4.69E-02  0.00E+00  1.36E+00  3.27E+01
 
 TH11
+       -4.83E+00 -2.77E+00  1.02E-01 -2.23E+01 -1.97E+00  5.37E-01  1.24E+00  0.00E+00  3.71E+00  8.20E+00  5.02E+01
 
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
 #CPUT: Total CPU Time in Seconds,      124.370
Stop Time:
Wed Sep 29 08:19:52 CDT 2021
