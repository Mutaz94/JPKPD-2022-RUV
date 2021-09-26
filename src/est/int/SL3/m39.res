Sat Sep 25 02:13:54 CDT 2021
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
$DATA ../../../../data/int/SL3/dat39.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF DATA RECS IN DATA SET:      979
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      879
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -633.816225554604        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.9246E+01  5.3580E+00  1.4804E+02  1.0323E+01 -1.4017E+01  1.9306E+01 -7.2449E+01 -2.4503E+02 -8.2082E+01 -4.9711E+01
            -5.8365E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2708.32850692167        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0212E+00  1.4911E+00  1.0889E+00  8.1502E-01  1.4468E+00  8.2979E-01  9.8135E-01  8.9132E-01  7.3285E-01  1.1755E+00
             2.6384E+00
 PARAMETER:  1.2103E-01  4.9950E-01  1.8512E-01 -1.0455E-01  4.6932E-01 -8.6579E-02  8.1170E-02 -1.5053E-02 -2.1082E-01  2.6169E-01
             1.0702E+00
 GRADIENT:   4.8111E+01  7.0632E+01 -1.9570E+01  9.4446E+01  2.6139E+01 -5.5487E+01  1.7928E+01  1.7453E+00 -1.2549E+01 -4.4084E+01
            -2.4773E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2715.80640158395        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0109E+00  1.4039E+00  2.1463E+00  8.7891E-01  1.6904E+00  9.3910E-01  9.3413E-01  6.1810E-01  4.8801E-01  1.5560E+00
             2.7046E+00
 PARAMETER:  1.1085E-01  4.3929E-01  8.6376E-01 -2.9071E-02  6.2499E-01  3.7169E-02  3.1857E-02 -3.8110E-01 -6.1741E-01  5.4211E-01
             1.0949E+00
 GRADIENT:   1.3520E+01  7.4450E+01 -2.0349E+01  1.4553E+02  7.3213E+00 -3.3123E+00 -1.9786E+01 -1.5857E+00 -1.4965E+01 -7.8675E+00
             4.8744E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2746.76939104323        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8691E-01  1.3550E+00  1.2368E+01  8.5969E-01  2.3776E+00  9.4606E-01  1.0657E+00  3.2968E+00  3.5939E-01  1.8959E+00
             2.5574E+00
 PARAMETER:  8.6823E-02  4.0377E-01  2.6151E+00 -5.1178E-02  9.6608E-01  4.4546E-02  1.6367E-01  1.2930E+00 -9.2335E-01  7.3970E-01
             1.0390E+00
 GRADIENT:  -3.5265E+01 -3.1023E+01 -1.0878E+00 -2.6990E+01  3.6943E+01 -1.5669E+00  1.9685E-01 -5.8953E+00  4.6812E-01 -1.7624E+01
             1.0173E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2750.72403607233        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0007E+00  1.3506E+00  1.5588E+01  8.7519E-01  2.2875E+00  9.5445E-01  1.0827E+00  4.4256E+00  3.2370E-01  1.9884E+00
             2.5476E+00
 PARAMETER:  1.0069E-01  4.0058E-01  2.8465E+00 -3.3311E-02  9.2746E-01  5.3381E-02  1.7950E-01  1.5874E+00 -1.0279E+00  7.8733E-01
             1.0352E+00
 GRADIENT:  -2.0753E+00 -8.7814E+00 -2.4871E+00 -5.5173E+00  9.8748E+00  2.4565E+00 -9.1764E-01 -3.3032E+00 -2.6152E-01  2.0162E+00
             7.4096E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2751.80033510451        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0012E+00  1.2639E+00  2.7091E+01  9.3676E-01  2.3579E+00  9.6239E-01  1.1617E+00  6.0235E+00  2.6149E-01  1.9844E+00
             2.5437E+00
 PARAMETER:  1.0115E-01  3.3419E-01  3.3992E+00  3.4671E-02  9.5779E-01  6.1666E-02  2.4984E-01  1.8957E+00 -1.2414E+00  7.8529E-01
             1.0336E+00
 GRADIENT:  -1.5981E+00 -6.9557E+00 -4.5462E+00  3.5718E+00  1.7352E+01  5.2983E+00 -2.4891E+00  3.4129E+00 -7.5099E-01 -8.4976E-02
             5.5933E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2752.39565568909        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0031E+00  1.2363E+00  3.2999E+01  9.5733E-01  2.3087E+00  9.5448E-01  1.1985E+00  6.5765E+00  2.4855E-01  1.9746E+00
             2.5350E+00
 PARAMETER:  1.0306E-01  3.1214E-01  3.5965E+00  5.6395E-02  9.3670E-01  5.3408E-02  2.8105E-01  1.9835E+00 -1.2921E+00  7.8036E-01
             1.0302E+00
 GRADIENT:   2.9776E+00 -2.6397E+00 -1.0743E+01  1.6268E+01  2.0273E+00  2.4343E+00 -3.4582E-01  1.4523E+01 -9.9059E-01 -8.1450E+00
            -3.7673E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2752.70049840894        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0033E+00  1.2398E+00  3.7900E+01  9.5624E-01  2.2484E+00  9.4338E-01  1.1983E+00  7.0091E+00  2.6365E-01  1.9496E+00
             2.5262E+00
 PARAMETER:  1.0326E-01  3.1493E-01  3.7349E+00  5.5256E-02  9.1022E-01  4.1710E-02  2.8087E-01  2.0472E+00 -1.2331E+00  7.6761E-01
             1.0267E+00
 GRADIENT:   3.8313E+00  1.1115E-01 -7.3811E+00  1.1878E+01 -8.9676E+00 -1.9653E+00  1.3167E+00  1.0422E+01 -8.6827E-01 -5.3565E+00
            -1.1951E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2753.39578238792        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      599
 NPARAMETR:  1.0026E+00  1.3025E+00  4.5544E+01  9.1861E-01  2.2699E+00  9.5115E-01  1.1120E+00  8.0627E+00  3.6502E-01  1.8982E+00
             2.5440E+00
 PARAMETER:  1.0260E-01  3.6427E-01  3.9187E+00  1.5106E-02  9.1974E-01  4.9915E-02  2.0616E-01  2.1872E+00 -9.0782E-01  7.4090E-01
             1.0337E+00
 GRADIENT:   1.6844E+00  4.4096E-01  8.0119E-01 -3.0654E+00  1.2668E+00  9.2156E-01  3.9631E-01 -1.6169E+00 -5.4165E-02  8.5360E-01
             6.2318E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2753.49416738874        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      777
 NPARAMETR:  1.0052E+00  1.3393E+00  4.9263E+01  8.9939E-01  2.2835E+00  9.5113E-01  1.0679E+00  8.6177E+00  4.2526E-01  1.8824E+00
             2.5441E+00
 PARAMETER:  1.0522E-01  3.9212E-01  3.9972E+00 -6.0408E-03  9.2572E-01  4.9897E-02  1.6573E-01  2.2538E+00 -7.5506E-01  7.3254E-01
             1.0338E+00
 GRADIENT:   6.5827E-01 -5.0093E-01 -1.1448E+00  1.4667E+00 -1.2463E+00  1.3149E-01  5.3684E-01  2.4368E+00  3.0808E-02  7.2992E-01
            -2.9941E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2753.50724118623        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      953
 NPARAMETR:  1.0049E+00  1.3902E+00  4.7783E+01  8.6350E-01  2.2876E+00  9.5077E-01  1.0418E+00  8.7580E+00  4.0303E-01  1.8772E+00
             2.5454E+00
 PARAMETER:  1.0488E-01  4.2941E-01  3.9667E+00 -4.6767E-02  9.2752E-01  4.9514E-02  1.4099E-01  2.2700E+00 -8.0873E-01  7.2979E-01
             1.0343E+00
 GRADIENT:  -2.3281E-01 -7.2719E-01 -8.7294E-01  7.5132E-01 -7.6105E-01 -1.2884E-02  2.8925E-01  1.9129E+00 -9.4586E-02 -6.8652E-01
            -1.4444E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2753.51549720171        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1128
 NPARAMETR:  1.0050E+00  1.4299E+00  4.4132E+01  8.3641E-01  2.2903E+00  9.5082E-01  1.0168E+00  8.6381E+00  4.0847E-01  1.8888E+00
             2.5445E+00
 PARAMETER:  1.0496E-01  4.5760E-01  3.8872E+00 -7.8638E-02  9.2870E-01  4.9574E-02  1.1663E-01  2.2562E+00 -7.9533E-01  7.3593E-01
             1.0339E+00
 GRADIENT:  -3.4142E-02 -2.0079E-01 -1.6520E-01  4.0398E-02 -8.3081E-02  4.4957E-02  1.0704E-01  2.1846E-01  2.1056E-02 -4.4181E-01
             1.9730E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2753.51597943118        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1305            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0050E+00  1.4449E+00  4.3289E+01  8.2611E-01  2.2915E+00  9.5073E-01  1.0084E+00  8.6393E+00  4.0587E-01  1.8907E+00
             2.5445E+00
 PARAMETER:  1.0498E-01  4.6806E-01  3.8679E+00 -9.1026E-02  9.2923E-01  4.9476E-02  1.0837E-01  2.2563E+00 -8.0173E-01  7.3697E-01
             1.0339E+00
 GRADIENT:   7.2336E+00  7.1844E+00  5.7531E-02  1.2161E+00  3.1604E+00  6.8115E-01  1.3801E-01  7.4289E-01  1.4812E-01  4.6900E-01
             2.3766E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2753.51648633717        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1481
 NPARAMETR:  1.0048E+00  1.4451E+00  4.3392E+01  8.2577E-01  2.2918E+00  9.5064E-01  1.0094E+00  8.6354E+00  4.0200E-01  1.8939E+00
             2.5441E+00
 PARAMETER:  1.0481E-01  4.6821E-01  3.8703E+00 -9.1435E-02  9.2936E-01  4.9385E-02  1.0935E-01  2.2559E+00 -8.1131E-01  7.3861E-01
             1.0338E+00
 GRADIENT:  -3.9150E-01 -9.4997E-02 -6.0332E-02 -6.8770E-03 -5.8245E-02 -3.2134E-02  2.4605E-02  2.0333E-02 -4.9927E-03  8.8033E-02
            -1.2656E-02

0ITERATION NO.:   67    OBJECTIVE VALUE:  -2753.51651567570        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1538
 NPARAMETR:  1.0049E+00  1.4451E+00  4.3394E+01  8.2578E-01  2.2917E+00  9.5069E-01  1.0093E+00  8.6356E+00  4.0227E-01  1.8936E+00
             2.5440E+00
 PARAMETER:  1.0489E-01  4.6815E-01  3.8703E+00 -9.1425E-02  9.2929E-01  4.9429E-02  1.0923E-01  2.2559E+00 -8.1064E-01  7.3850E-01
             1.0338E+00
 GRADIENT:  -1.9652E-01 -1.7214E-01 -3.6240E-02 -1.0016E-01 -5.2954E-02 -1.3221E-02  1.5370E-03 -3.9349E-02 -2.3037E-03  2.1357E-02
             5.0698E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1538
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4528E-03 -1.8385E-03 -2.3420E-02 -1.6768E-02 -2.2421E-02
 SE:             2.9368E-02  2.6532E-02  9.8142E-03  9.4666E-03  2.6636E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6055E-01  9.4476E-01  1.7019E-02  7.6513E-02  3.9993E-01

 ETASHRINKSD(%)  1.6118E+00  1.1113E+01  6.7121E+01  6.8286E+01  1.0765E+01
 ETASHRINKVR(%)  3.1977E+00  2.0991E+01  8.9190E+01  8.9942E+01  2.0372E+01
 EBVSHRINKSD(%)  1.6098E+00  1.1062E+01  7.6947E+01  6.9825E+01  9.2281E+00
 EBVSHRINKVR(%)  3.1937E+00  2.0900E+01  9.4686E+01  9.0894E+01  1.7605E+01
 RELATIVEINF(%)  9.6727E+01  4.6131E+00  4.6165E+00  5.2840E-01  7.4511E+01
 EPSSHRINKSD(%)  1.6406E+01
 EPSSHRINKVR(%)  3.0120E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          879
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1615.4939413738146     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2753.5165156756962     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1138.0225743018816     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.40
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2753.517       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.45E+00  4.34E+01  8.26E-01  2.29E+00  9.51E-01  1.01E+00  8.64E+00  4.02E-01  1.89E+00  2.54E+00
 


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
+       -1.95E+01  4.43E+02
 
 TH 3
+       -2.07E-01  1.60E+00  1.07E-01
 
 TH 4
+       -7.81E+01  5.21E+02 -4.24E+00  1.19E+03
 
 TH 5
+       -4.16E+00  2.04E+01  1.84E+00 -8.62E+01  8.40E+01
 
 TH 6
+        3.85E+01 -5.60E+00  2.02E-01  6.87E+00  4.18E+00  1.91E+02
 
 TH 7
+        3.89E+01 -1.19E+00 -1.35E+00  5.51E+01 -2.53E+01 -2.08E+01  1.81E+02
 
 TH 8
+        2.57E+00 -2.50E+01 -1.43E+01  6.41E+01 -2.94E+01 -3.42E+00  2.06E+01  1.30E+02
 
 TH 9
+       -9.74E+00  3.98E-01  7.20E-01 -9.86E+01  1.47E+01  3.52E+00  3.41E+01 -1.11E+01  2.31E+01
 
 TH10
+        3.89E+00 -4.89E+01 -2.91E+00  1.30E+02 -6.13E+01 -6.55E+00  4.12E+01  4.70E+01 -2.13E+01  1.29E+02
 
 TH11
+       -2.77E+01  9.71E+01  8.17E+00 -3.18E+02  1.33E+02  1.58E+01 -8.73E+01 -9.48E+02  5.45E+01 -2.08E+02  9.08E+02
 
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
 #CPUT: Total CPU Time in Seconds,       50.450
Stop Time:
Sat Sep 25 02:14:46 CDT 2021
