Sat Sep 18 06:48:58 CDT 2021
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
$DATA ../../../../data/int/D/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25001.8845741990        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5297E+02  2.8272E+02 -4.5853E+01  1.1287E+02  2.5742E+02 -1.8791E+03 -9.4186E+02 -4.7786E+01 -1.5048E+03 -6.8249E+02
            -5.2577E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -981.402809123739        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.6091E+00  1.8688E+00  1.0283E+00  2.0537E+00  8.3735E-01  4.2938E+00  4.7413E+00  9.7734E-01  3.0156E+00  1.9294E+00
             1.2485E+01
 PARAMETER:  5.7571E-01  7.2530E-01  1.2790E-01  8.1964E-01 -7.7512E-02  1.5572E+00  1.6563E+00  7.7084E-02  1.2038E+00  7.5719E-01
             2.6245E+00
 GRADIENT:   1.6067E+01  2.3275E+01 -3.5569E+01  7.0010E+01 -3.3214E+01  1.3772E+02  5.5037E+01  4.2175E+00  6.4545E+01  4.7724E+01
             5.1211E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1052.42939362908        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.4505E+00  1.3898E+00  5.9924E+01  4.1679E+00  2.4295E+00  2.8009E+00  1.2538E+01  6.9890E-01  3.0976E+00  2.4576E+00
             1.2231E+01
 PARAMETER:  4.7191E-01  4.2914E-01  4.1931E+00  1.5274E+00  9.8769E-01  1.1299E+00  2.6288E+00 -2.5825E-01  1.2306E+00  9.9918E-01
             2.6040E+00
 GRADIENT:   1.0379E+01  2.4131E+01 -2.5726E+00  9.0872E+01  5.3102E+00  6.7117E+01  2.4723E+01  1.5739E-02  1.4584E+01  7.0697E+01
             4.9425E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1300.95479986616        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0043E+00  1.1366E+00  9.4085E+00  1.2138E+00  1.9678E+00  1.9003E+00  5.6934E+00  7.2943E-01  1.7317E+00  4.6541E-01
             8.9145E+00
 PARAMETER:  1.0432E-01  2.2803E-01  2.3416E+00  2.9372E-01  7.7691E-01  7.4200E-01  1.8393E+00 -2.1549E-01  6.4911E-01 -6.6484E-01
             2.2877E+00
 GRADIENT:  -7.4413E+01 -1.9531E+00 -3.6309E+00 -4.7684E+01  1.0686E+01 -7.5555E+00  5.3309E+01  6.4996E-02  2.2746E+01  3.8652E+00
             2.5810E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1336.89871123131        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1363E+00  8.4443E-01  7.7838E+00  1.3622E+00  1.8037E+00  2.2971E+00  5.1647E+00  7.4995E-01  1.4161E+00  3.6218E-01
             7.6554E+00
 PARAMETER:  2.2776E-01 -6.9089E-02  2.1520E+00  4.0910E-01  6.8983E-01  9.3165E-01  1.7419E+00 -1.8775E-01  4.4793E-01 -9.1561E-01
             2.1354E+00
 GRADIENT:  -1.3028E+01 -4.1749E+00 -2.0186E+00 -1.0863E+00  5.4668E+00  2.1395E+00  3.9354E+00  9.4371E-02  2.9988E+00  1.2107E+00
             2.8449E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1337.09970069236        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.1697E+00  1.0611E+00  5.9750E+00  1.2321E+00  1.7656E+00  2.2712E+00  4.7362E+00  5.7101E-01  1.2434E+00  3.3188E-01
             7.6731E+00
 PARAMETER:  2.5673E-01  1.5933E-01  1.8876E+00  3.0872E-01  6.6849E-01  9.2033E-01  1.6552E+00 -4.6036E-01  3.1788E-01 -1.0030E+00
             2.1377E+00
 GRADIENT:  -1.8434E+00 -1.7681E+00 -2.3469E+00  2.9125E+00  1.3102E+00 -7.7756E-01  2.6682E+00  8.0082E-02 -3.0969E-01  9.4463E-01
             9.0364E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1337.10993022547        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  1.1744E+00  1.1549E+00  5.8949E+00  1.1808E+00  1.7889E+00  2.2695E+00  4.5644E+00  5.1846E-01  1.2017E+00  3.1679E-01
             7.6775E+00
 PARAMETER:  2.6071E-01  2.4405E-01  1.8741E+00  2.6622E-01  6.8160E-01  9.1954E-01  1.6183E+00 -5.5689E-01  2.8371E-01 -1.0495E+00
             2.1383E+00
 GRADIENT:  -2.6150E-01 -8.6483E-01 -1.4852E+00  2.3083E+00  4.9238E-01 -9.7200E-01  1.4577E+00  5.2861E-02 -5.5224E-01  8.4457E-01
            -7.1876E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1337.11145200755        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      511
 NPARAMETR:  1.1751E+00  1.2003E+00  5.9454E+00  1.1563E+00  1.8063E+00  2.2708E+00  4.4851E+00  4.8484E-01  1.1851E+00  3.0530E-01
             7.6797E+00
 PARAMETER:  2.6134E-01  2.8259E-01  1.8826E+00  2.4521E-01  6.9128E-01  9.2012E-01  1.6008E+00 -6.2393E-01  2.6981E-01 -1.0865E+00
             2.1386E+00
 GRADIENT:  -4.9904E-03 -4.9623E-01 -9.7266E-01  1.6473E+00  2.5761E-01 -7.5716E-01  8.4550E-01  3.8289E-02 -5.0904E-01  7.7040E-01
            -9.2735E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1337.11246644542        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.1752E+00  1.2246E+00  5.9980E+00  1.1433E+00  1.8173E+00  2.2718E+00  4.4440E+00  4.5758E-01  1.1771E+00  2.9564E-01
             7.6813E+00
 PARAMETER:  2.6143E-01  3.0260E-01  1.8914E+00  2.3393E-01  6.9733E-01  9.2058E-01  1.5915E+00 -6.8181E-01  2.6309E-01 -1.1186E+00
             2.1388E+00
 GRADIENT:   3.2080E-02 -3.0812E-01 -6.8139E-01  1.2322E+00  1.6633E-01 -5.9746E-01  5.1951E-01  2.9699E-02 -4.5187E-01  7.1050E-01
            -9.5553E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1337.11342173930        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  1.1752E+00  1.2530E+00  6.0809E+00  1.1283E+00  1.8314E+00  2.2732E+00  4.3969E+00  4.1244E-01  1.1685E+00  2.7925E-01
             7.6836E+00
 PARAMETER:  2.6147E-01  3.2554E-01  1.9052E+00  2.2068E-01  7.0509E-01  9.2119E-01  1.5809E+00 -7.8566E-01  2.5575E-01 -1.1757E+00
             2.1391E+00
 GRADIENT:   4.5372E-02 -9.1966E-02 -3.3215E-01  7.1408E-01  7.1365E-02 -3.9062E-01  1.4277E-01  1.9545E-02 -3.6266E-01  6.1581E-01
            -9.2910E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1337.11485616157        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.1752E+00  1.2758E+00  6.1673E+00  1.1164E+00  1.8439E+00  2.2743E+00  4.3602E+00  3.6059E-01  1.1623E+00  2.5997E-01
             7.6861E+00
 PARAMETER:  2.6144E-01  3.4357E-01  1.9193E+00  2.1010E-01  7.1191E-01  9.2168E-01  1.5725E+00 -9.2001E-01  2.5041E-01 -1.2472E+00
             2.1394E+00
 GRADIENT:   3.1689E-02  8.4153E-02 -4.9272E-02  2.8797E-01  1.0961E-02 -2.2498E-01 -1.5096E-01  1.1868E-02 -2.7432E-01  5.1616E-01
            -8.4068E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1337.11962917671        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      794
 NPARAMETR:  1.1752E+00  1.3011E+00  6.2975E+00  1.1033E+00  1.8598E+00  2.2758E+00  4.3200E+00  2.7430E-01  1.1564E+00  2.2644E-01
             7.6897E+00
 PARAMETER:  2.6140E-01  3.6324E-01  1.9402E+00  1.9827E-01  7.2048E-01  9.2232E-01  1.5633E+00 -1.1935E+00  2.4533E-01 -1.3853E+00
             2.1399E+00
 GRADIENT:   1.6122E-02  2.6375E-01  2.8283E-01 -2.4939E-01 -6.1799E-02 -9.5343E-03 -4.7539E-01  4.5965E-03 -1.5061E-01  3.7049E-01
            -6.7903E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1337.14400534063        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  1.1747E+00  1.3292E+00  6.5198E+00  1.0903E+00  1.8815E+00  2.2768E+00  4.2786E+00  1.0061E-01  1.1549E+00  1.4268E-01
             7.6967E+00
 PARAMETER:  2.6097E-01  3.8456E-01  1.9748E+00  1.8647E-01  7.3209E-01  9.2276E-01  1.5536E+00 -2.1965E+00  2.4399E-01 -1.8471E+00
             2.1408E+00
 GRADIENT:  -1.8416E-01  5.6139E-01  6.6694E-01 -7.6594E-01 -9.0723E-02  1.2610E-01 -7.9365E-01  1.6597E-04  4.4821E-02  1.3185E-01
            -1.7458E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1337.78250508758        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      987
 NPARAMETR:  1.1892E+00  1.1208E+00  7.0456E+00  1.2240E+00  1.8416E+00  2.3008E+00  4.8409E+00  1.0000E-02  1.2433E+00  1.0000E-02
             7.7190E+00
 PARAMETER:  2.7329E-01  2.1402E-01  2.0524E+00  3.0209E-01  7.1066E-01  9.3325E-01  1.6771E+00 -1.0108E+01  3.1780E-01 -5.1169E+00
             2.1437E+00
 GRADIENT:   2.2682E+00  1.8153E+00 -2.2820E-01 -1.3223E+00 -6.5264E-01 -1.2750E-01  2.6098E+00  0.0000E+00  8.8985E-02  0.0000E+00
             2.1999E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1337.94243430141        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1162
 NPARAMETR:  1.1818E+00  9.4152E-01  8.3443E+00  1.3287E+00  1.8542E+00  2.3035E+00  5.1461E+00  1.0000E-02  1.3310E+00  1.0000E-02
             7.7058E+00
 PARAMETER:  2.6706E-01  3.9738E-02  2.2216E+00  3.8417E-01  7.1744E-01  9.3443E-01  1.7382E+00 -1.2679E+01  3.8596E-01 -6.1746E+00
             2.1420E+00
 GRADIENT:   4.6254E-02 -9.5655E-03  7.8553E-03 -1.9695E-01 -1.3293E-02  2.1896E-02  6.7752E-02  0.0000E+00  1.4799E-01  0.0000E+00
            -1.7695E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1337.94265483781        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1338
 NPARAMETR:  1.1817E+00  9.4168E-01  8.3401E+00  1.3284E+00  1.8544E+00  2.3032E+00  5.1460E+00  1.0000E-02  1.3269E+00  1.0000E-02
             7.7069E+00
 PARAMETER:  2.6699E-01  3.9906E-02  2.2211E+00  3.8401E-01  7.1756E-01  9.3431E-01  1.7382E+00 -1.2622E+01  3.8287E-01 -6.1514E+00
             2.1421E+00
 GRADIENT:   4.2532E-03 -1.8498E-03 -1.0919E-03  3.4920E-03  6.1137E-04 -3.3417E-03  2.9183E-03  0.0000E+00 -5.1797E-03  0.0000E+00
            -2.2638E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1337.94265636439        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1519
 NPARAMETR:  1.1817E+00  9.4168E-01  8.3405E+00  1.3285E+00  1.8544E+00  2.3033E+00  5.1463E+00  1.0000E-02  1.3270E+00  1.0000E-02
             7.7069E+00
 PARAMETER:  2.6699E-01  3.9906E-02  2.2211E+00  3.8402E-01  7.1757E-01  9.3435E-01  1.7383E+00 -1.2620E+01  3.8291E-01 -6.1506E+00
             2.1421E+00
 GRADIENT:   5.4766E-03  1.6963E-03 -1.3816E-03 -1.4971E-03  2.4684E-03  1.0382E-02  1.5196E-02  0.0000E+00 -2.0205E-03  0.0000E+00
             3.1553E-03

0ITERATION NO.:   81    OBJECTIVE VALUE:  -1337.94265636439        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1541
 NPARAMETR:  1.1817E+00  9.4168E-01  8.3405E+00  1.3285E+00  1.8544E+00  2.3033E+00  5.1463E+00  1.0000E-02  1.3270E+00  1.0000E-02
             7.7069E+00
 PARAMETER:  2.6699E-01  3.9906E-02  2.2211E+00  3.8402E-01  7.1757E-01  9.3435E-01  1.7383E+00 -1.2620E+01  3.8291E-01 -6.1506E+00
             2.1421E+00
 GRADIENT:   5.4766E-03  1.6963E-03 -1.3816E-03 -1.4971E-03  2.4684E-03  1.0382E-02  1.5196E-02  0.0000E+00 -2.0205E-03  0.0000E+00
             3.1553E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1541
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1631E-03  3.5833E-02 -4.0027E-05 -6.9886E-02  6.0250E-05
 SE:             2.9118E-02  2.4170E-02  2.7387E-05  1.5346E-02  1.7493E-04
 N:                     100         100         100         100         100

 P VAL.:         7.7921E-01  1.3820E-01  1.4388E-01  5.2671E-06  7.3053E-01

 ETASHRINKSD(%)  2.4515E+00  1.9026E+01  9.9908E+01  4.8590E+01  9.9414E+01
 ETASHRINKVR(%)  4.8429E+00  3.4433E+01  1.0000E+02  7.3570E+01  9.9997E+01
 EBVSHRINKSD(%)  2.4609E+00  1.5044E+01  9.9901E+01  4.9104E+01  9.9408E+01
 EBVSHRINKVR(%)  4.8612E+00  2.7825E+01  1.0000E+02  7.4096E+01  9.9996E+01
 RELATIVEINF(%)  9.5003E+01  2.9604E+01  1.7081E-05  1.0573E+01  6.4971E-04
 EPSSHRINKSD(%)  7.6603E+00
 EPSSHRINKVR(%)  1.4734E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1337.9426563643901     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       316.14670340402063     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.55
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    15.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1337.943       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.18E+00  9.42E-01  8.34E+00  1.33E+00  1.85E+00  2.30E+00  5.15E+00  1.00E-02  1.33E+00  1.00E-02  7.71E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        6.77E+01
 
 TH 2
+        4.92E+00  1.13E+01
 
 TH 3
+        1.41E-01  2.49E-01  1.89E-01
 
 TH 4
+       -5.87E+00  2.82E+01 -1.48E+00  9.94E+01
 
 TH 5
+       -2.55E+00 -9.25E+00 -4.29E+00  2.31E+01  9.85E+01
 
 TH 6
+        1.29E+01 -3.07E+00  6.12E-03 -1.24E+01  1.35E+00  3.97E+00
 
 TH 7
+        2.04E+00 -2.20E+00 -1.36E-02 -6.70E+00  1.13E+00  1.26E+00  5.70E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.95E+00 -7.57E+00 -1.43E-01 -2.09E+01  5.89E+00  3.43E+00  1.76E+00  0.00E+00  5.68E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.03E+00 -2.27E+00  5.85E-02 -6.94E+00 -5.98E-01  7.30E-01  4.65E-01  0.00E+00  1.55E+00  0.00E+00  9.64E-01
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.40E+02
 
 TH 2
+       -6.08E-01  4.57E+01
 
 TH 3
+        1.45E-01  4.60E-01  2.38E-01
 
 TH 4
+       -2.54E+00  3.47E+01 -1.01E+00  1.14E+02
 
 TH 5
+       -2.16E+00 -1.51E+01 -4.34E+00  9.34E+00  1.02E+02
 
 TH 6
+        5.88E-01 -4.32E-01  2.68E-02  3.65E-01 -9.38E-01  3.36E+01
 
 TH 7
+        1.94E-01  4.73E+00 -1.18E-01 -8.29E+00  2.10E+00 -1.60E-01  4.20E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.63E-01 -1.98E+00 -3.27E-01 -2.42E+01  8.35E+00 -5.76E-01  2.49E+00  0.00E+00  1.71E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.62E+00 -3.13E+00 -6.21E-02 -7.92E+00  5.18E-01  1.42E+00  5.45E-01  0.00E+00  2.48E+00  0.00E+00  1.80E+01
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.43E+02
 
 TH 2
+        5.67E+01  3.51E+01
 
 TH 3
+        5.43E-01  3.77E-01  1.68E-01
 
 TH 4
+        6.00E+01  3.52E+01 -2.09E-01  1.17E+02
 
 TH 5
+       -2.35E+01 -1.17E+01 -3.85E+00 -3.70E+00  1.06E+02
 
 TH 6
+        5.03E+01  1.86E+01 -1.11E-01 -8.91E+00 -1.09E-01  5.73E+01
 
 TH 7
+        6.17E+00  4.26E+00 -1.33E-01 -8.77E+00  4.59E+00  8.96E+00  4.07E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.41E+00 -2.92E+00 -4.00E-01 -2.24E+01  1.21E+01  5.39E+00  2.59E+00  0.00E+00  1.64E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.07E+01 -2.94E+01 -2.61E-01 -3.18E+01  1.67E+01 -4.57E+00 -1.30E+00  0.00E+00  1.80E+00  0.00E+00  6.88E+02
 
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
 #CPUT: Total CPU Time in Seconds,       52.640
Stop Time:
Sat Sep 18 06:49:52 CDT 2021
