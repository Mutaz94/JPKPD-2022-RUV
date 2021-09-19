Sat Sep 18 12:45:55 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat30.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1626.98955541237        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8771E+01 -2.7110E+01 -3.9643E+01  3.6633E+01  6.9574E+01  3.1053E-01 -1.6264E+01 -1.1192E+00  4.9570E+00 -3.3667E+01
            -8.2400E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1641.95481114937        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0218E+00  1.2068E+00  1.1511E+00  8.8198E-01  1.1622E+00  1.0063E+00  1.2283E+00  9.8514E-01  9.0877E-01  1.2566E+00
             1.1850E+00
 PARAMETER:  1.2155E-01  2.8794E-01  2.4075E-01 -2.5587E-02  2.5029E-01  1.0625E-01  3.0563E-01  8.5033E-02  4.3367E-03  3.2842E-01
             2.6972E-01
 GRADIENT:   5.7937E+01  1.7472E+01 -7.5415E+00  2.5617E+01  2.8374E+01  2.0255E+00  8.4307E+00 -2.5905E+00  3.9188E+00 -4.9976E+00
             9.2860E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1642.27459041177        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0178E+00  1.2008E+00  1.1441E+00  8.7434E-01  1.1555E+00  1.0142E+00  1.2857E+00  1.3266E+00  7.6348E-01  1.2355E+00
             1.1734E+00
 PARAMETER:  1.1763E-01  2.8300E-01  2.3466E-01 -3.4289E-02  2.4449E-01  1.1411E-01  3.5134E-01  3.8259E-01 -1.6987E-01  3.1149E-01
             2.5993E-01
 GRADIENT:   4.8414E+01  1.0173E+01 -1.3285E+01  2.0435E+01  2.4459E+01  5.5360E+00  7.1859E+00  8.7047E-01 -3.2131E+00 -5.3051E+00
             1.1708E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1643.69436566659        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      275
 NPARAMETR:  1.0159E+00  1.3079E+00  1.1992E+00  8.1847E-01  1.1935E+00  1.0153E+00  1.1007E+00  1.4401E+00  9.2357E-01  1.2833E+00
             1.1724E+00
 PARAMETER:  1.1579E-01  3.6845E-01  2.8162E-01 -1.0032E-01  2.7689E-01  1.1518E-01  1.9598E-01  4.6474E-01  2.0488E-02  3.4941E-01
             2.5909E-01
 GRADIENT:   1.5780E+00  7.7145E+00 -3.3426E+00  1.5350E+01  4.8776E+00  1.8571E-01 -1.9700E+00 -6.3023E-01 -4.9917E-01 -1.1025E+00
            -9.0050E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1644.17267195186        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      457
 NPARAMETR:  1.0141E+00  1.4326E+00  1.2776E+00  7.2851E-01  1.2639E+00  1.0149E+00  1.0469E+00  1.8345E+00  9.5406E-01  1.3380E+00
             1.1752E+00
 PARAMETER:  1.1397E-01  4.5947E-01  3.4495E-01 -2.1676E-01  3.3422E-01  1.1482E-01  1.4587E-01  7.0675E-01  5.2967E-02  3.9118E-01
             2.6141E-01
 GRADIENT:  -2.6615E+00 -2.8029E+00  1.3883E+00 -2.5448E+00 -1.8849E+00 -1.4263E-01  1.9906E-01 -1.5776E-01 -1.8011E-01 -1.9992E-01
             3.1378E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1644.17584074933        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  1.0142E+00  1.4332E+00  1.2771E+00  7.2866E-01  1.2644E+00  1.0153E+00  1.0456E+00  1.8358E+00  9.5415E-01  1.3375E+00
             1.1749E+00
 PARAMETER:  1.1408E-01  4.5993E-01  3.4460E-01 -2.1654E-01  3.3456E-01  1.1521E-01  1.4464E-01  7.0746E-01  5.3067E-02  3.9078E-01
             2.6115E-01
 GRADIENT:  -2.4407E+00 -2.1612E+00  1.2693E+00 -1.8166E+00 -1.5522E+00  8.7567E-04  1.7714E-04 -1.2904E-01 -2.5318E-01 -3.0361E-01
             1.6958E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1644.18731512183        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  1.0152E+00  1.4351E+00  1.2661E+00  7.2935E-01  1.2655E+00  1.0153E+00  1.0449E+00  1.8294E+00  9.5770E-01  1.3396E+00
             1.1740E+00
 PARAMETER:  1.1513E-01  4.6125E-01  3.3595E-01 -2.1560E-01  3.3549E-01  1.1515E-01  1.4394E-01  7.0397E-01  5.6783E-02  3.9240E-01
             2.6043E-01
 GRADIENT:  -3.0888E-01 -2.7612E-01  2.3807E-01  1.0820E+00  3.5986E-01 -4.0310E-02  6.7938E-02  1.7747E-01 -4.1068E-02 -3.8190E-02
             3.0679E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1644.19928561707        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  1.0161E+00  1.4400E+00  1.2463E+00  7.2491E-01  1.2627E+00  1.0159E+00  1.0418E+00  1.7931E+00  9.6177E-01  1.3395E+00
             1.1736E+00
 PARAMETER:  1.1599E-01  4.6461E-01  3.2021E-01 -2.2171E-01  3.3324E-01  1.1578E-01  1.4091E-01  6.8393E-01  6.1017E-02  3.9233E-01
             2.6008E-01
 GRADIENT:   1.4765E+00 -1.3453E+00  6.1617E-01 -2.3911E-01 -1.9054E-01  1.9635E-01  4.4046E-02 -1.7720E-01 -6.2979E-02  2.0673E-01
            -9.8757E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1644.20020836773        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1126
 NPARAMETR:  1.0162E+00  1.4405E+00  1.2453E+00  7.2482E-01  1.2628E+00  1.0156E+00  1.0414E+00  1.7928E+00  9.6196E-01  1.3392E+00
             1.1734E+00
 PARAMETER:  1.1605E-01  4.6501E-01  3.1939E-01 -2.2183E-01  3.3332E-01  1.1545E-01  1.4053E-01  6.8378E-01  6.1216E-02  3.9206E-01
             2.5995E-01
 GRADIENT:   1.5913E+00 -9.8680E-01  5.5527E-01  1.2501E-01 -3.1600E-02  6.2117E-02  7.3990E-03 -1.6486E-01 -7.8842E-02  1.5041E-01
            -1.6746E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1644.20363329323        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1284
 NPARAMETR:  1.0154E+00  1.4419E+00  1.2401E+00  7.2435E-01  1.2623E+00  1.0154E+00  1.0407E+00  1.7959E+00  9.6158E-01  1.3376E+00
             1.1737E+00
 PARAMETER:  1.1532E-01  4.6596E-01  3.1517E-01 -2.2248E-01  3.3293E-01  1.1530E-01  1.3989E-01  6.8550E-01  6.0818E-02  3.9089E-01
             2.6017E-01
 GRADIENT:  -1.4153E-02 -5.0720E-01  1.9396E-01  9.3160E-01  2.5903E-01 -1.4245E-03 -6.4467E-02  3.7924E-02 -1.0752E-01  6.0656E-02
            -2.1470E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1644.24960387613        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1466
 NPARAMETR:  1.0169E+00  1.4830E+00  1.1357E+00  6.9536E-01  1.2489E+00  1.0165E+00  1.0208E+00  1.6998E+00  9.8681E-01  1.3194E+00
             1.1732E+00
 PARAMETER:  1.1673E-01  4.9405E-01  2.2723E-01 -2.6333E-01  3.2224E-01  1.1633E-01  1.2059E-01  6.3053E-01  8.6726E-02  3.7721E-01
             2.5972E-01
 GRADIENT:   2.3377E+00 -4.5404E-01  6.9177E-01 -3.1859E-01 -2.4295E+00  2.7746E-01  3.4975E-01 -2.0490E-01 -2.0450E-02 -6.5320E-01
            -1.6750E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1644.25122191810        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1659            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0168E+00  1.4830E+00  1.1347E+00  6.9536E-01  1.2490E+00  1.0157E+00  1.0189E+00  1.7006E+00  9.8712E-01  1.3196E+00
             1.1733E+00
 PARAMETER:  1.1661E-01  4.9406E-01  2.2637E-01 -2.6333E-01  3.2234E-01  1.1562E-01  1.1869E-01  6.3099E-01  8.7038E-02  3.7732E-01
             2.5983E-01
 GRADIENT:   4.4818E+01  3.1684E+01  6.2578E-01  9.3844E+00 -7.6645E-01  5.7793E+00  7.8579E-01 -5.5879E-02  1.9945E-01 -1.2165E-01
             9.2361E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1644.25198012276        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1800
 NPARAMETR:  1.0163E+00  1.4819E+00  1.1344E+00  6.9604E-01  1.2492E+00  1.0156E+00  1.0188E+00  1.7005E+00  9.8717E-01  1.3198E+00
             1.1735E+00
 PARAMETER:  1.1613E-01  4.9334E-01  2.2613E-01 -2.6235E-01  3.2250E-01  1.1547E-01  1.1867E-01  6.3093E-01  8.7091E-02  3.7751E-01
             2.5997E-01
 GRADIENT:   4.3593E+01  3.1348E+01  3.9496E-01  9.5565E+00 -3.0343E-01  5.7199E+00  7.1321E-01  1.3412E-02  2.2975E-01 -8.0006E-02
             2.0115E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1644.25744719548        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1978
 NPARAMETR:  1.0158E+00  1.4834E+00  1.1344E+00  6.9530E-01  1.2518E+00  1.0158E+00  1.0190E+00  1.7004E+00  9.8712E-01  1.3247E+00
             1.1732E+00
 PARAMETER:  1.1569E-01  4.9431E-01  2.2607E-01 -2.6341E-01  3.2459E-01  1.1565E-01  1.1887E-01  6.3084E-01  8.7037E-02  3.8121E-01
             2.5975E-01
 GRADIENT:   1.0925E-01 -8.5783E-01  5.0696E-02  4.6797E-01 -8.4194E-01  1.0716E-02  1.3561E-02 -1.0654E-01 -3.0805E-02 -1.9508E-01
            -2.8218E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1644.25763743691        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     2138
 NPARAMETR:  1.0159E+00  1.4837E+00  1.1342E+00  6.9539E-01  1.2520E+00  1.0157E+00  1.0186E+00  1.7009E+00  9.8813E-01  1.3245E+00
             1.1731E+00
 PARAMETER:  1.1575E-01  4.9456E-01  2.2596E-01 -2.6328E-01  3.2475E-01  1.1563E-01  1.1843E-01  6.3115E-01  8.8055E-02  3.8102E-01
             2.5962E-01
 GRADIENT:   2.2281E-01 -4.8322E-01 -1.3025E-02  8.4489E-01 -6.7930E-01 -7.2012E-04 -8.7773E-04 -8.8989E-02  2.1863E-04 -2.3997E-01
            -8.8362E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1644.26086890916        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2329
 NPARAMETR:  1.0158E+00  1.4849E+00  1.1344E+00  6.9392E-01  1.2543E+00  1.0158E+00  1.0187E+00  1.7097E+00  9.8779E-01  1.3279E+00
             1.1731E+00
 PARAMETER:  1.1564E-01  4.9532E-01  2.2613E-01 -2.6539E-01  3.2662E-01  1.1563E-01  1.1852E-01  6.3634E-01  8.7711E-02  3.8358E-01
             2.5969E-01
 GRADIENT:  -8.4902E-03 -1.7926E+00 -3.0561E-01  3.1858E-01 -7.8324E-03  2.3273E-03  8.3838E-02  2.2024E-02  5.6384E-02  8.4564E-03
             2.7838E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1644.26176817557        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2514             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0158E+00  1.4857E+00  1.1350E+00  6.9394E-01  1.2545E+00  1.0158E+00  1.0183E+00  1.7104E+00  9.8755E-01  1.3276E+00
             1.1732E+00
 PARAMETER:  1.1570E-01  4.9587E-01  2.2662E-01 -2.6536E-01  3.2677E-01  1.1563E-01  1.1814E-01  6.3673E-01  8.7471E-02  3.8338E-01
             2.5971E-01
 GRADIENT:   4.2509E+01  3.1562E+01 -1.6271E-01  1.0339E+01  1.3859E+00  5.7957E+00  8.0197E-01  8.6138E-02  2.9411E-01  4.4063E-01
             2.0185E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1644.26225792037        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2692
 NPARAMETR:  1.0158E+00  1.4859E+00  1.1356E+00  6.9376E-01  1.2546E+00  1.0158E+00  1.0181E+00  1.7105E+00  9.8772E-01  1.3277E+00
             1.1732E+00
 PARAMETER:  1.1569E-01  4.9604E-01  2.2717E-01 -2.6563E-01  3.2679E-01  1.1564E-01  1.1791E-01  6.3678E-01  8.7648E-02  3.8345E-01
             2.5973E-01
 GRADIENT:   8.6409E-02 -9.6172E-01 -1.5659E-01  6.9996E-01 -2.2112E-01  1.3164E-03  1.8983E-02 -4.3530E-02 -1.1334E-02 -5.6497E-02
            -1.5566E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1644.26305271952        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2878
 NPARAMETR:  1.0158E+00  1.4866E+00  1.1362E+00  6.9351E-01  1.2549E+00  1.0158E+00  1.0178E+00  1.7119E+00  9.8803E-01  1.3277E+00
             1.1731E+00
 PARAMETER:  1.1570E-01  4.9647E-01  2.2768E-01 -2.6599E-01  3.2703E-01  1.1564E-01  1.1766E-01  6.3762E-01  8.7960E-02  3.8347E-01
             2.5966E-01
 GRADIENT:   1.1599E-01 -7.1605E-01 -1.1093E-01  8.0691E-01 -2.5269E-01  9.0592E-04  2.5964E-02 -6.0393E-02 -1.6090E-02 -9.2326E-02
            -6.7743E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1644.26354993451        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     3059
 NPARAMETR:  1.0159E+00  1.4871E+00  1.1364E+00  6.9339E-01  1.2552E+00  1.0158E+00  1.0172E+00  1.7133E+00  9.8891E-01  1.3277E+00
             1.1730E+00
 PARAMETER:  1.1574E-01  4.9686E-01  2.2787E-01 -2.6616E-01  3.2729E-01  1.1565E-01  1.1710E-01  6.3844E-01  8.8847E-02  3.8343E-01
             2.5958E-01
 GRADIENT:   1.9619E-01 -4.1496E-01 -1.3577E-01  1.0742E+00 -1.6586E-01 -8.8497E-05  9.8677E-04 -5.2081E-02 -1.2401E-04 -1.2539E-01
            -1.1448E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1644.26528965589        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3246             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0159E+00  1.4883E+00  1.1371E+00  6.9267E-01  1.2559E+00  1.0158E+00  1.0167E+00  1.7173E+00  9.8892E-01  1.3282E+00
             1.1731E+00
 PARAMETER:  1.1573E-01  4.9762E-01  2.2852E-01 -2.6720E-01  3.2784E-01  1.1565E-01  1.1658E-01  6.4074E-01  8.8860E-02  3.8383E-01
             2.5964E-01
 GRADIENT:   4.2599E+01  3.2379E+01 -1.0147E-02  1.0592E+01  1.2467E+00  5.7974E+00  7.3955E-01  4.2129E-02  2.6468E-01  3.7630E-01
             1.0258E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1644.26571742513        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3428
 NPARAMETR:  1.0158E+00  1.4884E+00  1.1374E+00  6.9244E-01  1.2560E+00  1.0158E+00  1.0167E+00  1.7179E+00  9.8931E-01  1.3284E+00
             1.1732E+00
 PARAMETER:  1.1572E-01  4.9768E-01  2.2875E-01 -2.6753E-01  3.2790E-01  1.1565E-01  1.1656E-01  6.4109E-01  8.9257E-02  3.8399E-01
             2.5970E-01
 GRADIENT:   1.3570E-01 -6.1869E-01 -6.7226E-02  8.2742E-01 -2.6975E-01  1.8844E-03  1.4640E-02 -5.7989E-02 -2.6220E-03 -9.9902E-02
            -7.4996E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1644.26744217151        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3612
 NPARAMETR:  1.0158E+00  1.4897E+00  1.1378E+00  6.9173E-01  1.2568E+00  1.0158E+00  1.0159E+00  1.7221E+00  9.8939E-01  1.3289E+00
             1.1732E+00
 PARAMETER:  1.1572E-01  4.9860E-01  2.2912E-01 -2.6856E-01  3.2857E-01  1.1565E-01  1.1580E-01  6.4353E-01  8.9337E-02  3.8435E-01
             2.5970E-01
 GRADIENT:   1.4038E-01 -4.1283E-01 -8.1650E-02  9.9917E-01 -1.7804E-01 -1.2223E-03 -4.0526E-02 -4.8183E-02 -3.7979E-02 -1.1650E-01
            -9.5243E-02

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1644.28360559366        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3804             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0150E+00  1.5060E+00  1.1412E+00  6.8101E-01  1.2666E+00  1.0158E+00  1.0093E+00  1.7785E+00  9.9032E-01  1.3396E+00
             1.1754E+00
 PARAMETER:  1.1491E-01  5.0945E-01  2.3212E-01 -2.8418E-01  3.3635E-01  1.1571E-01  1.0930E-01  6.7575E-01  9.0269E-02  3.9240E-01
             2.6157E-01
 GRADIENT:   4.0288E+01  3.3226E+01 -1.7833E-01  1.0556E+01  1.5177E+00  5.7979E+00  5.8462E-01  3.1956E-01  8.1706E-02  9.6806E-01
             8.8157E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1644.28490463590        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3986            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0151E+00  1.5061E+00  1.1417E+00  6.8092E-01  1.2666E+00  1.0158E+00  1.0088E+00  1.7770E+00  9.9434E-01  1.3392E+00
             1.1751E+00
 PARAMETER:  1.1501E-01  5.0955E-01  2.3249E-01 -2.8431E-01  3.3633E-01  1.1572E-01  1.0878E-01  6.7490E-01  9.4322E-02  3.9206E-01
             2.6136E-01
 GRADIENT:   4.0573E+01  3.3237E+01 -1.1974E-01  1.0480E+01  1.4952E+00  5.8035E+00  7.2264E-01  2.9745E-01  2.6100E-01  9.2417E-01
             8.0495E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1644.28581436235        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     4145
 NPARAMETR:  1.0156E+00  1.5065E+00  1.1415E+00  6.8102E-01  1.2665E+00  1.0156E+00  1.0087E+00  1.7776E+00  9.9533E-01  1.3389E+00
             1.1732E+00
 PARAMETER:  1.1546E-01  5.0980E-01  2.3237E-01 -2.8416E-01  3.3624E-01  1.1550E-01  1.0867E-01  6.7524E-01  9.5322E-02  3.9187E-01
             2.5972E-01
 GRADIENT:  -5.0483E-01 -5.1737E-01 -1.8893E-01  9.8574E-01 -3.0758E-03 -8.9235E-02  3.1619E-02  1.9745E-01  2.2731E-03  3.5018E-01
            -1.2202E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1644.28665700973        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4326
 NPARAMETR:  1.0158E+00  1.5066E+00  1.1420E+00  6.8079E-01  1.2665E+00  1.0159E+00  1.0085E+00  1.7752E+00  9.9544E-01  1.3383E+00
             1.1735E+00
 PARAMETER:  1.1571E-01  5.0989E-01  2.3281E-01 -2.8449E-01  3.3625E-01  1.1573E-01  1.0846E-01  6.7389E-01  9.5425E-02  3.9138E-01
             2.6002E-01
 GRADIENT:   1.5059E-02 -6.2468E-01 -3.7009E-02  7.1974E-01 -7.6017E-02  1.8146E-03  1.9132E-02  1.1134E-01 -2.5686E-02  2.4546E-01
            -4.2341E-02

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1644.28718255103        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4507
 NPARAMETR:  1.0159E+00  1.5071E+00  1.1418E+00  6.8062E-01  1.2667E+00  1.0159E+00  1.0077E+00  1.7742E+00  9.9695E-01  1.3376E+00
             1.1735E+00
 PARAMETER:  1.1573E-01  5.1020E-01  2.3258E-01 -2.8475E-01  3.3638E-01  1.1574E-01  1.0771E-01  6.7335E-01  9.6944E-02  3.9091E-01
             2.6000E-01
 GRADIENT:   6.2594E-02 -4.7571E-01 -2.4929E-02  8.4471E-01  4.6544E-02 -1.1496E-04 -9.4135E-03  9.0062E-02  3.9196E-03  1.5173E-01
            -6.7556E-02

0ITERATION NO.:  137    OBJECTIVE VALUE:  -1644.28718360474        NO. OF FUNC. EVALS.:  86
 CUMULATIVE NO. OF FUNC. EVALS.:     4593
 NPARAMETR:  1.0159E+00  1.5072E+00  1.1418E+00  6.8043E-01  1.2666E+00  1.0159E+00  1.0079E+00  1.7732E+00  9.9679E-01  1.3374E+00
             1.1735E+00
 PARAMETER:  1.1573E-01  5.1020E-01  2.3258E-01 -2.8475E-01  3.3638E-01  1.1574E-01  1.0783E-01  6.7335E-01  9.6787E-02  3.9091E-01
             2.6000E-01
 GRADIENT:   4.1421E-02 -2.4764E-01 -1.6030E-02  7.9139E-01  4.4870E-02 -2.1262E-03  2.9171E-03  8.2546E-02  9.7670E-04  1.3120E-01
            -3.8955E-02
 NUMSIGDIG:         3.8         3.2         2.9         2.1         3.8         4.4         3.9         2.2         3.8         2.4
                    3.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4593
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.0729E-04 -1.1992E-02 -3.6274E-02  4.7342E-03 -4.2036E-02
 SE:             2.9797E-02  2.4082E-02  1.2133E-02  1.8217E-02  2.2249E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9713E-01  6.1851E-01  2.7925E-03  7.9496E-01  5.8848E-02

 ETASHRINKSD(%)  1.7641E-01  1.9324E+01  5.9353E+01  3.8970E+01  2.5463E+01
 ETASHRINKVR(%)  3.5251E-01  3.4913E+01  8.3479E+01  6.2754E+01  4.4443E+01
 EBVSHRINKSD(%)  5.6197E-01  1.8861E+01  6.6404E+01  4.2219E+01  2.0708E+01
 EBVSHRINKVR(%)  1.1208E+00  3.4165E+01  8.8713E+01  6.6613E+01  3.7128E+01
 RELATIVEINF(%)  9.8358E+01  1.6111E+00  1.0913E+00  7.5173E-01  2.2623E+01
 EPSSHRINKSD(%)  4.3635E+01
 EPSSHRINKVR(%)  6.8229E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1644.2871836047375     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -909.13635704099931     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    67.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1644.287       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.51E+00  1.14E+00  6.81E-01  1.27E+00  1.02E+00  1.01E+00  1.77E+00  9.97E-01  1.34E+00  1.17E+00
 


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
+        1.03E+03
 
 TH 2
+       -9.44E+00  3.29E+02
 
 TH 3
+        3.57E+00  1.97E+01  4.29E+01
 
 TH 4
+       -1.12E+01  1.06E+08 -7.03E+01  4.20E+08
 
 TH 5
+       -2.36E+00 -5.26E+01 -5.89E+01  8.64E+01  2.51E+02
 
 TH 6
+       -2.72E+00 -2.00E+00  5.65E-01 -2.81E+00  1.85E-01  1.90E+02
 
 TH 7
+        3.52E-01  1.68E+01  4.01E-01 -2.08E+01 -4.02E+00 -1.29E+00  8.98E+01
 
 TH 8
+       -4.06E-01  1.72E+07 -1.09E+01  6.82E+07  1.50E-02  2.17E-01  9.39E-01  1.11E+07
 
 TH 9
+       -2.23E+00 -9.05E+00 -4.49E+00  2.03E+00  3.90E+00  2.09E+00  3.35E+01  2.64E+00  3.08E+01
 
 TH10
+        5.99E-01 -3.78E+00 -5.24E+00 -1.56E+08 -3.86E+01  3.27E-01 -1.54E+00 -2.53E+07  2.82E+00  5.11E+01
 
 TH11
+       -7.50E+00 -1.04E+01 -8.35E+00  1.47E+00 -2.08E+00  2.94E+00  6.86E+00  2.55E+00  4.66E+00  9.51E+00  1.58E+02
 
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
 #CPUT: Total CPU Time in Seconds,       73.880
Stop Time:
Sat Sep 18 12:47:11 CDT 2021
