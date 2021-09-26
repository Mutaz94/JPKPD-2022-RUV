Sat Sep 25 02:04:20 CDT 2021
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
$DATA ../../../../data/int/SL3/dat25.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      983
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

 TOT. NO. OF OBS RECS:      883
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   253.248284231493        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.6988E+01 -5.6546E+00  7.1555E+01 -2.7043E+01  7.9403E+01  3.2616E+01 -7.7922E+01 -1.7152E+02 -1.0186E+02 -2.3786E+01
            -7.7693E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2348.79731345902        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0652E+00  1.2769E+00  1.1296E+00  1.0132E+00  1.2054E+00  6.9137E-01  9.6910E-01  9.4433E-01  9.1457E-01  7.4124E-01
             5.2391E+00
 PARAMETER:  1.6317E-01  3.4444E-01  2.2187E-01  1.1308E-01  2.8685E-01 -2.6909E-01  6.8614E-02  4.2718E-02  1.0697E-02 -1.9943E-01
             1.7561E+00
 GRADIENT:   1.7050E+01 -1.8693E+01 -3.1199E+01  3.9492E+01  5.3156E+01 -7.5689E+01  7.7207E+00  4.6767E+00  9.1132E+00  9.9101E+00
             7.7907E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2427.30596708091        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0198E+00  1.0106E+00  1.3271E+00  1.1050E+00  1.0563E+00  9.4089E-01  1.4953E+00  1.6577E+00  6.6227E-01  1.4358E-01
             4.7154E+00
 PARAMETER:  1.1962E-01  1.1050E-01  3.8303E-01  1.9981E-01  1.5476E-01  3.9073E-02  5.0233E-01  6.0546E-01 -3.1209E-01 -1.8408E+00
             1.6508E+00
 GRADIENT:  -4.6630E+01 -1.2525E+01 -1.0052E+01 -4.9082E+01 -9.3246E+00  2.5436E+01  3.4780E+01  1.6926E+01 -5.1443E+00  3.6243E-01
             6.9460E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2624.25696032776        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0027E+00  1.2521E+00  1.5808E+00  9.4195E-01  1.3383E+00  9.2339E-01  1.1004E+00  1.3823E+00  9.3936E-01  1.0643E+00
             2.9587E+00
 PARAMETER:  1.0269E-01  3.2481E-01  5.5796E-01  4.0196E-02  3.9142E-01  2.0300E-02  1.9567E-01  4.2372E-01  3.7448E-02  1.6233E-01
             1.1847E+00
 GRADIENT:  -1.7225E+01 -2.6589E+01 -8.1537E+00  5.7025E+00  2.8974E+01  1.0828E+01  1.6217E+01  4.7200E+00 -2.5331E+00  4.9995E+00
             2.6004E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2631.72832612584        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0071E+00  1.4757E+00  1.1219E+00  7.9504E-01  1.3246E+00  8.8854E-01  7.9250E-01  6.2736E-01  1.2087E+00  1.1536E+00
             2.9046E+00
 PARAMETER:  1.0711E-01  4.8915E-01  2.1502E-01 -1.2936E-01  3.8109E-01 -1.8171E-02 -1.3256E-01 -3.6623E-01  2.8957E-01  2.4293E-01
             1.1663E+00
 GRADIENT:  -6.0062E+00 -1.2903E+00 -1.4751E+00  1.6216E+00 -9.0879E+00 -3.5440E+00 -9.7292E-01  1.4985E+00  3.1763E+00  5.2966E+00
            -6.5335E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2633.23496102070        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0098E+00  1.6106E+00  1.1166E+00  7.1229E-01  1.4380E+00  8.9456E-01  7.6885E-01  2.7230E-01  1.2607E+00  1.1749E+00
             2.9169E+00
 PARAMETER:  1.0975E-01  5.7662E-01  2.1028E-01 -2.3927E-01  4.6325E-01 -1.1421E-02 -1.6285E-01 -1.2008E+00  3.3165E-01  2.6118E-01
             1.1705E+00
 GRADIENT:   1.6427E-01  7.5469E+00  8.5322E-01  1.7486E+00 -1.0204E+00 -1.0169E+00 -1.6830E+00  1.4732E-01 -1.3439E+00 -2.4212E+00
             2.7915E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2633.60825815809        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      485
 NPARAMETR:  1.0112E+00  1.6920E+00  1.0947E+00  6.6678E-01  1.4989E+00  8.9658E-01  7.5530E-01  1.6773E-01  1.3231E+00  1.2358E+00
             2.9145E+00
 PARAMETER:  1.1113E-01  6.2592E-01  1.9046E-01 -3.0530E-01  5.0472E-01 -9.1685E-03 -1.8064E-01 -1.6854E+00  3.7995E-01  3.1174E-01
             1.1697E+00
 GRADIENT:  -2.1039E+00  6.5462E+00  7.6786E-02  7.4280E+00 -1.2110E+00 -7.3979E-01 -9.0031E-01  4.5838E-02 -6.3688E-01  9.2590E-01
             1.2527E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2634.51288532521        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  1.0111E+00  1.9994E+00  8.8324E-01  4.5747E-01  1.7120E+00  8.9747E-01  7.1175E-01  1.0000E-02  1.6455E+00  1.3644E+00
             2.9019E+00
 PARAMETER:  1.1099E-01  7.9282E-01 -2.4159E-02 -6.8205E-01  6.3767E-01 -8.1775E-03 -2.4003E-01 -5.8987E+00  5.9805E-01  4.1074E-01
             1.1654E+00
 GRADIENT:  -2.2190E+00 -2.3527E+00 -2.4090E-01  3.3646E-01 -4.2854E-01 -5.4837E-01 -3.6845E-01  0.0000E+00 -4.7093E-01 -4.8713E-01
             4.7083E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2635.11997846572        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      837
 NPARAMETR:  1.0132E+00  2.3486E+00  5.6862E-01  2.3845E-01  1.9518E+00  9.0009E-01  6.8503E-01  1.0000E-02  2.3207E+00  1.5226E+00
             2.8921E+00
 PARAMETER:  1.1311E-01  9.5382E-01 -4.6454E-01 -1.3336E+00  7.6874E-01 -5.2625E-03 -2.7829E-01 -1.5051E+01  9.4189E-01  5.2041E-01
             1.1620E+00
 GRADIENT:   2.4579E+00  2.4496E+01  1.5803E-01  5.0715E+00 -2.0833E-01  1.7206E-01 -1.9337E+00  0.0000E+00  1.2198E+00 -2.9911E-01
            -9.4206E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2635.67254965156        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  1.0157E+00  2.5368E+00  3.3846E-01  1.1267E-01  2.1083E+00  9.0667E-01  6.8674E-01  1.0000E-02  3.1295E+00  1.6285E+00
             2.8798E+00
 PARAMETER:  1.1559E-01  1.0309E+00 -9.8334E-01 -2.0833E+00  8.4586E-01  2.0199E-03 -2.7580E-01 -2.6960E+01  1.2409E+00  5.8767E-01
             1.1577E+00
 GRADIENT:   9.3310E+00  2.5855E+01 -2.5461E+00  3.1152E+00  7.0540E+00  2.7242E+00  3.5192E+00  0.0000E+00 -3.1791E+00  1.3426E-01
            -6.4362E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2636.55600770446        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1190
 NPARAMETR:  1.0073E+00  2.5809E+00  2.9303E-01  6.5595E-02  2.1511E+00  8.9268E-01  6.6569E-01  1.0000E-02  4.1789E+00  1.6820E+00
             2.8847E+00
 PARAMETER:  1.0725E-01  1.0481E+00 -1.1275E+00 -2.6243E+00  8.6600E-01 -1.3530E-02 -3.0693E-01 -3.5664E+01  1.5301E+00  6.1995E-01
             1.1594E+00
 GRADIENT:  -1.2059E+01 -1.4156E+01 -2.4945E+00  1.2221E+00  4.5106E+00 -2.6037E+00 -1.0615E+00  0.0000E+00  1.0304E+00  1.3756E+00
             7.5407E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2637.18248454927        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1368
 NPARAMETR:  1.0107E+00  2.5799E+00  3.7139E-01  6.8093E-02  2.1429E+00  8.9764E-01  6.6510E-01  1.0000E-02  4.1262E+00  1.6796E+00
             2.8813E+00
 PARAMETER:  1.1068E-01  1.0477E+00 -8.9051E-01 -2.5869E+00  8.6214E-01 -7.9893E-03 -3.0782E-01 -3.4685E+01  1.5174E+00  6.1853E-01
             1.1582E+00
 GRADIENT:  -3.1303E+00 -4.7365E+00 -1.0535E+00  1.8303E+00  1.8128E+00 -7.3231E-01 -2.2589E+00  0.0000E+00 -1.9283E+00  2.1342E+00
            -3.1581E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2637.99411975087        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1547
 NPARAMETR:  1.0120E+00  2.6242E+00  2.8878E-01  3.9218E-02  2.1580E+00  9.0053E-01  6.7020E-01  1.0000E-02  5.4033E+00  1.6830E+00
             2.8773E+00
 PARAMETER:  1.1196E-01  1.0648E+00 -1.1421E+00 -3.1386E+00  8.6917E-01 -4.7695E-03 -3.0018E-01 -4.3804E+01  1.7870E+00  6.2059E-01
             1.1569E+00
 GRADIENT:   2.4390E-02 -9.8348E-01 -7.4058E-01  9.2353E-01  1.3129E-01  1.0174E-01 -3.8903E-01  0.0000E+00  6.7026E-01  1.5338E-01
            -9.5305E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2638.63917490284        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1731
 NPARAMETR:  1.0118E+00  2.6569E+00  3.4371E-01  2.2694E-02  2.1691E+00  8.9746E-01  6.7998E-01  1.0000E-02  6.9893E+00  1.6662E+00
             2.8807E+00
 PARAMETER:  1.1177E-01  1.0772E+00 -9.6797E-01 -3.6857E+00  8.7430E-01 -8.1851E-03 -2.8569E-01 -5.2147E+01  2.0444E+00  6.1055E-01
             1.1580E+00
 GRADIENT:  -1.1342E+00  3.1050E+01  2.6896E-01 -2.8220E+00  4.2547E+00 -1.2883E+00  9.5361E+00  0.0000E+00 -8.7965E+00  1.1404E+00
             3.8378E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2638.95040077467        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1898
 NPARAMETR:  1.0097E+00  2.6608E+00  3.5676E-01  1.8273E-02  2.1474E+00  8.9943E-01  6.7294E-01  1.0000E-02  7.7559E+00  1.6692E+00
             2.8791E+00
 PARAMETER:  1.0967E-01  1.0786E+00 -9.3069E-01 -3.9023E+00  8.6426E-01 -5.9945E-03 -2.9609E-01 -5.5553E+01  2.1484E+00  6.1233E-01
             1.1575E+00
 GRADIENT:  -7.7205E-01  7.1969E+01  5.4571E-01 -2.2555E+00  1.7112E-01 -6.2352E-02  7.2250E+00  0.0000E+00 -5.9355E+00  2.5391E+00
             4.8831E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2639.04430067616        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1974
 NPARAMETR:  1.0161E+00  2.6398E+00  3.5673E-01  1.8335E-02  2.1306E+00  9.0271E-01  6.7290E-01  1.0000E-02  7.7765E+00  1.6690E+00
             2.8773E+00
 PARAMETER:  1.1596E-01  1.0707E+00 -9.3077E-01 -3.8989E+00  8.5642E-01 -2.3504E-03 -2.9616E-01 -5.5553E+01  2.1511E+00  6.1222E-01
             1.1568E+00
 GRADIENT:   1.7003E+01  3.8645E+01  3.9806E-01 -2.5281E+00 -3.7835E+00  1.5292E+00  7.4659E+00  0.0000E+00 -5.7825E+00  2.8090E+00
             4.1022E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2639.09909501929        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     2131
 NPARAMETR:  1.0119E+00  2.6395E+00  3.5673E-01  1.8338E-02  2.1513E+00  9.0000E-01  6.7289E-01  1.0000E-02  7.7772E+00  1.6690E+00
             2.8772E+00
 PARAMETER:  1.1181E-01  1.0706E+00 -9.3078E-01 -3.8988E+00  8.6606E-01 -5.3646E-03 -2.9617E-01 -5.5553E+01  2.1512E+00  6.1221E-01
             1.1568E+00
 GRADIENT:  -5.4828E-02 -4.1033E+00  4.0797E-01 -2.8267E+00  1.1665E-02 -4.2735E-02  6.8251E+00  0.0000E+00 -7.7003E+00  2.0247E+00
             2.6004E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2639.37098114981        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     2226
 NPARAMETR:  1.0110E+00  2.6355E+00  3.4718E-01  1.7886E-02  2.1493E+00  8.9969E-01  6.6748E-01  1.0000E-02  8.1075E+00  1.6548E+00
             2.8743E+00
 PARAMETER:  1.1095E-01  1.0691E+00 -9.5791E-01 -3.9237E+00  8.6513E-01 -5.7001E-03 -3.0424E-01 -5.5553E+01  2.1928E+00  6.0367E-01
             1.1558E+00
 GRADIENT:   3.6896E+00  2.5440E+01 -7.8066E-03 -6.3793E-01  2.1611E+00  2.9077E-01  2.4577E+00  0.0000E+00 -1.5705E+00  4.1377E-01
             1.4200E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2639.47706252218        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2388
 NPARAMETR:  1.0107E+00  2.6357E+00  3.5748E-01  1.5470E-02  2.1443E+00  8.9941E-01  6.6557E-01  1.0000E-02  8.7047E+00  1.6544E+00
             2.8747E+00
 PARAMETER:  1.1067E-01  1.0692E+00 -9.2867E-01 -4.0689E+00  8.6281E-01 -6.0213E-03 -3.0712E-01 -5.5553E+01  2.2639E+00  6.0345E-01
             1.1559E+00
 GRADIENT:  -2.7457E+00 -1.6897E+01  2.6260E-01 -2.0411E+00 -1.3751E+00 -2.2395E-01  3.1070E+00  0.0000E+00 -5.4349E+00  1.1943E-01
             2.3640E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2639.48747392112        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     2556
 NPARAMETR:  1.0117E+00  2.6356E+00  3.5160E-01  1.5474E-02  2.1495E+00  8.9990E-01  6.6430E-01  1.0000E-02  8.7060E+00  1.6542E+00
             2.8749E+00
 PARAMETER:  1.1167E-01  1.0691E+00 -9.4526E-01 -4.0686E+00  8.6525E-01 -5.4745E-03 -3.0902E-01 -5.5553E+01  2.2640E+00  6.0332E-01
             1.1560E+00
 GRADIENT:  -1.7423E-01 -2.1000E+01  2.9415E-02 -8.3611E-01 -4.0748E-03 -1.1514E-02  6.8463E-01  0.0000E+00 -3.1163E+00 -2.6888E-02
            -1.7096E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2639.49209480439        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     2693
 NPARAMETR:  1.0119E+00  2.6371E+00  3.4531E-01  1.5444E-02  2.1489E+00  9.0010E-01  6.6359E-01  1.0000E-02  8.6968E+00  1.6547E+00
             2.8744E+00
 PARAMETER:  1.1183E-01  1.0697E+00 -9.6330E-01 -4.0706E+00  8.6497E-01 -5.2504E-03 -3.1009E-01 -5.5553E+01  2.2630E+00  6.0362E-01
             1.1558E+00
 GRADIENT:   2.2362E-01 -1.6198E+01 -4.8475E-02 -1.8315E+00 -1.7657E-01  6.3107E-02  2.0083E+00  0.0000E+00 -5.0577E+00 -1.6053E-02
            -3.9320E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2639.50251695863        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2870
 NPARAMETR:  1.0119E+00  2.6384E+00  3.5713E-01  1.5451E-02  2.1490E+00  8.9984E-01  6.6154E-01  1.0000E-02  8.7009E+00  1.6534E+00
             2.8751E+00
 PARAMETER:  1.1188E-01  1.0702E+00 -9.2966E-01 -4.0701E+00  8.6500E-01 -5.5330E-03 -3.1318E-01 -5.5553E+01  2.2634E+00  6.0283E-01
             1.1561E+00
 GRADIENT:   2.3425E-01 -1.6986E+01  1.7804E-01 -4.9644E-01 -2.9830E-01 -5.1771E-02 -1.0959E+00  0.0000E+00 -2.5714E+00 -7.4256E-02
            -3.2113E-01

0ITERATION NO.:  107    OBJECTIVE VALUE:  -2639.50377616909        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     2933
 NPARAMETR:  1.0119E+00  2.6380E+00  3.5764E-01  1.5466E-02  2.1496E+00  8.9992E-01  6.6155E-01  1.0000E-02  8.7061E+00  1.6537E+00
             2.8740E+00
 PARAMETER:  1.1188E-01  1.0703E+00 -9.2793E-01 -4.0700E+00  8.6489E-01 -5.5517E-03 -3.1326E-01 -5.5553E+01  2.2635E+00  6.0274E-01
             1.1561E+00
 GRADIENT:   1.5010E-01  2.7898E+03  1.9382E+03 -7.4123E+02 -2.2399E-01 -4.6715E-02 -9.3461E+02  0.0000E+00 -1.3329E+03 -5.1303E-02
             3.1011E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2933
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0090E-03 -6.9256E-03  1.4232E-06  1.2746E-02 -1.7427E-02
 SE:             2.9147E-02  2.7694E-02  3.9744E-06  7.5772E-03  2.6262E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7238E-01  8.0253E-01  7.2028E-01  9.2552E-02  5.0696E-01

 ETASHRINKSD(%)  2.3535E+00  7.2208E+00  9.9987E+01  7.4615E+01  1.2018E+01
 ETASHRINKVR(%)  4.6516E+00  1.3920E+01  1.0000E+02  9.3556E+01  2.2591E+01
 EBVSHRINKSD(%)  2.2765E+00  6.7910E+00  9.9939E+01  8.2270E+01  1.0610E+01
 EBVSHRINKVR(%)  4.5012E+00  1.3121E+01  1.0000E+02  9.6857E+01  2.0094E+01
 RELATIVEINF(%)  9.5404E+01  5.4065E+01  3.6117E-05  1.9825E+00  7.8885E+01
 EPSSHRINKSD(%)  1.5332E+01
 EPSSHRINKVR(%)  2.8313E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          883
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1622.8454496394520     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2639.5037761690924     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1016.6583265296404     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    72.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2639.504       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  2.64E+00  3.58E-01  1.55E-02  2.15E+00  9.00E-01  6.61E-01  1.00E-02  8.70E+00  1.65E+00  2.88E+00
 


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
+        1.30E+03
 
 TH 2
+        1.55E+02  1.31E+05
 
 TH 3
+        1.45E+03 -9.79E+01  8.16E+06
 
 TH 4
+       -7.75E+03  1.57E+04  4.25E+03  2.64E+08
 
 TH 5
+       -3.52E+00 -8.66E+01  1.46E+06  3.20E+03  7.76E+01
 
 TH 6
+        1.22E+01  9.33E+01  9.06E+02 -4.43E+03 -9.52E-01  2.55E+02
 
 TH 7
+        1.34E+02  2.61E+02  1.31E+07 -1.33E+04  3.80E+05 -6.03E+01  2.43E+07
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.47E+01  4.99E+01  1.40E+01  1.17E+05  1.04E+01 -1.42E+01 -4.12E+01  0.00E+00  2.69E+03
 
 TH10
+       -5.75E-01 -3.83E+02  2.72E+06  1.72E+04 -6.57E+00  1.19E+00 -2.47E+02  0.00E+00  5.49E+01  4.40E+01
 
 TH11
+        1.25E+02 -6.14E+02  8.12E+05  2.69E+04  1.69E+05  9.36E+01  1.51E+06  0.00E+00  8.60E+01  2.70E+05  9.43E+04
 
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
 #CPUT: Total CPU Time in Seconds,       86.301
Stop Time:
Sat Sep 25 02:05:48 CDT 2021
