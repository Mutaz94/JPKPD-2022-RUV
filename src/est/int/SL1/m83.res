Sat Sep 18 03:42:29 CDT 2021
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
$DATA ../../../../data/int/SL1/dat83.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3044.15261536690        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9540E+01  2.1466E+01  9.0820E+01  4.6530E+01 -4.3861E-01 -3.8982E+00 -3.9863E+01 -1.1760E+02 -4.3125E+01  2.7939E+01
            -1.4589E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3337.84817692088        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0104E+00  1.1039E+00  1.0279E+00  9.8517E-01  1.0786E+00  1.0176E+00  9.9699E-01  1.2619E+00  1.0745E+00  7.7059E-01
             1.4417E+00
 PARAMETER:  1.1035E-01  1.9883E-01  1.2756E-01  8.5062E-02  1.7565E-01  1.1749E-01  9.6990E-02  3.3263E-01  1.7188E-01 -1.6060E-01
             4.6581E-01
 GRADIENT:   2.3291E+01  5.0735E+01  2.1169E+00  8.2750E+01  1.4312E+01  2.4018E+00 -1.3080E+01 -1.8412E+01 -3.8153E+00 -1.3200E+01
            -2.0349E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3341.10809060413        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0144E+00  1.1243E+00  1.1588E+00  9.5659E-01  1.1336E+00  1.0320E+00  1.0914E+00  1.7875E+00  1.0108E+00  7.9658E-01
             1.4267E+00
 PARAMETER:  1.1430E-01  2.1712E-01  2.4737E-01  5.5623E-02  2.2535E-01  1.3147E-01  1.8744E-01  6.8081E-01  1.1074E-01 -1.2742E-01
             4.5540E-01
 GRADIENT:   3.3576E+01  4.2750E+01 -5.1076E+00  4.3786E+01  6.1570E-01  8.0913E+00  3.2852E+00 -8.2528E+00 -1.0959E+01 -1.1260E+01
            -1.9354E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3341.18698880198        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      326             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0145E+00  1.1252E+00  1.1596E+00  9.5618E-01  1.1347E+00  1.0048E+00  1.0898E+00  1.7908E+00  1.0104E+00  7.9776E-01
             1.4271E+00
 PARAMETER:  1.1439E-01  2.1796E-01  2.4809E-01  5.5195E-02  2.2635E-01  1.0481E-01  1.8596E-01  6.8264E-01  1.1030E-01 -1.2594E-01
             4.5565E-01
 GRADIENT:   3.4092E+01  4.2871E+01 -5.3089E+00  4.3958E+01  9.2081E-01 -2.9110E+00  3.2472E+00 -8.2391E+00 -1.0984E+01 -1.1338E+01
            -1.9305E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3341.22377188034        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      485
 NPARAMETR:  1.0145E+00  1.1252E+00  1.1596E+00  9.5618E-01  1.1347E+00  1.0215E+00  1.0898E+00  1.7908E+00  1.0104E+00  7.9776E-01
             1.4271E+00
 PARAMETER:  1.1439E-01  2.1796E-01  2.4809E-01  5.5195E-02  2.2635E-01  1.2127E-01  1.8596E-01  6.8264E-01  1.1030E-01 -1.2594E-01
             4.5567E-01
 GRADIENT:   7.8576E+00  3.2885E+01 -6.0272E+00  3.9623E+01 -4.2173E+00  1.2838E+00  2.4257E+00 -8.7599E+00 -1.1474E+01 -1.1506E+01
            -1.9381E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3341.57139740041        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      678             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0145E+00  1.1251E+00  1.1597E+00  9.5615E-01  1.1346E+00  1.0170E+00  1.0898E+00  1.7912E+00  1.0104E+00  7.9779E-01
             1.4297E+00
 PARAMETER:  1.1436E-01  2.1787E-01  2.4816E-01  5.5161E-02  2.2630E-01  1.1690E-01  1.8600E-01  6.8291E-01  1.1033E-01 -1.2591E-01
             4.5744E-01
 GRADIENT:   3.3733E+01  4.2537E+01 -5.3090E+00  4.3739E+01  6.7331E-01  2.1461E+00  3.3384E+00 -7.9957E+00 -1.0857E+01 -1.1232E+01
            -1.8874E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3341.57269729263        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      860
 NPARAMETR:  1.0145E+00  1.1251E+00  1.1597E+00  9.5615E-01  1.1346E+00  1.0162E+00  1.0898E+00  1.7912E+00  1.0104E+00  7.9779E-01
             1.4297E+00
 PARAMETER:  1.1436E-01  2.1787E-01  2.4816E-01  5.5160E-02  2.2630E-01  1.1606E-01  1.8600E-01  6.8291E-01  1.1033E-01 -1.2591E-01
             4.5745E-01
 GRADIENT:   7.8118E+00  3.2634E+01 -6.0493E+00  3.9480E+01 -4.4097E+00 -7.3599E-01  2.5277E+00 -8.5098E+00 -1.1345E+01 -1.1407E+01
            -1.8963E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3341.79217544731        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1030
 NPARAMETR:  1.0144E+00  1.1249E+00  1.1598E+00  9.5609E-01  1.1345E+00  1.0179E+00  1.0899E+00  1.7922E+00  1.0104E+00  7.9784E-01
             1.4312E+00
 PARAMETER:  1.1431E-01  2.1768E-01  2.4829E-01  5.5094E-02  2.2622E-01  1.1770E-01  1.8608E-01  6.8344E-01  1.1038E-01 -1.2585E-01
             4.5853E-01
 GRADIENT:   3.3517E+01  4.2213E+01 -5.3118E+00  4.3430E+01  4.3938E-01  2.4772E+00  3.3958E+00 -7.8219E+00 -1.0749E+01 -1.1152E+01
            -1.8607E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3341.89829176965        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1215
 NPARAMETR:  1.0144E+00  1.1248E+00  1.1599E+00  9.5605E-01  1.1345E+00  1.0410E+00  1.0899E+00  1.7928E+00  1.0105E+00  7.9786E-01
             1.4327E+00
 PARAMETER:  1.1428E-01  2.1756E-01  2.4835E-01  5.5058E-02  2.2619E-01  1.4017E-01  1.8610E-01  6.8376E-01  1.1041E-01 -1.2582E-01
             4.5957E-01
 GRADIENT:   7.3235E+00  3.2053E+01 -6.0419E+00  3.8938E+01 -4.7759E+00  8.5117E+00  2.6294E+00 -8.1925E+00 -1.1154E+01 -1.1255E+01
            -1.8438E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3342.02440097156        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1400
 NPARAMETR:  1.0107E+00  1.1247E+00  1.1600E+00  9.5603E-01  1.1383E+00  1.0184E+00  1.0704E+00  1.7931E+00  1.0105E+00  7.9789E-01
             1.4325E+00
 PARAMETER:  1.1068E-01  2.1751E-01  2.4841E-01  5.5033E-02  2.2954E-01  1.1820E-01  1.6805E-01  6.8393E-01  1.1043E-01 -1.2579E-01
             4.5945E-01
 GRADIENT:  -9.8174E-02  2.8919E+01 -7.1334E+00  4.0981E+01  4.8894E-02  1.2919E-01  1.4173E-01 -8.4019E+00 -1.1072E+01 -1.2544E+01
            -1.8587E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3342.36424672925        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1569
 NPARAMETR:  1.0107E+00  1.1247E+00  1.1599E+00  9.5605E-01  1.1373E+00  1.0180E+00  1.0691E+00  1.7934E+00  1.0658E+00  8.9271E-01
             1.4325E+00
 PARAMETER:  1.1062E-01  2.1755E-01  2.4835E-01  5.5055E-02  2.2864E-01  1.1780E-01  1.6680E-01  6.8412E-01  1.6375E-01 -1.3494E-02
             4.5939E-01
 GRADIENT:   2.4864E+01  3.0029E+01  6.1051E+00  4.6558E+01 -1.3868E+00  2.6119E+00  6.2279E+00 -4.9358E+00  3.4415E+00 -3.3570E-01
            -1.7876E+02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3342.42571626951        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1697
 NPARAMETR:  1.0107E+00  1.1246E+00  1.1601E+00  9.5600E-01  1.1429E+00  1.0161E+00  1.0683E+00  1.7928E+00  1.0656E+00  8.9329E-01
             1.4328E+00
 PARAMETER:  1.1067E-01  2.1745E-01  2.4848E-01  5.5005E-02  2.3359E-01  1.1592E-01  1.6610E-01  6.8378E-01  1.6352E-01 -1.2846E-02
             4.5962E-01
 GRADIENT:   2.4946E+01  2.6502E+01  4.2674E+00  4.8355E+01  6.0168E+00  1.8311E+00  6.3548E+00 -5.1591E+00  3.3145E+00 -5.4711E-01
            -1.7881E+02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3342.66332674477        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1884
 NPARAMETR:  1.0107E+00  1.1246E+00  1.1601E+00  9.5598E-01  1.1358E+00  1.0425E+00  1.0683E+00  1.7929E+00  1.0656E+00  8.9329E-01
             1.4358E+00
 PARAMETER:  1.1067E-01  2.1740E-01  2.4847E-01  5.4981E-02  2.2732E-01  1.4159E-01  1.6610E-01  6.8385E-01  1.6352E-01 -1.2845E-02
             4.6171E-01
 GRADIENT:   1.5957E-02  2.0349E+01  5.8081E+00  4.1582E+01 -8.7447E+00  9.1189E+00  5.5122E+00 -5.0914E+00  2.8119E+00 -2.1333E-01
            -1.7418E+02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3350.63573502580        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2067
 NPARAMETR:  1.0107E+00  1.1223E+00  1.1597E+00  9.5491E-01  1.1509E+00  1.0297E+00  1.0681E+00  1.7989E+00  1.0655E+00  8.9331E-01
             1.5799E+00
 PARAMETER:  1.1068E-01  2.1537E-01  2.4819E-01  5.3866E-02  2.4058E-01  1.2929E-01  1.6591E-01  6.8716E-01  1.6344E-01 -1.2822E-02
             5.5738E-01
 GRADIENT:  -2.4054E+00 -1.0688E-01 -1.5283E+01  4.1885E+01 -1.3740E+00  5.0817E+00  1.0145E+01  1.6029E+00  7.5578E+00  4.6675E+00
             3.3522E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3351.98044708464        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     2165
 NPARAMETR:  1.0046E+00  1.1486E+00  1.1598E+00  9.1666E-01  1.1511E+00  1.0041E+00  9.5205E-01  1.7984E+00  1.0436E+00  9.0719E-01
             1.5800E+00
 PARAMETER:  1.0462E-01  2.3854E-01  2.4829E-01  1.2979E-02  2.4069E-01  1.0411E-01  5.0867E-02  6.8689E-01  1.4265E-01  2.5966E-03
             5.5745E-01
 GRADIENT:   3.9515E+00  4.2290E+00 -5.7390E+00  3.4964E-01 -3.0912E+01 -2.7056E+00 -1.3135E+00  2.8102E+00  4.5117E-01 -8.9346E-01
             2.8240E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3352.08379772130        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     2323
 NPARAMETR:  1.0121E+00  1.1529E+00  1.1599E+00  9.1667E-01  1.1511E+00  1.0150E+00  9.5998E-01  1.7984E+00  1.0440E+00  9.1471E-01
             1.5800E+00
 PARAMETER:  1.1198E-01  2.4230E-01  2.4829E-01  1.2993E-02  2.4071E-01  1.1485E-01  5.9159E-02  6.8688E-01  1.4311E-01  1.0854E-02
             5.5741E-01
 GRADIENT:   4.7939E-01 -3.4218E-01 -5.7945E+00 -4.5513E-01 -3.8189E+01 -2.8506E-01  5.9352E-02  2.3242E+00  1.0205E-01  1.7048E-01
             2.7322E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3352.31509902255        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2504
 NPARAMETR:  1.0134E+00  1.1553E+00  1.1607E+00  9.1097E-01  1.1560E+00  1.0207E+00  9.4751E-01  1.7945E+00  1.0736E+00  9.2336E-01
             1.5674E+00
 PARAMETER:  1.1329E-01  2.4435E-01  2.4902E-01  6.7497E-03  2.4495E-01  1.2049E-01  4.6078E-02  6.8475E-01  1.7098E-01  2.0265E-02
             5.4939E-01
 GRADIENT:   3.5876E+00 -7.0116E+00 -5.0056E+00 -3.9341E+00 -3.5679E+01  1.8627E+00 -9.3617E-01  1.6108E+00  4.7972E+00  6.2246E-01
             1.1904E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -3352.35726416992        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2680
 NPARAMETR:  1.0132E+00  1.1560E+00  1.1608E+00  9.1076E-01  1.1566E+00  1.0204E+00  9.4772E-01  1.7940E+00  1.0723E+00  9.2384E-01
             1.5657E+00
 PARAMETER:  1.1314E-01  2.4497E-01  2.4911E-01  6.5256E-03  2.4551E-01  1.2015E-01  4.6303E-02  6.8446E-01  1.6979E-01  2.0781E-02
             5.4834E-01
 GRADIENT:   3.2794E+00 -6.6100E+00 -5.0760E+00 -3.7327E+00 -3.5395E+01  1.7242E+00 -8.7886E-01  1.4169E+00  4.4700E+00  5.8464E-01
             9.6326E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -3352.36771475048        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     2821
 NPARAMETR:  1.0132E+00  1.1562E+00  1.1610E+00  9.1072E-01  1.1568E+00  1.0146E+00  9.5359E-01  1.7934E+00  1.0724E+00  9.2379E-01
             1.5661E+00
 PARAMETER:  1.1308E-01  2.4510E-01  2.4924E-01  6.4752E-03  2.4565E-01  1.1449E-01  5.2478E-02  6.8411E-01  1.6987E-01  2.0731E-02
             5.4860E-01
 GRADIENT:   3.1598E+00 -6.1018E+00 -5.1131E+00 -3.6886E+00 -3.5223E+01 -4.7415E-01 -1.3341E-01  1.3792E+00  4.5181E+00  8.9823E-01
             1.0368E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -3354.05258393472        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     2956
 NPARAMETR:  9.9082E-01  1.1726E+00  1.2900E+00  9.0218E-01  1.2290E+00  1.0035E+00  9.5402E-01  1.9757E+00  1.0119E+00  9.2782E-01
             1.5513E+00
 PARAMETER:  9.0777E-02  2.5921E-01  3.5465E-01 -2.9435E-03  3.0617E-01  1.0351E-01  5.2925E-02  7.8094E-01  1.1188E-01  2.5078E-02
             5.3908E-01
 GRADIENT:  -2.6067E+01 -9.7259E+00 -4.9866E+00 -2.5230E+00  5.6892E+00 -3.8042E+00 -3.2507E-01 -6.5158E+00 -2.8011E+00 -3.0929E+00
            -3.6552E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -3354.42474712434        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:     3060
 NPARAMETR:  9.9105E-01  1.1757E+00  1.3137E+00  9.0127E-01  1.2365E+00  1.0040E+00  9.5348E-01  2.0371E+00  1.0079E+00  9.3381E-01
             1.5502E+00
 PARAMETER:  9.1010E-02  2.6186E-01  3.7284E-01 -3.9553E-03  3.1228E-01  1.0400E-01  5.2365E-02  8.1154E-01  1.0789E-01  3.1520E-02
             5.3836E-01
 GRADIENT:  -4.4869E+01 -1.9811E+01 -5.2834E+00 -6.5905E+00 -1.2371E+00 -5.6703E+00 -5.9793E-01 -6.2120E+00 -2.7768E+00 -3.0000E+00
            -3.9739E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -3354.46497109244        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     3223
 NPARAMETR:  9.9110E-01  1.1756E+00  1.3134E+00  9.0122E-01  1.2367E+00  1.0191E+00  9.5766E-01  2.0363E+00  1.0079E+00  9.3386E-01
             1.5506E+00
 PARAMETER:  9.1060E-02  2.6180E-01  3.7265E-01 -4.0054E-03  3.1244E-01  1.1891E-01  5.6737E-02  8.1113E-01  1.0790E-01  3.1570E-02
             5.3863E-01
 GRADIENT:  -2.4003E+01 -8.9253E+00 -4.5949E+00 -3.3584E+00  4.6414E+00  2.6096E+00  2.8224E-01 -5.6788E+00 -2.3540E+00 -2.6523E+00
            -2.2406E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -3355.16495774351        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     3358
 NPARAMETR:  1.0087E+00  1.1898E+00  1.3416E+00  9.0322E-01  1.2398E+00  1.0163E+00  9.5487E-01  2.0355E+00  1.0194E+00  9.5201E-01
             1.5530E+00
 PARAMETER:  1.0867E-01  2.7376E-01  3.9390E-01 -1.7931E-03  3.1498E-01  1.1620E-01  5.3821E-02  8.1073E-01  1.1918E-01  5.0822E-02
             5.4018E-01
 GRADIENT:   1.4701E+01  7.6874E+00  1.1222E+00  8.2551E+00 -5.9769E+00  2.2854E+00  1.6332E+00 -7.3578E+00  3.3181E-01 -3.2553E-01
             1.7010E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -3355.17727199930        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     3429
 NPARAMETR:  1.0036E+00  1.1911E+00  1.3373E+00  8.9787E-01  1.2451E+00  1.0119E+00  9.4236E-01  2.0355E+00  1.0205E+00  9.5942E-01
             1.5529E+00
 PARAMETER:  1.0361E-01  2.7486E-01  3.9063E-01 -7.7336E-03  3.1924E-01  1.1187E-01  4.0632E-02  8.1072E-01  1.2026E-01  5.8571E-02
             5.4015E-01
 GRADIENT:   2.8425E+00  2.4667E-01  2.9892E-01  2.6504E+00 -2.0082E+00  4.2837E-01  5.4293E-01 -7.2832E+00 -6.4892E-02 -2.3925E-01
            -3.1469E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -3355.32435179387        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     3501
 NPARAMETR:  1.0026E+00  1.2037E+00  1.3311E+00  8.8910E-01  1.2548E+00  1.0111E+00  9.2638E-01  2.0355E+00  1.0275E+00  9.7302E-01
             1.5536E+00
 PARAMETER:  1.0259E-01  2.8537E-01  3.8601E-01 -1.7545E-02  3.2699E-01  1.1104E-01  2.3525E-02  8.1072E-01  1.2714E-01  7.2653E-02
             5.4057E-01
 GRADIENT:   3.5000E-01 -3.7160E-01 -1.2607E-02  6.1470E-01 -3.2305E-01  4.7598E-02  1.0536E-01 -7.2013E+00 -3.6558E-02 -1.1766E-01
            -4.5624E-02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -3355.97642277877        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     3661
 NPARAMETR:  1.0122E+00  1.2931E+00  1.3073E+00  8.3959E-01  1.3192E+00  1.0164E+00  8.5653E-01  2.0355E+00  1.0729E+00  1.0432E+00
             1.5597E+00
 PARAMETER:  1.1215E-01  3.5706E-01  3.6798E-01 -7.4839E-02  3.7705E-01  1.1626E-01 -5.4863E-02  8.1072E-01  1.7038E-01  1.4225E-01
             5.4448E-01
 GRADIENT:   5.9708E-01  1.9229E-03  9.3904E-05  6.3241E-01  1.1827E-01  7.3414E-02 -2.6431E-01 -9.0412E+00 -5.1098E-02 -1.2767E-03
            -1.6437E-02

0ITERATION NO.:  130    OBJECTIVE VALUE:  -3356.31056956431        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     3821
 NPARAMETR:  1.0116E+00  1.2942E+00  1.3071E+00  8.3825E-01  1.3200E+00  1.0160E+00  8.5755E-01  2.1425E+00  1.0736E+00  1.0437E+00
             1.5596E+00
 PARAMETER:  1.1150E-01  3.5785E-01  3.6781E-01 -7.6440E-02  3.7766E-01  1.1588E-01 -5.3674E-02  8.6197E-01  1.7102E-01  1.4279E-01
             5.4443E-01
 GRADIENT:  -7.3233E-01 -1.4489E+00 -3.6720E+00 -1.0923E+00 -3.4011E+00 -8.1178E-02 -1.1387E-01 -3.9169E+00  1.0102E+00 -6.2157E-01
             3.6631E+00

0ITERATION NO.:  133    OBJECTIVE VALUE:  -3356.31525631489        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     3937
 NPARAMETR:  1.0116E+00  1.2942E+00  1.3076E+00  8.3826E-01  1.3201E+00  1.0160E+00  8.5767E-01  2.1440E+00  1.0735E+00  1.0438E+00
             1.5596E+00
 PARAMETER:  1.1151E-01  3.5787E-01  3.6818E-01 -7.6424E-02  3.7770E-01  1.1589E-01 -5.3639E-02  8.6270E-01  1.7090E-01  1.4288E-01
             5.4440E-01
 GRADIENT:  -1.7693E+05 -5.5134E+04  1.0718E+05  3.9457E+05 -5.2238E+04 -6.7672E-02 -1.1472E-01  4.5660E+04  1.1544E+05  2.7616E+05
            -7.2483E+04
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.6         1.6         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3937
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1104E-03 -3.2201E-02 -2.4769E-02  2.2506E-02 -3.5911E-02
 SE:             2.9809E-02  2.1140E-02  1.9984E-02  2.5122E-02  2.3562E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7029E-01  1.2771E-01  2.1518E-01  3.7031E-01  1.2749E-01

 ETASHRINKSD(%)  1.3474E-01  2.9177E+01  3.3051E+01  1.5840E+01  2.1065E+01
 ETASHRINKVR(%)  2.6930E-01  4.9841E+01  5.5179E+01  2.9170E+01  3.7692E+01
 EBVSHRINKSD(%)  5.6091E-01  2.9250E+01  3.6549E+01  1.7462E+01  2.0298E+01
 EBVSHRINKVR(%)  1.1187E+00  4.9945E+01  5.9740E+01  3.1875E+01  3.6477E+01
 RELATIVEINF(%)  9.8871E+01  1.4354E+01  3.0153E+01  2.3212E+01  2.4726E+01
 EPSSHRINKSD(%)  2.0002E+01
 EPSSHRINKVR(%)  3.6003E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3356.3152563148928     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1702.2258965464821     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   118.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3356.315       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.29E+00  1.31E+00  8.38E-01  1.32E+00  1.02E+00  8.58E-01  2.14E+00  1.07E+00  1.04E+00  1.56E+00
 


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
+        7.75E+08
 
 TH 2
+       -4.35E+02  4.60E+07
 
 TH 3
+        4.15E+02 -1.20E+04  4.26E+07
 
 TH 4
+        7.17E+04  1.78E+04 -1.68E+04  1.40E+09
 
 TH 5
+       -4.04E+02  6.18E+01 -3.03E+02  1.63E+04  3.97E+07
 
 TH 6
+       -7.51E+02 -1.88E+02  1.77E+02  1.02E+03 -1.71E+02  1.89E+02
 
 TH 7
+       -1.88E+04 -4.55E+03  4.40E+03  2.53E+04 -4.25E+03 -2.22E+00  7.38E+01
 
 TH 8
+        1.08E+02 -1.55E+03  2.99E+03 -4.37E+03 -7.35E+01  4.60E+01  1.15E+03  2.87E+06
 
 TH 9
+       -2.44E+03 -6.07E+02  5.75E+02 -4.40E+04 -5.45E+02  4.65E+02  1.16E+04  1.56E+02  2.93E+08
 
 TH10
+        1.75E+02  1.98E+01 -3.37E+01  7.89E+08  1.93E+01  5.73E+02  1.42E+04 -1.33E+01  3.61E+08  4.43E+08
 
 TH11
+       -2.44E+02  3.34E+03 -6.55E+03  9.53E+03  1.26E+02 -9.93E+01 -2.49E+03  2.04E+04 -3.14E+02  3.72E+01  1.37E+07
 
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
 #CPUT: Total CPU Time in Seconds,      133.328
Stop Time:
Sat Sep 18 03:44:44 CDT 2021
