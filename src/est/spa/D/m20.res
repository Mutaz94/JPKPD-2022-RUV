Wed Sep 29 19:49:35 CDT 2021
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
$DATA ../../../../data/spa/D/dat20.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -596.394138348061        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6734E+02 -1.2430E+02 -8.1959E+01 -1.9217E+02  3.4384E+02 -5.8241E+02 -3.4767E+02 -3.9332E+01 -6.5567E+02 -3.2478E+02
            -8.4603E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1318.70164537537        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0985E+00  8.4060E-01  1.4492E+00  1.2670E+00  9.5018E-01  1.9277E+00  1.6797E+00  1.1915E+00  2.7563E+00  1.9959E+00
             1.0605E+00
 PARAMETER:  1.9393E-01 -7.3641E-02  4.7105E-01  3.3663E-01  4.8897E-02  7.5633E-01  6.1861E-01  2.7522E-01  1.1139E+00  7.9110E-01
             1.5872E-01
 GRADIENT:  -1.1456E+01 -1.0354E+00 -6.2343E+01 -9.6708E+00 -8.1730E+00  9.4000E+00  4.4171E+00 -1.0149E+01  1.9213E+01  1.7710E+01
            -8.1193E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1373.01108021591        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      335
 NPARAMETR:  1.1426E+00  6.1938E-01  1.0379E+01  1.5868E+00  1.6109E+00  2.2362E+00  5.4162E-01  1.7821E+00  2.4521E+00  2.2851E+00
             1.2014E+00
 PARAMETER:  2.3334E-01 -3.7903E-01  2.4398E+00  5.6172E-01  5.7680E-01  9.0479E-01 -5.1318E-01  6.7780E-01  9.9694E-01  9.2643E-01
             2.8348E-01
 GRADIENT:   7.8371E+02  5.9526E+01  4.2392E+00  5.4539E+02  1.8006E+01  1.8006E+03  3.1143E+00 -2.5700E+00  4.6976E+02  6.0324E+01
             5.0601E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1388.19877836249        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      477             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1425E+00  6.1132E-01  1.0652E+01  1.5798E+00  1.5693E+00  2.2235E+00  2.5819E-01  1.8492E+00  2.4549E+00  1.6458E+00
             1.1905E+00
 PARAMETER:  2.3318E-01 -3.9214E-01  2.4657E+00  5.5733E-01  5.5061E-01  8.9908E-01 -1.2540E+00  7.1476E-01  9.9809E-01  5.9820E-01
             2.7438E-01
 GRADIENT:   7.9845E+02  5.9138E+01 -6.6729E-01  5.5175E+02  5.4408E+01  1.8582E+03  1.4973E+00 -5.4256E+00  4.9799E+02  2.4333E+00
             9.8124E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1388.80889786410        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      655
 NPARAMETR:  1.1427E+00  6.6895E-01  1.0677E+01  1.5814E+00  1.5676E+00  2.2075E+00  1.2829E-01  1.9368E+00  2.4836E+00  1.7002E+00
             1.1543E+00
 PARAMETER:  2.3342E-01 -3.0204E-01  2.4681E+00  5.5831E-01  5.4952E-01  8.9187E-01 -1.9534E+00  7.6106E-01  1.0097E+00  6.3076E-01
             2.4351E-01
 GRADIENT:   1.2506E+01  7.8163E+00 -1.6468E+00  1.6835E+01  2.7035E+01  1.2382E+02  1.0951E-01 -6.5357E+00 -4.0610E+01  7.1793E-01
            -4.2353E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1390.21610552282        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      821
 NPARAMETR:  1.1420E+00  6.6837E-01  1.0755E+01  1.5787E+00  1.5648E+00  2.1846E+00  1.0695E-01  1.9327E+00  2.6584E+00  1.6949E+00
             1.1625E+00
 PARAMETER:  2.3275E-01 -3.0291E-01  2.4754E+00  5.5661E-01  5.4777E-01  8.8145E-01 -2.1354E+00  7.5891E-01  1.0777E+00  6.2761E-01
             2.5059E-01
 GRADIENT:   1.2116E+01 -3.8739E+00 -1.7118E+00  1.3235E+01  2.6414E+01  1.1480E+02  9.4883E-02 -6.4367E+00 -1.1931E+01  8.5677E-01
            -6.3874E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1391.72791305485        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      977
 NPARAMETR:  1.1424E+00  6.9901E-01  1.0710E+01  1.4817E+00  1.5618E+00  2.2060E+00  2.4847E-02  1.9593E+00  2.6774E+00  1.6878E+00
             1.1641E+00
 PARAMETER:  2.3315E-01 -2.5809E-01  2.4712E+00  4.9322E-01  5.4586E-01  8.9119E-01 -3.5950E+00  7.7261E-01  1.0848E+00  6.2341E-01
             2.5198E-01
 GRADIENT:   8.3591E+02  4.6515E+01  4.9159E-01  4.6322E+02  4.6526E+01  1.9256E+03  4.9954E-02 -6.2287E+00  5.7632E+02  8.7732E+00
             2.4777E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1392.27951345447        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     1108
 NPARAMETR:  1.1401E+00  7.9122E-01  1.0936E+01  1.4705E+00  1.5543E+00  2.1442E+00  1.0000E-02  2.0121E+00  2.7453E+00  1.6800E+00
             1.1654E+00
 PARAMETER:  2.3115E-01 -1.3418E-01  2.4921E+00  4.8562E-01  5.4100E-01  8.6275E-01 -5.1655E+00  7.9918E-01  1.1099E+00  6.1877E-01
             2.5309E-01
 GRADIENT:   1.1104E+01  8.3462E+00 -9.7356E-01  1.8064E+01  2.0064E+01  1.0228E+02  0.0000E+00 -6.4982E+00 -3.1316E+01  1.0231E+00
            -9.1757E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1394.43875973278        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:     1259             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1400E+00  7.9009E-01  1.1079E+01  1.4665E+00  1.5529E+00  2.2414E+00  1.0000E-02  2.5623E+00  2.7485E+00  1.6746E+00
             1.1669E+00
 PARAMETER:  2.3103E-01 -1.3561E-01  2.5051E+00  4.8291E-01  5.4015E-01  9.0710E-01 -5.1655E+00  1.0409E+00  1.1111E+00  6.1558E-01
             2.5436E-01
 GRADIENT:   8.2286E+02  4.3865E+01  4.1759E+00  4.2435E+02  4.1550E+01  1.9727E+03  0.0000E+00 -1.0321E+01  5.5531E+02  1.1121E+01
             1.4065E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1396.70598780683        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1426             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1351E+00  7.6817E-01  1.1168E+01  1.3666E+00  1.5153E+00  2.2365E+00  1.0000E-02  4.0234E+00  2.8718E+00  1.6361E+00
             1.1355E+00
 PARAMETER:  2.2674E-01 -1.6374E-01  2.5130E+00  4.1235E-01  5.1565E-01  9.0490E-01 -5.1655E+00  1.4921E+00  1.1549E+00  5.9231E-01
             2.2705E-01
 GRADIENT:   8.5318E+02  2.8078E+01  1.4803E+00  3.7049E+02  2.7353E+01  2.1022E+03  0.0000E+00  2.7252E+00  6.7891E+02  1.0773E+01
             4.9055E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1398.31493802281        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  1.1631E+00  8.7344E-01  1.1191E+01  1.3485E+00  1.5106E+00  2.1360E+00  1.0000E-02  4.0230E+00  3.2039E+00  1.6194E+00
             1.1365E+00
 PARAMETER:  2.5106E-01 -3.5312E-02  2.5151E+00  3.9902E-01  5.1249E-01  8.5894E-01 -5.1655E+00  1.4920E+00  1.2644E+00  5.8207E-01
             2.2795E-01
 GRADIENT:   2.2957E+01 -1.3468E+01 -3.3530E+00  9.0370E+00  9.1334E-01  1.0210E+02  0.0000E+00  2.0151E+00 -1.6288E+00  6.2469E-01
            -5.0660E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1398.97119716701        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:     1730
 NPARAMETR:  1.1090E+00  8.7394E-01  1.1353E+01  1.2200E+00  1.5064E+00  2.2306E+00  1.0000E-02  3.9890E+00  3.2267E+00  1.6130E+00
             1.1376E+00
 PARAMETER:  2.0350E-01 -3.4743E-02  2.5295E+00  2.9887E-01  5.0970E-01  9.0229E-01 -5.1655E+00  1.4835E+00  1.2714E+00  5.7808E-01
             2.2888E-01
 GRADIENT:  -2.8721E+00 -2.9308E+01 -2.6149E+00 -1.0709E+01  5.7966E-01  1.4853E+02  0.0000E+00  1.4314E+00 -4.6444E+00  5.4005E-01
             8.3635E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1400.53338094244        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1898
 NPARAMETR:  1.1509E+00  9.4018E-01  1.5791E+01  1.1796E+00  1.5040E+00  2.2624E+00  1.0000E-02  4.0071E+00  3.2062E+00  1.6052E+00
             1.1435E+00
 PARAMETER:  2.4058E-01  3.8316E-02  2.8594E+00  2.6514E-01  5.0814E-01  9.1640E-01 -5.1655E+00  1.4881E+00  1.2651E+00  5.7324E-01
             2.3406E-01
 GRADIENT:   9.0201E+02  3.4329E+01  4.9060E-01  1.9306E+02  7.3856E+00  2.0943E+03  0.0000E+00  6.8032E+00  7.4678E+02  5.7125E+00
            -8.3621E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1400.80157079421        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     2069             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1499E+00  9.4659E-01  1.6755E+01  1.2072E+00  1.5161E+00  2.2490E+00  1.0000E-02  3.9857E+00  3.2200E+00  1.6102E+00
             1.1446E+00
 PARAMETER:  2.3972E-01  4.5116E-02  2.9187E+00  2.8831E-01  5.1613E-01  9.1049E-01 -5.1655E+00  1.4827E+00  1.2694E+00  5.7635E-01
             2.3510E-01
 GRADIENT:   8.9608E+02  4.0903E+01  2.9549E+00  2.1391E+02  1.4827E+01  2.0760E+03  0.0000E+00  5.0587E-02  7.2824E+02  7.2641E+00
             3.9225E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1400.90000235287        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     2264             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1326E+00  9.4532E-01  1.6894E+01  1.2063E+00  1.5333E+00  2.2545E+00  1.0000E-02  4.0752E+00  3.2674E+00  1.6126E+00
             1.1434E+00
 PARAMETER:  2.2454E-01  4.3771E-02  2.9270E+00  2.8756E-01  5.2741E-01  9.1291E-01 -5.1655E+00  1.5049E+00  1.2840E+00  5.7788E-01
             2.3398E-01
 GRADIENT:   8.3425E+02  3.5289E+01  1.3806E+00  2.1373E+02  2.0232E+01  2.1211E+03  0.0000E+00  2.7217E+00  7.4671E+02  5.7175E+00
             2.6519E+00

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1400.90185434556        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     2365
 NPARAMETR:  1.1344E+00  9.4573E-01  1.7063E+01  1.2069E+00  1.5368E+00  2.2474E+00  1.0000E-02  4.0715E+00  3.2495E+00  1.6223E+00
             1.1420E+00
 PARAMETER:  2.2708E-01  4.3771E-02  2.9269E+00  2.8927E-01  5.2745E-01  9.1171E-01 -5.1655E+00  1.5049E+00  1.2840E+00  5.7804E-01
             2.3395E-01
 GRADIENT:   1.4667E+04 -5.9221E+04 -1.0239E+00  5.7560E+03 -6.3108E+03  9.0358E-01  0.0000E+00  4.2196E-01  4.6062E+03 -7.6258E-01
             4.0421E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2365
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.9283E-05 -7.2387E-04 -9.7579E-03 -3.6531E-03 -4.0007E-02
 SE:             2.9863E-02  1.9927E-04  7.8594E-03  2.9400E-02  2.3244E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9842E-01  2.8070E-04  2.1440E-01  9.0112E-01  8.5218E-02

 ETASHRINKSD(%)  1.0000E-10  9.9332E+01  7.3670E+01  1.5046E+00  2.2130E+01
 ETASHRINKVR(%)  1.0000E-10  9.9996E+01  9.3067E+01  2.9866E+00  3.9362E+01
 EBVSHRINKSD(%)  1.1218E-01  9.9616E+01  7.9158E+01  6.1248E-01  1.8246E+01
 EBVSHRINKVR(%)  2.2423E-01  9.9999E+01  9.5656E+01  1.2212E+00  3.3163E+01
 RELATIVEINF(%)  9.9751E+01  8.1571E-04  3.6492E+00  5.4550E+01  5.6095E+01
 EPSSHRINKSD(%)  4.1584E+01
 EPSSHRINKVR(%)  6.5876E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1400.9018543455607     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -665.75102778182247     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    38.44
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1400.902       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  9.45E-01  1.69E+01  1.21E+00  1.53E+00  2.25E+00  1.00E-02  4.08E+00  3.27E+00  1.61E+00  1.14E+00
 


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
+        1.70E+02
 
 TH 2
+        3.41E+06  9.31E+06
 
 TH 3
+        6.62E+01  1.31E+01  1.07E-01
 
 TH 4
+       -1.85E+06 -5.04E+06 -9.63E+03  7.04E+01
 
 TH 5
+        2.92E+06 -5.38E+00 -1.95E+01  5.89E+05  1.27E+05
 
 TH 6
+        1.57E+05 -8.15E+02  6.30E-03  1.43E+01 -5.96E+00  3.91E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -2.29E+03  1.40E+05 -2.68E+02  3.80E+04  7.09E+02 -1.59E-01  0.00E+00  2.13E+03
 
 TH 9
+       -7.72E+04 -2.11E+05  2.69E+00  5.66E+04  6.83E-01  1.83E+01  0.00E+00 -3.15E+03  9.48E+03
 
 TH10
+       -1.05E+03  2.31E+00 -1.74E+01 -7.78E+02  3.25E+02  1.61E-01  0.00E+00 -1.43E+04  4.26E+04  4.33E+01
 
 TH11
+       -1.21E+06  7.98E+04  2.56E+00 -3.08E+02  1.76E+02  6.80E-01  0.00E+00  4.96E+04 -1.79E+03  3.36E+01  2.30E+02
 
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
 #CPUT: Total CPU Time in Seconds,       45.988
Stop Time:
Wed Sep 29 19:50:22 CDT 2021
