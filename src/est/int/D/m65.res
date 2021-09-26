Sat Sep 25 06:07:47 CDT 2021
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
$DATA ../../../../data/int/D/dat65.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   27620.1453641366        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.1625E+01  1.8090E+02 -1.7515E+01  1.5494E+02  3.3923E+02 -2.8794E+03 -2.5782E+03 -5.1165E+01 -2.1298E+03 -1.0589E+03
            -5.4551E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1073.56569792973        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  2.6947E+00  1.4916E+00  9.6055E-01  2.0867E+00  5.9710E-01  5.5912E+00  5.1448E+00  1.0405E+00  5.1926E+00  3.8685E+00
             9.5674E+00
 PARAMETER:  1.0913E+00  4.9987E-01  5.9749E-02  8.3558E-01 -4.1567E-01  1.8212E+00  1.7380E+00  1.3975E-01  1.7472E+00  1.4529E+00
             2.3584E+00
 GRADIENT:   4.8766E+01  2.0774E+01 -2.3695E+01  3.4777E+01 -3.1606E+01  1.0169E+02  8.4656E+01  5.0117E+00  1.0516E+02  6.7085E+01
             2.6139E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1175.55019666961        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.6141E+00  8.2235E-01  2.5300E+01  3.6647E+00  1.1025E+00  5.5043E+00  6.9628E+00  7.7735E-01  4.6721E+00  2.7989E+00
             9.6377E+00
 PARAMETER:  5.7878E-01 -9.5587E-02  3.3308E+00  1.3988E+00  1.9755E-01  1.8055E+00  2.0406E+00 -1.5187E-01  1.6416E+00  1.1292E+00
             2.3657E+00
 GRADIENT:   9.7628E+00  1.1714E+01  3.6319E+00  6.1910E+01 -1.0357E+02  1.2213E+02  8.5785E+00  6.1885E-01  4.7677E+01  4.8014E+01
             3.1605E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1340.68335455667        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.5286E+00  4.3457E-01  3.1031E+01  1.2427E+00  1.9229E+00  3.0042E+00  1.0442E+01  1.4640E-01  1.5006E+00  8.1011E-01
             9.3973E+00
 PARAMETER:  5.2434E-01 -7.3339E-01  3.5350E+00  3.1730E-01  7.5383E-01  1.2000E+00  2.4458E+00 -1.8214E+00  5.0584E-01 -1.1059E-01
             2.3404E+00
 GRADIENT:   3.7892E+01 -6.1733E+00  2.1786E+00 -8.0655E+01 -4.2078E+01  1.5540E+01  8.4290E+01  1.2476E-04  2.2166E+01  8.7870E+00
             2.9994E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1381.32514568007        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.1511E+00  7.1500E-01  3.3493E+01  1.2570E+00  2.1422E+00  2.9553E+00  6.5989E+00  3.0637E+00  1.4821E+00  8.0976E-01
             7.9765E+00
 PARAMETER:  2.4072E-01 -2.3547E-01  3.6113E+00  3.2869E-01  8.6184E-01  1.1836E+00  1.9869E+00  1.2196E+00  4.9348E-01 -1.1102E-01
             2.1765E+00
 GRADIENT:  -1.9462E+01 -1.4700E+01  1.3300E-02  7.8005E-01  1.1826E+00  9.3340E+00 -1.1287E+01  9.8833E-02  3.8012E+00  8.2707E+00
             2.3207E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1382.05664202194        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:      411
 NPARAMETR:  1.1587E+00  7.4750E-01  3.1214E+01  1.2347E+00  2.1317E+00  2.9433E+00  6.5479E+00  3.3386E+00  1.4537E+00  7.9724E-01
             7.9412E+00
 PARAMETER:  2.4728E-01 -1.9102E-01  3.5409E+00  3.1087E-01  8.5693E-01  1.1795E+00  1.9791E+00  1.3056E+00  4.7411E-01 -1.2660E-01
             2.1721E+00
 GRADIENT:  -1.9770E+01 -1.4485E+01  2.5129E-02 -4.7001E-01 -9.2510E-01  8.2976E-01 -3.8753E+01  1.0323E-01  3.3961E+00  7.8252E+00
             9.4343E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1384.70969006869        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      514
 NPARAMETR:  1.1585E+00  1.0529E+00  3.1158E+01  1.2349E+00  2.1327E+00  2.9415E+00  6.5550E+00  2.7428E+00  1.4533E+00  7.9729E-01
             7.9496E+00
 PARAMETER:  2.4715E-01  1.5159E-01  3.5391E+00  3.1102E-01  8.5737E-01  1.1789E+00  1.9802E+00  1.1090E+00  4.7387E-01 -1.2654E-01
             2.1731E+00
 GRADIENT:  -1.8130E+01 -4.0440E-01  1.9535E-01  1.7166E+01 -6.2378E+00  7.6379E+00  2.3981E+00  6.8840E-02  1.5661E+00  7.5090E+00
             5.5149E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1384.71162153274        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  1.1585E+00  1.0543E+00  3.1158E+01  1.2349E+00  2.1327E+00  2.9415E+00  6.5550E+00  2.6973E+00  1.4533E+00  7.9729E-01
             7.9496E+00
 PARAMETER:  2.4715E-01  1.5283E-01  3.5391E+00  3.1102E-01  8.5737E-01  1.1789E+00  1.9802E+00  1.0922E+00  4.7387E-01 -1.2654E-01
             2.1731E+00
 GRADIENT:  -2.0196E+01 -6.3410E-01  1.9252E-01  1.6050E+01 -6.7258E+00  4.7927E-01 -2.0950E+01  6.6928E-02  1.3108E+00  7.4995E+00
             1.5599E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1384.80309988890        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.1580E+00  1.0539E+00  3.0950E+01  1.2356E+00  2.1363E+00  2.9343E+00  6.5785E+00  2.9775E-01  1.4520E+00  7.9748E-01
             7.9756E+00
 PARAMETER:  2.4670E-01  1.5254E-01  3.5324E+00  3.1159E-01  8.5909E-01  1.1765E+00  1.9838E+00 -1.1115E+00  4.7297E-01 -1.2630E-01
             2.1764E+00
 GRADIENT:  -1.8566E+01 -3.0090E-01  2.3234E-01  1.6377E+01 -5.4055E+00  6.7772E+00  3.3981E+00  1.1420E-03  1.9235E+00  7.4915E+00
             1.1798E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1385.67372520223        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      914
 NPARAMETR:  1.1575E+00  1.0533E+00  3.0421E+01  1.2364E+00  2.1523E+00  2.9228E+00  7.1555E+00  4.1596E-01  1.4487E+00  7.9785E-01
             7.9194E+00
 PARAMETER:  2.4626E-01  1.5192E-01  3.5151E+00  3.1221E-01  8.6656E-01  1.1726E+00  2.0679E+00 -7.7716E-01  4.7067E-01 -1.2583E-01
             2.1693E+00
 GRADIENT:  -2.0334E+01  2.2632E+00 -7.7657E-02  6.6575E+00  1.6847E-01 -2.2107E+00 -4.8925E-01  2.0929E-03  4.2811E+00  7.4260E+00
            -3.3670E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1385.78833060084        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1102
 NPARAMETR:  1.1617E+00  1.0537E+00  3.0789E+01  1.2330E+00  2.1528E+00  2.9492E+00  7.1657E+00  1.0000E-02  1.4448E+00  7.9733E-01
             7.9298E+00
 PARAMETER:  2.4987E-01  1.5226E-01  3.5272E+00  3.0942E-01  8.6678E-01  1.1815E+00  2.0693E+00 -6.4840E+01  4.6794E-01 -1.2648E-01
             2.1706E+00
 GRADIENT:  -1.9197E+01  2.1760E+00 -5.9379E-02  5.7422E+00 -7.7069E-02  1.1621E+00  9.1111E-02  0.0000E+00  4.5272E+00  7.4219E+00
            -3.9231E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1385.83672840620        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1291
 NPARAMETR:  1.1642E+00  1.0539E+00  3.1001E+01  1.2309E+00  2.1535E+00  2.9649E+00  7.1870E+00  1.0000E-02  1.4423E+00  7.9703E-01
             7.9339E+00
 PARAMETER:  2.5204E-01  1.5245E-01  3.5340E+00  3.0775E-01  8.6711E-01  1.1869E+00  2.0723E+00 -1.0370E+02  4.6623E-01 -1.2686E-01
             2.1711E+00
 GRADIENT:  -1.8509E+01  2.1874E+00 -5.8118E-02  4.9843E+00 -6.2686E-02  3.1189E+00  8.8746E-01  0.0000E+00  4.7292E+00  7.4173E+00
             9.0893E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1385.90859254855        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     1423
 NPARAMETR:  1.1639E+00  9.8649E-01  3.2545E+01  1.2313E+00  2.1554E+00  2.9613E+00  7.2006E+00  1.0000E-02  1.4416E+00  7.9632E-01
             7.9505E+00
 PARAMETER:  2.5180E-01  8.6398E-02  3.5826E+00  3.0805E-01  8.6798E-01  1.1856E+00  2.0742E+00 -1.0380E+02  4.6576E-01 -1.2775E-01
             2.1732E+00
 GRADIENT:  -1.6583E+01  1.4539E-01 -2.0925E-02  3.0106E+00  6.5687E-01  9.9365E+00  2.5026E+01  0.0000E+00  5.4738E+00  7.4847E+00
             1.0624E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1385.91329207155        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1493
 NPARAMETR:  1.1639E+00  9.6086E-01  3.3781E+01  1.2313E+00  2.1553E+00  2.9610E+00  7.1959E+00  1.0000E-02  1.4416E+00  7.9468E-01
             7.9477E+00
 PARAMETER:  2.5181E-01  6.0074E-02  3.6199E+00  3.0803E-01  8.6794E-01  1.1855E+00  2.0735E+00 -1.0380E+02  4.6576E-01 -1.2982E-01
             2.1729E+00
 GRADIENT:  -1.6528E+01 -7.4796E-01  2.4625E-02  2.0312E+00  1.6643E-01  9.9138E+00  2.4335E+01  0.0000E+00  5.5543E+00  7.4428E+00
             1.0523E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1385.92438240816        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1563
 NPARAMETR:  1.1641E+00  9.1991E-01  3.6223E+01  1.2311E+00  2.1548E+00  2.9594E+00  7.1675E+00  1.0000E-02  1.4417E+00  7.8503E-01
             7.9305E+00
 PARAMETER:  2.5192E-01  1.6525E-02  3.6897E+00  3.0794E-01  8.6768E-01  1.1850E+00  2.0696E+00 -1.0380E+02  4.6582E-01 -1.4203E-01
             2.1707E+00
 GRADIENT:  -1.6313E+01 -2.2881E+00  1.0749E-01  9.3613E-01 -9.7006E-01  9.7289E+00  2.2395E+01  0.0000E+00  5.4805E+00  7.1611E+00
             6.8791E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1385.94141031767        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1635
 NPARAMETR:  1.1643E+00  8.9519E-01  3.8055E+01  1.2309E+00  2.1538E+00  2.9568E+00  7.1189E+00  1.0000E-02  1.4419E+00  7.6871E-01
             7.9010E+00
 PARAMETER:  2.5210E-01 -1.0714E-02  3.7390E+00  3.0776E-01  8.6723E-01  1.1841E+00  2.0627E+00 -1.0380E+02  4.6593E-01 -1.6304E-01
             2.1670E+00
 GRADIENT:  -1.5992E+01 -3.3649E+00  1.6952E-01  1.0302E+00 -2.0860E+00  9.3866E+00  2.0026E+01  0.0000E+00  5.1294E+00  6.6788E+00
            -6.5268E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1385.97796011774        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1831             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1637E+00  8.9533E-01  3.7777E+01  1.2317E+00  2.1574E+00  2.9498E+00  7.1470E+00  1.0000E-02  1.4405E+00  7.6837E-01
             7.9341E+00
 PARAMETER:  2.5162E-01 -1.0563E-02  3.7317E+00  3.0836E-01  8.6892E-01  1.1817E+00  2.0667E+00 -1.0380E+02  4.6502E-01 -1.6348E-01
             2.1712E+00
 GRADIENT:  -1.6501E+01 -3.3032E+00  1.2714E-01  2.4750E-01 -1.0730E+00  8.5933E+00  2.0946E+01  0.0000E+00  5.4476E+00  6.7748E+00
             7.7508E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1386.56565423876        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     1977
 NPARAMETR:  1.1656E+00  9.2951E-01  3.7601E+01  1.2315E+00  2.1577E+00  2.9491E+00  7.1439E+00  1.0000E-02  1.4381E+00  7.1438E-01
             7.9341E+00
 PARAMETER:  2.5323E-01  2.6903E-02  3.7270E+00  3.0820E-01  8.6905E-01  1.1815E+00  2.0663E+00 -1.0380E+02  4.6335E-01 -2.3634E-01
             2.1712E+00
 GRADIENT:  -1.8204E+01 -2.1767E+00  1.4804E-01  9.3108E-01 -2.8781E+00  1.4517E+00 -4.8603E+00  0.0000E+00  4.8029E+00  5.4765E+00
             1.2879E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1387.76968996402        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     2091
 NPARAMETR:  1.1654E+00  9.6495E-01  3.4241E+01  1.2317E+00  2.1588E+00  2.9471E+00  7.1521E+00  1.0000E-02  1.3360E+00  5.5916E-01
             7.9439E+00
 PARAMETER:  2.5309E-01  6.4324E-02  3.6334E+00  3.0837E-01  8.6955E-01  1.1808E+00  2.0674E+00 -1.0380E+02  3.8966E-01 -4.8131E-01
             2.1724E+00
 GRADIENT:  -1.6186E+01 -1.8181E-01  1.6127E-01  1.1767E+01 -5.3668E+00  8.7180E+00  2.0278E+01  0.0000E+00  4.1879E-01  2.5849E+00
            -2.2253E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1387.78604778831        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     2264
 NPARAMETR:  1.1654E+00  9.7682E-01  3.0768E+01  1.2318E+00  2.1593E+00  2.9463E+00  7.1558E+00  1.0000E-02  1.3358E+00  5.5923E-01
             7.9484E+00
 PARAMETER:  2.5302E-01  7.6550E-02  3.5265E+00  3.0845E-01  8.6977E-01  1.1805E+00  2.0679E+00 -1.0380E+02  3.8956E-01 -4.8119E-01
             2.1730E+00
 GRADIENT:  -1.8403E+01  2.3962E-02  3.4148E-03  1.0885E+01 -3.2626E+00  1.4701E+00 -4.9408E+00  0.0000E+00  2.7002E-01  2.6206E+00
            -5.3261E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1387.81265777151        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2444
 NPARAMETR:  1.1655E+00  9.7632E-01  3.1142E+01  1.2316E+00  2.1601E+00  2.9455E+00  7.1851E+00  1.0000E-02  1.3359E+00  5.5915E-01
             7.9866E+00
 PARAMETER:  2.5319E-01  7.6035E-02  3.5385E+00  3.0829E-01  8.7014E-01  1.1803E+00  2.0720E+00 -1.0380E+02  3.8962E-01 -4.8134E-01
             2.1778E+00
 GRADIENT:  -1.8684E+01  6.2674E-03  9.4346E-03  9.5713E+00 -3.3175E+00  1.5392E+00 -3.8198E+00  0.0000E+00  7.7704E-01  2.6875E+00
             4.6668E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1387.89515127500        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2625
 NPARAMETR:  1.1667E+00  9.5099E-01  3.3582E+01  1.2306E+00  2.1672E+00  2.9374E+00  7.3798E+00  1.0000E-02  1.3359E+00  5.5880E-01
             7.9619E+00
 PARAMETER:  2.5416E-01  4.9744E-02  3.6140E+00  3.0747E-01  8.7346E-01  1.1775E+00  2.0988E+00 -1.0380E+02  3.8962E-01 -4.8196E-01
             2.1747E+00
 GRADIENT:  -1.8287E+01 -9.6235E-04 -1.9888E-03  5.3431E+00 -2.2867E+00  4.2090E-01  1.2426E+00  0.0000E+00  1.7299E+00  2.6670E+00
             1.9376E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1388.26621675219        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     2758
 NPARAMETR:  1.1845E+00  9.5107E-01  3.3605E+01  1.2267E+00  2.1692E+00  2.9362E+00  7.3684E+00  1.0000E-02  1.3273E+00  5.3930E-01
             7.9624E+00
 PARAMETER:  2.6935E-01  4.9837E-02  3.6147E+00  3.0437E-01  8.7435E-01  1.1771E+00  2.0972E+00 -1.0380E+02  3.8314E-01 -5.1748E-01
             2.1747E+00
 GRADIENT:  -1.4648E+01 -1.5814E-01 -5.3470E-03  5.0922E+00 -2.3496E+00  7.8647E-01  1.1035E+00  0.0000E+00  1.6544E+00  2.3944E+00
            -1.4183E+00

0ITERATION NO.:  114    OBJECTIVE VALUE:  -1388.27294106700        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     2906
 NPARAMETR:  1.1845E+00  9.5118E-01  3.3612E+01  1.2268E+00  2.1692E+00  2.9361E+00  7.3686E+00  1.0000E-02  1.3266E+00  5.3798E-01
             7.9627E+00
 PARAMETER:  2.6935E-01  4.9844E-02  3.6147E+00  3.0437E-01  8.7435E-01  1.1771E+00  2.0972E+00 -1.0380E+02  3.8260E-01 -5.1993E-01
             2.1747E+00
 GRADIENT:   5.4841E+03 -1.6403E-01 -3.6142E-03 -4.8631E+03 -8.4920E+02  6.2756E+02 -7.0547E+02  0.0000E+00  1.9380E+03 -2.8473E+03
            -6.8012E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2906
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6417E-02  3.9194E-02 -7.1989E-07 -8.1856E-02 -1.0530E-03
 SE:             2.8929E-02  2.3960E-02  4.7222E-06  1.4203E-02  8.1633E-03
 N:                     100         100         100         100         100

 P VAL.:         3.6115E-01  1.0188E-01  8.7883E-01  8.2741E-09  8.9736E-01

 ETASHRINKSD(%)  3.0829E+00  1.9729E+01  9.9984E+01  5.2418E+01  7.2652E+01
 ETASHRINKVR(%)  6.0707E+00  3.5566E+01  1.0000E+02  7.7360E+01  9.2521E+01
 EBVSHRINKSD(%)  2.6782E+00  1.3063E+01  9.9981E+01  5.5108E+01  7.1030E+01
 EBVSHRINKVR(%)  5.2847E+00  2.4419E+01  1.0000E+02  7.9847E+01  9.1607E+01
 RELATIVEINF(%)  9.4625E+01  4.2782E+01  9.6166E-07  1.1460E+01  2.1473E+00
 EPSSHRINKSD(%)  7.6346E+00
 EPSSHRINKVR(%)  1.4686E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1388.2729410670015     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       265.81641870140925     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   102.90
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1388.273       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.18E+00  9.51E-01  3.36E+01  1.23E+00  2.17E+00  2.94E+00  7.37E+00  1.00E-02  1.33E+00  5.38E-01  7.96E+00
 


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
+        3.64E+06
 
 TH 2
+        1.66E+01 -3.20E+01
 
 TH 3
+       -2.90E-02  7.34E-03  6.09E-04
 
 TH 4
+       -4.42E+00 -1.18E+00  3.75E-02  2.66E+06
 
 TH 5
+        5.34E+02 -6.18E+00 -1.58E-01 -5.98E+00  1.03E+05
 
 TH 6
+       -2.90E+02  1.68E+00 -3.46E-03  8.43E-01  5.31E+01  3.08E+04
 
 TH 7
+        6.55E+01  1.90E+00 -1.67E-03 -6.60E+00 -1.11E+01  2.15E+00  1.55E+03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.85E-01  1.31E+01 -1.76E-02 -3.44E+01  3.75E+00 -1.44E-01  1.66E+00  0.00E+00  1.44E+06
 
 TH10
+        5.10E+02 -2.35E+01  3.08E-02 -4.33E+02 -8.12E+01  4.78E+01 -1.06E+01  0.00E+00 -2.61E+06  4.73E+06
 
 TH11
+        5.53E+01 -1.84E+00  8.68E-04 -7.84E+00 -1.11E+01  4.99E+01 -1.02E-01  0.00E+00  3.32E+00 -8.10E+00  1.24E+03
 
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
 #CPUT: Total CPU Time in Seconds,      121.811
Stop Time:
Sat Sep 25 06:09:51 CDT 2021
