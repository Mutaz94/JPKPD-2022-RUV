Sat Sep 25 05:42:01 CDT 2021
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
$DATA ../../../../data/int/D/dat36.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   28382.4958264153        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1993E+02  3.4973E+02 -3.3021E+01  1.7061E+02  1.5955E+02 -1.5696E+03 -7.6929E+02 -5.7838E+01 -1.3777E+03 -5.7439E+02
            -5.9897E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -939.730942252789        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3756E+00  1.8811E+00  9.0799E-01  2.1600E+00  9.0627E-01  4.0443E+00  4.0601E+00  9.8182E-01  2.4796E+00  1.6242E+00
             1.3359E+01
 PARAMETER:  4.1890E-01  7.3187E-01  3.4827E-03  8.7010E-01  1.5786E-03  1.4973E+00  1.5012E+00  8.1648E-02  1.0081E+00  5.8500E-01
             2.6922E+00
 GRADIENT:  -8.4901E+00  2.5965E+01 -4.0618E+01  1.1387E+02 -1.5735E+00  1.3991E+02  4.3417E+01  4.4801E+00  3.3425E+01  3.8924E+01
             4.8317E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1013.53270799301        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2097E+00  2.0954E+00  2.6123E+01  2.5989E+00  2.3082E+00  2.9807E+00  1.0827E+01  5.9191E-01  1.6896E+00  2.2023E+00
             1.3222E+01
 PARAMETER:  2.9035E-01  8.3975E-01  3.3628E+00  1.0551E+00  9.3645E-01  1.1922E+00  2.4820E+00 -4.2439E-01  6.2448E-01  8.8951E-01
             2.6818E+00
 GRADIENT:  -3.3050E+01  3.2313E+01 -2.9553E+00  5.7352E+01  2.8377E-01  1.0428E+02  2.8995E+01  4.5596E-02  9.5107E+00  5.6454E+01
             5.0530E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1245.33995812148        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1618E+00  4.5133E-01  9.8403E+01  1.7692E+00  3.1985E+00  2.1860E+00  7.9514E+00  1.4164E+00  1.0459E+00  2.0679E-01
             8.5276E+00
 PARAMETER:  2.5000E-01 -6.9556E-01  4.6891E+00  6.7054E-01  1.2627E+00  8.8209E-01  2.1733E+00  4.4815E-01  1.4486E-01 -1.4760E+00
             2.2433E+00
 GRADIENT:  -5.7681E+00  3.5688E+00 -2.6322E+00 -1.1525E+01  1.0903E+02  4.6435E+01  4.9820E+01  4.8699E-03 -2.2723E+01  4.1577E-01
             1.5790E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1276.97907980288        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.1798E+00  4.5347E-01  3.7253E+01  1.7209E+00  2.2072E+00  1.9622E+00  5.9342E+00  5.1813E-01  1.4454E+00  2.1341E-01
             8.2789E+00
 PARAMETER:  2.6535E-01 -6.9082E-01  3.7177E+00  6.4282E-01  8.9174E-01  7.7406E-01  1.8807E+00 -5.5753E-01  4.6842E-01 -1.4445E+00
             2.2137E+00
 GRADIENT:  -2.0171E+00 -2.7592E+00  4.5528E-01 -3.1642E+00 -1.6558E-01 -2.2004E-01  4.7924E+00 -4.1718E-03  2.6643E+00  1.0204E-01
             8.2266E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1277.29457503231        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.1938E+00  8.1307E-01  1.4811E+01  1.5198E+00  2.0459E+00  2.0034E+00  4.8423E+00  1.1882E+00  1.3950E+00  2.5522E-01
             8.1997E+00
 PARAMETER:  2.7716E-01 -1.0694E-01  2.7954E+00  5.1859E-01  8.1586E-01  7.9486E-01  1.6774E+00  2.7245E-01  4.3287E-01 -1.2656E+00
             2.2041E+00
 GRADIENT:   3.3893E+00  2.4059E+00  1.1361E+00  1.1568E+01 -6.0758E+00  2.7801E+00  1.0214E+00 -1.6324E-01  1.9356E+00  1.4698E-01
            -1.5086E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1277.64443198455        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.1887E+00  1.0774E+00  8.9343E+00  1.3266E+00  1.9463E+00  1.9782E+00  4.2320E+00  2.3386E+00  1.2598E+00  2.6299E-01
             8.2087E+00
 PARAMETER:  2.7290E-01  1.7457E-01  2.2899E+00  3.8259E-01  7.6590E-01  7.8219E-01  1.5427E+00  9.4953E-01  3.3093E-01 -1.2357E+00
             2.2052E+00
 GRADIENT:  -2.9447E-02  2.1899E+00 -1.6962E-01  1.2377E+01  5.2468E-01 -1.2055E+00 -4.4443E+00 -5.9325E-02  3.0202E-01  6.5832E-01
            -7.9858E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1277.81564766309        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      511
 NPARAMETR:  1.1872E+00  1.2971E+00  6.1037E+00  1.1412E+00  1.8685E+00  1.9850E+00  3.8965E+00  2.0428E+00  1.1063E+00  2.4156E-01
             8.2158E+00
 PARAMETER:  2.7158E-01  3.6011E-01  1.9089E+00  2.3211E-01  7.2514E-01  7.8561E-01  1.4601E+00  8.1435E-01  2.0099E-01 -1.3206E+00
             2.2061E+00
 GRADIENT:  -8.2373E-01  4.3714E-01  1.0564E+00  8.9093E-01 -2.4476E+00 -1.7216E-02 -1.4394E+00  1.4032E-01  4.1066E-02  7.0458E-01
            -1.4808E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1277.82346513177        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      581
 NPARAMETR:  1.1898E+00  1.4141E+00  4.7659E+00  1.0600E+00  1.8221E+00  1.9862E+00  3.7267E+00  1.7497E+00  1.0366E+00  2.2020E-01
             8.2165E+00
 PARAMETER:  2.7376E-01  4.4651E-01  1.6615E+00  1.5826E-01  6.9999E-01  7.8624E-01  1.4155E+00  6.5946E-01  1.3595E-01 -1.4132E+00
             2.2061E+00
 GRADIENT:   9.8014E-02  3.2647E-01  8.2475E-01 -8.1919E-01 -1.2213E+00  1.0443E-01 -8.5748E-01  1.1431E-01  1.6756E-01  6.3545E-01
            -8.0468E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1277.82970214140        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  1.1906E+00  1.4742E+00  4.0344E+00  1.0166E+00  1.7882E+00  1.9856E+00  3.6409E+00  1.5042E+00  9.9148E-01  2.0019E-01
             8.2194E+00
 PARAMETER:  2.7448E-01  4.8808E-01  1.4949E+00  1.1642E-01  6.8120E-01  7.8591E-01  1.3922E+00  5.0828E-01  9.1445E-02 -1.5085E+00
             2.2065E+00
 GRADIENT:   3.0251E-01 -6.6981E-02  2.0206E-01 -1.3500E+00  5.2972E-01 -4.3081E-03 -3.2180E-01  6.9481E-02  2.1926E-01  5.4507E-01
             1.7776E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1277.85983260679        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.1895E+00  1.4947E+00  3.4460E+00  9.9359E-01  1.7429E+00  1.9840E+00  3.6111E+00  1.1653E+00  9.4428E-01  1.6105E-01
             8.2246E+00
 PARAMETER:  2.7355E-01  5.0196E-01  1.3372E+00  9.3567E-02  6.5554E-01  7.8514E-01  1.3840E+00  2.5300E-01  4.2663E-02 -1.7260E+00
             2.2071E+00
 GRADIENT:  -2.9306E-01 -8.2809E-01 -7.2715E-01 -9.1890E-01  1.9191E+00 -1.3348E-01  6.2417E-01 -3.9915E-02  3.1293E-03  3.5367E-01
             8.4923E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1278.06630666024        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      797
 NPARAMETR:  1.1901E+00  1.4374E+00  3.8734E+00  1.0383E+00  1.7651E+00  1.9804E+00  3.6899E+00  1.4451E+00  9.8555E-01  9.1567E-02
             8.2247E+00
 PARAMETER:  2.7403E-01  4.6280E-01  1.4541E+00  1.3759E-01  6.6821E-01  7.8328E-01  1.4056E+00  4.6819E-01  8.5446E-02 -2.2907E+00
             2.2071E+00
 GRADIENT:  -1.1192E-01 -2.5606E-01 -1.1740E+00  2.4529E+00  3.3719E+00 -7.0945E-01 -3.5129E-01  1.1766E-01 -1.6106E-01  1.1529E-01
             1.9389E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1278.12573063523        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  1.1905E+00  1.4523E+00  3.8192E+00  1.0251E+00  1.7537E+00  1.9867E+00  3.6773E+00  1.4531E+00  9.8494E-01  3.2111E-02
             8.2199E+00
 PARAMETER:  2.7439E-01  4.7318E-01  1.4400E+00  1.2476E-01  6.6171E-01  7.8646E-01  1.4022E+00  4.7372E-01  8.4822E-02 -3.3386E+00
             2.2066E+00
 GRADIENT:   3.2077E-01 -6.4616E-02  1.4499E-01 -7.1794E-01 -1.0243E+00  3.4615E-01  5.7542E-01  2.0309E-02  5.7227E-02  1.4148E-02
            -2.1112E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1278.12747286701        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  1.1904E+00  1.4565E+00  3.8564E+00  1.0250E+00  1.7616E+00  1.9852E+00  3.6692E+00  1.4496E+00  9.8795E-01  1.9815E-02
             8.2217E+00
 PARAMETER:  2.7431E-01  4.7601E-01  1.4497E+00  1.2469E-01  6.6625E-01  7.8574E-01  1.4000E+00  4.7130E-01  8.7877E-02 -3.8213E+00
             2.2068E+00
 GRADIENT:   2.3280E-01  3.7248E-02  3.9960E-02 -2.3653E-01 -1.5359E-01  9.0979E-02  1.4682E-01 -9.0643E-03  9.2553E-02  5.3559E-03
            -5.8966E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1278.13313182711        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1037
 NPARAMETR:  1.1913E+00  1.4572E+00  3.8226E+00  1.0238E+00  1.7566E+00  1.9866E+00  3.6733E+00  1.4483E+00  9.8404E-01  1.2791E-02
             8.2202E+00
 PARAMETER:  2.7501E-01  4.7651E-01  1.4409E+00  1.2348E-01  6.6340E-01  7.8645E-01  1.4011E+00  4.7040E-01  8.3912E-02 -4.2590E+00
             2.2066E+00
 GRADIENT:  -1.7593E+00 -8.1311E-01  4.3473E-02 -7.7524E-01 -1.0848E+00 -2.5138E+00 -3.4657E+00  2.5467E-03  2.8920E-02  2.2047E-03
            -4.4080E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1278.37425669328        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1214
 NPARAMETR:  1.1950E+00  1.0364E+00  8.9542E+00  1.3215E+00  1.9398E+00  1.9987E+00  4.4654E+00  2.4351E+00  1.2155E+00  1.3364E-02
             8.2455E+00
 PARAMETER:  2.7812E-01  1.3579E-01  2.2921E+00  3.7880E-01  7.6259E-01  7.9249E-01  1.5964E+00  9.9001E-01  2.9516E-01 -4.2152E+00
             2.2097E+00
 GRADIENT:   2.1682E-01  8.1539E-01 -1.3830E-01  1.3402E+00  5.7145E-01 -1.7358E-01 -7.4579E-02  5.5315E-02  1.6971E-01  1.7564E-03
            -1.1703E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1278.39038751284        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1389
 NPARAMETR:  1.1940E+00  9.2969E-01  1.0499E+01  1.3862E+00  1.9650E+00  1.9986E+00  4.6791E+00  2.6127E+00  1.2468E+00  1.0000E-02
             8.2490E+00
 PARAMETER:  2.7728E-01  2.7092E-02  2.4513E+00  4.2656E-01  7.7547E-01  7.9244E-01  1.6431E+00  1.0604E+00  3.2060E-01 -4.9563E+00
             2.2101E+00
 GRADIENT:   9.0871E-03  8.8961E-03  6.9827E-02 -8.4452E-02 -1.3729E-01  4.0747E-02  2.7032E-02 -3.1319E-02 -9.4771E-02  0.0000E+00
            -2.6024E-02

0ITERATION NO.:   84    OBJECTIVE VALUE:  -1278.39054669239        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  1.1940E+00  9.3009E-01  1.0448E+01  1.3862E+00  1.9638E+00  1.9985E+00  4.6781E+00  2.6098E+00  1.2484E+00  1.0000E-02
             8.2489E+00
 PARAMETER:  2.7728E-01  2.7423E-02  2.4465E+00  4.2655E-01  7.7487E-01  7.9236E-01  1.6427E+00  1.0589E+00  3.2180E-01 -4.8729E+00
             2.2101E+00
 GRADIENT:  -2.4060E-03 -2.2134E-03  6.5334E-04 -5.7365E-03  1.4245E-03 -8.0419E-03 -1.6177E-02 -2.9727E-03 -2.9166E-03  0.0000E+00
             5.1600E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1522
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1505E-02  3.3044E-02 -1.8350E-02 -6.3225E-02 -1.1724E-05
 SE:             2.9228E-02  2.3923E-02  7.9864E-03  1.5051E-02  1.4870E-04
 N:                     100         100         100         100         100

 P VAL.:         6.9385E-01  1.6719E-01  2.1580E-02  2.6637E-05  9.3715E-01

 ETASHRINKSD(%)  2.0815E+00  1.9856E+01  7.3244E+01  4.9576E+01  9.9502E+01
 ETASHRINKVR(%)  4.1196E+00  3.5769E+01  9.2841E+01  7.4574E+01  9.9998E+01
 EBVSHRINKSD(%)  3.3036E+00  1.6977E+01  7.6634E+01  4.9180E+01  9.9460E+01
 EBVSHRINKVR(%)  6.4980E+00  3.1071E+01  9.4540E+01  7.4173E+01  9.9997E+01
 RELATIVEINF(%)  9.3299E+01  2.5543E+01  1.3092E+00  9.8676E+00  7.2166E-04
 EPSSHRINKSD(%)  7.5455E+00
 EPSSHRINKVR(%)  1.4522E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1278.3905466923882     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       375.69881307602259     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.50
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    17.22
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1278.391       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  9.30E-01  1.04E+01  1.39E+00  1.96E+00  2.00E+00  4.68E+00  2.61E+00  1.25E+00  1.00E-02  8.25E+00
 


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
+        1.11E+02
 
 TH 2
+       -3.14E+00  1.64E+01
 
 TH 3
+        4.88E-02  3.68E-01  1.09E-01
 
 TH 4
+       -2.30E+01  3.38E+01 -1.11E+00  1.06E+02
 
 TH 5
+        4.03E+00 -1.53E+01 -2.83E+00  1.40E+01  7.59E+01
 
 TH 6
+        1.34E+01 -4.67E+00  1.12E-01 -1.52E+01 -5.43E-01  3.17E+00
 
 TH 7
+        4.19E+00 -2.60E+00  8.73E-02 -8.51E+00 -1.04E+00  1.44E+00  7.33E-01
 
 TH 8
+       -7.83E-02 -2.60E-01 -9.12E-02  1.06E+00  2.35E+00 -1.12E-01 -8.34E-02  7.67E-02
 
 TH 9
+        7.59E+00 -8.53E+00  1.11E-01 -2.39E+01  7.00E-01  3.70E+00  1.96E+00 -1.23E-01  5.62E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.83E+00 -2.28E+00  7.42E-02 -6.73E+00 -1.11E+00  7.29E-01  4.66E-01 -6.61E-02  1.47E+00  0.00E+00  9.28E-01
 
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
+        1.74E+02
 
 TH 2
+       -1.80E+00  4.65E+01
 
 TH 3
+        1.44E-01  4.58E-01  1.44E-01
 
 TH 4
+       -4.33E+00  4.21E+01 -7.39E-01  1.15E+02
 
 TH 5
+       -1.49E+00 -1.38E+01 -2.59E+00  2.67E-01  7.24E+01
 
 TH 6
+       -3.27E+00 -5.27E-01  2.08E-02  2.34E-01 -1.55E-01  3.81E+01
 
 TH 7
+        6.67E-01  4.99E+00 -7.23E-02 -8.44E+00  1.66E+00  9.19E-02  4.83E+00
 
 TH 8
+       -9.13E-02 -6.80E-01 -2.44E-01  8.83E-01  2.03E+00  9.45E-02  9.88E-02  1.11E+00
 
 TH 9
+        1.09E+00 -5.15E+00 -2.85E-01 -2.13E+01  4.60E+00 -1.34E+00  2.34E+00  3.26E-01  2.34E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.69E+00 -4.10E+00 -4.05E-02 -7.49E+00 -6.32E-02  2.07E+00  1.78E-01  1.63E-01  2.27E+00  0.00E+00  1.52E+01
 
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
+        1.83E+02
 
 TH 2
+        6.39E+01  4.15E+01
 
 TH 3
+        2.70E-01  4.84E-01  8.96E-02
 
 TH 4
+        7.38E+01  4.08E+01  4.05E-01  1.16E+02
 
 TH 5
+       -8.56E+00 -1.75E+01 -2.31E+00 -1.71E+01  7.13E+01
 
 TH 6
+        4.00E+01  1.65E+01 -9.65E-02 -6.03E+00  6.83E+00  5.20E+01
 
 TH 7
+        6.57E+00  5.52E+00 -1.29E-01 -9.76E+00  2.72E+00  8.16E+00  5.27E+00
 
 TH 8
+       -7.50E-01 -4.80E-01 -8.68E-02 -9.40E-01  1.30E+00 -1.41E-01  1.56E-01  2.29E-01
 
 TH 9
+       -7.11E+00 -4.97E+00 -4.58E-02 -2.90E+01  6.01E+00  4.37E+00  3.17E+00 -2.06E-01  1.82E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.25E+01 -1.02E+01  1.04E-01 -3.01E+01  1.05E+01  7.21E-01 -1.15E-01 -6.71E-01  1.52E+01  0.00E+00  5.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       52.866
Stop Time:
Sat Sep 25 05:42:55 CDT 2021
