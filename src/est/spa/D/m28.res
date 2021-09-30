Wed Sep 29 19:54:58 CDT 2021
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
$DATA ../../../../data/spa/D/dat28.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   9790.11488652302        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5611E+02 -1.1342E+02 -2.8365E+01 -1.5017E+02  1.7152E+02 -8.9983E+02 -6.0319E+02 -4.7299E+01 -7.7326E+02 -1.9944E+02
            -2.0141E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -675.317827821779        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2907E+00  1.3044E+00  9.9271E-01  1.6526E+00  1.0694E+00  1.8046E+00  1.3112E+00  9.8471E-01  1.2610E+00  1.0331E+00
             1.4251E+01
 PARAMETER:  3.5520E-01  3.6573E-01  9.2687E-02  6.0236E-01  1.6707E-01  6.9033E-01  3.7098E-01  8.4589E-02  3.3191E-01  1.3255E-01
             2.7568E+00
 GRADIENT:  -3.6761E+01  2.3480E+01  2.0328E+00  3.1195E+01 -1.5309E+01  2.8189E+01  2.4079E+00  3.6512E+00  1.5216E+01  4.4440E+00
             2.2533E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -691.394782472402        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3069E+00  8.4315E-01  1.4631E+00  2.2921E+00  1.5838E+00  1.7747E+00  4.0439E+00  5.7501E-01  1.3975E+00  2.8538E+00
             1.2579E+01
 PARAMETER:  3.6769E-01 -7.0611E-02  4.8058E-01  9.2946E-01  5.5982E-01  6.7362E-01  1.4972E+00 -4.5337E-01  4.3466E-01  1.1487E+00
             2.6320E+00
 GRADIENT:  -1.7189E+01  2.7769E+01 -2.5539E+00  8.2572E+01 -1.6121E+01 -3.5336E+00  1.1992E+01  1.0690E+00  1.7514E+01  1.0259E+01
             1.6917E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -742.695971638541        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1174E+00  1.0855E+00  1.4174E+00  1.2447E+00  6.7308E+00  1.5176E+00  1.5915E+00  1.0000E-02  8.0167E-01  5.8656E+00
             9.1878E+00
 PARAMETER:  2.1105E-01  1.8207E-01  4.4882E-01  3.1893E-01  2.0067E+00  5.1716E-01  5.6469E-01 -5.3180E+00 -1.2106E-01  1.8691E+00
             2.3179E+00
 GRADIENT:   1.3994E+01  4.9793E+00  1.1842E+01 -1.7883E+01 -2.6143E-01  9.7921E-01  1.8420E+00  0.0000E+00 -4.0199E+00 -2.5876E+00
             3.1995E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -764.645548732737        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.9954E-01  3.9055E-01  6.9380E-01  1.5913E+00  7.4598E+00  1.3855E+00  2.1192E+00  1.0000E-02  4.8602E-01  7.2390E+00
             8.7003E+00
 PARAMETER:  9.9541E-02 -8.4019E-01 -2.6558E-01  5.6456E-01  2.1095E+00  4.2603E-01  8.5104E-01 -1.0387E+01 -6.2151E-01  2.0795E+00
             2.2634E+00
 GRADIENT:  -5.5098E+01  2.4187E+01  6.9633E+00  7.7950E+01  6.9899E+00 -3.4883E+01  1.2524E+00  0.0000E+00 -2.9577E-01 -7.8222E+00
            -4.7996E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -834.207067985257        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  6.8968E-01  1.0000E-02  8.8285E-02  8.0667E-01  5.6528E+00  1.3089E+00  2.2354E-01  1.0000E-02  1.0000E-02  1.3112E+00
             8.1931E+00
 PARAMETER: -2.7153E-01 -5.5694E+00 -2.3272E+00 -1.1483E-01  1.8321E+00  3.6920E-01 -1.3981E+00 -3.3849E+01 -5.3032E+00  3.7095E-01
             2.2033E+00
 GRADIENT:   2.8674E+01  0.0000E+00 -1.0658E+02  2.7193E+02  1.4039E+01 -4.8977E+01  1.1495E-04  0.0000E+00  0.0000E+00  3.9711E+00
            -9.4042E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -865.201598905201        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  5.2427E-01  1.0000E-02  3.6670E-02  4.1725E-01  6.9888E+00  1.4273E+00  2.9731E-02  1.0000E-02  1.0000E-02  7.1189E-01
             8.7699E+00
 PARAMETER: -5.4576E-01 -7.8274E+00 -3.2058E+00 -7.7407E-01  2.0443E+00  4.5576E-01 -3.4156E+00 -4.7064E+01 -8.3190E+00 -2.3984E-01
             2.2713E+00
 GRADIENT:   1.8836E+01  0.0000E+00 -2.0600E+01  2.1506E+01 -1.7957E+00  7.5289E+00  2.2829E-06  0.0000E+00  0.0000E+00  2.2443E-01
             5.9318E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -866.333367718551        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  4.8648E-01  1.0000E-02  3.1799E-02  3.7239E-01  8.1391E+00  1.3718E+00  2.1897E-02  1.0000E-02  1.0000E-02  3.8354E-01
             8.6412E+00
 PARAMETER: -6.2056E-01 -8.2582E+00 -3.3483E+00 -8.8781E-01  2.1967E+00  4.1609E-01 -3.7214E+00 -4.9341E+01 -8.8265E+00 -8.5830E-01
             2.2565E+00
 GRADIENT:   2.0859E+00  0.0000E+00 -7.9282E+00  9.3631E+00 -1.4368E-01 -1.4300E-01  7.1101E-07  0.0000E+00  0.0000E+00  9.7815E-03
            -1.5401E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -866.367126673581        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      866
 NPARAMETR:  4.8487E-01  1.0000E-02  3.1843E-02  3.7082E-01  9.4214E+00  1.3704E+00  1.9898E-02  1.0000E-02  1.0000E-02  1.4331E-01
             8.6433E+00
 PARAMETER: -6.2387E-01 -8.2582E+00 -3.3469E+00 -8.9204E-01  2.3430E+00  4.1513E-01 -3.8172E+00 -4.9341E+01 -8.8265E+00 -1.8428E+00
             2.2568E+00
 GRADIENT:  -2.2025E-01  0.0000E+00  2.6316E-01 -6.6954E-01 -2.8488E-02 -5.1856E-02  1.9421E-07  0.0000E+00  0.0000E+00  7.9385E-05
            -1.6343E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -866.375036122819        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1037
 NPARAMETR:  4.8395E-01  1.0000E-02  3.1703E-02  3.6995E-01  1.0271E+01  1.3698E+00  2.0237E-02  1.0000E-02  1.0000E-02  4.4699E-02
             8.6407E+00
 PARAMETER: -6.2576E-01 -8.2582E+00 -3.3514E+00 -8.9439E-01  2.4294E+00  4.1466E-01 -3.8002E+00 -4.9341E+01 -8.8265E+00 -3.0078E+00
             2.2565E+00
 GRADIENT:  -7.4312E-02  0.0000E+00 -2.5798E+00  3.0694E+00  1.3053E-02 -1.4452E-01  1.3876E-07  0.0000E+00  0.0000E+00  6.1864E-07
            -7.0350E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -866.377994531362        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1233             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8465E-01  1.0000E-02  3.1647E-02  3.6956E-01  1.0003E+01  1.3706E+00  2.1681E-02  1.0000E-02  1.0000E-02  2.6509E-02
             8.6457E+00
 PARAMETER: -6.2432E-01 -8.2582E+00 -3.3531E+00 -8.9545E-01  2.4029E+00  4.1528E-01 -3.7313E+00 -4.9341E+01 -8.8265E+00 -3.5303E+00
             2.2571E+00
 GRADIENT:   5.2513E+01  0.0000E+00  5.5644E+01  2.4203E+01  2.8029E-02  5.3081E+00  5.9381E-07  0.0000E+00  0.0000E+00  1.1635E-06
             1.7531E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -866.380710726270        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1429             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8405E-01  1.0000E-02  3.1634E-02  3.6907E-01  9.2836E+00  1.3704E+00  2.7476E-02  1.0000E-02  1.0000E-02  2.3989E-02
             8.6475E+00
 PARAMETER: -6.2557E-01 -8.2582E+00 -3.3535E+00 -8.9676E-01  2.3283E+00  4.1513E-01 -3.4944E+00 -4.9341E+01 -8.8265E+00 -3.6302E+00
             2.2573E+00
 GRADIENT:   5.1075E+01  0.0000E+00  5.8922E+01  2.0257E+01 -1.9054E-02  5.3450E+00  1.0727E-06  0.0000E+00  0.0000E+00  4.1147E-06
             1.8294E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -866.383599169678        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1629             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8392E-01  1.0000E-02  3.1598E-02  3.6890E-01  9.5405E+00  1.3703E+00  2.5123E-02  1.0000E-02  1.0000E-02  2.2469E-02
             8.6465E+00
 PARAMETER: -6.2583E-01 -8.2582E+00 -3.3547E+00 -8.9724E-01  2.3555E+00  4.1500E-01 -3.5840E+00 -4.9341E+01 -8.8265E+00 -3.6956E+00
             2.2571E+00
 GRADIENT:   5.1421E+01  0.0000E+00  5.7834E+01  2.1748E+01  1.9630E-03  5.3061E+00  8.5478E-07  0.0000E+00  0.0000E+00  2.1440E-06
             1.8031E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -866.389360843561        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1824             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8350E-01  1.0000E-02  3.1464E-02  3.6809E-01  9.8226E+00  1.3698E+00  2.2328E-02  1.0000E-02  1.0000E-02  1.4004E-02
             8.6440E+00
 PARAMETER: -6.2671E-01 -8.2582E+00 -3.3589E+00 -8.9942E-01  2.3847E+00  4.1464E-01 -3.7019E+00 -4.9341E+01 -8.8265E+00 -4.1684E+00
             2.2569E+00
 GRADIENT:   5.2360E+01  0.0000E+00  5.5548E+01  2.4846E+01  2.6494E-02  5.2094E+00  6.5727E-07  0.0000E+00  0.0000E+00  5.0521E-07
             1.7423E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -866.391687396555        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2009
 NPARAMETR:  4.8309E-01  1.0000E-02  3.1476E-02  3.6777E-01  9.8104E+00  1.3698E+00  2.1347E-02  1.0000E-02  1.0000E-02  1.4345E-02
             8.6461E+00
 PARAMETER: -6.2755E-01 -8.2582E+00 -3.3585E+00 -9.0029E-01  2.3834E+00  4.1466E-01 -3.7468E+00 -4.9341E+01 -8.8265E+00 -4.1444E+00
             2.2571E+00
 GRADIENT:   5.1323E-01  0.0000E+00 -9.3536E-01  3.1837E-01 -3.4486E-03  3.9503E-02  1.9983E-07  0.0000E+00  0.0000E+00  3.1869E-07
             2.6995E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -866.394303210700        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2198             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8299E-01  1.0000E-02  3.1401E-02  3.6745E-01  1.0017E+01  1.3696E+00  1.9489E-02  1.0000E-02  1.0000E-02  1.3619E-02
             8.6446E+00
 PARAMETER: -6.2775E-01 -8.2582E+00 -3.3609E+00 -9.0118E-01  2.4043E+00  4.1453E-01 -3.8379E+00 -4.9341E+01 -8.8265E+00 -4.1963E+00
             2.2569E+00
 GRADIENT:   5.2257E+01  0.0000E+00  5.6101E+01  2.4281E+01  2.3787E-02  5.2468E+00  4.9979E-07  0.0000E+00  0.0000E+00  3.3978E-07
             1.7583E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -866.396670312839        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2392             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8273E-01  1.0000E-02  3.1391E-02  3.6711E-01  9.5192E+00  1.3696E+00  1.9741E-02  1.0000E-02  1.0000E-02  1.2346E-02
             8.6459E+00
 PARAMETER: -6.2829E-01 -8.2582E+00 -3.3612E+00 -9.0210E-01  2.3533E+00  4.1453E-01 -3.8251E+00 -4.9341E+01 -8.8265E+00 -4.2944E+00
             2.2571E+00
 GRADIENT:   5.1613E+01  0.0000E+00  5.8203E+01  2.1622E+01 -1.3015E-03  5.3013E+00  5.6697E-07  0.0000E+00  0.0000E+00  7.2560E-07
             1.8049E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -866.398929180903        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2579
 NPARAMETR:  4.8253E-01  1.0000E-02  3.1356E-02  3.6683E-01  9.7518E+00  1.3695E+00  1.8080E-02  1.0000E-02  1.0000E-02  1.0277E-02
             8.6455E+00
 PARAMETER: -6.2871E-01 -8.2582E+00 -3.3624E+00 -9.0285E-01  2.3775E+00  4.1442E-01 -3.9129E+00 -4.9341E+01 -8.8265E+00 -4.4779E+00
             2.2570E+00
 GRADIENT:   7.2235E-01  0.0000E+00 -1.5096E+00  9.0634E-01 -1.4502E-03  2.7773E-02  1.5131E-07  0.0000E+00  0.0000E+00  1.9753E-07
            -1.1134E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -866.415471204580        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:     2685
 NPARAMETR:  4.7963E-01  1.0000E-02  3.0948E-02  3.6389E-01  1.0164E+01  1.3677E+00  2.3373E-02  1.0000E-02  1.0000E-02  1.0000E-02
             8.6338E+00
 PARAMETER: -6.3473E-01 -8.2582E+00 -3.3754E+00 -9.1091E-01  2.4189E+00  4.1311E-01 -3.6562E+00 -4.9341E+01 -8.8265E+00 -4.7513E+00
             2.2557E+00
 GRADIENT:   5.1364E+01  0.0000E+00  5.4955E+01  2.7580E+01  2.6095E-02  5.0733E+00  7.1042E-07  0.0000E+00  0.0000E+00  0.0000E+00
             1.6644E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -866.416464695227        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     2760
 NPARAMETR:  4.7710E-01  1.0000E-02  3.0726E-02  3.6211E-01  1.0480E+01  1.3660E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             8.6222E+00
 PARAMETER: -6.4003E-01 -8.2582E+00 -3.3826E+00 -9.1581E-01  2.4495E+00  4.1188E-01 -1.2260E+01 -4.9341E+01 -8.8265E+00 -4.7513E+00
             2.2543E+00
 GRADIENT:   4.9384E+01  0.0000E+00  5.4764E+01  2.9600E+01  1.8570E-02  4.8414E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
             1.5813E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -866.425219784384        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2956
 NPARAMETR:  4.7848E-01  1.0000E-02  3.0833E-02  3.6228E-01  9.2237E+00  1.3673E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             8.6338E+00
 PARAMETER: -6.3714E-01 -8.2582E+00 -3.3792E+00 -9.1535E-01  2.3218E+00  4.1280E-01 -9.2221E+00 -4.9341E+01 -8.8265E+00 -4.7513E+00
             2.2557E+00
 GRADIENT:  -1.7937E+00  0.0000E+00 -8.1388E-02 -1.8375E-01 -5.5110E-02 -7.8662E-02  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -6.5769E-01

0ITERATION NO.:  103    OBJECTIVE VALUE:  -866.427838783285        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     3048
 NPARAMETR:  4.7847E-01  1.0000E-02  3.0839E-02  3.6202E-01  1.0906E+01  1.3670E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             8.6390E+00
 PARAMETER: -6.3715E-01 -8.2582E+00 -3.3790E+00 -9.1605E-01  2.4893E+00  4.1265E-01 -9.2221E+00 -4.9341E+01 -8.8265E+00 -4.7513E+00
             2.2563E+00
 GRADIENT:  -9.7598E-01  0.0000E+00 -6.9871E-03 -6.6672E-01 -1.6716E-03 -4.8113E-02  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -2.3898E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3048
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0044E-03  1.4922E-06  9.6898E-05 -2.2978E-04  4.9416E-07
 SE:             2.8900E-02  1.1652E-06  3.2674E-04  4.1999E-04  1.4685E-06
 N:                     100         100         100         100         100

 P VAL.:         9.1720E-01  2.0031E-01  7.6680E-01  5.8431E-01  7.3649E-01

 ETASHRINKSD(%)  3.1810E+00  9.9996E+01  9.8905E+01  9.8593E+01  9.9995E+01
 ETASHRINKVR(%)  6.2607E+00  1.0000E+02  9.9988E+01  9.9980E+01  1.0000E+02
 EBVSHRINKSD(%)  3.5166E+00  9.9996E+01  9.8915E+01  9.8603E+01  9.9995E+01
 EBVSHRINKVR(%)  6.9095E+00  1.0000E+02  9.9988E+01  9.9980E+01  1.0000E+02
 RELATIVEINF(%)  3.1681E+00  7.5092E-09  7.3545E-05  9.7835E-05  1.8088E-08
 EPSSHRINKSD(%)  8.3471E+00
 EPSSHRINKVR(%)  1.5997E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -866.42783878328487     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -131.27701221954669     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.34
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -866.428       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.78E-01  1.00E-02  3.08E-02  3.62E-01  1.09E+01  1.37E+00  1.00E-02  1.00E-02  1.00E-02  1.00E-02  8.64E+00
 


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
+        7.83E+01
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.56E+03  0.00E+00  5.50E+05
 
 TH 4
+        7.08E+02  0.00E+00 -5.94E+04  6.41E+03
 
 TH 5
+        7.90E-02  0.00E+00 -6.62E+00  7.15E-01  7.97E-05
 
 TH 6
+       -2.71E+00  0.00E+00  2.27E+02 -2.45E+01 -2.73E-03  9.34E-02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.81E+00  0.00E+00  2.36E+02 -2.54E+01 -2.84E-03  9.72E-02  0.00E+00  0.00E+00  0.00E+00  0.00E+00  1.01E-01
 
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
+        2.37E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.04E+04  0.00E+00  8.63E+05
 
 TH 4
+       -3.47E+02  0.00E+00 -9.30E+04  1.11E+04
 
 TH 5
+        4.47E-01  0.00E+00 -1.04E+01  9.33E-01  5.01E-04
 
 TH 6
+       -1.16E+01  0.00E+00  3.53E+02 -6.26E+01  4.59E-03  8.29E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.46E+01  0.00E+00  3.71E+02 -2.87E+01 -7.77E-03  3.36E-01  0.00E+00  0.00E+00  0.00E+00  0.00E+00  6.48E+00
 
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
+        2.42E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.22E+04  0.00E+00  1.35E+06
 
 TH 4
+       -2.67E+02  0.00E+00 -1.43E+05  1.62E+04
 
 TH 5
+        4.79E-01  0.00E+00 -1.43E+01  1.31E+00  2.20E-04
 
 TH 6
+        9.32E+01  0.00E+00 -3.93E+02 -2.56E+01  3.07E-02  6.78E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        4.26E+01  0.00E+00  3.40E+02 -6.94E+01  2.05E-03  5.87E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  5.18E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.02
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       49.885
Stop Time:
Wed Sep 29 19:55:49 CDT 2021
