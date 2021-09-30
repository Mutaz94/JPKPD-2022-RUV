Wed Sep 29 20:11:22 CDT 2021
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
$DATA ../../../../data/spa/D/dat60.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13042.9210636267        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8200E+02  2.1924E+02 -5.4655E+01  4.2866E+01  2.8524E+02 -1.9137E+03 -8.8144E+02 -4.0796E+01 -1.5133E+03 -6.9505E+02
            -2.4068E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -586.914951231228        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3381E+00  1.0531E+00  8.8442E-01  2.0754E+00  1.6210E+00  2.5671E+00  1.4408E+00  9.5185E-01  1.8095E+00  1.2126E+00
             1.3116E+01
 PARAMETER:  3.9127E-01  1.5171E-01 -2.2823E-02  8.3017E-01  5.8303E-01  1.0428E+00  4.6518E-01  5.0656E-02  6.9305E-01  2.9274E-01
             2.6739E+00
 GRADIENT:  -4.5495E+00  3.9861E+01 -1.7269E+01  7.7927E+01 -6.9222E+00  6.1461E+01 -1.3923E+00  7.2495E+00 -1.9072E+01  7.9942E-01
             7.6338E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -616.337366221574        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3165E+00  1.0479E+00  1.1581E+00  2.1527E+00  4.0179E+00  2.3983E+00  3.4150E+00  4.0431E-01  3.2816E+00  6.9118E+00
             1.0331E+01
 PARAMETER:  3.7495E-01  1.4681E-01  2.4676E-01  8.6672E-01  1.4908E+00  9.7478E-01  1.3282E+00 -8.0557E-01  1.2883E+00  2.0332E+00
             2.4352E+00
 GRADIENT:   2.0311E+01  2.7621E+01 -1.2726E+01  6.2504E+01 -7.7436E+00  6.8252E+00  1.5090E+01  1.6029E-01  5.5894E+01 -1.3156E-01
             5.2402E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -649.490706128987        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  1.1330E+00  3.8934E-01  7.7861E-01  1.7196E+00  6.3855E+00  2.0989E+00  1.2384E+00  2.7881E-02  1.6381E+00  8.3178E+00
             1.0557E+01
 PARAMETER:  2.2484E-01 -8.4329E-01 -1.5024E-01  6.4209E-01  1.9540E+00  8.4141E-01  3.1381E-01 -3.4798E+00  5.9356E-01  2.2184E+00
             2.4567E+00
 GRADIENT:  -1.8044E+01  9.9427E+00  2.2693E+01  1.9760E+01 -5.7659E+00 -3.7144E+00  1.1955E+00 -5.4944E-03 -1.5906E+01  9.9450E+00
             5.1623E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -732.073956452425        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  7.9883E-01  7.0405E-02  6.9304E-02  6.2682E-01  1.0587E+02  1.8020E+00  1.8307E-02  1.0000E-02  9.2604E-01  9.8691E+00
             9.2563E+00
 PARAMETER: -1.2461E-01 -2.5535E+00 -2.5693E+00 -3.6709E-01  4.7622E+00  6.8892E-01 -3.9005E+00 -9.4943E+00  2.3167E-02  2.3894E+00
             2.3253E+00
 GRADIENT:   9.5977E+01  8.9671E+01 -3.8797E+01 -2.2539E+01 -8.3725E-01 -1.1120E+01  8.1818E-03  0.0000E+00 -3.1175E+01  1.6223E-01
            -2.6721E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -743.908993655858        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      417
 NPARAMETR:  7.1812E-01  5.1412E-02  5.5455E-02  5.4596E-01  1.3616E+02  1.8303E+00  1.0705E-02  1.0000E-02  9.2861E-01  9.0986E+00
             8.8662E+00
 PARAMETER: -2.3111E-01 -2.8679E+00 -2.7922E+00 -5.0521E-01  5.0138E+00  7.0450E-01 -4.4370E+00 -9.8374E+00  2.5935E-02  2.3081E+00
             2.2822E+00
 GRADIENT:   7.5673E+01  4.4777E+01 -5.8991E+01  1.1973E+01 -3.2179E-01 -1.5078E+01  1.7800E-03  0.0000E+00 -3.3755E+01  5.5124E-02
            -5.8377E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -767.326383063545        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      592
 NPARAMETR:  4.2538E-01  1.0000E-02  1.9728E-02  2.4582E-01  7.4076E+02  1.7445E+00  1.0000E-02  1.0000E-02  1.0409E+00  8.8020E+00
             9.0581E+00
 PARAMETER: -7.5478E-01 -4.6563E+00 -3.8257E+00 -1.3032E+00  6.7077E+00  6.5647E-01 -9.4189E+00 -1.2949E+01  1.4004E-01  2.2750E+00
             2.3037E+00
 GRADIENT:   1.5964E+00  0.0000E+00 -8.7571E+00  9.1209E+00  1.6156E-03 -2.6708E+00  0.0000E+00  0.0000E+00  6.8045E-01  6.3812E-07
             1.1503E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -767.372113558028        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:      791             RESET HESSIAN, TYPE I
 NPARAMETR:  4.3222E-01  1.0000E-02  2.0743E-02  2.5426E-01  5.1679E+02  1.7669E+00  1.0000E-02  1.0000E-02  1.0378E+00  8.9952E+00
             8.8992E+00
 PARAMETER: -7.3882E-01 -4.5730E+00 -3.7755E+00 -1.2694E+00  6.3476E+00  6.6921E-01 -9.1291E+00 -1.2739E+01  1.3707E-01  2.2967E+00
             2.2860E+00
 GRADIENT:   6.6813E+01  0.0000E+00  8.1820E+01  3.9286E+01  2.4535E-03  1.5142E+01  0.0000E+00  0.0000E+00  8.1136E-01  1.3934E-05
             1.7920E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -767.386148383151        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  4.2919E-01  1.0000E-02  2.0487E-02  2.5254E-01  2.8873E+02  1.7622E+00  1.0000E-02  1.0000E-02  1.0382E+00  8.1220E+00
             8.8728E+00
 PARAMETER: -7.4585E-01 -4.5730E+00 -3.7880E+00 -1.2762E+00  5.7655E+00  6.6658E-01 -9.1291E+00 -1.2739E+01  1.3746E-01  2.1946E+00
             2.2830E+00
 GRADIENT:   4.9216E-01  0.0000E+00 -5.6660E+00  6.4626E+00  4.3663E-03 -6.1041E-01  0.0000E+00  0.0000E+00 -5.0950E-01  1.2208E-06
            -2.3616E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -767.412421456937        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1125
 NPARAMETR:  4.2831E-01  1.0000E-02  2.0445E-02  2.5136E-01  1.9152E+02  1.7625E+00  1.0000E-02  1.0000E-02  1.0381E+00  8.1071E+00
             8.8880E+00
 PARAMETER: -7.4790E-01 -4.5730E+00 -3.7900E+00 -1.2809E+00  5.3550E+00  6.6674E-01 -9.1291E+00 -1.2739E+01  1.3738E-01  2.1927E+00
             2.2847E+00
 GRADIENT:  -5.4257E-01  0.0000E+00 -1.9772E+00  2.2903E+00  5.8895E-03 -3.9282E-01  0.0000E+00  0.0000E+00 -1.8091E-01  6.1372E-06
            -7.1144E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -767.429288200722        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1281
 NPARAMETR:  4.2629E-01  1.0000E-02  2.0254E-02  2.4979E-01  1.5393E+02  1.7608E+00  1.0000E-02  1.0000E-02  1.0380E+00  7.5609E+00
             8.8812E+00
 PARAMETER: -7.5263E-01 -4.5730E+00 -3.7994E+00 -1.2871E+00  5.1365E+00  6.6577E-01 -9.1291E+00 -1.2739E+01  1.3726E-01  2.1230E+00
             2.2839E+00
 GRADIENT:   6.5400E+01  0.0000E+00  8.2643E+01  4.1731E+01  7.2045E-03  1.4341E+01  0.0000E+00  0.0000E+00  7.5536E-01  1.1288E-04
             1.7165E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -767.449196912306        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1436
 NPARAMETR:  4.2736E-01  1.0000E-02  2.0200E-02  2.4793E-01  1.3885E+02  1.7644E+00  1.0000E-02  1.0000E-02  1.0410E+00  6.4619E+00
             8.8903E+00
 PARAMETER: -7.5014E-01 -4.5730E+00 -3.8021E+00 -1.2946E+00  5.0334E+00  6.6779E-01 -9.1291E+00 -1.2739E+01  1.4022E-01  1.9659E+00
             2.2850E+00
 GRADIENT:   1.1789E+00  0.0000E+00  1.8455E+00 -3.6124E+00  7.6824E-03  3.3031E-01  0.0000E+00  0.0000E+00  7.9541E-01  1.3106E-05
            -3.1695E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -767.474450150097        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1614
 NPARAMETR:  4.2639E-01  1.0000E-02  2.0240E-02  2.4852E-01  1.7587E+01  1.7588E+00  1.0000E-02  1.0000E-02  1.0452E+00  2.8621E+00
             8.8756E+00
 PARAMETER: -7.5241E-01 -4.5730E+00 -3.8001E+00 -1.2922E+00  2.9672E+00  6.6465E-01 -9.1291E+00 -1.2739E+01  1.4425E-01  1.1516E+00
             2.2833E+00
 GRADIENT:  -2.2951E+00  0.0000E+00  3.1275E+00 -3.6311E+00  2.1081E-02 -1.0209E+00  0.0000E+00  0.0000E+00  1.7614E+00  8.9172E-04
            -6.7001E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -767.508180267064        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1801             RESET HESSIAN, TYPE I
 NPARAMETR:  4.2846E-01  1.0000E-02  2.0148E-02  2.4838E-01  1.3960E+01  1.7643E+00  1.0000E-02  1.0000E-02  1.0371E+00  5.1934E-01
             8.8914E+00
 PARAMETER: -7.4756E-01 -4.5730E+00 -3.8047E+00 -1.2928E+00  2.7362E+00  6.6778E-01 -9.1291E+00 -1.2739E+01  1.3643E-01 -5.5519E-01
             2.2851E+00
 GRADIENT:   6.7507E+01  0.0000E+00  8.4726E+01  3.8220E+01  1.4094E-02  1.5152E+01  0.0000E+00  0.0000E+00  1.0117E+00  1.1702E-04
             1.7788E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -767.515007429979        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1995             RESET HESSIAN, TYPE I
 NPARAMETR:  4.2822E-01  1.0000E-02  2.0056E-02  2.4764E-01  1.3514E+01  1.7642E+00  1.0000E-02  1.0000E-02  1.0369E+00  4.6871E-01
             8.8893E+00
 PARAMETER: -7.4812E-01 -4.5730E+00 -3.8092E+00 -1.2958E+00  2.7037E+00  6.6772E-01 -9.1291E+00 -1.2739E+01  1.3621E-01 -6.5778E-01
             2.2848E+00
 GRADIENT:   6.8376E+01  0.0000E+00  8.3994E+01  3.9115E+01  4.9822E-03  1.5213E+01  0.0000E+00  0.0000E+00  9.1351E-01  1.2478E-04
             1.7467E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -767.531907286557        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2182
 NPARAMETR:  4.2635E-01  1.0000E-02  1.9923E-02  2.4593E-01  1.4484E+01  1.7628E+00  1.0000E-02  1.0000E-02  1.0369E+00  3.8278E-02
             8.9003E+00
 PARAMETER: -7.5250E-01 -4.5730E+00 -3.8159E+00 -1.3027E+00  2.7730E+00  6.6689E-01 -9.1291E+00 -1.2739E+01  1.3625E-01 -3.1629E+00
             2.2861E+00
 GRADIENT:   9.9776E-01  0.0000E+00 -8.2383E-01 -8.3847E-01  7.0391E-03  2.2080E-01  0.0000E+00  0.0000E+00  4.3910E-02  4.8527E-07
             4.0770E-01

0ITERATION NO.:   79    OBJECTIVE VALUE:  -767.535810977850        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     2318
 NPARAMETR:  4.2590E-01  1.0000E-02  1.9884E-02  2.4564E-01  1.4078E+01  1.7622E+00  1.0000E-02  1.0000E-02  1.0369E+00  3.7909E-02
             8.8902E+00
 PARAMETER: -7.5356E-01 -4.5730E+00 -3.8178E+00 -1.3039E+00  2.7446E+00  6.6654E-01 -9.1291E+00 -1.2739E+01  1.3622E-01 -3.1726E+00
             2.2849E+00
 GRADIENT:   8.7871E-01  0.0000E+00 -1.2146E+00 -3.3456E-01 -1.8799E-03  9.6969E-02  0.0000E+00  0.0000E+00 -2.1570E-02  5.8203E-07
            -3.7731E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2318
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.8712E-04  3.4945E-07  9.8358E-05 -2.2546E-02 -1.1305E-06
 SE:             2.9095E-02  1.8089E-06  2.1534E-04  2.4603E-02  2.8273E-06
 N:                     100         100         100         100         100

 P VAL.:         9.8116E-01  8.4682E-01  6.4784E-01  3.5946E-01  6.8927E-01

 ETASHRINKSD(%)  2.5276E+00  9.9994E+01  9.9279E+01  1.7577E+01  9.9991E+01
 ETASHRINKVR(%)  4.9912E+00  1.0000E+02  9.9995E+01  3.2065E+01  1.0000E+02
 EBVSHRINKSD(%)  2.3151E+00  9.9990E+01  9.9330E+01  1.7538E+01  9.9989E+01
 EBVSHRINKVR(%)  4.5766E+00  1.0000E+02  9.9996E+01  3.2000E+01  1.0000E+02
 RELATIVEINF(%)  3.7230E+00  1.6922E-07  3.1595E-05  4.7960E-01  6.9330E-08
 EPSSHRINKSD(%)  1.5095E+01
 EPSSHRINKVR(%)  2.7912E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -767.53581097785013     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -32.384984414111955     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -767.536       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.26E-01  1.00E-02  1.99E-02  2.46E-01  1.41E+01  1.76E+00  1.00E-02  1.00E-02  1.04E+00  3.79E-02  8.89E+00
 


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
+        1.87E+03
 
 TH 2
+        0.00E+00  9.45E+03
 
 TH 3
+       -1.43E+04  0.00E+00  1.20E+06
 
 TH 4
+       -1.46E+02  0.00E+00 -1.08E+05  1.11E+04
 
 TH 5
+        1.98E-01  0.00E+00 -3.25E+00  2.75E-01  8.44E-04
 
 TH 6
+        6.98E-01  0.00E+00 -1.01E+02 -2.40E+01  1.80E-03  5.59E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.77E-01  0.00E+00  1.55E+03 -1.60E+02 -9.47E-03 -1.66E+00  0.00E+00  0.00E+00  9.98E+01
 
 TH10
+       -1.14E-03  0.00E+00  2.94E-02 -1.86E-02 -2.94E-05  3.31E-03  0.00E+00  0.00E+00  7.35E-03 -6.01E-03
 
 TH11
+       -1.81E+01  0.00E+00  2.67E+02 -1.44E+01 -2.42E-03  1.04E+00  0.00E+00  0.00E+00  2.12E+00 -1.60E-04  4.58E+00
 
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
 #CPUT: Total CPU Time in Seconds,       45.233
Stop Time:
Wed Sep 29 20:12:09 CDT 2021
