Sat Sep 18 14:00:19 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1703.94310234915        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.8361E+01 -2.7966E+01 -4.4589E+01  1.1830E+01  6.7482E+01  1.1161E+01 -2.8024E+00  1.1621E+01  1.0565E+00 -1.4900E+01
             1.1300E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1711.57271147479        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  1.0644E+00  1.0163E+00  1.1176E+00  9.9716E-01  1.0086E+00  9.3085E-01  1.0160E+00  8.9179E-01  1.0070E+00  1.1035E+00
             9.4492E-01
 PARAMETER:  1.6237E-01  1.1620E-01  2.1116E-01  9.7153E-02  1.0857E-01  2.8342E-02  1.1585E-01 -1.4520E-02  1.0695E-01  1.9848E-01
             4.3343E-02
 GRADIENT:   1.7048E+01  2.1893E+00  1.4202E+00 -8.1741E+00 -1.0632E+01 -1.6358E+01  9.3701E-01  4.8008E+00 -9.5419E-01 -5.8722E+00
            -1.2729E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1713.07413516702        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      356
 NPARAMETR:  1.0589E+00  9.1835E-01  1.2274E+00  1.0667E+00  1.0205E+00  9.5084E-01  9.4040E-01  7.6476E-01  1.0028E+00  1.1991E+00
             9.5953E-01
 PARAMETER:  1.5721E-01  1.4825E-02  3.0492E-01  1.6455E-01  1.2032E-01  4.9593E-02  3.8552E-02 -1.6820E-01  1.0283E-01  2.8155E-01
             5.8687E-02
 GRADIENT:   6.3999E+00  5.8388E+00  4.8736E+00  1.2304E+00 -8.7464E+00 -6.7744E+00 -1.6075E-01 -8.3768E-01  3.2777E-01  3.0791E-01
            -5.5336E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1713.50798330184        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  1.0521E+00  7.5451E-01  1.3919E+00  1.1799E+00  1.0248E+00  9.6551E-01  1.0372E+00  8.9180E-01  9.3893E-01  1.2312E+00
             9.7649E-01
 PARAMETER:  1.5081E-01 -1.8168E-01  4.3069E-01  2.6541E-01  1.2446E-01  6.4896E-02  1.3657E-01 -1.4518E-02  3.6989E-02  3.0802E-01
             7.6211E-02
 GRADIENT:  -4.5753E+00  6.5017E+00  2.2276E+00  9.9800E+00 -3.2774E+00  2.7230E-01  1.0482E+00 -6.1611E-01  1.8585E+00  1.8631E+00
             2.6815E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1714.05528868982        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  1.0509E+00  5.2466E-01  1.4944E+00  1.3172E+00  9.8497E-01  9.6177E-01  1.0313E+00  9.6774E-01  8.6461E-01  1.2126E+00
             9.7147E-01
 PARAMETER:  1.4961E-01 -5.4501E-01  5.0172E-01  3.7553E-01  8.4854E-02  6.1016E-02  1.3083E-01  6.7205E-02 -4.5475E-02  2.9273E-01
             7.1052E-02
 GRADIENT:  -8.3865E-01  4.6754E-01 -1.1805E-01 -2.0952E+00  5.7077E-01 -1.4295E-01  7.5300E-02 -1.8705E-01 -1.1595E+00 -2.1978E-01
             9.5400E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1714.32459974039        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0490E+00  3.4180E-01  1.6314E+00  1.4370E+00  9.7114E-01  9.6064E-01  8.4041E-01  1.0991E+00  8.1807E-01  1.2179E+00
             9.7083E-01
 PARAMETER:  1.4783E-01 -9.7353E-01  5.8942E-01  4.6257E-01  7.0716E-02  5.9842E-02 -7.3865E-02  1.9451E-01 -1.0081E-01  2.9715E-01
             7.0397E-02
 GRADIENT:   1.0572E+00  1.3654E+00 -1.9056E-01  5.1494E+00 -8.8521E-01  3.6714E-01  2.6737E-02  3.4673E-01 -7.7432E-01 -1.0260E-01
            -3.5695E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1714.49165432901        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  1.0463E+00  1.9046E-01  1.7306E+00  1.5353E+00  9.5562E-01  9.5787E-01  5.1574E-01  1.1812E+00  7.7141E-01  1.2144E+00
             9.7140E-01
 PARAMETER:  1.4522E-01 -1.5583E+00  6.4846E-01  5.2872E-01  5.4610E-02  5.6959E-02 -5.6215E-01  2.6656E-01 -1.5954E-01  2.9425E-01
             7.0982E-02
 GRADIENT:   1.1167E-01  1.1151E+00  7.3226E-01  9.2697E+00 -2.1281E+00 -4.0276E-04  1.0864E-02 -8.8589E-02 -1.2331E+00 -1.5081E-02
            -2.8980E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1714.57254425291        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  1.0441E+00  1.0543E-01  1.7920E+00  1.5882E+00  9.4982E-01  9.5588E-01  2.5996E-01  1.2348E+00  7.4635E-01  1.2156E+00
             9.7231E-01
 PARAMETER:  1.4320E-01 -2.1497E+00  6.8332E-01  5.6260E-01  4.8517E-02  5.4875E-02 -1.2472E+00  3.1087E-01 -1.9256E-01  2.9527E-01
             7.1923E-02
 GRADIENT:  -1.5629E+00  3.0136E-01  6.1835E-01  4.2129E+00 -1.5969E+00 -4.0239E-01  2.0909E-03 -4.3250E-02 -7.7982E-02  4.5679E-01
             3.5112E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1714.60522710812        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  1.0437E+00  5.5405E-02  1.8150E+00  1.6191E+00  9.4338E-01  9.5606E-01  1.1242E-01  1.2537E+00  7.3052E-01  1.2098E+00
             9.7184E-01
 PARAMETER:  1.4277E-01 -2.7931E+00  6.9611E-01  5.8190E-01  4.1718E-02  5.5061E-02 -2.0855E+00  3.2611E-01 -2.1400E-01  2.9041E-01
             7.1439E-02
 GRADIENT:  -6.4479E-01  9.2378E-02  3.9408E-01  2.1492E+00 -6.0523E-01 -9.0369E-02  1.2972E-04 -2.0871E-01 -1.7455E-01  7.2357E-02
             1.1349E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1714.62109746647        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1585
 NPARAMETR:  1.0432E+00  2.3508E-02  1.8384E+00  1.6396E+00  9.4128E-01  9.5587E-01  3.3476E-02  1.2752E+00  7.2087E-01  1.2080E+00
             9.7170E-01
 PARAMETER:  1.4234E-01 -3.6504E+00  7.0887E-01  5.9444E-01  3.9485E-02  5.4870E-02 -3.2969E+00  3.4307E-01 -2.2729E-01  2.8896E-01
             7.1288E-02
 GRADIENT:  -4.8432E-01  2.9754E-02  2.5348E-01  1.4962E+00 -3.1797E-01 -1.3688E-02  2.3030E-06 -1.6105E-01 -1.1007E-01 -3.6028E-02
             7.9984E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1714.62720669576        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1761
 NPARAMETR:  1.0434E+00  1.0036E-02  1.8550E+00  1.6486E+00  9.4218E-01  9.5579E-01  1.0000E-02  1.2932E+00  7.1669E-01  1.2083E+00
             9.7113E-01
 PARAMETER:  1.4244E-01 -4.5016E+00  7.1790E-01  5.9990E-01  4.0442E-02  5.4781E-02 -4.5462E+00  3.5712E-01 -2.3312E-01  2.8919E-01
             7.0703E-02
 GRADIENT:   3.5770E-01  1.3908E-02 -2.9909E-03  1.1689E+00 -1.0185E-01  2.7942E-02  0.0000E+00  4.8988E-02 -1.5294E-01 -7.2367E-02
            -1.1630E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1714.62761273903        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1952             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0433E+00  1.0000E-02  1.8547E+00  1.6480E+00  9.4220E-01  9.5575E-01  1.0000E-02  1.2917E+00  7.1695E-01  1.2091E+00
             9.7134E-01
 PARAMETER:  1.4234E-01 -4.5109E+00  7.1771E-01  5.9955E-01  4.0460E-02  5.4745E-02 -4.5514E+00  3.5599E-01 -2.3275E-01  2.8987E-01
             7.0924E-02
 GRADIENT:   7.0060E+01  0.0000E+00  1.0469E+00  1.4021E+02  6.8605E-01  5.1362E+00  0.0000E+00  9.6151E-02  2.2022E+00  3.7080E-01
             8.6050E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1714.62763187988        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2145             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0433E+00  1.0000E-02  1.8540E+00  1.6479E+00  9.4201E-01  9.5579E-01  1.0000E-02  1.2912E+00  7.1697E-01  1.2090E+00
             9.7132E-01
 PARAMETER:  1.4235E-01 -4.5109E+00  7.1736E-01  5.9950E-01  4.0266E-02  5.4783E-02 -4.5514E+00  3.5560E-01 -2.3272E-01  2.8982E-01
             7.0898E-02
 GRADIENT:   7.0082E+01  0.0000E+00  1.0581E+00  1.4007E+02  6.5710E-01  5.1561E+00  0.0000E+00  9.7016E-02  2.2156E+00  3.8317E-01
             7.6832E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1714.62764421099        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     2313
 NPARAMETR:  1.0432E+00  1.0000E-02  1.8537E+00  1.6480E+00  9.4198E-01  9.5575E-01  1.0000E-02  1.2908E+00  7.1697E-01  1.2089E+00
             9.7134E-01
 PARAMETER:  1.4230E-01 -4.5109E+00  7.1717E-01  5.9955E-01  4.0228E-02  5.4738E-02 -4.5514E+00  3.5530E-01 -2.3272E-01  2.8971E-01
             7.0919E-02
 GRADIENT:  -9.0045E-02  0.0000E+00  3.9731E-02 -1.3695E-01 -4.4717E-03 -5.5534E-03  0.0000E+00  3.3404E-03 -6.2831E-03  1.6116E-02
             7.6822E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2313
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6183E-04 -3.7140E-06 -2.9365E-02 -7.4393E-03 -3.8310E-02
 SE:             2.9863E-02  1.7739E-06  1.5336E-02  2.9222E-02  2.2338E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9568E-01  3.6293E-02  5.5517E-02  7.9905E-01  8.6338E-02

 ETASHRINKSD(%)  1.0000E-10  9.9994E+01  4.8624E+01  2.1011E+00  2.5166E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  7.3605E+01  4.1580E+00  4.3999E+01
 EBVSHRINKSD(%)  4.1086E-01  9.9994E+01  5.2732E+01  2.3674E+00  2.0546E+01
 EBVSHRINKVR(%)  8.2004E-01  1.0000E+02  7.7657E+01  4.6787E+00  3.6871E+01
 RELATIVEINF(%)  9.7176E+01  1.6792E-08  5.5628E+00  5.7528E+00  1.1779E+01
 EPSSHRINKSD(%)  4.4485E+01
 EPSSHRINKVR(%)  6.9181E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.6276442109950     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -979.47681764725678     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.628       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  1.85E+00  1.65E+00  9.42E-01  9.56E-01  1.00E-02  1.29E+00  7.17E-01  1.21E+00  9.71E-01
 


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
+        1.11E+03
 
 TH 2
+        0.00E+00  1.14E+03
 
 TH 3
+       -1.97E+00  0.00E+00  4.79E+01
 
 TH 4
+       -8.60E+00  0.00E+00 -2.43E+01  7.60E+02
 
 TH 5
+        1.72E+01  0.00E+00 -1.25E+02 -4.36E+01  5.85E+02
 
 TH 6
+       -1.80E+01  0.00E+00  5.83E-01 -1.98E+00  8.14E+00  2.20E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.42E+00  0.00E+00 -1.56E+01 -3.91E+00 -9.18E+00  4.73E+00  0.00E+00  2.18E+01
 
 TH 9
+       -1.32E+00  0.00E+00  5.60E+00 -3.09E-01 -1.66E+00  1.76E+01  0.00E+00 -7.36E-01  3.63E+02
 
 TH10
+        1.22E-02  0.00E+00 -2.84E+00 -2.55E+00 -7.37E+01 -7.08E+00  0.00E+00  8.88E+00  4.64E+00  6.17E+01
 
 TH11
+       -1.87E+01  0.00E+00 -6.85E+00 -1.18E+01 -2.09E+01 -4.21E+00  0.00E+00  4.48E+00  1.43E+00  2.08E+01  2.91E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.359
Stop Time:
Sat Sep 18 14:00:56 CDT 2021
