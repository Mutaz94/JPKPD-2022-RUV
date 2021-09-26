Sat Sep 25 14:23:09 CDT 2021
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
$DATA ../../../../data/spa/D/dat48.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14069.9132915831        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6979E+02  2.4684E+02 -2.2666E+01  1.9998E+02  1.5800E+02 -1.7675E+03 -6.6775E+02 -5.6939E+01 -1.0290E+03 -4.9405E+02
            -2.7376E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -625.858822710172        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3861E+00  1.0897E+00  9.5006E-01  1.4519E+00  1.2276E+00  1.5393E+00  1.1562E+00  9.7028E-01  1.0965E+00  1.1210E+00
             1.5023E+01
 PARAMETER:  4.2648E-01  1.8594E-01  4.8773E-02  4.7289E-01  3.0502E-01  5.3130E-01  2.4515E-01  6.9834E-02  1.9214E-01  2.1419E-01
             2.8096E+00
 GRADIENT:   1.3647E+01  2.0184E+00 -5.6208E+00  1.7068E+00 -4.6700E+00  2.6854E+01  2.6996E+00  3.5061E+00  1.2492E+01  3.5023E+00
             1.3825E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -636.008571440516        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.3672E+00  9.4894E-01  1.6688E+00  1.6766E+00  2.6354E+00  1.3867E+00  2.2905E+00  4.4037E-01  9.7749E-01  5.2336E+00
             1.3788E+01
 PARAMETER:  4.1279E-01  4.7586E-02  6.1207E-01  6.1679E-01  1.0690E+00  4.2694E-01  9.2878E-01 -7.2014E-01  7.7235E-02  1.7551E+00
             2.7238E+00
 GRADIENT:   2.8456E+01  2.2892E+01  1.2550E+00  4.9852E+01 -5.1703E+00 -1.0896E+01  4.3462E+00  9.7459E-02  5.6360E+00  2.0397E+00
             6.8439E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -663.941641058925        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1297E+00  3.9190E-01  1.0585E+00  1.4621E+00  3.8230E+00  1.3652E+00  1.4348E+00  5.2285E-02  3.5508E-01  7.2568E+00
             1.1505E+01
 PARAMETER:  2.2199E-01 -8.3676E-01  1.5683E-01  4.7990E-01  1.4410E+00  4.1132E-01  4.6105E-01 -2.8511E+00 -9.3541E-01  2.0819E+00
             2.5428E+00
 GRADIENT:   2.1653E+01  2.2289E+00  7.8689E+00 -2.0213E+01 -8.6951E+00  1.5877E+01  2.1333E-01  5.4236E-04 -1.7428E-02  2.7959E+00
            -1.5597E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -676.767080036576        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.8842E-01  1.8153E-01  6.4176E-01  1.2954E+00  1.2023E+01  1.1522E+00  2.8380E-01  1.0000E-02  7.5340E-02  8.2607E+00
             1.1289E+01
 PARAMETER:  8.8350E-02 -1.6063E+00 -3.4354E-01  3.5880E-01  2.5868E+00  2.4170E-01 -1.1595E+00 -5.2492E+00 -2.4857E+00  2.2115E+00
             2.5239E+00
 GRADIENT:  -2.0423E+01 -4.9317E+00  4.0138E+01 -1.0815E+02  8.2956E+00 -2.1885E+01  1.5730E-02  0.0000E+00  2.8825E-01 -1.7626E+01
             1.4219E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -684.007320782510        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.0189E+00  1.5449E-01  6.1947E-01  1.4035E+00  1.2318E+01  1.2179E+00  3.3619E-01  1.0000E-02  6.8288E-02  9.4297E+00
             1.1222E+01
 PARAMETER:  1.1876E-01 -1.7677E+00 -3.7889E-01  4.3893E-01  2.6111E+00  2.9713E-01 -9.9007E-01 -5.4802E+00 -2.5840E+00  2.3439E+00
             2.5178E+00
 GRADIENT:  -1.3931E+01 -1.1371E+00  2.1427E+01 -4.5175E+01  8.0997E+00 -2.5372E+01  1.5537E-02  0.0000E+00  2.3043E-01 -7.0701E+00
            -6.7111E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -686.038564179577        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:      510             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0343E+00  1.3515E-01  5.9576E-01  1.4490E+00  1.4025E+01  1.2817E+00  3.0484E-01  1.0000E-02  5.4106E-02  1.0754E+01
             1.1267E+01
 PARAMETER:  1.3371E-01 -1.9014E+00 -4.1791E-01  4.7086E-01  2.7409E+00  3.4816E-01 -1.0880E+00 -5.9224E+00 -2.8168E+00  2.4753E+00
             2.5218E+00
 GRADIENT:  -5.8179E+00  3.7318E-01  1.2563E+01 -1.9155E+01  5.5309E+00 -1.6942E+01  7.4585E-03  0.0000E+00  1.2305E-01 -1.3532E+00
            -1.2222E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -686.227506105059        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      593
 NPARAMETR:  1.0449E+00  1.3521E-01  5.9562E-01  1.4494E+00  1.4021E+01  1.2819E+00  3.0346E-01  1.0000E-02  1.0000E-02  1.0751E+01
             1.1373E+01
 PARAMETER:  1.4388E-01 -1.9009E+00 -4.1816E-01  4.7113E-01  2.7406E+00  3.4837E-01 -1.0925E+00 -5.9224E+00 -6.6386E+00  2.4750E+00
             2.5312E+00
 GRADIENT:   2.4970E+00 -3.4193E-01  1.3976E+01 -2.9589E+01  7.4103E+00 -1.7522E+01  8.8347E-03  0.0000E+00  0.0000E+00 -2.5525E+00
            -2.7889E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -686.250904430501        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      668
 NPARAMETR:  1.0420E+00  1.3517E-01  5.9558E-01  1.4497E+00  1.4003E+01  1.2820E+00  2.6747E-01  1.0000E-02  1.0000E-02  1.0747E+01
             1.1430E+01
 PARAMETER:  1.4112E-01 -1.9013E+00 -4.1822E-01  4.7134E-01  2.7392E+00  3.4843E-01 -1.2187E+00 -5.9224E+00 -9.0234E+00  2.4746E+00
             2.5362E+00
 GRADIENT:  -2.3934E+00 -3.2459E-01  1.4470E+01 -2.9636E+01  7.6155E+00 -1.7364E+01  7.0265E-03  0.0000E+00  0.0000E+00 -2.7245E+00
             3.3011E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -686.351112715575        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      746
 NPARAMETR:  1.0351E+00  1.3325E-01  5.9640E-01  1.4525E+00  1.3629E+01  1.2806E+00  1.5441E-02  1.0000E-02  1.0000E-02  1.0669E+01
             1.1476E+01
 PARAMETER:  1.3451E-01 -1.9155E+00 -4.1684E-01  4.7330E-01  2.7122E+00  3.4730E-01 -4.0708E+00 -5.9224E+00 -8.4588E+00  2.4673E+00
             2.5403E+00
 GRADIENT:  -1.0938E+01  2.7034E-01  1.3663E+01 -2.2274E+01  5.6012E+00 -1.6468E+01  1.9775E-05  0.0000E+00  0.0000E+00 -1.4092E+00
             6.0229E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -686.403820551968        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  1.0418E+00  1.3270E-01  5.9656E-01  1.4537E+00  1.3511E+01  1.2803E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0642E+01
             1.1567E+01
 PARAMETER:  1.4092E-01 -1.9197E+00 -4.1658E-01  4.7411E-01  2.7035E+00  3.4710E-01 -4.9702E+00 -5.9224E+00 -1.0775E+01  2.4648E+00
             2.5482E+00
 GRADIENT:  -3.0316E+00  7.1219E-01  1.1264E+01 -2.1045E+01  9.8925E-01 -1.4190E+01  0.0000E+00  0.0000E+00  0.0000E+00  3.1744E+00
             6.5080E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -686.406718384898        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0427E+00  1.3282E-01  5.9650E-01  1.4535E+00  1.3537E+01  1.2804E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0645E+01
             1.1559E+01
 PARAMETER:  1.4179E-01 -1.9188E+00 -4.1667E-01  4.7399E-01  2.7054E+00  3.4717E-01 -4.7780E+00 -5.9224E+00 -1.0612E+01  2.4651E+00
             2.5474E+00
 GRADIENT:  -6.4979E-01  9.2404E-01  9.9233E+00 -1.8876E+01 -5.5177E-01 -1.2226E+01  0.0000E+00  0.0000E+00  0.0000E+00  4.3685E+00
             3.0392E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -686.407938893908        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     1011
 NPARAMETR:  1.0434E+00  1.3281E-01  5.9650E-01  1.4536E+00  1.3534E+01  1.2804E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0644E+01
             1.1564E+01
 PARAMETER:  1.4252E-01 -1.9189E+00 -4.1667E-01  4.7403E-01  2.7052E+00  3.4718E-01 -4.8034E+00 -5.9224E+00 -1.0759E+01  2.4650E+00
             2.5479E+00
 GRADIENT:   6.3051E+01  1.5777E+01 -5.5753E+01  1.2079E+02 -9.9902E+01  5.6908E+01  0.0000E+00  0.0000E+00  0.0000E+00  8.7563E+01
            -1.3323E+02

0ITERATION NO.:   64    OBJECTIVE VALUE:  -686.409821259033        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1151
 NPARAMETR:  1.0449E+00  1.3256E-01  5.9672E-01  1.4530E+00  1.3575E+01  1.2800E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0615E+01
             1.1597E+01
 PARAMETER:  1.4402E-01 -1.9188E+00 -4.1672E-01  4.7407E-01  2.7055E+00  3.4721E-01 -4.7804E+00 -5.9224E+00 -1.0938E+01  2.4647E+00
             2.5483E+00
 GRADIENT:   5.3392E+03  2.0044E+02 -1.8322E+03  1.5900E+03 -5.5867E+02  2.1987E+03  0.0000E+00  0.0000E+00  0.0000E+00  6.1745E+02
            -5.8825E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1151
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1643E-02 -1.9313E-05 -6.2338E-05 -4.3743E-04 -3.4326E-02
 SE:             2.8465E-02  9.8463E-06  6.8414E-05  2.8213E-04  8.3452E-03
 N:                     100         100         100         100         100

 P VAL.:         6.8253E-01  4.9822E-02  3.6220E-01  1.2103E-01  3.9043E-05

 ETASHRINKSD(%)  4.6389E+00  9.9967E+01  9.9771E+01  9.9055E+01  7.2042E+01
 ETASHRINKVR(%)  9.0627E+00  1.0000E+02  9.9999E+01  9.9991E+01  9.2184E+01
 EBVSHRINKSD(%)  7.0314E+00  9.9963E+01  9.9735E+01  9.9025E+01  7.3442E+01
 EBVSHRINKVR(%)  1.3568E+01  1.0000E+02  9.9999E+01  9.9990E+01  9.2947E+01
 RELATIVEINF(%)  3.0238E+01  3.6015E-07  1.0043E-04  1.8776E-04  5.4920E+00
 EPSSHRINKSD(%)  4.6483E+00
 EPSSHRINKVR(%)  9.0805E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -686.40982125903258     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       48.741005304705595     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -686.410       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.33E-01  5.96E-01  1.45E+00  1.35E+01  1.28E+00  1.00E-02  1.00E-02  1.00E-02  1.06E+01  1.16E+01
 


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
+        8.49E+06
 
 TH 2
+       -1.31E+03  2.96E+06
 
 TH 3
+        1.39E+03  1.54E+03  3.11E+06
 
 TH 4
+        1.83E+06  1.08E+06 -1.11E+06  8.05E+05
 
 TH 5
+        7.93E+00  1.07E+01  1.33E+01  2.10E+01  2.84E+02
 
 TH 6
+       -9.77E+02 -1.64E+03  1.69E+03 -6.21E+05  6.03E+00  9.72E+05
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -9.72E+00 -1.74E+01 -2.29E+01 -3.33E+01 -1.99E+02 -1.02E+01  0.00E+00  0.00E+00  0.00E+00  5.58E+02
 
 TH11
+       -1.32E+00  9.88E+00  1.97E+01  1.65E+01  1.76E+02  9.34E+00  0.00E+00  0.00E+00  0.00E+00 -2.47E+02  4.42E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.507
Stop Time:
Sat Sep 25 14:23:41 CDT 2021
