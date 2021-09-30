Wed Sep 29 20:21:37 CDT 2021
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
$DATA ../../../../data/spa/D/dat82.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m82.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   23358.7692145487        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.3636E+02  4.0055E+02 -7.4127E+01  5.1836E+02  1.1959E+02 -1.8649E+03 -1.0466E+03 -1.7552E+01 -1.3998E+03 -3.2453E+02
            -4.4796E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -526.516568054844        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1795E+00  1.0666E+00  9.2856E-01  1.3399E+00  1.6165E+00  1.1762E+00  8.9712E-01  9.2193E-01  5.6743E-01  8.3336E-01
             1.6011E+01
 PARAMETER:  2.6511E-01  1.6452E-01  2.5877E-02  3.9256E-01  5.8025E-01  2.6228E-01 -8.5639E-03  1.8710E-02 -4.6665E-01 -8.2292E-02
             2.8733E+00
 GRADIENT:  -9.2178E+01  1.6065E+01 -1.2730E+01  3.0547E+01 -7.8389E-01 -4.3070E+00  3.8639E+00  2.4453E+00  6.6584E+00  4.3858E-01
             1.0905E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -535.893239216484        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.2705E+00  9.3156E-01  1.3075E+00  1.3859E+00  2.0294E+00  1.1349E+00  7.7147E-01  3.6403E-01  4.1461E-01  1.6990E+00
             1.5948E+01
 PARAMETER:  3.3944E-01  2.9106E-02  3.6808E-01  4.2636E-01  8.0776E-01  2.2652E-01 -1.5946E-01 -9.1052E-01 -7.8043E-01  6.3003E-01
             2.8693E+00
 GRADIENT:   2.6551E+01  1.3719E+00  8.2558E-01 -6.5457E+00 -3.0848E+00  2.3390E+00  1.8674E+00  1.4238E-01  3.4354E+00  1.2718E-01
             6.7614E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -544.752187859501        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1186E+00  5.1826E-01  1.4913E+00  1.5463E+00  6.2830E+00  1.1639E+00  1.2799E+00  1.2764E-01  3.6142E-01  1.0272E+01
             1.3960E+01
 PARAMETER:  2.1206E-01 -5.5728E-01  4.9966E-01  5.3585E-01  1.9379E+00  2.5179E-01  3.4681E-01 -1.9585E+00 -9.1771E-01  2.4295E+00
             2.7362E+00
 GRADIENT:  -3.5139E+01  1.2582E+01  1.2860E+01  3.8500E+01  9.0628E-01  8.9552E-01  1.5947E+00  6.6866E-03  3.4684E+00  3.1818E+00
            -6.7522E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -558.726697139422        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0659E+00  1.8826E-01  7.3784E-01  1.4959E+00  5.9709E+00  1.0594E+00  1.1308E+00  1.2644E-02  7.5223E-02  7.7675E+00
             1.3305E+01
 PARAMETER:  1.6380E-01 -1.5700E+00 -2.0403E-01  5.0273E-01  1.8869E+00  1.5769E-01  2.2295E-01 -4.2706E+00 -2.4873E+00  2.1499E+00
             2.6881E+00
 GRADIENT:   3.6146E+00  2.4838E+00  9.1816E+00 -9.8802E+00  1.4399E+00 -1.8352E+01  3.0711E-01  2.8497E-04  4.2071E-01 -3.3325E+00
            -2.2897E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -564.139042329905        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0329E+00  3.4582E-02  5.1048E-01  1.5087E+00  6.1891E+00  1.1269E+00  2.5651E-01  1.0000E-02  1.0000E-02  6.3804E+00
             1.3853E+01
 PARAMETER:  1.3242E-01 -3.2644E+00 -5.7240E-01  5.1126E-01  1.9228E+00  2.1948E-01 -1.2606E+00 -6.7224E+00 -5.9111E+00  1.9532E+00
             2.7285E+00
 GRADIENT:  -1.3505E+01  9.5522E-01 -4.8514E+00  3.3178E+01  2.0590E+00 -9.9707E+00  4.2704E-04  0.0000E+00  0.0000E+00  7.7357E-01
             9.1115E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -588.902219036424        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  7.4145E-01  1.0000E-02  7.7193E-02  6.6385E-01  5.9540E+00  1.3682E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.2512E+00
             1.2704E+01
 PARAMETER: -1.9915E-01 -1.2223E+01 -2.4615E+00 -3.0970E-01  1.8841E+00  4.1352E-01 -6.5521E+00 -1.9065E+01 -2.5211E+01  3.2407E-01
             2.6419E+00
 GRADIENT:   1.1440E+02  0.0000E+00 -3.2060E+01  4.8689E+01  5.1497E-01  2.6942E+01  0.0000E+00  0.0000E+00  0.0000E+00  1.3759E+00
            -7.6973E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -602.733593360587        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  6.2935E-01  1.0000E-02  5.4950E-02  5.2300E-01  7.2522E+00  1.0847E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.1590E+00
             1.3780E+01
 PARAMETER: -3.6307E-01 -1.3906E+01 -2.8013E+00 -5.4817E-01  2.0813E+00  1.8126E-01 -8.6915E+00 -2.0901E+01 -2.9386E+01  2.4754E-01
             2.7232E+00
 GRADIENT:   4.4971E+01  0.0000E+00 -2.1415E+01  2.0775E+01 -3.2554E-01 -1.8887E+01  0.0000E+00  0.0000E+00  0.0000E+00  1.1197E-01
            -1.7642E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -607.379739294457        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      790             RESET HESSIAN, TYPE I
 NPARAMETR:  4.5892E-01  1.0000E-02  2.5499E-02  2.9655E-01  1.1137E+01  1.1035E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.5682E-01
             1.3582E+01
 PARAMETER: -6.7888E-01 -1.7800E+01 -3.5691E+00 -1.1155E+00  2.5103E+00  1.9849E-01 -1.1608E+01 -2.6169E+01 -3.8475E+01 -9.3051E-01
             2.7088E+00
 GRADIENT:   3.2299E+01  0.0000E+00  4.2884E+01  1.3015E+01 -2.8394E-02  9.6845E-01  0.0000E+00  0.0000E+00  0.0000E+00  7.3619E-05
             2.2214E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -607.388241588336        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  4.5383E-01  1.0000E-02  2.5063E-02  2.9314E-01  1.2575E+01  1.0999E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.4956E-01
             1.3501E+01
 PARAMETER: -6.9004E-01 -1.7800E+01 -3.5864E+00 -1.1271E+00  2.6317E+00  1.9520E-01 -1.1608E+01 -2.6169E+01 -3.8475E+01 -1.8000E+00
             2.7027E+00
 GRADIENT:   2.9198E+01  0.0000E+00  4.0914E+01  1.8378E+01 -2.0410E-02  1.0984E-01  0.0000E+00  0.0000E+00  0.0000E+00  4.6268E-06
             1.8888E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -607.405146886103        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1044             RESET HESSIAN, TYPE I
 NPARAMETR:  4.5655E-01  1.0000E-02  2.5135E-02  2.9397E-01  1.7661E+01  1.1052E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.8417E-01
             1.3552E+01
 PARAMETER: -6.8406E-01 -1.7800E+01 -3.5835E+00 -1.1243E+00  2.9714E+00  2.0002E-01 -1.1608E+01 -2.6169E+01 -3.8475E+01 -1.5919E+00
             2.7065E+00
 GRADIENT:   3.5849E+01  0.0000E+00  3.8228E+01  1.8288E+01  9.3997E-03  1.4510E+00  0.0000E+00  0.0000E+00  0.0000E+00  2.0023E-06
             1.9344E+01

0ITERATION NO.:   53    OBJECTIVE VALUE:  -607.407289623464        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1136
 NPARAMETR:  4.5619E-01  1.0000E-02  2.5147E-02  2.9385E-01  1.6072E+01  1.1033E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.8404E-01
             1.3568E+01
 PARAMETER: -6.8485E-01 -1.7800E+01 -3.5830E+00 -1.1247E+00  2.8771E+00  1.9834E-01 -1.1608E+01 -2.6169E+01 -3.8475E+01 -1.5926E+00
             2.7077E+00
 GRADIENT:  -6.5278E-01  0.0000E+00  2.0220E-01 -1.8845E-01  4.2679E-04  3.8128E-02  0.0000E+00  0.0000E+00  0.0000E+00  5.2472E-07
            -8.1656E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1136
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1636E-03  1.9081E-06  8.9217E-05 -1.9291E-04  7.3152E-06
 SE:             2.7524E-02  1.6401E-06  2.2913E-04  2.8031E-04  1.4363E-05
 N:                     100         100         100         100         100

 P VAL.:         9.0849E-01  2.4466E-01  6.9701E-01  4.9132E-01  6.1054E-01

 ETASHRINKSD(%)  7.7919E+00  9.9995E+01  9.9232E+01  9.9061E+01  9.9952E+01
 ETASHRINKVR(%)  1.4977E+01  1.0000E+02  9.9994E+01  9.9991E+01  1.0000E+02
 EBVSHRINKSD(%)  8.0834E+00  9.9994E+01  9.9110E+01  9.8888E+01  9.9952E+01
 EBVSHRINKVR(%)  1.5513E+01  1.0000E+02  9.9992E+01  9.9988E+01  1.0000E+02
 RELATIVEINF(%)  7.6161E-02  2.7086E-09  1.7875E-05  1.7312E-05  5.9802E-08
 EPSSHRINKSD(%)  5.0789E+00
 EPSSHRINKVR(%)  9.8998E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -607.40728962346361     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       127.74353694027457     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -607.407       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.56E-01  1.00E-02  2.51E-02  2.94E-01  1.61E+01  1.10E+00  1.00E-02  1.00E-02  1.00E-02  1.84E-01  1.36E+01
 


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
+        3.61E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.72E+04  0.00E+00  9.31E+05
 
 TH 4
+       -8.71E+02  0.00E+00 -9.13E+04  1.10E+04
 
 TH 5
+        2.86E-01  0.00E+00 -3.32E+00  1.55E-01  4.32E-05
 
 TH 6
+       -2.18E+01  0.00E+00  4.53E+02 -6.26E+01  2.27E-03  1.10E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        3.10E-04  0.00E+00  4.69E-03 -1.12E-04  1.67E-06 -5.83E-04  0.00E+00  0.00E+00  0.00E+00 -1.16E-03
 
 TH11
+       -3.64E+01  0.00E+00  4.03E+02 -1.82E+01 -3.42E-03  1.35E+00  0.00E+00  0.00E+00  0.00E+00 -2.39E-07  2.21E+00
 
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
 #CPUT: Total CPU Time in Seconds,       20.550
Stop Time:
Wed Sep 29 20:21:59 CDT 2021
