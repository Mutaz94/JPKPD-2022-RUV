Sat Sep 25 05:41:19 CDT 2021
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
$DATA ../../../../data/int/D/dat35.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   28307.6351507785        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.8666E+01  3.2593E+02 -6.8495E+01  2.0128E+02  1.7321E+02 -1.5183E+03 -8.2783E+02 -5.9452E+01 -1.5274E+03 -5.8068E+02
            -5.9713E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -976.158359981394        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.8130E+00  1.9501E+00  9.4747E-01  2.0861E+00  8.5020E-01  4.0894E+00  4.8968E+00  9.9896E-01  2.9319E+00  1.7588E+00
             1.3031E+01
 PARAMETER:  6.9501E-01  7.6786E-01  4.6035E-02  8.3530E-01 -6.2285E-02  1.5084E+00  1.6886E+00  9.8959E-02  1.1757E+00  6.6461E-01
             2.6673E+00
 GRADIENT:   1.9249E+01  2.4296E+01 -4.1561E+01  8.3409E+01 -2.6349E+01  1.3321E+02  5.9878E+01  4.2047E+00  5.1118E+01  3.7260E+01
             4.6911E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1063.19714262119        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.6329E+00  2.6184E+00  4.2214E+01  3.5728E+00  2.2456E+00  2.1948E+00  1.9946E+01  8.0014E-01  2.4606E+00  1.5739E+00
             1.3074E+01
 PARAMETER:  5.9037E-01  1.0626E+00  3.8428E+00  1.3733E+00  9.0898E-01  8.8609E-01  3.0931E+00 -1.2297E-01  1.0004E+00  5.5355E-01
             2.6706E+00
 GRADIENT:   3.7210E+01  2.0488E+01 -5.8258E+00  5.0673E+01  1.1053E+01  3.9039E+01  3.1493E+01 -1.8855E-03  9.0092E+00  2.8519E+01
             4.9422E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1268.49225243997        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0211E+00  7.0343E-01  3.4377E+01  1.6402E+00  2.3279E+00  2.0040E+00  2.9003E+00  1.5128E+00  3.9894E+00  6.0718E-01
             8.6608E+00
 PARAMETER:  1.2092E-01 -2.5179E-01  3.6374E+00  5.9483E-01  9.4499E-01  7.9515E-01  1.1648E+00  5.1396E-01  1.4837E+00 -3.9893E-01
             2.2588E+00
 GRADIENT:  -1.1545E+02 -2.4659E+01 -9.5472E+00 -7.6847E+00  6.1594E+01 -2.3519E+01  1.1228E+01  1.0763E+00  7.2962E+01  4.3400E+00
             1.3036E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1313.53808378451        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.3194E+00  8.0685E-01  5.2296E+01  1.5310E+00  2.1974E+00  2.0185E+00  5.1856E+00  6.2047E-01  1.7425E+00  5.4938E-01
             8.1122E+00
 PARAMETER:  3.7716E-01 -1.1462E-01  4.0569E+00  5.2591E-01  8.8729E-01  8.0236E-01  1.7459E+00 -3.7728E-01  6.5530E-01 -4.9896E-01
             2.1934E+00
 GRADIENT:   2.1652E+01 -9.1258E-01 -1.3448E-01  4.7186E+00  2.9400E+00 -3.6408E+00  1.2278E+00 -1.1334E-03  6.0165E+00  1.5067E+00
            -3.5875E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1315.33359575481        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.2652E+00  1.4710E+00  2.5118E+01  1.0294E+00  2.1991E+00  2.0057E+00  4.1988E+00  7.7637E-01  1.0041E+00  3.0062E-01
             8.3027E+00
 PARAMETER:  3.3522E-01  4.8592E-01  3.3236E+00  1.2899E-01  8.8804E-01  7.9602E-01  1.5348E+00 -1.5312E-01  1.0408E-01 -1.1019E+00
             2.2166E+00
 GRADIENT:  -2.4406E+00  1.7745E+00  7.0114E-01  4.8880E+00 -1.7945E+00 -3.9181E+00 -1.8335E-02 -1.2989E-03 -1.9066E+00  2.0061E-01
             9.6334E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1315.66134543399        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      513
 NPARAMETR:  1.2652E+00  1.4740E+00  1.8086E+01  1.0088E+00  2.1992E+00  2.0059E+00  4.1952E+00  7.7731E-01  1.0139E+00  1.0386E-01
             8.3041E+00
 PARAMETER:  3.3519E-01  4.8799E-01  2.9951E+00  1.0872E-01  8.8811E-01  7.9609E-01  1.5339E+00 -1.5192E-01  1.1381E-01 -2.1647E+00
             2.2167E+00
 GRADIENT:  -2.2946E+00 -4.0894E-01  4.7721E-01 -5.1938E-01  5.2015E+00 -4.1058E+00  3.0729E+00 -3.7657E-03 -3.0590E-01  2.2045E-02
             1.0671E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1315.72584461651        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  1.2652E+00  1.4740E+00  1.5122E+01  1.0123E+00  2.1991E+00  2.0059E+00  4.1950E+00  7.7950E-01  1.0265E+00  1.8121E-02
             8.3025E+00
 PARAMETER:  3.3519E-01  4.8800E-01  2.8162E+00  1.1218E-01  8.8804E-01  7.9609E-01  1.5339E+00 -1.4910E-01  1.2611E-01 -3.9107E+00
             2.2166E+00
 GRADIENT:  -2.4741E+00 -4.0358E-01  1.4188E-02  1.5683E-01  1.1221E+01 -4.1190E+00  3.2089E+00 -2.8446E-03  4.2368E-03  9.0403E-04
             9.9579E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1315.73620635311        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  1.2650E+00  1.4744E+00  1.5409E+01  1.0198E+00  2.1967E+00  2.0056E+00  4.1925E+00  8.2659E-01  1.0689E+00  1.0000E-02
             8.2767E+00
 PARAMETER:  3.3507E-01  4.8822E-01  2.8349E+00  1.1964E-01  8.8696E-01  7.9596E-01  1.5333E+00 -9.0446E-02  1.6664E-01 -2.2555E+01
             2.2134E+00
 GRADIENT:  -2.0701E+00  3.2740E-01  3.7841E-02  2.4121E-01  1.0371E+01 -4.3244E+00  3.1835E+00 -4.1154E-03  4.6670E-01  0.0000E+00
             3.6437E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1315.74037442809        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  1.2649E+00  1.4745E+00  1.5176E+01  1.0171E+00  2.1954E+00  2.0055E+00  4.1912E+00  8.5340E-01  1.0625E+00  1.0000E-02
             8.2631E+00
 PARAMETER:  3.3500E-01  4.8834E-01  2.8197E+00  1.1696E-01  8.8638E-01  7.9589E-01  1.5330E+00 -5.8528E-02  1.6065E-01 -3.2094E+01
             2.2118E+00
 GRADIENT:  -1.8496E+00  2.1849E-01  2.6421E-02  1.9396E-01  1.0583E+01 -4.3933E+00  3.2470E+00 -4.7434E-03  4.2154E-01  0.0000E+00
             3.4668E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1315.74037442809        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      808
 NPARAMETR:  1.2649E+00  1.4746E+00  1.5175E+01  1.0171E+00  2.1955E+00  2.0054E+00  4.1914E+00  8.5348E-01  1.0625E+00  1.0000E-02
             8.2638E+00
 PARAMETER:  3.3500E-01  4.8834E-01  2.8197E+00  1.1696E-01  8.8638E-01  7.9589E-01  1.5330E+00 -5.8528E-02  1.6065E-01 -3.2094E+01
             2.2118E+00
 GRADIENT:   6.1090E+03 -2.0985E+03  3.6308E+02  1.7337E+04 -2.3011E+03  2.5730E+03 -6.7234E+02 -3.4377E-03  6.3115E+03  0.0000E+00
            -4.6447E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      808
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3802E-02  1.9674E-02 -3.7465E-04 -5.4825E-02 -5.2224E-05
 SE:             2.9575E-02  2.6262E-02  7.7549E-04  1.1203E-02  1.6014E-04
 N:                     100         100         100         100         100

 P VAL.:         6.4073E-01  4.5379E-01  6.2901E-01  9.9114E-07  7.4434E-01

 ETASHRINKSD(%)  9.1938E-01  1.2018E+01  9.7402E+01  6.2468E+01  9.9463E+01
 ETASHRINKVR(%)  1.8303E+00  2.2591E+01  9.9933E+01  8.5913E+01  9.9997E+01
 EBVSHRINKSD(%)  3.2599E+00  8.4878E+00  9.7197E+01  6.7458E+01  9.9460E+01
 EBVSHRINKVR(%)  6.4136E+00  1.6255E+01  9.9921E+01  8.9410E+01  9.9997E+01
 RELATIVEINF(%)  9.3314E+01  3.3825E+01  1.4685E-02  4.3458E+00  5.4435E-04
 EPSSHRINKSD(%)  7.1233E+00
 EPSSHRINKVR(%)  1.3739E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1315.7403744280948     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       338.34898534031595     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.38
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1315.740       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.26E+00  1.47E+00  1.52E+01  1.02E+00  2.20E+00  2.01E+00  4.19E+00  8.53E-01  1.06E+00  1.00E-02  8.26E+00
 


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
+        2.85E+06
 
 TH 2
+       -2.41E+03  9.89E+05
 
 TH 3
+        4.06E+00 -2.34E+00  2.80E+02
 
 TH 4
+        2.53E+02 -1.44E+02  3.37E+00  3.58E+07
 
 TH 5
+       -7.05E+02  4.16E+02 -1.46E+00 -6.02E+01  1.35E+05
 
 TH 6
+       -8.15E+01 -3.26E+02  5.50E-01  2.50E+01 -6.34E+01  2.02E+05
 
 TH 7
+        6.35E+01  6.08E+01 -1.16E-01 -1.81E+01  7.54E+00 -9.43E+01  1.25E+04
 
 TH 8
+       -2.86E+01 -6.30E+00 -8.13E-02 -7.91E+00 -5.35E+00 -2.34E+00  1.88E-01  5.62E+00
 
 TH 9
+        1.79E+01  1.34E+01  6.24E-02  2.50E+07 -1.65E-01 -2.04E+00  3.12E+00 -3.48E+07  1.74E+07
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.05E+00  2.84E+01 -5.43E-02 -1.01E+01  6.66E+00  9.83E+01 -3.01E+01 -3.40E-02  1.48E+00  0.00E+00  1.54E+03
 
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
 #CPUT: Total CPU Time in Seconds,       39.352
Stop Time:
Sat Sep 25 05:42:00 CDT 2021
