Wed Sep 29 03:46:06 CDT 2021
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
$DATA ../../../../data/int/SL2/dat98.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2141.38951118682        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8616E+02  6.7206E+01  4.3027E+01  1.7993E+02  1.7833E+02  5.8699E+01 -6.3116E+01 -1.6200E+02 -8.7943E+01 -3.5688E+01
            -3.0857E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3042.56317419165        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.7789E-01  1.1986E+00  1.1349E+00  8.6613E-01  1.0849E+00  8.2748E-01  1.1246E+00  9.5771E-01  1.1522E+00  1.0026E+00
             2.0278E+00
 PARAMETER:  7.7641E-02  2.8117E-01  2.2654E-01 -4.3716E-02  1.8150E-01 -8.9370E-02  2.1744E-01  5.6787E-02  2.4167E-01  1.0262E-01
             8.0697E-01
 GRADIENT:   9.3170E+01  5.6177E+01 -7.2989E+00 -9.8714E-01 -3.4955E+00 -4.7069E+01  1.7787E+01  1.2359E+00  3.4228E+00 -9.3261E+00
             6.9781E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3052.67751354054        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      219
 NPARAMETR:  9.7316E-01  1.4993E+00  1.5592E+00  7.6238E-01  1.4110E+00  9.8666E-01  8.3348E-01  8.6612E-01  1.3129E+00  1.4152E+00
             2.0428E+00
 PARAMETER:  7.2792E-02  5.0500E-01  5.4420E-01 -1.7131E-01  4.4426E-01  8.6573E-02 -8.2142E-02 -4.3730E-02  3.7225E-01  4.4725E-01
             8.1432E-01
 GRADIENT:  -2.0999E+01  4.7199E+01  2.4233E+00  5.7106E+01  5.7910E+00  1.4739E+01  1.3279E+00 -3.4591E+00  1.8006E+00  1.4140E+01
             7.5055E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3059.02512613227        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  9.7582E-01  1.5343E+00  1.6098E+00  7.0774E-01  1.4423E+00  9.6846E-01  7.8537E-01  2.1843E+00  1.4524E+00  1.3207E+00
             1.9900E+00
 PARAMETER:  7.5518E-02  5.2810E-01  5.7613E-01 -2.4568E-01  4.6625E-01  6.7948E-02 -1.4160E-01  8.8129E-01  4.7323E-01  3.7819E-01
             7.8816E-01
 GRADIENT:  -1.3444E+01  6.0242E+00 -1.5409E+01  2.2731E+01  1.8567E+00  8.0064E+00  1.8614E+00  1.5452E+00  1.5246E+01 -3.7660E+00
            -1.4926E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3064.04671210002        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      575
 NPARAMETR:  9.7769E-01  1.4883E+00  2.5102E+00  7.4878E-01  1.4906E+00  9.7500E-01  8.3269E-01  3.1867E+00  1.2835E+00  1.3467E+00
             2.0031E+00
 PARAMETER:  7.7439E-02  4.9766E-01  1.0204E+00 -1.8930E-01  4.9921E-01  7.4682E-02 -8.3088E-02  1.2590E+00  3.4959E-01  3.9765E-01
             7.9470E-01
 GRADIENT:  -9.1553E+00  2.7236E+01 -3.3347E+00  7.9889E+00 -1.8812E+01  1.0742E+01  4.7613E+00  5.5642E+00  1.0664E+01 -3.0396E+00
             1.8253E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3066.74042647465        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  9.7961E-01  1.3111E+00  3.3174E+00  8.7573E-01  1.4555E+00  9.6053E-01  8.7453E-01  3.2684E+00  1.1072E+00  1.3101E+00
             1.9790E+00
 PARAMETER:  7.9403E-02  3.7087E-01  1.2992E+00 -3.2699E-02  4.7532E-01  5.9731E-02 -3.4069E-02  1.2843E+00  2.0180E-01  3.7012E-01
             7.8258E-01
 GRADIENT:  -4.2545E+00  3.0851E+01  4.7508E+00  1.4062E+01 -1.8324E+01  5.3409E+00 -5.4981E-01  2.5855E-01  5.1262E+00 -5.8823E-01
             3.0346E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3067.00920819808        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      950
 NPARAMETR:  9.8598E-01  1.3100E+00  3.2624E+00  8.6744E-01  1.4582E+00  9.5201E-01  8.7387E-01  3.2501E+00  1.0604E+00  1.3088E+00
             1.9790E+00
 PARAMETER:  8.5876E-02  3.7003E-01  1.2825E+00 -4.2209E-02  4.7718E-01  5.0823E-02 -3.4827E-02  1.2787E+00  1.5861E-01  3.6910E-01
             7.8260E-01
 GRADIENT:   1.1189E+01  2.1209E+01  3.7791E+00  1.5530E+00 -1.5525E+01  2.0760E+00 -3.0509E+00  6.2034E-02 -2.4253E+00 -1.1672E+00
             1.8822E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3067.33070032676        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1126            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8149E-01  1.2950E+00  3.1563E+00  8.6606E-01  1.4708E+00  9.4698E-01  8.9876E-01  3.1970E+00  1.0766E+00  1.3173E+00
             1.9767E+00
 PARAMETER:  8.1313E-02  3.5849E-01  1.2494E+00 -4.3799E-02  4.8577E-01  4.5523E-02 -6.7402E-03  1.2622E+00  1.7385E-01  3.7556E-01
             7.8144E-01
 GRADIENT:   1.1273E+02  9.0774E+01  9.0032E+00  1.1526E+01  4.6715E+01  1.1904E+01  1.7388E+00  8.0892E+00  5.4677E+00  6.5811E+00
             1.5092E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3067.65317705283        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1311             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8128E-01  1.2508E+00  3.3071E+00  9.0465E-01  1.4675E+00  9.4664E-01  9.6906E-01  3.1742E+00  1.0110E+00  1.2911E+00
             1.9787E+00
 PARAMETER:  8.1105E-02  3.2382E-01  1.2961E+00 -2.0784E-04  4.8358E-01  4.5168E-02  6.8569E-02  1.2551E+00  1.1098E-01  3.5553E-01
             7.8244E-01
 GRADIENT:   1.1174E+02  8.3386E+01  8.9482E+00  2.2903E+01  5.2006E+01  1.1744E+01  3.2738E+00  7.7000E+00  2.3508E+00  4.7161E+00
             1.8042E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3068.12818432084        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1487
 NPARAMETR:  9.8189E-01  1.2508E+00  3.3077E+00  9.0281E-01  1.4671E+00  9.4718E-01  9.2232E-01  3.1744E+00  1.0376E+00  1.2911E+00
             1.9783E+00
 PARAMETER:  8.1725E-02  3.2379E-01  1.2963E+00 -2.2446E-03  4.8327E-01  4.5735E-02  1.9133E-02  1.2551E+00  1.3693E-01  3.5550E-01
             7.8225E-01
 GRADIENT:   1.3876E+00  6.4107E+00 -1.9269E+00 -5.7632E-02  4.3917E+00  1.7077E-01  6.9593E-03 -2.7649E+00  1.0554E-01 -2.1351E+00
             6.2883E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3068.14471739877        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1644
 NPARAMETR:  9.8138E-01  1.2444E+00  3.3011E+00  9.0459E-01  1.4655E+00  9.4676E-01  9.2240E-01  3.1665E+00  1.0361E+00  1.3005E+00
             1.9811E+00
 PARAMETER:  8.1203E-02  3.1862E-01  1.2943E+00 -2.6953E-04  4.8221E-01  4.5291E-02  1.9224E-02  1.2526E+00  1.3544E-01  3.6274E-01
             7.8366E-01
 GRADIENT:   1.0915E-01  2.3079E+00 -2.2985E+00 -2.4789E+00  4.8204E+00  2.9278E-02 -6.6252E-02 -2.8364E+00  2.0563E-01 -2.8501E-01
             1.0136E+01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -3068.14471739877        NO. OF FUNC. EVALS.:  26
 CUMULATIVE NO. OF FUNC. EVALS.:     1670
 NPARAMETR:  9.8147E-01  1.2459E+00  3.2887E+00  9.0425E-01  1.4628E+00  9.4677E-01  9.2325E-01  3.1517E+00  1.0347E+00  1.2987E+00
             1.9866E+00
 PARAMETER:  8.1203E-02  3.1862E-01  1.2943E+00 -2.6953E-04  4.8221E-01  4.5291E-02  1.9224E-02  1.2526E+00  1.3544E-01  3.6274E-01
             7.8366E-01
 GRADIENT:  -3.0786E-01 -5.7412E+02  1.1857E+02  1.8380E+03  1.9160E+02 -5.3902E-03 -4.7979E-02  1.4564E+02  1.4529E-01  5.0296E+02
            -2.2915E+02
 NUMSIGDIG:         2.9         2.3         2.4         2.3         2.3         3.9         1.9         2.3         1.9         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1670
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3016E-03 -3.1495E-02 -3.0866E-02  1.8452E-02 -3.7742E-02
 SE:             2.9637E-02  1.9359E-02  1.9030E-02  2.3920E-02  2.3702E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6497E-01  1.0376E-01  1.0481E-01  4.4046E-01  1.1131E-01

 ETASHRINKSD(%)  7.1148E-01  3.5144E+01  3.6248E+01  1.9865E+01  2.0594E+01
 ETASHRINKVR(%)  1.4179E+00  5.7937E+01  5.9356E+01  3.5784E+01  3.6947E+01
 EBVSHRINKSD(%)  1.0298E+00  3.5543E+01  3.9818E+01  2.2038E+01  1.7132E+01
 EBVSHRINKVR(%)  2.0491E+00  5.8453E+01  6.3782E+01  3.9219E+01  3.1329E+01
 RELATIVEINF(%)  9.7923E+01  6.3563E+00  1.6510E+01  9.6238E+00  4.1406E+01
 EPSSHRINKSD(%)  1.8805E+01
 EPSSHRINKVR(%)  3.4073E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3068.1447173987726     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1414.0553576303619     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    52.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3068.145       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.24E+00  3.30E+00  9.05E-01  1.47E+00  9.47E-01  9.22E-01  3.17E+00  1.04E+00  1.30E+00  1.98E+00
 


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
+        1.26E+03
 
 TH 2
+       -4.54E+01  2.96E+04
 
 TH 3
+        2.90E+00 -3.23E+02  2.76E+02
 
 TH 4
+        3.08E+02  3.69E+02 -8.65E+00  5.63E+05
 
 TH 5
+        1.33E+01 -8.76E+01  1.42E+02  1.08E+02  1.04E+04
 
 TH 6
+        3.34E+00 -4.09E+00  1.75E-01  9.87E+00 -4.93E-01  2.14E+02
 
 TH 7
+        1.49E+00 -1.44E+01  4.24E+00 -5.51E+05 -1.68E+01 -5.82E-02  5.41E+05
 
 TH 8
+        3.25E+00  2.12E-01 -2.87E+01  8.59E+00  4.36E+01  2.26E-01  2.09E+00  3.05E+02
 
 TH 9
+        1.46E+00 -1.41E+01 -4.14E+00  3.62E+05  1.72E+01 -6.63E-01 -3.56E+05  1.23E+00  7.70E+01
 
 TH10
+        2.70E+01  4.29E+00 -2.66E+01  9.30E+01 -1.67E+02  1.45E+00  1.64E+00 -3.20E+01  8.15E+00  2.05E+04
 
 TH11
+       -2.10E+01 -2.09E+01 -2.36E+01 -3.50E+01  9.23E+01  1.69E+00  7.19E+00 -1.81E+01  2.57E+00  7.14E+01  2.27E+03
 
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
 #CPUT: Total CPU Time in Seconds,       68.769
Stop Time:
Wed Sep 29 03:47:16 CDT 2021
