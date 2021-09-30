Wed Sep 29 20:50:58 CDT 2021
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
$DATA ../../../../data/spa1/B/dat9.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2076.82606731277        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0818E+02 -2.3365E+01 -1.6307E+01  1.1792E+00  6.9841E+01  1.9437E+01  5.7780E+00  1.9543E+00  1.1410E+01 -5.3430E+00
             4.5387E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2083.11736182033        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      113
 NPARAMETR:  9.3033E-01  9.8620E-01  9.8218E-01  1.0075E+00  9.7405E-01  9.4723E-01  9.8826E-01  9.9697E-01  9.8661E-01  9.9900E-01
             9.6513E-01
 PARAMETER:  2.7783E-02  8.6100E-02  8.2016E-02  1.0746E-01  7.3707E-02  4.5789E-02  8.8186E-02  9.6967E-02  8.6519E-02  9.8998E-02
             6.4508E-02
 GRADIENT:  -1.6215E+01 -4.9511E+01 -1.4990E+01 -6.5890E+01  4.9188E+01 -2.0426E+01  2.4200E+00  3.0006E+00  3.1064E+00 -3.6029E+00
             2.0731E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2087.98951987587        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      291
 NPARAMETR:  9.3637E-01  8.0323E-01  9.2212E-01  1.1620E+00  8.0228E-01  1.0348E+00  9.4975E-01  8.0306E-01  9.0147E-01  9.6018E-01
             9.3707E-01
 PARAMETER:  3.4256E-02 -1.1911E-01  1.8915E-02  2.5015E-01 -1.2030E-01  1.3416E-01  4.8449E-02 -1.1932E-01 -3.7296E-03  5.9362E-02
             3.4998E-02
 GRADIENT:   2.0592E+00  2.5530E+01  2.4723E+01  6.7635E+00 -4.5750E+01  1.5385E+01 -2.5583E+00 -9.3068E-01 -4.1851E+00  3.1015E+00
             4.5920E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2089.66556376615        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  9.3388E-01  7.2001E-01  9.3509E-01  1.2149E+00  7.9663E-01  9.9219E-01  1.2132E+00  7.8201E-01  8.4895E-01  9.0856E-01
             9.3104E-01
 PARAMETER:  3.1591E-02 -2.2850E-01  3.2890E-02  2.9462E-01 -1.2737E-01  9.2157E-02  2.9323E-01 -1.4589E-01 -6.3758E-02  4.1100E-03
             2.8550E-02
 GRADIENT:  -1.4403E+00  1.2444E+01  7.1692E+00  1.4495E+01 -7.3771E+00 -2.7610E-01 -1.0112E+00 -2.1159E+00 -2.9218E+00 -2.4243E+00
            -9.6005E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2090.79065809079        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  9.2926E-01  4.0302E-01  1.1330E+00  1.4145E+00  7.9090E-01  9.8884E-01  1.4234E+00  9.7902E-01  7.9349E-01  9.6477E-01
             9.2991E-01
 PARAMETER:  2.6629E-02 -8.0876E-01  2.2491E-01  4.4681E-01 -1.3459E-01  8.8782E-02  4.5303E-01  7.8801E-02 -1.3132E-01  6.4134E-02
             2.7328E-02
 GRADIENT:   2.9191E-03  4.7466E+00  3.4285E+00  1.1588E+01 -7.1815E+00  5.1089E-01  3.4472E-01 -9.1894E-02  1.0620E+00  1.0034E+00
            -6.5433E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2090.90876797489        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      821
 NPARAMETR:  9.2728E-01  3.0934E-01  1.1958E+00  1.4744E+00  7.9266E-01  9.8651E-01  1.4685E+00  1.0470E+00  7.6969E-01  9.6985E-01
             9.3101E-01
 PARAMETER:  2.4501E-02 -1.0733E+00  2.7879E-01  4.8823E-01 -1.3236E-01  8.6423E-02  4.8424E-01  1.4591E-01 -1.6176E-01  6.9389E-02
             2.8516E-02
 GRADIENT:  -7.5117E-01  3.3886E+00  2.6179E+00  1.1854E+01 -4.7685E+00  1.7778E-01  9.8113E-02 -3.0138E-01 -6.3845E-01  2.6105E-01
            -1.6515E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2090.96399913007        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  9.2710E-01  2.5534E-01  1.2186E+00  1.4963E+00  7.9304E-01  9.8551E-01  1.4248E+00  1.0819E+00  7.6150E-01  9.7064E-01
             9.3133E-01
 PARAMETER:  2.4301E-02 -1.2652E+00  2.9774E-01  5.0300E-01 -1.3189E-01  8.5406E-02  4.5400E-01  1.7876E-01 -1.7246E-01  7.0199E-02
             2.8857E-02
 GRADIENT:   1.3036E+00 -9.4158E-01 -2.2852E+00 -1.5469E+01  5.2176E+00  1.5890E-01  1.0798E-01  8.5586E-02  2.3712E-01 -3.8382E-01
             2.7505E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2090.99789341533        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  9.2588E-01  2.3423E-01  1.2282E+00  1.5153E+00  7.8895E-01  9.8457E-01  6.5491E-01  1.0958E+00  7.6820E-01  9.7838E-01
             9.3114E-01
 PARAMETER:  2.2991E-02 -1.3514E+00  3.0552E-01  5.1561E-01 -1.3705E-01  8.4446E-02 -3.2326E-01  1.9150E-01 -1.6370E-01  7.8140E-02
             2.8659E-02
 GRADIENT:  -8.4748E-01  8.5267E-01 -8.8830E-01  1.8253E+00 -1.2585E+00 -8.7177E-02  5.5496E-02 -1.8112E-01  1.6377E+00  3.1534E-01
            -1.4981E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2091.02701452029        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1363
 NPARAMETR:  9.2703E-01  2.2568E-01  1.2318E+00  1.5132E+00  7.8956E-01  9.8487E-01  1.5848E-01  1.1012E+00  7.6195E-01  9.7467E-01
             9.3131E-01
 PARAMETER:  2.4225E-02 -1.3886E+00  3.0846E-01  5.1426E-01 -1.3628E-01  8.4751E-02 -1.7422E+00  1.9641E-01 -1.7188E-01  7.4342E-02
             2.8837E-02
 GRADIENT:   2.3779E+00 -7.2829E-01 -1.6064E+00 -1.5326E+01  2.1283E+00  1.2210E-01  5.0803E-03 -3.3308E-01 -1.2109E+00 -6.6147E-01
            -1.5058E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2091.04341148694        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1542
 NPARAMETR:  9.2639E-01  2.2724E-01  1.2350E+00  1.5159E+00  7.8940E-01  9.8468E-01  5.1497E-02  1.1062E+00  7.6513E-01  9.7837E-01
             9.3144E-01
 PARAMETER:  2.3542E-02 -1.3817E+00  3.1106E-01  5.1598E-01 -1.3649E-01  8.4564E-02 -2.8662E+00  2.0095E-01 -1.6771E-01  7.8129E-02
             2.8981E-02
 GRADIENT:   7.1133E-01  1.6443E-01 -2.0364E-01 -7.0446E+00 -1.6720E+00  2.3051E-02  8.3594E-04 -2.3870E-02  1.0893E-01  1.3360E-01
             3.7586E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2091.04787712640        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1726
 NPARAMETR:  9.2680E-01  2.2902E-01  1.2393E+00  1.5141E+00  7.9088E-01  9.8489E-01  1.0000E-02  1.1084E+00  7.6532E-01  9.8019E-01
             9.3151E-01
 PARAMETER:  2.3980E-02 -1.3740E+00  3.1452E-01  5.1481E-01 -1.3461E-01  8.4770E-02 -4.6487E+00  2.0293E-01 -1.6746E-01  7.9996E-02
             2.9049E-02
 GRADIENT:   1.6235E+00  3.8243E-02  8.6165E-01 -9.5454E+00 -2.6397E+00  9.1389E-02  0.0000E+00 -6.2449E-02 -2.5651E-02  2.3350E-01
             7.2141E-02

0ITERATION NO.:   52    OBJECTIVE VALUE:  -2091.04894656841        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1783
 NPARAMETR:  9.2659E-01  2.2895E-01  1.2387E+00  1.5153E+00  7.9136E-01  9.8481E-01  1.0000E-02  1.1086E+00  7.6534E-01  9.7921E-01
             9.3148E-01
 PARAMETER:  2.3756E-02 -1.3742E+00  3.1409E-01  5.1562E-01 -1.3401E-01  8.4697E-02 -4.6487E+00  2.0314E-01 -1.6744E-01  7.8996E-02
             2.9020E-02
 GRADIENT:   1.1102E+00  2.1235E-01 -9.2446E-03 -6.7099E+00 -1.4627E+00  5.7349E-02  0.0000E+00 -5.8181E-02 -1.9500E-02  3.4420E-02
             1.2279E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1783
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6320E-04 -7.7000E-05 -3.0676E-02 -5.4441E-03 -3.0122E-02
 SE:             2.9872E-02  4.3129E-05  1.8261E-02  2.9506E-02  2.2074E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8496E-01  7.4208E-02  9.2989E-02  8.5362E-01  1.7237E-01

 ETASHRINKSD(%)  1.0000E-10  9.9856E+01  3.8823E+01  1.1512E+00  2.6050E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.2574E+01  2.2891E+00  4.5314E+01
 EBVSHRINKSD(%)  3.0783E-01  9.9865E+01  4.1694E+01  1.6610E+00  2.3675E+01
 EBVSHRINKVR(%)  6.1472E-01  1.0000E+02  6.6004E+01  3.2945E+00  4.1745E+01
 RELATIVEINF(%)  9.7627E+01  1.4372E-05  8.1439E+00  9.1420E+00  1.1632E+01
 EPSSHRINKSD(%)  3.4834E+01
 EPSSHRINKVR(%)  5.7533E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2091.0489465684136     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1172.1104133637409     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2091.049       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.27E-01  2.29E-01  1.24E+00  1.52E+00  7.91E-01  9.85E-01  1.00E-02  1.11E+00  7.65E-01  9.79E-01  9.31E-01
 


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
+        1.33E+03
 
 TH 2
+       -3.07E+01  4.07E+02
 
 TH 3
+       -2.13E+00  7.78E+01  2.06E+02
 
 TH 4
+       -5.20E+00  4.84E+02 -4.21E+01  7.99E+02
 
 TH 5
+        5.33E+00 -3.09E+02 -3.87E+02 -6.37E+01  1.09E+03
 
 TH 6
+        1.56E+00 -4.24E+00 -2.68E-01 -1.95E+00 -1.19E+00  2.02E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -4.30E-01  3.30E+00 -3.57E+01 -5.06E+00 -2.32E-01 -5.08E-02  0.00E+00  3.87E+01
 
 TH 9
+        2.16E+00 -1.02E+02  6.91E+00 -1.39E+00  6.58E-01 -6.10E-01  0.00E+00 -1.20E+00  3.20E+02
 
 TH10
+        1.71E+00  9.22E+00 -6.16E+00  3.18E+00 -7.90E+01  3.50E-01  0.00E+00  1.97E+01  3.63E+00  7.62E+01
 
 TH11
+       -8.79E+00 -9.91E+00 -1.06E+01 -1.11E+01 -5.77E+00  2.26E+00  0.00E+00  8.28E+00  6.99E+00  1.37E+01  4.55E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.723
Stop Time:
Wed Sep 29 20:51:32 CDT 2021
