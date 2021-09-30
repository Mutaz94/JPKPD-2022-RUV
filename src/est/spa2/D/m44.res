Thu Sep 30 09:05:19 CDT 2021
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
$DATA ../../../../data/spa2/D/dat44.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   18193.8779918458        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.4100E+02  4.1268E+02  3.8294E+01  4.3429E+02  6.3436E+01 -1.8459E+03 -1.0876E+03 -9.5068E+01 -1.3395E+03 -3.2148E+02
            -3.5957E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -532.845461843909        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7457E-01  1.2044E+00  9.4586E-01  1.3239E+00  1.0981E+00  2.1002E+00  1.7716E+00  1.0067E+00  1.4403E+00  1.0199E+00
             1.3910E+01
 PARAMETER:  7.4239E-02  2.8600E-01  4.4334E-02  3.8061E-01  1.9361E-01  8.4202E-01  6.7189E-01  1.0672E-01  4.6483E-01  1.1967E-01
             2.7326E+00
 GRADIENT:  -9.6372E+01 -3.2112E+01 -2.7187E+01  2.1246E+01  3.5320E+01  4.1333E+01 -1.0612E+01  2.6090E+00  1.0910E+01  1.2546E+01
             3.4969E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -600.693299865368        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0188E+00  7.3077E-01  1.0640E+01  2.0671E+00  2.4510E+00  2.1075E+00  6.6501E+00  9.9709E-01  1.4738E+00  7.2856E-01
             1.2648E+01
 PARAMETER:  1.1866E-01 -2.1366E-01  2.4646E+00  8.2616E-01  9.9648E-01  8.4550E-01  1.9946E+00  9.7088E-02  4.8787E-01 -2.1668E-01
             2.6375E+00
 GRADIENT:  -5.2794E+01  2.0067E+01 -2.2319E+00  5.8811E+01  6.7475E+00  4.0079E+01  5.2269E+01  1.1685E-01  3.2775E+00  2.2964E+00
             3.0192E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -642.194866239511        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0764E+00  7.1085E-01  1.3843E+01  1.3683E+00  1.9551E+00  2.1984E+00  3.4003E+00  8.1240E+00  1.5740E+00  1.0197E+00
             9.9346E+00
 PARAMETER:  1.7362E-01 -2.4130E-01  2.7278E+00  4.1357E-01  7.7043E-01  8.8772E-01  1.3239E+00  2.1948E+00  5.5363E-01  1.1953E-01
             2.3960E+00
 GRADIENT:   8.5816E-01 -3.4705E+01 -3.5557E+00 -2.5055E+01  3.1929E+00  2.2650E+01 -3.9967E+01  7.1751E+00  8.5315E+00  1.0530E+01
             1.8009E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -671.302489074725        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0174E+00  1.2630E+00  5.1369E+00  9.4433E-01  1.7400E+00  2.1355E+00  3.4482E+00  4.0019E+00  8.5701E-01  2.5807E-01
             8.4229E+00
 PARAMETER:  1.1721E-01  3.3352E-01  1.7365E+00  4.2717E-02  6.5386E-01  8.5870E-01  1.3378E+00  1.4868E+00 -5.4303E-02 -1.2545E+00
             2.2309E+00
 GRADIENT:  -6.1447E+00  2.0635E+00  2.0970E+00  6.2524E+00 -9.6028E+00 -4.4069E-01  3.3562E+00 -2.5763E+00  4.5972E+00  8.9709E-01
             4.4242E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -673.083669963357        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0451E+00  1.5051E+00  5.4423E+00  7.5092E-01  1.8928E+00  2.1408E+00  3.0255E+00  4.9314E+00  3.6783E-01  1.3262E-01
             8.3312E+00
 PARAMETER:  1.4412E-01  5.0884E-01  1.7942E+00 -1.8645E-01  7.3803E-01  8.6116E-01  1.2071E+00  1.6956E+00 -9.0013E-01 -1.9203E+00
             2.2200E+00
 GRADIENT:   9.5692E+00 -2.8731E+00  4.8308E-01 -7.4333E+00  1.8124E+00 -1.8635E-01 -3.2663E+00  1.0907E-01  1.4708E+00  2.2173E-01
            -4.0578E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -677.867300491990        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      509
 NPARAMETR:  1.0415E+00  1.3758E+00  5.6387E+00  8.7097E-01  1.9075E+00  2.3411E+00  3.6112E+00  5.0696E+00  1.4312E-01  5.6334E-02
             8.5379E+00
 PARAMETER:  1.4064E-01  4.1901E-01  1.8297E+00 -3.8145E-02  7.4579E-01  9.5061E-01  1.3840E+00  1.7233E+00 -1.8441E+00 -2.7765E+00
             2.2445E+00
 GRADIENT:  -5.1788E+00 -7.2476E-01 -1.6952E+00 -4.9043E-01  1.5298E+00 -6.8508E-01 -8.0348E-01  1.6717E+00  2.3309E-01  3.7092E-02
             8.0129E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -678.314694129314        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      686
 NPARAMETR:  1.0584E+00  1.3188E+00  7.3249E+00  8.9965E-01  1.9586E+00  2.3478E+00  3.7147E+00  5.0549E+00  1.0097E-01  4.1940E-02
             8.5398E+00
 PARAMETER:  1.5671E-01  3.7671E-01  2.0913E+00 -5.7468E-03  7.7222E-01  9.5346E-01  1.4123E+00  1.7204E+00 -2.1930E+00 -3.0715E+00
             2.2447E+00
 GRADIENT:   8.5954E-01 -1.3563E+00 -1.2045E+00  3.4529E+00  1.1651E+00  2.0314E-01 -2.4779E+00  7.7710E-01  1.0289E-01  1.8629E-02
            -3.1077E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -679.531604649754        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      864
 NPARAMETR:  1.0530E+00  1.3750E+00  1.2459E+02  8.8268E-01  2.0165E+00  2.3263E+00  3.6271E+00  1.7296E+01  1.0000E-02  1.0000E-02
             8.5232E+00
 PARAMETER:  1.5161E-01  4.1848E-01  4.9250E+00 -2.4790E-02  8.0135E-01  9.4429E-01  1.3884E+00  2.9505E+00 -6.7108E+00 -6.5599E+00
             2.2428E+00
 GRADIENT:  -1.1020E+00  1.2201E-01  4.9858E-01  6.4626E-01 -1.4031E+00 -3.4035E+00 -2.2793E+00 -4.2658E-01  0.0000E+00  0.0000E+00
            -2.4464E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -679.643455909757        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1042
 NPARAMETR:  1.0595E+00  1.3505E+00  5.8369E+01  8.9547E-01  2.0360E+00  2.3623E+00  3.7218E+00  1.3500E+01  1.0000E-02  1.0000E-02
             8.5433E+00
 PARAMETER:  1.5778E-01  4.0048E-01  4.1668E+00 -1.0407E-02  8.1097E-01  9.5965E-01  1.4142E+00  2.7027E+00 -5.5310E+00 -5.6271E+00
             2.2451E+00
 GRADIENT:   1.0264E+00 -2.0780E-01  5.8695E-01 -1.5373E+00  1.4525E+00  1.8025E+00  1.0813E+00 -1.5258E+00  0.0000E+00  0.0000E+00
             1.7390E+00

0ITERATION NO.:   47    OBJECTIVE VALUE:  -679.655914033787        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1108
 NPARAMETR:  1.0578E+00  1.3517E+00  5.6643E+01  8.9592E-01  2.0388E+00  2.3617E+00  3.7144E+00  1.3372E+01  1.0000E-02  1.0000E-02
             8.5678E+00
 PARAMETER:  1.5774E-01  4.0094E-01  4.1294E+00 -1.0748E-02  8.1089E-01  9.5944E-01  1.4137E+00  2.6979E+00 -5.4757E+00 -5.5821E+00
             2.2451E+00
 GRADIENT:   8.0231E-01 -1.6204E-01 -2.2797E+01 -8.2480E-01 -1.1701E+02  4.7333E-02  1.0155E+00  3.4355E+01  0.0000E+00  0.0000E+00
            -4.1175E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1108
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6397E-03  7.7961E-03 -5.5128E-03 -7.3876E-04  1.2661E-05
 SE:             2.8696E-02  2.6303E-02  4.1673E-03  1.4944E-04  6.5317E-05
 N:                     100         100         100         100         100

 P VAL.:         7.3693E-01  7.6693E-01  1.8588E-01  7.6882E-07  8.4630E-01

 ETASHRINKSD(%)  3.8632E+00  1.1881E+01  8.6039E+01  9.9499E+01  9.9781E+01
 ETASHRINKVR(%)  7.5772E+00  2.2350E+01  9.8051E+01  9.9997E+01  1.0000E+02
 EBVSHRINKSD(%)  2.6940E+00  8.4859E+00  8.9524E+01  9.9602E+01  9.9681E+01
 EBVSHRINKVR(%)  5.3154E+00  1.6252E+01  9.8902E+01  9.9998E+01  9.9999E+01
 RELATIVEINF(%)  9.4224E+01  3.8186E+01  9.9101E-01  7.2418E-04  8.4100E-04
 EPSSHRINKSD(%)  8.9321E+00
 EPSSHRINKVR(%)  1.7066E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -679.65591403378744     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       423.07032581181966     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.54
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -679.656       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.35E+00  5.62E+01  8.95E-01  2.04E+00  2.36E+00  3.72E+00  1.34E+01  1.00E-02  1.00E-02  8.54E+00
 


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
+        1.66E+02
 
 TH 2
+       -2.58E+00  3.73E+01
 
 TH 3
+        1.44E-02  5.52E-02  4.43E-02
 
 TH 4
+       -5.11E+00  6.43E+01  1.65E-01  4.80E+02
 
 TH 5
+       -1.75E-01 -3.06E-01 -1.98E-01  5.55E-03  8.92E+02
 
 TH 6
+        1.73E+00 -5.61E-01 -6.54E-03  1.70E+00 -2.17E+00  2.98E+01
 
 TH 7
+        2.11E-01  3.04E+00 -4.60E-02 -5.14E+03 -4.74E+00 -8.18E-01  9.52E+01
 
 TH 8
+       -8.48E-02 -3.69E-01  4.57E-03 -9.62E-01  1.23E+00  3.35E-02  3.07E-01  1.87E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.52E+00 -3.30E+00  5.36E-01 -1.74E+01 -1.22E+00  1.49E+00  1.10E+00  6.07E-02  0.00E+00  0.00E+00  1.60E+01
 
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
 #CPUT: Total CPU Time in Seconds,       29.656
Stop Time:
Thu Sep 30 09:05:51 CDT 2021
