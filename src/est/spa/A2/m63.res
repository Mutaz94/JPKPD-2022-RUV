Wed Sep 29 12:57:11 CDT 2021
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
$DATA ../../../../data/spa/A2/dat63.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1192.44380041210        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3500E+02 -5.0487E+01 -1.6827E+01 -7.7782E+01  7.1149E+01  1.3328E+01 -1.4622E+01  8.2918E+00 -6.3087E+01 -1.0672E+00
            -8.7460E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1465.53105577221        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1526E+00  1.0407E+00  1.1778E+00  1.0648E+00  1.0237E+00  1.1874E+00  9.8700E-01  8.9439E-01  1.1781E+00  7.6344E-01
             2.2337E+00
 PARAMETER:  2.4200E-01  1.3988E-01  2.6365E-01  1.6280E-01  1.2339E-01  2.7175E-01  8.6919E-02 -1.1616E-02  2.6389E-01 -1.6992E-01
             9.0367E-01
 GRADIENT:   3.5396E+02 -4.8736E+00  9.9735E+00 -1.6117E+01 -1.3520E+01  4.9101E+01  5.3403E+00  3.7592E+00  6.6796E+00  8.3244E+00
             8.3220E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1472.77046650935        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.1007E+00  1.1680E+00  6.9983E-01  9.7929E-01  8.6167E-01  1.1034E+00  1.2492E+00  3.8882E-01  1.0984E+00  5.2181E-01
             2.1372E+00
 PARAMETER:  1.9596E-01  2.5531E-01 -2.5691E-01  7.9077E-02 -4.8883E-02  1.9842E-01  3.2253E-01 -8.4464E-01  1.9386E-01 -5.5045E-01
             8.5950E-01
 GRADIENT:   2.5856E+02  2.0677E+01 -1.4611E+01  3.5831E+00  2.2448E+01  3.0678E+01  2.3994E+01  1.3869E+00  1.4187E+00  6.1080E+00
            -1.0713E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1477.10839458679        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0276E+00  1.0585E+00  7.3652E-01  1.0302E+00  8.2061E-01  1.0258E+00  1.1262E+00  3.2545E-01  1.0939E+00  4.7839E-01
             2.1216E+00
 PARAMETER:  1.2726E-01  1.5687E-01 -2.0581E-01  1.2976E-01 -9.7712E-02  1.2542E-01  2.1887E-01 -1.0225E+00  1.8979E-01 -6.3733E-01
             8.5219E-01
 GRADIENT:   7.1679E+01  5.3884E+00  9.4296E-01 -3.4863E+00 -9.9923E-01  3.5369E+00  1.9102E+00  7.5506E-01  7.4645E-01  3.0410E+00
            -1.8532E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1479.68180209386        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      362
 NPARAMETR:  1.0526E+00  9.5964E-01  8.7311E-01  1.1316E+00  8.4978E-01  1.0398E+00  1.1628E+00  3.3267E-01  1.0582E+00  4.2815E-01
             2.2893E+00
 PARAMETER:  1.5126E-01  5.8800E-02 -3.5694E-02  2.2365E-01 -6.2783E-02  1.3907E-01  2.5083E-01 -1.0006E+00  1.5653E-01 -7.4829E-01
             9.2826E-01
 GRADIENT:   2.5924E+00  1.0405E+01  3.6639E+00  4.0203E+00 -1.0388E+01  5.7905E-01 -4.5025E-01  5.1979E-01 -1.9372E+00  1.2061E+00
             3.3014E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1481.01428724906        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  1.0447E+00  6.0632E-01  9.4254E-01  1.3365E+00  7.6406E-01  1.0276E+00  1.6710E+00  2.4906E-02  9.3733E-01  3.3997E-01
             2.3155E+00
 PARAMETER:  1.4374E-01 -4.0035E-01  4.0822E-02  3.9003E-01 -1.6911E-01  1.2725E-01  6.1344E-01 -3.5927E+00  3.5280E-02 -9.7891E-01
             9.3964E-01
 GRADIENT:  -3.8714E+00  4.2478E+00  7.6962E+00 -1.6490E+00 -9.2916E+00 -1.7753E+00  9.4070E-01  5.5913E-04 -3.8657E-01 -1.2223E+00
            -2.6922E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1481.87897909328        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.0415E+00  3.9947E-01  9.1826E-01  1.4588E+00  6.9413E-01  1.0316E+00  2.0913E+00  1.0000E-02  8.9023E-01  5.3414E-01
             2.2652E+00
 PARAMETER:  1.4070E-01 -8.1762E-01  1.4729E-02  4.7759E-01 -2.6509E-01  1.3115E-01  8.3778E-01 -7.4006E+00 -1.6276E-02 -5.2710E-01
             9.1765E-01
 GRADIENT:  -4.8928E+00  3.5775E+00 -1.5144E-01  7.2862E+00 -1.7359E+00 -1.1535E-01  1.3256E+00  0.0000E+00 -5.0542E-01  3.2536E-01
             2.0189E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1482.53912139499        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      893
 NPARAMETR:  1.0452E+00  1.8143E-01  8.5550E-01  1.5635E+00  6.0942E-01  1.0272E+00  3.0804E+00  1.0000E-02  8.5402E-01  5.7866E-01
             2.2341E+00
 PARAMETER:  1.4416E-01 -1.6069E+00 -5.6070E-02  5.4692E-01 -3.9524E-01  1.2683E-01  1.2250E+00 -1.3991E+01 -5.7803E-02 -4.4703E-01
             9.0384E-01
 GRADIENT:   1.1275E+01  3.6961E-01  1.3645E+01  9.0592E+00 -1.9777E+01 -9.1358E-01 -4.0538E+00  0.0000E+00  1.1932E+00  1.0474E+00
             2.5713E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1483.82969489220        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  1.0284E+00  4.2084E-02  7.4847E-01  1.6363E+00  5.3448E-01  1.0434E+00  6.3098E+00  1.0000E-02  8.6151E-01  5.3271E-01
             2.1992E+00
 PARAMETER:  1.2805E-01 -3.0681E+00 -1.8973E-01  5.9243E-01 -5.2645E-01  1.4252E-01  1.9421E+00 -2.6897E+01 -4.9071E-02 -5.2978E-01
             8.8809E-01
 GRADIENT:  -1.4676E+01 -2.2304E+00  1.2883E+01  6.8276E+01 -2.4927E+01  4.6960E+00 -8.8123E+00  0.0000E+00  9.0843E+00  7.5078E-01
             1.0854E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1485.98805319579        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  1.0323E+00  2.1742E-02  7.2155E-01  1.5989E+00  5.2190E-01  1.0252E+00  8.4543E+00  1.0000E-02  8.2599E-01  5.0957E-01
             2.2138E+00
 PARAMETER:  1.3181E-01 -3.7285E+00 -2.2635E-01  5.6932E-01 -5.5028E-01  1.2485E-01  2.2347E+00 -3.2886E+01 -9.1171E-02 -5.7419E-01
             8.9470E-01
 GRADIENT:  -3.6940E+00  7.9645E-01  1.1333E-01  2.0210E+00  1.2062E+00 -1.0040E+00  1.7615E+00  0.0000E+00 -1.6830E+00 -3.5479E-01
             6.6772E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1486.05096679325        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1428
 NPARAMETR:  1.0341E+00  2.2188E-02  6.6555E-01  1.5793E+00  4.9150E-01  1.0296E+00  8.2589E+00  1.0000E-02  8.3331E-01  5.2207E-01
             2.1865E+00
 PARAMETER:  1.3351E-01 -3.7082E+00 -3.0714E-01  5.5698E-01 -6.1030E-01  1.2919E-01  2.2113E+00 -3.2766E+01 -8.2347E-02 -5.4995E-01
             8.8230E-01
 GRADIENT:  -7.2439E-01 -2.3751E-01  4.0418E-01  3.8017E-02 -4.8253E-01  1.3383E-01 -4.7864E-01  0.0000E+00  5.2344E-01 -2.4667E-02
            -6.4748E-02

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1486.05096679325        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1456
 NPARAMETR:  1.0345E+00  2.1912E-02  6.6628E-01  1.5822E+00  4.9049E-01  1.0294E+00  8.2014E+00  1.0000E-02  8.3248E-01  5.2311E-01
             2.1932E+00
 PARAMETER:  1.3351E-01 -3.7082E+00 -3.0714E-01  5.5698E-01 -6.1030E-01  1.2919E-01  2.2113E+00 -3.2766E+01 -8.2347E-02 -5.4995E-01
             8.8230E-01
 GRADIENT:  -1.1660E+00  3.8820E+01 -2.9699E+02 -1.5618E+02  1.4640E+02  1.1204E-01  6.1160E+01  0.0000E+00  3.6202E-01 -1.6659E+02
            -1.0115E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1456
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2242E-04  1.6178E-02  3.7822E-05 -1.2467E-02 -7.8214E-03
 SE:             2.9358E-02  8.9076E-03  2.3420E-04  2.7431E-02  1.7532E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8309E-01  6.9337E-02  8.7170E-01  6.4948E-01  6.5550E-01

 ETASHRINKSD(%)  1.6456E+00  7.0158E+01  9.9215E+01  8.1037E+00  4.1267E+01
 ETASHRINKVR(%)  3.2641E+00  9.1095E+01  9.9994E+01  1.5551E+01  6.5504E+01
 EBVSHRINKSD(%)  1.6643E+00  8.0235E+01  9.9155E+01  7.7449E+00  4.0724E+01
 EBVSHRINKVR(%)  3.3008E+00  9.6093E+01  9.9993E+01  1.4890E+01  6.4863E+01
 RELATIVEINF(%)  9.6342E+01  2.9650E+00  2.7218E-04  4.3633E+01  1.3314E+00
 EPSSHRINKSD(%)  3.2188E+01
 EPSSHRINKVR(%)  5.4015E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1486.0509667932483     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -750.90014022951016     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1486.051       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.22E-02  6.66E-01  1.58E+00  4.91E-01  1.03E+00  8.26E+00  1.00E-02  8.33E-01  5.22E-01  2.19E+00
 


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
+        9.81E+02
 
 TH 2
+       -3.21E+02  9.04E+05
 
 TH 3
+        5.48E+01  4.06E+03  5.63E+04
 
 TH 4
+       -1.71E+01  2.29E+03 -2.21E+02  5.43E+03
 
 TH 5
+        1.69E+01 -7.92E+03 -2.70E+03 -5.69E+03  3.06E+04
 
 TH 6
+       -1.24E+01  8.07E+01 -2.49E+01 -1.17E+01  1.79E+01  1.80E+02
 
 TH 7
+       -1.02E+00  2.56E+03 -6.19E+02 -1.37E+02  3.97E+02  2.91E-01  1.86E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.76E+01  5.18E+02 -1.84E+01  3.63E+00 -1.16E+01  3.11E+01  1.69E+00  0.00E+00  4.24E+02
 
 TH10
+        4.27E+01  1.11E+03  3.89E+04 -5.62E+01  2.23E+02 -2.23E+01  4.45E+00  0.00E+00 -2.12E+01  2.78E+04
 
 TH11
+       -5.21E+00  5.87E+02 -5.54E+01  2.23E+03 -2.58E+03  1.78E+00 -6.45E+01  0.00E+00  1.59E+01  1.39E+01  1.09E+03
 
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
 #CPUT: Total CPU Time in Seconds,       25.233
Stop Time:
Wed Sep 29 12:57:40 CDT 2021
