Sat Sep 25 07:03:52 CDT 2021
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
$DATA ../../../../data/spa/B/dat2.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.18740710209        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.5192E+01 -4.9925E+01 -2.0877E+01 -5.2867E+01 -3.6933E+01  2.1795E+00  8.7191E+00  1.7518E+01  3.5420E+01  1.7916E+01
             5.4075E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1698.27018980628        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7992E-01  1.1502E+00  1.2802E+00  9.6365E-01  1.2622E+00  9.8662E-01  9.2797E-01  8.3606E-01  7.2537E-01  9.6287E-01
             1.0450E+00
 PARAMETER:  7.9717E-02  2.3994E-01  3.4703E-01  6.2972E-02  3.3286E-01  8.6531E-02  2.5244E-02 -7.9053E-02 -2.2108E-01  6.2159E-02
             1.4406E-01
 GRADIENT:   3.6481E+01  8.2049E+00 -3.9262E+00  1.3970E+01  5.7726E+01 -1.9501E+00 -1.0671E+01 -2.8337E+00 -2.3312E+01 -3.1315E+01
             1.8338E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1702.67381347493        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8858E-01  9.9766E-01  1.1665E+00  1.0611E+00  1.1468E+00  9.8350E-01  9.3051E-01  2.9532E-01  8.5155E-01  1.0721E+00
             1.0077E+00
 PARAMETER:  8.8515E-02  9.7657E-02  2.5401E-01  1.5932E-01  2.3695E-01  8.3358E-02  2.7977E-02 -1.1197E+00 -6.0698E-02  1.6964E-01
             1.0763E-01
 GRADIENT:   6.0902E+01 -4.1234E+00 -1.6628E+01  3.5233E+01  4.5352E+01 -3.8630E+00 -2.5668E+00 -5.0700E-01  3.3668E+00 -6.3218E+00
            -2.9412E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1703.86972776452        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      250
 NPARAMETR:  9.6853E-01  1.0325E+00  1.0660E+00  1.0257E+00  1.0838E+00  9.9034E-01  9.9852E-01  3.2475E-01  8.2280E-01  1.0189E+00
             1.0026E+00
 PARAMETER:  6.8029E-02  1.3195E-01  1.6391E-01  1.2538E-01  1.8049E-01  9.0295E-02  9.8516E-02 -1.0247E+00 -9.5044E-02  1.1870E-01
             1.0261E-01
 GRADIENT:  -2.6239E+01 -4.9573E+00 -5.9998E+00  1.9597E+00  1.0565E+01 -4.9516E+00 -6.1650E-01  1.7968E-02 -4.0367E-01 -1.2665E+00
            -9.1475E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1704.18839126411        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      427
 NPARAMETR:  9.7986E-01  1.1394E+00  1.0469E+00  9.5980E-01  1.1212E+00  1.0010E+00  9.1965E-01  3.7362E-01  8.7123E-01  1.0422E+00
             1.0060E+00
 PARAMETER:  7.9659E-02  2.3052E-01  1.4580E-01  5.8965E-02  2.1439E-01  1.0097E-01  1.6241E-02 -8.8452E-01 -3.7846E-02  1.4137E-01
             1.0601E-01
 GRADIENT:  -1.3512E+00  2.1091E-01 -3.5674E-01  1.3653E+00  1.1584E-02 -5.0444E-01 -4.4236E-02  3.4481E-02 -4.2409E-01  8.0965E-02
             2.8306E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1704.37749724883        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.8149E-01  1.4100E+00  9.2365E-01  7.8532E-01  1.2125E+00  1.0022E+00  7.7370E-01  3.4565E-01  1.0316E+00  1.0803E+00
             1.0066E+00
 PARAMETER:  8.1320E-02  4.4358E-01  2.0581E-02 -1.4166E-01  2.9268E-01  1.0223E-01 -1.5658E-01 -9.6233E-01  1.3115E-01  1.7723E-01
             1.0659E-01
 GRADIENT:  -5.4178E-01  3.5003E+00  9.8016E-01  2.3693E+00 -1.8937E+00 -3.4666E-01  1.7789E-02  8.2518E-02  9.5845E-02 -3.6627E-02
             1.4365E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1704.49247009584        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  9.8264E-01  1.6230E+00  7.8213E-01  6.4441E-01  1.2818E+00  1.0036E+00  7.0344E-01  1.8738E-01  1.1913E+00  1.1052E+00
             1.0065E+00
 PARAMETER:  8.2487E-02  5.8430E-01 -1.4573E-01 -3.3943E-01  3.4824E-01  1.0362E-01 -2.5178E-01 -1.5746E+00  2.7504E-01  2.0007E-01
             1.0646E-01
 GRADIENT:   2.0556E-01  3.7264E+00 -7.2580E-01  2.8379E+00  9.3731E-01  2.8582E-03 -3.8933E-01  6.0938E-02 -3.4890E-01 -1.4143E-01
             3.9653E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1704.50459804872        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  9.8260E-01  1.6426E+00  7.7154E-01  6.2815E-01  1.2898E+00  1.0034E+00  6.9767E-01  1.5114E-01  1.2207E+00  1.1103E+00
             1.0067E+00
 PARAMETER:  8.2443E-02  5.9628E-01 -1.5936E-01 -3.6497E-01  3.5448E-01  1.0343E-01 -2.6001E-01 -1.7896E+00  2.9939E-01  2.0463E-01
             1.0672E-01
 GRADIENT:   1.3826E-01 -1.6825E+00  1.8370E-02 -1.0181E+00 -4.5447E-01 -6.2078E-02  1.9548E-01  3.9545E-02  4.4073E-01  9.4290E-02
             1.9351E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1704.52651625791        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  9.8253E-01  1.6218E+00  7.7946E-01  6.4257E-01  1.2799E+00  1.0036E+00  7.0450E-01  3.0427E-02  1.1960E+00  1.1050E+00
             1.0063E+00
 PARAMETER:  8.2372E-02  5.8356E-01 -1.4915E-01 -3.4228E-01  3.4678E-01  1.0357E-01 -2.5027E-01 -3.3924E+00  2.7901E-01  1.9988E-01
             1.0630E-01
 GRADIENT:   2.6618E-02 -4.3137E-01 -8.5235E-02 -1.3699E-01  1.5427E-01 -1.7599E-02 -7.5670E-03  1.4899E-03  2.3849E-02 -2.3694E-03
             1.4039E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1704.52723368486        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1304
 NPARAMETR:  9.8252E-01  1.6208E+00  7.8031E-01  6.4338E-01  1.2795E+00  1.0036E+00  7.0475E-01  1.0000E-02  1.1947E+00  1.1050E+00
             1.0063E+00
 PARAMETER:  8.2362E-02  5.8296E-01 -1.4800E-01 -3.4102E-01  3.4649E-01  1.0360E-01 -2.4966E-01 -4.5387E+00  2.7811E-01  1.9989E-01
             1.0633E-01
 GRADIENT:  -5.0011E-03  1.3294E-02  5.8479E-03 -2.4722E-03  2.2909E-02 -4.4750E-03  1.1896E-02  0.0000E+00  1.1627E-02  4.9305E-03
             5.9006E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1304
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.4610E-04 -2.7458E-02 -2.5822E-04  1.6871E-02 -3.6441E-02
 SE:             2.9822E-02  2.1234E-02  1.0473E-04  2.2665E-02  2.3122E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9342E-01  1.9597E-01  1.3682E-02  4.5666E-01  1.1503E-01

 ETASHRINKSD(%)  9.1202E-02  2.8864E+01  9.9649E+01  2.4069E+01  2.2537E+01
 ETASHRINKVR(%)  1.8232E-01  4.9397E+01  9.9999E+01  4.2345E+01  3.9995E+01
 EBVSHRINKSD(%)  4.2653E-01  2.7441E+01  9.9686E+01  2.6205E+01  2.0318E+01
 EBVSHRINKVR(%)  8.5124E-01  4.7352E+01  9.9999E+01  4.5542E+01  3.6508E+01
 RELATIVEINF(%)  9.8960E+01  1.9273E+00  1.1552E-04  2.1559E+00  1.3296E+01
 EPSSHRINKSD(%)  4.2210E+01
 EPSSHRINKVR(%)  6.6604E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1704.5272336848589     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -969.37640712112068     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1704.527       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.62E+00  7.80E-01  6.43E-01  1.28E+00  1.00E+00  7.05E-01  1.00E-02  1.19E+00  1.11E+00  1.01E+00
 


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
+        1.13E+03
 
 TH 2
+       -7.75E+00  4.70E+02
 
 TH 3
+        8.98E+00  1.16E+02  2.06E+02
 
 TH 4
+       -1.50E+01  4.89E+02 -1.57E+02  1.04E+03
 
 TH 5
+       -2.34E+00 -1.51E+02 -1.89E+02  1.45E+02  3.54E+02
 
 TH 6
+       -2.39E+00 -1.58E+00  2.37E+00 -4.47E+00 -1.24E+00  1.95E+02
 
 TH 7
+       -2.91E-01  3.21E+00  1.09E+01 -1.67E+01 -1.44E+01 -1.81E-01  1.18E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.41E+00 -1.84E+01 -2.26E+01  4.89E+01  4.74E+00 -7.48E-01  3.35E+01  0.00E+00  4.92E+01
 
 TH10
+        1.46E+00 -1.03E+01 -2.15E+01 -3.71E+00 -4.76E+01  1.18E+00  9.49E+00  0.00E+00  4.45E+00  7.29E+01
 
 TH11
+       -2.77E+00 -2.59E+01 -3.47E+01  3.93E+00  2.70E+00  4.70E+00  9.64E+00  0.00E+00  9.71E+00  1.62E+01  2.24E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.349
Stop Time:
Sat Sep 25 07:04:15 CDT 2021
