Fri Sep 24 23:15:58 CDT 2021
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
$DATA ../../../../data/int/S1/dat36.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3382.43933114987        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0936E+01 -6.6128E+01  1.3070E+02 -3.3955E+01  8.5381E+01 -2.9011E+01 -5.2850E+01 -4.6676E+02 -9.5264E+01 -3.2356E+01
            -3.7447E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3627.65436301507        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0211E+00  1.0856E+00  9.1611E-01  1.0006E+00  1.0473E+00  1.0686E+00  1.2369E+00  3.2696E+00  8.6504E-01  1.2304E+00
             1.1493E+00
 PARAMETER:  1.2086E-01  1.8211E-01  1.2386E-02  1.0063E-01  1.4622E-01  1.6632E-01  3.1263E-01  1.2847E+00 -4.4976E-02  3.0734E-01
             2.3915E-01
 GRADIENT:   5.9005E+01  6.1368E+00 -2.0465E+01  6.3974E+01  1.0923E+01  1.4894E+00  7.1990E+00  3.5121E+01  8.6520E+00 -1.1745E+01
             5.6756E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3629.97174215385        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0177E+00  1.1545E+00  9.6819E-01  9.7350E-01  1.0391E+00  1.0955E+00  1.1856E+00  3.3135E+00  8.0937E-01  1.4850E+00
             1.1437E+00
 PARAMETER:  1.1751E-01  2.4364E-01  6.7676E-02  7.3140E-02  1.3833E-01  1.9123E-01  2.7028E-01  1.2980E+00 -1.1150E-01  4.9544E-01
             2.3423E-01
 GRADIENT:   5.1315E+01  4.2008E+01 -4.7821E+00  7.0660E+01 -3.2342E+01  1.3234E+01 -2.0864E+00  3.5389E+01  8.9699E+00  1.3355E+01
             4.8393E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3630.63701602968        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      245
 NPARAMETR:  1.0152E+00  1.1486E+00  9.7277E-01  9.7226E-01  1.0391E+00  1.0874E+00  1.1732E+00  3.2886E+00  7.9738E-01  1.4726E+00
             1.1409E+00
 PARAMETER:  1.1509E-01  2.3851E-01  7.2391E-02  7.1869E-02  1.3836E-01  1.8383E-01  2.5977E-01  1.2905E+00 -1.2643E-01  4.8704E-01
             2.3185E-01
 GRADIENT:   4.6293E+01  3.6414E+01 -3.6641E+00  6.4872E+01 -3.1111E+01  1.0048E+01 -5.0791E+00  3.3658E+01  6.1679E+00  1.1348E+01
             4.3103E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3630.68560754161        NO. OF FUNC. EVALS.: 220
 CUMULATIVE NO. OF FUNC. EVALS.:      465             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0152E+00  1.1486E+00  9.7277E-01  9.7194E-01  1.0391E+00  1.0875E+00  1.1757E+00  3.2885E+00  7.9738E-01  1.4726E+00
             1.1406E+00
 PARAMETER:  1.1506E-01  2.3851E-01  7.2392E-02  7.1542E-02  1.3836E-01  1.8389E-01  2.6190E-01  1.2904E+00 -1.2643E-01  4.8704E-01
             2.3159E-01
 GRADIENT:   4.6250E+01  3.6246E+01 -3.6598E+00  6.4173E+01 -3.1218E+01  1.0083E+01 -4.6860E+00  3.3638E+01  6.3424E+00  1.1283E+01
             4.2716E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3630.69614821403        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  1.0152E+00  1.1486E+00  9.7279E-01  9.7197E-01  1.0391E+00  1.0888E+00  1.1759E+00  3.2874E+00  7.9740E-01  1.4724E+00
             1.1406E+00
 PARAMETER:  1.1509E-01  2.3857E-01  7.2417E-02  7.1567E-02  1.3833E-01  1.8509E-01  2.6203E-01  1.2901E+00 -1.2639E-01  4.8691E-01
             2.3153E-01
 GRADIENT:   5.4692E+00  1.7668E+01 -4.1161E+00  5.7351E+01 -3.6579E+01 -1.4340E-01 -8.4256E+00  2.9026E+01  5.8676E+00  7.4630E+00
             4.2196E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3631.11309375583        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:      806
 NPARAMETR:  1.0152E+00  1.1485E+00  9.7278E-01  9.7193E-01  1.0391E+00  1.0718E+00  1.1976E+00  3.2733E+00  7.9738E-01  1.4723E+00
             1.1404E+00
 PARAMETER:  1.1507E-01  2.3848E-01  7.2406E-02  7.1525E-02  1.3838E-01  1.6930E-01  2.8032E-01  1.2858E+00 -1.2642E-01  4.8683E-01
             2.3136E-01
 GRADIENT:   4.6444E+01  3.5853E+01 -6.2040E+00  6.5450E+01 -3.3329E+01  3.3654E+00 -2.6485E+00  3.0613E+01  7.4098E+00  1.3322E+01
             4.2635E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3631.13453849681        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      956
 NPARAMETR:  1.0148E+00  1.1485E+00  9.7279E-01  9.7193E-01  1.0391E+00  1.0731E+00  1.1998E+00  3.2734E+00  7.9739E-01  1.4723E+00
             1.1404E+00
 PARAMETER:  1.1469E-01  2.3848E-01  7.2409E-02  7.1528E-02  1.3838E-01  1.7060E-01  2.8214E-01  1.2858E+00 -1.2641E-01  4.8683E-01
             2.3135E-01
 GRADIENT:   4.8632E+00  1.7254E+01 -6.6740E+00  5.8508E+01 -3.8648E+01 -5.9559E+00 -6.4120E+00  2.6082E+01  7.0466E+00  9.5637E+00
             4.2325E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3631.17295804077        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1123
 NPARAMETR:  1.0148E+00  1.1487E+00  9.7283E-01  9.7198E-01  1.0391E+00  1.0890E+00  1.1996E+00  3.2734E+00  7.9744E-01  1.4723E+00
             1.1404E+00
 PARAMETER:  1.1464E-01  2.3860E-01  7.2459E-02  7.1578E-02  1.3838E-01  1.8526E-01  2.8200E-01  1.2858E+00 -1.2635E-01  4.8683E-01
             2.3135E-01
 GRADIENT:   4.6683E+00  1.7372E+01 -6.6623E+00  5.8673E+01 -3.8680E+01 -8.3985E-02 -6.4414E+00  2.6078E+01  7.0376E+00  9.5734E+00
             4.2355E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3631.20261071509        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1301             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0147E+00  1.1487E+00  9.7283E-01  9.7180E-01  1.0391E+00  1.0893E+00  1.1998E+00  3.2724E+00  7.9741E-01  1.4723E+00
             1.1401E+00
 PARAMETER:  1.1459E-01  2.3860E-01  7.2459E-02  7.1390E-02  1.3838E-01  1.8553E-01  2.8212E-01  1.2855E+00 -1.2638E-01  4.8682E-01
             2.3113E-01
 GRADIENT:   4.5275E+01  3.5931E+01 -6.1811E+00  6.5199E+01 -3.3392E+01  1.0836E+01 -2.3089E+00  3.0525E+01  7.5300E+00  1.3280E+01
             4.2347E+01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -3631.20261071509        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1360
 NPARAMETR:  1.0147E+00  1.1487E+00  9.7283E-01  9.7180E-01  1.0391E+00  1.0893E+00  1.1998E+00  3.2724E+00  7.9741E-01  1.4723E+00
             1.1401E+00
 PARAMETER:  1.1459E-01  2.3860E-01  7.2459E-02  7.1390E-02  1.3838E-01  1.8553E-01  2.8212E-01  1.2855E+00 -1.2638E-01  4.8682E-01
             2.3113E-01
 GRADIENT:   4.5775E+00  1.7250E+01 -6.6474E+00  5.8290E+01 -3.8724E+01  2.1161E-02 -6.4296E+00  2.5996E+01  7.0577E+00  9.5483E+00
             4.1941E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1360
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.8141E-03 -4.5315E-02 -2.3585E-02  7.3749E-03 -2.2051E-02
 SE:             2.9886E-02  2.4518E-02  2.3166E-02  2.3894E-02  2.3945E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5160E-01  6.4566E-02  3.0862E-01  7.5759E-01  3.5709E-01

 ETASHRINKSD(%)  1.0000E-10  1.7863E+01  2.2392E+01  1.9953E+01  1.9782E+01
 ETASHRINKVR(%)  1.0000E-10  3.2535E+01  3.9770E+01  3.5925E+01  3.5651E+01
 EBVSHRINKSD(%)  2.7861E-01  2.0416E+01  1.3618E+01  2.1518E+01  1.5548E+01
 EBVSHRINKVR(%)  5.5644E-01  3.6664E+01  2.5382E+01  3.8406E+01  2.8678E+01
 RELATIVEINF(%)  9.9442E+01  3.1465E+01  7.0195E+01  3.2324E+01  4.4289E+01
 EPSSHRINKSD(%)  2.4938E+01
 EPSSHRINKVR(%)  4.3656E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3631.2026107150887     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1977.1132509466779     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    46.07
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3631.203       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.15E+00  9.73E-01  9.72E-01  1.04E+00  1.09E+00  1.20E+00  3.27E+00  7.97E-01  1.47E+00  1.14E+00
 


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
+        4.08E+08
 
 TH 2
+        1.72E+08  1.47E+08
 
 TH 3
+       -4.86E+08 -2.24E+03  1.16E+09
 
 TH 4
+        4.86E+08  2.06E+08 -5.78E+08  5.84E+08
 
 TH 5
+        1.57E+06 -2.20E+03 -5.89E+03 -5.70E+03  5.33E+08
 
 TH 6
+       -1.28E+01  4.48E+02 -1.25E+03  9.69E+00  9.04E+05  1.67E+02
 
 TH 7
+       -1.40E+08 -5.92E+07  1.67E+08 -1.30E+03  1.69E+03 -3.03E+00  2.30E+05
 
 TH 8
+        2.25E+07 -5.35E+01 -2.68E+07 -1.40E+02 -1.43E+02  2.93E+01 -5.01E+02  6.22E+05
 
 TH 9
+        2.83E+01  1.99E+08  1.12E+09  5.64E+08 -5.68E+03  1.22E+03 -1.77E+04 -2.59E+07  1.08E+09
 
 TH10
+       -2.15E+00  2.85E+02  9.02E+02 -5.75E+02  1.07E+08 -1.12E+00  1.08E+05  8.72E+01  4.03E+01  2.14E+07
 
 TH11
+       -1.78E+08 -3.63E+05  2.14E+08 -2.15E+08  6.94E+05 -4.67E+02  7.16E+03  4.96E+06 -2.07E+08  1.74E+02  1.59E+08
 
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
 #CPUT: Total CPU Time in Seconds,       62.274
Stop Time:
Fri Sep 24 23:17:02 CDT 2021
