Sat Sep 25 09:53:04 CDT 2021
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
$DATA ../../../../data/spa/S1/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1696.64016837550        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3447E+01 -1.0503E+01 -3.4569E+01  3.1759E+01  4.4433E+01  5.4840E+01 -6.2828E+00  7.5470E+00 -7.8788E+00  1.1939E+01
             2.8457E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1702.33063126884        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0236E+00  1.0294E+00  1.0439E+00  9.6578E-01  1.0077E+00  8.3291E-01  1.0583E+00  9.4463E-01  1.0623E+00  8.9099E-01
             1.0087E+00
 PARAMETER:  1.2328E-01  1.2893E-01  1.4292E-01  6.5179E-02  1.0771E-01 -8.2824E-02  1.5668E-01  4.3036E-02  1.6045E-01 -1.5427E-02
             1.0864E-01
 GRADIENT:   1.1336E+02 -9.2735E+00 -2.6915E+00  3.0029E+00  2.5109E+01 -1.0433E+01 -1.5009E+00  9.9461E-01  5.1335E+00 -4.0663E+00
             1.3463E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1703.69371833978        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0188E+00  1.0144E+00  7.4769E-01  9.6573E-01  8.3824E-01  8.4103E-01  1.1911E+00  5.8887E-01  9.9443E-01  6.7888E-01
             9.8764E-01
 PARAMETER:  1.1859E-01  1.1428E-01 -1.9077E-01  6.5134E-02 -7.6456E-02 -7.3123E-02  2.7490E-01 -4.2955E-01  9.4415E-02 -2.8731E-01
             8.7562E-02
 GRADIENT:   9.0004E+01  8.7172E+00  2.1566E+00  2.2975E+01 -3.2599E+00 -6.5715E+00 -3.7745E+00  1.0308E+00  6.5929E+00 -7.5948E+00
            -5.8942E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1705.06422036194        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.9630E-01  1.0788E+00  7.0734E-01  9.1248E-01  8.5529E-01  8.5322E-01  1.1609E+00  3.4208E-01  9.9093E-01  7.7479E-01
             9.9177E-01
 PARAMETER:  9.6295E-02  1.7587E-01 -2.4624E-01  8.4087E-03 -5.6314E-02 -5.8737E-02  2.4921E-01 -9.7270E-01  9.0892E-02 -1.5517E-01
             9.1739E-02
 GRADIENT:   1.1599E+01  3.6470E+00 -2.1863E+00  2.3470E+00 -3.9414E+00 -6.1629E-01  5.7087E-01  8.0597E-01  8.7886E-01  2.5238E+00
            -4.2136E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1705.07085295658        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      295
 NPARAMETR:  9.9397E-01  1.0640E+00  7.1665E-01  9.2080E-01  8.5541E-01  8.5428E-01  1.1752E+00  3.0435E-01  9.8307E-01  7.7699E-01
             9.9317E-01
 PARAMETER:  9.3956E-02  1.6200E-01 -2.3317E-01  1.7483E-02 -5.6172E-02 -5.7494E-02  2.6140E-01 -1.0896E+00  8.2927E-02 -1.5233E-01
             9.3145E-02
 GRADIENT:   4.8394E+00  1.3284E+00 -1.6454E+00  1.0414E+00 -1.5451E+00 -2.3377E-01  4.4033E-01  5.5742E-01  3.6775E-01  1.5524E+00
            -4.9389E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1705.08005861385        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      365
 NPARAMETR:  9.9218E-01  1.0528E+00  7.1952E-01  9.2645E-01  8.5362E-01  8.5520E-01  1.1865E+00  2.3343E-01  9.7705E-01  7.7785E-01
             9.9423E-01
 PARAMETER:  9.2153E-02  1.5145E-01 -2.2917E-01  2.3603E-02 -5.8273E-02 -5.6424E-02  2.7104E-01 -1.3549E+00  7.6782E-02 -1.5122E-01
             9.4211E-02
 GRADIENT:  -4.4334E-01 -5.7848E-01 -1.0950E+00  1.2662E-03  4.8963E-01  6.1435E-02  2.6878E-01  2.8140E-01 -6.9432E-03  6.3366E-01
             2.1836E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1705.17330778023        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  9.9213E-01  1.0563E+00  7.0731E-01  9.2282E-01  8.4904E-01  8.5542E-01  1.1854E+00  4.4559E-02  9.7819E-01  7.7389E-01
             9.9403E-01
 PARAMETER:  9.2103E-02  1.5473E-01 -2.4629E-01  1.9678E-02 -6.3653E-02 -5.6163E-02  2.7007E-01 -3.0109E+00  7.7947E-02 -1.5632E-01
             9.4008E-02
 GRADIENT:  -8.8609E-01 -4.7771E-01 -3.7728E-01 -1.1402E-01  6.1741E-01  9.1891E-02  1.0814E-01  9.8296E-03 -3.8913E-03  3.0223E-02
             8.7209E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1705.66463646697        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  1.0086E+00  1.1223E+00  7.0891E-01  8.9127E-01  8.8071E-01  8.6267E-01  1.1344E+00  1.0000E-02  1.0160E+00  7.9985E-01
             9.9552E-01
 PARAMETER:  1.0860E-01  2.1542E-01 -2.4403E-01 -1.5106E-02 -2.7027E-02 -4.7728E-02  2.2612E-01 -1.0120E+01  1.1588E-01 -1.2333E-01
             9.5506E-02
 GRADIENT:   3.2119E+00 -1.0977E-02 -6.1518E-01  1.6444E+00  2.0406E-01  1.2286E-01  5.2647E-04  0.0000E+00 -1.1232E-01  4.3586E-01
             1.1246E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1705.70713774972        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      775
 NPARAMETR:  1.0077E+00  1.2137E+00  6.8267E-01  8.3351E-01  9.1169E-01  8.6250E-01  1.0613E+00  1.0000E-02  1.0735E+00  8.0544E-01
             9.9619E-01
 PARAMETER:  1.0770E-01  2.9368E-01 -2.8175E-01 -8.2110E-02  7.5437E-03 -4.7915E-02  1.5952E-01 -1.1912E+01  1.7088E-01 -1.1637E-01
             9.6186E-02
 GRADIENT:  -1.0366E-01  2.3448E-01  1.3159E-01  1.0587E-01 -1.8861E-01 -6.2710E-04 -1.7274E-02  0.0000E+00 -2.7245E-02 -3.7763E-02
            -1.4188E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1705.70720811005        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  1.0078E+00  1.2168E+00  6.8147E-01  8.3143E-01  9.1274E-01  8.6251E-01  1.0590E+00  1.0000E-02  1.0756E+00  8.0583E-01
             9.9620E-01
 PARAMETER:  1.0774E-01  2.9621E-01 -2.8350E-01 -8.4610E-02  8.6975E-03 -4.7907E-02  1.5736E-01 -1.1960E+01  1.7287E-01 -1.1589E-01
             9.6197E-02
 GRADIENT:   6.6686E-03  9.0410E-03  2.4287E-03  4.2919E-03 -6.1992E-03  5.6724E-04  6.7208E-04  0.0000E+00 -6.8176E-04  5.2483E-05
            -7.0738E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      867
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.9226E-05 -1.1587E-02 -3.8102E-04  6.9684E-03 -2.0399E-02
 SE:             2.9803E-02  2.2923E-02  1.5360E-04  2.4119E-02  2.1968E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9761E-01  6.1321E-01  1.3115E-02  7.7265E-01  3.5310E-01

 ETASHRINKSD(%)  1.5509E-01  2.3206E+01  9.9485E+01  1.9199E+01  2.6405E+01
 ETASHRINKVR(%)  3.0994E-01  4.1026E+01  9.9997E+01  3.4711E+01  4.5838E+01
 EBVSHRINKSD(%)  5.6429E-01  2.2791E+01  9.9531E+01  1.9657E+01  2.5974E+01
 EBVSHRINKVR(%)  1.1254E+00  4.0388E+01  9.9998E+01  3.5450E+01  4.5202E+01
 RELATIVEINF(%)  9.8639E+01  3.7257E+00  2.2906E-04  4.5495E+00  6.1035E+00
 EPSSHRINKSD(%)  4.3966E+01
 EPSSHRINKVR(%)  6.8602E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1705.7072081100503     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -970.55638154631208     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1705.707       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.22E+00  6.81E-01  8.31E-01  9.13E-01  8.63E-01  1.06E+00  1.00E-02  1.08E+00  8.06E-01  9.96E-01
 


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
+        1.45E+03
 
 TH 2
+       -6.84E+00  4.16E+02
 
 TH 3
+        1.88E+01  2.21E+02  5.70E+02
 
 TH 4
+       -1.93E+01  3.01E+02 -3.27E+02  9.10E+02
 
 TH 5
+       -7.74E+00 -3.48E+02 -6.69E+02  3.78E+02  1.12E+03
 
 TH 6
+        1.77E+00 -1.27E+00  2.80E+00 -1.63E+00 -3.19E+00  2.64E+02
 
 TH 7
+        9.96E-01  2.61E+01 -1.26E+01 -1.20E+01 -3.50E+00 -3.63E-01  6.76E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.17E+00 -2.33E+01 -3.41E+01  4.44E+01 -4.85E+00  8.18E-01  1.57E+01  0.00E+00  8.09E+01
 
 TH10
+        3.86E-02 -1.27E+01 -5.66E+01 -2.17E+01 -6.57E+01  2.47E+00  1.86E+01  0.00E+00  9.91E+00  1.00E+02
 
 TH11
+       -8.31E+00 -1.57E+01 -3.36E+01 -2.07E+00  7.01E-01  9.80E-01  7.75E+00  0.00E+00  7.48E+00  2.26E+01  2.16E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.118
Stop Time:
Sat Sep 25 09:53:19 CDT 2021
