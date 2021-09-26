Sat Sep 25 14:52:10 CDT 2021
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
$DATA ../../../../data/spa/D/dat100.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   31197.5524623488        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.7706E+02  7.5089E+02 -2.8440E+01  8.0345E+02  4.3412E+01 -2.9351E+03 -1.2421E+03 -4.4362E+01 -1.8178E+03 -5.7567E+02
            -5.8507E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -357.494228709462        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0392E+00  9.4032E-01  8.8222E-01  1.1906E+00  1.6610E+00  1.7846E+00  9.0063E-01  9.5119E-01  7.1519E-01  8.5396E-01
             1.5147E+01
 PARAMETER:  1.3844E-01  3.8470E-02 -2.5312E-02  2.7442E-01  6.0744E-01  6.7922E-01 -4.6589E-03  4.9956E-02 -2.3520E-01 -5.7867E-02
             2.8178E+00
 GRADIENT:  -2.8581E+01 -6.8001E-01 -3.4052E+00 -1.0280E+01 -2.5152E+00  4.3948E+01  1.2408E+00  2.4467E+00  4.4455E+00  3.6304E-01
            -5.3501E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -375.781644968029        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.1509E+00  6.4669E-01  6.8679E-01  1.4198E+00  3.1206E+00  1.4910E+00  3.3765E-01  1.1903E-01  1.5736E-01  1.2223E+00
             1.7227E+01
 PARAMETER:  2.4051E-01 -3.3589E-01 -2.7573E-01  4.5049E-01  1.2380E+00  4.9943E-01 -9.8574E-01 -2.0284E+00 -1.7492E+00  3.0073E-01
             2.9464E+00
 GRADIENT:   1.3080E+01  1.0838E+01 -9.5171E+00  1.9729E+00 -1.2550E+01 -4.7599E+00  6.3898E-02  6.8748E-02  3.7351E-01  1.4706E+00
            -7.7148E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -393.496245655607        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0496E+00  2.6711E-01  2.8062E-01  1.5646E+00  3.4814E+01  1.7526E+00  1.0133E-02  5.4399E-02  6.1723E-02  3.5636E+00
             1.9032E+01
 PARAMETER:  1.4840E-01 -1.2201E+00 -1.1707E+00  5.4760E-01  3.6500E+00  6.6109E-01 -4.4920E+00 -2.8114E+00 -2.6851E+00  1.3708E+00
             3.0461E+00
 GRADIENT:  -1.1997E+01 -1.1482E+01  1.2533E+01  1.0625E+02  7.6741E-01 -1.4353E+01  3.4568E-04  4.7434E-02  1.1653E-01  7.0853E-03
            -4.6023E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -414.266380452707        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  8.9055E-01  1.3562E-01  1.6868E-01  1.0863E+00  4.6365E+02  1.6378E+00  1.0000E-02  1.6064E-01  6.9366E-02  2.5350E+01
             1.7315E+01
 PARAMETER: -1.5919E-02 -1.8979E+00 -1.6797E+00  1.8274E-01  6.2391E+00  5.9336E-01 -5.6333E+00 -1.7286E+00 -2.5684E+00  3.3328E+00
             2.9516E+00
 GRADIENT:   2.1886E+00 -9.0832E+00  1.3825E+01  3.1666E+01  2.3350E-02 -5.5621E+00  0.0000E+00  4.2436E-01  2.1119E-01  1.1724E-03
            -1.0336E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -426.110580862474        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  4.8642E-01  1.6126E-02  3.0503E-02  3.6385E-01  2.2490E+06  1.5143E+00  1.0000E-02  3.0485E-02  1.3989E-02  2.3548E+03
             1.5895E+01
 PARAMETER: -6.2069E-01 -4.0273E+00 -3.3899E+00 -9.1102E-01  1.4726E+01  5.1493E-01 -1.0771E+01 -3.3905E+00 -4.1695E+00  7.8642E+00
             2.8660E+00
 GRADIENT:  -5.2948E+00  4.0488E-01 -1.9902E+01  3.6809E+01  1.5866E-07 -6.2218E+00  0.0000E+00  2.5024E-02  9.5625E-03 -1.3300E-08
            -1.5933E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -426.909544785595        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  4.8402E-01  1.4880E-02  2.8711E-02  3.4054E-01  3.9841E+06  1.5401E+00  1.0000E-02  2.1624E-02  1.0000E-02  3.2936E+03
             1.6337E+01
 PARAMETER: -6.2563E-01 -4.1077E+00 -3.4505E+00 -9.7721E-01  1.5298E+01  5.3187E-01 -1.1464E+01 -3.7340E+00 -4.5583E+00  8.1997E+00
             2.8935E+00
 GRADIENT:  -2.1569E-02  4.0659E-01 -4.5625E-01 -1.8063E-01 -3.1844E-09  1.8732E-01  0.0000E+00  1.3101E-02  0.0000E+00 -4.3258E-09
            -7.2531E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -426.959166764104        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  4.8998E-01  1.1687E-02  2.9768E-02  3.4990E-01  3.6118E+06  1.5425E+00  1.0000E-02  1.0000E-02  1.0261E-02  2.4550E+03
             1.6374E+01
 PARAMETER: -6.1339E-01 -4.3493E+00 -3.4143E+00 -9.5010E-01  1.5200E+01  5.3343E-01 -1.1770E+01 -4.9731E+00 -4.4794E+00  7.9059E+00
             2.8957E+00
 GRADIENT:  -3.2366E-02  1.1248E-02 -6.4810E-02  8.1122E-02  9.3943E-08 -1.3945E-01  0.0000E+00  0.0000E+00  4.9617E-03 -1.0857E-08
             2.3995E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -426.959567319567        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      888
 NPARAMETR:  4.8991E-01  1.1307E-02  2.9781E-02  3.5000E-01  3.6496E+06  1.5436E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.4009E+03
             1.6367E+01
 PARAMETER: -6.1354E-01 -4.3823E+00 -3.4139E+00 -9.4983E-01  1.5210E+01  5.3411E-01 -1.1808E+01 -5.1336E+00 -5.0627E+00  7.8836E+00
             2.8953E+00
 GRADIENT:   4.2844E-03  2.0552E-05 -2.1518E-04 -1.5069E-02  9.5691E-08  2.0065E-03  0.0000E+00  0.0000E+00  0.0000E+00 -1.0368E-08
             1.8189E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -426.959567812839        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:      990
 NPARAMETR:  4.8990E-01  1.1302E-02  2.9781E-02  3.5000E-01  3.5945E+06  1.5436E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.3875E+03
             1.6366E+01
 PARAMETER: -6.1356E-01 -4.3824E+00 -3.4139E+00 -9.4982E-01  1.5210E+01  5.3410E-01 -1.1808E+01 -5.1336E+00 -5.0830E+00  7.8836E+00
             2.8952E+00
 GRADIENT:  -2.9113E-03  3.3148E-05  4.3351E-03 -1.7832E-03  7.9012E-07 -1.1597E-03  0.0000E+00  0.0000E+00  0.0000E+00  1.8684E-06
            -2.9252E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      990
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0685E-03  1.6956E-06  1.0359E-04 -2.1728E-04  2.3389E-07
 SE:             2.8277E-02  1.0176E-06  2.0770E-04  2.6111E-04  6.5477E-07
 N:                     100         100         100         100         100

 P VAL.:         9.1359E-01  9.5668E-02  6.1795E-01  4.0532E-01  7.2093E-01

 ETASHRINKSD(%)  5.2696E+00  9.9997E+01  9.9304E+01  9.9125E+01  9.9998E+01
 ETASHRINKVR(%)  1.0261E+01  1.0000E+02  9.9995E+01  9.9992E+01  1.0000E+02
 EBVSHRINKSD(%)  5.4802E+00  9.9996E+01  9.9226E+01  9.9020E+01  9.9998E+01
 EBVSHRINKVR(%)  1.0660E+01  1.0000E+02  9.9994E+01  9.9990E+01  1.0000E+02
 RELATIVEINF(%)  1.1800E-01  1.8777E-09  2.3687E-05  1.9699E-05  1.3255E-10
 EPSSHRINKSD(%)  4.0314E+00
 EPSSHRINKVR(%)  7.9002E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -426.95956781283945     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       308.19125875089873     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.35
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -426.960       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.90E-01  1.13E-02  2.98E-02  3.50E-01  3.65E+06  1.54E+00  1.00E-02  1.00E-02  1.00E-02  2.40E+03  1.64E+01
 


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
+        1.69E+03
 
 TH 2
+       -4.83E+02  1.26E+03
 
 TH 3
+       -7.80E+03  4.67E+02  4.51E+05
 
 TH 4
+       -2.81E+02  2.36E+02 -4.64E+04  5.58E+03
 
 TH 5
+       -2.98E-11  4.20E-09  1.87E-09  5.64E-10  5.08E-20
 
 TH 6
+       -1.65E+01  2.27E+01  4.54E+02 -5.76E+01 -1.18E-10  5.28E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.30E-07 -8.02E-06 -4.79E-06 -3.93E-07  8.31E-15 -5.46E-07  0.00E+00  0.00E+00  0.00E+00  1.00E-11
 
 TH11
+       -1.92E+01  4.96E+00  2.21E+02 -1.29E+01 -2.13E-12  4.30E-01  0.00E+00  0.00E+00  0.00E+00 -2.79E-09  1.25E+00
 
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
 #CPUT: Total CPU Time in Seconds,       19.175
Stop Time:
Sat Sep 25 14:52:35 CDT 2021
