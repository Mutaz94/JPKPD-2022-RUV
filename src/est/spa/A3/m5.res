Sat Sep 25 09:02:16 CDT 2021
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
$DATA ../../../../data/spa/A3/dat5.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -508.300513632428        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1016E+02  4.0667E+01  1.0786E+02 -7.4639E+01  4.7806E+01  1.7096E+00 -5.8092E+01 -2.6543E+01 -8.3800E+01 -1.4757E+02
            -1.8948E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1194.08401571502        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.7143E-01  9.7901E-01  8.7987E-01  1.1972E+00  9.3728E-01  8.7420E-01  1.0870E+00  9.9903E-01  1.0973E+00  1.2513E+00
             5.3119E+00
 PARAMETER:  7.1012E-02  7.8788E-02 -2.7978E-02  2.8002E-01  3.5222E-02 -3.4444E-02  1.8340E-01  9.9025E-02  1.9288E-01  3.2419E-01
             1.7699E+00
 GRADIENT:  -5.5877E+01  7.0100E+00 -1.7802E+01  3.1062E+01 -2.1953E+01 -5.2696E+00  1.3250E+01  7.1783E+00  3.3023E+01  3.1597E+01
             2.5098E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1243.58772932920        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.4736E-01  6.9075E-01  1.7114E+00  1.4233E+00  1.3256E+00  8.8184E-01  2.0099E+00  2.4083E-01  8.4247E-01  1.7902E+00
             4.1095E+00
 PARAMETER:  4.5926E-02 -2.6998E-01  6.3731E-01  4.5299E-01  3.8186E-01 -2.5743E-02  7.9806E-01 -1.3237E+00 -7.1417E-02  6.8231E-01
             1.5133E+00
 GRADIENT:  -4.5231E+01  3.0606E+01 -2.5600E+00  5.7731E+01 -2.2596E+01 -1.0675E+01  2.0522E+01  1.0125E-01  2.4242E+01  3.5945E+01
             1.2026E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1279.29191414120        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.5975E-01  7.3384E-01  1.1571E+00  1.2485E+00  1.0188E+00  9.0913E-01  1.0269E+00  9.0871E-02  8.0274E-01  5.6881E-01
             3.4902E+00
 PARAMETER:  5.8917E-02 -2.0947E-01  2.4595E-01  3.2194E-01  1.1862E-01  4.7365E-03  1.2650E-01 -2.2983E+00 -1.1972E-01 -4.6421E-01
             1.3499E+00
 GRADIENT:   2.9640E+01 -1.0981E+01 -6.4858E+00 -2.2380E+01  8.8224E+00 -1.3480E+00  9.8720E-01  3.6237E-02  3.1924E-01  3.6250E+00
            -1.3765E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1281.09908745828        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      295
 NPARAMETR:  9.4995E-01  9.4959E-01  1.3871E+00  1.1266E+00  1.2217E+00  9.1460E-01  4.3904E-01  3.0093E-02  9.3768E-01  1.7442E-01
             3.5605E+00
 PARAMETER:  4.8652E-02  4.8273E-02  4.2721E-01  2.1924E-01  3.0022E-01  1.0732E-02 -7.2316E-01 -3.4035E+00  3.5653E-02 -1.6463E+00
             1.3699E+00
 GRADIENT:   1.1987E+00 -3.4517E+00 -2.0529E-01 -5.7046E+00  1.3196E+00  9.7213E-01  1.4670E-01  2.5152E-03  4.8387E-01  1.5788E-01
             1.3360E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1281.14056192976        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  9.5013E-01  1.0028E+00  1.4152E+00  1.0937E+00  1.2606E+00  9.1328E-01  2.7325E-01  1.2492E-02  9.8220E-01  7.0861E-02
             3.5663E+00
 PARAMETER:  4.8843E-02  1.0284E-01  4.4728E-01  1.8959E-01  3.3162E-01  9.2862E-03 -1.1974E+00 -4.2827E+00  8.2042E-02 -2.5470E+00
             1.3715E+00
 GRADIENT:   5.0016E-01 -1.2455E+00 -1.5288E-01 -2.7079E+00  8.1770E-01  2.6964E-01  2.9631E-02  4.0356E-04  1.1690E-01  2.3688E-02
             6.1676E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1281.25757584045        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  9.4911E-01  8.5939E-01  1.4093E+00  1.1824E+00  1.1921E+00  9.1272E-01  1.2843E-01  1.0000E-02  9.2822E-01  1.0000E-02
             3.5594E+00
 PARAMETER:  4.7765E-02 -5.1529E-02  4.4308E-01  2.6758E-01  2.7569E-01  8.6710E-03 -1.9523E+00 -6.0127E+00  2.5509E-02 -4.7581E+00
             1.3696E+00
 GRADIENT:   1.2163E+00 -2.1418E+00 -2.5283E-01 -4.1862E+00  7.7141E-01  3.6308E-01  1.4841E-02  0.0000E+00  1.5613E-01  0.0000E+00
             2.3416E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1281.31672215536        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      596
 NPARAMETR:  9.4942E-01  7.5789E-01  1.4241E+00  1.2535E+00  1.1557E+00  9.1241E-01  6.5642E-02  1.0000E-02  8.8190E-01  1.0000E-02
             3.5643E+00
 PARAMETER:  4.8100E-02 -1.7722E-01  4.5356E-01  3.2595E-01  2.4472E-01  8.3290E-03 -2.6235E+00 -7.5264E+00 -2.5680E-02 -6.6400E+00
             1.3710E+00
 GRADIENT:  -7.1560E-02  1.0539E+00 -9.8856E-02  2.2354E+00 -6.7617E-02 -5.5135E-04  2.5408E-03  0.0000E+00  6.5160E-02  0.0000E+00
             4.1286E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1281.33187167150        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      771
 NPARAMETR:  9.4843E-01  6.7260E-01  1.4246E+00  1.3052E+00  1.1197E+00  9.1168E-01  2.7322E-02  1.0000E-02  8.4784E-01  1.0000E-02
             3.5587E+00
 PARAMETER:  4.7055E-02 -2.9661E-01  4.5391E-01  3.6636E-01  2.1305E-01  7.5392E-03 -3.5001E+00 -9.4892E+00 -6.5061E-02 -8.9229E+00
             1.3694E+00
 GRADIENT:   2.2942E-02  2.3495E-02 -8.5091E-03  5.2722E-02  1.0481E-02  2.3763E-03  4.8867E-04  0.0000E+00 -6.5813E-03  0.0000E+00
            -3.9018E-04

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1281.33191399000        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  9.4836E-01  6.6742E-01  1.4243E+00  1.3084E+00  1.1174E+00  9.1164E-01  2.5578E-02  1.0000E-02  8.4595E-01  1.0000E-02
             3.5585E+00
 PARAMETER:  4.6983E-02 -3.0434E-01  4.5368E-01  3.6878E-01  2.1100E-01  7.4863E-03 -3.5660E+00 -9.6358E+00 -6.7290E-02 -9.0907E+00
             1.3693E+00
 GRADIENT:   9.1906E-04 -8.7731E-03  1.1301E-03 -1.0270E-02  1.6675E-03  1.7190E-04  4.2987E-04  0.0000E+00  4.8727E-03  0.0000E+00
             4.5116E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      898
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0506E-03 -5.9227E-04  8.1073E-05 -9.8834E-03  9.9740E-06
 SE:             2.8475E-02  2.8820E-04  6.9803E-05  2.3607E-02  1.3407E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4259E-01  3.9871E-02  2.4546E-01  6.7546E-01  9.4070E-01

 ETASHRINKSD(%)  4.6061E+00  9.9035E+01  9.9766E+01  2.0913E+01  9.9551E+01
 ETASHRINKVR(%)  9.0001E+00  9.9991E+01  9.9999E+01  3.7453E+01  9.9998E+01
 EBVSHRINKSD(%)  4.5062E+00  9.9078E+01  9.9735E+01  2.0438E+01  9.9527E+01
 EBVSHRINKVR(%)  8.8093E+00  9.9991E+01  9.9999E+01  3.6699E+01  9.9998E+01
 RELATIVEINF(%)  8.7772E+01  2.4048E-04  6.8785E-05  2.5235E+00  1.2850E-04
 EPSSHRINKSD(%)  1.9601E+01
 EPSSHRINKVR(%)  3.5361E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1281.3319139900002     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -546.18108742626202     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1281.332       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.48E-01  6.67E-01  1.42E+00  1.31E+00  1.12E+00  9.12E-01  2.56E-02  1.00E-02  8.46E-01  1.00E-02  3.56E+00
 


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
+        1.36E+03
 
 TH 2
+       -8.94E+01  3.33E+02
 
 TH 3
+        4.54E+00  4.16E+01  3.25E+01
 
 TH 4
+       -1.08E+02  3.99E+02  2.05E+01  5.65E+02
 
 TH 5
+        8.53E+00 -1.54E+02 -7.49E+01 -1.08E+02  1.98E+02
 
 TH 6
+        5.12E+00 -1.75E+01  2.89E+00 -2.47E+01 -2.11E+00  1.96E+02
 
 TH 7
+       -1.12E+00  1.45E+00 -2.96E-01 -5.08E-01 -4.24E-01 -1.74E+00 -3.91E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.13E+01 -4.66E+01  2.20E+00 -3.40E+00  1.62E+01  6.76E+00 -2.15E+00  0.00E+00  9.65E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.99E+01 -1.65E+01  6.09E-02 -1.33E+01  3.58E+00  4.09E+00  1.33E-01  0.00E+00  1.45E+01  0.00E+00  3.53E+01
 
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
 #CPUT: Total CPU Time in Seconds,       15.301
Stop Time:
Sat Sep 25 09:02:33 CDT 2021
