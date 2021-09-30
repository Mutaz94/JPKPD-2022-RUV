Wed Sep 29 21:56:45 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat25.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1665.88627004855        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3113E+02  4.9659E+00 -7.0839E+00  6.2578E+00  5.0745E+01  6.0048E+01 -3.6081E+01  1.6443E+01 -4.3042E+01 -2.2858E+01
            -8.1029E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1851.72258278188        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      150
 NPARAMETR:  1.0051E+00  9.4043E-01  1.1223E+00  1.1365E+00  1.0087E+00  8.6317E-01  1.3039E+00  8.2288E-01  1.3474E+00  1.1725E+00
             1.9004E+00
 PARAMETER:  1.0513E-01  3.8579E-02  2.1534E-01  2.2799E-01  1.0871E-01 -4.7145E-02  3.6539E-01 -9.4948E-02  3.9815E-01  2.5913E-01
             7.4204E-01
 GRADIENT:   2.0050E+01  1.2996E+01 -1.4147E+01  3.6212E+01 -3.4446E+00 -2.2117E+01  4.4110E+00  8.3336E+00  4.2893E+01  4.0352E+00
             9.2437E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1858.49495891253        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  1.0107E+00  6.2044E-01  1.3557E+00  1.3757E+00  1.0036E+00  8.5368E-01  1.7894E+00  4.6587E-01  1.0159E+00  1.3082E+00
             1.8834E+00
 PARAMETER:  1.1060E-01 -3.7732E-01  4.0435E-01  4.1896E-01  1.0356E-01 -5.8198E-02  6.8188E-01 -6.6386E-01  1.1576E-01  3.6868E-01
             7.3307E-01
 GRADIENT:   4.8543E+01  2.7985E+01  1.4176E+00  8.1202E+01  9.0309E-01 -2.5811E+01  3.4299E+00 -8.4198E-02  5.3872E+00  1.5826E+01
             7.7418E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1872.04077747638        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  9.9198E-01  4.8977E-01  9.0419E-01  1.3659E+00  7.2241E-01  9.2110E-01  1.7714E+00  1.6264E-01  9.5766E-01  9.4100E-01
             1.6679E+00
 PARAMETER:  9.1953E-02 -6.1382E-01 -7.1237E-04  4.1178E-01 -2.2516E-01  1.7813E-02  6.7175E-01 -1.7162E+00  5.6733E-02  3.9193E-02
             6.1159E-01
 GRADIENT:   2.8470E-01  1.4614E+01  7.8041E-02  3.4846E+01 -3.6763E+00  3.6468E+00 -1.8985E+00  4.1437E-01 -7.3119E+00 -6.2081E-01
            -1.8088E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1876.33801256683        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  9.8879E-01  1.8296E-01  8.5723E-01  1.5092E+00  6.2322E-01  8.9685E-01  2.7882E+00  8.5690E-02  9.1763E-01  8.8255E-01
             1.7300E+00
 PARAMETER:  8.8728E-02 -1.5985E+00 -5.4046E-02  5.1157E-01 -3.7286E-01 -8.8612E-03  1.1254E+00 -2.3570E+00  1.4043E-02 -2.4936E-02
             6.4811E-01
 GRADIENT:   9.0151E+00  2.6510E+00 -1.7162E+00  1.4118E+01 -6.0733E-01 -3.9860E+00 -5.7114E-01  9.3271E-02 -5.4101E-01 -4.4459E-02
             7.6856E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1877.44385124040        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      856
 NPARAMETR:  9.8028E-01  3.9317E-02  8.8163E-01  1.5870E+00  6.0551E-01  9.0404E-01  6.2655E+00  3.5360E-02  8.8921E-01  8.9247E-01
             1.7229E+00
 PARAMETER:  8.0080E-02 -3.1361E+00 -2.5980E-02  5.6185E-01 -4.0168E-01 -8.7683E-04  1.9351E+00 -3.2422E+00 -1.7417E-02 -1.3757E-02
             6.4403E-01
 GRADIENT:  -2.3561E+00  7.0038E-01  3.7938E+00  1.1666E+01 -7.2772E+00  2.1626E-01  2.2496E-01  8.1389E-03 -4.9958E-01  5.6645E-01
             6.3378E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1877.65985223793        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1032
 NPARAMETR:  9.8008E-01  1.0000E-02  8.8242E-01  1.5948E+00  6.0290E-01  9.0300E-01  1.2341E+01  1.2273E-02  8.8337E-01  8.9054E-01
             1.7226E+00
 PARAMETER:  7.9881E-02 -4.5422E+00 -2.5090E-02  5.6677E-01 -4.0600E-01 -2.0298E-03  2.6129E+00 -4.3004E+00 -2.4017E-02 -1.5930E-02
             6.4386E-01
 GRADIENT:  -1.1403E-01  2.1896E-02 -5.2657E-01 -1.2012E+00  6.9875E-01  7.1232E-02  6.7901E-02  8.6398E-04 -5.5124E-01  3.5656E-03
            -2.7045E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1877.66492492235        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1226             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8026E-01  1.0000E-02  8.8247E-01  1.5921E+00  6.0237E-01  9.0283E-01  1.1991E+01  1.0000E-02  8.8498E-01  8.9030E-01
             1.7229E+00
 PARAMETER:  8.0059E-02 -4.5649E+00 -2.5033E-02  5.6504E-01 -4.0689E-01 -2.2210E-03  2.5842E+00 -5.0979E+00 -2.2190E-02 -1.6201E-02
             6.4403E-01
 GRADIENT:   1.3102E+02  0.0000E+00  2.7170E+00  3.5110E+02  2.6654E+01  9.7932E+00  3.6279E-01  0.0000E+00  8.6726E+00  1.9324E+00
             6.3970E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1877.66533929227        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1359
 NPARAMETR:  9.8026E-01  1.0000E-02  8.8159E-01  1.5921E+00  6.0252E-01  9.0283E-01  1.1990E+01  1.0000E-02  8.8487E-01  8.9019E-01
             1.7230E+00
 PARAMETER:  8.0066E-02 -4.5649E+00 -2.6033E-02  5.6506E-01 -4.0663E-01 -2.2254E-03  2.5841E+00 -5.0979E+00 -2.2312E-02 -1.6326E-02
             6.4405E-01
 GRADIENT:   4.5209E-01  0.0000E+00 -2.5553E-01 -5.4833E+00  9.5644E-01  1.8175E-02  5.7083E-03  0.0000E+00  4.3116E-02  3.2785E-02
            -5.0098E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1359
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.8713E-04  8.5100E-04 -6.7836E-05 -5.0787E-03 -1.7240E-02
 SE:             2.9553E-02  1.8912E-03  1.9754E-04  2.8882E-02  2.3501E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7875E-01  6.5272E-01  7.3129E-01  8.6042E-01  4.6320E-01

 ETASHRINKSD(%)  9.9381E-01  9.3664E+01  9.9338E+01  3.2432E+00  2.1268E+01
 ETASHRINKVR(%)  1.9777E+00  9.9599E+01  9.9996E+01  6.3813E+00  3.8012E+01
 EBVSHRINKSD(%)  1.1487E+00  9.4164E+01  9.9322E+01  2.9458E+00  2.0081E+01
 EBVSHRINKVR(%)  2.2843E+00  9.9659E+01  9.9995E+01  5.8048E+00  3.6130E+01
 RELATIVEINF(%)  8.8841E+01  1.7626E-02  4.0319E-04  7.5361E+00  4.2378E+00
 EPSSHRINKSD(%)  2.9014E+01
 EPSSHRINKVR(%)  4.9609E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1877.6653392922731     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -958.72680608760038     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.17
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1877.665       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.00E-02  8.82E-01  1.59E+00  6.03E-01  9.03E-01  1.20E+01  1.00E-02  8.85E-01  8.90E-01  1.72E+00
 


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
+        1.39E+03
 
 TH 2
+        0.00E+00  1.05E+03
 
 TH 3
+       -1.65E+01  0.00E+00  6.66E+02
 
 TH 4
+       -1.21E+01  0.00E+00 -8.40E+01  5.29E+02
 
 TH 5
+        2.80E+01  0.00E+00 -1.18E+03 -5.32E+01  2.49E+03
 
 TH 6
+        2.79E+00  0.00E+00  1.11E+00 -4.33E+00 -2.00E+00  2.32E+02
 
 TH 7
+        3.58E-03  0.00E+00 -6.61E-03 -2.06E-02  1.61E-02  1.72E-03  6.71E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.67E+00  0.00E+00  1.38E+01 -1.43E+00 -2.47E+00 -5.03E-01  1.78E-02  0.00E+00  2.24E+02
 
 TH10
+        6.24E-01  0.00E+00 -2.47E+01  2.33E+00 -8.44E+01  1.03E+00 -1.29E-02  0.00E+00 -1.40E+00  1.17E+02
 
 TH11
+       -1.46E+01  0.00E+00 -8.37E+00 -9.09E+00  1.37E+01  2.94E+00 -1.82E-03  0.00E+00  7.03E+00  1.96E+01  1.53E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.213
Stop Time:
Wed Sep 29 21:57:18 CDT 2021
