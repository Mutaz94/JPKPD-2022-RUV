Wed Sep 29 23:38:54 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1265.69452207474        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5106E+02  3.6329E+01  5.9008E+01  1.8524E+01  1.1599E+02  5.5837E+01 -5.6591E+01 -2.1567E+01 -5.1629E+01 -6.4451E+01
            -1.4772E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1718.68810632448        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0513E+00  9.2012E-01  8.7288E-01  1.1410E+00  8.1616E-01  9.4203E-01  1.4727E+00  7.6499E-01  1.1933E+00  1.2341E+00
             2.0333E+00
 PARAMETER:  1.5002E-01  1.6754E-02 -3.5958E-02  2.3193E-01 -1.0315E-01  4.0283E-02  4.8712E-01 -1.6790E-01  2.7669E-01  3.1034E-01
             8.0966E-01
 GRADIENT:   2.9169E+02  5.1290E+01  1.2026E+01  9.0659E+01 -3.2825E+01  5.4061E+00  1.0476E+01  7.0951E+00  2.9020E+01  2.7431E+01
            -3.3088E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1731.56435470722        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0256E+00  6.1336E-01  5.1889E-01  1.3092E+00  4.9064E-01  1.0525E+00  1.9576E+00  1.4749E-01  1.0527E+00  7.4444E-01
             2.0019E+00
 PARAMETER:  1.2523E-01 -3.8881E-01 -5.5607E-01  3.6938E-01 -6.1204E-01  1.5115E-01  7.7174E-01 -1.8140E+00  1.5140E-01 -1.9512E-01
             7.9410E-01
 GRADIENT:   1.8211E+02  8.6107E+01  5.3188E+01  2.1511E+02 -3.4428E+01  5.0150E+01  1.9240E+01  2.3578E-01  9.8386E+00  2.2838E-01
            -3.3354E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1747.95156066411        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      337
 NPARAMETR:  9.8491E-01  4.8303E-01  4.4544E-01  1.2659E+00  4.3030E-01  9.0035E-01  1.5711E+00  2.2958E-02  1.0262E+00  6.2088E-01
             2.0714E+00
 PARAMETER:  8.4791E-02 -6.2768E-01 -7.0870E-01  3.3581E-01 -7.4327E-01 -4.9765E-03  5.5181E-01 -3.6741E+00  1.2584E-01 -3.7662E-01
             8.2822E-01
 GRADIENT:  -2.3107E+01  2.5284E+01  3.4191E+01  4.1078E+01 -4.2887E+01 -1.3237E+01 -8.1746E+00 -1.0997E-03 -3.6423E+00 -1.0998E+01
            -1.0403E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1753.19081502551        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  9.8838E-01  2.7413E-01  4.5578E-01  1.3276E+00  3.9518E-01  9.1638E-01  2.4416E+00  1.0000E-02  9.8724E-01  6.9063E-01
             2.0860E+00
 PARAMETER:  8.8311E-02 -1.1941E+00 -6.8574E-01  3.8340E-01 -8.2842E-01  1.2673E-02  9.9265E-01 -5.2750E+00  8.7158E-02 -2.7015E-01
             8.3527E-01
 GRADIENT:   3.7077E+00  1.1759E+01  2.9451E+01 -4.8401E+00 -5.0107E+01 -2.9495E+00  2.9807E+00  0.0000E+00 -2.6446E+00  8.0875E-01
            -6.1757E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1757.25523516215        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  9.7543E-01  1.2905E-01  5.0203E-01  1.4205E+00  4.1263E-01  9.1906E-01  3.7953E+00  1.0000E-02  9.7376E-01  6.7904E-01
             2.1117E+00
 PARAMETER:  7.5126E-02 -1.9476E+00 -5.8910E-01  4.5100E-01 -7.8520E-01  1.5601E-02  1.4338E+00 -7.4793E+00  7.3413E-02 -2.8708E-01
             8.4748E-01
 GRADIENT:  -1.2032E+01 -1.7778E+00  1.2705E+01  1.8464E+01 -1.3462E+01 -4.0567E-03 -9.0404E+00  0.0000E+00  6.7782E+00 -5.1660E-01
             3.4436E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1758.55971300622        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      872             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7755E-01  7.0551E-02  4.8261E-01  1.4351E+00  3.9647E-01  9.3340E-01  5.3077E+00  1.0000E-02  9.5643E-01  6.7618E-01
             2.1254E+00
 PARAMETER:  7.7298E-02 -2.5514E+00 -6.2854E-01  4.6120E-01 -8.2516E-01  3.1077E-02  1.7692E+00 -9.7932E+00  5.5448E-02 -2.9130E-01
             8.5396E-01
 GRADIENT:   8.0000E+01  4.7667E+00  1.6153E+01  1.8129E+02  4.6395E+01  1.2679E+01  1.9334E+01  0.0000E+00  7.6359E+00  4.7387E+00
             1.7667E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1758.62055928770        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  9.7827E-01  7.0474E-02  4.8278E-01  1.4347E+00  3.9653E-01  9.1766E-01  5.3109E+00  1.0000E-02  9.5649E-01  6.7630E-01
             2.1246E+00
 PARAMETER:  7.8031E-02 -2.5525E+00 -6.2819E-01  4.6096E-01 -8.2500E-01  1.4072E-02  1.7698E+00 -9.7932E+00  5.5516E-02 -2.9113E-01
             8.5360E-01
 GRADIENT:   2.6584E+00  2.7259E-01  3.3281E-01  1.7990E+01 -4.6981E+00 -2.5511E-01 -1.5761E+00  0.0000E+00  1.6021E+00  2.1909E+00
             8.0690E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1758.74919319990        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1187
 NPARAMETR:  9.7774E-01  6.8289E-02  4.8621E-01  1.4275E+00  3.9904E-01  9.1947E-01  5.3557E+00  1.0000E-02  9.5647E-01  6.7688E-01
             2.0998E+00
 PARAMETER:  7.7605E-02 -2.5830E+00 -6.2139E-01  4.5574E-01 -8.1837E-01  1.6413E-02  1.7787E+00 -9.7932E+00  5.6492E-02 -2.8881E-01
             8.4166E-01
 GRADIENT:   1.7126E+00  2.2685E+01 -1.0986E+02 -1.4298E+02  8.2628E+01  7.8170E-01  2.8925E+01  0.0000E+00  3.9941E+00  2.1436E+00
            -1.7962E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1187
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1755E-04  2.5562E-02 -3.5354E-05 -1.4449E-02 -2.6536E-03
 SE:             2.9276E-02  1.1585E-02  2.4339E-04  2.7525E-02  2.1419E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9135E-01  2.7351E-02  8.8451E-01  5.9963E-01  9.0140E-01

 ETASHRINKSD(%)  1.9225E+00  6.1188E+01  9.9185E+01  7.7882E+00  2.8243E+01
 ETASHRINKVR(%)  3.8080E+00  8.4936E+01  9.9993E+01  1.4970E+01  4.8510E+01
 EBVSHRINKSD(%)  1.6559E+00  7.3041E+01  9.9163E+01  5.3939E+00  2.5940E+01
 EBVSHRINKVR(%)  3.2844E+00  9.2732E+01  9.9993E+01  1.0497E+01  4.5151E+01
 RELATIVEINF(%)  9.5670E+01  3.9238E+00  4.4134E-04  4.3741E+01  3.4428E+00
 EPSSHRINKSD(%)  2.7889E+01
 EPSSHRINKVR(%)  4.8000E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1758.7491931998964     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -839.81065999522366     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1758.749       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  6.84E-02  4.86E-01  1.43E+00  3.99E-01  9.20E-01  5.36E+00  1.00E-02  9.57E-01  6.78E-01  2.10E+00
 


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
+        1.34E+03
 
 TH 2
+       -9.33E+01  5.63E+04
 
 TH 3
+        1.03E+01 -1.76E+02  2.24E+04
 
 TH 4
+        2.73E+00  1.67E+02  8.63E+03  4.62E+03
 
 TH 5
+        3.89E+01 -1.82E+01 -5.00E+03 -3.17E+01  2.49E+04
 
 TH 6
+       -6.81E+00  2.24E+00 -1.17E+01 -1.48E+01  9.96E+00  2.22E+02
 
 TH 7
+        9.86E-02  4.24E+01 -1.36E+01 -2.13E+00  2.03E+01 -4.14E-02  2.22E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.66E+01  1.94E+01 -6.17E+01 -5.99E+01  8.35E+01  3.10E+01  3.47E-01  0.00E+00  4.15E+02
 
 TH10
+       -4.72E+01 -1.70E+02 -7.33E+01  4.99E+00 -5.57E+01  2.66E+01 -4.09E+00  0.00E+00  2.01E+02  3.40E+02
 
 TH11
+       -2.72E+01 -1.86E+01 -2.64E+01 -1.36E+01  2.82E+01  8.77E+00 -5.16E-01  0.00E+00  1.03E+04  5.05E+03  6.62E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.05
 #CPUT: Total CPU Time in Seconds,       29.147
Stop Time:
Wed Sep 29 23:39:26 CDT 2021
