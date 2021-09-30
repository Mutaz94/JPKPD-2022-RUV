Wed Sep 29 03:03:22 CDT 2021
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
$DATA ../../../../data/int/SL2/dat34.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      997
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      897
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1622.26562551161        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0493E+02  1.7162E+01  1.9161E+02  4.7402E+01  2.6575E+02 -4.4809E+00 -7.3308E+01 -2.5542E+02 -1.2379E+01 -5.5174E+01
            -4.0176E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2894.20525961057        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      107
 NPARAMETR:  9.1630E-01  1.1695E+00  8.2257E-01  9.8163E-01  8.9861E-01  1.0737E+00  1.0077E+00  1.0802E+00  8.1617E-01  1.0195E+00
             2.6545E+00
 PARAMETER:  1.2588E-02  2.5661E-01 -9.5319E-02  8.1455E-02 -6.9022E-03  1.7110E-01  1.0763E-01  1.7711E-01 -1.0314E-01  1.1926E-01
             1.0763E+00
 GRADIENT:  -9.6568E+01  2.2535E+01 -1.4715E+01 -8.6321E+00 -5.1539E+01  4.2087E+00  9.5181E+00  7.3779E+00  1.0499E+00 -2.5908E+00
             2.9346E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2900.83501116393        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      282
 NPARAMETR:  9.4548E-01  1.6112E+00  1.1138E+00  7.7560E-01  1.2652E+00  1.0085E+00  8.4890E-01  1.3684E+00  8.4160E-01  1.3369E+00
             2.6772E+00
 PARAMETER:  4.3933E-02  5.7700E-01  2.0780E-01 -1.5412E-01  3.3520E-01  1.0843E-01 -6.3812E-02  4.1366E-01 -7.2451E-02  3.9034E-01
             1.0848E+00
 GRADIENT:  -4.6456E+01  9.6934E+01  6.5926E+00  6.5828E+01 -3.2471E+01 -1.5267E+01  3.2675E+00  5.5126E-01  1.0380E+00  9.0467E+00
             3.0631E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2932.27361893556        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      460
 NPARAMETR:  9.4896E-01  1.5292E+00  1.7577E+00  7.7947E-01  1.4435E+00  1.1247E+00  8.3856E-01  2.9832E+00  6.6245E-01  1.1722E+00
             2.2783E+00
 PARAMETER:  4.7613E-02  5.2477E-01  6.6399E-01 -1.4914E-01  4.6707E-01  2.1751E-01 -7.6064E-02  1.1930E+00 -3.1181E-01  2.5889E-01
             9.2342E-01
 GRADIENT:  -1.8660E+01 -8.5391E+00 -1.9090E+01  1.1667E+01  1.8896E+01  2.2923E+01 -1.4730E+01 -9.2042E+00 -5.5561E+00 -2.1318E+01
             2.0726E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2934.65070506748        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      652             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4958E-01  1.5376E+00  1.8207E+00  7.7388E-01  1.4631E+00  1.1159E+00  8.5588E-01  3.1952E+00  6.5728E-01  1.1971E+00
             2.2651E+00
 PARAMETER:  4.8262E-02  5.3020E-01  6.9920E-01 -1.5634E-01  4.8053E-01  2.0963E-01 -5.5630E-02  1.2616E+00 -3.1965E-01  2.7991E-01
             9.1760E-01
 GRADIENT:   5.9225E+01  1.1252E+02 -1.9114E+01  2.5926E+01  5.0910E+01  3.7097E+01 -9.1002E+00 -1.4138E+00 -2.3040E+00 -1.5305E+01
             2.8931E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2935.28949866999        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  9.4934E-01  1.5355E+00  1.8175E+00  7.7358E-01  1.4648E+00  1.0518E+00  8.5566E-01  3.1851E+00  7.3761E-01  1.1963E+00
             2.2702E+00
 PARAMETER:  4.8012E-02  5.2888E-01  6.9747E-01 -1.5673E-01  4.8173E-01  1.5046E-01 -5.5880E-02  1.2585E+00 -2.0435E-01  2.7921E-01
             9.1987E-01
 GRADIENT:  -2.0555E+01 -1.4354E+01 -2.2687E+01  7.2432E+00  1.9873E+01 -6.9124E-01 -4.6075E+00 -3.8638E+00  2.9401E-01 -1.7462E+01
             1.8817E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2938.87009688848        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      954
 NPARAMETR:  9.4924E-01  1.5314E+00  2.3171E+00  7.7095E-01  1.4415E+00  1.0449E+00  8.8424E-01  3.3916E+00  7.1522E-01  1.3009E+00
             2.2802E+00
 PARAMETER:  4.7905E-02  5.2616E-01  9.4031E-01 -1.6013E-01  4.6567E-01  1.4394E-01 -2.3031E-02  1.3213E+00 -2.3516E-01  3.6303E-01
             9.2425E-01
 GRADIENT:   5.5872E+01  1.0687E+02 -8.8432E-01 -6.0251E+00 -4.1791E+00  8.9161E+00  1.2393E+00  2.2753E+00  2.8722E+00  3.0488E+00
             5.1920E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2941.83405490658        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1141             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5917E-01  1.5284E+00  2.5612E+00  7.8503E-01  1.5237E+00  1.0516E+00  9.0761E-01  3.5525E+00  6.4652E-01  1.3284E+00
             2.2841E+00
 PARAMETER:  5.8315E-02  5.2424E-01  1.0405E+00 -1.4203E-01  5.2111E-01  1.5033E-01  3.0582E-03  1.3677E+00 -3.3615E-01  3.8398E-01
             9.2599E-01
 GRADIENT:   7.5314E+01  1.1380E+02 -5.6385E+00  1.5811E+01  3.4031E+01  1.1474E+01  5.4391E-01  5.4223E+00  2.3391E+00  3.6568E+00
             5.6428E+01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2941.86334728631        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1277
 NPARAMETR:  9.5927E-01  1.5283E+00  2.5607E+00  7.8490E-01  1.5254E+00  1.0517E+00  9.1073E-01  3.5511E+00  6.3501E-01  1.3311E+00
             2.2842E+00
 PARAMETER:  5.8383E-02  5.2429E-01  1.0406E+00 -1.4224E-01  5.2214E-01  1.5042E-01  5.9058E-03  1.3677E+00 -3.5060E-01  3.8611E-01
             9.2579E-01
 GRADIENT:  -1.9347E+03  3.6911E+02  1.7380E+02 -1.3563E+03 -1.4828E+00 -8.6735E-02 -1.0357E+00  1.4270E+02  8.0247E-01  4.9870E+02
            -1.7446E+02
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.5         2.8         1.0         2.3         0.8         2.3
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1277
 NO. OF SIG. DIGITS IN FINAL EST.:  0.8
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7248E-03 -8.8586E-03 -1.7648E-02  3.9387E-03 -2.5337E-02
 SE:             2.9481E-02  2.5286E-02  1.7920E-02  1.5125E-02  2.4163E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5335E-01  7.2609E-01  3.2473E-01  7.9454E-01  2.9438E-01

 ETASHRINKSD(%)  1.2339E+00  1.5288E+01  3.9964E+01  4.9331E+01  1.9051E+01
 ETASHRINKVR(%)  2.4527E+00  2.8239E+01  6.3957E+01  7.4326E+01  3.4472E+01
 EBVSHRINKSD(%)  1.2735E+00  1.6099E+01  4.3909E+01  5.1502E+01  1.6561E+01
 EBVSHRINKVR(%)  2.5308E+00  2.9607E+01  6.8537E+01  7.6480E+01  3.0379E+01
 RELATIVEINF(%)  9.7440E+01  7.5504E+00  1.3595E+01  2.3453E+00  3.7033E+01
 EPSSHRINKSD(%)  1.8851E+01
 EPSSHRINKVR(%)  3.4149E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2941.8633472863107     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1293.2876187171280     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.43
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2941.863       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.59E-01  1.53E+00  2.56E+00  7.85E-01  1.53E+00  1.05E+00  9.10E-01  3.55E+00  6.37E-01  1.33E+00  2.28E+00
 


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
+        5.27E+05
 
 TH 2
+       -6.29E+04  7.95E+03
 
 TH 3
+        2.01E+00  4.58E+00  6.29E+02
 
 TH 4
+       -5.70E+01  5.17E+02  4.16E+01  3.88E+05
 
 TH 5
+        6.38E+04 -7.68E+03 -1.83E+02  8.81E+02  1.87E+02
 
 TH 6
+        1.32E+01 -5.47E+00 -1.48E-01 -1.41E+00 -1.93E+00  1.67E+02
 
 TH 7
+       -4.28E+01 -6.63E+04  1.33E+00 -8.03E+01  8.33E-01 -5.73E-01  1.26E+02
 
 TH 8
+        6.91E-01 -1.08E+00  6.07E+00  3.59E+01 -1.66E+01 -1.83E-01  1.12E+00  1.63E+03
 
 TH 9
+       -2.36E+02  1.65E+01  8.75E+00 -2.22E+02  7.63E+00  5.61E-02  3.71E+01  6.69E+00  3.35E+01
 
 TH10
+        9.23E+00  3.48E+00 -2.13E+01  3.56E+02 -1.88E+02 -9.54E-01  6.55E+00 -1.07E+01  4.87E+01  1.83E+04
 
 TH11
+       -1.61E+01 -1.68E+01 -4.39E+01 -8.73E+01  3.62E+01  3.71E+00  6.44E+00 -4.08E+00 -5.93E+00  3.38E+01  1.34E+03
 
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
 #CPUT: Total CPU Time in Seconds,       55.149
Stop Time:
Wed Sep 29 03:04:19 CDT 2021
