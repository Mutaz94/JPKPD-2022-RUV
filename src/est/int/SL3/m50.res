Sat Sep 25 02:21:42 CDT 2021
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
$DATA ../../../../data/int/SL3/dat50.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      982
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

 TOT. NO. OF OBS RECS:      882
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -621.549746854606        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6383E+01  6.3358E+01  4.8871E+01  5.3492E+00  1.2582E+02  5.5298E-01 -9.7230E+01 -1.5923E+02 -1.2379E+02 -8.8967E-01
            -5.9687E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2722.76225388509        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0112E+00  1.1290E+00  1.1911E+00  9.9456E-01  1.0933E+00  9.6821E-01  8.1143E-01  9.2419E-01  1.0594E+00  7.6802E-01
             2.6535E+00
 PARAMETER:  1.1109E-01  2.2131E-01  2.7489E-01  9.4544E-02  1.8923E-01  6.7691E-02 -1.0895E-01  2.1163E-02  1.5769E-01 -1.6394E-01
             1.0759E+00
 GRADIENT:  -5.5370E+00  1.1132E+01  5.6387E+00  2.4124E+01  1.5773E+01 -8.5946E+00 -2.0641E+00 -3.7941E-01 -5.5080E+00 -3.0484E+01
            -1.6458E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2726.47460398363        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0245E+00  1.2550E+00  1.3652E+00  9.4601E-01  1.2608E+00  1.0405E+00  5.8712E-01  6.5694E-01  1.1487E+00  1.0542E+00
             2.6775E+00
 PARAMETER:  1.2420E-01  3.2715E-01  4.1127E-01  4.4498E-02  3.3176E-01  1.3970E-01 -4.3253E-01 -3.2016E-01  2.3867E-01  1.5278E-01
             1.0849E+00
 GRADIENT:   2.2836E+01  2.0675E+01 -4.8295E+00  6.0477E+01  2.9787E+01  1.8080E+01 -6.1469E+00  6.2538E-01 -6.5768E+00 -9.0502E+00
             5.1664E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2730.68731634444        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0129E+00  1.3728E+00  1.2718E+00  8.3124E-01  1.2985E+00  9.9628E-01  6.4752E-01  3.4569E-01  1.2537E+00  1.1202E+00
             2.6622E+00
 PARAMETER:  1.1285E-01  4.1687E-01  3.4046E-01 -8.4833E-02  3.6122E-01  9.6273E-02 -3.3460E-01 -9.6221E-01  3.2608E-01  2.1347E-01
             1.0792E+00
 GRADIENT:   6.4150E-01  4.9778E+00 -9.0316E-01  5.1872E+00 -7.8792E-02  2.5830E+00  1.3820E+00  2.6646E-01  1.3330E+00 -4.9549E-01
             3.3002E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2731.02359776453        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0126E+00  1.4114E+00  1.3028E+00  8.0649E-01  1.3382E+00  9.9517E-01  6.2532E-01  2.1511E-01  1.2897E+00  1.1453E+00
             2.6619E+00
 PARAMETER:  1.1249E-01  4.4460E-01  3.6449E-01 -1.1506E-01  3.9133E-01  9.5154E-02 -3.6950E-01 -1.4366E+00  3.5444E-01  2.3571E-01
             1.0790E+00
 GRADIENT:  -9.2908E-02  5.3344E+00 -8.0295E-01  4.9401E+00  1.4251E+00  2.1814E+00  7.4665E-01  1.0555E-01  1.8061E+00 -7.8099E-01
             3.2357E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2731.46244383987        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      438
 NPARAMETR:  1.0165E+00  1.5464E+00  1.4152E+00  7.2571E-01  1.4734E+00  9.9012E-01  5.8974E-01  1.0962E-01  1.3849E+00  1.2321E+00
             2.6569E+00
 PARAMETER:  1.1641E-01  5.3595E-01  4.4724E-01 -2.2061E-01  4.8758E-01  9.0073E-02 -4.2807E-01 -2.1107E+00  4.2561E-01  3.0868E-01
             1.0771E+00
 GRADIENT:   1.1455E+00  8.8203E+00 -6.0583E-01  7.8445E+00  1.6841E+00 -3.5387E-01 -5.9274E-01  2.7329E-02 -1.1288E+00  3.9885E-01
            -2.4892E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2731.53223885682        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      614
 NPARAMETR:  1.0159E+00  1.5767E+00  1.4216E+00  6.9841E-01  1.4975E+00  9.9080E-01  5.8643E-01  5.3311E-02  1.4253E+00  1.2443E+00
             2.6582E+00
 PARAMETER:  1.1579E-01  5.5533E-01  4.5179E-01 -2.5895E-01  5.0382E-01  9.0756E-02 -4.3371E-01 -2.8316E+00  4.5438E-01  3.1855E-01
             1.0776E+00
 GRADIENT:   8.9689E-02 -4.2117E-01 -1.1210E-01 -2.4642E-01  2.7844E-01 -3.2739E-02  2.2897E-02  6.4868E-03 -2.7960E-02  4.3708E-02
            -2.2325E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2731.53561414464        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      790
 NPARAMETR:  1.0159E+00  1.5740E+00  1.4226E+00  7.0035E-01  1.4954E+00  9.9087E-01  5.8660E-01  1.0000E-02  1.4228E+00  1.2431E+00
             2.6582E+00
 PARAMETER:  1.1576E-01  5.5365E-01  4.5248E-01 -2.5618E-01  5.0238E-01  9.0823E-02 -4.3341E-01 -4.6329E+00  4.5260E-01  3.1761E-01
             1.0776E+00
 GRADIENT:   1.7241E-02 -3.8531E-02 -7.4778E-03 -3.7721E-02  3.1590E-02 -7.0310E-03  5.7825E-03  0.0000E+00 -3.2360E-02  5.7566E-03
            -9.5789E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2731.53561414464        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  1.0159E+00  1.5740E+00  1.4226E+00  7.0035E-01  1.4954E+00  9.9087E-01  5.8660E-01  1.0000E-02  1.4228E+00  1.2431E+00
             2.6582E+00
 PARAMETER:  1.1576E-01  5.5365E-01  4.5248E-01 -2.5618E-01  5.0238E-01  9.0823E-02 -4.3341E-01 -4.6329E+00  4.5260E-01  3.1761E-01
             1.0776E+00
 GRADIENT:   1.7241E-02 -3.8531E-02 -7.4778E-03 -3.7721E-02  3.1590E-02 -7.0310E-03  5.7825E-03  0.0000E+00 -3.2360E-02  5.7566E-03
            -9.5789E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      812
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8280E-03 -3.3257E-02 -7.4512E-05  1.5283E-02 -1.8859E-02
 SE:             2.9400E-02  1.6770E-02  6.2753E-05  2.4703E-02  2.5659E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5042E-01  4.7350E-02  2.3508E-01  5.3613E-01  4.6235E-01

 ETASHRINKSD(%)  1.5075E+00  4.3820E+01  9.9790E+01  1.7242E+01  1.4040E+01
 ETASHRINKVR(%)  2.9923E+00  6.8438E+01  1.0000E+02  3.1512E+01  2.6109E+01
 EBVSHRINKSD(%)  1.6397E+00  4.4071E+01  9.9779E+01  1.7503E+01  1.3228E+01
 EBVSHRINKVR(%)  3.2526E+00  6.8719E+01  1.0000E+02  3.1943E+01  2.4707E+01
 RELATIVEINF(%)  9.6669E+01  4.0011E+00  2.0023E-04  1.0892E+01  1.7951E+01
 EPSSHRINKSD(%)  1.5950E+01
 EPSSHRINKVR(%)  2.9356E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          882
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1621.0075725730426     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2731.5356141446350     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1110.5280415715924     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2731.536       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.57E+00  1.42E+00  7.00E-01  1.50E+00  9.91E-01  5.87E-01  1.00E-02  1.42E+00  1.24E+00  2.66E+00
 


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
+        1.06E+03
 
 TH 2
+       -1.35E+01  4.17E+02
 
 TH 3
+        1.24E+00  3.27E+01  2.38E+01
 
 TH 4
+       -1.96E+01  4.16E+02 -2.19E+01  7.61E+02
 
 TH 5
+       -2.66E+00 -1.48E+02 -5.91E+01  4.80E+01  2.41E+02
 
 TH 6
+        5.19E+00 -3.42E+00  1.06E+00 -3.31E+00 -1.39E+00  1.83E+02
 
 TH 7
+        2.17E+00 -2.57E+01  2.23E+00 -1.98E+01 -5.33E+00 -1.02E+00  7.81E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.23E+00 -1.90E+01 -4.54E-01  3.49E+01  1.30E+00 -1.16E+00  2.49E+01  0.00E+00  4.67E+01
 
 TH10
+        1.75E+00 -1.51E+01 -9.51E+00  7.99E+00 -1.37E+01  6.73E-01  7.07E+00  0.00E+00  8.31E-01  7.45E+01
 
 TH11
+       -1.55E+01 -1.84E+01 -1.29E+00 -1.70E+01  2.98E+00  2.96E+00  8.04E+00  0.00E+00  4.92E+00  7.13E+00  1.63E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.863
Stop Time:
Sat Sep 25 02:22:11 CDT 2021
