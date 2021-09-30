Thu Sep 30 02:26:09 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat87.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1624.51064321008        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0313E+02 -8.3446E+01 -1.2686E+01 -4.8544E+00  1.3234E+02  3.1126E+01 -1.0090E+01 -1.9802E+02 -2.5944E+01 -1.4124E+00
            -7.5199E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2121.89062880570        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.1821E-01  9.8809E-01  1.1354E+00  1.0463E+00  9.8982E-01  1.0383E+00  9.6752E-01  9.4829E-01  1.0284E+00  9.4250E-01
             1.1368E+00
 PARAMETER:  1.4670E-02  8.8017E-02  2.2695E-01  1.4528E-01  8.9768E-02  1.3762E-01  6.6981E-02  4.6906E-02  1.2799E-01  4.0775E-02
             2.2819E-01
 GRADIENT:   1.4804E+02 -9.7349E+00 -2.9717E+01  7.9935E+01  6.6630E+01  7.1453E+01  1.3123E+00  4.3619E+00  1.2536E+01 -8.3639E+00
             8.9125E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2126.56693629024        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      264
 NPARAMETR:  9.5702E-01  8.4397E-01  1.3306E+00  1.1693E+00  9.9855E-01  9.2731E-01  1.2977E+00  9.1160E-01  9.2404E-01  1.0927E+00
             1.1174E+00
 PARAMETER:  5.6070E-02 -6.9636E-02  3.8563E-01  2.5645E-01  9.8548E-02  2.4529E-02  3.6057E-01  7.4423E-03  2.1002E-02  1.8869E-01
             2.1103E-01
 GRADIENT:  -1.4952E+02 -2.0886E+00 -1.5652E+01  4.9582E+00  2.7258E+01 -6.1276E+01  9.4821E+00 -4.9989E-01  1.7229E+00  8.6218E+00
             7.6473E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2138.44674563196        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  9.9546E-01  8.6114E-01  1.2048E+00  1.1441E+00  9.3290E-01  1.0160E+00  1.1685E+00  8.0733E-01  9.4136E-01  9.8900E-01
             1.0286E+00
 PARAMETER:  9.5447E-02 -4.9502E-02  2.8631E-01  2.3464E-01  3.0542E-02  1.1590E-01  2.5570E-01 -1.1403E-01  3.9569E-02  8.8935E-02
             1.2819E-01
 GRADIENT:  -4.0033E+01 -4.3952E+00 -9.6753E-01 -7.1632E+00 -1.7803E-01 -1.2243E+01  3.7336E+00 -1.6713E+00 -1.2946E+00 -6.6658E-01
             1.3179E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2139.99908759878        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  1.0140E+00  9.6411E-01  1.3014E+00  1.0898E+00  1.0026E+00  1.0494E+00  9.5456E-01  1.0103E+00  1.0084E+00  1.0421E+00
             1.0075E+00
 PARAMETER:  1.1388E-01  6.3447E-02  3.6345E-01  1.8596E-01  1.0255E-01  1.4823E-01  5.3497E-02  1.1028E-01  1.0837E-01  1.4124E-01
             1.0744E-01
 GRADIENT:  -9.1404E-01 -1.3885E+00  4.5615E-01  1.5703E-01  8.1604E-01  1.4557E+00 -4.5540E-01 -3.3031E-01 -2.3042E-02 -7.5632E-01
            -4.8683E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2140.03777923506        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:      811             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0154E+00  9.8377E-01  1.2857E+00  1.0755E+00  1.0038E+00  1.0471E+00  9.5397E-01  1.0094E+00  1.0187E+00  1.0422E+00
             1.0083E+00
 PARAMETER:  1.1529E-01  8.3632E-02  3.5134E-01  1.7280E-01  1.0380E-01  1.4603E-01  5.2872E-02  1.0933E-01  1.1851E-01  1.4134E-01
             1.0822E-01
 GRADIENT:   5.2658E+02  3.8114E+01  5.9465E+00  1.5859E+02  6.8730E+00  9.7225E+01  3.0526E+00  1.8382E-01  1.2158E+01  8.8713E-01
             1.4315E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2140.04012015016        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      990             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0160E+00  9.8426E-01  1.2811E+00  1.0755E+00  1.0033E+00  1.0483E+00  9.5609E-01  1.0019E+00  1.0193E+00  1.0415E+00
             1.0078E+00
 PARAMETER:  1.1587E-01  8.4132E-02  3.4771E-01  1.7279E-01  1.0328E-01  1.4717E-01  5.5095E-02  1.0190E-01  1.1910E-01  1.4064E-01
             1.0781E-01
 GRADIENT:   5.3074E+02  3.8205E+01  5.4793E+00  1.5954E+02  7.8810E+00  9.8338E+01  3.1410E+00  1.0965E-01  1.2431E+01  8.1114E-01
             1.0500E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2140.04309507782        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1151
 NPARAMETR:  1.0158E+00  9.8773E-01  1.2801E+00  1.0764E+00  1.0027E+00  1.0475E+00  9.5543E-01  1.0017E+00  1.0150E+00  1.0439E+00
             1.0078E+00
 PARAMETER:  1.1572E-01  8.7653E-02  3.4696E-01  1.7365E-01  1.0271E-01  1.4643E-01  5.4402E-02  1.0170E-01  1.1493E-01  1.4295E-01
             1.0781E-01
 GRADIENT:   2.4035E+00  1.1951E+00  8.1759E-01  3.1944E+00 -1.6277E+00  6.3338E-01 -1.9496E-01 -6.8205E-02 -1.0497E-01  3.4683E-02
            -7.1566E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2140.04484235364        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1322
 NPARAMETR:  1.0159E+00  9.8636E-01  1.2781E+00  1.0772E+00  1.0033E+00  1.0475E+00  9.6008E-01  1.0012E+00  1.0142E+00  1.0434E+00
             1.0079E+00
 PARAMETER:  1.1560E-01  8.7074E-02  3.4600E-01  1.7365E-01  1.0302E-01  1.4620E-01  5.9502E-02  1.0109E-01  1.1528E-01  1.4234E-01
             1.0784E-01
 GRADIENT:  -3.0663E-01  3.3809E-01  2.9349E-01 -1.4756E+05 -4.6172E-01 -1.0297E-01  1.2125E-02 -4.6523E-03  1.6894E-01 -3.8827E-02
            -1.8090E-02
 NUMSIGDIG:         2.9         2.0         2.7         2.3         2.5         2.7         2.5         2.9         1.9         2.8
                    3.7

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1322
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.3991E-04 -2.0635E-02 -2.8098E-02  5.6687E-03 -3.1232E-02
 SE:             2.9885E-02  1.5359E-02  1.3060E-02  2.6302E-02  2.3118E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8559E-01  1.7912E-01  3.1435E-02  8.2936E-01  1.7671E-01

 ETASHRINKSD(%)  1.0000E-10  4.8544E+01  5.6248E+01  1.1884E+01  2.2551E+01
 ETASHRINKVR(%)  1.0000E-10  7.3523E+01  8.0858E+01  2.2356E+01  4.0017E+01
 EBVSHRINKSD(%)  3.2018E-01  4.8169E+01  5.9972E+01  1.2419E+01  1.9844E+01
 EBVSHRINKVR(%)  6.3933E-01  7.3135E+01  8.3977E+01  2.3296E+01  3.5750E+01
 RELATIVEINF(%)  9.8582E+01  9.6488E-01  2.4289E+00  3.1942E+00  1.3608E+01
 EPSSHRINKSD(%)  3.3893E+01
 EPSSHRINKVR(%)  5.6299E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2140.0448423536427     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1221.1063091489700     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.14
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2140.045       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  9.87E-01  1.28E+00  1.08E+00  1.00E+00  1.05E+00  9.60E-01  1.00E+00  1.02E+00  1.04E+00  1.01E+00
 


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
+        9.76E+02
 
 TH 2
+       -9.68E+00  6.70E+07
 
 TH 3
+        1.71E+00  7.26E+01  1.18E+02
 
 TH 4
+       -2.85E+02 -3.54E+07 -7.22E+03  1.83E+07
 
 TH 5
+        2.14E+00 -2.02E+02 -2.09E+02 -3.51E+02  6.60E+02
 
 TH 6
+        8.55E-01 -2.29E+00  4.41E-01 -4.86E+02 -2.49E-01  1.80E+02
 
 TH 7
+        8.06E-01 -1.56E+00  3.08E+00  1.38E+03 -6.78E+00  1.26E-01  2.25E+01
 
 TH 8
+       -3.23E-01 -1.36E+01 -2.80E+01  2.98E+04  5.70E+00  1.90E-02  1.05E+00  1.91E+01
 
 TH 9
+        4.75E+07  5.54E+07 -1.44E+00 -2.93E+07 -5.40E+07 -2.38E-01  2.53E+01  2.24E+00  4.68E+07
 
 TH10
+       -3.74E+07 -4.37E+07 -1.49E+01  2.31E+07  4.26E+07  3.44E-01  4.67E+00  1.07E+01 -7.37E+07  5.81E+07
 
 TH11
+       -7.49E+00 -9.60E+00 -7.81E+00  9.90E+04  1.68E+00  1.18E+00  2.24E+00  6.04E+00  5.93E+00  1.54E+01  3.97E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.625
Stop Time:
Thu Sep 30 02:26:39 CDT 2021
