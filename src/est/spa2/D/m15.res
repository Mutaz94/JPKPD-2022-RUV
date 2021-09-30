Thu Sep 30 08:37:56 CDT 2021
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
$DATA ../../../../data/spa2/D/dat15.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1741.49559488680        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2280E+02 -1.3054E+02 -2.7100E+01 -2.4530E+02  2.3322E+02 -3.7019E+02 -3.0553E+02 -6.6739E+01 -5.3240E+02 -1.8220E+02
            -1.0940E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2193.13265663943        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      170
 NPARAMETR:  9.4950E-01  1.1158E+00  1.2247E+00  1.1503E+00  1.0578E+00  1.6085E+00  3.0741E+00  1.6912E+00  2.1618E+00  1.6402E+00
             1.1218E+00
 PARAMETER:  4.8182E-02  2.0959E-01  3.0272E-01  2.4000E-01  1.5621E-01  5.7531E-01  1.2230E+00  6.2543E-01  8.7092E-01  5.9480E-01
             2.1496E-01
 GRADIENT:  -1.6155E+02 -5.5615E+00 -2.7582E+01  2.5024E+01 -2.2975E+01 -6.8267E+01  1.2310E+01 -3.0965E+00  6.3203E+01  5.4098E+01
             5.6470E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2201.23658737625        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      355
 NPARAMETR:  9.6295E-01  1.0849E+00  1.4104E+00  1.2059E+00  1.1022E+00  1.6081E+00  3.1374E+00  2.0508E+00  2.1236E+00  1.5656E+00
             1.1140E+00
 PARAMETER:  6.2251E-02  1.8144E-01  4.4386E-01  2.8727E-01  1.9727E-01  5.7504E-01  1.2434E+00  8.1825E-01  8.5313E-01  5.4827E-01
             2.0795E-01
 GRADIENT:  -1.4930E+02 -2.7009E+00 -1.9944E+01  3.1122E+01 -1.5490E+01 -6.4018E+01  1.1177E+01  6.7801E+00  5.9351E+01  4.8775E+01
             5.2986E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2241.42259609133        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  1.1482E+00  1.2140E+00  1.5238E+00  1.0148E+00  1.2000E+00  1.7295E+00  3.4207E+00  2.0679E+00  1.3868E+00  9.8515E-01
             1.1091E+00
 PARAMETER:  2.3817E-01  2.9390E-01  5.2120E-01  1.1467E-01  2.8231E-01  6.4781E-01  1.3298E+00  8.2656E-01  4.2698E-01  8.5040E-02
             2.0352E-01
 GRADIENT:   7.3932E+00  2.3968E+00  7.7490E+00 -2.3652E+01  9.3997E+00 -5.9286E+00  2.1209E+01 -9.4590E-01 -1.5996E+00 -1.0003E+01
            -2.7076E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2244.39006017099        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      672
 NPARAMETR:  1.1447E+00  1.1692E+00  1.4581E+00  1.0469E+00  1.1769E+00  1.8160E+00  3.1182E+00  2.0018E+00  1.4033E+00  1.0374E+00
             1.1068E+00
 PARAMETER:  2.3518E-01  2.5634E-01  4.7711E-01  1.4579E-01  2.6286E-01  6.9666E-01  1.2373E+00  7.9406E-01  4.3880E-01  1.3672E-01
             2.0145E-01
 GRADIENT:   8.4452E+02  1.3199E+02  1.6389E+01  1.0151E+02  4.5593E+01  7.6471E+02  7.1373E+02  9.5997E+00  4.6965E+01  1.4924E+00
             2.6723E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2244.49499062193        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  1.1452E+00  1.1944E+00  1.4411E+00  1.0493E+00  1.1738E+00  1.8326E+00  3.0941E+00  2.0108E+00  1.4109E+00  1.0451E+00
             1.1064E+00
 PARAMETER:  2.3554E-01  2.7761E-01  4.6542E-01  1.4809E-01  2.6025E-01  7.0576E-01  1.2295E+00  7.9854E-01  4.4423E-01  1.4415E-01
             2.0115E-01
 GRADIENT:   5.4870E+00 -1.0109E+00  2.2703E-02  2.5344E+00  3.5957E-01  2.2607E+01 -3.8437E+00  8.6595E-01 -1.9871E-02  2.8414E-01
             3.6938E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2244.52414977033        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      997
 NPARAMETR:  1.1481E+00  1.2037E+00  1.4335E+00  1.0413E+00  1.1744E+00  1.8347E+00  3.0864E+00  1.9945E+00  1.4146E+00  1.0420E+00
             1.1062E+00
 PARAMETER:  2.3809E-01  2.8543E-01  4.6010E-01  1.4050E-01  2.6075E-01  7.0687E-01  1.2270E+00  7.9038E-01  4.4684E-01  1.4117E-01
             2.0094E-01
 GRADIENT:   7.3238E+00 -8.2390E-01  9.4296E-01  2.1789E-01 -7.5718E-01  2.2989E+01 -3.3631E+00  1.8328E-01  8.5930E-02 -2.7495E-01
            -4.5383E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2244.53876737120        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1179             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1486E+00  1.2104E+00  1.4193E+00  1.0359E+00  1.1753E+00  1.8371E+00  3.0800E+00  1.9860E+00  1.4172E+00  1.0459E+00
             1.1067E+00
 PARAMETER:  2.3857E-01  2.9098E-01  4.5018E-01  1.3527E-01  2.6154E-01  7.0820E-01  1.2249E+00  7.8611E-01  4.4868E-01  1.4485E-01
             2.0141E-01
 GRADIENT:   8.5908E+02  1.5800E+02  1.3894E+01  9.4168E+01  4.1347E+01  7.8043E+02  6.9226E+02  1.0133E+01  4.8463E+01  2.5220E+00
             2.2691E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2244.53876737120        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  1.1486E+00  1.2104E+00  1.4193E+00  1.0359E+00  1.1753E+00  1.8371E+00  3.0800E+00  1.9860E+00  1.4172E+00  1.0459E+00
             1.1067E+00
 PARAMETER:  2.3857E-01  2.9098E-01  4.5018E-01  1.3527E-01  2.6154E-01  7.0820E-01  1.2249E+00  7.8611E-01  4.4868E-01  1.4485E-01
             2.0141E-01
 GRADIENT:  -1.4100E+04 -5.7816E+03  2.5076E+01 -2.4871E+04  1.2858E+04  2.2812E-04 -2.3432E+00  3.0407E-01  1.5406E-02 -2.8923E-02
             1.6693E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1238
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.7649E-04  6.3004E-03 -4.7495E-02 -2.2703E-03 -3.6425E-02
 SE:             2.9946E-02  2.5892E-02  1.8208E-02  2.3661E-02  2.0204E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8997E-01  8.0775E-01  9.0933E-03  9.2356E-01  7.1405E-02

 ETASHRINKSD(%)  1.0000E-10  1.3259E+01  3.9002E+01  2.0733E+01  3.2315E+01
 ETASHRINKVR(%)  1.0000E-10  2.4760E+01  6.2793E+01  3.7167E+01  5.4187E+01
 EBVSHRINKSD(%)  1.2116E-01  1.0386E+01  4.0878E+01  2.5166E+01  2.9906E+01
 EBVSHRINKVR(%)  2.4217E-01  1.9693E+01  6.5046E+01  4.3999E+01  5.0868E+01
 RELATIVEINF(%)  9.9752E+01  4.9329E+01  2.2455E+01  2.8974E+01  2.7226E+01
 EPSSHRINKSD(%)  3.1795E+01
 EPSSHRINKVR(%)  5.3481E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2244.5387673711998     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1141.8125275255927     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2244.539       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.15E+00  1.21E+00  1.42E+00  1.04E+00  1.18E+00  1.84E+00  3.08E+00  1.99E+00  1.42E+00  1.05E+00  1.11E+00
 


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
+        1.12E+06
 
 TH 2
+       -1.46E+01  1.36E+06
 
 TH 3
+       -4.80E+05  1.26E+03  2.05E+05
 
 TH 4
+        2.19E+06  3.09E+02 -9.38E+05  4.28E+06
 
 TH 5
+       -9.98E+05  4.40E+02  4.27E+05 -1.95E+06  8.90E+05
 
 TH 6
+       -2.97E+01 -2.31E+01  9.27E-03 -5.82E+01  2.61E+01  5.90E+01
 
 TH 7
+       -9.23E+01 -6.50E+01 -3.48E+04 -1.59E+05  7.25E+04 -1.06E-01  5.92E+03
 
 TH 8
+        4.49E+03  3.49E+03 -8.52E+04  8.78E+03 -4.01E+03  4.84E-02  5.57E-01  1.55E+01
 
 TH 9
+       -4.82E+05  3.30E+02 -6.94E+02 -9.44E+05  4.30E+05 -5.97E-02  4.45E+01  5.14E-01  2.08E+05
 
 TH10
+        2.03E+06  5.18E+02 -8.68E+05  3.96E+06 -6.35E+02 -3.61E-02 -4.24E-01  1.70E+00 -8.73E+05  5.66E+01
 
 TH11
+       -1.38E+06  1.18E+03  5.89E+05 -2.69E+06  1.23E+06  3.70E+01  1.00E+05 -5.51E+03  5.93E+05 -2.49E+06  1.69E+06
 
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
 #CPUT: Total CPU Time in Seconds,       40.603
Stop Time:
Thu Sep 30 08:38:39 CDT 2021
