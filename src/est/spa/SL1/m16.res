Wed Sep 29 14:52:56 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat16.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1717.92477114617        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5474E+02 -6.5750E+01 -5.8045E+01  2.1741E+01  1.1235E+02  7.6510E+01 -5.1790E+00  5.1504E+00  3.3242E+01 -7.8659E+00
             2.2565E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1730.51187925367        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0170E+00  1.1181E+00  1.0613E+00  9.7587E-01  9.7368E-01  8.5204E-01  1.0268E+00  9.8523E-01  8.6308E-01  1.0026E+00
             9.2998E-01
 PARAMETER:  1.1688E-01  2.1163E-01  1.5947E-01  7.5576E-02  7.3326E-02 -6.0118E-02  1.2642E-01  8.5122E-02 -4.7242E-02  1.0256E-01
             2.7413E-02
 GRADIENT:  -5.0677E+01  4.8389E+00  1.7151E+01 -8.9448E+00 -2.6239E+01 -2.0853E+01 -7.2089E+00 -2.2467E+00  9.0235E-01 -3.2294E+00
            -1.3305E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1731.04117367099        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0308E+00  1.1117E+00  9.4361E-01  9.6618E-01  9.2540E-01  8.8328E-01  1.2498E+00  9.6454E-01  7.1002E-01  9.1250E-01
             9.6254E-01
 PARAMETER:  1.3038E-01  2.0587E-01  4.1962E-02  6.5599E-02  2.2470E-02 -2.4108E-02  3.2301E-01  6.3901E-02 -2.4246E-01  8.4280E-03
             6.1820E-02
 GRADIENT:  -1.0907E+01  5.1666E+00  1.5856E+01 -3.0364E+01 -2.5271E+01 -5.0590E+00  8.2505E+00 -1.9514E-01 -8.9181E+00 -4.8076E+00
             7.6476E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1732.76760666602        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0365E+00  1.1623E+00  8.3846E-01  9.4127E-01  9.1656E-01  8.9617E-01  1.0995E+00  6.6293E-01  8.2291E-01  9.6343E-01
             9.5613E-01
 PARAMETER:  1.3583E-01  2.5038E-01 -7.6183E-02  3.9478E-02  1.2872E-02 -9.6298E-03  1.9489E-01 -3.1108E-01 -9.4909E-02  6.2742E-02
             5.5141E-02
 GRADIENT:   1.7445E+00  3.2873E+00 -4.7239E-01  5.6820E+00 -1.6613E+00  1.9473E-01 -5.6818E-01  4.1214E-01  2.6764E-01  1.2274E+00
             1.4784E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1732.93366680493        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0370E+00  1.3619E+00  6.4180E-01  8.0308E-01  9.1264E-01  8.9654E-01  9.8230E-01  2.8871E-01  8.9767E-01  9.2617E-01
             9.5671E-01
 PARAMETER:  1.3637E-01  4.0891E-01 -3.4348E-01 -1.1930E-01  8.5839E-03 -9.2137E-03  8.2146E-02 -1.1423E+00 -7.9486E-03  2.3303E-02
             5.5748E-02
 GRADIENT:  -9.8207E-01  1.0747E+00  1.2690E-01  1.9158E-01 -1.5696E+00 -4.2804E-01  8.0978E-01  9.3802E-02 -5.2097E-02 -3.1395E-02
             2.3584E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1732.95599102359        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      905             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0395E+00  1.3968E+00  6.2331E-01  7.7876E-01  9.2344E-01  8.9837E-01  9.5478E-01  1.4982E-01  9.1764E-01  9.3125E-01
             9.5637E-01
 PARAMETER:  1.3877E-01  4.3422E-01 -3.7271E-01 -1.5005E-01  2.0349E-02 -7.1732E-03  5.3730E-02 -1.7983E+00  1.4054E-02  2.8767E-02
             5.5391E-02
 GRADIENT:   6.5440E+02  2.9493E+02  6.4898E+00  6.8819E+01  1.0126E+01  3.8998E+01  4.6983E+00  6.4189E-02  2.8829E+00 -5.4343E-02
             4.5967E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1732.96302223329        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  1.0379E+00  1.3982E+00  6.2238E-01  7.7968E-01  9.2279E-01  8.9737E-01  9.5745E-01  1.2284E-01  9.1932E-01  9.3606E-01
             9.5680E-01
 PARAMETER:  1.3718E-01  4.3518E-01 -3.7420E-01 -1.4887E-01  1.9651E-02 -8.2871E-03  5.6521E-02 -1.9969E+00  1.5884E-02  3.3920E-02
             5.5835E-02
 GRADIENT:   1.0576E+00  2.9188E-01  1.9003E-01  5.2939E-01 -2.6645E-01 -1.1263E-01  7.6900E-02  1.1287E-02  1.2153E-02  4.2502E-02
             5.3170E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1732.96440535560        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  1.0371E+00  1.4003E+00  6.1972E-01  7.7781E-01  9.2289E-01  8.9751E-01  9.5612E-01  6.3623E-02  9.2086E-01  9.3673E-01
             9.5669E-01
 PARAMETER:  1.3640E-01  4.3671E-01 -3.7848E-01 -1.5127E-01  1.9756E-02 -8.1300E-03  5.5127E-02 -2.6548E+00  1.7556E-02  3.4640E-02
             5.5729E-02
 GRADIENT:  -1.1231E+00 -4.0282E-01  4.9407E-02  2.2311E-01  6.4409E-02 -6.1521E-02  7.7061E-02  3.1753E-03  5.0650E-02  1.3381E-01
             3.3883E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1732.96515313028        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  1.0371E+00  1.4056E+00  6.1512E-01  7.7420E-01  9.2300E-01  8.9768E-01  9.5306E-01  2.4683E-02  9.2301E-01  9.3509E-01
             9.5661E-01
 PARAMETER:  1.3646E-01  4.4047E-01 -3.8594E-01 -1.5592E-01  1.9870E-02 -7.9396E-03  5.1922E-02 -3.6016E+00  1.9885E-02  3.2889E-02
             5.5641E-02
 GRADIENT:  -1.0345E+00 -4.8270E-01  1.3753E-02  1.6183E-01  5.1798E-02 -4.6026E-04  2.2920E-02  5.6469E-04  9.3236E-03  5.9229E-02
             8.5027E-03

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1732.96799681641        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:     1524
 NPARAMETR:  1.0382E+00  1.4042E+00  6.1505E-01  7.7391E-01  9.2294E-01  8.9763E-01  9.5289E-01  1.0000E-02  9.2290E-01  9.3474E-01
             9.5658E-01
 PARAMETER:  1.3752E-01  4.3947E-01 -3.8605E-01 -1.5629E-01  1.9812E-02 -7.9967E-03  5.1740E-02 -4.6529E+00  1.9768E-02  3.2513E-02
             5.5614E-02
 GRADIENT:   1.9536E+00 -2.4035E+00 -3.6250E-02 -1.2391E+00  4.6800E-01 -1.3708E-02 -7.0978E-02  0.0000E+00 -1.1280E-03  3.5194E-02
             1.0118E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1524
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3991E-04 -1.2871E-02 -3.5600E-04  9.9006E-03 -2.1000E-02
 SE:             2.9844E-02  2.4008E-02  1.4066E-04  2.2199E-02  2.3221E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9091E-01  5.9189E-01  1.1374E-02  6.5561E-01  3.6582E-01

 ETASHRINKSD(%)  1.8566E-02  1.9569E+01  9.9529E+01  2.5629E+01  2.2205E+01
 ETASHRINKVR(%)  3.7128E-02  3.5308E+01  9.9998E+01  4.4690E+01  3.9480E+01
 EBVSHRINKSD(%)  4.6961E-01  1.9501E+01  9.9593E+01  2.7226E+01  2.0559E+01
 EBVSHRINKVR(%)  9.3702E-01  3.5199E+01  9.9998E+01  4.7040E+01  3.6891E+01
 RELATIVEINF(%)  9.8895E+01  3.2084E+00  1.3201E-04  2.3234E+00  7.5436E+00
 EPSSHRINKSD(%)  4.4531E+01
 EPSSHRINKVR(%)  6.9232E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1732.9679968164075     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -997.81717025266937     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1732.968       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.40E+00  6.15E-01  7.74E-01  9.23E-01  8.98E-01  9.53E-01  1.00E-02  9.23E-01  9.35E-01  9.57E-01
 


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
+        1.27E+03
 
 TH 2
+       -6.81E+00  4.15E+02
 
 TH 3
+        1.58E+01  1.56E+02  5.34E+02
 
 TH 4
+       -2.05E+01  3.91E+02 -4.23E+02  1.17E+03
 
 TH 5
+       -5.20E+00 -2.53E+02 -5.92E+02  4.40E+02  9.42E+02
 
 TH 6
+        9.74E-01 -1.19E+00  3.08E+00 -3.79E+00 -1.13E+00  2.43E+02
 
 TH 7
+        7.39E-01  2.18E+01 -1.73E+01 -1.55E+01  1.40E+00 -2.85E-01  9.76E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.45E+00 -1.78E+01 -4.36E+01  4.96E+01  7.57E+00 -6.32E-01  2.46E+01  0.00E+00  7.74E+01
 
 TH10
+       -8.45E-01 -1.56E+01 -4.89E+01 -6.36E+00 -6.51E+01  1.65E-01  1.30E+01  0.00E+00  1.44E+01  9.25E+01
 
 TH11
+       -7.54E+00 -1.50E+01 -3.11E+01  3.88E+00  1.20E+00  2.71E+00  7.43E+00  0.00E+00  1.16E+01  1.79E+01  2.30E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.174
Stop Time:
Wed Sep 29 14:53:23 CDT 2021
