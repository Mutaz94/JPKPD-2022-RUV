Wed Sep 29 14:45:38 CDT 2021
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
$DATA ../../../../data/spa/S1/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1697.88944382532        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1233E+02 -4.0151E+01 -3.2299E+01  2.0728E+01  5.4963E+01  7.4462E+01 -1.3261E+01 -2.4077E+00  1.2534E+01  1.8938E+01
             3.1061E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1704.89705679604        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  1.0052E+00  1.1160E+00  1.0552E+00  9.6582E-01  1.0238E+00  8.6095E-01  1.1013E+00  1.0394E+00  9.6196E-01  8.7184E-01
             8.6619E-01
 PARAMETER:  1.0521E-01  2.0977E-01  1.5377E-01  6.5221E-02  1.2347E-01 -4.9720E-02  1.9647E-01  1.3862E-01  6.1218E-02 -3.7149E-02
            -4.3646E-02
 GRADIENT:   3.6539E+00 -4.7295E+00  1.3398E+01 -4.0073E+00  5.9036E+00 -2.9840E+01 -9.9281E+00 -1.3539E+01 -2.2236E+00  2.0805E+00
            -3.5391E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1707.85167614589        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0068E+00  1.1879E+00  9.9695E-01  9.1535E-01  9.8956E-01  8.7877E-01  1.2335E+00  1.3713E+00  9.1136E-01  6.4087E-01
             9.0886E-01
 PARAMETER:  1.0676E-01  2.7218E-01  9.6948E-02  1.1547E-02  8.9503E-02 -2.9230E-02  3.0989E-01  4.1576E-01  7.1810E-03 -3.4493E-01
             4.4398E-03
 GRADIENT:   5.8331E+00  1.3544E+01  2.1002E+01 -2.1105E+01 -2.3557E+01 -2.0832E+01  7.3814E+00 -5.9187E+00 -5.8541E+00 -8.4396E+00
            -1.5485E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1710.11986561224        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0053E+00  1.2379E+00  1.0947E+00  9.0013E-01  1.0634E+00  9.2293E-01  1.0903E+00  1.6006E+00  9.9242E-01  7.8061E-01
             9.3461E-01
 PARAMETER:  1.0527E-01  3.1345E-01  1.9051E-01 -5.2149E-03  1.6151E-01  1.9793E-02  1.8649E-01  5.7038E-01  9.2388E-02 -1.4769E-01
             3.2370E-02
 GRADIENT:   1.2906E+00  7.0689E+00  3.4252E+00  4.3183E+00 -6.7762E+00 -4.7288E-01  3.7140E-01 -7.2277E-01  1.1150E-01 -4.0197E-01
            -1.1739E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1710.23932042787        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  1.0056E+00  1.3909E+00  9.5460E-01  8.0127E-01  1.0895E+00  9.2480E-01  9.8385E-01  1.6013E+00  1.0738E+00  7.9985E-01
             9.3252E-01
 PARAMETER:  1.0556E-01  4.2998E-01  5.3535E-02 -1.2156E-01  1.8576E-01  2.1823E-02  8.3722E-02  5.7080E-01  1.7124E-01 -1.2334E-01
             3.0132E-02
 GRADIENT:  -2.5791E-01  7.5275E+00  3.1460E-01  9.3721E+00  1.4790E+00 -9.3404E-03 -2.3154E+00 -7.7805E-01 -1.0243E+00 -9.9061E-03
            -9.7476E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1710.34447306148        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  1.0062E+00  1.5356E+00  7.8806E-01  7.0760E-01  1.0943E+00  9.2572E-01  9.2557E-01  1.5592E+00  1.1618E+00  7.8111E-01
             9.3251E-01
 PARAMETER:  1.0614E-01  5.2894E-01 -1.3818E-01 -2.4588E-01  1.9012E-01  2.2817E-02  2.2659E-02  5.4416E-01  2.4994E-01 -1.4703E-01
             3.0128E-02
 GRADIENT:  -7.3861E-01  1.4479E+01 -2.0984E-01  1.4188E+01  2.4981E+00  1.2006E-01 -2.2679E+00 -5.0462E-01 -7.7486E-01 -1.1477E+00
            -1.2186E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1710.62090814577        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1052
 NPARAMETR:  1.0072E+00  1.7462E+00  5.1507E-01  5.6051E-01  1.0782E+00  9.2632E-01  8.6844E-01  1.3625E+00  1.3226E+00  7.3233E-01
             9.3490E-01
 PARAMETER:  1.0717E-01  6.5745E-01 -5.6345E-01 -4.7892E-01  1.7530E-01  2.3469E-02 -4.1059E-02  4.0933E-01  3.7957E-01 -2.1152E-01
             3.2681E-02
 GRADIENT:   1.2335E-02  1.9760E+01  8.6886E-01  1.1885E+01 -4.6118E+00  2.4267E-01  5.7403E-01  3.7574E-01  5.0992E-01 -1.8903E+00
            -3.2518E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1710.75483383344        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  1.0073E+00  1.8557E+00  4.0418E-01  4.8111E-01  1.0858E+00  9.2615E-01  8.3913E-01  1.2855E+00  1.4409E+00  7.2088E-01
             9.3738E-01
 PARAMETER:  1.0731E-01  7.1824E-01 -8.0590E-01 -6.3166E-01  1.8228E-01  2.3286E-02 -7.5386E-02  3.5111E-01  4.6530E-01 -2.2728E-01
             3.5329E-02
 GRADIENT:   5.5947E-01  1.5455E+01  7.0253E-01  8.2750E+00 -6.1164E+00  3.6745E-01  2.5736E+00  8.5253E-01  1.2760E+00 -2.3230E+00
             4.6214E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1711.06594298578        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1415
 NPARAMETR:  1.0068E+00  1.8961E+00  3.6920E-01  4.4472E-01  1.1103E+00  9.2533E-01  8.1598E-01  1.1054E+00  1.4925E+00  7.6927E-01
             9.3532E-01
 PARAMETER:  1.0681E-01  7.3978E-01 -8.9643E-01 -7.1031E-01  2.0464E-01  2.2391E-02 -1.0337E-01  2.0024E-01  5.0044E-01 -1.6231E-01
             3.3135E-02
 GRADIENT:  -9.0057E-01 -5.5174E+00 -8.6014E-01  3.4228E+00  4.6038E+00 -3.7219E-02  4.0096E-01  2.8807E-01 -2.0367E+00  4.4889E-02
             9.9372E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1711.07379986553        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1513
 NPARAMETR:  1.0069E+00  1.8961E+00  3.6880E-01  4.4416E-01  1.1100E+00  9.2534E-01  8.1569E-01  1.0939E+00  1.4925E+00  7.6904E-01
             9.3524E-01
 PARAMETER:  1.0690E-01  7.3978E-01 -8.9749E-01 -7.1157E-01  2.0437E-01  2.2411E-02 -1.0372E-01  1.8975E-01  5.0044E-01 -1.6262E-01
             3.3049E-02
 GRADIENT:   6.4114E+04 -9.2598E+03  7.6658E+03 -9.6437E+03  9.9852E+00 -3.4270E+04 -6.6084E+04 -3.6214E+04 -8.6961E+01 -4.2148E+04
             6.8537E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1513
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3820E-03 -2.6937E-02 -2.1058E-02  2.7423E-02 -4.3294E-02
 SE:             2.9872E-02  2.6314E-02  8.4315E-03  2.1740E-02  1.9555E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6310E-01  3.0599E-01  1.2505E-02  2.0717E-01  2.6832E-02

 ETASHRINKSD(%)  1.0000E-10  1.1845E+01  7.1754E+01  2.7168E+01  3.4488E+01
 ETASHRINKVR(%)  1.0000E-10  2.2287E+01  9.2021E+01  4.6954E+01  5.7082E+01
 EBVSHRINKSD(%)  4.3623E-01  1.2348E+01  7.3334E+01  2.8579E+01  3.2847E+01
 EBVSHRINKVR(%)  8.7055E-01  2.3171E+01  9.2889E+01  4.8990E+01  5.4905E+01
 RELATIVEINF(%)  9.8987E+01  6.1288E+00  4.1788E-01  2.7899E+00  1.0157E+01
 EPSSHRINKSD(%)  4.5415E+01
 EPSSHRINKVR(%)  7.0205E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1711.0737998655281     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -975.92297330178997     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.54
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1711.074       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.90E+00  3.69E-01  4.44E-01  1.11E+00  9.25E-01  8.16E-01  1.09E+00  1.49E+00  7.69E-01  9.35E-01
 


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
+        1.48E+07
 
 TH 2
+       -5.09E+01  8.75E+04
 
 TH 3
+        3.69E+04 -2.68E+03  1.58E+06
 
 TH 4
+       -2.15E+02  1.43E+03 -1.31E+04  1.72E+06
 
 TH 5
+        1.23E+03 -2.75E+02 -5.33E+01  2.77E+01  6.66E+06
 
 TH 6
+        1.14E+03 -8.86E+01  3.71E+02 -3.92E+02 -2.28E+00  2.00E+07
 
 TH 7
+       -1.96E+02  2.41E+01 -8.88E+01  4.49E+01  4.46E+00 -1.45E+03  2.39E+07
 
 TH 8
+       -3.03E+02  2.98E+03 -1.92E+04  1.33E+04 -6.42E+02 -5.93E+02  1.05E+02  4.00E+06
 
 TH 9
+        2.12E+06  9.98E+02 -5.37E+03 -7.22E+05 -1.90E+02  1.71E+00  1.14E+01 -1.10E+06  3.09E+05
 
 TH10
+        2.60E+02 -3.24E+01  6.97E+01 -1.06E+02  6.04E+06 -9.82E+02  1.85E+02 -1.31E+02  7.82E+00  1.10E+07
 
 TH11
+       -1.70E+07  1.37E+02 -5.51E+06  6.60E+02  1.41E+03  1.32E+03 -2.21E+02  1.01E+03  2.45E+06  3.19E+02  1.96E+07
 
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
 #CPUT: Total CPU Time in Seconds,       25.577
Stop Time:
Wed Sep 29 14:46:05 CDT 2021
