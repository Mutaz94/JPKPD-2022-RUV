Fri Sep 24 19:51:28 CDT 2021
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
$DATA ../../../../data/int/B/dat87.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3405.88997308220        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4294E+01 -4.6992E+01 -3.6521E+00 -5.6721E+01  1.3163E+02 -6.5020E+00 -1.0152E+01 -2.0803E+02 -3.1211E+01 -2.1327E+01
            -7.5723E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3921.72398856306        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.9954E-01  1.0662E+00  1.0880E+00  1.0057E+00  1.0001E+00  1.0194E+00  8.9034E-01  9.8032E-01  1.0402E+00  9.1890E-01
             1.1632E+00
 PARAMETER:  9.9535E-02  1.6409E-01  1.8438E-01  1.0573E-01  1.0007E-01  1.1917E-01 -1.6151E-02  8.0122E-02  1.3938E-01  1.5422E-02
             2.5117E-01
 GRADIENT:  -1.6745E-01  5.2414E+00 -1.1068E+01  7.1556E+00  2.3454E+01 -7.0299E-01  5.1443E+00  6.7823E+00  1.8390E+00 -4.6333E+00
             2.6720E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3922.80807697976        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0078E+00  1.0302E+00  9.8780E-01  1.0087E+00  9.4917E-01  9.6225E-01  8.6124E-01  5.4088E-01  9.9310E-01  1.0589E+00
             1.1314E+00
 PARAMETER:  1.0774E-01  1.2974E-01  8.7721E-02  1.0867E-01  4.7832E-02  6.1517E-02 -4.9381E-02 -5.1457E-01  9.3074E-02  1.5723E-01
             2.2347E-01
 GRADIENT:   2.0425E+01 -2.1770E+00 -2.5963E+01 -1.5811E+01  2.2572E+01 -2.5885E+01  6.6241E+00 -8.3318E+00 -1.7411E+01  2.7810E+01
             1.9614E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3941.63499221634        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9634E-01  1.0214E+00  1.0445E+00  1.0186E+00  9.5048E-01  1.0210E+00  8.7556E-01  8.6773E-01  1.0335E+00  9.2352E-01
             1.0133E+00
 PARAMETER:  9.6329E-02  1.2122E-01  1.4352E-01  1.1842E-01  4.9213E-02  1.2075E-01 -3.2887E-02 -4.1877E-02  1.3298E-01  2.0437E-02
             1.1319E-01
 GRADIENT:   6.5601E+00 -1.4097E+00 -1.5455E+00 -6.3269E-01  2.9641E+00  1.5248E+00  2.4550E-01 -1.2224E+00 -2.8908E+00 -6.2972E-01
             8.1963E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3941.66161080521        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      405             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0338E+00  1.0363E+00  1.0515E+00  1.0184E+00  9.5100E-01  1.0573E+00  8.7926E-01  8.7002E-01  1.0342E+00  9.3346E-01
             1.0130E+00
 PARAMETER:  1.3326E-01  1.3570E-01  1.5022E-01  1.1827E-01  4.9756E-02  1.5576E-01 -2.8673E-02 -3.9239E-02  1.3366E-01  3.1142E-02
             1.1290E-01
 GRADIENT:   9.7683E+01  2.0741E+01  3.1822E+00  6.3812E+00 -1.3114E+01  1.6836E+01  1.5321E+00 -1.7711E+00 -2.7239E+00  7.2370E-01
             6.5750E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3942.11080747961        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      540
 NPARAMETR:  1.0150E+00  1.0297E+00  1.0485E+00  1.0184E+00  9.5100E-01  1.0364E+00  8.7394E-01  8.7002E-01  1.0342E+00  9.3115E-01
             1.0130E+00
 PARAMETER:  1.1490E-01  1.2928E-01  1.4732E-01  1.1827E-01  4.9756E-02  1.3575E-01 -3.4747E-02 -3.9239E-02  1.3366E-01  2.8668E-02
             1.1290E-01
 GRADIENT:   8.1000E-03  2.3146E-02 -2.8369E-02 -8.5929E+00 -1.1018E+01 -3.3751E-02  2.6728E-02 -1.4982E+00 -4.2253E+00  2.4853E-02
             6.8148E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3942.18080834346        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      708             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0150E+00  1.0296E+00  1.0511E+00  1.0201E+00  9.5321E-01  1.0365E+00  8.7372E-01  8.7769E-01  1.0399E+00  9.3102E-01
             1.0118E+00
 PARAMETER:  1.1487E-01  1.2917E-01  1.4981E-01  1.1994E-01  5.2075E-02  1.3584E-01 -3.4995E-02 -3.0467E-02  1.3916E-01  2.8525E-02
             1.1172E-01
 GRADIENT:   5.2930E+01  9.0044E+00  1.0153E+00  7.1887E+00 -2.7847E+00  8.9371E+00  5.3363E-01 -1.1325E+00 -8.4276E-01  3.2707E-01
             4.9212E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3942.18277903134        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      893
 NPARAMETR:  1.0149E+00  1.0296E+00  1.0511E+00  1.0202E+00  9.5329E-01  1.0364E+00  8.7370E-01  8.7806E-01  1.0401E+00  9.3102E-01
             1.0117E+00
 PARAMETER:  1.1482E-01  1.2916E-01  1.4981E-01  1.1998E-01  5.2168E-02  1.3580E-01 -3.5014E-02 -3.0042E-02  1.3932E-01  2.8529E-02
             1.1168E-01
 GRADIENT:   5.2809E+01  8.9345E+00  9.2421E-01  7.3032E+00 -2.6123E+00  8.9181E+00  5.3302E-01 -1.1048E+00 -7.8788E-01  3.2400E-01
             4.8687E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -3942.18277903134        NO. OF FUNC. EVALS.:  68
 CUMULATIVE NO. OF FUNC. EVALS.:      961
 NPARAMETR:  1.0150E+00  1.0296E+00  1.0512E+00  1.0201E+00  9.5332E-01  1.0365E+00  8.7368E-01  8.7803E-01  1.0401E+00  9.3100E-01
             1.0118E+00
 PARAMETER:  1.1482E-01  1.2916E-01  1.4981E-01  1.1998E-01  5.2168E-02  1.3580E-01 -3.5014E-02 -3.0042E-02  1.3932E-01  2.8529E-02
             1.1168E-01
 GRADIENT:  -1.0037E-01  2.1492E+06 -1.0008E-01  2.3136E+06 -1.3879E+06 -1.5579E-02  2.7759E+06  2.7758E+06  1.9924E+06  2.7759E+06
            -2.4859E+06
 NUMSIGDIG:         3.1         3.3         2.8         3.3         3.3         3.2         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      961
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.0590E-04 -2.6733E-02 -1.9555E-02  1.4104E-02 -2.1752E-02
 SE:             2.9903E-02  2.1768E-02  1.6466E-02  2.8181E-02  2.4856E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8917E-01  2.1941E-01  2.3497E-01  6.1673E-01  3.8150E-01

 ETASHRINKSD(%)  1.0000E-10  2.7076E+01  4.4838E+01  5.5900E+00  1.6731E+01
 ETASHRINKVR(%)  1.0000E-10  4.6820E+01  6.9572E+01  1.0867E+01  3.0662E+01
 EBVSHRINKSD(%)  2.5170E-01  2.7079E+01  4.6278E+01  6.7050E+00  1.6402E+01
 EBVSHRINKVR(%)  5.0277E-01  4.6825E+01  7.1139E+01  1.2960E+01  3.0114E+01
 RELATIVEINF(%)  9.9495E+01  2.6190E+01  1.8072E+01  6.3312E+01  2.8007E+01
 EPSSHRINKSD(%)  2.1087E+01
 EPSSHRINKVR(%)  3.7727E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3942.1827790313437     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2288.0934192629329     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3942.183       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.03E+00  1.05E+00  1.02E+00  9.53E-01  1.04E+00  8.74E-01  8.78E-01  1.04E+00  9.31E-01  1.01E+00
 


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
+        9.98E+02
 
 TH 2
+        7.00E+03  3.92E+09
 
 TH 3
+        8.54E-01  1.43E+05  2.80E+09
 
 TH 4
+        7.61E+03  6.21E+04  1.56E+05  4.63E+09
 
 TH 5
+       -9.77E+03 -7.98E+04 -2.00E+05  7.89E+04  7.64E+09
 
 TH 6
+        3.31E+00  1.15E+04  3.13E+09  1.25E+04 -1.61E+04  1.84E+02
 
 TH 7
+        1.07E+04  5.97E+09  2.18E+05 -7.63E+04  9.79E+04  1.76E+04  9.09E+09
 
 TH 8
+        1.06E+04  5.94E+09  2.17E+05 -5.30E+05 -8.29E+09  1.75E+04  9.05E+09  9.00E+09
 
 TH 9
+        6.43E+03  5.22E+04  1.31E+05 -9.41E+04  6.65E+04  1.06E+04 -6.44E+04 -4.47E+05  3.30E+09
 
 TH10
+        1.00E+04  5.61E+09  2.05E+05 -7.53E+04 -7.82E+09  1.65E+04  8.53E+09  8.49E+09 -6.36E+04  8.01E+09
 
 TH11
+       -8.25E+03 -6.70E+04 -1.69E+05  8.32E+04 -8.53E+04 -1.36E+04  8.26E+04  5.74E+05  1.02E+05  8.16E+04  5.44E+09
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       37.049
Stop Time:
Fri Sep 24 19:52:08 CDT 2021
