Wed Sep 29 12:27:44 CDT 2021
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
$DATA ../../../../data/spa/A1/dat93.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1280.32702784695        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5765E+02  4.6217E+01  2.3979E+01  6.6631E+01  1.3841E+02  5.9422E+01 -1.0311E+01 -7.1503E+00  2.3039E+00 -7.6354E+01
            -6.4588E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1446.91680159981        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1191E+00  9.1991E-01  8.9093E-01  1.0740E+00  8.1002E-01  1.1393E+00  9.3644E-01  9.3402E-01  9.4621E-01  1.0727E+00
             2.3951E+00
 PARAMETER:  2.1251E-01  1.6521E-02 -1.5495E-02  1.7135E-01 -1.1070E-01  2.3044E-01  3.4334E-02  3.1748E-02  4.4709E-02  1.7019E-01
             9.7344E-01
 GRADIENT:   3.4745E+02  2.4745E+01  1.4122E+01  2.3301E+01 -3.4370E+01  5.5018E+01  6.0678E+00  8.7088E+00  7.9192E+00  2.0261E+01
             8.2726E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1461.24096596096        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  1.1037E+00  6.1209E-01  4.1554E-01  1.2619E+00  4.4751E-01  1.0229E+00  3.6773E-01  3.8024E-01  9.9003E-01  5.9448E-01
             2.2368E+00
 PARAMETER:  1.9869E-01 -3.9088E-01 -7.7817E-01  3.3260E-01 -7.0406E-01  1.2268E-01 -9.0039E-01 -8.6694E-01  8.9980E-02 -4.2007E-01
             9.0504E-01
 GRADIENT:   1.9170E+02  5.7821E+01 -2.7440E+01  1.8462E+02  3.8308E+00  8.0315E+00 -2.5644E+00  1.1829E+00  3.2678E+00 -9.0605E+00
             4.4300E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1490.46027128143        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      359
 NPARAMETR:  9.8937E-01  5.1593E-01  5.0575E-01  1.2669E+00  4.8082E-01  9.4555E-01  1.0765E+00  9.2886E-02  8.4651E-01  7.5470E-01
             1.8517E+00
 PARAMETER:  8.9308E-02 -5.6178E-01 -5.8171E-01  3.3661E-01 -6.3226E-01  4.4013E-02  1.7374E-01 -2.2764E+00 -6.6633E-02 -1.8144E-01
             7.1610E-01
 GRADIENT:  -1.9402E+01  2.7290E+01 -1.5107E+01  9.9212E+01  1.5633E+01 -4.6179E-01 -2.8687E+00  1.3922E-01 -1.0148E+01  2.2310E+00
            -1.0163E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1496.73844353212        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  9.9723E-01  3.4758E-01  4.1179E-01  1.2511E+00  3.7882E-01  9.4476E-01  1.3424E+00  1.0000E-02  8.5181E-01  6.8764E-01
             1.8442E+00
 PARAMETER:  9.7225E-02 -9.5676E-01 -7.8724E-01  3.2405E-01 -8.7071E-01  4.3171E-02  3.9448E-01 -5.0898E+00 -6.0390E-02 -2.7449E-01
             7.1204E-01
 GRADIENT:   1.9434E+00  5.0863E+00  9.7014E+00 -6.0076E+00 -1.4903E+01 -3.2987E-01 -2.1824E-01  0.0000E+00  8.7609E-01 -5.5462E-01
            -1.2755E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1498.55140003818        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  9.8524E-01  1.4566E-01  4.8579E-01  1.3869E+00  3.9482E-01  9.2878E-01  2.6149E+00  1.0000E-02  7.9411E-01  7.4383E-01
             1.8586E+00
 PARAMETER:  8.5125E-02 -1.8265E+00 -6.2197E-01  4.2705E-01 -8.2934E-01  2.6122E-02  1.0612E+00 -7.4087E+00 -1.3053E-01 -1.9594E-01
             7.1983E-01
 GRADIENT:  -4.9080E+00  2.5702E+00  5.0313E+00  1.9372E+01 -1.1086E+01 -3.7024E+00  7.5673E-02  0.0000E+00 -3.6271E+00 -1.1787E+00
             8.0994E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1499.06832258885        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.8096E-01  3.8550E-02  5.0848E-01  1.4345E+00  3.9674E-01  9.3532E-01  5.2118E+00  1.0000E-02  7.8439E-01  7.7268E-01
             1.8502E+00
 PARAMETER:  8.0781E-02 -3.1558E+00 -5.7633E-01  4.6079E-01 -8.2446E-01  3.3138E-02  1.7509E+00 -1.2674E+01 -1.4285E-01 -1.5789E-01
             7.1532E-01
 GRADIENT:  -1.2538E+00  1.3850E-01 -2.3512E+00 -1.1972E+00  3.4634E+00  2.3784E-01  7.0433E-02  0.0000E+00 -4.6825E-01 -1.1038E-01
             1.5000E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1499.10732146473        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  9.8021E-01  1.0559E-02  5.1688E-01  1.4520E+00  3.9703E-01  9.3369E-01  9.3582E+00  1.0000E-02  7.8120E-01  7.8046E-01
             1.8498E+00
 PARAMETER:  8.0013E-02 -4.4508E+00 -5.5994E-01  4.7292E-01 -8.2375E-01  3.1389E-02  2.3363E+00 -1.8037E+01 -1.4693E-01 -1.4787E-01
             7.1509E-01
 GRADIENT:   2.1775E-01  5.6925E-02  2.0114E+00  2.3096E+00 -3.2618E+00 -1.2971E-01 -3.8681E-03  0.0000E+00  1.3753E-01  9.4161E-02
            -7.5509E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1499.11491435243        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1254             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8017E-01  1.0000E-02  5.1652E-01  1.4497E+00  3.9732E-01  9.3401E-01  1.0044E+01  1.0000E-02  7.8082E-01  7.8023E-01
             1.8498E+00
 PARAMETER:  7.9972E-02 -4.5781E+00 -5.6063E-01  4.7136E-01 -8.2302E-01  3.1732E-02  2.4070E+00 -1.8299E+01 -1.4741E-01 -1.4817E-01
             7.1508E-01
 GRADIENT:   1.0404E+02  0.0000E+00  1.4838E+01  2.0469E+02  7.6087E+01  9.3646E+00  1.8052E-01  0.0000E+00  7.0573E+00  1.0139E+00
             5.7568E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1499.11515904417        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1390
 NPARAMETR:  9.8015E-01  1.0000E-02  5.1647E-01  1.4497E+00  3.9721E-01  9.3398E-01  9.7718E+00  1.0000E-02  7.8077E-01  7.8016E-01
             1.8496E+00
 PARAMETER:  7.9948E-02 -4.5781E+00 -5.6075E-01  4.7133E-01 -8.2328E-01  3.1694E-02  2.3795E+00 -1.8299E+01 -1.4748E-01 -1.4825E-01
             7.1497E-01
 GRADIENT:   3.0837E-01  0.0000E+00  4.7363E-01 -3.0421E+00  9.8950E-02  1.5967E-02  1.2373E-03  0.0000E+00  1.3276E-02  3.3592E-03
            -3.3438E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1390
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.6485E-05  4.7847E-04 -4.5193E-05 -7.2468E-03 -8.6344E-03
 SE:             2.9460E-02  1.5964E-03  2.3996E-04  2.7829E-02  2.4377E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9928E-01  7.6439E-01  8.5061E-01  7.9455E-01  7.2319E-01

 ETASHRINKSD(%)  1.3043E+00  9.4652E+01  9.9196E+01  6.7687E+00  1.8332E+01
 ETASHRINKVR(%)  2.5917E+00  9.9714E+01  9.9994E+01  1.3079E+01  3.3304E+01
 EBVSHRINKSD(%)  1.4535E+00  9.4894E+01  9.9195E+01  6.2810E+00  1.7333E+01
 EBVSHRINKVR(%)  2.8860E+00  9.9739E+01  9.9994E+01  1.2168E+01  3.1662E+01
 RELATIVEINF(%)  8.3452E+01  1.4662E-02  2.9509E-04  7.7814E+00  2.6031E+00
 EPSSHRINKSD(%)  3.7657E+01
 EPSSHRINKVR(%)  6.1133E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1499.1151590441732     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -763.96433248043502     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.19
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1499.115       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.00E-02  5.16E-01  1.45E+00  3.97E-01  9.34E-01  9.77E+00  1.00E-02  7.81E-01  7.80E-01  1.85E+00
 


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
+        1.29E+03
 
 TH 2
+        0.00E+00  1.92E+02
 
 TH 3
+       -5.91E-02  0.00E+00  2.73E+03
 
 TH 4
+       -3.05E+01  0.00E+00 -2.98E+02  7.63E+02
 
 TH 5
+        5.63E+01  0.00E+00 -4.33E+03 -1.56E+02  7.94E+03
 
 TH 6
+       -7.15E-01  0.00E+00  8.22E+00 -7.33E+00 -7.24E-01  2.15E+02
 
 TH 7
+        2.77E-03  0.00E+00  1.30E-02 -3.82E-03 -1.56E-02  4.05E-03  1.79E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.80E+00  0.00E+00  4.91E+01 -9.42E+00  2.40E+00  3.72E+00  2.55E-02  0.00E+00  2.55E+02
 
 TH10
+       -3.40E+00  0.00E+00 -4.86E+01  3.00E+00 -6.88E+01 -3.51E-01  1.14E-02  0.00E+00 -2.30E+00  1.51E+02
 
 TH11
+       -1.21E+01  0.00E+00 -2.33E+01 -9.61E+00  1.70E+01  1.82E+00  4.41E-03  0.00E+00  1.01E+01  2.82E+01  7.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       23.186
Stop Time:
Wed Sep 29 12:28:11 CDT 2021
