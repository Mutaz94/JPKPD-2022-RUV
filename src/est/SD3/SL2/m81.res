Sun Oct 24 00:03:31 CDT 2021
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
$DATA ../../../../data/SD3/SL2/dat81.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1599.58946279504        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9172E+02 -5.4814E+01  1.2938E+01  2.9763E+01  1.1203E+02  2.5952E+01 -1.2355E+01 -1.8790E+02  3.3714E+00 -1.5727E+01
            -8.1204E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2063.63501279988        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:       92
 NPARAMETR:  9.3863E-01  1.0537E+00  9.9691E-01  9.9456E-01  9.2952E-01  1.1463E+00  1.0239E+00  1.1291E+00  9.6440E-01  1.0197E+00
             1.4716E+00
 PARAMETER:  3.6671E-02  1.5233E-01  9.6903E-02  9.4550E-02  2.6916E-02  2.3654E-01  1.2360E-01  2.2144E-01  6.3746E-02  1.1955E-01
             4.8632E-01
 GRADIENT:   6.0865E+01  1.5532E+01 -7.6531E+00  3.0016E+01 -1.9175E+01  1.2247E+02  6.8742E-01  9.3230E+00  2.1029E+01  1.1140E+01
             2.0555E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2080.00041588788        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  9.2243E-01  9.6305E-01  8.3524E-01  1.0421E+00  8.4133E-01  1.0827E+00  1.2720E+00  4.4067E-01  6.7906E-01  1.1527E+00
             1.2860E+00
 PARAMETER:  1.9253E-02  6.2355E-02 -8.0037E-02  1.4122E-01 -7.2766E-02  1.7943E-01  3.4062E-01 -7.1946E-01 -2.8704E-01  2.4213E-01
             3.5154E-01
 GRADIENT:   7.2575E+01  1.9979E+01 -1.7974E+01  1.0597E+02 -8.9018E+00  9.1179E+01  6.1788E-01  1.7693E+00 -1.2126E+01  3.0512E+01
             1.3503E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2101.63870106342        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      283
 NPARAMETR:  9.7654E-01  8.5835E-01  1.0692E+00  1.1324E+00  9.0382E-01  1.0619E+00  1.5381E+00  8.4049E-01  7.7662E-01  1.0433E+00
             1.1834E+00
 PARAMETER:  7.6265E-02 -5.2744E-02  1.6694E-01  2.2434E-01 -1.1209E-03  1.6001E-01  5.3055E-01 -7.3773E-02 -1.5280E-01  1.4239E-01
             2.6842E-01
 GRADIENT:  -9.8856E+01  1.2891E+01 -7.2512E+00  2.2337E+01 -1.4145E+00 -1.9726E+01  1.1322E+01  1.0771E+00  5.5132E+00  4.9466E+00
             8.0065E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2109.66998339461        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  1.0245E+00  9.0216E-01  1.0787E+00  1.0924E+00  9.2846E-01  1.0731E+00  1.3451E+00  8.7163E-01  7.8095E-01  1.0503E+00
             1.0681E+00
 PARAMETER:  1.2421E-01 -2.9582E-03  1.7576E-01  1.8835E-01  2.5769E-02  1.7058E-01  3.9646E-01 -3.7386E-02 -1.4725E-01  1.4903E-01
             1.6586E-01
 GRADIENT:  -2.9940E+00 -2.8492E+00 -5.2984E-01  3.1182E-01 -3.5417E-02 -1.0572E+01 -6.5312E-01 -7.6936E-01 -1.2440E+00 -8.1525E-02
             3.9504E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2109.87005217423        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      645
 NPARAMETR:  1.0245E+00  9.0516E-01  1.0867E+00  1.0907E+00  9.3288E-01  1.1031E+00  1.3471E+00  8.9052E-01  7.8770E-01  1.0531E+00
             1.0673E+00
 PARAMETER:  1.2418E-01  3.5197E-04  1.8311E-01  1.8683E-01  3.0524E-02  1.9812E-01  3.9794E-01 -1.5954E-02 -1.3864E-01  1.5175E-01
             1.6511E-01
 GRADIENT:  -2.6946E+00 -3.1447E+00 -6.0845E-01 -8.6309E-01  1.8496E-01  9.4237E-01  3.3256E-01 -5.6564E-01  1.9001E-01  1.2669E-01
             3.7160E+00

0ITERATION NO.:   28    OBJECTIVE VALUE:  -2109.87595299484        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      744
 NPARAMETR:  1.0248E+00  9.0563E-01  1.0868E+00  1.0907E+00  9.3287E-01  1.1032E+00  1.3469E+00  8.9235E-01  7.8763E-01  1.0531E+00
             1.0668E+00
 PARAMETER:  1.2447E-01  8.7991E-04  1.8320E-01  1.8678E-01  3.0516E-02  1.9818E-01  3.9781E-01 -1.3900E-02 -1.3873E-01  1.5171E-01
             1.6465E-01
 GRADIENT:   1.9523E+05  2.4299E+05  1.3347E+05  1.3008E+05  2.4300E+05 -1.2261E+05  3.0542E+04  2.4432E+05  8.7576E+04 -1.6017E+05
            -1.4792E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      744
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8638E-03  3.1295E-03 -2.8887E-02 -6.2053E-03 -2.3669E-02
 SE:             2.9919E-02  2.0807E-02  1.2835E-02  2.2965E-02  2.2637E-02
 N:                     100         100         100         100         100

 P VAL.:         9.2375E-01  8.8044E-01  2.4406E-02  7.8700E-01  2.9574E-01

 ETASHRINKSD(%)  1.0000E-10  3.0295E+01  5.7001E+01  2.3064E+01  2.4164E+01
 ETASHRINKVR(%)  1.0000E-10  5.1411E+01  8.1511E+01  4.0808E+01  4.2489E+01
 EBVSHRINKSD(%)  3.1692E-01  3.0841E+01  6.1177E+01  2.3297E+01  2.1144E+01
 EBVSHRINKVR(%)  6.3283E-01  5.2171E+01  8.4928E+01  4.1167E+01  3.7817E+01
 RELATIVEINF(%)  9.8721E+01  2.6412E+00  2.0446E+00  3.2176E+00  1.3986E+01
 EPSSHRINKSD(%)  3.4073E+01
 EPSSHRINKVR(%)  5.6537E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2109.8759529948375     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1190.9374197901648     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     3.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2109.876       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  9.06E-01  1.09E+00  1.09E+00  9.33E-01  1.10E+00  1.35E+00  8.92E-01  7.88E-01  1.05E+00  1.07E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       24.596
Stop Time:
Sun Oct 24 00:03:38 CDT 2021
