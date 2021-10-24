Sun Oct 24 00:30:03 CDT 2021
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
$DATA ../../../../data/SD3/TD1/dat55.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2100.77216742561        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1454E+02 -6.9130E+01 -3.9839E+01 -1.2595E+01  1.0592E+02  7.8116E+01 -1.1503E+01  3.1612E+00  2.1219E+00 -2.0824E+01
            -2.2440E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2112.85809372251        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  1.0083E+00  1.0681E+00  1.0068E+00  1.0325E+00  9.6955E-01  8.3176E-01  1.0379E+00  9.9516E-01  1.0072E+00  1.0502E+00
             1.0413E+00
 PARAMETER:  1.0826E-01  1.6585E-01  1.0677E-01  1.3201E-01  6.9078E-02 -8.4212E-02  1.3721E-01  9.5149E-02  1.0721E-01  1.4895E-01
             1.4043E-01
 GRADIENT:   4.1136E+01 -7.9646E+00 -1.2502E+01  1.4624E+01  1.9865E+01 -1.0881E+01 -4.1335E+00  2.5957E+00  1.3540E+00 -4.1612E+00
             9.1202E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2113.37934105814        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0020E+00  1.2079E+00  9.6954E-01  9.5234E-01  1.0176E+00  8.3044E-01  1.0196E+00  8.7912E-01  1.0701E+00  1.1531E+00
             1.0121E+00
 PARAMETER:  1.0203E-01  2.8888E-01  6.9068E-02  5.1162E-02  1.1740E-01 -8.5799E-02  1.1939E-01 -2.8832E-02  1.6780E-01  2.4243E-01
             1.1207E-01
 GRADIENT:   1.9799E+01  8.9768E+00 -2.9063E+00  2.1089E+01  9.5414E+00 -1.1716E+01  1.5257E+00 -8.0908E-01  2.4868E+00  4.5131E+00
            -1.2578E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2114.98645714358        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.9800E-01  1.3527E+00  6.3015E-01  8.3891E-01  8.8693E-01  8.5905E-01  9.7058E-01  4.3433E-01  1.1261E+00  9.6896E-01
             1.0256E+00
 PARAMETER:  9.7997E-02  4.0210E-01 -3.6180E-01 -7.5657E-02 -1.9995E-02 -5.1923E-02  7.0134E-02 -7.3395E-01  2.1874E-01  6.8472E-02
             1.2530E-01
 GRADIENT:  -1.6706E+00  1.9450E+01  1.8122E+00  1.2494E+01 -2.0103E+01  9.2842E-01  2.1628E+00  9.1065E-01  3.5815E-01  3.3351E+00
             2.1931E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2115.89080400775        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.9883E-01  1.5828E+00  5.6545E-01  6.8551E-01  1.0006E+00  8.5948E-01  8.4884E-01  2.6700E-01  1.3046E+00  1.0371E+00
             1.0236E+00
 PARAMETER:  9.8834E-02  5.5920E-01 -4.7014E-01 -2.7759E-01  1.0059E-01 -5.1431E-02 -6.3889E-02 -1.2205E+00  3.6588E-01  1.3645E-01
             1.2328E-01
 GRADIENT:  -5.5115E-01  3.8314E+00 -6.2562E-01  3.9830E+00 -2.6882E-01  7.1441E-01 -7.9777E-01  2.0260E-01 -3.3695E-01 -9.0105E-01
             1.4726E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2115.95771374605        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.9898E-01  1.6185E+00  5.5512E-01  6.6205E-01  1.0191E+00  8.5792E-01  8.4108E-01  1.0569E-01  1.3405E+00  1.0576E+00
             1.0235E+00
 PARAMETER:  9.8975E-02  5.8149E-01 -4.8857E-01 -3.1241E-01  1.1891E-01 -5.3248E-02 -7.3072E-02 -2.1472E+00  3.9304E-01  1.5596E-01
             1.2318E-01
 GRADIENT:  -1.9691E-01  3.1386E+00  1.6439E-01  2.9954E+00  6.8264E-02 -3.3324E-02  8.5677E-03  2.4446E-02  1.0751E-01 -2.2414E-01
            -1.7963E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2115.97790477151        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1044
 NPARAMETR:  9.9956E-01  1.6189E+00  5.5196E-01  6.5789E-01  1.0201E+00  8.5810E-01  8.3974E-01  2.1109E-02  1.3445E+00  1.0586E+00
             1.0235E+00
 PARAMETER:  9.9557E-02  5.8174E-01 -4.9428E-01 -3.1871E-01  1.1990E-01 -5.3032E-02 -7.4664E-02 -3.7580E+00  3.9606E-01  1.5696E-01
             1.2323E-01
 GRADIENT:   3.8650E+02  4.8868E+02  4.6899E+00  8.7411E+01  1.0309E+01  2.4374E+01  5.4981E+00  2.5669E-03  1.8663E+01  1.4891E+00
             1.2054E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2115.98031711649        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1118
 NPARAMETR:  9.9952E-01  1.6129E+00  5.5372E-01  6.6382E-01  1.0159E+00  8.5841E-01  8.4139E-01  1.0000E-02  1.3396E+00  1.0549E+00
             1.0232E+00
 PARAMETER:  9.9522E-02  5.7804E-01 -4.9110E-01 -3.0974E-01  1.1579E-01 -5.2671E-02 -7.2694E-02 -1.5331E+01  3.9240E-01  1.5346E-01
             1.2296E-01
 GRADIENT:   3.8645E+02  4.8461E+02  4.7482E+00  8.8313E+01  9.7025E+00  2.4518E+01  5.3527E+00  0.0000E+00  1.8911E+01  1.3262E+00
             9.4023E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2115.98053387377        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:     1197
 NPARAMETR:  9.9888E-01  1.6094E+00  5.5452E-01  6.6663E-01  1.0136E+00  8.5807E-01  8.4320E-01  1.0000E-02  1.3357E+00  1.0531E+00
             1.0232E+00
 PARAMETER:  9.8875E-02  5.7587E-01 -4.8965E-01 -3.0552E-01  1.1353E-01 -5.3070E-02 -7.0549E-02 -2.5722E+01  3.8947E-01  1.5174E-01
             1.2293E-01
 GRADIENT:   3.8446E+02  4.8143E+02  4.7784E+00  8.8336E+01  9.3067E+00  2.4410E+01  5.4056E+00  0.0000E+00  1.8816E+01  1.3138E+00
             9.3854E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2115.98058281506        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:     1301
 NPARAMETR:  9.9850E-01  1.6079E+00  5.5478E-01  6.6759E-01  1.0127E+00  8.5786E-01  8.4377E-01  1.0000E-02  1.3342E+00  1.0524E+00
             1.0232E+00
 PARAMETER:  9.8501E-02  5.7494E-01 -4.8918E-01 -3.0409E-01  1.1266E-01 -5.3315E-02 -6.9870E-02 -3.0893E+01  3.8830E-01  1.5104E-01
             1.2291E-01
 GRADIENT:  -1.6072E+00  1.1076E+00 -3.8226E-01  2.6500E+00  4.3619E-01 -5.9426E-02 -8.1072E-02  0.0000E+00  5.1704E-01 -1.3734E-01
            -2.9487E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -2115.98586748065        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1367
 NPARAMETR:  9.9931E-01  1.6020E+00  5.5399E-01  6.6369E-01  1.0143E+00  8.5755E-01  8.4375E-01  1.0000E-02  1.3320E+00  1.0540E+00
             1.0237E+00
 PARAMETER:  9.9624E-02  5.7541E-01 -4.8889E-01 -3.0687E-01  1.1317E-01 -5.2726E-02 -6.9126E-02 -2.6463E+01  3.8756E-01  1.5205E-01
             1.2321E-01
 GRADIENT:   3.5434E-01  3.0297E+00  1.9754E-01  9.1310E-01 -5.0363E-01  1.3877E-01  4.3082E-02  0.0000E+00  6.3049E-02 -3.2200E-02
            -6.0652E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1367
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9020E-04 -2.9654E-02 -2.8980E-04  2.2579E-02 -3.0389E-02
 SE:             2.9814E-02  2.2990E-02  1.2493E-04  2.4053E-02  2.3223E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9223E-01  1.9711E-01  2.0362E-02  3.4788E-01  1.9069E-01

 ETASHRINKSD(%)  1.2019E-01  2.2980E+01  9.9581E+01  1.9418E+01  2.2198E+01
 ETASHRINKVR(%)  2.4023E-01  4.0679E+01  9.9998E+01  3.5066E+01  3.9469E+01
 EBVSHRINKSD(%)  4.7274E-01  2.2872E+01  9.9639E+01  2.0347E+01  2.0611E+01
 EBVSHRINKVR(%)  9.4324E-01  4.0512E+01  9.9999E+01  3.6554E+01  3.6974E+01
 RELATIVEINF(%)  9.8989E+01  5.2737E+00  2.9782E-04  6.1238E+00  1.2571E+01
 EPSSHRINKSD(%)  3.4062E+01
 EPSSHRINKVR(%)  5.6522E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2115.9858674806496     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1197.0473342759769     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2115.986       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.61E+00  5.55E-01  6.66E-01  1.01E+00  8.58E-01  8.44E-01  1.00E-02  1.33E+00  1.05E+00  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,       41.629
Stop Time:
Sun Oct 24 00:30:12 CDT 2021
