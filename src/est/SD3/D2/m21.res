Sun Oct 24 01:12:07 CDT 2021
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
$DATA ../../../../data/SD3/D2/dat21.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2098.02315135906        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5269E+02  7.2855E+00 -1.3944E+01  1.9562E+01 -2.5406E+01  1.7481E+01 -6.8008E+01  1.6554E+01 -8.5176E+01 -1.3448E+01
             7.4609E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2128.82212476598        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.1326E+00  1.1470E+00  1.3106E+00  1.0280E+00  1.2714E+00  1.0471E+00  1.7754E+00  8.2760E-01  1.4169E+00  1.2157E+00
             1.0452E+00
 PARAMETER:  2.2452E-01  2.3718E-01  3.7051E-01  1.2759E-01  3.4011E-01  1.4603E-01  6.7400E-01 -8.9225E-02  4.4845E-01  2.9531E-01
             1.4421E-01
 GRADIENT:   1.5500E+02  2.6254E+01 -1.2217E+01  4.8225E+01  1.0609E+01 -1.9403E+01  1.4169E+01  5.0381E+00  2.9516E+01 -5.0816E+00
             1.0344E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2135.71580481511        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.1250E+00  8.6579E-01  1.5939E+00  1.2163E+00  1.2427E+00  1.0292E+00  2.1165E+00  5.5279E-01  1.1810E+00  1.3849E+00
             9.8874E-01
 PARAMETER:  2.1776E-01 -4.4115E-02  5.6615E-01  2.9579E-01  3.1730E-01  1.2879E-01  8.4975E-01 -4.9277E-01  2.6639E-01  4.2560E-01
             8.8679E-02
 GRADIENT:   1.5468E+02  3.3486E+01  1.9843E+01  5.0338E+01 -9.7573E+00 -2.5012E+01  1.0557E+01 -1.6479E+00  1.0354E+01  2.9321E+00
             5.8066E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2150.85127094671        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0483E+00  8.9630E-01  1.1008E+00  1.0878E+00  1.0765E+00  1.0703E+00  2.0346E+00  1.9624E-01  1.0342E+00  1.1637E+00
             8.9184E-01
 PARAMETER:  1.4718E-01 -9.4791E-03  1.9608E-01  1.8413E-01  1.7375E-01  1.6791E-01  8.1028E-01 -1.5284E+00  1.3358E-01  2.5158E-01
            -1.4468E-02
 GRADIENT:   4.2564E-01 -1.0500E+00  5.4366E-01 -1.7732E+00 -1.4035E+00 -1.7230E-01  1.7863E-01  2.2246E-01  2.0177E-01 -5.7223E-02
            -3.5614E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2150.94778913714        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0494E+00  9.4963E-01  1.0758E+00  1.0555E+00  1.0899E+00  1.0713E+00  1.9594E+00  8.5700E-02  1.0471E+00  1.1674E+00
             8.9246E-01
 PARAMETER:  1.4825E-01  4.8314E-02  1.7303E-01  1.5405E-01  1.8609E-01  1.6890E-01  7.7265E-01 -2.3569E+00  1.4602E-01  2.5479E-01
            -1.3773E-02
 GRADIENT:   1.4227E+00 -3.4458E-01  2.8480E-01 -1.6039E+00 -7.2447E-01  8.6162E-02  4.5511E-02  4.3543E-02  1.7230E-01  3.7242E-01
             4.2696E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2150.96971612903        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  1.0486E+00  9.3998E-01  1.0752E+00  1.0624E+00  1.0850E+00  1.0708E+00  1.9746E+00  2.0335E-02  1.0422E+00  1.1628E+00
             8.9204E-01
 PARAMETER:  1.4744E-01  3.8108E-02  1.7252E-01  1.6056E-01  1.8161E-01  1.6843E-01  7.8038E-01 -3.7954E+00  1.4133E-01  2.5084E-01
            -1.4241E-02
 GRADIENT:  -3.5512E-02  9.2409E-02  4.5436E-03  1.7227E-01  3.7228E-02 -9.4426E-02  1.3753E-03  2.3965E-03 -4.4053E-02 -3.5581E-02
             3.1272E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2150.97985898305        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1084             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0513E+00  9.3976E-01  1.0744E+00  1.0617E+00  1.0853E+00  1.0720E+00  1.9840E+00  1.0000E-02  1.0402E+00  1.1626E+00
             8.9197E-01
 PARAMETER:  1.5000E-01  3.7865E-02  1.7178E-01  1.5983E-01  1.8181E-01  1.6953E-01  7.8514E-01 -4.8927E+00  1.3943E-01  2.5063E-01
            -1.4322E-02
 GRADIENT:   8.8740E+02  4.5809E+01  3.3190E+00  2.1369E+02  1.4686E+01  9.6679E+01  1.1830E+02  0.0000E+00  1.6942E+01  3.1574E+00
             1.0231E+00

0ITERATION NO.:   32    OBJECTIVE VALUE:  -2150.97985898305        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1149
 NPARAMETR:  1.0513E+00  9.3976E-01  1.0746E+00  1.0617E+00  1.0844E+00  1.0720E+00  1.9886E+00  1.0000E-02  1.0417E+00  1.1632E+00
             8.9192E-01
 PARAMETER:  1.5000E-01  3.7865E-02  1.7178E-01  1.5983E-01  1.8181E-01  1.6953E-01  7.8514E-01 -4.8927E+00  1.3943E-01  2.5063E-01
            -1.4322E-02
 GRADIENT:   5.7125E-04 -1.3879E-03 -5.4393E-02 -2.1909E-02  3.8291E-01  3.8176E-04 -1.8287E-01  0.0000E+00 -1.0742E-01 -4.5127E-02
             2.3653E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1149
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.8102E-04  7.8624E-03 -5.1531E-04 -8.9932E-03 -1.6506E-02
 SE:             2.9919E-02  2.2976E-02  1.5975E-04  2.3151E-02  2.3084E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8184E-01  7.3220E-01  1.2564E-03  6.9767E-01  4.7458E-01

 ETASHRINKSD(%)  1.0000E-10  2.3028E+01  9.9465E+01  2.2442E+01  2.2665E+01
 ETASHRINKVR(%)  1.0000E-10  4.0753E+01  9.9997E+01  3.9847E+01  4.0192E+01
 EBVSHRINKSD(%)  2.4330E-01  2.2140E+01  9.9522E+01  2.3449E+01  1.9180E+01
 EBVSHRINKVR(%)  4.8600E-01  3.9378E+01  9.9998E+01  4.1399E+01  3.4681E+01
 RELATIVEINF(%)  9.9097E+01  9.3703E+00  5.1026E-04  8.5308E+00  1.6187E+01
 EPSSHRINKSD(%)  3.3046E+01
 EPSSHRINKVR(%)  5.5171E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2150.9798589830457     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1232.0413257783730     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2150.980       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  9.40E-01  1.07E+00  1.06E+00  1.09E+00  1.07E+00  1.98E+00  1.00E-02  1.04E+00  1.16E+00  8.92E-01
 


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
 #CPUT: Total CPU Time in Seconds,       39.656
Stop Time:
Sun Oct 24 01:12:15 CDT 2021
