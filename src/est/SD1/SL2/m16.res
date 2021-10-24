Sat Oct 23 15:01:48 CDT 2021
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
$DATA ../../../../data/SD1/SL2/dat16.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      999
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

 TOT. NO. OF OBS RECS:      899
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2656.28971872824        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7476E+02  1.4884E+01 -8.0924E-01  1.5843E+02  2.0539E+02  8.0676E+01 -1.0186E+02 -1.1068E+02 -6.6731E+01 -2.7403E+01
            -2.2574E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3239.73901417857        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  1.0002E+00  1.1185E+00  1.0198E+00  9.1847E-01  9.9503E-01  7.6413E-01  1.2895E+00  9.9068E-01  9.6767E-01  1.0753E+00
             1.8249E+00
 PARAMETER:  1.0022E-01  2.1203E-01  1.1964E-01  1.4949E-02  9.5018E-02 -1.6902E-01  3.5424E-01  9.0641E-02  6.7133E-02  1.7260E-01
             7.0154E-01
 GRADIENT:  -1.2551E+01  2.6161E+01 -1.8993E+01 -1.1483E+01  8.5202E+00 -4.5964E+01  2.2450E+01  5.0081E+00 -5.5172E+00  4.5231E+00
             9.8867E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3243.37263592921        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.9633E-01  1.1375E+00  1.0253E+00  9.1240E-01  1.0340E+00  8.2003E-01  1.3012E+00  4.5387E-01  9.7286E-01  1.0785E+00
             1.8767E+00
 PARAMETER:  9.6328E-02  2.2886E-01  1.2499E-01  8.3258E-03  1.3341E-01 -9.8414E-02  3.6329E-01 -6.8995E-01  7.2489E-02  1.7559E-01
             7.2952E-01
 GRADIENT:  -1.4984E+01  2.2034E+01 -2.1977E+01  6.1699E+00  4.7737E+01 -2.2212E+01  2.5565E+01  1.7956E-01 -2.4202E+00  4.0097E-01
             5.3164E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3261.65584292798        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      341
 NPARAMETR:  1.0470E+00  1.5383E+00  1.5894E+00  7.4871E-01  1.4336E+00  8.9747E-01  9.3244E-01  3.2314E+00  1.0799E+00  1.5173E+00
             1.7954E+00
 PARAMETER:  1.4589E-01  5.3069E-01  5.6337E-01 -1.8941E-01  4.6021E-01 -8.1803E-03  3.0046E-02  1.2729E+00  1.7685E-01  5.1694E-01
             6.8523E-01
 GRADIENT:   1.5415E+01  4.8307E+00 -2.9556E+01  3.4218E+01  1.3801E+00  8.5984E+00  5.8535E+00  1.1037E+01  2.1336E+01  6.7369E+00
             8.8291E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3274.43985072303        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  1.0504E+00  1.6086E+00  2.2924E+00  6.9124E-01  1.5945E+00  8.8405E-01  9.1705E-01  3.9681E+00  8.2812E-01  1.5480E+00
             1.7878E+00
 PARAMETER:  1.4915E-01  5.7535E-01  9.2959E-01 -2.6927E-01  5.6655E-01 -2.3246E-02  1.3404E-02  1.4783E+00 -8.8592E-02  5.3694E-01
             6.8099E-01
 GRADIENT:   2.5005E+01 -2.3387E+00 -1.3715E+01 -4.3355E+00  1.4994E+01  2.7641E+00 -5.3731E+00  5.3554E+00  4.4056E+00 -2.5761E+00
             5.2173E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3275.31557712335        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      689
 NPARAMETR:  1.0398E+00  1.5992E+00  2.3289E+00  6.9948E-01  1.5911E+00  8.7704E-01  9.4695E-01  3.9499E+00  7.4125E-01  1.5453E+00
             1.7909E+00
 PARAMETER:  1.3905E-01  5.6952E-01  9.4541E-01 -2.5742E-01  5.6443E-01 -3.1201E-02  4.5488E-02  1.4737E+00 -1.9941E-01  5.3525E-01
             6.8273E-01
 GRADIENT:  -3.6212E+00  4.9278E+00 -1.3306E+01 -2.7673E+00  1.4621E+01 -1.3354E-01 -4.3081E+00  3.9071E+00  1.3432E+00 -2.5905E+00
             9.9516E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3275.36040821070        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      872            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0406E+00  1.5959E+00  2.3217E+00  7.0135E-01  1.5941E+00  8.7756E-01  9.6581E-01  3.9291E+00  7.1790E-01  1.5540E+00
             1.7950E+00
 PARAMETER:  1.3979E-01  5.6745E-01  9.4231E-01 -2.5475E-01  5.6631E-01 -3.0615E-02  6.5208E-02  1.4684E+00 -2.3143E-01  5.4083E-01
             6.8500E-01
 GRADIENT:   1.9802E+02  2.4480E+02 -8.5455E+00  3.8421E+01  8.0912E+01  1.1351E+01  2.3547E+00  1.0914E+01  2.8371E+00  1.2343E+01
             2.7410E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3275.36816457759        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      947
 NPARAMETR:  1.0370E+00  1.5808E+00  2.3205E+00  7.0553E-01  1.5946E+00  8.7560E-01  9.6953E-01  3.9256E+00  7.0742E-01  1.5478E+00
             1.7957E+00
 PARAMETER:  1.3634E-01  5.5795E-01  9.4177E-01 -2.4880E-01  5.6661E-01 -3.2846E-02  6.9052E-02  1.4675E+00 -2.4613E-01  5.3686E-01
             6.8538E-01
 GRADIENT:   1.8316E+02  2.2818E+02 -9.2036E+00  3.2771E+01  8.3821E+01  1.0566E+01  1.6948E+00  1.1053E+01  2.6873E+00  1.1667E+01
             2.9263E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3275.48558180522        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  1.0385E+00  1.5686E+00  2.3202E+00  7.1482E-01  1.5947E+00  8.7703E-01  9.8297E-01  3.9239E+00  6.7299E-01  1.5604E+00
             1.7958E+00
 PARAMETER:  1.3774E-01  5.5017E-01  9.4165E-01 -2.3573E-01  5.6668E-01 -3.1213E-02  8.2820E-02  1.4671E+00 -2.9602E-01  5.4497E-01
             6.8548E-01
 GRADIENT:  -7.3593E+00 -2.5850E+00 -1.5797E+01 -4.0816E+00  2.1454E+01 -1.1868E-01 -2.0001E+00  3.1070E+00 -2.0628E-03  8.2361E-01
             1.8512E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3277.42614030863        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1243            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0414E+00  1.5487E+00  2.5997E+00  7.2974E-01  1.5653E+00  8.7743E-01  1.0007E+00  3.8428E+00  6.6451E-01  1.5482E+00
             1.7873E+00
 PARAMETER:  1.4057E-01  5.3745E-01  1.0554E+00 -2.1506E-01  5.4807E-01 -3.0755E-02  1.0067E-01  1.4462E+00 -3.0870E-01  5.3707E-01
             6.8072E-01
 GRADIENT:   2.0376E+02  2.2068E+02 -1.2163E+00  2.4356E+01  6.9966E+01  1.1446E+01  2.8536E+00  8.0685E+00  1.8215E+00  1.4665E+01
             2.2147E+01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -3277.59042011049        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:     1322
 NPARAMETR:  1.0412E+00  1.5483E+00  2.6557E+00  7.2950E-01  1.5588E+00  8.7730E-01  1.0016E+00  3.8242E+00  6.6592E-01  1.5472E+00
             1.7863E+00
 PARAMETER:  1.4034E-01  5.3691E-01  1.0773E+00 -2.1528E-01  5.4418E-01 -3.0877E-02  1.0089E-01  1.4421E+00 -3.0968E-01  5.3671E-01
             6.7979E-01
 GRADIENT:  -1.3626E+03 -3.5127E+02  1.6971E+02  8.7344E+02  1.7605E+02  8.6451E-02 -8.2506E-01  1.3236E+02 -4.2609E-01  1.7825E+02
            -2.7871E+02
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         2.6         1.2         2.3         1.0         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1322
 NO. OF SIG. DIGITS IN FINAL EST.:  1.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.7661E-04 -1.0444E-02 -2.4655E-02  1.2130E-02 -2.6166E-02
 SE:             2.9676E-02  2.6506E-02  1.8527E-02  1.4935E-02  2.4778E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7375E-01  6.9355E-01  1.8327E-01  4.1670E-01  2.9096E-01

 ETASHRINKSD(%)  5.8031E-01  1.1201E+01  3.7932E+01  4.9965E+01  1.6992E+01
 ETASHRINKVR(%)  1.1573E+00  2.1147E+01  6.1476E+01  7.4965E+01  3.1097E+01
 EBVSHRINKSD(%)  9.2836E-01  1.2079E+01  4.1648E+01  5.3666E+01  1.1850E+01
 EBVSHRINKVR(%)  1.8481E+00  2.2699E+01  6.5950E+01  7.8532E+01  2.2296E+01
 RELATIVEINF(%)  9.8127E+01  1.0909E+01  1.6451E+01  2.6961E+00  5.1929E+01
 EPSSHRINKSD(%)  1.8871E+01
 EPSSHRINKVR(%)  3.4181E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          899
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1652.2514827020016     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3277.5904201104872     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1625.3389374084857     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3277.590       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.55E+00  2.66E+00  7.30E-01  1.56E+00  8.77E-01  1.00E+00  3.83E+00  6.64E-01  1.55E+00  1.79E+00
 


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
 #CPUT: Total CPU Time in Seconds,       93.469
Stop Time:
Sat Oct 23 15:02:04 CDT 2021
