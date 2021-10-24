Sun Oct 24 03:28:37 CDT 2021
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
$DATA ../../../../data/SD4/SL3/dat24.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      499
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

 TOT. NO. OF OBS RECS:      399
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1615.10515187447        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3572E+02 -5.2248E+01  3.0896E+01 -1.0759E+02 -3.1678E+00  3.1094E+01 -1.9741E+00 -4.2869E+00 -1.3081E+01 -7.6334E+00
            -2.4099E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1625.74044524247        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7970E-01  1.0298E+00  9.6754E-01  1.0918E+00  9.8816E-01  1.0190E+00  1.0036E+00  1.0080E+00  1.0098E+00  1.0258E+00
             1.0516E+00
 PARAMETER:  7.9495E-02  1.2939E-01  6.6999E-02  1.8780E-01  8.8092E-02  1.1884E-01  1.0358E-01  1.0794E-01  1.0975E-01  1.2545E-01
             1.5032E-01
 GRADIENT:   2.5242E+00  5.7163E+00  6.6590E+00 -2.2069E+00 -5.0445E+00 -1.0796E+00  2.4053E-01 -9.8368E-02  1.0004E+00 -1.9298E+00
            -9.9902E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1625.90002496498        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  9.7967E-01  9.1847E-01  9.6640E-01  1.1714E+00  9.3278E-01  1.0123E+00  1.0114E+00  9.2301E-01  9.7237E-01  1.0208E+00
             1.0420E+00
 PARAMETER:  7.9463E-02  1.4951E-02  6.5821E-02  2.5821E-01  3.0417E-02  1.1224E-01  1.1133E-01  1.9888E-02  7.1982E-02  1.2062E-01
             1.4112E-01
 GRADIENT:   3.5646E+00  1.3945E+01  6.1371E+00  1.5721E+01 -9.3482E+00 -3.5727E+00 -1.3743E+00 -1.3575E+00 -2.7220E-01  2.2718E-01
            -4.1554E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1626.90210923540        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  9.7857E-01  7.8430E-01  8.5440E-01  1.2477E+00  8.1055E-01  1.0204E+00  1.3545E+00  7.4444E-01  8.6550E-01  8.8331E-01
             1.0513E+00
 PARAMETER:  7.8340E-02 -1.4296E-01 -5.7353E-02  3.2126E-01 -1.1005E-01  1.2015E-01  4.0340E-01 -1.9513E-01 -4.4443E-02 -2.4079E-02
             1.5000E-01
 GRADIENT:   5.8319E-01  1.4898E+01  3.9975E+00  2.8553E+01 -4.8060E+00 -4.3198E-01 -2.5936E+00 -7.9691E-01 -1.0410E+00 -1.4430E+00
             2.2042E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1627.91706623552        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  9.7685E-01  5.3809E-01  7.7246E-01  1.3651E+00  6.7788E-01  1.0181E+00  1.8349E+00  5.6236E-01  7.9571E-01  8.3753E-01
             1.0434E+00
 PARAMETER:  7.6578E-02 -5.1973E-01 -1.5817E-01  4.1121E-01 -2.8878E-01  1.1792E-01  7.0699E-01 -4.7561E-01 -1.2852E-01 -7.7293E-02
             1.4250E-01
 GRADIENT:   9.5215E-01  6.4227E+00  5.4922E+00  1.3443E+01 -9.2926E+00 -5.1068E-01 -1.0429E+00 -1.8327E+00  6.1236E-01  7.3038E-01
            -1.2619E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1628.52721862826        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  9.7033E-01  3.4301E-01  9.5214E-01  1.4973E+00  7.3000E-01  1.0136E+00  2.5053E+00  8.4618E-01  7.5122E-01  8.9218E-01
             1.0428E+00
 PARAMETER:  6.9882E-02 -9.7000E-01  5.0954E-02  5.0369E-01 -2.1471E-01  1.1356E-01  1.0184E+00 -6.7023E-02 -1.8606E-01 -1.4090E-02
             1.4195E-01
 GRADIENT:  -8.9667E-01  4.6330E+00  2.6309E+00  1.2317E+01 -6.3019E+00  1.2038E-02  9.5919E-01  8.8114E-01 -1.3236E+00 -3.3092E-01
            -2.8055E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1628.63571059809        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  9.6831E-01  2.6523E-01  1.0256E+00  1.5465E+00  7.5251E-01  1.0118E+00  2.8614E+00  9.1651E-01  7.3962E-01  9.2789E-01
             1.0434E+00
 PARAMETER:  6.7794E-02 -1.2272E+00  1.2529E-01  5.3596E-01 -1.8435E-01  1.1168E-01  1.1513E+00  1.2819E-02 -2.0162E-01  2.5153E-02
             1.4251E-01
 GRADIENT:  -5.5841E-01  2.3455E+00  2.1706E+00  6.4293E+00 -2.7935E+00  3.4737E-02  2.1376E-01  1.0557E-01 -1.0320E+00 -7.1509E-01
            -4.1910E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1628.64452547115        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  9.6743E-01  2.2888E-01  1.0515E+00  1.5694E+00  7.5762E-01  1.0109E+00  3.0781E+00  9.4566E-01  7.3536E-01  9.4127E-01
             1.0437E+00
 PARAMETER:  6.6886E-02 -1.3746E+00  1.5026E-01  5.5072E-01 -1.7757E-01  1.1087E-01  1.2243E+00  4.4129E-02 -2.0739E-01  3.9472E-02
             1.4276E-01
 GRADIENT:  -4.2935E-01  1.8074E+00  1.7406E+00  5.5592E+00 -2.2928E+00  1.3756E-02  2.1945E-01  4.2250E-02 -8.3367E-01 -4.9727E-01
            -3.0416E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1628.69623087779        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1426             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6799E-01  2.2181E-01  1.0502E+00  1.5649E+00  7.5777E-01  1.0110E+00  3.1357E+00  9.4180E-01  7.3687E-01  9.4697E-01
             1.0440E+00
 PARAMETER:  6.7467E-02 -1.4060E+00  1.4899E-01  5.4783E-01 -1.7737E-01  1.1097E-01  1.2428E+00  4.0033E-02 -2.0534E-01  4.5511E-02
             1.4305E-01
 GRADIENT:   3.3815E+02  2.8837E+01  4.9723E+00  6.9114E+02  1.5057E+01  4.0756E+01  2.6626E+01  4.5034E-01  1.6099E+01  7.2371E-01
             1.2311E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1628.69804629464        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1593
 NPARAMETR:  9.6794E-01  2.2176E-01  1.0499E+00  1.5657E+00  7.5684E-01  1.0110E+00  3.1338E+00  9.3921E-01  7.3649E-01  9.4747E-01
             1.0439E+00
 PARAMETER:  6.7414E-02 -1.4062E+00  1.4869E-01  5.4832E-01 -1.7860E-01  1.1090E-01  1.2423E+00  3.7287E-02 -2.0587E-01  4.6045E-02
             1.4292E-01
 GRADIENT:   1.1603E+00  4.2765E-01  1.1349E+00 -1.1191E+01  2.0100E-03  1.1514E-01  7.5663E-01  2.5217E-03  1.0212E-01  5.4543E-02
            -9.3785E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1593
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6130E-04  2.7886E-02 -3.1703E-02 -2.1320E-02 -2.3397E-02
 SE:             2.9824E-02  1.5047E-02  1.6933E-02  2.5928E-02  2.0446E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7429E-01  6.3841E-02  6.1171E-02  4.1092E-01  2.5248E-01

 ETASHRINKSD(%)  8.4285E-02  4.9591E+01  4.3272E+01  1.3138E+01  3.1503E+01
 ETASHRINKVR(%)  1.6850E-01  7.4590E+01  6.7819E+01  2.4550E+01  5.3081E+01
 EBVSHRINKSD(%)  4.7723E-01  5.9035E+01  4.5284E+01  9.5242E+00  2.7719E+01
 EBVSHRINKVR(%)  9.5217E-01  8.3218E+01  7.0061E+01  1.8141E+01  4.7755E+01
 RELATIVEINF(%)  9.8194E+01  3.2928E+00  4.8026E+00  1.9136E+01  7.4640E+00
 EPSSHRINKSD(%)  4.5004E+01
 EPSSHRINKVR(%)  6.9755E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          399
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    733.31294949732876     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1628.6980462946385     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -895.38509679730976     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1628.698       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  2.22E-01  1.05E+00  1.57E+00  7.57E-01  1.01E+00  3.13E+00  9.39E-01  7.36E-01  9.47E-01  1.04E+00
 


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
 #CPUT: Total CPU Time in Seconds,       48.229
Stop Time:
Sun Oct 24 03:28:47 CDT 2021
