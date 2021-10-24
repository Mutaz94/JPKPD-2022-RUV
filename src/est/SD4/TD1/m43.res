Sun Oct 24 03:46:59 CDT 2021
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
$DATA ../../../../data/SD4/TD1/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1688.71735426303        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4009E+02 -3.9071E+00 -1.4490E+01  4.0512E+01  2.5367E+01  4.0289E+01  1.1654E+01  9.8648E+00  4.8445E+01 -4.5776E+00
             4.0637E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1694.84047394922        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7398E-01  9.8814E-01  9.5545E-01  1.0247E+00  9.6846E-01  9.4076E-01  9.4882E-01  9.2318E-01  8.0815E-01  1.0571E+00
             9.8695E-01
 PARAMETER:  7.3640E-02  8.8074E-02  5.4426E-02  1.2437E-01  6.7956E-02  3.8935E-02  4.7468E-02  2.0071E-02 -1.1301E-01  1.5549E-01
             8.6869E-02
 GRADIENT:   3.8620E+02  1.8828E+01 -1.8005E+01  8.3689E+01  1.6951E+01  1.8885E+01 -3.6632E-01  1.0225E+01  1.3355E+01  7.3067E+00
             3.5829E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1700.06942188146        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      204
 NPARAMETR:  9.6187E-01  9.9360E-01  8.8068E-01  1.0132E+00  9.4696E-01  1.0035E+00  1.0216E+00  4.7072E-01  7.5000E-01  1.0496E+00
             8.9932E-01
 PARAMETER:  6.1126E-02  9.3581E-02 -2.7060E-02  1.1315E-01  4.5505E-02  1.0353E-01  1.2135E-01 -6.5349E-01 -1.8768E-01  1.4842E-01
            -6.1109E-03
 GRADIENT:  -4.7634E+01 -1.1931E+01 -1.5786E+01 -3.9053E+00  2.3709E+01 -2.5406E+00 -5.0227E+00  2.1042E+00  3.1806E-01  6.0264E-01
            -2.8207E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1703.02089051765        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      381
 NPARAMETR:  9.8289E-01  8.3269E-01  8.3488E-01  1.1187E+00  8.2652E-01  1.0092E+00  1.3249E+00  1.4314E-01  6.6351E-01  9.8656E-01
             8.9928E-01
 PARAMETER:  8.2744E-02 -8.3090E-02 -8.0468E-02  2.1215E-01 -9.0527E-02  1.0915E-01  3.8134E-01 -1.8440E+00 -3.1021E-01  8.6466E-02
            -6.1660E-03
 GRADIENT:   5.5858E-02  2.0090E+01 -2.6783E+00  1.6548E+01 -5.0362E+00  5.3825E-01  4.5679E+00  2.2517E-01 -1.1488E+00  5.2580E+00
            -2.5351E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1706.38242965058        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  9.8107E-01  4.4974E-01  8.4460E-01  1.3310E+00  7.0404E-01  1.0122E+00  1.9139E+00  1.0000E-02  5.9793E-01  9.1005E-01
             9.0907E-01
 PARAMETER:  8.0889E-02 -6.9909E-01 -6.8888E-02  3.8595E-01 -2.5092E-01  1.1216E-01  7.4913E-01 -5.9945E+00 -4.1428E-01  5.7441E-03
             4.6639E-03
 GRADIENT:   4.5708E+00  3.3891E+00 -8.6928E+00  2.1481E+01  9.0543E+00  2.6558E+00 -4.4568E+00  0.0000E+00 -4.0270E+00 -1.9647E+00
             5.6697E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1707.73010906145        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      733
 NPARAMETR:  9.6850E-01  1.7694E-01  9.3199E-01  1.4922E+00  6.7743E-01  9.8792E-01  3.3425E+00  1.0000E-02  5.9539E-01  9.7119E-01
             9.1758E-01
 PARAMETER:  6.7996E-02 -1.6320E+00  2.9568E-02  5.0024E-01 -2.8946E-01  8.7843E-02  1.3067E+00 -1.2782E+01 -4.1853E-01  7.0771E-02
             1.3985E-02
 GRADIENT:  -9.3354E+00  4.6437E+00  4.3540E+00  1.5208E+01 -7.7531E+00 -4.7378E+00  3.7175E+00  0.0000E+00  4.8987E-01 -7.4156E-01
             2.0317E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1708.02710064087        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      917
 NPARAMETR:  9.7182E-01  1.2961E-01  9.3816E-01  1.5106E+00  6.7355E-01  9.9880E-01  3.7670E+00  1.0000E-02  5.8786E-01  9.7781E-01
             9.1426E-01
 PARAMETER:  7.1414E-02 -1.9432E+00  3.6165E-02  5.1247E-01 -2.9520E-01  9.8803E-02  1.4263E+00 -1.5120E+01 -4.3127E-01  7.7556E-02
             1.0365E-02
 GRADIENT:   1.0860E+00  3.3050E-02  7.7313E-01 -5.4426E+00 -4.4731E-01  9.2973E-02 -4.6492E-01  0.0000E+00 -1.5155E-01  9.1282E-02
             3.3225E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1708.03979207881        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1111             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7196E-01  1.2995E-01  9.3678E-01  1.5086E+00  6.7311E-01  9.9870E-01  3.7875E+00  1.0000E-02  5.8933E-01  9.7697E-01
             9.1331E-01
 PARAMETER:  7.1561E-02 -1.9406E+00  3.4697E-02  5.1119E-01 -2.9585E-01  9.8702E-02  1.4317E+00 -1.5120E+01 -4.2877E-01  7.6703E-02
             9.3197E-03
 GRADIENT:   4.5273E+02  3.1073E+01  5.0834E+00  8.9915E+02  3.7380E+01  4.8173E+01  6.4819E+01  0.0000E+00  2.7662E+01  8.1270E-01
             5.1526E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1708.04055980205        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     1243
 NPARAMETR:  9.7196E-01  1.2975E-01  9.3645E-01  1.5086E+00  6.7293E-01  9.9866E-01  3.7850E+00  1.0000E-02  5.8958E-01  9.7815E-01
             9.1367E-01
 PARAMETER:  7.1554E-02 -1.9422E+00  3.4346E-02  5.1121E-01 -2.9611E-01  9.8654E-02  1.4310E+00 -1.5120E+01 -4.2835E-01  7.7906E-02
             9.7133E-03
 GRADIENT:   1.4960E+00  3.6794E-01  5.1454E-01 -1.3309E+01  6.6860E-01  1.1310E-01  1.0784E+00  0.0000E+00  2.7241E-01  1.1745E-01
             5.1890E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1243
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1306E-04  3.3916E-02 -3.6014E-04 -2.4923E-02 -3.2701E-03
 SE:             2.9925E-02  1.5478E-02  2.1969E-04  2.5657E-02  2.5595E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8099E-01  2.8436E-02  1.0116E-01  3.3134E-01  8.9834E-01

 ETASHRINKSD(%)  1.0000E-10  4.8147E+01  9.9264E+01  1.4047E+01  1.4253E+01
 ETASHRINKVR(%)  1.0000E-10  7.3112E+01  9.9995E+01  2.6120E+01  2.6474E+01
 EBVSHRINKSD(%)  3.4181E-01  6.0549E+01  9.9296E+01  8.9012E+00  1.0038E+01
 EBVSHRINKVR(%)  6.8246E-01  8.4436E+01  9.9995E+01  1.7010E+01  1.9068E+01
 RELATIVEINF(%)  9.8887E+01  4.3455E+00  4.4463E-04  2.3372E+01  7.4388E+00
 EPSSHRINKSD(%)  4.3550E+01
 EPSSHRINKVR(%)  6.8134E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1708.0405598020463     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -972.88973323830817     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1708.041       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.72E-01  1.30E-01  9.36E-01  1.51E+00  6.73E-01  9.99E-01  3.79E+00  1.00E-02  5.90E-01  9.78E-01  9.14E-01
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       36.532
Stop Time:
Sun Oct 24 03:47:08 CDT 2021
