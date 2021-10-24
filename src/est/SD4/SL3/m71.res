Sun Oct 24 03:35:47 CDT 2021
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
$DATA ../../../../data/SD4/SL3/dat71.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1556.99358303154        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7153E+02 -5.1024E+01  1.2345E+01 -6.7205E+01  1.3899E+01  3.5902E+01 -1.8053E+01 -7.7829E+00 -1.3095E+01 -4.2897E+00
            -9.5443E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1573.12610535539        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.6419E-01  1.1487E+00  9.2048E-01  1.0133E+00  1.0363E+00  9.3765E-01  1.2611E+00  1.0744E+00  1.0109E+00  9.2058E-01
             1.2871E+00
 PARAMETER:  6.3538E-02  2.3866E-01  1.7138E-02  1.1322E-01  1.3566E-01  3.5618E-02  3.3200E-01  1.7172E-01  1.1081E-01  1.7244E-02
             3.5240E-01
 GRADIENT:   2.3599E+02  5.7053E+01 -1.8529E+01  7.2615E+01  4.1941E+01  1.4523E+00  1.7096E+01  8.1080E-01  1.0582E+01 -2.4395E+00
             2.6066E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1573.64544461534        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      242
 NPARAMETR:  9.6411E-01  1.2157E+00  8.6058E-01  9.7809E-01  9.8427E-01  9.7864E-01  1.1905E+00  1.1408E+00  9.9181E-01  9.1888E-01
             1.3017E+00
 PARAMETER:  6.3453E-02  2.9530E-01 -5.0150E-02  7.7843E-02  8.4141E-02  7.8406E-02  2.7440E-01  2.3173E-01  9.1774E-02  1.5399E-02
             3.6364E-01
 GRADIENT:   1.0195E+01  3.4134E+01  2.3363E+00  3.1630E+01 -2.2143E+01 -5.7722E+00  4.3643E+00  3.7887E+00  2.0734E+00  6.6391E+00
             2.9554E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1575.95206837710        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      421
 NPARAMETR:  9.5863E-01  1.1660E+00  7.4569E-01  9.6390E-01  9.1706E-01  9.8932E-01  1.1799E+00  9.0561E-01  9.7001E-01  8.1093E-01
             1.2075E+00
 PARAMETER:  5.7750E-02  2.5360E-01 -1.9345E-01  6.3233E-02  1.3415E-02  8.9261E-02  2.6545E-01  8.5756E-04  6.9552E-02 -1.0957E-01
             2.8854E-01
 GRADIENT:  -2.1365E+00 -4.2922E+00 -1.3307E+00 -3.8014E-01  3.1028E+00 -2.0229E+00 -1.4900E+00 -9.2579E-02 -2.7719E-01  1.1978E-01
            -4.3050E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1576.15956975127        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      601
 NPARAMETR:  9.6017E-01  1.3812E+00  4.6729E-01  8.0127E-01  8.3953E-01  9.9346E-01  1.0486E+00  6.5820E-01  1.0341E+00  6.6102E-01
             1.2096E+00
 PARAMETER:  5.9350E-02  4.2294E-01 -6.6080E-01 -1.2156E-01 -7.4909E-02  9.3437E-02  1.4743E-01 -3.1825E-01  1.3353E-01 -3.1397E-01
             2.9026E-01
 GRADIENT:  -7.1030E-01  2.3256E+00  2.7084E+00 -4.0481E+00 -6.7597E+00 -6.1881E-01  1.6356E+00  2.4712E-01  1.8251E-01  6.6903E-01
             7.6372E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1576.28961140288        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  9.6058E-01  1.4993E+00  3.1015E-01  7.0162E-01  7.8276E-01  9.9811E-01  9.5615E-01  4.4058E-01  1.0818E+00  5.5222E-01
             1.2005E+00
 PARAMETER:  5.9782E-02  5.0499E-01 -1.0707E+00 -2.5437E-01 -1.4493E-01  9.8111E-02  5.5158E-02 -7.1966E-01  1.7863E-01 -4.9380E-01
             2.8275E-01
 GRADIENT:   6.9091E+00  1.0461E+01  2.6348E+00 -5.6608E-01 -1.0331E+01  2.0518E+00  9.4298E-01  3.0012E-01 -7.2048E-01  7.7400E-01
            -7.5335E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1576.37986916806        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      956
 NPARAMETR:  9.5637E-01  1.5470E+00  2.8041E-01  6.6306E-01  7.9271E-01  9.9155E-01  9.2693E-01  3.6187E-01  1.1259E+00  5.3630E-01
             1.2041E+00
 PARAMETER:  5.5394E-02  5.3630E-01 -1.1715E+00 -3.1089E-01 -1.3230E-01  9.1512E-02  2.4125E-02 -9.1647E-01  2.1857E-01 -5.2306E-01
             2.8570E-01
 GRADIENT:  -5.3844E-01  5.4770E-01 -7.6273E-02  2.4693E-01 -1.9089E-01 -1.5750E-01  3.7582E-01  2.0545E-01  3.0203E-01 -2.6251E-01
             3.3841E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1576.44316483003        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1136
 NPARAMETR:  9.5709E-01  1.5781E+00  2.7073E-01  6.4205E-01  8.0636E-01  9.9192E-01  9.0889E-01  8.9580E-02  1.1454E+00  5.5973E-01
             1.2026E+00
 PARAMETER:  5.6138E-02  5.5622E-01 -1.2066E+00 -3.4309E-01 -1.1523E-01  9.1886E-02  4.4689E-03 -2.3126E+00  2.3578E-01 -4.8031E-01
             2.8452E-01
 GRADIENT:   1.4074E+00 -2.3222E+00  6.2153E-01 -1.7336E+00  8.3267E-01  4.5512E-02 -1.1673E+00  6.6436E-03 -7.3420E-01 -3.9133E-01
            -3.2048E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1576.45187637520        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1312
 NPARAMETR:  9.5655E-01  1.5804E+00  2.6802E-01  6.4142E-01  8.0486E-01  9.9147E-01  9.1222E-01  4.5549E-02  1.1499E+00  5.5919E-01
             1.2022E+00
 PARAMETER:  5.5574E-02  5.5770E-01 -1.2167E+00 -3.4408E-01 -1.1708E-01  9.1436E-02  8.1281E-03 -2.9890E+00  2.3970E-01 -4.8127E-01
             2.8412E-01
 GRADIENT:   4.2503E-01 -5.3582E-01 -3.2160E-01  6.6488E-01  1.1991E+00 -1.1903E-01  2.9674E-01  2.2507E-03  1.2139E-01 -9.5381E-02
            -1.3659E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1576.45387523395        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1491
 NPARAMETR:  9.5657E-01  1.5811E+00  2.6506E-01  6.3928E-01  8.0294E-01  9.9186E-01  9.0954E-01  1.0000E-02  1.1505E+00  5.5694E-01
             1.2022E+00
 PARAMETER:  5.5595E-02  5.5815E-01 -1.2278E+00 -3.4742E-01 -1.1947E-01  9.1822E-02  5.1877E-03 -4.9027E+00  2.4017E-01 -4.8530E-01
             2.8413E-01
 GRADIENT:   8.8792E-01 -1.5420E+00 -2.0450E-01 -3.2158E-01  6.0627E-01  8.8764E-02 -7.1938E-03  0.0000E+00  1.9684E-02  1.7489E-02
             1.6277E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1576.45387523395        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1513
 NPARAMETR:  9.5657E-01  1.5811E+00  2.6506E-01  6.3928E-01  8.0294E-01  9.9186E-01  9.0954E-01  1.0000E-02  1.1505E+00  5.5694E-01
             1.2022E+00
 PARAMETER:  5.5595E-02  5.5815E-01 -1.2278E+00 -3.4742E-01 -1.1947E-01  9.1822E-02  5.1877E-03 -4.9027E+00  2.4017E-01 -4.8530E-01
             2.8413E-01
 GRADIENT:   8.8792E-01 -1.5420E+00 -2.0450E-01 -3.2158E-01  6.0627E-01  8.8764E-02 -7.1938E-03  0.0000E+00  1.9684E-02  1.7489E-02
             1.6277E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1513
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8727E-04 -1.2954E-02 -1.9730E-04  1.2714E-02 -2.2640E-02
 SE:             2.9799E-02  2.6654E-02  1.1161E-04  2.4740E-02  1.5977E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9231E-01  6.2695E-01  7.7107E-02  6.0731E-01  1.5647E-01

 ETASHRINKSD(%)  1.6937E-01  1.0706E+01  9.9626E+01  1.7119E+01  4.6475E+01
 ETASHRINKVR(%)  3.3846E-01  2.0265E+01  9.9999E+01  3.1307E+01  7.1350E+01
 EBVSHRINKSD(%)  5.8675E-01  1.0940E+01  9.9644E+01  1.6861E+01  4.7341E+01
 EBVSHRINKVR(%)  1.1701E+00  2.0683E+01  9.9999E+01  3.0880E+01  7.2270E+01
 RELATIVEINF(%)  9.7151E+01  6.9484E+00  1.1127E-04  4.6644E+00  2.7890E+00
 EPSSHRINKSD(%)  4.3602E+01
 EPSSHRINKVR(%)  6.8193E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1576.4538752339495     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -841.30304867021130     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1576.454       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.58E+00  2.65E-01  6.39E-01  8.03E-01  9.92E-01  9.10E-01  1.00E-02  1.15E+00  5.57E-01  1.20E+00
 


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
 #CPUT: Total CPU Time in Seconds,       41.556
Stop Time:
Sun Oct 24 03:35:56 CDT 2021
