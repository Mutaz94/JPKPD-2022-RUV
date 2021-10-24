Sun Oct 24 02:39:10 CDT 2021
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
$DATA ../../../../data/SD4/S1/dat2.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1676.08274157010        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4707E+02 -2.8575E+01 -2.2775E+01 -3.6984E+00 -5.7859E+00  3.8356E+01  1.3142E+01  1.5295E+01  3.5487E+01  1.7998E+01
            -2.4081E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1681.03239187677        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.7680E-01  1.1342E+00  1.2391E+00  9.6930E-01  1.1559E+00  9.7089E-01  9.0289E-01  8.7581E-01  7.6161E-01  8.8436E-01
             1.1911E+00
 PARAMETER:  7.6523E-02  2.2594E-01  3.1441E-01  6.8818E-02  2.4486E-01  7.0458E-02 -2.1596E-03 -3.2610E-02 -1.7232E-01 -2.2888E-02
             2.7489E-01
 GRADIENT:  -6.8606E+00  1.5111E+00  1.7059E+01 -2.0231E+01  2.8972E+00 -5.1990E+00 -7.5145E+00 -3.1292E+00 -1.8672E+01 -2.2890E+01
             2.4638E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1682.12905935610        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.7134E-01  1.0921E+00  1.3785E+00  1.0001E+00  1.2299E+00  1.0279E+00  9.1371E-01  6.1518E-01  7.8501E-01  1.1386E+00
             1.1832E+00
 PARAMETER:  7.0926E-02  1.8810E-01  4.2097E-01  1.0008E-01  3.0690E-01  1.2751E-01  9.7545E-03 -3.8585E-01 -1.4206E-01  2.2977E-01
             2.6826E-01
 GRADIENT:  -1.5868E+01 -5.1792E+00  9.9466E+00 -1.8287E+01  1.0034E+01  1.6508E+01 -2.3925E+00 -1.8829E+00 -1.2399E+01 -2.1746E-02
             2.3408E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1687.03356610670        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  9.8169E-01  9.6200E-01  9.9734E-01  1.0690E+00  9.8796E-01  9.8805E-01  1.0759E+00  3.3844E-01  7.8586E-01  9.5254E-01
             1.0716E+00
 PARAMETER:  8.1521E-02  6.1258E-02  9.7336E-02  1.6672E-01  8.7887E-02  8.7977E-02  1.7314E-01 -9.8341E-01 -1.4098E-01  5.1373E-02
             1.6915E-01
 GRADIENT:   5.8294E+00 -4.8172E+00 -2.9467E+00 -5.8918E+00  6.1011E+00  6.8689E-01  8.8544E-01  1.0491E-01 -1.9558E+00  4.5826E-01
            -9.8545E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1687.12697018650        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  9.7909E-01  1.0265E+00  9.7717E-01  1.0325E+00  1.0031E+00  9.8630E-01  9.9422E-01  3.2655E-01  8.2836E-01  9.5726E-01
             1.0731E+00
 PARAMETER:  7.8866E-02  1.2612E-01  7.6909E-02  1.3197E-01  1.0307E-01  8.6210E-02  9.4202E-02 -1.0192E+00 -8.8302E-02  5.6318E-02
             1.7051E-01
 GRADIENT:  -1.2568E+00  6.8601E-02 -7.5983E-01  1.4198E+00  1.0518E+00 -2.0571E-01  1.1136E-01  8.2825E-02  1.0662E-01  2.6268E-01
            -3.3416E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1687.16985339191        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.8049E-01  1.1803E+00  8.9789E-01  9.3391E-01  1.0380E+00  9.8799E-01  8.9177E-01  2.0728E-01  8.9339E-01  9.6216E-01
             1.0742E+00
 PARAMETER:  8.0301E-02  2.6581E-01 -7.7032E-03  3.1622E-02  1.3727E-01  8.7921E-02 -1.4549E-02 -1.4737E+00 -1.2733E-02  6.1422E-02
             1.7153E-01
 GRADIENT:  -4.5372E-01  1.6539E+00 -1.7409E-01  1.7696E+00 -3.0635E-01  2.5575E-02 -1.9893E-02  4.9089E-02 -2.9295E-01  1.4307E-01
            -1.4869E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1687.18098410039        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1062             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8164E-01  1.2166E+00  8.8078E-01  9.0914E-01  1.0486E+00  9.8833E-01  8.6932E-01  9.0871E-02  9.1361E-01  9.6472E-01
             1.0742E+00
 PARAMETER:  8.1470E-02  2.9602E-01 -2.6943E-02  4.7422E-03  1.4746E-01  8.8260E-02 -4.0048E-02 -2.2983E+00  9.6477E-03  6.4083E-02
             1.7157E-01
 GRADIENT:   3.4735E+02  1.1596E+02  2.3038E+00  4.2047E+01  7.6372E+00  3.1967E+01  2.3787E+00  2.2856E-02  2.9901E+00  1.7853E-01
             9.0333E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1687.18457772342        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1219
 NPARAMETR:  9.8054E-01  1.2165E+00  8.7653E-01  9.0949E-01  1.0472E+00  9.8791E-01  8.6933E-01  4.8087E-02  9.1366E-01  9.6576E-01
             1.0752E+00
 PARAMETER:  8.0352E-02  2.9597E-01 -3.1783E-02  5.1266E-03  1.4616E-01  8.7835E-02 -4.0033E-02 -2.9347E+00  9.7008E-03  6.5157E-02
             1.7250E-01
 GRADIENT:  -8.2500E-01  4.3425E-02 -2.8149E-01  6.4679E-01  7.7088E-01 -9.4215E-02 -1.3779E-01  2.0945E-03 -3.2769E-02  6.4654E-02
             2.3093E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1687.18766886635        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1393
 NPARAMETR:  9.8147E-01  1.2161E+00  8.7228E-01  9.0907E-01  1.0441E+00  9.8838E-01  8.7199E-01  1.0000E-02  9.1272E-01  9.6196E-01
             1.0744E+00
 PARAMETER:  8.1295E-02  2.9564E-01 -3.6641E-02  4.6635E-03  1.4314E-01  8.8312E-02 -3.6974E-02 -4.5849E+00  8.6725E-03  6.1213E-02
             1.7177E-01
 GRADIENT:   1.2599E+00 -2.4311E-01 -5.0982E-02 -2.5292E-01  5.8618E-02  8.1168E-02  2.5576E-03  0.0000E+00  1.9726E-02  1.6692E-03
            -1.6340E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1687.18766886635        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1415
 NPARAMETR:  9.8147E-01  1.2161E+00  8.7228E-01  9.0907E-01  1.0441E+00  9.8838E-01  8.7199E-01  1.0000E-02  9.1272E-01  9.6196E-01
             1.0744E+00
 PARAMETER:  8.1295E-02  2.9564E-01 -3.6641E-02  4.6635E-03  1.4314E-01  8.8312E-02 -3.6974E-02 -4.5849E+00  8.6725E-03  6.1213E-02
             1.7177E-01
 GRADIENT:   1.2599E+00 -2.4311E-01 -5.0982E-02 -2.5292E-01  5.8618E-02  8.1168E-02  2.5576E-03  0.0000E+00  1.9726E-02  1.6692E-03
            -1.6340E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1415
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.3982E-04 -1.3764E-02 -3.0878E-04  3.5887E-03 -2.6073E-02
 SE:             2.9789E-02  2.0184E-02  1.4530E-04  2.3500E-02  2.3173E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9358E-01  4.9529E-01  3.3577E-02  8.7863E-01  2.6052E-01

 ETASHRINKSD(%)  2.0225E-01  3.2381E+01  9.9513E+01  2.1272E+01  2.2369E+01
 ETASHRINKVR(%)  4.0408E-01  5.4277E+01  9.9998E+01  3.8019E+01  3.9734E+01
 EBVSHRINKSD(%)  4.9974E-01  3.1741E+01  9.9548E+01  2.1937E+01  2.0831E+01
 EBVSHRINKVR(%)  9.9698E-01  5.3407E+01  9.9998E+01  3.9061E+01  3.7322E+01
 RELATIVEINF(%)  9.8555E+01  1.2385E+00  1.8277E-04  1.8765E+00  6.6582E+00
 EPSSHRINKSD(%)  4.1710E+01
 EPSSHRINKVR(%)  6.6023E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1687.1876688663544     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -952.03684230261626     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1687.188       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.22E+00  8.72E-01  9.09E-01  1.04E+00  9.88E-01  8.72E-01  1.00E-02  9.13E-01  9.62E-01  1.07E+00
 


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
 #CPUT: Total CPU Time in Seconds,       39.452
Stop Time:
Sun Oct 24 02:39:18 CDT 2021
