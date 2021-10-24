Sat Oct 23 14:11:36 CDT 2021
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
$DATA ../../../../data/SD1/A3/dat20.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17.9695671914301        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7885E+02  2.1592E+02  2.0585E+02  7.1297E+01  2.9389E+02  4.0498E+01 -2.3310E+02 -2.4068E+02 -4.8646E+01 -1.3149E+02
            -6.9982E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2345.36121353229        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9253E-01  9.7663E-01  8.0834E-01  1.1721E+00  8.3163E-01  7.4707E-01  8.9481E-01  1.4799E+00  8.2126E-01  7.2689E-01
             5.1741E+00
 PARAMETER:  9.2504E-02  7.6350E-02 -1.1277E-01  2.5884E-01 -8.4365E-02 -1.9160E-01 -1.1144E-02  4.9197E-01 -9.6920E-02 -2.1898E-01
             1.7437E+00
 GRADIENT:  -9.5694E+01  3.1833E+01 -3.4684E+01  1.0808E+02  9.7894E+00 -6.0727E+01  3.2335E+00  2.2618E+01  4.3949E+00  1.4269E+01
             8.0387E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2411.31147317985        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.2811E-01  7.4346E-01  1.2330E+00  1.2843E+00  7.9784E-01  1.0763E+00  1.3476E+00  1.3304E+00  8.6546E-01  6.9199E-02
             4.7198E+00
 PARAMETER:  2.5390E-02 -1.9644E-01  3.0949E-01  3.5024E-01 -1.2585E-01  1.7355E-01  3.9833E-01  3.8548E-01 -4.4499E-02 -2.5708E+00
             1.6518E+00
 GRADIENT:  -1.3285E+02  2.8043E+01  4.6380E+01  8.0051E+01 -8.1271E+01  4.3433E+01  6.1266E+00  4.5466E-01  2.8296E+01  6.2649E-02
             7.3024E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2619.77951788648        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      256
 NPARAMETR:  9.1213E-01  5.2182E-01  4.1417E-01  1.2477E+00  4.4771E-01  9.8616E-01  1.6952E+00  9.4310E-01  9.4573E-01  2.3817E-01
             2.9124E+00
 PARAMETER:  8.0278E-03 -5.5043E-01 -7.8148E-01  3.2133E-01 -7.0361E-01  8.6059E-02  6.2782E-01  4.1412E-02  4.4203E-02 -1.3348E+00
             1.1690E+00
 GRADIENT:  -1.5897E+02  7.1289E+00 -6.9415E+01  9.0130E+01  8.7164E+01  2.9282E+00  4.6419E+01 -2.5297E+00  2.9098E+00  1.9373E+00
             7.4494E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2635.59991049850        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  9.7713E-01  5.3384E-01  4.6028E-01  1.2037E+00  4.6609E-01  9.5012E-01  1.3751E+00  1.0556E+00  9.1513E-01  2.4604E-01
             2.8137E+00
 PARAMETER:  7.6866E-02 -5.2766E-01 -6.7593E-01  2.8536E-01 -6.6337E-01  4.8829E-02  4.1854E-01  1.5409E-01  1.1306E-02 -1.3022E+00
             1.1345E+00
 GRADIENT:  -4.0372E+00  4.0525E+00 -1.1299E+01  5.1797E+00  1.9130E+01  2.2123E+00 -8.7140E-02  9.5988E-02 -1.0523E+00  1.1395E+00
            -5.3112E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2637.59888819072        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      609
 NPARAMETR:  9.8390E-01  4.5542E-01  3.7482E-01  1.2182E+00  3.8734E-01  9.1776E-01  1.4310E+00  1.0249E+00  9.5208E-01  1.0000E-02
             2.7621E+00
 PARAMETER:  8.3765E-02 -6.8653E-01 -8.8131E-01  2.9738E-01 -8.4845E-01  1.4176E-02  4.5840E-01  1.2463E-01  5.0897E-02 -6.2241E+00
             1.1160E+00
 GRADIENT:   1.3280E+01  2.2289E+01  4.0575E+00  3.9964E+01 -1.1184E+01 -1.1473E+01 -7.5580E-01 -3.0669E+00 -6.3067E+00  0.0000E+00
            -2.5257E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2644.68296846763        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  9.7766E-01  3.2555E-01  2.4679E-01  1.1346E+00  2.7901E-01  9.4504E-01  1.4734E+00  1.0797E+00  1.1468E+00  1.0000E-02
             2.7060E+00
 PARAMETER:  7.7408E-02 -1.0223E+00 -1.2992E+00  2.2629E-01 -1.1765E+00  4.3477E-02  4.8755E-01  1.7667E-01  2.3701E-01 -1.8464E+01
             1.0955E+00
 GRADIENT:  -1.4111E+00  2.1558E-01 -6.4116E-01  9.5821E-01  2.2132E+00  2.8630E-01  1.5869E-01  1.0370E+00  1.8472E-01  0.0000E+00
             4.3836E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2644.70569720156        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      951
 NPARAMETR:  9.7823E-01  3.2189E-01  2.4307E-01  1.1303E+00  2.7582E-01  9.4419E-01  1.4747E+00  1.0773E+00  1.1550E+00  1.0000E-02
             2.7037E+00
 PARAMETER:  7.7988E-02 -1.0336E+00 -1.3144E+00  2.2247E-01 -1.1880E+00  4.2573E-02  4.8847E-01  1.7450E-01  2.4406E-01 -1.8878E+01
             1.0946E+00
 GRADIENT:   4.5382E-03  1.0143E-01 -9.3077E-02 -1.3042E-01 -1.2102E+00  5.3783E-03  9.0033E-02  1.8211E-01  6.2125E-02  0.0000E+00
             1.9099E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      951
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7632E-03  5.5532E-03 -7.2651E-03 -5.6318E-03  3.2388E-04
 SE:             2.9305E-02  2.6656E-02  2.3642E-02  2.8137E-02  5.0051E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5202E-01  8.3497E-01  7.5862E-01  8.4136E-01  5.1757E-01

 ETASHRINKSD(%)  1.8243E+00  1.0700E+01  2.0795E+01  5.7363E+00  9.8323E+01
 ETASHRINKVR(%)  3.6152E+00  2.0254E+01  3.7265E+01  1.1144E+01  9.9972E+01
 EBVSHRINKSD(%)  1.9972E+00  9.8449E+00  2.0587E+01  5.8118E+00  9.8717E+01
 EBVSHRINKVR(%)  3.9546E+00  1.8721E+01  3.6937E+01  1.1286E+01  9.9984E+01
 RELATIVEINF(%)  9.6019E+01  2.2744E+01  9.2527E+00  6.8160E+01  1.4914E-03
 EPSSHRINKSD(%)  1.9245E+01
 EPSSHRINKVR(%)  3.4787E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2644.7056972015589     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -990.61633743314815     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2644.706       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  3.22E-01  2.43E-01  1.13E+00  2.76E-01  9.44E-01  1.47E+00  1.08E+00  1.15E+00  1.00E-02  2.70E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       71.412
Stop Time:
Sat Oct 23 14:11:49 CDT 2021
