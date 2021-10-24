Sun Oct 24 02:14:04 CDT 2021
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
$DATA ../../../../data/SD4/A2/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1135.59733043072        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4533E+02  7.8237E+01  2.7707E+01  1.1294E+02  9.7844E+01  9.0750E+01 -3.2804E+01 -1.5457E+01 -3.4523E+01 -9.0540E+01
            -9.2183E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1439.49405694115        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      108
 NPARAMETR:  1.0012E+00  9.6034E-01  9.5797E-01  9.9907E-01  9.2069E-01  8.6722E-01  1.0468E+00  1.0005E+00  1.0557E+00  1.1487E+00
             2.3794E+00
 PARAMETER:  1.0118E-01  5.9531E-02  5.7064E-02  9.9067E-02  1.7372E-02 -4.2461E-02  1.4572E-01  1.0049E-01  1.5422E-01  2.3865E-01
             9.6683E-01
 GRADIENT:  -4.1227E+01 -6.9931E+00 -5.2339E+00 -1.6964E+01 -2.7667E+00  1.1242E+01 -1.8173E+00  6.5557E+00 -1.9731E-01 -7.1366E+00
             2.8588E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1442.45263439820        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      285
 NPARAMETR:  1.0069E+00  7.8921E-01  1.3106E+00  1.1581E+00  9.9064E-01  8.8417E-01  1.1719E+00  7.2764E-01  9.7431E-01  1.3185E+00
             2.3943E+00
 PARAMETER:  1.0691E-01 -1.3672E-01  3.7047E-01  2.4678E-01  9.0596E-02 -2.3110E-02  2.5862E-01 -2.1795E-01  7.3978E-02  3.7652E-01
             9.7308E-01
 GRADIENT:  -1.9780E+01  2.1228E+01  6.7334E+00  2.7375E+01 -1.6787E+01  1.8406E+01 -6.6020E-01  1.4606E+00 -3.1845E+00 -1.3274E+00
             2.7699E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1445.82712218121        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      465
 NPARAMETR:  1.0106E+00  6.3449E-01  1.0824E+00  1.2346E+00  8.5985E-01  8.3405E-01  1.4210E+00  3.1066E-01  9.2506E-01  1.2508E+00
             2.2076E+00
 PARAMETER:  1.1055E-01 -3.5493E-01  1.7922E-01  3.1074E-01 -5.1002E-02 -8.1462E-02  4.5139E-01 -1.0690E+00  2.2107E-02  3.2377E-01
             8.9191E-01
 GRADIENT:  -3.7807E+00  1.3616E+01  2.4971E+00  2.5146E+01 -7.0729E+00 -2.3067E+00 -1.8482E+00  4.2315E-01 -5.3598E+00  3.3105E-02
            -9.6654E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1448.13067213851        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.0117E+00  3.2021E-01  1.0412E+00  1.3871E+00  7.5084E-01  8.2420E-01  1.6705E+00  1.4922E-02  8.8883E-01  1.2047E+00
             2.2632E+00
 PARAMETER:  1.1159E-01 -1.0388E+00  1.4037E-01  4.2725E-01 -1.8656E-01 -9.3344E-02  6.1310E-01 -4.1049E+00 -1.7852E-02  2.8620E-01
             9.1680E-01
 GRADIENT:   1.3047E+01 -8.8646E-02 -1.9150E+00 -1.5038E+01  2.3942E-01 -5.0070E+00 -1.2220E-01  1.4783E-03 -7.3687E-02  3.2099E+00
             7.1683E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1449.18340909677        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.0017E+00  1.3775E-01  1.2117E+00  1.5165E+00  7.8196E-01  8.3799E-01  1.7071E+00  1.0000E-02  8.4476E-01  1.2361E+00
             2.2311E+00
 PARAMETER:  1.0170E-01 -1.8823E+00  2.9202E-01  5.1637E-01 -1.4595E-01 -7.6755E-02  6.3482E-01 -8.2846E+00 -6.8706E-02  3.1197E-01
             9.0249E-01
 GRADIENT:  -4.1651E+00  7.6685E-01 -1.1713E+00  1.0743E+00  4.3231E+00  1.4481E+00 -3.1272E-02  0.0000E+00 -8.0282E-01 -1.4494E+00
            -3.3666E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1449.64145245535        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      990
 NPARAMETR:  1.0004E+00  3.2686E-02  1.1405E+00  1.5709E+00  7.2063E-01  8.3521E-01  1.2469E+00  1.0000E-02  8.2095E-01  1.1991E+00
             2.2356E+00
 PARAMETER:  1.0040E-01 -3.3208E+00  2.3149E-01  5.5167E-01 -2.2764E-01 -8.0066E-02  3.2068E-01 -1.6775E+01 -9.7289E-02  2.8159E-01
             9.0450E-01
 GRADIENT:  -2.1872E+00  2.8629E-01  1.8245E+00  5.8614E+00 -4.3026E+00  6.8605E-01 -3.4917E-05  0.0000E+00 -5.4749E-01  7.3343E-02
            -1.3537E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1449.73633007248        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  1.0007E+00  1.0000E-02  1.1695E+00  1.5820E+00  7.3148E-01  8.3316E-01  8.3598E-01  1.0000E-02  8.1561E-01  1.2114E+00
             2.2406E+00
 PARAMETER:  1.0068E-01 -4.5915E+00  2.5654E-01  5.5866E-01 -2.1268E-01 -8.2531E-02 -7.9153E-02 -2.4568E+01 -1.0382E-01  2.9174E-01
             9.0676E-01
 GRADIENT:   7.1285E-01  0.0000E+00 -4.4999E-01 -4.8126E+00  1.2095E+00  8.2760E-02  4.8486E-05  0.0000E+00  4.2200E-01 -2.7321E-02
             2.1464E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1449.73692376718        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1302
 NPARAMETR:  1.0003E+00  1.0000E-02  1.1693E+00  1.5823E+00  7.3118E-01  8.3291E-01  8.3526E-01  1.0000E-02  8.1453E-01  1.2108E+00
             2.2399E+00
 PARAMETER:  1.0029E-01 -4.5915E+00  2.5638E-01  5.5891E-01 -2.1310E-01 -8.2833E-02 -8.0008E-02 -2.4568E+01 -1.0514E-01  2.9125E-01
             9.0645E-01
 GRADIENT:  -5.2033E-01  0.0000E+00 -3.4763E-01 -3.9073E+00  9.3274E-01 -4.5717E-02  3.9295E-05  0.0000E+00  1.9067E-02 -7.7740E-02
            -3.5139E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1302
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.6493E-04 -2.5637E-04 -2.4354E-06 -1.2837E-02 -3.1100E-02
 SE:             2.9097E-02  1.2806E-04  1.0935E-04  2.7412E-02  2.2332E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8725E-01  4.5291E-02  9.8223E-01  6.3958E-01  1.6372E-01

 ETASHRINKSD(%)  2.5223E+00  9.9571E+01  9.9634E+01  8.1654E+00  2.5186E+01
 ETASHRINKVR(%)  4.9810E+00  9.9998E+01  9.9999E+01  1.5664E+01  4.4029E+01
 EBVSHRINKSD(%)  2.5857E+00  9.9594E+01  9.9595E+01  7.3493E+00  2.3365E+01
 EBVSHRINKVR(%)  5.1045E+00  9.9998E+01  9.9998E+01  1.4158E+01  4.1271E+01
 RELATIVEINF(%)  8.8263E+01  6.8925E-05  1.3634E-04  4.7517E+00  3.4031E+00
 EPSSHRINKSD(%)  3.4434E+01
 EPSSHRINKVR(%)  5.7011E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1449.7369237671794     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -714.58609720344123     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1449.737       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.17E+00  1.58E+00  7.31E-01  8.33E-01  8.35E-01  1.00E-02  8.15E-01  1.21E+00  2.24E+00
 


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
 #CPUT: Total CPU Time in Seconds,       38.798
Stop Time:
Sun Oct 24 02:14:13 CDT 2021
