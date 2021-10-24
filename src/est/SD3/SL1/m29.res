Sat Oct 23 23:24:48 CDT 2021
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
$DATA ../../../../data/SD3/SL1/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2026.85601292034        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8051E+02 -1.7123E+01 -3.1643E+01  3.8984E+01  7.1497E+01  1.7943E+01 -2.7062E+01  4.7606E+00 -1.8216E+01 -9.6323E+00
            -8.0596E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2038.19841860119        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8670E-01  1.0469E+00  1.0337E+00  1.0244E+00  9.9984E-01  1.0827E+00  1.1219E+00  9.7881E-01  1.0874E+00  1.0086E+00
             1.0889E+00
 PARAMETER:  8.6611E-02  1.4581E-01  1.3312E-01  1.2409E-01  9.9841E-02  1.7947E-01  2.1499E-01  7.8580E-02  1.8375E-01  1.0856E-01
             1.8513E-01
 GRADIENT:   5.4669E+00 -1.7668E+00 -1.5795E+01  2.6242E+01  1.9894E+01  2.1313E+00 -1.4980E+01  4.0540E+00  6.3829E-01 -3.3293E+00
            -7.1078E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2040.42051569273        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.8778E-01  1.0850E+00  1.0427E+00  9.9188E-01  1.0219E+00  1.0671E+00  1.3767E+00  8.7532E-01  1.0338E+00  1.0721E+00
             1.0773E+00
 PARAMETER:  8.7705E-02  1.8159E-01  1.4181E-01  9.1844E-02  1.2167E-01  1.6495E-01  4.1966E-01 -3.3163E-02  1.3327E-01  1.6966E-01
             1.7442E-01
 GRADIENT:   8.3029E+00  4.5755E+00 -3.2070E+00  9.6530E+00  2.7703E+00 -3.9321E+00  3.0447E+00  2.4624E+00  3.4823E+00  2.9511E+00
            -1.5015E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2041.45279550056        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.8564E-01  1.2872E+00  7.4822E-01  8.5111E-01  9.6431E-01  1.0851E+00  1.2418E+00  4.3409E-01  1.0635E+00  9.6379E-01
             1.1038E+00
 PARAMETER:  8.5538E-02  3.5250E-01 -1.9006E-01 -6.1210E-02  6.3658E-02  1.8168E-01  3.1653E-01 -7.3450E-01  1.6152E-01  6.3121E-02
             1.9875E-01
 GRADIENT:  -1.2806E+00  1.1474E+01  3.3449E+00  5.5703E+00 -1.1361E+01  1.8340E+00  4.0412E-01  3.5415E-01 -1.2025E+00  1.5601E+00
             6.2855E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2041.67512645276        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.8706E-01  1.3939E+00  6.7802E-01  7.7323E-01  9.9776E-01  1.0807E+00  1.1728E+00  2.6702E-01  1.1265E+00  9.6737E-01
             1.0942E+00
 PARAMETER:  8.6972E-02  4.3211E-01 -2.8858E-01 -1.5718E-01  9.7762E-02  1.7761E-01  2.5935E-01 -1.2204E+00  2.1908E-01  6.6828E-02
             1.9005E-01
 GRADIENT:   3.4180E-01 -3.9914E-01 -6.3375E-01 -8.1060E-02  7.1320E-01 -1.5545E-01  1.4335E-01  1.7824E-01  3.2316E-01  1.4897E-01
            -4.2464E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2041.72994578123        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.8668E-01  1.4071E+00  6.6248E-01  7.6472E-01  9.9541E-01  1.0822E+00  1.1658E+00  7.3838E-02  1.1302E+00  9.6823E-01
             1.0964E+00
 PARAMETER:  8.6587E-02  4.4152E-01 -3.1177E-01 -1.6824E-01  9.5397E-02  1.7898E-01  2.5340E-01 -2.5059E+00  2.2235E-01  6.7710E-02
             1.9201E-01
 GRADIENT:  -6.4064E-01  8.0735E-01 -1.5914E-01 -2.1367E-01 -1.2169E+00  3.4622E-01 -1.3078E-02  1.3285E-02 -9.9651E-02  5.1109E-01
             1.2705E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2041.76506930835        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1067
 NPARAMETR:  9.8780E-01  1.3795E+00  6.7707E-01  7.8060E-01  9.8926E-01  1.0828E+00  1.1807E+00  1.8558E-02  1.1190E+00  9.6372E-01
             1.0948E+00
 PARAMETER:  8.7725E-02  4.2175E-01 -2.8998E-01 -1.4769E-01  8.9202E-02  1.7959E-01  2.6607E-01 -3.8868E+00  2.1245E-01  6.3050E-02
             1.9053E-01
 GRADIENT:   1.8259E+00 -8.6883E-01  4.7779E-01 -1.2388E+00 -3.9704E-01  6.4362E-01 -1.1316E-01  6.6944E-04 -1.2445E-01 -2.8623E-01
            -2.2250E-01

0ITERATION NO.:   34    OBJECTIVE VALUE:  -2041.76626276268        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     1199
 NPARAMETR:  9.8762E-01  1.3792E+00  6.7663E-01  7.8104E-01  9.8945E-01  1.0823E+00  1.1813E+00  1.0000E-02  1.1195E+00  9.6562E-01
             1.0949E+00
 PARAMETER:  8.7544E-02  4.2148E-01 -2.9063E-01 -1.4713E-01  8.9397E-02  1.7909E-01  2.6658E-01 -4.8549E+00  2.1293E-01  6.5011E-02
             1.9070E-01
 GRADIENT:   1.4752E+00 -1.1544E+00 -2.5176E-01 -4.7829E-01  4.5355E-01  4.5102E-01  2.4765E-02  0.0000E+00  9.6321E-02  1.8648E-02
             4.8967E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1199
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9278E-04 -1.1871E-02 -3.4723E-04  1.0559E-02 -2.2615E-02
 SE:             2.9830E-02  2.4333E-02  1.3368E-04  2.3149E-02  2.2126E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8949E-01  6.2565E-01  9.3912E-03  6.4830E-01  3.0675E-01

 ETASHRINKSD(%)  6.4770E-02  1.8480E+01  9.9552E+01  2.2447E+01  2.5874E+01
 ETASHRINKVR(%)  1.2950E-01  3.3545E+01  9.9998E+01  3.9855E+01  4.5053E+01
 EBVSHRINKSD(%)  3.7429E-01  1.7907E+01  9.9610E+01  2.4020E+01  2.4516E+01
 EBVSHRINKVR(%)  7.4717E-01  3.2608E+01  9.9998E+01  4.2270E+01  4.3022E+01
 RELATIVEINF(%)  9.9086E+01  5.9162E+00  2.0086E-04  4.5688E+00  9.7771E+00
 EPSSHRINKSD(%)  3.3307E+01
 EPSSHRINKVR(%)  5.5520E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2041.7662627626833     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1122.8277295580106     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2041.766       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  1.38E+00  6.77E-01  7.81E-01  9.89E-01  1.08E+00  1.18E+00  1.00E-02  1.12E+00  9.66E-01  1.09E+00
 


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
 #CPUT: Total CPU Time in Seconds,       96.847
Stop Time:
Sat Oct 23 23:25:03 CDT 2021
