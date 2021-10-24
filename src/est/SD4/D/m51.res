Sun Oct 24 04:20:17 CDT 2021
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
$DATA ../../../../data/SD4/D/dat51.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1081.76801839708        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6144E+02 -1.1528E+01  4.7107E+01  1.2842E+01  1.1565E+02  3.0322E+01 -1.9860E+01 -2.0289E+02 -6.2028E+01 -3.4029E+01
            -7.6920E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1570.90501067497        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0814E+00  9.8526E-01  9.8521E-01  9.4723E-01  9.1764E-01  1.1649E+00  1.0577E+00  1.1507E+00  1.0980E+00  1.0116E+00
             1.2096E+00
 PARAMETER:  1.7823E-01  8.5151E-02  8.5103E-02  4.5790E-02  1.4055E-02  2.5263E-01  1.5611E-01  2.4039E-01  1.9345E-01  1.1149E-01
             2.9030E-01
 GRADIENT:   6.9457E+02 -3.7079E+01  1.1086E+01 -5.6387E+01 -1.2040E+01  1.1981E+02 -7.4629E+00  4.2654E+00 -7.8456E+00  2.6903E+00
             6.2797E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1580.32796436397        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0446E+00  9.9378E-01  6.7520E-01  9.2766E-01  8.1724E-01  1.1123E+00  1.5035E+00  7.2820E-01  9.8742E-01  1.0514E+00
             9.6580E-01
 PARAMETER:  1.4362E-01  9.3766E-02 -2.9274E-01  2.4912E-02 -1.0182E-01  2.0639E-01  5.0782E-01 -2.1718E-01  8.7343E-02  1.5014E-01
             6.5205E-02
 GRADIENT:   7.6402E+02 -1.3379E+01 -2.1622E+01  1.1678E+01  2.6354E+01  1.3979E+02  5.3180E+01  7.5510E+00 -4.0939E+00  2.8868E+01
            -1.5719E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1590.00713818487        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      342
 NPARAMETR:  1.0100E+00  8.0011E-01  9.1125E-01  1.1064E+00  8.2209E-01  1.0822E+00  1.5655E+00  8.6523E-01  1.0183E+00  9.1291E-01
             9.5123E-01
 PARAMETER:  1.0991E-01 -1.2301E-01  7.0606E-03  2.0108E-01 -9.5900E-02  1.7898E-01  5.4822E-01 -4.4764E-02  1.1816E-01  8.8822E-03
             5.0004E-02
 GRADIENT:   6.7392E+01 -1.4333E+01  5.0170E+00 -2.5717E+01 -1.1167E+00  7.0183E+00 -2.8611E-01  3.9144E-01  3.2522E+00 -6.1706E+00
            -2.7756E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1592.71270473308        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.7181E-01  8.5619E-01  8.8437E-01  1.0910E+00  8.3279E-01  1.0545E+00  1.5080E+00  7.8685E-01  1.0142E+00  9.7910E-01
             1.0116E+00
 PARAMETER:  7.1402E-02 -5.5259E-02 -2.2878E-02  1.8710E-01 -8.2972E-02  1.5305E-01  5.1075E-01 -1.3972E-01  1.1406E-01  7.8877E-02
             1.1156E-01
 GRADIENT:  -8.2871E+00 -3.3293E+00 -1.8974E+00 -2.3234E+00 -7.2041E-01 -9.9975E-01  1.1134E+00  9.4427E-01 -5.0461E-02  2.0018E+00
             5.3420E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1593.74508307124        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  9.8152E-01  1.1848E+00  6.5658E-01  8.8411E-01  8.5130E-01  1.0617E+00  1.2003E+00  4.3942E-01  1.1659E+00  8.9524E-01
             1.0123E+00
 PARAMETER:  8.1345E-02  2.6955E-01 -3.2072E-01 -2.3168E-02 -6.0987E-02  1.5989E-01  2.8259E-01 -7.2230E-01  2.5351E-01 -1.0668E-02
             1.1222E-01
 GRADIENT:   5.8765E+00  1.0028E+01  3.3649E+00  1.1239E+01 -2.3798E+00  6.1037E-01 -2.0640E-01 -1.9784E-01 -3.8155E-01 -3.0642E+00
             1.4476E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1594.17094748668        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  9.7815E-01  1.4318E+00  5.4437E-01  7.1640E-01  9.2155E-01  1.0616E+00  1.0205E+00  2.4002E-01  1.3396E+00  9.4678E-01
             1.0083E+00
 PARAMETER:  7.7904E-02  4.5892E-01 -5.0813E-01 -2.3351E-01  1.8302E-02  1.5982E-01  1.2030E-01 -1.3270E+00  3.9236E-01  4.5314E-02
             1.0822E-01
 GRADIENT:  -1.4026E+00  2.6228E+00  2.7565E-01  1.8794E+00 -1.8222E+00  4.2813E-01 -1.5876E+00  1.2857E-01 -1.2424E+00  5.6557E-01
            -4.7294E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1594.23411764488        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1051
 NPARAMETR:  9.7910E-01  1.4463E+00  5.3579E-01  7.0483E-01  9.2760E-01  1.0605E+00  1.0242E+00  1.0180E-01  1.3620E+00  9.4603E-01
             1.0098E+00
 PARAMETER:  7.8875E-02  4.6901E-01 -5.2401E-01 -2.4980E-01  2.4844E-02  1.5873E-01  1.2387E-01 -2.1847E+00  4.0894E-01  4.4514E-02
             1.0970E-01
 GRADIENT:   5.1109E-01 -6.5086E-02 -8.1948E-02  1.3044E+00  9.9762E-01 -1.8259E-02  1.9854E-01  1.8586E-02  1.3217E-01  2.7769E-01
             2.2383E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1594.24912810870        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1208
 NPARAMETR:  9.7939E-01  1.4496E+00  5.2967E-01  6.9925E-01  9.2586E-01  1.0629E+00  1.0201E+00  1.5467E-02  1.3656E+00  9.4001E-01
             1.0091E+00
 PARAMETER:  7.9170E-02  4.7130E-01 -5.3551E-01 -2.5775E-01  2.2966E-02  1.6105E-01  1.1992E-01 -4.0690E+00  4.1158E-01  3.8135E-02
             1.0909E-01
 GRADIENT:   4.2257E+02  3.5711E+02  8.8839E+00  9.8561E+01  8.9477E+00  1.1254E+02  8.1214E+00  1.8320E-03  2.1783E+01  7.8636E-01
             8.4268E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1594.25013754679        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1391             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7968E-01  1.4500E+00  5.2924E-01  6.9985E-01  9.2489E-01  1.0627E+00  1.0208E+00  1.0000E-02  1.3651E+00  9.3913E-01
             1.0091E+00
 PARAMETER:  7.9466E-02  4.7156E-01 -5.3631E-01 -2.5689E-01  2.1916E-02  1.6084E-01  1.2062E-01 -5.0806E+00  4.1124E-01  3.7201E-02
             1.0908E-01
 GRADIENT:   4.2309E+02  3.5835E+02  9.0926E+00  9.8865E+01  8.2000E+00  1.1221E+02  8.2783E+00  0.0000E+00  2.1827E+01  8.0661E-01
             8.2367E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1594.25013754679        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     1452
 NPARAMETR:  9.7967E-01  1.4482E+00  5.2750E-01  7.0073E-01  9.2581E-01  1.0627E+00  1.0215E+00  1.0000E-02  1.3657E+00  9.3869E-01
             1.0093E+00
 PARAMETER:  7.9466E-02  4.7156E-01 -5.3631E-01 -2.5689E-01  2.1916E-02  1.6084E-01  1.2062E-01 -5.0806E+00  4.1124E-01  3.7201E-02
             1.0908E-01
 GRADIENT:   6.7602E-03  1.0199E+00  6.7470E-01 -5.7774E-01 -4.1089E-01  2.1177E-03 -6.9741E-02  0.0000E+00 -4.8780E-02  3.7638E-02
            -4.8515E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1452
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3250E-04 -2.3508E-02 -3.5672E-04  1.8654E-02 -2.9495E-02
 SE:             2.9865E-02  2.3840E-02  1.3312E-04  2.4173E-02  2.1976E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9646E-01  3.2409E-01  7.3694E-03  4.4030E-01  1.7956E-01

 ETASHRINKSD(%)  1.0000E-10  2.0133E+01  9.9554E+01  1.9016E+01  2.6377E+01
 ETASHRINKVR(%)  1.0000E-10  3.6213E+01  9.9998E+01  3.4416E+01  4.5796E+01
 EBVSHRINKSD(%)  3.8674E-01  1.9915E+01  9.9631E+01  1.9534E+01  2.5579E+01
 EBVSHRINKVR(%)  7.7199E-01  3.5864E+01  9.9999E+01  3.5253E+01  4.4615E+01
 RELATIVEINF(%)  9.9210E+01  5.6694E+00  1.9600E-04  6.1707E+00  8.8792E+00
 EPSSHRINKSD(%)  4.5472E+01
 EPSSHRINKVR(%)  7.0266E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1594.2501375467916     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -859.09931098305344     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1594.250       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.45E+00  5.29E-01  7.00E-01  9.25E-01  1.06E+00  1.02E+00  1.00E-02  1.37E+00  9.39E-01  1.01E+00
 


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
 #CPUT: Total CPU Time in Seconds,       41.392
Stop Time:
Sun Oct 24 04:20:27 CDT 2021
