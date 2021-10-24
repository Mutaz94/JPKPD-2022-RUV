Sat Oct 23 18:11:55 CDT 2021
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
$DATA ../../../../data/SD2/S1/dat7.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2402.01432501758        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9334E+02  7.6850E+01  5.2768E+01  1.0763E+02  1.2035E+02  2.2661E+01 -1.7436E+01 -2.1241E+02  1.1003E+01 -3.7934E+01
            -7.9553E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2753.76262749103        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      135
 NPARAMETR:  7.8077E-01  9.3999E-01  9.7415E-01  9.2596E-01  9.1587E-01  1.0377E+00  1.0049E+00  1.1534E+00  9.7364E-01  1.0289E+00
             1.6035E+00
 PARAMETER: -1.4748E-01  3.8116E-02  7.3809E-02  2.3073E-02  1.2122E-02  1.3703E-01  1.0494E-01  2.4269E-01  7.3288E-02  1.2845E-01
             5.7218E-01
 GRADIENT:  -4.8173E+02 -8.2668E+01  3.9511E+00 -1.5732E+02 -3.6844E+01 -1.3461E+02  4.8012E+00  1.3009E+01  9.2441E+00  1.2041E+01
             4.6904E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2806.93430822314        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      310
 NPARAMETR:  8.3577E-01  8.5905E-01  6.1163E-01  9.7750E-01  6.9307E-01  9.9506E-01  1.2535E+00  1.3351E-01  9.6242E-01  1.3620E+00
             1.3912E+00
 PARAMETER: -7.9402E-02 -5.1932E-02 -3.9162E-01  7.7244E-02 -2.6662E-01  9.5050E-02  3.2597E-01 -1.9136E+00  6.1691E-02  4.0896E-01
             4.3018E-01
 GRADIENT:  -3.7254E+02 -1.1764E+01 -1.0582E+02 -9.5024E+01 -3.5489E+01 -1.0401E+02  3.4431E+01  3.9296E-01  1.9557E+01  7.9575E+01
             3.3073E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2929.68200575172        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      485
 NPARAMETR:  9.5888E-01  9.2470E-01  7.9800E-01  1.0205E+00  8.3904E-01  9.8969E-01  1.1109E+00  5.9710E-01  7.6437E-01  1.0053E+00
             9.9092E-01
 PARAMETER:  5.8007E-02  2.1710E-02 -1.2564E-01  1.2031E-01 -7.5501E-02  8.9638E-02  2.0515E-01 -4.1567E-01 -1.6870E-01  1.0527E-01
             9.0875E-02
 GRADIENT:  -6.2317E+01  2.5255E+01 -4.9188E+01  2.5109E+01  1.5306E-01 -4.7269E+01  2.8849E+00 -1.2282E+00 -4.2135E+01  1.9411E+00
            -5.4157E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2941.16058175004        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      662
 NPARAMETR:  9.7945E-01  9.5505E-01  9.7578E-01  1.0021E+00  9.3611E-01  1.0739E+00  1.0247E+00  9.4458E-01  9.0219E-01  1.0264E+00
             1.0214E+00
 PARAMETER:  7.9231E-02  5.4004E-02  7.5480E-02  1.0211E-01  3.3974E-02  1.7129E-01  1.2438E-01  4.2987E-02 -2.9255E-03  1.2610E-01
             1.2122E-01
 GRADIENT:  -1.1219E+01 -7.6052E-01 -3.6272E+00 -1.6171E+00  2.2375E+00 -8.1883E+00  2.0700E+00  1.1208E+00  1.7944E+00 -1.4311E+00
             4.9207E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2941.27886261974        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  9.7960E-01  9.5498E-01  9.7618E-01  1.0017E+00  9.3677E-01  1.0972E+00  1.0232E+00  9.4420E-01  8.9552E-01  1.0280E+00
             1.0181E+00
 PARAMETER:  7.9392E-02  5.3930E-02  7.5895E-02  1.0174E-01  3.4678E-02  1.9278E-01  1.2291E-01  4.2580E-02 -1.0352E-02  1.2765E-01
             1.1797E-01
 GRADIENT:  -1.0301E+01 -1.2521E+00 -3.7285E+00 -2.7604E+00  2.7230E+00  6.9335E-01  1.6606E+00  9.9455E-01 -1.2221E-02 -1.4144E+00
             3.7264E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2941.31672026244        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  9.8036E-01  9.5468E-01  9.7716E-01  1.0017E+00  9.3669E-01  1.0974E+00  1.0041E+00  9.2784E-01  8.9576E-01  1.0298E+00
             1.0182E+00
 PARAMETER:  8.0163E-02  5.3618E-02  7.6895E-02  1.0170E-01  3.4597E-02  1.9292E-01  1.0409E-01  2.5101E-02 -1.0084E-02  1.2940E-01
             1.1800E-01
 GRADIENT:  -8.8999E+00 -2.0689E+00 -1.4234E+00 -3.2263E+00  2.3574E+00  7.7021E-01 -1.1188E-02 -3.3349E-02 -7.0458E-01 -2.0796E+00
            -4.7713E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2941.33291377734        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1145
 NPARAMETR:  9.8086E-01  9.5505E-01  9.7809E-01  1.0019E+00  9.3600E-01  1.0944E+00  1.0020E+00  9.2833E-01  8.9825E-01  1.0353E+00
             1.0185E+00
 PARAMETER:  8.0672E-02  5.4012E-02  7.7849E-02  1.0192E-01  3.3859E-02  1.9018E-01  1.0201E-01  2.5632E-02 -7.3092E-03  1.3464E-01
             1.1829E-01
 GRADIENT:   4.5448E+02  5.0475E+01  4.7207E+00  1.0182E+02  2.5170E+01  1.3136E+02  4.7987E+00  4.4028E-01  8.0580E+00  2.2246E+00
             1.6720E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2941.34217490992        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1326
 NPARAMETR:  9.8136E-01  9.5486E-01  9.7802E-01  1.0024E+00  9.3683E-01  1.0971E+00  9.9931E-01  9.2722E-01  8.9846E-01  1.0384E+00
             1.0185E+00
 PARAMETER:  8.1188E-02  5.3804E-02  7.7772E-02  1.0236E-01  3.4750E-02  1.9272E-01  9.9307E-02  2.4432E-02 -7.0762E-03  1.3769E-01
             1.1830E-01
 GRADIENT:  -7.0300E+00 -1.8725E+00 -1.1202E+00 -1.7143E+00  1.5691E+00  7.0602E-01 -2.6114E-03 -1.9816E-04 -2.6763E-02 -7.0737E-01
             1.8247E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -2941.34217490992        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  9.8235E-01  9.5518E-01  9.7842E-01  1.0025E+00  9.3666E-01  1.0972E+00  9.9932E-01  9.2722E-01  8.9847E-01  1.0392E+00
             1.0185E+00
 PARAMETER:  8.1188E-02  5.3804E-02  7.7772E-02  1.0236E-01  3.4750E-02  1.9272E-01  9.9307E-02  2.4432E-02 -7.0762E-03  1.3769E-01
             1.1830E-01
 GRADIENT:  -8.4474E+00 -1.6083E+00 -9.8862E-01 -9.9373E-01  1.2409E+00 -1.9996E-01 -7.9937E-03 -5.8510E-04 -4.1501E-02 -6.4863E-01
             1.8928E-02
 NUMSIGDIG:         1.0         1.5         1.4         2.0         1.8         2.6         3.0         3.9         2.8         1.3
                    4.0

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1355
 NO. OF SIG. DIGITS IN FINAL EST.:  1.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.7263E-03 -1.5933E-02 -2.4072E-02  8.3110E-03 -2.2431E-02
 SE:             2.9923E-02  2.0044E-02  1.4304E-02  2.6304E-02  2.4572E-02
 N:                     100         100         100         100         100

 P VAL.:         8.7450E-01  4.2667E-01  9.2384E-02  7.5203E-01  3.6132E-01

 ETASHRINKSD(%)  1.0000E-10  3.2849E+01  5.2081E+01  1.1880E+01  1.7680E+01
 ETASHRINKVR(%)  1.0000E-10  5.4907E+01  7.7038E+01  2.2348E+01  3.2234E+01
 EBVSHRINKSD(%)  2.3613E-01  3.2673E+01  5.4616E+01  1.2903E+01  1.6586E+01
 EBVSHRINKVR(%)  4.7170E-01  5.4671E+01  7.9403E+01  2.4141E+01  3.0420E+01
 RELATIVEINF(%)  9.9523E+01  1.5149E+01  1.1004E+01  3.1596E+01  2.3515E+01
 EPSSHRINKSD(%)  2.5371E+01
 EPSSHRINKVR(%)  4.4304E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2941.3421749099152     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1654.8282284233735     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2941.342       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  9.55E-01  9.78E-01  1.00E+00  9.37E-01  1.10E+00  9.99E-01  9.27E-01  8.98E-01  1.04E+00  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,       86.325
Stop Time:
Sat Oct 23 18:12:10 CDT 2021
