Sun Oct 24 01:34:23 CDT 2021
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
$DATA ../../../../data/SD4/B/dat15.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1688.03931462483        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3732E+02 -5.7458E+01 -1.4979E+01 -4.8398E+01  5.8671E+01  6.0038E+01 -8.0599E+00 -2.2406E+00 -1.9285E+01 -3.2487E+00
             2.3688E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1698.33582560039        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0420E+00  1.0442E+00  9.9048E-01  1.0556E+00  9.7203E-01  9.2946E-01  1.0349E+00  1.0256E+00  1.0736E+00  1.0006E+00
             9.4078E-01
 PARAMETER:  1.4117E-01  1.4328E-01  9.0432E-02  1.5415E-01  7.1634E-02  2.6844E-02  1.3429E-01  1.2532E-01  1.7098E-01  1.0059E-01
             3.8956E-02
 GRADIENT:   4.9623E+00 -6.0935E+00 -5.3007E+00  4.1020E+00  9.8908E+00 -1.3010E+00 -1.7731E+00 -1.8068E+00  1.8743E-01  5.8617E-01
            -7.4941E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1698.76115537359        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0434E+00  1.2345E+00  9.4230E-01  9.4715E-01  1.0161E+00  9.3356E-01  9.6532E-01  1.1925E+00  1.1476E+00  9.7172E-01
             9.4197E-01
 PARAMETER:  1.4251E-01  3.1070E-01  4.0568E-02  4.5704E-02  1.1593E-01  3.1247E-02  6.4702E-02  2.7604E-01  2.3766E-01  7.1308E-02
             4.0214E-02
 GRADIENT:   6.0081E+00  1.5301E+01  2.5495E+00  1.4831E+01 -6.7808E+00  9.7958E-02  7.3763E-02  2.0990E-01 -2.7781E+00 -2.6499E+00
            -9.9531E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1699.29684161602        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0413E+00  1.4851E+00  8.6782E-01  7.7864E-01  1.1072E+00  9.3655E-01  8.2570E-01  1.3417E+00  1.3507E+00  1.0406E+00
             9.4230E-01
 PARAMETER:  1.4049E-01  4.9545E-01 -4.1772E-02 -1.5020E-01  2.0180E-01  3.4447E-02 -9.1522E-02  3.9394E-01  4.0061E-01  1.3982E-01
             4.0565E-02
 GRADIENT:  -5.9169E-01  1.5897E+01  7.2704E+00  5.8921E+00 -1.2769E+01  1.1332E+00 -1.0168E+00 -1.1841E+00 -7.1742E-01 -1.1987E+00
            -1.3910E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1699.58199992384        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0432E+00  1.6864E+00  7.2195E-01  6.4386E-01  1.1583E+00  9.3494E-01  7.8614E-01  1.3742E+00  1.5259E+00  1.0613E+00
             9.4567E-01
 PARAMETER:  1.4225E-01  6.2257E-01 -2.2580E-01 -3.4027E-01  2.4696E-01  3.2729E-02 -1.4062E-01  4.1789E-01  5.2256E-01  1.5954E-01
             4.4136E-02
 GRADIENT:   2.4789E+00  1.6290E+01  5.4551E+00  5.0475E+00 -9.4817E+00  2.3008E-01 -1.4430E-01 -6.1600E-01 -1.0986E+00 -8.2601E-01
            -5.6278E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1699.64092192489        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0432E+00  1.7843E+00  6.3183E-01  5.7901E-01  1.1772E+00  9.3438E-01  7.6762E-01  1.3236E+00  1.6295E+00  1.0697E+00
             9.4723E-01
 PARAMETER:  1.4226E-01  6.7903E-01 -3.5913E-01 -4.4643E-01  2.6317E-01  3.2127E-02 -1.6446E-01  3.8035E-01  5.8826E-01  1.6738E-01
             4.5782E-02
 GRADIENT:   1.8267E+00  1.8041E+01  3.9426E+00  6.4719E+00 -7.9370E+00 -1.3579E-01 -6.5910E-01 -4.6624E-01 -1.1311E+00 -3.4412E-01
            -2.9970E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1699.65662959569        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  1.0431E+00  1.8374E+00  5.8037E-01  5.4228E-01  1.1883E+00  9.3433E-01  7.5900E-01  1.2774E+00  1.6951E+00  1.0738E+00
             9.4814E-01
 PARAMETER:  1.4222E-01  7.0833E-01 -4.4408E-01 -5.1197E-01  2.7256E-01  3.2069E-02 -1.7575E-01  3.4484E-01  6.2777E-01  1.7125E-01
             4.6745E-02
 GRADIENT:   1.4684E+00  1.6089E+01  2.7239E+00  6.4841E+00 -6.0752E+00 -2.1984E-01 -8.1504E-01 -3.4880E-01 -9.3485E-01 -1.4746E-01
            -1.6620E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1699.69964489758        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  1.0428E+00  1.8356E+00  5.7217E-01  5.4264E-01  1.1894E+00  9.3459E-01  7.6275E-01  1.2613E+00  1.6974E+00  1.0736E+00
             9.4808E-01
 PARAMETER:  1.4195E-01  7.0739E-01 -4.5832E-01 -5.1132E-01  2.7341E-01  3.2348E-02 -1.7083E-01  3.3213E-01  6.2909E-01  1.7103E-01
             4.6687E-02
 GRADIENT:   6.7953E-01  1.0482E+01  3.0383E-01  7.9087E+00 -7.3400E-01 -1.4617E-01 -5.1647E-02  1.7935E-02 -8.3579E-02  6.3000E-02
             2.1212E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1699.75066141802        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1403             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0437E+00  1.8389E+00  5.6686E-01  5.3621E-01  1.1910E+00  9.3495E-01  7.6182E-01  1.2667E+00  1.6996E+00  1.0730E+00
             9.4794E-01
 PARAMETER:  1.4275E-01  7.0915E-01 -4.6765E-01 -5.2322E-01  2.7481E-01  3.2736E-02 -1.7205E-01  3.3639E-01  6.3038E-01  1.7045E-01
             4.6533E-02
 GRADIENT:   6.9266E+02  8.7538E+02  3.6565E+00  1.1119E+02  1.8452E+01  3.7921E+01  1.2584E+01  4.8174E-01  2.9069E+01  1.5308E+00
             7.6689E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1699.78421507580        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1584
 NPARAMETR:  1.0435E+00  1.8422E+00  5.6627E-01  5.2790E-01  1.1942E+00  9.3507E-01  7.6093E-01  1.2444E+00  1.7286E+00  1.0738E+00
             9.4818E-01
 PARAMETER:  1.4257E-01  7.1093E-01 -4.6869E-01 -5.3884E-01  2.7747E-01  3.2865E-02 -1.7321E-01  3.1863E-01  6.4733E-01  1.7117E-01
             4.6791E-02
 GRADIENT:   2.5390E+00 -3.5919E+00  1.5158E+00  7.8746E-01 -1.4797E+00  1.0286E-01  2.9278E-01 -2.1191E-01  6.0089E-01 -2.8862E-01
             7.0010E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1699.79444302564        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:     1686
 NPARAMETR:  1.0438E+00  1.8442E+00  5.6308E-01  5.2716E-01  1.1943E+00  9.3506E-01  7.5929E-01  1.2501E+00  1.7286E+00  1.0740E+00
             9.4792E-01
 PARAMETER:  1.4269E-01  7.1012E-01 -4.7312E-01 -5.3884E-01  2.7778E-01  3.2858E-02 -1.7379E-01  3.2650E-01  6.4733E-01  1.7187E-01
             4.6748E-02
 GRADIENT:  -9.3351E-01 -8.1249E+03  3.3556E-01  1.1620E+00  5.3552E-01  2.3755E-03  2.3480E-01  3.9712E-02 -7.0971E+01  3.3569E+04
             1.8073E-01
 NUMSIGDIG:         2.6         2.3         2.3         2.3         2.8         4.2         1.8         1.7         4.7         2.3
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1686
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1853E-04 -4.2246E-02 -2.6665E-02  3.2540E-02 -4.9484E-02
 SE:             2.9864E-02  2.2581E-02  9.9051E-03  2.3129E-02  2.1749E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7813E-01  6.1365E-02  7.1014E-03  1.5947E-01  2.2893E-02

 ETASHRINKSD(%)  1.0000E-10  2.4351E+01  6.6817E+01  2.2514E+01  2.7138E+01
 ETASHRINKVR(%)  1.0000E-10  4.2772E+01  8.8989E+01  3.9959E+01  4.6911E+01
 EBVSHRINKSD(%)  4.5262E-01  2.3762E+01  6.9484E+01  2.3047E+01  2.4470E+01
 EBVSHRINKVR(%)  9.0320E-01  4.1878E+01  9.0688E+01  4.0782E+01  4.2952E+01
 RELATIVEINF(%)  9.9070E+01  4.2703E+00  1.0741E+00  4.6607E+00  1.7382E+01
 EPSSHRINKSD(%)  4.6229E+01
 EPSSHRINKVR(%)  7.1087E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1699.7944430256416     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -964.64361646190343     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1699.794       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.84E+00  5.64E-01  5.28E-01  1.19E+00  9.35E-01  7.60E-01  1.25E+00  1.73E+00  1.07E+00  9.48E-01
 


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
 #CPUT: Total CPU Time in Seconds,       51.011
Stop Time:
Sun Oct 24 01:34:34 CDT 2021
