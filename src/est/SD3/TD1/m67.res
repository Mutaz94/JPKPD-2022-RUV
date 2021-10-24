Sun Oct 24 00:31:46 CDT 2021
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
$DATA ../../../../data/SD3/TD1/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2112.68745842998        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1769E+02 -8.9681E+01 -2.7621E+01 -6.9422E+01  5.0633E+01  4.0481E+01 -1.5561E+01  8.2479E+00  8.5212E+00  5.4405E+00
            -2.2090E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2117.47549373066        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      143
 NPARAMETR:  9.7677E-01  1.0058E+00  1.0364E+00  1.0044E+00  9.9679E-01  1.0102E+00  1.3045E+00  8.7649E-01  9.8737E-01  9.9977E-01
             1.0020E+00
 PARAMETER:  7.6497E-02  1.0579E-01  1.3571E-01  1.0437E-01  9.6788E-02  1.1019E-01  3.6583E-01 -3.1831E-02  8.7286E-02  9.9765E-02
             1.0205E-01
 GRADIENT:   3.6521E+02 -3.6408E+01  9.2458E+00 -4.1704E+01  1.7220E+00  5.2667E+01  2.3436E+01  2.9164E+00  2.0606E+01  6.8938E+00
             5.4901E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2119.05900156026        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      321
 NPARAMETR:  9.9628E-01  1.0058E+00  9.7986E-01  1.0044E+00  9.9678E-01  1.0344E+00  1.3190E+00  6.9553E-01  9.0324E-01  9.9975E-01
             1.0029E+00
 PARAMETER:  9.6275E-02  1.0581E-01  7.9656E-02  1.0438E-01  9.6778E-02  1.3377E-01  3.7690E-01 -2.6308E-01 -1.7616E-03  9.9755E-02
             1.0292E-01
 GRADIENT:   2.2872E-01 -8.5039E+01 -9.5563E-01 -1.2510E+02  2.2972E+01  1.3257E+00 -1.9786E-01  1.5871E-01 -2.1309E-02  2.8184E+00
            -2.9029E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2123.11286890340        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      457
 NPARAMETR:  9.7783E-01  1.1206E+00  9.7924E-01  1.0773E+00  9.7928E-01  1.0153E+00  1.2962E+00  6.9877E-01  8.9856E-01  9.8039E-01
             1.0038E+00
 PARAMETER:  7.7578E-02  2.1386E-01  7.9017E-02  1.7444E-01  7.9063E-02  1.1520E-01  3.5947E-01 -2.5843E-01 -6.9667E-03  8.0195E-02
             1.0375E-01
 GRADIENT:   3.6175E+02  1.2983E+02  9.5759E+00  1.9751E+02 -1.2262E+00  5.6463E+01  2.5796E+01 -1.4688E+00  8.9326E+00  9.3888E-01
            -2.3284E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2125.49815861443        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      594
 NPARAMETR:  9.9811E-01  1.0761E+00  9.7532E-01  1.0526E+00  9.7591E-01  1.0328E+00  1.2377E+00  7.3133E-01  8.7947E-01  9.6972E-01
             1.0057E+00
 PARAMETER:  9.8106E-02  1.7332E-01  7.5008E-02  1.5123E-01  7.5611E-02  1.3224E-01  3.1325E-01 -2.1289E-01 -2.8440E-02  6.9251E-02
             1.0572E-01
 GRADIENT:   2.5742E+00  2.7206E-01  6.8048E+00 -7.0851E+00 -2.7462E+00  4.3067E-01 -7.9859E-01 -7.6983E-01 -2.2010E+00  2.2973E-01
            -9.7832E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2125.77529416709        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      770
 NPARAMETR:  9.9697E-01  1.0541E+00  8.6679E-01  1.0615E+00  9.0507E-01  1.0312E+00  1.2530E+00  5.8423E-01  8.8388E-01  8.9766E-01
             1.0045E+00
 PARAMETER:  9.6968E-02  1.5272E-01 -4.2954E-02  1.5965E-01  2.5723E-04  1.3074E-01  3.2557E-01 -4.3746E-01 -2.3429E-02 -7.9612E-03
             1.0444E-01
 GRADIENT:  -1.1077E+00  2.1506E+00  2.8133E+00  2.3326E+00 -3.5143E+00 -3.9060E-01 -5.5929E-01  2.7855E-01  4.7416E-01 -2.1008E+00
            -3.7149E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2125.97544613765        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      950
 NPARAMETR:  9.9855E-01  1.0753E+00  7.8476E-01  1.0375E+00  8.7498E-01  1.0332E+00  1.2361E+00  2.9237E-01  8.8667E-01  9.0675E-01
             1.0029E+00
 PARAMETER:  9.8547E-02  1.7257E-01 -1.4237E-01  1.3684E-01 -3.3556E-02  1.3267E-01  3.1199E-01 -1.1297E+00 -2.0285E-02  2.1144E-03
             1.0291E-01
 GRADIENT:   8.0336E-01 -2.8626E+00  1.5357E+00 -1.3589E+00 -1.0974E+00  6.5362E-02 -5.6198E-01 -2.0603E-02 -8.1062E-01  6.8791E-01
            -1.2074E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2126.22521886132        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1125
 NPARAMETR:  9.9938E-01  1.2194E+00  7.1439E-01  9.4869E-01  9.0389E-01  1.0350E+00  1.1326E+00  7.2833E-02  9.4062E-01  8.9939E-01
             1.0029E+00
 PARAMETER:  9.9380E-02  2.9839E-01 -2.3633E-01  4.7327E-02 -1.0448E-03  1.3436E-01  2.2449E-01 -2.5196E+00  3.8785E-02 -6.0405E-03
             1.0286E-01
 GRADIENT:   1.2759E-01  5.1451E-01  5.5265E-01 -1.5947E-01  3.6851E-02  1.0177E-01  5.7800E-02  1.1240E-02 -4.5080E-02 -5.0004E-01
            -9.2686E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2126.23309526910        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1283
 NPARAMETR:  1.0003E+00  1.2149E+00  7.1347E-01  9.5144E-01  9.0111E-01  1.0354E+00  1.1349E+00  3.0177E-02  9.3908E-01  9.0080E-01
             1.0028E+00
 PARAMETER:  1.0026E-01  2.9466E-01 -2.3762E-01  5.0221E-02 -4.1282E-03  1.3479E-01  2.2654E-01 -3.4007E+00  3.7146E-02 -4.4694E-03
             1.0276E-01
 GRADIENT:   1.9605E+00  3.3121E-01 -9.3008E-03  3.5605E-01  1.4166E-01  2.7839E-01  5.3965E-02  2.2495E-03  3.6168E-02 -9.1046E-03
             9.9198E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2126.23419081835        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1462
 NPARAMETR:  9.9997E-01  1.2146E+00  7.1334E-01  9.5124E-01  9.0099E-01  1.0354E+00  1.1348E+00  1.0082E-02  9.3900E-01  9.0083E-01
             1.0027E+00
 PARAMETER:  9.9975E-02  2.9441E-01 -2.3779E-01  5.0007E-02 -4.2572E-03  1.3477E-01  2.2643E-01 -4.4970E+00  3.7058E-02 -4.4412E-03
             1.0273E-01
 GRADIENT:   1.3715E+00 -2.5986E-02  7.3247E-02 -1.7819E-01  7.0703E-02  2.7022E-01  1.8006E-02  5.7591E-04  7.3852E-03 -3.4878E-03
            -1.5285E-03

0ITERATION NO.:   46    OBJECTIVE VALUE:  -2126.23419081835        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     1486
 NPARAMETR:  9.9997E-01  1.2146E+00  7.1334E-01  9.5124E-01  9.0099E-01  1.0354E+00  1.1348E+00  1.0082E-02  9.3900E-01  9.0083E-01
             1.0027E+00
 PARAMETER:  9.9975E-02  2.9441E-01 -2.3779E-01  5.0007E-02 -4.2572E-03  1.3477E-01  2.2643E-01 -4.4970E+00  3.7058E-02 -4.4412E-03
             1.0273E-01
 GRADIENT:  -1.4385E-01  1.6370E-01  7.9718E-02 -8.0981E-02  8.2299E-02 -3.8133E-02 -1.9031E-02  9.1851E-05 -2.1681E-03 -4.4321E-03
            -1.1065E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1486
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1509E-04 -7.4955E-03 -4.1533E-04  4.5135E-03 -1.5976E-02
 SE:             2.9876E-02  2.3044E-02  1.7484E-04  2.4196E-02  2.2701E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8891E-01  7.4497E-01  1.7525E-02  8.5203E-01  4.8160E-01

 ETASHRINKSD(%)  1.0000E-10  2.2801E+01  9.9414E+01  1.8939E+01  2.3948E+01
 ETASHRINKVR(%)  1.0000E-10  4.0403E+01  9.9997E+01  3.4291E+01  4.2161E+01
 EBVSHRINKSD(%)  3.2055E-01  2.2441E+01  9.9486E+01  1.9639E+01  2.3178E+01
 EBVSHRINKVR(%)  6.4006E-01  3.9846E+01  9.9997E+01  3.5421E+01  4.0984E+01
 RELATIVEINF(%)  9.9153E+01  4.7395E+00  3.6849E-04  5.3534E+00  7.6253E+00
 EPSSHRINKSD(%)  3.3546E+01
 EPSSHRINKVR(%)  5.5839E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2126.2341908183548     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1207.2956576136821     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2126.234       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.21E+00  7.13E-01  9.51E-01  9.01E-01  1.04E+00  1.13E+00  1.01E-02  9.39E-01  9.01E-01  1.00E+00
 


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
 #CPUT: Total CPU Time in Seconds,       48.701
Stop Time:
Sun Oct 24 00:31:56 CDT 2021
