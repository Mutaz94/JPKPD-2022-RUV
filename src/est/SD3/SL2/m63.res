Sun Oct 24 00:00:57 CDT 2021
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
$DATA ../../../../data/SD3/SL2/dat63.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2129.66065126851        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2677E+02 -8.2781E+01 -3.8448E+01 -6.4357E+01  5.3012E+01  1.9794E+01 -1.2771E+01  1.2516E+01 -1.6734E+00  8.0395E+00
            -1.5544E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2144.45581819690        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0468E+00  1.0906E+00  1.1111E+00  1.0550E+00  1.0426E+00  1.0702E+00  1.0776E+00  9.3183E-01  9.8826E-01  9.6059E-01
             1.0123E+00
 PARAMETER:  1.4576E-01  1.8675E-01  2.0531E-01  1.5351E-01  1.4170E-01  1.6786E-01  1.7470E-01  2.9393E-02  8.8195E-02  5.9797E-02
             1.1227E-01
 GRADIENT:   3.9780E-01 -3.6913E-01 -2.2316E+00  2.4274E+00  6.2281E+00  4.2395E+00 -7.0173E+00  2.9934E+00  1.1593E+00 -7.1191E+00
            -1.2891E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2145.57843085675        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0436E+00  9.2314E-01  1.1821E+00  1.1729E+00  1.0101E+00  1.0532E+00  1.3923E+00  7.0309E-01  8.4609E-01  1.0649E+00
             1.0256E+00
 PARAMETER:  1.4268E-01  2.0027E-02  2.6729E-01  2.5952E-01  1.1005E-01  1.5179E-01  4.3095E-01 -2.5227E-01 -6.7129E-02  1.6286E-01
             1.2532E-01
 GRADIENT:  -2.5118E+00  1.3895E+01  2.7457E+00  2.0490E+01 -2.0058E+00 -1.1956E+00 -1.0109E+00 -1.2840E+00 -6.7041E+00  5.1611E+00
            -2.7008E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2146.28827901422        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0449E+00  8.5137E-01  1.0992E+00  1.1961E+00  9.3582E-01  1.0596E+00  1.4540E+00  6.6620E-01  8.6289E-01  9.4424E-01
             1.0284E+00
 PARAMETER:  1.4396E-01 -6.0906E-02  1.9459E-01  2.7910E-01  3.3667E-02  1.5792E-01  4.7434E-01 -3.0617E-01 -4.7465E-02  4.2621E-02
             1.2804E-01
 GRADIENT:   3.7476E-01  2.8492E+00  3.8603E+00 -3.6474E-01 -4.7177E+00  1.1590E+00  1.0629E-01 -3.5919E-01 -2.3968E-01 -1.0802E+00
             1.5411E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2146.31085278539        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      726             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0446E+00  7.5575E-01  1.1910E+00  1.2570E+00  9.4699E-01  1.0560E+00  1.5321E+00  7.6812E-01  8.5010E-01  9.7218E-01
             1.0281E+00
 PARAMETER:  1.4359E-01 -1.8005E-01  2.7478E-01  3.2876E-01  4.5535E-02  1.5449E-01  5.2664E-01 -1.6381E-01 -6.2401E-02  7.1782E-02
             1.2767E-01
 GRADIENT:   5.8123E+02  3.2358E+01  2.9067E+00  3.4899E+02  1.0548E+01  6.0556E+01  1.9142E+01  5.6053E-01  7.6600E+00  5.1905E-01
             1.6753E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2146.31923052714        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      905             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0453E+00  7.5516E-01  1.1947E+00  1.2587E+00  9.4569E-01  1.0555E+00  1.5350E+00  7.6030E-01  8.5018E-01  9.7472E-01
             1.0278E+00
 PARAMETER:  1.4430E-01 -1.8083E-01  2.7791E-01  3.3011E-01  4.4156E-02  1.5398E-01  5.2851E-01 -1.7404E-01 -6.2305E-02  7.4397E-02
             1.2744E-01
 GRADIENT:   5.8593E+02  3.4133E+01  5.9340E+00  3.5224E+02  6.3876E+00  6.0088E+01  1.9349E+01  1.5024E-01  7.7060E+00  6.4313E-01
             1.2087E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2146.31958068868        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  1.0448E+00  7.5523E-01  1.1944E+00  1.2589E+00  9.4582E-01  1.0554E+00  1.5353E+00  7.6032E-01  8.4972E-01  9.7494E-01
             1.0279E+00
 PARAMETER:  1.4387E-01 -1.8074E-01  2.7768E-01  3.3026E-01  4.4295E-02  1.5397E-01  5.2874E-01 -1.7401E-01 -6.2849E-02  7.4617E-02
             1.2747E-01
 GRADIENT:   3.0660E+00  8.2885E-02  9.7996E-02 -2.2166E+00 -1.3230E-01  2.1977E-01  9.2691E-02 -2.2648E-02  1.7519E-02 -3.0127E-02
            -2.8013E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -2146.31958068868        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  1.0448E+00  7.5523E-01  1.1944E+00  1.2589E+00  9.4582E-01  1.0554E+00  1.5353E+00  7.6032E-01  8.4972E-01  9.7494E-01
             1.0279E+00
 PARAMETER:  1.4387E-01 -1.8074E-01  2.7768E-01  3.3026E-01  4.4295E-02  1.5397E-01  5.2874E-01 -1.7401E-01 -6.2849E-02  7.4617E-02
             1.2747E-01
 GRADIENT:  -1.2521E-04  7.8699E-02  1.1109E-01  3.0617E-02 -1.2291E-01  4.9288E-05 -1.6024E-02 -1.1616E-02  6.4298E-03 -3.1054E-02
            -2.7299E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1112
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9676E-04  5.0210E-03 -2.7645E-02 -7.4128E-03 -2.3633E-02
 SE:             2.9879E-02  1.8658E-02  1.2645E-02  2.4880E-02  2.2183E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8407E-01  7.8785E-01  2.8807E-02  7.6575E-01  2.8672E-01

 ETASHRINKSD(%)  1.0000E-10  3.7494E+01  5.7636E+01  1.6649E+01  2.5683E+01
 ETASHRINKVR(%)  1.0000E-10  6.0930E+01  8.2053E+01  3.0526E+01  4.4770E+01
 EBVSHRINKSD(%)  3.3101E-01  3.9008E+01  6.0876E+01  1.6275E+01  2.3032E+01
 EBVSHRINKVR(%)  6.6093E-01  6.2800E+01  8.4693E+01  2.9901E+01  4.0760E+01
 RELATIVEINF(%)  9.8341E+01  2.1496E+00  2.3166E+00  4.5587E+00  1.0124E+01
 EPSSHRINKSD(%)  3.3501E+01
 EPSSHRINKVR(%)  5.5779E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2146.3195806886811     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1227.3810474840084     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2146.320       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  7.55E-01  1.19E+00  1.26E+00  9.46E-01  1.06E+00  1.54E+00  7.60E-01  8.50E-01  9.75E-01  1.03E+00
 


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
 #CPUT: Total CPU Time in Seconds,       37.722
Stop Time:
Sun Oct 24 00:01:06 CDT 2021
