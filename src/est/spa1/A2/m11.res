Wed Sep 29 23:03:15 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat11.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 RAW OUTPUT FILE (FILE): m11.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -998.143168011204        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4648E+02  5.0279E+01  6.6206E+01  5.0723E+01  1.9528E+02  5.4277E+01 -7.0547E+01 -1.6638E+02 -9.8884E+01 -8.3998E+01
            -1.6717E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1599.47755334889        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0481E+00  8.9910E-01  1.1532E+00  1.1498E+00  9.5260E-01  8.6896E-01  1.2310E+00  4.4695E-01  1.1089E+00  1.0969E+00
             2.3516E+00
 PARAMETER:  1.4700E-01 -6.3557E-03  2.4256E-01  2.3957E-01  5.1442E-02 -4.0459E-02  3.0784E-01 -7.0530E-01  2.0332E-01  1.9245E-01
             9.5510E-01
 GRADIENT:   4.1612E+02  4.3988E+01  5.5371E+00  9.3357E+01  2.1546E+01 -3.3900E+01 -2.9795E+00  9.6384E-01  2.9642E-02 -1.6857E+00
            -5.8704E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1614.09463140929        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0057E+00  9.1460E-01  3.5224E-01  1.0101E+00  5.1374E-01  9.3870E-01  9.6781E-01  8.8248E-02  9.8222E-01  5.8807E-01
             2.6220E+00
 PARAMETER:  1.0565E-01  1.0727E-02 -9.4344E-01  1.1001E-01 -5.6605E-01  3.6737E-02  6.7280E-02 -2.3276E+00  8.2065E-02 -4.3091E-01
             1.0639E+00
 GRADIENT:   1.9504E+02  1.9790E+01 -4.7145E+01  4.5103E+01  7.5356E+01  1.2472E+01 -1.3249E+01  1.7446E-01 -3.3924E+01  7.2932E+00
             4.8657E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1628.78595637654        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      269
 NPARAMETR:  9.5015E-01  7.3362E-01  3.3060E-01  1.1051E+00  4.3813E-01  8.5563E-01  1.1398E+00  9.2373E-02  1.0246E+00  5.3038E-01
             2.4726E+00
 PARAMETER:  4.8866E-02 -2.0977E-01 -1.0069E+00  1.9995E-01 -7.2524E-01 -5.5922E-02  2.3089E-01 -2.2819E+00  1.2433E-01 -5.3416E-01
             1.0053E+00
 GRADIENT:   1.3378E+01  1.1760E+01 -8.3160E+01  7.6440E+01  9.3013E+01 -1.8729E+01 -4.2587E+00  2.0292E-01 -2.3179E+01  3.2306E+00
             1.3623E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1649.90613229691        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.4194E-01  3.2451E-01  4.3039E-01  1.2672E+00  3.7763E-01  8.8436E-01  1.7508E+00  6.9161E-02  1.0496E+00  6.5040E-01
             2.4086E+00
 PARAMETER:  4.0185E-02 -1.0254E+00 -7.4307E-01  3.3684E-01 -8.7385E-01 -2.2895E-02  6.6009E-01 -2.5713E+00  1.4841E-01 -3.3017E-01
             9.7903E-01
 GRADIENT:   1.3384E+01  1.1305E+01  2.1905E+01 -7.0889E+00 -3.5984E+01 -3.5160E-01 -3.4532E+00  1.7180E-01  7.5959E+00  9.8450E+00
             9.7482E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1654.64363923920        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  9.3650E-01  1.7864E-01  4.2870E-01  1.3283E+00  3.5784E-01  9.0061E-01  2.9328E+00  4.7610E-02  9.9868E-01  5.8313E-01
             2.3939E+00
 PARAMETER:  3.4390E-02 -1.6224E+00 -7.4701E-01  3.8391E-01 -9.2768E-01 -4.6825E-03  1.1760E+00 -2.9447E+00  9.8680E-02 -4.3935E-01
             9.7293E-01
 GRADIENT:   1.4719E+01  5.0447E+00  1.7711E+01 -3.1953E+00 -3.1823E+01  8.3925E+00  1.4496E+00  6.4306E-02  4.4862E+00 -6.8804E-01
            -1.6203E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1656.04757501322        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      795
 NPARAMETR:  9.2797E-01  1.1947E-01  4.4945E-01  1.3669E+00  3.6922E-01  8.7843E-01  4.0118E+00  3.6478E-02  9.6802E-01  5.5986E-01
             2.3924E+00
 PARAMETER:  2.5245E-02 -2.0247E+00 -6.9974E-01  4.1252E-01 -8.9635E-01 -2.9615E-02  1.4892E+00 -3.2110E+00  6.7494E-02 -4.8007E-01
             9.7231E-01
 GRADIENT:  -1.3296E+00 -3.0223E+00  9.3109E+00  1.1013E+01 -1.1314E+01 -7.7766E-02 -7.4212E+00  4.7529E-02  3.0447E+00  2.2853E+00
             2.0405E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1656.47563853818        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  9.2507E-01  7.7467E-02  4.3807E-01  1.3727E+00  3.5968E-01  8.7645E-01  5.1496E+00  2.4125E-02  9.6009E-01  5.4547E-01
             2.3952E+00
 PARAMETER:  2.2110E-02 -2.4579E+00 -7.2537E-01  4.1676E-01 -9.2253E-01 -3.1880E-02  1.7389E+00 -3.6245E+00  5.9275E-02 -5.0611E-01
             9.7348E-01
 GRADIENT:   5.0079E+00  2.5365E+01 -3.5782E+01 -3.6415E+01  4.5613E+01 -2.1764E+00  4.7286E+01  7.5095E-03 -1.8000E+01 -1.0391E+01
            -1.2015E+01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1656.52932730358        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     1096
 NPARAMETR:  9.2573E-01  7.3203E-02  4.3308E-01  1.3690E+00  3.5702E-01  8.7829E-01  5.2468E+00  2.2120E-02  9.6166E-01  5.3134E-01
             2.3949E+00
 PARAMETER:  2.2812E-02 -2.5142E+00 -7.3696E-01  4.1402E-01 -9.2983E-01 -2.9771E-02  1.7578E+00 -3.6745E+00  6.1064E-02 -5.3245E-01
             9.7319E-01
 GRADIENT:  -1.3813E+00  2.0124E+01 -8.3960E+01 -1.5302E+02  6.2219E+01  1.2234E-01  2.6007E+01  1.8261E-02  1.8408E+00 -1.2400E+02
            -6.5533E+01
 NUMSIGDIG:         2.3         2.4         2.4         2.3         2.4         2.5         2.5         0.6         1.4         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1096
 NO. OF SIG. DIGITS IN FINAL EST.:  0.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3617E-03  3.4056E-02 -5.9266E-05 -1.7600E-02  3.7890E-03
 SE:             2.9178E-02  1.3383E-02  5.0853E-04  2.7185E-02  1.7905E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3549E-01  1.0935E-02  9.0722E-01  5.1736E-01  8.3241E-01

 ETASHRINKSD(%)  2.2490E+00  5.5167E+01  9.8296E+01  8.9265E+00  4.0015E+01
 ETASHRINKVR(%)  4.4474E+00  7.9900E+01  9.9971E+01  1.7056E+01  6.4018E+01
 EBVSHRINKSD(%)  2.2295E+00  6.7062E+01  9.8002E+01  6.8626E+00  3.7255E+01
 EBVSHRINKVR(%)  4.4093E+00  8.9151E+01  9.9960E+01  1.3254E+01  6.0631E+01
 RELATIVEINF(%)  9.4798E+01  6.4203E+00  2.1149E-03  4.2084E+01  2.0571E+00
 EPSSHRINKSD(%)  2.6084E+01
 EPSSHRINKVR(%)  4.5364E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1656.5293273035791     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -737.59079409890637     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.38
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1656.529       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.26E-01  7.32E-02  4.33E-01  1.37E+00  3.57E-01  8.78E-01  5.25E+00  2.29E-02  9.62E-01  5.31E-01  2.39E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.68E+03
 
 TH 2
+       -1.12E+02  4.97E+04
 
 TH 3
+        1.04E+01 -1.88E+02  2.09E+04
 
 TH 4
+        9.53E+00  1.85E+02 -3.47E+02  5.72E+03
 
 TH 5
+        6.89E+01 -1.85E+02 -6.91E+03 -9.72E+01  2.74E+04
 
 TH 6
+       -1.05E+01 -1.57E+00  1.03E+00 -1.22E+01  5.58E+00  2.38E+02
 
 TH 7
+       -1.26E-01  2.81E+01 -1.13E+01 -7.62E-01  1.34E+01  1.07E-01  2.15E+01
 
 TH 8
+       -3.98E+00 -3.95E+00 -3.29E+00  1.23E-01  4.40E+00  6.86E-01 -7.12E-02  3.55E+01
 
 TH 9
+       -1.56E+02  2.89E+01 -3.59E+01 -6.68E+01  9.27E+01  3.05E+01  1.30E+00  8.86E+00  5.38E+02
 
 TH10
+        8.18E+01 -1.62E+02 -1.27E+02 -6.18E+00  6.50E+01 -1.68E+01 -2.54E+00  6.11E-02  6.13E+04  2.11E+04
 
 TH11
+       -1.11E+01 -2.81E+01 -3.78E+01  1.25E+03  4.46E+01  1.49E+00 -4.78E-01  4.20E-01 -8.23E+00  2.50E+01  3.91E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     5.02
 #CPUT: Total CPU Time in Seconds,       25.891
Stop Time:
Wed Sep 29 23:03:47 CDT 2021
