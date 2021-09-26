Sat Sep 25 11:31:25 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat10.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m10.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1598.30904192240        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.2918E+01 -1.9602E+01 -2.9571E+01  4.2361E+00  4.3349E+01 -1.6452E+01  1.0638E+01  7.2221E+00  1.1661E+01  4.8852E+00
            -1.1583E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1613.32672877092        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.8370E-01  1.0075E+00  1.0625E+00  1.0087E+00  9.8167E-01  1.0370E+00  9.2745E-01  9.6203E-01  9.4525E-01  9.4358E-01
             1.2354E+00
 PARAMETER:  8.3566E-02  1.0749E-01  1.6066E-01  1.0862E-01  8.1503E-02  1.3629E-01  2.4680E-02  6.1292E-02  4.3698E-02  4.1925E-02
             3.1142E-01
 GRADIENT:   2.8662E+01  1.3785E+01  5.9347E+00  5.2889E+00 -1.2797E+01  9.7354E-01  5.5355E+00  1.6745E+00 -2.3893E-01  2.3882E+00
            -4.6711E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1614.32040458720        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.7668E-01  9.1121E-01  1.0958E+00  1.0626E+00  9.6236E-01  1.0485E+00  6.8617E-01  8.7194E-01  9.9278E-01  9.8734E-01
             1.2120E+00
 PARAMETER:  7.6399E-02  7.0215E-03  1.9152E-01  1.6072E-01  6.1636E-02  1.4739E-01 -2.7662E-01 -3.7033E-02  9.2751E-02  8.7257E-02
             2.9227E-01
 GRADIENT:   1.7258E+01  4.3958E+00  2.9417E+00  7.2056E+00 -9.0844E+00  5.9837E+00  2.4845E+00 -1.0579E+00  4.8714E+00  3.5660E+00
            -1.1033E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1615.40728771536        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.6722E-01  1.0589E+00  1.2892E+00  9.6117E-01  1.1052E+00  1.0288E+00  2.5119E-01  1.2318E+00  1.1474E+00  1.0670E+00
             1.2589E+00
 PARAMETER:  6.6672E-02  1.5724E-01  3.5401E-01  6.0394E-02  1.9999E-01  1.2842E-01 -1.2815E+00  3.0850E-01  2.3746E-01  1.6486E-01
             3.3027E-01
 GRADIENT:  -5.2621E+00 -2.0307E+00  3.5159E-01 -4.2153E+00  3.7544E+00 -1.9948E+00  7.6150E-01 -6.4381E-01  2.8910E+00 -2.6109E-01
             4.8664E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1615.82214519106        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  9.6899E-01  9.5818E-01  1.3139E+00  1.0290E+00  1.0695E+00  1.0334E+00  1.0594E-01  1.1983E+00  1.0794E+00  1.0531E+00
             1.2422E+00
 PARAMETER:  6.8497E-02  5.7278E-02  3.7299E-01  1.2863E-01  1.6723E-01  1.3281E-01 -2.1449E+00  2.8094E-01  1.7637E-01  1.5176E-01
             3.1688E-01
 GRADIENT:   8.0613E-01  1.2602E+00  6.9791E-01  1.9233E+00 -8.4901E-01  2.0767E-01  1.3139E-01 -8.4687E-02  1.3016E-01 -4.7105E-02
            -7.2175E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1615.83426149035        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  9.6865E-01  9.6289E-01  1.3063E+00  1.0245E+00  1.0700E+00  1.0328E+00  4.8284E-02  1.1960E+00  1.0852E+00  1.0538E+00
             1.2442E+00
 PARAMETER:  6.8151E-02  6.2187E-02  3.6718E-01  1.2416E-01  1.6762E-01  1.3224E-01 -2.9306E+00  2.7900E-01  1.8175E-01  1.5238E-01
             3.1851E-01
 GRADIENT:  -1.2496E-01 -2.8610E-01 -1.3515E-01 -2.1432E-02  2.9979E-01 -6.1731E-02  3.1786E-02  2.0524E-02  1.4839E-01  8.3917E-04
             1.4182E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1615.84448735745        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  9.6872E-01  9.6431E-01  1.3061E+00  1.0235E+00  1.0701E+00  1.0329E+00  1.0000E-02  1.1947E+00  1.0867E+00  1.0544E+00
             1.2438E+00
 PARAMETER:  6.8219E-02  6.3654E-02  3.6707E-01  1.2325E-01  1.6778E-01  1.3240E-01 -4.6145E+00  2.7791E-01  1.8314E-01  1.5300E-01
             3.1820E-01
 GRADIENT:  -2.7716E+01 -2.2691E+00 -1.4951E-01 -5.3819E+00 -8.2076E-01 -4.9051E+00  0.0000E+00 -1.1415E-01 -6.1942E-01 -5.8819E-02
            -2.6660E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1616.07872384317        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      649
 NPARAMETR:  9.8164E-01  9.3880E-01  1.3366E+00  1.0449E+00  1.0711E+00  1.0448E+00  1.0000E-02  1.2096E+00  1.0689E+00  1.0558E+00
             1.2450E+00
 PARAMETER:  8.1466E-02  3.6849E-02  3.9016E-01  1.4392E-01  1.6867E-01  1.4387E-01 -4.7179E+00  2.9031E-01  1.6663E-01  1.5431E-01
             3.1915E-01
 GRADIENT:  -7.3208E-02  1.8758E-01  1.5248E-01  4.6344E-02 -3.0331E-01  2.4288E-02  0.0000E+00 -6.6142E-03 -2.7237E-02  1.2635E-02
            -1.4422E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1616.07886744767        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      741
 NPARAMETR:  9.8164E-01  9.3471E-01  1.3368E+00  1.0476E+00  1.0698E+00  1.0447E+00  1.0000E-02  1.2069E+00  1.0661E+00  1.0549E+00
             1.2454E+00
 PARAMETER:  8.1473E-02  3.2482E-02  3.9026E-01  1.4650E-01  1.6749E-01  1.4377E-01 -4.6946E+00  2.8807E-01  1.6398E-01  1.5347E-01
             3.1945E-01
 GRADIENT:  -1.4089E-02  4.2911E-03  9.2551E-03 -4.3734E-03 -1.1841E-02 -5.7893E-03  0.0000E+00 -1.1013E-03 -2.3209E-03 -1.9386E-03
             4.4384E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      741
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.8964E-04 -6.4289E-04 -2.5775E-02 -3.6621E-03 -3.1772E-02
 SE:             2.9803E-02  2.0651E-04  1.3548E-02  2.8801E-02  2.2207E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7886E-01  1.8519E-03  5.7101E-02  8.9882E-01  1.5250E-01

 ETASHRINKSD(%)  1.5604E-01  9.9308E+01  5.4613E+01  3.5141E+00  2.5605E+01
 ETASHRINKVR(%)  3.1184E-01  9.9995E+01  7.9400E+01  6.9047E+00  4.4654E+01
 EBVSHRINKSD(%)  5.8459E-01  9.9373E+01  5.7663E+01  4.0620E+00  2.3609E+01
 EBVSHRINKVR(%)  1.1658E+00  9.9996E+01  8.2075E+01  7.9590E+00  4.1645E+01
 RELATIVEINF(%)  9.8577E+01  3.3304E-04  5.2210E+00  9.4258E+00  1.5102E+01
 EPSSHRINKSD(%)  4.1927E+01
 EPSSHRINKVR(%)  6.6276E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1616.0788674476712     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -880.92804088393302     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1616.079       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  9.35E-01  1.34E+00  1.05E+00  1.07E+00  1.04E+00  1.00E-02  1.21E+00  1.07E+00  1.05E+00  1.25E+00
 


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
+        1.05E+03
 
 TH 2
+       -1.94E+01  5.37E+02
 
 TH 3
+        3.16E+00  7.11E+01  8.35E+01
 
 TH 4
+       -1.73E+01  5.60E+02 -6.47E+00  8.20E+02
 
 TH 5
+        4.23E-01 -2.19E+02 -1.55E+02 -4.70E+01  5.17E+02
 
 TH 6
+        1.22E+00 -2.87E+00  8.30E-01 -3.11E+00 -2.36E+00  1.78E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -4.93E-01 -1.61E+01 -2.08E+01 -1.79E+00  6.82E-01 -2.25E-02  0.00E+00  1.83E+01
 
 TH 9
+        3.12E+00 -1.03E+02  2.70E+00  5.79E+00  4.11E+00 -4.27E-01  0.00E+00 -9.11E-01  1.47E+02
 
 TH10
+       -1.76E-02 -2.65E+00 -4.72E+00 -2.23E+00 -5.88E+01  9.59E-01  0.00E+00  8.22E+00  2.24E+00  7.02E+01
 
 TH11
+       -7.63E+00 -2.36E+01 -8.43E+00 -1.26E+01 -1.16E+00  3.22E+00  0.00E+00  3.42E+00  1.23E+01  1.54E+01  1.44E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       12.859
Stop Time:
Sat Sep 25 11:31:39 CDT 2021
