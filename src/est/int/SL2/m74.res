Sat Sep 25 01:33:58 CDT 2021
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
$DATA ../../../../data/int/SL2/dat74.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      997
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

 TOT. NO. OF OBS RECS:      897
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1479.24436440086        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2634E+01 -3.4342E+01  1.9404E+02 -1.7195E+01  6.6692E+01  1.9929E+01 -5.3306E+01 -2.2044E+02 -7.3457E+01 -1.1901E-01
            -4.4496E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2896.75645822668        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0061E+00  1.2926E+00  9.0097E-01  8.6539E-01  1.1547E+00  8.7842E-01  9.8697E-01  8.9168E-01  1.1315E+00  8.6297E-01
             2.0199E+00
 PARAMETER:  1.0603E-01  3.5666E-01 -4.2885E-03 -4.4578E-02  2.4382E-01 -2.9628E-02  8.6887E-02 -1.4648E-02  2.2354E-01 -4.7376E-02
             8.0307E-01
 GRADIENT:   1.6289E+00  2.5706E+00 -1.8080E+01 -7.9318E+00  2.2150E+01 -3.1211E+01  8.9423E+00  3.2967E+00 -1.9298E+01 -2.5573E+01
            -3.7690E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2914.09569943750        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0194E+00  1.5624E+00  6.9812E-01  7.6017E-01  1.2611E+00  1.0382E+00  7.3283E-01  9.8187E-02  1.4076E+00  1.2504E+00
             2.1649E+00
 PARAMETER:  1.1917E-01  5.4622E-01 -2.5937E-01 -1.7421E-01  3.3200E-01  1.3750E-01 -2.1084E-01 -2.2209E+00  4.4190E-01  3.2345E-01
             8.7237E-01
 GRADIENT:   2.7899E+01  4.3508E+01 -5.9394E+01  6.4315E+01 -2.1572E+01  3.2078E+01 -6.9868E-01  1.1071E-02  9.8272E+00  1.3313E+01
            -1.9576E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2932.87709037552        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0103E+00  1.4970E+00  9.5359E-01  7.5067E-01  1.3237E+00  9.4409E-01  7.8017E-01  2.9232E-02  1.3253E+00  1.1251E+00
             2.3270E+00
 PARAMETER:  1.1020E-01  5.0349E-01  5.2480E-02 -1.8679E-01  3.8041E-01  4.2461E-02 -1.4825E-01 -3.4325E+00  3.8167E-01  2.1785E-01
             9.4459E-01
 GRADIENT:   2.4493E+00 -1.1263E+00  3.5549E+00 -2.9330E+00 -2.9517E+00  2.0780E-02 -3.0499E-01  3.9520E-03 -4.9609E-01 -2.1823E+00
            -1.8834E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2933.22919545518        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0091E+00  1.5287E+00  9.2652E-01  7.3429E-01  1.3442E+00  9.4463E-01  7.6815E-01  2.7593E-02  1.3504E+00  1.1570E+00
             2.3285E+00
 PARAMETER:  1.0904E-01  5.2443E-01  2.3676E-02 -2.0885E-01  3.9583E-01  4.3034E-02 -1.6377E-01 -3.4902E+00  4.0037E-01  2.4580E-01
             9.4524E-01
 GRADIENT:  -7.0620E-01  7.0496E-01 -1.2780E+00  2.4639E+00  5.0820E-01  1.3136E-01  2.0685E-01  3.6764E-03  5.7129E-01  9.6728E-01
             1.0853E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2933.33413973277        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0093E+00  1.5385E+00  9.3085E-01  7.2627E-01  1.3537E+00  9.4428E-01  7.6552E-01  2.7831E-02  1.3563E+00  1.1575E+00
             2.3273E+00
 PARAMETER:  1.0931E-01  5.3084E-01  2.8344E-02 -2.1983E-01  4.0284E-01  4.2672E-02 -1.6721E-01 -3.4816E+00  4.0478E-01  2.4623E-01
             9.4472E-01
 GRADIENT:   3.9409E-02  1.3056E-01  3.2398E-02  1.7512E-01 -6.1302E-02  2.8011E-02  5.9833E-03  3.5943E-03  1.6470E-02  2.1693E-02
             1.2626E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2935.21445509481        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  1.0150E+00  1.8948E+00  8.1094E-01  5.0713E-01  1.6333E+00  9.4948E-01  7.0929E-01  1.1316E-02  1.6982E+00  1.3304E+00
             2.3232E+00
 PARAMETER:  1.1489E-01  7.3909E-01 -1.0956E-01 -5.7899E-01  5.9057E-01  4.8160E-02 -2.4349E-01 -4.3815E+00  6.2955E-01  3.8552E-01
             9.4295E-01
 GRADIENT:   3.9848E+00  8.8743E+00  1.0506E+00  4.8695E+00  3.8522E-01  1.2428E+00  3.7108E-01  3.3331E-04  1.7758E-01  1.2148E+00
             3.4658E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2935.34236393967        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  1.0132E+00  1.9799E+00  7.4239E-01  4.4466E-01  1.6976E+00  9.4639E-01  7.0054E-01  1.0000E-02  1.8170E+00  1.3639E+00
             2.3171E+00
 PARAMETER:  1.1313E-01  7.8306E-01 -1.9788E-01 -7.1045E-01  6.2919E-01  4.4905E-02 -2.5590E-01 -4.9134E+00  6.9721E-01  4.1032E-01
             9.4032E-01
 GRADIENT:   1.8342E-02  3.6774E-02 -1.2533E-03  3.5844E-02 -2.2579E-03  4.7983E-03  4.2675E-03  0.0000E+00  4.5714E-03  1.7658E-03
             2.0385E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2935.34237555014        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:      772
 NPARAMETR:  1.0132E+00  1.9812E+00  7.4164E-01  4.4380E-01  1.6985E+00  9.4640E-01  7.0046E-01  1.0000E-02  1.8188E+00  1.3645E+00
             2.3171E+00
 PARAMETER:  1.1312E-01  7.8368E-01 -1.9909E-01 -7.1237E-01  6.2977E-01  4.4887E-02 -2.5602E-01 -4.9214E+00  6.9814E-01  4.1074E-01
             9.4030E-01
 GRADIENT:  -9.7793E-03 -1.3111E-02 -7.8986E-03  5.7191E-04  1.5655E-02 -3.4607E-03 -6.3775E-04  0.0000E+00 -2.1165E-03 -4.0594E-03
            -7.6701E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      772
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4649E-03 -3.4178E-02 -6.1412E-05  2.9682E-02 -2.3347E-02
 SE:             2.9501E-02  2.2392E-02  5.7788E-05  2.2074E-02  2.5696E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6040E-01  1.2692E-01  2.8791E-01  1.7874E-01  3.6356E-01

 ETASHRINKSD(%)  1.1680E+00  2.4985E+01  9.9806E+01  2.6049E+01  1.3916E+01
 ETASHRINKVR(%)  2.3224E+00  4.3728E+01  1.0000E+02  4.5313E+01  2.5895E+01
 EBVSHRINKSD(%)  1.4059E+00  2.3917E+01  9.9797E+01  2.9042E+01  1.2061E+01
 EBVSHRINKVR(%)  2.7919E+00  4.2113E+01  1.0000E+02  4.9650E+01  2.2667E+01
 RELATIVEINF(%)  9.7168E+01  7.7662E+00  2.1330E-04  7.0629E+00  2.6440E+01
 EPSSHRINKSD(%)  1.6682E+01
 EPSSHRINKVR(%)  3.0581E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2935.3423755501435     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1286.7666469809608     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.51
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2935.342       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.98E+00  7.41E-01  4.44E-01  1.70E+00  9.46E-01  7.00E-01  1.00E-02  1.82E+00  1.36E+00  2.32E+00
 


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
+        1.18E+03
 
 TH 2
+       -9.84E+00  3.46E+02
 
 TH 3
+        2.51E+00  4.29E+01  7.00E+01
 
 TH 4
+       -1.68E+01  3.59E+02 -6.89E+01  8.55E+02
 
 TH 5
+       -2.63E+00 -7.74E+01 -4.25E+01  1.07E+02  1.59E+02
 
 TH 6
+        4.15E+00 -2.07E+00 -1.12E+00 -6.79E+00 -6.94E-01  2.13E+02
 
 TH 7
+        2.25E+00 -1.52E+01 -1.48E+00 -4.63E+00 -7.84E+00  5.85E-01  1.57E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.22E+00 -8.57E+00 -6.70E+00  4.15E+01  9.99E-01 -3.01E-01  1.78E+01  0.00E+00  2.15E+01
 
 TH10
+        2.69E-01 -1.41E+01 -7.25E-02  1.56E+01 -1.26E+01  5.17E-01  3.89E+00  0.00E+00  3.76E+00  6.26E+01
 
 TH11
+       -1.55E+01 -1.27E+01 -4.14E+00 -1.38E+01  1.42E+00  2.77E+00  7.19E+00  0.00E+00  1.80E+00  5.77E+00  2.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.775
Stop Time:
Sat Sep 25 01:34:27 CDT 2021
