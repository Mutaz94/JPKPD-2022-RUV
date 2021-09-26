Sat Sep 25 14:10:08 CDT 2021
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
$DATA ../../../../data/spa/D/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   2877.80278413876        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.4953E+02  6.9077E+00 -9.9705E+01 -9.7464E+01  2.5949E+02 -1.0429E+03 -3.4912E+02 -1.9664E+01 -7.2304E+02 -4.1027E+02
            -6.6985E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -875.761397186206        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2450E+00  7.4528E-01  2.7318E+00  2.3048E+00  2.9926E+00  2.9934E+00  2.5637E+00  8.5306E-01  5.1779E+00  2.1464E+00
             5.5228E+00
 PARAMETER:  3.1910E-01 -1.9399E-01  1.1050E+00  9.3499E-01  1.1962E+00  1.1964E+00  1.0414E+00 -5.8928E-02  1.7444E+00  8.6379E-01
             1.8089E+00
 GRADIENT:   9.3503E+00 -2.5746E+01 -6.0025E+00  4.0147E+01 -2.6924E-01  9.1774E+01  1.2696E+01  1.7005E-01  1.2477E+02  7.2336E-01
             9.8868E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -899.581055578195        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.8679E+00  3.9089E-01  1.8395E+00  1.6000E+00  1.9710E+01  2.5296E+00  9.1407E+00  5.9649E-01  2.5256E+00  1.4337E+01
             5.4595E+00
 PARAMETER:  7.2483E-01 -8.3934E-01  7.0949E-01  5.6997E-01  3.0811E+00  1.0281E+00  2.3127E+00 -4.1670E-01  1.0265E+00  2.7629E+00
             1.7974E+00
 GRADIENT:   1.7909E+02  1.6386E+01  1.8486E+01 -7.9141E+01 -9.0974E-01 -2.0788E+01  6.4449E+01 -6.6174E-01  2.2682E+01  1.0870E+01
             7.4172E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -971.064528080686        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.0900E+00  3.2288E-01  1.2793E+00  1.7680E+00  1.0679E+01  2.1780E+00  7.8643E+00  1.9633E+00  1.9043E+00  1.3202E+01
             5.3498E+00
 PARAMETER:  1.8622E-01 -1.0305E+00  3.4629E-01  6.6985E-01  2.4683E+00  8.7839E-01  2.1623E+00  7.7461E-01  7.4411E-01  2.6804E+00
             1.7771E+00
 GRADIENT:  -9.9593E+00  1.0066E+01  1.4908E+01  4.2046E+01 -1.9179E+00  3.1541E+01  6.1401E+00 -1.4583E+01  3.8512E+01  2.0210E+01
             6.7390E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1058.65350015211        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  7.7461E-01  2.9916E-02  1.3459E-01  8.3682E-01  1.1214E+01  1.4074E+00  4.8194E+00  2.6518E+00  6.0313E-01  8.6693E+00
             4.0509E+00
 PARAMETER: -1.5539E-01 -3.4094E+00 -1.9055E+00 -7.8150E-02  2.5172E+00  4.4173E-01  1.6726E+00  1.0752E+00 -4.0562E-01  2.2598E+00
             1.4989E+00
 GRADIENT:   3.0304E+01  3.8081E-01 -4.3076E+01  6.5099E+01 -3.6549E-01 -1.7706E+01  2.1880E-01  3.3546E+01  1.7891E+01  1.7232E+01
            -1.3884E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1096.59798943474        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  7.0736E-01  3.9113E-02  1.1703E-01  7.3733E-01  1.7253E+01  1.4744E+00  3.7147E+00  1.6664E+00  5.5271E-01  9.8531E+00
             4.6963E+00
 PARAMETER: -2.4621E-01 -3.1413E+00 -2.0454E+00 -2.0472E-01  2.9480E+00  4.8827E-01  1.4123E+00  6.1065E-01 -4.9293E-01  2.3878E+00
             1.6468E+00
 GRADIENT:  -1.8499E+01  3.9124E+00  1.0537E+01  2.4961E+01  8.1235E-01  1.0395E+00  3.5902E-01 -6.1008E+00  9.0943E+00  9.5882E-01
            -1.7627E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1097.54655954071        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      484
 NPARAMETR:  7.0554E-01  3.6948E-02  1.1288E-01  7.2397E-01  1.6249E+01  1.4627E+00  3.6449E+00  1.6488E+00  5.3519E-01  9.6569E+00
             4.6790E+00
 PARAMETER: -2.4879E-01 -3.1982E+00 -2.0814E+00 -2.2301E-01  2.8881E+00  4.8027E-01  1.3933E+00  6.0007E-01 -5.2512E-01  2.3677E+00
             1.6431E+00
 GRADIENT:  -1.7727E+01  2.5129E+00  4.2086E+00  2.2464E+01  8.8326E-01 -3.5641E+00  1.7287E-01 -6.7938E+00  8.2318E+00  7.8094E-01
            -6.1299E+00

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1097.54655954071        NO. OF FUNC. EVALS.:  40
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  7.0554E-01  3.6948E-02  1.1288E-01  7.2397E-01  1.6249E+01  1.4627E+00  3.6448E+00  1.6490E+00  5.3491E-01  9.6569E+00
             4.6790E+00
 PARAMETER: -2.4879E-01 -3.1982E+00 -2.0814E+00 -2.2301E-01  2.8881E+00  4.8027E-01  1.3933E+00  6.0007E-01 -5.2512E-01  2.3677E+00
             1.6431E+00
 GRADIENT:   6.6097E+03  2.6020E+02  4.0022E+02 -7.3722E+03 -5.6989E+02 -3.5232E+00  5.9187E+02 -6.7719E+00  8.2276E+00  3.4936E+02
            -5.0753E+02
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         1.4         3.3         0.8         0.0         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      524
 NO. OF SIG. DIGITS IN FINAL EST.:  0.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2774E-02 -8.5660E-03 -1.3617E-02 -2.5363E-02 -7.2469E-03
 SE:             2.9431E-02  1.3269E-03  2.3114E-02  1.3636E-02  4.1215E-03
 N:                     100         100         100         100         100

 P VAL.:         6.6426E-01  1.0831E-10  5.5579E-01  6.2884E-02  7.8692E-02

 ETASHRINKSD(%)  1.4036E+00  9.5555E+01  2.2566E+01  5.4318E+01  8.6193E+01
 ETASHRINKVR(%)  2.7874E+00  9.9802E+01  4.0039E+01  7.9132E+01  9.8094E+01
 EBVSHRINKSD(%)  2.6159E+00  9.3639E+01  2.2618E+01  5.0885E+01  8.5706E+01
 EBVSHRINKVR(%)  5.1633E+00  9.9595E+01  4.0121E+01  7.5877E+01  9.7957E+01
 RELATIVEINF(%)  3.4912E+01  2.4111E-01  4.1089E+00  1.4595E+00  1.9195E+00
 EPSSHRINKSD(%)  2.1319E+01
 EPSSHRINKVR(%)  3.8093E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1097.5465595407118     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -362.39573297697359     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.79
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1097.547       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         7.06E-01  3.69E-02  1.13E-01  7.24E-01  1.62E+01  1.46E+00  3.64E+00  1.65E+00  5.35E-01  9.66E+00  4.68E+00
 


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
+        1.34E+07
 
 TH 2
+       -1.01E+03  2.95E+07
 
 TH 3
+       -8.28E+02  1.11E+04  7.49E+06
 
 TH 4
+        2.25E+02 -5.34E+03 -1.09E+07  1.58E+07
 
 TH 5
+        2.03E+00  3.92E+01  2.01E+01 -2.97E+01  1.87E+02
 
 TH 6
+        1.48E+02  1.99E+02  1.17E+02 -1.74E+02 -4.84E-01  8.52E+01
 
 TH 7
+       -1.29E+01 -1.44E+01 -1.02E+01  1.39E+01  4.73E-02  4.52E+00  1.60E+04
 
 TH 8
+        1.86E+02  2.75E+02 -4.89E+01 -2.10E+02 -7.05E-01  3.46E+00  6.33E+00  4.21E+05
 
 TH 9
+       -5.66E+02 -9.45E+02 -3.21E+02  5.51E+02  2.13E+00  2.09E+06 -1.96E+01 -1.48E+06  5.22E+06
 
 TH10
+        2.25E+02  3.34E+02  1.69E+02 -2.45E+02 -8.51E-01  1.01E+00  7.77E+00  1.46E+00 -4.35E+00  7.90E+02
 
 TH11
+        1.68E+01  1.97E+01  7.56E+01 -4.82E+01 -1.04E-01 -2.27E+00  2.96E-01 -1.87E+00  1.65E+01 -5.14E+00  6.99E+03
 
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
 #CPUT: Total CPU Time in Seconds,       13.308
Stop Time:
Sat Sep 25 14:10:26 CDT 2021
