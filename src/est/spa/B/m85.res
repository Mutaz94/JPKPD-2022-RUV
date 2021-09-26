Sat Sep 25 07:35:27 CDT 2021
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
$DATA ../../../../data/spa/B/dat85.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1639.98048671835        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.5477E+01 -4.9757E+01 -1.1352E+01 -6.1907E+01 -3.9364E+00  3.1913E+01 -1.5193E+01  9.7379E+00 -1.4451E+01  8.9637E+00
            -3.0822E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1645.88032028807        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8378E-01  1.1203E+00  1.1140E+00  9.8758E-01  1.1127E+00  8.9256E-01  1.2089E+00  8.8103E-01  1.0771E+00  9.1867E-01
             1.0948E+00
 PARAMETER:  8.3645E-02  2.1361E-01  2.0797E-01  8.7498E-02  2.0680E-01 -1.3659E-02  2.8972E-01 -2.6668E-02  1.7425E-01  1.5174E-02
             1.9057E-01
 GRADIENT:   5.6267E+01  1.9844E+01  1.2353E+01  1.8364E+01  2.3708E+01 -9.7435E+00  5.2001E+00 -1.6499E+00  6.7841E+00 -1.3313E+01
             1.9611E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1648.51969623597        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8082E-01  9.0754E-01  8.5479E-01  1.1064E+00  9.0656E-01  9.2511E-01  1.5916E+00  2.9787E-01  9.2133E-01  8.5305E-01
             1.1006E+00
 PARAMETER:  8.0638E-02  2.9849E-03 -5.6895E-02  2.0113E-01  1.9018E-03  2.2154E-02  5.6473E-01 -1.1111E+00  1.8058E-02 -5.8943E-02
             1.9590E-01
 GRADIENT:   4.3348E+01  5.5594E+00 -3.6737E+01  5.4392E+01  5.4483E+01  4.3719E+00  1.5874E+01  8.2997E-01  5.5792E+00  5.1466E-01
             1.2682E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1650.25577928449        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.6642E-01  9.9238E-01  8.1291E-01  1.0302E+00  8.9958E-01  9.1832E-01  1.3498E+00  3.1043E-01  9.4909E-01  8.2789E-01
             1.0750E+00
 PARAMETER:  6.5841E-02  9.2352E-02 -1.0713E-01  1.2976E-01 -5.8299E-03  1.4786E-02  3.9999E-01 -1.0698E+00  4.7750E-02 -8.8872E-02
             1.7234E-01
 GRADIENT:   5.3904E+00 -1.3960E+00 -6.8588E+00  7.7255E+00  9.0566E+00  1.2188E+00  2.1014E+00  8.3578E-01  6.9897E-01  8.7808E-01
             2.2843E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1650.27369865915        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.6458E-01  1.0152E+00  7.9074E-01  1.0117E+00  8.9391E-01  9.1707E-01  1.3106E+00  2.4936E-01  9.5770E-01  8.1854E-01
             1.0710E+00
 PARAMETER:  6.3941E-02  1.1507E-01 -1.3478E-01  1.1159E-01 -1.2148E-02  1.3431E-02  3.7050E-01 -1.2889E+00  5.6775E-02 -1.0024E-01
             1.6856E-01
 GRADIENT:   1.5252E-02 -1.3240E+00 -2.3345E+00  5.5940E-01  1.5042E+00  4.7795E-01  2.2793E-01  5.3870E-01 -4.5499E-02  7.8548E-01
             6.3780E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1650.34587083050        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.6367E-01  1.0296E+00  7.7038E-01  9.9920E-01  8.8717E-01  9.1578E-01  1.2895E+00  1.1590E-01  9.6392E-01  8.0948E-01
             1.0682E+00
 PARAMETER:  6.2995E-02  1.2920E-01 -1.6087E-01  9.9195E-02 -1.9724E-02  1.2023E-02  3.5425E-01 -2.0550E+00  6.3251E-02 -1.1137E-01
             1.6595E-01
 GRADIENT:  -2.8720E+00 -9.5960E-01  8.3409E-01 -3.7363E+00 -2.8951E+00 -2.2971E-01 -7.9652E-01  1.1314E-01 -3.6147E-01  4.1766E-01
            -5.2571E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1650.77535995331        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  9.7738E-01  9.8786E-01  7.9001E-01  1.0324E+00  8.8001E-01  9.2367E-01  1.3528E+00  1.0000E-02  9.4504E-01  8.1915E-01
             1.0698E+00
 PARAMETER:  7.7120E-02  8.7784E-02 -1.3571E-01  1.3185E-01 -2.7819E-02  2.0595E-02  4.0217E-01 -5.5620E+00  4.3467E-02 -9.9492E-02
             1.6748E-01
 GRADIENT:  -6.9515E-02  6.5248E-02  5.7564E-02 -4.3812E-01 -6.8352E-01 -1.4173E-01 -3.8797E-01  0.0000E+00 -2.4822E-01  2.7049E-01
            -3.7139E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1650.77961071442        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      690
 NPARAMETR:  9.7728E-01  9.6778E-01  7.9485E-01  1.0446E+00  8.7402E-01  9.2384E-01  1.3798E+00  1.0000E-02  9.3720E-01  8.1612E-01
             1.0710E+00
 PARAMETER:  7.7022E-02  6.7250E-02 -1.2960E-01  1.4365E-01 -3.4649E-02  2.0787E-02  4.2195E-01 -5.9956E+00  3.5138E-02 -1.0319E-01
             1.6856E-01
 GRADIENT:  -1.5574E-03 -5.9083E-03 -5.1181E-03 -3.8153E-03  6.2039E-03 -2.8591E-04  5.9596E-04  0.0000E+00  1.9672E-03  5.2585E-04
             6.6653E-04

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1650.77961071442        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  9.7728E-01  9.6778E-01  7.9485E-01  1.0446E+00  8.7402E-01  9.2384E-01  1.3798E+00  1.0000E-02  9.3720E-01  8.1612E-01
             1.0710E+00
 PARAMETER:  7.7022E-02  6.7250E-02 -1.2960E-01  1.4365E-01 -3.4649E-02  2.0787E-02  4.2195E-01 -5.9956E+00  3.5138E-02 -1.0319E-01
             1.6856E-01
 GRADIENT:  -1.5574E-03 -5.9083E-03 -5.1181E-03 -3.8153E-03  6.2039E-03 -2.8591E-04  5.9596E-04  0.0000E+00  1.9672E-03  5.2585E-04
             6.6653E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      712
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.1586E-05  2.8023E-03 -4.5150E-04 -6.3759E-03 -1.2830E-02
 SE:             2.9798E-02  2.2253E-02  1.8367E-04  2.4065E-02  2.1655E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9889E-01  8.9979E-01  1.3963E-02  7.9106E-01  5.5355E-01

 ETASHRINKSD(%)  1.7311E-01  2.5451E+01  9.9385E+01  1.9378E+01  2.7452E+01
 ETASHRINKVR(%)  3.4593E-01  4.4424E+01  9.9996E+01  3.5000E+01  4.7369E+01
 EBVSHRINKSD(%)  5.5994E-01  2.5009E+01  9.9425E+01  1.9399E+01  2.6535E+01
 EBVSHRINKVR(%)  1.1167E+00  4.3763E+01  9.9997E+01  3.5034E+01  4.6029E+01
 RELATIVEINF(%)  9.8519E+01  4.7085E+00  4.3281E-04  6.4423E+00  5.4603E+00
 EPSSHRINKSD(%)  4.2669E+01
 EPSSHRINKVR(%)  6.7132E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1650.7796107144241     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -915.62878415068587     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.93
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1650.780       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.77E-01  9.68E-01  7.95E-01  1.04E+00  8.74E-01  9.24E-01  1.38E+00  1.00E-02  9.37E-01  8.16E-01  1.07E+00
 


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
+        1.35E+03
 
 TH 2
+       -8.94E+00  3.68E+02
 
 TH 3
+        2.16E+01  2.10E+02  6.49E+02
 
 TH 4
+       -1.27E+01  2.68E+02 -2.82E+02  7.66E+02
 
 TH 5
+       -6.66E+00 -3.37E+02 -7.41E+02  3.26E+02  1.17E+03
 
 TH 6
+        8.56E-01  2.17E-01  3.45E+00 -3.61E+00 -1.12E+00  2.28E+02
 
 TH 7
+        1.15E+00  3.08E+01 -9.18E+00 -1.24E+01 -1.65E-02  1.74E-02  3.91E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.13E+00 -2.19E+01 -3.12E+01  3.13E+01  2.05E+00 -2.58E-01  1.22E+01  0.00E+00  1.06E+02
 
 TH10
+       -1.04E+00 -4.18E+00 -5.82E+01 -2.67E+01 -5.34E+01  1.32E-01  1.12E+01  0.00E+00  1.35E+01  8.98E+01
 
 TH11
+       -9.04E+00 -1.23E+01 -4.02E+01 -4.96E+00  8.44E+00  2.58E+00  4.40E+00  0.00E+00  1.01E+01  2.36E+01  1.93E+02
 
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
 #CPUT: Total CPU Time in Seconds,       12.821
Stop Time:
Sat Sep 25 07:35:41 CDT 2021
