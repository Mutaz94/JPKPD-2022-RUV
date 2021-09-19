Sat Sep 18 15:56:43 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat64.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA
0.9 FIX ;     IIV CL
0.9 FIX  ;     IIV V
0.9 FIX ;      IIV KA
$SIGMA  1  FIX;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E+00
 0.0000E+00   0.9000E+00
 0.0000E+00   0.0000E+00   0.9000E+00
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

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   9081.90875798555        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   4.6027E+02 -3.1079E+02  1.7095E+02 -1.6548E+03  1.7175E+03 -3.3463E+02 -2.0121E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -547.188718012112        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  1.1816E+00  2.5710E+00  4.7319E+00  8.1741E-01  1.2797E+00  3.8135E-01  1.4896E+01
 PARAMETER:  2.6689E-01  1.0443E+00  1.6543E+00 -1.0161E-01  3.4661E-01 -8.6404E-01  2.8011E+00
 GRADIENT:  -3.7425E+01  2.3547E+01 -1.8192E-01  4.4664E+01  3.0822E+01  5.3819E-03  2.2247E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -598.131895470471        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      111
 NPARAMETR:  1.2431E+00  2.0346E+00  1.4162E+01  6.1811E-01  1.0532E+00  2.5175E+00  1.0924E+01
 PARAMETER:  3.1759E-01  8.1032E-01  2.7506E+00 -3.8108E-01  1.5185E-01  1.0233E+00  2.4910E+00
 GRADIENT:   1.1949E+01  2.2469E+01  4.0503E-01 -4.2111E+01 -2.5546E+01 -1.1043E-02  8.2604E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -616.890293580496        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.1985E+00  1.6406E+00  3.0900E+01  8.3515E-01  1.1369E+00  4.3239E+00  8.7525E+00
 PARAMETER:  2.8110E-01  5.9504E-01  3.5308E+00 -8.0149E-02  2.2832E-01  1.5642E+00  2.2693E+00
 GRADIENT:   5.0387E-01  1.2839E-02  9.5032E-02 -1.7415E-01 -1.1883E+00  3.2475E-02 -2.2324E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -616.932422127879        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  1.2084E+00  1.6710E+00  2.8816E+01  8.3767E-01  1.1442E+00  4.1561E+00  8.8058E+00
 PARAMETER:  2.8928E-01  6.1343E-01  3.4609E+00 -7.7134E-02  2.3469E-01  1.5246E+00  2.2754E+00
 GRADIENT:   1.3219E-01  1.8698E-01  1.1425E-01  4.2218E-01  2.1397E-01  2.3013E-02  1.2054E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -617.770766685945        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      392
 NPARAMETR:  1.2305E+00  1.6798E+00  3.2079E+00  8.3731E-01  1.1369E+00  4.7284E-01  8.7690E+00
 PARAMETER:  3.0743E-01  6.1866E-01  1.2656E+00 -7.7567E-02  2.2833E-01 -6.4899E-01  2.2712E+00
 GRADIENT:   2.5961E-01 -7.7544E-01 -1.4851E-01  1.1385E-01 -6.2701E-01  2.6744E-01 -3.2581E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -617.990255809909        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      509
 NPARAMETR:  1.2319E+00  1.6819E+00  2.6137E+00  8.3495E-01  1.1354E+00  7.4827E-02  8.7528E+00
 PARAMETER:  3.0853E-01  6.1992E-01  1.0608E+00 -8.0382E-02  2.2695E-01 -2.4926E+00  2.2694E+00
 GRADIENT:   1.2896E-01  2.8083E-01 -1.7128E-01 -5.4658E-01 -6.1933E-01  2.4701E-02 -1.2328E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -618.004950481683        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  1.2315E+00  1.6807E+00  2.6457E+00  8.3633E-01  1.1375E+00  1.0914E-02  8.7651E+00
 PARAMETER:  3.0823E-01  6.1918E-01  1.0729E+00 -7.8728E-02  2.2881E-01 -4.4177E+00  2.2708E+00
 GRADIENT:  -2.6716E-02  1.6802E-02  3.6317E-03 -2.2539E-02  9.0726E-03  4.8722E-04 -2.6801E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -618.004993319853        NO. OF FUNC. EVALS.:  37
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  1.2316E+00  1.6806E+00  2.6436E+00  8.3635E-01  1.1374E+00  1.0000E-02  8.7652E+00
 PARAMETER:  3.0832E-01  6.1913E-01  1.0721E+00 -7.8704E-02  2.2876E-01 -4.7367E+00  2.2708E+00
 GRADIENT:  -4.6192E-03  1.1402E-02 -2.7617E-03 -1.2462E-02 -4.4983E-03  0.0000E+00 -2.4720E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      661
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.2634E-03 -4.2870E-02 -2.3271E-04
 SE:             8.8997E-02  9.1694E-02  1.0457E-04
 N:                     100         100         100

 P VAL.:         9.5284E-01  6.4012E-01  2.6048E-02

 ETASHRINKSD(%)  5.7165E+00  2.8589E+00  9.9889E+01
 ETASHRINKVR(%)  1.1106E+01  5.6360E+00  1.0000E+02
 EBVSHRINKSD(%)  5.6138E+00  3.3323E+00  9.9820E+01
 EBVSHRINKVR(%)  1.0912E+01  6.5536E+00  1.0000E+02
 RELATIVEINF(%)  8.8843E+01  9.1451E+01  3.1844E-04
 EPSSHRINKSD(%)  1.7972E+01
 EPSSHRINKVR(%)  3.2714E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -618.00499331985350     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       117.14583324388468     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     5.79
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     2.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -618.005       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.23E+00  1.68E+00  2.64E+00  8.36E-01  1.14E+00  1.00E-02  8.77E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-01
 
 ETA2
+        0.00E+00  9.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.49E-01
 
 ETA2
+        0.00E+00  9.49E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.49E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        9.93E+01
 
 TH 2
+        3.34E+00  1.73E+00
 
 TH 3
+        3.50E-01 -1.41E-02  2.63E-03
 
 TH 4
+       -6.86E+01 -5.75E+00 -2.97E-01  2.71E+02
 
 TH 5
+       -9.27E+00  3.22E+00 -3.10E-01 -5.63E+01  8.95E+01
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -3.22E+00 -9.69E-01 -7.71E-03  3.64E+00  1.31E+00  0.00E+00  6.81E-01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        9.12E+01
 
 TH 2
+       -1.82E+00  2.81E+01
 
 TH 3
+        3.32E-01 -6.10E-01  6.41E-01
 
 TH 4
+       -4.37E+00 -1.54E+00 -2.93E-01  2.28E+02
 
 TH 5
+        2.85E+00 -1.02E-01 -4.00E-01 -1.05E+01  1.32E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -3.69E+00 -1.81E+00 -6.71E-02  2.71E+00  2.04E+00  0.00E+00  4.77E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        9.34E+01
 
 TH 2
+        3.48E+01  2.85E+01
 
 TH 3
+        9.20E-01  5.57E-01  2.02E-01
 
 TH 4
+        6.52E+01  3.15E+01  3.23E-01  2.67E+02
 
 TH 5
+        2.78E+01  5.15E+01 -2.58E-01  8.00E+01  2.14E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -9.74E+00 -9.15E+00 -3.37E-01  1.08E+01  1.17E+01  0.00E+00  4.16E+01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,        7.906
Stop Time:
Sat Sep 18 15:56:53 CDT 2021
