Sat Sep 25 10:28:25 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat34.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1582.46872738134        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5870E+02 -1.1371E+02 -3.5942E+01 -7.4946E+01  1.1233E+02 -3.2118E+01 -8.2124E+00 -1.1234E+01  3.3840E+01 -3.3885E+01
            -3.6798E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1601.95943439614        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6089E-01  1.2578E+00  9.1656E-01  9.3283E-01  9.8006E-01  1.0380E+00  1.0998E+00  1.2124E+00  5.9916E-01  1.1691E+00
             1.0514E+00
 PARAMETER:  6.0100E-02  3.2939E-01  1.2876E-02  3.0471E-02  7.9854E-02  1.3733E-01  1.9509E-01  2.9257E-01 -4.1223E-01  2.5623E-01
             1.5014E-01
 GRADIENT:   5.9042E+01  6.1834E+01  5.7094E-01  5.5335E+01 -1.8450E+01 -8.5858E+00 -6.4038E+00 -4.0395E+00 -1.2278E+01  7.4702E+00
            -1.4379E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1602.32758317321        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      199
 NPARAMETR:  9.6008E-01  1.2584E+00  9.2688E-01  9.3205E-01  9.8657E-01  1.0396E+00  1.1290E+00  1.2396E+00  5.9361E-01  1.1578E+00
             1.0545E+00
 PARAMETER:  5.9266E-02  3.2986E-01  2.4073E-02  2.9632E-02  8.6480E-02  1.3887E-01  2.2137E-01  3.1476E-01 -4.2153E-01  2.4652E-01
             1.5304E-01
 GRADIENT:   2.4032E+01  4.9274E+01 -4.9917E-02  4.6337E+01 -1.6182E+01 -1.2453E+01 -2.0913E+00 -4.0185E+00 -1.1602E+01  4.7805E+00
            -1.3026E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1603.47644670378        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.3384E-01  1.2582E+00  9.3038E-01  9.0187E-01  1.0127E+00  1.0585E+00  1.1132E+00  1.2394E+00  5.9374E-01  1.1414E+00
             1.0870E+00
 PARAMETER:  3.1553E-02  3.2969E-01  2.7841E-02 -3.2860E-03  1.1265E-01  1.5681E-01  2.0721E-01  3.1460E-01 -4.2131E-01  2.3227E-01
             1.8346E-01
 GRADIENT:   4.1129E-01  1.4508E+01  1.3616E-01 -2.1221E+00 -1.4683E+00  5.3810E-01  3.9852E-01 -5.7119E+00 -9.5996E+00 -3.7235E-01
             2.2035E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1603.77768108210        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      480
 NPARAMETR:  9.4952E-01  1.2582E+00  9.4889E-01  9.0678E-01  1.0227E+00  1.0705E+00  1.1188E+00  1.2394E+00  5.9374E-01  1.1584E+00
             1.0878E+00
 PARAMETER:  4.8202E-02  3.2969E-01  4.7538E-02  2.1466E-03  1.2244E-01  1.6817E-01  2.1229E-01  3.1460E-01 -4.2131E-01  2.4703E-01
             1.8418E-01
 GRADIENT:   1.8387E+00  8.6290E+00 -4.1983E-01  7.0394E-01 -7.4399E-02  1.0622E-01  1.1960E-01 -6.3057E+00 -1.0040E+01  1.0486E-01
             2.1004E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1605.19253962846        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  9.3946E-01  1.2480E+00  9.8168E-01  9.0649E-01  1.0225E+00  1.0619E+00  1.0825E+00  1.3961E+00  7.1288E-01  1.1403E+00
             1.0783E+00
 PARAMETER:  3.7550E-02  3.2152E-01  8.1508E-02  1.8250E-03  1.2228E-01  1.6005E-01  1.7931E-01  4.3366E-01 -2.3844E-01  2.3128E-01
             1.7541E-01
 GRADIENT:   1.3947E+01  3.2660E+00 -1.7231E+00 -8.9210E+00 -6.0416E+00  2.2052E+00  3.2715E+00 -3.0729E-01  1.2667E-01  1.1751E+00
             3.7482E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1605.56704446053        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  9.5098E-01  1.2531E+00  1.0298E+00  9.1280E-01  1.0423E+00  1.0734E+00  1.0449E+00  1.4633E+00  7.4187E-01  1.1416E+00
             1.0789E+00
 PARAMETER:  4.9742E-02  3.2559E-01  1.2937E-01  8.7578E-03  1.4145E-01  1.7083E-01  1.4389E-01  4.8067E-01 -1.9859E-01  2.3241E-01
             1.7596E-01
 GRADIENT:   5.8437E+00 -9.3805E-01 -1.9771E+00  1.7485E-01 -4.1738E+00  1.0772E+00 -1.2789E+00 -2.2682E-01  6.8838E-01 -1.4094E+00
            -1.0319E-01

0ITERATION NO.:   34    OBJECTIVE VALUE:  -1605.58352691203        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      976
 NPARAMETR:  9.5090E-01  1.2531E+00  1.0336E+00  9.1291E-01  1.0440E+00  1.0733E+00  1.0447E+00  1.4683E+00  7.4187E-01  1.1434E+00
             1.0790E+00
 PARAMETER:  4.9656E-02  3.2565E-01  1.3291E-01  8.8777E-03  1.4299E-01  1.7074E-01  1.4361E-01  4.8414E-01 -1.9858E-01  2.3388E-01
             1.7601E-01
 GRADIENT:  -6.6697E+05  2.0481E+05 -1.9786E+00  3.3348E+05 -4.0863E+00  1.0332E+00 -1.3208E+00  1.3765E+05 -3.3588E+05 -1.3841E+00
            -1.8947E+05
 NUMSIGDIG:         3.3         3.3         0.9         3.3         1.2         1.5         0.9         3.3         3.3         1.1
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      976
 NO. OF SIG. DIGITS IN FINAL EST.:  0.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.2209E-03 -6.6610E-03 -3.9308E-02 -9.9292E-04 -3.5996E-02
 SE:             2.9721E-02  2.3240E-02  1.4837E-02  1.9286E-02  2.1685E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4043E-01  7.7441E-01  8.0659E-03  9.5894E-01  9.6929E-02

 ETASHRINKSD(%)  4.3001E-01  2.2144E+01  5.0294E+01  3.5391E+01  2.7352E+01
 ETASHRINKVR(%)  8.5818E-01  3.9384E+01  7.5293E+01  5.8257E+01  4.7222E+01
 EBVSHRINKSD(%)  4.6763E-01  2.2534E+01  5.4939E+01  3.6525E+01  2.3965E+01
 EBVSHRINKVR(%)  9.3307E-01  3.9990E+01  7.9695E+01  5.9710E+01  4.2186E+01
 RELATIVEINF(%)  9.8600E+01  1.7729E+00  1.6748E+00  1.1088E+00  1.3118E+01
 EPSSHRINKSD(%)  4.4530E+01
 EPSSHRINKVR(%)  6.9231E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1605.5835269120253     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -870.43270034828709     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1605.584       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.51E-01  1.25E+00  1.03E+00  9.13E-01  1.04E+00  1.07E+00  1.04E+00  1.47E+00  7.42E-01  1.14E+00  1.08E+00
 


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
+        1.84E+09
 
 TH 2
+        8.57E+02  1.00E+08
 
 TH 3
+       -1.28E+05  3.00E+04  8.84E+08
 
 TH 4
+        8.44E+04 -1.92E+04  1.34E+05  2.00E+09
 
 TH 5
+       -1.00E+04  2.23E+03 -1.72E+02  1.06E+04  4.60E+02
 
 TH 6
+       -4.33E+03  1.01E+03  1.12E+00  4.50E+03 -2.53E+00  1.68E+02
 
 TH 7
+        1.52E+03 -3.24E+02  8.09E+08 -1.60E+03 -3.24E+00 -4.33E-01  7.41E+08
 
 TH 8
+        4.97E+02  5.74E+07  1.71E+04 -1.13E+04  1.34E+03  5.79E+02 -2.01E+02  3.29E+07
 
 TH 9
+        1.19E+09 -1.01E+04 -8.29E+04  5.45E+04 -6.48E+03 -2.79E+03  1.01E+03 -5.80E+03  7.68E+08
 
 TH10
+       -4.85E+02  1.08E+02 -1.14E+01  5.08E+02 -5.71E+01 -1.82E-01 -5.08E-01  7.18E+01 -3.12E+02  6.29E+01
 
 TH11
+        3.01E+04 -7.03E+03 -6.39E+08 -9.62E+08 -5.03E+03 -2.16E+03 -5.85E+08 -4.02E+03  1.95E+04 -2.34E+02  4.62E+08
 
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
 #CPUT: Total CPU Time in Seconds,       17.658
Stop Time:
Sat Sep 25 10:28:44 CDT 2021
