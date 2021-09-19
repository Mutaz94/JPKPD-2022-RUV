Sat Sep 18 15:57:10 CDT 2021
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
$DATA ../../../../data/spa/All/dat67.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12147.6341857733        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   6.2483E+01  1.7324E+01 -1.1148E+02 -1.9813E+01 -1.2881E+02 -1.7611E+02 -2.6021E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -493.501562552325        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  1.5462E+00  1.9244E+00  3.4515E+00  7.3134E-01  6.9814E-01  6.1686E-01  1.5132E+01
 PARAMETER:  5.3579E-01  7.5460E-01  1.3388E+00 -2.1288E-01 -2.5934E-01 -3.8312E-01  2.8168E+00
 GRADIENT:   1.4987E+01  3.3035E+00 -2.2335E+00  2.4831E+00  1.5612E+01  1.9213E-01  1.6115E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -521.879648686504        NO. OF FUNC. EVALS.:  53
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  1.2011E+00  1.4734E+00  1.3784E+01  6.1829E-01  5.2698E-01  2.8162E+00  1.1639E+01
 PARAMETER:  2.8321E-01  4.8756E-01  2.7235E+00 -3.8079E-01 -5.4059E-01  1.1354E+00  2.5544E+00
 GRADIENT:  -1.5607E+01 -5.5644E+00 -2.6122E-01 -1.4771E+01  4.7249E+00 -1.5595E-01  2.1981E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -525.736813213758        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  1.2435E+00  1.6695E+00  2.1530E+01  7.6811E-01  6.4894E-01  6.2683E+00  1.1006E+01
 PARAMETER:  3.1795E-01  6.1251E-01  3.1695E+00 -1.6383E-01 -3.3242E-01  1.9355E+00  2.4984E+00
 GRADIENT:  -9.8529E+00  1.0641E+01 -1.6824E+00  9.0385E+00 -2.3300E+01  5.4588E+00  2.0113E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -529.324301490965        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      221
 NPARAMETR:  1.2743E+00  1.8670E+00  2.8946E+01  9.6990E-01  8.8587E-01  8.4168E+00  1.0435E+01
 PARAMETER:  3.4239E-01  7.2433E-01  3.4654E+00  6.9442E-02 -2.1188E-02  2.2302E+00  2.4452E+00
 GRADIENT:   1.8096E+00  9.4399E+00 -1.4664E+00  5.1585E+01 -2.0792E+01  9.2318E+00  2.5359E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -535.482529215793        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:      330             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2741E+00  1.8665E+00  3.1163E+01  7.8882E-01  8.8499E-01  7.7517E+00  1.0427E+01
 PARAMETER:  3.4223E-01  7.2407E-01  3.5392E+00 -1.3721E-01 -2.2182E-02  2.1479E+00  2.4444E+00
 GRADIENT:   5.5583E+00  1.4046E+01 -7.6739E-01  2.3371E+00 -9.0999E+00  6.5870E+00  1.5701E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -537.239642721013        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      387
 NPARAMETR:  1.2740E+00  1.8656E+00  2.4983E+01  7.7833E-01  8.8500E-01  5.0871E+00  1.0302E+01
 PARAMETER:  3.4212E-01  7.2359E-01  3.3182E+00 -1.5061E-01 -2.2167E-02  1.7267E+00  2.4323E+00
 GRADIENT:   6.8480E+00  1.5988E+01  9.2339E-02 -5.5430E-02 -1.1707E+01  8.5371E-02  1.0672E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -537.527913970746        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  1.2735E+00  1.8624E+00  6.0488E+00  7.9999E-01  8.8504E-01  2.6530E+00  9.8987E+00
 PARAMETER:  3.4177E-01  7.2185E-01  1.8999E+00 -1.2316E-01 -2.2117E-02  1.0757E+00  2.3924E+00
 GRADIENT:   5.3357E+00  1.6462E+01 -3.8002E-02  2.1832E+00 -1.5140E+01  3.4340E-01 -1.3014E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -537.704663601510        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      497
 NPARAMETR:  1.2736E+00  1.8634E+00  5.7209E+00  7.7979E-01  8.8504E-01  2.5791E+00  1.0132E+01
 PARAMETER:  3.4187E-01  7.2238E-01  1.8441E+00 -1.4873E-01 -2.2126E-02  1.0474E+00  2.4157E+00
 GRADIENT:   5.2363E+00  1.5073E+01 -3.9153E-01 -1.0392E+00 -1.2292E+01  4.8115E-01  1.1905E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -537.771523002815        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.2733E+00  1.8611E+00  4.5207E+00  7.8568E-01  8.8509E-01  9.0019E-01  1.0084E+01
 PARAMETER:  3.4160E-01  7.2117E-01  1.6087E+00 -1.4121E-01 -2.2071E-02 -5.1461E-03  2.4109E+00
 GRADIENT:   4.0535E+00  1.4494E+01 -1.9036E-01 -4.4294E-02 -1.2972E+01  2.9245E-02 -5.5763E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -538.182857817136        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  1.2578E+00  1.7591E+00  4.8318E+00  7.6985E-01  8.8745E-01  1.0000E-02  1.0075E+01
 PARAMETER:  3.2936E-01  6.6483E-01  1.6752E+00 -1.6156E-01 -1.9403E-02 -4.9417E+01  2.4101E+00
 GRADIENT:   3.2107E+00  2.1549E+00  4.5174E-02 -2.8332E+00 -8.0210E+00  0.0000E+00  2.2875E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -538.195586275749        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  1.2539E+00  1.7345E+00  4.8081E+00  7.7905E-01  8.8804E-01  1.0000E-02  1.0009E+01
 PARAMETER:  3.2630E-01  6.5071E-01  1.6703E+00 -1.4968E-01 -1.8738E-02 -6.1666E+01  2.4035E+00
 GRADIENT:   7.0745E-01 -2.3445E+00  1.1423E-02 -4.9222E-01 -8.4170E+00  0.0000E+00 -1.6917E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -543.639971404430        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      764
 NPARAMETR:  1.2320E+00  1.7714E+00  4.7142E+00  7.8813E-01  9.7244E-01  1.0000E-02  9.8944E+00
 PARAMETER:  3.0868E-01  6.7179E-01  1.6506E+00 -1.3809E-01  7.2049E-02 -5.3622E+01  2.3920E+00
 GRADIENT:   1.6613E+00  5.4029E-01  6.5325E-03  1.8807E+00  9.0972E-01  0.0000E+00  1.7246E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -543.647763556093        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      859
 NPARAMETR:  1.2317E+00  1.7796E+00  4.6900E+00  7.8433E-01  9.6895E-01  1.0000E-02  9.9059E+00
 PARAMETER:  3.0836E-01  6.7638E-01  1.6454E+00 -1.4293E-01  6.8454E-02 -5.3622E+01  2.3931E+00
 GRADIENT:  -3.9751E-02 -1.9331E-02  3.7546E-04  2.3951E-02 -1.0598E-02  0.0000E+00  9.3406E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      859
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2887E-02 -6.2832E-02 -2.3350E-05
 SE:             8.9116E-02  8.8292E-02  2.3679E-05
 N:                     100         100         100

 P VAL.:         8.8502E-01  4.7669E-01  3.2409E-01

 ETASHRINKSD(%)  5.5903E+00  6.4628E+00  9.9975E+01
 ETASHRINKVR(%)  1.0868E+01  1.2508E+01  1.0000E+02
 EBVSHRINKSD(%)  5.6510E+00  5.3752E+00  9.9975E+01
 EBVSHRINKVR(%)  1.0983E+01  1.0462E+01  1.0000E+02
 RELATIVEINF(%)  4.8278E+01  7.1137E+01  3.0388E-06
 EPSSHRINKSD(%)  1.6086E+01
 EPSSHRINKVR(%)  2.9585E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -543.64776355609320     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       191.50306300764498     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     7.96
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
 





 #OBJV:********************************************     -543.648       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.23E+00  1.78E+00  4.69E+00  7.84E-01  9.69E-01  1.00E-02  9.91E+00
 


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
+        9.04E+01
 
 TH 2
+        3.07E+00  3.26E+00
 
 TH 3
+        3.05E-01  2.22E-02  1.08E-03
 
 TH 4
+       -1.01E+01 -1.18E+01 -7.30E-02  5.23E+02
 
 TH 5
+       -2.40E+01  2.18E+00 -7.39E-02 -1.37E+02  1.14E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -3.10E+00 -1.42E+00 -1.55E-02  8.42E-01  2.42E+00  0.00E+00  7.34E-01
 
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
+        1.06E+02
 
 TH 2
+       -3.57E+00  3.27E+01
 
 TH 3
+        3.36E-01  9.66E-02  4.08E-02
 
 TH 4
+        6.52E+00 -7.27E+00 -5.59E-03  2.89E+02
 
 TH 5
+       -3.04E+00 -3.85E+00 -2.79E-02 -4.10E+01  1.58E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -4.38E+00 -2.05E+00 -2.54E-02  1.42E+00  2.99E+00  0.00E+00  3.84E+00
 
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
+        1.07E+02
 
 TH 2
+        4.23E+01  3.37E+01
 
 TH 3
+        4.71E-01  2.11E-01  2.84E-03
 
 TH 4
+        2.31E+01  1.74E+01  2.06E-01  1.82E+02
 
 TH 5
+        2.14E+01  6.17E+01  1.17E-02  4.93E+01  2.74E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -5.01E+00 -3.07E+00 -4.77E-02  9.99E+00  1.56E+01  0.00E+00  2.58E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       10.081
Stop Time:
Sat Sep 18 15:57:21 CDT 2021
