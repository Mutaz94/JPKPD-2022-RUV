Wed Sep 29 10:41:00 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/28/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/int/All/dat63.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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
$OMEGA (0.09 FIX)x3
$SIGMA  1  FIX;        [P]
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26858.1947490097        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   3.2590E+02  6.6871E+00 -1.5785E+03 -1.6282E+03 -2.3606E+03 -1.7060E+03 -5.4582E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1226.20408056007        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:       61
 NPARAMETR:  1.2203E+00  2.2146E+00  6.1450E+00  4.3174E+00  5.0008E+00  1.8326E+00  9.9291E+00
 PARAMETER:  2.9907E-01  8.9506E-01  1.9156E+00  1.5627E+00  1.7096E+00  7.0576E-01  2.3955E+00
 GRADIENT:   5.5922E+00  5.6097E+01  5.8601E+01  1.5399E+02  1.8636E+02  2.6380E+01  3.7906E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1265.40558704897        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:      116
 NPARAMETR:  1.5126E+00  1.9026E+00  2.3498E+01  2.1757E+00  3.5029E+00  6.2991E+00  9.1101E+00
 PARAMETER:  5.1382E-01  7.4323E-01  3.2569E+00  8.7733E-01  1.3536E+00  1.9404E+00  2.3094E+00
 GRADIENT:   1.1260E+02  4.6297E+01  1.0642E+00 -4.4795E+01  6.6398E+01 -1.7833E-03  1.5873E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1275.90567711046        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      168
 NPARAMETR:  1.3665E+00  1.7162E+00  2.3690E+01  2.2550E+00  3.2132E+00  7.8766E+00  8.7710E+00
 PARAMETER:  4.1227E-01  6.4013E-01  3.2650E+00  9.1316E-01  1.2673E+00  2.1639E+00  2.2714E+00
 GRADIENT:   6.6413E+01  2.5911E+01  7.0627E-01 -2.1923E+01  2.4551E+01  8.5986E-02  1.1189E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1278.21387642451        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      218
 NPARAMETR:  1.2232E+00  1.5840E+00  2.4806E+01  2.3011E+00  3.1033E+00  8.7171E+00  8.3864E+00
 PARAMETER:  3.0147E-01  5.5995E-01  3.3111E+00  9.3337E-01  1.2325E+00  2.2653E+00  2.2266E+00
 GRADIENT:   2.3037E+01  8.6124E+00  4.0433E-01 -1.0085E+01  5.1650E+00  1.2144E-01  4.1734E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1290.83163785659        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  1.1662E+00  1.6900E+00  2.2812E+01  2.3799E+00  3.5279E+00  7.8280E+00  8.4908E+00
 PARAMETER:  2.5378E-01  6.2471E-01  3.2273E+00  9.6707E-01  1.3607E+00  2.1577E+00  2.2390E+00
 GRADIENT:   1.3409E+01  2.6570E+01  3.9906E-01  3.1787E+01  7.4236E+01  2.3318E-01  6.2356E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1291.07838170699        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.1676E+00  1.7501E+00  2.2799E+01  2.3773E+00  3.5493E+00  7.8270E+00  8.4176E+00
 PARAMETER:  2.5496E-01  6.5965E-01  3.2267E+00  9.6595E-01  1.3668E+00  2.1576E+00  2.2303E+00
 GRADIENT:   4.1206E-01  8.2238E-01  3.3337E-01  9.9473E+00 -1.7509E+00  2.3878E-01  2.9138E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1291.29374991469        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      560
 NPARAMETR:  1.1681E+00  1.7573E+00  2.2023E+01  2.2971E+00  3.5803E+00  6.2811E+00  8.4501E+00
 PARAMETER:  2.5536E-01  6.6377E-01  3.1921E+00  9.3165E-01  1.3755E+00  1.9375E+00  2.2342E+00
 GRADIENT:   1.1138E+00  1.1906E+00  1.0348E+00  3.5248E-01  1.8080E+00  7.4972E-02  8.4304E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1291.59165574930        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:      688
 NPARAMETR:  1.1613E+00  1.7190E+00  1.7814E+01  2.3741E+00  3.5932E+00  5.8822E-01  8.2847E+00
 PARAMETER:  2.4953E-01  6.4176E-01  2.9800E+00  9.6461E-01  1.3790E+00 -4.3065E-01  2.2144E+00
 GRADIENT:  -2.9073E-01 -1.3939E+00  4.6454E+00  8.5264E+00  2.6898E+00  1.9088E-03 -2.7494E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1380.03431684206        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:      770
 NPARAMETR:  1.1725E+00  1.7351E+00  2.7010E+00  2.2035E+00  3.1021E+00  1.0000E-02  6.6732E+00
 PARAMETER:  2.5916E-01  6.5107E-01  1.0936E+00  8.9005E-01  1.2321E+00 -1.0754E+01  1.9981E+00
 GRADIENT:   1.8484E+01  8.5187E+01  4.7378E+01 -6.2423E+01  4.8733E+01  0.0000E+00 -2.8802E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1410.06607501650        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:      861
 NPARAMETR:  1.1717E+00  1.7278E+00  2.5670E+00  2.2593E+00  3.1610E+00  1.0000E-02  7.7217E+00
 PARAMETER:  2.5848E-01  6.4684E-01  1.0427E+00  9.1505E-01  1.2509E+00 -9.9533E+00  2.1440E+00
 GRADIENT:  -7.7944E+00  2.0049E+01  8.8478E-01 -4.7121E+01 -4.4612E+01  0.0000E+00  9.7194E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1415.21646267940        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  1.1927E+00  1.4943E+00  2.5453E+00  2.5923E+00  3.4325E+00  1.0000E-02  7.7733E+00
 PARAMETER:  2.7621E-01  5.0162E-01  1.0343E+00  1.0525E+00  1.3333E+00 -9.8537E+00  2.1507E+00
 GRADIENT:   1.2282E+01  2.5535E+01  1.0052E+01  1.9658E+01  8.1889E+01  0.0000E+00  9.5858E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1417.12646091382        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:     1085
 NPARAMETR:  1.2229E+00  1.5524E+00  2.5257E+00  2.7393E+00  3.4941E+00  1.0000E-02  7.5654E+00
 PARAMETER:  3.0122E-01  5.3981E-01  1.0265E+00  1.1077E+00  1.3511E+00 -9.8537E+00  2.1236E+00
 GRADIENT:   1.8192E-01 -3.6111E-01 -2.1665E-01 -4.3290E-01  1.4313E-01  0.0000E+00  2.2716E+00

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1417.12646091382        NO. OF FUNC. EVALS.:  14
 CUMULATIVE NO. OF FUNC. EVALS.:     1099
 NPARAMETR:  1.2229E+00  1.5524E+00  2.5257E+00  2.7393E+00  3.4941E+00  1.0000E-02  7.5654E+00
 PARAMETER:  3.0122E-01  5.3981E-01  1.0265E+00  1.1077E+00  1.3511E+00 -9.8537E+00  2.1236E+00
 GRADIENT:   1.8192E-01 -3.6111E-01 -2.1665E-01 -4.3290E-01  1.4313E-01  0.0000E+00  2.2716E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1099
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.6684E-03 -4.0053E-03 -1.0203E-04
 SE:             2.9532E-02  3.0009E-02  1.5057E-04
 N:                     100         100         100

 P VAL.:         7.9512E-01  8.9382E-01  4.9802E-01

 ETASHRINKSD(%)  1.0644E+00  1.0000E-10  9.9496E+01
 ETASHRINKVR(%)  2.1174E+00  1.0000E-10  9.9997E+01
 EBVSHRINKSD(%)  2.9084E+00  8.7419E-01  9.9458E+01
 EBVSHRINKVR(%)  5.7323E+00  1.7407E+00  9.9997E+01
 RELATIVEINF(%)  9.4225E+01  9.7483E+01  2.9115E-03
 EPSSHRINKSD(%)  8.9179E+00
 EPSSHRINKVR(%)  1.7040E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1417.1264609138229     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       236.96289885458782     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:    20.82
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     3.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1417.126       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.22E+00  1.55E+00  2.53E+00  2.74E+00  3.49E+00  1.00E-02  7.57E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 


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
+        1.10E+02
 
 TH 2
+        1.13E+00  3.80E+00
 
 TH 3
+        7.04E-01 -7.09E+00  4.65E+01
 
 TH 4
+       -2.06E+01  3.56E+00 -8.46E-01  8.93E+00
 
 TH 5
+       -9.50E+00  1.87E+00 -1.48E+00  4.22E+00  2.02E+00
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -1.40E+00 -4.87E-01 -1.07E-01 -1.45E-01 -5.74E-02  0.00E+00  7.38E-01
 
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
+        9.07E+01
 
 TH 2
+       -3.64E-01  3.70E+01
 
 TH 3
+       -9.06E-01 -3.86E+00  4.57E+01
 
 TH 4
+       -2.74E-01  9.42E-02 -1.99E-01  2.37E+01
 
 TH 5
+        6.05E-01  1.81E-01 -3.91E-01 -4.26E-01  1.58E+01
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -3.53E+00 -2.25E+00 -3.22E-01  1.14E+00  8.07E-01  0.00E+00  1.85E+01
 
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
+        9.49E+01
 
 TH 2
+        4.62E+01  3.79E+01
 
 TH 3
+        3.07E+00  6.45E-01  4.55E+01
 
 TH 4
+        2.51E+01  1.37E+01  3.40E+00  4.48E+01
 
 TH 5
+        1.20E+01  1.85E+01  1.63E+00  1.87E+01  2.69E+01
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -3.84E+01 -2.52E+01 -1.86E+00  5.37E+01  4.58E+01  0.00E+00  5.00E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.776
Stop Time:
Wed Sep 29 10:41:27 CDT 2021
