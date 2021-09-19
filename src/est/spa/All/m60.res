Sat Sep 18 15:56:08 CDT 2021
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
$DATA ../../../../data/spa/All/dat60.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   10517.9311258406        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   6.5208E+01  2.2199E+01 -1.2577E+02  8.8174E-01 -1.0115E+02 -7.2116E+01 -2.2914E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:   1897.99217707854        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:       61
 NPARAMETR:  4.7179E-01  1.1479E-01  3.8962E+00  1.5309E+00  1.5966E+00  1.0517E+00  2.1367E+00
 PARAMETER: -6.5122E-01 -2.0647E+00  1.4600E+00  5.2582E-01  5.6786E-01  1.5043E-01  8.5928E-01
 GRADIENT:   1.6327E+00 -8.1539E+01  6.7045E+01  1.3831E+02 -2.1115E+02 -6.0055E+01 -4.9271E+03

0ITERATION NO.:   10    OBJECTIVE VALUE:  -384.358431317379        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      113
 NPARAMETR:  1.0747E+00  9.3537E-01  2.2194E+00  1.2821E+00  1.4598E+00  3.2394E-01  5.4488E+00
 PARAMETER:  1.7207E-01  3.3191E-02  8.9726E-01  3.4850E-01  4.7829E-01 -1.0272E+00  1.7954E+00
 GRADIENT:   2.3557E+00 -3.3230E+01 -4.3892E+00  1.2044E+02  1.0457E+02  1.7647E+00 -5.6457E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -592.427175302029        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.2230E+00  1.5725E+00  2.7288E+00  6.8360E-01  8.0048E-01  4.8410E-01  9.3156E+00
 PARAMETER:  3.0130E-01  5.5268E-01  1.1038E+00 -2.8038E-01 -1.2255E-01 -6.2546E-01  2.3317E+00
 GRADIENT:   4.0298E+00  1.1842E+01 -3.7380E+00 -2.7287E-01 -4.5739E-01  9.4673E-01 -1.0550E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -593.033117761557        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      216
 NPARAMETR:  1.1998E+00  1.5087E+00  3.7915E+00  6.8042E-01  7.9848E-01  5.0820E-01  9.2713E+00
 PARAMETER:  2.8216E-01  5.1128E-01  1.4328E+00 -2.8504E-01 -1.2504E-01 -5.7688E-01  2.3269E+00
 GRADIENT:  -1.5566E-01 -2.4106E-01 -1.6061E-01 -6.4850E-02 -4.6330E-01  8.3260E-02 -5.1410E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -593.042172868873        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  1.1998E+00  1.5103E+00  3.9948E+00  6.8095E-01  7.9993E-01  5.0772E-01  9.2791E+00
 PARAMETER:  2.8215E-01  5.1231E-01  1.4850E+00 -2.8426E-01 -1.2323E-01 -5.7783E-01  2.3278E+00
 GRADIENT:  -1.6891E+00 -1.5942E+00 -7.1309E-03 -5.7488E-01 -1.5990E-01  4.6571E-02 -1.8259E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -593.068802878107        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:      384
 NPARAMETR:  1.2055E+00  1.5186E+00  3.9807E+00  6.8176E-01  7.9830E-01  3.2453E-01  9.2988E+00
 PARAMETER:  2.8692E-01  5.1776E-01  1.4815E+00 -2.8308E-01 -1.2528E-01 -1.0254E+00  2.3299E+00
 GRADIENT:   1.5269E+00  1.2402E+00  1.8080E-01  5.0751E-01 -6.2018E-01  1.9502E-02  6.3437E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -593.104342291843        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      438
 NPARAMETR:  1.2027E+00  1.5090E+00  3.2515E+00  6.8144E-01  8.0050E-01  4.9323E-02  9.2678E+00
 PARAMETER:  2.8453E-01  5.1144E-01  1.2791E+00 -2.8355E-01 -1.2251E-01 -2.9094E+00  2.3265E+00
 GRADIENT:  -4.3299E-02 -1.5505E-01  4.3467E-02  7.8037E-02  1.8095E-01  2.8238E-03 -8.0151E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -593.118888877729        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      558
 NPARAMETR:  1.2089E+00  1.5188E+00  3.2136E+00  6.8211E-01  8.0040E-01  7.5562E-02  9.3006E+00
 PARAMETER:  2.8975E-01  5.1795E-01  1.2674E+00 -2.8257E-01 -1.2264E-01 -2.4828E+00  2.3301E+00
 GRADIENT:  -2.9157E-01 -3.0284E-01 -1.1943E-02 -2.8299E-01  7.9607E-02  7.1586E-03 -2.8905E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -593.123014278064        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  1.2103E+00  1.5210E+00  3.1997E+00  6.8279E-01  8.0032E-01  1.0000E-02  9.3063E+00
 PARAMETER:  2.9088E-01  5.1936E-01  1.2631E+00 -2.8157E-01 -1.2275E-01 -4.6213E+00  2.3307E+00
 GRADIENT:   7.6936E-03  1.8468E-03 -1.0540E-04  2.4967E-02 -9.2443E-03  0.0000E+00  2.6983E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -593.123014278064        NO. OF FUNC. EVALS.:  14
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  1.2103E+00  1.5210E+00  3.1997E+00  6.8279E-01  8.0032E-01  1.0000E-02  9.3063E+00
 PARAMETER:  2.9088E-01  5.1936E-01  1.2631E+00 -2.8157E-01 -1.2275E-01 -4.6213E+00  2.3307E+00
 GRADIENT:   7.6936E-03  1.8468E-03 -1.0540E-04  2.4967E-02 -9.2443E-03  0.0000E+00  2.6983E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      687
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9536E-02 -7.2485E-02 -1.1876E-04
 SE:             9.0169E-02  8.4366E-02  4.8780E-05
 N:                     100         100         100

 P VAL.:         8.2847E-01  3.9025E-01  1.4908E-02

 ETASHRINKSD(%)  4.4748E+00  1.0622E+01  9.9948E+01
 ETASHRINKVR(%)  8.7494E+00  2.0116E+01  1.0000E+02
 EBVSHRINKSD(%)  5.3201E+00  6.9981E+00  9.9905E+01
 EBVSHRINKVR(%)  1.0357E+01  1.3506E+01  1.0000E+02
 RELATIVEINF(%)  8.6988E+01  8.5051E+01  8.6111E-05
 EPSSHRINKSD(%)  1.6130E+01
 EPSSHRINKVR(%)  2.9658E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -593.12301427806437     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       142.02781228567380     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     6.56
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     2.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -593.123       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.21E+00  1.52E+00  3.20E+00  6.83E-01  8.00E-01  1.00E-02  9.31E+00
 


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
+        1.22E+02
 
 TH 2
+       -1.16E+01  2.14E+01
 
 TH 3
+        6.07E-01 -1.04E-01  3.19E-03
 
 TH 4
+       -4.17E+01 -1.14E+01 -1.77E-01  4.99E+02
 
 TH 5
+       -1.47E+01 -9.33E+00 -1.25E-01 -2.12E+02  2.18E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -1.79E+00 -3.47E+00 -3.24E-03  4.62E+00  4.98E+00  0.00E+00  7.96E-01
 
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
+        1.35E+02
 
 TH 2
+       -8.07E+00  6.33E+01
 
 TH 3
+        6.71E-01 -3.07E-01  1.50E-01
 
 TH 4
+       -1.50E+01 -1.76E+01 -8.11E-02  3.49E+02
 
 TH 5
+       -1.81E+00 -4.55E+00 -1.15E-01 -8.08E+01  2.26E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -4.28E+00 -3.01E+00 -6.39E-02  6.47E+00  5.20E+00  0.00E+00  4.23E+00
 
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
+        1.49E+02
 
 TH 2
+        6.86E+01  6.56E+01
 
 TH 3
+        4.42E-01  2.84E-01  2.94E-02
 
 TH 4
+        1.32E+01  1.46E+01 -2.63E-01  2.88E+02
 
 TH 5
+        3.84E+01  8.83E+01  4.73E-02  5.33E+01  3.34E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -7.59E+00 -1.27E+01 -6.91E-02  2.04E+01  9.34E+00  0.00E+00  3.46E+01
 
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
 #CPUT: Total CPU Time in Seconds,        8.869
Stop Time:
Sat Sep 18 15:56:18 CDT 2021
