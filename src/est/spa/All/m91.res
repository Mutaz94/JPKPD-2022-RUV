Sat Sep 25 15:15:33 CDT 2021
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
$DATA ../../../../data/spa/All/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   19286.1548979326        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   8.8192E+02 -4.9798E+02  4.2146E+02 -3.5127E+03  3.6080E+03 -4.0139E+02 -3.9143E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -360.675335754706        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  1.4834E+00  1.5511E+00  3.0568E+00  4.7073E-01  2.8454E-01  4.0613E-01  1.5770E+01
 PARAMETER:  4.9436E-01  5.3898E-01  1.2174E+00 -6.5347E-01 -1.1569E+00 -8.0107E-01  2.8581E+00
 GRADIENT:   1.2117E+02  1.3567E+02 -2.1737E+00 -6.0098E+01 -2.3827E+01  2.0005E-01 -1.3745E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -390.980460930050        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      111
 NPARAMETR:  1.1395E+00  1.1637E+00  2.5426E+00  4.5537E-01  1.8579E-01  7.1687E-02  1.5794E+01
 PARAMETER:  2.3059E-01  2.5160E-01  1.0332E+00 -6.8665E-01 -1.5831E+00 -2.5354E+00  2.8596E+00
 GRADIENT:  -3.5466E+00  2.4584E+00 -1.0045E+00 -3.7857E+00  1.0572E+00  1.8813E-02 -1.7666E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -393.592550016461        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.1179E+00  1.1391E+00  8.5786E+01  4.4949E-01  1.6680E-01  2.5908E-01  1.5830E+01
 PARAMETER:  2.1148E-01  2.3025E-01  4.5519E+00 -6.9965E-01 -1.6910E+00 -1.2506E+00  2.8619E+00
 GRADIENT:   3.8607E+00 -3.2169E+00 -1.0251E-01 -3.3254E+00  8.7114E-01 -3.6707E-05 -6.1272E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -393.742156523189        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      214
 NPARAMETR:  1.1115E+00  1.1388E+00  1.5836E+03  4.5428E-01  1.5481E-01  7.4625E-01  1.5892E+01
 PARAMETER:  2.0575E-01  2.2996E-01  7.4674E+00 -6.8904E-01 -1.7656E+00 -1.9270E-01  2.8658E+00
 GRADIENT:  -4.1866E-01  6.3472E-01 -5.5956E-03  3.5131E-01 -1.5667E-01 -8.7848E-07  1.2725E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -393.744882565116        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      266
 NPARAMETR:  1.1118E+00  1.1381E+00  3.4953E+03  4.5372E-01  1.5753E-01  1.0495E+00  1.5866E+01
 PARAMETER:  2.0602E-01  2.2932E-01  8.2592E+00 -6.9027E-01 -1.7481E+00  1.4828E-01  2.8642E+00
 GRADIENT:   2.4741E-01 -3.9115E-01 -2.5347E-03 -1.4815E-01  5.6195E-02 -3.5651E-07  6.6718E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -393.744985363403        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      316
 NPARAMETR:  1.1113E+00  1.1380E+00  7.0642E+03  4.5397E-01  1.5751E-01  1.3958E+00  1.5856E+01
 PARAMETER:  2.0557E-01  2.2924E-01  8.9628E+00 -6.8972E-01 -1.7483E+00  4.3346E-01  2.8635E+00
 GRADIENT:   8.4978E-02 -1.4114E-01 -1.2531E-03 -5.7751E-02  2.4925E-02 -1.5448E-07  2.2356E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -393.745074689192        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  1.1111E+00  1.1379E+00  1.7117E+04  4.5408E-01  1.5743E-01  1.9922E+00  1.5852E+01
 PARAMETER:  2.0537E-01  2.2918E-01  9.8478E+00 -6.8949E-01 -1.7488E+00  7.8926E-01  2.8633E+00
 GRADIENT:   1.5091E-02 -2.9076E-02 -5.1681E-04 -1.4272E-02  7.5633E-03 -5.3440E-08  3.6847E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -393.746458801406        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:      424
 NPARAMETR:  1.1114E+00  1.1379E+00  9.2063E+12  4.5391E-01  1.5729E-01  6.2815E+03  1.5859E+01
 PARAMETER:  2.0562E-01  2.2919E-01  2.9951E+01 -6.8986E-01 -1.7497E+00  8.8454E+00  2.8637E+00
 GRADIENT:   1.1635E-01 -1.8034E-01  0.0000E+00 -5.6915E-02  1.8623E-02  0.0000E+00  3.3654E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -393.756229766175        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.1153E+00  1.1431E+00  9.6548E+17  4.5595E-01  1.6000E-01  6.3945E+05  1.5897E+01
 PARAMETER:  2.0916E-01  2.3377E-01  4.1511E+01 -6.8538E-01 -1.7326E+00  1.3468E+01  2.8661E+00
 GRADIENT:  -2.3821E-01  9.4432E-02  0.0000E+00 -1.2688E-01 -5.5040E-02  0.0000E+00 -9.6579E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -393.756928160302        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  1.1161E+00  1.1434E+00  1.0490E+18  4.5627E-01  1.6023E-01  6.5658E+05  1.5911E+01
 PARAMETER:  2.0985E-01  2.3399E-01  4.1594E+01 -6.8467E-01 -1.7311E+00  1.3495E+01  2.8670E+00
 GRADIENT:   3.5486E-02 -4.8126E-01  0.0000E+00  1.1155E-01 -1.2024E-02  0.0000E+00 -2.1506E-01

0ITERATION NO.:   53    OBJECTIVE VALUE:  -393.757280058362        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.1171E+00  1.1447E+00  1.0490E+18  4.5622E-01  1.6073E-01  6.5658E+05  1.5923E+01
 PARAMETER:  2.1073E-01  2.3513E-01  4.1594E+01 -6.8479E-01 -1.7280E+00  1.3495E+01  2.8678E+00
 GRADIENT:   1.1108E-01 -7.2120E-02  0.0000E+00 -3.9378E-02 -4.0261E-03  0.0000E+00  1.7480E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      725
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4503E-02 -5.8776E-02  0.0000E+00
 SE:             8.5325E-02  3.6300E-02  0.0000E+00
 N:                     100         100         100

 P VAL.:         6.8594E-01  1.0541E-01  1.0000E+00

 ETASHRINKSD(%)  9.6066E+00  6.1544E+01  1.0000E+02
 ETASHRINKVR(%)  1.8290E+01  8.5212E+01  1.0000E+02
 EBVSHRINKSD(%)  8.9210E+00  6.2865E+01  1.0000E+02
 EBVSHRINKVR(%)  1.7046E+01  8.6210E+01  1.0000E+02
 RELATIVEINF(%)  7.0747E+01  1.1761E+01  1.0000E-10
 EPSSHRINKSD(%)  3.1169E+00
 EPSSHRINKVR(%)  6.1366E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -393.75728005836237     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       341.39354650537581     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           200
  
 #TERE:
 Elapsed estimation  time in seconds:     5.14
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     1.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -393.757       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.12E+00  1.14E+00  1.05E+18  4.56E-01  1.61E-01  6.57E+05  1.59E+01
 


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
+        5.45E+02
 
 TH 2
+       -4.50E+02  5.64E+02
 
 TH 3
+        0.00E+00  0.00E+00  0.00E+00
 
 TH 4
+        2.08E+02 -1.74E+02  0.00E+00  7.87E+02
 
 TH 5
+        1.50E+01 -1.81E+02  0.00E+00 -2.69E+02  3.73E+02
 
 TH 6
+        8.53E-58 -1.35E-57  0.00E+00  2.21E-57 -5.11E-57  2.80E-99
 
 TH 7
+       -3.46E+00 -5.36E+00  0.00E+00 -1.67E+00  9.18E+00 -7.49E-43  4.04E-01
 
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
+        3.56E+02
 
 TH 2
+       -1.55E+02  4.04E+02
 
 TH 3
+        0.00E+00  0.00E+00  0.00E+00
 
 TH 4
+        2.91E+01 -1.57E+02  0.00E+00  6.32E+02
 
 TH 5
+       -4.21E+01 -1.69E+02  0.00E+00 -5.96E+01  3.08E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  1.38E-06
 
 TH 7
+       -7.64E+00 -1.01E+01  0.00E+00  6.15E+00  1.09E+01  0.00E+00  1.90E+00
 
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
+        3.56E+02
 
 TH 2
+        1.23E+02  4.66E+02
 
 TH 3
+        0.00E+00  0.00E+00  0.00E+00
 
 TH 4
+       -1.53E+02 -2.60E+02  0.00E+00  6.72E+02
 
 TH 5
+       -1.41E+02 -2.37E+02  0.00E+00  1.62E+02  3.31E+02
 
 TH 6
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       -1.47E+01 -2.90E+01  0.00E+00  2.93E+01  1.37E+01  0.00E+00  7.21E+00
 
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
 #CPUT: Total CPU Time in Seconds,        6.898
Stop Time:
Sat Sep 25 15:15:46 CDT 2021
