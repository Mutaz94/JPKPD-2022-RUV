Sat Sep 18 15:57:22 CDT 2021
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
$DATA ../../../../data/spa/All/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   11134.2571319982        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   5.3037E+01  1.4144E+01 -1.2721E+02 -7.6400E+01 -1.9208E+02 -6.8976E+01 -2.4023E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -510.042962992229        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  1.6391E+00  1.8788E+00  4.7938E+00  5.6741E-01  7.2131E-01  3.9032E-01  1.5154E+01
 PARAMETER:  5.9413E-01  7.3066E-01  1.6673E+00 -4.6668E-01 -2.2668E-01 -8.4079E-01  2.8182E+00
 GRADIENT:   4.7776E+01  3.0654E+01  7.7612E-01 -4.1293E+01  3.4513E+01 -1.1385E-02  1.1140E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -536.607298943940        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:      115
 NPARAMETR:  1.5180E+00  1.3776E+00  6.3632E+00  5.9480E-01  4.9142E-01  2.1736E-01  1.3459E+01
 PARAMETER:  5.1741E-01  4.2036E-01  1.9505E+00 -4.1953E-01 -6.1046E-01 -1.4262E+00  2.6997E+00
 GRADIENT:   7.4671E+01 -3.1705E+01  6.4235E-02  1.7822E+01  1.6953E+01 -4.1588E-03  5.7405E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -561.075648706681        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.1412E+00  1.0892E+00  7.8862E+00  3.3959E-01  2.3172E-01  1.0916E-01  1.3157E+01
 PARAMETER:  2.3206E-01  1.8548E-01  2.1651E+00 -9.8003E-01 -1.3622E+00 -2.1149E+00  2.6770E+00
 GRADIENT:  -5.6683E+00 -1.3958E+01 -1.3662E+00 -3.2108E+01  2.4057E+01 -1.3811E-03  1.3098E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -572.113369982187        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      215
 NPARAMETR:  1.0876E+00  1.0345E+00  7.0109E+01  3.4849E-01  5.3115E-02  1.1751E-01  1.3151E+01
 PARAMETER:  1.8393E-01  1.3394E-01  4.3501E+00 -9.5415E-01 -2.8353E+00 -2.0413E+00  2.6765E+00
 GRADIENT:  -1.3789E+01 -9.0220E-01 -1.5386E-01 -5.9563E+00  8.5897E-01 -2.0157E-05  9.7396E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -572.790724891473        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      265
 NPARAMETR:  1.1012E+00  1.0418E+00  1.1499E+03  3.5510E-01  1.0000E-02  1.9288E-01  1.3266E+01
 PARAMETER:  1.9642E-01  1.4091E-01  7.1474E+00 -9.3536E-01 -4.6038E+00 -1.5457E+00  2.6852E+00
 GRADIENT:   1.4333E-01 -2.1934E-01 -8.2817E-03  1.1799E-01  0.0000E+00 -2.0057E-07  1.5081E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -572.804539311905        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      346
 NPARAMETR:  1.1032E+00  1.0438E+00  7.6922E+04  3.5546E-01  1.0000E-02  4.2475E-01  1.3297E+01
 PARAMETER:  1.9822E-01  1.4291E-01  1.1351E+01 -9.3435E-01 -7.2412E+00 -7.5624E-01  2.6875E+00
 GRADIENT:  -5.7773E-01  1.7566E-01 -1.2271E-04 -3.5852E-01  0.0000E+00 -2.2550E-10 -8.5568E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -572.805715918473        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  1.1044E+00  1.0448E+00  2.3761E+06  3.5591E-01  1.0000E-02  8.0874E-01  1.3315E+01
 PARAMETER:  1.9933E-01  1.4384E-01  1.4781E+01 -9.3309E-01 -9.3930E+00 -1.1228E-01  2.6889E+00
 GRADIENT:  -2.5661E-02  7.7553E-03 -3.9517E-06 -7.9955E-03  0.0000E+00  0.0000E+00  1.9049E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -572.805720387899        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      576
 NPARAMETR:  1.1044E+00  1.0448E+00  3.3170E+08  3.5592E-01  1.0000E-02  2.0577E+00  1.3315E+01
 PARAMETER:  1.9935E-01  1.4384E-01  1.9720E+01 -9.3306E-01 -1.2489E+01  8.2160E-01  2.6889E+00
 GRADIENT:   1.6289E-03 -1.2545E-03 -2.8304E-08 -4.2925E-04  0.0000E+00  0.0000E+00  1.7875E-04

0ITERATION NO.:   42    OBJECTIVE VALUE:  -572.805720454543        NO. OF FUNC. EVALS.:  47
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.1044E+00  1.0448E+00  5.3815E+08  3.5594E-01  1.0000E-02  2.2499E+00  1.3315E+01
 PARAMETER:  1.9935E-01  1.4384E-01  2.0187E+01 -9.3305E-01 -1.2782E+01  9.0998E-01  2.6889E+00
 GRADIENT:   1.0988E-04 -5.2174E-04 -5.4934E-07 -1.5648E-03  0.0000E+00 -4.2137E-06  9.9669E-04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      623
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.5684E-02 -3.8412E-03  1.9890E-10
 SE:             8.6178E-02  2.6023E-03  1.3245E-10
 N:                     100         100         100

 P VAL.:         6.7882E-01  1.3993E-01  1.3318E-01

 ETASHRINKSD(%)  8.7026E+00  9.7243E+01  1.0000E+02
 ETASHRINKVR(%)  1.6648E+01  9.9924E+01  1.0000E+02
 EBVSHRINKSD(%)  8.4806E+00  9.7166E+01  1.0000E+02
 EBVSHRINKVR(%)  1.6242E+01  9.9920E+01  1.0000E+02
 RELATIVEINF(%)  6.3581E+01  6.0998E-02  1.0000E-10
 EPSSHRINKSD(%)  3.8016E+00
 EPSSHRINKVR(%)  7.4588E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -572.80572045454301     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       162.34510610919517     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     4.16
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     1.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -572.806       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.10E+00  1.04E+00  5.29E+08  3.56E-01  1.00E-02  2.25E+00  1.33E+01
 


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
+        8.75E+02
 
 TH 2
+       -8.62E+02  1.21E+03
 
 TH 3
+       -1.44E-12  5.39E-12  4.59E-26
 
 TH 4
+        1.72E+02 -2.72E+02 -1.81E-12  8.57E+02
 
 TH 5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 6
+        1.08E-03 -3.90E-04  5.86E-18 -5.41E-04  0.00E+00  2.98E-09
 
 TH 7
+       -1.46E+00 -1.26E+01 -1.51E-13  4.97E+00  0.00E+00 -2.86E-05  5.44E-01
 
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
+        6.11E+02
 
 TH 2
+       -3.48E+02  8.04E+02
 
 TH 3
+        1.75E-12  5.63E-12 -1.28E-23
 
 TH 4
+        7.02E+01 -2.58E+02 -2.70E-12  1.09E+03
 
 TH 5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 6
+        1.21E-03  4.64E-04 -1.58E-14 -9.90E-04  0.00E+00 -2.13E-05
 
 TH 7
+       -1.07E+01 -1.72E+01 -3.40E-14  8.78E+00  0.00E+00 -5.30E-05  3.19E+00
 
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
+        6.03E+02
 
 TH 2
+        7.54E+01  7.25E+02
 
 TH 3
+        2.57E-16 -3.56E-15  1.70E-30
 
 TH 4
+       -5.36E+01 -3.42E+02  6.58E-16  1.43E+03
 
 TH 5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 6
+       -2.39E-06  3.22E-06 -6.30E-23  4.74E-07  0.00E+00  1.24E-12
 
 TH 7
+       -1.77E+01 -3.77E+01  1.14E-15  4.39E+01  0.00E+00 -9.79E-07  1.13E+01
 
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
 #CPUT: Total CPU Time in Seconds,        5.624
Stop Time:
Sat Sep 18 15:57:29 CDT 2021
