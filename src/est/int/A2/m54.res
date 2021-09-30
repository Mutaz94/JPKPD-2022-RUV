Tue Sep 28 23:01:25 CDT 2021
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
$DATA ../../../../data/int/A2/dat54.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       28 SEP 2021
Days until program expires : 201
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2103.56911045000        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8615E+02  2.1726E+02  1.4333E+02  9.7591E+01  1.6582E+02  2.9383E+01 -7.9359E+01 -2.4033E+02 -5.6186E+01 -1.1382E+01
            -3.1160E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3031.92661818446        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1457E+00  9.1407E-01  8.7643E-01  1.0284E+00  8.8501E-01  1.0493E+00  9.6936E-01  9.9799E-01  9.5726E-01  6.5773E-01
             2.2400E+00
 PARAMETER:  2.3601E-01  1.0150E-02 -3.1900E-02  1.2800E-01 -2.2153E-02  1.4816E-01  6.8878E-02  9.7985E-02  5.6322E-02 -3.1896E-01
             9.0647E-01
 GRADIENT:   3.9849E+02  2.5622E+01  3.3594E+00  2.5440E+01  2.5357E+00  9.8690E+00 -1.9623E+01  1.1423E+01 -5.3791E+00 -1.1006E-01
             1.1375E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3037.48877594363        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      178
 NPARAMETR:  1.1320E+00  6.3340E-01  5.4933E-01  1.1787E+00  5.8130E-01  9.1844E-01  1.4863E+00  8.9496E-01  1.0122E+00  3.7471E-01
             2.1480E+00
 PARAMETER:  2.2403E-01 -3.5666E-01 -4.9905E-01  2.6438E-01 -4.4249E-01  1.4927E-02  4.9629E-01 -1.0979E-02  1.1213E-01 -8.8161E-01
             8.6453E-01
 GRADIENT:   1.8963E+02  2.8570E+01 -1.7898E+01  9.0699E+01  1.6775E+01 -5.8704E+01  2.9780E+01  1.0633E+01  1.7374E+01  8.3214E+00
             7.3255E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3052.52159902058        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  1.0807E+00  6.3831E-01  5.7788E-01  1.1420E+00  6.0564E-01  9.9438E-01  1.3504E+00  4.6308E-01  9.5276E-01  4.5816E-01
             2.1531E+00
 PARAMETER:  1.7760E-01 -3.4894E-01 -4.4838E-01  2.3279E-01 -4.0148E-01  9.4366E-02  4.0037E-01 -6.6986E-01  5.1612E-02 -6.8054E-01
             8.6690E-01
 GRADIENT:   6.2710E+01  2.9342E+00 -1.0383E+01  3.3750E+01  1.5339E+01 -1.1028E+01  5.7643E+00  2.2657E+00  2.2577E+00  3.5079E+00
             4.2090E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3054.58777669482        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  1.0495E+00  6.2150E-01  5.6362E-01  1.1290E+00  5.9005E-01  1.0204E+00  1.3270E+00  3.0978E-01  9.4496E-01  4.6755E-01
             2.1165E+00
 PARAMETER:  1.4836E-01 -3.7562E-01 -4.7338E-01  2.2131E-01 -4.2755E-01  1.2020E-01  3.8289E-01 -1.0719E+00  4.3384E-02 -6.6025E-01
             8.4978E-01
 GRADIENT:   2.1030E-01  9.1438E-02 -8.6563E-01  1.7898E-01  1.0623E+00  6.7360E-01 -7.5251E-02  8.1938E-03 -7.1768E-01 -4.4192E-01
             7.8562E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3054.60262398712        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  1.0495E+00  6.1673E-01  5.5980E-01  1.1306E+00  5.8557E-01  1.0186E+00  1.3208E+00  2.8433E-01  9.4704E-01  4.8342E-01
             2.1152E+00
 PARAMETER:  1.4835E-01 -3.8332E-01 -4.8017E-01  2.2276E-01 -4.3517E-01  1.1844E-01  3.7827E-01 -1.1576E+00  4.5589E-02 -6.2687E-01
             8.4916E-01
 GRADIENT:   3.0262E-01 -1.6743E-01 -9.6761E-01 -1.1005E-01  6.9934E-01  2.2447E-02  2.8212E-01  5.4547E-02  9.8912E-03 -6.2611E-02
            -4.7211E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3054.65779745626        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  1.0475E+00  6.1830E-01  5.6460E-01  1.1320E+00  5.8780E-01  1.0170E+00  1.2877E+00  8.7210E-02  9.4736E-01  5.3831E-01
             2.1180E+00
 PARAMETER:  1.4643E-01 -3.8077E-01 -4.7163E-01  2.2395E-01 -4.3137E-01  1.1690E-01  3.5283E-01 -2.3394E+00  4.5926E-02 -5.1933E-01
             8.5047E-01
 GRADIENT:  -3.5472E+00  2.3521E-02 -1.5260E-01 -9.3760E-01 -2.2742E+00 -2.9743E-01 -2.5213E-01  2.7881E-02  7.4885E-02  5.9481E-01
            -1.4340E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3054.70995016447        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  1.0451E+00  6.3548E-01  5.8322E-01  1.1262E+00  6.0620E-01  1.0160E+00  1.2932E+00  1.0000E-02  9.4670E-01  5.3003E-01
             2.1200E+00
 PARAMETER:  1.4414E-01 -3.5338E-01 -4.3918E-01  2.1889E-01 -4.0055E-01  1.1589E-01  3.5711E-01 -1.3156E+01  4.5224E-02 -5.3482E-01
             8.5141E-01
 GRADIENT:  -8.4726E+00 -9.6129E-01  2.0474E+00 -2.4754E+00 -2.8618E+00 -5.8177E-01 -5.4596E-01  0.0000E+00 -9.9418E-02 -2.6316E-01
            -2.0572E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -3054.74995722027        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  1.0498E+00  6.5148E-01  5.9712E-01  1.1227E+00  6.2172E-01  1.0174E+00  1.2984E+00  1.0000E-02  9.4759E-01  5.2954E-01
             2.1241E+00
 PARAMETER:  1.4856E-01 -3.2850E-01 -4.1564E-01  2.1570E-01 -3.7526E-01  1.1723E-01  3.6112E-01 -2.2824E+01  4.6165E-02 -5.3574E-01
             8.5335E-01
 GRADIENT:   7.1610E-01 -1.2233E-01 -8.9710E-01  7.1496E-02  1.0009E+00  7.4773E-02  4.6740E-01  0.0000E+00  2.9481E-02  1.7951E-01
            -1.0652E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1185
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.4742E-04 -1.5235E-04 -1.6964E-04 -2.3802E-03 -2.7842E-03
 SE:             2.9622E-02  2.4868E-02  2.5789E-04  2.8012E-02  1.8856E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7448E-01  9.9511E-01  5.1068E-01  9.3229E-01  8.8261E-01

 ETASHRINKSD(%)  7.6312E-01  1.6690E+01  9.9136E+01  6.1552E+00  3.6831E+01
 ETASHRINKVR(%)  1.5204E+00  3.0594E+01  9.9993E+01  1.1932E+01  6.0097E+01
 EBVSHRINKSD(%)  9.9560E-01  1.5886E+01  9.9146E+01  6.2068E+00  3.8243E+01
 EBVSHRINKVR(%)  1.9813E+00  2.9248E+01  9.9993E+01  1.2028E+01  6.1860E+01
 RELATIVEINF(%)  9.7980E+01  1.4637E+01  1.0123E-03  6.3692E+01  2.8437E+00
 EPSSHRINKSD(%)  1.7423E+01
 EPSSHRINKVR(%)  3.1810E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3054.7499572202714     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1400.6605974518607     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3054.750       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  6.51E-01  5.97E-01  1.12E+00  6.22E-01  1.02E+00  1.30E+00  1.00E-02  9.48E-01  5.30E-01  2.12E+00
 


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
+        9.55E+02
 
 TH 2
+       -3.74E+00  1.10E+03
 
 TH 3
+        6.39E+00  3.43E+02  2.03E+03
 
 TH 4
+       -1.29E+01  1.29E+02 -2.16E+02  8.64E+02
 
 TH 5
+       -8.54E+00 -1.37E+03 -2.36E+03  3.44E+02  3.80E+03
 
 TH 6
+        2.67E+00 -1.09E+00  6.51E+00 -4.11E+00 -3.64E+00  1.84E+02
 
 TH 7
+        1.98E-01  2.43E+01 -1.95E+01  4.30E+00  5.41E-02 -2.03E-01  6.05E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.40E+00 -1.53E+01  3.56E+01  9.47E+00 -2.49E+01 -1.24E+00  3.19E+00  0.00E+00  1.72E+02
 
 TH10
+       -4.26E-01 -2.35E+00 -1.11E+02 -1.97E+01  1.04E+02  7.38E-01  4.33E+01  0.00E+00 -2.01E-01  1.03E+02
 
 TH11
+       -1.21E+01 -1.77E+01 -2.12E+01 -1.15E+01 -6.20E-01  2.33E+00  5.72E-01  0.00E+00  7.85E+00  1.80E+01  2.58E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       41.240
Stop Time:
Tue Sep 28 23:02:08 CDT 2021
