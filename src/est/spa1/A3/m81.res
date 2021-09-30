Thu Sep 30 00:35:08 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat81.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -46.8076295058501        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2397E+02  8.8944E+01  1.6679E+02  4.3774E+01  2.8680E+02  2.0634E+01 -8.5969E+01 -1.1664E+02 -4.1414E+01 -2.2095E+02
            -3.6068E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1388.52734093854        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1301E+00  9.5936E-01  8.3535E-01  1.2457E+00  7.7680E-01  1.0285E+00  9.9700E-01  9.3689E-01  9.2815E-01  1.1551E+00
             6.1800E+00
 PARAMETER:  2.2233E-01  5.8510E-02 -7.9904E-02  3.1969E-01 -1.5257E-01  1.2813E-01  9.6994E-02  3.4806E-02  2.5437E-02  2.4417E-01
             1.9213E+00
 GRADIENT:   7.0549E+01  2.0546E+01 -1.0874E+01  5.6987E+01 -2.3932E+01  8.5811E+00  9.4641E+00  7.3051E+00  2.3927E+01  3.2333E+01
             3.9461E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1476.06835603292        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0519E+00  5.4186E-01  1.7860E-01  1.2777E+00  2.5409E-01  1.1402E+00  3.1384E-01  7.2926E-02  1.6723E+00  9.7999E-01
             4.2033E+00
 PARAMETER:  1.5056E-01 -5.1274E-01 -1.6226E+00  3.4508E-01 -1.2701E+00  2.3121E-01 -1.0589E+00 -2.5183E+00  6.1420E-01  7.9791E-02
             1.5359E+00
 GRADIENT:   1.2885E+01  7.3183E+01 -2.1068E+01  1.5450E+02 -7.2515E+00  2.2248E+01  1.8083E+00  8.9602E-02  3.2959E+01  5.7776E+01
             2.6853E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1501.47671738509        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  9.9916E-01  4.5053E-01  1.2691E-01  1.2182E+00  2.0695E-01  1.1998E+00  1.1102E-01  2.0629E-02  2.5156E+00  6.1171E-01
             3.4652E+00
 PARAMETER:  9.9162E-02 -6.9733E-01 -1.9643E+00  2.9737E-01 -1.4753E+00  2.8212E-01 -2.0981E+00 -3.7810E+00  1.0225E+00 -3.9150E-01
             1.3428E+00
 GRADIENT:  -4.9221E+01  1.0302E+02 -3.5965E+01  7.7788E+01 -5.5829E+01  2.5072E+01  7.2507E-02 -4.6151E-03  7.7404E+01 -1.4358E+01
             1.3179E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1505.97318583079        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.9972E-01  4.4364E-01  1.2694E-01  1.2141E+00  1.9981E-01  1.2011E+00  1.0000E-02  4.5931E-01  2.4413E+00  6.1318E-01
             3.3839E+00
 PARAMETER:  9.9715E-02 -7.1274E-01 -1.9640E+00  2.9398E-01 -1.5104E+00  2.8320E-01 -2.8149E+01 -6.7803E-01  9.9251E-01 -3.8910E-01
             1.3190E+00
 GRADIENT:  -2.5243E+01  1.5533E+02  3.0775E+01  9.7117E+01 -2.4666E+01  3.4582E+01  0.0000E+00 -1.3174E+00  9.3109E+01 -1.2530E+01
             1.3274E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1506.40850065689        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:      482
 NPARAMETR:  9.9972E-01  4.4344E-01  1.2696E-01  1.2140E+00  1.9956E-01  1.2011E+00  1.0000E-02  5.9675E-01  2.4391E+00  6.1321E-01
             3.3809E+00
 PARAMETER:  9.9725E-02 -7.1319E-01 -1.9639E+00  2.9390E-01 -1.5116E+00  2.8321E-01 -2.8963E+01 -4.1626E-01  9.9162E-01 -3.8906E-01
             1.3181E+00
 GRADIENT:  -5.0641E+01  1.4691E+02  8.9384E+00  8.3152E+01 -1.6381E+02  2.6481E+01  0.0000E+00 -1.4480E+00  6.8354E+01 -9.0498E+00
             1.2549E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1578.38491066329        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.0221E+00  3.6699E-01  1.3675E-01  9.2616E-01  2.0418E-01  1.0581E+00  1.0000E-02  1.0715E+00  1.4344E+00  6.3562E-01
             2.5804E+00
 PARAMETER:  1.2189E-01 -9.0241E-01 -1.8896E+00  2.3292E-02 -1.4888E+00  1.5649E-01 -2.8963E+01  1.6901E-01  4.6078E-01 -3.5315E-01
             1.0479E+00
 GRADIENT:   4.1423E+01  5.2869E+01  6.4858E+00 -2.5414E+01  2.0162E+02 -7.4213E-01  0.0000E+00  2.0682E+00 -6.5153E+00  3.6244E+00
            -1.4499E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1583.01787905698        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      760
 NPARAMETR:  1.0321E+00  3.3743E-01  1.4605E-01  9.7966E-01  2.0355E-01  1.0854E+00  1.0000E-02  1.1523E+00  1.5221E+00  6.1047E-01
             2.6371E+00
 PARAMETER:  1.3164E-01 -9.8640E-01 -1.8238E+00  7.9453E-02 -1.4918E+00  1.8196E-01 -2.8963E+01  2.4172E-01  5.2010E-01 -3.9353E-01
             1.0697E+00
 GRADIENT:   4.5375E+00  9.5302E+00 -7.2541E+00  4.6182E-01 -1.2329E+01  2.2178E+00  0.0000E+00  6.5769E+00  7.5395E+00  5.5928E+00
             1.6549E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1584.31974098173        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      935
 NPARAMETR:  1.0280E+00  3.4858E-01  1.6341E-01  1.0144E+00  2.1779E-01  1.0838E+00  1.0000E-02  1.1994E+00  1.4488E+00  5.0349E-01
             2.6361E+00
 PARAMETER:  1.2758E-01 -9.5389E-01 -1.7115E+00  1.1430E-01 -1.4242E+00  1.8047E-01 -2.8963E+01  2.8179E-01  4.7073E-01 -5.8619E-01
             1.0693E+00
 GRADIENT:  -3.4190E+00  1.8388E+00 -3.4119E+00  1.4126E-01 -2.5233E+00  1.4233E+00  0.0000E+00  2.0135E+00  2.9289E+00 -1.8626E+00
            -2.4599E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1584.43932082101        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1062
 NPARAMETR:  1.0300E+00  3.5402E-01  1.6807E-01  1.0238E+00  2.2169E-01  1.0788E+00  1.0000E-02  1.2025E+00  1.4202E+00  5.1548E-01
             2.6342E+00
 PARAMETER:  1.2958E-01 -9.3841E-01 -1.6834E+00  1.2356E-01 -1.4065E+00  1.7586E-01 -2.8963E+01  2.8440E-01  4.5079E-01 -5.6266E-01
             1.0686E+00
 GRADIENT:   1.4809E-01 -3.9028E-01  3.5315E-02 -3.0018E-01  6.1092E-01  4.7003E-03  0.0000E+00 -3.4506E-01 -2.1059E-01 -5.5432E-02
             5.3242E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1062
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8119E-03 -2.2709E-04  1.7167E-02 -5.6113E-03  1.9874E-02
 SE:             2.9284E-02  1.3707E-04  1.9959E-02  2.6696E-02  1.7716E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5066E-01  9.7585E-02  3.8973E-01  8.3352E-01  2.6193E-01

 ETASHRINKSD(%)  1.8944E+00  9.9541E+01  3.3135E+01  1.0566E+01  4.0649E+01
 ETASHRINKVR(%)  3.7528E+00  9.9998E+01  5.5290E+01  2.0015E+01  6.4775E+01
 EBVSHRINKSD(%)  1.9342E+00  9.9551E+01  3.1627E+01  7.6471E+00  4.1284E+01
 EBVSHRINKVR(%)  3.8309E+00  9.9998E+01  5.3251E+01  1.4709E+01  6.5524E+01
 RELATIVEINF(%)  9.6098E+01  4.3122E-04  6.6833E+00  4.8485E+01  2.3511E+00
 EPSSHRINKSD(%)  2.8682E+01
 EPSSHRINKVR(%)  4.9137E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1584.4393208210056     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -665.50078761633290     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1584.439       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  3.54E-01  1.68E-01  1.02E+00  2.22E-01  1.08E+00  1.00E-02  1.20E+00  1.42E+00  5.15E-01  2.63E+00
 


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
+        8.69E+02
 
 TH 2
+       -8.05E+01  3.52E+03
 
 TH 3
+        1.42E+02 -1.54E+03  1.85E+04
 
 TH 4
+       -2.42E+01  5.01E+02 -3.22E+03  5.62E+02
 
 TH 5
+        1.93E+02 -8.70E+03 -1.10E+04 -1.19E+03  3.40E+04
 
 TH 6
+        3.58E+00 -3.99E+01  1.09E+02 -1.36E+01  5.16E+01  1.58E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -2.46E+01  5.96E+02 -3.27E+03  1.67E+02 -8.88E+02 -1.10E+01  0.00E+00  2.85E+02
 
 TH 9
+       -2.32E-01  1.29E+02 -1.05E+03  4.87E+01 -5.44E+01 -3.54E+00  0.00E+00  8.27E+01  8.96E+01
 
 TH10
+        9.10E-01 -1.20E+02  4.27E+02 -9.16E+00  3.59E+02  3.66E+00  0.00E+00  1.45E+01  1.04E+01  9.65E+01
 
 TH11
+       -1.32E+01 -2.98E+02  4.66E+02 -6.79E+01  4.08E+02  2.59E+00  0.00E+00 -1.63E+02 -1.06E+01  2.63E+01  1.09E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.956
Stop Time:
Thu Sep 30 00:35:36 CDT 2021
