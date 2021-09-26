Sat Sep 25 14:10:53 CDT 2021
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
$DATA ../../../../data/spa/D/dat25.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   6549.38392476200        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.3172E+02 -6.5890E+01 -8.2619E+01 -1.6817E+02  1.9371E+02 -1.3577E+03 -4.8764E+02 -3.8078E+00 -7.8161E+02 -3.7598E+02
            -1.3519E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -714.192857168952        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.7688E+00  1.2076E+00  1.0127E+00  1.8018E+00  1.1737E+00  2.0350E+00  1.5300E+00  1.0279E+00  1.8175E+00  1.1537E+00
             1.3300E+01
 PARAMETER:  6.7028E-01  2.8867E-01  1.1260E-01  6.8880E-01  2.6013E-01  8.1051E-01  5.2527E-01  1.2751E-01  6.9747E-01  2.4295E-01
             2.6877E+00
 GRADIENT:   9.9264E+01  4.4144E+00 -1.3358E+01  8.6104E+00 -6.8476E+00  3.7519E+01  2.2567E+00  5.8021E+00  1.8877E+01  4.2988E+00
             1.9703E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -747.043965862638        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.5772E+00  1.0155E+00  5.8598E+00  2.2743E+00  6.1885E+00  2.0586E+00  4.3645E+00  1.4063E+00  2.1379E+00  1.1606E+01
             1.0713E+01
 PARAMETER:  5.5568E-01  1.1537E-01  1.8681E+00  9.2165E-01  1.9227E+00  8.2202E-01  1.5735E+00  4.4095E-01  8.5985E-01  2.5515E+00
             2.4715E+00
 GRADIENT:   7.5234E+01  1.7032E+01  7.6599E-01  4.3416E+01 -3.2904E+00  2.0592E+01  1.6543E+01 -6.2122E-02  3.6116E+01  1.3754E+01
             1.5066E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -795.427633734737        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1806E+00  5.1686E-01  4.7509E+00  2.0427E+00  2.7137E+00  1.7963E+00  4.3279E+00  5.8164E+00  1.5000E+00  6.6754E+00
             9.3534E+00
 PARAMETER:  2.6599E-01 -5.5998E-01  1.6583E+00  8.1428E-01  1.0983E+00  6.8574E-01  1.5651E+00  1.8607E+00  5.0545E-01  1.9984E+00
             2.3357E+00
 GRADIENT:  -3.6690E+01  1.6602E+01  1.5871E+01 -9.5171E-02 -5.5096E+01  9.8793E+00  5.5197E+00  3.6836E+00  1.1348E+01  4.4422E+01
             1.4429E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -853.298868101972        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.1655E+00  4.6134E-01  1.5440E+00  1.6918E+00  4.0517E+00  1.6923E+00  2.7698E+00  3.1092E+00  1.4196E+00  6.2923E+00
             6.6828E+00
 PARAMETER:  2.5311E-01 -6.7362E-01  5.3440E-01  6.2579E-01  1.4991E+00  6.2612E-01  1.1188E+00  1.2344E+00  4.5035E-01  1.9393E+00
             1.9995E+00
 GRADIENT:   2.0632E+01  2.5099E+01  1.1023E+01  4.1056E+00 -1.3943E+01 -2.7407E+01  3.6719E-01 -1.2909E+01 -8.1285E+00  4.0688E+00
            -3.0583E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -896.337195005533        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.6679E-01  8.3595E-02  6.6148E-01  1.4799E+00  7.5341E+00  1.6186E+00  9.9941E-01  3.6399E+00  7.0770E-01  9.2945E+00
             6.6387E+00
 PARAMETER:  6.6222E-02 -2.3818E+00 -3.1328E-01  4.9199E-01  2.1194E+00  5.8157E-01  9.9407E-02  1.3920E+00 -2.4574E-01  2.3294E+00
             1.9929E+00
 GRADIENT:  -2.4613E+01  2.2348E+00  9.9358E+00  2.1174E+01 -1.9402E+00  1.0413E+01  1.7217E-01 -6.7846E+00  1.2160E+00  1.2695E+01
            -1.0461E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -902.567088739894        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  9.4009E-01  6.0649E-02  4.1120E-01  1.3187E+00  7.5484E+00  1.5519E+00  7.2326E-01  3.2735E+00  5.8447E-01  8.5638E+00
             6.4778E+00
 PARAMETER:  3.8219E-02 -2.7026E+00 -7.8869E-01  3.7666E-01  2.1213E+00  5.3947E-01 -2.2399E-01  1.2859E+00 -4.3705E-01  2.2475E+00
             1.9684E+00
 GRADIENT:   5.9500E+00  3.0947E+00 -1.0184E+01  5.8639E+01 -1.6632E+01 -3.0625E+00  8.8244E-02 -8.5539E+00 -1.5643E+00  3.3827E+01
            -2.4880E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -905.803527461215        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  9.0952E-01  4.2452E-02  2.8307E-01  1.1827E+00  7.9555E+00  1.5235E+00  5.1481E-01  3.0065E+00  4.6221E-01  8.0604E+00
             6.3845E+00
 PARAMETER:  5.1670E-03 -3.0594E+00 -1.1620E+00  2.6779E-01  2.1739E+00  5.2101E-01 -5.6395E-01  1.2008E+00 -6.7174E-01  2.1870E+00
             1.9539E+00
 GRADIENT:   3.1841E+01  3.6299E+00 -2.8477E+01  3.8338E+01 -4.0452E+01 -4.8889E+00  1.3856E-02 -7.1420E-03  3.5183E+00  4.2542E+01
            -3.8906E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -906.208082873931        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      652
 NPARAMETR:  9.0919E-01  4.2222E-02  2.8305E-01  1.1818E+00  7.9740E+00  1.5223E+00  5.0992E-01  3.0024E+00  4.3903E-01  8.0353E+00
             6.3991E+00
 PARAMETER:  4.7938E-03 -3.0648E+00 -1.1621E+00  2.6708E-01  2.1762E+00  5.2020E-01 -5.7350E-01  1.1994E+00 -7.2320E-01  2.1838E+00
             1.9562E+00
 GRADIENT:   2.6420E+01  1.1066E+00 -1.9501E+01  3.7873E+01 -1.2956E+01 -6.3608E+00  1.3287E-02 -5.6186E+00  2.7574E+00  3.7275E+01
            -2.9102E+01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -906.208082873931        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  9.0919E-01  4.2227E-02  2.8303E-01  1.1818E+00  7.9734E+00  1.5222E+00  5.0963E-01  3.0025E+00  4.3904E-01  8.0359E+00
             6.3986E+00
 PARAMETER:  4.7938E-03 -3.0648E+00 -1.1621E+00  2.6708E-01  2.1762E+00  5.2020E-01 -5.7350E-01  1.1994E+00 -7.2320E-01  2.1838E+00
             1.9562E+00
 GRADIENT:  -2.0954E+04 -6.8329E+02  1.7841E+03  3.9648E+03  9.5159E+02  4.0259E+03  1.2587E-02 -8.8092E+02 -2.8978E+03 -9.2601E+02
             1.0414E+03
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         3.3         1.9         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      724
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.5828E-03 -1.0905E-03 -3.5012E-02 -2.7482E-02 -3.3218E-02
 SE:             2.9279E-02  2.0485E-04  2.0445E-02  9.7770E-03  7.7049E-03
 N:                     100         100         100         100         100

 P VAL.:         7.4345E-01  1.0211E-07  8.6801E-02  4.9411E-03  1.6248E-05

 ETASHRINKSD(%)  1.9126E+00  9.9314E+01  3.1508E+01  6.7246E+01  7.4187E+01
 ETASHRINKVR(%)  3.7886E+00  9.9995E+01  5.3088E+01  8.9272E+01  9.3337E+01
 EBVSHRINKSD(%)  3.8020E+00  9.9116E+01  2.1772E+01  6.6564E+01  5.9870E+01
 EBVSHRINKVR(%)  7.4594E+00  9.9992E+01  3.8803E+01  8.8820E+01  8.3896E+01
 RELATIVEINF(%)  5.2691E+01  4.8962E-03  1.4823E+01  1.8753E+00  1.2766E+01
 EPSSHRINKSD(%)  1.5663E+01
 EPSSHRINKVR(%)  2.8873E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -906.20808287393061     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -171.05725631019243     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.34
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -906.208       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.09E-01  4.22E-02  2.83E-01  1.18E+00  7.97E+00  1.52E+00  5.10E-01  3.00E+00  4.39E-01  8.04E+00  6.40E+00
 


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
+        6.35E+07
 
 TH 2
+        8.19E+03  3.13E+07
 
 TH 3
+       -3.32E+03  2.86E+03  4.85E+06
 
 TH 4
+        1.83E+07 -7.45E+03  2.62E+03  2.62E+02
 
 TH 5
+       -6.54E+01  1.84E+01 -2.27E+01  5.68E+01  1.75E+03
 
 TH 6
+       -1.36E+03 -2.87E+03  1.13E+03  2.10E+06  2.11E+01  8.37E+05
 
 TH 7
+        1.09E+00  1.40E+01 -1.07E+00 -4.40E+00  2.39E-02  5.17E-01 -1.69E+00
 
 TH 8
+        2.99E+02 -2.99E+02  1.91E+02 -2.72E+02 -1.45E+01 -1.05E+02 -3.73E-03  4.04E+04
 
 TH 9
+        3.37E+03 -1.58E+03  6.14E+02  6.17E+02  1.16E+01  2.52E+02 -1.66E+00 -5.20E+01  5.20E+06
 
 TH10
+        5.94E+01 -5.51E+01  4.36E+01 -6.32E+01 -9.14E+00 -2.00E+01  6.72E-02  4.16E+01 -1.08E+01  1.71E+03
 
 TH11
+       -9.31E+01  6.75E+01  1.27E+05 -1.33E+05 -5.65E-01  3.03E+01  3.67E-02  5.97E+00  1.94E+01 -1.10E+01  9.70E+00
 
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
 #CPUT: Total CPU Time in Seconds,       21.678
Stop Time:
Sat Sep 25 14:11:18 CDT 2021
