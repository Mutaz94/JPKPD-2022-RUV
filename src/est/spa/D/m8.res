Sat Sep 18 15:05:53 CDT 2021
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
$DATA ../../../../data/spa/D/dat8.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m8.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1286.90065461848        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.4238E+01 -4.4278E+01 -5.4688E+00 -6.4777E+01  5.5327E+01 -3.1443E+02 -1.6879E+02 -2.1250E+01 -2.1131E+02 -3.8487E+01
            -5.2609E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1468.32370255376        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.5812E-01  1.3054E+00  9.7870E-01  8.7210E-01  1.0674E+00  1.4237E+00  2.6235E+00  1.1357E+00  1.1985E+00  9.5773E-01
             8.9463E-01
 PARAMETER:  5.7219E-02  3.6649E-01  7.8466E-02 -3.6846E-02  1.6526E-01  4.5328E-01  1.0645E+00  2.2721E-01  2.8108E-01  5.6816E-02
            -1.1345E-02
 GRADIENT:  -1.5146E+00  5.8267E+01  2.7633E+01 -1.5307E+01 -2.2873E+01  1.1626E+00  6.8695E+01 -6.1412E+00 -1.6568E-01 -2.4067E+00
            -1.6522E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1472.26958267091        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  1.0345E+00  9.6042E-01  1.1885E+00  1.0663E+00  1.0036E+00  1.4828E+00  3.1666E+00  1.3778E+00  1.1869E+00  8.8329E-01
             9.1865E-01
 PARAMETER:  1.3393E-01  5.9618E-02  2.7265E-01  1.6421E-01  1.0355E-01  4.9395E-01  1.2526E+00  4.2047E-01  2.7136E-01 -2.4098E-02
             1.5155E-02
 GRADIENT:   2.8985E+01  2.2117E+01  2.4057E+01 -1.4302E+01 -3.3806E+01 -4.5183E+01  3.0568E+01 -3.6332E+00  1.0139E+01 -2.1860E+00
            -4.9260E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1479.96164207386        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  1.0068E+00  9.9903E-01  7.7109E-01  9.8550E-01  9.0706E-01  1.6034E+00  2.5449E+00  9.8767E-01  1.0491E+00  7.7565E-01
             9.2081E-01
 PARAMETER:  1.0682E-01  9.9034E-02 -1.5995E-01  8.5396E-02  2.4520E-03  5.7212E-01  1.0341E+00  8.7598E-02  1.4793E-01 -1.5405E-01
             1.7498E-02
 GRADIENT:  -1.4100E+00 -3.6248E+00 -2.1200E+01  2.2388E+01  2.9847E+01 -9.0368E+00  2.2435E+00  4.0848E+00 -3.1661E+00 -1.2720E+00
            -2.8894E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1481.15504912343        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      629
 NPARAMETR:  1.0103E+00  1.1533E+00  6.8319E-01  8.8183E-01  9.1126E-01  1.6455E+00  2.2853E+00  7.1602E-01  1.1572E+00  7.7064E-01
             9.3056E-01
 PARAMETER:  1.1023E-01  2.4265E-01 -2.8099E-01 -2.5752E-02  7.0738E-03  5.9807E-01  9.2648E-01 -2.3404E-01  2.4597E-01 -1.6054E-01
             2.8027E-02
 GRADIENT:   5.3313E-01 -9.3490E-01 -8.3757E-03 -2.1802E+00  7.2002E-01  1.6821E+00  1.4898E+00  4.4588E-01 -2.9147E-02 -1.4394E-01
             1.5034E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1481.29363812756        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      804
 NPARAMETR:  1.0104E+00  1.2801E+00  5.5998E-01  7.9784E-01  8.9594E-01  1.6376E+00  2.0886E+00  3.4803E-01  1.2227E+00  7.3451E-01
             9.2808E-01
 PARAMETER:  1.1034E-01  3.4694E-01 -4.7985E-01 -1.2585E-01 -9.8873E-03  5.9320E-01  8.3649E-01 -9.5547E-01  3.0103E-01 -2.0856E-01
             2.5359E-02
 GRADIENT:  -1.0565E-01 -4.8923E-01 -1.3011E+00 -1.7437E+00 -5.0335E-01  1.2115E-01  4.5635E-01  3.3252E-01  6.7054E-01  3.6967E-01
             7.6365E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1481.37945508844        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      980
 NPARAMETR:  1.0102E+00  1.2881E+00  5.3956E-01  7.9211E-01  8.8668E-01  1.6345E+00  2.0715E+00  1.2769E-01  1.2192E+00  7.2992E-01
             9.2721E-01
 PARAMETER:  1.1011E-01  3.5318E-01 -5.1699E-01 -1.3306E-01 -2.0275E-02  5.9136E-01  8.2827E-01 -1.9581E+00  2.9822E-01 -2.1482E-01
             2.4424E-02
 GRADIENT:  -3.4406E-01  3.0042E-02  8.1325E-02 -2.2857E-01 -8.4577E-01 -5.1438E-01 -4.6405E-01  2.8733E-02 -9.3559E-02  2.7365E-01
             3.8653E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1481.39386947475        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  1.0105E+00  1.2918E+00  5.3626E-01  7.8941E-01  8.8689E-01  1.6363E+00  2.0699E+00  2.9355E-02  1.2230E+00  7.2807E-01
             9.2650E-01
 PARAMETER:  1.1043E-01  3.5605E-01 -5.2314E-01 -1.3647E-01 -2.0031E-02  5.9243E-01  8.2752E-01 -3.4283E+00  3.0133E-01 -2.1736E-01
             2.3662E-02
 GRADIENT:  -7.4651E-02 -4.1895E-02  1.9972E-02  3.6127E-02 -9.3798E-02 -7.6651E-02 -1.7784E-02  1.4459E-03  6.5366E-02  6.4283E-02
             7.4589E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1481.39463576392        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1331
 NPARAMETR:  1.0106E+00  1.2927E+00  5.3580E-01  7.8885E-01  8.8702E-01  1.6366E+00  2.0690E+00  1.0000E-02  1.2233E+00  7.2747E-01
             9.2637E-01
 PARAMETER:  1.1052E-01  3.5674E-01 -5.2400E-01 -1.3718E-01 -1.9882E-02  5.9263E-01  8.2706E-01 -4.5804E+00  3.0156E-01 -2.1819E-01
             2.3519E-02
 GRADIENT:  -3.3362E-03 -1.4891E-02  1.2528E-02  6.5803E-03 -3.1089E-03  4.3955E-03 -1.7537E-03  0.0000E+00  9.4183E-03  8.3946E-03
             6.1789E-03

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1481.39463652200        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1388
 NPARAMETR:  1.0106E+00  1.2929E+00  5.3565E-01  7.8872E-01  8.8701E-01  1.6366E+00  2.0687E+00  1.0000E-02  1.2234E+00  7.2724E-01
             9.2636E-01
 PARAMETER:  1.1053E-01  3.5691E-01 -5.2427E-01 -1.3734E-01 -1.9901E-02  5.9262E-01  8.2693E-01 -4.5776E+00  3.0162E-01 -2.1850E-01
             2.3506E-02
 GRADIENT:  -7.2275E-04 -3.7494E-04  8.6449E-04 -8.3803E-04  8.1846E-04 -3.3663E-04  4.7679E-04  0.0000E+00  7.8949E-05 -4.3969E-04
            -6.5966E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1388
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4065E-04 -1.1176E-03 -5.2695E-04  5.6667E-03 -1.5412E-02
 SE:             2.9958E-02  2.7067E-02  1.7167E-04  2.4157E-02  1.8340E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9625E-01  9.6706E-01  2.1439E-03  8.1453E-01  4.0071E-01

 ETASHRINKSD(%)  1.0000E-10  9.3234E+00  9.9425E+01  1.9072E+01  3.8559E+01
 ETASHRINKVR(%)  1.0000E-10  1.7778E+01  9.9997E+01  3.4507E+01  6.2250E+01
 EBVSHRINKSD(%)  1.4095E-01  7.3746E+00  9.9550E+01  2.0402E+01  3.9544E+01
 EBVSHRINKVR(%)  2.8169E-01  1.4205E+01  9.9998E+01  3.6642E+01  6.3450E+01
 RELATIVEINF(%)  9.9701E+01  2.3881E+01  4.2285E-04  1.1075E+01  6.2653E+00
 EPSSHRINKSD(%)  4.5711E+01
 EPSSHRINKVR(%)  7.0528E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1481.3946365219983     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -746.24380995826016     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1481.395       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.29E+00  5.36E-01  7.89E-01  8.87E-01  1.64E+00  2.07E+00  1.00E-02  1.22E+00  7.27E-01  9.26E-01
 


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
+        4.08E+02
 
 TH 2
+        1.08E-01  1.36E+02
 
 TH 3
+        3.65E+00  1.14E+02  9.84E+02
 
 TH 4
+       -2.56E+00  7.27E+01 -5.30E+02  7.81E+02
 
 TH 5
+       -7.85E-01 -1.37E+02 -7.55E+02  5.61E+02  9.94E+02
 
 TH 6
+        9.92E-02 -6.49E-02 -8.53E-01 -6.99E-01 -9.32E-01  7.43E+01
 
 TH 7
+        1.04E-01  8.65E+00 -4.87E+01 -1.12E+01  2.42E+01 -1.67E-01  3.28E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.37E-01 -7.76E+00 -3.70E+01  4.71E+01 -1.88E+01 -1.86E-01  5.19E+00  0.00E+00  6.09E+01
 
 TH10
+        5.51E-01 -7.94E+00 -5.60E+01 -2.32E+01 -9.04E+01 -1.02E-01  2.73E+00  0.00E+00  2.00E+01  7.81E+01
 
 TH11
+       -7.64E-01 -5.90E+00 -3.48E+01 -4.88E+00 -1.33E+01  1.41E+00  3.26E+00  0.00E+00  7.93E+00  2.13E+01  2.44E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.748
Stop Time:
Sat Sep 18 15:06:20 CDT 2021
