Wed Sep 29 13:38:46 CDT 2021
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
$DATA ../../../../data/spa/A3/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -43.4946077165983        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8960E+02  4.8112E+01  7.9329E+01 -2.4199E+01  1.8540E+02  7.3002E+01 -5.3075E+01 -3.4083E+01 -1.1793E+02 -1.2983E+02
            -2.9193E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1090.51095108430        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0222E+00  7.7324E-01  6.4387E-01  1.1563E+00  5.9073E-01  8.9174E-01  9.0058E-01  7.9735E-01  1.2410E+00  9.5071E-01
             1.8184E+00
 PARAMETER:  1.2198E-01 -1.5717E-01 -3.4026E-01  2.4522E-01 -4.2640E-01 -1.4577E-02 -4.7158E-03 -1.2646E-01  3.1588E-01  4.9449E-02
             6.9796E-01
 GRADIENT:   1.4410E+02  8.1162E+01  8.5648E+01  5.7930E+01 -5.5575E+01  1.6262E+01 -5.3603E+00  3.5987E+00 -7.9803E-01  9.9197E+00
            -7.1839E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1098.39731866471        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      221
 NPARAMETR:  1.0277E+00  5.3888E-01  7.4350E-01  1.3754E+00  5.5673E-01  9.0076E-01  6.1890E-01  5.7576E-01  1.2367E+00  9.5780E-01
             1.8512E+00
 PARAMETER:  1.2736E-01 -5.1826E-01 -1.9639E-01  4.1871E-01 -4.8568E-01 -4.5173E-03 -3.7981E-01 -4.5206E-01  3.1248E-01  5.6880E-02
             7.1582E-01
 GRADIENT:   1.6374E+01  7.7026E+01  1.5613E+02  8.2397E+01 -1.6239E+02  1.2818E+01 -2.5708E+00 -7.0133E+00 -1.9508E+00 -9.0440E+00
            -7.0637E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1262.89602723124        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  1.0346E+00  6.7420E-01  6.0677E-01  1.1071E+00  5.7648E-01  1.0633E+00  2.4130E-01  1.2618E-01  1.1495E+00  1.0989E+00
             2.8576E+00
 PARAMETER:  1.3401E-01 -2.9423E-01 -3.9960E-01  2.0177E-01 -4.5082E-01  1.6140E-01 -1.3217E+00 -1.9700E+00  2.3933E-01  1.9430E-01
             1.1500E+00
 GRADIENT:   2.6689E+00 -1.8989E+01  1.8192E+01 -1.3003E+02 -7.8634E+00  6.5922E+01 -1.4032E-01  2.5718E-01 -2.9772E+01  4.7744E+01
            -1.4802E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1314.34469529733        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      577
 NPARAMETR:  1.0691E+00  4.0484E-01  3.1223E-01  1.2635E+00  3.1742E-01  8.5955E-01  1.0000E-02  2.4499E-01  9.6704E-01  5.1214E-01
             3.6263E+00
 PARAMETER:  1.6680E-01 -8.0427E-01 -1.0640E+00  3.3386E-01 -1.0475E+00 -5.1341E-02 -8.1956E+00 -1.3065E+00  6.6486E-02 -5.6915E-01
             1.3882E+00
 GRADIENT:   2.1798E+01  2.5073E+01  2.5097E+01  5.2792E+01 -6.1860E+01 -2.0372E+00  0.0000E+00  8.8644E-01 -3.1053E+01  1.1622E+01
             2.7621E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1326.37196182320        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      752
 NPARAMETR:  1.0372E+00  5.2471E-01  1.9230E-01  1.0639E+00  2.8415E-01  8.8509E-01  1.0000E-02  2.2826E-01  1.2694E+00  2.9420E-01
             3.2766E+00
 PARAMETER:  1.3654E-01 -5.4492E-01 -1.5487E+00  1.6195E-01 -1.1583E+00 -2.2064E-02 -1.2885E+01 -1.3773E+00  3.3851E-01 -1.1235E+00
             1.2868E+00
 GRADIENT:  -1.6999E+00 -5.0774E+00  8.9465E+00  1.2820E+01  2.9666E+00 -2.1064E+00  0.0000E+00 -1.0824E+00 -2.7259E+00 -3.9612E+00
            -3.1158E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1332.50411953227        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      930
 NPARAMETR:  1.0095E+00  7.7078E-01  1.4875E-01  9.0687E-01  3.2443E-01  9.0078E-01  1.0000E-02  8.5470E-01  1.4457E+00  2.5873E-01
             3.2967E+00
 PARAMETER:  1.0946E-01 -1.6035E-01 -1.8055E+00  2.2402E-03 -1.0257E+00 -4.4936E-03 -1.5710E+01 -5.7003E-02  4.6860E-01 -1.2520E+00
             1.2929E+00
 GRADIENT:  -2.3452E+01  6.3971E+01  2.9576E+01 -5.7059E+00 -6.0011E+01  6.0297E+00  0.0000E+00 -1.1188E+01 -3.2261E+00 -2.2070E-01
             1.2069E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1345.49368922496        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1107
 NPARAMETR:  1.0201E+00  5.4727E-01  1.2587E-01  9.5453E-01  2.5193E-01  8.7816E-01  1.0000E-02  1.3084E+00  1.5509E+00  2.2227E-01
             2.9660E+00
 PARAMETER:  1.1988E-01 -5.0282E-01 -1.9725E+00  5.3462E-02 -1.2786E+00 -2.9924E-02 -1.7200E+01  3.6879E-01  5.3882E-01 -1.4039E+00
             1.1872E+00
 GRADIENT:   1.7437E+01 -3.5891E+01 -2.5906E+01 -4.0426E+00  5.4713E+01 -2.4427E-02  0.0000E+00 -6.5634E+00 -7.4740E+00 -6.9033E-01
             1.8008E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1347.55750257519        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1282
 NPARAMETR:  1.0134E+00  5.6365E-01  1.3722E-01  9.6519E-01  2.5886E-01  8.6976E-01  1.0000E-02  1.4824E+00  1.6020E+00  2.4693E-01
             2.7925E+00
 PARAMETER:  1.1334E-01 -4.7332E-01 -1.8862E+00  6.4571E-02 -1.2515E+00 -3.9533E-02 -1.5832E+01  4.9368E-01  5.7128E-01 -1.2986E+00
             1.1269E+00
 GRADIENT:  -9.6256E-03 -8.6230E-02  3.0303E+00 -3.2681E+00 -2.2945E+00  3.5928E-01  0.0000E+00 -3.4828E-01  6.6120E-02  1.6392E-01
             1.2024E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1347.57823180360        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:     1408
 NPARAMETR:  1.0127E+00  5.6710E-01  1.3588E-01  9.6688E-01  2.5944E-01  8.6902E-01  1.0000E-02  1.4864E+00  1.6095E+00  2.4773E-01
             2.7880E+00
 PARAMETER:  1.1258E-01 -4.6722E-01 -1.8959E+00  6.6319E-02 -1.2492E+00 -4.0386E-02 -1.5974E+01  4.9635E-01  5.7590E-01 -1.2954E+00
             1.1253E+00
 GRADIENT:  -1.0265E-01 -3.6572E-01 -5.3195E-01 -2.9454E-01  9.6582E-01  6.5170E-02  0.0000E+00  1.1332E-01  3.9245E-01  8.8638E-02
             2.7001E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1408
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.8847E-03 -5.3483E-04  1.7520E-02 -4.0707E-03  1.5298E-02
 SE:             2.8624E-02  2.6592E-04  2.1425E-02  2.6813E-02  1.0337E-02
 N:                     100         100         100         100         100

 P VAL.:         8.6450E-01  4.4302E-02  4.1351E-01  8.7933E-01  1.3888E-01

 ETASHRINKSD(%)  4.1054E+00  9.9109E+01  2.8222E+01  1.0172E+01  6.5371E+01
 ETASHRINKVR(%)  8.0422E+00  9.9992E+01  4.8479E+01  1.9309E+01  8.8008E+01
 EBVSHRINKSD(%)  4.0577E+00  9.9112E+01  2.9496E+01  8.6190E+00  6.6034E+01
 EBVSHRINKVR(%)  7.9508E+00  9.9992E+01  5.0292E+01  1.6495E+01  8.8463E+01
 RELATIVEINF(%)  8.8504E+01  2.6881E-04  1.3588E+01  6.8422E+01  2.9893E-01
 EPSSHRINKSD(%)  3.3765E+01
 EPSSHRINKVR(%)  5.6129E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1347.5782318035958     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -612.42740523985765     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1347.578       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  5.67E-01  1.36E-01  9.67E-01  2.59E-01  8.69E-01  1.00E-02  1.49E+00  1.61E+00  2.48E-01  2.79E+00
 


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
+        1.33E+03
 
 TH 2
+       -3.44E+00  3.03E+03
 
 TH 3
+       -4.88E+02  3.78E+03  1.31E+04
 
 TH 4
+       -2.70E+01  1.98E+02 -3.35E+02  3.83E+02
 
 TH 5
+        3.01E+02 -9.66E+03 -1.48E+04 -3.82E+02  3.32E+04
 
 TH 6
+        1.48E+00 -2.60E+01  1.20E+01 -8.67E+00  1.08E+02  2.23E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        5.41E-01  8.81E+00 -6.60E+01  1.59E+00 -3.85E+01  6.23E+00  0.00E+00  2.46E+01
 
 TH 9
+        1.56E+01 -6.14E+01  7.56E+01 -8.07E+00  2.02E+02  3.66E+00  0.00E+00 -1.62E+00  4.72E+01
 
 TH10
+        1.30E+00 -1.79E+02 -3.64E+01 -3.02E+00  5.58E+02  3.53E+00  0.00E+00  1.58E+01  7.91E+00  5.73E+01
 
 TH11
+       -2.17E+01 -5.11E+01  1.53E+01 -1.07E+00  1.10E+02  2.23E-02  0.00E+00  8.58E+00  7.08E+00  1.61E+01  3.48E+01
 
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
 #CPUT: Total CPU Time in Seconds,       25.605
Stop Time:
Wed Sep 29 13:39:13 CDT 2021
