Sat Sep 18 12:34:11 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat95.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1601.61501377418        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0160E+01 -7.4357E+01 -3.9677E+01 -1.0063E+02  7.8915E+01  1.1467E+01 -1.8175E+01  1.6223E+01 -8.1241E+01 -1.8796E-01
            -3.4286E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1620.03564634850        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0200E+00  9.7140E-01  1.0830E+00  1.1030E+00  9.7717E-01  9.7608E-01  1.0110E+00  8.3555E-01  1.3168E+00  9.3294E-01
             1.0396E+00
 PARAMETER:  1.1978E-01  7.0987E-02  1.7978E-01  1.9803E-01  7.6909E-02  7.5791E-02  1.1092E-01 -7.9670E-02  3.7522E-01  3.0581E-02
             1.3883E-01
 GRADIENT:   9.1342E+01  2.2037E+00 -7.0783E+00  2.9982E+01  3.4316E+01  7.1005E-01  2.7446E+00  8.5276E+00  1.0955E+01 -1.3009E+01
            -1.5802E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1624.76026603919        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0145E+00  8.6781E-01  1.0541E+00  1.1651E+00  9.4740E-01  9.5699E-01  8.1773E-01  3.0833E-01  1.2625E+00  1.1443E+00
             1.1157E+00
 PARAMETER:  1.1437E-01 -4.1786E-02  1.5272E-01  2.5280E-01  4.5969E-02  5.6040E-02 -1.0123E-01 -1.0766E+00  3.3310E-01  2.3476E-01
             2.0946E-01
 GRADIENT:   7.2347E+01 -1.6618E+00 -2.8197E+01  3.2452E+01  3.1691E+01 -5.8312E+00  2.0334E+00  1.6978E+00  5.5656E+00  1.0185E+01
             2.2671E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1626.46919869409        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9041E-01  8.2296E-01  1.0316E+00  1.1797E+00  8.9188E-01  9.6276E-01  9.0403E-01  3.0270E-01  1.2138E+00  1.0598E+00
             1.0700E+00
 PARAMETER:  9.0363E-02 -9.4848E-02  1.3114E-01  2.6524E-01 -1.4423E-02  6.2044E-02 -8.9633E-04 -1.0950E+00  2.9372E-01  1.5807E-01
             1.6763E-01
 GRADIENT:   1.4341E+01  5.2257E+00 -4.3188E+00  2.1625E+01  1.6816E+00 -2.8253E+00  6.2422E-01  1.4291E+00  1.8915E+00  2.4749E+00
             4.1388E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1626.52372071080        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.8343E-01  9.0314E-01  9.5511E-01  1.1154E+00  8.8016E-01  9.7362E-01  1.1777E+00  2.0246E-01  1.1972E+00  9.7440E-01
             1.0573E+00
 PARAMETER:  8.3289E-02 -1.8822E-03  5.4075E-02  2.0918E-01 -2.7649E-02  7.3262E-02  2.6353E-01 -1.4972E+00  2.7996E-01  7.4062E-02
             1.5569E-01
 GRADIENT:  -2.8717E+00  4.1367E+00  2.5330E+00  2.0311E+00 -4.6419E+00  6.6265E-01  1.0685E+00  7.7743E-01 -1.6340E+00  4.6596E+00
            -5.3611E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1626.80581891502        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.8345E-01  8.9260E-01  8.7062E-01  1.1043E+00  8.3256E-01  9.7414E-01  1.2615E+00  9.6389E-02  1.1803E+00  8.5382E-01
             1.0586E+00
 PARAMETER:  8.3312E-02 -1.3613E-02 -3.8554E-02  1.9919E-01 -8.3252E-02  7.3805E-02  3.3228E-01 -2.2394E+00  2.6574E-01 -5.8041E-02
             1.5691E-01
 GRADIENT:  -4.5176E+00 -3.3926E+00 -1.7151E+00 -6.8179E+00  2.2775E+00  4.7419E-01 -4.8603E-01  1.9982E-01 -8.2276E-01 -5.2403E-01
            -5.1099E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1628.05121383250        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      495
 NPARAMETR:  9.9544E-01  7.4771E-01  8.7441E-01  1.2147E+00  7.7594E-01  9.7467E-01  1.4834E+00  4.2913E-02  1.1044E+00  8.3236E-01
             1.0590E+00
 PARAMETER:  9.5428E-02 -1.9074E-01 -3.4204E-02  2.9446E-01 -1.5369E-01  7.4342E-02  4.9434E-01 -3.0486E+00  1.9933E-01 -8.3487E-02
             1.5732E-01
 GRADIENT:  -9.3885E+00  6.2099E+00  1.4417E+00  6.0580E+00 -4.6160E+00 -1.0775E+00 -1.2683E+00  3.8546E-02 -4.2756E-01 -1.6917E+00
             6.5501E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1628.91455118276        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  9.9568E-01  5.0185E-01  9.0365E-01  1.3537E+00  7.1861E-01  9.7216E-01  1.9900E+00  1.0000E-02  1.0065E+00  8.7212E-01
             1.0517E+00
 PARAMETER:  9.5675E-02 -5.8946E-01 -1.3102E-03  4.0284E-01 -2.3043E-01  7.1766E-02  7.8813E-01 -4.5189E+00  1.0643E-01 -3.6827E-02
             1.5045E-01
 GRADIENT:  -1.5419E+00  5.7963E-01  2.7959E-01 -3.7094E-01 -3.9738E-01 -1.6948E-01 -3.5992E-01  0.0000E+00 -3.0977E-01 -1.6624E-01
            -1.8213E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1628.99065420955        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      845
 NPARAMETR:  9.9466E-01  4.0581E-01  9.1014E-01  1.4092E+00  6.9699E-01  9.7027E-01  2.2870E+00  1.0000E-02  9.7633E-01  8.9139E-01
             1.0493E+00
 PARAMETER:  9.4647E-02 -8.0188E-01  5.8455E-03  4.4299E-01 -2.6098E-01  6.9821E-02  9.2723E-01 -5.4839E+00  7.6043E-02 -1.4973E-02
             1.4812E-01
 GRADIENT:  -1.3368E-01  1.5199E-01  5.5369E-02  5.2153E-01 -2.8229E-01 -4.0887E-02  5.9710E-02  0.0000E+00  9.0276E-02  2.5584E-02
             3.9295E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1628.99082874691        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  9.9471E-01  4.0574E-01  9.1076E-01  1.4089E+00  6.9743E-01  9.7036E-01  2.2856E+00  1.0000E-02  9.7612E-01  8.9184E-01
             1.0492E+00
 PARAMETER:  9.4693E-02 -8.0205E-01  6.5290E-03  4.4279E-01 -2.6035E-01  6.9912E-02  9.2664E-01 -5.4814E+00  7.5834E-02 -1.4471E-02
             1.4806E-01
 GRADIENT:  -9.4440E-03 -5.1318E-04 -1.8624E-03 -1.9828E-02 -7.4579E-03 -1.9256E-03  3.3421E-03  0.0000E+00  1.2992E-02  7.4307E-03
             6.8207E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      937
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6219E-04  1.8916E-02 -4.2444E-04 -1.3469E-02 -5.9132E-03
 SE:             2.9810E-02  1.6434E-02  2.0472E-04  2.6928E-02  2.3868E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9298E-01  2.4974E-01  3.8144E-02  6.1694E-01  8.0433E-01

 ETASHRINKSD(%)  1.3108E-01  4.4943E+01  9.9314E+01  9.7869E+00  2.0038E+01
 ETASHRINKVR(%)  2.6199E-01  6.9687E+01  9.9995E+01  1.8616E+01  3.6061E+01
 EBVSHRINKSD(%)  5.0908E-01  5.0652E+01  9.9273E+01  8.5611E+00  1.5843E+01
 EBVSHRINKVR(%)  1.0156E+00  7.5647E+01  9.9995E+01  1.6389E+01  2.9175E+01
 RELATIVEINF(%)  9.8196E+01  3.5630E+00  7.0854E-04  1.8908E+01  6.2354E+00
 EPSSHRINKSD(%)  4.3322E+01
 EPSSHRINKVR(%)  6.7876E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1628.9908287469145     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -893.84000218317635     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1628.991       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  4.06E-01  9.11E-01  1.41E+00  6.97E-01  9.70E-01  2.29E+00  1.00E-02  9.76E-01  8.92E-01  1.05E+00
 


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
+        1.18E+03
 
 TH 2
+       -1.89E+01  3.72E+02
 
 TH 3
+        1.36E+01  2.05E+02  7.58E+02
 
 TH 4
+       -4.39E+00  2.34E+02 -1.15E+02  4.96E+02
 
 TH 5
+       -5.49E+00 -4.60E+02 -1.11E+03  1.09E+02  2.03E+03
 
 TH 6
+       -1.90E+00 -4.67E+00  9.75E+00 -9.16E-01 -1.27E-01  2.06E+02
 
 TH 7
+        6.92E-01  2.86E+01 -6.70E-01 -1.98E+00 -3.65E+00  5.92E-02  8.16E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.41E+00 -3.68E+01  1.25E+00  1.06E+01 -1.55E+01 -1.32E+00  1.31E+00  0.00E+00  1.60E+02
 
 TH10
+        6.90E-01  2.71E+01 -7.24E+01 -2.54E+01 -1.98E+01  8.54E-01  5.61E+00  0.00E+00 -6.57E-01  1.12E+02
 
 TH11
+       -1.08E+01 -7.58E+00 -4.04E+01 -7.03E+00  2.17E+01  1.69E+00  1.29E+00  0.00E+00  6.70E+00  2.79E+01  2.00E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.120
Stop Time:
Sat Sep 18 12:34:29 CDT 2021
