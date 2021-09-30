Wed Sep 29 16:03:43 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat68.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1733.38426275356        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0752E+02  4.6718E+01 -2.4660E+01  1.2422E+02  3.7445E+01  4.7401E+01  2.5324E+01  9.4222E+00  5.0547E+01 -2.4510E+00
             6.2448E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1747.55129264294        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0360E+00  9.7408E-01  1.0185E+00  9.8385E-01  1.0078E+00  1.0027E+00  8.3166E-01  9.4065E-01  7.6248E-01  1.0718E+00
             8.4687E-01
 PARAMETER:  1.3536E-01  7.3735E-02  1.1837E-01  8.3721E-02  1.0778E-01  1.0268E-01 -8.4328E-02  3.8818E-02 -1.7117E-01  1.6934E-01
            -6.6210E-02
 GRADIENT:   4.5407E+01 -1.6686E+01 -1.9147E+01 -1.3767E+01  2.9936E+01 -3.1436E-02  7.6843E-01  5.6855E+00 -1.0374E+01 -7.6603E-02
            -9.5982E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1749.16432481545        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0262E+00  7.7561E-01  1.1214E+00  1.1313E+00  9.4939E-01  1.0105E+00  7.1223E-01  8.3136E-01  7.9084E-01  1.0861E+00
             8.2839E-01
 PARAMETER:  1.2589E-01 -1.5411E-01  2.1458E-01  2.2336E-01  4.8068E-02  1.1042E-01 -2.3936E-01 -8.4696E-02 -1.3466E-01  1.8263E-01
            -8.8273E-02
 GRADIENT:   2.8133E+01  1.9045E+01 -3.7051E+00  4.5575E+01  6.2716E-01  4.2149E+00 -2.9427E-02  7.6721E-01  5.1570E+00  9.9518E-01
            -1.0725E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1749.98334257038        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0112E+00  7.3500E-01  1.0540E+00  1.1371E+00  8.9710E-01  9.9869E-01  8.6605E-01  6.9750E-01  7.3935E-01  1.0376E+00
             8.5281E-01
 PARAMETER:  1.1116E-01 -2.0789E-01  1.5260E-01  2.2852E-01 -8.5834E-03  9.8687E-02 -4.3813E-02 -2.6025E-01 -2.0199E-01  1.3691E-01
            -5.9217E-02
 GRADIENT:  -5.3780E+00  5.3687E+00  5.0983E+00  1.3694E+00 -8.0159E+00 -7.8643E-02  5.9245E-02 -8.5935E-01  3.2000E-02 -2.2804E-02
             7.9287E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1750.48572256547        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0136E+00  4.8316E-01  1.2614E+00  1.3077E+00  8.9955E-01  9.9875E-01  1.1688E+00  9.2987E-01  6.6292E-01  1.0575E+00
             8.4592E-01
 PARAMETER:  1.1353E-01 -6.2740E-01  3.3223E-01  3.6828E-01 -5.8585E-03  9.8753E-02  2.5599E-01  2.7285E-02 -3.1109E-01  1.5593E-01
            -6.7332E-02
 GRADIENT:   9.5664E+00  8.7961E+00  6.8489E-01  2.2317E+01 -4.0334E+00  1.6890E+00  9.4480E-01  5.8732E-01  1.6657E+00 -4.7758E-01
            -2.3374E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1750.87964714174        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0079E+00  3.0729E-01  1.4900E+00  1.4301E+00  9.2931E-01  9.9005E-01  1.4140E+00  1.1450E+00  6.1801E-01  1.1005E+00
             8.4837E-01
 PARAMETER:  1.0789E-01 -1.0799E+00  4.9878E-01  4.5777E-01  2.6683E-02  9.0002E-02  4.4640E-01  2.3542E-01 -3.8125E-01  1.9572E-01
            -6.4433E-02
 GRADIENT:   4.0120E+00  8.6674E+00 -1.8073E-01  3.8264E+01 -2.1666E+00 -1.8554E-01  2.5498E-01  7.0931E-02 -8.2880E-01 -1.1464E-01
            -1.4135E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1751.57973928697        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.0014E+00  1.3960E-01  1.7810E+00  1.5394E+00  9.6574E-01  9.8329E-01  1.5290E+00  1.3947E+00  5.8842E-01  1.1394E+00
             8.5294E-01
 PARAMETER:  1.0142E-01 -1.8690E+00  6.7717E-01  5.3137E-01  6.5142E-02  8.3148E-02  5.2464E-01  4.3269E-01 -4.3031E-01  2.3047E-01
            -5.9070E-02
 GRADIENT:  -4.7323E+00  3.2766E+00  2.0162E+00  1.6254E+01 -1.9896E+00 -1.4517E+00  4.2635E-03 -8.7109E-01 -1.1663E-01  1.0369E-01
             9.3024E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1751.76133295284        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1253
 NPARAMETR:  1.0016E+00  7.9214E-02  1.8894E+00  1.5823E+00  9.7359E-01  9.8461E-01  1.5004E+00  1.4955E+00  5.7507E-01  1.1430E+00
             8.5153E-01
 PARAMETER:  1.0161E-01 -2.4356E+00  7.3625E-01  5.5888E-01  7.3237E-02  8.4495E-02  5.0571E-01  5.0247E-01 -4.5326E-01  2.3364E-01
            -6.0724E-02
 GRADIENT:  -2.3032E+00  2.2375E+00  2.4328E+00  2.1356E+01 -4.6106E+00 -4.5163E-01 -5.5992E-03 -4.8033E-01 -8.5434E-01 -1.1236E-01
             2.2034E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1752.13838841886        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1432
 NPARAMETR:  1.0017E+00  2.0222E-02  1.9332E+00  1.6174E+00  9.7044E-01  9.8517E-01  1.2247E+00  1.5393E+00  5.6141E-01  1.1402E+00
             8.5008E-01
 PARAMETER:  1.0172E-01 -3.8010E+00  7.5918E-01  5.8085E-01  6.9992E-02  8.5062E-02  3.0269E-01  5.3134E-01 -4.7731E-01  2.3117E-01
            -6.2424E-02
 GRADIENT:   1.3535E-01  4.6284E-01  1.0721E+00  7.7475E+00 -2.3132E+00  1.8946E-01  4.2222E-04 -1.1675E-01 -1.0487E+00 -2.3014E-01
            -2.6001E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1752.24231459660        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1526
 NPARAMETR:  1.0016E+00  1.0000E-02  1.9124E+00  1.6217E+00  9.6506E-01  9.8474E-01  1.0260E+00  1.5250E+00  5.6023E-01  1.1380E+00
             8.4992E-01
 PARAMETER:  1.0162E-01 -4.6401E+00  7.4834E-01  5.8347E-01  6.4432E-02  8.4626E-02  1.2564E-01  5.2197E-01 -4.7941E-01  2.2930E-01
            -6.2607E-02
 GRADIENT:   3.9860E-01  0.0000E+00 -6.1456E-02  2.1228E+00 -3.5543E-01  7.0944E-02  1.4959E-04  4.5234E-02 -3.4341E-01 -4.9085E-02
            -1.7434E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1526
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.0680E-04 -4.7720E-04 -3.6741E-02 -1.7692E-02 -4.2881E-02
 SE:             2.9909E-02  2.2762E-04  1.7794E-02  2.9052E-02  2.1497E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8115E-01  3.6034E-02  3.8940E-02  5.4255E-01  4.6068E-02

 ETASHRINKSD(%)  1.0000E-10  9.9237E+01  4.0389E+01  2.6728E+00  2.7983E+01
 ETASHRINKVR(%)  1.0000E-10  9.9994E+01  6.4465E+01  5.2741E+00  4.8135E+01
 EBVSHRINKSD(%)  3.0053E-01  9.9279E+01  4.5017E+01  3.1664E+00  2.3040E+01
 EBVSHRINKVR(%)  6.0015E-01  9.9995E+01  6.9769E+01  6.2326E+00  4.0771E+01
 RELATIVEINF(%)  9.7655E+01  2.1126E-04  8.0855E+00  4.1432E+00  1.2354E+01
 EPSSHRINKSD(%)  4.5677E+01
 EPSSHRINKVR(%)  7.0491E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1752.2423145966000     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1017.0914880328618     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1752.242       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.91E+00  1.62E+00  9.65E-01  9.85E-01  1.03E+00  1.52E+00  5.60E-01  1.14E+00  8.50E-01
 


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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.50E+00  0.00E+00  4.32E+01
 
 TH 4
+       -1.41E+01  0.00E+00 -3.25E+01  1.26E+03
 
 TH 5
+       -6.52E-01  0.00E+00 -1.12E+02 -7.88E+01  5.68E+02
 
 TH 6
+        1.76E+00  0.00E+00  2.99E-02 -2.62E+00  3.62E-02  2.03E+02
 
 TH 7
+        6.32E-02  0.00E+00  6.68E-03 -2.58E-02  3.21E-02  4.70E-02  2.85E-02
 
 TH 8
+        9.36E-02  0.00E+00 -1.66E+01 -5.17E+00 -6.25E+00  6.77E-02 -1.93E-02  2.12E+01
 
 TH 9
+        2.53E+00  0.00E+00  6.22E+00 -3.48E+01  2.76E+00 -1.90E+00  1.54E-02  6.26E-01  5.66E+02
 
 TH10
+        1.04E+00  0.00E+00 -1.73E+00 -3.37E+00 -6.72E+01  5.13E-01 -3.62E-02  8.93E+00  1.64E+00  6.29E+01
 
 TH11
+       -8.24E+00  0.00E+00 -5.15E+00 -1.94E+01 -2.24E+00  1.64E+00 -1.77E-01  6.22E+00  2.09E+01  8.45E+00  2.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.480
Stop Time:
Wed Sep 29 16:04:10 CDT 2021
