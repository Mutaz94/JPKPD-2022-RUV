Wed Sep 29 11:01:55 CDT 2021
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
$DATA ../../../../data/spa/B/dat20.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1656.35385342521        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2055E+02  1.3528E+01 -3.4211E+01  8.6139E+01  7.0091E+01  6.8700E+01 -6.6008E-01  8.4079E+00  1.3416E+01 -2.9937E+00
             1.7569E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1659.20423627069        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.7222E-01  1.0289E+00  1.0428E+00  9.8865E-01  9.7804E-01  1.0320E+00  1.0250E+00  9.5512E-01  9.9308E-01  9.9357E-01
             9.9698E-01
 PARAMETER:  7.1828E-02  1.2850E-01  1.4196E-01  8.8586E-02  7.7798E-02  1.3145E-01  1.2472E-01  5.4086E-02  9.3058E-02  9.3547E-02
             9.6971E-02
 GRADIENT:   1.5620E+00  5.5081E+00  6.2860E+00  5.7592E-01 -7.9609E+00  6.9512E+00 -2.6265E+00  3.6484E+00 -9.8827E-03 -2.4551E+00
            -1.5572E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1660.73260498161        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.7275E-01  8.9421E-01  9.2116E-01  1.0754E+00  8.7797E-01  9.9163E-01  1.2068E+00  4.7857E-01  9.1386E-01  9.8737E-01
             1.0021E+00
 PARAMETER:  7.2375E-02 -1.1819E-02  1.7878E-02  1.7272E-01 -3.0146E-02  9.1598E-02  2.8797E-01 -6.3696E-01  9.9223E-03  8.7294E-02
             1.0214E-01
 GRADIENT:   1.8196E+00  9.6703E+00 -3.8715E+00  2.2048E+01  5.9269E-01 -9.3703E+00 -1.9271E+00  7.3638E-01 -1.1907E+00  7.5358E+00
             6.4595E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1661.61482620734        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.7212E-01  8.2265E-01  8.0879E-01  1.0947E+00  7.8893E-01  1.0210E+00  1.3255E+00  3.5044E-01  8.9015E-01  8.2683E-01
             1.0036E+00
 PARAMETER:  7.1727E-02 -9.5219E-02 -1.1222E-01  1.9046E-01 -1.3708E-01  1.2076E-01  3.8178E-01 -9.4857E-01 -1.6364E-02 -9.0151E-02
             1.0355E-01
 GRADIENT:  -1.0662E+00 -1.2884E+00 -4.9399E+00  5.7384E+00  5.7336E+00  2.2469E+00 -2.2032E+00  3.9925E-01  1.2435E+00 -9.8568E-01
            -1.5362E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1661.72880239060        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.7187E-01  7.4819E-01  8.3683E-01  1.1391E+00  7.7492E-01  1.0139E+00  1.4644E+00  2.8648E-01  8.5785E-01  8.4349E-01
             1.0065E+00
 PARAMETER:  7.1469E-02 -1.9010E-01 -7.8130E-02  2.3027E-01 -1.5499E-01  1.1376E-01  4.8145E-01 -1.1501E+00 -5.3329E-02 -7.0210E-02
             1.0651E-01
 GRADIENT:   6.5569E-02  1.3573E+00  1.1019E+00  2.2207E-01 -2.0450E+00 -2.9932E-01  3.1576E-01  1.0783E-01 -2.3750E-01 -4.6042E-02
             6.1945E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1661.74790672100        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.7257E-01  7.3880E-01  8.3638E-01  1.1438E+00  7.7314E-01  1.0159E+00  1.4768E+00  2.1137E-01  8.5722E-01  8.4678E-01
             1.0074E+00
 PARAMETER:  7.2186E-02 -2.0272E-01 -7.8677E-02  2.3438E-01 -1.5730E-01  1.1578E-01  4.8986E-01 -1.4542E+00 -5.4065E-02 -6.6309E-02
             1.0740E-01
 GRADIENT:   1.7647E+00  1.5608E-01  5.7844E-01 -7.3609E-02  3.7827E-01  5.3474E-01  1.3966E-01  2.5933E-02  1.6928E-02 -7.0748E-01
            -3.8464E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1661.76036859127        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  9.7235E-01  7.3595E-01  8.3423E-01  1.1460E+00  7.7182E-01  1.0156E+00  1.4771E+00  9.2530E-02  8.5797E-01  8.5855E-01
             1.0095E+00
 PARAMETER:  7.1965E-02 -2.0659E-01 -8.1247E-02  2.3628E-01 -1.5901E-01  1.1546E-01  4.9007E-01 -2.2802E+00 -5.3186E-02 -5.2510E-02
             1.0942E-01
 GRADIENT:   1.2475E+00  2.3222E-01 -2.8358E-01  1.2055E+00  3.8905E-01  3.9966E-01  5.4102E-02  9.1101E-03  2.2204E-01  3.2516E-01
             1.4173E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1661.76400231913        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  9.7237E-01  7.3337E-01  8.3356E-01  1.1460E+00  7.7064E-01  1.0153E+00  1.4825E+00  2.0244E-02  8.5660E-01  8.5723E-01
             1.0096E+00
 PARAMETER:  7.1981E-02 -2.1011E-01 -8.2047E-02  2.3631E-01 -1.6054E-01  1.1523E-01  4.9375E-01 -3.7999E+00 -5.4787E-02 -5.4048E-02
             1.0953E-01
 GRADIENT:   1.3387E+00 -3.3939E-01  6.7145E-01 -1.2835E+00 -2.4660E-01  3.2118E-01  9.8769E-02  4.2730E-04 -1.7866E-02 -2.6901E-02
            -5.7265E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1661.76580809238        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1405             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7259E-01  7.3511E-01  8.3268E-01  1.1450E+00  7.7036E-01  1.0161E+00  1.4813E+00  1.0000E-02  8.5685E-01  8.5617E-01
             1.0095E+00
 PARAMETER:  7.2207E-02 -2.0774E-01 -8.3103E-02  2.3539E-01 -1.6090E-01  1.1594E-01  4.9291E-01 -5.0913E+00 -5.4495E-02 -5.5286E-02
             1.0945E-01
 GRADIENT:   4.3219E+02  4.6411E+01  4.9807E+00  2.6019E+02  1.6355E+01  8.5225E+01  2.1358E+01  0.0000E+00  1.3208E+01  6.3833E-01
             7.4832E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1661.76667606446        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1547
 NPARAMETR:  9.7261E-01  7.3658E-01  8.3193E-01  1.1441E+00  7.7011E-01  1.0161E+00  1.4797E+00  1.0000E-02  8.5725E-01  8.5592E-01
             1.0095E+00
 PARAMETER:  7.2205E-02 -2.0781E-01 -8.4278E-02  2.3529E-01 -1.6059E-01  1.1594E-01  4.9259E-01 -5.0913E+00 -5.4222E-02 -5.4888E-02
             1.0959E-01
 GRADIENT:  -3.2591E-02 -3.0710E-01 -1.6356E-01  8.6661E-01  7.8952E-01 -5.7680E-03  5.4941E-02  0.0000E+00 -2.7918E-02  7.2956E-02
             4.5511E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1547
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1190E-05  8.0556E-03 -4.3631E-04 -9.9946E-03 -8.5134E-03
 SE:             2.9832E-02  1.9885E-02  2.0260E-04  2.4924E-02  2.3579E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9917E-01  6.8540E-01  3.1273E-02  6.8842E-01  7.1805E-01

 ETASHRINKSD(%)  5.9018E-02  3.3381E+01  9.9321E+01  1.6502E+01  2.1009E+01
 ETASHRINKVR(%)  1.1800E-01  5.5619E+01  9.9995E+01  3.0281E+01  3.7604E+01
 EBVSHRINKSD(%)  4.2547E-01  3.4107E+01  9.9371E+01  1.6080E+01  1.9326E+01
 EBVSHRINKVR(%)  8.4912E-01  5.6581E+01  9.9996E+01  2.9574E+01  3.4917E+01
 RELATIVEINF(%)  9.8657E+01  3.3020E+00  4.7891E-04  6.8693E+00  5.8522E+00
 EPSSHRINKSD(%)  4.3496E+01
 EPSSHRINKVR(%)  6.8073E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1661.7666760644574     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -926.61584950071926     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1661.767       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.73E-01  7.35E-01  8.32E-01  1.14E+00  7.71E-01  1.02E+00  1.48E+00  1.00E-02  8.57E-01  8.57E-01  1.01E+00
 


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
+        1.13E+03
 
 TH 2
+       -9.49E+00  4.33E+02
 
 TH 3
+        1.82E+01  2.31E+02  6.88E+02
 
 TH 4
+       -8.50E+00  3.49E+02 -2.40E+02  8.21E+02
 
 TH 5
+       -4.03E+00 -4.44E+02 -9.46E+02  2.82E+02  1.66E+03
 
 TH 6
+        9.23E-02 -1.66E+00  3.61E+00 -2.74E+00 -1.42E+00  1.90E+02
 
 TH 7
+        1.22E+00  3.43E+01 -2.52E+00 -7.11E+00 -5.30E+00  2.97E-01  2.56E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.83E+00 -2.53E+01 -1.74E+01  1.61E+01  4.71E+00  4.93E-01  1.13E+01  0.00E+00  1.51E+02
 
 TH10
+       -1.32E+00 -5.71E+00 -7.44E+01 -2.59E+01 -4.65E+01 -2.38E-01  1.00E+01  0.00E+00  8.51E+00  1.14E+02
 
 TH11
+       -7.75E+00 -1.05E+01 -3.96E+01 -5.90E+00  2.01E+01  1.32E+00  2.88E+00  0.00E+00  6.96E+00  2.53E+01  2.17E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.492
Stop Time:
Wed Sep 29 11:02:24 CDT 2021
