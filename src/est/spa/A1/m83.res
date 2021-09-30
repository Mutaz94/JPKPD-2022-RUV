Wed Sep 29 12:23:55 CDT 2021
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
$DATA ../../../../data/spa/A1/dat83.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1364.82710944080        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4908E+02  8.6618E+01  7.3797E+01  5.5634E+01 -1.0521E+01  2.8157E+01  3.3269E+00 -2.9031E+01  3.4800E+00 -2.0051E+01
            -5.2447E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1465.14989593990        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1961E+00  9.5311E-01  8.9072E-01  1.0783E+00  9.3371E-01  1.2249E+00  9.2346E-01  1.0446E+00  9.2746E-01  9.5160E-01
             2.3729E+00
 PARAMETER:  2.7911E-01  5.1973E-02 -1.5730E-02  1.7536E-01  3.1410E-02  3.0289E-01  2.0376E-02  1.4364E-01  2.4697E-02  5.0391E-02
             9.6411E-01
 GRADIENT:   4.5899E+02  3.1329E+01 -1.6003E+00  6.3093E+01 -1.6337E+01  3.6329E+01  8.6003E+00  7.3676E+00  1.0676E+01  2.0392E+01
             9.7554E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1476.38774461062        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      179
 NPARAMETR:  1.1759E+00  8.8085E-01  3.9243E-01  1.0518E+00  5.6505E-01  1.2754E+00  5.4017E-01  8.0086E-01  1.0317E+00  4.3593E-01
             2.0605E+00
 PARAMETER:  2.6203E-01 -2.6870E-02 -8.3539E-01  1.5050E-01 -4.7084E-01  3.4325E-01 -5.1588E-01 -1.2206E-01  1.3118E-01 -7.3027E-01
             8.2294E-01
 GRADIENT:   2.0324E+02  3.7748E+01 -2.5852E+01  9.2417E+01  5.8284E+00  3.5609E+01 -9.1419E+00  4.7683E+00  1.5734E+01 -2.8012E+00
             3.4141E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1506.79736578644        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      357
 NPARAMETR:  1.0249E+00  8.1384E-01  6.2329E-01  1.0956E+00  7.0515E-01  1.0022E+00  5.0730E-01  6.6470E-01  1.0057E+00  7.3497E-01
             1.7728E+00
 PARAMETER:  1.2459E-01 -1.0600E-01 -3.7275E-01  1.9127E-01 -2.4935E-01  1.0216E-01 -5.7865E-01 -3.0842E-01  1.0565E-01 -2.0792E-01
             6.7256E-01
 GRADIENT:   3.7622E+01  1.8033E+01 -6.3017E+00  2.6504E+01  7.2243E+00 -1.2724E+01 -3.0199E+00  2.5053E-01  6.0623E+00  9.3072E-01
            -2.3987E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1514.03286918982        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.9702E-01  4.0755E-01  4.9945E-01  1.2950E+00  4.7132E-01  1.0481E+00  2.8151E-01  3.3398E-01  8.4657E-01  7.0832E-01
             1.7440E+00
 PARAMETER:  9.7013E-02 -7.9760E-01 -5.9426E-01  3.5852E-01 -6.5222E-01  1.4702E-01 -1.1676E+00 -9.9668E-01 -6.6565E-02 -2.4486E-01
             6.5619E-01
 GRADIENT:  -2.0876E+01  3.5423E+01  4.0729E+01  7.2254E+01 -6.7931E+01  3.0312E+00 -4.3353E-01  3.8417E-01 -1.2753E+00  2.0505E+00
            -9.4872E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1519.89682551361        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  1.0076E+00  2.3871E-01  3.9940E-01  1.2774E+00  3.7519E-01  1.0337E+00  1.9105E-01  1.5496E-01  8.1932E-01  7.0254E-01
             1.7261E+00
 PARAMETER:  1.0757E-01 -1.3325E+00 -8.1778E-01  3.4486E-01 -8.8033E-01  1.3315E-01 -1.5552E+00 -1.7646E+00 -9.9286E-02 -2.5305E-01
             6.4584E-01
 GRADIENT:   6.6642E+00  6.5871E+00  3.0851E+01 -2.4664E+01 -4.1341E+01 -2.0278E+00 -3.5577E-02  3.7188E-01  4.1253E+00  6.2592E+00
             3.1036E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1523.26903495964        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.9796E-01  6.5675E-02  4.2016E-01  1.3641E+00  3.7163E-01  1.0355E+00  9.8021E-02  2.8606E-02  7.5950E-01  7.0308E-01
             1.7224E+00
 PARAMETER:  9.7962E-02 -2.6230E+00 -7.6713E-01  4.1051E-01 -8.8987E-01  1.3485E-01 -2.2226E+00 -3.4541E+00 -1.7509E-01 -2.5228E-01
             6.4369E-01
 GRADIENT:   8.3867E+00  9.2229E-01  1.3104E+01 -1.2841E+01 -1.6548E+01  8.4950E-01 -4.5565E-04  4.6471E-03 -3.2740E+00  6.5430E-01
             4.2022E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1523.88912348794        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.9105E-01  1.0000E-02  4.3389E-01  1.4007E+00  3.7541E-01  1.0305E+00  3.5694E-02  1.0000E-02  7.5157E-01  7.1083E-01
             1.7204E+00
 PARAMETER:  9.1006E-02 -4.6293E+00 -7.3496E-01  4.3694E-01 -8.7975E-01  1.3003E-01 -3.2328E+00 -6.2428E+00 -1.8559E-01 -2.4132E-01
             6.4255E-01
 GRADIENT:   1.6863E+00  0.0000E+00  6.2411E+00  1.7303E+00 -9.1260E+00 -1.7024E-01 -1.0381E-06  0.0000E+00 -1.1283E+00  4.2810E-01
            -5.1959E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1523.90739206124        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  9.9051E-01  1.0000E-02  4.3329E-01  1.3997E+00  3.7612E-01  1.0312E+00  3.6081E-02  1.0000E-02  7.5426E-01  7.0899E-01
             1.7229E+00
 PARAMETER:  9.0466E-02 -4.6121E+00 -7.3635E-01  4.3624E-01 -8.7786E-01  1.3068E-01 -3.2220E+00 -6.2272E+00 -1.8202E-01 -2.4391E-01
             6.4404E-01
 GRADIENT:   7.1167E-01  0.0000E+00  4.3308E-02 -4.1847E-01 -4.4579E-02  1.0334E-01 -9.9472E-07  0.0000E+00  7.0639E-04  7.5278E-02
             2.7324E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1150
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.6759E-05 -7.6030E-06 -4.6203E-05 -7.3973E-03 -7.0607E-03
 SE:             2.9589E-02  6.3766E-06  2.7420E-04  2.8074E-02  2.4641E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9793E-01  2.3313E-01  8.6619E-01  7.9217E-01  7.7446E-01

 ETASHRINKSD(%)  8.7248E-01  9.9979E+01  9.9081E+01  5.9472E+00  1.7449E+01
 ETASHRINKVR(%)  1.7373E+00  1.0000E+02  9.9992E+01  1.1541E+01  3.1854E+01
 EBVSHRINKSD(%)  1.0648E+00  9.9981E+01  9.9089E+01  5.7151E+00  1.6825E+01
 EBVSHRINKVR(%)  2.1183E+00  1.0000E+02  9.9992E+01  1.1104E+01  3.0819E+01
 RELATIVEINF(%)  8.4198E+01  2.8479E-07  3.4042E-04  1.1418E+01  2.2558E+00
 EPSSHRINKSD(%)  3.8511E+01
 EPSSHRINKVR(%)  6.2191E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1523.9073920612391     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -788.75656549750090     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.02
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1523.907       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  1.00E-02  4.33E-01  1.40E+00  3.76E-01  1.03E+00  3.61E-02  1.00E-02  7.54E-01  7.09E-01  1.72E+00
 


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
+        1.04E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.00E+01  0.00E+00  4.93E+03
 
 TH 4
+       -2.34E+01  0.00E+00 -4.77E+02  8.88E+02
 
 TH 5
+        6.14E+01  0.00E+00 -6.85E+03 -1.68E+02  1.08E+04
 
 TH 6
+       -2.61E-01  0.00E+00  8.51E+00 -6.30E+00  1.66E+00  1.79E+02
 
 TH 7
+       -1.74E-02  0.00E+00  1.45E-02  9.87E-03  2.60E-03 -3.05E-03 -6.73E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.83E+00  0.00E+00  7.37E+01 -1.66E+01  1.55E+00  1.99E+00  7.23E-03  0.00E+00  2.78E+02
 
 TH10
+       -2.46E+00  0.00E+00 -4.91E+01  4.71E+00 -5.32E+01 -2.28E-01 -4.92E-02  0.00E+00  5.02E-01  1.84E+02
 
 TH11
+       -9.88E+00  0.00E+00 -4.24E+01 -8.52E+00  1.30E+01  1.61E+00  9.19E-03  0.00E+00  1.15E+01  3.39E+01  8.14E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.490
Stop Time:
Wed Sep 29 12:24:16 CDT 2021
