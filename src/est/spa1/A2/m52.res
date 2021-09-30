Wed Sep 29 23:26:22 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat52.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m52.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1190.75971408619        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2558E+02  1.0311E+01  8.8768E+01 -6.7267E+01  4.9535E+01  3.8117E+01 -5.8788E+01 -2.0328E+01 -4.7257E+01 -8.0400E+01
            -1.6149E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1718.41647524842        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0182E+00  1.0698E+00  1.0750E+00  1.0959E+00  1.0810E+00  8.8401E-01  1.0552E+00  8.8956E-01  1.0194E+00  9.0045E-01
             2.4245E+00
 PARAMETER:  1.1803E-01  1.6744E-01  1.7235E-01  1.9162E-01  1.7792E-01 -2.3290E-02  1.5377E-01 -1.7028E-02  1.1918E-01 -4.8602E-03
             9.8562E-01
 GRADIENT:   1.1927E+02  2.0000E+01 -8.7730E+00  3.9941E+01 -2.2451E+00 -2.3189E+01 -1.9906E-01  6.6465E+00  3.3470E+00  8.9222E+00
             1.8463E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1722.67957306073        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9960E-01  8.6244E-01  1.1701E+00  1.2350E+00  1.0309E+00  9.4785E-01  1.2158E+00  3.7404E-01  9.0397E-01  8.9787E-01
             2.3817E+00
 PARAMETER:  9.9596E-02 -4.7987E-02  2.5709E-01  3.1105E-01  1.3040E-01  4.6437E-02  2.9537E-01 -8.8340E-01 -9.6328E-04 -7.7332E-03
             9.6783E-01
 GRADIENT:   6.6364E+01  2.3878E+01  4.1754E+00  7.9009E+01 -6.1821E+00  5.0598E+00 -3.1725E+00  8.8617E-01 -6.4177E+00  4.9608E+00
            -1.7218E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1724.47668608587        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.8171E-01  9.2879E-01  6.1693E-01  1.1235E+00  7.2720E-01  9.5976E-01  1.3289E+00  1.9091E-01  9.0387E-01  5.4640E-01
             2.3533E+00
 PARAMETER:  8.1541E-02  2.6127E-02 -3.8301E-01  2.1642E-01 -2.1856E-01  5.8933E-02  3.8439E-01 -1.5560E+00 -1.0702E-03 -5.0440E-01
             9.5582E-01
 GRADIENT:   1.3218E+01  1.1567E+01 -1.2084E+01  3.9869E+01  1.9640E+01  5.7431E+00  3.3787E+00  6.5651E-01  1.2494E+00  6.0341E+00
            -5.4735E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1727.21285364646        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      394
 NPARAMETR:  1.0046E+00  8.2487E-01  6.2027E-01  1.1601E+00  6.8685E-01  9.5328E-01  1.4585E+00  2.1923E-01  8.6161E-01  3.7278E-01
             2.4734E+00
 PARAMETER:  1.0464E-01 -9.2535E-02 -3.7761E-01  2.4853E-01 -2.7564E-01  5.2156E-02  4.7741E-01 -1.4176E+00 -4.8952E-02 -8.8677E-01
             1.0056E+00
 GRADIENT:  -7.3718E-01 -3.1220E+00 -4.4593E+00 -2.4348E+01  4.6340E+00  8.6844E-01  7.1161E-01  7.2242E-01  1.4970E-01  2.0207E+00
             2.1574E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1729.18778431885        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      569
 NPARAMETR:  9.9929E-01  5.7881E-01  6.4491E-01  1.3023E+00  6.1893E-01  9.4510E-01  1.9420E+00  5.9699E-02  8.0325E-01  1.7666E-01
             2.4834E+00
 PARAMETER:  9.9285E-02 -4.4678E-01 -3.3864E-01  3.6410E-01 -3.7976E-01  4.3537E-02  7.6373E-01 -2.7184E+00 -1.1909E-01 -1.6335E+00
             1.0096E+00
 GRADIENT:  -4.8956E+00  3.7804E+00  6.1531E+00 -2.2142E+00 -9.0504E+00 -9.0610E-01  3.9123E-02  3.5221E-03 -9.5301E-01 -8.7993E-01
             6.9919E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1730.58568979965        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      748
 NPARAMETR:  9.6913E-01  3.7777E-01  7.2490E-01  1.4541E+00  6.0827E-01  9.3667E-01  2.5552E+00  1.0000E-02  7.8425E-01  7.7762E-01
             2.3859E+00
 PARAMETER:  6.8646E-02 -8.7347E-01 -2.2172E-01  4.7438E-01 -3.9713E-01  3.4571E-02  1.0381E+00 -8.1168E+00 -1.4303E-01 -1.5151E-01
             9.6957E-01
 GRADIENT:  -6.6552E+01  1.3967E+01 -1.1675E+01  6.6635E+01  1.6017E+01 -5.9234E+00 -1.0178E+00  0.0000E+00 -2.4570E+01 -1.9993E-01
             5.7779E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1734.55360396275        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      927
 NPARAMETR:  9.9156E-01  3.3588E-01  5.0633E-01  1.3931E+00  4.6468E-01  9.4219E-01  2.4561E+00  1.0000E-02  8.4029E-01  6.4904E-01
             2.4128E+00
 PARAMETER:  9.1525E-02 -9.9101E-01 -5.8057E-01  4.3150E-01 -6.6640E-01  4.0450E-02  9.9858E-01 -9.6802E+00 -7.4009E-02 -3.3226E-01
             9.8078E-01
 GRADIENT:  -1.4975E+01  1.0739E+00 -4.2684E+01  3.6313E+01  4.8127E+01 -2.8553E+00 -3.7756E+00  0.0000E+00 -1.5217E+01  4.9175E+00
             7.8694E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1741.23719912941        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1103
 NPARAMETR:  9.8921E-01  2.5169E-01  4.6548E-01  1.3814E+00  4.1576E-01  9.5101E-01  3.0154E+00  1.0000E-02  9.0033E-01  6.3196E-01
             2.1740E+00
 PARAMETER:  8.9156E-02 -1.2796E+00 -6.6469E-01  4.2311E-01 -7.7764E-01  4.9768E-02  1.2037E+00 -1.2177E+01 -4.9973E-03 -3.5892E-01
             8.7659E-01
 GRADIENT:  -5.4553E+00  4.2425E-01  2.0533E+00 -7.3428E+00 -1.5776E+00 -1.7123E-01  7.9921E-01  0.0000E+00 -1.1009E+00  2.1059E-01
             1.5465E+00

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1741.25759273681        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1195
 NPARAMETR:  9.9073E-01  2.4714E-01  4.6203E-01  1.3850E+00  4.1316E-01  9.5121E-01  3.0351E+00  1.0000E-02  9.0429E-01  6.2899E-01
             2.1728E+00
 PARAMETER:  9.0689E-02 -1.2978E+00 -6.7213E-01  4.2573E-01 -7.8392E-01  4.9984E-02  1.2102E+00 -1.2405E+01 -5.9974E-04 -3.6364E-01
             8.7603E-01
 GRADIENT:  -1.5681E+00 -1.7968E-01 -8.8453E-03 -1.0393E+00  1.1753E-01 -1.4906E-01  2.3701E-01  0.0000E+00  6.6240E-02 -6.7860E-02
             1.0862E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1195
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3141E-03  4.3731E-02 -3.0251E-04 -2.4365E-02  1.5052E-02
 SE:             2.9370E-02  1.7701E-02  1.9630E-04  2.6554E-02  1.8383E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1016E-01  1.3488E-02  1.2331E-01  3.5885E-01  4.1288E-01

 ETASHRINKSD(%)  1.6069E+00  4.0701E+01  9.9342E+01  1.1040E+01  3.8415E+01
 ETASHRINKVR(%)  3.1880E+00  6.4836E+01  9.9996E+01  2.0861E+01  6.2073E+01
 EBVSHRINKSD(%)  1.7347E+00  4.9876E+01  9.9176E+01  9.0144E+00  3.3305E+01
 EBVSHRINKVR(%)  3.4394E+00  7.4876E+01  9.9993E+01  1.7216E+01  5.5517E+01
 RELATIVEINF(%)  9.5684E+01  9.6723E+00  4.0931E-04  4.1431E+01  2.5020E+00
 EPSSHRINKSD(%)  2.8350E+01
 EPSSHRINKVR(%)  4.8662E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1741.2575927368127     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -822.31905953214005     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1741.258       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  2.47E-01  4.62E-01  1.39E+00  4.13E-01  9.51E-01  3.04E+00  1.00E-02  9.04E-01  6.29E-01  2.17E+00
 


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
+        1.21E+03
 
 TH 2
+       -3.90E+01  8.28E+02
 
 TH 3
+       -5.95E+00  4.04E+02  3.68E+03
 
 TH 4
+       -1.31E+01  1.65E+02 -3.17E+02  6.00E+02
 
 TH 5
+        4.52E+01 -9.69E+02 -4.95E+03  6.84E+01  7.43E+03
 
 TH 6
+        5.48E+00 -4.10E+00  7.75E+00 -7.99E+00  1.03E+00  2.04E+02
 
 TH 7
+        2.22E+00  5.29E+01 -5.88E+00 -4.26E+00 -1.04E+01  5.25E-01  5.91E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.26E+00 -3.08E+01  2.28E+01 -5.83E+00  6.06E+01 -1.99E+00  2.58E+00  0.00E+00  1.68E+02
 
 TH10
+       -1.13E+00  2.44E+01 -1.83E+02 -1.68E+01  1.37E+02 -3.94E-01  4.80E+00  0.00E+00 -1.25E+01  1.18E+02
 
 TH11
+       -1.64E+01 -1.29E+01 -3.18E+01 -1.05E+01  3.31E+01  3.48E+00  1.53E-01  0.00E+00  1.09E+01  2.00E+01  9.53E+01
 
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
 #CPUT: Total CPU Time in Seconds,       26.519
Stop Time:
Wed Sep 29 23:26:50 CDT 2021
