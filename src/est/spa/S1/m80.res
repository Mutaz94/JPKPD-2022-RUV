Sat Sep 25 10:06:47 CDT 2021
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
$DATA ../../../../data/spa/S1/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1723.56991466374        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.5332E-01 -1.1461E+01  1.2087E+01 -1.7637E+01 -4.0247E+01  4.1206E+01  5.2016E+00  1.8939E+00  3.8254E+01  1.0186E+01
            -2.8532E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1730.50289866838        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  1.0090E+00  1.0627E+00  1.0083E+00  9.8062E-01  1.0611E+00  8.5290E-01  9.8102E-01  9.9243E-01  8.1771E-01  9.6164E-01
             1.0172E+00
 PARAMETER:  1.0901E-01  1.6081E-01  1.0830E-01  8.0432E-02  1.5928E-01 -5.9110E-02  8.0842E-02  9.2403E-02 -1.0125E-01  6.0888E-02
             1.1701E-01
 GRADIENT:  -3.7264E+01  1.1223E+00  3.2751E+00  2.4952E+00  2.3884E+00 -2.2331E+01 -4.2178E+00 -2.4529E+00 -7.6225E-01 -1.8774E+00
            -1.3433E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1731.16923612847        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      292
 NPARAMETR:  1.0207E+00  1.0924E+00  9.6149E-01  9.6959E-01  1.0386E+00  8.8206E-01  1.0273E+00  1.0227E+00  7.8069E-01  9.2349E-01
             1.0169E+00
 PARAMETER:  1.2052E-01  1.8838E-01  6.0728E-02  6.9117E-02  1.3792E-01 -2.5499E-02  1.2696E-01  1.2249E-01 -1.4758E-01  2.0403E-02
             1.1671E-01
 GRADIENT:  -4.0668E+00  1.9116E+01  5.2012E+00  1.7022E+01 -8.3498E+00 -8.0167E+00 -7.8251E-01 -5.6218E-01 -5.6487E+00 -1.1806E+00
            -1.2096E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1731.59220256053        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  1.0237E+00  1.1689E+00  8.2020E-01  9.0392E-01  1.0134E+00  9.0301E-01  9.8104E-01  8.4857E-01  8.3600E-01  9.0632E-01
             1.0150E+00
 PARAMETER:  1.2340E-01  2.5610E-01 -9.8209E-02 -1.0191E-03  1.1329E-01 -2.0194E-03  8.0855E-02 -6.4208E-02 -7.9123E-02  1.6327E-03
             1.1490E-01
 GRADIENT:   5.8918E-01 -1.5443E-01 -5.1548E-01  1.0431E+00 -6.4648E-02  7.0953E-01  5.2561E-01  1.6453E-01  2.6896E-01  4.1867E-01
            -8.7958E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1731.71232285674        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  1.0253E+00  1.3623E+00  6.0989E-01  7.7083E-01  9.9394E-01  9.0322E-01  8.8522E-01  5.7917E-01  9.0723E-01  8.6367E-01
             1.0150E+00
 PARAMETER:  1.2501E-01  4.0916E-01 -3.9448E-01 -1.6029E-01  9.3923E-02 -1.7925E-03 -2.1919E-02 -4.4616E-01  2.6404E-03 -4.6562E-02
             1.1485E-01
 GRADIENT:  -5.7535E-02  1.3488E+00 -3.9742E-01  1.1591E+00 -1.7500E+00 -1.2690E-01  2.2204E-01  3.3239E-01  5.3395E-03  6.1241E-01
             1.9788E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1731.75437192572        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  1.0257E+00  1.4660E+00  5.3036E-01  6.9987E-01  1.0081E+00  9.0374E-01  8.3919E-01  4.1675E-01  9.6216E-01  8.6062E-01
             1.0154E+00
 PARAMETER:  1.2541E-01  4.8252E-01 -5.3419E-01 -2.5686E-01  1.0804E-01 -1.2121E-03 -7.5316E-02 -7.7526E-01  6.1423E-02 -5.0099E-02
             1.1524E-01
 GRADIENT:   2.6577E-02 -3.6472E-01 -4.2906E-01 -4.0086E-01 -2.4701E-01 -1.2683E-01  1.7348E-01  1.9044E-01  1.3490E-01  3.0890E-01
             1.6449E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1731.80983330553        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  1.0259E+00  1.4467E+00  5.1105E-01  7.0910E-01  9.8111E-01  9.0457E-01  8.4884E-01  1.9259E-01  9.4946E-01  8.4056E-01
             1.0144E+00
 PARAMETER:  1.2553E-01  4.6932E-01 -5.7128E-01 -2.4376E-01  8.0929E-02 -2.9371E-04 -6.3887E-02 -1.5472E+00  4.8142E-02 -7.3681E-02
             1.1433E-01
 GRADIENT:   7.9304E-02 -8.3810E-01 -1.1546E+00  1.0113E+00  1.2377E+00  1.4297E-01  9.1382E-02  4.2484E-02  1.4920E-01  2.2163E-01
             1.3442E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1731.82553316100        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  1.0258E+00  1.4567E+00  5.0323E-01  7.0191E-01  9.8157E-01  9.0389E-01  8.4338E-01  3.8820E-02  9.5584E-01  8.3950E-01
             1.0144E+00
 PARAMETER:  1.2546E-01  4.7616E-01 -5.8671E-01 -2.5395E-01  8.1395E-02 -1.0448E-03 -7.0339E-02 -3.1488E+00  5.4837E-02 -7.4946E-02
             1.1434E-01
 GRADIENT:  -1.3533E-01 -3.9617E-01  6.8706E-02 -2.4040E-01 -7.6493E-02 -1.6435E-01 -3.1293E-01  1.2362E-03  2.0661E-02 -6.7080E-02
            -2.6224E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1731.82653206970        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1353
 NPARAMETR:  1.0258E+00  1.4552E+00  5.0392E-01  7.0307E-01  9.8105E-01  9.0427E-01  8.4567E-01  1.0000E-02  9.5441E-01  8.3978E-01
             1.0144E+00
 PARAMETER:  1.2552E-01  4.7513E-01 -5.8534E-01 -2.5230E-01  8.0869E-02 -6.2228E-04 -6.7628E-02 -4.6170E+00  5.3334E-02 -7.4613E-02
             1.1427E-01
 GRADIENT:  -7.8264E-04 -5.3014E-03 -1.8911E-02  1.3149E-02  9.1542E-03  1.1129E-04  1.0729E-02  0.0000E+00  3.3415E-03  7.6341E-03
             2.7836E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1731.82653206970        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1375
 NPARAMETR:  1.0258E+00  1.4552E+00  5.0392E-01  7.0307E-01  9.8105E-01  9.0427E-01  8.4567E-01  1.0000E-02  9.5441E-01  8.3978E-01
             1.0144E+00
 PARAMETER:  1.2552E-01  4.7513E-01 -5.8534E-01 -2.5230E-01  8.0869E-02 -6.2228E-04 -6.7628E-02 -4.6170E+00  5.3334E-02 -7.4613E-02
             1.1427E-01
 GRADIENT:  -7.8264E-04 -5.3014E-03 -1.8911E-02  1.3149E-02  9.1542E-03  1.1129E-04  1.0729E-02  0.0000E+00  3.3415E-03  7.6341E-03
             2.7836E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1375
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.4432E-05 -1.6321E-02 -3.4372E-04  1.2857E-02 -2.4848E-02
 SE:             2.9826E-02  2.4295E-02  1.4267E-04  2.2717E-02  2.1699E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9801E-01  5.0172E-01  1.5990E-02  5.7142E-01  2.5217E-01

 ETASHRINKSD(%)  8.0735E-02  1.8608E+01  9.9522E+01  2.3896E+01  2.7305E+01
 ETASHRINKVR(%)  1.6140E-01  3.3753E+01  9.9998E+01  4.2082E+01  4.7154E+01
 EBVSHRINKSD(%)  5.1359E-01  1.8533E+01  9.9578E+01  2.4995E+01  2.6421E+01
 EBVSHRINKVR(%)  1.0245E+00  3.3632E+01  9.9998E+01  4.3743E+01  4.5861E+01
 RELATIVEINF(%)  9.8923E+01  3.2046E+00  1.2032E-04  2.4125E+00  6.0870E+00
 EPSSHRINKSD(%)  4.3792E+01
 EPSSHRINKVR(%)  6.8406E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1731.8265320697010     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -996.67570550596281     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1731.827       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.46E+00  5.04E-01  7.03E-01  9.81E-01  9.04E-01  8.46E-01  1.00E-02  9.54E-01  8.40E-01  1.01E+00
 


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
+        1.28E+03
 
 TH 2
+       -9.04E+00  5.01E+02
 
 TH 3
+        1.23E+01  2.47E+02  8.62E+02
 
 TH 4
+       -2.51E+01  4.10E+02 -6.19E+02  1.40E+03
 
 TH 5
+       -5.90E+00 -3.05E+02 -7.27E+02  4.94E+02  8.86E+02
 
 TH 6
+        1.30E-01 -1.19E+00  3.31E+00 -3.83E+00 -6.91E-01  2.37E+02
 
 TH 7
+        8.79E-01  1.94E+01 -3.75E+01 -1.18E+01  3.42E-01  1.50E+00  1.24E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.31E+00 -2.08E+01 -4.85E+01  6.01E+01  1.25E-01 -1.13E+00  2.03E+01  0.00E+00  7.99E+01
 
 TH10
+       -1.57E+00 -1.55E+01 -5.17E+01 -1.32E+01 -7.14E+01 -2.01E-01  2.07E+01  0.00E+00  1.77E+01  8.81E+01
 
 TH11
+       -7.21E+00 -1.72E+01 -3.42E+01  4.49E+00 -3.25E+00  2.24E+00  1.09E+01  0.00E+00  1.23E+01  2.08E+01  2.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.069
Stop Time:
Sat Sep 25 10:07:11 CDT 2021
