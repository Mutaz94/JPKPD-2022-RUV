Thu Sep 30 01:57:20 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2163.01644399526        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6225E+02 -2.7090E+00 -1.0182E+01  2.7751E+01 -2.8382E+01  5.5950E+01  2.5418E+00  1.5883E+01  2.0422E+01  1.9597E+01
            -3.3643E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2173.06712317404        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0467E+00  1.1089E+00  1.1545E+00  9.9750E-01  1.1410E+00  9.8395E-01  1.0186E+00  8.9440E-01  9.1893E-01  9.2754E-01
             1.0358E+00
 PARAMETER:  1.4560E-01  2.0335E-01  2.4365E-01  9.7502E-02  2.3187E-01  8.3818E-02  1.1842E-01 -1.1597E-02  1.5451E-02  2.4779E-02
             1.3514E-01
 GRADIENT:  -2.1004E+00  1.7304E+01  6.9040E-01  2.3626E+01  2.5428E+00 -1.6387E+00 -5.1412E+00  4.1001E+00 -3.7471E+00 -1.6134E+01
            -1.7078E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2175.19416409416        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0451E+00  1.0352E+00  1.3185E+00  1.0404E+00  1.2181E+00  9.8683E-01  1.0824E+00  5.9243E-01  9.3579E-01  1.1587E+00
             1.0415E+00
 PARAMETER:  1.4407E-01  1.3461E-01  3.7646E-01  1.3958E-01  2.9733E-01  8.6741E-02  1.7916E-01 -4.2352E-01  3.3637E-02  2.4734E-01
             1.4070E-01
 GRADIENT:  -1.4209E+00  3.9128E+00 -1.0206E-02  1.1301E+01  1.0460E+01  2.0682E-01  7.6933E-01 -8.4211E-01  4.0458E+00  1.0236E+00
            -1.3042E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2175.91217307982        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      544
 NPARAMETR:  1.0477E+00  1.1454E+00  1.0222E+00  9.5889E-01  1.1137E+00  9.8791E-01  1.0729E+00  3.5497E-01  9.2108E-01  1.0317E+00
             1.0504E+00
 PARAMETER:  1.4661E-01  2.3573E-01  1.2200E-01  5.8025E-02  2.0767E-01  8.7831E-02  1.7033E-01 -9.3571E-01  1.7786E-02  1.3119E-01
             1.4921E-01
 GRADIENT:  -8.6295E-01  7.1040E+00 -2.1540E-01  8.0763E+00 -4.4874E+00 -2.2636E-01  3.1489E-02  4.1412E-01 -3.0115E-01  8.8504E-01
            -1.0201E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2176.10234027046        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      720
 NPARAMETR:  1.0493E+00  1.2505E+00  9.4545E-01  8.8474E-01  1.1384E+00  9.8971E-01  1.0082E+00  1.1897E-01  9.7006E-01  1.0338E+00
             1.0510E+00
 PARAMETER:  1.4811E-01  3.2355E-01  4.3908E-02 -2.2457E-02  2.2960E-01  8.9654E-02  1.0821E-01 -2.0289E+00  6.9607E-02  1.3320E-01
             1.4977E-01
 GRADIENT:   4.8956E-01 -5.3634E-01 -3.5361E-01 -2.8856E-01  1.3890E-01  7.4876E-02  1.0852E-01  4.3759E-02  1.2255E-01  1.8869E-01
             1.5190E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2176.10276427395        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  1.0493E+00  1.2644E+00  9.3808E-01  8.7594E-01  1.1424E+00  9.8972E-01  1.0001E+00  1.0159E-01  9.7637E-01  1.0346E+00
             1.0509E+00
 PARAMETER:  1.4809E-01  3.3462E-01  3.6080E-02 -3.2456E-02  2.3317E-01  8.9670E-02  1.0010E-01 -2.1868E+00  7.6091E-02  1.3398E-01
             1.4962E-01
 GRADIENT:   2.3202E-01 -3.6095E-01 -2.4752E-01 -2.7535E-01  1.5206E-01  3.4277E-02  5.5536E-02  3.2671E-02  7.9117E-02  9.4964E-02
             8.2630E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2176.10334962974        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  1.0492E+00  1.2770E+00  9.3118E-01  8.6785E-01  1.1462E+00  9.8972E-01  9.9291E-01  7.5552E-02  9.8216E-01  1.0354E+00
             1.0507E+00
 PARAMETER:  1.4806E-01  3.4448E-01  2.8696E-02 -4.1735E-02  2.3644E-01  8.9670E-02  9.2882E-02 -2.4829E+00  8.1996E-02  1.3477E-01
             1.4949E-01
 GRADIENT:  -4.7752E-02 -4.2741E-01 -1.7878E-01 -4.4310E-01  2.6530E-01 -7.8210E-03  4.7861E-03  1.8418E-02  3.5874E-02  1.0238E-02
             3.2937E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2176.12286348279        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  1.0491E+00  1.2415E+00  9.4786E-01  8.9159E-01  1.1333E+00  9.8958E-01  1.0131E+00  1.0000E-02  9.6523E-01  1.0308E+00
             1.0510E+00
 PARAMETER:  1.4797E-01  3.1636E-01  4.6455E-02 -1.4746E-02  2.2511E-01  8.9522E-02  1.1297E-01 -6.4480E+00  6.4615E-02  1.3036E-01
             1.4979E-01
 GRADIENT:   2.9930E-01  1.4615E+00  4.0687E-01  1.2536E+00 -6.7802E-01  4.3461E-02 -1.3578E-02  0.0000E+00 -9.1768E-02 -6.6283E-02
            -1.2256E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2176.12609147103        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  1.0511E+00  1.2421E+00  9.4694E-01  8.9066E-01  1.1338E+00  9.9000E-01  1.0128E+00  1.0000E-02  9.6603E-01  1.0311E+00
             1.0511E+00
 PARAMETER:  1.4979E-01  3.1683E-01  4.5481E-02 -1.5795E-02  2.2559E-01  8.9946E-02  1.1271E-01 -6.3772E+00  6.5435E-02  1.3064E-01
             1.4984E-01
 GRADIENT:   4.4707E+00  5.1218E-01  9.1066E-02  6.9825E-01 -1.0159E-01  2.1291E-01  3.1309E-02  0.0000E+00  1.1884E-02 -6.8335E-03
            -3.0415E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1420
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6823E-05 -1.1957E-02 -3.4954E-04  6.3601E-03 -2.5963E-02
 SE:             2.9851E-02  2.1310E-02  1.4160E-04  2.3592E-02  2.3085E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9955E-01  5.7473E-01  1.3568E-02  7.8748E-01  2.6074E-01

 ETASHRINKSD(%)  1.0000E-10  2.8609E+01  9.9526E+01  2.0964E+01  2.2661E+01
 ETASHRINKVR(%)  1.0000E-10  4.9034E+01  9.9998E+01  3.7534E+01  4.0187E+01
 EBVSHRINKSD(%)  3.7550E-01  2.7633E+01  9.9567E+01  2.2381E+01  2.0482E+01
 EBVSHRINKVR(%)  7.4958E-01  4.7630E+01  9.9998E+01  3.9753E+01  3.6769E+01
 RELATIVEINF(%)  9.8666E+01  2.1423E+00  2.0024E-04  2.5465E+00  9.9317E+00
 EPSSHRINKSD(%)  3.2146E+01
 EPSSHRINKVR(%)  5.3958E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2176.1260914710324     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1257.1875582663597     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.69
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.19
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2176.126       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.24E+00  9.47E-01  8.91E-01  1.13E+00  9.90E-01  1.01E+00  1.00E-02  9.66E-01  1.03E+00  1.05E+00
 


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
+        1.02E+03
 
 TH 2
+       -7.29E+00  3.83E+02
 
 TH 3
+        6.76E+00  1.07E+02  2.46E+02
 
 TH 4
+       -3.75E+00  3.99E+02 -1.86E+02  9.06E+02
 
 TH 5
+        1.32E+00 -1.78E+02 -2.89E+02  2.06E+02  5.13E+02
 
 TH 6
+       -6.20E-02 -1.54E+00  1.90E+00 -1.72E+00 -3.15E-01  2.00E+02
 
 TH 7
+        1.24E+00  1.89E+01  6.34E+00 -1.03E+01 -1.25E+01 -1.20E-01  5.81E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.31E-01 -1.62E+01 -2.33E+01  2.70E+01  1.31E+01 -6.17E-01  3.32E+01  0.00E+00  8.17E+01
 
 TH10
+        1.33E+00 -8.96E+00 -3.08E+01  2.51E-01 -5.09E+01  2.27E-01  1.49E+00  0.00E+00  8.50E+00  8.10E+01
 
 TH11
+       -8.51E+00 -1.79E+01 -3.51E+01 -3.75E+00  4.96E+00  2.26E+00  7.01E+00  0.00E+00  7.03E+00  2.17E+01  3.78E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.950
Stop Time:
Thu Sep 30 01:57:50 CDT 2021
