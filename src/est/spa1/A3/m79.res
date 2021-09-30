Thu Sep 30 00:34:10 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat79.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   711.718451353895        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1368E+02  4.3232E+01  3.3096E+02 -8.7455E+01  2.5010E+02  6.6606E+01 -5.9891E+01 -4.8090E+02 -9.6287E+01 -1.4624E+02
            -4.8143E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1440.89998068021        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0425E+00  1.0484E+00  8.8985E-01  1.1663E+00  9.5770E-01  6.6575E-01  9.1831E-01  1.0267E+00  8.4561E-01  9.4075E-01
             5.3276E+00
 PARAMETER:  1.4159E-01  1.4730E-01 -1.6706E-02  2.5383E-01  5.6783E-02 -3.0684E-01  1.4779E-02  1.2634E-01 -6.7700E-02  3.8926E-02
             1.7729E+00
 GRADIENT:  -2.8294E+01  1.4715E+01 -2.0369E+01  5.0393E+01 -4.3814E+00 -5.0378E+01  9.0459E+00  7.8243E+00  1.6502E+01  1.8988E+01
             2.8386E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1470.93626471600        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0276E+00  8.0310E-01  3.5223E-01  1.2058E+00  4.4922E-01  7.6750E-01  5.9665E-01  1.8075E-01  1.0715E+00  2.0616E-01
             4.8259E+00
 PARAMETER:  1.2727E-01 -1.1927E-01 -9.4348E-01  2.8714E-01 -7.0025E-01 -1.6462E-01 -4.1643E-01 -1.6106E+00  1.6902E-01 -1.4791E+00
             1.6740E+00
 GRADIENT:  -6.1542E+01  7.0492E+01  2.6477E+01  8.5304E+01 -7.4613E+01 -2.5123E+01  1.7276E+00  9.6279E-02  2.0536E+01  8.7807E-01
             2.2709E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1505.31567652120        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.0101E+00  6.6624E-01  3.0116E-01  1.1446E+00  3.9950E-01  8.7794E-01  4.0073E-01  2.1899E-01  1.0952E+00  2.8718E-01
             3.6985E+00
 PARAMETER:  1.1010E-01 -3.0610E-01 -1.1001E+00  2.3502E-01 -8.1754E-01 -3.0175E-02 -8.1447E-01 -1.4187E+00  1.9089E-01 -1.1476E+00
             1.4079E+00
 GRADIENT:  -5.7190E+00 -3.3253E+00 -1.7341E+01  1.1690E+01  5.3065E+01  1.1157E+01 -1.1950E+00 -1.3079E+00  1.1972E+01 -3.6965E+00
             1.8671E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1516.32694432436        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      416
 NPARAMETR:  1.0454E+00  4.5677E-01  2.6826E-01  1.2885E+00  2.9341E-01  8.2682E-01  5.4368E-01  1.6494E+00  1.0880E+00  2.6981E-01
             3.3230E+00
 PARAMETER:  1.4437E-01 -6.8358E-01 -1.2158E+00  3.5351E-01 -1.1262E+00 -9.0166E-02 -5.0940E-01  6.0043E-01  1.8432E-01 -1.2101E+00
             1.3009E+00
 GRADIENT:   5.0604E+01  3.7222E+01  5.5653E+01  1.3533E+02 -8.0990E+01 -1.3734E+01  3.1129E+00  2.1420E+01 -1.7729E+01  5.0159E+00
             2.6334E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1544.60643603507        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      592
 NPARAMETR:  1.0074E+00  3.7834E-01  1.5758E-01  1.0979E+00  2.2487E-01  8.5808E-01  1.4448E-01  1.5743E+00  1.4104E+00  1.9310E-01
             2.9789E+00
 PARAMETER:  1.0738E-01 -8.7196E-01 -1.7478E+00  1.9343E-01 -1.3922E+00 -5.3062E-02 -1.8346E+00  5.5380E-01  4.4384E-01 -1.5445E+00
             1.1916E+00
 GRADIENT:  -1.7408E+01  5.4335E+00  2.6563E-01  2.3820E+00  1.1549E+01 -1.4785E+00  2.4461E-01 -3.6783E+00  1.3525E+00  5.9107E-01
            -3.6849E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1545.28600864075        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      771
 NPARAMETR:  1.0131E+00  3.5363E-01  1.5057E-01  1.0856E+00  2.1454E-01  8.6037E-01  6.4989E-02  1.6045E+00  1.4419E+00  1.5781E-01
             2.9832E+00
 PARAMETER:  1.1304E-01 -9.3949E-01 -1.7933E+00  1.8211E-01 -1.4393E+00 -5.0388E-02 -2.6335E+00  5.7282E-01  4.6596E-01 -1.7464E+00
             1.1930E+00
 GRADIENT:  -6.7384E-01  6.3896E+00  5.2291E+00  2.0606E+00 -1.6300E+01 -1.5240E-01  3.9671E-02 -4.9600E-01 -1.4140E+00 -1.9186E-01
            -2.0146E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1545.35092027855        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      949
 NPARAMETR:  1.0134E+00  3.5258E-01  1.4979E-01  1.0837E+00  2.1476E-01  8.6047E-01  1.0023E-02  1.6043E+00  1.4455E+00  1.6971E-01
             2.9849E+00
 PARAMETER:  1.1330E-01 -9.4247E-01 -1.7985E+00  1.8042E-01 -1.4382E+00 -5.0276E-02 -4.5029E+00  5.7271E-01  4.6848E-01 -1.6737E+00
             1.1936E+00
 GRADIENT:   4.9936E-01  3.4338E-01  8.5973E-02 -9.6013E-02  4.0371E-01 -6.5641E-03  7.2366E-04 -5.7622E-03 -3.9345E-02 -1.1167E-02
            -2.6745E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1545.35337189628        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1126
 NPARAMETR:  1.0132E+00  3.5090E-01  1.4924E-01  1.0828E+00  2.1407E-01  8.6036E-01  1.0000E-02  1.6045E+00  1.4487E+00  1.7182E-01
             2.9846E+00
 PARAMETER:  1.1312E-01 -9.4724E-01 -1.8022E+00  1.7957E-01 -1.4415E+00 -5.0399E-02 -1.1177E+01  5.7284E-01  4.7066E-01 -1.6613E+00
             1.1935E+00
 GRADIENT:   9.5336E-02 -1.3703E-01  1.1451E-01 -9.1126E-02 -6.9416E-02 -3.1037E-02  0.0000E+00  7.8525E-02 -3.3837E-02 -8.0257E-03
            -4.9972E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1545.35337189628        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  1.0132E+00  3.5090E-01  1.4924E-01  1.0828E+00  2.1407E-01  8.6036E-01  1.0000E-02  1.6045E+00  1.4487E+00  1.7182E-01
             2.9846E+00
 PARAMETER:  1.1312E-01 -9.4724E-01 -1.8022E+00  1.7957E-01 -1.4415E+00 -5.0399E-02 -1.1177E+01  5.7284E-01  4.7066E-01 -1.6613E+00
             1.1935E+00
 GRADIENT:   9.5336E-02 -1.3703E-01  1.1451E-01 -9.1126E-02 -6.9416E-02 -3.1037E-02  0.0000E+00  7.8525E-02 -3.3837E-02 -8.0257E-03
            -4.9972E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1148
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4873E-03 -3.8548E-04  1.3394E-02 -8.0886E-03  1.5252E-02
 SE:             2.8752E-02  1.4030E-04  2.4682E-02  2.6698E-02  6.7454E-03
 N:                     100         100         100         100         100

 P VAL.:         9.0346E-01  6.0046E-03  5.8738E-01  7.6192E-01  2.3753E-02

 ETASHRINKSD(%)  3.6778E+00  9.9530E+01  1.7311E+01  1.0557E+01  7.7402E+01
 ETASHRINKVR(%)  7.2204E+00  9.9998E+01  3.1625E+01  2.0000E+01  9.4893E+01
 EBVSHRINKSD(%)  3.8322E+00  9.9492E+01  1.6286E+01  9.2671E+00  7.9183E+01
 EBVSHRINKVR(%)  7.5176E+00  9.9997E+01  2.9920E+01  1.7675E+01  9.5666E+01
 RELATIVEINF(%)  9.2303E+01  5.5536E-04  1.7734E+01  6.2177E+01  4.2364E-01
 EPSSHRINKSD(%)  2.6812E+01
 EPSSHRINKVR(%)  4.6436E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1545.3533718962756     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -626.41483869160288     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1545.353       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  3.51E-01  1.49E-01  1.08E+00  2.14E-01  8.60E-01  1.00E-02  1.60E+00  1.45E+00  1.72E-01  2.98E+00
 


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
+        1.36E+03
 
 TH 2
+       -4.04E+01  2.59E+03
 
 TH 3
+       -1.52E+02  2.82E+03  1.34E+04
 
 TH 4
+       -5.64E+00  1.42E+02 -5.91E+01  3.62E+02
 
 TH 5
+        2.27E+02 -8.37E+03 -1.74E+04 -8.37E+02  4.01E+04
 
 TH 6
+        4.69E+00 -1.98E+01  1.62E+00 -7.14E+00  5.11E+01  2.28E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        2.75E+00 -1.49E+01 -2.28E+00 -4.94E-01  4.46E+01  3.58E+00  0.00E+00  4.04E+01
 
 TH 9
+        1.62E+01 -7.55E+01  5.74E+01 -1.48E+01  4.31E+02  2.04E+00  0.00E+00  3.64E-01  5.67E+01
 
 TH10
+       -2.75E+00 -1.07E+02 -2.57E+01  2.11E+00  4.07E+02  5.64E-01  0.00E+00  1.04E+01  8.73E+00  3.35E+01
 
 TH11
+       -2.68E+01 -4.16E+01 -2.72E+01  3.67E-01  1.17E+02  2.40E+00  0.00E+00  4.48E+00  6.88E+00  8.63E+00  5.20E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.284
Stop Time:
Thu Sep 30 00:34:44 CDT 2021
