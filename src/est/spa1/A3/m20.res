Thu Sep 30 00:03:35 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat20.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   259.037196290924        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9202E+02  6.6773E+01  1.2785E+02  7.4230E+01  2.4181E+02  4.7370E+01 -5.8428E+01 -1.9428E+02 -5.0999E+01 -1.2787E+02
            -4.2007E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1431.94394684858        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.9664E-01  1.0336E+00  9.1587E-01  1.0617E+00  8.8299E-01  7.6536E-01  9.6539E-01  1.1895E+00  9.4932E-01  9.2709E-01
             5.0173E+00
 PARAMETER:  9.6630E-02  1.3305E-01  1.2120E-02  1.5991E-01 -2.4445E-02 -1.6740E-01  6.4779E-02  2.7356E-01  4.7996E-02  2.4299E-02
             1.7129E+00
 GRADIENT:  -2.3374E+01 -2.5264E+00 -9.5403E+00  3.6845E+00 -1.5406E+01 -1.8338E+01  1.0215E+01  9.6764E+00  1.6770E+01  2.2958E+01
             2.9246E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1450.97773263737        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.8778E-01  1.3985E+00  3.7656E-01  8.2548E-01  6.4239E-01  8.2303E-01  6.8379E-01  1.5130E+00  1.3072E+00  4.3090E-01
             4.5366E+00
 PARAMETER:  8.7702E-02  4.3541E-01 -8.7668E-01 -9.1787E-02 -3.4256E-01 -9.4762E-02 -2.8010E-01  5.1411E-01  3.6787E-01 -7.4187E-01
             1.6122E+00
 GRADIENT:  -5.8509E+01  1.5267E+02  4.1572E+01  3.4135E+01 -1.3163E+02 -1.7702E+01  3.6255E+00  4.6124E+00  1.6746E+01  6.5813E+00
             2.3667E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1503.74721257580        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.6322E-01  7.1353E-01  2.6447E-01  1.0889E+00  3.6678E-01  8.9110E-01  2.7358E-01  1.0001E+00  1.2071E+00  4.5653E-01
             3.0870E+00
 PARAMETER:  6.2529E-02 -2.3753E-01 -1.2300E+00  1.8516E-01 -9.0300E-01 -1.5301E-02 -1.1962E+00  1.0006E-01  2.8824E-01 -6.8410E-01
             1.2272E+00
 GRADIENT:  -2.8127E+01  1.1781E+02  3.5320E+01  8.4822E+01 -4.1630E+01 -5.7247E+00 -1.9362E+00 -8.1051E+00 -1.3770E+01 -7.7817E-01
            -2.1311E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1504.05129033950        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.6004E-01  6.1840E-01  2.3763E-01  1.1259E+00  3.2134E-01  8.9970E-01  2.5760E-01  9.2502E-01  1.1557E+00  4.2743E-01
             2.9321E+00
 PARAMETER:  5.9222E-02 -3.8062E-01 -1.3370E+00  2.1858E-01 -1.0353E+00 -5.6961E-03 -1.2563E+00  2.2060E-02  2.4468E-01 -7.4996E-01
             1.1757E+00
 GRADIENT:  -2.9603E+01  1.4473E+02  4.4294E+01  1.2880E+02 -6.1035E+01 -4.3117E+00 -2.4482E+00 -1.5132E+01 -3.9794E+01 -7.8279E+00
            -5.5292E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1512.41290453030        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      456
 NPARAMETR:  9.9420E-01  5.4024E-01  2.7119E-01  1.1328E+00  3.2915E-01  9.1136E-01  1.9650E-01  1.1207E+00  1.0962E+00  2.4857E-01
             3.0057E+00
 PARAMETER:  9.4182E-02 -5.1575E-01 -1.2049E+00  2.2467E-01 -1.0113E+00  7.1791E-03 -1.5271E+00  2.1399E-01  1.9188E-01 -1.2920E+00
             1.2005E+00
 GRADIENT:   2.5756E+01  4.5587E+01  2.3397E+01  5.4730E+01 -2.4311E+01  3.1160E+00 -1.1196E+00 -1.1662E+01 -3.2064E+01 -3.6614E+00
            -4.5146E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1525.17313724946        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      631
 NPARAMETR:  9.7803E-01  3.8286E-01  1.9983E-01  1.0584E+00  2.5021E-01  8.8872E-01  2.2214E-01  1.3207E+00  1.3705E+00  1.9192E-02
             2.9632E+00
 PARAMETER:  7.7783E-02 -8.6009E-01 -1.5103E+00  1.5676E-01 -1.2854E+00 -1.7977E-02 -1.4045E+00  3.7816E-01  4.1516E-01 -3.8533E+00
             1.1863E+00
 GRADIENT:  -3.3251E+00 -4.7692E+00  9.6631E+00 -6.5091E+00 -6.4664E+00 -1.7723E+00 -8.4028E-01 -5.2588E-01  2.2687E+00 -3.0851E-02
             4.8419E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1526.83349244636        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      808
 NPARAMETR:  9.7899E-01  4.1361E-01  1.9176E-01  1.0524E+00  2.5169E-01  9.0118E-01  7.1347E-01  1.3001E+00  1.3426E+00  2.5539E-02
             2.9182E+00
 PARAMETER:  7.8768E-02 -7.8284E-01 -1.5515E+00  1.5108E-01 -1.2795E+00 -4.0520E-03 -2.3762E-01  3.6241E-01  3.9461E-01 -3.5675E+00
             1.1710E+00
 GRADIENT:  -8.5704E-01  6.8202E-01 -2.1863E+00 -4.9256E-01  1.9222E+00  1.6225E+00  4.8122E-01  1.9609E-01 -2.3640E+00 -2.7432E-02
            -2.1369E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1526.87703370421        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  9.7924E-01  4.1044E-01  1.9242E-01  1.0544E+00  2.5131E-01  8.9655E-01  6.6499E-01  1.3015E+00  1.3567E+00  2.4207E-02
             2.9263E+00
 PARAMETER:  7.9026E-02 -7.9052E-01 -1.5481E+00  1.5293E-01 -1.2811E+00 -9.2064E-03 -3.0798E-01  3.6349E-01  4.0507E-01 -3.6211E+00
             1.1737E+00
 GRADIENT:  -2.2797E-01  3.0581E-01  4.8246E-02  4.2683E-01 -4.3433E-01  1.1333E-02 -2.8083E-02 -3.2309E-02 -1.2916E-01 -2.6291E-02
            -1.4274E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      900
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1394E-03 -1.5758E-02  7.5202E-03 -5.3366E-03  1.1442E-03
 SE:             2.8916E-02  9.7368E-03  2.2925E-02  2.7136E-02  9.7963E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4102E-01  1.0557E-01  7.4289E-01  8.4409E-01  2.4283E-01

 ETASHRINKSD(%)  3.1281E+00  6.7381E+01  2.3197E+01  9.0922E+00  9.6718E+01
 ETASHRINKVR(%)  6.1584E+00  8.9360E+01  4.1013E+01  1.7358E+01  9.9892E+01
 EBVSHRINKSD(%)  3.1676E+00  6.6910E+01  2.3331E+01  7.8398E+00  9.7078E+01
 EBVSHRINKVR(%)  6.2349E+00  8.9050E+01  4.1219E+01  1.5065E+01  9.9915E+01
 RELATIVEINF(%)  9.3589E+01  1.6623E+00  8.4567E+00  6.0811E+01  4.8303E-03
 EPSSHRINKSD(%)  2.6959E+01
 EPSSHRINKVR(%)  4.6650E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1526.8770337042058     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -607.93850049953312     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.32
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1526.877       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  4.10E-01  1.92E-01  1.05E+00  2.51E-01  8.97E-01  6.65E-01  1.30E+00  1.36E+00  2.42E-02  2.93E+00
 


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
+       -6.04E+01  1.73E+03
 
 TH 3
+       -9.73E+01  2.61E+03  1.03E+04
 
 TH 4
+       -2.21E+01  2.10E+02 -1.17E+02  4.59E+02
 
 TH 5
+        2.48E+02 -6.00E+03 -1.41E+04 -7.02E+02  2.87E+04
 
 TH 6
+        3.55E+00 -1.98E+01  1.82E+01 -8.61E+00  3.54E+01  2.16E+02
 
 TH 7
+       -3.17E-01 -3.11E+01 -3.18E+01 -2.58E+00  1.20E+02  1.07E-01  1.05E+01
 
 TH 8
+       -5.03E-01 -1.13E+01 -3.06E+01  1.69E+00  4.08E+01  2.28E+00  4.07E+00  4.38E+01
 
 TH 9
+        1.22E+01 -5.76E+01  6.51E+01 -9.11E+00  2.39E+02  2.45E+00  6.67E+00 -9.31E-01  7.05E+01
 
 TH10
+       -3.66E-01 -7.03E+00 -3.83E+00 -2.33E-01  3.02E+01 -1.54E-01  1.59E+00  2.42E+00  1.46E+00 -4.51E+01
 
 TH11
+       -2.18E+01 -1.84E+01 -3.02E+01 -3.67E+00  5.42E+01  2.28E+00  3.69E+00  9.19E+00  6.19E+00  1.32E+00  5.26E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.524
Stop Time:
Thu Sep 30 00:04:04 CDT 2021
