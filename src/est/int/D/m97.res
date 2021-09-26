Sat Sep 25 06:39:02 CDT 2021
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
$DATA ../../../../data/int/D/dat97.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E22.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   52357.4406645423        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8094E+02  7.9739E+02  4.9475E+02  6.3820E+02  4.1260E+02 -4.9444E+03 -2.0534E+03 -1.6711E+03 -2.8669E+03 -5.9484E+02
            -1.0076E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -472.879492118737        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1856E+00  1.5059E+00  8.9664E-01  1.5166E+00  1.0652E+00  3.8482E+00  2.5035E+00  1.0016E+00  1.4592E+00  1.1837E+00
             1.3492E+01
 PARAMETER:  2.7022E-01  5.0938E-01 -9.1048E-03  5.1650E-01  1.6318E-01  1.4476E+00  1.0177E+00  1.0159E-01  4.7789E-01  2.6865E-01
             2.7021E+00
 GRADIENT:  -8.4956E+00 -6.5267E+00 -5.1715E+00  1.5198E+02 -1.0306E+01  5.7357E+01 -2.0605E+02 -5.9394E+00 -5.4327E+01  1.5903E+01
             6.9996E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -550.284660531634        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  8.0159E-01  3.9246E+00  3.8649E+01  2.2942E+00  5.3776E+00  4.2769E+00  3.3702E+00  5.8440E-01  1.0682E+01  2.0601E+00
             1.2191E+01
 PARAMETER: -1.2115E-01  1.4673E+00  3.7545E+00  9.3040E-01  1.7822E+00  1.5532E+00  1.3150E+00 -4.3717E-01  2.4686E+00  8.2274E-01
             2.6007E+00
 GRADIENT:  -2.8897E+01  4.9706E+01 -1.5181E+01  3.3595E+01  2.9005E+01  1.1173E+02  3.1437E+01 -1.9537E-01  1.2003E+01  2.0911E+01
             2.1663E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -637.533476629946        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  8.5540E-01  2.0714E+00  1.1646E+01  2.7470E+00  2.4206E+00  1.7156E+00  2.3235E+00  5.1668E+00  1.2884E+01  5.9053E-01
             1.3310E+01
 PARAMETER: -5.6183E-02  8.2822E-01  2.5550E+00  1.1105E+00  9.8402E-01  6.3975E-01  9.4306E-01  1.7422E+00  2.6560E+00 -4.2674E-01
             2.6885E+00
 GRADIENT:  -4.5800E+01 -7.5438E+01 -2.6207E+01  2.5403E+01  4.9154E+00  7.9874E+00  5.0722E+01 -1.9864E+01  7.2046E+01  5.1693E+00
             1.7379E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -733.120689366281        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0110E+00  2.0420E+00  2.0392E+02  1.1755E+00  3.5674E+00  1.8633E+00  9.9934E-01  1.4209E+01  7.8655E+00  4.3570E-01
             1.2584E+01
 PARAMETER:  1.1092E-01  8.1391E-01  5.4177E+00  2.6168E-01  1.3718E+00  7.2234E-01  9.9340E-02  2.7539E+00  2.1625E+00 -7.3080E-01
             2.6324E+00
 GRADIENT:   2.3226E+01 -1.0076E+00 -3.6119E+00  2.0239E+01  5.0470E+01  4.3026E+00  6.2577E+00  1.9610E+00  1.4801E+01  2.0772E+00
             2.6112E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -748.368066672732        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.1659E-01  1.8218E+00  2.1795E+02  6.6098E-01  2.7179E+00  1.7369E+00  6.3246E-01  1.1186E+01  7.1015E+00  2.4873E-01
             1.2164E+01
 PARAMETER:  1.2903E-02  6.9985E-01  5.4843E+00 -3.1403E-01  1.0999E+00  6.5209E-01 -3.5814E-01  2.5147E+00  2.0603E+00 -1.2914E+00
             2.5985E+00
 GRADIENT:  -1.6073E+01 -7.2891E+00 -3.7427E+00 -3.2030E-01  8.9342E+00 -1.2991E+01 -6.1037E+00  3.9536E+00 -8.4433E+00  8.5851E-01
             8.0563E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -763.353323100882        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.6360E-01  1.9700E+00  7.9903E+04  6.2301E-01  2.4943E+00  1.8093E+00  9.6649E-01  2.9063E+00  7.7368E+00  4.4561E-02
             1.2039E+01
 PARAMETER:  6.2918E-02  7.7802E-01  1.1389E+01 -3.7319E-01  1.0140E+00  6.9293E-01  6.5911E-02  1.1669E+00  2.1460E+00 -3.0109E+00
             2.5882E+00
 GRADIENT:   7.5791E+00 -4.6651E+00  2.6011E-02  5.7241E+00 -1.4767E+01 -4.6534E+00  5.5031E-01  7.6830E-02  4.8457E+00  3.1086E-02
            -5.1179E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -763.942326837161        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      542
 NPARAMETR:  9.5617E-01  2.0008E+00  1.3048E+05  5.7878E-01  2.6110E+00  1.8277E+00  9.3145E-01  2.7935E+00  7.6776E+00  3.9729E-02
             1.1981E+01
 PARAMETER:  5.5185E-02  7.9356E-01  1.1879E+01 -4.4684E-01  1.0597E+00  7.0307E-01  2.8986E-02  1.1273E+00  2.1383E+00 -3.1257E+00
             2.5834E+00
 GRADIENT:  -4.0374E+01  1.8316E+01  4.0957E-02  4.2055E+00 -1.5320E+00  1.5049E+01 -8.2762E+00  2.7088E-02 -2.3162E+01  2.4687E-02
            -9.9135E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -766.332800954111        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      719
 NPARAMETR:  9.7604E-01  2.3193E+00  7.2000E+07  4.1044E-01  2.6497E+00  1.8585E+00  1.1652E+00  6.7837E-01  9.6049E+00  1.0000E-02
             1.1881E+01
 PARAMETER:  7.5750E-02  9.4128E-01  1.8192E+01 -7.9053E-01  1.0745E+00  7.1979E-01  2.5290E-01 -2.8807E-01  2.3623E+00 -5.0038E+00
             2.5749E+00
 GRADIENT:  -3.2167E+00  3.4316E+01 -1.3098E-03  1.2328E+01 -2.4328E-01  8.0036E+00 -4.3641E+00 -1.6605E-02  1.1904E+01  0.0000E+00
            -4.2599E+01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -766.432021652831        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  9.7565E-01  2.3292E+00  9.9611E+07  4.0451E-01  2.6436E+00  1.8574E+00  1.1895E+00  6.3745E-01  9.5066E+00  1.0000E-02
             1.1884E+01
 PARAMETER:  7.5250E-02  9.4457E-01  1.8498E+01 -8.0429E-01  1.0722E+00  7.1987E-01  2.7325E-01 -3.5051E-01  2.3516E+00 -5.0858E+00
             2.5752E+00
 GRADIENT:  -1.0004E+05 -8.7959E+03 -5.4253E+02  1.5749E+04  1.7103E+03  1.3951E+04 -1.3522E+04 -1.6858E+02 -5.1162E+00  0.0000E+00
            -1.4676E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      822
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.4316E-03 -1.0224E-01 -4.3792E-08  3.1135E-02  6.6237E-06
 SE:             2.7478E-02  1.5658E-02  4.0576E-08  2.1872E-02  8.1880E-05
 N:                     100         100         100         100         100

 P VAL.:         7.8681E-01  6.6394E-11  2.8048E-01  1.5458E-01  9.3553E-01

 ETASHRINKSD(%)  7.9437E+00  4.7543E+01  1.0000E+02  2.6727E+01  9.9726E+01
 ETASHRINKVR(%)  1.5256E+01  7.2482E+01  1.0000E+02  4.6311E+01  9.9999E+01
 EBVSHRINKSD(%)  1.0039E+01  4.6289E+01  1.0000E+02  2.2025E+01  9.9606E+01
 EBVSHRINKVR(%)  1.9071E+01  7.1151E+01  1.0000E+02  3.9198E+01  9.9998E+01
 RELATIVEINF(%)  8.0675E+01  1.6763E+01  1.4024E-09  3.6116E+01  1.4651E-03
 EPSSHRINKSD(%)  2.7872E+00
 EPSSHRINKVR(%)  5.4968E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -766.43202165283083     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       887.65733811557993     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    29.19
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -766.432       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  2.33E+00  9.78E+07  4.05E-01  2.64E+00  1.86E+00  1.19E+00  6.37E-01  9.50E+00  1.00E-02  1.19E+01
 


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
+        3.75E+06
 
 TH 2
+       -1.01E+07  3.32E+03
 
 TH 3
+        8.51E-03 -1.21E-03  7.75E-13
 
 TH 4
+        1.86E+05 -2.35E+06 -1.60E-03  9.69E+06
 
 TH 5
+        1.66E+07 -2.06E+05 -2.63E-06  3.71E+06  2.57E+05
 
 TH 6
+        1.48E+06 -1.55E+06 -3.50E-06 -1.42E+04 -8.46E+05  1.01E+04
 
 TH 7
+        1.02E+08 -3.19E+04 -2.38E-05  5.19E+04 -1.31E+06 -3.13E+06  2.00E+05
 
 TH 8
+        5.55E+07 -5.77E+03  5.44E-05  1.35E+05  3.98E+06  8.52E+06  4.84E+05  5.06E+07
 
 TH 9
+        7.68E+05  5.11E+04  6.21E-05  8.66E+05 -8.75E+04 -8.39E+04  8.29E+04 -1.87E+03  2.68E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        6.06E+03  7.95E+01  2.15E-05 -3.45E+05  2.00E+04  9.40E+04 -1.18E+05 -5.84E+05 -1.23E+01  0.00E+00  1.95E+01
 
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
 #CPUT: Total CPU Time in Seconds,       55.327
Stop Time:
Sat Sep 25 06:39:59 CDT 2021
