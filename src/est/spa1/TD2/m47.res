Thu Sep 30 02:08:53 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat47.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2079.66309556366        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3816E+02 -3.6269E+01 -5.9092E+00 -3.3110E+01  5.6393E+01  6.0436E+01 -2.9818E+00 -6.8576E+00 -1.9544E+01  7.2601E+00
            -1.9108E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2087.32203204686        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      178
 NPARAMETR:  9.8746E-01  1.0032E+00  9.7090E-01  1.0637E+00  9.4941E-01  9.2560E-01  1.0130E+00  1.0403E+00  1.0881E+00  9.4668E-01
             1.0293E+00
 PARAMETER:  8.7385E-02  1.0315E-01  7.0471E-02  1.6180E-01  4.8085E-02  2.2686E-02  1.1287E-01  1.3949E-01  1.8446E-01  4.5210E-02
             1.2883E-01
 GRADIENT:   1.0137E+01 -6.7754E+00 -5.7637E+00 -8.4690E-01  1.8597E+01 -5.9396E+00  1.5895E+00 -3.2120E+00  2.3027E+00  7.8315E+00
             5.4199E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2088.02408887104        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.8501E-01  9.8704E-01  9.5181E-01  1.0671E+00  9.0905E-01  9.3573E-01  1.0059E+00  1.0917E+00  1.0880E+00  8.8145E-01
             1.0049E+00
 PARAMETER:  8.4894E-02  8.6960E-02  5.0611E-02  1.6494E-01  4.6499E-03  3.3567E-02  1.0593E-01  1.8771E-01  1.8430E-01 -2.6192E-02
             1.0487E-01
 GRADIENT:   3.9792E+02  3.4831E+01  1.0396E+01  1.1083E+02 -2.0046E+00  3.7123E+01  4.5008E+00 -2.4410E-01  2.1492E+01  5.0446E+00
            -1.0976E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2088.04712554324        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  9.8550E-01  9.8704E-01  9.5181E-01  1.0671E+00  9.0905E-01  9.3497E-01  1.0054E+00  1.0917E+00  1.0597E+00  8.8145E-01
             1.0332E+00
 PARAMETER:  8.5395E-02  8.6960E-02  5.0611E-02  1.6494E-01  4.6499E-03  3.2762E-02  1.0540E-01  1.8771E-01  1.5798E-01 -2.6192E-02
             1.3266E-01
 GRADIENT:   4.4861E+00  2.0402E-01  7.4104E+00 -1.3260E+01 -1.0670E+01 -1.8080E+00 -3.6610E-01 -3.3916E-01 -3.2865E+00  4.9684E+00
             9.7626E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2088.28109434015        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      628
 NPARAMETR:  9.8454E-01  9.8746E-01  9.4508E-01  1.0708E+00  9.1158E-01  9.3917E-01  1.0032E+00  1.0964E+00  1.0747E+00  8.6453E-01
             1.0204E+00
 PARAMETER:  8.4422E-02  8.7383E-02  4.3510E-02  1.6842E-01  7.4220E-03  3.7240E-02  1.0316E-01  1.9200E-01  1.7201E-01 -4.5572E-02
             1.2018E-01
 GRADIENT:   3.8420E+02  3.2661E+01  4.0182E+00  1.1625E+02  9.3899E+00  3.7500E+01  3.2197E+00  3.3412E-01  1.6831E+01  2.8143E+00
             1.0651E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2088.29310506853        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      791
 NPARAMETR:  9.8444E-01  9.8726E-01  9.4351E-01  1.0740E+00  9.1159E-01  9.3933E-01  1.0031E+00  1.0971E+00  1.0728E+00  8.6251E-01
             1.0205E+00
 PARAMETER:  8.4319E-02  8.7182E-02  4.1853E-02  1.7138E-01  7.4362E-03  3.7410E-02  1.0307E-01  1.9271E-01  1.7030E-01 -4.7911E-02
             1.2034E-01
 GRADIENT:   1.8764E+00  2.6632E-01  2.1139E-01  1.6295E-01  3.6496E+00 -3.8908E-02 -3.4987E-01 -1.1848E-01  2.1372E-02  1.8944E+00
            -6.3014E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2088.36275032442        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      978             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8346E-01  9.8817E-01  9.2984E-01  1.0722E+00  9.0309E-01  9.3904E-01  1.0256E+00  1.0976E+00  1.0646E+00  8.4309E-01
             1.0208E+00
 PARAMETER:  8.3318E-02  8.8095E-02  2.7261E-02  1.6972E-01 -1.9277E-03  3.7102E-02  1.2526E-01  1.9309E-01  1.6257E-01 -7.0680E-02
             1.2057E-01
 GRADIENT:   3.8041E+02  3.4376E+01  3.3744E+00  1.2093E+02  9.9406E+00  3.7444E+01  4.0518E+00  1.1673E+00  1.4890E+01  1.6020E+00
             1.5332E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2088.36637216727        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1161            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8379E-01  9.8785E-01  9.2922E-01  1.0727E+00  9.0311E-01  9.3954E-01  1.0256E+00  1.0966E+00  1.0673E+00  8.4244E-01
             1.0210E+00
 PARAMETER:  8.3653E-02  8.7779E-02  2.6592E-02  1.7014E-01 -1.9120E-03  3.7637E-02  1.2529E-01  1.9223E-01  1.6517E-01 -7.1452E-02
             1.2078E-01
 GRADIENT:   3.8108E+02  3.4135E+01  2.8933E+00  1.2187E+02  1.0795E+01  3.7604E+01  4.1568E+00  1.1599E+00  1.5843E+01  1.5005E+00
             1.6959E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2088.38596685644        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1346
 NPARAMETR:  9.8428E-01  9.8727E-01  9.2235E-01  1.0736E+00  9.0013E-01  9.3976E-01  1.0250E+00  1.0773E+00  1.0674E+00  8.3070E-01
             1.0206E+00
 PARAMETER:  8.4153E-02  8.7188E-02  1.9174E-02  1.7099E-01 -5.2117E-03  3.7869E-02  1.2466E-01  1.7447E-01  1.6523E-01 -8.5488E-02
             1.2044E-01
 GRADIENT:   1.1005E+00  7.0612E-01  4.2879E-01  1.6847E+00  5.4015E+00  8.6935E-02 -6.5286E-01 -3.4986E-01  1.1098E-01 -7.5875E-01
            -1.6510E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -2088.38613665731        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1403
 NPARAMETR:  9.8428E-01  9.8727E-01  9.2143E-01  1.0736E+00  9.0013E-01  9.3964E-01  1.0250E+00  1.0773E+00  1.0671E+00  8.3070E-01
             1.0206E+00
 PARAMETER:  8.4153E-02  8.7188E-02  1.8174E-02  1.7099E-01 -5.2118E-03  3.7738E-02  1.2466E-01  1.7447E-01  1.6496E-01 -8.5488E-02
             1.2044E-01
 GRADIENT:   1.0892E+00  3.9915E-01 -9.2700E-02  1.8921E+00  6.1657E+00  3.3272E-02 -6.7552E-01 -2.7704E-01  5.0714E-02 -7.3809E-01
            -1.3985E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1403
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6550E-04 -1.6567E-02 -2.8091E-02  5.2310E-03 -3.4743E-02
 SE:             2.9875E-02  1.7326E-02  1.6224E-02  2.6437E-02  2.0634E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8490E-01  3.3897E-01  8.3377E-02  8.4315E-01  9.2229E-02

 ETASHRINKSD(%)  1.0000E-10  4.1956E+01  4.5646E+01  1.1433E+01  3.0872E+01
 ETASHRINKVR(%)  1.0000E-10  6.6309E+01  7.0457E+01  2.1559E+01  5.2214E+01
 EBVSHRINKSD(%)  3.8362E-01  4.2362E+01  4.7822E+01  1.1704E+01  3.0242E+01
 EBVSHRINKVR(%)  7.6576E-01  6.6779E+01  7.2775E+01  2.2037E+01  5.1339E+01
 RELATIVEINF(%)  9.8572E+01  1.7203E+00  3.7456E+00  5.5620E+00  7.9749E+00
 EPSSHRINKSD(%)  3.4675E+01
 EPSSHRINKVR(%)  5.7326E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2088.3861366573092     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1169.4476034526365     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2088.386       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  9.87E-01  9.21E-01  1.07E+00  9.00E-01  9.40E-01  1.02E+00  1.08E+00  1.07E+00  8.31E-01  1.02E+00
 


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
+        6.89E+07
 
 TH 2
+       -1.16E+01  6.84E+07
 
 TH 3
+        5.65E+00  1.75E+02  3.08E+02
 
 TH 4
+       -1.75E+00  3.29E+02 -1.14E+02  1.98E+07
 
 TH 5
+        3.44E+00 -7.51E+07 -4.02E+07  1.15E+02  8.23E+07
 
 TH 6
+       -9.75E-01 -1.80E+00  9.61E-01 -6.34E-01 -2.66E+02  2.23E+02
 
 TH 7
+        4.71E-01  1.13E+01  4.48E+00 -3.66E-01 -9.49E+00 -4.77E-01  4.09E+07
 
 TH 8
+       -1.67E+00 -3.60E+07 -1.05E+02  2.98E+01  3.94E+07 -6.38E-01  5.43E+00  1.89E+07
 
 TH 9
+        2.48E+00 -1.50E+01 -3.56E+00  3.85E+01 -2.10E+07 -9.47E-02  1.95E+01  1.32E+02  1.08E+02
 
 TH10
+        1.30E+00 -8.13E+07 -1.35E+01 -6.76E+00  4.46E+07  7.67E-02  1.15E+01  4.27E+07  9.62E-01  9.67E+07
 
 TH11
+       -8.10E+00 -7.32E+00 -1.47E+01 -5.06E+00 -3.01E+07  1.89E+00  2.14E+00  4.29E+01  3.59E+00  1.49E+01  3.79E+02
 
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
 #CPUT: Total CPU Time in Seconds,       29.881
Stop Time:
Thu Sep 30 02:09:24 CDT 2021
