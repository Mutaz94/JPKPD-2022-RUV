Wed Sep 29 11:53:14 CDT 2021
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
$DATA ../../../../data/spa/A1/dat6.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1282.57869529627        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3872E+02 -1.9041E+01 -5.1277E+01  6.0558E+01  2.2652E+02  5.3770E+01 -1.6731E+01  3.4788E+00 -1.8410E+01 -6.5065E+01
            -6.6150E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1461.57014144286        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.4476E-01  9.6578E-01  1.0255E+00  1.0337E+00  8.8828E-01  8.2015E-01  9.8076E-01  9.4210E-01  1.0744E+00  9.4829E-01
             1.6747E+00
 PARAMETER:  4.3181E-02  6.5182E-02  1.2518E-01  1.3314E-01 -1.8469E-02 -9.8272E-02  8.0572E-02  4.0354E-02  1.7181E-01  4.6905E-02
             6.1566E-01
 GRADIENT:  -4.2448E+01 -1.5546E+00 -1.6348E+01  2.6258E+01  5.9422E+01 -6.5827E+01  6.3761E-01  5.0933E+00  2.6457E+00 -1.0867E+01
            -9.8544E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1471.07295975688        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.4822E-01  1.0937E+00  5.2776E-01  9.3099E-01  6.7518E-01  8.6977E-01  9.0549E-01  2.0012E-01  1.1091E+00  7.2201E-01
             1.7171E+00
 PARAMETER:  4.6830E-02  1.8955E-01 -5.3911E-01  2.8494E-02 -2.9277E-01 -3.9530E-02  7.2467E-04 -1.5089E+00  2.0353E-01 -2.2572E-01
             6.4062E-01
 GRADIENT:  -3.7418E+01  3.5641E+01 -1.0187E+01  4.2137E+01  4.5119E+01 -4.0584E+01 -6.2768E+00  4.5277E-01  7.9834E+00 -5.6963E+00
            -7.2461E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1475.99388320832        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.6086E-01  1.0848E+00  3.3900E-01  8.9386E-01  5.3396E-01  9.2939E-01  1.0952E+00  8.3757E-02  1.0459E+00  6.0787E-01
             1.7963E+00
 PARAMETER:  6.0075E-02  1.8143E-01 -9.8176E-01 -1.2202E-02 -5.2744E-01  2.6769E-02  1.9094E-01 -2.3798E+00  1.4487E-01 -3.9779E-01
             6.8575E-01
 GRADIENT:   3.7268E-01  8.1309E+01  1.3866E+01  3.0832E+01  1.2732E+00 -1.1086E+01  2.1324E+01  9.3155E-02  1.2590E+01  7.1590E+00
            -2.6449E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1489.19334501661        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      409
 NPARAMETR:  1.0126E+00  7.7829E-01  6.6127E-01  1.1397E+00  6.1467E-01  9.9090E-01  1.2537E+00  1.3752E-01  9.2257E-01  7.4110E-01
             1.9823E+00
 PARAMETER:  1.1257E-01 -1.5066E-01 -3.1359E-01  2.3074E-01 -3.8668E-01  9.0855E-02  3.2612E-01 -1.8840E+00  1.9411E-02 -1.9962E-01
             7.8427E-01
 GRADIENT:  -1.7686E+00  2.0991E+01  6.6620E+00  1.1749E+01 -1.3075E+01  7.5801E+00  1.1251E+00  2.0109E-01 -4.9544E+00  5.7004E-01
             5.4831E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1491.33278392188        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      587
 NPARAMETR:  1.0213E+00  5.5412E-01  6.5672E-01  1.2581E+00  5.4589E-01  1.0060E+00  1.6935E+00  1.2635E-02  8.9323E-01  7.1597E-01
             1.9785E+00
 PARAMETER:  1.2107E-01 -4.9038E-01 -3.2050E-01  3.2963E-01 -5.0534E-01  1.0598E-01  6.2679E-01 -4.2712E+00 -1.2913E-02 -2.3411E-01
             7.8232E-01
 GRADIENT:   1.9183E+01  1.5387E+01  1.1654E+01  6.0068E+00 -1.6323E+01  1.2676E+01  6.8020E+00  1.6494E-03  3.9504E+00 -2.3497E+00
             1.1202E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1494.44430450359        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  1.0067E+00  3.3050E-01  6.0706E-01  1.3643E+00  4.6958E-01  9.6948E-01  1.9134E+00  1.0000E-02  8.6432E-01  7.5547E-01
             1.8928E+00
 PARAMETER:  1.0670E-01 -1.0071E+00 -3.9913E-01  4.1068E-01 -6.5592E-01  6.9002E-02  7.4886E-01 -9.7936E+00 -4.5812E-02 -1.8041E-01
             7.3805E-01
 GRADIENT:  -6.8779E+00  6.9635E+00  3.5420E+00  1.8869E+01 -1.2626E+01 -1.3672E+00 -1.0307E+00  0.0000E+00  9.7261E-01 -1.6105E+00
            -3.1875E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1496.80011160201        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  1.0007E+00  1.0211E-01  7.0443E-01  1.5043E+00  4.8145E-01  9.7424E-01  3.2685E+00  1.0000E-02  8.1002E-01  8.1464E-01
             1.9083E+00
 PARAMETER:  1.0071E-01 -2.1817E+00 -2.5037E-01  5.0830E-01 -6.3096E-01  7.3906E-02  1.2843E+00 -2.2256E+01 -1.1069E-01 -1.0500E-01
             7.4619E-01
 GRADIENT:  -2.6666E+00  1.8121E+00  6.0570E+00  2.1323E+01 -8.9547E+00  2.2782E+00 -7.8195E-01  0.0000E+00 -4.8994E+00 -8.1576E-01
            -9.9434E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1497.75676839317        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  9.9830E-01  2.0048E-02  6.5163E-01  1.5194E+00  4.4382E-01  9.6352E-01  7.4244E+00  1.0000E-02  8.1160E-01  8.0824E-01
             1.8907E+00
 PARAMETER:  9.8303E-02 -3.8096E+00 -3.2827E-01  5.1830E-01 -7.1233E-01  6.2833E-02  2.1048E+00 -4.0363E+01 -1.0875E-01 -1.1289E-01
             7.3695E-01
 GRADIENT:  -2.6389E+00  2.4968E-01  5.5468E+00  6.5858E+00 -1.0731E+01 -1.8736E+00 -1.5179E-01  0.0000E+00 -1.6856E+00  4.3406E-01
            -1.7013E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1497.87590652622        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1288
 NPARAMETR:  9.9879E-01  1.0000E-02  6.6378E-01  1.5250E+00  4.5037E-01  9.6794E-01  1.1186E+01  1.0000E-02  8.1406E-01  8.1300E-01
             1.8923E+00
 PARAMETER:  9.8786E-02 -4.5764E+00 -3.0981E-01  5.2200E-01 -6.9768E-01  6.7418E-02  2.5147E+00 -4.8919E+01 -1.0572E-01 -1.0702E-01
             7.3779E-01
 GRADIENT:  -4.1733E-02  0.0000E+00 -3.2377E-01 -1.0879E+00  7.2409E-01  4.5109E-02  3.2347E-02  0.0000E+00  3.6849E-01  2.4435E-01
             2.9369E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1497.87636073913        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1345
 NPARAMETR:  9.9879E-01  1.0000E-02  6.6263E-01  1.5253E+00  4.4967E-01  9.6791E-01  1.1167E+01  1.0000E-02  8.1341E-01  8.1125E-01
             1.8920E+00
 PARAMETER:  9.8789E-02 -4.5736E+00 -3.1155E-01  5.2216E-01 -6.9924E-01  6.7388E-02  2.5130E+00 -4.8882E+01 -1.0651E-01 -1.0918E-01
             7.3761E-01
 GRADIENT:  -9.6870E-02  0.0000E+00  1.4552E-02  2.5026E-01 -6.0482E-02  1.5077E-02  2.2254E-02  0.0000E+00  5.4248E-02 -3.1609E-02
             6.6742E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1345
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2126E-04  8.8410E-04 -3.6312E-05 -9.6575E-03 -1.1604E-02
 SE:             2.9443E-02  1.9838E-03  2.0697E-04  2.7995E-02  2.3967E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9400E-01  6.5584E-01  8.6073E-01  7.3011E-01  6.2826E-01

 ETASHRINKSD(%)  1.3626E+00  9.3354E+01  9.9307E+01  6.2144E+00  1.9707E+01
 ETASHRINKVR(%)  2.7066E+00  9.9558E+01  9.9995E+01  1.2043E+01  3.5531E+01
 EBVSHRINKSD(%)  1.4771E+00  9.4188E+01  9.9278E+01  5.8948E+00  1.8735E+01
 EBVSHRINKVR(%)  2.9323E+00  9.9662E+01  9.9995E+01  1.1442E+01  3.3960E+01
 RELATIVEINF(%)  8.8251E+01  1.6084E-02  2.8121E-04  6.6955E+00  2.9131E+00
 EPSSHRINKSD(%)  3.7429E+01
 EPSSHRINKVR(%)  6.0849E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1497.8763607391325     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -762.72553417539427     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1497.876       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  1.00E-02  6.63E-01  1.53E+00  4.50E-01  9.68E-01  1.12E+01  1.00E-02  8.13E-01  8.11E-01  1.89E+00
 


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
+        1.15E+03
 
 TH 2
+        0.00E+00  6.96E+04
 
 TH 3
+        3.96E+00  0.00E+00  1.33E+03
 
 TH 4
+       -2.61E+01  0.00E+00 -1.56E+02  6.41E+02
 
 TH 5
+        3.70E+01  0.00E+00 -2.55E+03 -1.29E+02  5.52E+03
 
 TH 6
+        5.58E-01  0.00E+00  6.63E+00 -7.05E+00 -6.35E-01  1.99E+02
 
 TH 7
+       -3.04E-01  0.00E+00  8.92E+00  5.78E+00 -2.10E+01  4.06E-02  4.74E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.35E+00  0.00E+00  2.97E+01 -1.25E+01 -6.15E+00  1.10E+00  2.19E+00  0.00E+00  2.39E+02
 
 TH10
+       -3.94E+00  0.00E+00 -1.64E+01  1.28E+00 -6.68E+01 -5.52E-01  5.06E+00  0.00E+00  3.03E+00  1.31E+02
 
 TH11
+       -1.23E+01  0.00E+00 -2.04E+01 -7.68E+00  9.59E+00  3.17E+00  1.08E+00  0.00E+00  9.03E+00  2.63E+01  7.03E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.641
Stop Time:
Wed Sep 29 11:53:38 CDT 2021
