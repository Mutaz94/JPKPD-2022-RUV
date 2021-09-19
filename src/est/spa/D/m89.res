Sat Sep 18 15:42:33 CDT 2021
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
$DATA ../../../../data/spa/D/dat89.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   24441.9268672017        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3024E+02  5.0177E+02  4.3144E-01  5.5140E+02  2.7277E+02 -3.2234E+03 -1.1110E+03 -6.6411E+01 -1.4580E+03 -9.6246E+02
            -4.5110E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -466.988618877870        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2590E+00  1.1448E+00  8.1277E-01  1.6562E+00  1.6011E+00  1.7134E+00  9.6009E-01  9.5352E-01  7.2144E-01  9.1774E-01
             1.5509E+01
 PARAMETER:  3.3029E-01  2.3527E-01 -1.0731E-01  6.0455E-01  5.7067E-01  6.3848E-01  5.9275E-02  5.2407E-02 -2.2650E-01  1.4156E-02
             2.8414E+00
 GRADIENT:  -2.6882E+01  3.7505E+01 -5.2572E+00  7.4845E+01 -4.1701E+00  4.0434E+00 -8.7375E-01  4.0709E+00  3.7724E+00  2.6272E-01
            -5.5816E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -496.590045622977        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.3754E+00  6.5151E-01  3.3106E-01  1.5515E+00  4.2509E-01  1.7956E+00  1.7559E+00  1.0000E-02  1.1620E-01  6.0804E-01
             1.6596E+01
 PARAMETER:  4.1871E-01 -3.2846E-01 -1.0055E+00  5.3924E-01 -7.5545E-01  6.8533E-01  6.6298E-01 -5.0328E+00 -2.0524E+00 -3.9752E-01
             2.9092E+00
 GRADIENT:   1.1490E+01  1.4756E+01  5.4291E+00 -1.3027E+00  2.1061E+01  7.7821E+00  9.3069E+00  0.0000E+00  5.0873E-01  6.9060E+00
             1.9063E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -571.446104694764        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      237
 NPARAMETR:  2.0880E+00  1.5837E-01  2.1849E-02  1.8200E+00  8.9529E-02  3.2466E+00  1.8429E+00  1.0000E-02  6.6403E+00  1.0000E-02
             1.6516E+01
 PARAMETER:  8.3621E-01 -1.7428E+00 -3.7236E+00  6.9884E-01 -2.3132E+00  1.2776E+00  7.1132E-01 -2.5533E+01  1.9932E+00 -6.9518E+00
             2.9044E+00
 GRADIENT:   6.2922E+01  1.3666E+01 -3.5935E+00  4.5887E+00  3.3682E+01 -1.0864E+01  6.8716E+00  0.0000E+00  5.6750E+00  0.0000E+00
             5.6348E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -606.219109771446        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.2643E+00  4.9683E-02  2.9533E-02  2.3490E+00  8.9827E-02  3.4719E+00  2.4369E+00  1.0000E-02  5.8093E-01  9.3187E-02
             1.3297E+01
 PARAMETER:  3.3452E-01 -2.9021E+00 -3.4222E+00  9.5400E-01 -2.3099E+00  1.3447E+00  9.9071E-01 -2.3210E+01 -4.4313E-01 -2.2731E+00
             2.6875E+00
 GRADIENT:   6.4871E+00 -7.2529E-01 -1.8882E+00 -6.3628E+00  1.5639E+01  7.7276E+00  3.8483E-01  0.0000E+00  3.8808E-01  8.0925E-01
             6.5122E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -618.223148234316        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  7.9441E-01  1.0000E-02  4.5889E-02  9.3593E+00  8.7019E-02  3.9053E+00  4.4553E+00  1.0000E-02  5.9697E-01  1.0000E-02
             1.2926E+01
 PARAMETER: -1.3016E-01 -6.3454E+00 -2.9815E+00  2.3364E+00 -2.3416E+00  1.4623E+00  1.5941E+00 -3.0332E+01 -4.1589E-01 -6.1664E+00
             2.6592E+00
 GRADIENT:  -1.3709E+01  0.0000E+00 -1.4712E+01 -2.6597E+00  1.6549E+01 -2.4503E+01  1.2136E-03  0.0000E+00  3.0526E-01  0.0000E+00
            -8.3793E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -623.821716509212        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      456
 NPARAMETR:  6.4784E-01  1.0000E-02  5.7977E-02  1.7707E+01  8.4053E-02  4.7641E+00  5.0500E+00  1.0000E-02  1.6276E+00  1.0000E-02
             1.3367E+01
 PARAMETER: -3.3412E-01 -7.8484E+00 -2.7477E+00  2.9739E+00 -2.3763E+00  1.6611E+00  1.7194E+00 -3.4776E+01  5.8712E-01 -9.6671E+00
             2.6928E+00
 GRADIENT:  -6.7894E+00  0.0000E+00 -7.1639E+00 -3.6889E+00 -9.3738E+00 -4.1478E+00 -8.4471E-03  0.0000E+00  4.0599E-01  0.0000E+00
             1.2242E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -625.142735182941        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  5.7189E-01  1.0000E-02  6.4644E-02  2.6443E+01  8.3549E-02  5.2921E+00  5.5597E+00  1.0000E-02  2.4720E+00  1.0000E-02
             1.3340E+01
 PARAMETER: -4.5881E-01 -8.7889E+00 -2.6389E+00  3.3750E+00 -2.3823E+00  1.7662E+00  1.8155E+00 -3.7315E+01  1.0050E+00 -1.1294E+01
             2.6907E+00
 GRADIENT:  -1.0987E+00  0.0000E+00 -6.7416E-01 -9.5784E-01 -2.3430E+00 -2.3562E+00 -1.5853E-02  0.0000E+00  3.7432E-01  0.0000E+00
            -1.4701E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -625.239923021805        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      733            RESET HESSIAN, TYPE II
 NPARAMETR:  5.6612E-01  1.0000E-02  6.5883E-02  2.9810E+01  8.3304E-02  5.4401E+00  9.8419E+00  1.0000E-02  2.6469E+00  1.0000E-02
             1.3405E+01
 PARAMETER: -4.6895E-01 -9.0238E+00 -2.6199E+00  3.4948E+00 -2.3853E+00  1.7938E+00  2.3867E+00 -3.7999E+01  1.0734E+00 -1.1743E+01
             2.6957E+00
 GRADIENT:   1.3527E+00  0.0000E+00  2.0139E+00  8.8645E-01  5.7945E+00  3.2326E+00 -2.4533E-01  0.0000E+00 -1.2793E-01  0.0000E+00
             1.0415E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -625.294272309582        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  5.6028E-01  1.0000E-02  6.5757E-02  2.9471E+01  8.3155E-02  5.4340E+00  1.1473E+01  1.0000E-02  2.6364E+00  1.0000E-02
             1.3398E+01
 PARAMETER: -4.7932E-01 -9.0238E+00 -2.6218E+00  3.4834E+00 -2.3870E+00  1.7927E+00  2.5400E+00 -3.7999E+01  1.0694E+00 -1.1743E+01
             2.6951E+00
 GRADIENT:   2.1025E-01  0.0000E+00 -1.2333E-01  6.2497E-02 -3.3676E-01 -1.5156E-01  1.1619E-02  0.0000E+00 -1.8102E-02  0.0000E+00
            -3.2742E-01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -625.295979163470        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      996
 NPARAMETR:  5.5571E-01  1.0000E-02  6.6009E-02  2.9413E+01  8.3251E-02  5.4530E+00  1.1474E+01  1.0000E-02  2.6303E+00  1.0000E-02
             1.3403E+01
 PARAMETER: -4.8752E-01 -9.0238E+00 -2.6180E+00  3.4814E+00 -2.3859E+00  1.7962E+00  2.5401E+00 -3.7999E+01  1.0671E+00 -1.1743E+01
             2.6955E+00
 GRADIENT:   7.6197E-03  0.0000E+00  1.2217E-02 -3.9485E-02 -8.0103E-02  3.8150E-02  1.2622E-03  0.0000E+00  8.4175E-03  0.0000E+00
            -9.5606E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      996
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8079E-02  5.2138E-03 -4.9564E-04 -3.6255E-02  2.3525E-04
 SE:             2.5699E-02  1.6857E-03  1.3944E-04  1.1320E-02  2.2413E-04
 N:                     100         100         100         100         100

 P VAL.:         2.7457E-01  1.9824E-03  3.7861E-04  1.3609E-03  2.9389E-01

 ETASHRINKSD(%)  1.3906E+01  9.4353E+01  9.9533E+01  6.2077E+01  9.9249E+01
 ETASHRINKVR(%)  2.5878E+01  9.9681E+01  9.9998E+01  8.5619E+01  9.9994E+01
 EBVSHRINKSD(%)  8.7343E+00  9.6976E+01  9.9643E+01  6.8562E+01  9.9018E+01
 EBVSHRINKVR(%)  1.6706E+01  9.9909E+01  9.9999E+01  9.0117E+01  9.9990E+01
 RELATIVEINF(%)  4.6032E+01  2.1133E-02  2.1124E-04  4.4138E+00  3.8544E-03
 EPSSHRINKSD(%)  1.8263E+00
 EPSSHRINKVR(%)  3.6193E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -625.29597916346984     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       109.85484740026834     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.82
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
 





 #OBJV:********************************************     -625.296       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.56E-01  1.00E-02  6.60E-02  2.94E+01  8.33E-02  5.45E+00  1.15E+01  1.00E-02  2.63E+00  1.00E-02  1.34E+01
 


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
+        9.96E+01
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        8.92E+02  0.00E+00  4.94E+04
 
 TH 4
+        3.89E-01  0.00E+00 -1.46E+00  2.51E-02
 
 TH 5
+        8.34E+01  0.00E+00 -3.10E+04  2.59E+01  1.09E+05
 
 TH 6
+        4.35E+00  0.00E+00 -2.04E+02 -8.79E-02  2.31E+01  4.34E+00
 
 TH 7
+        1.42E-02  0.00E+00 -4.78E+00  7.36E-03  2.55E+00 -3.75E-03  7.92E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.33E+00  0.00E+00 -3.02E+01 -1.95E-01 -9.78E+00  6.11E-01 -1.67E-01  0.00E+00  2.76E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.07E+00  0.00E+00 -1.49E+02 -3.90E-02  4.77E+01  5.35E-01 -2.31E-02  0.00E+00  4.09E-01  0.00E+00  2.35E+00
 
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
 #CPUT: Total CPU Time in Seconds,       19.080
Stop Time:
Sat Sep 18 15:42:54 CDT 2021
