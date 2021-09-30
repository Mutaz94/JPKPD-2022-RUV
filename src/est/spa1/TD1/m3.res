Thu Sep 30 01:05:53 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2148.98383924513        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0479E+02 -3.7821E+01 -4.2012E+01  1.5032E+01  4.8233E+01  6.9610E+01  7.2176E+00  1.3962E+01  1.8053E+01  2.9213E+01
            -4.0058E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2157.75216269093        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0227E+00  1.0794E+00  1.1655E+00  1.0140E+00  1.0582E+00  1.0110E+00  9.5248E-01  9.1351E-01  9.2919E-01  7.8779E-01
             1.1056E+00
 PARAMETER:  1.2243E-01  1.7641E-01  2.5317E-01  1.1388E-01  1.5658E-01  1.1095E-01  5.1309E-02  9.5376E-03  2.6559E-02 -1.3853E-01
             2.0040E-01
 GRADIENT:  -3.2856E+00 -8.7705E-01  1.0538E+01 -8.5316E+00  2.0457E+01  2.2078E+00 -2.8143E+00 -5.8569E+00 -1.0731E+01 -7.4991E+00
             2.5120E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2158.79902031901        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0260E+00  1.0593E+00  9.8081E-01  1.0172E+00  9.5825E-01  1.0084E+00  1.0646E+00  8.6398E-01  9.1975E-01  6.3662E-01
             1.0990E+00
 PARAMETER:  1.2567E-01  1.5763E-01  8.0624E-02  1.1706E-01  5.7350E-02  1.0837E-01  1.6262E-01 -4.6211E-02  1.6348E-02 -3.5158E-01
             1.9442E-01
 GRADIENT:   2.3169E+00  9.5511E-01  4.1487E+00 -2.3778E+00  1.0534E+01  7.7712E-01  2.8367E-01 -1.0017E+00 -5.1389E+00 -9.4223E+00
             2.1844E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2160.40426447892        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0257E+00  1.1639E+00  8.3358E-01  9.4378E-01  9.4109E-01  1.0078E+00  9.6543E-01  5.8998E-01  9.9605E-01  7.3200E-01
             1.0594E+00
 PARAMETER:  1.2536E-01  2.5178E-01 -8.2024E-02  4.2134E-02  3.9289E-02  1.0780E-01  6.4816E-02 -4.2767E-01  9.6041E-02 -2.1197E-01
             1.5769E-01
 GRADIENT:  -2.2597E-01 -1.4343E+00 -3.9350E-01  1.7211E-01  4.7443E-01 -8.2678E-02 -3.6581E-01  2.7488E-01  4.4857E-01 -1.0169E-01
            -4.5184E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2160.47212385519        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0272E+00  1.3299E+00  7.3368E-01  8.3817E-01  9.7345E-01  1.0096E+00  8.9088E-01  4.2558E-01  1.0791E+00  7.3994E-01
             1.0624E+00
 PARAMETER:  1.2679E-01  3.8509E-01 -2.0969E-01 -7.6540E-02  7.3095E-02  1.0952E-01 -1.5543E-02 -7.5431E-01  1.7617E-01 -2.0119E-01
             1.6057E-01
 GRADIENT:   2.9384E-01  3.8596E+00  2.2054E-02  1.9299E+00 -2.4753E+00 -4.3627E-02  4.2593E-01  3.4463E-01 -5.2468E-01  5.2745E-01
             2.6765E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2160.48361445555        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0274E+00  1.3852E+00  7.0390E-01  8.0082E-01  9.8979E-01  1.0100E+00  8.6622E-01  3.5023E-01  1.1169E+00  7.4484E-01
             1.0632E+00
 PARAMETER:  1.2703E-01  4.2587E-01 -2.5112E-01 -1.2212E-01  8.9738E-02  1.0993E-01 -4.3616E-02 -9.4917E-01  2.1060E-01 -1.9458E-01
             1.6132E-01
 GRADIENT:   1.8784E-01  1.0556E+00 -5.7348E-01  5.0719E-01 -5.3034E-01 -4.4117E-02  3.7977E-01  2.8270E-01 -2.4613E-01  4.0392E-01
             2.0856E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2160.52456834312        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.0281E+00  1.3864E+00  6.9794E-01  7.9995E-01  9.8853E-01  1.0106E+00  8.6688E-01  2.5961E-01  1.1216E+00  7.4902E-01
             1.0629E+00
 PARAMETER:  1.2770E-01  4.2674E-01 -2.5962E-01 -1.2320E-01  8.8462E-02  1.1058E-01 -4.2860E-02 -1.2486E+00  2.1475E-01 -1.8899E-01
             1.6098E-01
 GRADIENT:   1.6289E+00  9.2883E-01  4.4746E-01  1.6956E+00  7.9404E-01  2.2183E-01  3.5746E-01  5.2604E-02  3.1498E-01  4.7853E-01
             1.2051E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2160.56032358709        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1253
 NPARAMETR:  1.0262E+00  1.3969E+00  6.7324E-01  7.8940E-01  9.8057E-01  1.0094E+00  8.6352E-01  1.1424E-01  1.1255E+00  7.3251E-01
             1.0630E+00
 PARAMETER:  1.2583E-01  4.3429E-01 -2.9565E-01 -1.3649E-01  8.0383E-02  1.0933E-01 -4.6743E-02 -2.0695E+00  2.1820E-01 -2.1128E-01
             1.6114E-01
 GRADIENT:  -2.8803E+00 -1.7063E+00 -7.4531E-02 -7.7363E-01 -4.4483E-02 -3.9529E-01 -3.7252E-01  2.1700E-02 -3.3795E-01 -5.0191E-01
            -1.7628E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2160.56784981788        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1429
 NPARAMETR:  1.0256E+00  1.3992E+00  6.7205E-01  7.8806E-01  9.8177E-01  1.0090E+00  8.6388E-01  3.2906E-02  1.1286E+00  7.3857E-01
             1.0627E+00
 PARAMETER:  1.2532E-01  4.3587E-01 -2.9742E-01 -1.3818E-01  8.1601E-02  1.0897E-01 -4.6316E-02 -3.3141E+00  2.2101E-01 -2.0303E-01
             1.6079E-01
 GRADIENT:  -4.0196E+00 -1.9510E+00 -1.4022E-01 -4.4001E-01  2.9180E-01 -5.3997E-01  3.3595E-02  1.7133E-03  9.7243E-03  1.3111E-01
            -4.0613E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2160.57851773707        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1615             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0298E+00  1.3989E+00  6.7154E-01  7.8835E-01  9.8110E-01  1.0114E+00  8.6385E-01  1.0000E-02  1.1284E+00  7.3697E-01
             1.0627E+00
 PARAMETER:  1.2936E-01  4.3572E-01 -2.9818E-01 -1.3782E-01  8.0917E-02  1.1135E-01 -4.6351E-02 -1.1602E+01  2.2079E-01 -2.0520E-01
             1.6083E-01
 GRADIENT:   5.2727E+02  2.9108E+02  4.6627E+00  6.6520E+01  7.0531E+00  6.9350E+01  5.1237E+00  0.0000E+00  1.2900E+01  7.1728E-01
             1.5102E+00

0ITERATION NO.:   48    OBJECTIVE VALUE:  -2160.57909616397        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1687
 NPARAMETR:  1.0285E+00  1.3983E+00  6.7153E-01  7.8823E-01  9.8109E-01  1.0105E+00  8.6377E-01  1.0000E-02  1.1282E+00  7.3695E-01
             1.0627E+00
 PARAMETER:  1.2807E-01  4.3529E-01 -2.9820E-01 -1.3796E-01  8.0908E-02  1.1042E-01 -4.6446E-02 -1.1602E+01  2.2063E-01 -2.0524E-01
             1.6083E-01
 GRADIENT:   2.0536E+00 -2.4389E+00 -1.2837E-01 -7.3138E-01  4.2837E-01  5.4822E-02 -1.0051E-01  0.0000E+00 -7.8384E-03 -2.5305E-02
            -6.1329E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1687
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8517E-04 -1.6003E-02 -3.0730E-04  1.2193E-02 -2.2868E-02
 SE:             2.9898E-02  2.2675E-02  1.6016E-04  2.4719E-02  2.1263E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8438E-01  4.8036E-01  5.5012E-02  6.2183E-01  2.8215E-01

 ETASHRINKSD(%)  1.0000E-10  2.4034E+01  9.9463E+01  1.7189E+01  2.8767E+01
 ETASHRINKVR(%)  1.0000E-10  4.2292E+01  9.9997E+01  3.1423E+01  4.9258E+01
 EBVSHRINKSD(%)  3.5965E-01  2.3814E+01  9.9519E+01  1.7377E+01  2.8934E+01
 EBVSHRINKVR(%)  7.1800E-01  4.1958E+01  9.9998E+01  3.1734E+01  4.9496E+01
 RELATIVEINF(%)  9.9173E+01  4.1101E+00  3.9647E-04  5.8468E+00  6.1766E+00
 EPSSHRINKSD(%)  3.2704E+01
 EPSSHRINKVR(%)  5.4713E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2160.5790961639655     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1241.6405629592928     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2160.579       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.40E+00  6.72E-01  7.88E-01  9.81E-01  1.01E+00  8.64E-01  1.00E-02  1.13E+00  7.37E-01  1.06E+00
 


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
+       -6.08E+00  4.59E+02
 
 TH 3
+        6.39E+00  2.52E+02  6.84E+02
 
 TH 4
+       -5.25E+00  3.21E+02 -3.29E+02  9.67E+02
 
 TH 5
+        1.18E+00 -3.75E+02 -6.37E+02  3.73E+02  1.06E+03
 
 TH 6
+        8.16E-01 -1.36E+00  2.28E+00 -1.68E+00  2.94E-01  1.93E+02
 
 TH 7
+        1.03E+00  2.37E+01 -1.76E+01 -7.87E+00 -9.27E+00  5.83E-01  9.60E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.72E-01 -2.50E+01 -2.34E+01  4.78E+01 -7.19E+00 -2.93E-01  1.59E+01  0.00E+00  8.37E+01
 
 TH10
+        9.40E-01 -7.42E+00 -5.19E+01 -1.96E+01 -5.88E+01  5.51E-01  3.09E+01  0.00E+00  8.20E+00  9.70E+01
 
 TH11
+       -7.49E+00 -1.94E+01 -2.53E+01 -5.68E+00 -6.88E+00  9.59E-01  7.95E+00  0.00E+00  6.26E+00  2.81E+01  3.66E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.744
Stop Time:
Thu Sep 30 01:06:26 CDT 2021
