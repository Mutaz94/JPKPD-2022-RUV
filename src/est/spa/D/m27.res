Sat Sep 18 15:14:28 CDT 2021
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
$DATA ../../../../data/spa/D/dat27.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m27.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   6454.42760979643        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.0411E+02 -2.5343E+01 -7.0306E+01 -1.1245E+02  3.2794E+02 -1.1091E+03 -4.3586E+02 -5.3525E+01 -7.6153E+02 -4.8303E+02
            -1.3465E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -712.969913430352        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.5759E+00  1.1512E+00  1.0483E+00  1.7395E+00  1.0597E+00  2.0920E+00  1.4991E+00  9.9959E-01  1.8086E+00  1.3398E+00
             1.3214E+01
 PARAMETER:  5.5480E-01  2.4085E-01  1.4712E-01  6.5358E-01  1.5802E-01  8.3810E-01  5.0488E-01  9.9591E-02  6.9257E-01  3.9253E-01
             2.6813E+00
 GRADIENT:   4.3796E+01  1.1744E+01 -7.9886E+00  1.4524E+01 -1.4825E+01  4.6680E+01  1.5639E+00  4.7999E+00  1.7104E+01  8.3271E+00
             2.1230E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -741.468934030928        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.5111E+00  6.6674E-01  5.4814E+00  2.7347E+00  3.4561E+00  2.1993E+00  7.7830E+00  7.0609E-01  2.4431E+00  5.6204E+00
             9.9684E+00
 PARAMETER:  5.1287E-01 -3.0535E-01  1.8014E+00  1.1060E+00  1.3401E+00  8.8816E-01  2.1519E+00 -2.4801E-01  9.9326E-01  1.8264E+00
             2.3994E+00
 GRADIENT:   4.7922E+01  1.8035E+01  4.2514E+00  5.1982E+01 -1.0254E+01  8.2280E+00  2.2154E+01 -8.4249E-02  3.7217E+01  1.2272E+01
             1.5337E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -798.238810582398        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.1341E+00  1.8134E+00  7.7984E+00  9.7824E-01  7.1781E+00  1.9987E+00  1.1329E+00  6.6529E+00  2.7537E+00  7.7185E+00
             8.2420E+00
 PARAMETER:  2.2580E-01  6.9519E-01  2.1539E+00  7.8005E-02  2.0710E+00  7.9252E-01  2.2477E-01  1.9951E+00  1.1129E+00  2.1436E+00
             2.2092E+00
 GRADIENT:  -4.4134E+01  2.7358E+01  4.2849E+00  2.6483E-02 -4.7744E+00  1.1569E+01 -1.1324E+00 -1.7180E+00 -1.8000E+01  9.3310E+00
             8.5334E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -820.276744731906        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.2111E+00  1.5493E+00  3.2372E+00  9.4355E-01  4.4060E+00  2.1148E+00  1.1028E+00  3.2890E+00  3.4114E+00  4.5758E+00
             6.6722E+00
 PARAMETER:  2.9151E-01  5.3781E-01  1.2747E+00  4.1893E-02  1.5830E+00  8.4896E-01  1.9789E-01  1.2906E+00  1.3271E+00  1.6208E+00
             1.9979E+00
 GRADIENT:   4.8438E+00  1.1153E-02  1.6416E+00 -1.1888E+01 -8.2669E+00 -2.0710E+00  6.2569E+00  3.6634E+00 -1.4643E+00  3.4551E+00
             1.7360E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -841.126171903452        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.1136E+00  8.0943E-01  8.3867E-01  1.4691E+00  3.0031E+01  2.0506E+00  1.6060E-01  3.0949E-02  2.3184E+00  7.2597E+00
             6.8558E+00
 PARAMETER:  2.0761E-01 -1.1142E-01 -7.5944E-02  4.8466E-01  3.5022E+00  8.1814E-01 -1.7288E+00 -3.3754E+00  9.4088E-01  2.0823E+00
             2.0251E+00
 GRADIENT:  -8.7668E+00  2.8993E+01 -8.7866E+00  2.4598E+01 -9.7512E-01  7.9822E+00  1.5834E-01  4.3942E-03  7.9535E+00  5.2134E-03
             7.9374E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -884.520669872425        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.1448E-01  1.5099E-01  2.2051E-01  1.1181E+00  3.7695E+02  1.6552E+00  1.0000E-02  1.0000E-02  1.0372E+00  8.8693E+00
             6.0993E+00
 PARAMETER:  1.0604E-02 -1.7905E+00 -1.4118E+00  2.1160E-01  6.0321E+00  6.0391E-01 -6.1992E+00 -9.7196E+00  1.3656E-01  2.2826E+00
             1.9082E+00
 GRADIENT:   2.9729E+01  2.6263E+01  1.5530E+01  1.3752E+01 -3.6162E-02 -1.6341E+01  0.0000E+00  0.0000E+00 -4.8762E+01  1.1509E-03
            -1.0866E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -931.126274179303        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  4.9642E-01  1.5423E-02  3.4519E-02  3.8267E-01  7.4671E+03  1.4705E+00  1.0000E-02  1.0000E-02  9.1357E-01  9.4738E+00
             6.6996E+00
 PARAMETER: -6.0034E-01 -4.0719E+00 -3.2662E+00 -8.6057E-01  9.0183E+00  4.8562E-01 -1.4191E+01 -1.8453E+01  9.6038E-03  2.3485E+00
             2.0021E+00
 GRADIENT:  -7.3883E+00  3.9003E-01 -3.0713E+01  6.7070E+01  2.2525E-04  4.2191E+00  0.0000E+00  0.0000E+00 -1.1120E+01 -1.5006E-08
             1.3538E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -931.904188718426        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  4.3578E-01  1.0000E-02  2.5106E-02  2.9859E-01  1.3953E+04  1.4441E+00  1.0000E-02  1.0000E-02  9.1240E-01  9.4521E+00
             6.6368E+00
 PARAMETER: -7.3061E-01 -4.5602E+00 -3.5846E+00 -1.1087E+00  9.6435E+00  4.6746E-01 -1.5862E+01 -2.0224E+01  8.3237E-03  2.3462E+00
             1.9926E+00
 GRADIENT:  -2.2651E+01  0.0000E+00 -3.1584E+01  5.1068E+01  4.2105E-05  4.8016E+00  0.0000E+00  0.0000E+00 -1.0651E+01  1.2114E-09
             8.8917E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -933.094677678225        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  4.4835E-01  1.0000E-02  2.5725E-02  2.9632E-01  1.3420E+04  1.4250E+00  1.0000E-02  1.0000E-02  9.4930E-01  9.3063E+00
             6.5478E+00
 PARAMETER: -7.0218E-01 -4.5602E+00 -3.5603E+00 -1.1163E+00  9.6045E+00  4.5416E-01 -1.5824E+01 -2.0144E+01  4.7965E-02  2.3307E+00
             1.9791E+00
 GRADIENT:  -3.2227E-02  0.0000E+00  8.4716E-02 -8.2331E-02  4.8803E-05 -7.5463E-03  0.0000E+00  0.0000E+00  2.1336E-02  2.6828E-09
            -3.0014E-04

0ITERATION NO.:   46    OBJECTIVE VALUE:  -933.094677678225        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  4.4835E-01  1.0000E-02  2.5725E-02  2.9632E-01  1.3420E+04  1.4250E+00  1.0000E-02  1.0000E-02  9.4930E-01  9.3063E+00
             6.5478E+00
 PARAMETER: -7.0218E-01 -4.5602E+00 -3.5603E+00 -1.1163E+00  9.6045E+00  4.5416E-01 -1.5824E+01 -2.0144E+01  4.7965E-02  2.3307E+00
             1.9791E+00
 GRADIENT:  -3.2227E-02  0.0000E+00  8.4716E-02 -8.2331E-02  4.8803E-05 -7.5463E-03  0.0000E+00  0.0000E+00  2.1336E-02  2.6828E-09
            -3.0014E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      820
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4641E-04  1.3923E-06  7.1245E-05 -1.7817E-02 -5.3380E-07
 SE:             2.9167E-02  1.2254E-06  2.3136E-04  2.5554E-02  8.1636E-07
 N:                     100         100         100         100         100

 P VAL.:         9.8779E-01  2.5586E-01  7.5813E-01  4.8566E-01  5.1319E-01

 ETASHRINKSD(%)  2.2853E+00  9.9996E+01  9.9225E+01  1.4391E+01  9.9997E+01
 ETASHRINKVR(%)  4.5183E+00  1.0000E+02  9.9994E+01  2.6711E+01  1.0000E+02
 EBVSHRINKSD(%)  2.3400E+00  9.9995E+01  9.9261E+01  1.4221E+01  9.9997E+01
 EBVSHRINKVR(%)  4.6253E+00  1.0000E+02  9.9995E+01  2.6420E+01  1.0000E+02
 RELATIVEINF(%)  2.8576E+00  2.0301E-08  4.0793E-05  5.5434E-01  4.4026E-09
 EPSSHRINKSD(%)  1.8586E+01
 EPSSHRINKVR(%)  3.3717E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -933.09467767822457     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -197.94385111448639     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -933.095       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.48E-01  1.00E-02  2.57E-02  2.96E-01  1.34E+04  1.42E+00  1.00E-02  1.00E-02  9.49E-01  9.31E+00  6.55E+00
 


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
+        2.58E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.54E+04  0.00E+00  8.75E+05
 
 TH 4
+       -1.79E+02  0.00E+00 -8.63E+04  9.86E+03
 
 TH 5
+       -8.84E-07  0.00E+00 -1.57E-06  1.92E-07 -1.48E-12
 
 TH 6
+       -2.22E+00  0.00E+00 -1.73E+02 -2.05E+01  2.22E-07  8.74E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.36E+00  0.00E+00  1.20E+03 -1.38E+02  7.51E-07 -8.65E-01  0.00E+00  0.00E+00  1.16E+02
 
 TH10
+        1.28E-03  0.00E+00 -2.96E-02  2.65E-03 -2.17E-09 -1.28E-03  0.00E+00  0.00E+00  2.29E-02  4.90E-05
 
 TH11
+       -2.31E+01  0.00E+00  2.58E+02 -1.36E+01 -4.75E-09  9.84E-01  0.00E+00  0.00E+00  5.04E+00 -3.49E-05  9.09E+00
 
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
 #CPUT: Total CPU Time in Seconds,       17.024
Stop Time:
Sat Sep 18 15:14:47 CDT 2021
