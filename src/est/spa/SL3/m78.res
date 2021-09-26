Sat Sep 25 11:54:35 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1606.36932860959        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.1061E+01 -1.2334E+02 -6.8299E+01 -9.8165E+01  5.9188E+01 -2.1134E+00 -2.0003E+01  1.5465E+01 -8.0582E+00  1.1187E+01
            -9.4763E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1632.71900188365        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7863E-01  1.1869E+00  1.3339E+00  9.7425E-01  1.1443E+00  1.0016E+00  1.0993E+00  8.8075E-01  1.0131E+00  9.1868E-01
             1.1896E+00
 PARAMETER:  7.8395E-02  2.7131E-01  3.8811E-01  7.3909E-02  2.3482E-01  1.0161E-01  1.9469E-01 -2.6983E-02  1.1298E-01  1.5184E-02
             2.7365E-01
 GRADIENT:   3.0162E+01  3.8294E+00 -9.4918E-01 -3.4559E+00  3.1226E+00 -8.0789E-01  9.5535E-01  3.2119E+00 -3.7206E+00 -1.8541E+01
            -1.1772E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1635.50269233929        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.7603E-01  1.1809E+00  1.6457E+00  9.8716E-01  1.2473E+00  1.0472E+00  8.7443E-01  7.0777E-01  1.1573E+00  1.1375E+00
             1.1776E+00
 PARAMETER:  7.5735E-02  2.6630E-01  5.9816E-01  8.7079E-02  3.2100E-01  1.4607E-01 -3.4183E-02 -2.4564E-01  2.4605E-01  2.2880E-01
             2.6346E-01
 GRADIENT:   2.6572E+01  3.3990E+00  2.1843E-01  2.8043E+00  6.8734E+00  1.7215E+01  1.6488E+00 -7.4225E-02  3.2683E+00 -3.8420E+00
            -1.2912E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1635.98052173736        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      255
 NPARAMETR:  9.6673E-01  1.2678E+00  1.5783E+00  9.3673E-01  1.2589E+00  1.0031E+00  7.9741E-01  6.5861E-01  1.2220E+00  1.1602E+00
             1.1917E+00
 PARAMETER:  6.6165E-02  3.3730E-01  5.5636E-01  3.4641E-02  3.3025E-01  1.0306E-01 -1.2639E-01 -3.1763E-01  3.0047E-01  2.4857E-01
             2.7537E-01
 GRADIENT:  -2.5321E+01  4.6579E+00 -5.0316E-01  9.2330E+00  6.9404E-01 -3.1811E+00  5.3405E-01  2.5548E-02  1.1611E+00 -8.2035E-01
            -6.8948E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1636.73216181031        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      432
 NPARAMETR:  9.7761E-01  1.5274E+00  1.5234E+00  7.5410E-01  1.3420E+00  1.0098E+00  6.0622E-01  7.1532E-01  1.5042E+00  1.2075E+00
             1.2102E+00
 PARAMETER:  7.7356E-02  5.2359E-01  5.2093E-01 -1.8223E-01  3.9417E-01  1.0972E-01 -4.0051E-01 -2.3502E-01  5.0824E-01  2.8852E-01
             2.9079E-01
 GRADIENT:  -3.5377E+00  6.7162E+00  3.8414E-01  5.0608E+00 -6.6671E-01 -4.6290E-01 -1.0762E+00 -1.0493E-01 -1.0029E+00 -2.8299E-01
            -1.4732E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1637.13050136687        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      609
 NPARAMETR:  9.8341E-01  1.7750E+00  1.2261E+00  5.9040E-01  1.3782E+00  1.0131E+00  7.2563E-01  6.1435E-01  1.6761E+00  1.2053E+00
             1.2192E+00
 PARAMETER:  8.3269E-02  6.7383E-01  3.0384E-01 -4.2695E-01  4.2076E-01  1.1300E-01 -2.2071E-01 -3.8720E-01  6.1650E-01  2.8674E-01
             2.9816E-01
 GRADIENT:   7.2589E+00  1.8788E+00  1.1953E+00  1.5163E+00 -7.2660E+00  6.0397E-01  2.5293E+00  6.5178E-01  2.6410E+00  1.6544E+00
             2.5165E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1637.61048259899        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      784
 NPARAMETR:  9.8184E-01  1.9145E+00  1.0303E+00  4.9578E-01  1.4186E+00  1.0125E+00  7.1378E-01  2.3178E-01  1.8052E+00  1.2077E+00
             1.2147E+00
 PARAMETER:  8.1672E-02  7.4945E-01  1.2984E-01 -6.0162E-01  4.4967E-01  1.1241E-01 -2.3718E-01 -1.3620E+00  6.9068E-01  2.8869E-01
             2.9447E-01
 GRADIENT:   2.6801E+00 -7.5071E-01 -1.6177E+00  8.3102E-01  5.2733E+00  1.7548E-01 -9.2914E-01  1.4687E-01 -1.3875E+00 -5.3740E-03
             6.9829E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1637.68672797495        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  9.8071E-01  1.9734E+00  1.0042E+00  4.5286E-01  1.4318E+00  1.0120E+00  7.0201E-01  1.3694E-01  1.9396E+00  1.2144E+00
             1.2150E+00
 PARAMETER:  8.0520E-02  7.7973E-01  1.0420E-01 -6.9217E-01  4.5895E-01  1.1198E-01 -2.5381E-01 -1.8882E+00  7.6246E-01  2.9423E-01
             2.9473E-01
 GRADIENT:   9.2449E-02 -1.9501E+00 -1.8875E-01 -8.1877E-01  4.2727E-01  2.9028E-02 -1.2469E-01  5.0010E-02 -6.1934E-02  8.7685E-02
             1.1588E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1637.71902495307        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1137
 NPARAMETR:  9.8059E-01  1.9446E+00  1.0191E+00  4.7360E-01  1.4208E+00  1.0119E+00  7.0710E-01  2.7079E-02  1.8865E+00  1.2081E+00
             1.2139E+00
 PARAMETER:  8.0398E-02  7.6503E-01  1.1891E-01 -6.4739E-01  4.5124E-01  1.1188E-01 -2.4659E-01 -3.5090E+00  7.3473E-01  2.8903E-01
             2.9387E-01
 GRADIENT:  -1.5730E-02 -7.4466E-02 -2.0725E-02 -6.9263E-03 -6.2209E-02 -8.7424E-03  5.9339E-02  1.9212E-03  3.5749E-02 -1.0396E-02
            -1.3526E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1637.71985108410        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1319
 NPARAMETR:  9.8061E-01  1.9452E+00  1.0194E+00  4.7307E-01  1.4213E+00  1.0120E+00  7.0654E-01  1.0000E-02  1.8882E+00  1.2085E+00
             1.2140E+00
 PARAMETER:  8.0418E-02  7.6537E-01  1.1923E-01 -6.4850E-01  4.5157E-01  1.1190E-01 -2.4738E-01 -4.5598E+00  7.3563E-01  2.8939E-01
             2.9392E-01
 GRADIENT:   2.2555E-02 -1.4615E-01  6.1622E-03 -5.0357E-02 -4.2789E-03  2.6305E-03  1.3863E-02  0.0000E+00  8.9406E-03  2.6103E-03
             1.4441E-03

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1637.71985244237        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1376
 NPARAMETR:  9.8060E-01  1.9453E+00  1.0193E+00  4.7312E-01  1.4213E+00  1.0120E+00  7.0648E-01  1.0000E-02  1.8881E+00  1.2085E+00
             1.2140E+00
 PARAMETER:  8.0411E-02  7.6542E-01  1.1911E-01 -6.4840E-01  4.5158E-01  1.1190E-01 -2.4746E-01 -4.5598E+00  7.3559E-01  2.8938E-01
             2.9392E-01
 GRADIENT:   3.8989E-03  1.6974E-02 -3.5983E-03  1.3515E-02  7.2227E-03  8.4230E-05 -4.2025E-03  0.0000E+00  5.2675E-03 -4.1283E-05
            -2.5498E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1376
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.0276E-04 -3.6292E-02 -9.0668E-05  2.2883E-02 -4.1823E-02
 SE:             2.9704E-02  2.0219E-02  5.0695E-05  2.1912E-02  2.2807E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9187E-01  7.2666E-02  7.3692E-02  2.9632E-01  6.6687E-02

 ETASHRINKSD(%)  4.8850E-01  3.2263E+01  9.9830E+01  2.6593E+01  2.3594E+01
 ETASHRINKVR(%)  9.7461E-01  5.4118E+01  1.0000E+02  4.6115E+01  4.1621E+01
 EBVSHRINKSD(%)  6.3925E-01  2.8963E+01  9.9822E+01  2.9780E+01  2.1965E+01
 EBVSHRINKVR(%)  1.2744E+00  4.9538E+01  1.0000E+02  5.0691E+01  3.9105E+01
 RELATIVEINF(%)  9.8638E+01  4.1843E+00  1.2990E-04  4.4106E+00  2.7618E+01
 EPSSHRINKSD(%)  4.0064E+01
 EPSSHRINKVR(%)  6.4077E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1637.7198524423707     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -902.56902587863249     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1637.720       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.95E+00  1.02E+00  4.73E-01  1.42E+00  1.01E+00  7.06E-01  1.00E-02  1.89E+00  1.21E+00  1.21E+00
 


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
+        1.11E+03
 
 TH 2
+       -9.60E+00  3.12E+02
 
 TH 3
+       -1.31E+00  3.06E+01  3.76E+01
 
 TH 4
+       -1.18E+01  3.71E+02 -3.46E+01  7.05E+02
 
 TH 5
+       -3.79E+00 -7.35E+01 -4.31E+01  4.48E+01  2.31E+02
 
 TH 6
+        3.26E+00 -2.15E+00 -1.02E+00 -4.38E+00 -2.54E+00  1.86E+02
 
 TH 7
+        2.81E+00 -2.27E+01  1.17E+01 -1.45E+01 -1.86E+01  2.59E+00  1.08E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.27E+00 -9.86E+00 -6.60E+00  3.83E+01  3.14E+00 -2.85E-01  2.29E+01  0.00E+00  1.94E+01
 
 TH10
+        4.06E-02 -4.54E+00 -5.45E+00 -3.67E+00 -3.87E+01  3.37E-03  3.40E+00  0.00E+00  1.20E+00  5.89E+01
 
 TH11
+       -1.04E+01 -1.80E+01 -8.61E+00 -1.77E+00  3.26E+00  1.68E+00  6.56E+00  0.00E+00  2.35E+00  1.82E+01  1.66E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.461
Stop Time:
Sat Sep 25 11:55:00 CDT 2021
