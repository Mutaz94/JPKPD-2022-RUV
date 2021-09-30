Thu Sep 30 03:37:19 CDT 2021
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
$DATA ../../../../data/spa1/D/dat77.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   24324.0705818989        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5174E+02  4.7709E+02 -4.9065E+01  2.9146E+02  2.3232E+02 -2.3349E+03 -9.6301E+02 -9.1731E+01 -1.7218E+03 -5.2500E+02
            -4.6696E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -575.860474395569        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1465E+00  9.2885E-01  7.9483E-01  1.8611E+00  1.6977E+00  2.8192E+00  1.2920E+00  9.3219E-01  1.5221E+00  9.8390E-01
             1.3685E+01
 PARAMETER:  2.3670E-01  2.6195E-02 -1.2963E-01  7.2116E-01  6.2926E-01  1.1365E+00  3.5620E-01  2.9783E-02  5.2009E-01  8.3767E-02
             2.7163E+00
 GRADIENT:  -4.8534E+01  4.7155E+01 -1.6882E+01  7.3502E+01 -7.8597E+00  8.3169E+01 -3.3388E+00  6.0119E+00 -4.7154E+01  4.4981E-01
             7.3909E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -624.868333637254        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.2365E+00  7.3833E-01  6.9937E-01  1.8923E+00  5.2483E+00  2.0651E+00  1.0790E+00  2.1991E-01  2.2425E+00  1.9692E+00
             1.3489E+01
 PARAMETER:  3.1229E-01 -2.0337E-01 -2.5757E-01  7.3777E-01  1.7579E+00  8.2519E-01  1.7605E-01 -1.4145E+00  9.0759E-01  7.7761E-01
             2.7019E+00
 GRADIENT:  -4.0144E+01  4.2511E+01 -1.6968E+01  6.5378E+01 -2.1477E+01 -3.7030E+01  4.1566E+00  3.0826E-01  2.7802E+01  1.0646E+00
             1.1772E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -710.315774705310        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0917E+00  1.5713E-01  4.5380E-01  1.4667E+00  2.3664E+01  1.8303E+00  4.4529E-01  2.5262E+00  1.2029E+00  5.1325E+00
             1.2387E+01
 PARAMETER:  1.8773E-01 -1.7507E+00 -6.9010E-01  4.8300E-01  3.2639E+00  7.0450E-01 -7.0904E-01  1.0267E+00  2.8473E-01  1.7356E+00
             2.6167E+00
 GRADIENT:  -1.3634E+01  9.9121E+00  2.6228E+01  5.6301E+01 -1.6633E-01 -1.9306E+01  1.4281E-01 -1.4528E+01  2.7990E+01  1.0223E-02
             8.4779E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -760.337573694215        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  7.0102E-01  1.5674E-02  6.2546E-02  5.9371E-01  2.2662E+02  1.3818E+00  5.7062E-01  1.9091E+00  1.8193E-01  2.1379E+00
             1.1064E+01
 PARAMETER: -2.5522E-01 -4.0558E+00 -2.6719E+00 -4.2136E-01  5.5233E+00  4.2339E-01 -4.6103E-01  7.4661E-01 -1.6041E+00  8.5981E-01
             2.5037E+00
 GRADIENT:   2.3512E+01 -3.4772E-01 -6.1354E+01  8.1926E+01  1.4079E-02 -6.2336E+01  1.3416E-04  1.9705E+01  1.2576E+00 -8.7585E-06
            -1.3574E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -785.544962966940        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      502
 NPARAMETR:  5.6402E-01  1.3009E-02  3.6560E-02  3.7741E-01  4.4025E+02  1.7613E+00  4.4512E-01  9.9316E-01  1.1145E-01  2.1157E+00
             1.0568E+01
 PARAMETER: -4.7267E-01 -4.2421E+00 -3.2088E+00 -8.7443E-01  6.1873E+00  6.6603E-01 -7.0940E-01  9.3134E-02 -2.0942E+00  8.4937E-01
             2.4578E+00
 GRADIENT:   1.3613E+01 -1.0756E-01 -3.3661E+01  4.3016E+01  6.3672E-03 -2.7215E+00  1.5655E-04 -5.0893E+00  5.1134E-02 -1.1079E-06
            -3.8915E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -786.943039668655        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      686             RESET HESSIAN, TYPE I
 NPARAMETR:  5.4328E-01  1.4463E-02  3.4372E-02  3.5516E-01  4.1449E+02  1.7622E+00  2.0786E-01  1.0065E+00  4.6115E-02  2.2634E+00
             1.0905E+01
 PARAMETER: -5.1013E-01 -4.1362E+00 -3.2705E+00 -9.3517E-01  6.1271E+00  6.6654E-01 -1.4709E+00  1.0652E-01 -2.9766E+00  9.1688E-01
             2.4893E+00
 GRADIENT:   3.7938E+01 -1.2925E-01  5.2482E+01  2.9176E+01  2.6868E-03  1.2644E+01  1.1064E-04  1.2752E-01  4.2918E-02 -2.1865E-07
             2.4101E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -787.017638986181        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      863
 NPARAMETR:  5.4443E-01  1.9393E-02  3.4513E-02  3.5446E-01  2.0654E+02  1.7587E+00  1.2839E-02  1.0154E+00  2.7550E-02  2.2312E+00
             1.0895E+01
 PARAMETER: -5.0801E-01 -3.8429E+00 -3.2664E+00 -9.3715E-01  5.4305E+00  6.6460E-01 -4.2553E+00  1.1531E-01 -3.4918E+00  9.0253E-01
             2.4883E+00
 GRADIENT:  -6.5061E-02  2.9024E-02  2.9408E-01 -1.3162E+00  3.7179E-03  4.3113E-02  2.0417E-05  3.6295E-02  9.8030E-03  8.5054E-06
            -3.8631E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -787.064263827691        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1042
 NPARAMETR:  5.4590E-01  1.8076E-02  3.4596E-02  3.5534E-01  1.2654E+01  1.7613E+00  1.0000E-02  1.0172E+00  1.5577E-02  3.5125E-01
             1.0896E+01
 PARAMETER: -5.0532E-01 -3.9132E+00 -3.2640E+00 -9.3468E-01  2.6380E+00  6.6605E-01 -5.7985E+00  1.1707E-01 -4.0619E+00 -9.4626E-01
             2.4884E+00
 GRADIENT:  -8.2774E-01  1.6610E-01  1.7824E+00 -3.4193E+00 -1.2118E-02  1.8739E-01  0.0000E+00  3.3852E-01  3.2663E-03  5.2218E-04
             7.5705E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -787.073097720177        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  5.4689E-01  1.7361E-02  3.4521E-02  3.5593E-01  1.1923E+01  1.7620E+00  1.0000E-02  1.0124E+00  1.0000E-02  3.1046E-02
             1.0891E+01
 PARAMETER: -5.0351E-01 -3.9535E+00 -3.2662E+00 -9.3302E-01  2.5785E+00  6.6648E-01 -5.5413E+00  1.1228E-01 -5.0506E+00 -3.3723E+00
             2.4879E+00
 GRADIENT:   5.0683E-01 -5.8058E-03 -2.7710E+00  2.4481E+00  1.3944E-02  9.4894E-02  0.0000E+00  1.3334E-02  0.0000E+00  4.3500E-06
            -9.6400E-01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -787.073097720177        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1246
 NPARAMETR:  5.4689E-01  1.7361E-02  3.4521E-02  3.5593E-01  1.1923E+01  1.7620E+00  1.0000E-02  1.0124E+00  1.0000E-02  3.1046E-02
             1.0891E+01
 PARAMETER: -5.0351E-01 -3.9535E+00 -3.2662E+00 -9.3302E-01  2.5785E+00  6.6648E-01 -5.5413E+00  1.1228E-01 -5.0506E+00 -3.3723E+00
             2.4879E+00
 GRADIENT:   5.0683E-01 -5.8058E-03 -2.7710E+00  2.4481E+00  1.3944E-02  9.4894E-02  0.0000E+00  1.3334E-02  0.0000E+00  4.3500E-06
            -9.6400E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1246
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2049E-03 -1.1642E-05  1.2036E-02 -3.7767E-04  5.1109E-06
 SE:             2.9038E-02  3.6227E-06  2.1918E-02  2.9493E-04  4.3106E-06
 N:                     100         100         100         100         100

 P VAL.:         8.3080E-01  1.3114E-03  5.8291E-01  2.0035E-01  2.3575E-01

 ETASHRINKSD(%)  2.7178E+00  9.9988E+01  2.6572E+01  9.9012E+01  9.9986E+01
 ETASHRINKVR(%)  5.3617E+00  1.0000E+02  4.6084E+01  9.9990E+01  1.0000E+02
 EBVSHRINKSD(%)  2.9143E+00  9.9976E+01  2.7310E+01  9.9016E+01  9.9980E+01
 EBVSHRINKVR(%)  5.7436E+00  1.0000E+02  4.7161E+01  9.9990E+01  1.0000E+02
 RELATIVEINF(%)  8.1943E+00  3.7024E-07  7.0556E-01  9.7933E-05  1.4614E-07
 EPSSHRINKSD(%)  8.4329E+00
 EPSSHRINKVR(%)  1.6155E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -787.07309772017675     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       131.86543548449595     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -787.073       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.47E-01  1.74E-02  3.45E-02  3.56E-01  1.19E+01  1.76E+00  1.00E-02  1.01E+00  1.00E-02  3.10E-02  1.09E+01
 


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
+       -3.64E+02  6.74E+03
 
 TH 3
+       -2.59E+03  2.27E+03  4.38E+05
 
 TH 4
+       -4.82E+02 -4.26E+02 -5.49E+04  7.69E+03
 
 TH 5
+        2.09E-01 -2.05E+00 -4.80E+00  6.36E-01  2.40E-03
 
 TH 6
+        1.57E+00  2.80E+01  1.46E+02 -3.49E+01  3.65E-03  5.29E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -7.70E+00 -5.02E+01 -3.58E+02 -3.04E+00 -3.03E-04  2.12E+00  0.00E+00  4.10E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        8.48E-03  7.15E-02  1.63E-02 -3.73E-03 -4.12E-05  1.10E-03  0.00E+00 -3.62E-03  0.00E+00 -1.12E-02
 
 TH11
+       -1.37E+01  1.40E+00  2.56E+02 -2.74E+01 -4.63E-03  5.60E-01  0.00E+00  3.16E+00  0.00E+00 -7.69E-05  3.88E+00
 
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
 #CPUT: Total CPU Time in Seconds,       34.331
Stop Time:
Thu Sep 30 03:37:55 CDT 2021
