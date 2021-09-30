Wed Sep 29 23:19:20 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat39.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1043.10801769770        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5993E+02  7.7242E+01  8.0787E+01  5.5758E+01  8.0862E+01  5.5295E+01 -1.5364E+01 -2.7166E+01 -1.3143E+01 -6.1473E+01
            -1.9959E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1706.71512791362        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0215E+00  1.0212E+00  1.1699E+00  1.0358E+00  1.1107E+00  7.7614E-01  8.5307E-01  8.2571E-01  8.3547E-01  6.7990E-01
             2.5554E+00
 PARAMETER:  1.2123E-01  1.2101E-01  2.5696E-01  1.3519E-01  2.0498E-01 -1.5343E-01 -5.8909E-02 -9.1515E-02 -7.9764E-02 -2.8582E-01
             1.0382E+00
 GRADIENT:   1.5117E+02  8.0782E+00 -6.8423E+00  9.0405E+00  2.1842E+01 -6.6043E+01 -3.0956E+00  4.8480E+00 -2.4164E+01  4.1841E+00
            -2.3887E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1719.17902229499        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0314E+00  9.5477E-01  9.3792E-01  1.0809E+00  9.2021E-01  9.6179E-01  6.5074E-01  3.4865E-01  1.1218E+00  2.7564E-01
             2.6955E+00
 PARAMETER:  1.3091E-01  5.3716E-02  3.5908E-02  1.7781E-01  1.6844E-02  6.1046E-02 -3.2965E-01 -9.5369E-01  2.1494E-01 -1.1887E+00
             1.0916E+00
 GRADIENT:   1.4273E+02  3.3257E+01  2.4250E+01  4.3016E+01 -3.3915E+01  1.7938E+01  2.5700E+00  1.2749E+00  2.5019E+01  5.1051E-01
             3.2016E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1725.47262254735        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.7588E-01  8.0106E-01  5.2494E-01  1.1142E+00  6.1341E-01  9.0724E-01  1.0481E+00  7.1205E-02  9.1029E-01  2.8138E-01
             2.5108E+00
 PARAMETER:  7.5587E-02 -1.2182E-01 -5.4447E-01  2.0816E-01 -3.8872E-01  2.6544E-03  1.4696E-01 -2.5422E+00  6.0087E-03 -1.1681E+00
             1.0206E+00
 GRADIENT:  -9.4998E+00  1.9986E+01 -1.6141E+01  6.2480E+01  3.7587E+01 -4.7985E+00  2.0765E+00  9.5335E-02 -1.3417E+00  1.3740E+00
             7.4222E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1730.93920768515        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      405
 NPARAMETR:  1.0043E+00  5.3402E-01  4.1950E-01  1.2164E+00  4.4652E-01  9.1315E-01  1.2256E+00  3.3001E-02  8.5553E-01  2.0706E-01
             2.4797E+00
 PARAMETER:  1.0426E-01 -5.2733E-01 -7.6868E-01  2.9590E-01 -7.0627E-01  9.1418E-03  3.0347E-01 -3.3112E+00 -5.6030E-02 -1.4748E+00
             1.0081E+00
 GRADIENT:   5.0528E+00  1.5163E+01 -4.4200E+00  5.8249E+01 -4.5746E+00 -6.3345E+00 -9.5484E+00 -3.7296E-03 -1.3526E+01 -2.6414E+00
            -3.8679E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1742.38328962096        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      581
 NPARAMETR:  1.0115E+00  3.4403E-01  2.9841E-01  1.1431E+00  3.1821E-01  9.2111E-01  1.7876E+00  1.0000E-02  1.0528E+00  4.9812E-01
             2.2036E+00
 PARAMETER:  1.1146E-01 -9.6702E-01 -1.1093E+00  2.3378E-01 -1.0450E+00  1.7826E-02  6.8090E-01 -5.9774E+00  1.5145E-01 -5.9692E-01
             8.9009E-01
 GRADIENT:   3.8816E+01  6.0459E+00 -2.8148E+00 -5.0604E+01  1.0659E+01 -3.0872E+00  1.4865E+01  0.0000E+00  1.8597E+01  6.0723E+00
            -1.0630E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1746.00393579649        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      756
 NPARAMETR:  9.9907E-01  3.2722E-01  3.2379E-01  1.2042E+00  3.3027E-01  9.2435E-01  1.4555E+00  1.0000E-02  9.7337E-01  5.0773E-01
             2.2614E+00
 PARAMETER:  9.9072E-02 -1.0171E+00 -1.0277E+00  2.8579E-01 -1.0078E+00  2.1332E-02  4.7536E-01 -5.5808E+00  7.3013E-02 -5.7781E-01
             9.1598E-01
 GRADIENT:   5.5743E+00  3.0649E+00  8.7268E-01  3.8256E+00 -3.7924E+00 -1.1909E+00 -2.0145E-01  0.0000E+00 -2.4893E-02  1.4073E+00
             1.1849E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1746.61422144878        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      933
 NPARAMETR:  9.8471E-01  1.8187E-01  3.9914E-01  1.3164E+00  3.5579E-01  9.2090E-01  2.1356E+00  1.0000E-02  9.1323E-01  5.6916E-01
             2.3061E+00
 PARAMETER:  8.4595E-02 -1.6045E+00 -8.1844E-01  3.7491E-01 -9.3343E-01  1.7601E-02  8.5876E-01 -6.3098E+00  9.2364E-03 -4.6360E-01
             9.3556E-01
 GRADIENT:  -1.0191E+01  7.8074E+00  1.0886E+01  1.0050E+01 -2.0429E+01  5.7904E-01  7.4404E+00  0.0000E+00 -1.8368E+00  1.1432E+00
             5.4947E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1748.65695565564        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1109
 NPARAMETR:  9.8131E-01  6.1232E-02  3.9651E-01  1.3474E+00  3.4519E-01  9.2210E-01  3.5051E+00  1.0000E-02  8.9757E-01  5.4154E-01
             2.2827E+00
 PARAMETER:  8.1135E-02 -2.6931E+00 -8.2506E-01  3.9818E-01 -9.6365E-01  1.8896E-02  1.3542E+00 -8.3033E+00 -8.0674E-03 -5.1334E-01
             9.2534E-01
 GRADIENT:   1.8901E+00 -1.9754E-02 -4.2621E+00 -1.3041E+01  9.6566E+00  2.4617E+00 -4.2416E-01  0.0000E+00  2.5826E-02 -1.3875E+00
            -5.8849E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1749.12736254403        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1284
 NPARAMETR:  9.7733E-01  1.0000E-02  4.0635E-01  1.3821E+00  3.4504E-01  9.1531E-01  9.6260E+00  1.0000E-02  8.8458E-01  5.5037E-01
             2.2978E+00
 PARAMETER:  7.7066E-02 -4.7513E+00 -8.0054E-01  4.2361E-01 -9.6409E-01  1.1509E-02  2.3645E+00 -1.2036E+01 -2.2641E-02 -4.9716E-01
             9.3193E-01
 GRADIENT:  -1.0065E+00  0.0000E+00  1.7625E+00  2.6178E+00 -3.2198E+00  3.6243E-01 -5.3152E-02  0.0000E+00  2.4215E-02 -7.9943E-02
            -1.3345E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1749.13593284420        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1450             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7816E-01  1.0000E-02  4.0605E-01  1.3795E+00  3.4497E-01  9.1429E-01  1.0519E+01  1.0000E-02  8.8461E-01  5.5091E-01
             2.2977E+00
 PARAMETER:  7.7917E-02 -4.7513E+00 -8.0129E-01  4.2173E-01 -9.6430E-01  1.0395E-02  2.4532E+00 -1.2036E+01 -2.2612E-02 -4.9618E-01
             9.3191E-01
 GRADIENT:   6.9329E+01  0.0000E+00  1.8152E+01  1.2563E+02  5.7410E+01  5.0361E+00  1.4157E-01  0.0000E+00  3.4201E+00  1.0748E+00
             9.3161E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1749.13593284420        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:     1508
 NPARAMETR:  9.7816E-01  1.0000E-02  4.0605E-01  1.3795E+00  3.4497E-01  9.1429E-01  1.0519E+01  1.0000E-02  8.8461E-01  5.5091E-01
             2.2977E+00
 PARAMETER:  7.7917E-02 -4.7513E+00 -8.0129E-01  4.2173E-01 -9.6430E-01  1.0395E-02  2.4532E+00 -1.2036E+01 -2.2612E-02 -4.9618E-01
             9.3191E-01
 GRADIENT:   1.2938E+00  0.0000E+00  1.4707E+00 -1.8570E+00 -1.5081E+00 -1.0944E-03 -1.0893E-02  0.0000E+00  8.7503E-02  2.3327E-02
            -2.7592E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1508
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3283E-05  9.6243E-04  9.3282E-05 -6.8930E-03 -4.6816E-03
 SE:             2.9308E-02  1.8581E-03  2.6683E-04  2.8063E-02  2.0446E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9964E-01  6.0448E-01  7.2664E-01  8.0597E-01  8.1889E-01

 ETASHRINKSD(%)  1.8156E+00  9.3775E+01  9.9106E+01  5.9849E+00  3.1503E+01
 ETASHRINKVR(%)  3.5982E+00  9.9613E+01  9.9992E+01  1.1612E+01  5.3081E+01
 EBVSHRINKSD(%)  1.8586E+00  9.4427E+01  9.9081E+01  5.4993E+00  3.1464E+01
 EBVSHRINKVR(%)  3.6826E+00  9.9689E+01  9.9992E+01  1.0696E+01  5.3029E+01
 RELATIVEINF(%)  8.0725E+01  2.7897E-02  4.0387E-04  1.4395E+01  1.9520E+00
 EPSSHRINKSD(%)  2.5721E+01
 EPSSHRINKVR(%)  4.4827E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1749.1359328442020     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -830.19739963952929     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1749.136       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  1.00E-02  4.06E-01  1.38E+00  3.45E-01  9.14E-01  1.05E+01  1.00E-02  8.85E-01  5.51E-01  2.30E+00
 


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
+        1.34E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.07E+01  0.00E+00  5.55E+03
 
 TH 4
+       -2.27E+01  0.00E+00 -3.09E+02  6.67E+02
 
 TH 5
+        1.21E+02  0.00E+00 -8.44E+03 -3.13E+02  1.45E+04
 
 TH 6
+        8.81E-01  0.00E+00  1.19E+01 -7.95E+00 -4.38E+00  2.20E+02
 
 TH 7
+        5.69E-03  0.00E+00 -2.14E-02 -2.69E-02  9.14E-02  2.66E-03  2.04E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.77E+00  0.00E+00  8.44E+01 -8.18E+00  7.35E+00  2.20E+00  1.52E-02  0.00E+00  2.02E+02
 
 TH10
+       -6.16E+00  0.00E+00 -7.05E+01  5.70E+00  2.69E+01  1.31E+00 -4.15E-03  0.00E+00  3.94E+00  1.53E+02
 
 TH11
+       -1.61E+01  0.00E+00 -1.79E+01 -8.99E+00  1.12E+01  2.04E+00 -6.91E-03  0.00E+00  5.69E+00  3.31E+01  9.19E+01
 
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
 #CPUT: Total CPU Time in Seconds,       30.198
Stop Time:
Wed Sep 29 23:20:01 CDT 2021
