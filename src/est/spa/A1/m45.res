Wed Sep 29 12:08:50 CDT 2021
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
$DATA ../../../../data/spa/A1/dat45.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1218.02254992957        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1196E+02  1.2654E+01  8.8197E+00  2.3790E+01  6.2757E+01 -3.4569E+00 -7.7379E+00 -3.2852E+00 -1.3060E+01 -1.3849E+01
            -7.4649E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1425.87319933585        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0931E+00  1.0220E+00  1.0655E+00  1.0485E+00  9.8936E-01  1.4437E+00  9.6138E-01  9.3862E-01  1.0268E+00  8.5063E-01
             2.2201E+00
 PARAMETER:  1.8906E-01  1.2178E-01  1.6345E-01  1.4732E-01  8.9306E-02  4.6721E-01  6.0617E-02  3.6652E-02  1.2642E-01 -6.1774E-02
             8.9753E-01
 GRADIENT:   2.8879E+02  1.6506E+01  3.7766E+00  2.5572E+01 -6.0629E+00  1.1395E+02  5.5537E+00  2.7142E+00  6.7058E+00  4.6731E+00
             3.2175E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1437.72867211826        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0272E+00  1.0467E+00  6.2123E-01  9.8917E-01  7.6603E-01  1.2191E+00  1.0146E+00  4.0783E-01  1.0997E+00  5.3485E-01
             2.1061E+00
 PARAMETER:  1.2681E-01  1.4561E-01 -3.7605E-01  8.9109E-02 -1.6653E-01  2.9811E-01  1.1451E-01 -7.9689E-01  1.9503E-01 -5.2577E-01
             8.4486E-01
 GRADIENT:   1.9629E+02  8.7505E+00 -7.2774E+00  2.2352E+01  1.5775E+01  5.9118E+01  7.7736E+00  3.3326E-01  2.2907E+01 -4.8142E+00
             3.1403E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1443.80453250459        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.7833E-01  9.2847E-01  9.5761E-01  1.0877E+00  9.0096E-01  1.1524E+00  8.7372E-01  3.6583E-01  1.0035E+00  7.9840E-01
             2.1161E+00
 PARAMETER:  7.8092E-02  2.5780E-02  5.6688E-02  1.8407E-01 -4.2924E-03  2.4186E-01 -3.5000E-02 -9.0558E-01  1.0352E-01 -1.2514E-01
             8.4957E-01
 GRADIENT:   1.1379E+01  5.9281E+00  2.1330E+00 -5.5831E-01 -2.7191E+00  1.4998E+00 -1.4954E+00  3.1051E-01 -2.2701E+00 -1.3336E+00
            -5.7742E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1445.74374134477        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  9.6793E-01  5.8660E-01  9.7884E-01  1.2972E+00  7.7640E-01  1.1581E+00  1.2728E+00  6.0125E-02  8.8048E-01  7.9449E-01
             2.1235E+00
 PARAMETER:  6.7409E-02 -4.3341E-01  7.8609E-02  3.6019E-01 -1.5308E-01  2.4679E-01  3.4125E-01 -2.7113E+00 -2.7283E-02 -1.3006E-01
             8.5307E-01
 GRADIENT:  -1.6695E+00  9.6900E+00  1.3273E+01  3.8498E+00 -2.2216E+01  4.3901E+00  6.2362E-01  1.9492E-02  1.5573E+00  9.0611E-01
             2.2925E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1447.28076257577        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  9.6441E-01  2.6729E-01  1.0222E+00  1.4941E+00  7.1775E-01  1.1392E+00  1.7818E+00  1.0000E-02  8.0014E-01  8.0252E-01
             2.1150E+00
 PARAMETER:  6.3762E-02 -1.2194E+00  1.2199E-01  5.0151E-01 -2.3164E-01  2.3034E-01  6.7764E-01 -6.2229E+00 -1.2297E-01 -1.2000E-01
             8.4907E-01
 GRADIENT:   8.3574E-01  4.1873E+00  4.0924E+00  2.3609E+01 -7.6704E+00 -6.9472E-01 -5.8679E-01  0.0000E+00 -6.5790E-01 -7.0370E-01
             9.5967E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1447.94996260919        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  9.6143E-01  1.0612E-01  9.6713E-01  1.5686E+00  6.5503E-01  1.1356E+00  3.2444E+00  1.0000E-02  7.7167E-01  8.1395E-01
             2.0761E+00
 PARAMETER:  6.0669E-02 -2.1432E+00  6.6580E-02  5.5019E-01 -3.2308E-01  2.2719E-01  1.2769E+00 -1.1063E+01 -1.5919E-01 -1.0586E-01
             8.3048E-01
 GRADIENT:   1.8166E+00  1.1926E+00  4.9063E+00  1.0290E+01 -1.0275E+01 -1.2824E+00 -7.3901E-02  0.0000E+00  5.4442E-02  1.1938E+00
            -2.5742E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1448.10153907080        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1000
 NPARAMETR:  9.5873E-01  5.1280E-02  1.0414E+00  1.6076E+00  6.8077E-01  1.1382E+00  4.8943E+00  1.0000E-02  7.5924E-01  8.2220E-01
             2.0989E+00
 PARAMETER:  5.7856E-02 -2.8705E+00  1.4061E-01  5.7474E-01 -2.8453E-01  2.2941E-01  1.6881E+00 -1.4878E+01 -1.7543E-01 -9.5768E-02
             8.4141E-01
 GRADIENT:  -4.2827E-01  6.2171E-01  6.7839E-01  5.1133E+00 -1.0766E-01  1.7549E-01  5.6167E-01  0.0000E+00  5.7681E-01 -2.3739E-01
             4.7671E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1448.22216410677        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1180
 NPARAMETR:  9.5783E-01  1.0662E-02  1.0141E+00  1.6263E+00  6.5924E-01  1.1379E+00  9.9765E+00  1.0000E-02  7.5219E-01  8.2204E-01
             2.0866E+00
 PARAMETER:  5.6912E-02 -4.4411E+00  1.1399E-01  5.8630E-01 -3.1666E-01  2.2915E-01  2.4002E+00 -2.3236E+01 -1.8476E-01 -9.5967E-02
             8.3552E-01
 GRADIENT:  -3.7901E-01  8.6239E-02  1.3055E+00  6.8704E+00 -2.9397E+00  1.4645E-01 -1.1876E-02  0.0000E+00 -5.8497E-01  2.3901E-01
            -1.0551E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1448.25207375928        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  9.5807E-01  1.0000E-02  1.0065E+00  1.6198E+00  6.5625E-01  1.1379E+00  1.0208E+01  1.0000E-02  7.5382E-01  8.1768E-01
             2.0892E+00
 PARAMETER:  5.7166E-02 -4.5833E+00  1.0651E-01  5.8228E-01 -3.2121E-01  2.2920E-01  2.4232E+00 -2.3397E+01 -1.8261E-01 -1.0129E-01
             8.3678E-01
 GRADIENT:   3.2544E-01  0.0000E+00  7.6285E-01 -5.2962E+00 -7.3308E-01  2.4609E-01 -1.3631E-02  0.0000E+00  8.7801E-02  4.0268E-02
            -4.6816E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1448.25262491887        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:     1469
 NPARAMETR:  9.5807E-01  1.0000E-02  1.0059E+00  1.6196E+00  6.5595E-01  1.1380E+00  1.0301E+01  1.0000E-02  7.5360E-01  8.1751E-01
             2.0892E+00
 PARAMETER:  5.7168E-02 -4.5833E+00  1.0577E-01  5.8225E-01 -3.2099E-01  2.2921E-01  2.4568E+00 -2.3397E+01 -1.8271E-01 -1.0148E-01
             8.3679E-01
 GRADIENT:   5.2147E-03  0.0000E+00 -2.8630E-02  1.2118E-01  4.4347E-01 -8.6603E-03  1.8680E-02  0.0000E+00  2.2062E-02  2.1643E-04
             4.0302E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1469
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5287E-04  3.2853E-04  4.0688E-05 -1.0677E-02 -2.1651E-02
 SE:             2.9494E-02  1.8032E-03  1.6372E-04  2.7451E-02  2.0742E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9586E-01  8.5543E-01  8.0373E-01  6.9731E-01  2.9658E-01

 ETASHRINKSD(%)  1.1900E+00  9.3959E+01  9.9452E+01  8.0346E+00  3.0510E+01
 ETASHRINKVR(%)  2.3659E+00  9.9635E+01  9.9997E+01  1.5424E+01  5.1712E+01
 EBVSHRINKSD(%)  1.2830E+00  9.4411E+01  9.9419E+01  7.6668E+00  3.0502E+01
 EBVSHRINKVR(%)  2.5495E+00  9.9688E+01  9.9997E+01  1.4746E+01  5.1700E+01
 RELATIVEINF(%)  9.2372E+01  9.3266E-03  1.9583E-04  3.6394E+00  2.1560E+00
 EPSSHRINKSD(%)  3.3690E+01
 EPSSHRINKVR(%)  5.6030E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1448.2526249188654     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -713.10179835512724     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1448.253       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.58E-01  1.00E-02  1.01E+00  1.62E+00  6.56E-01  1.14E+00  1.06E+01  1.00E-02  7.54E-01  8.18E-01  2.09E+00
 


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
+        9.11E+02
 
 TH 2
+        0.00E+00  2.51E+02
 
 TH 3
+       -6.46E-01  0.00E+00  3.72E+02
 
 TH 4
+       -2.46E+01  0.00E+00 -4.28E+01  6.37E+02
 
 TH 5
+        1.65E+01  0.00E+00 -7.86E+02 -1.17E+02  1.84E+03
 
 TH 6
+       -1.33E-01  0.00E+00  3.32E+00 -7.11E+00 -1.88E+00  1.46E+02
 
 TH 7
+        2.67E-03  0.00E+00 -1.70E-02 -5.36E-02  4.87E-02  3.22E-04  2.12E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.26E+00  0.00E+00  2.02E+01 -1.18E+01 -9.14E+00  1.41E+00  2.06E-03  0.00E+00  2.55E+02
 
 TH10
+       -3.87E+00  0.00E+00  5.30E+00 -2.65E+00 -5.80E+01  3.95E-01 -3.81E-02  0.00E+00  1.70E+00  7.63E+01
 
 TH11
+       -9.38E+00  0.00E+00 -9.60E+00 -9.84E+00 -5.46E+00  1.28E+00 -6.41E-03  0.00E+00  1.23E+01  2.48E+01  6.50E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.268
Stop Time:
Wed Sep 29 12:09:15 CDT 2021
