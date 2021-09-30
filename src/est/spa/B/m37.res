Wed Sep 29 11:11:01 CDT 2021
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
$DATA ../../../../data/spa/B/dat37.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1651.73712820860        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3756E+02  1.1399E+01  2.3665E+01  4.1928E+01 -1.1438E+01  6.4625E+01 -2.6997E+01 -1.0483E+01  1.0868E+01 -9.6677E+00
            -8.5605E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1659.40920497055        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0130E+00  1.1896E+00  8.6209E-01  9.1461E-01  1.0603E+00  9.3696E-01  1.5616E+00  1.1414E+00  9.0057E-01  1.1342E+00
             1.0638E+00
 PARAMETER:  1.1287E-01  2.7359E-01 -4.8395E-02  1.0747E-02  1.5857E-01  3.4883E-02  5.4574E-01  2.3228E-01 -4.7309E-03  2.2592E-01
             1.6184E-01
 GRADIENT:   1.2991E+01  1.6106E+01 -1.9384E+01  1.9868E+01  5.2075E+00 -1.8703E+01  3.0586E+01  7.2984E+00  1.8348E+01  1.6708E+01
             2.7145E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1660.73634348429        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0102E+00  1.0591E+00  1.0533E+00  9.9664E-01  1.0943E+00  9.4615E-01  1.7040E+00  1.2778E+00  7.8678E-01  1.0697E+00
             1.1096E+00
 PARAMETER:  1.1011E-01  1.5745E-01  1.5190E-01  9.6637E-02  1.9007E-01  4.4649E-02  6.3298E-01  3.4518E-01 -1.3981E-01  1.6733E-01
             2.0398E-01
 GRADIENT:   9.7494E+00  1.6328E+01 -8.6034E+00  8.5457E+00  1.5925E+01 -1.3522E+01  2.6620E+01  3.7960E+00  9.7591E+00  3.3659E+00
             3.7836E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.49774734803        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      561             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0056E+00  9.1344E-01  1.1774E+00  1.0844E+00  1.0736E+00  9.6355E-01  1.7326E+00  1.2836E+00  7.1440E-01  1.0677E+00
             1.0437E+00
 PARAMETER:  1.0557E-01  9.4588E-03  2.6335E-01  1.8101E-01  1.7103E-01  6.2867E-02  6.4963E-01  3.4965E-01 -2.3632E-01  1.6548E-01
             1.4276E-01
 GRADIENT:   4.3826E+02  3.7117E+01  8.5795E+00  1.8172E+02  1.3718E+01  4.4574E+01  4.8559E+01  3.6071E-01  1.4682E+01  2.2071E+00
             1.3945E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.79071126645        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      719
 NPARAMETR:  1.0001E+00  9.0215E-01  1.1775E+00  1.0844E+00  1.0736E+00  9.6331E-01  1.6917E+00  1.2836E+00  7.0795E-01  1.0641E+00
             1.0283E+00
 PARAMETER:  1.0012E-01 -2.9758E-03  2.6336E-01  1.8102E-01  1.7102E-01  6.2619E-02  6.2573E-01  3.4964E-01 -2.4539E-01  1.6217E-01
             1.2791E-01
 GRADIENT:  -9.3391E+00  4.3428E+00  4.4011E+00  5.5794E+00  6.7380E+00 -5.8197E+00  3.2133E-01 -9.6505E-01 -6.2298E-01  2.7884E-01
             6.8459E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.97274298781        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      896            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0040E+00  8.9416E-01  1.1741E+00  1.0841E+00  1.0715E+00  9.7791E-01  1.6945E+00  1.2892E+00  7.1190E-01  1.0658E+00
             1.0106E+00
 PARAMETER:  1.0395E-01 -1.1876E-02  2.6051E-01  1.8078E-01  1.6902E-01  7.7665E-02  6.2739E-01  3.5404E-01 -2.3982E-01  1.6374E-01
             1.1058E-01
 GRADIENT:   4.5688E+02  2.6764E+01  8.0050E+00  1.8094E+02  1.7727E+01  5.3920E+01  4.3042E+01  5.1893E-01  1.4783E+01  1.7963E+00
             1.2459E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.57670229966        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0045E+00  8.8806E-01  1.0524E+00  1.0841E+00  9.9680E-01  9.7744E-01  1.6966E+00  1.2073E+00  7.1181E-01  9.6997E-01
             1.0106E+00
 PARAMETER:  1.0453E-01 -1.8721E-02  1.5111E-01  1.8077E-01  9.6791E-02  7.7185E-02  6.2865E-01  2.8839E-01 -2.3995E-01  6.9507E-02
             1.1059E-01
 GRADIENT:  -2.6891E-01  1.9584E+00  2.4622E+00  2.5934E+00 -2.3593E+00 -1.9527E-01 -3.0934E-01  1.7035E+00 -2.2171E-01  1.4337E-01
             1.5803E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.81379468704        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  1.0053E+00  8.7441E-01  9.5115E-01  1.0841E+00  9.4465E-01  9.7845E-01  1.7161E+00  1.0102E+00  7.1367E-01  9.3764E-01
             1.0106E+00
 PARAMETER:  1.0527E-01 -3.4202E-02  4.9912E-02  1.8077E-01  4.3061E-02  7.8212E-02  6.4007E-01  1.1014E-01 -2.3734E-01  3.5608E-02
             1.1059E-01
 GRADIENT:   4.1538E-03 -7.9611E-02 -3.1000E-03  2.2735E+00  3.1410E-02  7.3439E-06  5.3613E-02  4.7332E-03 -2.0064E-02 -1.3915E-02
             1.8724E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1665.82870145864        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1427
 NPARAMETR:  1.0075E+00  8.8473E-01  9.4552E-01  1.0787E+00  9.4678E-01  9.8038E-01  1.7001E+00  1.0112E+00  7.1798E-01  9.4037E-01
             1.0054E+00
 PARAMETER:  1.0748E-01 -2.2471E-02  4.3975E-02  1.7576E-01  4.5310E-02  8.0186E-02  6.3067E-01  1.1111E-01 -2.3132E-01  3.8520E-02
             1.0542E-01
 GRADIENT:   4.9469E+00  2.0043E-01 -9.9968E-01  4.2412E+00  1.1492E+00  6.9178E-01  1.1701E-03  2.1151E-01  2.9602E-01  2.0841E-01
            -1.1682E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1665.84006593326        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  1.0058E+00  8.8649E-01  9.4489E-01  1.0758E+00  9.4685E-01  9.7889E-01  1.6997E+00  1.0090E+00  7.1612E-01  9.3880E-01
             1.0059E+00
 PARAMETER:  1.0578E-01 -2.0485E-02  4.3312E-02  1.7307E-01  4.5387E-02  7.8665E-02  6.3045E-01  1.0899E-01 -2.3391E-01  3.6846E-02
             1.0585E-01
 GRADIENT:   9.7786E-01 -5.8124E-01  2.1878E-01  5.5126E-01  1.0313E-03  1.0367E-01  2.4619E-01 -1.6567E-03 -4.4650E-03  2.6761E-03
             3.6999E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1665.88279470007        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1783
 NPARAMETR:  1.0035E+00  9.6373E-01  9.0209E-01  1.0325E+00  9.5515E-01  9.7792E-01  1.5875E+00  9.9902E-01  7.3267E-01  9.3255E-01
             1.0052E+00
 PARAMETER:  1.0345E-01  6.3053E-02 -3.0436E-03  1.3202E-01  5.4110E-02  7.7669E-02  5.6219E-01  9.9021E-02 -2.1107E-01  3.0163E-02
             1.0516E-01
 GRADIENT:  -6.6253E+00  3.1640E+00  4.1526E-01  6.7653E+00 -5.2863E-01 -6.9995E-01 -9.4321E-01  5.8120E-03 -4.0294E-01 -3.3842E-01
            -2.0602E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1665.93537878081        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1958
 NPARAMETR:  1.0070E+00  1.0220E+00  8.6212E-01  9.9181E-01  9.6303E-01  9.8034E-01  1.5193E+00  9.7415E-01  7.4898E-01  9.3618E-01
             1.0050E+00
 PARAMETER:  1.0700E-01  1.2179E-01 -4.8364E-02  9.1777E-02  6.2334E-02  8.0147E-02  5.1824E-01  7.3807E-02 -1.8905E-01  3.4050E-02
             1.0494E-01
 GRADIENT:   2.1983E-01 -2.8039E-04  2.5210E-01 -2.4373E-01  3.2433E-01  2.4560E-02  7.9415E-03  3.8492E-02  6.5641E-03  2.3834E-02
             8.0519E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1665.93976996585        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2137
 NPARAMETR:  1.0072E+00  1.0220E+00  8.4874E-01  9.9090E-01  9.5489E-01  9.8046E-01  1.5187E+00  9.5255E-01  7.4968E-01  9.2867E-01
             1.0049E+00
 PARAMETER:  1.0714E-01  1.2172E-01 -6.4008E-02  9.0860E-02  5.3842E-02  8.0267E-02  5.1784E-01  5.1388E-02 -1.8810E-01  2.6002E-02
             1.0488E-01
 GRADIENT:   2.3577E-01  8.2134E-03  7.3657E-02 -1.4215E-01 -1.6728E-01  2.4014E-02  6.0358E-02  1.7310E-02  2.6210E-02  3.8085E-02
             3.7288E-02

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1665.93976996585        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     2162
 NPARAMETR:  1.0072E+00  1.0220E+00  8.4874E-01  9.9090E-01  9.5489E-01  9.8046E-01  1.5187E+00  9.5255E-01  7.4968E-01  9.2867E-01
             1.0049E+00
 PARAMETER:  1.0714E-01  1.2172E-01 -6.4008E-02  9.0860E-02  5.3842E-02  8.0267E-02  5.1784E-01  5.1388E-02 -1.8810E-01  2.6002E-02
             1.0488E-01
 GRADIENT:  -1.1277E+00  2.7005E-02  7.5043E-02  6.9148E-02 -1.5708E-01 -1.7592E-01 -1.8547E-01  9.1579E-03 -1.0045E-02  3.7262E-02
             3.7705E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2162
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2980E-03  5.4341E-03 -3.4989E-02 -1.2176E-02 -2.5135E-02
 SE:             2.9857E-02  2.4212E-02  1.4088E-02  2.0301E-02  2.0411E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6532E-01  8.2242E-01  1.3002E-02  5.4864E-01  2.1816E-01

 ETASHRINKSD(%)  1.0000E-10  1.8886E+01  5.2805E+01  3.1991E+01  3.1622E+01
 ETASHRINKVR(%)  1.0000E-10  3.4205E+01  7.7726E+01  5.3747E+01  5.3244E+01
 EBVSHRINKSD(%)  4.6470E-01  1.8425E+01  5.6842E+01  3.2869E+01  2.8960E+01
 EBVSHRINKVR(%)  9.2724E-01  3.3454E+01  8.1374E+01  5.4934E+01  4.9534E+01
 RELATIVEINF(%)  9.8563E+01  6.4602E+00  2.2530E+00  3.8033E+00  8.5538E+00
 EPSSHRINKSD(%)  4.4944E+01
 EPSSHRINKVR(%)  6.9689E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.9397699658512     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.78894340211298     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.77
 Elapsed covariance  time in seconds:     6.00
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.940       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.02E+00  8.49E-01  9.91E-01  9.55E-01  9.80E-01  1.52E+00  9.53E-01  7.50E-01  9.29E-01  1.00E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.00E-02  4.62E-01  2.97E-01  2.91E-01  1.24E-01  7.37E-02  5.32E-01  3.63E-01  1.45E-01  2.06E-01  6.74E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.98E-04
 
 TH 2
+        1.86E-03  2.13E-01
 
 TH 3
+       -1.15E-03 -1.19E-01  8.85E-02
 
 TH 4
+       -1.19E-03 -1.33E-01  7.67E-02  8.45E-02
 
 TH 5
+        2.71E-04  3.83E-02 -9.53E-03 -2.30E-02  1.53E-02
 
 TH 6
+        1.27E-04  6.90E-03 -2.78E-03 -4.19E-03  2.17E-03  5.43E-03
 
 TH 7
+       -2.25E-03 -2.35E-01  1.38E-01  1.48E-01 -3.78E-02 -7.96E-03  2.83E-01
 
 TH 8
+       -5.10E-04 -6.10E-02  6.26E-02  4.12E-02 -1.11E-03 -2.48E-03  7.28E-02  1.32E-01
 
 TH 9
+        8.19E-04  5.52E-02 -3.25E-02 -3.52E-02  8.87E-03  1.86E-03 -6.29E-02 -2.06E-02  2.10E-02
 
 TH10
+       -2.19E-04  1.90E-02  2.87E-03 -1.11E-02  1.54E-02  1.47E-03 -1.46E-02 -2.02E-02  4.22E-03  4.23E-02
 
 TH11
+       -2.77E-04  9.88E-04  1.22E-04 -5.38E-04  5.03E-04  1.13E-04 -1.92E-03  3.11E-03 -1.03E-04 -8.80E-04  4.54E-03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.00E-02
 
 TH 2
+        1.35E-01  4.62E-01
 
 TH 3
+       -1.29E-01 -8.65E-01  2.97E-01
 
 TH 4
+       -1.36E-01 -9.91E-01  8.88E-01  2.91E-01
 
 TH 5
+        7.30E-02  6.71E-01 -2.59E-01 -6.39E-01  1.24E-01
 
 TH 6
+        5.73E-02  2.03E-01 -1.27E-01 -1.96E-01  2.38E-01  7.37E-02
 
 TH 7
+       -1.41E-01 -9.56E-01  8.74E-01  9.59E-01 -5.74E-01 -2.03E-01  5.32E-01
 
 TH 8
+       -4.69E-02 -3.64E-01  5.80E-01  3.91E-01 -2.46E-02 -9.26E-02  3.77E-01  3.63E-01
 
 TH 9
+        1.88E-01  8.24E-01 -7.55E-01 -8.36E-01  4.95E-01  1.74E-01 -8.15E-01 -3.92E-01  1.45E-01
 
 TH10
+       -3.56E-02  2.00E-01  4.69E-02 -1.86E-01  6.07E-01  9.72E-02 -1.33E-01 -2.70E-01  1.41E-01  2.06E-01
 
 TH11
+       -1.37E-01  3.18E-02  6.07E-03 -2.75E-02  6.03E-02  2.27E-02 -5.37E-02  1.27E-01 -1.06E-02 -6.35E-02  6.74E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.19E+03
 
 TH 2
+        1.66E+01  3.49E+02
 
 TH 3
+        2.22E+01  8.89E+01  3.78E+02
 
 TH 4
+       -3.05E+01  3.63E+02 -3.01E+02  1.03E+03
 
 TH 5
+       -6.05E+01 -2.20E+02 -4.12E+02  2.83E+02  8.18E+02
 
 TH 6
+       -1.18E+01  3.86E+00  1.85E-01 -1.53E+01 -5.00E+01  2.01E+02
 
 TH 7
+        9.16E+00  2.97E+01 -3.99E+00 -3.97E+01 -2.09E+01  8.30E+00  4.98E+01
 
 TH 8
+       -2.36E+00 -1.02E+01 -5.49E+01  3.12E+01  1.56E+01  6.62E+00  1.66E+00  2.36E+01
 
 TH 9
+       -5.91E+01 -1.34E+00 -4.30E+01  1.02E+02  4.83E+01 -2.72E+00  7.47E+00  1.10E+01  1.71E+02
 
 TH10
+        2.13E+01  1.96E+01 -1.71E+01  1.73E+01 -1.02E+02  1.17E+01  1.35E+00  2.12E+01  3.56E+00  6.78E+01
 
 TH11
+        8.01E+01  1.35E+01  1.34E+01 -1.39E+01 -3.88E+01 -1.65E+00  1.20E+01 -5.76E+00  4.92E+00  9.81E+00  2.36E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       34.829
Stop Time:
Wed Sep 29 11:11:37 CDT 2021
