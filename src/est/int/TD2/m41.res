Wed Sep 29 07:23:02 CDT 2021
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
$DATA ../../../../data/int/TD2/dat41.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3197.66682133497        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7258E+02 -8.2163E+01  1.4790E+02  1.4754E+02  2.0637E+02  5.0234E+01 -5.2901E+01 -1.7340E+02 -5.0255E+01 -3.8771E+01
            -1.1524E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3443.25672369622        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0467E+00  1.0795E+00  7.6041E-01  9.4844E-01  9.0051E-01  9.6012E-01  1.0646E+00  1.2022E+00  1.0673E+00  1.0761E+00
             1.3952E+00
 PARAMETER:  1.4566E-01  1.7654E-01 -1.7390E-01  4.7068E-02 -4.7980E-03  5.9302E-02  1.6262E-01  2.8419E-01  1.6509E-01  1.7333E-01
             4.3306E-01
 GRADIENT:   3.5968E+02  3.5427E+01 -1.5596E+01  1.2859E+01  4.9535E+01  1.4299E+01  4.4244E+00 -3.7187E+01 -1.4357E+00 -8.9981E+00
            -8.7327E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3464.14320527562        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0304E+00  1.1818E+00  8.2772E-01  9.2953E-01  1.0526E+00  1.0192E+00  8.6157E-01  2.3673E+00  1.0790E+00  1.6467E+00
             1.3500E+00
 PARAMETER:  1.2999E-01  2.6700E-01 -8.9084E-02  2.6922E-02  1.5129E-01  1.1900E-01 -4.8998E-02  9.6175E-01  1.7605E-01  5.9878E-01
             4.0008E-01
 GRADIENT:   3.0917E+02  4.7989E+01 -2.3136E+01  5.9272E+01  6.9489E+01  4.3087E+01  1.0704E+01  1.2903E+01  2.2548E+01  6.1507E+01
            -5.3910E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3482.76093604974        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      293
 NPARAMETR:  1.0366E+00  1.3891E+00  9.3653E-01  8.4253E-01  1.1966E+00  9.7649E-01  6.3793E-01  2.4451E+00  1.1037E+00  1.6009E+00
             1.4064E+00
 PARAMETER:  1.3595E-01  4.2863E-01  3.4424E-02 -7.1342E-02  2.7949E-01  7.6204E-02 -3.4953E-01  9.9410E-01  1.9867E-01  5.7057E-01
             4.4103E-01
 GRADIENT:   7.8810E+00 -2.3045E+01 -1.0599E+01  4.6067E+01  2.5238E+01  2.7028E+00 -4.3305E+00 -8.1012E+00 -1.6510E+00  6.2095E-01
             1.2396E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3486.25076982723        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      470
 NPARAMETR:  1.0350E+00  1.5143E+00  1.0141E+00  7.5943E-01  1.2851E+00  9.7098E-01  6.3111E-01  2.5905E+00  1.1648E+00  1.6449E+00
             1.4170E+00
 PARAMETER:  1.3441E-01  5.1494E-01  1.1403E-01 -1.7519E-01  3.5082E-01  7.0551E-02 -3.6027E-01  1.0519E+00  2.5257E-01  5.9768E-01
             4.4853E-01
 GRADIENT:   3.7958E+00  4.3585E+00  1.1962E+01  2.1252E+01  1.8621E+01  4.6316E-01 -4.4893E+00 -1.2776E+01  2.6383E-01  1.1544E+00
             7.8701E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3486.87595720925        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  1.0341E+00  1.5191E+00  9.5335E-01  7.5678E-01  1.2778E+00  9.6924E-01  6.3284E-01  2.6158E+00  1.1699E+00  1.6146E+00
             1.4160E+00
 PARAMETER:  1.3351E-01  5.1809E-01  5.2226E-02 -1.7869E-01  3.4516E-01  6.8757E-02 -3.5754E-01  1.0616E+00  2.5692E-01  5.7907E-01
             4.4781E-01
 GRADIENT:   1.6506E+00  3.2177E+00  2.3498E-01  2.3768E+01  1.7089E+01 -2.7472E-01 -4.0037E+00 -6.8839E+00  2.4655E-01 -4.4783E+00
             6.4450E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3487.40087707084        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  1.0309E+00  1.5190E+00  9.5245E-01  7.5180E-01  1.2497E+00  9.6887E-01  6.3740E-01  2.7109E+00  1.1689E+00  1.6427E+00
             1.4072E+00
 PARAMETER:  1.3047E-01  5.1806E-01  5.1285E-02 -1.8529E-01  3.2291E-01  6.8376E-02 -3.5036E-01  1.0973E+00  2.5608E-01  5.9635E-01
             4.4160E-01
 GRADIENT:   2.8626E+02  3.6529E+02  1.0034E+00  5.8212E+01  7.2071E+01  1.8954E+01  5.5284E+00  8.8471E+00  8.2773E+00  3.3380E+01
             5.1519E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3487.63440147022        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  1.0307E+00  1.5229E+00  9.5462E-01  7.3851E-01  1.2588E+00  9.7048E-01  6.8811E-01  2.7004E+00  1.1443E+00  1.6414E+00
             1.4090E+00
 PARAMETER:  1.3024E-01  5.2062E-01  5.3558E-02 -2.0312E-01  3.3017E-01  7.0040E-02 -2.7380E-01  1.0934E+00  2.3477E-01  5.9554E-01
             4.4286E-01
 GRADIENT:  -5.4879E+00 -1.6553E+01  3.9976E-01 -4.6368E+00 -1.4174E+00  2.5161E-01  3.1170E-01 -3.5295E+00  1.5393E+00  5.0056E-01
            -4.1225E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3487.67804781454        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1169
 NPARAMETR:  1.0308E+00  1.5241E+00  9.5998E-01  7.4310E-01  1.2614E+00  9.7039E-01  7.1105E-01  2.7103E+00  1.1138E+00  1.6397E+00
             1.4088E+00
 PARAMETER:  1.3038E-01  5.2141E-01  5.9160E-02 -1.9693E-01  3.3226E-01  6.9947E-02 -2.4101E-01  1.0971E+00  2.0780E-01  5.9449E-01
             4.4273E-01
 GRADIENT:  -5.2574E+00 -8.7142E+00  4.8430E-01  1.0626E+00 -5.0102E-01  1.9507E-01  5.3507E-01 -4.0996E+00  1.3970E-01  1.5129E-01
            -1.1046E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3487.90168487934        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1346
 NPARAMETR:  1.0327E+00  1.5517E+00  9.5868E-01  7.2567E-01  1.2876E+00  9.6935E-01  7.1270E-01  2.8680E+00  1.1142E+00  1.6578E+00
             1.4090E+00
 PARAMETER:  1.3217E-01  5.3937E-01  5.7806E-02 -2.2067E-01  3.5277E-01  6.8867E-02 -2.3869E-01  1.1536E+00  2.0810E-01  6.0551E-01
             4.4286E-01
 GRADIENT:  -1.0869E+00 -9.5740E+00 -4.6636E-01 -1.0308E+00  1.2946E+00 -2.4619E-01 -2.7606E-01  1.7708E+00  5.3971E-03  3.0849E-01
             9.7683E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3487.93566457517        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  1.0342E+00  1.5595E+00  9.5622E-01  7.2352E-01  1.2842E+00  9.7033E-01  7.2157E-01  2.8445E+00  1.1078E+00  1.6608E+00
             1.4098E+00
 PARAMETER:  1.3358E-01  5.4439E-01  5.5237E-02 -2.2363E-01  3.5015E-01  6.9882E-02 -2.2633E-01  1.1454E+00  2.0236E-01  6.0731E-01
             4.4342E-01
 GRADIENT:   2.2035E+00 -1.9336E+00  1.5416E-01  9.3408E-01 -2.5594E+00  1.5107E-01  2.2390E-01  4.2637E-01 -1.6376E-01  7.8338E-01
            -4.3862E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3487.94147483723        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1688
 NPARAMETR:  1.0336E+00  1.5605E+00  9.5060E-01  7.2157E-01  1.2842E+00  9.7011E-01  7.1083E-01  2.8216E+00  1.1213E+00  1.6610E+00
             1.4105E+00
 PARAMETER:  1.3309E-01  5.4503E-01  4.9337E-02 -2.2633E-01  3.5017E-01  6.9657E-02 -2.4132E-01  1.1373E+00  2.1448E-01  6.0741E-01
             4.4394E-01
 GRADIENT:   1.0517E+03 -2.6253E+02 -8.5458E-02 -6.1673E+02  3.9629E+02  2.4132E-02  3.2610E-02 -2.4205E-01  5.3760E-03 -2.2670E+02
             3.1130E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1688
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.1399E-04 -4.7637E-02 -3.0154E-02  3.5561E-02 -3.1351E-02
 SE:             2.9829E-02  1.9827E-02  2.0609E-02  2.4708E-02  2.6568E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7556E-01  1.6275E-02  1.4343E-01  1.5008E-01  2.3798E-01

 ETASHRINKSD(%)  6.9594E-02  3.3578E+01  3.0958E+01  1.7225E+01  1.0995E+01
 ETASHRINKVR(%)  1.3914E-01  5.5881E+01  5.2331E+01  3.1483E+01  2.0782E+01
 EBVSHRINKSD(%)  5.1128E-01  3.2348E+01  3.0613E+01  2.0730E+01  8.8212E+00
 EBVSHRINKVR(%)  1.0199E+00  5.4232E+01  5.1854E+01  3.7162E+01  1.6864E+01
 RELATIVEINF(%)  9.8972E+01  1.0462E+01  4.2229E+01  1.5445E+01  4.2893E+01
 EPSSHRINKSD(%)  2.1014E+01
 EPSSHRINKVR(%)  3.7611E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3487.9414748372287     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1833.8521150688180     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    49.93
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3487.941       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.56E+00  9.51E-01  7.22E-01  1.28E+00  9.70E-01  7.11E-01  2.82E+00  1.12E+00  1.66E+00  1.41E+00
 


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
+        1.86E+05
 
 TH 2
+       -7.35E+00  5.37E+03
 
 TH 3
+       -1.67E+03  3.00E+02  1.00E+02
 
 TH 4
+        9.93E+02  3.43E+02  1.35E+03  1.31E+05
 
 TH 5
+       -4.70E+02 -1.03E+01  8.13E+04  4.75E+02  1.74E+04
 
 TH 6
+        2.88E-01 -1.41E+00  1.92E-01 -8.74E-01 -9.80E-01  2.07E+02
 
 TH 7
+       -8.03E+02  1.14E+02 -2.61E+00  6.74E+02 -4.53E+04 -3.27E-02  8.25E+01
 
 TH 8
+       -1.10E+02  1.25E+01 -9.58E+00  9.20E+01 -2.44E+03 -6.85E-02 -1.27E+00  8.57E+00
 
 TH 9
+       -9.03E+01  9.61E+00  3.08E+00  1.10E+02 -3.22E+04 -1.31E-01  5.33E+01  2.72E+00  5.97E+01
 
 TH10
+        6.73E+02 -1.18E+02  2.28E+02 -5.58E+02 -7.59E+03  3.43E-01  1.10E+02  1.27E+01  1.54E+01  3.40E+03
 
 TH11
+       -1.12E+03  1.41E+02 -3.66E+02  9.62E+02  1.22E+04  1.80E+00 -1.66E+02 -1.55E+01 -1.55E+01 -5.38E+03  9.21E+03
 
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
 #CPUT: Total CPU Time in Seconds,       65.590
Stop Time:
Wed Sep 29 07:24:10 CDT 2021
