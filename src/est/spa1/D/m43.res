Thu Sep 30 03:05:19 CDT 2021
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
$DATA ../../../../data/spa1/D/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17130.4116893522        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.7206E+02  3.3802E+02 -5.6080E+01  3.5141E+02  6.2735E+01 -1.5356E+03 -8.7765E+02 -1.7306E+01 -1.2381E+03 -1.7567E+02
            -3.4090E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -563.783061836404        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1936E+00  1.1591E+00  9.9311E-01  1.7140E+00  1.2163E+00  2.1029E+00  1.4360E+00  9.4793E-01  1.6950E+00  9.5407E-01
             1.4140E+01
 PARAMETER:  2.7699E-01  2.4762E-01  9.3091E-02  6.3882E-01  2.9581E-01  8.4334E-01  4.6184E-01  4.6523E-02  6.2771E-01  5.2981E-02
             2.7490E+00
 GRADIENT:  -3.0509E+01  3.8247E+01 -1.5619E+01  7.4838E+01 -4.1771E+00  5.4134E+01 -3.7020E+00  5.0022E+00 -4.5926E+00  3.0722E+00
             2.1416E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -583.180330021972        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.2206E+00  1.2405E+00  1.9826E+00  1.6531E+00  3.6904E+00  1.6915E+00  2.5670E+00  5.0187E-01  1.8453E+00  8.7270E-01
             1.3643E+01
 PARAMETER:  2.9934E-01  3.1553E-01  7.8440E-01  6.0267E-01  1.4057E+00  6.2564E-01  1.0427E+00 -5.8940E-01  7.1264E-01 -3.6166E-02
             2.7132E+00
 GRADIENT:  -1.3747E+01  2.8798E+01  2.3196E+00  5.2679E+01 -1.2872E+00 -8.7526E+00  8.2947E+00  7.3645E-02  1.9156E+01  3.7354E-02
             2.1273E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -634.826670565913        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0585E+00  1.0850E+00  1.1228E+00  1.0970E+00  1.4981E+01  1.8040E+00  1.9757E+00  2.2642E-02  1.8574E+00  6.3516E+00
             9.3001E+00
 PARAMETER:  1.5683E-01  1.8158E-01  2.1580E-01  1.9262E-01  2.8068E+00  6.8999E-01  7.8093E-01 -3.6880E+00  7.1920E-01  1.9487E+00
             2.3300E+00
 GRADIENT:   1.2649E+01  1.4078E+01  5.5486E+00 -2.1296E+01 -1.0978E+00  2.8622E+00  7.3888E+00  5.0086E-04 -1.0387E+01 -1.5032E-01
            -6.1811E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -670.802254620283        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  8.7147E-01  2.2606E-01  3.0374E-01  1.1782E+00  9.0919E+00  1.3803E+00  1.1966E+00  1.0000E-02  9.8513E-01  3.3041E+00
             8.9669E+00
 PARAMETER: -3.7578E-02 -1.3870E+00 -1.0916E+00  2.6398E-01  2.3074E+00  4.2230E-01  2.7951E-01 -8.0039E+00  8.5018E-02  1.2952E+00
             2.2935E+00
 GRADIENT:  -8.7920E+00  5.9745E+01 -9.6469E+00  3.2269E+01 -1.3526E+01 -6.4284E+01  4.5242E+00  0.0000E+00 -2.2746E+01  1.7566E+00
            -7.6546E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -724.854401290178        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  5.6204E-01  4.1274E-02  5.1930E-02  5.5192E-01  1.4523E+01  1.4565E+00  7.0296E-02  1.0000E-02  5.4420E-01  3.6596E+00
             9.6641E+00
 PARAMETER: -4.7618E-01 -3.0875E+00 -2.8579E+00 -4.9435E-01  2.7757E+00  4.7606E-01 -2.5550E+00 -1.5885E+01 -5.0844E-01  1.3973E+00
             2.3684E+00
 GRADIENT:   4.3949E+01  2.6334E+01 -1.3706E+02  2.6926E+02 -2.7418E+00  5.1330E+00  6.5189E-02  0.0000E+00 -3.6616E+00  2.1036E+00
            -1.3194E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -739.398340843118        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      495
 NPARAMETR:  5.3960E-01  3.6577E-02  4.8976E-02  4.9881E-01  1.6132E+01  1.4740E+00  5.1661E-02  1.0000E-02  5.2334E-01  4.5741E+00
             9.6935E+00
 PARAMETER: -5.1693E-01 -3.2083E+00 -2.9164E+00 -5.9552E-01  2.8808E+00  4.8800E-01 -2.8631E+00 -1.6389E+01 -5.4752E-01  1.6204E+00
             2.3715E+00
 GRADIENT:  -8.4944E+00  1.8025E+01 -9.8820E+01  1.4498E+02 -2.3746E+00  8.5280E+00  1.8885E-02  0.0000E+00 -2.3333E+00  1.7394E+00
            -1.4458E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -756.555429418243        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  3.8798E-01  1.1036E-02  2.1139E-02  2.5060E-01  2.2001E+01  1.2898E+00  1.0000E-02  1.0000E-02  5.6062E-01  5.2282E+00
             9.1220E+00
 PARAMETER: -8.4681E-01 -4.4066E+00 -3.7566E+00 -1.2839E+00  3.1911E+00  3.5449E-01 -5.1113E+00 -2.0355E+01 -4.7872E-01  1.7541E+00
             2.3107E+00
 GRADIENT:  -9.4355E+00  7.1540E-01 -1.8044E+01  2.5185E+01 -2.6973E-02 -1.1470E+01  0.0000E+00  0.0000E+00 -4.9234E+00  1.7361E-03
            -1.5566E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -757.169457418432        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      857
 NPARAMETR:  4.1485E-01  1.0354E-02  2.4141E-02  2.7616E-01  2.8162E+01  1.3544E+00  1.0000E-02  1.0000E-02  6.6593E-01  1.9806E+00
             9.0777E+00
 PARAMETER: -7.7983E-01 -4.4704E+00 -3.6238E+00 -1.1868E+00  3.4380E+00  4.0337E-01 -4.8898E+00 -1.9598E+01 -3.0657E-01  7.8342E-01
             2.3058E+00
 GRADIENT:   8.3026E+00  6.1562E-02 -1.6814E+01  1.6387E+01  3.9152E-02 -5.5503E-01  0.0000E+00  0.0000E+00 -3.2451E-01 -7.5579E-05
            -4.8391E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -757.244104330099        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1033
 NPARAMETR:  4.1279E-01  1.0000E-02  2.4326E-02  2.7590E-01  1.5263E+01  1.3567E+00  1.0000E-02  1.0000E-02  6.6383E-01  2.6132E+00
             9.0808E+00
 PARAMETER: -7.8481E-01 -4.5691E+00 -3.6162E+00 -1.1877E+00  2.8254E+00  4.0509E-01 -4.8898E+00 -1.9598E+01 -3.0973E-01  1.0606E+00
             2.3062E+00
 GRADIENT:  -1.5542E+00  0.0000E+00  6.8510E-02 -6.5061E-01  7.7780E-03 -7.4475E-02  0.0000E+00  0.0000E+00 -3.5089E-01 -6.6942E-05
            -1.0651E+00

0ITERATION NO.:   49    OBJECTIVE VALUE:  -757.250656221401        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1160
 NPARAMETR:  4.1321E-01  1.0000E-02  2.4306E-02  2.7580E-01  1.4322E+01  1.3568E+00  1.0000E-02  1.0000E-02  6.6647E-01  2.6627E+00
             9.0884E+00
 PARAMETER: -7.8381E-01 -4.5593E+00 -3.6170E+00 -1.1881E+00  2.7618E+00  4.0511E-01 -4.8898E+00 -1.9598E+01 -3.0576E-01  1.0793E+00
             2.3070E+00
 GRADIENT:  -6.0696E-01  0.0000E+00 -4.8205E-01 -4.7504E-01  1.9417E-03 -3.1396E-02  0.0000E+00  0.0000E+00  1.6702E-03  1.5156E-04
            -1.3990E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1160
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.5258E-03  1.0277E-06  1.2499E-04 -1.8987E-02 -8.4773E-05
 SE:             2.8874E-02  1.4471E-06  2.7587E-04  2.2099E-02  2.3991E-04
 N:                     100         100         100         100         100

 P VAL.:         9.0281E-01  4.7758E-01  6.5050E-01  3.9022E-01  7.2382E-01

 ETASHRINKSD(%)  3.2690E+00  9.9995E+01  9.9076E+01  2.5967E+01  9.9196E+01
 ETASHRINKVR(%)  6.4311E+00  1.0000E+02  9.9991E+01  4.5191E+01  9.9994E+01
 EBVSHRINKSD(%)  3.2489E+00  9.9994E+01  9.9125E+01  2.7693E+01  9.9217E+01
 EBVSHRINKVR(%)  6.3922E+00  1.0000E+02  9.9992E+01  4.7717E+01  9.9994E+01
 RELATIVEINF(%)  2.1574E+00  4.0503E-08  4.8150E-05  2.8720E-01  2.1815E-04
 EPSSHRINKSD(%)  1.0440E+01
 EPSSHRINKVR(%)  1.9791E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -757.25065622140119     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       161.68787698327151     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.39
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -757.251       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.13E-01  1.00E-02  2.43E-02  2.76E-01  1.43E+01  1.36E+00  1.00E-02  1.00E-02  6.66E-01  2.66E+00  9.09E+00
 


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
+        3.26E+03
 
 TH 2
+        0.00E+00  4.00E+03
 
 TH 3
+       -1.94E+04  0.00E+00  1.29E+06
 
 TH 4
+       -1.72E+02  0.00E+00 -1.33E+05  1.55E+04
 
 TH 5
+        3.22E-01  0.00E+00 -4.95E+00  3.77E-01  4.99E-04
 
 TH 6
+        3.66E+00  0.00E+00 -7.20E+01 -3.83E+01  4.24E-03  8.25E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.83E+01  0.00E+00 -9.67E+00  6.18E+00 -6.98E-03 -2.37E+00  0.00E+00  0.00E+00  8.25E+01
 
 TH10
+       -1.00E-03  0.00E+00 -2.61E-02 -4.90E-03 -5.67E-05  1.68E-03  0.00E+00  0.00E+00  1.19E-03  1.83E-05
 
 TH11
+       -3.02E+01  0.00E+00  5.32E+02 -4.24E+01 -4.43E-03  2.20E+00  0.00E+00  0.00E+00  7.14E+00  1.15E-04  6.07E+00
 
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
 #CPUT: Total CPU Time in Seconds,       29.085
Stop Time:
Thu Sep 30 03:05:49 CDT 2021
