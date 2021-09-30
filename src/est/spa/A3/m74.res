Wed Sep 29 13:46:35 CDT 2021
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
$DATA ../../../../data/spa/A3/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -621.560264539465        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5003E+02  2.4568E+00  8.7164E+01 -1.0695E+02  5.5712E+01  4.5008E+01 -4.5857E+01 -3.8842E+01 -1.0467E+02 -4.0949E+01
            -1.8010E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1322.09589014962        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6847E-01  1.0104E+00  9.6488E-01  1.0895E+00  1.0014E+00  8.3747E-01  1.0201E+00  1.0178E+00  1.1034E+00  8.5432E-01
             2.8012E+00
 PARAMETER:  6.7961E-02  1.1031E-01  6.4252E-02  1.8575E-01  1.0140E-01 -7.7365E-02  1.1990E-01  1.1765E-01  1.9842E-01 -5.7455E-02
             1.1300E+00
 GRADIENT:  -8.4249E+00 -1.4843E+01 -4.5252E+00 -1.7312E+01  1.3619E+01 -2.4323E+01  4.5694E-01 -1.2656E+00 -4.0063E+00  1.5361E+01
            -5.0031E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1331.66708613754        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.7005E-01  1.1533E+00  8.2629E-01  1.0441E+00  9.9997E-01  9.6240E-01  6.7257E-01  1.1676E+00  1.1853E+00  3.0249E-01
             3.0136E+00
 PARAMETER:  6.9592E-02  2.4266E-01 -9.0813E-02  1.4312E-01  9.9975E-02  6.1679E-02 -2.9665E-01  2.5496E-01  2.7001E-01 -1.0957E+00
             1.2031E+00
 GRADIENT:  -1.4084E+01  9.8380E+00 -1.8159E+01  2.3615E+01  2.6513E+01  2.2164E+01 -3.1202E+00  1.9635E+00 -8.4381E+00  1.8817E+00
            -1.3177E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1334.54519330751        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.7871E-01  1.0276E+00  6.3349E-01  1.0724E+00  7.7029E-01  8.9316E-01  1.0287E+00  7.1541E-01  1.1056E+00  2.2742E-01
             3.0001E+00
 PARAMETER:  7.8477E-02  1.2723E-01 -3.5651E-01  1.6990E-01 -1.6098E-01 -1.2995E-02  1.2832E-01 -2.3490E-01  2.0038E-01 -1.3810E+00
             1.1987E+00
 GRADIENT:   2.8729E+00  5.4271E+00  2.5494E+00  4.5720E+00 -2.5617E+00 -8.8089E-01 -3.6380E-01 -1.0912E-01  9.8069E-01  1.4644E+00
             1.8724E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1335.85035090138        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  9.8981E-01  8.6908E-01  7.3990E-01  1.1909E+00  7.7625E-01  9.0318E-01  1.1218E+00  7.6716E-01  1.0374E+00  1.3082E-01
             3.0761E+00
 PARAMETER:  8.9758E-02 -4.0323E-02 -2.0124E-01  2.7473E-01 -1.5328E-01 -1.8302E-03  2.1495E-01 -1.6506E-01  1.3668E-01 -1.9339E+00
             1.2237E+00
 GRADIENT:  -7.8497E+00  3.4869E-01 -2.8116E+00  1.7650E+00  3.3853E+00  1.2422E+00 -7.9708E-01 -1.2287E-01 -1.6555E-01  4.0197E-01
             1.6381E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1337.35217989007        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  9.9381E-01  5.6601E-01  6.8916E-01  1.3571E+00  6.2895E-01  8.9875E-01  1.8143E+00  7.2100E-01  9.2921E-01  1.7609E-02
             3.0126E+00
 PARAMETER:  9.3787E-02 -4.6914E-01 -2.7229E-01  4.0535E-01 -3.6371E-01 -6.7462E-03  6.9572E-01 -2.2712E-01  2.6580E-02 -3.9393E+00
             1.2028E+00
 GRADIENT:   8.7982E+00  7.8777E+00 -5.4047E-01  8.2147E+00 -1.1747E+00 -5.9250E-01  3.2348E+00 -2.9539E-01  1.3628E+00  7.1164E-03
             3.8587E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1339.22004392166        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  9.8726E-01  3.1048E-01  5.7101E-01  1.4464E+00  4.8821E-01  8.9797E-01  2.5676E+00  7.5111E-01  8.9071E-01  1.0000E-02
             2.8793E+00
 PARAMETER:  8.7174E-02 -1.0696E+00 -4.6034E-01  4.6910E-01 -6.1702E-01 -7.6220E-03  1.0430E+00 -1.8620E-01 -1.5737E-02 -7.4690E+00
             1.1575E+00
 GRADIENT:   4.5253E+00  4.7036E+00  9.9290E-01  2.2108E+00 -3.3747E+00 -1.8874E+00  1.7364E+00  4.3717E-01 -4.7088E+00  0.0000E+00
            -7.6657E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1341.62202255524        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  9.8131E-01  1.6365E-01  4.6961E-01  1.4667E+00  4.0135E-01  9.0789E-01  3.5975E+00  7.4699E-01  9.3337E-01  1.0000E-02
             2.8529E+00
 PARAMETER:  8.1138E-02 -1.7100E+00 -6.5585E-01  4.8301E-01 -8.1292E-01  3.3721E-03  1.3802E+00 -1.9170E-01  3.1050E-02 -1.1909E+01
             1.1483E+00
 GRADIENT:  -1.0138E+00  5.8257E-01  9.6399E+00  6.1817E+00 -1.0532E+01  1.2174E+00 -4.7092E+00  2.3742E+00  3.8242E+00  0.0000E+00
             6.0484E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1344.89324921106        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1117
 NPARAMETR:  9.7673E-01  5.0139E-02  4.0035E-01  1.4678E+00  3.5034E-01  8.9688E-01  6.6802E+00  7.0602E-01  9.0786E-01  1.0000E-02
             2.8006E+00
 PARAMETER:  7.6453E-02 -2.8930E+00 -8.1541E-01  4.8373E-01 -9.4886E-01 -8.8283E-03  1.9992E+00 -2.4811E-01  3.3366E-03 -2.0137E+01
             1.1298E+00
 GRADIENT:   6.4026E+00  7.8153E+00 -1.3715E+00 -2.2987E+00  4.0725E+00 -2.1989E+00  1.4461E+01 -8.5967E-01 -1.1752E+01  0.0000E+00
            -1.1970E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1345.14126468365        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1292
 NPARAMETR:  9.7470E-01  4.4251E-02  3.8981E-01  1.4563E+00  3.4190E-01  9.0469E-01  7.1384E+00  7.2526E-01  9.2997E-01  1.0000E-02
             2.7936E+00
 PARAMETER:  7.4182E-02 -3.0159E+00 -8.4196E-01  4.7559E-01 -9.7269E-01 -2.2368E-04  2.0668E+00 -2.2346E-01  2.8242E-02 -2.1081E+01
             1.1266E+00
 GRADIENT:  -1.7299E+00  3.8541E+01  1.8001E+00 -9.9763E+01  7.4836E+01 -1.6637E-01  5.5030E+01 -5.6674E-01  1.0350E+00  0.0000E+00
            -4.1726E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1292
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2476E-03  2.4701E-02 -2.6151E-03 -1.5449E-02  1.3686E-04
 SE:             2.8794E-02  1.0134E-02  1.7660E-02  2.6266E-02  3.6127E-04
 N:                     100         100         100         100         100

 P VAL.:         9.3778E-01  1.4790E-02  8.8228E-01  5.5641E-01  7.0480E-01

 ETASHRINKSD(%)  3.5366E+00  6.6051E+01  4.0835E+01  1.2004E+01  9.8790E+01
 ETASHRINKVR(%)  6.9480E+00  8.8474E+01  6.4996E+01  2.2567E+01  9.9985E+01
 EBVSHRINKSD(%)  3.3220E+00  7.4662E+01  4.0519E+01  1.0619E+01  9.8872E+01
 EBVSHRINKVR(%)  6.5336E+00  9.3580E+01  6.4620E+01  2.0111E+01  9.9987E+01
 RELATIVEINF(%)  9.1899E+01  4.1950E+00  1.4897E+00  3.1955E+01  5.2072E-04
 EPSSHRINKSD(%)  3.0450E+01
 EPSSHRINKVR(%)  5.1628E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1345.1412646836507     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -609.99043811991248     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.25
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
 





 #OBJV:********************************************    -1345.141       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  4.43E-02  3.90E-01  1.46E+00  3.42E-01  9.05E-01  7.15E+00  7.24E-01  9.31E-01  1.00E-02  2.79E+00
 


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
+        1.35E+03
 
 TH 2
+       -2.01E+02  4.31E+05
 
 TH 3
+       -2.16E+02 -3.91E+04  2.29E+04
 
 TH 4
+       -1.81E+04 -1.87E+04  2.62E+03  6.86E+03
 
 TH 5
+        2.45E+02  3.75E+04 -2.52E+04 -1.36E+04  3.95E+04
 
 TH 6
+       -1.72E+00 -9.57E+01  1.34E+02  3.24E+01  3.39E+01  2.17E+02
 
 TH 7
+        2.58E+00  4.36E+03 -4.90E+02 -2.31E+02  4.77E+02 -2.31E+00  4.64E+01
 
 TH 8
+       -1.39E+01 -1.52E+02  6.43E+01  1.46E+01  1.00E+02  9.56E+00 -2.99E+00  6.75E+01
 
 TH 9
+       -1.19E+01  1.36E+02  1.62E+02  1.89E+04  8.15E+01  6.35E+00 -1.72E+03  8.60E+00  1.98E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.49E+01  1.19E+02  1.71E+03  1.39E+03 -1.64E+03  1.03E+01 -5.01E+01  1.78E+01 -2.04E+00  0.00E+00  1.64E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.866
Stop Time:
Wed Sep 29 13:47:01 CDT 2021
