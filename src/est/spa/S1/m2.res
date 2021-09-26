Sat Sep 25 09:39:17 CDT 2021
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
$DATA ../../../../data/spa/S1/dat2.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1676.08274169437        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.7048E+01 -5.8580E+01 -2.4303E+01 -6.0900E+01 -1.1824E+01  7.4280E+00  1.0873E+01  1.5158E+01  3.1347E+01  1.7346E+01
            -2.4823E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1681.10886918814        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8850E-01  1.1415E+00  1.2520E+00  9.8372E-01  1.1748E+00  9.7271E-01  8.9443E-01  8.4709E-01  7.3495E-01  9.2012E-01
             1.1538E+00
 PARAMETER:  8.8433E-02  2.3233E-01  3.2478E-01  8.3583E-02  2.6111E-01  7.2328E-02 -1.1565E-02 -6.5949E-02 -2.0796E-01  1.6746E-02
             2.4306E-01
 GRADIENT:   4.9805E+01  3.3269E+01  1.3672E+01  2.3560E+01  1.4533E+01 -2.3127E+00 -1.1526E+01 -3.8668E+00 -2.4239E+01 -2.3150E+01
             1.1084E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1684.95765702239        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9096E-01  9.4400E-01  1.0742E+00  1.1049E+00  1.0311E+00  9.8851E-01  8.7448E-01  2.5505E-01  8.7096E-01  9.8882E-01
             1.1390E+00
 PARAMETER:  9.0919E-02  4.2372E-02  1.7154E-01  1.9979E-01  1.3059E-01  8.8440E-02 -3.4131E-02 -1.2663E+00 -3.8154E-02  8.8758E-02
             2.3013E-01
 GRADIENT:   5.6088E+01  3.7419E+00 -1.5934E+01  4.8683E+01  3.0988E+01  3.5242E+00 -3.0709E+00 -1.9230E-01  7.1507E+00 -5.7567E+00
             2.1543E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1686.88763085135        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.6755E-01  1.0139E+00  9.6187E-01  1.0388E+00  9.9212E-01  9.7967E-01  1.0149E+00  3.1167E-01  8.0972E-01  9.4606E-01
             1.0743E+00
 PARAMETER:  6.7013E-02  1.1379E-01  6.1121E-02  1.3806E-01  9.2092E-02  7.9460E-02  1.1475E-01 -1.0658E+00 -1.1107E-01  4.4555E-02
             1.7163E-01
 GRADIENT:   4.3275E+00  1.0696E+00 -5.0386E+00  1.0070E+01  7.0725E+00 -7.7886E-02 -7.5510E-02  2.0021E-01 -1.3827E+00  3.8808E-01
             8.5637E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1686.88814859514        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.6716E-01  1.0184E+00  9.6051E-01  1.0352E+00  9.9277E-01  9.7978E-01  1.0091E+00  3.0939E-01  8.1400E-01  9.4620E-01
             1.0738E+00
 PARAMETER:  6.6613E-02  1.1820E-01  5.9714E-02  1.3462E-01  9.2746E-02  7.9571E-02  1.0906E-01 -1.0732E+00 -1.0580E-01  4.4697E-02
             1.7124E-01
 GRADIENT:   3.4310E+00  7.7933E-01 -4.1249E+00  8.1133E+00  5.6745E+00 -5.6689E-02 -6.7472E-02  1.8837E-01 -1.0698E+00  3.7839E-01
             6.9267E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1686.88817928990        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.6707E-01  1.0196E+00  9.5986E-01  1.0342E+00  9.9283E-01  9.7981E-01  1.0077E+00  3.0748E-01  8.1508E-01  9.4611E-01
             1.0738E+00
 PARAMETER:  6.6520E-02  1.1941E-01  5.9033E-02  1.3366E-01  9.2801E-02  7.9602E-02  1.0772E-01 -1.0793E+00 -1.0447E-01  4.4606E-02
             1.7116E-01
 GRADIENT:   3.2158E+00  7.0874E-01 -3.8889E+00  7.6106E+00  5.3285E+00 -5.1639E-02 -5.9730E-02  1.8324E-01 -9.8989E-01  3.6919E-01
             6.5543E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1687.13411599752        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  9.7981E-01  1.0292E+00  9.6956E-01  1.0288E+00  9.9944E-01  9.8783E-01  9.8961E-01  3.0207E-01  8.2928E-01  9.5589E-01
             1.0734E+00
 PARAMETER:  7.9600E-02  1.2878E-01  6.9090E-02  1.2842E-01  9.9442E-02  8.7756E-02  8.9552E-02 -1.0971E+00 -8.7194E-02  5.4888E-02
             1.7084E-01
 GRADIENT:  -9.4808E-01 -1.1558E+00  3.0686E-01 -1.7401E+00 -9.5671E-01  3.1497E-01 -1.4241E-01  6.7985E-02 -8.4246E-02  4.9438E-01
            -1.8233E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1687.17430565812        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  9.8137E-01  1.1376E+00  9.1249E-01  9.5959E-01  1.0237E+00  9.8771E-01  9.1647E-01  2.0228E-01  8.7333E-01  9.5597E-01
             1.0737E+00
 PARAMETER:  8.1197E-02  2.2892E-01  8.4221E-03  5.8748E-02  1.2341E-01  8.7638E-02  1.2769E-02 -1.4981E+00 -3.5438E-02  5.4975E-02
             1.7112E-01
 GRADIENT:   8.4254E-01 -3.8999E-01  6.7050E-02 -6.1179E-01 -1.0563E-01 -6.6083E-02 -2.8290E-01  3.6862E-02 -4.0486E-01 -1.1041E-01
            -2.6392E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1687.18153523119        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      866
 NPARAMETR:  9.8124E-01  1.1935E+00  8.8732E-01  9.2378E-01  1.0399E+00  9.8820E-01  8.8355E-01  1.4240E-01  9.0250E-01  9.6194E-01
             1.0743E+00
 PARAMETER:  8.1066E-02  2.7690E-01 -1.9554E-02  2.0723E-02  1.3915E-01  8.8132E-02 -2.3813E-02 -1.8491E+00 -2.5881E-03  6.1198E-02
             1.7169E-01
 GRADIENT:  -1.5463E-01 -1.7972E-01 -1.0038E-02 -3.5560E-01 -4.0403E-02 -3.7256E-03  1.2949E-02  1.9738E-02  8.9889E-02  4.8837E-02
             1.1324E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1687.19166423815        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1045
 NPARAMETR:  9.8118E-01  1.1593E+00  8.9156E-01  9.4539E-01  1.0240E+00  9.8805E-01  9.0764E-01  4.9627E-02  8.8345E-01  9.5424E-01
             1.0743E+00
 PARAMETER:  8.0998E-02  2.4783E-01 -1.4788E-02  4.3841E-02  1.2367E-01  8.7973E-02  3.0908E-03 -2.9032E+00 -2.3918E-02  5.3156E-02
             1.7170E-01
 GRADIENT:  -9.4351E-02  4.3050E-02 -2.2936E-01  2.6288E-01  1.4373E-01 -2.9566E-02 -7.0042E-04  2.2536E-03 -6.5391E-03  4.2474E-02
             4.6491E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1687.19273299780        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1220
 NPARAMETR:  9.8122E-01  1.1612E+00  8.9160E-01  9.4415E-01  1.0250E+00  9.8811E-01  9.0612E-01  1.0341E-02  8.8466E-01  9.5491E-01
             1.0744E+00
 PARAMETER:  8.1041E-02  2.4945E-01 -1.4738E-02  4.2525E-02  1.2466E-01  8.8040E-02  1.4139E-03 -4.4716E+00 -2.2548E-02  5.3857E-02
             1.7173E-01
 GRADIENT:   3.7347E-04  4.0280E-02  5.9893E-03  6.3419E-02  1.7263E-02 -3.2867E-03 -1.0049E-02  8.5295E-05 -8.0364E-04 -1.6969E-02
            -1.4328E-02

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1687.19273884372        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1347
 NPARAMETR:  9.8122E-01  1.1610E+00  8.9154E-01  9.4423E-01  1.0248E+00  9.8812E-01  9.0637E-01  1.0000E-02  8.8453E-01  9.5488E-01
             1.0744E+00
 PARAMETER:  8.1040E-02  2.4928E-01 -1.4803E-02  4.2612E-02  1.2453E-01  8.8047E-02  1.6907E-03 -4.6383E+00 -2.2702E-02  5.3825E-02
             1.7174E-01
 GRADIENT:  -1.8198E-04 -1.0882E-02 -1.9343E-03 -1.0623E-02  1.3536E-03 -1.4349E-04  3.3880E-04  0.0000E+00  1.2292E-03  1.4504E-03
             8.4667E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1347
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.3841E-04 -1.1750E-02 -3.1349E-04  1.7805E-03 -2.4936E-02
 SE:             2.9789E-02  2.0033E-02  1.5022E-04  2.3568E-02  2.3218E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9361E-01  5.5751E-01  3.6891E-02  9.3978E-01  2.8282E-01

 ETASHRINKSD(%)  2.0402E-01  3.2888E+01  9.9497E+01  2.1045E+01  2.2218E+01
 ETASHRINKVR(%)  4.0761E-01  5.4959E+01  9.9997E+01  3.7662E+01  3.9499E+01
 EBVSHRINKSD(%)  4.9977E-01  3.2341E+01  9.9531E+01  2.1597E+01  2.0646E+01
 EBVSHRINKVR(%)  9.9704E-01  5.4223E+01  9.9998E+01  3.8530E+01  3.7030E+01
 RELATIVEINF(%)  9.8502E+01  1.2174E+00  2.0003E-04  1.9029E+00  6.3430E+00
 EPSSHRINKSD(%)  4.1662E+01
 EPSSHRINKVR(%)  6.5967E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1687.1927388437207     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -952.04191227998251     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.11
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1687.193       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.16E+00  8.92E-01  9.44E-01  1.02E+00  9.88E-01  9.06E-01  1.00E-02  8.85E-01  9.55E-01  1.07E+00
 


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
+        1.17E+03
 
 TH 2
+       -9.77E+00  4.77E+02
 
 TH 3
+        1.54E+01  1.60E+02  3.34E+02
 
 TH 4
+       -1.80E+01  4.78E+02 -1.77E+02  9.81E+02
 
 TH 5
+       -4.19E+00 -2.73E+02 -4.14E+02  1.91E+02  7.29E+02
 
 TH 6
+        1.80E+00 -2.48E+00  4.43E+00 -6.47E+00 -2.65E+00  2.00E+02
 
 TH 7
+        5.00E+00  1.76E+01  8.04E+00 -6.77E+00 -1.36E+01 -1.22E+00  5.71E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.67E+00 -1.97E+01 -2.04E+01  3.18E+01  1.00E+01  8.60E-01  3.20E+01  0.00E+00  1.04E+02
 
 TH10
+       -1.70E+00 -6.56E+00 -3.12E+01 -1.21E+01 -5.44E+01 -4.69E-01  7.97E+00  0.00E+00  4.96E+00  9.31E+01
 
 TH11
+       -8.82E+00 -1.99E+01 -3.84E+01  1.18E-02  7.23E+00  2.92E+00  5.86E+00  0.00E+00  1.36E+01  2.20E+01  1.98E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.795
Stop Time:
Sat Sep 25 09:39:38 CDT 2021
