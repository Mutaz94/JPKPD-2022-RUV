Sat Sep 25 07:50:45 CDT 2021
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
$DATA ../../../../data/spa/A1/dat6.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1282.57869529627        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3100E+01 -5.4981E+01 -5.2872E+01 -5.5156E+00  2.1704E+02  1.2347E+01 -2.0354E+01  3.3332E+00 -2.6859E+01 -6.7897E+01
            -6.6281E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1479.36288429287        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0570E+00  9.7483E-01  1.0496E+00  1.0286E+00  8.8601E-01  9.5292E-01  1.0009E+00  9.5486E-01  1.0681E+00  9.9206E-01
             2.1428E+00
 PARAMETER:  1.5543E-01  7.4505E-02  1.4837E-01  1.2824E-01 -2.1032E-02  5.1776E-02  1.0094E-01  5.3810E-02  1.6586E-01  9.2030E-02
             8.6211E-01
 GRADIENT:   1.2116E+02 -2.2785E+01 -1.0805E+01 -2.7560E+01  2.0911E+01 -6.6943E+00  5.1145E+00  6.0213E+00  4.3602E+00  9.9971E+00
             3.0739E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1484.54899225635        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0511E+00  7.3313E-01  8.3254E-01  1.2000E+00  7.0311E-01  9.4196E-01  8.9031E-01  3.3367E-01  1.1026E+00  8.7131E-01
             2.0309E+00
 PARAMETER:  1.4979E-01 -2.1044E-01 -8.3271E-02  2.8234E-01 -2.5224E-01  4.0210E-02 -1.6185E-02 -9.9761E-01  1.9763E-01 -3.7761E-02
             8.0849E-01
 GRADIENT:   1.0921E+02  3.6743E-01 -2.0038E+01  3.4981E+01  2.9682E+01 -1.2316E+01 -5.0371E-01  1.0398E+00  2.4056E+01  2.4863E+00
             1.1172E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1491.28977968029        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.9618E-01  5.8745E-01  5.9912E-01  1.2265E+00  5.2713E-01  9.7606E-01  1.5324E+00  1.0583E-01  8.6205E-01  7.2033E-01
             1.9047E+00
 PARAMETER:  9.6176E-02 -4.3197E-01 -4.1230E-01  3.0415E-01 -5.4031E-01  7.5769E-02  5.2684E-01 -2.1459E+00 -4.8439E-02 -2.2805E-01
             7.4431E-01
 GRADIENT:  -2.7825E+01  1.4741E+01  7.5175E-01  2.0988E+01 -1.6381E+00  1.2289E+00  1.5185E+00  1.6906E-01 -7.4526E+00  1.7957E+00
            -2.4373E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1492.51539720637        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  1.0163E+00  5.1687E-01  5.9197E-01  1.2647E+00  5.0184E-01  9.7599E-01  1.5494E+00  6.9982E-02  8.7987E-01  7.1653E-01
             1.9138E+00
 PARAMETER:  1.1616E-01 -5.5997E-01 -4.2431E-01  3.3485E-01 -5.8948E-01  7.5695E-02  5.3785E-01 -2.5595E+00 -2.7976E-02 -2.3334E-01
             7.4911E-01
 GRADIENT:   6.6615E+00  1.1973E+01  6.6275E+00  1.5652E+01 -1.4973E+01  5.9027E-01 -1.2949E+00  6.6907E-02 -1.3219E+00 -5.3115E-01
            -8.6117E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1495.48971145656        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      542
 NPARAMETR:  1.0054E+00  2.4176E-01  6.3526E-01  1.4048E+00  4.6906E-01  9.7385E-01  2.3995E+00  1.0000E-02  8.4249E-01  7.6881E-01
             1.9019E+00
 PARAMETER:  1.0543E-01 -1.3198E+00 -3.5373E-01  4.3991E-01 -6.5703E-01  7.3503E-02  9.7528E-01 -6.3182E+00 -7.1388E-02 -1.6291E-01
             7.4284E-01
 GRADIENT:  -3.1968E+00  3.3229E+00  7.2151E+00  1.1442E+00 -1.0963E+01  1.2165E+00  6.3328E-01  0.0000E+00 -3.3697E-01 -2.2583E+00
            -4.0428E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1497.25074613624        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      717
 NPARAMETR:  1.0002E+00  6.2063E-02  6.4454E-01  1.5013E+00  4.4746E-01  9.6916E-01  3.1866E+00  1.0000E-02  8.1853E-01  8.0418E-01
             1.8856E+00
 PARAMETER:  1.0024E-01 -2.6796E+00 -3.3923E-01  5.0636E-01 -7.0416E-01  6.8675E-02  1.2590E+00 -1.4298E+01 -1.0025E-01 -1.1793E-01
             7.3423E-01
 GRADIENT:  -2.0690E+00  9.3126E-01  5.1487E-01  1.8386E+01 -5.3386E+00 -1.0153E-01 -4.9256E-01  0.0000E+00 -2.6104E+00 -1.3034E-01
            -2.2085E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1497.83261385035        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      892
 NPARAMETR:  9.9963E-01  1.0000E-02  6.4570E-01  1.5177E+00  4.4111E-01  9.6728E-01  4.4358E+00  1.0000E-02  8.1571E-01  8.0431E-01
             1.8918E+00
 PARAMETER:  9.9629E-02 -4.5189E+00 -3.3742E-01  5.1719E-01 -7.1846E-01  6.6735E-02  1.5897E+00 -2.5863E+01 -1.0370E-01 -1.1777E-01
             7.3755E-01
 GRADIENT:   1.1753E+00  0.0000E+00  8.8320E-01  1.7296E+00 -2.0347E+00 -3.7107E-01 -2.9371E-02  0.0000E+00 -3.3559E-02 -2.8438E-01
             4.5464E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1497.87459376467        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  9.9843E-01  1.0000E-02  6.5039E-01  1.5167E+00  4.4392E-01  9.6807E-01  9.4370E+00  1.0000E-02  8.1501E-01  8.0806E-01
             1.8898E+00
 PARAMETER:  9.8433E-02 -4.5156E+00 -3.3018E-01  5.1652E-01 -7.1212E-01  6.7551E-02  2.3446E+00 -2.5848E+01 -1.0456E-01 -1.1312E-01
             7.3645E-01
 GRADIENT:   9.4794E+00  0.0000E+00  5.3843E-01  2.1508E+01  5.5700E+00  1.0723E+00 -6.6604E-02  0.0000E+00  6.2776E-01  1.7695E-01
             8.0489E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1497.88439774300        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1188
 NPARAMETR:  9.9809E-01  1.0000E-02  6.5000E-01  1.5159E+00  4.4357E-01  9.6770E-01  1.1115E+01  1.0000E-02  8.1457E-01  8.0766E-01
             1.8889E+00
 PARAMETER:  9.8085E-02 -4.5156E+00 -3.3078E-01  5.1600E-01 -7.1290E-01  6.7162E-02  2.5083E+00 -2.5848E+01 -1.0510E-01 -1.1361E-01
             7.3599E-01
 GRADIENT:  -2.0203E+00  0.0000E+00  5.5405E-01 -4.8777E+00  2.6739E-01 -1.3894E-01  4.6312E-03  0.0000E+00  2.8144E-02  7.1372E-02
             1.2661E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1497.88928072669        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  9.9902E-01  1.0000E-02  6.5080E-01  1.5187E+00  4.4400E-01  9.6811E-01  1.1098E+01  1.0000E-02  8.1443E-01  8.0756E-01
             1.8895E+00
 PARAMETER:  9.9017E-02 -4.5156E+00 -3.2955E-01  5.1788E-01 -7.1194E-01  6.7586E-02  2.5068E+00 -2.5848E+01 -1.0527E-01 -1.1374E-01
             7.3634E-01
 GRADIENT:   1.9480E-03  0.0000E+00 -9.2004E-02  7.4695E-02  2.4446E-01 -3.7798E-03 -9.9046E-03  0.0000E+00 -3.3216E-03  1.1153E-02
             1.0191E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1315
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9740E-05  9.3497E-04 -3.5935E-05 -7.8357E-03 -1.1300E-02
 SE:             2.9444E-02  1.9718E-03  2.0870E-04  2.8003E-02  2.4008E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9919E-01  6.3538E-01  8.6329E-01  7.7962E-01  6.3787E-01

 ETASHRINKSD(%)  1.3592E+00  9.3394E+01  9.9301E+01  6.1853E+00  1.9571E+01
 ETASHRINKVR(%)  2.7000E+00  9.9564E+01  9.9995E+01  1.1988E+01  3.5311E+01
 EBVSHRINKSD(%)  1.4732E+00  9.4227E+01  9.9272E+01  5.8783E+00  1.8601E+01
 EBVSHRINKVR(%)  2.9247E+00  9.9667E+01  9.9995E+01  1.1411E+01  3.3742E+01
 RELATIVEINF(%)  8.7948E+01  1.5883E-02  2.8299E-04  6.7343E+00  2.8936E+00
 EPSSHRINKSD(%)  3.7467E+01
 EPSSHRINKVR(%)  6.0896E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1497.8892807266907     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -762.73845416295251     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1497.889       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  1.00E-02  6.51E-01  1.52E+00  4.44E-01  9.68E-01  1.11E+01  1.00E-02  8.14E-01  8.08E-01  1.89E+00
 


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
+        1.13E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        5.21E+00  0.00E+00  1.39E+03
 
 TH 4
+       -2.50E+01  0.00E+00 -1.64E+02  6.45E+02
 
 TH 5
+        3.72E+01  0.00E+00 -2.67E+03 -1.31E+02  5.73E+03
 
 TH 6
+       -1.94E+00  0.00E+00  4.10E+00 -6.65E+00 -5.88E+00  1.81E+02
 
 TH 7
+        2.61E-02  0.00E+00 -6.77E-02 -9.48E-02  1.90E-01  3.01E-02  2.12E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.06E+00  0.00E+00  2.70E+01 -9.09E+00 -4.88E+00  1.05E+01 -4.26E-03  0.00E+00  2.14E+02
 
 TH10
+        2.57E+00  0.00E+00 -1.60E+01  2.24E+00 -6.25E+01  5.18E+00 -6.75E-02  0.00E+00 -7.91E+00  1.18E+02
 
 TH11
+       -1.23E+01  0.00E+00 -2.07E+01 -7.56E+00  1.03E+01  2.89E+00 -2.15E-02  0.00E+00  8.66E+00  2.64E+01  7.04E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.632
Stop Time:
Sat Sep 25 07:51:07 CDT 2021
