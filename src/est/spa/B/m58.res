Sat Sep 25 07:24:18 CDT 2021
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
$DATA ../../../../data/spa/B/dat58.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m58.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1628.75141285650        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6072E+02 -1.4229E+01  3.0017E+01 -4.0757E+01 -1.7427E+01  1.3265E+01 -1.7298E+00 -4.7902E+00  1.2320E+01 -8.9951E+00
             1.7749E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1635.23711088879        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.4492E-01  1.0384E+00  8.4362E-01  9.9773E-01  9.5679E-01  9.8963E-01  1.1028E+00  9.9452E-01  8.5768E-01  1.0544E+00
             9.6592E-01
 PARAMETER:  4.3340E-02  1.3769E-01 -7.0055E-02  9.7726E-02  5.5824E-02  8.9576E-02  1.9788E-01  9.4503E-02 -5.3530E-02  1.5300E-01
             6.5327E-02
 GRADIENT:   3.4368E+01  2.7182E+00 -9.2327E+00  1.0672E+01  2.3637E+00  1.6404E+01 -3.6932E+00  6.8057E+00 -4.0639E+00  1.1813E+01
             9.5763E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1636.54494262483        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.4344E-01  9.4204E-01  6.7991E-01  1.0486E+00  8.0275E-01  9.8940E-01  1.2716E+00  7.3844E-01  7.6030E-01  8.4499E-01
             9.5343E-01
 PARAMETER:  4.1776E-02  4.0288E-02 -2.8580E-01  1.4748E-01 -1.1971E-01  8.9340E-02  3.4031E-01 -2.0322E-01 -1.7404E-01 -6.8428E-02
             5.2307E-02
 GRADIENT:   2.7719E+01  1.2588E+01 -2.2035E+01  4.1207E+01  2.0503E+01  1.5633E+01  8.8963E-01  5.7594E+00 -1.0312E+01  5.5796E+00
             3.6427E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1638.06322621083        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      257
 NPARAMETR:  9.3091E-01  8.8990E-01  5.5058E-01  1.0433E+00  6.7905E-01  9.4435E-01  1.2649E+00  4.1056E-01  8.0618E-01  6.9676E-01
             9.4010E-01
 PARAMETER:  2.8409E-02 -1.6649E-02 -4.9679E-01  1.4241E-01 -2.8706E-01  4.2743E-02  3.3501E-01 -7.9023E-01 -1.1545E-01 -2.6132E-01
             3.8231E-02
 GRADIENT:  -4.8296E+01  6.9170E+00  1.2403E-01 -9.4364E-01 -1.2277E+01 -9.2758E+00 -2.2730E+00  1.7594E+00 -9.2407E-01  2.1104E+00
            -1.5217E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1639.97116525040        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      433
 NPARAMETR:  9.4873E-01  7.1283E-01  6.3691E-01  1.1543E+00  6.7718E-01  9.6700E-01  1.5786E+00  2.2271E-01  7.5608E-01  7.8548E-01
             9.4652E-01
 PARAMETER:  4.7369E-02 -2.3851E-01 -3.5112E-01  2.4348E-01 -2.8982E-01  6.6446E-02  5.5653E-01 -1.4019E+00 -1.7960E-01 -1.4146E-01
             4.5033E-02
 GRADIENT:   2.6198E+00  1.2610E+00  6.2707E-01 -1.6120E-01 -9.8919E-01  2.0391E+00  2.6996E-01  7.8113E-02  5.7632E-01  1.0193E-01
            -4.9593E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1640.01409478655        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      608
 NPARAMETR:  9.4740E-01  6.7234E-01  6.3638E-01  1.1753E+00  6.6334E-01  9.6146E-01  1.6542E+00  1.7869E-01  7.4161E-01  7.8735E-01
             9.4880E-01
 PARAMETER:  4.5961E-02 -2.9699E-01 -3.5196E-01  2.6155E-01 -3.1047E-01  6.0694E-02  6.0331E-01 -1.6221E+00 -1.9893E-01 -1.3908E-01
             4.7445E-02
 GRADIENT:   1.1582E-02 -3.2439E-02 -4.1440E-02 -1.5402E-02  1.9134E-02  6.2527E-03  6.5657E-03  6.2648E-03  1.3305E-02  2.2487E-02
             6.9657E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1640.01430263306        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  9.4739E-01  6.7250E-01  6.3608E-01  1.1752E+00  6.6338E-01  9.6144E-01  1.6537E+00  1.6586E-01  7.4174E-01  7.8883E-01
             9.4886E-01
 PARAMETER:  4.5956E-02 -2.9675E-01 -3.5243E-01  2.6146E-01 -3.1040E-01  6.0681E-02  6.0302E-01 -1.6966E+00 -1.9876E-01 -1.3720E-01
             4.7504E-02
 GRADIENT:  -1.7681E-02 -8.8570E-02 -1.0672E-01  1.4518E-01  2.2934E-01 -3.1966E-03 -1.3122E-03  1.4936E-03 -1.7206E-03  2.5347E-02
             9.3946E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1640.01537904288        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  9.4751E-01  6.7883E-01  6.3121E-01  1.1712E+00  6.6250E-01  9.6161E-01  1.6405E+00  1.0965E-01  7.4424E-01  7.9029E-01
             9.4884E-01
 PARAMETER:  4.6083E-02 -2.8738E-01 -3.6011E-01  2.5799E-01 -3.1173E-01  6.0855E-02  5.9502E-01 -2.1104E+00 -1.9539E-01 -1.3535E-01
             4.7487E-02
 GRADIENT:  -9.6750E-03  1.0276E-01 -1.3983E-01  5.2625E-01  9.8377E-02  7.6846E-04  8.9910E-02  1.1165E-03  2.8562E-02  1.7763E-01
             4.1356E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1640.01705856085        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1145
 NPARAMETR:  9.4766E-01  6.8720E-01  6.2537E-01  1.1651E+00  6.6150E-01  9.6180E-01  1.6208E+00  7.0905E-02  7.4733E-01  7.8561E-01
             9.4863E-01
 PARAMETER:  4.6236E-02 -2.7513E-01 -3.6942E-01  2.5282E-01 -3.1324E-01  6.1054E-02  5.8290E-01 -2.5464E+00 -1.9125E-01 -1.4129E-01
             4.7268E-02
 GRADIENT:   7.8386E-03 -1.7850E-01 -1.4187E-02 -4.2493E-01 -8.5595E-03  2.5915E-03 -4.2150E-02  1.6296E-03 -2.1925E-04  2.2866E-02
             1.3524E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1640.01762584170        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1320
 NPARAMETR:  9.4770E-01  6.9182E-01  6.2462E-01  1.1627E+00  6.6271E-01  9.6187E-01  1.6122E+00  3.0000E-02  7.4890E-01  7.8668E-01
             9.4856E-01
 PARAMETER:  4.6279E-02 -2.6843E-01 -3.7062E-01  2.5075E-01 -3.1142E-01  6.1122E-02  5.7758E-01 -3.4066E+00 -1.8915E-01 -1.3994E-01
             4.7192E-02
 GRADIENT:  -2.7450E-03  4.7226E-03  5.4371E-03 -4.2159E-02  6.4351E-02 -2.9352E-04 -5.1443E-03  2.8326E-04  2.9343E-03  6.9884E-03
             2.9247E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1640.01781787879        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1499
 NPARAMETR:  9.4769E-01  6.9018E-01  6.2389E-01  1.1634E+00  6.6161E-01  9.6186E-01  1.6151E+00  1.0000E-02  7.4841E-01  7.8595E-01
             9.4859E-01
 PARAMETER:  4.6276E-02 -2.7081E-01 -3.7178E-01  2.5138E-01 -3.1308E-01  6.1111E-02  5.7943E-01 -4.6147E+00 -1.8981E-01 -1.4086E-01
             4.7224E-02
 GRADIENT:  -1.7859E-03 -1.4986E-02 -2.1298E-02 -2.0956E-02  9.4382E-03 -2.0636E-04 -8.9267E-04  0.0000E+00  1.3852E-03  7.5843E-03
             3.8007E-03

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1640.01781928609        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1591
 NPARAMETR:  9.4769E-01  6.9028E-01  6.2398E-01  1.1634E+00  6.6171E-01  9.6186E-01  1.6150E+00  1.0000E-02  7.4844E-01  7.8604E-01
             9.4859E-01
 PARAMETER:  4.6276E-02 -2.7065E-01 -3.7163E-01  2.5136E-01 -3.1293E-01  6.1111E-02  5.7933E-01 -4.5821E+00 -1.8977E-01 -1.4074E-01
             4.7220E-02
 GRADIENT:   2.6964E-04 -4.4723E-04 -7.8861E-04 -7.8910E-04  6.6817E-04 -1.4261E-04 -1.0372E-05  0.0000E+00  7.5720E-05  3.0492E-04
             1.8009E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1591
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6965E-04  1.5462E-02 -5.9098E-04 -1.4333E-02  2.8273E-03
 SE:             2.9862E-02  2.1490E-02  2.4633E-04  2.5024E-02  2.2909E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9280E-01  4.7182E-01  1.6435E-02  5.6678E-01  9.0178E-01

 ETASHRINKSD(%)  1.0000E-10  2.8006E+01  9.9175E+01  1.6167E+01  2.3253E+01
 ETASHRINKVR(%)  1.0000E-10  4.8169E+01  9.9993E+01  2.9720E+01  4.1099E+01
 EBVSHRINKSD(%)  4.0905E-01  2.7805E+01  9.9262E+01  1.6273E+01  2.2193E+01
 EBVSHRINKVR(%)  8.1643E-01  4.7879E+01  9.9995E+01  2.9899E+01  3.9460E+01
 RELATIVEINF(%)  9.8859E+01  6.1313E+00  5.2831E-04  1.0615E+01  4.0831E+00
 EPSSHRINKSD(%)  4.4595E+01
 EPSSHRINKVR(%)  6.9303E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1640.0178192860924     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -904.86699272235421     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1640.018       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.48E-01  6.90E-01  6.24E-01  1.16E+00  6.62E-01  9.62E-01  1.61E+00  1.00E-02  7.48E-01  7.86E-01  9.49E-01
 


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
+        1.33E+03
 
 TH 2
+       -1.11E+01  4.91E+02
 
 TH 3
+        2.69E+01  3.68E+02  1.71E+03
 
 TH 4
+       -1.02E+01  3.38E+02 -5.37E+02  1.04E+03
 
 TH 5
+       -5.15E+00 -5.99E+02 -1.87E+03  5.12E+02  2.49E+03
 
 TH 6
+       -1.85E+00 -2.75E+00  5.24E+00 -2.48E+00 -2.32E+00  2.14E+02
 
 TH 7
+        1.82E+00  4.08E+01 -2.05E+01 -1.15E+01  5.67E+00  7.91E-02  2.77E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.39E+00 -2.93E+01 -3.84E+01  2.17E+01  1.33E+01  1.02E+00  1.22E+01  0.00E+00  1.82E+02
 
 TH10
+       -2.51E+00 -1.06E+01 -1.40E+02 -3.42E+01 -2.57E+01 -1.92E+00  8.03E+00  0.00E+00  1.74E+01  1.29E+02
 
 TH11
+       -9.61E+00 -8.50E+00 -4.91E+01 -1.04E+01  1.98E+01  4.08E+00  3.69E+00  0.00E+00  1.73E+01  2.39E+01  2.33E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.203
Stop Time:
Sat Sep 25 07:24:44 CDT 2021
