Wed Sep 29 13:02:31 CDT 2021
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
$DATA ../../../../data/spa/A2/dat76.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1017.60567858861        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3653E+02  1.0235E+02  8.1174E+01  5.6547E+01  6.1513E+01 -6.4740E+00 -3.0496E+01 -2.0532E+01 -5.6689E+01 -9.6631E+01
            -9.8745E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1314.86108772949        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1246E+00  9.0015E-01  8.8994E-01  1.0988E+00  8.4606E-01  1.4551E+00  1.0588E+00  9.8977E-01  1.1893E+00  1.1106E+00
             3.0652E+00
 PARAMETER:  2.1745E-01 -5.1944E-03 -1.6606E-02  1.9422E-01 -6.7170E-02  4.7509E-01  1.5711E-01  8.9721E-02  2.7336E-01  2.0493E-01
             1.2201E+00
 GRADIENT:   2.3993E+02  1.6613E+01  1.1080E+01  2.0889E+01 -3.5098E+01  8.6817E+01  6.7160E+00  6.7082E+00  2.2379E+01  2.4124E+01
             1.2132E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1328.54674034628        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.1015E+00  5.0336E-01  3.7932E-01  1.3847E+00  4.0445E-01  1.3669E+00  1.1157E+00  5.0629E-01  1.1098E+00  6.5585E-01
             2.7219E+00
 PARAMETER:  1.9669E-01 -5.8644E-01 -8.6937E-01  4.2545E-01 -8.0523E-01  4.1255E-01  2.0946E-01 -5.8065E-01  2.0417E-01 -3.2182E-01
             1.1013E+00
 GRADIENT:   2.2435E+02  5.6734E+01 -2.5362E+01  2.7931E+02  2.5561E+01  7.1650E+01 -1.2418E+00  4.4802E+00  1.9914E-01  1.5746E+01
             8.8108E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1362.39816814403        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0124E+00  5.0950E-01  4.8164E-01  1.3295E+00  4.8217E-01  1.1756E+00  1.4273E+00  4.4306E-01  1.0956E+00  5.9524E-01
             2.4583E+00
 PARAMETER:  1.1235E-01 -5.7433E-01 -6.3056E-01  3.8482E-01 -6.2946E-01  2.6179E-01  4.5580E-01 -7.1406E-01  1.9132E-01 -4.1880E-01
             9.9947E-01
 GRADIENT:   5.5423E+01  2.6359E+01 -2.5189E+01  1.2772E+02  2.4964E+01  2.6710E+01  2.2159E+00  1.8464E-01  1.2815E+01 -2.1460E+00
             3.3595E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1380.05418893669        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      482
 NPARAMETR:  9.7266E-01  2.9360E-01  3.2764E-01  1.2382E+00  3.2774E-01  1.0791E+00  1.4942E+00  1.0999E-02  9.9577E-01  7.5314E-01
             2.0650E+00
 PARAMETER:  7.2278E-02 -1.1255E+00 -1.0158E+00  3.1370E-01 -1.0155E+00  1.7615E-01  5.0158E-01 -4.4099E+00  9.5762E-02 -1.8350E-01
             8.2512E-01
 GRADIENT:  -7.2555E-01  5.6997E+00 -8.1890E+00  2.7034E+01  5.6883E+00 -2.5175E+00 -2.3427E+00  4.5117E-04 -8.6018E+00  7.7751E+00
            -8.7139E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1381.46716562367        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  9.7213E-01  2.2230E-01  3.0015E-01  1.2202E+00  2.9714E-01  1.0828E+00  1.7354E+00  1.0000E-02  1.0391E+00  7.2405E-01
             2.0839E+00
 PARAMETER:  7.1729E-02 -1.4037E+00 -1.1035E+00  2.9903E-01 -1.1136E+00  1.7951E-01  6.5123E-01 -6.4786E+00  1.3838E-01 -2.2290E-01
             8.3423E-01
 GRADIENT:   8.8043E-01  1.2891E+00 -4.3942E-01 -1.1995E+00 -1.3668E+00 -3.2330E-01 -7.0431E-01  0.0000E+00  3.9586E+00 -1.7917E+00
             2.7037E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1382.78444524154        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      834
 NPARAMETR:  9.6177E-01  1.1032E-01  3.5231E-01  1.3097E+00  3.1786E-01  1.0751E+00  2.9656E+00  1.0000E-02  9.6481E-01  7.6780E-01
             2.1196E+00
 PARAMETER:  6.1015E-02 -2.1044E+00 -9.4324E-01  3.6978E-01 -1.0462E+00  1.7245E-01  1.1871E+00 -8.0010E+00  6.4172E-02 -1.6422E-01
             8.5124E-01
 GRADIENT:  -4.2551E+00  1.3894E+00  1.6925E+01  3.9464E+00 -2.4266E+01 -1.4086E-01 -1.3040E+00  0.0000E+00 -4.7667E+00  1.8240E+00
             2.2330E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1383.39054043898        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1010
 NPARAMETR:  9.5899E-01  4.8937E-02  3.4469E-01  1.3249E+00  3.0966E-01  1.0732E+00  5.1702E+00  1.0000E-02  9.7398E-01  7.6388E-01
             2.1036E+00
 PARAMETER:  5.8126E-02 -2.9172E+00 -9.6511E-01  3.8135E-01 -1.0723E+00  1.7067E-01  1.7429E+00 -1.1475E+01  7.3634E-02 -1.6935E-01
             8.4367E-01
 GRADIENT:  -1.2148E+00  1.2328E+00  1.9616E+00  3.0337E+00 -4.2428E+00 -1.6444E-01  1.8217E+00  0.0000E+00 -3.2952E-02 -2.0864E+00
            -1.5785E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1383.50027594385        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  9.5689E-01  1.1404E-02  3.5365E-01  1.3447E+00  3.1249E-01  1.0721E+00  1.0460E+01  1.0000E-02  9.6494E-01  7.7667E-01
             2.1055E+00
 PARAMETER:  5.5938E-02 -4.3738E+00 -9.3946E-01  3.9619E-01 -1.0632E+00  1.6966E-01  2.4475E+00 -1.7373E+01  6.4307E-02 -1.5274E-01
             8.4457E-01
 GRADIENT:  -1.4471E-01  6.5225E-02  1.4165E+00  1.8998E+00 -2.3843E+00  1.7887E-02  1.7245E-02  0.0000E+00 -4.6792E-01 -1.3815E-01
            -5.8633E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1383.50720013266        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1354
 NPARAMETR:  9.5689E-01  1.0000E-02  3.5346E-01  1.3435E+00  3.1255E-01  1.0720E+00  1.1068E+01  1.0000E-02  9.6600E-01  7.7722E-01
             2.1073E+00
 PARAMETER:  5.5938E-02 -4.5656E+00 -9.4000E-01  3.9529E-01 -1.0630E+00  1.6953E-01  2.5041E+00 -1.7901E+01  6.5405E-02 -1.5203E-01
             8.4543E-01
 GRADIENT:   1.0062E-01  0.0000E+00  1.0535E-01 -6.5053E-01  2.4518E-01  2.3098E-02 -2.0051E-02  0.0000E+00 -4.2060E-02  1.6717E-02
             1.6317E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1354
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9224E-04  1.7336E-03  1.3704E-05 -7.8446E-03 -4.0722E-03
 SE:             2.9354E-02  1.8636E-03  2.5372E-04  2.7866E-02  2.4510E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8934E-01  3.5227E-01  9.5693E-01  7.7832E-01  8.6804E-01

 ETASHRINKSD(%)  1.6612E+00  9.3757E+01  9.9150E+01  6.6462E+00  1.7887E+01
 ETASHRINKVR(%)  3.2947E+00  9.9610E+01  9.9993E+01  1.2851E+01  3.2575E+01
 EBVSHRINKSD(%)  1.5803E+00  9.4929E+01  9.9191E+01  5.8325E+00  1.7244E+01
 EBVSHRINKVR(%)  3.1355E+00  9.9743E+01  9.9993E+01  1.1325E+01  3.1514E+01
 RELATIVEINF(%)  8.4647E+01  3.4292E-02  3.0278E-04  1.7544E+01  2.7688E+00
 EPSSHRINKSD(%)  3.7498E+01
 EPSSHRINKVR(%)  6.0935E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1383.5072001326632     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -648.35637356892505     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.44
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1383.507       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.00E-02  3.53E-01  1.34E+00  3.13E-01  1.07E+00  1.11E+01  1.00E-02  9.66E-01  7.77E-01  2.11E+00
 


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
+        1.02E+03
 
 TH 2
+        0.00E+00  5.35E+02
 
 TH 3
+       -1.85E+01  0.00E+00  5.84E+03
 
 TH 4
+       -1.68E+01  0.00E+00 -4.41E+02  5.87E+02
 
 TH 5
+        8.92E+01  0.00E+00 -8.09E+03 -1.69E+02  1.30E+04
 
 TH 6
+        2.29E+00  0.00E+00  1.01E+01 -7.81E+00 -1.55E+00  1.60E+02
 
 TH 7
+        1.72E-02  0.00E+00 -1.25E-01 -7.87E-02  3.40E-01  4.76E-03  4.58E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.54E+00  0.00E+00  8.10E+01 -9.98E+00  1.54E+01 -7.19E-01  3.90E-02  0.00E+00  1.70E+02
 
 TH10
+       -3.39E+00  0.00E+00 -4.87E+01  7.60E+00 -7.09E+00  2.22E+00 -1.07E-01  0.00E+00  6.20E-01  1.50E+02
 
 TH11
+       -1.27E+01  0.00E+00 -3.39E+01 -6.68E+00  9.15E+00  2.82E+00 -4.96E-02  0.00E+00  6.63E+00  2.44E+01  5.65E+01
 
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
 #CPUT: Total CPU Time in Seconds,       23.869
Stop Time:
Wed Sep 29 13:02:56 CDT 2021
