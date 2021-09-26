Fri Sep 24 22:48:10 CDT 2021
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
$DATA ../../../../data/int/A3/dat88.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -376.676391929845        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3991E+01  2.6108E+02  2.5645E+02 -1.1798E+01  3.6307E+02  9.5074E+00 -3.7094E+02 -1.6531E+02 -4.5013E+01 -2.9320E+02
            -6.1331E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2520.96655688050        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0416E+00  6.0089E-01  7.5765E-01  1.4339E+00  5.6047E-01  8.3709E-01  1.2241E+00  8.6261E-01  8.4746E-01  9.6448E-01
             4.3261E+00
 PARAMETER:  1.4077E-01 -4.0935E-01 -1.7753E-01  4.6037E-01 -4.7899E-01 -7.7829E-02  3.0220E-01 -4.7788E-02 -6.5517E-02  6.3838E-02
             1.5647E+00
 GRADIENT:  -1.9885E+01  4.7597E+01  2.1391E+01  2.6529E+02 -6.9599E+01 -5.1573E+01  1.8614E+01  1.5536E+01 -2.3270E+01  4.2150E+01
             6.5620E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2553.14687593536        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8658E-01  3.2586E-01  1.9379E-01  1.7059E+00  2.1495E-01  9.5490E-01  1.0297E+00  3.6839E-01  2.0025E+00  7.1101E-01
             3.7631E+00
 PARAMETER:  8.6492E-02 -1.0213E+00 -1.5410E+00  6.3409E-01 -1.4373E+00  5.3851E-02  1.2923E-01 -8.9862E-01  7.9442E-01 -2.4107E-01
             1.4252E+00
 GRADIENT:  -9.7795E+01  2.0436E+02  1.8847E+01  2.4166E+02 -2.7447E+02 -1.6571E+01 -9.4924E+00  2.9461E+00 -1.8155E+01  5.2527E+00
             6.1086E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2807.05830414026        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0079E+00  2.9101E-01  2.3016E-01  1.1835E+00  2.5379E-01  9.4827E-01  1.1652E+00  2.8352E-01  1.3225E+00  6.9352E-01
             2.3589E+00
 PARAMETER:  1.0788E-01 -1.1344E+00 -1.3690E+00  2.6851E-01 -1.2712E+00  4.6885E-02  2.5288E-01 -1.1605E+00  3.7949E-01 -2.6597E-01
             9.5818E-01
 GRADIENT:   3.1203E+00 -1.3865E+01 -1.6948E+00  1.0107E+02  3.7020E+01 -1.5898E+01 -3.2324E+00 -2.1758E+00 -5.8681E+00  8.3678E+00
            -8.9501E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2816.01199080674        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0084E+00  2.4303E-01  1.5732E-01  1.0208E+00  1.9625E-01  1.0346E+00  1.2215E+00  1.4083E+00  1.3791E+00  5.5071E-01
             2.2550E+00
 PARAMETER:  1.0841E-01 -1.3146E+00 -1.7495E+00  1.2060E-01 -1.5284E+00  1.3402E-01  3.0012E-01  4.4237E-01  4.2141E-01 -4.9655E-01
             9.1315E-01
 GRADIENT:   2.7848E+00  5.5956E+01  2.9625E+01  3.0048E+01 -4.3580E+01  1.4886E+01 -5.0693E+00  5.5582E-01 -3.4279E+01  1.0963E+01
            -3.5875E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2819.32767813106        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0103E+00  2.1597E-01  1.3067E-01  9.5028E-01  1.7874E-01  1.0123E+00  1.2600E+00  1.5608E+00  1.5973E+00  4.8075E-01
             2.2693E+00
 PARAMETER:  1.1020E-01 -1.4326E+00 -1.9351E+00  4.9005E-02 -1.6218E+00  1.1225E-01  3.3110E-01  5.4518E-01  5.6831E-01 -6.3242E-01
             9.1947E-01
 GRADIENT:   8.9719E+00  2.4832E+01  9.7541E+00  4.6798E+00 -3.6514E+01  7.2411E+00 -2.0806E+00  4.7885E+00 -1.7192E+01 -4.8611E+00
             3.8568E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2823.06880092320        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  1.0097E+00  2.3043E-01  1.4757E-01  9.7398E-01  1.9491E-01  1.0024E+00  1.2962E+00  1.4951E+00  1.5836E+00  4.5273E-01
             2.2919E+00
 PARAMETER:  1.0967E-01 -1.3678E+00 -1.8134E+00  7.3634E-02 -1.5352E+00  1.0237E-01  3.5946E-01  5.0218E-01  5.5971E-01 -6.9245E-01
             9.2939E-01
 GRADIENT:   1.8324E-01  3.2403E+00 -2.4056E+00 -1.0945E+01  2.6803E+00  2.8184E+00  2.0401E+00  1.8716E+00 -1.1838E+00  5.0381E+00
             1.7143E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2823.10448503679        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  1.0096E+00  2.3026E-01  1.4749E-01  9.7425E-01  1.9464E-01  9.9439E-01  1.2960E+00  1.4953E+00  1.5831E+00  4.5212E-01
             2.2901E+00
 PARAMETER:  1.0960E-01 -1.3685E+00 -1.8140E+00  7.3911E-02 -1.5366E+00  9.4372E-02  3.5926E-01  5.0232E-01  5.5936E-01 -6.9382E-01
             9.2859E-01
 GRADIENT:   6.4507E+00  1.1458E+01  6.9700E+00 -9.1046E+00  5.5146E+01  3.7898E-01  2.2588E+00  2.2208E+00  7.4381E-02  5.2765E+00
             1.7158E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2823.16037701789        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  1.0096E+00  2.2913E-01  1.4746E-01  9.7432E-01  1.9404E-01  9.8718E-01  1.2957E+00  1.4953E+00  1.5831E+00  4.5153E-01
             2.2759E+00
 PARAMETER:  1.0957E-01 -1.3735E+00 -1.8142E+00  7.3981E-02 -1.5397E+00  8.7094E-02  3.5905E-01  5.0231E-01  5.5939E-01 -6.9511E-01
             9.2235E-01
 GRADIENT:   4.9051E-01  2.0573E+00  1.6285E+00 -9.4252E+00 -7.0758E+00 -2.9230E+00  1.5443E+00  1.0578E+00 -2.4577E+00  4.3691E+00
             3.0728E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2823.34682487275        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1014
 NPARAMETR:  1.0096E+00  2.2884E-01  1.4750E-01  9.7433E-01  1.9443E-01  9.9571E-01  1.2957E+00  1.4952E+00  1.5833E+00  4.5134E-01
             2.2714E+00
 PARAMETER:  1.0958E-01 -1.3747E+00 -1.8139E+00  7.3998E-02 -1.5377E+00  9.5701E-02  3.5901E-01  5.0224E-01  5.5948E-01 -6.9555E-01
             9.2040E-01
 GRADIENT:   7.2936E-01 -1.0324E+00  1.6400E-02 -9.9879E+00  2.2818E-01  4.7280E-01 -2.3553E+02  8.3488E-01 -2.2766E+00  5.3837E+00
            -6.0674E-01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -2823.34682487275        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1042
 NPARAMETR:  1.0096E+00  2.2884E-01  1.4750E-01  9.7433E-01  1.9443E-01  9.9571E-01  1.2957E+00  1.4952E+00  1.5833E+00  4.5134E-01
             2.2714E+00
 PARAMETER:  1.0958E-01 -1.3747E+00 -1.8139E+00  7.3998E-02 -1.5377E+00  9.5701E-02  3.5901E-01  5.0224E-01  5.5948E-01 -6.9555E-01
             9.2040E-01
 GRADIENT:   9.8704E-01 -9.6686E-01 -7.3845E-03 -9.9389E+00  4.4033E-01 -4.2290E+00  8.4153E+00  8.9730E-01 -2.1480E+00  4.5492E+00
            -2.0974E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1042
 NO. OF SIG. DIGITS IN FINAL EST.:  5.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.6301E-04  1.2449E-02  9.9682E-03  7.5928E-03  2.4871E-02
 SE:             2.9409E-02  2.5021E-02  2.4583E-02  2.7745E-02  1.5768E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7659E-01  6.1881E-01  6.8512E-01  7.8434E-01  1.1473E-01

 ETASHRINKSD(%)  1.4769E+00  1.6176E+01  1.7644E+01  7.0509E+00  4.7175E+01
 ETASHRINKVR(%)  2.9319E+00  2.9735E+01  3.2175E+01  1.3605E+01  7.2095E+01
 EBVSHRINKSD(%)  1.4888E+00  1.4643E+01  1.6064E+01  6.0250E+00  4.7722E+01
 EBVSHRINKVR(%)  2.9554E+00  2.7141E+01  2.9547E+01  1.1687E+01  7.2670E+01
 RELATIVEINF(%)  9.7007E+01  2.5294E+01  1.8108E+01  5.7479E+01  4.5355E+00
 EPSSHRINKSD(%)  2.1183E+01
 EPSSHRINKVR(%)  3.7879E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2823.3468248727531     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1169.2574651043424     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.45
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2823.347       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  2.29E-01  1.48E-01  9.74E-01  1.94E-01  9.96E-01  1.30E+00  1.50E+00  1.58E+00  4.51E-01  2.27E+00
 


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
+        3.31E+06
 
 TH 2
+       -2.13E+00  4.21E+05
 
 TH 3
+        1.37E+01  1.23E+01  5.85E+05
 
 TH 4
+       -1.04E+04 -5.18E+01 -1.66E+02  4.30E+06
 
 TH 5
+       -3.16E+01 -9.97E+03 -1.65E+04 -7.53E+02  5.00E+05
 
 TH 6
+        8.08E+08 -3.45E+01  1.12E+01 -6.36E+04 -9.06E+00  3.62E+06
 
 TH 7
+       -9.36E+05  6.32E+01 -3.99E+00  7.34E+04 -6.94E+01  1.04E+04  1.89E+05
 
 TH 8
+        4.26E-01  3.57E+01  3.78E+01  2.01E-01  2.43E+01  3.69E+00 -2.46E+00  7.24E+04
 
 TH 9
+        2.95E+01 -7.97E+00  1.13E+02 -3.47E+00  2.93E+02 -1.24E+01  2.15E+00  3.25E+01  5.20E+04
 
 TH10
+       -3.12E+01  5.27E+01  2.14E+01 -1.51E+01  5.69E+02  1.50E+01  2.50E+01  3.54E+03 -3.51E+03  4.14E+05
 
 TH11
+       -1.38E+01 -4.15E+01 -8.51E+01  4.20E+00  4.83E+01  4.80E+00  7.71E+00  1.69E+02  8.38E+01 -2.81E+01  9.52E+03
 
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
 #CPUT: Total CPU Time in Seconds,       42.792
Stop Time:
Fri Sep 24 22:48:54 CDT 2021
