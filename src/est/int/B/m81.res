Fri Sep 24 19:48:26 CDT 2021
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
$DATA ../../../../data/int/B/dat81.csv ignore=@
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3382.02968683753        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.1699E+00 -1.6438E+01  8.2607E+01 -4.1967E+01  8.8085E+01 -3.0649E+01 -4.4057E+01 -2.2397E+02 -1.0380E+01 -3.9373E+01
            -7.1597E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3889.29162868449        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8990E-01  9.9676E-01  9.3599E-01  1.0242E+00  9.3257E-01  1.1071E+00  8.8820E-01  8.8783E-01  8.5477E-01  1.0688E+00
             1.2596E+00
 PARAMETER:  8.9853E-02  9.6753E-02  3.3850E-02  1.2395E-01  3.0192E-02  2.0173E-01 -1.8555E-02 -1.8979E-02 -5.6928E-02  1.6653E-01
             3.3078E-01
 GRADIENT:  -3.5294E+01 -7.0889E+00 -8.9543E+00  1.0386E+00 -1.3689E+01  1.3982E+01 -3.0807E+00  1.0926E+01 -1.6808E+01  5.1969E+00
             4.2877E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3898.59855512236        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.9140E-01  1.0175E+00  9.3650E-01  1.0308E+00  9.6255E-01  1.1067E+00  8.2529E-01  5.1461E-01  1.0031E+00  1.0089E+00
             1.2267E+00
 PARAMETER:  9.1365E-02  1.1734E-01  3.4398E-02  1.3037E-01  6.1831E-02  2.0140E-01 -9.2022E-02 -5.6435E-01  1.0309E-01  1.0882E-01
             3.0432E-01
 GRADIENT:  -3.0547E+01  5.8828E+00  6.7801E-01  4.3945E+01  2.7700E+01  1.4652E+01 -8.9197E+00 -5.7855E+00  2.2851E+01 -2.1202E+01
             3.7033E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3944.62866293376        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9951E-01  9.1172E-01  8.7081E-01  1.0690E+00  8.5600E-01  1.0721E+00  1.0286E+00  6.9367E-01  8.8333E-01  9.1665E-01
             9.9758E-01
 PARAMETER:  9.9511E-02  7.5803E-03 -3.8329E-02  1.6675E-01 -5.5483E-02  1.6962E-01  1.2818E-01 -2.6575E-01 -2.4055E-02  1.2968E-02
             9.7580E-02
 GRADIENT:   2.4006E-01 -3.6573E+00  6.4879E+00  4.3030E+01 -1.2728E+01  5.1041E+00  3.8032E+00 -6.8999E+00 -1.0564E+01 -1.1847E+01
             3.7464E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3944.68090792748        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:      427             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9969E-01  9.1193E-01  8.7093E-01  1.0688E+00  8.5620E-01  1.1273E+00  1.0283E+00  6.9381E-01  8.8350E-01  9.1709E-01
             9.9760E-01
 PARAMETER:  9.9686E-02  7.8041E-03 -3.8192E-02  1.6656E-01 -5.5246E-02  2.1985E-01  1.2789E-01 -2.6556E-01 -2.3859E-02  1.3451E-02
             9.7594E-02
 GRADIENT:   5.1794E+00 -3.7769E+00  6.4850E+00  4.2573E+01 -1.2738E+01  2.9743E+01  3.7752E+00 -6.9124E+00 -1.0508E+01 -1.1770E+01
             3.7570E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3944.81786869498        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      590            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9982E-01  9.1197E-01  8.7091E-01  1.0688E+00  8.5625E-01  1.1004E+00  1.0282E+00  6.9429E-01  8.8366E-01  9.1738E-01
             9.9750E-01
 PARAMETER:  9.9818E-02  7.8489E-03 -3.8217E-02  1.6651E-01 -5.5191E-02  1.9566E-01  1.2779E-01 -2.6487E-01 -2.3684E-02  1.3766E-02
             9.7495E-02
 GRADIENT:   3.2743E+00 -3.7683E+00  6.3535E+00  4.2559E+01 -1.2664E+01  1.8138E+01  3.7957E+00 -6.8681E+00 -1.0455E+01 -1.1729E+01
             3.7393E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3946.22301441880        NO. OF FUNC. EVALS.: 153
 CUMULATIVE NO. OF FUNC. EVALS.:      743
 NPARAMETR:  1.0142E+00  9.2289E-01  8.6940E-01  1.0687E+00  8.6692E-01  1.1004E+00  9.8009E-01  6.9451E-01  8.8372E-01  9.9468E-01
             9.9747E-01
 PARAMETER:  1.1408E-01  1.9757E-02 -3.9947E-02  1.6641E-01 -4.2804E-02  1.9565E-01  7.9892E-02 -2.6455E-01 -2.3615E-02  9.4662E-02
             9.7467E-02
 GRADIENT:  -1.8919E+01 -1.0992E+01 -7.4966E+00  3.2790E+01 -8.0747E+00  9.6667E-01  2.0152E+00 -6.3569E+00 -1.2142E+01 -8.0064E-01
             3.7426E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3946.74693745592        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      920
 NPARAMETR:  1.0248E+00  9.4833E-01  8.9323E-01  1.0687E+00  8.9213E-01  1.1004E+00  9.4728E-01  6.9451E-01  8.8372E-01  1.0240E+00
             9.9747E-01
 PARAMETER:  1.2452E-01  4.6947E-02 -1.2909E-02  1.6641E-01 -1.4140E-02  1.9565E-01  4.5840E-02 -2.6455E-01 -2.3615E-02  1.2374E-01
             9.7467E-02
 GRADIENT:  -5.1668E-03 -3.2466E-01  4.6217E-01  5.6889E+01 -5.4004E-01  1.1641E+00 -4.0925E-03 -8.2022E+00 -1.2969E+01  9.3206E-02
             3.3617E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3947.24391772370        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  1.0243E+00  9.4873E-01  8.9324E-01  1.0650E+00  8.9264E-01  1.0991E+00  9.4713E-01  7.1037E-01  8.8841E-01  1.0237E+00
             9.9514E-01
 PARAMETER:  1.2397E-01  4.7367E-02 -1.2898E-02  1.6296E-01 -1.3576E-02  1.9445E-01  4.5678E-02 -2.4198E-01 -1.8322E-02  1.2342E-01
             9.5129E-02
 GRADIENT:   5.9368E+01  5.0889E+00  2.2047E-01  6.8069E+01  4.2371E+00  1.7890E+01  6.2556E-01 -7.1029E+00 -1.0315E+01  6.9773E-01
             3.0663E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3947.27726005800        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  1.0241E+00  9.4876E-01  8.9326E-01  1.0645E+00  8.9266E-01  1.0987E+00  9.4709E-01  7.1073E-01  8.8897E-01  1.0237E+00
             9.9508E-01
 PARAMETER:  1.2379E-01  4.7406E-02 -1.2873E-02  1.6254E-01 -1.3551E-02  1.9416E-01  4.5644E-02 -2.4146E-01 -1.7696E-02  1.2340E-01
             9.5064E-02
 GRADIENT:  -1.2857E+00 -2.6676E+00 -6.2404E-01  4.7653E+01 -1.5951E+00  5.6706E-01  2.6313E-02 -7.1755E+00 -1.0957E+01  2.8170E-02
             3.0416E+01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -3947.27726005800        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1258
 NPARAMETR:  1.0241E+00  9.4876E-01  8.9326E-01  1.0645E+00  8.9266E-01  1.0987E+00  9.4709E-01  7.1073E-01  8.8897E-01  1.0237E+00
             9.9508E-01
 PARAMETER:  1.2379E-01  4.7406E-02 -1.2873E-02  1.6254E-01 -1.3551E-02  1.9416E-01  4.5644E-02 -2.4146E-01 -1.7696E-02  1.2340E-01
             9.5064E-02
 GRADIENT:  -1.3037E+06 -1.6139E+06 -1.6138E+06  9.9296E+05 -1.6139E+06 -1.6624E+06  1.6138E+06  1.3366E+06  3.2277E+06  1.3078E+06
            -3.2281E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1258
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.8949E-04 -2.4079E-02 -1.7283E-02 -7.2418E-03 -2.0816E-02
 SE:             2.9874E-02  2.2481E-02  1.6634E-02  2.8630E-02  2.5560E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7625E-01  2.8414E-01  2.9880E-01  8.0031E-01  4.1543E-01

 ETASHRINKSD(%)  1.0000E-10  2.4684E+01  4.4274E+01  4.0856E+00  1.4369E+01
 ETASHRINKVR(%)  1.0000E-10  4.3276E+01  6.8947E+01  8.0042E+00  2.6673E+01
 EBVSHRINKSD(%)  2.0963E-01  2.4863E+01  4.8489E+01  7.4597E+00  1.4319E+01
 EBVSHRINKVR(%)  4.1881E-01  4.3544E+01  7.3466E+01  1.4363E+01  2.6588E+01
 RELATIVEINF(%)  9.9578E+01  2.6717E+01  1.5337E+01  6.0231E+01  2.9410E+01
 EPSSHRINKSD(%)  2.2240E+01
 EPSSHRINKVR(%)  3.9533E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3947.2772600579983     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2293.1879002895876     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.48
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3947.277       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  9.49E-01  8.93E-01  1.06E+00  8.93E-01  1.10E+00  9.47E-01  7.11E-01  8.89E-01  1.02E+00  9.95E-01
 


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
+        5.02E+09
 
 TH 2
+       -5.95E+04  8.96E+09
 
 TH 3
+       -2.51E+05 -8.46E+04  1.01E+10
 
 TH 4
+        5.17E+03  4.39E+04  1.84E+05  2.70E+09
 
 TH 5
+       -6.94E+04  2.83E+09 -9.89E+04  5.10E+04  1.01E+10
 
 TH 6
+       -8.85E+08 -1.50E+04  1.26E+09  8.22E+03 -1.59E+04  1.77E+09
 
 TH 7
+       -6.72E+09  7.97E+04  3.35E+05  1.46E+09  9.29E+04  1.50E+04  9.00E+09
 
 TH 8
+        5.21E+03  4.40E+04  1.85E+05 -7.64E+04  5.13E+04  8.29E+03 -7.24E+04  2.74E+09
 
 TH 9
+        1.01E+04  8.49E+04  3.57E+05 -1.41E+05  9.90E+04  1.60E+04 -1.40E+05 -1.42E+05  1.02E+10
 
 TH10
+        4.17E+04  5.97E+04  5.92E+04 -3.06E+04  6.96E+04  1.13E+04 -5.58E+04 -3.08E+04 -5.95E+04  5.06E+09
 
 TH11
+       -9.00E+03 -7.59E+04 -3.19E+05  1.32E+05 -8.84E+04 -1.43E+04  1.25E+05 -1.21E+06  2.45E+05  5.32E+04  8.15E+09
 
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
 #CPUT: Total CPU Time in Seconds,       46.248
Stop Time:
Fri Sep 24 19:49:14 CDT 2021
