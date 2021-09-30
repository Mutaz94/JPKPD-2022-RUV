Thu Sep 30 02:51:07 CDT 2021
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
$DATA ../../../../data/spa1/D/dat26.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   10195.4215107934        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0427E+02  5.1558E+01 -4.4021E+01 -4.5411E+01  2.4060E+02 -9.2345E+02 -4.0919E+02 -2.5615E+02 -8.4262E+02 -2.8591E+02
            -2.1698E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -783.351196889603        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.3139E+00  1.1686E+00  1.1579E+00  1.9779E+00  1.0552E+00  2.3009E+00  1.6345E+00  9.4758E-01  2.6084E+00  1.1951E+00
             1.2932E+01
 PARAMETER:  3.7301E-01  2.5584E-01  2.4660E-01  7.8204E-01  1.5372E-01  9.3329E-01  5.9135E-01  4.6156E-02  1.0587E+00  2.7825E-01
             2.6597E+00
 GRADIENT:  -3.1130E+01  1.4110E+01 -1.8249E+01  5.1044E+01 -1.0567E+01  6.3896E+01  4.7121E+00  5.0900E+00  3.9780E+01  7.5183E+00
             3.5766E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -808.697734140377        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  1.3815E+00  1.2545E+00  8.3857E+00  2.4512E+00  1.0732E+01  1.7340E+00  6.1465E+00  5.6254E-01  2.7034E+00  4.7572E+00
             1.1929E+01
 PARAMETER:  4.2319E-01  3.2673E-01  2.2265E+00  9.9657E-01  2.4732E+00  6.5043E-01  1.9159E+00 -4.7529E-01  1.0945E+00  1.6597E+00
             2.5789E+00
 GRADIENT:   2.4871E+01  2.9473E+01  2.2159E+00  6.1901E+01  3.4695E-01 -2.0299E+00  4.0021E+01 -1.5823E-02  5.6428E+01 -4.0444E-02
             3.1915E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -838.928421943943        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      239
 NPARAMETR:  1.2738E+00  1.1577E+00  1.0609E+01  1.9561E+00  8.4399E+00  1.8271E+00  4.3026E+00  1.4199E+00  2.4365E+00  1.3051E+01
             1.0534E+01
 PARAMETER:  3.4202E-01  2.4641E-01  2.4617E+00  7.7095E-01  2.2330E+00  7.0273E-01  1.5592E+00  4.5057E-01  9.9056E-01  2.6689E+00
             2.4546E+00
 GRADIENT:  -4.7741E+00  1.8258E+01  1.3401E+00  4.3208E+01 -7.5693E+00  2.8813E+00  2.8404E+01 -3.2995E-02  5.4355E+01  3.2787E+01
             2.7143E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -938.163637925718        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  1.1097E+00  1.1267E+00  3.0864E+01  1.1538E+00  7.4093E+00  1.9295E+00  1.6100E+00  5.8188E+00  2.8249E+00  7.3605E+00
             6.9850E+00
 PARAMETER:  2.0405E-01  2.1926E-01  3.5296E+00  2.4305E-01  2.1027E+00  7.5726E-01  5.7621E-01  1.8611E+00  1.1385E+00  2.0961E+00
             2.0438E+00
 GRADIENT:  -1.8835E+01 -1.8349E+01  4.9794E+00 -2.5838E+01 -1.8035E+00  2.4793E+01  5.9682E+00  2.3513E+00  3.1295E+01  5.8796E-01
             6.7867E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -938.778062199638        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      434
 NPARAMETR:  1.1146E+00  1.1150E+00  2.8859E+01  1.1688E+00  7.3244E+00  1.9111E+00  1.6541E+00  5.7742E+00  2.7751E+00  7.4141E+00
             6.9669E+00
 PARAMETER:  2.0849E-01  2.0882E-01  3.4624E+00  2.5597E-01  2.0912E+00  7.4769E-01  6.0324E-01  1.8534E+00  1.1207E+00  2.1034E+00
             2.0412E+00
 GRADIENT:  -3.7401E+01 -1.8315E+01  3.9522E+00 -3.0217E+01 -7.2195E+00 -1.3575E+01  3.6402E+00 -9.4869E+00  1.0403E+01  5.8278E+00
             4.1681E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -950.864411506241        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      615
 NPARAMETR:  1.2170E+00  5.5793E-01  2.7506E+00  1.8675E+00  5.3100E+00  1.6322E+00  3.9938E+00  4.9171E+00  1.6415E+00  5.8134E+00
             6.7392E+00
 PARAMETER:  2.9637E-01 -4.8352E-01  1.1118E+00  7.2462E-01  1.7696E+00  5.8993E-01  1.4847E+00  1.6927E+00  5.9564E-01  1.8602E+00
             2.0079E+00
 GRADIENT:   4.0445E+01  3.1345E+01  6.4223E-01  6.9819E+01 -2.9678E+00 -3.6787E+01 -1.6437E+00  1.3254E+01  3.1481E+00  4.6783E-01
            -4.7242E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -993.739738256861        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      794
 NPARAMETR:  1.0700E+00  1.6418E-01  6.0436E-01  1.4508E+00  6.7149E+00  1.7511E+00  2.8628E+00  2.1375E+00  1.0004E+00  9.9333E+00
             7.1207E+00
 PARAMETER:  1.6762E-01 -1.7068E+00 -4.0359E-01  4.7215E-01  2.0043E+00  6.6024E-01  1.1518E+00  8.5963E-01  1.0042E-01  2.3959E+00
             2.0630E+00
 GRADIENT:   4.0456E+01  1.3412E+01  5.9788E+00  2.6129E+01 -2.0472E+01  2.7997E+01 -1.0366E+00 -2.8192E+01 -1.4178E+01  3.7103E+01
             1.9174E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1018.11452978465        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      971
 NPARAMETR:  8.7667E-01  8.8707E-02  3.5441E-01  1.1603E+00  7.0821E+00  1.4974E+00  1.4985E+00  3.2801E+00  7.1018E-01  9.0521E+00
             6.4974E+00
 PARAMETER: -3.1623E-02 -2.3224E+00 -9.3731E-01  2.4871E-01  2.0576E+00  5.0373E-01  5.0443E-01  1.2879E+00 -2.4224E-01  2.3030E+00
             1.9714E+00
 GRADIENT:  -1.5737E+01  1.2580E+01 -9.5645E+00 -5.7652E+00 -2.9186E+01  2.0008E+01 -3.8586E-01  2.0736E+01  4.3213E+00  3.4593E+01
            -4.5411E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1042.23612213202        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1146
 NPARAMETR:  8.6569E-01  2.7278E-02  2.1953E-01  9.8479E-01  6.9814E+00  1.2862E+00  5.5456E-01  2.2559E+00  3.2516E-01  6.0515E+00
             6.9484E+00
 PARAMETER: -4.4227E-02 -3.5017E+00 -1.4163E+00  8.4671E-02  2.0432E+00  3.5171E-01 -4.8958E-01  9.1356E-01 -1.0234E+00  1.9003E+00
             2.0385E+00
 GRADIENT:   5.5569E+01  1.1167E+00  2.2649E+01 -5.1162E+01 -1.8319E+01 -2.3591E+01 -2.2842E-03 -5.5866E+00 -1.7856E+00  1.8385E+01
             1.2129E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1064.27641670283        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1324
 NPARAMETR:  7.0358E-01  1.0000E-02  1.0751E-01  7.4256E-01  8.2213E+00  1.3078E+00  1.6376E-01  1.7812E+00  9.9390E-02  4.3368E+00
             6.7694E+00
 PARAMETER: -2.5158E-01 -5.2859E+00 -2.1301E+00 -1.9765E-01  2.2067E+00  3.6834E-01 -1.7093E+00  6.7729E-01 -2.2087E+00  1.5671E+00
             2.0124E+00
 GRADIENT:  -1.5419E+01  0.0000E+00 -7.5831E+00  2.4751E+01 -3.1920E+00 -1.2975E+01 -1.7384E-05  1.6675E+00  1.0956E-01 -5.8228E+00
             7.6580E+00

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1064.27641670283        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  7.0279E-01  1.0000E-02  1.0869E-01  7.4168E-01  8.3079E+00  1.3108E+00  1.6659E-01  1.7868E+00  1.0049E-01  4.3074E+00
             6.7906E+00
 PARAMETER: -2.5158E-01 -5.2859E+00 -2.1301E+00 -1.9765E-01  2.2067E+00  3.6834E-01 -1.7093E+00  6.7729E-01 -2.2087E+00  1.5671E+00
             2.0124E+00
 GRADIENT:   1.3936E+02  0.0000E+00 -1.2204E+01  1.1969E+02 -6.7215E+01 -6.5631E+01 -1.1395E-04 -5.4798E+01 -1.7558E+01  5.5954E+01
            -8.7863E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1355
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5572E-02 -5.9458E-05  8.8514E-04 -6.8971E-03 -1.0813E-02
 SE:             2.9889E-02  1.8875E-05  2.4928E-02  3.0442E-03  6.2109E-03
 N:                     100         100         100         100         100

 P VAL.:         6.0236E-01  1.6322E-03  9.7167E-01  2.3471E-02  8.1691E-02

 ETASHRINKSD(%)  1.0000E-10  9.9937E+01  1.6489E+01  8.9802E+01  7.9193E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  3.0259E+01  9.8960E+01  9.5671E+01
 EBVSHRINKSD(%)  3.5342E+00  9.9941E+01  1.4448E+01  9.0487E+01  8.3419E+01
 EBVSHRINKVR(%)  6.9435E+00  1.0000E+02  2.6808E+01  9.9095E+01  9.7251E+01
 RELATIVEINF(%)  2.5585E+01  1.2194E-06  1.1869E+00  7.2841E-03  8.2526E-01
 EPSSHRINKSD(%)  1.3416E+01
 EPSSHRINKVR(%)  2.5032E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1064.2764167028326     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -145.33788349815995     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.89
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1064.276       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         7.04E-01  1.00E-02  1.08E-01  7.43E-01  8.22E+00  1.31E+00  1.64E-01  1.78E+00  9.94E-02  4.34E+00  6.77E+00
 


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
+        3.23E+04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.83E+02  0.00E+00  6.93E+04
 
 TH 4
+       -5.77E+02  0.00E+00 -5.43E+03  4.71E+04
 
 TH 5
+        2.49E+00  0.00E+00  3.67E+02  1.27E+00  3.92E+00
 
 TH 6
+       -1.14E+04  0.00E+00  4.37E+01 -1.76E+01 -1.50E-01  4.30E+03
 
 TH 7
+        4.70E-04  0.00E+00 -2.76E-02  2.33E-02  3.36E-04  3.68E-03 -3.47E-02
 
 TH 8
+        1.50E+01  0.00E+00  4.39E+00 -5.34E+01  1.23E+00  2.79E+00  1.11E-03  7.11E+02
 
 TH 9
+       -2.51E+04  0.00E+00  1.25E+01  5.00E-02 -8.56E-01  9.23E+03 -2.79E-02  4.96E+00  2.03E+04
 
 TH10
+        8.19E-01  0.00E+00 -1.06E+03  1.61E+03 -1.38E+01 -1.03E+00 -4.63E-04 -2.00E+02 -2.28E+00  5.62E+01
 
 TH11
+       -5.70E+01  0.00E+00  5.85E+02 -5.72E+01  6.33E+00 -1.43E+00 -1.66E-03  8.62E+00  3.29E+01 -1.65E+01  2.55E+01
 
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
 #CPUT: Total CPU Time in Seconds,       34.322
Stop Time:
Thu Sep 30 02:51:43 CDT 2021
