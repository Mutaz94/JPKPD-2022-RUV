Thu Sep 30 03:06:55 CDT 2021
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
$DATA ../../../../data/spa1/D/dat45.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   15228.2964881289        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.2642E+02  3.4952E+02 -8.7289E+01  3.1135E+02  1.9493E+02 -1.6592E+03 -8.3573E+02 -5.2341E+01 -1.3928E+03 -3.0175E+02
            -2.9873E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -562.768487053478        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0899E+00  1.0641E+00  1.0412E+00  2.0884E+00  1.3365E+00  3.1845E+00  1.7850E+00  9.3926E-01  3.2938E+00  1.0422E+00
             1.2353E+01
 PARAMETER:  1.8606E-01  1.6216E-01  1.4041E-01  8.3639E-01  3.9007E-01  1.2583E+00  6.7942E-01  3.7338E-02  1.2920E+00  1.4131E-01
             2.6139E+00
 GRADIENT:  -1.8002E+01  8.7359E+00 -3.8843E+01  7.3640E+01  2.5388E+00  1.2942E+02  8.1154E+00  5.3254E+00  5.6667E+01  3.0682E+00
             2.1939E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -590.050435623557        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.1579E+00  1.0985E+00  1.8744E+00  2.2129E+00  4.2050E+00  2.3392E+00  5.9270E+00  7.4443E-01  2.7612E+00  1.7842E+00
             1.2233E+01
 PARAMETER:  2.4657E-01  1.9394E-01  7.2827E-01  8.9433E-01  1.5363E+00  9.4982E-01  1.8795E+00 -1.9514E-01  1.1157E+00  6.7897E-01
             2.6041E+00
 GRADIENT:  -5.7678E+00  2.8099E+01 -1.7268E+00  7.7172E+01 -3.9431E+00  3.9669E+01  1.6259E+01 -8.6351E-02  5.6230E+01  3.1163E-01
             2.1861E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -633.966634020949        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  9.7530E-01  8.0373E-01  1.9348E+00  1.1600E+00  6.8680E+00  1.9019E+00  2.7978E+00  5.2884E+00  1.4606E+00  6.7708E+00
             1.0835E+01
 PARAMETER:  7.4985E-02 -1.1849E-01  7.6000E-01  2.4841E-01  2.0269E+00  7.4285E-01  1.1288E+00  1.7655E+00  4.7885E-01  2.0126E+00
             2.4828E+00
 GRADIENT:  -4.8256E+01  2.1818E+00  7.7134E+00 -7.1490E+01 -6.5590E+00  2.5424E+01 -3.0499E+01  2.6329E+01 -1.3390E+00  3.2971E+00
             2.0782E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -681.390516919018        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.9925E-01  3.5925E-01  9.8531E-01  1.5483E+00  8.0338E+00  1.5792E+00  5.9849E+00  2.6026E+00  1.2413E+00  8.3336E+00
             8.9815E+00
 PARAMETER:  9.9247E-02 -9.2374E-01  8.5201E-02  5.3718E-01  2.1837E+00  5.5694E-01  1.8892E+00  1.0565E+00  3.1613E-01  2.2203E+00
             2.2952E+00
 GRADIENT:   1.4852E+01  3.2277E+01  2.2386E+01  6.4097E+01 -2.6781E+00 -6.0176E+01  5.2887E+00 -2.3631E+01 -1.1455E+01  8.1215E+00
             9.8670E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -706.712445468210        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  7.8619E-01  1.0326E-01  2.1730E-01  1.0136E+00  1.7962E+01  1.5640E+00  6.0463E+00  2.3608E+00  6.8120E-01  6.9775E+00
             7.1402E+00
 PARAMETER: -1.4055E-01 -2.1705E+00 -1.4265E+00  1.1352E-01  2.9883E+00  5.4727E-01  1.8994E+00  9.5901E-01 -2.8390E-01  2.0427E+00
             2.0657E+00
 GRADIENT:   4.2080E+01  3.8073E+01 -1.4485E+01  1.1178E+02  1.0828E+00 -1.9971E+01  1.2440E+00  8.8964E+00 -1.8748E+01 -3.8460E+00
            -2.3695E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -719.354938108729        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  6.5793E-01  4.7409E-02  9.9597E-02  7.3615E-01  3.1114E+01  1.5914E+00  4.8548E+00  2.4372E+00  4.2411E-01  7.2130E+00
             6.9161E+00
 PARAMETER: -3.1866E-01 -2.9489E+00 -2.2066E+00 -2.0633E-01  3.5377E+00  5.6461E-01  1.6800E+00  9.9083E-01 -7.5777E-01  2.0759E+00
             2.0339E+00
 GRADIENT:   5.7732E+01  3.6402E+01 -2.3669E+01  1.0824E+02 -6.5422E-01 -3.7470E+00  1.2736E+01  2.9581E+01 -7.3922E+00 -4.6851E-03
            -2.4943E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -760.004244472220        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  7.1513E-01  3.5982E-02  1.2240E-01  7.7228E-01  2.3020E+01  1.7122E+00  2.5269E+00  2.3077E+00  4.7938E-01  4.5014E+00
             8.4335E+00
 PARAMETER: -2.3530E-01 -3.2247E+00 -2.0005E+00 -1.5841E-01  3.2364E+00  6.3779E-01  1.0270E+00  9.3624E-01 -6.3527E-01  1.6044E+00
             2.2322E+00
 GRADIENT:   3.5995E+01  2.4469E+00 -6.4833E+00 -3.5240E+00  9.3499E-02  9.5922E+00  2.2075E-01  1.7242E+01 -2.4387E+00 -3.5602E-02
            -7.2797E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -765.440233313155        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      799
 NPARAMETR:  6.0279E-01  2.0322E-02  7.5550E-02  6.0139E-01  1.0200E+01  1.6304E+00  1.9956E+00  1.5057E+00  6.1736E-01  2.2675E+00
             8.3792E+00
 PARAMETER: -4.0619E-01 -3.7960E+00 -2.4830E+00 -4.0852E-01  2.4224E+00  5.8882E-01  7.9095E-01  5.0922E-01 -3.8230E-01  9.1866E-01
             2.2258E+00
 GRADIENT:   3.3800E+00  2.5644E-02 -2.6530E+00  4.6723E+00 -1.5733E-01 -2.8296E+00  2.7850E-03 -2.9924E+00 -3.2410E+00 -1.0264E-02
            -1.0793E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -765.661683123770        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      978             RESET HESSIAN, TYPE I
 NPARAMETR:  5.9742E-01  2.3202E-02  7.4502E-02  5.9591E-01  1.1293E+01  1.6419E+00  4.2951E-01  1.4682E+00  6.8425E-01  2.3506E+00
             8.3556E+00
 PARAMETER: -4.1514E-01 -3.6635E+00 -2.4969E+00 -4.1767E-01  2.5242E+00  5.9587E-01 -7.4511E-01  4.8401E-01 -2.7943E-01  9.5467E-01
             2.2229E+00
 GRADIENT:   3.8203E+01 -3.3546E-02  4.0714E+01  1.6998E+01  2.3435E-02  1.4489E+01  1.9473E-03  8.7649E-01 -7.0618E-02 -5.7095E-03
             2.1289E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -765.681932583790        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1156
 NPARAMETR:  5.9807E-01  2.6537E-02  7.4571E-02  5.9525E-01  1.2578E+01  1.6391E+00  9.4636E-02  1.4753E+00  6.9419E-01  2.3813E+00
             8.3611E+00
 PARAMETER: -4.1405E-01 -3.5292E+00 -2.4960E+00 -4.1878E-01  2.6320E+00  5.9412E-01 -2.2577E+00  4.8885E-01 -2.6502E-01  9.6766E-01
             2.2236E+00
 GRADIENT:   4.7407E-02  3.7982E-02 -8.5600E-02  2.1399E-01 -7.6850E-03 -1.0781E-01  4.9317E-04 -1.3112E-02 -1.1950E-02 -1.2148E-03
            -4.4624E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -765.683029366746        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1322
 NPARAMETR:  5.9758E-01  2.5935E-02  7.4461E-02  5.9485E-01  1.2586E+01  1.6401E+00  1.3500E-02  1.4764E+00  6.9360E-01  2.6182E+00
             8.3596E+00
 PARAMETER: -4.1486E-01 -3.5522E+00 -2.4975E+00 -4.1945E-01  2.6326E+00  5.9476E-01 -4.2050E+00  4.8960E-01 -2.6586E-01  1.0625E+00
             2.2234E+00
 GRADIENT:  -3.8934E-02 -4.2764E-02 -2.2537E-01  4.2372E-01  1.3161E-02  9.9014E-02  7.5038E-06  1.0824E-01  3.8498E-02 -4.2370E-03
            -6.9753E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1322
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4078E-03 -2.2014E-05 -2.6776E-03 -2.2007E-02  2.4758E-05
 SE:             2.8568E-02  4.5177E-06  2.0117E-02  1.6429E-02  4.2569E-04
 N:                     100         100         100         100         100

 P VAL.:         8.4986E-01  1.1015E-06  8.9411E-01  1.8040E-01  9.5362E-01

 ETASHRINKSD(%)  4.2946E+00  9.9985E+01  3.2606E+01  4.4960E+01  9.8574E+01
 ETASHRINKVR(%)  8.4048E+00  1.0000E+02  5.4580E+01  6.9706E+01  9.9980E+01
 EBVSHRINKSD(%)  4.0167E+00  9.9978E+01  3.1440E+01  4.5595E+01  9.8806E+01
 EBVSHRINKVR(%)  7.8721E+00  1.0000E+02  5.2995E+01  7.0401E+01  9.9986E+01
 RELATIVEINF(%)  9.5700E+00  8.8494E-07  1.4442E+00  5.8474E-01  1.3294E-03
 EPSSHRINKSD(%)  1.1873E+01
 EPSSHRINKVR(%)  2.2336E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -765.68302936674593     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       153.25550383792677     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.60
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -765.683       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.98E-01  2.59E-02  7.45E-02  5.95E-01  1.26E+01  1.64E+00  1.35E-02  1.48E+00  6.94E-01  2.62E+00  8.36E+00
 


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
+        1.05E+03
 
 TH 2
+       -2.78E+02  2.62E+03
 
 TH 3
+       -1.21E+03 -4.30E+02  4.26E+04
 
 TH 4
+       -2.59E+02  1.53E+02 -8.01E+03  1.87E+03
 
 TH 5
+        2.05E-01 -1.63E+00 -1.11E+00  2.61E-01  3.01E-03
 
 TH 6
+        6.53E+00  4.32E+00  4.73E+01 -2.80E+01  1.04E-02  6.05E+01
 
 TH 7
+       -2.14E-02  6.11E-02 -4.68E-02 -1.95E-02 -1.61E-04 -2.77E-03  9.54E-02
 
 TH 8
+        4.51E+00 -2.18E+01 -1.55E+02 -5.51E+00 -1.28E-02  1.95E+00  9.83E-03  2.01E+01
 
 TH 9
+        3.42E+00 -4.86E+01  7.20E+01 -4.43E+01 -1.50E-02  2.76E+00 -5.90E-03  2.07E+01  4.49E+01
 
 TH10
+       -6.95E-03  6.93E-01 -8.90E-02 -3.05E-02 -6.54E-04  2.39E-03 -3.45E-03  7.39E-03  6.39E-03 -9.80E-04
 
 TH11
+       -1.55E+01 -4.56E+00  7.71E+01 -1.38E+01 -4.04E-03  1.40E+00  1.96E-04  1.41E+00  2.08E+00  3.00E-04  7.38E+00
 
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
 #CPUT: Total CPU Time in Seconds,       32.194
Stop Time:
Thu Sep 30 03:07:28 CDT 2021
