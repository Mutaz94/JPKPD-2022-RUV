Wed Sep 29 13:20:39 CDT 2021
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
$DATA ../../../../data/spa/A3/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   33.1982556083512        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8652E+02  5.0016E+01  9.2047E+01 -5.1866E+01  1.3491E+02  4.3287E+01 -6.1011E+01 -2.9298E+01 -1.5745E+02 -1.2126E+02
            -3.0181E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1254.64079492292        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0631E+00  1.0082E+00  9.8823E-01  1.1763E+00  1.0750E+00  8.0284E-01  9.9358E-01  9.4720E-01  1.1427E+00  8.5189E-01
             3.5814E+00
 PARAMETER:  1.6115E-01  1.0820E-01  8.8160E-02  2.6239E-01  1.7228E-01 -1.1960E-01  9.3561E-02  4.5755E-02  2.3337E-01 -6.0299E-02
             1.3758E+00
 GRADIENT:   9.6260E+01  7.1403E+00 -3.1878E+01  5.6804E+01  3.8246E+01 -3.9780E+01  2.5400E+00  5.3261E+00  6.9542E-02  9.2989E+00
            -6.7200E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1266.50397173804        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0616E+00  8.3798E-01  8.1719E-01  1.2762E+00  7.6278E-01  8.9415E-01  6.9426E-01  5.6538E-01  1.1591E+00  2.5358E-01
             3.9226E+00
 PARAMETER:  1.5978E-01 -7.6756E-02 -1.0188E-01  3.4388E-01 -1.7078E-01 -1.1880E-02 -2.6491E-01 -4.7026E-01  2.4767E-01 -1.2721E+00
             1.4667E+00
 GRADIENT:   5.8311E+01  5.6486E+01  3.1849E+01  8.6404E+01 -6.4405E+01 -2.9066E-01 -2.4576E+00  2.4213E+00  7.2442E-01  7.6507E-01
             4.9074E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1271.76694416101        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0300E+00  8.0638E-01  7.6865E-01  1.1935E+00  7.6687E-01  8.9131E-01  8.5878E-01  1.7875E-01  1.1062E+00  2.4898E-01
             3.8731E+00
 PARAMETER:  1.2951E-01 -1.1520E-01 -1.6312E-01  2.7685E-01 -1.6544E-01 -1.5065E-02 -5.2237E-02 -1.6218E+00  2.0093E-01 -1.2904E+00
             1.4541E+00
 GRADIENT:  -1.8847E+00  4.1565E-01 -1.3555E+00 -5.5629E-01  1.1914E-01  1.5364E+00  1.9853E-01  3.6750E-01 -1.0628E-01  1.3222E+00
             6.6103E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1272.62978351939        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0282E+00  6.2694E-01  6.6771E-01  1.2756E+00  6.3468E-01  8.9423E-01  1.0966E+00  4.6210E-02  1.0475E+00  8.3038E-02
             3.8184E+00
 PARAMETER:  1.2785E-01 -3.6691E-01 -3.0390E-01  3.4340E-01 -3.5464E-01 -1.1795E-02  1.9224E-01 -2.9746E+00  1.4642E-01 -2.3885E+00
             1.4398E+00
 GRADIENT:  -3.7658E+01  2.4220E+00  2.6892E+00 -1.5748E+01 -7.8769E+00 -7.9681E-01 -2.6384E-01  2.5290E-02 -1.1546E+00  1.1143E-01
            -1.3551E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1274.55198855935        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  1.0418E+00  4.3525E-01  7.6089E-01  1.4267E+00  6.4257E-01  8.8765E-01  1.3874E+00  1.1389E-02  9.5111E-01  1.9484E-02
             3.9190E+00
 PARAMETER:  1.4096E-01 -7.3183E-01 -1.7326E-01  4.5537E-01 -3.4228E-01 -1.9175E-02  4.2746E-01 -4.3751E+00  4.9879E-02 -3.8382E+00
             1.4658E+00
 GRADIENT:  -1.3087E+00  4.5784E+00  3.6957E+00  1.0510E+01 -6.0558E+00 -8.8960E-01 -3.0470E-01  1.7720E-03 -1.7591E+00  6.4409E-03
            -2.1637E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1276.09059077627        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  1.0342E+00  1.4539E-01  5.0075E-01  1.5003E+00  4.1878E-01  8.9841E-01  3.0426E+00  1.0000E-02  9.5350E-01  1.0000E-02
             3.8560E+00
 PARAMETER:  1.3363E-01 -1.8284E+00 -5.9165E-01  5.0569E-01 -7.7041E-01 -7.1262E-03  1.2127E+00 -1.1433E+01  5.2385E-02 -1.1071E+01
             1.4496E+00
 GRADIENT:  -8.1036E+00  3.3485E+00  1.8795E+01  3.6942E+01 -3.1170E+01 -1.1694E+00 -1.0614E+00  0.0000E+00 -1.6223E+00  0.0000E+00
            -2.4738E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1277.39340552452        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0328E+00  1.1240E-01  3.7671E-01  1.4161E+00  3.3897E-01  9.0598E-01  3.8085E+00  1.0000E-02  1.0196E+00  1.0000E-02
             3.7959E+00
 PARAMETER:  1.3226E-01 -2.0857E+00 -8.7628E-01  4.4788E-01 -9.8184E-01  1.2579E-03  1.4372E+00 -1.3816E+01  1.1946E-01 -1.3547E+01
             1.4339E+00
 GRADIENT:  -4.0301E-01  3.4026E+01 -5.5021E-01 -6.5477E+00 -1.9542E+01 -7.6084E+00  3.7538E+01  0.0000E+00 -1.1464E+01  0.0000E+00
            -1.8356E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1278.08797697830        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  1.0322E+00  8.7274E-02  3.4425E-01  1.3901E+00  3.1625E-01  9.1375E-01  4.4934E+00  1.0000E-02  1.0281E+00  1.0000E-02
             3.7825E+00
 PARAMETER:  1.3172E-01 -2.3387E+00 -9.6638E-01  4.2938E-01 -1.0512E+00  9.8058E-03  1.6026E+00 -1.5708E+01  1.2775E-01 -1.5553E+01
             1.4304E+00
 GRADIENT:   8.8641E-01  1.4483E+01  6.2047E+00 -2.5502E+00 -1.7567E+01 -2.6238E+00  1.6754E+01  0.0000E+00  1.7073E+00  0.0000E+00
            -1.5371E+01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1278.08797697830        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:     1094
 NPARAMETR:  1.0309E+00  8.5333E-02  3.4695E-01  1.3955E+00  3.1959E-01  9.1463E-01  4.4294E+00  1.0000E-02  1.0281E+00  1.0000E-02
             3.8284E+00
 PARAMETER:  1.3172E-01 -2.3387E+00 -9.6638E-01  4.2938E-01 -1.0512E+00  9.8058E-03  1.6026E+00 -1.5708E+01  1.2775E-01 -1.5553E+01
             1.4304E+00
 GRADIENT:   3.5650E+02  1.8119E+01 -8.0431E+01 -1.9418E+02 -5.4449E+01 -2.3865E+02  4.9015E+01  0.0000E+00  3.9620E+01  0.0000E+00
            -5.6234E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1094
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2565E-03  4.6895E-03  9.5666E-05 -2.2617E-02  1.4848E-04
 SE:             2.8125E-02  7.9337E-03  2.4725E-04  2.5721E-02  3.7308E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6437E-01  5.5446E-01  6.9881E-01  3.7923E-01  6.9064E-01

 ETASHRINKSD(%)  5.7769E+00  7.3421E+01  9.9172E+01  1.3831E+01  9.8750E+01
 ETASHRINKVR(%)  1.1220E+01  9.2936E+01  9.9993E+01  2.5748E+01  9.9984E+01
 EBVSHRINKSD(%)  5.2031E+00  8.0547E+01  9.9172E+01  1.3291E+01  9.8867E+01
 EBVSHRINKVR(%)  1.0135E+01  9.6216E+01  9.9993E+01  2.4816E+01  9.9987E+01
 RELATIVEINF(%)  8.3344E+01  1.2519E+00  2.2945E-04  2.3022E+01  4.0712E-04
 EPSSHRINKSD(%)  2.2050E+01
 EPSSHRINKVR(%)  3.9239E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1278.0879769782960     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -542.93715041455778     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.33
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1278.088       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  8.73E-02  3.44E-01  1.39E+00  3.16E-01  9.14E-01  4.49E+00  1.00E-02  1.03E+00  1.00E-02  3.78E+00
 


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
+        1.22E+05
 
 TH 2
+       -4.16E+02  6.41E+04
 
 TH 3
+        1.75E+02  3.61E+03  2.41E+04
 
 TH 4
+        6.52E+01  7.73E+02 -1.25E+03  6.42E+03
 
 TH 5
+       -2.66E+04 -4.40E+03  1.84E+03  7.27E+02  3.60E+04
 
 TH 6
+       -8.47E+04  1.72E+02 -1.30E+02 -9.57E+01 -2.06E+02  1.42E+05
 
 TH 7
+       -8.25E+00  2.72E+02  8.66E+01  1.76E+01 -9.17E+01  5.31E+00  5.08E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.96E+04 -1.06E+02  8.49E+01  6.07E+00 -2.79E+04 -1.01E+04 -1.95E+00  0.00E+00  1.29E+05
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.92E+00  2.79E+02 -1.50E+02 -3.24E+01  1.76E+02 -3.11E+00  8.31E+00  0.00E+00  1.02E+01  0.00E+00  9.43E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.716
Stop Time:
Wed Sep 29 13:21:03 CDT 2021
