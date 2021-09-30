Wed Sep 29 11:20:31 CDT 2021
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
$DATA ../../../../data/spa/B/dat54.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1714.89921749156        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0199E+02  6.5382E+01 -2.4017E+01  1.4107E+02 -4.2988E+00  4.8540E+01 -8.2073E-01  1.2885E+01  2.6012E+00  1.2459E+01
             1.3571E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1723.92443692200        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      196
 NPARAMETR:  1.0636E+00  1.0434E+00  1.1762E+00  9.4386E-01  1.1007E+00  1.0835E+00  1.0455E+00  8.8857E-01  1.0643E+00  9.6965E-01
             9.9022E-01
 PARAMETER:  1.6168E-01  1.4249E-01  2.6231E-01  4.2222E-02  1.9595E-01  1.8022E-01  1.4450E-01 -1.8137E-02  1.6232E-01  6.9175E-02
             9.0172E-02
 GRADIENT:   1.6173E+01 -5.5504E+00  8.0371E+00 -1.4284E+01 -1.9922E+00 -1.7834E+00  2.4945E+00  1.5066E+00  9.3272E-01 -1.0271E+01
             1.5615E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1724.93751324273        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0641E+00  1.0828E+00  1.1708E+00  9.2692E-01  1.1231E+00  1.0894E+00  8.1564E-01  6.3013E-01  1.1658E+00  1.1399E+00
             9.9185E-01
 PARAMETER:  1.6210E-01  1.7957E-01  2.5770E-01  2.4117E-02  2.1611E-01  1.8560E-01 -1.0379E-01 -3.6183E-01  2.5344E-01  2.3093E-01
             9.1820E-02
 GRADIENT:   1.6429E+01  3.0844E+00  7.5624E+00  2.2775E+00 -1.0676E+01  4.1397E-01  2.6348E+00 -1.9124E+00  4.5785E+00  2.1347E+00
             5.0357E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1725.56127419748        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0516E+00  1.3478E+00  9.0301E-01  7.6282E-01  1.1431E+00  1.0897E+00  8.5969E-01  5.5082E-01  1.2329E+00  1.0554E+00
             9.7043E-01
 PARAMETER:  1.5035E-01  3.9849E-01 -2.0253E-03 -1.7073E-01  2.3372E-01  1.8589E-01 -5.1180E-02 -4.9635E-01  3.0935E-01  1.5390E-01
             6.9984E-02
 GRADIENT:  -1.0172E+01  1.5677E+01  1.0720E+00  1.2779E+01 -4.8304E+00 -4.0637E-01 -5.7869E-01  8.9052E-01 -3.0976E+00  3.7469E-02
            -3.6806E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1726.17938492719        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  1.0597E+00  1.6174E+00  7.3383E-01  5.8704E-01  1.2374E+00  1.0928E+00  7.7036E-01  2.8751E-01  1.5160E+00  1.0869E+00
             9.8819E-01
 PARAMETER:  1.5802E-01  5.8082E-01 -2.0948E-01 -4.3267E-01  3.1302E-01  1.8870E-01 -1.6090E-01 -1.1465E+00  5.1609E-01  1.8334E-01
             8.8120E-02
 GRADIENT:   2.5724E+00  1.0365E+01 -3.9222E+00  1.5097E+01  8.3659E+00  4.8518E-01  1.1322E+00  3.1994E-01  4.5843E-01 -1.0941E+00
             2.4403E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1726.53975138356        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      911
 NPARAMETR:  1.0565E+00  1.6238E+00  7.2935E-01  5.6901E-01  1.2366E+00  1.0893E+00  7.4849E-01  1.0106E-01  1.5463E+00  1.1046E+00
             9.8207E-01
 PARAMETER:  1.5493E-01  5.8477E-01 -2.1561E-01 -4.6386E-01  3.1239E-01  1.8556E-01 -1.8969E-01 -2.1920E+00  5.3586E-01  1.9945E-01
             8.1912E-02
 GRADIENT:  -2.9430E+00 -5.2967E+00 -4.1980E-01  2.3404E+00 -4.0184E-01 -7.5737E-01 -3.3001E-02  3.3324E-02 -4.8373E-01  4.2729E-01
            -3.5020E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1726.60185875178        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1100             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0612E+00  1.6508E+00  7.1482E-01  5.4613E-01  1.2484E+00  1.0952E+00  7.3651E-01  2.4001E-02  1.5947E+00  1.1095E+00
             9.8315E-01
 PARAMETER:  1.5939E-01  6.0129E-01 -2.3573E-01 -5.0490E-01  3.2185E-01  1.9098E-01 -2.0583E-01 -3.6296E+00  5.6666E-01  2.0387E-01
             8.3002E-02
 GRADIENT:   9.0508E+02  8.3761E+02  9.5747E-01  1.5686E+02  1.8276E+01  1.7104E+02  1.3853E+01  3.0514E-03  3.4344E+01  1.7013E+00
             8.1100E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1726.60574579219        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1284
 NPARAMETR:  1.0615E+00  1.6509E+00  7.1700E-01  5.4661E-01  1.2500E+00  1.0956E+00  7.3514E-01  1.0000E-02  1.5995E+00  1.1107E+00
             9.8374E-01
 PARAMETER:  1.5972E-01  6.0134E-01 -2.3268E-01 -5.0402E-01  3.2312E-01  1.9132E-01 -2.0769E-01 -4.6912E+00  5.6966E-01  2.0500E-01
             8.3608E-02
 GRADIENT:   6.3691E+00 -1.3681E+01 -5.4730E-01 -8.5067E-01  1.1897E-01  1.6406E+00  1.5133E-01  0.0000E+00  2.2110E-01  3.0821E-02
             2.0279E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1726.60737424039        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1474             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0615E+00  1.6517E+00  7.2007E-01  5.4677E-01  1.2506E+00  1.0956E+00  7.3377E-01  1.0000E-02  1.6009E+00  1.1119E+00
             9.8352E-01
 PARAMETER:  1.5971E-01  6.0181E-01 -2.2841E-01 -5.0373E-01  3.2360E-01  1.9129E-01 -2.0956E-01 -4.6912E+00  5.7056E-01  2.0607E-01
             8.3379E-02
 GRADIENT:   9.0703E+02  8.4090E+02  1.5079E+00  1.5717E+02  1.7955E+01  1.7124E+02  1.3976E+01  0.0000E+00  3.4972E+01  1.4642E+00
             7.1260E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1726.60737424039        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     1535
 NPARAMETR:  1.0615E+00  1.6508E+00  7.1843E-01  5.4700E-01  1.2531E+00  1.0955E+00  7.3244E-01  1.0000E-02  1.6054E+00  1.1140E+00
             9.8410E-01
 PARAMETER:  1.5971E-01  6.0181E-01 -2.2841E-01 -5.0373E-01  3.2360E-01  1.9129E-01 -2.0956E-01 -4.6912E+00  5.7056E-01  2.0607E-01
             8.3379E-02
 GRADIENT:   2.0931E-02  3.9959E-01  1.5646E-01 -7.3523E-02 -8.1921E-01  4.9862E-03  5.9553E-02  0.0000E+00 -1.7670E-01 -1.0848E-01
            -8.7094E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1535
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.9160E-05 -3.8047E-02 -2.5173E-04  2.5884E-02 -3.8216E-02
 SE:             2.9877E-02  2.0822E-02  9.8400E-05  2.3582E-02  2.3490E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9735E-01  6.7660E-02  1.0521E-02  2.7236E-01  1.0376E-01

 ETASHRINKSD(%)  1.0000E-10  3.0244E+01  9.9670E+01  2.0998E+01  2.1304E+01
 ETASHRINKVR(%)  1.0000E-10  5.1341E+01  9.9999E+01  3.7587E+01  3.8070E+01
 EBVSHRINKSD(%)  3.5652E-01  2.8649E+01  9.9710E+01  2.2660E+01  1.9293E+01
 EBVSHRINKVR(%)  7.1177E-01  4.9090E+01  9.9999E+01  4.0186E+01  3.4864E+01
 RELATIVEINF(%)  9.9204E+01  3.6997E+00  1.7280E-04  4.9598E+00  1.8734E+01
 EPSSHRINKSD(%)  4.3572E+01
 EPSSHRINKVR(%)  6.8159E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1726.6073742403935     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -991.45654767665530     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1726.607       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.65E+00  7.20E-01  5.47E-01  1.25E+00  1.10E+00  7.34E-01  1.00E-02  1.60E+00  1.11E+00  9.84E-01
 


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
+        8.16E+02
 
 TH 2
+       -4.95E+00  4.14E+02
 
 TH 3
+        6.19E+00  1.14E+02  1.94E+02
 
 TH 4
+       -1.05E+01  3.88E+02 -1.29E+02  8.80E+02
 
 TH 5
+       -2.39E+00 -1.51E+02 -1.68E+02  1.31E+02  3.77E+02
 
 TH 6
+       -1.21E-01 -5.93E-01  1.54E+00 -2.51E+00 -2.97E-01  1.64E+02
 
 TH 7
+        6.09E-01 -3.69E+00  1.78E+01 -2.58E+01 -2.29E+01 -3.27E-01  9.00E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.09E-01 -2.25E+01 -2.23E+01  5.19E+01  2.61E+00 -4.24E-01  2.26E+01  0.00E+00  3.59E+01
 
 TH10
+        1.50E-02 -8.14E+00 -2.42E+01 -5.45E+00 -4.50E+01  2.01E-01  2.09E+01  0.00E+00  1.37E+00  6.70E+01
 
 TH11
+       -5.87E+00 -2.62E+01 -3.12E+01  5.43E+00  3.48E+00  2.12E+00  1.08E+01  0.00E+00  3.21E+00  2.13E+01  2.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       29.646
Stop Time:
Wed Sep 29 11:21:02 CDT 2021
