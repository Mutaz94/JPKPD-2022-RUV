Wed Sep 29 12:08:03 CDT 2021
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
$DATA ../../../../data/spa/A1/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1460.23341238595        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3239E+02  2.7480E+01  2.7390E+01  3.8715E+01  6.3984E+01  5.2477E+01  6.6152E+00 -9.5084E-01  2.1732E+01 -5.0587E+01
            -3.2954E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1497.89097972158        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.5162E-01  9.4705E-01  8.6229E-01  1.0431E+00  8.7982E-01  8.8454E-01  9.1728E-01  9.6002E-01  8.3244E-01  1.1795E+00
             2.3934E+00
 PARAMETER:  5.0408E-02  4.5592E-02 -4.8169E-02  1.4222E-01 -2.8041E-02 -2.2690E-02  1.3656E-02  5.9199E-02 -8.3399E-02  2.6508E-01
             9.7270E-01
 GRADIENT:  -4.9380E+01 -9.2615E+00 -1.7524E+01 -3.9096E+00 -7.0795E+00 -1.7097E+01  9.5138E+00  1.0806E+01  1.3689E+01  2.4343E+01
             1.5523E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1509.40160637583        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.4482E-01  6.5570E-01  5.0489E-01  1.2333E+00  5.2914E-01  9.3888E-01  6.2320E-01  3.2172E-01  7.7669E-01  7.4401E-01
             2.2104E+00
 PARAMETER:  4.3237E-02 -3.2205E-01 -5.8342E-01  3.0972E-01 -5.3651E-01  3.6936E-02 -3.7289E-01 -1.0341E+00 -1.5272E-01 -1.9570E-01
             8.9317E-01
 GRADIENT:  -6.8614E+01  6.4153E+01 -1.0509E+01  2.1131E+02  5.5580E+00 -3.8907E+00 -1.3624E+00  2.3334E+00 -2.5354E-01  3.2274E+00
             1.1668E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1540.89040764428        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      343
 NPARAMETR:  9.7623E-01  6.3837E-01  7.3171E-01  1.2150E+00  6.6934E-01  9.3572E-01  1.6597E+00  2.5725E-01  6.6328E-01  9.4848E-01
             1.6484E+00
 PARAMETER:  7.5940E-02 -3.4883E-01 -2.1237E-01  2.9472E-01 -3.0147E-01  3.3561E-02  6.0663E-01 -1.2577E+00 -3.1056E-01  4.7107E-02
             5.9979E-01
 GRADIENT:  -2.9724E+01  2.7139E+01  1.8689E+01 -1.9866E+01 -2.9574E+01 -7.1142E+00  1.3861E+01  6.6030E-01 -1.0731E+01  6.8002E+00
             2.0911E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1547.75575021408        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.8709E-01  4.4165E-01  5.5699E-01  1.2913E+00  5.0808E-01  9.5783E-01  1.8510E+00  4.6452E-02  7.0480E-01  7.3170E-01
             1.5901E+00
 PARAMETER:  8.7002E-02 -7.1725E-01 -4.8521E-01  3.5562E-01 -5.7711E-01  5.6911E-02  7.1571E-01 -2.9693E+00 -2.4984E-01 -2.1238E-01
             5.6381E-01
 GRADIENT:  -3.8069E+00  1.6227E+01  1.0535E+01  2.9608E+01 -1.9122E+01  1.1986E+00  1.9758E+00  1.4316E-02  3.8774E+00 -2.7385E+00
             4.2481E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1550.02677335702        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  9.8486E-01  2.2408E-01  6.0112E-01  1.3915E+00  4.9537E-01  9.5875E-01  2.8667E+00  1.0000E-02  6.5121E-01  7.9065E-01
             1.5876E+00
 PARAMETER:  8.4747E-02 -1.3958E+00 -4.0895E-01  4.3035E-01 -6.0246E-01  5.7872E-02  1.1531E+00 -6.0414E+00 -3.2893E-01 -1.3490E-01
             5.6222E-01
 GRADIENT:   7.2679E+00  1.6066E+00  7.2348E+00 -8.6248E+00 -1.0336E+01  3.2955E+00  1.4915E+00  0.0000E+00 -1.9969E+00 -1.5815E-01
            -1.2356E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1550.70778249768        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      873
 NPARAMETR:  9.7488E-01  1.2233E-01  6.7715E-01  1.4708E+00  5.2894E-01  9.4200E-01  3.8824E+00  1.0000E-02  6.4596E-01  8.5486E-01
             1.5923E+00
 PARAMETER:  7.4564E-02 -2.0010E+00 -2.8986E-01  4.8581E-01 -5.3689E-01  4.0246E-02  1.4565E+00 -8.3828E+00 -3.3702E-01 -5.6819E-02
             5.6517E-01
 GRADIENT:  -6.4092E+00  1.3383E+00 -5.9003E+00  2.1759E+01  8.2066E+00 -2.3570E+00 -4.8755E-01  0.0000E+00  1.8309E+00 -8.2026E-01
            -8.1812E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1551.29697723780        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  9.7470E-01  7.1975E-02  6.7493E-01  1.4884E+00  5.1756E-01  9.4533E-01  4.9463E+00  1.0000E-02  6.2773E-01  8.5654E-01
             1.6062E+00
 PARAMETER:  7.4369E-02 -2.5314E+00 -2.9314E-01  4.9771E-01 -5.5862E-01  4.3777E-02  1.6986E+00 -1.1435E+01 -3.6564E-01 -5.4857E-02
             5.7385E-01
 GRADIENT:  -2.5337E+00  2.2231E+00 -2.5963E+00  2.5130E+00  3.3879E+00 -8.7314E-01  4.0357E+00  0.0000E+00 -3.9207E+00 -1.2004E-01
             1.0197E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1551.79003209090        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  9.7511E-01  4.7647E-02  6.6541E-01  1.4917E+00  5.0872E-01  9.4702E-01  5.8792E+00  1.0000E-02  6.2601E-01  8.4005E-01
             1.6049E+00
 PARAMETER:  7.4794E-02 -2.9439E+00 -3.0735E-01  4.9994E-01 -5.7586E-01  4.5561E-02  1.8714E+00 -1.3842E+01 -3.6839E-01 -7.4296E-02
             5.7304E-01
 GRADIENT:   3.4742E-01  1.8074E-01  1.0834E+00  1.3348E-01 -8.7055E-01  1.0884E-01 -3.9703E-01  0.0000E+00  6.7265E-01 -3.0758E-01
            -2.0887E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1552.50005548449        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1416
 NPARAMETR:  9.7262E-01  2.5935E-02  6.4102E-01  1.4875E+00  4.9172E-01  9.4450E-01  7.5457E+00  1.0000E-02  6.0828E-01  8.0008E-01
             1.6145E+00
 PARAMETER:  7.2239E-02 -3.5522E+00 -3.4470E-01  4.9713E-01 -6.0985E-01  4.2895E-02  2.1210E+00 -1.7615E+01 -3.9712E-01 -1.2304E-01
             5.7902E-01
 GRADIENT:  -4.9737E+00 -3.2325E+00  5.7534E+00 -4.9454E-01 -8.0185E+00 -6.1419E-01 -5.1298E+00  0.0000E+00 -1.1437E+00 -8.0690E-01
             1.3686E-01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1552.50005548449        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1445
 NPARAMETR:  9.7345E-01  2.5844E-02  6.4015E-01  1.4883E+00  4.9142E-01  9.4508E-01  7.5311E+00  1.0000E-02  6.0893E-01  8.0107E-01
             1.6144E+00
 PARAMETER:  7.2239E-02 -3.5522E+00 -3.4470E-01  4.9713E-01 -6.0985E-01  4.2895E-02  2.1210E+00 -1.7615E+01 -3.9712E-01 -1.2304E-01
             5.7902E-01
 GRADIENT:  -4.7745E+00  4.4717E+01  4.5044E+00 -1.7182E+02  1.3603E+02 -5.4711E-01  7.0224E+01  0.0000E+00 -8.8614E-01 -7.3076E-01
             1.4417E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1445
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4874E-03  2.5078E-02 -1.4328E-04 -1.2970E-02 -6.6626E-03
 SE:             2.9664E-02  1.0722E-02  2.5830E-04  2.7034E-02  2.3941E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3317E-01  1.9340E-02  5.7909E-01  6.3139E-01  7.8078E-01

 ETASHRINKSD(%)  6.2122E-01  6.4080E+01  9.9135E+01  9.4339E+00  1.9796E+01
 ETASHRINKVR(%)  1.2386E+00  8.7098E+01  9.9993E+01  1.7978E+01  3.5673E+01
 EBVSHRINKSD(%)  1.0670E+00  7.6511E+01  9.9206E+01  9.1638E+00  1.8008E+01
 EBVSHRINKVR(%)  2.1226E+00  9.4483E+01  9.9994E+01  1.7488E+01  3.2774E+01
 RELATIVEINF(%)  9.7559E+01  3.9025E+00  3.4233E-04  3.2899E+01  3.7858E+00
 EPSSHRINKSD(%)  3.7970E+01
 EPSSHRINKVR(%)  6.1523E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1552.5000554844905     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -817.34922892075235     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1552.500       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.73E-01  2.59E-02  6.41E-01  1.49E+00  4.92E-01  9.44E-01  7.55E+00  1.00E-02  6.08E-01  8.00E-01  1.61E+00
 


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
+        1.29E+03
 
 TH 2
+       -2.03E+02  5.08E+05
 
 TH 3
+       -2.50E+00 -1.18E+02  1.77E+03
 
 TH 4
+       -3.63E+01  1.42E+03 -2.14E+02  9.06E+03
 
 TH 5
+        7.74E+00 -1.51E+03 -3.58E+04 -1.98E+04  5.37E+04
 
 TH 6
+        4.99E+00  1.91E+01  4.65E+00 -5.16E+00 -1.54E+00  2.14E+02
 
 TH 7
+       -7.89E-01  1.47E+03 -1.12E+00  6.19E+00 -6.79E+00  1.59E-01  2.63E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.77E+00  7.96E+02  1.09E+02  1.22E+04 -3.03E+04  2.89E+00  5.04E+00  0.00E+00  5.02E+02
 
 TH10
+       -1.07E+01  1.27E+02 -2.46E+01 -1.90E+00 -7.42E+04  3.57E-01  9.36E-01  0.00E+00  6.15E+01  2.02E+02
 
 TH11
+       -1.22E+01  6.14E+01 -7.83E+00 -1.86E+01 -7.79E+03  2.54E+00  3.74E-01  0.00E+00  4.04E+01  4.12E+01  9.71E+01
 
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
 #CPUT: Total CPU Time in Seconds,       25.539
Stop Time:
Wed Sep 29 12:08:30 CDT 2021
