Sat Sep 18 13:02:52 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat77.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1632.43525124697        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.6753E+01 -1.6210E+01  2.8577E+00 -3.3826E+00 -1.6752E+01 -1.5680E+01 -4.0162E+00  5.2771E+00  4.4527E+01 -7.7663E+00
            -1.2028E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1653.74958057817        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0283E+00  1.1289E+00  1.0402E+00  8.8166E-01  1.1296E+00  1.0281E+00  1.1375E+00  9.0596E-01  6.4626E-01  1.0646E+00
             1.2483E+00
 PARAMETER:  1.2787E-01  2.2124E-01  1.3939E-01 -2.5944E-02  2.2182E-01  1.2773E-01  2.2887E-01  1.2417E-03 -3.3655E-01  1.6257E-01
             3.2179E-01
 GRADIENT:   2.6185E+01 -6.0871E+01  1.1248E+01 -9.2104E+01  3.3995E+00 -1.5325E+00 -1.5497E+00  2.9386E-01  2.8697E-02 -5.8643E+00
            -5.5377E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1656.52921383832        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0165E+00  9.3395E-01  9.7048E-01  1.0340E+00  9.8353E-01  1.0129E+00  1.3603E+00  7.3834E-01  4.8685E-01  1.0340E+00
             1.2607E+00
 PARAMETER:  1.1635E-01  3.1664E-02  7.0037E-02  1.3345E-01  8.3397E-02  1.1286E-01  4.0774E-01 -2.0335E-01 -6.1980E-01  1.3347E-01
             3.3165E-01
 GRADIENT:  -5.2994E+00 -6.6914E+00  1.3114E+00 -5.1336E+00 -5.5669E+00 -8.8639E+00 -7.9630E+00  2.2998E+00 -8.0027E+00  5.2483E+00
             1.0386E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1659.08847791210        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0208E+00  9.8385E-01  6.5792E-01  9.8901E-01  7.9732E-01  1.0443E+00  1.2858E+00  2.4348E-01  6.1679E-01  7.7909E-01
             1.2504E+00
 PARAMETER:  1.2056E-01  8.3721E-02 -3.1868E-01  8.8947E-02 -1.2650E-01  1.4334E-01  3.5134E-01 -1.3127E+00 -3.8323E-01 -1.4962E-01
             3.2344E-01
 GRADIENT:  -1.3973E+00  1.0583E+01  9.5745E+00 -4.9433E+00 -1.4849E+01  2.9816E+00  9.9044E-01  4.3462E-01  1.5809E+00 -2.8365E-01
             2.4426E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1659.47811278490        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0212E+00  8.6590E-01  5.5812E-01  1.0415E+00  6.8069E-01  1.0348E+00  1.3989E+00  1.3617E-01  5.8121E-01  6.6704E-01
             1.2342E+00
 PARAMETER:  1.2097E-01 -4.3981E-02 -4.8317E-01  1.4069E-01 -2.8465E-01  1.3419E-01  4.3569E-01 -1.8939E+00 -4.4264E-01 -3.0490E-01
             3.1043E-01
 GRADIENT:  -3.3608E+00  9.1075E+00 -8.8947E+00  2.6899E+01  9.5461E+00 -1.6978E+00  9.5503E-02  2.3901E-01 -8.6292E-01  1.3274E+00
            -2.1368E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1660.40482500799        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  1.0430E+00  7.9419E-01  6.6872E-01  1.0929E+00  7.2804E-01  1.0535E+00  1.5310E+00  1.4242E-01  5.7498E-01  7.5418E-01
             1.2439E+00
 PARAMETER:  1.4213E-01 -1.3043E-01 -3.0239E-01  1.8886E-01 -2.1740E-01  1.5208E-01  5.2589E-01 -1.8490E+00 -4.5342E-01 -1.8213E-01
             3.1826E-01
 GRADIENT:   1.0052E+01  2.1684E+00  3.3561E+00 -6.2954E+00 -4.9181E+00  1.8618E+00  2.1148E-01  1.7344E-01 -3.5290E-01  1.8304E-01
             1.0581E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1660.71839623004        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      636
 NPARAMETR:  1.0368E+00  6.6180E-01  6.9022E-01  1.1674E+00  6.9407E-01  1.0483E+00  1.7619E+00  7.0777E-02  5.5990E-01  7.6254E-01
             1.2387E+00
 PARAMETER:  1.3618E-01 -3.1279E-01 -2.7075E-01  2.5474E-01 -2.6518E-01  1.4712E-01  6.6640E-01 -2.5482E+00 -4.8000E-01 -1.7111E-01
             3.1409E-01
 GRADIENT:   9.4458E-01  1.0274E+00  1.7980E+00 -6.6859E-01 -3.2747E+00  1.4653E-01  1.7409E-01  4.5655E-02 -5.7639E-01  3.4510E-01
             2.3842E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1660.73016975265        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  1.0359E+00  6.4286E-01  7.0429E-01  1.1790E+00  6.9791E-01  1.0486E+00  1.8013E+00  5.1031E-02  5.6042E-01  7.7168E-01
             1.2388E+00
 PARAMETER:  1.3526E-01 -3.4183E-01 -2.5056E-01  2.6464E-01 -2.5967E-01  1.4749E-01  6.8848E-01 -2.8753E+00 -4.7907E-01 -1.5918E-01
             3.1415E-01
 GRADIENT:  -3.9169E-02 -1.9539E-01  6.3644E-01 -1.1288E+00 -5.7200E-01  3.7963E-01  2.4799E-01  2.2489E-02 -1.2554E-01 -1.3245E-01
             2.9973E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1660.74575437516        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      989
 NPARAMETR:  1.0362E+00  6.5967E-01  6.9640E-01  1.1693E+00  6.9857E-01  1.0478E+00  1.7626E+00  1.0000E-02  5.6281E-01  7.6788E-01
             1.2378E+00
 PARAMETER:  1.3558E-01 -3.1602E-01 -2.6183E-01  2.5641E-01 -2.5872E-01  1.4666E-01  6.6681E-01 -4.5227E+00 -4.7481E-01 -1.6412E-01
             3.1335E-01
 GRADIENT:  -1.0235E-02  3.3700E-02  6.7191E-03  1.6922E-01 -1.4741E-03 -1.7209E-02 -4.0571E-02  0.0000E+00 -2.6986E-02 -7.2870E-03
            -3.9391E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1660.74577136038        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1046
 NPARAMETR:  1.0362E+00  6.5915E-01  6.9579E-01  1.1694E+00  6.9801E-01  1.0478E+00  1.7640E+00  1.0000E-02  5.6284E-01  7.6722E-01
             1.2379E+00
 PARAMETER:  1.3559E-01 -3.1680E-01 -2.6271E-01  2.5648E-01 -2.5952E-01  1.4672E-01  6.6757E-01 -4.5158E+00 -4.7476E-01 -1.6498E-01
             3.1343E-01
 GRADIENT:   7.5037E-03 -2.8808E-02 -1.7644E-02 -7.4054E-02  3.4286E-02  6.2596E-03  1.1630E-02  0.0000E+00  7.5975E-03 -1.5340E-03
             1.4454E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1046
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4894E-04  2.1962E-02 -4.1795E-04 -2.7955E-02  2.2085E-03
 SE:             2.9815E-02  2.2693E-02  2.0882E-04  2.1413E-02  2.2067E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8799E-01  3.3315E-01  4.5344E-02  1.9172E-01  9.2028E-01

 ETASHRINKSD(%)  1.1665E-01  2.3977E+01  9.9300E+01  2.8264E+01  2.6072E+01
 ETASHRINKVR(%)  2.3316E-01  4.2205E+01  9.9995E+01  4.8540E+01  4.5346E+01
 EBVSHRINKSD(%)  5.6943E-01  2.4053E+01  9.9302E+01  2.7752E+01  2.4494E+01
 EBVSHRINKVR(%)  1.1356E+00  4.2321E+01  9.9995E+01  4.7803E+01  4.2989E+01
 RELATIVEINF(%)  9.8245E+01  7.3918E+00  2.8827E-04  6.7023E+00  3.1393E+00
 EPSSHRINKSD(%)  4.0803E+01
 EPSSHRINKVR(%)  6.4957E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1660.7457713603824     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -925.59494479664420     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.23
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1660.746       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  6.59E-01  6.96E-01  1.17E+00  6.98E-01  1.05E+00  1.76E+00  1.00E-02  5.63E-01  7.67E-01  1.24E+00
 


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
+        9.31E+02
 
 TH 2
+       -1.49E+01  4.88E+02
 
 TH 3
+        2.54E+01  2.21E+02  1.20E+03
 
 TH 4
+       -2.04E+01  4.73E+02 -5.89E+02  1.37E+03
 
 TH 5
+       -1.01E+01 -4.47E+02 -1.53E+03  5.94E+02  2.21E+03
 
 TH 6
+        1.90E-01 -3.06E+00  3.77E+00 -3.78E+00 -3.03E+00  1.78E+02
 
 TH 7
+        1.32E+00  3.73E+01 -9.11E+00 -2.08E+01  4.23E+00 -3.73E-01  2.63E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.94E+00 -2.15E+01 -5.16E+01 -3.66E+01  5.32E+01 -1.25E-02  1.89E+01  0.00E+00  1.96E+02
 
 TH10
+       -4.52E-01 -4.00E+00 -6.39E+01 -4.22E+01 -4.87E+01 -5.84E-01  6.14E+00  0.00E+00  1.80E+01  1.11E+02
 
 TH11
+       -6.70E+00 -1.10E+01 -3.64E+01 -9.91E+00  1.52E+01  2.31E+00  3.37E+00  0.00E+00  2.21E+01  2.69E+01  1.49E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.087
Stop Time:
Sat Sep 18 13:03:11 CDT 2021
