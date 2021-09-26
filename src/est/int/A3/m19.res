Fri Sep 24 22:07:18 CDT 2021
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
$DATA ../../../../data/int/A3/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -398.135976508415        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3926E+01  2.1410E+02  2.9457E+02 -5.0965E+01  2.3864E+02  3.8512E+01 -2.4806E+02 -2.4343E+02 -5.0211E+01 -2.5953E+02
            -6.1726E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2758.55533643801        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0416E+00  1.0007E+00  8.8884E-01  1.0585E+00  9.3836E-01  8.0163E-01  9.2080E-01  9.5772E-01  6.9586E-01  1.0450E+00
             2.9006E+00
 PARAMETER:  1.4074E-01  1.0073E-01 -1.7843E-02  1.5683E-01  3.6381E-02 -1.2110E-01  1.7490E-02  5.6801E-02 -2.6261E-01  1.4398E-01
             1.1649E+00
 GRADIENT:   1.8666E+01  1.0928E+01 -9.9855E+00  6.9725E+00  1.9315E+00 -3.7227E+01 -6.1834E-01  9.4306E+00 -9.9318E-01 -3.6562E-01
             9.0268E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2764.77013796820        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0352E+00  8.0762E-01  7.1462E-01  1.1731E+00  7.3637E-01  8.6154E-01  9.7503E-01  3.2723E-01  6.6034E-01  9.2291E-01
             2.8740E+00
 PARAMETER:  1.3458E-01 -1.1366E-01 -2.3600E-01  2.5962E-01 -2.0602E-01 -4.9037E-02  7.4711E-02 -1.0171E+00 -3.1500E-01  1.9774E-02
             1.1557E+00
 GRADIENT:  -1.9454E+00  4.2095E+01 -1.0720E+01  8.6433E+01  6.3366E+00 -1.0066E+01 -5.9109E+00  1.4677E+00 -2.2255E+01 -1.5256E+01
             8.3792E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2779.20363381135        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0193E+00  5.1085E-01  4.7389E-01  1.2816E+00  4.6804E-01  8.9756E-01  1.1323E+00  9.6610E-02  8.1945E-01  8.5918E-01
             2.6048E+00
 PARAMETER:  1.1913E-01 -5.7169E-01 -6.4678E-01  3.4813E-01 -6.5920E-01 -8.0796E-03  2.2423E-01 -2.2371E+00 -9.9126E-02 -5.1781E-02
             1.0574E+00
 GRADIENT:  -3.0886E+01  2.5962E+01  2.6558E+01  1.6616E+02 -5.6385E-01  2.9169E+00 -2.7224E+00  2.6525E-01 -1.3836E+01 -4.8993E-01
            -5.0333E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2791.69433035945        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0166E+00  2.7077E-01  2.2683E-01  1.2027E+00  2.4814E-01  9.6498E-01  1.3990E+00  1.0000E-02  1.0021E+00  7.6642E-01
             2.2462E+00
 PARAMETER:  1.1648E-01 -1.2065E+00 -1.3836E+00  2.8457E-01 -1.2938E+00  6.4349E-02  4.3576E-01 -5.5327E+00  1.0210E-01 -1.6602E-01
             9.0924E-01
 GRADIENT:  -2.2805E+01 -2.6461E+01  1.2822E+02  1.9938E+02 -1.1159E+02  2.5728E+01 -2.9276E+00  0.0000E+00 -5.2041E+01  4.1880E+00
            -2.3759E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2818.95297329915        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0352E+00  1.8305E-01  1.1715E-01  8.7977E-01  1.6601E-01  9.2288E-01  1.4733E+00  1.0000E-02  1.4482E+00  8.1821E-01
             2.2480E+00
 PARAMETER:  1.3458E-01 -1.5980E+00 -2.0443E+00 -2.8100E-02 -1.6957E+00  1.9743E-02  4.8747E-01 -1.1868E+01  4.7033E-01 -1.0064E-01
             9.1004E-01
 GRADIENT:   2.1713E+01 -3.3358E+01  4.8680E+01  3.4218E+00 -5.0643E+01  1.4617E+01  1.0783E-01  0.0000E+00 -2.2514E+01  5.4576E+00
            -1.1263E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2825.19778778128        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      475
 NPARAMETR:  1.0281E+00  1.9021E-01  1.1340E-01  8.8444E-01  1.6774E-01  8.9367E-01  1.4258E+00  1.0000E-02  1.5436E+00  7.9358E-01
             2.3497E+00
 PARAMETER:  1.2767E-01 -1.5596E+00 -2.0768E+00 -2.2796E-02 -1.6854E+00 -1.2419E-02  4.5477E-01 -1.2135E+01  5.3413E-01 -1.3120E-01
             9.5429E-01
 GRADIENT:  -6.9034E+00 -8.1088E+00 -1.9288E+01  1.0740E+01 -5.0004E+01  2.8794E+00 -2.4106E+00  0.0000E+00 -7.8701E+00  1.2520E+00
            -1.3169E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2830.05036384425        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  1.0313E+00  2.1552E-01  1.3844E-01  9.5112E-01  1.9111E-01  8.8915E-01  1.4207E+00  1.0000E-02  1.4208E+00  7.6205E-01
             2.3755E+00
 PARAMETER:  1.3085E-01 -1.4347E+00 -1.8773E+00  4.9883E-02 -1.5549E+00 -1.7487E-02  4.5112E-01 -1.0533E+01  4.5120E-01 -1.7175E-01
             9.6519E-01
 GRADIENT:  -3.1246E-02  2.8338E-02 -1.1735E-02  1.0331E-02 -7.6121E-03 -1.1036E-02  2.4398E-03  0.0000E+00  4.4118E-02  5.8006E-03
            -7.8258E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2830.05036384425        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  1.0313E+00  2.1552E-01  1.3844E-01  9.5112E-01  1.9111E-01  8.8915E-01  1.4207E+00  1.0000E-02  1.4208E+00  7.6205E-01
             2.3755E+00
 PARAMETER:  1.3085E-01 -1.4347E+00 -1.8773E+00  4.9883E-02 -1.5549E+00 -1.7487E-02  4.5112E-01 -1.0533E+01  4.5120E-01 -1.7175E-01
             9.6519E-01
 GRADIENT:  -3.1246E-02  2.8338E-02 -1.1735E-02  1.0331E-02 -7.6121E-03 -1.1036E-02  2.4398E-03  0.0000E+00  4.4118E-02  5.8006E-03
            -7.8258E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      673
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6098E-03  9.8786E-03  1.8297E-04 -4.9248E-03  7.4952E-03
 SE:             2.9322E-02  2.4965E-02  2.6625E-04  2.7029E-02  2.6800E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5622E-01  6.9233E-01  4.9196E-01  8.5542E-01  7.7973E-01

 ETASHRINKSD(%)  1.7662E+00  1.6363E+01  9.9108E+01  9.4488E+00  1.0215E+01
 ETASHRINKVR(%)  3.5012E+00  3.0049E+01  9.9992E+01  1.8005E+01  1.9387E+01
 EBVSHRINKSD(%)  1.8515E+00  1.5285E+01  9.9191E+01  6.2690E+00  1.0908E+01
 EBVSHRINKVR(%)  3.6688E+00  2.8233E+01  9.9993E+01  1.2145E+01  2.0626E+01
 RELATIVEINF(%)  9.6313E+01  2.4988E+01  7.8681E-04  4.4312E+01  1.0055E+01
 EPSSHRINKSD(%)  1.9577E+01
 EPSSHRINKVR(%)  3.5322E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2830.0503638442465     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1175.9610040758357     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2830.050       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.16E-01  1.38E-01  9.51E-01  1.91E-01  8.89E-01  1.42E+00  1.00E-02  1.42E+00  7.62E-01  2.38E+00
 


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
+        1.27E+03
 
 TH 2
+        5.51E+00  8.62E+03
 
 TH 3
+        4.73E+01 -3.58E+03  3.59E+04
 
 TH 4
+       -2.45E+00 -1.39E+02 -9.52E+02  5.37E+02
 
 TH 5
+       -1.83E+01 -4.41E+03 -2.86E+04 -3.94E+02  4.17E+04
 
 TH 6
+        4.71E+00 -9.42E+00  2.46E+01  6.73E-01 -2.01E+01  2.38E+02
 
 TH 7
+       -3.33E-01  4.43E+01  5.95E+01 -3.10E+00 -7.01E+01  7.01E-01  4.87E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.13E+01  4.61E+00  3.84E+02 -1.21E+01  4.81E+01 -6.47E-01  9.33E-01  0.00E+00  6.16E+01
 
 TH10
+       -7.39E+00  7.07E+00  1.16E+02  2.39E+01  4.02E+01  1.65E-01  4.05E+00  0.00E+00  8.57E-01  2.24E+02
 
 TH11
+       -2.07E+01 -1.19E+01 -2.18E+02 -9.75E-02  1.13E+02  2.84E+00  9.59E+00  0.00E+00  8.68E+00  1.19E+01  1.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.994
Stop Time:
Fri Sep 24 22:07:49 CDT 2021
