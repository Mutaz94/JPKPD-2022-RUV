Sat Sep 25 09:45:46 CDT 2021
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
$DATA ../../../../data/spa/S1/dat21.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1746.11389139992        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.8199E+01 -2.9151E+01 -1.7883E+01 -3.6851E+01 -4.2126E+01  4.7349E+00  2.5952E+00  1.8256E+01  1.7790E+01  1.6722E+01
             3.5402E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1755.00130632347        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0545E+00  1.0721E+00  1.2414E+00  9.5828E-01  1.1900E+00  9.6356E-01  1.0108E+00  8.4676E-01  8.7875E-01  9.8074E-01
             9.0087E-01
 PARAMETER:  1.5305E-01  1.6966E-01  3.1622E-01  5.7383E-02  2.7399E-01  6.2881E-02  1.1075E-01 -6.6337E-02 -2.9257E-02  8.0554E-02
            -4.3892E-03
 GRADIENT:   1.1720E+02 -3.1391E+01  1.2136E+01 -5.7949E+01  1.2204E+01 -5.9800E+00 -8.5407E+00  7.5366E-01 -1.1160E+01 -2.3335E+01
            -1.9899E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1758.77579342961        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0327E+00  8.3846E-01  1.0852E+00  1.1320E+00  1.0139E+00  9.7716E-01  1.4183E+00  3.8402E-01  7.1629E-01  9.7623E-01
             8.7666E-01
 PARAMETER:  1.3215E-01 -7.6192E-02  1.8173E-01  2.2398E-01  1.1384E-01  7.6897E-02  4.4946E-01 -8.5705E-01 -2.3367E-01  7.5946E-02
            -3.1641E-02
 GRADIENT:   5.9840E+01  7.9992E+00 -1.7133E+00  3.8381E+01  8.0786E+00  5.2478E-01 -1.0437E+00  4.6224E-01 -1.1183E+01 -7.7199E+00
            -2.6169E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1759.50039861670        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      251
 NPARAMETR:  1.0225E+00  8.9822E-01  1.0579E+00  1.0857E+00  1.0293E+00  9.7847E-01  1.2979E+00  3.5011E-01  7.9440E-01  9.9798E-01
             9.0448E-01
 PARAMETER:  1.2220E-01 -7.3443E-03  1.5631E-01  1.8223E-01  1.2884E-01  7.8232E-02  3.6076E-01 -9.4951E-01 -1.3017E-01  9.7982E-02
            -3.9389E-04
 GRADIENT:  -4.3716E+01 -8.0134E+00 -8.1473E+00 -1.0418E+01  8.4205E+00 -6.2078E+00 -1.7959E+00  7.8941E-01 -2.1181E+00 -1.6736E+00
            -9.3443E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1760.57421280712        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  1.0396E+00  7.4141E-01  1.1223E+00  1.1908E+00  9.9014E-01  9.8982E-01  1.5130E+00  2.0565E-01  7.5708E-01  1.0182E+00
             9.2484E-01
 PARAMETER:  1.3882E-01 -1.9920E-01  2.1534E-01  2.7462E-01  9.0086E-02  8.9766E-02  5.1408E-01 -1.4816E+00 -1.7828E-01  1.1800E-01
             2.1867E-02
 GRADIENT:  -1.2840E+00  8.4157E-01 -1.2349E+00  2.2715E+00  5.2427E-01  9.7503E-02 -2.7652E-01  9.7692E-02 -2.7303E-01  1.8505E-01
            -8.7259E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1760.61390606735        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      604
 NPARAMETR:  1.0392E+00  6.6720E-01  1.1396E+00  1.2340E+00  9.7020E-01  9.8850E-01  1.6405E+00  1.3705E-01  7.3894E-01  1.0193E+00
             9.2770E-01
 PARAMETER:  1.3846E-01 -3.0467E-01  2.3068E-01  3.1027E-01  6.9751E-02  8.8436E-02  5.9499E-01 -1.8874E+00 -2.0254E-01  1.1908E-01
             2.4951E-02
 GRADIENT:  -3.3329E-03 -1.5153E-01 -8.6548E-02 -1.5179E-01  5.5242E-02 -1.0768E-05 -1.3281E-03  2.0817E-02 -4.8884E-03  5.4579E-02
             1.9230E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1760.62445669927        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      785
 NPARAMETR:  1.0395E+00  6.8729E-01  1.1259E+00  1.2214E+00  9.7101E-01  9.8880E-01  1.6062E+00  3.8656E-02  7.4363E-01  1.0166E+00
             9.2786E-01
 PARAMETER:  1.3875E-01 -2.7500E-01  2.1855E-01  3.0003E-01  7.0579E-02  8.8741E-02  5.7385E-01 -3.1530E+00 -1.9621E-01  1.1651E-01
             2.5130E-02
 GRADIENT:  -8.7104E-02  3.3286E-03 -5.7957E-02 -3.0164E-02 -5.0500E-02 -1.6446E-02  1.4431E-02  1.6551E-03  1.5604E-02  3.7408E-02
             7.6169E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1760.62520201366        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      961
 NPARAMETR:  1.0395E+00  6.8663E-01  1.1263E+00  1.2219E+00  9.7103E-01  9.8883E-01  1.6071E+00  1.0000E-02  7.4347E-01  1.0169E+00
             9.2779E-01
 PARAMETER:  1.3876E-01 -2.7596E-01  2.1894E-01  3.0037E-01  7.0604E-02  8.8765E-02  5.7443E-01 -4.6798E+00 -1.9643E-01  1.1676E-01
             2.5050E-02
 GRADIENT:  -3.4856E-02  2.6603E-03  7.9506E-03  5.4110E-03 -9.7222E-03 -2.5567E-03  3.0981E-04  0.0000E+00  1.8081E-04 -3.1970E-03
            -5.5149E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1760.62521808656        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  1.0395E+00  6.8776E-01  1.1261E+00  1.2212E+00  9.7131E-01  9.8884E-01  1.6053E+00  1.0000E-02  7.4377E-01  1.0168E+00
             9.2774E-01
 PARAMETER:  1.3882E-01 -2.7436E-01  2.1870E-01  2.9982E-01  7.0922E-02  8.8797E-02  5.7324E-01 -4.6755E+00 -1.9609E-01  1.1681E-01
             2.5068E-02
 GRADIENT:   3.0115E-02 -4.0945E-03 -8.1843E-03  1.8300E-02  1.1408E-02  1.9352E-03 -1.8493E-03  0.0000E+00 -2.4235E-03  5.3482E-03
             7.3382E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1026
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.4909E-04  1.0549E-02 -3.6215E-04 -1.4696E-02 -1.6929E-02
 SE:             2.9847E-02  1.9690E-02  1.7049E-04  2.3639E-02  2.3878E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9601E-01  5.9212E-01  3.3657E-02  5.3416E-01  4.7833E-01

 ETASHRINKSD(%)  7.7626E-03  3.4036E+01  9.9429E+01  2.0807E+01  2.0006E+01
 ETASHRINKVR(%)  1.5525E-02  5.6487E+01  9.9997E+01  3.7284E+01  3.6009E+01
 EBVSHRINKSD(%)  3.8206E-01  3.5333E+01  9.9451E+01  1.9958E+01  1.6729E+01
 EBVSHRINKVR(%)  7.6265E-01  5.8181E+01  9.9997E+01  3.5933E+01  3.0659E+01
 RELATIVEINF(%)  9.8324E+01  2.6647E+00  4.0162E-04  4.4684E+00  8.5119E+00
 EPSSHRINKSD(%)  4.2298E+01
 EPSSHRINKVR(%)  6.6705E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1760.6252180865636     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1025.4743915228255     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1760.625       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  6.88E-01  1.13E+00  1.22E+00  9.71E-01  9.89E-01  1.61E+00  1.00E-02  7.44E-01  1.02E+00  9.28E-01
 


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
+        1.04E+03
 
 TH 2
+       -1.26E+01  4.03E+02
 
 TH 3
+        1.24E+01  1.03E+02  2.83E+02
 
 TH 4
+       -5.29E+00  4.34E+02 -1.46E+02  8.68E+02
 
 TH 5
+       -2.01E+00 -2.23E+02 -4.08E+02  1.67E+02  7.83E+02
 
 TH 6
+        5.30E+00 -2.54E+00  2.01E+00 -3.64E-01 -4.26E-01  2.01E+02
 
 TH 7
+        7.47E-01  2.86E+01  5.75E+00 -4.57E+00 -8.19E+00  6.41E-01  1.95E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.36E+00 -1.54E+01 -2.14E+01 -1.47E+01  2.53E+01 -2.02E+00  2.45E+01  0.00E+00  1.52E+02
 
 TH10
+        6.02E-01 -9.86E-01 -2.93E+01 -1.55E+01 -5.35E+01 -3.91E+00  1.71E+00  0.00E+00  4.69E+00  9.08E+01
 
 TH11
+       -7.41E+00 -1.80E+01 -4.33E+01 -5.07E+00  1.06E+01  3.22E+00  2.57E+00  0.00E+00  1.08E+01  2.42E+01  2.63E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.606
Stop Time:
Sat Sep 25 09:46:06 CDT 2021
