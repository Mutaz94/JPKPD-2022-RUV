Sat Sep 25 01:35:21 CDT 2021
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
$DATA ../../../../data/int/SL2/dat77.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      999
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

 TOT. NO. OF OBS RECS:      899
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1343.28408155599        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.7257E+01  3.7787E+01  1.9264E+02  7.7104E+01  1.7317E+01 -1.1063E+01 -1.0202E+02 -2.0878E+02 -3.3393E+01 -3.6635E+01
            -4.6766E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2856.92893027331        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0361E+00  1.1958E+00  9.1171E-01  8.5918E-01  1.1519E+00  1.0266E+00  1.0580E+00  9.6260E-01  7.8995E-01  1.0907E+00
             2.0228E+00
 PARAMETER:  1.3550E-01  2.7882E-01  7.5683E-03 -5.1777E-02  2.4141E-01  1.2625E-01  1.5634E-01  6.1886E-02 -1.3578E-01  1.8678E-01
             8.0448E-01
 GRADIENT:   3.9415E+00 -4.1763E+01  2.5216E+00 -3.9726E+01 -7.8781E-01  2.3258E+00  3.7516E+00 -4.0042E+00 -8.1459E+00 -1.7140E+01
            -4.5033E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2880.73296810652        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0459E+00  1.4549E+00  7.9732E-01  7.3548E-01  1.3013E+00  9.5007E-01  8.0787E-01  2.5177E-01  8.9979E-01  1.5125E+00
             2.1518E+00
 PARAMETER:  1.4491E-01  4.7493E-01 -1.2650E-01 -2.0724E-01  3.6335E-01  4.8780E-02 -1.1335E-01 -1.2792E+00 -5.5920E-03  5.1373E-01
             8.6628E-01
 GRADIENT:   1.9175E+01 -4.7412E+00 -5.8803E+00  2.3741E+01 -1.3268E+01 -2.9241E+01 -5.2742E+00 -1.0552E+00  1.5508E+00  2.5646E+01
            -2.8449E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2900.74450554941        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0422E+00  1.5778E+00  7.8103E-01  6.6250E-01  1.4192E+00  1.0169E+00  8.0548E-01  7.2624E-02  9.0285E-01  1.4108E+00
             2.4090E+00
 PARAMETER:  1.4130E-01  5.5600E-01 -1.4715E-01 -3.1173E-01  4.5010E-01  1.1676E-01 -1.1632E-01 -2.5225E+00 -2.1939E-03  4.4417E-01
             9.7920E-01
 GRADIENT:   1.4252E+00  3.8979E+00 -6.3412E-01  1.2080E+01  2.2923E-01 -1.5682E-01 -2.1875E+00 -4.0957E-02  1.1837E+00  2.4011E+00
             1.4588E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2903.02565666793        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0408E+00  1.7400E+00  6.8104E-01  5.4711E-01  1.5504E+00  1.0163E+00  7.7002E-01  4.5136E-02  9.9052E-01  1.4734E+00
             2.3775E+00
 PARAMETER:  1.4003E-01  6.5390E-01 -2.8414E-01 -5.0311E-01  5.3853E-01  1.1619E-01 -1.6134E-01 -2.9981E+00  9.0478E-02  4.8756E-01
             9.6606E-01
 GRADIENT:   2.1831E-01 -8.4949E+00 -7.0873E-01 -2.7030E+00  9.6222E-01 -4.9052E-01  1.2308E+00 -1.3318E-02 -4.7934E-01 -1.4168E+00
            -4.9933E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2904.49020862024        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0410E+00  1.8661E+00  6.1990E-01  4.8069E-01  1.6433E+00  1.0194E+00  7.2303E-01  2.7579E-02  1.1164E+00  1.5405E+00
             2.3788E+00
 PARAMETER:  1.4016E-01  7.2386E-01 -3.7819E-01 -6.3254E-01  5.9673E-01  1.1921E-01 -2.2431E-01 -3.4907E+00  2.1014E-01  5.3214E-01
             9.6661E-01
 GRADIENT:  -1.3389E-01  2.9369E+01  4.8703E-01  1.1365E+01  1.5789E+00  5.0509E-01 -9.4140E-01 -4.1752E-03 -4.7606E-01  8.2769E-01
             1.8172E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2907.53699475059        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  1.0478E+00  2.1714E+00  3.6866E-01  2.7739E-01  1.9048E+00  1.0229E+00  6.5453E-01  1.0000E-02  1.6572E+00  1.6823E+00
             2.3609E+00
 PARAMETER:  1.4670E-01  8.7536E-01 -8.9788E-01 -1.1823E+00  7.4440E-01  1.2268E-01 -3.2383E-01 -6.2621E+00  6.0514E-01  6.2015E-01
             9.5904E-01
 GRADIENT:   2.0448E+00  3.7240E-01 -2.3950E+00  5.8428E+00  3.3840E+00  5.8998E-01 -3.4989E-01  0.0000E+00  1.7557E-01  3.5396E-02
            -1.4413E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2908.03725650687        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  1.0467E+00  2.2869E+00  2.8606E-01  1.9711E-01  2.0082E+00  1.0212E+00  6.3537E-01  1.0000E-02  2.0365E+00  1.7371E+00
             2.3597E+00
 PARAMETER:  1.4561E-01  9.2721E-01 -1.1515E+00 -1.5240E+00  7.9723E-01  1.2093E-01 -3.5355E-01 -7.7873E+00  8.1125E-01  6.5224E-01
             9.5854E-01
 GRADIENT:   2.3179E-01 -1.7861E+00 -5.4874E-01  1.0582E-01  4.2047E-01  6.6664E-02  3.0024E-02  0.0000E+00 -6.8161E-02 -1.3269E-01
            -2.6919E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2908.04352388372        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      856
 NPARAMETR:  1.0466E+00  2.2833E+00  2.9466E-01  2.0001E-01  2.0049E+00  1.0210E+00  6.3581E-01  1.0000E-02  2.0295E+00  1.7358E+00
             2.3601E+00
 PARAMETER:  1.4550E-01  9.2560E-01 -1.1219E+00 -1.5094E+00  7.9557E-01  1.2077E-01 -3.5286E-01 -7.6515E+00  8.0777E-01  6.5149E-01
             9.5869E-01
 GRADIENT:  -1.1986E-02  3.3361E-02  7.0375E-04  4.9256E-03  1.3741E-02 -2.3738E-03  1.5273E-02  0.0000E+00  1.0966E-03  6.8955E-03
             2.2515E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      856
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.9101E-04 -1.9827E-02 -4.3648E-05  2.6325E-02 -1.6781E-02
 SE:             2.9588E-02  2.6277E-02  3.0614E-05  1.5673E-02  2.6711E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7598E-01  4.5052E-01  1.5394E-01  9.3020E-02  5.2984E-01

 ETASHRINKSD(%)  8.7767E-01  1.1970E+01  9.9897E+01  4.7494E+01  1.0514E+01
 ETASHRINKVR(%)  1.7476E+00  2.2508E+01  1.0000E+02  7.2431E+01  1.9923E+01
 EBVSHRINKSD(%)  1.1364E+00  1.2022E+01  9.9896E+01  5.3680E+01  7.5597E+00
 EBVSHRINKVR(%)  2.2599E+00  2.2599E+01  1.0000E+02  7.8545E+01  1.4548E+01
 RELATIVEINF(%)  9.7686E+01  1.2524E+01  4.3230E-05  2.8688E+00  4.6293E+01
 EPSSHRINKSD(%)  1.6604E+01
 EPSSHRINKVR(%)  3.0450E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          899
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1652.2514827020016     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2908.0435238837240     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1255.7920411817224     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.68
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2908.044       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.28E+00  2.95E-01  2.00E-01  2.00E+00  1.02E+00  6.36E-01  1.00E-02  2.03E+00  1.74E+00  2.36E+00
 


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
+        9.52E+02
 
 TH 2
+       -1.21E+01  4.17E+02
 
 TH 3
+        1.72E+00  3.98E+01  1.98E+02
 
 TH 4
+       -2.83E+01  5.07E+02 -2.79E+02  1.66E+03
 
 TH 5
+       -2.03E+00 -2.39E+01 -2.68E+01  1.01E+02  7.86E+01
 
 TH 6
+        3.82E+00 -3.15E+00  6.59E-01 -8.24E+00 -4.48E-01  1.82E+02
 
 TH 7
+       -7.95E-01 -2.12E+00 -3.36E+00 -2.69E+01 -5.29E-01 -2.80E+00  3.05E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -8.50E-02 -4.34E+00 -1.21E+01  5.67E+01 -3.94E-02 -1.24E-01  7.07E+00  0.00E+00  8.67E+00
 
 TH10
+        3.04E-01 -3.84E+00 -1.22E-01  2.20E+01 -6.22E+00 -3.34E-01  5.19E+00  0.00E+00  8.90E-01  4.37E+01
 
 TH11
+       -1.18E+01 -1.62E+01 -3.78E+00 -8.73E+00  1.13E+00  2.56E+00  1.11E+01  0.00E+00  1.78E+00  4.27E+00  2.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.150
Stop Time:
Sat Sep 25 01:35:52 CDT 2021
