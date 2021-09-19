Sat Sep 18 07:25:32 CDT 2021
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
$DATA ../../../../data/int/D/dat71.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26152.6266872596        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2488E+01  3.1403E+02 -3.6625E+01 -1.0456E+02  2.2006E+02 -2.8835E+03 -1.2480E+03 -9.3335E+01 -2.3294E+03 -8.6795E+02
            -5.2630E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -911.088548038189        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  3.5119E+00  1.7659E+00  8.5955E-01  3.2963E+00  1.3163E+00  8.6798E+00  7.6037E+00  1.0383E+00  7.1480E+00  2.7672E+00
             9.0420E+00
 PARAMETER:  1.3562E+00  6.6868E-01 -5.1346E-02  1.2928E+00  3.7484E-01  2.2610E+00  2.1286E+00  1.3755E-01  2.0668E+00  1.1178E+00
             2.3019E+00
 GRADIENT:   3.6742E+01 -6.5606E+00 -4.2061E+01  5.7416E+01 -4.1488E+01  1.6464E+02  1.2009E+02  3.6908E+00  8.2141E+01  6.1768E+01
             2.4413E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1053.94193725431        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  2.0892E+00  5.5392E+00  1.7698E+01  6.0264E+00  3.8366E+00  3.5919E+00  6.0763E+00  8.3851E-01  3.7614E+01  3.9264E+00
             8.6699E+00
 PARAMETER:  8.3678E-01  1.8118E+00  2.9734E+00  1.8962E+00  1.4446E+00  1.3787E+00  1.9044E+00 -7.6123E-02  3.7274E+00  1.4677E+00
             2.2599E+00
 GRADIENT:   9.3599E+01 -7.1237E+00  2.8191E+01  1.8195E+01 -1.3613E+01  1.4175E+01  8.5604E+01 -8.2582E+00  1.0059E+02  1.3568E+02
             1.7185E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1210.58695281601        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.2630E+00  4.2863E+00  6.2459E+00  1.6948E+00  2.6212E+00  3.3907E+00  3.5038E+00  3.6355E-01  1.8981E+01  7.2081E-01
             8.8386E+00
 PARAMETER:  3.3346E-01  1.5554E+00  1.9319E+00  6.2757E-01  1.0636E+00  1.3210E+00  1.3538E+00 -9.1185E-01  3.0435E+00 -2.2737E-01
             2.2791E+00
 GRADIENT:   1.5054E+01  6.4374E-01 -2.2125E+01  2.2044E+01 -5.2748E+01  5.7363E+01  4.8983E+00 -5.7530E-01  2.2806E+01  2.5458E+00
             2.3356E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1268.48355788053        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  8.7383E-01  3.0828E+00  6.5126E+00  2.0014E-01  2.7242E+00  3.2381E+00  3.5857E+00  5.3736E-01  1.2796E+01  1.0233E-01
             7.4658E+00
 PARAMETER: -3.4867E-02  1.2258E+00  1.9737E+00 -1.5087E+00  1.1022E+00  1.2750E+00  1.3770E+00 -5.2109E-01  2.6492E+00 -2.1796E+00
             2.1103E+00
 GRADIENT:  -5.5339E+01  3.1897E+00 -2.1837E+01  1.2416E+00  4.4238E+00  1.0333E+01  2.2162E+01 -3.7114E-01  4.4262E+00 -5.0616E-03
            -3.6755E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1284.16376780107        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.1569E+00  3.2736E+00  9.7011E+00  1.1452E-01  2.7585E+00  3.0205E+00  3.2343E+00  3.9929E-01  1.4830E+01  4.6667E-02
             7.7081E+00
 PARAMETER:  2.4575E-01  1.2859E+00  2.3722E+00 -2.0670E+00  1.1147E+00  1.2054E+00  1.2738E+00 -8.1807E-01  2.7966E+00 -2.9647E+00
             2.1423E+00
 GRADIENT:  -2.1196E-01  7.4357E+00  6.7485E+00 -1.5494E+00  1.3331E+00 -2.3389E-01 -9.1416E+00  3.1681E-01  3.4371E+00  3.3024E-03
             6.2916E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1285.16729683762        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.1436E+00  3.0954E+00  8.7598E+00  1.2996E-01  2.7402E+00  3.0447E+00  3.3405E+00  2.8881E-01  1.4016E+01  4.9076E-02
             7.6776E+00
 PARAMETER:  2.3414E-01  1.2299E+00  2.2702E+00 -1.9406E+00  1.1080E+00  1.2134E+00  1.3061E+00 -1.1420E+00  2.7402E+00 -2.9144E+00
             2.1383E+00
 GRADIENT:  -2.5658E+00  2.1353E+00 -3.5721E+00 -2.4877E+00  9.0125E+00  2.4952E+00  1.0039E+01  6.8952E-01 -3.2037E+00  1.6955E-03
             2.7206E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1285.58414806237        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      578
 NPARAMETR:  1.1558E+00  3.0965E+00  8.7568E+00  1.3019E-01  2.7401E+00  3.0452E+00  3.3401E+00  1.0000E-02  1.3996E+01  4.7850E-02
             7.6715E+00
 PARAMETER:  2.4478E-01  1.2303E+00  2.2698E+00 -1.9388E+00  1.1080E+00  1.2136E+00  1.3060E+00 -4.7667E+00  2.7388E+00 -2.9397E+00
             2.1375E+00
 GRADIENT:  -2.6952E-02  2.2713E+00 -1.1174E+00 -2.4880E+00  8.8933E+00  2.5042E+00  9.9192E+00  0.0000E+00 -3.1821E+00  1.4349E-03
             6.2786E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1285.60892133505        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:      683
 NPARAMETR:  1.1630E+00  3.0958E+00  8.7650E+00  1.3023E-01  2.7351E+00  3.0409E+00  3.3381E+00  1.0000E-02  1.4008E+01  1.6074E-02
             7.6673E+00
 PARAMETER:  2.5101E-01  1.2300E+00  2.2708E+00 -1.9384E+00  1.1062E+00  1.2122E+00  1.3054E+00 -5.0457E+00  2.7396E+00 -4.0306E+00
             2.1370E+00
 GRADIENT:  -4.4483E-01 -5.4620E+00 -1.4280E+00 -3.2733E+00  7.5879E+00 -4.5432E+00  1.4657E+00  0.0000E+00 -9.8977E+00  7.4208E-05
            -4.7952E+00

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1285.60892133505        NO. OF FUNC. EVALS.:  35
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  1.1633E+00  3.0959E+00  8.7643E+00  1.3022E-01  2.7352E+00  3.0410E+00  3.3382E+00  1.0000E-02  1.4006E+01  1.6074E-02
             7.6668E+00
 PARAMETER:  2.5101E-01  1.2300E+00  2.2708E+00 -1.9384E+00  1.1062E+00  1.2122E+00  1.3054E+00 -5.0457E+00  2.7396E+00 -4.0306E+00
             2.1370E+00
 GRADIENT:  -4.3075E-01 -7.3255E+01  7.4364E+01  8.1572E+01 -1.4873E+02 -6.7451E+01 -1.4280E+02  0.0000E+00  2.7456E+02 -6.8070E-06
             7.4752E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      718
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.4014E-03  1.1494E-03 -1.0540E-05  1.8525E-02 -2.4382E-04
 SE:             2.9769E-02  2.8027E-02  3.6231E-05  1.2961E-02  2.4269E-04
 N:                     100         100         100         100         100

 P VAL.:         8.0365E-01  9.6729E-01  7.7112E-01  1.5294E-01  3.1506E-01

 ETASHRINKSD(%)  2.7096E-01  6.1064E+00  9.9879E+01  5.6577E+01  9.9187E+01
 ETASHRINKVR(%)  5.4118E-01  1.1840E+01  1.0000E+02  8.1145E+01  9.9993E+01
 EBVSHRINKSD(%)  2.2695E+00  5.2381E+00  9.9730E+01  6.2327E+01  9.9182E+01
 EBVSHRINKVR(%)  4.4874E+00  1.0202E+01  9.9999E+01  8.5807E+01  9.9993E+01
 RELATIVEINF(%)  9.5467E+01  6.2131E+01  7.2941E-04  1.0020E+01  6.4853E-03
 EPSSHRINKSD(%)  8.7500E+00
 EPSSHRINKVR(%)  1.6734E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1285.6089213350506     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       368.48043843336018     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    19.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1285.609       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.16E+00  3.10E+00  8.76E+00  1.30E-01  2.74E+00  3.04E+00  3.34E+00  1.00E-02  1.40E+01  1.61E-02  7.67E+00
 


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
+        4.33E+05
 
 TH 2
+        5.12E+01  3.05E+03
 
 TH 3
+       -8.20E+00 -6.47E+01  1.17E+02
 
 TH 4
+       -8.31E+02 -5.10E+03  7.38E+02  2.94E+06
 
 TH 5
+        5.25E+01  4.13E+02 -1.42E+02 -4.76E+03  5.09E+03
 
 TH 6
+        1.49E+02  4.86E+01 -7.85E+00 -7.97E+02  5.11E+01  2.89E+03
 
 TH 7
+        2.55E+01  8.94E+03 -5.45E+01 -1.34E+05  1.14E+04  2.40E+01  1.04E+04
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.82E+00 -2.25E+01  6.38E+00  2.15E+02 -4.85E+01 -1.62E+00 -9.35E+02  0.00E+00  1.37E+02
 
 TH10
+       -6.27E+00 -6.25E+00  1.10E+00  1.14E+02 -4.33E+00 -1.50E+01 -2.97E+00  0.00E+00  3.06E-01  7.90E+00
 
 TH11
+       -1.39E+01 -9.06E+01  2.03E+01  1.00E+03 -1.32E+02 -9.13E+00 -5.63E+01  0.00E+00  6.55E+00  6.00E-01  1.86E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       42.867
Stop Time:
Sat Sep 18 07:26:17 CDT 2021
