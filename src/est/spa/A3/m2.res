Wed Sep 29 13:14:34 CDT 2021
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
$DATA ../../../../data/spa/A3/dat2.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   224.216985294171        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8369E+02  1.1290E+02  9.2354E+01  4.1888E+01  1.6672E+02  3.2750E+01 -6.9944E+01 -1.9809E+01 -1.2101E+02 -1.6990E+02
            -3.3110E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1189.14642649254        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0305E+00  9.5629E-01  9.1992E-01  1.1468E+00  9.3799E-01  8.5955E-01  9.9544E-01  9.6979E-01  1.0595E+00  1.0628E+00
             5.4292E+00
 PARAMETER:  1.3008E-01  5.5309E-02  1.6529E-02  2.3702E-01  3.5983E-02 -5.1347E-02  9.5432E-02  6.9328E-02  1.5780E-01  1.6093E-01
             1.7918E+00
 GRADIENT:   1.4828E+01 -1.1643E+01 -1.6672E+01  2.6623E+00 -6.6714E+00 -5.6064E+00  1.3233E+01  6.7857E+00  2.9663E+01  2.5762E+01
             1.9891E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1215.75091079073        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.0110E+00  6.4508E-01  3.6494E-01  1.2461E+00  4.1360E-01  9.0405E-01  6.2819E-01  9.0775E-02  1.0466E+00  4.0618E-01
             4.8068E+00
 PARAMETER:  1.1096E-01 -3.3838E-01 -9.0803E-01  3.2000E-01 -7.8285E-01 -8.7023E-04 -3.6491E-01 -2.2994E+00  1.4550E-01 -8.0095E-01
             1.6700E+00
 GRADIENT:  -3.8230E+01  5.8706E+01  1.7595E+01  9.9805E+01 -5.5701E+01 -9.9414E+00  2.4810E+00  1.6536E-01  1.2314E+01  8.5768E+00
             1.4387E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1236.09537532276        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.7359E-01  5.5332E-01  2.0447E-01  1.0733E+00  2.9021E-01  9.5068E-01  2.6252E-01  1.0000E-02  1.5049E+00  2.3215E-01
             3.6460E+00
 PARAMETER:  7.3234E-02 -4.9181E-01 -1.4873E+00  1.7077E-01 -1.1371E+00  4.9424E-02 -1.2374E+00 -4.6961E+00  5.0872E-01 -1.3604E+00
             1.3936E+00
 GRADIENT:  -2.9714E+01  1.2738E+02  9.3378E+01  4.4646E+01 -1.5185E+02 -4.0728E+00 -2.7390E+00  0.0000E+00  3.7074E+01 -6.3221E+00
            -9.7323E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1246.72964364594        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      391
 NPARAMETR:  1.0002E+00  4.7067E-01  2.4349E-01  1.1269E+00  3.0588E-01  9.4312E-01  5.0917E-01  1.0000E-02  1.1869E+00  1.5621E-01
             3.6201E+00
 PARAMETER:  1.0021E-01 -6.5360E-01 -1.3127E+00  2.1943E-01 -1.0846E+00  4.1434E-02 -5.7498E-01 -4.8342E+00  2.7134E-01 -1.7566E+00
             1.3865E+00
 GRADIENT:   2.2254E+00  2.2515E+01  1.0422E+01  2.6591E+01 -1.8586E+01 -5.4405E+00 -2.6017E+00  0.0000E+00 -7.2971E-01 -1.4438E+00
            -2.7158E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1249.53126590786        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      566
 NPARAMETR:  1.0031E+00  3.5675E-01  2.3646E-01  1.1289E+00  2.7601E-01  9.5649E-01  9.3052E-01  1.0000E-02  1.1601E+00  8.0380E-02
             3.7026E+00
 PARAMETER:  1.0308E-01 -9.3073E-01 -1.3420E+00  2.2126E-01 -1.1873E+00  5.5520E-02  2.7986E-02 -5.9470E+00  2.4847E-01 -2.4210E+00
             1.4090E+00
 GRADIENT:   3.6606E+00  1.8681E+00 -2.3568E+00 -4.3496E+00  2.5032E+00 -3.8408E-01  7.0015E-01  0.0000E+00  1.0855E-01 -2.5719E-01
             9.4723E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1250.53041296298        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      743
 NPARAMETR:  9.9742E-01  2.4714E-01  2.8344E-01  1.2255E+00  2.8828E-01  9.4546E-01  1.0404E+00  1.0000E-02  1.0910E+00  3.8475E-02
             3.7660E+00
 PARAMETER:  9.7417E-02 -1.2978E+00 -1.1608E+00  3.0332E-01 -1.1438E+00  4.3913E-02  1.3958E-01 -6.6912E+00  1.8711E-01 -3.1577E+00
             1.4260E+00
 GRADIENT:  -5.4126E+00  3.0040E+00  2.4942E+00  5.2176E+00 -4.9815E+00 -5.8531E-01 -7.9000E-01  0.0000E+00 -4.3201E-01 -6.1128E-02
             6.8849E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1251.85918446379        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      918
 NPARAMETR:  9.9135E-01  1.2242E-01  2.9765E-01  1.2767E+00  2.8292E-01  9.3875E-01  2.8334E+00  1.0000E-02  1.0457E+00  1.0000E-02
             3.7599E+00
 PARAMETER:  9.1315E-02 -2.0003E+00 -1.1119E+00  3.4425E-01 -1.1626E+00  3.6789E-02  1.1415E+00 -9.7918E+00  1.4472E-01 -5.0375E+00
             1.4244E+00
 GRADIENT:  -1.3664E+00  1.4831E+00  1.1251E+01  1.7695E+00 -1.6084E+01 -2.3863E-01 -1.6699E-01  0.0000E+00  3.0072E-03  0.0000E+00
            -3.0031E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1252.18683605462        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1093
 NPARAMETR:  9.8691E-01  6.5671E-02  2.9491E-01  1.2899E+00  2.7802E-01  9.3797E-01  4.8612E+00  1.0000E-02  1.0333E+00  1.0000E-02
             3.7538E+00
 PARAMETER:  8.6827E-02 -2.6231E+00 -1.1211E+00  3.5456E-01 -1.1801E+00  3.5960E-02  1.6813E+00 -1.2799E+01  1.3279E-01 -6.7315E+00
             1.4228E+00
 GRADIENT:  -1.6659E-01  1.2516E+00  8.8212E-01 -1.9866E+00 -1.4980E+00 -2.5657E-01  1.6073E+00  0.0000E+00 -5.8905E-01  0.0000E+00
            -1.1290E+00

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1252.19053204296        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1186
 NPARAMETR:  9.8697E-01  6.5153E-02  2.9481E-01  1.2912E+00  2.7810E-01  9.3838E-01  4.8368E+00  1.0000E-02  1.0347E+00  1.0000E-02
             3.7663E+00
 PARAMETER:  8.6883E-02 -2.6310E+00 -1.1214E+00  3.5558E-01 -1.1798E+00  3.6400E-02  1.6763E+00 -1.2799E+01  1.3411E-01 -6.7315E+00
             1.4261E+00
 GRADIENT:  -6.9897E-01  2.1800E-01 -6.7962E-01 -8.5629E-02  7.4497E-01  1.7144E-01  1.5511E-01  0.0000E+00  2.0939E-01  0.0000E+00
             2.5985E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1186
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.3169E-04  6.0589E-05  1.5862E-04 -1.6329E-02  1.2002E-04
 SE:             2.8245E-02  4.5141E-03  2.4711E-04  2.5668E-02  3.7161E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8498E-01  9.8929E-01  5.2094E-01  5.2468E-01  7.4672E-01

 ETASHRINKSD(%)  5.3768E+00  8.4877E+01  9.9172E+01  1.4008E+01  9.8755E+01
 ETASHRINKVR(%)  1.0465E+01  9.7713E+01  9.9993E+01  2.6054E+01  9.9985E+01
 EBVSHRINKSD(%)  4.9292E+00  8.5808E+01  9.9154E+01  1.3169E+01  9.8799E+01
 EBVSHRINKVR(%)  9.6154E+00  9.7986E+01  9.9993E+01  2.4603E+01  9.9986E+01
 RELATIVEINF(%)  6.9195E+01  2.1017E-01  2.5737E-04  1.4643E+01  4.3390E-04
 EPSSHRINKSD(%)  2.1737E+01
 EPSSHRINKVR(%)  3.8748E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1252.1905320429562     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -517.03970547921801     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1252.191       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  6.52E-02  2.95E-01  1.29E+00  2.78E-01  9.38E-01  4.84E+00  1.00E-02  1.03E+00  1.00E-02  3.77E+00
 


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
+        1.17E+03
 
 TH 2
+       -2.67E+02  8.71E+04
 
 TH 3
+       -1.59E+02  3.60E+03  8.49E+03
 
 TH 4
+       -6.27E+01  5.17E+02 -2.25E+02  4.71E+02
 
 TH 5
+        4.35E+02 -4.86E+03 -1.23E+04 -5.39E+02  2.09E+04
 
 TH 6
+       -3.20E+00  1.04E+02  4.67E+01 -1.70E+01  1.49E+01  1.83E+02
 
 TH 7
+       -2.59E+00  1.84E+03  5.89E+01 -2.19E+02 -6.72E+01  2.49E+00  3.92E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.31E+01  8.63E+01  3.64E+01 -1.06E+01  1.45E+02  3.68E+00  2.53E+00  0.00E+00  1.02E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.12E+01 -8.52E+02 -2.87E+01 -6.19E+00  5.37E+01  3.68E+00 -1.77E+01  0.00E+00  7.30E+00  0.00E+00  5.70E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.435
Stop Time:
Wed Sep 29 13:14:58 CDT 2021
