Sat Oct 23 23:06:56 CDT 2021
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
$DATA ../../../../data/SD3/S1/dat58.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m58.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2061.06731389546        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1979E+02  5.9209E+00  4.0179E+00  4.7631E+01  2.4445E+01  5.4411E+01 -3.5346E+00 -1.6359E+00  2.3234E+01 -4.6964E+00
            -7.5474E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2067.29918399406        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.3341E-01  1.0769E+00  9.0644E-01  9.7487E-01  9.9445E-01  8.9139E-01  1.1052E+00  1.0193E+00  8.1485E-01  1.0349E+00
             1.0178E+00
 PARAMETER:  3.1089E-02  1.7409E-01  1.7742E-03  7.4550E-02  9.4434E-02 -1.4978E-02  2.0004E-01  1.1911E-01 -1.0475E-01  1.3433E-01
             1.1760E-01
 GRADIENT:  -2.4532E+01  3.1751E+00 -1.1890E+01  2.5713E+01  1.5266E+01 -2.9263E+01 -1.2615E+01  3.5732E+00 -5.5793E+00  6.4253E+00
             7.6594E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2068.33076047866        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.4079E-01  8.8980E-01  8.6823E-01  1.0894E+00  8.7142E-01  9.1065E-01  1.4625E+00  9.2558E-01  6.6674E-01  8.7270E-01
             1.0200E+00
 PARAMETER:  3.8962E-02 -1.6764E-02 -4.1304E-02  1.8560E-01 -3.7635E-02  6.4039E-03  4.8017E-01  2.2661E-02 -3.0535E-01 -3.6167E-02
             1.1978E-01
 GRADIENT:   6.1103E-03  2.0663E+01 -1.6003E+01  4.4589E+01  8.7453E+00 -1.9031E+01 -8.9246E-01  2.2477E+00 -1.3576E+01 -1.1127E+00
             5.1788E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2070.40255962248        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.3968E-01  8.3974E-01  9.4436E-01  1.1083E+00  8.8925E-01  9.5264E-01  1.4536E+00  8.6875E-01  7.5651E-01  9.0888E-01
             1.0099E+00
 PARAMETER:  3.7785E-02 -7.4667E-02  4.2753E-02  2.0282E-01 -1.7379E-02  5.1482E-02  4.7406E-01 -4.0701E-02 -1.7903E-01  4.4568E-03
             1.0983E-01
 GRADIENT:   6.2226E-02  3.6799E+00  3.1731E+00  4.4143E+00 -2.1068E+00 -2.6087E-01 -8.0771E-01 -8.2861E-01  1.9930E+00 -1.6864E+00
            -9.8644E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2070.67795377652        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.3687E-01  6.3889E-01  1.0804E+00  1.2316E+00  8.8205E-01  9.5055E-01  1.8384E+00  1.0100E+00  6.7222E-01  9.2586E-01
             1.0121E+00
 PARAMETER:  3.4793E-02 -3.4802E-01  1.7731E-01  3.0830E-01 -2.5503E-02  4.9284E-02  7.0887E-01  1.0996E-01 -2.9717E-01  2.2965E-02
             1.1200E-01
 GRADIENT:  -1.2574E-01  1.3356E+00 -6.2837E-01  3.1392E+00  3.5899E-02  3.9803E-02 -2.6647E-01 -7.4419E-03 -3.3808E-01  9.3211E-02
             1.5686E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2070.68736111634        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      885             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3665E-01  6.0016E-01  1.1114E+00  1.2531E+00  8.8324E-01  9.5018E-01  1.9267E+00  1.0402E+00  6.6400E-01  9.2850E-01
             1.0120E+00
 PARAMETER:  3.4553E-02 -4.1056E-01  2.0563E-01  3.2564E-01 -2.4164E-02  4.8899E-02  7.5583E-01  1.3937E-01 -3.0947E-01  2.5818E-02
             1.1195E-01
 GRADIENT:   3.6977E+02  5.0973E+01  5.0546E+00  3.6092E+02  7.9052E+00  4.7542E+01  3.3173E+01  4.4478E-01  1.1426E+01  7.6469E-01
             1.2419E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2070.68848222378        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1067
 NPARAMETR:  9.3684E-01  6.0085E-01  1.1112E+00  1.2531E+00  8.8328E-01  9.5035E-01  1.9314E+00  1.0399E+00  6.6296E-01  9.2886E-01
             1.0119E+00
 PARAMETER:  3.4757E-02 -4.0940E-01  2.0543E-01  3.2564E-01 -2.4109E-02  4.9078E-02  7.5823E-01  1.3910E-01 -3.1104E-01  2.6204E-02
             1.1183E-01
 GRADIENT:   1.2985E+00  2.0902E-02 -9.5150E-02 -2.1728E+00  6.0864E-02  1.7919E-01  2.9015E-01  4.2680E-02  1.2108E-01  9.6684E-03
             1.6904E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -2070.68848222378        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  9.3684E-01  6.0085E-01  1.1112E+00  1.2531E+00  8.8328E-01  9.5035E-01  1.9314E+00  1.0399E+00  6.6296E-01  9.2886E-01
             1.0119E+00
 PARAMETER:  3.4757E-02 -4.0940E-01  2.0543E-01  3.2564E-01 -2.4109E-02  4.9078E-02  7.5823E-01  1.3910E-01 -3.1104E-01  2.6204E-02
             1.1183E-01
 GRADIENT:  -1.1232E-03 -5.2296E-02 -8.6948E-02  2.4152E-01  7.1161E-02  4.9545E-04 -3.4490E-02  2.2034E-02  7.9663E-02  8.6435E-03
             1.7359E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1091
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0664E-03  2.0554E-02 -3.5878E-02 -2.3707E-02 -2.2078E-02
 SE:             2.9876E-02  2.0969E-02  1.6762E-02  2.2349E-02  2.0701E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7153E-01  3.2697E-01  3.2317E-02  2.8880E-01  2.8620E-01

 ETASHRINKSD(%)  1.0000E-10  2.9753E+01  4.3846E+01  2.5127E+01  3.0649E+01
 ETASHRINKVR(%)  1.0000E-10  5.0653E+01  6.8467E+01  4.3940E+01  5.1904E+01
 EBVSHRINKSD(%)  3.8931E-01  3.1106E+01  4.7556E+01  2.3847E+01  2.7851E+01
 EBVSHRINKVR(%)  7.7710E-01  5.2537E+01  7.2496E+01  4.2007E+01  4.7946E+01
 RELATIVEINF(%)  9.8522E+01  5.2778E+00  4.7924E+00  6.3991E+00  1.1248E+01
 EPSSHRINKSD(%)  3.4227E+01
 EPSSHRINKVR(%)  5.6739E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2070.6884822237762     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1151.7499490191035     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2070.688       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.37E-01  6.01E-01  1.11E+00  1.25E+00  8.83E-01  9.50E-01  1.93E+00  1.04E+00  6.63E-01  9.29E-01  1.01E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       88.276
Stop Time:
Sat Oct 23 23:07:11 CDT 2021
