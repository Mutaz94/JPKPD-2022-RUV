Sat Oct 23 15:21:23 CDT 2021
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
$DATA ../../../../data/SD1/SL3/dat9.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      980
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      880
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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -367.859422084206        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3466E+02 -3.2033E+00  1.0884E+02  1.8601E+02  2.1839E+02  1.9363E+01 -9.7212E+01 -2.9229E+02 -1.4003E+02 -3.1683E+01
            -6.1361E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2540.80114177062        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.7036E-01  1.2953E+00  9.7042E-01  8.7892E-01  1.0048E+00  9.0035E-01  1.0912E+00  1.0863E+00  1.0885E+00  1.0562E+00
             3.4548E+00
 PARAMETER:  6.9908E-02  3.5873E-01  6.9978E-02 -2.9061E-02  1.0482E-01 -4.9729E-03  1.8727E-01  1.8280E-01  1.8482E-01  1.5469E-01
             1.3398E+00
 GRADIENT:   5.8545E+01  5.3818E+01  8.8414E+00 -9.2589E+00 -6.3190E+01 -1.3464E+01  2.2243E+01  2.9564E+00  7.5281E+00 -5.3676E+00
             4.4771E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2575.53127995579        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  9.4143E-01  1.8236E+00  1.8870E+00  6.6705E-01  1.4620E+00  8.8280E-01  7.9236E-01  4.4237E+00  2.1219E+00  1.5763E+00
             3.1099E+00
 PARAMETER:  3.9646E-02  7.0082E-01  7.3498E-01 -3.0489E-01  4.7983E-01 -2.4653E-02 -1.3274E-01  1.5870E+00  8.5231E-01  5.5508E-01
             1.2346E+00
 GRADIENT:  -4.1385E-01  1.9097E+02 -5.3065E+00  8.3145E+01 -4.1887E+01 -2.6313E+01  2.9819E+01  2.5294E+01  5.6409E+01  1.2244E+01
             3.4185E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2626.29285022662        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.2118E-01  1.4507E+00  2.2913E+00  7.7197E-01  1.4282E+00  9.1929E-01  8.0382E-01  3.3465E+00  1.0829E+00  1.3797E+00
             2.6454E+00
 PARAMETER:  1.7901E-02  4.7207E-01  9.2913E-01 -1.5881E-01  4.5641E-01  1.5851E-02 -1.1838E-01  1.3079E+00  1.7966E-01  4.2183E-01
             1.0728E+00
 GRADIENT:  -1.7268E+01  3.8857E+01 -7.5588E+00 -8.8589E+00  1.2539E+01 -1.2762E+01  1.2890E+00  4.2392E+00 -6.8346E+00 -2.0829E+00
             3.9544E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2631.61133837552        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  9.4148E-01  1.5027E+00  2.8097E+00  7.5366E-01  1.5345E+00  9.5378E-01  7.6053E-01  3.7902E+00  1.2081E+00  1.4449E+00
             2.6239E+00
 PARAMETER:  3.9701E-02  5.0724E-01  1.1331E+00 -1.8281E-01  5.2821E-01  5.2675E-02 -1.7374E-01  1.4324E+00  2.8908E-01  4.6805E-01
             1.0647E+00
 GRADIENT:  -1.6313E+01 -2.1242E+01 -7.4620E+00 -4.5105E+00  8.6628E+00 -2.4672E+00  4.7480E+00  2.8895E+00  1.0045E+00 -6.0996E+00
             9.0409E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2632.07242645191        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      555
 NPARAMETR:  9.4283E-01  1.5285E+00  2.8281E+00  7.3894E-01  1.5497E+00  9.5449E-01  7.3395E-01  3.8139E+00  1.2449E+00  1.4706E+00
             2.6246E+00
 PARAMETER:  4.1135E-02  5.2426E-01  1.1396E+00 -2.0254E-01  5.3805E-01  5.3422E-02 -2.0931E-01  1.4386E+00  3.1907E-01  4.8569E-01
             1.0649E+00
 GRADIENT:  -1.3068E+01 -1.6642E+01 -6.4608E+00 -1.2616E+00  8.5618E+00 -2.2083E+00  4.5848E+00  2.7715E+00  9.4862E-01 -4.8090E+00
             9.6849E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2632.24928923064        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      744            RESET HESSIAN, TYPE II
 NPARAMETR:  9.4811E-01  1.5438E+00  2.8114E+00  7.3950E-01  1.5445E+00  9.6030E-01  6.6864E-01  3.7805E+00  1.2367E+00  1.4667E+00
             2.6382E+00
 PARAMETER:  4.6711E-02  5.3426E-01  1.1337E+00 -2.0179E-01  5.3473E-01  5.9493E-02 -3.0251E-01  1.4299E+00  3.1243E-01  4.8301E-01
             1.0701E+00
 GRADIENT:   5.5129E+01  1.0004E+02 -2.6773E+00  2.3036E+01  2.8952E+01  4.9506E+00  1.4294E+00  6.4974E+00 -4.2828E+00 -1.0403E+00
             3.5863E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2632.98605322675        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.4675E-01  1.5259E+00  2.8101E+00  7.3273E-01  1.5439E+00  9.5730E-01  5.9272E-01  3.7759E+00  1.3702E+00  1.4662E+00
             2.6386E+00
 PARAMETER:  4.5275E-02  5.2258E-01  1.1332E+00 -2.1098E-01  5.3428E-01  5.6363E-02 -4.2303E-01  1.4286E+00  4.1495E-01  4.8269E-01
             1.0703E+00
 GRADIENT:  -3.4902E+00 -1.7679E+01 -1.0090E+01 -1.9820E+00  2.2201E+00 -7.0629E-01  4.5807E+00  1.3336E+00  2.7600E+00 -5.0416E+00
             1.5574E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2634.63146535270        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  9.4672E-01  1.5484E+00  2.8174E+00  6.8784E-01  1.5435E+00  9.5718E-01  1.8055E-01  3.7717E+00  1.6504E+00  1.4673E+00
             2.6237E+00
 PARAMETER:  4.5253E-02  5.3719E-01  1.1358E+00 -2.7419E-01  5.3406E-01  5.6231E-02 -1.6118E+00  1.4275E+00  6.0101E-01  4.8339E-01
             1.0646E+00
 GRADIENT:  -2.5252E+00 -6.8178E+00 -8.6417E+00 -5.9153E+00 -6.6454E+00 -7.7801E-01  3.0465E-01 -3.1485E+00  1.1328E+00 -1.4369E+01
            -4.9367E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2634.80994311577        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1225
 NPARAMETR:  9.4905E-01  1.5206E+00  2.8256E+00  7.1256E-01  1.5442E+00  9.5977E-01  9.5368E-02  3.7738E+00  1.6251E+00  1.4688E+00
             2.6157E+00
 PARAMETER:  4.7805E-02  5.1901E-01  1.1390E+00 -2.3895E-01  5.3461E-01  5.8956E-02 -2.2277E+00  1.4284E+00  5.8542E-01  4.8455E-01
             1.0613E+00
 GRADIENT:   2.6155E+00 -1.0818E+03  4.8065E+02 -2.3529E+03  1.0478E+03  1.2726E-01  4.9671E-02  3.9114E+02 -9.6087E+02  1.1461E+03
            -5.5144E+02
 NUMSIGDIG:         1.6         2.3         2.3         2.3         2.3         2.5         0.6         2.3         2.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1225
 NO. OF SIG. DIGITS IN FINAL EST.:  0.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4113E-03 -1.5286E-02 -1.6774E-02  5.6558E-03 -2.1564E-02
 SE:             2.9405E-02  3.7458E-03  1.7059E-02  2.7736E-02  2.5880E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3464E-01  4.4886E-05  3.2547E-01  8.3842E-01  4.0471E-01

 ETASHRINKSD(%)  1.4904E+00  8.7451E+01  4.2849E+01  7.0807E+00  1.3300E+01
 ETASHRINKVR(%)  2.9586E+00  9.8425E+01  6.7338E+01  1.3660E+01  2.4831E+01
 EBVSHRINKSD(%)  1.8226E+00  8.9816E+01  4.9518E+01  6.5354E+00  1.4148E+01
 EBVSHRINKVR(%)  3.6120E+00  9.8963E+01  7.4516E+01  1.2644E+01  2.6294E+01
 RELATIVEINF(%)  9.6305E+01  2.0770E-01  1.4963E+01  1.9375E+01  4.3771E+01
 EPSSHRINKSD(%)  1.6300E+01
 EPSSHRINKVR(%)  2.9944E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          880
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1617.3318184402240     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2634.8099431157652     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1017.4781246755413     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2634.810       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  1.52E+00  2.83E+00  7.13E-01  1.54E+00  9.60E-01  9.75E-02  3.77E+00  1.62E+00  1.47E+00  2.62E+00
 


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
 #CPUT: Total CPU Time in Seconds,       85.329
Stop Time:
Sat Oct 23 15:21:37 CDT 2021
