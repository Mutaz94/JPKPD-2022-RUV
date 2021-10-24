Sun Oct 24 01:35:21 CDT 2021
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
$DATA ../../../../data/SD4/B/dat22.csv ignore=@
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
Current Date:       24 OCT 2021
Days until program expires : 175
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
 RAW OUTPUT FILE (FILE): m22.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1701.76839656284        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9799E+02  1.9640E+01 -3.5731E+01  9.4147E+01  5.6093E+01  8.2005E+01  8.0856E+00  4.9055E+00  2.0685E+01  1.6293E+01
             5.3580E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1711.67763698191        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      133
 NPARAMETR:  9.5871E-01  1.0262E+00  1.0411E+00  9.6227E-01  9.5998E-01  7.6512E-01  9.7633E-01  9.8699E-01  9.4644E-01  9.3359E-01
             8.4605E-01
 PARAMETER:  5.7830E-02  1.2586E-01  1.4031E-01  6.1544E-02  5.9159E-02 -1.6772E-01  7.6043E-02  8.6901E-02  4.4951E-02  3.1281E-02
            -6.7180E-02
 GRADIENT:  -6.1847E+00  7.8619E+00  2.4366E+01 -2.4152E+01 -4.5412E+01 -4.3934E+01  3.0726E-01 -3.2542E+00 -3.3227E+00  6.1357E+00
            -1.3854E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1713.39837742119        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  9.6496E-01  1.1763E+00  9.5138E-01  8.7619E-01  9.8778E-01  8.1752E-01  9.4523E-01  1.2314E+00  9.5371E-01  8.2957E-01
             8.4353E-01
 PARAMETER:  6.4335E-02  2.6238E-01  5.0159E-02 -3.2173E-02  8.7701E-02 -1.0148E-01  4.3677E-02  3.0814E-01  5.2599E-02 -8.6853E-02
            -7.0163E-02
 GRADIENT:   1.2676E+01  1.9414E+01  1.4376E+01 -5.4423E+00 -2.6318E+01 -1.4780E+01  6.4865E-01  1.7573E+00 -9.4631E+00 -4.8011E+00
            -1.6058E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1715.18642852267        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      484
 NPARAMETR:  9.6393E-01  1.0978E+00  8.8691E-01  9.1846E-01  9.4693E-01  8.4577E-01  9.5880E-01  9.2544E-01  9.7936E-01  8.4961E-01
             8.7070E-01
 PARAMETER:  6.3268E-02  1.9335E-01 -2.0011E-02  1.4944E-02  4.5469E-02 -6.7512E-02  5.7927E-02  2.2519E-02  7.9143E-02 -6.2979E-02
            -3.8458E-02
 GRADIENT:   7.4720E+00 -1.0003E+00  2.1469E-02 -3.0647E-01 -1.0854E+00 -9.7675E-01  1.5615E-01  4.5129E-01  5.9803E-01  3.3839E-01
            -8.5284E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1715.28113972719        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  9.6292E-01  1.2729E+00  7.2476E-01  8.0367E-01  9.5406E-01  8.4955E-01  8.8575E-01  7.4834E-01  1.0585E+00  8.3360E-01
             8.7298E-01
 PARAMETER:  6.2216E-02  3.4132E-01 -2.2192E-01 -1.1857E-01  5.2973E-02 -6.3044E-02 -2.1320E-02 -1.8989E-01  1.5683E-01 -8.2000E-02
            -3.5840E-02
 GRADIENT:   1.7078E-01  5.8199E+00  1.1924E+00  4.8933E+00 -1.9464E+00  2.6247E-01 -4.5552E-02 -3.0761E-01 -3.5383E-01 -6.0197E-01
             7.9683E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1715.29684609632        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  9.6389E-01  1.3547E+00  6.5810E-01  7.4824E-01  9.6379E-01  8.4907E-01  8.5283E-01  6.6721E-01  1.1073E+00  8.3334E-01
             8.7242E-01
 PARAMETER:  6.3219E-02  4.0360E-01 -3.1841E-01 -1.9003E-01  6.3114E-02 -6.3616E-02 -5.9197E-02 -3.0465E-01  2.0190E-01 -8.2311E-02
            -3.6489E-02
 GRADIENT:   2.1105E+00  5.8419E+00  1.0861E+00  4.4663E+00 -1.9745E+00 -1.2721E-01 -1.1389E-01 -2.5350E-01 -3.0789E-01 -5.5003E-01
            -1.6814E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1715.31209405879        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1015             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6398E-01  1.3848E+00  6.3104E-01  7.2255E-01  9.7022E-01  8.4947E-01  8.4036E-01  6.5652E-01  1.1309E+00  8.3705E-01
             8.7251E-01
 PARAMETER:  6.3313E-02  4.2553E-01 -3.6038E-01 -2.2496E-01  6.9767E-02 -6.3140E-02 -7.3921E-02 -3.2080E-01  2.2303E-01 -7.7876E-02
            -3.6382E-02
 GRADIENT:   5.4446E+02  3.8002E+02  4.3762E+00  1.0356E+02  1.2793E+01  4.3547E+01  6.2693E+00  3.6654E-01  1.3293E+01  1.2900E+00
             1.0757E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1715.31690546754        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  9.6388E-01  1.3862E+00  6.3203E-01  7.2278E-01  9.7010E-01  8.4954E-01  8.3974E-01  6.5095E-01  1.1308E+00  8.3566E-01
             8.7234E-01
 PARAMETER:  6.3215E-02  4.2659E-01 -3.5883E-01 -2.2465E-01  6.9640E-02 -6.3058E-02 -7.4659E-02 -3.2932E-01  2.2294E-01 -7.9538E-02
            -3.6575E-02
 GRADIENT:  -2.0198E-01 -1.5284E-01 -3.1553E-01  1.6258E-01  9.3615E-02 -4.6703E-02  6.8469E-02  9.1261E-03 -4.7209E-02  3.7537E-02
             4.6742E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1185
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2172E-04 -2.1567E-02 -1.9631E-02  1.5494E-02 -3.0613E-02
 SE:             2.9856E-02  2.2491E-02  8.5222E-03  2.3975E-02  2.1977E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9407E-01  3.3761E-01  2.1254E-02  5.1810E-01  1.6363E-01

 ETASHRINKSD(%)  1.0000E-10  2.4652E+01  7.1449E+01  1.9681E+01  2.6375E+01
 ETASHRINKVR(%)  1.0000E-10  4.3227E+01  9.1849E+01  3.5489E+01  4.5794E+01
 EBVSHRINKSD(%)  4.4210E-01  2.4374E+01  7.4199E+01  2.0459E+01  2.4952E+01
 EBVSHRINKVR(%)  8.8225E-01  4.2807E+01  9.3343E+01  3.6732E+01  4.3677E+01
 RELATIVEINF(%)  9.8992E+01  2.5398E+00  5.0879E-01  3.0755E+00  7.9512E+00
 EPSSHRINKSD(%)  4.5660E+01
 EPSSHRINKVR(%)  7.0472E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1715.3169054675366     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -980.16607890379839     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1715.317       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.64E-01  1.39E+00  6.32E-01  7.23E-01  9.70E-01  8.50E-01  8.40E-01  6.51E-01  1.13E+00  8.36E-01  8.72E-01
 


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
 #CPUT: Total CPU Time in Seconds,       33.477
Stop Time:
Sun Oct 24 01:35:29 CDT 2021
