Sat Oct 23 22:01:51 CDT 2021
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
$DATA ../../../../data/SD3/A2/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1170.87036066278        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5294E+02  8.0731E+01  1.3983E+02  3.9301E+01  8.7669E+01  4.8484E+01 -4.4405E+01 -1.5889E+02 -5.2965E+01 -7.3294E+01
            -1.5212E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1700.51070752608        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0258E+00  7.8775E-01  8.6223E-01  1.1557E+00  8.7802E-01  1.0599E+00  1.0382E+00  4.6754E-01  1.0217E+00  1.0609E+00
             2.2525E+00
 PARAMETER:  1.2548E-01 -1.3857E-01 -4.8232E-02  2.4472E-01 -3.0090E-02  1.5820E-01  1.3753E-01 -6.6028E-01  1.2147E-01  1.5912E-01
             9.1203E-01
 GRADIENT:   1.5412E+02  8.1178E-01 -2.6067E+01  7.1183E+01  6.7902E+01  4.6813E+01 -3.7996E+00  8.3626E-01  8.9071E+00  4.3741E-01
            -3.0198E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1714.93722455673        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0126E+00  8.0896E-01  4.1461E-01  1.0701E+00  5.3900E-01  9.1566E-01  1.2055E+00  1.3404E-01  9.4169E-01  6.1921E-01
             2.4246E+00
 PARAMETER:  1.1253E-01 -1.1201E-01 -7.8042E-01  1.6777E-01 -5.1803E-01  1.1894E-02  2.8692E-01 -1.9097E+00  3.9923E-02 -3.7931E-01
             9.8566E-01
 GRADIENT:   8.2784E+01  2.8847E+01 -1.3610E+01  5.1234E+01  4.7312E+01 -1.1300E+01  4.4808E+00  2.2496E-01 -8.2866E+00  3.0095E+00
             6.2038E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1718.15638581776        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      275
 NPARAMETR:  9.9086E-01  7.0742E-01  3.9128E-01  1.1186E+00  4.8424E-01  9.0716E-01  1.4591E+00  1.3198E-01  9.3033E-01  5.5732E-01
             2.3615E+00
 PARAMETER:  9.0816E-02 -2.4613E-01 -8.3833E-01  2.1207E-01 -6.2518E-01  2.5633E-03  4.7779E-01 -1.9251E+00  2.7788E-02 -4.8461E-01
             9.5929E-01
 GRADIENT:  -4.6234E+01  3.3179E+01 -2.7445E+01  4.4296E+01  3.0069E+01 -2.2054E+01  1.2519E+01  1.6144E-01 -7.0807E+00 -4.0269E-01
             4.2465E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1730.02524859346        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  1.0073E+00  3.5978E-01  3.4582E-01  1.2069E+00  3.5511E-01  9.6423E-01  1.9646E+00  1.0012E-01  9.3986E-01  5.7669E-01
             2.1915E+00
 PARAMETER:  1.0727E-01 -9.2227E-01 -9.6184E-01  2.8808E-01 -9.3532E-01  6.3572E-02  7.7529E-01 -2.2014E+00  3.7978E-02 -4.5045E-01
             8.8457E-01
 GRADIENT:   1.2831E+01  6.1174E+00  2.6896E+01  7.8368E+00 -3.9400E+01  3.2413E+00  2.0207E+00 -6.9655E-02 -1.5638E+00 -3.3169E+00
            -2.6227E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1730.64326450173        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      627            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0031E+00  3.9174E-01  3.2624E-01  1.1775E+00  3.5309E-01  9.5813E-01  1.7706E+00  1.0552E-01  9.5855E-01  5.7452E-01
             2.1833E+00
 PARAMETER:  1.0309E-01 -8.3716E-01 -1.0201E+00  2.6342E-01 -9.4104E-01  5.7224E-02  6.7132E-01 -2.1489E+00  5.7669E-02 -4.5423E-01
             8.8085E-01
 GRADIENT:   7.8999E+01  1.2373E+01  2.1372E+01  6.2679E+01  7.7932E+01  6.1008E+00  3.3949E+00  1.8664E-02  2.6304E+00  1.2940E+00
             8.9854E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1730.66374725661        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      802
 NPARAMETR:  1.0023E+00  3.9271E-01  3.2638E-01  1.1768E+00  3.5350E-01  9.5807E-01  1.7574E+00  1.8206E-01  9.5939E-01  5.7231E-01
             2.1817E+00
 PARAMETER:  1.0233E-01 -8.3468E-01 -1.0197E+00  2.6276E-01 -9.3986E-01  5.7161E-02  6.6381E-01 -1.6034E+00  5.8544E-02 -4.5807E-01
             8.8012E-01
 GRADIENT:  -1.5317E+00 -1.8969E+00 -4.3062E+00 -5.6643E-01  6.8632E+00 -1.1157E-01 -4.3321E-01 -3.8791E-02 -1.7826E-02  7.8011E-01
             5.4696E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1730.72472131148        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      981
 NPARAMETR:  1.0033E+00  3.9590E-01  3.2649E-01  1.1763E+00  3.5308E-01  9.5889E-01  1.7669E+00  3.2402E-01  9.5974E-01  5.4256E-01
             2.1764E+00
 PARAMETER:  1.0327E-01 -8.2659E-01 -1.0194E+00  2.6242E-01 -9.4105E-01  5.8025E-02  6.6923E-01 -1.0269E+00  5.8907E-02 -5.1145E-01
             8.7767E-01
 GRADIENT:   3.4844E-01  4.6696E-01  1.4731E+00  6.3277E-01 -3.7820E-01  3.1080E-02  1.0581E-01 -1.6179E-02  9.3562E-02  7.0200E-02
            -7.9659E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1730.72868286079        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1159            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0032E+00  3.9604E-01  3.2586E-01  1.1753E+00  3.5255E-01  9.5891E-01  1.7656E+00  3.3361E-01  9.5975E-01  5.3850E-01
             2.1752E+00
 PARAMETER:  1.0323E-01 -8.2623E-01 -1.0213E+00  2.6152E-01 -9.4257E-01  5.8039E-02  6.6848E-01 -9.9777E-01  5.8917E-02 -5.1898E-01
             8.7713E-01
 GRADIENT:   7.9552E+01  1.4309E+01  2.7133E+01  6.1926E+01  7.2901E+01  6.1977E+00  3.5696E+00  1.9303E-01  2.7482E+00  1.0600E+00
             8.4681E+00

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1730.72992689326        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  1.0032E+00  3.9558E-01  3.2554E-01  1.1754E+00  3.5265E-01  9.5890E-01  1.7650E+00  3.3696E-01  9.5978E-01  5.3916E-01
             2.1756E+00
 PARAMETER:  1.0318E-01 -8.2740E-01 -1.0223E+00  2.6159E-01 -9.4229E-01  5.8030E-02  6.6817E-01 -9.8780E-01  5.8953E-02 -5.1775E-01
             8.7732E-01
 GRADIENT:   1.1328E-01 -2.5515E-01 -3.8162E-02  8.6714E-02  2.1529E+00 -1.1800E-03 -1.0468E-03  1.0794E-02 -3.8075E-02  9.6202E-02
             1.4611E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1252
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3864E-03  2.9864E-02 -9.9578E-03 -1.4891E-02  1.1020E-02
 SE:             2.9414E-02  1.8393E-02  7.5569E-03  2.7098E-02  1.7518E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6241E-01  1.0445E-01  1.8760E-01  5.8265E-01  5.2930E-01

 ETASHRINKSD(%)  1.4608E+00  3.8381E+01  7.4684E+01  9.2186E+00  4.1314E+01
 ETASHRINKVR(%)  2.9003E+00  6.2030E+01  9.3591E+01  1.7587E+01  6.5559E+01
 EBVSHRINKSD(%)  1.6216E+00  4.1440E+01  7.4084E+01  8.3652E+00  3.9419E+01
 EBVSHRINKVR(%)  3.2168E+00  6.5707E+01  9.3284E+01  1.6031E+01  6.3299E+01
 RELATIVEINF(%)  9.6471E+01  8.2788E+00  5.2315E-01  4.7002E+01  2.0700E+00
 EPSSHRINKSD(%)  2.9392E+01
 EPSSHRINKVR(%)  5.0146E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1730.7299268932636     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -811.79139368859092     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1730.730       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  3.96E-01  3.26E-01  1.18E+00  3.53E-01  9.59E-01  1.77E+00  3.37E-01  9.60E-01  5.39E-01  2.18E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,      101.095
Stop Time:
Sat Oct 23 22:02:07 CDT 2021
