Sat Oct 23 15:36:22 CDT 2021
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
$DATA ../../../../data/SD1/SL3/dat81.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      982
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

 TOT. NO. OF OBS RECS:      882
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1479.54155578875        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2440E+02  3.9974E+01  2.0060E+02  1.0131E+02  2.0586E+02  2.6522E+01 -1.1388E+02 -4.1231E+02 -5.6842E+01 -5.7221E+01
            -9.9268E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2324.11655108319        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1258E+00  1.2486E+00  9.7828E-01  1.0191E+00  9.6250E-01  1.1829E+00  1.0481E+00  1.0314E+00  7.7240E-01  9.3021E-01
             5.3080E+00
 PARAMETER:  2.1848E-01  3.2199E-01  7.8037E-02  1.1889E-01  6.1782E-02  2.6799E-01  1.4699E-01  1.3091E-01 -1.5825E-01  2.7654E-02
             1.7692E+00
 GRADIENT:   8.4395E+01  7.2035E+01  2.6418E+00  4.0893E+01 -5.6962E+01  3.6463E+01  1.0599E+01  6.2117E+00  1.4026E+01  1.1996E+01
             7.5752E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2484.35550053043        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0245E+00  1.1139E+00  2.0814E+00  1.0182E+00  1.2731E+00  1.0633E+00  1.4061E+00  2.8314E+00  5.8465E-01  4.6886E-01
             4.0489E+00
 PARAMETER:  1.2418E-01  2.0790E-01  8.3306E-01  1.1803E-01  3.4145E-01  1.6136E-01  4.4080E-01  1.1408E+00 -4.3674E-01 -6.5746E-01
             1.4984E+00
 GRADIENT:  -2.6634E+01 -1.7003E+01 -1.6784E+01 -3.2208E+01  4.3898E+01  2.3949E-01  3.2884E+01  1.6946E+01  4.9011E+00 -3.7178E+00
             4.5094E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2557.31436951555        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0064E+00  1.2359E+00  1.0475E+00  8.9210E-01  1.1099E+00  1.0645E+00  1.0231E+00  1.6400E-01  7.4018E-01  1.1425E+00
             3.0948E+00
 PARAMETER:  1.0640E-01  3.1180E-01  1.4638E-01 -1.4179E-02  2.0425E-01  1.6246E-01  1.2280E-01 -1.7079E+00 -2.0086E-01  2.3323E-01
             1.2297E+00
 GRADIENT:  -1.6045E+01 -8.9385E+00 -6.0725E+00 -2.2396E+01  2.1798E+01  3.7129E+00  5.0099E+00  7.7707E-02 -6.4520E+00  4.3050E+00
            -1.7225E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2558.17615202470        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      321
 NPARAMETR:  1.0151E+00  1.1738E+00  9.8221E-01  9.3244E-01  1.0305E+00  1.0598E+00  1.0123E+00  1.9068E-01  8.2681E-01  1.0238E+00
             3.1207E+00
 PARAMETER:  1.1498E-01  2.6021E-01  8.2053E-02  3.0050E-02  1.3006E-01  1.5811E-01  1.1219E-01 -1.5572E+00 -9.0178E-02  1.2349E-01
             1.2381E+00
 GRADIENT:  -5.7303E+01 -3.1400E+01 -5.5998E+00 -2.5084E+01  6.0890E+00 -1.4542E+01 -2.3798E-01  1.2378E-01  4.4657E-01 -3.7764E+00
            -2.0368E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2561.69513499292        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      498
 NPARAMETR:  1.0466E+00  1.4930E+00  1.1096E+00  7.7002E-01  1.2794E+00  1.0945E+00  8.4400E-01  1.0259E-01  9.1845E-01  1.2415E+00
             3.1470E+00
 PARAMETER:  1.4551E-01  5.0080E-01  2.0402E-01 -1.6134E-01  3.4639E-01  1.9026E-01 -6.9606E-02 -2.1770E+00  1.4934E-02  3.1633E-01
             1.2465E+00
 GRADIENT:  -1.1291E+00  3.3483E+00  1.4643E-01  4.0163E+00  4.3050E-01 -4.4424E-01  4.4569E-01  6.6020E-03  1.8335E-01  2.1457E-01
            -5.1080E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2561.71868541698        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  1.0473E+00  1.5293E+00  1.0919E+00  7.4447E-01  1.3021E+00  1.0959E+00  8.2460E-01  5.2616E-02  9.3849E-01  1.2570E+00
             3.1463E+00
 PARAMETER:  1.4617E-01  5.2479E-01  1.8792E-01 -1.9509E-01  3.6398E-01  1.9162E-01 -9.2851E-02 -2.8447E+00  3.6518E-02  3.2872E-01
             1.2462E+00
 GRADIENT:   7.6175E+01  7.6683E+01  3.5005E-01  1.2105E+01  1.4283E+01  1.8987E+01  8.6148E-01  3.3540E-03  4.6313E-01  1.7898E+00
             2.1791E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2561.71964809047        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      830
 NPARAMETR:  1.0472E+00  1.5291E+00  1.0900E+00  7.4448E-01  1.3006E+00  1.0957E+00  8.2472E-01  1.5740E-02  9.3858E-01  1.2563E+00
             3.1466E+00
 PARAMETER:  1.4608E-01  5.2470E-01  1.8619E-01 -1.9508E-01  3.6283E-01  1.9138E-01 -9.2705E-02 -4.0515E+00  3.6617E-02  3.2817E-01
             1.2463E+00
 GRADIENT:  -9.0842E-03  3.4265E-02 -3.5442E-02  3.1147E-02 -9.1792E-03  4.9829E-03  7.7393E-03  1.5600E-04  5.5953E-03  1.7272E-04
            -1.7444E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2561.71973972428        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:      888
 NPARAMETR:  1.0472E+00  1.5290E+00  1.0902E+00  7.4455E-01  1.3005E+00  1.0957E+00  8.2470E-01  1.0000E-02  9.3867E-01  1.2563E+00
             3.1466E+00
 PARAMETER:  1.4613E-01  5.2463E-01  1.8635E-01 -1.9497E-01  3.6271E-01  1.9143E-01 -9.2740E-02 -4.6316E+00  3.6704E-02  3.2816E-01
             1.2463E+00
 GRADIENT:   8.3137E-02  7.8152E-02 -5.2180E-03  1.4355E-02 -1.0097E-01  2.4834E-02  2.4183E-03  0.0000E+00  1.2502E-02  3.9537E-03
            -2.8943E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      888
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6522E-03 -1.2261E-02 -7.1566E-05  1.5364E-03 -1.3548E-02
 SE:             2.9324E-02  2.2104E-02  6.4810E-05  1.9356E-02  2.4971E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5507E-01  5.7911E-01  2.6948E-01  9.3673E-01  5.8745E-01

 ETASHRINKSD(%)  1.7616E+00  2.5950E+01  9.9783E+01  3.5156E+01  1.6345E+01
 ETASHRINKVR(%)  3.4922E+00  4.5166E+01  1.0000E+02  5.7953E+01  3.0018E+01
 EBVSHRINKSD(%)  1.8047E+00  2.6232E+01  9.9787E+01  3.5868E+01  1.5963E+01
 EBVSHRINKVR(%)  3.5769E+00  4.5583E+01  1.0000E+02  5.8870E+01  2.9379E+01
 RELATIVEINF(%)  9.6309E+01  4.1923E+00  1.7148E-04  3.3750E+00  1.0222E+01
 EPSSHRINKSD(%)  1.5202E+01
 EPSSHRINKVR(%)  2.8092E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          882
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1621.0075725730426     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2561.7197397242849     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -940.71216715124228     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2561.720       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.53E+00  1.09E+00  7.45E-01  1.30E+00  1.10E+00  8.25E-01  1.00E-02  9.39E-01  1.26E+00  3.15E+00
 


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
 #CPUT: Total CPU Time in Seconds,       52.815
Stop Time:
Sat Oct 23 15:36:32 CDT 2021
