Sat Oct 23 22:42:26 CDT 2021
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
$DATA ../../../../data/SD3/A3/dat64.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   150.634875811999        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2333E+02  4.7550E+01  1.4119E+02  1.4283E+01  2.5652E+02  3.5535E+01 -6.8037E+01 -2.0015E+02 -9.8263E+01 -1.6426E+02
            -3.9605E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1493.90241193408        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0565E+00  1.0550E+00  9.5829E-01  1.1023E+00  9.4604E-01  8.8684E-01  9.3860E-01  9.6065E-01  7.6880E-01  8.9139E-01
             4.6834E+00
 PARAMETER:  1.5494E-01  1.5354E-01  5.7399E-02  1.9742E-01  4.4531E-02 -2.0085E-02  3.6632E-02  5.9860E-02 -1.6293E-01 -1.4979E-02
             1.6440E+00
 GRADIENT:   1.9661E+01  1.2258E+01 -1.8851E+01  4.1374E+01  1.5891E+00 -2.8623E+01  3.9228E+00  5.6625E+00  4.2121E+00  1.8290E+01
             2.3620E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1518.06355964530        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0323E+00  7.5993E-01  3.7467E-01  1.1820E+00  4.6636E-01  9.5823E-01  1.1857E+00  1.4682E-01  9.1777E-01  2.9569E-01
             4.1324E+00
 PARAMETER:  1.3183E-01 -1.7453E-01 -8.8170E-01  2.6722E-01 -6.6281E-01  5.7336E-02  2.7035E-01 -1.8185E+00  1.4190E-02 -1.1184E+00
             1.5188E+00
 GRADIENT:  -2.4917E+01  2.1086E+01 -3.6099E+01  1.1656E+02  3.7661E+01 -1.2556E+01  8.0942E+00  1.8913E-01  4.1942E+00  3.2283E+00
             1.7121E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1538.12675272647        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      237
 NPARAMETR:  1.0001E+00  6.6220E-01  2.7914E-01  1.1013E+00  3.7349E-01  1.0185E+00  1.0611E+00  2.4282E-02  1.0794E+00  3.5538E-01
             3.1358E+00
 PARAMETER:  1.0011E-01 -3.1219E-01 -1.1761E+00  1.9653E-01 -8.8487E-01  1.1831E-01  1.5933E-01 -3.6180E+00  1.7637E-01 -9.3456E-01
             1.2429E+00
 GRADIENT:  -3.3981E+01  3.0219E+00 -1.7367E+01  6.5467E+01  8.5694E+01  5.5826E+00 -8.7477E+00 -1.0396E-02 -9.5455E-01 -5.7763E+00
            -3.5237E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1540.85034277267        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0270E+00  6.1822E-01  2.7339E-01  1.1265E+00  3.4465E-01  1.0091E+00  1.1343E+00  1.5535E-02  1.0766E+00  4.4256E-01
             3.0900E+00
 PARAMETER:  1.2667E-01 -3.8092E-01 -1.1968E+00  2.1914E-01 -9.6523E-01  1.0910E-01  2.2604E-01 -4.0647E+00  1.7382E-01 -7.1517E-01
             1.2282E+00
 GRADIENT:  -2.0666E+01  3.9477E+01  1.1024E+01  6.4256E+01 -3.9262E+01 -9.7024E-01 -4.6049E+00 -4.6408E-03 -6.8126E+00 -6.0238E+00
            -5.0917E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1545.15217278828        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0368E+00  5.0899E-01  2.6904E-01  1.1132E+00  3.1978E-01  1.0087E+00  1.2319E+00  1.0000E-02  1.0809E+00  5.2376E-01
             3.2287E+00
 PARAMETER:  1.3618E-01 -5.7534E-01 -1.2129E+00  2.0722E-01 -1.0401E+00  1.0868E-01  3.0854E-01 -5.0016E+00  1.7776E-01 -5.4672E-01
             1.2721E+00
 GRADIENT:  -1.4608E+00  1.7959E+00 -3.6651E-01  3.0821E+00 -1.6047E+00  1.5634E+00  5.4958E-01  0.0000E+00  1.1519E-01 -1.1285E-01
             8.5207E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1545.22386981687        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0371E+00  4.9953E-01  2.6721E-01  1.1115E+00  3.1640E-01  1.0047E+00  1.1931E+00  1.0000E-02  1.0870E+00  5.4486E-01
             3.1936E+00
 PARAMETER:  1.3639E-01 -5.9409E-01 -1.2197E+00  2.0570E-01 -1.0507E+00  1.0468E-01  2.7659E-01 -5.1509E+00  1.8343E-01 -5.0722E-01
             1.2611E+00
 GRADIENT:   1.3991E-01 -8.8305E-02  6.1775E-01  9.9189E-03 -4.7244E-01 -7.0188E-05 -4.1650E-02  0.0000E+00 -1.2066E-01 -1.3942E-01
            -3.5615E-01

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1545.22386981687        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      749
 NPARAMETR:  1.0371E+00  4.9953E-01  2.6721E-01  1.1115E+00  3.1640E-01  1.0047E+00  1.1931E+00  1.0000E-02  1.0870E+00  5.4486E-01
             3.1936E+00
 PARAMETER:  1.3639E-01 -5.9409E-01 -1.2197E+00  2.0570E-01 -1.0507E+00  1.0468E-01  2.7659E-01 -5.1509E+00  1.8343E-01 -5.0722E-01
             1.2611E+00
 GRADIENT:   1.3991E-01 -8.8305E-02  6.1775E-01  9.9189E-03 -4.7244E-01 -7.0188E-05 -4.1650E-02  0.0000E+00 -1.2066E-01 -1.3942E-01
            -3.5615E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      749
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.0073E-04  1.0355E-02 -6.2422E-05 -1.2798E-02  3.8435E-03
 SE:             2.8940E-02  1.5718E-02  1.9474E-04  2.6026E-02  1.7024E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7517E-01  5.1004E-01  7.4857E-01  6.2290E-01  8.2138E-01

 ETASHRINKSD(%)  3.0485E+00  4.7342E+01  9.9348E+01  1.2810E+01  4.2967E+01
 ETASHRINKVR(%)  6.0040E+00  7.2271E+01  9.9996E+01  2.3978E+01  6.7472E+01
 EBVSHRINKSD(%)  2.8951E+00  4.8354E+01  9.9377E+01  1.1559E+01  4.2797E+01
 EBVSHRINKVR(%)  5.7063E+00  7.3326E+01  9.9996E+01  2.1781E+01  6.7278E+01
 RELATIVEINF(%)  9.4032E+01  3.4308E+00  3.1201E-04  4.5061E+01  1.3079E+00
 EPSSHRINKSD(%)  2.4714E+01
 EPSSHRINKVR(%)  4.3320E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1545.2238698168740     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -626.28533661220126     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1545.224       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  5.00E-01  2.67E-01  1.11E+00  3.16E-01  1.00E+00  1.19E+00  1.00E-02  1.09E+00  5.45E-01  3.19E+00
 


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
 #CPUT: Total CPU Time in Seconds,       59.047
Stop Time:
Sat Oct 23 22:42:37 CDT 2021
