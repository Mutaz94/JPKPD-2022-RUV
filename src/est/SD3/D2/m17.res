Sun Oct 24 01:11:31 CDT 2021
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
$DATA ../../../../data/SD3/D2/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1926.32416166158        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.7018E+02 -6.0760E+01  7.0830E+00 -1.3739E+02  2.6094E+01 -8.4444E+01 -6.4748E+01 -1.1482E+01 -1.3517E+02 -4.6697E+01
            -3.0465E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1990.43775678313        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      132
 NPARAMETR:  1.1143E+00  9.7032E-01  9.9713E-01  9.7573E-01  9.7728E-01  1.2764E+00  1.3341E+00  1.0601E+00  1.5973E+00  1.2635E+00
             1.0522E+00
 PARAMETER:  2.0827E-01  6.9870E-02  9.7123E-02  7.5434E-02  7.7023E-02  3.4401E-01  3.8827E-01  1.5834E-01  5.6830E-01  3.3386E-01
             1.5087E-01
 GRADIENT:   5.7140E+01 -6.8458E+01 -1.4710E+01 -1.1438E+02 -2.4491E+01 -3.6344E+01 -6.4669E+00  5.3053E-01  2.5013E+01  1.6951E+00
             2.5671E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2003.76226115979        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      310
 NPARAMETR:  1.1275E+00  1.2341E+00  1.3661E+00  9.2861E-01  1.2681E+00  1.6239E+00  1.7306E+00  1.4141E+00  1.1367E+00  1.6777E+00
             1.0841E+00
 PARAMETER:  2.1996E-01  3.1035E-01  4.1195E-01  2.5931E-02  3.3748E-01  5.8484E-01  6.4846E-01  4.4648E-01  2.2812E-01  6.1745E-01
             1.8078E-01
 GRADIENT:   4.8902E+01 -2.8638E+01  3.1322E+01 -7.6861E+01 -5.1834E+01  6.4809E+01 -1.0489E+01 -8.9496E+00 -1.0513E+01  2.6061E+01
             4.1997E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2017.57582308462        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      489
 NPARAMETR:  1.0881E+00  1.0140E+00  1.2848E+00  1.1098E+00  1.1841E+00  1.3757E+00  2.0040E+00  1.2095E+00  1.2034E+00  1.2948E+00
             1.0424E+00
 PARAMETER:  1.8447E-01  1.1394E-01  3.5057E-01  2.0422E-01  2.6899E-01  4.1895E-01  7.9514E-01  2.9020E-01  2.8515E-01  3.5834E-01
             1.4148E-01
 GRADIENT:   2.3446E+01 -1.0015E+01  3.5289E+00 -1.9663E+01  3.8902E+00  1.9333E+00 -3.3903E-01 -6.4726E+00  1.2409E+01 -2.1695E+00
             1.2498E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2019.54218725408        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.0686E+00  8.3069E-01  1.4274E+00  1.2794E+00  1.1547E+00  1.3535E+00  2.3770E+00  1.3500E+00  1.0589E+00  1.3392E+00
             1.0248E+00
 PARAMETER:  1.6638E-01 -8.5501E-02  4.5584E-01  3.4640E-01  2.4380E-01  4.0269E-01  9.6583E-01  4.0008E-01  1.5720E-01  3.9206E-01
             1.2448E-01
 GRADIENT:   4.5152E+00  1.1393E+01 -1.9499E+00  1.8892E+01  2.5677E+00 -3.9341E+00  7.6533E-01 -3.2223E+00  3.1757E+00  3.4347E-01
             1.3640E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2019.62118753239        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      834
 NPARAMETR:  1.0691E+00  7.3794E-01  1.4389E+00  1.2800E+00  1.1575E+00  1.4108E+00  2.3736E+00  1.3558E+00  1.0571E+00  1.3438E+00
             1.0253E+00
 PARAMETER:  1.6684E-01 -2.0390E-01  4.6389E-01  3.4683E-01  2.4630E-01  4.4413E-01  9.6441E-01  4.0436E-01  1.5557E-01  3.9549E-01
             1.2498E-01
 GRADIENT:   5.9007E+00 -1.2855E+01 -7.2761E+00 -2.2880E+01  1.4374E+01  1.4887E+01 -8.0461E+00 -2.2925E+00  4.8807E+00 -2.0176E-01
             3.2946E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2020.14803354915        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1025             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0688E+00  7.8908E-01  1.4441E+00  1.2787E+00  1.1534E+00  1.3929E+00  2.3798E+00  1.3658E+00  1.0456E+00  1.3425E+00
             1.0241E+00
 PARAMETER:  1.6657E-01 -1.3689E-01  4.6746E-01  3.4584E-01  2.4275E-01  4.3141E-01  9.6703E-01  4.1173E-01  1.4454E-01  3.9452E-01
             1.2381E-01
 GRADIENT:   6.9104E+02  3.3952E+01  6.5558E+00  3.8951E+02  1.5763E+01  4.5211E+02  1.8847E+02 -1.2862E+00  1.6980E+01  5.5959E+00
             2.6117E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2020.16916936869        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1212
 NPARAMETR:  1.0691E+00  7.8462E-01  1.4449E+00  1.2779E+00  1.1526E+00  1.3941E+00  2.3960E+00  1.3640E+00  1.0413E+00  1.3385E+00
             1.0240E+00
 PARAMETER:  1.6678E-01 -1.4256E-01  4.6801E-01  3.4519E-01  2.4198E-01  4.3225E-01  9.7381E-01  4.1045E-01  1.4048E-01  3.9157E-01
             1.2369E-01
 GRADIENT:   5.4185E+00  4.6492E-01 -7.0147E-01 -4.8182E+00  3.0433E+00  9.5879E+00 -2.4907E+00 -2.9848E+00  2.0089E+00  1.5521E-01
             1.2058E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2020.17493121534        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1271
 NPARAMETR:  1.0691E+00  7.8458E-01  1.4452E+00  1.2783E+00  1.1520E+00  1.3929E+00  2.3992E+00  1.3640E+00  1.0399E+00  1.3385E+00
             1.0238E+00
 PARAMETER:  1.6678E-01 -1.4261E-01  4.6826E-01  3.4554E-01  2.4152E-01  4.3141E-01  9.7513E-01  4.1045E-01  1.3908E-01  3.9153E-01
             1.2352E-01
 GRADIENT:  -2.9546E+04  1.7278E+04  7.4545E+03  3.8540E+04 -5.5162E+04 -5.9432E-01  1.3656E+04  1.1923E+04 -3.5436E+04  1.3842E-01
            -3.9902E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1271
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.2003E-04  1.6679E-02 -5.1794E-02 -1.9866E-02 -4.1441E-02
 SE:             3.0016E-02  2.1377E-02  1.6573E-02  2.2480E-02  2.0567E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8618E-01  4.3526E-01  1.7766E-03  3.7685E-01  4.3909E-02

 ETASHRINKSD(%)  1.0000E-10  2.8383E+01  4.4479E+01  2.4689E+01  3.1099E+01
 ETASHRINKVR(%)  1.0000E-10  4.8711E+01  6.9174E+01  4.3282E+01  5.2527E+01
 EBVSHRINKSD(%)  2.1180E-01  2.9579E+01  5.0385E+01  2.3375E+01  2.6015E+01
 EBVSHRINKVR(%)  4.2315E-01  5.0409E+01  7.5384E+01  4.1286E+01  4.5263E+01
 RELATIVEINF(%)  9.9099E+01  7.1715E+00  7.5817E+00  8.4782E+00  1.8502E+01
 EPSSHRINKSD(%)  3.4844E+01
 EPSSHRINKVR(%)  5.7546E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2020.1749312153386     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1101.2363980106659     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2020.175       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  7.85E-01  1.45E+00  1.28E+00  1.15E+00  1.39E+00  2.40E+00  1.36E+00  1.04E+00  1.34E+00  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,       47.474
Stop Time:
Sun Oct 24 01:11:41 CDT 2021
