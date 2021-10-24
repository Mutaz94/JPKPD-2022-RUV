Sun Oct 24 01:46:09 CDT 2021
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
$DATA ../../../../data/SD4/B/dat89.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1660.93887182181        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9777E+02 -2.8138E+01 -8.7129E+00 -1.4494E+00  6.7195E+01  6.0144E+01 -2.7616E-01 -8.6762E+00  1.0793E+01 -1.0332E+01
             2.0511E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1666.02272688209        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9766E-01  1.0184E+00  9.4965E-01  1.0275E+00  9.3343E-01  8.6496E-01  1.0058E+00  1.0585E+00  9.6444E-01  1.0321E+00
             9.2186E-01
 PARAMETER:  9.7661E-02  1.1826E-01  4.8339E-02  1.2715E-01  3.1111E-02 -4.5067E-02  1.0583E-01  1.5688E-01  6.3793E-02  1.3155E-01
             1.8635E-02
 GRADIENT:   1.0162E-02 -2.2012E+00  2.2007E+00 -7.1435E+00 -5.0843E+00 -2.1790E+01 -9.2921E-01 -2.6395E+00  1.3818E+00  5.8363E+00
            -1.2551E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1666.56814885124        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.9785E-01  1.0077E+00  9.1923E-01  1.0321E+00  9.0616E-01  8.7765E-01  1.0971E+00  1.1447E+00  9.2737E-01  9.5925E-01
             9.2532E-01
 PARAMETER:  9.7845E-02  1.0769E-01  1.5779E-02  1.3161E-01  1.4644E-03 -3.0504E-02  1.9271E-01  2.3514E-01  2.4593E-02  5.8399E-02
             2.2386E-02
 GRADIENT:   1.0184E-01 -5.0591E-01 -3.7311E-03 -8.0119E+00 -9.7453E+00 -1.5494E+01  8.3593E-01  3.2743E+00  9.4052E-01  3.8095E+00
            -1.0838E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1667.71733238159        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.9747E-01  8.2556E-01  9.2674E-01  1.1608E+00  8.3283E-01  9.1457E-01  1.3472E+00  9.7504E-01  8.2596E-01  8.7560E-01
             9.5333E-01
 PARAMETER:  9.7470E-02 -9.1690E-02  2.3916E-02  2.4912E-01 -8.2931E-02  1.0698E-02  3.9801E-01  7.4725E-02 -9.1209E-02 -3.2852E-02
             5.2204E-02
 GRADIENT:   6.6202E-01  1.7584E+01  4.1989E+00  2.3686E+01 -5.1434E+00  1.2294E+00 -3.5132E-01 -9.8159E-01 -1.9472E+00 -8.4958E-01
             3.4817E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1669.63140515230        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.9827E-01  5.7027E-01  6.4314E-01  1.2514E+00  5.9100E-01  9.0786E-01  1.8046E+00  5.8818E-01  7.5082E-01  6.2966E-01
             9.5958E-01
 PARAMETER:  9.8271E-02 -4.6164E-01 -3.4140E-01  3.2428E-01 -4.2594E-01  3.3395E-03  6.9035E-01 -4.3072E-01 -1.8659E-01 -3.6258E-01
             5.8738E-02
 GRADIENT:  -9.4700E-01  7.8728E+00  7.4660E+00  7.2951E+00 -1.2709E+01 -2.5439E+00 -1.4586E+00 -1.8595E+00 -2.0965E+00 -7.9337E-01
            -1.7846E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1670.94559162481        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.9171E-01  3.4564E-01  8.6655E-01  1.4157E+00  6.5753E-01  9.1315E-01  2.6141E+00  8.9556E-01  7.0052E-01  7.2934E-01
             9.5605E-01
 PARAMETER:  9.1676E-02 -9.6236E-01 -4.3232E-02  4.4765E-01 -3.1927E-01  9.1401E-03  1.0609E+00 -1.0306E-02 -2.5594E-01 -2.1561E-01
             5.5054E-02
 GRADIENT:   2.1213E-01  6.2461E+00 -3.3629E+00  1.7232E+01 -1.0814E+00  1.2393E+00  2.3682E+00  1.7491E+00 -9.1007E-01  2.2823E-01
            -5.8867E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1671.45803340990        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.8629E-01  1.9395E-01  1.2455E+00  1.5418E+00  7.8881E-01  9.0607E-01  3.4976E+00  1.2494E+00  6.7233E-01  8.8399E-01
             9.5041E-01
 PARAMETER:  8.6197E-02 -1.5402E+00  3.1950E-01  5.3297E-01 -1.3723E-01  1.3653E-03  1.3521E+00  3.2263E-01 -2.9701E-01 -2.3307E-02
             4.9143E-02
 GRADIENT:  -1.1436E+00  2.8512E+00 -2.2353E-01  1.3821E+01  1.9650E+00 -6.7747E-01  4.5770E-01 -1.4570E+00 -1.2248E-01 -3.0899E-01
            -1.1972E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1671.64075440682        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1237
 NPARAMETR:  9.8499E-01  1.1075E-01  1.4974E+00  1.6045E+00  8.4938E-01  9.0674E-01  4.4201E+00  1.5195E+00  6.6014E-01  9.3289E-01
             9.5203E-01
 PARAMETER:  8.4871E-02 -2.1004E+00  5.0375E-01  5.7279E-01 -6.3254E-02  2.1052E-03  1.5862E+00  5.1841E-01 -3.1531E-01  3.0529E-02
             5.0840E-02
 GRADIENT:   2.8669E-01  8.0923E-01  3.1185E+00  5.3841E+00 -3.7714E+00 -6.6714E-02 -2.7272E-01 -5.0517E-01  2.5684E-01 -4.3266E-01
             1.4645E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1672.00065054921        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  9.8407E-01  6.4382E-02  1.6454E+00  1.6348E+00  8.8335E-01  9.0607E-01  5.4697E+00  1.6887E+00  6.4607E-01  9.5452E-01
             9.5169E-01
 PARAMETER:  8.3947E-02 -2.6429E+00  5.9797E-01  5.9152E-01 -2.4031E-02  1.3646E-03  1.7992E+00  6.2394E-01 -3.3685E-01  5.3453E-02
             5.0485E-02
 GRADIENT:   7.7892E-01  1.5584E+00  1.3719E+00 -7.3273E+00 -2.9821E-02  1.9670E-02  3.4985E+00 -4.3277E-01 -3.3389E+00 -3.9099E-01
            -9.7251E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1672.34279742017        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1595
 NPARAMETR:  9.8332E-01  5.8499E-02  1.4070E+00  1.6178E+00  8.1276E-01  9.0489E-01  5.6559E+00  1.5103E+00  6.4720E-01  8.7910E-01
             9.5185E-01
 PARAMETER:  8.3181E-02 -2.7388E+00  4.4145E-01  5.8104E-01 -1.0731E-01  6.1557E-05  1.8327E+00  5.1229E-01 -3.3509E-01 -2.8858E-02
             5.0657E-02
 GRADIENT:   4.1106E-01  6.7123E+00 -4.5794E+00 -3.5423E+01  8.6197E+00 -2.1162E-01  1.6204E+01  2.2782E+00 -1.5932E+01 -8.3326E-01
            -5.6622E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1672.35070500612        NO. OF FUNC. EVALS.:  62
 CUMULATIVE NO. OF FUNC. EVALS.:     1657
 NPARAMETR:  9.8360E-01  5.8042E-02  1.4117E+00  1.6194E+00  8.1341E-01  9.0548E-01  5.6567E+00  1.5027E+00  6.4735E-01  8.8042E-01
             9.5205E-01
 PARAMETER:  8.3183E-02 -2.7435E+00  4.4362E-01  5.8154E-01 -1.0634E-01  2.4623E-04  1.8341E+00  5.1238E-01 -3.3535E-01 -2.7279E-02
             5.0778E-02
 GRADIENT:  -1.3024E+00  1.1464E+01 -1.1582E+00 -5.4344E+01  6.4665E-01 -3.0791E-01  1.1895E+01  1.8119E+00 -1.2041E+02  3.3493E-02
            -1.1469E-01
 NUMSIGDIG:         2.0         2.4         2.0         2.5         2.2         1.8         2.6         1.5         2.3         2.6
                    2.6

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1657
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6486E-03  3.1162E-02 -3.3882E-02 -2.2428E-02 -3.6977E-02
 SE:             2.9899E-02  1.3131E-02  1.9781E-02  2.7147E-02  1.8660E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5603E-01  1.7639E-02  8.6738E-02  4.0870E-01  4.7527E-02

 ETASHRINKSD(%)  1.0000E-10  5.6008E+01  3.3732E+01  9.0548E+00  3.7486E+01
 ETASHRINKVR(%)  1.0000E-10  8.0647E+01  5.6085E+01  1.7290E+01  6.0920E+01
 EBVSHRINKSD(%)  4.4773E-01  7.0419E+01  3.3311E+01  6.7566E+00  3.3755E+01
 EBVSHRINKVR(%)  8.9346E-01  9.1250E+01  5.5526E+01  1.3057E+01  5.6116E+01
 RELATIVEINF(%)  9.8920E+01  4.5937E+00  9.4255E+00  4.5587E+01  9.0895E+00
 EPSSHRINKSD(%)  4.5910E+01
 EPSSHRINKVR(%)  7.0742E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1672.3507050061201     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -937.19987844238187     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1672.351       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  5.82E-02  1.41E+00  1.62E+00  8.14E-01  9.05E-01  5.66E+00  1.51E+00  6.47E-01  8.80E-01  9.52E-01
 


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
 #CPUT: Total CPU Time in Seconds,       51.337
Stop Time:
Sun Oct 24 01:46:20 CDT 2021
