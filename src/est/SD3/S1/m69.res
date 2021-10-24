Sat Oct 23 23:09:39 CDT 2021
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
$DATA ../../../../data/SD3/S1/dat69.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2091.68104399462        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3266E+02  6.3786E+01 -5.1534E+01  1.8138E+02  1.3537E+02  4.3776E+01 -4.1124E-02  7.5702E+00 -2.0659E+00 -1.4435E+01
            -4.5638E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2101.81156009673        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0315E+00  9.7246E-01  1.0148E+00  9.6955E-01  9.2072E-01  1.0228E+00  1.0051E+00  9.8040E-01  1.0263E+00  1.0035E+00
             1.0637E+00
 PARAMETER:  1.3099E-01  7.2078E-02  1.1471E-01  6.9081E-02  1.7398E-02  1.2258E-01  1.0507E-01  8.0209E-02  1.2593E-01  1.0351E-01
             1.6171E-01
 GRADIENT:   4.4651E+00  5.0296E+00  9.5508E-01 -2.8523E+00  7.9378E+00 -7.4331E-01 -7.4333E-02  6.9339E+00 -2.2875E+00  3.4113E+00
             9.2138E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2105.42500205020        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      355
 NPARAMETR:  1.0102E+00  6.6124E-01  8.7664E-01  1.1506E+00  7.4504E-01  1.0488E+00  1.2291E+00  5.1901E-01  8.8672E-01  9.4911E-01
             1.0295E+00
 PARAMETER:  1.1019E-01 -3.1364E-01 -3.1662E-02  2.4025E-01 -1.9431E-01  1.4768E-01  3.0631E-01 -5.5584E-01 -2.0224E-02  4.7773E-02
             1.2907E-01
 GRADIENT:  -3.2706E+01  3.7075E+00 -1.1043E+01  2.4885E+01  2.6970E+00  9.5303E+00 -5.8077E+00  1.8129E+00 -8.9996E+00  6.7709E+00
            -9.6383E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2107.51736459759        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  1.0295E+00  7.2556E-01  7.9067E-01  1.0863E+00  7.2715E-01  1.0236E+00  1.4113E+00  3.2582E-01  9.1418E-01  8.5442E-01
             1.0435E+00
 PARAMETER:  1.2907E-01 -2.2082E-01 -1.3487E-01  1.8282E-01 -2.1862E-01  1.2335E-01  4.4450E-01 -1.0214E+00  1.0274E-02 -5.7334E-02
             1.4259E-01
 GRADIENT:   3.2675E+00 -2.1539E+00  8.3080E-01 -5.3921E+00 -4.0120E-01  7.8322E-02  1.0341E+00  4.1169E-01  1.0205E+00  1.0528E+00
             6.8400E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2107.76329985456        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  1.0299E+00  7.8017E-01  7.4458E-01  1.0519E+00  7.2471E-01  1.0259E+00  1.3262E+00  1.2668E-01  9.3505E-01  8.4596E-01
             1.0421E+00
 PARAMETER:  1.2946E-01 -1.4825E-01 -1.9494E-01  1.5060E-01 -2.2198E-01  1.2554E-01  3.8231E-01 -1.9661E+00  3.2847E-02 -6.7282E-02
             1.4120E-01
 GRADIENT:   2.0376E+00 -1.7761E+00 -1.8809E+00 -7.8653E-01  2.8935E+00  5.4994E-01  2.3789E-01  6.7566E-02  4.9874E-01  7.4494E-01
            -7.4120E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2107.80130596340        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  1.0293E+00  7.9322E-01  7.3549E-01  1.0441E+00  7.2361E-01  1.0249E+00  1.3097E+00  2.3577E-02  9.3863E-01  8.4016E-01
             1.0423E+00
 PARAMETER:  1.2889E-01 -1.3166E-01 -2.0722E-01  1.4319E-01 -2.2350E-01  1.2460E-01  3.6981E-01 -3.6475E+00  3.6664E-02 -7.4167E-02
             1.4143E-01
 GRADIENT:   3.4609E-01 -5.1895E-02  1.4348E-01 -2.9834E-01 -8.7397E-02  8.6376E-02  2.4170E-02  2.2411E-03  6.1043E-02  5.6188E-04
            -3.0386E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2107.80164828273        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0291E+00  7.9292E-01  7.3491E-01  1.0444E+00  7.2321E-01  1.0247E+00  1.3100E+00  1.0000E-02  9.3825E-01  8.3977E-01
             1.0423E+00
 PARAMETER:  1.2871E-01 -1.3203E-01 -2.0801E-01  1.4341E-01 -2.2405E-01  1.2436E-01  3.7001E-01 -4.5733E+00  3.6260E-02 -7.4626E-02
             1.4147E-01
 GRADIENT:  -4.1204E-02 -4.6399E-03 -4.0974E-02  2.0297E-02  2.4524E-02 -9.2401E-03 -3.1426E-03  0.0000E+00 -4.4536E-03  8.5391E-03
             8.1468E-03

0ITERATION NO.:   33    OBJECTIVE VALUE:  -2107.80371440755        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     1151
 NPARAMETR:  1.0305E+00  7.9291E-01  7.3493E-01  1.0441E+00  7.2318E-01  1.0252E+00  1.3109E+00  1.0000E-02  9.3832E-01  8.3974E-01
             1.0423E+00
 PARAMETER:  1.3007E-01 -1.3204E-01 -2.0797E-01  1.4317E-01 -2.2410E-01  1.2488E-01  3.7073E-01 -4.5733E+00  3.6341E-02 -7.4669E-02
             1.4146E-01
 GRADIENT:   2.8909E+00 -9.0725E-02  1.7982E-01 -5.0755E-01 -2.3871E-01  1.9818E-01  7.1866E-02  0.0000E+00  2.8531E-02  3.6279E-02
            -1.5105E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1151
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7823E-04  1.3439E-03 -4.0703E-04 -2.6200E-03 -7.8072E-03
 SE:             2.9868E-02  1.9642E-02  1.9946E-04  2.6146E-02  2.4016E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9257E-01  9.4545E-01  4.1283E-02  9.2018E-01  7.4512E-01

 ETASHRINKSD(%)  1.0000E-10  3.4197E+01  9.9332E+01  1.2408E+01  1.9544E+01
 ETASHRINKVR(%)  1.0000E-10  5.6700E+01  9.9996E+01  2.3276E+01  3.5268E+01
 EBVSHRINKSD(%)  3.4613E-01  3.4424E+01  9.9384E+01  1.2368E+01  1.9089E+01
 EBVSHRINKVR(%)  6.9106E-01  5.6998E+01  9.9996E+01  2.3206E+01  3.4534E+01
 RELATIVEINF(%)  9.8888E+01  4.8794E+00  5.5086E-04  1.1367E+01  9.8266E+00
 EPSSHRINKSD(%)  3.3723E+01
 EPSSHRINKVR(%)  5.6074E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2107.8037144075506     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1188.8651812028779     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2107.804       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  7.93E-01  7.35E-01  1.04E+00  7.23E-01  1.03E+00  1.31E+00  1.00E-02  9.38E-01  8.40E-01  1.04E+00
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       92.459
Stop Time:
Sat Oct 23 23:09:54 CDT 2021
