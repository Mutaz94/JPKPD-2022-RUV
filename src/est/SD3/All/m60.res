Sun Oct 24 01:29:43 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/28/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/SD3/All/dat60.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA (0.09 FIX)x3
$SIGMA  1  FIX;        [P]

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
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1354.62360691321        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   5.0069E+02  9.1629E+01 -6.7815E+02  3.5660E+01  4.6209E+01 -3.3535E+02 -6.5288E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1402.89865945148        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  9.9154E-01  1.0837E+00  2.8985E+00  7.8557E-01  6.2219E-01  1.0751E+00  4.6846E+00
 PARAMETER:  9.1500E-02  1.8035E-01  1.1642E+00 -1.4134E-01 -3.7452E-01  1.7242E-01  1.6443E+00
 GRADIENT:   2.1508E+01 -1.5552E+02  1.5117E+01 -3.5448E+01 -1.8978E+01  1.0417E+00  1.8929E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1423.98607511726        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      113
 NPARAMETR:  9.6721E-01  1.1116E+00  4.8684E+00  9.1196E-01  7.5946E-01  3.1459E+00  3.9434E+00
 PARAMETER:  6.6661E-02  2.0582E-01  1.6828E+00  7.8379E-03 -1.7515E-01  1.2461E+00  1.4720E+00
 GRADIENT:   6.2878E+00 -1.6921E+01  5.3456E-01  2.6905E+00  3.8428E-01  1.5740E+00  3.0812E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1444.71492919994        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  9.7666E-01  1.0997E+00  1.7233E+00  9.0745E-01  7.3782E-01  5.9080E-01  3.7411E+00
 PARAMETER:  7.6381E-02  1.9503E-01  6.4426E-01  2.8877E-03 -2.0405E-01 -4.2628E-01  1.4194E+00
 GRADIENT:   8.2333E+00 -8.3017E+00  1.1973E+01 -7.2536E-01 -5.7424E+00  3.9473E+00  8.3275E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1446.90071563259        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      221
 NPARAMETR:  9.7287E-01  1.0962E+00  1.6138E+00  9.0692E-01  7.4824E-01  9.8628E-02  3.7308E+00
 PARAMETER:  7.2499E-02  1.9187E-01  5.7861E-01  2.2972E-03 -1.9004E-01 -2.2164E+00  1.4166E+00
 GRADIENT:   5.4795E-01 -1.8149E-01  2.2486E-01 -3.6461E-01 -1.5026E+00  1.1923E-01  5.6249E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1447.32941675285        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      337
 NPARAMETR:  9.8513E-01  1.1084E+00  1.6333E+00  9.1558E-01  7.5575E-01  3.7222E-02  3.7643E+00
 PARAMETER:  8.5017E-02  2.0295E-01  5.9060E-01  1.1802E-02 -1.8005E-01 -3.1909E+00  1.4256E+00
 GRADIENT:   5.2809E-01 -3.1994E-01  2.2051E-01  7.4124E-01 -1.1389E-01  1.5439E-02 -1.5238E-01

0ITERATION NO.:   29    OBJECTIVE VALUE:  -1447.33753364095        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      420
 NPARAMETR:  9.8496E-01  1.1084E+00  1.6322E+00  9.1365E-01  7.5585E-01  1.0000E-02  3.7652E+00
 PARAMETER:  8.4849E-02  2.0294E-01  5.8996E-01  9.6947E-03 -1.7992E-01 -4.7741E+00  1.4258E+00
 GRADIENT:   7.6537E-02 -1.0164E-01  6.0597E-02  6.0907E-02 -4.6223E-02  0.0000E+00  6.2236E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      420
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0866E-03 -8.8416E-03 -4.1895E-05
 SE:             2.8493E-02  2.6565E-02  8.8836E-05
 N:                     100         100         100

 P VAL.:         9.4162E-01  7.3927E-01  6.3721E-01

 ETASHRINKSD(%)  4.5440E+00  1.1003E+01  9.9702E+01
 ETASHRINKVR(%)  8.8816E+00  2.0795E+01  9.9999E+01
 EBVSHRINKSD(%)  4.3871E+00  1.1082E+01  9.9620E+01
 EBVSHRINKVR(%)  8.5818E+00  2.0936E+01  9.9999E+01
 RELATIVEINF(%)  9.0871E+01  7.4702E+01  1.3709E-03
 EPSSHRINKSD(%)  1.7151E+01
 EPSSHRINKVR(%)  3.1360E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1447.3375336409524     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -528.39900043627972     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     1.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1447.338       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         9.85E-01  1.11E+00  1.63E+00  9.14E-01  7.56E-01  1.00E-02  3.77E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,        7.943
Stop Time:
Sun Oct 24 01:29:47 CDT 2021
