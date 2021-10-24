Sun Oct 24 03:19:33 CDT 2021
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
$DATA ../../../../data/SD4/SL2/dat61.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.48534126793        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6220E+02 -2.1206E+01  1.8079E+00 -1.9079E+01 -1.4691E+01  5.7248E+01  1.0386E+01  1.1880E+01  2.1143E+01  7.0750E+00
            -1.5378E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1705.19254070531        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0307E+00  1.0529E+00  1.0485E+00  1.0416E+00  1.0475E+00  9.8597E-01  9.4980E-01  9.2938E-01  9.0756E-01  9.9289E-01
             1.0589E+00
 PARAMETER:  1.3020E-01  1.5155E-01  1.4736E-01  1.4077E-01  1.4641E-01  8.5873E-02  4.8498E-02  2.6762E-02  3.0096E-03  9.2865E-02
             1.5727E-01
 GRADIENT:  -1.9953E+00  1.4272E+01  4.1129E+00  1.0985E+01 -8.3301E+00 -6.5055E-01  2.5512E+00  5.4008E+00 -1.6410E+00 -4.6945E+00
             2.6049E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1706.63744715930        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0365E+00  9.8298E-01  1.0076E+00  1.0778E+00  1.0137E+00  9.8421E-01  8.3813E-01  6.6146E-01  9.2487E-01  1.0508E+00
             1.0409E+00
 PARAMETER:  1.3584E-01  8.2833E-02  1.0759E-01  1.7497E-01  1.1361E-01  8.4085E-02 -7.6585E-02 -3.1331E-01  2.1899E-02  1.4957E-01
             1.4010E-01
 GRADIENT:   1.1603E+01 -1.0119E-01 -7.0661E+00  9.4385E+00  4.3991E+00 -1.1180E+00 -2.3129E+00  2.5847E+00  4.7526E-02 -1.4321E+00
            -3.3952E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1707.46234501206        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0290E+00  8.8050E-01  9.5722E-01  1.1377E+00  9.3378E-01  9.8694E-01  1.0847E+00  3.6523E-01  8.4111E-01  1.0069E+00
             1.0539E+00
 PARAMETER:  1.2861E-01 -2.7270E-02  5.6281E-02  2.2901E-01  3.1487E-02  8.6853E-02  1.8128E-01 -9.0724E-01 -7.3036E-02  1.0691E-01
             1.5254E-01
 GRADIENT:  -4.7550E+00  6.3028E+00  3.3246E+00  8.3122E+00 -6.4109E+00  2.8381E-01 -6.3184E-01  5.0581E-01 -9.1662E-01 -4.9421E-02
             1.6878E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1707.70580526046        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0313E+00  8.3651E-01  9.1663E-01  1.1559E+00  8.9670E-01  9.8617E-01  1.1772E+00  1.5865E-01  8.2030E-01  9.8308E-01
             1.0486E+00
 PARAMETER:  1.3086E-01 -7.8517E-02  1.2952E-02  2.4487E-01 -9.0290E-03  8.6070E-02  2.6316E-01 -1.7410E+00 -9.8089E-02  8.2935E-02
             1.4748E-01
 GRADIENT:   4.0313E-01  4.9028E-01 -4.9058E-01  1.3485E+00 -1.4132E-01 -6.9690E-03  1.6687E-01  1.0130E-01  5.9242E-02  1.7775E-01
            -7.0155E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1707.73109813836        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0318E+00  8.1633E-01  9.1960E-01  1.1667E+00  8.9018E-01  9.8604E-01  1.2001E+00  4.8227E-02  8.1270E-01  9.8258E-01
             1.0492E+00
 PARAMETER:  1.3131E-01 -1.0294E-01  1.6180E-02  2.5416E-01 -1.6337E-02  8.5945E-02  2.8236E-01 -2.9318E+00 -1.0739E-01  8.2423E-02
             1.4803E-01
 GRADIENT:   5.0326E+02  2.0093E+01  3.7791E+00  2.1098E+02  6.2893E+00  4.4255E+01  6.2756E+00  2.1581E-02  4.5077E+00  1.8799E-01
             7.7110E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1707.73427652393        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0310E+00  8.1783E-01  9.1829E-01  1.1666E+00  8.9055E-01  9.8572E-01  1.1984E+00  3.0359E-02  8.1383E-01  9.8578E-01
             1.0498E+00
 PARAMETER:  1.3049E-01 -1.0110E-01  1.4757E-02  2.5411E-01 -1.5914E-02  8.5620E-02  2.8097E-01 -3.3947E+00 -1.0601E-01  8.5683E-02
             1.4862E-01
 GRADIENT:  -9.0559E-02  9.7823E-02 -2.9126E-02  6.9065E-01 -4.6890E-02 -5.6702E-02  3.7581E-02  3.3336E-03  2.2317E-02  1.0592E-01
             7.0460E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1707.75180321475        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1239
 NPARAMETR:  1.0314E+00  8.3743E-01  9.1020E-01  1.1540E+00  8.9290E-01  9.8611E-01  1.1775E+00  1.0000E-02  8.2211E-01  9.8378E-01
             1.0496E+00
 PARAMETER:  1.3095E-01 -7.7418E-02  5.9138E-03  2.4321E-01 -1.3280E-02  8.6017E-02  2.6336E-01 -8.7777E+00 -9.5878E-02  8.3642E-02
             1.4840E-01
 GRADIENT:   4.7202E-01  2.8063E-01  6.5579E-01 -4.4107E-01 -1.7938E+00 -3.9121E-02  2.6371E-01  0.0000E+00  4.5360E-01  4.5031E-01
             2.1147E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1707.78250675505        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1416
 NPARAMETR:  1.0318E+00  8.8767E-01  8.9192E-01  1.1233E+00  9.0322E-01  9.8701E-01  1.1214E+00  1.0000E-02  8.3925E-01  9.7675E-01
             1.0487E+00
 PARAMETER:  1.3126E-01 -1.9157E-02 -1.4384E-02  2.1625E-01 -1.7873E-03  8.6927E-02  2.1460E-01 -2.1763E+01 -7.5246E-02  7.6477E-02
             1.4752E-01
 GRADIENT:  -6.3231E-02  7.4202E-01  1.5363E-01  6.9562E-01 -1.4931E+00 -9.8765E-03  1.0307E-01  0.0000E+00  1.2658E-01  3.6960E-02
             8.4743E-03

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1707.78494616331        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     1547
 NPARAMETR:  1.0332E+00  8.8711E-01  8.9205E-01  1.1228E+00  9.0417E-01  9.8742E-01  1.1206E+00  1.0000E-02  8.3875E-01  9.7668E-01
             1.0486E+00
 PARAMETER:  1.3270E-01 -1.9792E-02 -1.4232E-02  2.1580E-01 -7.3712E-04  8.7343E-02  2.1387E-01 -2.1763E+01 -7.5847E-02  7.6406E-02
             1.4748E-01
 GRADIENT:   3.2638E+00 -6.4606E-01 -7.0892E-01 -5.9554E-01  1.9197E-01  1.5916E-01 -1.5551E-02  0.0000E+00 -1.7138E-02 -9.2949E-02
            -1.4877E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1547
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.2821E-04 -4.6379E-03 -4.0179E-04 -2.9506E-03 -1.9985E-02
 SE:             2.9826E-02  1.8337E-02  1.8456E-04  2.4981E-02  2.3880E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8320E-01  8.0032E-01  2.9480E-02  9.0598E-01  4.0265E-01

 ETASHRINKSD(%)  7.9062E-02  3.8570E+01  9.9382E+01  1.6310E+01  1.9999E+01
 ETASHRINKVR(%)  1.5806E-01  6.2263E+01  9.9996E+01  2.9960E+01  3.5998E+01
 EBVSHRINKSD(%)  4.6863E-01  3.8697E+01  9.9415E+01  1.6429E+01  1.7916E+01
 EBVSHRINKVR(%)  9.3506E-01  6.2420E+01  9.9997E+01  3.0160E+01  3.2622E+01
 RELATIVEINF(%)  9.8274E+01  1.4077E+00  4.1758E-04  3.3235E+00  5.9310E+00
 EPSSHRINKSD(%)  4.2460E+01
 EPSSHRINKVR(%)  6.6891E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1707.7849461633084     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -972.63411959957023     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1707.785       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  8.87E-01  8.92E-01  1.12E+00  9.04E-01  9.87E-01  1.12E+00  1.00E-02  8.39E-01  9.77E-01  1.05E+00
 


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
 #CPUT: Total CPU Time in Seconds,       43.113
Stop Time:
Sun Oct 24 03:19:42 CDT 2021
