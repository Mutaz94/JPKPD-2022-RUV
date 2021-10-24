Sat Oct 23 15:30:56 CDT 2021
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
$DATA ../../../../data/SD1/SL3/dat54.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      986
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

 TOT. NO. OF OBS RECS:      886
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1944.75997291885        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3562E+02  1.6593E+02  1.6444E+01  2.3363E+02  1.0632E+02  4.3331E+01 -8.1918E+01 -1.1180E+02 -8.8668E+01  1.8988E+01
            -3.5008E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2982.04623469751        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0543E+00  1.1241E+00  1.2996E+00  8.5698E-01  1.1950E+00  1.2367E+00  1.1249E+00  8.6680E-01  1.1802E+00  7.7385E-01
             1.9000E+00
 PARAMETER:  1.5286E-01  2.1696E-01  3.6209E-01 -5.4339E-02  2.7811E-01  3.1242E-01  2.1772E-01 -4.2951E-02  2.6565E-01 -1.5637E-01
             7.4183E-01
 GRADIENT:   2.3361E+02  1.5112E+01 -3.2581E+00 -5.4824E-02  7.2791E+01  1.1595E+02 -3.3901E+00  6.1669E-02  6.9674E+00 -2.8428E+01
            -2.7130E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3003.10778606807        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.0409E+00  1.1623E+00  1.3460E+00  8.7547E-01  1.2407E+00  1.1658E+00  1.1624E+00  1.1227E-01  1.1532E+00  9.8780E-01
             2.0701E+00
 PARAMETER:  1.4005E-01  2.5040E-01  3.9713E-01 -3.2996E-02  3.1567E-01  2.5339E-01  2.5048E-01 -2.0868E+00  2.4256E-01  8.7725E-02
             8.2760E-01
 GRADIENT:   1.5420E+02  3.4600E+01 -1.4223E+01  3.8380E+01  6.3664E+01  7.5765E+01  1.5821E+01  3.2643E-02  9.6851E+00  4.6358E+00
            -5.9264E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3003.36346943734        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      248
 NPARAMETR:  1.0312E+00  1.1417E+00  1.2946E+00  8.8128E-01  1.2034E+00  1.1333E+00  1.1575E+00  9.5395E-02  1.1351E+00  9.5463E-01
             2.0795E+00
 PARAMETER:  1.3068E-01  2.3253E-01  3.5822E-01 -2.6384E-02  2.8515E-01  2.2515E-01  2.4630E-01 -2.2497E+00  2.2671E-01  5.3566E-02
             8.3214E-01
 GRADIENT:  -4.8833E+01 -3.1058E+01 -1.5820E+01  6.3030E+00  2.5512E+01  1.3133E+01  5.5221E+00  7.3573E-03  1.5524E+00  2.2323E+00
            -6.5042E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3010.84666583164        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  1.0640E+00  1.5014E+00  1.9860E+00  6.7710E-01  1.5975E+00  1.0910E+00  8.5924E-01  1.1304E-01  1.3925E+00  1.2365E+00
             2.1275E+00
 PARAMETER:  1.6203E-01  5.0638E-01  7.8613E-01 -2.8994E-01  5.6842E-01  1.8713E-01 -5.1712E-02 -2.0800E+00  4.3111E-01  3.1227E-01
             8.5493E-01
 GRADIENT:   3.6125E+00 -9.8284E+00  1.0629E+00 -2.5298E+00  1.8837E+00  9.3435E-01  1.5576E+00  5.3522E-03  8.2036E-01 -2.7490E+00
            -8.0286E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3011.11222121506        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      601
 NPARAMETR:  1.0618E+00  1.6299E+00  1.9444E+00  5.9494E-01  1.6688E+00  1.0887E+00  7.5628E-01  9.6492E-02  1.5708E+00  1.3006E+00
             2.1363E+00
 PARAMETER:  1.5998E-01  5.8850E-01  7.6494E-01 -4.1930E-01  6.1213E-01  1.8501E-01 -1.7935E-01 -2.2383E+00  5.5157E-01  3.6284E-01
             8.5909E-01
 GRADIENT:  -7.4946E-01  3.1289E-01  1.1086E+00 -3.6324E-01 -2.3205E+00  6.2702E-02 -2.5353E-01  5.5558E-03  4.5712E-02 -3.6563E-01
             8.1469E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3011.11765847907        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      776
 NPARAMETR:  1.0626E+00  1.6454E+00  1.8500E+00  5.8412E-01  1.6684E+00  1.0887E+00  7.5453E-01  8.8174E-02  1.5870E+00  1.2989E+00
             2.1358E+00
 PARAMETER:  1.6076E-01  5.9796E-01  7.1516E-01 -4.3765E-01  6.1184E-01  1.8497E-01 -1.8166E-01 -2.3284E+00  5.6185E-01  3.6152E-01
             8.5883E-01
 GRADIENT:   6.6369E-01 -4.1267E-01 -8.0874E-02 -1.5664E-01  2.5344E-01  1.0067E-02 -2.0055E-02  4.8572E-03  8.1268E-02  1.2427E-01
             2.0336E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3011.13972816101        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      961             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0631E+00  1.6264E+00  1.8529E+00  5.9623E-01  1.6558E+00  1.0895E+00  7.6620E-01  1.1838E-02  1.5578E+00  1.2906E+00
             2.1355E+00
 PARAMETER:  1.6116E-01  5.8639E-01  7.1677E-01 -4.1713E-01  6.0430E-01  1.8571E-01 -1.6632E-01 -4.3365E+00  5.4326E-01  3.5509E-01
             8.5871E-01
 GRADIENT:   2.0428E+02  2.4531E+02  1.9297E+00  4.1875E+01  4.1056E+01  3.5245E+01  4.0686E+00  2.3191E-04  1.1836E+01  4.0975E+00
             1.5300E+01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -3011.14019043743        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  1.0630E+00  1.6264E+00  1.8536E+00  5.9687E-01  1.6553E+00  1.0895E+00  7.6482E-01  1.0000E-02  1.5577E+00  1.2904E+00
             2.1354E+00
 PARAMETER:  1.6114E-01  5.8637E-01  7.1711E-01 -4.1606E-01  6.0400E-01  1.8569E-01 -1.6811E-01 -4.9340E+00  5.4321E-01  3.5497E-01
             8.5866E-01
 GRADIENT:   1.4230E+00 -8.1780E-01  5.9796E-03 -1.0030E-01 -1.7168E-01  3.0201E-01 -4.4916E-02  0.0000E+00  3.8276E-02 -7.7238E-03
            -1.1175E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1091
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1212E-04 -3.2875E-02 -9.1408E-05  2.1638E-02 -1.9560E-02
 SE:             2.9651E-02  2.0176E-02  5.5714E-05  2.3428E-02  2.6413E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8084E-01  1.0323E-01  1.0086E-01  3.5570E-01  4.5896E-01

 ETASHRINKSD(%)  6.6706E-01  3.2407E+01  9.9813E+01  2.1513E+01  1.1515E+01
 ETASHRINKVR(%)  1.3297E+00  5.4312E+01  1.0000E+02  3.8398E+01  2.1703E+01
 EBVSHRINKSD(%)  8.7918E-01  3.1869E+01  9.9826E+01  2.3184E+01  1.0543E+01
 EBVSHRINKVR(%)  1.7506E+00  5.3582E+01  1.0000E+02  4.0994E+01  1.9974E+01
 RELATIVEINF(%)  9.8214E+01  7.3870E+00  1.0090E-04  9.9262E+00  2.8545E+01
 EPSSHRINKSD(%)  1.7042E+01
 EPSSHRINKVR(%)  3.1180E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          886
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1628.3590808386800     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3011.1401904374302     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1382.7811095987502     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3011.140       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.63E+00  1.85E+00  5.97E-01  1.66E+00  1.09E+00  7.65E-01  1.00E-02  1.56E+00  1.29E+00  2.14E+00
 


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
 #CPUT: Total CPU Time in Seconds,       69.187
Stop Time:
Sat Oct 23 15:31:08 CDT 2021
