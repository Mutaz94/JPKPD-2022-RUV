Thu Sep 30 23:08:45 CDT 2021
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
$DATA ../../../../data/SD1/TD2/dat2.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3106.51796802373        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6481E+02  3.0105E+01  1.8645E+01  1.5469E+02  2.5208E+01  3.6376E+01 -2.9160E+01 -7.6533E+01 -1.0061E+02 -2.7409E+01
            -1.2425E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3383.53875137536        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.3392E-01  1.5755E+00  1.3926E+00  6.1932E-01  1.7505E+00  9.2749E-01  9.7816E-01  9.1390E-01  1.8044E+00  1.7226E+00
             1.6324E+00
 PARAMETER:  3.1635E-02  5.5456E-01  4.3119E-01 -3.7914E-01  6.5989E-01  2.4731E-02  7.7921E-02  9.9613E-03  6.9025E-01  6.4385E-01
             5.9003E-01
 GRADIENT:   2.8880E+01  2.0031E+02 -1.3419E+01  6.2711E+01  1.2826E+02 -1.7838E+01  5.2558E+01 -1.0643E+00  6.0395E+01  5.0581E+01
             2.2535E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3401.78004157245        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.6947E-01  1.6763E+00  1.6372E+00  5.7858E-01  1.5076E+00  8.8655E-01  5.5826E-01  8.3225E-01  1.5741E+00  1.8274E+00
             1.6408E+00
 PARAMETER:  6.8993E-02  6.1660E-01  5.9299E-01 -4.4718E-01  5.1054E-01 -2.0413E-02 -4.8293E-01 -8.3620E-02  5.5368E-01  7.0287E-01
             5.9519E-01
 GRADIENT:   1.1953E+02  3.1038E+02  3.1370E+01  3.5180E+01 -2.4379E+01 -3.3513E+01  1.5259E+01 -7.5380E+00  2.2976E+01  4.2581E+01
             2.1095E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3434.64376343395        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      270
 NPARAMETR:  9.8812E-01  1.8502E+00  1.9928E+00  5.0427E-01  1.8627E+00  9.5752E-01  5.0515E-01  3.0832E+00  1.6694E+00  1.7173E+00
             1.4998E+00
 PARAMETER:  8.8047E-02  7.1531E-01  7.8956E-01 -5.8465E-01  7.2205E-01  5.6588E-02 -5.8290E-01  1.2260E+00  6.1244E-01  6.4076E-01
             5.0534E-01
 GRADIENT:   1.2909E+01  2.0850E+01  5.4071E+00  7.0948E+00 -9.7023E+00 -1.5765E+01 -1.0426E+00 -2.4313E-01  7.6689E+00 -1.5258E+00
             5.8078E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3436.87337612317        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      411
 NPARAMETR:  9.8187E-01  1.8334E+00  2.0118E+00  5.0932E-01  1.8868E+00  9.9251E-01  5.3144E-01  3.6128E+00  1.5233E+00  1.7398E+00
             1.4557E+00
 PARAMETER:  8.1702E-02  7.0617E-01  7.9904E-01 -5.7468E-01  7.3486E-01  9.2480E-02 -5.3216E-01  1.3845E+00  5.2087E-01  6.5377E-01
             4.7549E-01
 GRADIENT:   2.0042E+02  5.8386E+02 -3.7308E+00  6.8094E+01  1.0180E+02  1.8126E+01  9.7784E+00  1.4866E+01  1.3245E+01  3.3374E+01
             9.0911E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3437.08309637006        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      596            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8259E-01  1.8306E+00  2.0318E+00  5.1587E-01  1.8893E+00  9.9590E-01  5.5334E-01  3.5170E+00  1.4908E+00  1.7375E+00
             1.4535E+00
 PARAMETER:  8.2441E-02  7.0467E-01  8.0890E-01 -5.6190E-01  7.3618E-01  9.5896E-02 -4.9177E-01  1.3576E+00  4.9929E-01  6.5243E-01
             4.7400E-01
 GRADIENT:   2.0207E+02  5.9565E+02  3.7768E+00  6.8612E+01  1.1011E+02  1.9290E+01  1.1208E+01  6.5481E+00  1.2643E+01  3.5548E+01
             8.8438E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3437.15835468933        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  9.8259E-01  1.8189E+00  2.0329E+00  5.2704E-01  1.8784E+00  9.9617E-01  5.5823E-01  3.4453E+00  1.4723E+00  1.7225E+00
             1.4533E+00
 PARAMETER:  8.2434E-02  6.9824E-01  8.0947E-01 -5.4048E-01  7.3044E-01  9.6161E-02 -4.8298E-01  1.3370E+00  4.8681E-01  6.4377E-01
             4.7385E-01
 GRADIENT:   6.8957E-01  8.2595E+00 -1.5104E+00  1.9139E+00  2.6509E-01  1.2712E-01  1.5061E-02  3.4563E-01  2.7264E-01 -5.1156E-01
            -3.7249E-01

0ITERATION NO.:   31    OBJECTIVE VALUE:  -3437.15835468933        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  9.8259E-01  1.8189E+00  2.0329E+00  5.2704E-01  1.8784E+00  9.9617E-01  5.5823E-01  3.4453E+00  1.4723E+00  1.7225E+00
             1.4533E+00
 PARAMETER:  8.2434E-02  6.9824E-01  8.0947E-01 -5.4048E-01  7.3044E-01  9.6161E-02 -4.8298E-01  1.3370E+00  4.8681E-01  6.4377E-01
             4.7385E-01
 GRADIENT:  -5.0518E-02  1.0402E+01 -2.0010E+01  1.0064E+00 -4.0811E+01  7.7735E-02 -4.8174E-02 -2.1361E+01 -6.2148E+01  4.6895E+01
            -6.6120E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      778
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1249E-04 -5.2673E-02 -2.5112E-02  3.4317E-02 -2.8470E-02
 SE:             2.9810E-02  1.9185E-02  1.5964E-02  2.4171E-02  2.6455E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7826E-01  6.0410E-03  1.1571E-01  1.5568E-01  2.8185E-01

 ETASHRINKSD(%)  1.3207E-01  3.5729E+01  4.6519E+01  1.9026E+01  1.1373E+01
 ETASHRINKVR(%)  2.6397E-01  5.8692E+01  7.1398E+01  3.4431E+01  2.1452E+01
 EBVSHRINKSD(%)  5.3161E-01  3.4602E+01  5.0721E+01  2.2917E+01  7.2497E+00
 EBVSHRINKVR(%)  1.0604E+00  5.7231E+01  7.5716E+01  4.0582E+01  1.3974E+01
 RELATIVEINF(%)  9.8930E+01  6.8543E+00  1.4476E+01  9.7337E+00  6.0256E+01
 EPSSHRINKSD(%)  1.9176E+01
 EPSSHRINKVR(%)  3.4675E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3437.1583546893294     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1783.0689949209186     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3437.158       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.82E+00  2.03E+00  5.27E-01  1.88E+00  9.96E-01  5.58E-01  3.45E+00  1.47E+00  1.72E+00  1.45E+00
 


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
 #CPUT: Total CPU Time in Seconds,       13.573
Stop Time:
Thu Sep 30 23:09:00 CDT 2021
