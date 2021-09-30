Wed Sep 29 20:12:45 CDT 2021
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
$DATA ../../../../data/spa/D/dat62.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14451.8795347326        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5489E+02  3.3873E+02 -9.6223E+01  1.0969E+02  3.7322E+02 -1.7652E+03 -6.9455E+02 -6.6840E+01 -1.5891E+03 -6.1520E+02
            -2.7218E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -629.532683717306        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4192E+00  1.0066E+00  9.5470E-01  1.5362E+00  1.2769E+00  1.9664E+00  1.1964E+00  9.7060E-01  1.3488E+00  1.0697E+00
             1.4437E+01
 PARAMETER:  4.5006E-01  1.0661E-01  5.3640E-02  5.2930E-01  3.4443E-01  7.7621E-01  2.7935E-01  7.0164E-02  3.9923E-01  1.6739E-01
             2.7698E+00
 GRADIENT:   1.7675E+01  1.1148E+01 -3.6452E+00  1.3223E+01 -6.1583E+00  5.1970E+01  5.5682E-01  4.4976E+00  7.9003E+00  2.4678E+00
             1.0987E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -637.379916183666        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.4098E+00  8.8128E-01  1.0521E+00  1.7113E+00  2.4433E+00  1.7282E+00  2.5887E+00  5.8073E-01  1.3963E+00  2.6375E+00
             1.3501E+01
 PARAMETER:  4.4346E-01 -2.6380E-02  1.5082E-01  6.3726E-01  9.9335E-01  6.4709E-01  1.0512E+00 -4.4347E-01  4.3383E-01  1.0698E+00
             2.7028E+00
 GRADIENT:   2.8305E+01  1.6528E+01 -1.1449E+01  4.6822E+01 -5.0180E+00  5.7862E-01  5.3964E+00  1.0325E+00  1.6175E+01  3.3746E+00
             6.8300E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -662.901893022601        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1004E+00  3.3247E-01  8.7421E-01  1.4075E+00  6.0580E+00  1.4123E+00  1.1414E+00  5.2551E-02  6.1387E-01  5.2672E+00
             1.3013E+01
 PARAMETER:  1.9569E-01 -1.0012E+00 -3.4434E-02  4.4183E-01  1.9014E+00  4.4525E-01  2.3226E-01 -2.8460E+00 -3.8798E-01  1.7615E+00
             2.6659E+00
 GRADIENT:  -3.9895E+01  2.9719E+00  2.2610E+01 -1.9490E+01 -2.9073E+00 -6.6621E+00  4.8173E-01  2.4099E-05  8.4595E+00  1.7436E-01
             7.7405E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -744.374289969285        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  6.8754E-01  1.7690E-02  6.5224E-02  6.0470E-01  2.5086E+02  1.4876E+00  1.8798E-01  1.0000E-02  4.1826E-02  2.2340E+00
             1.0080E+01
 PARAMETER: -2.7464E-01 -3.9347E+00 -2.6299E+00 -4.0303E-01  5.6249E+00  4.9713E-01 -1.5714E+00 -1.8145E+01 -3.0742E+00  9.0379E-01
             2.4105E+00
 GRADIENT:   5.8062E+01  4.8828E-01 -6.5645E+01  2.0244E+02  2.7494E-02 -5.0718E+01  7.3509E-05  0.0000E+00 -4.7694E-02 -2.8222E-05
            -1.1717E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -747.070575655878        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  5.4455E-01  1.0000E-02  3.2465E-02  3.9389E-01  7.5349E+02  1.5355E+00  1.1621E-01  1.0000E-02  1.5001E-02  1.8336E+00
             9.6450E+00
 PARAMETER: -5.0779E-01 -4.5252E+00 -3.3276E+00 -8.3169E-01  6.7247E+00  5.2885E-01 -2.0524E+00 -2.2957E+01 -4.0996E+00  7.0630E-01
             2.3664E+00
 GRADIENT:   8.4421E+01  2.5437E-02 -1.7067E+02  3.3352E+02  1.0548E-02 -4.3826E+01  1.3042E-05  0.0000E+00 -7.6241E-03 -1.1302E-06
            -1.4171E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -761.505951140961        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  5.6209E-01  1.3736E-02  3.1503E-02  3.8154E-01  6.3859E+02  1.6618E+00  7.0694E-02  1.0000E-02  1.3884E-02  1.7314E+00
             1.0569E+01
 PARAMETER: -4.7610E-01 -4.1877E+00 -3.3577E+00 -8.6354E-01  6.5593E+00  6.0792E-01 -2.5494E+00 -2.3488E+01 -4.1770E+00  6.4893E-01
             2.4579E+00
 GRADIENT:   5.5196E+01 -2.4448E-01 -1.9066E+02  2.4807E+02  1.1930E-02 -2.4825E+01  8.5182E-06  0.0000E+00 -1.5284E-03 -1.0176E-06
            -6.7675E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -785.934194175543        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      698
 NPARAMETR:  3.7907E-01  1.0000E-02  1.4179E-02  1.8213E-01  2.6495E+03  1.6465E+00  3.4471E-02  1.0000E-02  1.0000E-02  1.3008E+00
             1.0615E+01
 PARAMETER: -8.7004E-01 -5.0025E+00 -4.1560E+00 -1.6031E+00  7.9821E+00  5.9868E-01 -3.2676E+00 -3.0174E+01 -6.0433E+00  3.6295E-01
             2.4623E+00
 GRADIENT:  -5.8352E-01  0.0000E+00 -2.8482E+00  2.4090E+00 -2.5957E-04  3.9495E-01  1.0110E-04  0.0000E+00  0.0000E+00  1.1934E-08
             1.2104E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -785.940950152187        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  3.8033E-01  1.0000E-02  1.4159E-02  1.8168E-01  9.4538E+04  1.6426E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.3000E+00
             1.0611E+01
 PARAMETER: -8.6671E-01 -5.0025E+00 -4.1574E+00 -1.6055E+00  1.1557E+01  5.9627E-01 -4.7409E+00 -3.0174E+01 -6.0433E+00  3.6233E-01
             2.4619E+00
 GRADIENT:   2.6699E+00  0.0000E+00 -1.2763E+00 -1.1024E+00 -6.3640E-06 -1.7958E-01  0.0000E+00  0.0000E+00  0.0000E+00  3.1377E-11
            -8.1314E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -785.949047203309        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1029
 NPARAMETR:  3.7931E-01  1.0000E-02  1.4128E-02  1.8122E-01  1.3372E+21  1.6446E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.3237E+00
             1.0614E+01
 PARAMETER: -8.6940E-01 -5.0025E+00 -4.1596E+00 -1.6081E+00  4.8745E+01  5.9747E-01 -4.7409E+00 -3.0174E+01 -6.0433E+00  3.8044E-01
             2.4622E+00
 GRADIENT:   1.2962E+00  0.0000E+00  7.6491E-01 -2.9796E+00  0.0000E+00  2.4142E-01  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -8.8629E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1029
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6113E-04  3.2965E-06  6.7808E-05 -1.5938E-04  0.0000E+00
 SE:             2.9099E-02  4.3759E-06  3.3574E-04  3.8274E-04  0.0000E+00
 N:                     100         100         100         100         100

 P VAL.:         9.7365E-01  4.5126E-01  8.3994E-01  6.7709E-01  1.0000E+00

 ETASHRINKSD(%)  2.5139E+00  9.9985E+01  9.8875E+01  9.8718E+01  1.0000E+02
 ETASHRINKVR(%)  4.9647E+00  1.0000E+02  9.9987E+01  9.9984E+01  1.0000E+02
 EBVSHRINKSD(%)  2.7853E+00  9.9975E+01  9.8880E+01  9.8682E+01  1.0000E+02
 EBVSHRINKVR(%)  5.4930E+00  1.0000E+02  9.9987E+01  9.9983E+01  1.0000E+02
 RELATIVEINF(%)  6.4681E+00  3.6496E-06  4.8568E-05  6.6769E-05  0.0000E+00
 EPSSHRINKSD(%)  7.4710E+00
 EPSSHRINKVR(%)  1.4384E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -785.94904720330931     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -50.798220639571127     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           400
  
 #TERE:
 Elapsed estimation  time in seconds:    13.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -785.949       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.79E-01  1.00E-02  1.41E-02  1.81E-01  1.34E+21  1.64E+00  1.00E-02  1.00E-02  1.00E-02  1.32E+00  1.06E+01
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.65E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.03E+04  0.00E+00  4.47E+06
 
 TH 4
+       -7.97E+02  0.00E+00 -4.01E+05  3.90E+04
 
 TH 5
+        3.43E-28  0.00E+00 -2.51E-25  2.47E-27 -3.71E-51
 
 TH 6
+       -2.45E+01  0.00E+00  1.36E+03 -1.23E+02  2.23E-26  4.24E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.18E-05  0.00E+00 -2.34E-03 -3.57E-04  1.86E-28  7.30E-05  0.00E+00  0.00E+00  0.00E+00 -6.74E-05
 
 TH11
+       -2.17E+01  0.00E+00  4.48E+02 -2.45E+01 -3.22E-28  7.92E-01  0.00E+00  0.00E+00  0.00E+00  2.81E-05  3.69E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       20.532
Stop Time:
Wed Sep 29 20:13:07 CDT 2021
