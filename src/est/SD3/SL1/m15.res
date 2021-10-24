Sat Oct 23 23:21:12 CDT 2021
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
$DATA ../../../../data/SD3/SL1/dat15.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2110.75968372550        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3950E+02 -6.0508E+01 -1.2011E+01 -4.0814E+01  6.1321E+01  6.4814E+01 -7.5608E+00 -1.7910E+01 -2.1494E+01 -6.6283E+00
            -2.8858E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2123.61498175186        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0478E+00  1.0867E+00  9.8152E-01  1.0346E+00  9.9340E-01  8.7543E-01  1.0305E+00  1.1686E+00  1.1161E+00  1.0019E+00
             1.0210E+00
 PARAMETER:  1.4672E-01  1.8318E-01  8.1348E-02  1.3406E-01  9.3377E-02 -3.3035E-02  1.3003E-01  2.5583E-01  2.0985E-01  1.0187E-01
             1.2082E-01
 GRADIENT:   1.3688E+01 -1.6880E+01 -1.9059E+01  1.4192E+01  2.1315E+01 -1.9728E+01  2.8207E+00 -8.1399E+00  5.3447E+00  1.0639E+00
            -1.0282E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2127.05651465993        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      385             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0625E+00  1.3064E+00  1.0674E+00  9.0659E-01  1.0900E+00  9.2643E-01  9.6342E-01  1.7527E+00  1.2534E+00  1.0058E+00
             1.0002E+00
 PARAMETER:  1.6063E-01  3.6725E-01  1.6520E-01  1.9298E-03  1.8620E-01  2.3586E-02  6.2735E-02  6.6117E-01  3.2587E-01  1.0583E-01
             1.0022E-01
 GRADIENT:   7.4695E+02  2.2983E+02 -1.9176E+00  6.5157E+01 -4.0875E+00  3.3880E+01  1.6515E+01  7.7880E+00  3.3701E+01 -2.5689E+00
            -2.1922E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2127.69079068625        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0625E+00  1.3071E+00  1.0674E+00  9.0658E-01  1.0900E+00  9.2399E-01  9.2109E-01  1.7527E+00  1.2534E+00  1.0139E+00
             1.0077E+00
 PARAMETER:  1.6063E-01  3.6781E-01  1.6520E-01  1.9282E-03  1.8620E-01  2.0948E-02  1.7804E-02  6.6117E-01  3.2587E-01  1.1384E-01
             1.0767E-01
 GRADIENT:   4.6959E+01  3.5045E+00 -4.1997E+00  9.8606E+00 -1.7339E+01  1.2725E+00  8.7254E+00  4.3050E+00  8.2071E+00 -2.7652E+00
            -1.6676E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2128.06884005890        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0625E+00  1.3062E+00  1.0674E+00  9.0658E-01  1.0900E+00  9.2316E-01  8.8801E-01  1.7527E+00  1.2534E+00  1.0227E+00
             1.0131E+00
 PARAMETER:  1.6063E-01  3.6712E-01  1.6520E-01  1.9241E-03  1.8619E-01  2.0047E-02 -1.8767E-02  6.6117E-01  3.2587E-01  1.2243E-01
             1.1300E-01
 GRADIENT:   4.6890E+01  2.7197E+00 -4.4643E+00  9.8018E+00 -1.8024E+01  9.6799E-01  6.4908E+00  4.0195E+00  6.9854E+00 -1.9747E+00
            -1.2198E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2128.44351400784        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  1.0625E+00  1.3061E+00  1.0671E+00  9.0658E-01  1.0900E+00  9.2321E-01  8.1440E-01  1.7506E+00  1.2303E+00  1.0234E+00
             1.0135E+00
 PARAMETER:  1.6063E-01  3.6702E-01  1.6490E-01  1.9234E-03  1.8619E-01  2.0099E-02 -1.0530E-01  6.5998E-01  3.0729E-01  1.2316E-01
             1.1343E-01
 GRADIENT:   4.6705E+01  5.6889E+00 -4.4822E+00  9.9676E+00 -1.7353E+01  1.0351E+00  7.8306E-01  2.6320E+00 -9.3541E-01 -3.3013E+00
            -1.2322E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2128.77186713010        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  1.0508E+00  1.2930E+00  1.0669E+00  9.0658E-01  1.0898E+00  9.2312E-01  7.9156E-01  1.7518E+00  1.2423E+00  1.0236E+00
             1.0135E+00
 PARAMETER:  1.4956E-01  3.5700E-01  1.6474E-01  1.9233E-03  1.8600E-01  1.9999E-02 -1.3375E-01  6.6064E-01  3.1695E-01  1.2328E-01
             1.1343E-01
 GRADIENT:   6.5276E+02  2.0880E+02 -4.3939E+00  5.6235E+01 -1.2666E+00  3.3480E+01  4.6598E+00  5.8219E+00  2.4248E+01 -2.7437E+00
            -1.0783E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2128.81436340637        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  1.0512E+00  1.2942E+00  1.0673E+00  9.0635E-01  1.0893E+00  9.1808E-01  7.8588E-01  1.7547E+00  1.2433E+00  1.0487E+00
             1.0138E+00
 PARAMETER:  1.4994E-01  3.5790E-01  1.6515E-01  1.6730E-03  1.8554E-01  1.4528E-02 -1.4095E-01  6.6229E-01  3.1775E-01  1.4753E-01
             1.1371E-01
 GRADIENT:   1.8974E+01 -4.7542E+00 -6.7499E+00  2.9923E+00 -1.7187E+01 -4.9048E-01  2.1016E-01  3.4065E+00  7.8156E-01 -9.7415E-02
            -1.1352E+01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2128.81436340637        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  1.0512E+00  1.2942E+00  1.0673E+00  9.0635E-01  1.0893E+00  9.1808E-01  7.8588E-01  1.7547E+00  1.2433E+00  1.0487E+00
             1.0138E+00
 PARAMETER:  1.4994E-01  3.5790E-01  1.6515E-01  1.6730E-03  1.8554E-01  1.4528E-02 -1.4095E-01  6.6229E-01  3.1775E-01  1.4753E-01
             1.1371E-01
 GRADIENT:   1.5484E+04  1.2947E+04  2.8098E+04 -4.6394E+04 -1.2517E+04 -2.3194E+04 -3.2914E+04  3.4676E+03  7.2955E+03  3.1444E+04
             4.0779E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1249
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.6696E-03 -3.3834E-02 -3.5141E-02  1.4353E-02 -3.9891E-02
 SE:             2.9901E-02  1.6392E-02  1.6670E-02  2.5960E-02  2.1054E-02
 N:                     100         100         100         100         100

 P VAL.:         8.4961E-01  3.9006E-02  3.5026E-02  5.8034E-01  5.8135E-02

 ETASHRINKSD(%)  1.0000E-10  4.5086E+01  4.4153E+01  1.3031E+01  2.9466E+01
 ETASHRINKVR(%)  1.0000E-10  6.9844E+01  6.8811E+01  2.4364E+01  5.0249E+01
 EBVSHRINKSD(%)  4.1097E-01  4.5190E+01  4.4995E+01  1.3172E+01  2.6918E+01
 EBVSHRINKVR(%)  8.2025E-01  6.9959E+01  6.9745E+01  2.4608E+01  4.6591E+01
 RELATIVEINF(%)  9.8846E+01  2.0073E+00  6.5783E+00  5.9922E+00  1.5562E+01
 EPSSHRINKSD(%)  3.4152E+01
 EPSSHRINKVR(%)  5.6640E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2128.8143634063736     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1209.8758302017009     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2128.814       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.29E+00  1.07E+00  9.06E-01  1.09E+00  9.18E-01  7.86E-01  1.75E+00  1.24E+00  1.05E+00  1.01E+00
 


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
 #CPUT: Total CPU Time in Seconds,       98.411
Stop Time:
Sat Oct 23 23:21:27 CDT 2021
