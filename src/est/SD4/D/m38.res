Sun Oct 24 04:18:17 CDT 2021
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
$DATA ../../../../data/SD4/D/dat38.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1627.55805205580        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5417E+02  5.3274E+00 -3.1916E+01  3.8731E+01 -1.2434E+01  2.1502E+01 -2.9635E+01  1.3476E+01 -5.0445E+01  2.0842E+01
             5.8588E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1640.10707047821        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9918E-01  1.1046E+00  1.2763E+00  1.0120E+00  1.1461E+00  1.2306E+00  1.2629E+00  9.1319E-01  1.4508E+00  8.8236E-01
             9.8000E-01
 PARAMETER:  9.9176E-02  1.9945E-01  3.4400E-01  1.1193E-01  2.3639E-01  3.0747E-01  3.3342E-01  9.1891E-03  4.7210E-01 -2.5150E-02
             7.9801E-02
 GRADIENT:   7.0067E+00  1.7967E+01  1.1102E+01  4.0732E+01 -1.1523E+01  4.2358E+01  2.9959E+00 -2.7377E+00  3.1424E+01 -1.0644E+01
            -1.0127E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1642.20072539008        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  1.0017E+00  1.2493E+00  1.1498E+00  9.3310E-01  1.1844E+00  1.1339E+00  1.3136E+00  7.2757E-01  1.5268E+00  9.5463E-01
             9.7439E-01
 PARAMETER:  1.0166E-01  3.2258E-01  2.3962E-01  3.0762E-02  2.6925E-01  2.2562E-01  3.7275E-01 -2.1804E-01  5.2316E-01  5.3567E-02
             7.4052E-02
 GRADIENT:   1.0524E+01  2.8546E+01  2.7813E+00  5.0279E+01  1.9971E+00  1.2877E+01  1.2497E+01 -6.1097E-01  2.9838E+01 -1.4419E+00
            -9.7469E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1647.86955898458        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.9780E-01  1.2866E+00  9.3530E-01  8.4710E-01  1.1186E+00  1.0917E+00  1.2232E+00  3.8587E-01  1.2790E+00  9.1566E-01
             9.9411E-01
 PARAMETER:  9.7801E-02  3.5197E-01  3.3114E-02 -6.5940E-02  2.1210E-01  1.8774E-01  3.0148E-01 -8.5224E-01  3.4607E-01  1.1887E-02
             9.4094E-02
 GRADIENT:   1.5072E+00  7.2367E+00  2.4639E-01  9.2906E+00 -8.5850E+00 -2.0382E+00 -6.6770E-01  3.5423E-01 -2.6046E+00  2.3860E+00
             1.4792E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1648.97282577347        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.9857E-01  1.7088E+00  7.4977E-01  5.7679E-01  1.2741E+00  1.0993E+00  9.9229E-01  2.2122E-01  1.7098E+00  9.5319E-01
             9.9469E-01
 PARAMETER:  9.8564E-02  6.3578E-01 -1.8799E-01 -4.5028E-01  3.4227E-01  1.9468E-01  9.2263E-02 -1.4086E+00  6.3639E-01  5.2061E-02
             9.4677E-02
 GRADIENT:   4.7317E-01  1.6584E+01  3.3186E+00  7.9473E+00 -4.8179E+00  6.2137E-01 -7.7491E-01  1.0210E-01 -1.9999E+00 -2.0261E+00
            -1.0409E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1649.00258249589        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  9.9870E-01  1.7785E+00  7.0705E-01  5.2974E-01  1.3011E+00  1.0991E+00  9.6566E-01  1.7946E-01  1.8150E+00  9.6613E-01
             9.9513E-01
 PARAMETER:  9.8703E-02  6.7575E-01 -2.4666E-01 -5.3537E-01  3.6322E-01  1.9449E-01  6.5056E-02 -1.6178E+00  6.9609E-01  6.5541E-02
             9.5119E-02
 GRADIENT:   4.0656E-01  1.4079E+01  2.0087E+00  7.3562E+00 -2.8157E+00  5.4915E-01 -4.3106E-01  7.8930E-02 -1.2577E+00 -1.7990E+00
            -9.0215E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1649.20697171625        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  9.9851E-01  1.7730E+00  7.0064E-01  5.2086E-01  1.3032E+00  1.0978E+00  9.6516E-01  8.7381E-02  1.8309E+00  9.7521E-01
             9.9464E-01
 PARAMETER:  9.8508E-02  6.7267E-01 -2.5577E-01 -5.5227E-01  3.6480E-01  1.9328E-01  6.4543E-02 -2.3375E+00  7.0479E-01  7.4893E-02
             9.4629E-02
 GRADIENT:   1.7031E-01 -1.4416E+00  9.2739E-01  1.7957E+00 -2.0996E+00  8.6957E-02  1.7209E-01  2.1305E-02 -4.6436E-01 -2.5057E-01
            -3.0628E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1649.24455778755        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1245
 NPARAMETR:  9.9944E-01  1.7673E+00  7.0085E-01  5.1800E-01  1.3036E+00  1.1002E+00  9.6388E-01  1.8527E-02  1.8376E+00  9.7516E-01
             9.9406E-01
 PARAMETER:  9.9441E-02  6.6946E-01 -2.5547E-01 -5.5777E-01  3.6510E-01  1.9547E-01  6.3214E-02 -3.8885E+00  7.0845E-01  7.4850E-02
             9.4042E-02
 GRADIENT:   1.9781E+00 -9.3649E+00  7.5392E-01 -1.1640E+00 -1.4499E+00  9.8116E-01  8.6240E-03  9.8916E-04 -1.4293E-01 -1.7951E-01
            -4.1667E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1649.24798234065        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1435             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9942E-01  1.7649E+00  6.9901E-01  5.1934E-01  1.3041E+00  1.1002E+00  9.6411E-01  1.0000E-02  1.8388E+00  9.7601E-01
             9.9473E-01
 PARAMETER:  9.9417E-02  6.6810E-01 -2.5808E-01 -5.5520E-01  3.6552E-01  1.9546E-01  6.3455E-02 -5.2754E+00  7.0909E-01  7.5721E-02
             9.4717E-02
 GRADIENT:   4.6623E+02  8.0130E+02  8.9888E-01  1.2341E+02  1.9139E+01  1.2932E+02  9.0897E+00  0.0000E+00  3.8933E+01  5.4699E-01
             9.0385E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1649.24875352136        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1577
 NPARAMETR:  9.9940E-01  1.7621E+00  6.9903E-01  5.2130E-01  1.3044E+00  1.1002E+00  9.6482E-01  1.0000E-02  1.8380E+00  9.7610E-01
             9.9497E-01
 PARAMETER:  9.9423E-02  6.6810E-01 -2.5550E-01 -5.5531E-01  3.6505E-01  1.9547E-01  6.3626E-02 -5.2754E+00  7.0808E-01  7.5644E-02
             9.4579E-02
 GRADIENT:   1.5123E-02  9.0249E-01  1.5109E-01 -4.8679E-01 -2.8155E-01  2.9176E-03 -3.5499E-02  0.0000E+00 -3.2417E-02 -7.7590E-03
            -5.4364E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1577
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.4244E-05 -3.1006E-02 -2.6053E-04  2.8946E-02 -4.2357E-02
 SE:             2.9870E-02  2.3981E-02  9.4446E-05  2.2517E-02  2.1517E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9828E-01  1.9603E-01  5.8070E-03  1.9863E-01  4.9003E-02

 ETASHRINKSD(%)  1.0000E-10  1.9660E+01  9.9684E+01  2.4564E+01  2.7917E+01
 ETASHRINKVR(%)  1.0000E-10  3.5455E+01  9.9999E+01  4.3094E+01  4.8040E+01
 EBVSHRINKSD(%)  3.6195E-01  1.8386E+01  9.9727E+01  2.6636E+01  2.6645E+01
 EBVSHRINKVR(%)  7.2259E-01  3.3392E+01  9.9999E+01  4.6177E+01  4.6191E+01
 RELATIVEINF(%)  9.9216E+01  7.1592E+00  1.4750E-04  6.0436E+00  1.7299E+01
 EPSSHRINKSD(%)  4.3688E+01
 EPSSHRINKVR(%)  6.8290E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1649.2487535213556     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -914.09792695761746     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1649.249       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  1.76E+00  7.01E-01  5.19E-01  1.30E+00  1.10E+00  9.64E-01  1.00E-02  1.84E+00  9.76E-01  9.95E-01
 


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
 #CPUT: Total CPU Time in Seconds,       49.471
Stop Time:
Sun Oct 24 04:18:27 CDT 2021
