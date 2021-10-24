Sat Oct 23 22:00:14 CDT 2021
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
$DATA ../../../../data/SD3/A1/dat95.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1751.12494597021        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8857E+02 -3.1404E+01  2.2585E+01 -4.6751E+01  1.3955E+02  4.4243E+01 -2.8062E+01 -5.6397E+01 -7.6573E+01 -7.4789E+00
            -5.2900E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1871.74778170914        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1338E+00  9.1299E-01  9.7584E-01  1.0780E+00  8.7872E-01  1.0939E+00  1.0827E+00  9.2980E-01  1.3197E+00  8.3149E-01
             1.6752E+00
 PARAMETER:  2.2557E-01  8.9711E-03  7.5546E-02  1.7512E-01 -2.9291E-02  1.8972E-01  1.7941E-01  2.7209E-02  3.7737E-01 -8.4533E-02
             6.1595E-01
 GRADIENT:   5.3973E+02 -1.7496E+01 -7.4296E-01 -3.6103E+00  2.1806E+01  3.9255E+01  4.0730E+00  9.8363E+00  3.2811E+01  9.6593E+00
             5.4159E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1887.62049739509        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:      187
 NPARAMETR:  1.1200E+00  6.7580E-01  5.5136E-01  1.2717E+00  5.4719E-01  1.1453E+00  1.3389E+00  4.8650E-01  1.0374E+00  3.5883E-01
             1.5781E+00
 PARAMETER:  2.1337E-01 -2.9186E-01 -4.9536E-01  3.4037E-01 -5.0295E-01  2.3562E-01  3.9182E-01 -6.2053E-01  1.3669E-01 -9.2491E-01
             5.5621E-01
 GRADIENT:   1.8036E+02  2.1238E+01 -2.0062E+01  9.7684E+01  2.7843E+01  3.9714E+01 -2.1251E+01 -5.3389E+00 -1.6210E+01 -6.3283E+00
             9.8401E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1909.92433575116        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  1.0220E+00  4.6082E-01  5.9443E-01  1.3417E+00  5.1765E-01  9.5635E-01  2.0946E+00  2.5484E-01  9.9285E-01  5.7538E-01
             1.4751E+00
 PARAMETER:  1.2177E-01 -6.7476E-01 -4.2015E-01  3.9394E-01 -5.5846E-01  5.5374E-02  8.3936E-01 -1.2671E+00  9.2828E-02 -4.5273E-01
             4.8869E-01
 GRADIENT:   5.0449E+01  1.2773E+01 -1.1645E+01  3.0982E+01  9.7701E+00  2.1839E-01  3.9518E+00 -1.4694E+00 -5.9370E+00  2.0516E+00
            -3.3145E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1913.06610531426        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      541
 NPARAMETR:  9.9346E-01  2.5243E-01  6.8332E-01  1.4621E+00  5.2020E-01  9.5148E-01  2.7866E+00  2.3883E-01  9.6965E-01  6.6068E-01
             1.5467E+00
 PARAMETER:  9.3439E-02 -1.2766E+00 -2.8080E-01  4.7985E-01 -5.5353E-01  5.0266E-02  1.1248E+00 -1.3320E+00  6.9177E-02 -3.1449E-01
             5.3610E-01
 GRADIENT:  -3.4036E+00  6.3604E+00  1.7404E+01  1.7130E+01 -2.2577E+01  2.5000E+00  9.3231E-01 -9.9672E-01 -2.4100E-01  1.1004E+00
            -9.7787E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1914.41365104830        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  1.0012E+00  1.2996E-01  6.7884E-01  1.5321E+00  5.0310E-01  9.9835E-01  3.9357E+00  6.5911E-01  9.3711E-01  5.5059E-01
             1.5189E+00
 PARAMETER:  1.0123E-01 -1.9405E+00 -2.8736E-01  5.2667E-01 -5.8696E-01  9.8346E-02  1.4701E+00 -3.1686E-01  3.5045E-02 -4.9677E-01
             5.1801E-01
 GRADIENT:   2.3256E+01 -1.6777E-01 -2.2120E+01  3.6653E+01  2.7339E+01  2.0363E+01  7.5086E-01 -6.4458E-01 -1.0942E+00  2.6270E+00
            -1.9666E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1916.97775426992        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  9.8793E-01  1.3504E-01  6.4596E-01  1.4920E+00  4.7976E-01  9.4052E-01  3.6714E+00  8.3743E-01  9.4579E-01  4.4623E-01
             1.5096E+00
 PARAMETER:  8.7855E-02 -1.9022E+00 -3.3701E-01  5.0011E-01 -6.3448E-01  3.8674E-02  1.4006E+00 -7.7415E-02  4.4266E-02 -7.0692E-01
             5.1187E-01
 GRADIENT:  -6.5087E+00  6.7493E+00 -8.1274E+00 -1.2650E+01  1.5080E+01 -1.1504E+00  9.8467E+00  2.8304E+00 -3.1283E+00  4.3801E-01
            -9.4745E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1919.45018694466        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  9.8320E-01  4.3982E-02  6.3702E-01  1.5369E+00  4.6174E-01  9.4246E-01  5.8804E+00  7.9092E-01  9.2950E-01  4.4037E-01
             1.5094E+00
 PARAMETER:  8.3056E-02 -3.0240E+00 -3.5095E-01  5.2976E-01 -6.7275E-01  4.0735E-02  1.8716E+00 -1.3456E-01  2.6888E-02 -7.2014E-01
             5.1173E-01
 GRADIENT:  -1.1117E+01 -2.4063E+00  1.5745E+01  2.1818E+01 -2.3144E+01  2.3858E-01 -8.0008E+00  3.6957E+00 -8.3747E-01  3.0956E+00
             3.7689E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1921.40227078995        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1269             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9039E-01  1.7264E-02  5.7340E-01  1.5078E+00  4.2731E-01  9.4016E-01  8.7643E+00  8.2077E-01  9.3445E-01  3.5598E-01
             1.5006E+00
 PARAMETER:  9.0340E-02 -3.9592E+00 -4.5617E-01  5.1062E-01 -7.5024E-01  3.8297E-02  2.2707E+00 -9.7509E-02  3.2199E-02 -9.3288E-01
             5.0587E-01
 GRADIENT:   1.6809E+02  1.4105E+01  1.5714E+01  3.5164E+02  6.2290E+01  1.0310E+01  1.3055E+02  2.2048E-01  5.7878E+00  1.1425E+00
             5.8799E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1921.42955891850        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1439
 NPARAMETR:  9.8574E-01  1.7207E-02  5.7326E-01  1.5105E+00  4.2721E-01  9.4238E-01  8.7330E+00  8.1539E-01  9.4011E-01  3.5202E-01
             1.4986E+00
 PARAMETER:  8.5516E-02 -3.9601E+00 -4.5670E-01  5.1214E-01 -7.5007E-01  4.0820E-02  2.2683E+00 -1.0305E-01  3.8763E-02 -9.4469E-01
             5.0418E-01
 GRADIENT:  -1.1157E+00  2.1177E+01 -2.0641E+02 -1.7843E+02  1.1427E+02  2.5052E-01  3.2822E+01  1.5876E+00  7.9007E-01 -1.0147E+02
            -1.8990E+02
 NUMSIGDIG:         2.0         2.4         2.3         2.3         2.4         1.9         2.4         1.1         1.4         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1439
 NO. OF SIG. DIGITS IN FINAL EST.:  1.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2451E-03  1.6832E-02 -7.3283E-03 -7.4622E-03 -2.0207E-03
 SE:             2.9675E-02  8.8870E-03  2.0989E-02  2.9033E-02  1.2922E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6653E-01  5.8221E-02  7.2698E-01  7.9716E-01  8.7574E-01

 ETASHRINKSD(%)  5.8611E-01  7.0227E+01  2.9683E+01  2.7350E+00  5.6709E+01
 ETASHRINKVR(%)  1.1688E+00  9.1136E+01  5.0555E+01  5.3953E+00  8.1259E+01
 EBVSHRINKSD(%)  7.9945E-01  8.0673E+01  2.8591E+01  2.8149E+00  5.6264E+01
 EBVSHRINKVR(%)  1.5925E+00  9.6265E+01  4.9008E+01  5.5506E+00  8.0872E+01
 RELATIVEINF(%)  9.8193E+01  3.0463E+00  4.7932E+00  6.5351E+01  1.7585E+00
 EPSSHRINKSD(%)  3.1324E+01
 EPSSHRINKVR(%)  5.2836E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1921.4295589184999     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1002.4910257138272     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1921.430       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.72E-02  5.73E-01  1.51E+00  4.27E-01  9.43E-01  8.74E+00  8.16E-01  9.41E-01  3.52E-01  1.50E+00
 


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
 #CPUT: Total CPU Time in Seconds,      116.052
Stop Time:
Sat Oct 23 22:00:32 CDT 2021
