Sat Sep 25 05:53:30 CDT 2021
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
$DATA ../../../../data/int/D/dat49.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   24913.0551145071        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8035E+02  3.0988E+02  8.5468E+00  2.5948E+02  2.2539E+02 -2.2634E+03 -8.6667E+02 -8.6063E+01 -1.1585E+03 -7.4560E+02
            -5.2342E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -968.518776813102        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.9830E+00  1.6580E+00  9.9283E-01  1.4565E+00  9.5990E-01  5.7577E+00  6.2045E+00  1.0080E+00  2.8601E+00  2.0523E+00
             1.1699E+01
 PARAMETER:  7.8460E-01  6.0564E-01  9.2808E-02  4.7602E-01  5.9073E-02  1.8505E+00  1.9253E+00  1.0800E-01  1.1509E+00  8.1895E-01
             2.5595E+00
 GRADIENT:   2.7202E+01  4.1651E+00 -3.9671E+01  3.5519E+01 -2.9371E+01  1.4411E+02  7.9250E+01  3.5564E+00  7.6781E+01  5.0152E+01
             5.0950E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1052.32770810316        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.4424E+00  9.3287E-01  1.1768E+02  3.9964E+00  2.8115E+00  4.2806E+00  1.6892E+01  9.4813E-01  2.8673E+00  2.3026E+00
             1.1559E+01
 PARAMETER:  4.6631E-01  3.0507E-02  4.8680E+00  1.4854E+00  1.1337E+00  1.5541E+00  2.9268E+00  4.6738E-02  1.1534E+00  9.3404E-01
             2.5474E+00
 GRADIENT:   1.7969E+01  1.7315E+01 -2.4485E+00  8.5834E+01  1.6706E+01  1.3413E+02  3.9530E+01  7.9453E-03 -5.9684E+00  6.6330E+01
             4.6478E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1248.35381302806        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0040E+00  7.9243E-01  4.0766E+01  1.3982E+00  2.3442E+00  1.5216E+00  5.1062E+00  3.6958E+00  2.1109E+00  7.5430E-01
             9.7324E+00
 PARAMETER:  1.0402E-01 -1.3265E-01  3.8079E+00  4.3522E-01  9.5194E-01  5.1979E-01  1.7305E+00  1.4072E+00  8.4709E-01 -1.8197E-01
             2.3755E+00
 GRADIENT:  -3.3842E+01 -2.7198E+01 -1.8688E+00 -3.5930E+01  1.0423E+01 -4.3746E+01 -4.9656E+01  2.6922E-01  1.7927E+01  8.6222E+00
             3.6780E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1305.07413547068        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.0010E+00  4.2653E-01  9.3097E+01  1.7677E+00  2.3319E+00  1.7760E+00  7.0084E+00  4.2554E+00  2.0396E+00  5.0757E-01
             7.9466E+00
 PARAMETER:  1.0099E-01 -7.5206E-01  4.6336E+00  6.6969E-01  9.4670E-01  6.7435E-01  2.0471E+00  1.5482E+00  8.1276E-01 -5.7812E-01
             2.1727E+00
 GRADIENT:   1.8656E+00 -4.8808E+00 -1.1881E+00 -3.4882E+00 -2.8348E+00 -1.0091E-01 -5.3365E-01  8.4007E-02  5.5142E+00  2.0958E+00
             1.0688E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1307.02616201521        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      435
 NPARAMETR:  9.9824E-01  3.7991E-01  1.4101E+02  1.8981E+00  2.3774E+00  1.7991E+00  8.4417E+00  4.3169E+00  2.0515E+00  4.5879E-01
             7.9051E+00
 PARAMETER:  9.8241E-02 -8.6781E-01  5.0488E+00  7.4085E-01  9.6603E-01  6.8730E-01  2.2332E+00  1.5625E+00  8.1858E-01 -6.7916E-01
             2.1675E+00
 GRADIENT:   2.2173E-01  2.5201E+00 -9.9678E-01 -8.4117E-01 -4.9697E-01  2.2439E+00  7.5522E+00  2.1799E-02  2.1679E+00  1.4882E+00
            -1.0237E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1308.24163515132        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      610
 NPARAMETR:  9.9466E-01  2.2217E-01  4.1317E+02  1.9697E+00  2.4265E+00  1.7841E+00  9.3206E+00  6.2433E+00  2.0090E+00  3.0677E-01
             7.9537E+00
 PARAMETER:  9.4650E-02 -1.4043E+00  6.1239E+00  7.7787E-01  9.8645E-01  6.7893E-01  2.3322E+00  1.9315E+00  7.9764E-01 -1.0816E+00
             2.1736E+00
 GRADIENT:  -2.7484E+00 -9.4422E-01 -3.8521E-01  1.6666E-01  1.2213E+00 -2.7014E-01  9.2893E-01  8.8078E-03 -1.8156E-01  5.2885E-01
             3.7309E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1308.92948671921        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  1.0002E+00  2.9005E-01  3.3102E+03  1.9393E+00  2.4533E+00  1.7853E+00  8.6798E+00  9.6559E+00  2.0259E+00  2.6569E-02
             7.9654E+00
 PARAMETER:  1.0025E-01 -1.1377E+00  8.2048E+00  7.6234E-01  9.9745E-01  6.7961E-01  2.2610E+00  2.3676E+00  8.0599E-01 -3.5280E+00
             2.1751E+00
 GRADIENT:   4.6205E-01 -3.5582E-02 -4.4824E-02 -2.3641E-01  1.1234E+00 -1.4191E-01  6.9522E-02  2.1289E-04 -3.7597E-01  2.9546E-03
             8.9870E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1308.97305377195        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      966
 NPARAMETR:  9.9905E-01  2.9325E-01  5.1929E+04  1.9388E+00  2.4499E+00  1.7864E+00  8.6612E+00  1.7746E+01  2.0325E+00  1.0000E-02
             7.9605E+00
 PARAMETER:  9.9048E-02 -1.1267E+00  1.0958E+01  7.6206E-01  9.9604E-01  6.8021E-01  2.2589E+00  2.9761E+00  8.0925E-01 -6.3375E+00
             2.1745E+00
 GRADIENT:  -2.7077E+00  1.0375E-01 -2.9777E-03  1.3619E+00 -2.0945E-01  4.6012E-01 -8.9145E-02  2.2898E-04 -5.1337E-01  0.0000E+00
             1.9854E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1308.97407980234        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1155             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9943E-01  2.9221E-01  6.3375E+04  1.9389E+00  2.4504E+00  1.7868E+00  8.6757E+00  1.8441E+01  2.0286E+00  1.0000E-02
             7.9610E+00
 PARAMETER:  9.9433E-02 -1.1303E+00  1.1157E+01  7.6212E-01  9.9627E-01  6.8042E-01  2.2605E+00  3.0146E+00  8.0737E-01 -6.5209E+00
             2.1746E+00
 GRADIENT:   6.8575E-01  7.6921E-01 -2.2069E-03  4.2971E+00  4.0452E-01  1.4640E+00  2.6066E+01 -8.5120E-06  1.1632E+00  0.0000E+00
             3.8710E+00

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1308.97419734391        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  9.9944E-01  2.9208E-01  6.3180E+04  1.9391E+00  2.4510E+00  1.7856E+00  8.6926E+00  1.8433E+01  2.0305E+00  1.0000E-02
             7.9610E+00
 PARAMETER:  9.9433E-02 -1.1308E+00  1.1157E+01  7.6212E-01  9.9630E-01  6.8035E-01  2.2602E+00  3.0146E+00  8.0898E-01 -6.5209E+00
             2.1746E+00
 GRADIENT:  -1.0047E-02 -6.9020E-03  3.1633E-03 -2.5977E-02 -5.4656E-02  1.6422E-01 -2.4091E-01  3.9381E-03  1.4580E-01  0.0000E+00
             2.7977E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1254
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0563E-02  6.6313E-02  1.2970E-05 -8.5719E-02 -4.5368E-05
 SE:             2.8146E-02  2.0342E-02  8.3598E-06  1.8649E-02  1.5282E-04
 N:                     100         100         100         100         100

 P VAL.:         7.0744E-01  1.1148E-03  1.2078E-01  4.3011E-06  7.6657E-01

 ETASHRINKSD(%)  5.7078E+00  3.1851E+01  9.9972E+01  3.7524E+01  9.9488E+01
 ETASHRINKVR(%)  1.1090E+01  5.3557E+01  1.0000E+02  6.0968E+01  9.9997E+01
 EBVSHRINKSD(%)  6.7917E+00  3.9768E+01  9.9960E+01  2.2895E+01  9.9467E+01
 EBVSHRINKVR(%)  1.3122E+01  6.3721E+01  1.0000E+02  4.0548E+01  9.9997E+01
 RELATIVEINF(%)  8.6772E+01  2.1436E+01  3.6264E-06  3.4687E+01  6.4842E-04
 EPSSHRINKSD(%)  7.7261E+00
 EPSSHRINKVR(%)  1.4855E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1308.9741973439100     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       345.11516242450080     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    42.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    19.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1308.974       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  2.92E-01  6.34E+04  1.94E+00  2.45E+00  1.79E+00  8.67E+00  1.84E+01  2.03E+00  1.00E-02  7.96E+00
 


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
+        2.04E+03
 
 TH 2
+        8.31E+01  1.83E+02
 
 TH 3
+        1.37E-04  1.03E-04  4.43E-11
 
 TH 4
+        8.15E+01 -5.37E+00  6.24E-06  7.38E+01
 
 TH 5
+       -1.96E+01  2.63E+01 -1.21E-05 -9.11E+00  5.26E+01
 
 TH 6
+       -3.46E+01 -2.43E+01  3.79E-05  8.89E+00  6.59E+00  6.40E+01
 
 TH 7
+        5.52E+00  1.05E+01  9.74E-07 -4.20E+00  2.19E-01 -3.99E-01  1.58E+00
 
 TH 8
+       -6.59E-01  4.11E-01 -5.92E-07 -5.09E-02  1.74E-01 -1.82E-01 -5.78E-03 -1.43E-03
 
 TH 9
+       -4.77E+01 -4.41E+01 -5.15E-08 -1.39E+00  1.12E+01 -5.34E+00 -2.79E-01  2.80E-02  2.72E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.19E+01 -2.09E+00 -1.88E-07 -4.63E+00  4.80E-02  9.74E-01  1.28E-01  3.85E-03  2.05E+00  0.00E+00  1.68E+01
 
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
 #CPUT: Total CPU Time in Seconds,       61.748
Stop Time:
Sat Sep 25 05:54:34 CDT 2021
