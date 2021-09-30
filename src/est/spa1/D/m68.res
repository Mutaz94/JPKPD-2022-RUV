Thu Sep 30 03:28:43 CDT 2021
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
$DATA ../../../../data/spa1/D/dat68.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16237.0313669737        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0006E+02  4.3257E+02 -5.6401E+01  1.3300E+02  4.9385E+02 -2.6554E+03 -8.4419E+02 -8.5972E+01 -1.8495E+03 -9.9280E+02
            -2.9974E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -643.051501586302        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.7402E+00  9.6231E-01  8.0113E-01  2.5366E+00  1.3091E+00  3.1026E+00  1.2863E+00  9.7740E-01  2.4470E+00  1.3560E+00
             1.2255E+01
 PARAMETER:  6.5400E-01  6.1585E-02 -1.2173E-01  1.0308E+00  3.6932E-01  1.2322E+00  3.5176E-01  7.7141E-02  9.9485E-01  4.0453E-01
             2.6059E+00
 GRADIENT:   5.9896E+01  5.3714E+01 -4.7743E+01  1.4082E+02 -1.4900E+00  7.4743E+01  3.4584E-01  9.1194E+00 -5.2886E+01  5.2416E+00
             8.4591E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -706.801988189940        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.5100E+00  9.5935E-01  4.4903E+00  2.6360E+00  9.6035E+00  4.2439E+00  2.4788E+00  2.8367E-01  4.1468E+00  9.6384E+00
             9.1688E+00
 PARAMETER:  5.1209E-01  5.8499E-02  1.6019E+00  1.0693E+00  2.3621E+00  1.5455E+00  1.0078E+00 -1.1600E+00  1.5223E+00  2.3658E+00
             2.3158E+00
 GRADIENT:   4.5920E+01  1.1192E+01  5.9317E+00  7.8062E+01  6.7911E-01  1.8323E+02  1.2542E+01 -1.8620E-02  4.6069E+01  1.7231E-01
             5.9846E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -725.741196758115        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.2326E+00  6.7473E-01  2.3292E+00  1.5469E+00  8.8890E+00  2.9188E+00  8.2753E+00  2.0454E-02  1.6616E+00  1.0058E+01
             1.0059E+01
 PARAMETER:  3.0916E-01 -2.9345E-01  9.4550E-01  5.3625E-01  2.2848E+00  1.1712E+00  2.2133E+00 -3.7896E+00  6.0778E-01  2.4084E+00
             2.4085E+00
 GRADIENT:   8.1882E+00  3.8687E+00 -1.9874E-02  1.3883E+01 -2.8266E-01  6.1399E+01  1.1315E+02  2.6102E-05  2.7375E+01  5.4192E-01
             1.0955E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -730.152498936019        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.2277E+00  6.7723E-01  3.0918E+00  1.5567E+00  3.5170E+02  2.9430E+00  8.1851E+00  1.0000E-02  1.6457E+00  6.5009E+00
             9.0714E+00
 PARAMETER:  3.0517E-01 -2.8975E-01  1.2287E+00  5.4257E-01  5.9628E+00  1.1794E+00  2.2023E+00 -4.8979E+00  5.9816E-01  1.9719E+00
             2.3051E+00
 GRADIENT:  -2.9241E+00  3.6850E+00 -1.8097E-01  2.0209E+01  4.4562E-03  2.1005E+01  2.0775E+00  0.0000E+00  1.8648E+01 -4.3487E-05
            -6.7292E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -733.328356021626        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      570
 NPARAMETR:  1.2275E+00  6.7729E-01  7.5333E+00  1.5515E+00  1.0713E+02  2.8824E+00  7.9249E+00  1.0000E-02  1.1988E+00  2.3426E+01
             9.0588E+00
 PARAMETER:  3.0498E-01 -2.8965E-01  2.1193E+00  5.3925E-01  4.7740E+00  1.1586E+00  2.1700E+00 -4.7784E+00  2.8136E-01  3.2538E+00
             2.3037E+00
 GRADIENT:  -6.3789E+00  3.5232E+00 -2.6809E-01  5.2670E+01  4.3326E-02  1.5834E+01 -9.6005E+00  0.0000E+00 -9.7886E-01 -6.1372E-03
            -2.5123E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -734.106761152692        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      753
 NPARAMETR:  1.2275E+00  6.7731E-01  8.9612E+00  1.5498E+00  5.4225E+00  2.8785E+00  7.9968E+00  1.0000E-02  1.2441E+00  5.1391E+00
             9.1348E+00
 PARAMETER:  3.0498E-01 -2.8963E-01  2.2929E+00  5.3815E-01  1.7906E+00  1.1573E+00  2.1790E+00 -4.7784E+00  3.1844E-01  1.7369E+00
             2.3121E+00
 GRADIENT:  -1.0495E+01  2.5643E+00  8.8919E-02  3.3846E+01  3.1580E-01  1.4961E+01 -4.5576E+00  0.0000E+00  4.8525E+00 -6.4948E-02
            -7.7816E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -735.096047491720        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      931
 NPARAMETR:  1.2279E+00  6.7721E-01  1.1324E+01  1.5449E+00  2.4463E+00  2.8646E+00  8.1060E+00  1.0000E-02  1.0208E+00  6.1814E-01
             9.2926E+00
 PARAMETER:  3.0527E-01 -2.8977E-01  2.5269E+00  5.3495E-01  9.9458E-01  1.1524E+00  2.1926E+00 -4.7784E+00  1.2062E-01 -3.8103E-01
             2.3292E+00
 GRADIENT:  -1.5692E+01  2.6813E+00 -8.1551E-02  5.0375E+01 -2.0887E-01  1.3517E+01 -4.3602E+00  0.0000E+00 -3.8878E+00  3.0636E-01
             9.9898E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -735.779970065654        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1109
 NPARAMETR:  1.2288E+00  6.7718E-01  1.5994E+01  1.5322E+00  2.4333E+00  2.8321E+00  8.3582E+00  1.0000E-02  1.0224E+00  3.2710E-02
             9.3253E+00
 PARAMETER:  3.0602E-01 -2.8981E-01  2.8722E+00  5.2669E-01  9.8925E-01  1.1410E+00  2.2232E+00 -4.7784E+00  1.2214E-01 -3.3201E+00
             2.3327E+00
 GRADIENT:  -1.6174E+01  2.9938E+00  1.8655E-02  3.8444E+01  5.1628E-01  9.6558E+00  3.0848E+00  0.0000E+00 -1.2467E+00  9.2902E-04
             9.2552E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -735.871672277997        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  1.2265E+00  6.7750E-01  1.6038E+01  1.5325E+00  2.4327E+00  2.8060E+00  8.3059E+00  1.0000E-02  1.0535E+00  1.5931E-02
             9.2614E+00
 PARAMETER:  3.0419E-01 -2.8935E-01  2.8749E+00  5.2690E-01  9.8902E-01  1.1318E+00  2.2170E+00 -4.7784E+00  1.5210E-01 -4.0395E+00
             2.3259E+00
 GRADIENT:  -1.6376E+01  2.9489E+00  1.7061E-02  3.6733E+01  6.2638E-01  6.4402E+00  2.3967E+00  0.0000E+00 -3.2846E-01  2.2309E-04
             4.1218E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -736.022250762132        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1451
 NPARAMETR:  1.2313E+00  6.5862E-01  1.5069E+01  1.5325E+00  2.4013E+00  2.8033E+00  8.3043E+00  1.0000E-02  1.0577E+00  1.3988E-02
             9.2514E+00
 PARAMETER:  3.0804E-01 -3.1760E-01  2.8126E+00  5.2691E-01  9.7564E-01  1.1311E+00  2.2168E+00 -4.7784E+00  1.5563E-01 -4.1283E+00
             2.3254E+00
 GRADIENT:  -2.3397E+01  2.4678E+00  2.5330E-02  3.5694E+01 -3.3152E+02  2.4258E+02  1.8604E+00  0.0000E+00 -1.6102E-01  2.0087E-04
             1.9793E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1451
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0286E-02  3.8840E-02  7.3408E-06 -9.1137E-02 -2.0802E-05
 SE:             2.8177E-02  2.1929E-02  1.0520E-06  1.3598E-02  1.6212E-05
 N:                     100         100         100         100         100

 P VAL.:         2.8244E-01  7.6533E-02  3.0161E-12  2.0638E-11  1.9945E-01

 ETASHRINKSD(%)  5.6040E+00  2.6535E+01  9.9996E+01  5.4445E+01  9.9946E+01
 ETASHRINKVR(%)  1.0894E+01  4.6029E+01  1.0000E+02  7.9247E+01  1.0000E+02
 EBVSHRINKSD(%)  4.0244E+00  2.1212E+01  9.9995E+01  5.4482E+01  9.9888E+01
 EBVSHRINKVR(%)  7.8868E+00  3.7924E+01  1.0000E+02  7.9281E+01  1.0000E+02
 RELATIVEINF(%)  7.7163E+01  3.7067E+01  3.0290E-08  4.0654E+00  3.9319E-05
 EPSSHRINKSD(%)  1.0794E+01
 EPSSHRINKVR(%)  2.0423E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -736.02225076213222     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       182.91628244254048     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -736.022       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.23E+00  6.59E-01  1.51E+01  1.53E+00  2.40E+00  2.80E+00  8.30E+00  1.00E-02  1.06E+00  1.46E-02  9.26E+00
 


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
+        1.21E+05
 
 TH 2
+       -1.86E+00  3.70E+05
 
 TH 3
+        5.02E+02  9.14E+02  4.51E+00
 
 TH 4
+        2.64E+04  4.80E+04 -2.37E+02  2.50E+04
 
 TH 5
+        7.92E-01 -2.55E+00 -3.45E-03 -1.28E+00  2.96E+03
 
 TH 6
+       -3.12E+01  1.71E+00 -5.79E+01 -6.33E+03  1.91E+00  1.04E+02
 
 TH 7
+       -1.66E+02  2.50E+00  3.71E-03  4.28E-01  1.08E-02  2.59E+00  2.84E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.97E-01  2.35E+05  2.21E-02 -4.55E+01  6.44E-01 -1.61E+04  1.25E+00  0.00E+00  2.44E+01
 
 TH10
+        9.66E-02 -6.36E-02  3.50E-04  1.85E-03 -3.00E-03  5.50E-03  1.78E-03  0.00E+00 -9.89E-02  9.38E-01
 
 TH11
+       -1.98E+03 -1.86E+03  8.35E-03 -9.48E+00 -1.66E+02  1.00E+02  2.00E+01  0.00E+00  3.97E+00 -1.74E-03  2.58E+01
 
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
 #CPUT: Total CPU Time in Seconds,       39.664
Stop Time:
Thu Sep 30 03:29:24 CDT 2021
