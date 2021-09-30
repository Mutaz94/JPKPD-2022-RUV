Wed Sep 29 13:46:08 CDT 2021
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
$DATA ../../../../data/spa/A3/dat73.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m73.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -180.780962162326        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2386E+02  1.0353E+02  1.3341E+02  1.1164E+01  1.2863E+02  7.6001E+01 -4.0840E+01 -5.2581E+01 -1.0826E+02 -1.3433E+02
            -2.6302E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1234.79494473278        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0674E+00  9.8727E-01  9.4000E-01  1.0562E+00  9.6758E-01  8.0803E-01  9.5172E-01  1.0051E+00  1.0075E+00  9.6460E-01
             2.4906E+00
 PARAMETER:  1.6522E-01  8.7193E-02  3.8123E-02  1.5472E-01  6.7040E-02 -1.1316E-01  5.0520E-02  1.0509E-01  1.0743E-01  6.3962E-02
             1.0125E+00
 GRADIENT:   2.7632E+02  8.7889E+00  1.8939E+01 -2.8953E-01  3.0739E+01 -1.5521E+01  4.6142E+00 -7.0112E+00 -1.5953E+01 -1.1136E+01
            -3.3408E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1296.25918520296        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0216E+00  1.0556E+00  5.0460E-01  9.5702E-01  6.8505E-01  1.0081E+00  2.3190E-01  1.2125E+00  1.0214E+00  1.9945E-02
             3.4572E+00
 PARAMETER:  1.2138E-01  1.5413E-01 -5.8399E-01  5.6065E-02 -2.7826E-01  1.0805E-01 -1.3614E+00  2.9270E-01  1.2113E-01 -3.8148E+00
             1.3405E+00
 GRADIENT:   2.7003E+01  5.5351E+01  3.1650E+01 -1.2869E+01 -6.6411E+01  6.0236E+01 -1.9694E+00 -8.9643E-02 -1.8328E+01 -3.2221E-03
            -5.1528E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1307.82311255659        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      274
 NPARAMETR:  1.0018E+00  8.5486E-01  5.8416E-01  1.0787E+00  6.8484E-01  8.6751E-01  7.4739E-01  5.6096E-01  8.4240E-01  1.5727E-02
             3.7585E+00
 PARAMETER:  1.0175E-01 -5.6812E-02 -4.3758E-01  1.7573E-01 -2.7857E-01 -4.2131E-02 -1.9116E-01 -4.7811E-01 -7.1497E-02 -4.0524E+00
             1.4240E+00
 GRADIENT:  -6.7242E+01 -1.4600E+01 -1.8992E+00 -3.2501E+01  2.4851E+00  1.9658E+01  4.6332E-01  1.6510E+00 -7.4088E+00  4.1258E-03
            -1.7264E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1312.04630562276        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  1.0291E+00  1.0057E+00  8.5222E-01  1.0441E+00  9.5349E-01  7.9798E-01  3.0226E-01  4.1649E-01  1.0076E+00  3.9321E-02
             3.8973E+00
 PARAMETER:  1.2864E-01  1.0570E-01 -5.9911E-02  1.4313E-01  5.2378E-02 -1.2568E-01 -1.0965E+00 -7.7588E-01  1.0759E-01 -3.1360E+00
             1.4603E+00
 GRADIENT:   5.9465E+00 -7.8789E+00 -2.0670E+00 -5.9340E+00  6.4499E+00  1.2170E-01  5.8927E-01  7.1136E-01 -3.7917E-01  1.3616E-02
             7.5208E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1312.64373073753        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      625
 NPARAMETR:  1.0307E+00  1.2164E+00  7.6303E-01  9.2060E-01  1.0052E+00  7.9956E-01  9.0484E-02  2.9992E-01  1.1810E+00  5.3058E-02
             3.8983E+00
 PARAMETER:  1.3026E-01  2.9593E-01 -1.7045E-01  1.7275E-02  1.0523E-01 -1.2369E-01 -2.3026E+00 -1.1042E+00  2.6633E-01 -2.8364E+00
             1.4605E+00
 GRADIENT:   2.1757E+00  1.4682E+01  2.3351E+00  1.1421E+01 -7.1351E+00  8.2543E-02  3.0807E-02  2.8254E-01 -1.2912E+00  1.3965E-02
             2.5395E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1312.89143754340        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  1.0284E+00  1.3297E+00  7.3973E-01  8.3588E-01  1.0677E+00  7.9744E-01  2.3870E-02  1.4669E-01  1.2987E+00  8.9210E-02
             3.8864E+00
 PARAMETER:  1.2799E-01  3.8498E-01 -2.0146E-01 -7.9273E-02  1.6550E-01 -1.2635E-01 -3.6351E+00 -1.8194E+00  3.6137E-01 -2.3168E+00
             1.4575E+00
 GRADIENT:  -8.1585E-02  5.9907E-01 -1.8698E-01  3.2605E-01  1.3991E-01  9.2088E-02  4.3893E-03  6.4360E-02  8.1641E-03  4.3332E-02
             6.7915E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1312.93344727987        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  1.0281E+00  1.2526E+00  7.5827E-01  8.8565E-01  1.0329E+00  7.9729E-01  3.6411E-02  5.3334E-02  1.2293E+00  7.3653E-02
             3.8825E+00
 PARAMETER:  1.2768E-01  3.2522E-01 -1.7671E-01 -2.1432E-02  1.3238E-01 -1.2653E-01 -3.2129E+00 -2.8312E+00  3.0646E-01 -2.5084E+00
             1.4565E+00
 GRADIENT:  -9.5424E-02 -4.5498E-01 -4.8630E-01  3.7241E-01  7.5806E-01 -1.1290E-02  9.8818E-03  9.4337E-03 -8.6213E-03  3.0352E-02
             2.5189E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1312.94314876940        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  1.0280E+00  1.2571E+00  7.6168E-01  8.8256E-01  1.0365E+00  7.9739E-01  1.5636E-02  1.0000E-02  1.2333E+00  7.3472E-02
             3.8836E+00
 PARAMETER:  1.2762E-01  3.2884E-01 -1.7223E-01 -2.4924E-02  1.3586E-01 -1.2641E-01 -4.0582E+00 -6.4735E+00  3.0973E-01 -2.5108E+00
             1.4568E+00
 GRADIENT:  -1.1168E-01  1.1541E-01  5.5222E-02  2.2210E-01 -1.2826E-01  1.0483E-01  1.9033E-03  0.0000E+00  3.3040E-02  2.9713E-02
             3.3394E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1312.95863410143        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1332
 NPARAMETR:  1.0280E+00  1.2735E+00  7.5809E-01  8.7174E-01  1.0446E+00  7.9702E-01  1.0000E-02  1.0000E-02  1.2475E+00  1.0046E-02
             3.8839E+00
 PARAMETER:  1.2764E-01  3.4173E-01 -1.7695E-01 -3.7269E-02  1.4363E-01 -1.2687E-01 -2.8770E+01 -1.1883E+02  3.2116E-01 -4.5006E+00
             1.4568E+00
 GRADIENT:   4.9976E-03 -5.4587E-02  8.2071E-03 -3.5184E-02  6.2283E-03  4.0879E-03  0.0000E+00  0.0000E+00  1.4523E-02  3.9115E-04
             2.6778E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1312.95864257502        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1465
 NPARAMETR:  1.0281E+00  1.2735E+00  7.5798E-01  8.7174E-01  1.0446E+00  7.9703E-01  1.0000E-02  1.0000E-02  1.2474E+00  1.0000E-02
             3.8832E+00
 PARAMETER:  1.2774E-01  3.4180E-01 -1.7710E-01 -3.7261E-02  1.4365E-01 -1.2687E-01 -2.9069E+01 -1.2020E+02  3.2105E-01 -4.5705E+00
             1.4567E+00
 GRADIENT:   3.1182E-01  6.3828E-03  3.3260E-03 -1.1934E-03  5.4588E-03 -6.6245E-03  0.0000E+00  0.0000E+00 -1.7812E-02  0.0000E+00
            -1.7716E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1465
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1686E-03 -4.8099E-04  5.1995E-05 -2.9892E-03  1.1084E-05
 SE:             2.8074E-02  2.4154E-04  8.2435E-05  2.2297E-02  1.6625E-04
 N:                     100         100         100         100         100

 P VAL.:         8.8196E-01  4.6440E-02  5.2821E-01  8.9335E-01  9.4684E-01

 ETASHRINKSD(%)  5.9502E+00  9.9191E+01  9.9724E+01  2.5302E+01  9.9443E+01
 ETASHRINKVR(%)  1.1546E+01  9.9993E+01  9.9999E+01  4.4203E+01  9.9997E+01
 EBVSHRINKSD(%)  6.2970E+00  9.9209E+01  9.9697E+01  2.4291E+01  9.9424E+01
 EBVSHRINKVR(%)  1.2198E+01  9.9994E+01  9.9999E+01  4.2682E+01  9.9997E+01
 RELATIVEINF(%)  8.3958E+01  2.4680E-04  1.1648E-04  4.1431E+00  2.4040E-04
 EPSSHRINKSD(%)  1.7804E+01
 EPSSHRINKVR(%)  3.2437E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1312.9586425750231     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -577.80781601128490     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1312.959       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.27E+00  7.58E-01  8.72E-01  1.04E+00  7.97E-01  1.00E-02  1.00E-02  1.25E+00  1.00E-02  3.88E+00
 


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
+        1.46E+03
 
 TH 2
+       -1.27E+02  4.40E+02
 
 TH 3
+        1.36E+01  1.41E+02  1.44E+02
 
 TH 4
+       -1.73E+02  3.86E+02  2.58E+01  5.32E+02
 
 TH 5
+        1.83E+01 -2.72E+02 -2.00E+02 -9.49E+01  3.34E+02
 
 TH 6
+       -1.80E+00 -2.34E+01  9.44E+00 -3.19E+01 -2.81E+00  2.40E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.65E+00 -3.20E+01  2.76E+00  7.68E+00  9.69E+00  7.58E+00  0.00E+00  0.00E+00  3.76E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.59E+00
 
 TH11
+       -2.05E+01 -2.10E+01 -4.89E+00 -1.19E+01  8.50E+00  4.90E+00  0.00E+00  0.00E+00  8.14E+00  0.00E+00  3.14E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.935
Stop Time:
Wed Sep 29 13:46:35 CDT 2021
