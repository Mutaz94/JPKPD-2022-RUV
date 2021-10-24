Sun Oct 24 04:27:03 CDT 2021
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
$DATA ../../../../data/SD4/D/dat90.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m90.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1649.61785730599        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5333E+02 -2.5651E+00 -1.4025E+01  4.1687E+01  3.0718E+01  3.0367E+01 -7.9413E+00 -1.5970E-01  3.6497E+00 -9.4498E+00
             2.7369E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1651.14001358318        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  9.8235E-01  1.0204E+00  1.0071E+00  1.0099E+00  9.8857E-01  9.9290E-01  1.0089E+00  1.0005E+00  1.0017E+00  1.0074E+00
             9.8137E-01
 PARAMETER:  8.2197E-02  1.2018E-01  1.0707E-01  1.0985E-01  8.8502E-02  9.2870E-02  1.0883E-01  1.0055E-01  1.0172E-01  1.0735E-01
             8.1192E-02
 GRADIENT:  -4.5199E+00 -5.9435E+00 -2.7218E+00  1.7128E+00 -4.0807E+00 -1.1804E+01 -9.9278E+00 -1.3632E+00 -3.4500E+00 -7.8971E+00
             1.9528E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1652.51702427884        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  9.8814E-01  1.3522E+00  1.1381E+00  8.3703E-01  1.2223E+00  1.0761E+00  9.7145E-01  1.5091E+00  1.1697E+00  1.1642E+00
             9.0003E-01
 PARAMETER:  8.8065E-02  4.0170E-01  2.2935E-01 -7.7894E-02  3.0071E-01  1.7332E-01  7.1030E-02  5.1151E-01  2.5674E-01  2.5208E-01
            -5.3252E-03
 GRADIENT:   7.6545E+00  3.1151E+01 -3.9446E-01  3.4536E+01  1.8803E+01  1.8666E+01  2.8853E+00 -2.5772E+00  1.0425E+00 -6.6265E+00
            -1.5416E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1654.15197119975        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      487
 NPARAMETR:  9.8663E-01  1.3108E+00  1.0856E+00  8.3519E-01  1.1608E+00  1.0260E+00  9.5594E-01  1.4412E+00  1.1458E+00  1.1474E+00
             9.3382E-01
 PARAMETER:  8.6544E-02  3.7064E-01  1.8209E-01 -8.0091E-02  2.4907E-01  1.2568E-01  5.4937E-02  4.6545E-01  2.3612E-01  2.3750E-01
             3.1530E-02
 GRADIENT:   4.5791E+00  5.6051E+00  2.8824E+00  3.6590E+00 -7.2696E+00  8.4574E-01  3.6246E-01 -1.5899E-01 -4.0338E-01 -4.5605E-01
             1.4031E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1654.36395740045        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      664
 NPARAMETR:  9.8525E-01  1.5506E+00  9.5004E-01  6.8136E-01  1.2387E+00  1.0229E+00  8.6596E-01  1.5504E+00  1.3039E+00  1.1956E+00
             9.2519E-01
 PARAMETER:  8.5137E-02  5.3864E-01  4.8750E-02 -2.8366E-01  3.1408E-01  1.2269E-01 -4.3915E-02  5.3854E-01  3.6534E-01  2.7864E-01
             2.2248E-02
 GRADIENT:  -2.7858E-01  1.3361E+01  5.5719E+00  5.2709E+00 -4.2758E+00 -6.9923E-01 -7.5411E-01 -1.8959E+00 -1.4411E+00  7.0819E-02
            -3.3938E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1654.51279945547        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      840
 NPARAMETR:  9.8443E-01  1.7236E+00  7.7699E-01  5.6734E-01  1.2749E+00  1.0231E+00  8.2491E-01  1.5255E+00  1.4508E+00  1.2114E+00
             9.2380E-01
 PARAMETER:  8.4312E-02  6.4444E-01 -1.5233E-01 -4.6680E-01  3.4287E-01  1.2287E-01 -9.2475E-02  5.2231E-01  4.7209E-01  2.9178E-01
             2.0739E-02
 GRADIENT:  -3.3643E+00  1.6815E+01  5.6320E+00  5.6893E+00 -2.1299E+00 -7.8424E-01 -1.0569E+00 -2.2368E+00 -1.5429E+00  3.7420E-01
            -4.5837E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1654.52122522031        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  9.8465E-01  1.8215E+00  6.6689E-01  5.0656E-01  1.2846E+00  1.0241E+00  8.0756E-01  1.4641E+00  1.5436E+00  1.2094E+00
             9.2552E-01
 PARAMETER:  8.4529E-02  6.9966E-01 -3.0513E-01 -5.8012E-01  3.5046E-01  1.2386E-01 -1.1374E-01  4.8123E-01  5.3411E-01  2.9011E-01
             2.2605E-02
 GRADIENT:  -3.7025E+00  2.6676E+01  5.9953E+00  8.1172E+00 -4.2311E+00 -4.4031E-01 -1.1850E+00 -2.2214E+00 -1.6604E+00  1.6377E-01
            -4.2040E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1654.90455462181        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1209
 NPARAMETR:  9.8506E-01  1.7995E+00  6.4770E-01  5.0550E-01  1.2846E+00  1.0242E+00  8.0963E-01  1.4694E+00  1.5374E+00  1.2040E+00
             9.3248E-01
 PARAMETER:  8.4948E-02  6.8749E-01 -3.3432E-01 -5.8221E-01  3.5041E-01  1.2388E-01 -1.1118E-01  4.8487E-01  5.3009E-01  2.8569E-01
             3.0088E-02
 GRADIENT:  -2.8957E+00 -7.1806E+00  3.9114E-01  2.3866E+00  6.8234E+00 -3.7858E-01 -3.5496E-01 -6.0571E-01 -3.9165E-01  6.3986E-01
            -1.8024E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1654.95355617571        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1389
 NPARAMETR:  9.8555E-01  1.7934E+00  6.4877E-01  5.0810E-01  1.2702E+00  1.0253E+00  8.1201E-01  1.4659E+00  1.5376E+00  1.1876E+00
             9.3201E-01
 PARAMETER:  8.5442E-02  6.8410E-01 -3.3267E-01 -5.7707E-01  3.3917E-01  1.2499E-01 -1.0824E-01  4.8248E-01  5.3020E-01  2.7192E-01
             2.9593E-02
 GRADIENT:  -1.7952E+00 -5.7215E+00  2.1426E+00  2.1882E-01  1.9689E-02  6.5748E-02  2.5504E-01 -5.7877E-01  1.2717E-01 -6.6448E-02
            -2.9205E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1654.98970589194        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1566
 NPARAMETR:  9.8604E-01  1.7873E+00  6.3544E-01  5.1118E-01  1.2580E+00  1.0247E+00  8.1144E-01  1.4624E+00  1.5218E+00  1.1766E+00
             9.3133E-01
 PARAMETER:  8.5941E-02  6.8072E-01 -3.5343E-01 -5.7103E-01  3.2954E-01  1.2443E-01 -1.0894E-01  4.8007E-01  5.1986E-01  2.6263E-01
             2.8855E-02
 GRADIENT:  -8.6083E-01 -7.7119E+00  9.0247E-01  6.3467E-01 -8.9651E-01 -1.6007E-01 -2.1263E-01  6.3736E-02 -1.4541E-01  1.4055E-01
            -3.3057E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1655.01052421902        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1746             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8713E-01  1.7882E+00  6.1926E-01  5.0863E-01  1.2546E+00  1.0257E+00  8.1256E-01  1.4219E+00  1.5234E+00  1.1692E+00
             9.3140E-01
 PARAMETER:  8.7047E-02  6.8122E-01 -3.7922E-01 -5.7604E-01  3.2685E-01  1.2535E-01 -1.0756E-01  4.5203E-01  5.2094E-01  2.5634E-01
             2.8937E-02
 GRADIENT:   4.9475E+02  8.3239E+02  2.2027E+00  1.2984E+02  2.2849E+01  6.2277E+01  1.0539E+01  7.5052E-01  2.5150E+01  2.6847E+00
             5.5825E-01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1655.01052421902        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:     1806
 NPARAMETR:  9.8713E-01  1.7882E+00  6.1926E-01  5.0863E-01  1.2546E+00  1.0257E+00  8.1256E-01  1.4219E+00  1.5234E+00  1.1692E+00
             9.3140E-01
 PARAMETER:  8.7047E-02  6.8122E-01 -3.7922E-01 -5.7604E-01  3.2685E-01  1.2535E-01 -1.0756E-01  4.5203E-01  5.2094E-01  2.5634E-01
             2.8937E-02
 GRADIENT:  -1.5504E+05 -8.4591E+03 -2.5673E+04  9.6701E-01 -2.9797E+04 -9.9618E-03 -7.3175E-03  6.3276E+03 -1.1079E+04  2.2468E+04
            -1.3777E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1806
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.9646E-04 -3.3294E-02 -3.1438E-02  3.0590E-02 -4.9927E-02
 SE:             2.9879E-02  2.3469E-02  1.0498E-02  2.1831E-02  2.2008E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8674E-01  1.5602E-01  2.7469E-03  1.6115E-01  2.3293E-02

 ETASHRINKSD(%)  1.0000E-10  2.1375E+01  6.4831E+01  2.6864E+01  2.6271E+01
 ETASHRINKVR(%)  1.0000E-10  3.8181E+01  8.7632E+01  4.6511E+01  4.5641E+01
 EBVSHRINKSD(%)  3.5629E-01  2.0369E+01  6.8287E+01  2.9798E+01  2.2441E+01
 EBVSHRINKVR(%)  7.1131E-01  3.6589E+01  8.9943E+01  5.0716E+01  3.9847E+01
 RELATIVEINF(%)  9.9224E+01  3.3674E+00  9.0393E-01  2.4524E+00  1.9848E+01
 EPSSHRINKSD(%)  4.6121E+01
 EPSSHRINKVR(%)  7.0971E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1655.0105242190227     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -919.85969765528455     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1655.011       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  1.79E+00  6.19E-01  5.09E-01  1.25E+00  1.03E+00  8.13E-01  1.42E+00  1.52E+00  1.17E+00  9.31E-01
 


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
 #CPUT: Total CPU Time in Seconds,       55.408
Stop Time:
Sun Oct 24 04:27:14 CDT 2021
