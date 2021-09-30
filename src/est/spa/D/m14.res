Wed Sep 29 19:45:24 CDT 2021
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
$DATA ../../../../data/spa/D/dat14.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1124.33599626464        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0598E+02 -1.8255E+02 -1.1841E+02 -1.9677E+02  2.2045E+02 -3.1235E+02 -2.3107E+02  4.1603E-01 -4.0660E+02 -6.6119E+01
            -6.8754E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1451.21290308576        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.2120E+00  1.1482E+00  2.3142E+00  1.1080E+00  1.4347E+00  1.2534E+00  2.4668E+00  1.0142E+00  1.5316E+00  1.4367E+00
             1.1352E+00
 PARAMETER:  2.9227E-01  2.3822E-01  9.3908E-01  2.0256E-01  4.6096E-01  3.2583E-01  1.0029E+00  1.1413E-01  5.2629E-01  4.6234E-01
             2.2683E-01
 GRADIENT:   1.1376E+03  7.2553E+01  1.2153E+00  1.4317E+02  7.3574E+01  3.6761E+01  1.9670E+02 -2.8013E-01  5.3026E+01  2.2106E+01
             5.3463E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1466.10312841037        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      264
 NPARAMETR:  1.2654E+00  1.2746E+00  1.1728E+00  7.2557E-01  1.2019E+00  1.7136E+00  3.1184E+00  1.0414E+00  1.4099E+00  9.4105E-01
             1.1050E+00
 PARAMETER:  3.3540E-01  3.4260E-01  2.5936E-01 -2.2080E-01  2.8391E-01  6.3863E-01  1.2373E+00  1.4057E-01  4.4353E-01  3.9243E-02
             1.9980E-01
 GRADIENT:   1.1677E+02 -4.3948E+01  4.3189E+00 -1.1186E+02  1.2039E+01 -1.9307E+01  6.1335E+01  2.9482E+00 -2.0644E+01 -8.3358E+00
            -7.9102E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1494.90738646852        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.2003E+00  8.0242E-01  1.0445E+00  1.2067E+00  8.9829E-01  1.5100E+00  3.6697E+00  2.8714E-01  1.2312E+00  7.5961E-01
             1.0805E+00
 PARAMETER:  2.8258E-01 -1.2012E-01  1.4358E-01  2.8791E-01 -7.2654E-03  5.1209E-01  1.4001E+00 -1.1478E+00  3.0799E-01 -1.7496E-01
             1.7746E-01
 GRADIENT:   9.1620E+01 -4.6072E+00 -3.7615E+01 -5.4339E+00  4.0855E+01 -7.4858E+01  1.6416E+01  6.4571E-01  1.2514E+01 -7.6705E-01
            -1.5359E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1500.37277655707        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      631             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1681E+00  8.2978E-01  1.1022E+00  1.2065E+00  9.2116E-01  1.5579E+00  3.4935E+00  2.7350E-01  1.2301E+00  8.0542E-01
             1.0843E+00
 PARAMETER:  2.5539E-01 -8.6591E-02  1.9726E-01  2.8770E-01  1.7883E-02  5.4333E-01  1.3509E+00 -1.1964E+00  3.0706E-01 -1.1639E-01
             1.8094E-01
 GRADIENT:   9.7563E+02  2.2090E+01 -1.9735E+01  2.5403E+02  2.8115E+01  3.8740E+02  4.1671E+02  5.1696E-01  3.3887E+01  1.4973E+00
             1.2631E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1504.56866487152        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      814
 NPARAMETR:  1.1681E+00  8.2978E-01  1.1507E+00  1.2056E+00  9.2116E-01  1.7746E+00  3.4068E+00  2.7597E-01  1.2301E+00  8.0542E-01
             1.0853E+00
 PARAMETER:  2.5538E-01 -8.6591E-02  2.4038E-01  2.8698E-01  1.7883E-02  6.7357E-01  1.3258E+00 -1.1875E+00  3.0706E-01 -1.1639E-01
             1.8184E-01
 GRADIENT:   4.8572E+01 -1.2925E+00  3.3828E-01 -1.7274E+01 -1.3022E+01  9.0975E+00 -1.7989E+00  1.6918E-01  5.7126E+00  1.5922E-01
            -1.2466E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1506.09693839580        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  1.0985E+00  8.2257E-01  1.1627E+00  1.2281E+00  9.2680E-01  1.7341E+00  3.4674E+00  1.9587E-01  1.1781E+00  8.2271E-01
             1.0904E+00
 PARAMETER:  1.9396E-01 -9.5326E-02  2.5076E-01  3.0548E-01  2.3987E-02  6.5051E-01  1.3434E+00 -1.5303E+00  2.6394E-01 -9.5152E-02
             1.8652E-01
 GRADIENT:   1.4390E+00  2.0365E+00 -3.5785E+00 -5.1940E+00 -2.1619E+00  3.9961E+00 -2.1125E+00  4.2869E-02  8.7529E-01 -1.1858E+00
            -1.1370E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1506.26901434162        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1185             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1001E+00  8.0745E-01  1.1735E+00  1.2349E+00  9.2805E-01  1.7472E+00  3.4975E+00  1.8653E-01  1.1729E+00  8.3233E-01
             1.0907E+00
 PARAMETER:  1.9541E-01 -1.1388E-01  2.5995E-01  3.1098E-01  2.5326E-02  6.5803E-01  1.3520E+00 -1.5792E+00  2.5946E-01 -8.3523E-02
             1.8680E-01
 GRADIENT:   6.9649E+02  2.7596E+01  1.2522E-01  2.8250E+02  4.4926E+00  5.7484E+02  3.9982E+02  9.1587E-02  2.1784E+01 -3.8092E-01
             1.5803E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1506.30294897012        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1367
 NPARAMETR:  1.0999E+00  8.0676E-01  1.1814E+00  1.2363E+00  9.2815E-01  1.7616E+00  3.4910E+00  1.8461E-01  1.1720E+00  8.3277E-01
             1.0912E+00
 PARAMETER:  1.9522E-01 -1.1473E-01  2.6673E-01  3.1214E-01  2.5443E-02  6.6625E-01  1.3502E+00 -1.5895E+00  2.5871E-01 -8.2995E-02
             1.8729E-01
 GRADIENT:   2.6401E+00  1.8106E+00 -6.9663E-01 -6.9159E+00 -6.1131E+00  1.0862E+01 -3.4861E+00  1.7546E-02  3.7469E-01 -9.6490E-01
            -5.5317E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1506.31074611268        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1467
 NPARAMETR:  1.0998E+00  8.0679E-01  1.1816E+00  1.2364E+00  9.2813E-01  1.7618E+00  3.4971E+00  1.8471E-01  1.1719E+00  8.3376E-01
             1.0913E+00
 PARAMETER:  1.9516E-01 -1.1469E-01  2.6683E-01  3.1224E-01  2.5412E-02  6.6635E-01  1.3519E+00 -1.5890E+00  2.5863E-01 -8.1815E-02
             1.8735E-01
 GRADIENT:   6.1974E+03 -1.0545E+04 -2.2739E+03 -1.9444E+03  1.2083E+04 -5.7965E-01  8.9096E+02 -7.6181E+02  2.3347E+03  1.2089E+04
            -6.4648E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1467
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3065E-03  2.4615E-02 -8.0488E-03 -2.7854E-02  3.0845E-03
 SE:             2.9976E-02  2.3992E-02  2.5716E-03  2.3450E-02  2.0050E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3867E-01  3.0492E-01  1.7489E-03  2.3490E-01  8.7774E-01

 ETASHRINKSD(%)  1.0000E-10  1.9623E+01  9.1385E+01  2.1441E+01  3.2829E+01
 ETASHRINKVR(%)  1.0000E-10  3.5395E+01  9.9258E+01  3.8285E+01  5.4880E+01
 EBVSHRINKSD(%)  1.6928E-01  1.7317E+01  9.2001E+01  2.4203E+01  3.0872E+01
 EBVSHRINKVR(%)  3.3828E-01  3.1635E+01  9.9360E+01  4.2548E+01  5.2214E+01
 RELATIVEINF(%)  9.9574E+01  2.0968E+01  1.1112E-01  1.5489E+01  8.2001E+00
 EPSSHRINKSD(%)  4.2487E+01
 EPSSHRINKVR(%)  6.6922E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1506.3107461126804     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -771.15991954894218     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1506.311       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.10E+00  8.07E-01  1.18E+00  1.24E+00  9.28E-01  1.76E+00  3.50E+00  1.85E-01  1.17E+00  8.34E-01  1.09E+00
 


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
+        6.57E+05
 
 TH 2
+        7.15E+02  3.53E+06
 
 TH 3
+       -2.48E+03 -1.03E+06  3.06E+05
 
 TH 4
+       -1.11E+03 -3.31E+02  1.27E+03  2.04E+05
 
 TH 5
+       -1.52E+06  1.56E+03 -6.06E+03 -2.35E+03  3.51E+06
 
 TH 6
+       -1.15E-01  2.78E+00  1.62E+00  3.58E-01 -3.93E+00  6.42E+01
 
 TH 7
+       -3.00E+01  6.92E+04 -1.18E+02  4.61E+00 -5.70E+01 -1.55E-01  1.36E+03
 
 TH 8
+        4.80E+05 -5.25E+02  1.80E+03  8.12E+02  1.83E+03  9.03E-01  2.24E+01  3.51E+05
 
 TH 9
+       -1.91E+03  4.95E+02 -1.78E+03 -7.70E+02 -1.08E+06 -1.08E+00 -2.12E+04  3.41E+05  3.28E+05
 
 TH10
+       -2.57E+03  3.92E+06 -6.38E+03  1.42E+03 -6.01E+03 -3.41E+00 -7.90E+01  1.89E+03 -1.20E+06  4.34E+06
 
 TH11
+        6.89E+05 -7.56E+02  2.58E+03  1.16E+03 -3.27E+03  2.27E+00  3.30E+01  1.04E+03  4.89E+05  2.74E+03  7.25E+05
 
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
 #CPUT: Total CPU Time in Seconds,       30.599
Stop Time:
Wed Sep 29 19:45:56 CDT 2021
