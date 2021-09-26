Sat Sep 25 08:16:46 CDT 2021
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
$DATA ../../../../data/spa/A1/dat82.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m82.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1274.46768226971        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1914E+01 -3.9650E+01 -5.7817E+01  3.2056E+01  1.6850E+02  3.6332E+01 -3.7281E+01  7.0939E-01 -4.2666E+01 -6.4217E+01
            -6.4637E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1465.19667286887        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0147E+00  1.0266E+00  1.1012E+00  9.8754E-01  9.5879E-01  8.5163E-01  1.1510E+00  9.4944E-01  1.1171E+00  1.0803E+00
             1.9129E+00
 PARAMETER:  1.1462E-01  1.2625E-01  1.9640E-01  8.7459E-02  5.7921E-02 -6.0607E-02  2.4061E-01  4.8117E-02  2.1070E-01  1.7727E-01
             7.4861E-01
 GRADIENT:   5.1682E+01 -1.3236E+01 -1.2020E+01 -5.5679E+00  1.7791E+01 -1.8460E+01 -3.6796E+00  4.2260E+00  4.9045E+00 -5.3175E+00
            -1.4035E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1467.63854958871        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0113E+00  8.5723E-01  1.2720E+00  1.1388E+00  9.4342E-01  8.8417E-01  1.3874E+00  5.6657E-01  9.9628E-01  1.1473E+00
             1.9214E+00
 PARAMETER:  1.1120E-01 -5.4055E-02  3.4058E-01  2.2996E-01  4.1752E-02 -2.3107E-02  4.2741E-01 -4.6815E-01  9.6278E-02  2.3743E-01
             7.5307E-01
 GRADIENT:   4.0928E+01  1.9286E+01 -1.5329E-01  3.8133E+01  1.6423E+00 -2.9877E+00  4.9120E-01  4.4326E-01 -1.0319E+00 -1.1807E+00
            -1.4273E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1469.07821508753        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.9790E-01  7.2833E-01  1.1264E+00  1.1788E+00  8.4506E-01  8.8874E-01  1.5296E+00  1.9332E-01  9.5185E-01  1.0762E+00
             1.9717E+00
 PARAMETER:  9.7899E-02 -2.1701E-01  2.1900E-01  2.6452E-01 -6.8352E-02 -1.7951E-02  5.2502E-01 -1.5434E+00  5.0657E-02  1.7348E-01
             7.7889E-01
 GRADIENT:   3.4631E+00  2.4277E+00  1.4409E+00  1.5230E+00 -3.4430E+00 -1.7416E-01  1.6599E-01  1.0065E-01 -2.6117E-01  2.2015E-01
             1.9654E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1469.63242178328        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      341
 NPARAMETR:  9.9682E-01  5.1004E-01  1.2019E+00  1.3277E+00  7.9922E-01  8.8922E-01  1.8593E+00  3.7882E-02  8.9511E-01  1.0824E+00
             1.9633E+00
 PARAMETER:  9.6818E-02 -5.7327E-01  2.8389E-01  3.8344E-01 -1.2412E-01 -1.7410E-02  7.2019E-01 -3.1733E+00 -1.0806E-02  1.7917E-01
             7.7464E-01
 GRADIENT:  -4.1762E+00  6.9183E+00  5.4703E+00  1.4287E+01 -1.2486E+01 -1.3369E-02 -1.2172E-01  2.0030E-03 -1.7973E+00 -5.3709E-01
            -3.6302E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1470.29487026564        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  9.9447E-01  2.9649E-01  1.3697E+00  1.4571E+00  8.0766E-01  8.8495E-01  2.3355E+00  1.0000E-02  8.5687E-01  1.1197E+00
             1.9856E+00
 PARAMETER:  9.4452E-02 -1.1158E+00  4.1458E-01  4.7643E-01 -1.1361E-01 -2.2227E-02  9.4824E-01 -7.0730E+00 -5.4467E-02  2.1309E-01
             7.8594E-01
 GRADIENT:  -1.9666E+00  9.1037E-01  5.4083E-01  4.7280E-02  6.8746E-01 -1.8173E-01  2.2552E-01  0.0000E+00 -1.7104E+00 -1.6189E+00
            -2.3730E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1470.65503375186        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  9.9249E-01  1.3471E-01  1.3720E+00  1.5506E+00  7.6890E-01  8.8299E-01  3.3462E+00  1.0000E-02  8.3953E-01  1.1284E+00
             1.9866E+00
 PARAMETER:  9.2466E-02 -1.9047E+00  4.1624E-01  5.3867E-01 -1.6280E-01 -2.4443E-02  1.3078E+00 -1.3906E+01 -7.4909E-02  2.2083E-01
             7.8643E-01
 GRADIENT:  -8.4979E-01 -1.2749E-02  3.4516E-01 -3.3208E+00 -6.2395E-01 -2.7328E-01 -1.2930E-01  0.0000E+00  2.2085E+00  1.3011E-01
             7.9677E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1470.79903089459        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  9.9119E-01  6.1484E-02  1.4214E+00  1.6005E+00  7.6930E-01  8.8255E-01  4.9251E+00  1.0000E-02  8.2117E-01  1.1388E+00
             1.9869E+00
 PARAMETER:  9.1152E-02 -2.6890E+00  4.5167E-01  5.7029E-01 -1.6227E-01 -2.4936E-02  1.6943E+00 -2.1097E+01 -9.7027E-02  2.2994E-01
             7.8656E-01
 GRADIENT:  -1.5312E+00  8.6845E-02 -1.0394E+00  1.7026E+00  1.9556E+00 -8.1640E-02 -3.9844E-02  0.0000E+00  6.9152E-01 -1.4371E-01
            -1.1641E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1470.87754965240        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1042
 NPARAMETR:  9.9080E-01  1.4850E-02  1.4224E+00  1.6274E+00  7.5693E-01  8.8186E-01  9.3127E+00  1.0000E-02  8.0989E-01  1.1381E+00
             1.9882E+00
 PARAMETER:  9.0762E-02 -4.1097E+00  4.5236E-01  5.8698E-01 -1.7849E-01 -2.5724E-02  2.3314E+00 -3.4808E+01 -1.1086E-01  2.2938E-01
             7.8723E-01
 GRADIENT:  -5.6293E-01  5.7051E-03  3.1888E-01  1.5472E+00 -1.0892E+00 -1.1486E-01 -7.4532E-02  0.0000E+00 -4.7131E-01  1.3889E-01
            -1.3741E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1471.13788617038        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1221
 NPARAMETR:  9.9065E-01  1.0000E-02  1.4160E+00  1.6265E+00  7.5668E-01  8.8197E-01  1.2579E+01  1.0000E-02  8.0104E-01  1.1349E+00
             1.9888E+00
 PARAMETER:  9.0607E-02 -4.8080E+00  4.4781E-01  5.8643E-01 -1.7882E-01 -2.5601E-02  2.6320E+00 -4.1709E+01 -1.2185E-01  2.2651E-01
             7.8755E-01
 GRADIENT:  -5.8772E-01  0.0000E+00  1.3363E+00  8.1967E+00 -2.6567E+00 -3.4080E-02 -1.2460E+00  0.0000E+00  2.5581E+00  1.5669E+00
             5.7990E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1471.16950537618        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1384
 NPARAMETR:  9.9078E-01  1.0000E-02  1.3968E+00  1.6212E+00  7.5198E-01  8.8192E-01  1.2545E+01  1.0000E-02  7.9454E-01  1.1268E+00
             1.9884E+00
 PARAMETER:  9.0734E-02 -4.8212E+00  4.3417E-01  5.8318E-01 -1.8505E-01 -2.5659E-02  2.6293E+00 -4.1896E+01 -1.2999E-01  2.1936E-01
             7.8732E-01
 GRADIENT:   1.6808E-02  0.0000E+00 -1.0770E-01 -4.1262E-01  1.3417E-01  4.2247E-03  2.5671E-01  0.0000E+00 -3.5415E-02 -2.0875E-01
            -6.6300E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1384
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.7584E-04  5.0721E-03 -1.7391E-06 -1.2637E-02 -2.7341E-02
 SE:             2.9324E-02  4.5694E-03  1.1104E-04  2.7680E-02  2.2601E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7889E-01  2.6699E-01  9.8750E-01  6.4799E-01  2.2638E-01

 ETASHRINKSD(%)  1.7607E+00  8.4692E+01  9.9628E+01  7.2694E+00  2.4285E+01
 ETASHRINKVR(%)  3.4904E+00  9.7657E+01  9.9999E+01  1.4010E+01  4.2672E+01
 EBVSHRINKSD(%)  1.8667E+00  8.8146E+01  9.9615E+01  6.8715E+00  2.2383E+01
 EBVSHRINKVR(%)  3.6985E+00  9.8595E+01  9.9999E+01  1.3271E+01  3.9756E+01
 RELATIVEINF(%)  9.6064E+01  1.0009E+00  1.1133E-04  4.3698E+01  4.5437E+00
 EPSSHRINKSD(%)  3.5923E+01
 EPSSHRINKVR(%)  5.8941E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1471.1695053761803     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -736.01867881244209     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1471.170       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  1.00E-02  1.40E+00  1.62E+00  7.52E-01  8.82E-01  1.25E+01  1.00E-02  7.95E-01  1.13E+00  1.99E+00
 


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
+        1.43E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        6.46E+00  0.00E+00  1.71E+02
 
 TH 4
+       -4.48E+01  0.00E+00  3.91E+05  2.51E+05
 
 TH 5
+        1.10E+01  0.00E+00 -4.20E+02  3.71E+02  1.24E+03
 
 TH 6
+        2.64E+01  0.00E+00 -8.93E+00  8.93E+00 -9.39E+00  2.24E+02
 
 TH 7
+        1.74E-01  0.00E+00  4.34E+00  1.06E+01 -9.97E+00 -3.32E-01  2.06E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.68E+01  0.00E+00  2.78E+00  1.49E+02 -1.30E+01  3.35E+01 -4.95E+00  0.00E+00  2.34E+02
 
 TH10
+       -5.98E+00  0.00E+00  1.96E+02  9.59E+05 -4.20E+02 -8.18E+00  9.51E+00  0.00E+00  1.02E+01  4.93E+02
 
 TH11
+       -1.27E+01  0.00E+00  2.58E+01  1.51E+05 -5.69E+01  4.76E+00  1.71E+00  0.00E+00  1.31E+01  8.54E+01  7.83E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.921
Stop Time:
Sat Sep 25 08:17:11 CDT 2021
