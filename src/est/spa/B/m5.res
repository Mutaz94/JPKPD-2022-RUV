Sat Sep 25 07:05:01 CDT 2021
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
$DATA ../../../../data/spa/B/dat5.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1656.72992715986        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5594E+02 -3.7838E+01  1.9414E+00 -5.4796E+01 -3.8604E+01  1.4918E+01 -1.2142E+01  7.8784E+00  9.1917E+00  8.4072E+00
             3.7878E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1666.24951379073        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.6358E-01  1.1574E+00  1.0906E+00  9.3922E-01  1.1719E+00  9.2924E-01  1.2872E+00  9.2634E-01  8.8803E-01  1.0260E+00
             8.9372E-01
 PARAMETER:  6.2904E-02  2.4617E-01  1.8670E-01  3.7297E-02  2.5864E-01  2.6610E-02  3.5251E-01  2.3481E-02 -1.8745E-02  1.2571E-01
            -1.2359E-02
 GRADIENT:   9.1471E+01  1.7343E+01  1.1250E+01 -5.0076E+00  1.7294E+01 -7.4920E+00  1.2907E+01 -4.2457E+00 -1.7968E+00 -1.2212E+01
            -1.0537E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1666.88452804302        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.6219E-01  1.0414E+00  1.0981E+00  1.0219E+00  1.1070E+00  9.1681E-01  1.3971E+00  8.4349E-01  8.3623E-01  1.0465E+00
             8.9332E-01
 PARAMETER:  6.1455E-02  1.4053E-01  1.9356E-01  1.2162E-01  2.0162E-01  1.3143E-02  4.3439E-01 -7.0204E-02 -7.8852E-02  1.4543E-01
            -1.2811E-02
 GRADIENT:   9.0038E+01  1.9628E+01  1.0373E+01  1.2414E+01  4.0974E+00 -1.2884E+01  1.1796E+01 -3.5900E+00 -1.8133E+00 -2.5015E+00
            -9.1173E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1668.46839721241        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.3798E-01  1.0424E+00  8.4803E-01  9.9725E-01  9.7249E-01  9.4304E-01  1.3127E+00  5.8565E-01  8.4901E-01  9.1010E-01
             9.0366E-01
 PARAMETER:  3.5975E-02  1.4152E-01 -6.4839E-02  9.7245E-02  7.2107E-02  4.1353E-02  3.7205E-01 -4.3503E-01 -6.3683E-02  5.7967E-03
            -1.3016E-03
 GRADIENT:   1.7692E+01  1.4370E+00 -6.9614E+00  8.2996E+00  9.7967E+00 -1.4158E+00  3.0323E+00  1.3874E+00  4.1701E-01  2.2260E-01
            -1.2303E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1668.48138720566        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.3524E-01  1.0409E+00  7.9663E-01  9.9218E-01  9.3696E-01  9.4546E-01  1.3056E+00  4.8791E-01  8.4668E-01  8.7467E-01
             9.0504E-01
 PARAMETER:  3.3050E-02  1.4005E-01 -1.2736E-01  9.2150E-02  3.4880E-02  4.3914E-02  3.6668E-01 -6.1762E-01 -6.6430E-02 -3.3912E-02
             2.2534E-04
 GRADIENT:   9.1463E+00  3.9692E-01 -5.0209E+00  5.0870E+00  5.9836E+00 -7.4147E-01  1.7763E+00  1.1425E+00  4.0460E-01  4.1666E-01
            -4.0909E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1668.48292566957        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.3403E-01  1.0394E+00  7.6476E-01  9.8924E-01  9.1446E-01  9.4671E-01  1.3030E+00  4.1934E-01  8.4432E-01  8.5152E-01
             9.0556E-01
 PARAMETER:  3.1751E-02  1.3865E-01 -1.6819E-01  8.9181E-02  1.0582E-02  4.5239E-02  3.6466E-01 -7.6906E-01 -6.9222E-02 -6.0737E-02
             8.0192E-04
 GRADIENT:   5.1501E+00 -1.0627E-02 -3.8137E+00  3.3874E+00  3.9354E+00 -4.2339E-01  1.1337E+00  9.4183E-01  3.4790E-01  4.3925E-01
            -7.9396E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1668.48387255757        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  9.3324E-01  1.0377E+00  7.4243E-01  9.8739E-01  8.9818E-01  9.4756E-01  1.3017E+00  3.6304E-01  8.4227E-01  8.3458E-01
             9.0586E-01
 PARAMETER:  3.0907E-02  1.3705E-01 -1.9782E-01  8.7313E-02 -7.3820E-03  4.6134E-02  3.6365E-01 -9.1324E-01 -7.1649E-02 -8.0825E-02
             1.1242E-03
 GRADIENT:   2.5260E+00 -2.1394E-01 -2.6898E+00  2.0876E+00  2.3682E+00 -2.1541E-01  6.5625E-01  7.4607E-01  2.6721E-01  3.7556E-01
             9.0020E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1668.48507917220        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      513
 NPARAMETR:  9.3256E-01  1.0357E+00  7.2078E-01  9.8572E-01  8.8216E-01  9.4833E-01  1.3009E+00  2.9986E-01  8.4007E-01  8.1784E-01
             9.0607E-01
 PARAMETER:  3.0183E-02  1.3512E-01 -2.2742E-01  8.5618E-02 -2.5383E-02  4.6949E-02  3.6305E-01 -1.1045E+00 -7.4270E-02 -1.0109E-01
             1.3622E-03
 GRADIENT:   2.2438E-01 -3.8248E-01 -1.6022E+00  8.9156E-01  9.4257E-01 -3.2202E-02  2.2031E-01  5.3745E-01  1.8386E-01  2.9671E-01
             2.1993E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1668.48840538931        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  9.3217E-01  1.0342E+00  7.0423E-01  9.8442E-01  8.6993E-01  9.4885E-01  1.3008E+00  2.4393E-01  8.3831E-01  8.0498E-01
             9.0617E-01
 PARAMETER:  2.9764E-02  1.3362E-01 -2.5065E-01  8.4299E-02 -3.9341E-02  4.7494E-02  3.6297E-01 -1.3109E+00 -7.6364E-02 -1.1693E-01
             1.4731E-03
 GRADIENT:  -1.2013E+00 -5.1013E-01 -1.0134E+00  1.6823E-01  1.2578E-01  8.4615E-02 -2.4291E-02  3.7850E-01  1.4092E-01  2.7919E-01
             3.0726E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1668.49593149215        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  9.3180E-01  1.0324E+00  6.8714E-01  9.8309E-01  8.5709E-01  9.4935E-01  1.3007E+00  1.6528E-01  8.3643E-01  7.9130E-01
             9.0622E-01
 PARAMETER:  2.9361E-02  1.3184E-01 -2.7521E-01  8.2948E-02 -5.4207E-02  4.8025E-02  3.6288E-01 -1.7001E+00 -7.8608E-02 -1.3408E-01
             1.5287E-03
 GRADIENT:  -2.5892E+00 -5.8744E-01 -1.2813E-01 -6.8582E-01 -8.4951E-01  1.9516E-01 -3.2059E-01  1.8172E-01  6.6235E-02  1.7714E-01
             3.4659E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1668.56745027395        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  9.3282E-01  1.0359E+00  6.6885E-01  9.7918E-01  8.4701E-01  9.4903E-01  1.3013E+00  2.1251E-02  8.3569E-01  7.8059E-01
             9.0574E-01
 PARAMETER:  3.0461E-02  1.3524E-01 -3.0220E-01  7.8957E-02 -6.6043E-02  4.7685E-02  3.6338E-01 -3.7513E+00 -7.9499E-02 -1.4771E-01
             9.9620E-04
 GRADIENT:  -3.9261E-01 -4.3831E-01 -2.1359E+00  6.0678E-01  1.0256E+00  6.2015E-02  4.1197E-01  3.9148E-03  2.0093E-01  6.4063E-01
             3.5321E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1668.99380313242        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  9.4261E-01  1.0520E+00  6.8678E-01  9.7545E-01  8.6832E-01  9.5301E-01  1.2973E+00  1.0000E-02  8.4651E-01  8.0289E-01
             9.0649E-01
 PARAMETER:  4.0894E-02  1.5068E-01 -2.7574E-01  7.5147E-02 -4.1194E-02  5.1867E-02  3.6031E-01 -1.3083E+01 -6.6631E-02 -1.1954E-01
             1.8269E-03
 GRADIENT:  -1.9282E+01 -2.4551E+00 -1.9998E+00 -6.1494E-01  2.3865E+00 -1.3847E+00 -4.0767E-01  0.0000E+00  5.9074E-01 -1.3727E-01
             1.6332E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1669.07932318536        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  9.5008E-01  1.0524E+00  6.9125E-01  9.7714E-01  8.7026E-01  9.5596E-01  1.3024E+00  1.0000E-02  8.4303E-01  8.0837E-01
             9.0663E-01
 PARAMETER:  4.8791E-02  1.5110E-01 -2.6925E-01  7.6877E-02 -3.8962E-02  5.4962E-02  3.6419E-01 -1.7949E+01 -7.0754E-02 -1.1273E-01
             1.9823E-03
 GRADIENT:  -2.3609E-02  2.1485E-02  1.7982E-02 -2.3175E-02  4.3124E-03  6.4952E-03  1.3300E-02  0.0000E+00 -3.0141E-02 -1.6948E-02
            -2.0248E-03

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1669.07932870296        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  9.5009E-01  1.0522E+00  6.9127E-01  9.7729E-01  8.7014E-01  9.5594E-01  1.3025E+00  1.0000E-02  8.4309E-01  8.0838E-01
             9.0663E-01
 PARAMETER:  4.8800E-02  1.5087E-01 -2.6923E-01  7.7027E-02 -3.9096E-02  5.4944E-02  3.6425E-01 -1.7936E+01 -7.0676E-02 -1.1273E-01
             1.9754E-03
 GRADIENT:   1.0853E-03 -4.5205E-04  2.3711E-04  5.9031E-04  2.3640E-03 -4.2744E-04 -2.0404E-04  0.0000E+00  3.2747E-04 -2.0537E-03
            -3.3302E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1084
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6192E-05 -2.0248E-04 -5.1023E-04 -3.0638E-03 -1.1685E-02
 SE:             2.9846E-02  2.3763E-02  2.0115E-04  2.3567E-02  2.1966E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9930E-01  9.9320E-01  1.1196E-02  8.9656E-01  5.9477E-01

 ETASHRINKSD(%)  1.1846E-02  2.0392E+01  9.9326E+01  2.1048E+01  2.6411E+01
 ETASHRINKVR(%)  2.3691E-02  3.6625E+01  9.9995E+01  3.7666E+01  4.5846E+01
 EBVSHRINKSD(%)  3.8903E-01  1.9634E+01  9.9412E+01  2.1750E+01  2.5707E+01
 EBVSHRINKVR(%)  7.7655E-01  3.5414E+01  9.9997E+01  3.8769E+01  4.4806E+01
 RELATIVEINF(%)  9.9018E+01  5.6542E+00  4.3024E-04  5.4346E+00  5.3280E+00
 EPSSHRINKSD(%)  4.4284E+01
 EPSSHRINKVR(%)  6.8957E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1669.0793287029649     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -933.92850213922668     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.47
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1669.079       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.50E-01  1.05E+00  6.91E-01  9.77E-01  8.70E-01  9.56E-01  1.30E+00  1.00E-02  8.43E-01  8.08E-01  9.07E-01
 


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
+        1.34E+03
 
 TH 2
+       -5.51E+00  3.95E+02
 
 TH 3
+        1.85E+01  2.32E+02  9.19E+02
 
 TH 4
+       -1.07E+01  3.00E+02 -4.65E+02  1.02E+03
 
 TH 5
+       -4.33E+00 -3.33E+02 -8.99E+02  4.91E+02  1.23E+03
 
 TH 6
+        5.83E-01 -1.51E+00  4.12E+00 -3.29E+00 -2.28E+00  2.14E+02
 
 TH 7
+        7.94E-01  2.89E+01 -2.43E+01 -1.51E+01  6.32E+00 -8.63E-02  5.18E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.37E+00 -2.11E+01 -4.45E+01  3.52E+01  5.86E+00 -1.58E+00  1.74E+01  0.00E+00  1.22E+02
 
 TH10
+       -1.26E+00 -1.44E+01 -7.49E+01 -2.23E+01 -6.59E+01  2.95E-01  1.10E+01  0.00E+00  1.50E+01  1.08E+02
 
 TH11
+       -7.99E+00 -1.30E+01 -4.34E+01 -8.33E-01  7.86E+00  1.70E+00  6.24E+00  0.00E+00  1.06E+01  2.06E+01  2.61E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.278
Stop Time:
Sat Sep 25 07:05:18 CDT 2021
