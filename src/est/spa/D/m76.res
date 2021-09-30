Wed Sep 29 20:19:16 CDT 2021
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
$DATA ../../../../data/spa/D/dat76.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17174.0072116104        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6259E+02  2.2112E+02  3.3446E-01  1.6525E+02  2.4742E+02 -2.6226E+03 -1.2462E+03 -4.9806E+01 -1.4521E+03 -5.5492E+02
            -3.1130E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -559.481027099078        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1873E+00  1.0909E+00  8.0334E-01  1.7498E+00  1.4140E+00  2.4731E+00  1.1216E+00  9.7498E-01  1.0271E+00  1.0144E+00
             1.4241E+01
 PARAMETER:  2.7169E-01  1.8703E-01 -1.1898E-01  6.5948E-01  4.4639E-01  1.0055E+00  2.1478E-01  7.4659E-02  1.2676E-01  1.1428E-01
             2.7561E+00
 GRADIENT:  -2.9034E+01  3.3617E+01 -1.9471E+00  6.2150E+01 -6.7932E+00  4.9287E+01 -5.2443E-01  6.1468E+00  4.8074E+00  9.1778E-01
             4.9767E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -577.704486994131        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2217E+00  7.7397E-01  9.9840E-01  1.6381E+00  2.1898E+00  1.9119E+00  2.3572E+00  3.1441E-01  6.4699E-01  2.3378E+00
             1.4253E+01
 PARAMETER:  3.0025E-01 -1.5623E-01  9.8396E-02  5.9354E-01  8.8380E-01  7.4809E-01  9.5746E-01 -1.0570E+00 -3.3543E-01  9.4922E-01
             2.7570E+00
 GRADIENT:  -5.8730E+00  1.7997E+01  6.2898E+00  2.5664E+01 -5.9050E+00 -1.5960E+00  4.0387E+00  3.2837E-01  6.5268E+00  2.5758E+00
             6.3612E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -597.547691076826        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1095E+00  1.2843E-01  6.2337E-01  1.3883E+00  3.7046E+00  1.7833E+00  5.7910E+00  1.0000E-02  3.1120E-01  3.5171E+00
             1.1408E+01
 PARAMETER:  2.0395E-01 -1.9524E+00 -3.7261E-01  4.2808E-01  1.4096E+00  6.7845E-01  1.8563E+00 -6.7275E+00 -1.0673E+00  1.3576E+00
             2.5343E+00
 GRADIENT:   7.8985E+01  2.5262E+00  5.2649E+01 -1.3954E+02 -2.5900E+01 -1.4922E+00  4.4608E+00  0.0000E+00  2.9628E+00  2.9625E+00
            -6.4250E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -652.506408451093        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  8.1296E-01  1.3471E-02  7.5100E-02  7.4322E-01  3.5322E+01  1.9811E+00  4.6146E+00  1.0000E-02  4.6743E-02  4.3869E-01
             9.5955E+00
 PARAMETER: -1.0707E-01 -4.2072E+00 -2.4889E+00 -1.9677E-01  3.6645E+00  7.8364E-01  1.6292E+00 -1.9525E+01 -2.9631E+00 -7.2397E-01
             2.3613E+00
 GRADIENT:   1.2013E+02  3.0016E-01 -2.0694E+02  3.7698E+02  3.8603E-01 -4.8525E+01  1.1336E-02  0.0000E+00 -1.4866E-01 -3.1819E-05
            -2.4538E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -712.254111766547        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      415
 NPARAMETR:  6.0050E-01  1.1822E-02  5.5157E-02  5.2670E-01  8.8881E+01  1.9954E+00  5.2507E-01  1.0000E-02  1.0000E-02  3.8482E-01
             1.1257E+01
 PARAMETER: -4.0999E-01 -4.3378E+00 -2.7976E+00 -5.4113E-01  4.5873E+00  7.9084E-01 -5.4423E-01 -1.9959E+01 -4.5431E+00 -8.5499E-01
             2.5210E+00
 GRADIENT:  -1.5850E+01  2.6933E-02  7.7206E+00  3.7962E+00 -1.4343E-03 -8.1233E+00  1.8778E-05  0.0000E+00  4.2034E-04 -7.8406E-07
             3.0339E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -713.045226149847        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:      569
 NPARAMETR:  5.8986E-01  1.7191E-02  4.6859E-02  4.7082E-01  7.8784E+01  2.0235E+00  3.3794E-01  1.0000E-02  1.0000E-02  4.9281E-01
             1.1207E+01
 PARAMETER: -4.2787E-01 -3.9634E+00 -2.9606E+00 -6.5328E-01  4.4667E+00  8.0481E-01 -9.8488E-01 -2.0435E+01 -4.7742E+00 -6.0763E-01
             2.5166E+00
 GRADIENT:   2.6368E+01 -1.9678E-02  3.1399E+01  1.7956E+01  8.3104E-03  1.6415E+01  1.0412E-04  0.0000E+00  0.0000E+00 -4.9166E-07
             1.6576E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -713.067043242911        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  5.9313E-01  2.0766E-02  4.6676E-02  4.6950E-01  7.8907E+00  2.0196E+00  4.4450E-02  1.0000E-02  1.0000E-02  3.9565E-01
             1.1204E+01
 PARAMETER: -4.2234E-01 -3.7744E+00 -2.9645E+00 -6.5609E-01  2.1657E+00  8.0292E-01 -3.0134E+00 -2.0435E+01 -4.7742E+00 -8.2723E-01
             2.5163E+00
 GRADIENT:   1.4180E+00  7.6526E-01  4.9178E+00 -1.0676E+01  6.4553E-02 -1.0224E+00  2.4147E-04  0.0000E+00  0.0000E+00  2.5363E-02
            -1.5485E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -713.249855054968        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  5.9263E-01  1.7646E-02  4.6386E-02  4.7283E-01  6.7192E+00  2.0284E+00  1.0000E-02  1.0000E-02  1.0000E-02  6.6726E-02
             1.1230E+01
 PARAMETER: -4.2319E-01 -3.9372E+00 -2.9708E+00 -6.4902E-01  2.0050E+00  8.0724E-01 -5.9872E+00 -2.0435E+01 -4.7742E+00 -2.6072E+00
             2.5186E+00
 GRADIENT:   2.5061E-01  1.2936E-01 -3.4251E+00  2.1024E+00  4.4716E-01 -8.4265E-01  0.0000E+00  0.0000E+00  0.0000E+00  4.2920E-03
            -3.6823E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -713.269992791411        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1032
 NPARAMETR:  5.9226E-01  1.6131E-02  4.6492E-02  4.7347E-01  6.4972E+00  2.0337E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0383E-02
             1.1227E+01
 PARAMETER: -4.2380E-01 -4.0270E+00 -2.9685E+00 -6.4766E-01  1.9714E+00  8.0986E-01 -6.8170E+00 -2.0435E+01 -4.7742E+00 -4.4676E+00
             2.5184E+00
 GRADIENT:  -3.6393E-01  7.8450E-03 -5.6878E-01 -1.3274E+00 -6.1987E-02 -1.0783E-01  0.0000E+00  0.0000E+00  0.0000E+00  1.3202E-04
            -1.7287E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -713.270028241930        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1099
 NPARAMETR:  5.9363E-01  1.6096E-02  4.6503E-02  4.7395E-01  6.5112E+00  2.0393E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.1225E+01
 PARAMETER: -4.2387E-01 -4.0278E+00 -2.9683E+00 -6.4746E-01  1.9718E+00  8.0956E-01 -6.8170E+00 -2.0435E+01 -4.7742E+00 -4.5123E+00
             2.5183E+00
 GRADIENT:  -6.0125E-01  1.3548E-03 -4.0431E-02 -8.6065E-01 -2.3902E-02 -5.3935E-01  0.0000E+00  0.0000E+00  0.0000E+00  2.9015E-05
             9.0240E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1099
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2734E-03  4.1260E-06  1.3004E-04 -2.7858E-04  1.5317E-06
 SE:             2.9139E-02  3.7971E-06  2.5501E-04  3.4080E-04  3.0034E-05
 N:                     100         100         100         100         100

 P VAL.:         9.1055E-01  2.7721E-01  6.1010E-01  4.1368E-01  9.5933E-01

 ETASHRINKSD(%)  2.3822E+00  9.9987E+01  9.9146E+01  9.8858E+01  9.9899E+01
 ETASHRINKVR(%)  4.7076E+00  1.0000E+02  9.9993E+01  9.9987E+01  1.0000E+02
 EBVSHRINKSD(%)  2.4654E+00  9.9981E+01  9.9094E+01  9.8781E+01  9.9873E+01
 EBVSHRINKVR(%)  4.8701E+00  1.0000E+02  9.9992E+01  9.9985E+01  1.0000E+02
 RELATIVEINF(%)  2.2376E+01  6.9425E-07  1.1682E-04  2.0978E-04  3.2022E-05
 EPSSHRINKSD(%)  6.6037E+00
 EPSSHRINKVR(%)  1.2771E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -713.27002824192982     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       21.880798321808356     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -713.270       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.92E-01  1.61E-02  4.65E-02  4.74E-01  6.50E+00  2.03E+00  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.12E+01
 


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
+        7.22E+02
 
 TH 2
+       -2.72E+02  3.89E+03
 
 TH 3
+       -2.08E+03  8.27E+03  2.50E+05
 
 TH 4
+       -1.35E+02 -1.31E+03 -3.30E+04  4.75E+03
 
 TH 5
+        1.96E+00 -1.44E+01 -1.74E+02  2.46E+01  3.31E-01
 
 TH 6
+        3.18E+00  2.12E+01  2.83E+02 -5.20E+01  4.36E-02  4.29E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        9.64E-03  2.55E-01 -9.83E-02 -1.53E-02 -5.74E-03 -1.95E-03  0.00E+00  0.00E+00  0.00E+00  7.33E+00
 
 TH11
+       -1.07E+01  1.35E+01  1.51E+02 -1.81E+01 -1.53E-01  7.29E-01  0.00E+00  0.00E+00  0.00E+00  7.24E-05  3.20E+00
 
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
 #CPUT: Total CPU Time in Seconds,       21.346
Stop Time:
Wed Sep 29 20:19:39 CDT 2021
