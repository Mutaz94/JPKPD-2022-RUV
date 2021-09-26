Sat Sep 25 10:33:37 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat49.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1584.76357570073        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4385E+02 -5.0281E+01 -1.2809E+01 -4.3630E+01  1.9114E+01 -2.1071E+01 -3.7344E+01  1.7599E+00 -2.5898E+01 -1.5510E+01
            -1.1788E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1592.27031181313        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.6492E-01  1.1351E+00  1.0643E+00  9.2387E-01  1.0930E+00  1.0995E+00  1.5614E+00  9.7218E-01  1.2382E+00  1.1723E+00
             1.0670E+00
 PARAMETER:  6.4292E-02  2.2671E-01  1.6233E-01  2.0814E-02  1.8894E-01  1.9487E-01  5.4561E-01  7.1781E-02  3.1364E-01  2.5897E-01
             1.6481E-01
 GRADIENT:   5.6785E+01 -2.5587E+00  3.3784E+00 -1.9622E+01 -1.3278E+01  2.7040E+01  2.9280E+01  2.7642E+00  3.0918E+01  1.0420E+01
             2.0826E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1594.78504799560        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  9.6857E-01  1.0358E+00  1.2004E+00  1.0317E+00  1.1172E+00  1.1088E+00  1.5654E+00  8.5190E-01  1.0917E+00  1.2500E+00
             1.0448E+00
 PARAMETER:  6.8067E-02  1.3514E-01  2.8263E-01  1.3125E-01  2.1080E-01  2.0330E-01  5.4817E-01 -6.0287E-02  1.8773E-01  3.2316E-01
             1.4379E-01
 GRADIENT:   6.6661E+01  1.4152E+01  1.3260E+01  1.6737E+01  1.9565E+00  3.1108E+01  1.6260E+01 -5.3279E+00  1.8283E+01  7.6693E+00
             7.8025E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1598.72072117427        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.3559E-01  1.1029E+00  7.8063E-01  9.4628E-01  9.3625E-01  1.0390E+00  1.4145E+00  4.5453E-01  9.9377E-01  9.8046E-01
             1.0135E+00
 PARAMETER:  3.3418E-02  1.9798E-01 -1.4765E-01  4.4786E-02  3.4124E-02  1.3826E-01  4.4675E-01 -6.8850E-01  9.3749E-02  8.0265E-02
             1.1336E-01
 GRADIENT:  -5.0672E+00 -1.1095E+00 -5.5816E+00  4.4085E+00  3.0754E+00  1.1908E+00  1.3193E+00  1.2292E+00  2.4093E+00  2.4461E+00
             1.1218E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1599.14508921913        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.3877E-01  1.1150E+00  7.0653E-01  9.2937E-01  8.9354E-01  1.0363E+00  1.4077E+00  1.5618E-01  9.7628E-01  9.1994E-01
             1.0098E+00
 PARAMETER:  3.6817E-02  2.0889E-01 -2.4739E-01  2.6748E-02 -1.2565E-02  1.3565E-01  4.4195E-01 -1.7567E+00  7.5994E-02  1.6553E-02
             1.0973E-01
 GRADIENT:   4.3946E-01  1.0380E+00  4.9651E-02 -9.3574E-01 -2.5727E+00 -2.1288E-01  1.6087E-01  1.7692E-01 -4.7824E-01  5.6226E-01
            -8.0848E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1599.16362051844        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  9.3868E-01  1.1143E+00  7.0647E-01  9.2968E-01  8.9379E-01  1.0363E+00  1.4077E+00  1.3482E-01  9.7682E-01  9.2017E-01
             1.0099E+00
 PARAMETER:  3.6721E-02  2.0825E-01 -2.4747E-01  2.7085E-02 -1.2279E-02  1.3569E-01  4.4199E-01 -1.9038E+00  7.6544E-02  1.6806E-02
             1.0984E-01
 GRADIENT:   2.2430E+01 -2.5156E+01  3.8911E+00 -3.7255E+01  1.8149E+00  5.1758E+01  5.9119E+00  1.2449E-01  1.2630E+01  3.3146E-01
            -1.4178E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1599.54978459538        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      596            RESET HESSIAN, TYPE II
 NPARAMETR:  9.5734E-01  1.1183E+00  7.0915E-01  9.3416E-01  8.9330E-01  1.0356E+00  1.4001E+00  1.3568E-01  9.7577E-01  9.1891E-01
             1.0102E+00
 PARAMETER:  5.6403E-02  2.1178E-01 -2.4369E-01  3.1889E-02 -1.2833E-02  1.3501E-01  4.3652E-01 -1.8974E+00  7.5475E-02  1.5428E-02
             1.1017E-01
 GRADIENT:   4.0551E+01  6.5261E+00  2.1136E+00  5.6293E+00 -4.2269E+00 -2.0524E-01 -1.0100E+00  1.1229E-01 -7.7393E-01 -4.1744E-01
            -5.8028E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1599.64342058863        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  9.5705E-01  1.1181E+00  7.0913E-01  9.3404E-01  8.9882E-01  1.0480E+00  1.4012E+00  1.3290E-01  9.7624E-01  9.2556E-01
             1.0140E+00
 PARAMETER:  5.6103E-02  2.1167E-01 -2.4372E-01  3.1759E-02 -6.6719E-03  1.4691E-01  4.3735E-01 -1.9182E+00  7.5958E-02  2.2644E-02
             1.1392E-01
 GRADIENT:   4.4569E-02 -5.4154E+00 -3.8564E+00  2.2422E+00  3.2118E+00 -4.4523E+00 -3.7218E+00  1.1452E-01 -1.1540E+00  6.7565E-02
             1.2997E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1599.67463997338        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      925
 NPARAMETR:  9.5705E-01  1.1181E+00  7.0913E-01  9.3402E-01  8.9669E-01  1.0599E+00  1.4012E+00  1.3290E-01  9.7665E-01  9.2406E-01
             1.0110E+00
 PARAMETER:  5.6103E-02  2.1167E-01 -2.4372E-01  3.1744E-02 -9.0416E-03  1.5818E-01  4.3735E-01 -1.9182E+00  7.6376E-02  2.1017E-02
             1.1095E-01
 GRADIENT:   1.0949E-01 -4.2303E+00 -1.7257E+00  7.6636E-01 -1.3148E-01  1.1468E-02 -3.8120E+00  1.0983E-01 -1.1912E+00 -1.8639E-02
             1.4858E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1599.67910100368        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1105
 NPARAMETR:  9.5705E-01  1.1181E+00  7.0913E-01  9.3382E-01  8.9693E-01  1.0601E+00  1.4012E+00  1.3290E-01  9.8321E-01  9.2281E-01
             1.0108E+00
 PARAMETER:  5.6103E-02  2.1167E-01 -2.4372E-01  3.1527E-02 -8.7740E-03  1.5837E-01  4.3735E-01 -1.9182E+00  8.3069E-02  1.9673E-02
             1.1077E-01
 GRADIENT:   1.3360E-01 -4.6872E+00 -2.0903E+00  1.0786E+00  3.0274E-01  7.7288E-02 -3.5905E+00  1.1140E-01 -1.0285E-01 -8.6481E-02
            -1.3338E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1599.67934727629        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1169
 NPARAMETR:  9.5705E-01  1.1181E+00  7.0913E-01  9.3381E-01  8.9690E-01  1.0601E+00  1.4012E+00  1.3290E-01  9.8340E-01  9.2290E-01
             1.0108E+00
 PARAMETER:  5.6103E-02  2.1167E-01 -2.4372E-01  3.1519E-02 -8.8163E-03  1.5832E-01  4.3735E-01 -1.9182E+00  8.3263E-02  1.9769E-02
             1.1076E-01
 GRADIENT:   4.1154E+04  1.1436E+05  2.2830E+05 -3.9719E+05 -2.5600E+06 -1.7030E+06 -2.8659E+05 -9.6132E+03 -1.7237E+06  2.0512E+06
             2.0652E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1169
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.7904E-05  1.3059E-04 -6.0422E-03 -1.3422E-03 -1.6226E-02
 SE:             2.9833E-02  2.4162E-02  2.2030E-03  2.3226E-02  2.2057E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9925E-01  9.9569E-01  6.0949E-03  9.5392E-01  4.6195E-01

 ETASHRINKSD(%)  5.4623E-02  1.9056E+01  9.2620E+01  2.2190E+01  2.6106E+01
 ETASHRINKVR(%)  1.0922E-01  3.4480E+01  9.9455E+01  3.9456E+01  4.5397E+01
 EBVSHRINKSD(%)  3.9922E-01  1.9407E+01  9.3483E+01  2.3067E+01  2.4987E+01
 EBVSHRINKVR(%)  7.9684E-01  3.5047E+01  9.9575E+01  4.0813E+01  4.3730E+01
 RELATIVEINF(%)  9.9032E+01  6.4461E+00  6.2537E-02  5.7970E+00  7.1211E+00
 EPSSHRINKSD(%)  4.4136E+01
 EPSSHRINKVR(%)  6.8792E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1599.6793472762874     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -864.52852071254927     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1599.679       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.12E+00  7.09E-01  9.34E-01  8.97E-01  1.06E+00  1.40E+00  1.33E-01  9.83E-01  9.23E-01  1.01E+00
 


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
+       -8.88E+03
 
 TH 2
+       -4.50E+00  1.62E+08
 
 TH 3
+        1.59E+08  1.59E+08  2.98E+09
 
 TH 4
+       -5.40E+09 -1.47E+08  2.76E+09 -1.04E+04
 
 TH 5
+       -8.03E+09 -1.99E+09  9.88E+08 -2.38E+08  7.08E+03
 
 TH 6
+        1.23E+01  8.64E+08  7.79E+06  7.33E+08 -7.11E+08  1.42E+09
 
 TH 7
+       -5.68E+08 -7.35E+07  4.30E+08 -3.55E+08  1.95E+08  7.04E+08  7.42E+07
 
 TH 8
+       -2.51E+07 -1.01E+09  1.59E+09 -4.08E+08 -5.46E+08  1.48E+00 -2.21E+08  1.21E+09
 
 TH 9
+       -6.93E+09 -4.59E+09 -4.92E+08 -6.38E+08 -1.12E+09  1.64E+07 -7.98E+08 -5.95E+08  5.22E+03
 
 TH10
+       -3.08E+02 -2.49E+09  2.59E+08 -4.15E+09 -2.21E+09 -4.22E+09  5.80E+07 -2.37E+09 -3.32E+09  8.27E+08
 
 TH11
+        2.24E+08  4.32E+08  1.87E+09  2.78E+08  4.75E+09 -4.14E+03  3.98E+07 -1.07E+08 -3.64E+08  5.37E+09  1.19E+08
 
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
 #CPUT: Total CPU Time in Seconds,       26.782
Stop Time:
Sat Sep 25 10:34:05 CDT 2021
