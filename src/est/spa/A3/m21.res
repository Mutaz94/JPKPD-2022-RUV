Sat Sep 18 10:19:43 CDT 2021
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
$DATA ../../../../data/spa/A3/dat21.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -120.156912792149        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.2683E+01  3.4436E+01  1.3990E+02 -1.5572E+02  9.2480E+01  1.7644E+01 -6.9652E+01 -5.7935E+01 -1.7140E+02 -1.4322E+02
            -2.6863E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1252.74837273607        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0973E+00  9.8556E-01  8.8253E-01  1.2387E+00  9.9432E-01  8.3478E-01  9.7230E-01  1.0147E+00  1.0605E+00  9.9976E-01
             5.2840E+00
 PARAMETER:  1.9281E-01  8.5456E-02 -2.4961E-02  3.1406E-01  9.4300E-02 -8.0583E-02  7.1910E-02  1.1462E-01  1.5874E-01  9.9759E-02
             1.7647E+00
 GRADIENT:   7.8952E+00 -7.2586E+00 -2.4578E+01  1.8092E+01  1.7193E+00 -8.8745E+00  1.0870E+01  8.1993E+00  2.6066E+01  1.9323E+01
             1.9489E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1275.52765691516        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0894E+00  9.4950E-01  4.2987E-01  1.1545E+00  5.4192E-01  9.1284E-01  2.8590E-01  2.6931E-01  1.3418E+00  4.3452E-01
             4.6436E+00
 PARAMETER:  1.8559E-01  4.8182E-02 -7.4427E-01  2.4363E-01 -5.1263E-01  8.8027E-03 -1.1521E+00 -1.2119E+00  3.9404E-01 -7.3350E-01
             1.6355E+00
 GRADIENT:   4.3901E+00  8.9693E+01  5.4849E+01  3.6929E+01 -1.2737E+02  6.4335E+00  3.5608E-01  4.5538E-01  3.3053E+01  3.8806E+00
             1.2702E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1305.01924236901        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0586E+00  5.8214E-01  3.6736E-01  1.2604E+00  4.1893E-01  9.0827E-01  4.3748E-01  1.4073E-01  1.0077E+00  2.4022E-01
             3.6908E+00
 PARAMETER:  1.5697E-01 -4.4104E-01 -9.0140E-01  3.3145E-01 -7.7006E-01  3.7919E-03 -7.2673E-01 -1.8609E+00  1.0771E-01 -1.3262E+00
             1.4058E+00
 GRADIENT:  -1.3309E+01  3.1517E+01  2.6902E+01  3.3317E+01 -4.3194E+01 -4.7090E+00 -1.5144E+00 -2.0686E-01 -4.8209E+00 -1.8260E+00
            -1.4545E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1308.75888570177        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0652E+00  4.3130E-01  1.8469E-01  1.1172E+00  2.6206E-01  9.2917E-01  8.9796E-01  4.8284E-02  1.1658E+00  1.2783E-01
             3.5483E+00
 PARAMETER:  1.6316E-01 -7.4095E-01 -1.5891E+00  2.1081E-01 -1.2392E+00  2.6538E-02 -7.6349E-03 -2.9307E+00  2.5341E-01 -1.9571E+00
             1.3665E+00
 GRADIENT:   2.3858E+01 -9.5644E+00 -1.5891E+01  2.2140E+01  1.4232E+01 -8.3238E+00 -6.5893E-03 -7.9102E-02 -1.2369E+01 -1.1807E+00
             1.2548E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1311.50262482747        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0390E+00  4.9783E-01  1.6898E-01  1.0585E+00  2.6559E-01  9.5354E-01  8.9643E-01  6.4537E-02  1.2949E+00  1.9792E-01
             3.3208E+00
 PARAMETER:  1.3831E-01 -5.9750E-01 -1.6780E+00  1.5688E-01 -1.2258E+00  5.2430E-02 -9.3298E-03 -2.6405E+00  3.5844E-01 -1.5199E+00
             1.3002E+00
 GRADIENT:  -1.2134E+01  2.7131E+00  8.7290E-02  5.5424E+00  5.5424E-01 -6.9706E-02  8.9503E-01 -1.2382E-01 -1.9953E+00 -1.9096E+00
             9.7602E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1313.74823017665        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0513E+00  4.4069E-01  1.5763E-01  1.0458E+00  2.4231E-01  9.5299E-01  5.7907E-01  4.9792E-02  1.3959E+00  5.0190E-01
             3.1923E+00
 PARAMETER:  1.5006E-01 -7.1940E-01 -1.7475E+00  1.4482E-01 -1.3176E+00  5.1850E-02 -4.4633E-01 -2.8999E+00  4.3357E-01 -5.8935E-01
             1.2607E+00
 GRADIENT:   1.4248E+01 -3.0202E+00 -5.9409E+00 -1.8413E+00  2.2528E+01 -6.2579E-01  1.7502E+00 -3.3364E-02  4.6053E+00  2.0462E+00
             1.8104E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1313.80130328007        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      519
 NPARAMETR:  1.0502E+00  4.1793E-01  1.5247E-01  1.0415E+00  2.3299E-01  9.5445E-01  5.3239E-01  4.5420E-02  1.4084E+00  5.1958E-01
             3.1868E+00
 PARAMETER:  1.4899E-01 -7.7244E-01 -1.7808E+00  1.4071E-01 -1.3568E+00  5.3376E-02 -5.3038E-01 -2.9918E+00  4.4244E-01 -5.5473E-01
             1.2590E+00
 GRADIENT:   9.6522E+00 -2.1773E+00 -3.9915E+00 -6.3598E-01  1.4636E+01 -3.8251E-01  1.5219E+00 -2.8916E-02  3.5835E+00  1.2526E+00
             2.9242E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1313.80222965787        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      594
 NPARAMETR:  1.0498E+00  4.1406E-01  1.5112E-01  1.0395E+00  2.3110E-01  9.5472E-01  5.1306E-01  4.4994E-02  1.4119E+00  5.2567E-01
             3.1829E+00
 PARAMETER:  1.4863E-01 -7.8174E-01 -1.7897E+00  1.3873E-01 -1.3649E+00  5.3662E-02 -5.6737E-01 -3.0012E+00  4.4497E-01 -5.4308E-01
             1.2578E+00
 GRADIENT:   8.7312E+00 -2.0692E+00 -3.6792E+00 -6.0418E-01  1.3472E+01 -3.2440E-01  1.4356E+00 -2.8445E-02  3.2821E+00  1.1828E+00
             2.8723E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1313.80259341256        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  1.0496E+00  4.1193E-01  1.5025E-01  1.0381E+00  2.2998E-01  9.5486E-01  4.9959E-01  4.4912E-02  1.4144E+00  5.2947E-01
             3.1802E+00
 PARAMETER:  1.4840E-01 -7.8689E-01 -1.7954E+00  1.3737E-01 -1.3698E+00  5.3811E-02 -5.9397E-01 -3.0030E+00  4.4673E-01 -5.3587E-01
             1.2569E+00
 GRADIENT:   8.2388E+00 -1.9799E+00 -3.4898E+00 -5.8246E-01  1.2778E+01 -2.9918E-01  1.3742E+00 -2.8407E-02  3.1092E+00  1.1318E+00
             2.7832E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1313.80275064287        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      744
 NPARAMETR:  1.0494E+00  4.1065E-01  1.4969E-01  1.0371E+00  2.2928E-01  9.5494E-01  4.9046E-01  4.4966E-02  1.4161E+00  5.3192E-01
             3.1783E+00
 PARAMETER:  1.4826E-01 -7.9000E-01 -1.7992E+00  1.3647E-01 -1.3728E+00  5.3896E-02 -6.1241E-01 -3.0019E+00  4.4792E-01 -5.3127E-01
             1.2564E+00
 GRADIENT:   7.9479E+00 -1.9165E+00 -3.3701E+00 -5.6748E-01  1.2345E+01 -2.8616E-01  1.3327E+00 -2.8525E-02  3.0034E+00  1.0977E+00
             2.7101E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1313.80280840945        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  1.0493E+00  4.0973E-01  1.4926E-01  1.0364E+00  2.2876E-01  9.5500E-01  4.8337E-01  4.5089E-02  1.4174E+00  5.3375E-01
             3.1769E+00
 PARAMETER:  1.4815E-01 -7.9227E-01 -1.8020E+00  1.3576E-01 -1.3751E+00  5.3956E-02 -6.2698E-01 -2.9991E+00  4.4883E-01 -5.2782E-01
             1.2559E+00
 GRADIENT:   7.7386E+00 -1.8661E+00 -3.2802E+00 -5.5599E-01  1.2023E+01 -2.7752E-01  1.3007E+00 -2.8712E-02  2.9255E+00  1.0714E+00
             2.6488E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1313.80282440378        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  1.0493E+00  4.0920E-01  1.4902E-01  1.0360E+00  2.2846E-01  9.5503E-01  4.7919E-01  4.5210E-02  1.4182E+00  5.3482E-01
             3.1761E+00
 PARAMETER:  1.4809E-01 -7.9355E-01 -1.8037E+00  1.3535E-01 -1.3764E+00  5.3990E-02 -6.3567E-01 -2.9964E+00  4.4937E-01 -5.2583E-01
             1.2556E+00
 GRADIENT:   7.6203E+00 -1.8360E+00 -3.2280E+00 -5.4946E-01  1.1837E+01 -2.7292E-01  1.2819E+00 -2.8891E-02  2.8809E+00  1.0558E+00
             2.6112E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1313.80283590959        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  1.0492E+00  4.0876E-01  1.4881E-01  1.0356E+00  2.2821E-01  9.5506E-01  4.7561E-01  4.5357E-02  1.4188E+00  5.3571E-01
             3.1754E+00
 PARAMETER:  1.4804E-01 -7.9462E-01 -1.8051E+00  1.3500E-01 -1.3775E+00  5.4017E-02 -6.4316E-01 -2.9932E+00  4.4983E-01 -5.2416E-01
             1.2554E+00
 GRADIENT:   7.5218E+00 -1.8099E+00 -3.1836E+00 -5.4402E-01  1.1681E+01 -2.6909E-01  1.2659E+00 -2.9104E-02  2.8435E+00  1.0426E+00
             2.5783E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1314.04224469084        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:     1103
 NPARAMETR:  1.0492E+00  4.1474E-01  1.5540E-01  1.0481E+00  2.3413E-01  9.5235E-01  4.5302E-01  4.9511E-02  1.3813E+00  5.2076E-01
             3.2168E+00
 PARAMETER:  1.4799E-01 -7.8011E-01 -1.7617E+00  1.4697E-01 -1.3519E+00  5.1177E-02 -6.9181E-01 -2.9056E+00  4.2303E-01 -5.5247E-01
             1.2684E+00
 GRADIENT:  -2.1178E-01 -2.9229E+00 -4.7164E+00 -6.0666E-01  4.2414E+00 -9.8638E-01  8.1826E-01 -3.3163E-02  1.7118E-01  6.9073E-01
             5.0595E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1315.11175911264        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1281
 NPARAMETR:  1.0496E+00  4.2887E-01  1.6286E-01  1.0543E+00  2.4095E-01  9.6087E-01  8.5020E-02  4.9102E-01  1.3655E+00  5.6545E-01
             3.1489E+00
 PARAMETER:  1.4845E-01 -7.4661E-01 -1.7149E+00  1.5289E-01 -1.3232E+00  6.0081E-02 -2.3649E+00 -6.1128E-01  4.1155E-01 -4.7013E-01
             1.2471E+00
 GRADIENT:   1.4004E-01  2.3899E+00  6.6389E+00 -3.1650E+00  2.8348E-01  3.8762E+00  2.6067E-02 -1.1246E+00 -6.5623E-01  7.4168E+00
             2.6845E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1318.19085397539        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1461
 NPARAMETR:  1.0384E+00  4.3415E-01  1.3912E-01  1.0159E+00  2.2829E-01  9.5202E-01  7.3028E-02  1.1196E+00  1.4852E+00  4.2530E-01
             3.0065E+00
 PARAMETER:  1.3765E-01 -7.3437E-01 -1.8724E+00  1.1576E-01 -1.3772E+00  5.0828E-02 -2.5169E+00  2.1299E-01  4.9552E-01 -7.5495E-01
             1.2008E+00
 GRADIENT:  -2.4945E+00  1.0602E+00  6.8526E-01  6.0087E-01 -2.4678E+00  3.9465E-01 -9.0266E-03  4.6661E-02 -2.9518E-01  2.6916E-01
             1.0867E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1318.19797327408        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1636
 NPARAMETR:  1.0393E+00  4.3533E-01  1.3976E-01  1.0161E+00  2.2918E-01  9.5087E-01  7.8750E-02  1.1321E+00  1.4861E+00  4.1531E-01
             3.0018E+00
 PARAMETER:  1.3854E-01 -7.3164E-01 -1.8678E+00  1.1599E-01 -1.3733E+00  4.9621E-02 -2.4415E+00  2.2409E-01  4.9618E-01 -7.7872E-01
             1.1992E+00
 GRADIENT:  -1.0679E-01  1.9687E-01  2.8811E-01 -2.0539E-01 -5.5201E-01  6.1865E-02 -1.5573E-02 -1.3335E-03 -5.0228E-02 -7.8296E-02
            -3.3458E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1318.44123082312        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1827
 NPARAMETR:  1.0411E+00  4.6465E-01  1.4532E-01  1.0173E+00  2.4042E-01  9.5116E-01  4.4478E-01  1.1625E+00  1.4615E+00  3.3027E-01
             3.0399E+00
 PARAMETER:  1.4027E-01 -6.6647E-01 -1.8288E+00  1.1718E-01 -1.3254E+00  4.9923E-02 -7.1017E-01  2.5058E-01  4.7947E-01 -1.0078E+00
             1.2118E+00
 GRADIENT:   5.1327E+00 -2.5022E+00 -1.7344E+00 -4.5374E+00  7.6857E+00  4.7160E-01 -3.6621E-01  4.5919E-01  1.0512E+00  9.6001E-01
             4.4412E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1319.08863980631        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2003
 NPARAMETR:  1.0386E+00  4.8368E-01  1.5248E-01  1.0371E+00  2.4938E-01  9.5114E-01  7.4388E-01  1.1843E+00  1.4179E+00  1.5580E-01
             3.0201E+00
 PARAMETER:  1.3789E-01 -6.2632E-01 -1.7807E+00  1.3640E-01 -1.2888E+00  4.9909E-02 -1.9588E-01  2.6916E-01  4.4919E-01 -1.7592E+00
             1.2053E+00
 GRADIENT:  -2.0984E+00  9.0450E-01 -7.9476E-01  5.4361E+00  2.6238E-01  3.0078E-02  7.2928E-01 -6.1022E-02  2.0320E-01  4.8345E-01
             6.1233E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1319.30048421382        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2180
 NPARAMETR:  1.0385E+00  4.8801E-01  1.5245E-01  1.0306E+00  2.5074E-01  9.5062E-01  7.5148E-01  1.2226E+00  1.4196E+00  4.9221E-02
             3.0164E+00
 PARAMETER:  1.3779E-01 -6.1742E-01 -1.7809E+00  1.3014E-01 -1.2833E+00  4.9364E-02 -1.8571E-01  3.0098E-01  4.5034E-01 -2.9114E+00
             1.2041E+00
 GRADIENT:  -7.2664E-01 -5.2174E-01 -4.6473E-01  8.0782E-01  1.0868E+00 -2.2313E-02  1.2523E-01  2.2344E-01  3.2628E-02  4.2490E-02
             3.9334E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1319.32103851683        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2355
 NPARAMETR:  1.0388E+00  4.8887E-01  1.5254E-01  1.0295E+00  2.5095E-01  9.5069E-01  7.5266E-01  1.2217E+00  1.4191E+00  1.0000E-02
             3.0158E+00
 PARAMETER:  1.3805E-01 -6.1566E-01 -1.7804E+00  1.2907E-01 -1.2825E+00  4.9431E-02 -1.8415E-01  3.0028E-01  4.5004E-01 -4.5412E+00
             1.2039E+00
 GRADIENT:   5.7103E-03 -2.7029E-03  1.2270E-02 -2.6518E-02 -2.7107E-02 -1.8017E-03  1.8419E-03  1.2208E-02 -1.1810E-02  0.0000E+00
            -1.5270E-02

0ITERATION NO.:  107    OBJECTIVE VALUE:  -1319.32104187890        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2412
 NPARAMETR:  1.0388E+00  4.8892E-01  1.5255E-01  1.0295E+00  2.5098E-01  9.5070E-01  7.5261E-01  1.2213E+00  1.4191E+00  1.0000E-02
             3.0161E+00
 PARAMETER:  1.3805E-01 -6.1555E-01 -1.7803E+00  1.2912E-01 -1.2824E+00  4.9439E-02 -1.8421E-01  2.9991E-01  4.5003E-01 -4.5391E+00
             1.2040E+00
 GRADIENT:  -2.8871E-03  7.3905E-03  6.9087E-03  2.6452E-03 -8.8442E-03  1.2785E-03 -4.6656E-04 -3.5906E-03  2.7816E-04  0.0000E+00
            -8.4848E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2412
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1231E-03 -1.1737E-02  1.0515E-02 -8.5362E-03  5.0094E-04
 SE:             2.8754E-02  1.4075E-02  1.8350E-02  2.6222E-02  3.7952E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4114E-01  4.0434E-01  5.6663E-01  7.4478E-01  1.8685E-01

 ETASHRINKSD(%)  3.6694E+00  5.2846E+01  3.8526E+01  1.2152E+01  9.8729E+01
 ETASHRINKVR(%)  7.2042E+00  7.7765E+01  6.2210E+01  2.2827E+01  9.9984E+01
 EBVSHRINKSD(%)  3.6118E+00  5.2065E+01  3.8198E+01  1.0512E+01  9.8768E+01
 EBVSHRINKVR(%)  7.0932E+00  7.7022E+01  6.1805E+01  1.9920E+01  9.9985E+01
 RELATIVEINF(%)  8.7126E+01  1.1340E+00  6.6690E+00  5.3100E+01  4.9817E-04
 EPSSHRINKSD(%)  3.2507E+01
 EPSSHRINKVR(%)  5.4447E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1319.3210418789049     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -584.17021531516673     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.82
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1319.321       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  4.89E-01  1.53E-01  1.03E+00  2.51E-01  9.51E-01  7.53E-01  1.22E+00  1.42E+00  1.00E-02  3.02E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.41E+01
 
 TH 2
+       -9.65E+01  1.92E+03
 
 TH 3
+       -3.98E+02  3.62E+03  1.19E+04
 
 TH 4
+        2.35E+01  1.69E+02 -4.31E+02  1.26E+02
 
 TH 5
+        4.28E+02 -7.87E+03 -1.56E+04 -5.84E+02  3.24E+04
 
 TH 6
+       -1.61E+00 -2.75E+01  1.82E+01 -1.28E+01  1.03E+02  1.37E+00
 
 TH 7
+        2.28E+00 -4.21E+01 -8.34E+01 -3.16E+00  1.74E+02  5.53E-01  9.29E-01
 
 TH 8
+       -4.76E-01 -2.64E+01 -7.55E+00 -8.62E+00  1.02E+02  9.67E-01  5.49E-01  7.19E-01
 
 TH 9
+       -4.19E+00 -6.76E+01  5.04E+01 -3.24E+01  2.52E+02  3.44E+00  1.35E+00  2.42E+00  8.66E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        4.74E-01 -9.19E+00 -1.76E+01 -7.71E-01  3.78E+01  1.28E-01  2.02E-01  1.24E-01  3.15E-01  0.00E+00  4.41E-02
 
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
+        1.06E+03
 
 TH 2
+       -1.96E+01  1.96E+03
 
 TH 3
+       -4.08E+02  3.47E+03  1.20E+04
 
 TH 4
+       -3.03E+01  1.70E+02 -4.48E+02  4.16E+02
 
 TH 5
+        3.72E+02 -7.04E+03 -1.49E+04 -3.44E+02  2.95E+04
 
 TH 6
+       -1.57E+00 -1.51E+01  2.40E+01 -1.12E+01  8.08E+01  1.91E+02
 
 TH 7
+       -1.77E+00 -2.48E+01 -7.95E+01 -3.65E+00  1.61E+02 -8.67E-01  2.27E+01
 
 TH 8
+        1.30E+00 -2.65E+01 -2.34E-01 -9.13E-01  8.54E+01  1.81E+00  6.30E+00  1.97E+01
 
 TH 9
+        1.11E+01 -3.46E+01  7.51E+01 -6.12E+00  2.07E+02  2.06E+00  7.98E+00 -2.43E+00  5.75E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.80E+01 -1.01E+01 -1.95E+01 -3.11E+00  3.27E+01  2.18E+00  7.24E+00  1.03E+01  6.21E+00  0.00E+00  3.00E+01
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.06E+03
 
 TH 2
+       -3.36E+02  1.90E+03
 
 TH 3
+       -9.23E+02  3.33E+03  1.23E+04
 
 TH 4
+       -1.83E+02  5.10E+02  8.37E+00  4.18E+02
 
 TH 5
+        1.55E+03 -6.69E+03 -1.45E+04 -1.44E+03  2.60E+04
 
 TH 6
+        4.10E+00  2.31E+01 -1.84E+01  2.98E-02  2.32E+02  1.93E+02
 
 TH 7
+        2.04E+00 -3.09E+01 -2.46E+01 -4.32E+00  1.35E+02 -3.78E+00  1.62E+01
 
 TH 8
+        1.77E+01 -4.44E+01 -1.86E+02 -3.58E+00  2.22E+02  4.07E+00  7.62E+00  1.92E+01
 
 TH 9
+        6.12E+01 -6.58E+01  2.74E+02 -7.71E+01  1.06E+02 -9.82E-02  6.91E+00  6.24E+00  9.20E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        5.62E+01 -9.37E+01  2.04E+01 -8.07E+01  2.92E+02  2.05E+01  7.63E+00  1.72E+01  5.90E+01  0.00E+00  8.54E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       33.384
Stop Time:
Sat Sep 18 10:20:18 CDT 2021
