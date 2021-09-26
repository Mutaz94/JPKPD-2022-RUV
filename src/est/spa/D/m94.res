Sat Sep 25 14:47:56 CDT 2021
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
$DATA ../../../../data/spa/D/dat94.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26349.2528545582        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.4097E+02  6.2068E+02  7.1332E+01  5.2514E+02 -7.3872E+01 -3.3234E+03 -1.2548E+03 -1.7235E+02 -1.8074E+03 -5.9777E+02
            -4.8542E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -434.440561679989        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1408E+00  6.7334E-01  6.4034E-01  1.5778E+00  2.9421E+00  2.1652E+00  8.8472E-01  8.8642E-01  5.4752E-01  7.5990E-01
             1.5188E+01
 PARAMETER:  2.3173E-01 -2.9551E-01 -3.4576E-01  5.5604E-01  1.1791E+00  8.7250E-01 -2.2487E-02 -2.0567E-02 -5.0236E-01 -1.7457E-01
             2.8205E+00
 GRADIENT:  -4.7990E+00  2.4958E+01 -1.2062E+01  1.9670E+01 -1.9686E+01  2.0640E+01  2.3716E-01  6.4939E+00  4.1663E+00  9.8687E-01
            -5.6660E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -459.554832217226        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0088E+00  2.9054E-01  3.2735E-01  1.9870E+00  5.4075E+00  2.2628E+00  8.8673E-02  3.0116E-02  5.8368E-02  4.9693E-01
             1.8637E+01
 PARAMETER:  1.0880E-01 -1.1360E+00 -1.0167E+00  7.8663E-01  1.7878E+00  9.1659E-01 -2.3228E+00 -3.4027E+00 -2.7410E+00 -5.9931E-01
             3.0251E+00
 GRADIENT:  -6.0278E+01  1.7332E+01 -1.4321E+01  1.5476E+02 -1.1096E+01 -6.0020E+01  8.9849E-02  2.6729E-02  1.3534E-01  1.0654E+00
             1.7070E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -518.135251768605        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.1296E-01  1.0916E-01  1.2790E-01  9.5913E-01  4.0491E+01  2.1904E+00  2.0595E-02  1.2224E-01  2.9233E-01  1.4238E-02
             1.6672E+01
 PARAMETER:  8.9410E-03 -2.1149E+00 -1.9565E+00  5.8271E-02  3.8011E+00  8.8408E-01 -3.7827E+00 -2.0018E+00 -1.1299E+00 -4.1518E+00
             2.9137E+00
 GRADIENT:   9.7280E+00  1.0842E+01  1.3562E+01 -1.7708E+01 -2.6275E-01 -2.3323E+00  3.5803E-03  2.7003E-01  5.0219E+00  1.2419E-06
             5.0788E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -533.081701178277        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  5.5996E-01  2.2689E-02  3.7243E-02  4.0918E-01  2.4236E+02  2.2341E+00  2.9222E-01  1.3760E-01  3.6753E-02  2.0630E-02
             1.3876E+01
 PARAMETER: -4.7989E-01 -3.6859E+00 -3.1903E+00 -7.9359E-01  5.5904E+00  9.0383E-01 -1.1302E+00 -1.8834E+00 -3.2035E+00 -3.7810E+00
             2.7302E+00
 GRADIENT:   9.1105E+00  1.1654E+00 -1.2073E+01  1.4999E+01  1.2999E-03  3.4082E-01  1.9723E-02  6.7278E-01  8.6603E-02  1.3831E-09
            -2.2261E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -533.253881149310        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  5.2655E-01  1.9218E-02  3.2983E-02  3.7171E-01  2.9700E+02  2.2084E+00  3.3688E-01  9.7256E-02  2.2581E-02  2.0734E-02
             1.3943E+01
 PARAMETER: -5.4140E-01 -3.8519E+00 -3.3118E+00 -8.8965E-01  5.7937E+00  8.9228E-01 -9.8802E-01 -2.2304E+00 -3.6906E+00 -3.7760E+00
             2.7350E+00
 GRADIENT:   2.1899E+00  1.2531E+00 -5.3201E+00  7.3199E+00 -9.6890E-04 -2.1744E+00  1.1496E-02  3.8260E-01  3.4076E-02  3.9140E-10
            -1.4319E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -533.254062058694        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  5.2504E-01  1.9004E-02  3.2738E-02  3.6958E-01  3.0020E+02  2.2084E+00  3.3785E-01  9.5130E-02  2.1943E-02  2.0582E-02
             1.3951E+01
 PARAMETER: -5.4429E-01 -3.8631E+00 -3.3192E+00 -8.9538E-01  5.8044E+00  8.9226E-01 -9.8514E-01 -2.2525E+00 -3.7193E+00 -3.7834E+00
             2.7356E+00
 GRADIENT:   2.0977E+00  1.2321E+00 -5.1273E+00  7.0211E+00 -9.4377E-04 -2.1025E+00  1.0816E-02  3.6829E-01  3.2219E-02  3.9064E-10
            -1.3810E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -533.254132347047        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  5.2398E-01  1.8846E-02  3.2567E-02  3.6810E-01  3.0229E+02  2.2084E+00  3.3855E-01  9.3665E-02  2.1508E-02  2.0430E-02
             1.3957E+01
 PARAMETER: -5.4630E-01 -3.8714E+00 -3.3244E+00 -8.9941E-01  5.8114E+00  8.9227E-01 -9.8308E-01 -2.2680E+00 -3.7393E+00 -3.7907E+00
             2.7360E+00
 GRADIENT:   2.0403E+00  1.2122E+00 -4.9942E+00  6.8185E+00 -9.1094E-04 -2.0460E+00  1.0291E-02  3.5854E-01  3.0981E-02  3.2990E-10
            -1.3450E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -533.254143717087        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      599
 NPARAMETR:  5.2314E-01  1.8714E-02  3.2431E-02  3.6691E-01  3.0388E+02  2.2084E+00  3.3919E-01  9.2482E-02  2.1164E-02  2.0281E-02
             1.3961E+01
 PARAMETER: -5.4791E-01 -3.8785E+00 -3.3286E+00 -9.0263E-01  5.8166E+00  8.9227E-01 -9.8121E-01 -2.2807E+00 -3.7554E+00 -3.7981E+00
             2.7363E+00
 GRADIENT:   1.9947E+00  1.1934E+00 -4.8864E+00  6.6588E+00 -8.7734E-04 -2.0006E+00  9.8529E-03  3.5073E-01  3.0022E-02  3.2926E-10
            -1.3160E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -534.120185446376        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  5.7324E-01  2.0981E-02  3.9985E-02  4.2776E-01  1.6078E+02  2.2397E+00  3.3624E-01  2.1069E-02  2.5513E-02  1.0000E-02
             1.4340E+01
 PARAMETER: -4.5646E-01 -3.7642E+00 -3.1193E+00 -7.4919E-01  5.1800E+00  9.0635E-01 -9.8992E-01 -3.7599E+00 -3.5686E+00 -4.5218E+00
             2.7630E+00
 GRADIENT:   2.9445E+00  4.4825E-01 -5.7538E-01 -1.7114E+00  3.3319E-03  1.0156E+00  4.3320E-03  1.6397E-02  4.3391E-02  0.0000E+00
            -3.2540E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -534.216388806393        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      953
 NPARAMETR:  5.7184E-01  1.2820E-02  4.0597E-02  4.3317E-01  6.1467E+01  2.2370E+00  9.9008E-01  1.0000E-02  1.4664E-02  1.0000E-02
             1.4423E+01
 PARAMETER: -4.5890E-01 -4.2567E+00 -3.1041E+00 -7.3662E-01  4.2185E+00  9.0513E-01  9.0028E-02 -8.0655E+00 -4.1224E+00 -6.7576E+00
             2.7688E+00
 GRADIENT:   3.5059E-01 -1.4702E-02 -3.2985E-01  2.6598E-01  1.4247E-02  4.6576E-03  6.8937E-04  0.0000E+00  1.4478E-02  0.0000E+00
             3.2066E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -534.267651623927        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1132
 NPARAMETR:  5.7488E-01  1.0000E-02  4.0976E-02  4.3683E-01  1.0489E+01  2.2404E+00  7.8495E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.4417E+01
 PARAMETER: -4.5359E-01 -5.1760E+00 -3.0948E+00 -7.2821E-01  2.4503E+00  9.0665E-01  2.1605E+00 -1.5974E+01 -5.0773E+00 -1.0866E+01
             2.7684E+00
 GRADIENT:   2.1987E-01  0.0000E+00 -2.1469E-02  2.3670E-02  2.7816E-02 -4.2775E-02  2.8140E-02  0.0000E+00  0.0000E+00  0.0000E+00
             9.3169E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -534.280866085584        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1310
 NPARAMETR:  5.7362E-01  1.0000E-02  4.0778E-02  4.3530E-01  1.0154E+01  2.2405E+00  2.2971E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.4410E+01
 PARAMETER: -4.5578E-01 -5.2077E+00 -3.0996E+00 -7.3173E-01  2.4179E+00  9.0670E-01  9.3164E-01 -1.6023E+01 -4.5110E+00 -1.0710E+01
             2.7679E+00
 GRADIENT:   2.1044E-02  0.0000E+00  1.6619E-02 -3.3402E-02  6.7574E-04 -2.5751E-02  2.2576E-03  0.0000E+00  0.0000E+00  0.0000E+00
            -1.5045E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -534.281009827550        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1492
 NPARAMETR:  5.7354E-01  1.0000E-02  4.0761E-02  4.3518E-01  1.0170E+01  2.2407E+00  1.9580E+00  1.0000E-02  1.0238E-02  1.0000E-02
             1.4412E+01
 PARAMETER: -4.5592E-01 -5.2063E+00 -3.1000E+00 -7.3200E-01  2.4194E+00  9.0679E-01  7.7192E-01 -1.5980E+01 -4.4816E+00 -1.0569E+01
             2.7680E+00
 GRADIENT:   1.1113E-02  0.0000E+00 -1.5684E-02  2.7658E-03  2.1661E-03  6.5598E-03  1.6383E-03  0.0000E+00  7.0750E-03  0.0000E+00
             5.1211E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -534.281958740424        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1669
 NPARAMETR:  5.7367E-01  1.0000E-02  4.0774E-02  4.3529E-01  1.0136E+01  2.2407E+00  3.7830E-01  1.0000E-02  1.0000E-02  1.0000E-02
             1.4412E+01
 PARAMETER: -4.5571E-01 -5.1730E+00 -3.0997E+00 -7.3175E-01  2.4161E+00  9.0681E-01 -8.7207E-01 -1.5369E+01 -4.5370E+00 -8.1671E+00
             2.7680E+00
 GRADIENT:   4.4957E-02  0.0000E+00 -2.4014E-02  6.4861E-03 -7.0542E-04  8.9722E-03  6.0767E-05  0.0000E+00  0.0000E+00  0.0000E+00
             3.6233E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -534.281993285631        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1844
 NPARAMETR:  5.7360E-01  1.0000E-02  4.0778E-02  4.3531E-01  1.0146E+01  2.2407E+00  6.8340E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.4410E+01
 PARAMETER: -4.5582E-01 -5.1366E+00 -3.0996E+00 -7.3170E-01  2.4171E+00  9.0678E-01 -2.5833E+00 -1.4714E+01 -4.5971E+00 -5.6507E+00
             2.7679E+00
 GRADIENT:   4.1871E-04  0.0000E+00 -3.9754E-04  2.3826E-04 -1.1677E-05 -1.1740E-04  1.9984E-06  0.0000E+00  0.0000E+00  0.0000E+00
             5.8791E-04

0ITERATION NO.:   80    OBJECTIVE VALUE:  -534.284163722465        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2025
 NPARAMETR:  5.7359E-01  1.0000E-02  4.0779E-02  4.3531E-01  1.0149E+01  2.2407E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.6226E+00
             1.4410E+01
 PARAMETER: -4.5584E-01 -5.0364E+00 -3.0996E+00 -7.3169E-01  2.4174E+00  9.0678E-01 -7.3735E+00 -1.2893E+01 -4.7689E+00  1.3872E+00
             2.7679E+00
 GRADIENT:  -1.1054E-02  0.0000E+00  1.3752E-02 -2.6688E-02 -1.1581E-02  2.2495E-02  0.0000E+00  0.0000E+00  0.0000E+00 -1.7939E-03
             6.7733E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -534.284428547754        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2187
 NPARAMETR:  5.7365E-01  1.0000E-02  4.0797E-02  4.3544E-01  1.0323E+01  2.2405E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.3142E+00
             1.4410E+01
 PARAMETER: -4.5573E-01 -5.0245E+00 -3.0991E+00 -7.3141E-01  2.4344E+00  9.0668E-01 -7.4717E+00 -1.2780E+01 -4.7688E+00  1.5619E+00
             2.7680E+00
 GRADIENT:   2.8303E-03  0.0000E+00 -1.5736E-03  6.8254E-04 -1.1347E-05 -1.1738E-04  0.0000E+00  0.0000E+00  0.0000E+00  1.2600E-05
             2.0653E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2187
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8870E-03  9.8038E-07  1.2614E-04 -2.6620E-04 -3.0418E-04
 SE:             2.9081E-02  4.9711E-07  2.1270E-04  2.7489E-04  5.4092E-04
 N:                     100         100         100         100         100

 P VAL.:         9.2092E-01  4.8594E-02  5.5317E-01  3.3286E-01  5.7388E-01

 ETASHRINKSD(%)  2.5763E+00  9.9998E+01  9.9287E+01  9.9079E+01  9.8188E+01
 ETASHRINKVR(%)  5.0863E+00  1.0000E+02  9.9995E+01  9.9992E+01  9.9967E+01
 EBVSHRINKSD(%)  2.5912E+00  9.9998E+01  9.9182E+01  9.8935E+01  9.8273E+01
 EBVSHRINKVR(%)  5.1153E+00  1.0000E+02  9.9993E+01  9.9989E+01  9.9970E+01
 RELATIVEINF(%)  6.9017E+00  3.3581E-09  6.0501E-05  7.7610E-05  3.4859E-03
 EPSSHRINKSD(%)  4.8662E+00
 EPSSHRINKVR(%)  9.4956E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -534.28442854775437     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       200.86639801598380     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.88
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
 





 #OBJV:********************************************     -534.284       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.74E-01  1.00E-02  4.08E-02  4.35E-01  1.03E+01  2.24E+00  1.00E-02  1.00E-02  1.00E-02  4.31E+00  1.44E+01
 


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
+        6.32E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.15E+03  0.00E+00  2.66E+05
 
 TH 4
+       -1.18E+02  0.00E+00 -3.24E+04  4.28E+03
 
 TH 5
+        1.75E-01  0.00E+00 -4.30E+00  5.06E-01  8.24E-03
 
 TH 6
+        1.37E+00  0.00E+00  3.90E+02 -5.48E+01  2.77E-03  3.48E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.96E-03  0.00E+00  3.32E-02 -1.91E-02 -1.69E-03  1.31E-03  0.00E+00  0.00E+00  0.00E+00  5.93E-04
 
 TH11
+       -9.05E+00  0.00E+00  1.39E+02 -1.46E+01 -5.12E-03  4.68E-01  0.00E+00  0.00E+00  0.00E+00  1.51E-04  1.60E+00
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       33.820
Stop Time:
Sat Sep 25 14:48:44 CDT 2021
