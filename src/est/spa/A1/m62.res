Sat Sep 25 08:09:28 CDT 2021
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
$DATA ../../../../data/spa/A1/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1367.35844100540        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.6007E+01  1.7776E+01  2.5845E+01  1.0552E+01  9.3184E+01  2.4541E+01 -3.1159E+01 -1.8397E+01 -3.5765E+01 -5.8157E+01
            -4.9617E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1501.64863210280        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0434E+00  9.2967E-01  8.9792E-01  1.0533E+00  8.8363E-01  8.8257E-01  1.1150E+00  1.0074E+00  1.0820E+00  1.0655E+00
             1.8183E+00
 PARAMETER:  1.4245E-01  2.7069E-02 -7.6701E-03  1.5196E-01 -2.3722E-02 -2.4914E-02  2.0883E-01  1.0735E-01  1.7881E-01  1.6348E-01
             6.9789E-01
 GRADIENT:   4.0173E+01  5.5966E+00 -5.7639E+00  1.5600E+01  1.3505E+01 -2.0459E+01 -3.3756E+00  4.7498E+00  8.3469E+00  4.1182E+00
             1.0455E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1507.24411112395        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0489E+00  7.3079E-01  5.2640E-01  1.1611E+00  5.8889E-01  9.2405E-01  1.4229E+00  5.3492E-01  9.6185E-01  7.1881E-01
             1.7642E+00
 PARAMETER:  1.4775E-01 -2.1364E-01 -5.4169E-01  2.4938E-01 -4.2952E-01  2.1010E-02  4.5273E-01 -5.2565E-01  6.1101E-02 -2.3015E-01
             6.6770E-01
 GRADIENT:   4.0501E+01  1.5620E+01 -5.0751E+01  1.0787E+02  6.7704E+01 -5.2326E+00  9.0738E-02  4.5724E+00  8.9289E+00 -1.6642E-01
            -2.4112E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1512.68006952065        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0272E+00  6.1491E-01  3.5823E-01  1.0967E+00  4.2560E-01  9.5486E-01  1.4917E+00  2.8734E-01  8.5532E-01  5.0339E-01
             1.7490E+00
 PARAMETER:  1.2687E-01 -3.8628E-01 -9.2658E-01  1.9233E-01 -7.5426E-01  5.3807E-02  4.9993E-01 -1.1471E+00 -5.6285E-02 -5.8638E-01
             6.5905E-01
 GRADIENT:  -1.4387E+01  1.7924E+01  9.7439E+00  1.7068E+01 -1.5100E+01  5.7456E+00 -2.3714E-01  1.5324E-01 -7.0508E+00 -5.9281E+00
            -4.5899E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1515.50700658241        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      343
 NPARAMETR:  1.0358E+00  4.4883E-01  4.3500E-01  1.2036E+00  4.2669E-01  9.3661E-01  1.8251E+00  1.4392E-01  8.5587E-01  7.1746E-01
             1.7404E+00
 PARAMETER:  1.3519E-01 -7.0112E-01 -7.3241E-01  2.8529E-01 -7.5170E-01  3.4509E-02  7.0163E-01 -1.8385E+00 -5.5633E-02 -2.3204E-01
             6.5410E-01
 GRADIENT:  -8.6333E+00  9.1422E+00  1.9312E+01 -1.2351E+01 -3.8043E+01 -1.5552E-01 -1.0373E+00  2.1026E-01 -1.2480E+00  3.1302E+00
             7.6591E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1516.80041425733        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.0353E+00  3.5155E-01  5.0518E-01  1.2811E+00  4.5898E-01  9.3015E-01  2.2840E+00  8.4576E-02  8.3777E-01  7.7780E-01
             1.7470E+00
 PARAMETER:  1.3467E-01 -9.4541E-01 -5.8285E-01  3.4769E-01 -6.7876E-01  2.7587E-02  9.2592E-01 -2.3701E+00 -7.7016E-02 -1.5129E-01
             6.5792E-01
 GRADIENT:  -7.5537E-01  8.9663E-02 -4.7114E+00  3.1724E-01  4.8787E+00 -9.0506E-01  1.0569E+00  3.4796E-02  6.5822E-01 -8.2490E-01
             1.7018E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1517.19363310840        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      695
 NPARAMETR:  1.0322E+00  2.5472E-01  5.4120E-01  1.3432E+00  4.5909E-01  9.3104E-01  2.5064E+00  3.5635E-02  8.3646E-01  8.6073E-01
             1.7399E+00
 PARAMETER:  1.3172E-01 -1.2676E+00 -5.1397E-01  3.9508E-01 -6.7851E-01  2.8543E-02  1.0189E+00 -3.2344E+00 -7.8580E-02 -4.9971E-02
             6.5382E-01
 GRADIENT:   4.1613E-01  1.7543E+00  2.2653E+00  2.3983E+00 -2.3034E+00  8.3698E-01  4.3344E-01  6.5007E-03 -1.2218E+00 -1.5097E+00
            -3.4746E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1519.34955610449        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      871
 NPARAMETR:  1.0225E+00  1.0970E-01  5.5076E-01  1.4149E+00  4.3951E-01  9.2205E-01  4.0218E+00  1.0000E-02  8.1659E-01  8.9458E-01
             1.7293E+00
 PARAMETER:  1.2226E-01 -2.1100E+00 -4.9646E-01  4.4707E-01 -7.2209E-01  1.8843E-02  1.4917E+00 -5.7658E+00 -1.0262E-01 -1.1402E-02
             6.4769E-01
 GRADIENT:  -7.9824E+00  3.9651E+00  1.2039E+01  1.0436E+01 -1.9751E+01 -1.3713E+00  2.8331E+00  0.0000E+00  1.8387E-01 -1.0350E+00
            -3.4039E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1520.42574260332        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  1.0230E+00  4.8193E-02  5.5668E-01  1.4460E+00  4.4038E-01  9.2669E-01  5.7858E+00  1.0000E-02  8.0233E-01  8.6615E-01
             1.7483E+00
 PARAMETER:  1.2278E-01 -2.9325E+00 -4.8576E-01  4.6878E-01 -7.2011E-01  2.3864E-02  1.8554E+00 -8.4161E+00 -1.2024E-01 -4.3691E-02
             6.5864E-01
 GRADIENT:   5.3341E-02  5.2285E-01  1.3831E-01  2.3989E+01 -5.9519E-01  1.3249E+00 -2.2464E+00  0.0000E+00  4.9307E-01 -5.9199E+00
            -4.2701E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1522.00881925102        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1226
 NPARAMETR:  1.0207E+00  1.0000E-02  5.5770E-01  1.4509E+00  4.3671E-01  9.2140E-01  1.2106E+01  1.0000E-02  7.8956E-01  8.8608E-01
             1.7420E+00
 PARAMETER:  1.2053E-01 -4.5916E+00 -4.8394E-01  4.7221E-01 -7.2848E-01  1.8135E-02  2.5937E+00 -1.3985E+01 -1.3628E-01 -2.0951E-02
             6.5501E-01
 GRADIENT:  -9.0991E-01  0.0000E+00  4.1518E+00  1.0263E+01 -5.5987E+00  8.9990E-02 -3.1840E+00  0.0000E+00  5.1963E-01  8.5758E-01
             7.5859E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1522.09523784260        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1404
 NPARAMETR:  1.0210E+00  1.0000E-02  5.4741E-01  1.4436E+00  4.3156E-01  9.2165E-01  1.2285E+01  1.0000E-02  7.9100E-01  8.8297E-01
             1.7385E+00
 PARAMETER:  1.2082E-01 -4.6469E+00 -5.0256E-01  4.6714E-01 -7.4036E-01  1.8414E-02  2.6084E+00 -1.4199E+01 -1.3445E-01 -2.4463E-02
             6.5301E-01
 GRADIENT:   1.6131E-01  0.0000E+00 -2.2384E-01  5.3015E-01  5.6064E-01 -2.0122E-02  8.2799E-01  0.0000E+00 -3.3260E-01 -5.4164E-01
            -1.1388E-01

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1522.09617475713        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1496
 NPARAMETR:  1.0210E+00  1.0000E-02  5.4514E-01  1.4420E+00  4.3026E-01  9.2162E-01  1.2275E+01  1.0000E-02  7.9141E-01  8.8265E-01
             1.7376E+00
 PARAMETER:  1.2078E-01 -4.6485E+00 -5.0672E-01  4.6601E-01 -7.4336E-01  1.8380E-02  2.6076E+00 -1.4210E+01 -1.3393E-01 -2.4832E-02
             6.5252E-01
 GRADIENT:   9.1380E-03  0.0000E+00 -7.0091E-02 -6.9407E-02  1.0580E-01  7.8786E-04  4.8420E-02  0.0000E+00 -2.3362E-02 -2.3830E-02
            -1.7802E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1496
 NO. OF SIG. DIGITS IN FINAL EST.:  4.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9952E-05  1.0255E-02 -1.1258E-04 -8.3225E-03 -8.5957E-03
 SE:             2.9524E-02  6.0767E-03  2.2508E-04  2.7923E-02  2.4949E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9919E-01  9.1480E-02  6.1696E-01  7.6567E-01  7.3045E-01

 ETASHRINKSD(%)  1.0917E+00  7.9642E+01  9.9246E+01  6.4529E+00  1.6417E+01
 ETASHRINKVR(%)  2.1715E+00  9.5856E+01  9.9994E+01  1.2489E+01  3.0139E+01
 EBVSHRINKSD(%)  1.3385E+00  8.4473E+01  9.9255E+01  5.7760E+00  1.4759E+01
 EBVSHRINKVR(%)  2.6591E+00  9.7589E+01  9.9994E+01  1.1218E+01  2.7339E+01
 RELATIVEINF(%)  9.6979E+01  2.0122E+00  3.0020E-04  3.5313E+01  4.1334E+00
 EPSSHRINKSD(%)  3.9346E+01
 EPSSHRINKVR(%)  6.3211E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1522.0961747571300     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -786.94534819339185     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1522.096       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E-02  5.45E-01  1.44E+00  4.30E-01  9.22E-01  1.23E+01  1.00E-02  7.91E-01  8.83E-01  1.74E+00
 


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
+        1.22E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -7.57E+06  0.00E+00  3.38E+06
 
 TH 4
+        1.09E+02  0.00E+00 -2.59E+03  5.71E+05
 
 TH 5
+       -2.44E+02  0.00E+00  1.59E+03 -1.20E+06  2.52E+06
 
 TH 6
+       -2.30E+01  0.00E+00  1.15E+03 -1.21E+02  2.34E+02  1.69E+02
 
 TH 7
+       -2.32E+00  0.00E+00  4.26E+01  3.20E+01 -9.17E+01  1.99E+00  2.51E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.17E+02  0.00E+00  8.81E+06 -4.55E+02  9.03E+02  1.16E+02  6.42E+00  0.00E+00  8.64E+02
 
 TH10
+       -1.60E+02  0.00E+00  1.06E+07 -9.58E+02  1.92E+03  9.02E+01  1.62E+01  0.00E+00  8.21E+02  1.28E+03
 
 TH11
+       -4.76E+01  0.00E+00  8.24E+05 -2.24E+02  4.59E+02  3.54E+01  3.69E+00  0.00E+00  2.38E+02  3.43E+02  2.01E+05
 
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
 #CPUT: Total CPU Time in Seconds,       24.472
Stop Time:
Sat Sep 25 08:09:54 CDT 2021
