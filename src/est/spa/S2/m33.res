Sat Sep 18 13:23:09 CDT 2021
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
$DATA ../../../../data/spa/S2/dat33.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1610.34103698614        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6473E+02 -1.1814E+01  7.3135E+00  1.0507E+00 -5.5693E+00  1.5423E+01 -8.5684E+00 -3.0506E+00  2.6728E+01 -1.0286E+01
            -3.9447E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1620.19606988352        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.5848E-01  1.0957E+00  9.7566E-01  9.4127E-01  1.0493E+00  9.2401E-01  1.1192E+00  1.0123E+00  8.0385E-01  1.0629E+00
             1.1176E+00
 PARAMETER:  5.7589E-02  1.9136E-01  7.5360E-02  3.9479E-02  1.4813E-01  2.0971E-02  2.1259E-01  1.1223E-01 -1.1835E-01  1.6096E-01
             2.1119E-01
 GRADIENT:   6.5539E+01 -4.1613E-01  3.2136E+00 -2.0870E+00 -1.3964E+00 -8.0283E+00 -5.8495E+00 -3.2005E-01 -8.9292E-01  1.1344E+00
             9.7654E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1620.64491274346        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.5348E-01  9.7080E-01  9.5815E-01  1.0264E+00  9.7492E-01  9.3577E-01  1.2809E+00  8.9546E-01  7.3384E-01  1.0147E+00
             1.1176E+00
 PARAMETER:  5.2360E-02  7.0370E-02  5.7246E-02  1.2601E-01  7.4596E-02  3.3617E-02  3.4753E-01 -1.0421E-02 -2.0946E-01  1.1454E-01
             2.1118E-01
 GRADIENT:   5.1952E+01  9.4675E+00  7.3580E-01  2.0117E+01 -3.7029E+00 -2.4887E+00 -2.3943E+00  6.4757E-01 -3.2225E+00  1.1036E+00
             1.0818E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1621.02045755691        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.3964E-01  1.0119E+00  8.4349E-01  9.8632E-01  9.3932E-01  9.3823E-01  1.2619E+00  7.0257E-01  7.5357E-01  9.8358E-01
             1.0961E+00
 PARAMETER:  3.7737E-02  1.1183E-01 -7.0204E-02  8.6227E-02  3.7400E-02  3.6237E-02  3.3264E-01 -2.5301E-01 -1.8293E-01  8.3442E-02
             1.9179E-01
 GRADIENT:   1.4150E+01  2.2648E+00 -1.4929E+00  6.0204E+00  1.5154E+00 -1.7281E+00  3.8905E-01  5.4464E-01 -1.1513E+00  6.0842E-01
             4.0235E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1621.02108586825        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.3841E-01  1.0189E+00  8.1905E-01  9.7946E-01  9.2883E-01  9.3947E-01  1.2544E+00  6.5497E-01  7.5767E-01  9.7483E-01
             1.0933E+00
 PARAMETER:  3.6428E-02  1.1869E-01 -9.9614E-02  7.9250E-02  2.6170E-02  3.7561E-02  3.2667E-01 -3.2316E-01 -1.7751E-01  7.4512E-02
             1.8917E-01
 GRADIENT:   1.0499E+01  1.5519E+00 -1.4790E+00  4.6113E+00  1.3904E+00 -1.3492E+00  3.6673E-01  5.2546E-01 -8.8005E-01  5.6585E-01
             3.1211E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1621.02135384475        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.3783E-01  1.0237E+00  8.0011E-01  9.7461E-01  9.2017E-01  9.4023E-01  1.2493E+00  6.1469E-01  7.6024E-01  9.6765E-01
             1.0916E+00
 PARAMETER:  3.5814E-02  1.2344E-01 -1.2301E-01  7.4286E-02  1.6803E-02  3.8371E-02  3.2258E-01 -3.8663E-01 -1.7412E-01  6.7115E-02
             1.8769E-01
 GRADIENT:   8.6150E+00  1.2076E+00 -1.4035E+00  3.8540E+00  1.2625E+00 -1.1369E+00  3.4062E-01  4.9537E-01 -7.3036E-01  5.2694E-01
             2.6282E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1621.20379916320        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      555
 NPARAMETR:  9.4798E-01  1.0725E+00  7.4387E-01  9.4201E-01  9.0974E-01  9.5250E-01  1.2093E+00  4.9848E-01  7.8646E-01  9.5467E-01
             1.0833E+00
 PARAMETER:  4.6581E-02  1.6996E-01 -1.9589E-01  4.0263E-02  5.4084E-03  5.1338E-02  2.9005E-01 -5.9619E-01 -1.4021E-01  5.3606E-02
             1.8001E-01
 GRADIENT:   1.3808E+00 -4.8589E-01 -4.5715E-01 -8.5444E-01 -8.1569E-01  7.6069E-01  2.2700E-01  3.0557E-01  4.8199E-01  2.6841E-01
            -6.4971E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1621.24860021460        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  9.4816E-01  1.1254E+00  6.7726E-01  9.0463E-01  8.9766E-01  9.5160E-01  1.1633E+00  2.5762E-01  8.0300E-01  9.4199E-01
             1.0853E+00
 PARAMETER:  4.6764E-02  2.1816E-01 -2.8970E-01 -2.2767E-04 -7.9628E-03  5.0395E-02  2.5125E-01 -1.2563E+00 -1.1940E-01  4.0241E-02
             1.8188E-01
 GRADIENT:   6.3719E-02 -8.6158E-01 -1.1313E+00  2.8269E-02  1.5305E+00  6.6613E-02  1.4343E-01  5.6475E-02  1.2428E-01  1.6666E-01
             6.1073E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1621.26574460406        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  9.4804E-01  1.1143E+00  6.6169E-01  9.0960E-01  8.8067E-01  9.5167E-01  1.1727E+00  7.9503E-02  7.9771E-01  9.3124E-01
             1.0851E+00
 PARAMETER:  4.6641E-02  2.0820E-01 -3.1296E-01  5.2505E-03 -2.7074E-02  5.0463E-02  2.5932E-01 -2.4320E+00 -1.2600E-01  2.8764E-02
             1.8164E-01
 GRADIENT:  -5.5477E-01 -2.6269E-01 -3.8163E-01 -1.6826E-01  2.1908E-01  3.7654E-02  1.2734E-01  3.8543E-03 -8.6860E-02  1.9837E-01
             6.0673E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1621.26699232730        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1086
 NPARAMETR:  9.4834E-01  1.1150E+00  6.6060E-01  9.0911E-01  8.8024E-01  9.5159E-01  1.1715E+00  4.4489E-02  7.9888E-01  9.3012E-01
             1.0852E+00
 PARAMETER:  4.6958E-02  2.0887E-01 -3.1461E-01  4.7160E-03 -2.7555E-02  5.0376E-02  2.5825E-01 -3.0125E+00 -1.2454E-01  2.7559E-02
             1.8172E-01
 GRADIENT:   1.9086E-01 -1.5630E-01 -1.7943E-01 -6.7875E-02  1.2983E-01 -2.6205E-03  4.8537E-02  9.8985E-04  2.1561E-02  4.4727E-02
             5.0617E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1621.26744033180        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  9.4826E-01  1.1145E+00  6.6084E-01  9.0955E-01  8.8003E-01  9.5159E-01  1.1718E+00  1.0000E-02  7.9870E-01  9.3021E-01
             1.0851E+00
 PARAMETER:  4.6877E-02  2.0836E-01 -3.1424E-01  5.1915E-03 -2.7800E-02  5.0377E-02  2.5851E-01 -4.7358E+00 -1.2477E-01  2.7650E-02
             1.8169E-01
 GRADIENT:   2.6519E-03  4.7966E-03 -8.1777E-03  6.5914E-03 -4.0459E-03 -1.9552E-03  2.0172E-03  0.0000E+00 -8.6345E-04  3.6862E-03
             5.9679E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1621.26744033180        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1283
 NPARAMETR:  9.4826E-01  1.1145E+00  6.6084E-01  9.0955E-01  8.8003E-01  9.5159E-01  1.1718E+00  1.0000E-02  7.9870E-01  9.3021E-01
             1.0851E+00
 PARAMETER:  4.6877E-02  2.0836E-01 -3.1424E-01  5.1915E-03 -2.7800E-02  5.0377E-02  2.5851E-01 -4.7358E+00 -1.2477E-01  2.7650E-02
             1.8169E-01
 GRADIENT:   2.6519E-03  4.7966E-03 -8.1777E-03  6.5914E-03 -4.0459E-03 -1.9552E-03  2.0172E-03  0.0000E+00 -8.6345E-04  3.6862E-03
             5.9679E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1283
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2161E-04 -3.4017E-03 -4.1155E-04 -2.6585E-03 -1.5621E-02
 SE:             2.9799E-02  2.3558E-02  1.6843E-04  2.2015E-02  2.2743E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9674E-01  8.8519E-01  1.4546E-02  9.0388E-01  4.9217E-01

 ETASHRINKSD(%)  1.7061E-01  2.1078E+01  9.9436E+01  2.6247E+01  2.3808E+01
 ETASHRINKVR(%)  3.4093E-01  3.7713E+01  9.9997E+01  4.5606E+01  4.1948E+01
 EBVSHRINKSD(%)  5.3330E-01  2.0705E+01  9.9504E+01  2.7037E+01  2.2565E+01
 EBVSHRINKVR(%)  1.0637E+00  3.7123E+01  9.9998E+01  4.6764E+01  4.0038E+01
 RELATIVEINF(%)  9.8678E+01  3.8588E+00  2.4500E-04  3.1010E+00  6.0058E+00
 EPSSHRINKSD(%)  4.3004E+01
 EPSSHRINKVR(%)  6.7514E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1621.2674403317976     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -886.11661376805944     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1621.267       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.48E-01  1.11E+00  6.61E-01  9.10E-01  8.80E-01  9.52E-01  1.17E+00  1.00E-02  7.99E-01  9.30E-01  1.09E+00
 


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
+        1.35E+03
 
 TH 2
+       -7.79E+00  4.22E+02
 
 TH 3
+        1.99E+01  1.82E+02  7.05E+02
 
 TH 4
+       -2.27E+01  3.92E+02 -4.36E+02  1.13E+03
 
 TH 5
+       -4.56E+00 -2.77E+02 -7.13E+02  4.33E+02  9.92E+02
 
 TH 6
+       -2.86E+00 -2.74E+00  3.76E+00 -5.64E+00 -2.75E+00  2.15E+02
 
 TH 7
+        6.39E-01  2.75E+01 -1.96E+01 -1.68E+01  3.60E+00 -1.32E+00  6.21E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.92E+00 -1.85E+01 -3.69E+01  3.00E+01  9.76E+00  2.06E+00  1.98E+01  0.00E+00  1.07E+02
 
 TH10
+       -1.80E+00 -1.20E+01 -5.35E+01 -1.40E+01 -5.72E+01  2.22E-01  1.03E+01  0.00E+00  9.93E+00  8.61E+01
 
 TH11
+       -8.84E+00 -1.55E+01 -3.80E+01 -1.09E+00  1.89E+00  4.51E+00  6.68E+00  0.00E+00  1.35E+01  1.93E+01  1.84E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.386
Stop Time:
Sat Sep 18 13:23:30 CDT 2021
