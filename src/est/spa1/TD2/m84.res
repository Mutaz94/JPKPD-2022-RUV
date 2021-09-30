Thu Sep 30 02:24:49 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat84.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2143.61814126889        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6801E+02 -9.8283E+00 -5.3428E+01  6.7965E+01  3.8654E+01  8.3962E+01 -6.6829E+00  1.4805E+01  1.8647E+00  9.4793E+00
            -3.3287E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2155.55158083046        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0236E+00  1.1271E+00  1.2647E+00  9.6242E-01  1.1254E+00  7.6894E-01  1.0607E+00  8.9535E-01  1.0533E+00  9.6348E-01
             1.0540E+00
 PARAMETER:  1.2330E-01  2.1965E-01  3.3481E-01  6.1692E-02  2.1816E-01 -1.6274E-01  1.5888E-01 -1.0539E-02  1.5198E-01  6.2797E-02
             1.5256E-01
 GRADIENT:  -4.4232E+01  1.6470E+01  5.4249E+00  1.2608E+01  5.3458E-01 -4.1812E+01  3.0283E+00  1.6093E-01 -3.0441E+00 -2.1761E+01
             3.7886E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2157.07915360692        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0244E+00  1.2780E+00  1.4702E+00  8.7120E-01  1.2873E+00  8.1898E-01  8.8821E-01  7.9652E-01  1.1956E+00  1.1957E+00
             1.0704E+00
 PARAMETER:  1.2415E-01  3.4526E-01  4.8539E-01 -3.7881E-02  3.5253E-01 -9.9700E-02 -1.8542E-02 -1.2751E-01  2.7864E-01  2.7876E-01
             1.6799E-01
 GRADIENT:  -3.7136E+01  1.9697E+01  1.8556E+01  1.1398E+01  1.2367E+01 -1.3331E+01  2.0237E+00 -8.6805E+00 -1.7313E+00 -8.1648E+00
             1.6685E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2160.15406347056        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  1.0375E+00  1.3665E+00  1.2188E+00  8.0080E-01  1.2378E+00  8.4728E-01  8.2352E-01  9.0110E-01  1.2701E+00  1.1799E+00
             1.0416E+00
 PARAMETER:  1.3677E-01  4.1222E-01  2.9786E-01 -1.2214E-01  3.1334E-01 -6.5723E-02 -9.4173E-02 -4.1426E-03  3.3907E-01  2.6540E-01
             1.4076E-01
 GRADIENT:   2.4751E+00  5.8098E+00  8.2893E-01  6.2062E+00 -2.2167E+00  1.0812E-01 -1.2561E+00 -8.9822E-02 -3.2301E-01  4.7890E-01
             5.0896E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2160.21530056878        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0385E+00  1.4846E+00  1.0718E+00  7.2686E-01  1.2439E+00  8.4735E-01  8.2248E-01  7.7731E-01  1.3438E+00  1.1692E+00
             1.0421E+00
 PARAMETER:  1.3778E-01  4.9514E-01  1.6935E-01 -2.1902E-01  3.1823E-01 -6.5642E-02 -9.5426E-02 -1.5192E-01  3.9550E-01  2.5631E-01
             1.4128E-01
 GRADIENT:   3.0202E+00  1.3719E+01  1.3067E+00  1.0447E+01 -4.5276E+00 -3.4798E-01 -9.8933E-02  6.6492E-02  1.0807E-01  8.3268E-01
             4.2578E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2160.33982105934        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  1.0387E+00  1.6400E+00  8.9022E-01  6.2152E-01  1.2646E+00  8.4812E-01  8.0431E-01  5.8430E-01  1.4697E+00  1.1594E+00
             1.0419E+00
 PARAMETER:  1.3795E-01  5.9470E-01 -1.6289E-02 -3.7559E-01  3.3472E-01 -6.4728E-02 -1.1777E-01 -4.3734E-01  4.8508E-01  2.4789E-01
             1.4103E-01
 GRADIENT:   6.4903E-01  1.2311E+01  2.7542E-02  8.4847E+00 -2.7736E+00 -5.3464E-01  9.4991E-01  3.7983E-01  2.9020E-01  3.1126E-01
             2.7651E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2160.47223702775        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1078
 NPARAMETR:  1.0382E+00  1.6726E+00  8.4625E-01  5.9081E-01  1.2737E+00  8.4947E-01  7.9171E-01  4.4678E-01  1.5117E+00  1.1608E+00
             1.0416E+00
 PARAMETER:  1.3744E-01  6.1440E-01 -6.6946E-02 -4.2626E-01  3.4193E-01 -6.3143E-02 -1.3356E-01 -7.0568E-01  5.1326E-01  2.4907E-01
             1.4074E-01
 GRADIENT:  -1.3865E+00 -8.0291E-01  6.6849E-01  1.7998E+00  4.7330E-01  4.5189E-03  1.0313E-01  5.2323E-02 -1.1909E-01  1.8628E-01
             1.0787E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2160.48829462419        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.0378E+00  1.7124E+00  7.8153E-01  5.6786E-01  1.2684E+00  8.4993E-01  7.8682E-01  2.7266E-01  1.5530E+00  1.1485E+00
             1.0415E+00
 PARAMETER:  1.3714E-01  6.3791E-01 -1.4651E-01 -4.6588E-01  3.3776E-01 -6.2606E-02 -1.3975E-01 -1.1995E+00  5.4020E-01  2.3845E-01
             1.4063E-01
 GRADIENT:  -3.3027E+00  3.8791E+00 -1.2652E+00  6.2865E+00  5.2319E-02  8.8544E-04  2.9358E-02  1.1322E-01  1.2273E+00  3.3852E-01
             5.3163E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2160.55124955680        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1435
 NPARAMETR:  1.0391E+00  1.7204E+00  7.7069E-01  5.5858E-01  1.2699E+00  8.4997E-01  7.8582E-01  1.7133E-01  1.5587E+00  1.1472E+00
             1.0415E+00
 PARAMETER:  1.3833E-01  6.4255E-01 -1.6046E-01 -4.8236E-01  3.3893E-01 -6.2556E-02 -1.4103E-01 -1.6641E+00  5.4388E-01  2.3728E-01
             1.4066E-01
 GRADIENT:   3.0841E-01 -7.8175E-01  6.8146E-02  2.5855E+00 -4.6170E-01  2.1976E-02 -7.8666E-02  3.2984E-02  3.5126E-01  1.2176E-02
            -1.2396E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2160.57493153996        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1617
 NPARAMETR:  1.0405E+00  1.7309E+00  7.5676E-01  5.4674E-01  1.2748E+00  8.5008E-01  7.8420E-01  4.0145E-02  1.5718E+00  1.1494E+00
             1.0413E+00
 PARAMETER:  1.3970E-01  6.4864E-01 -1.7870E-01 -5.0378E-01  3.4278E-01 -6.2425E-02 -1.4309E-01 -3.1152E+00  5.5224E-01  2.3921E-01
             1.4049E-01
 GRADIENT:   4.4527E+00 -9.0170E+00 -3.5097E-01 -9.4796E-01  1.0718E+00  6.0782E-02 -2.4711E-02  2.4910E-03  4.7089E-02  9.2503E-02
             6.7296E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2160.57793612538        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1797
 NPARAMETR:  1.0401E+00  1.7318E+00  7.5699E-01  5.4724E-01  1.2737E+00  8.5007E-01  7.8367E-01  1.0000E-02  1.5728E+00  1.1482E+00
             1.0413E+00
 PARAMETER:  1.3927E-01  6.4915E-01 -1.7840E-01 -5.0288E-01  3.4190E-01 -6.2442E-02 -1.4376E-01 -5.3288E+00  5.5285E-01  2.3818E-01
             1.4048E-01
 GRADIENT:   3.0985E+00 -6.7236E+00  3.9547E-02 -3.0624E-01  1.0356E-01  4.9628E-02 -7.1518E-02  0.0000E+00  7.6640E-02 -1.3402E-02
            -4.4309E-02

0ITERATION NO.:   51    OBJECTIVE VALUE:  -2160.57793612538        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1819
 NPARAMETR:  1.0401E+00  1.7318E+00  7.5699E-01  5.4724E-01  1.2737E+00  8.5007E-01  7.8367E-01  1.0000E-02  1.5728E+00  1.1482E+00
             1.0413E+00
 PARAMETER:  1.3927E-01  6.4915E-01 -1.7840E-01 -5.0288E-01  3.4190E-01 -6.2442E-02 -1.4376E-01 -5.3288E+00  5.5285E-01  2.3818E-01
             1.4048E-01
 GRADIENT:   3.0985E+00 -6.7236E+00  3.9547E-02 -3.0624E-01  1.0356E-01  4.9628E-02 -7.1518E-02  0.0000E+00  7.6640E-02 -1.3402E-02
            -4.4309E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1819
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9163E-04 -3.4996E-02 -2.5063E-04  2.8362E-02 -3.8171E-02
 SE:             2.9850E-02  2.1752E-02  1.0318E-04  2.3338E-02  2.3332E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8953E-01  1.0764E-01  1.5143E-02  2.2426E-01  1.0184E-01

 ETASHRINKSD(%)  1.0000E-10  2.7129E+01  9.9654E+01  2.1815E+01  2.1835E+01
 ETASHRINKVR(%)  1.0000E-10  4.6898E+01  9.9999E+01  3.8871E+01  3.8902E+01
 EBVSHRINKSD(%)  4.7183E-01  2.5151E+01  9.9712E+01  2.4708E+01  1.8941E+01
 EBVSHRINKVR(%)  9.4143E-01  4.3976E+01  9.9999E+01  4.3311E+01  3.4294E+01
 RELATIVEINF(%)  9.8932E+01  5.0084E+00  2.1266E-04  5.3533E+00  2.0646E+01
 EPSSHRINKSD(%)  3.2851E+01
 EPSSHRINKVR(%)  5.4911E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2160.5779361253844     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1241.6394029207117     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.08
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2160.578       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.73E+00  7.57E-01  5.47E-01  1.27E+00  8.50E-01  7.84E-01  1.00E-02  1.57E+00  1.15E+00  1.04E+00
 


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
+        1.41E+03
 
 TH 2
+       -6.88E+00  3.61E+02
 
 TH 3
+        5.70E+00  8.98E+01  1.93E+02
 
 TH 4
+       -7.36E+00  3.58E+02 -1.43E+02  8.66E+02
 
 TH 5
+        1.50E+00 -1.18E+02 -1.43E+02  1.51E+02  3.41E+02
 
 TH 6
+        6.87E-01 -9.27E-01  1.84E+00 -2.08E+00  3.28E-01  2.71E+02
 
 TH 7
+        2.05E+00  1.85E-01  7.65E+00 -1.82E+01 -1.98E+01 -5.18E-01  1.15E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.83E-02 -1.69E+01 -2.13E+01  4.99E+01  4.78E+00 -3.88E-01  2.33E+01  0.00E+00  3.33E+01
 
 TH10
+        1.59E+00 -1.36E+01 -2.26E+01  7.00E+00 -4.50E+01  8.11E-01  8.29E-01  0.00E+00  5.09E+00  7.03E+01
 
 TH11
+       -1.05E+01 -1.83E+01 -2.18E+01 -2.47E+00  2.47E+00  2.41E+00  1.01E+01  0.00E+00  1.86E+00  1.74E+01  3.83E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.635
Stop Time:
Thu Sep 30 02:25:27 CDT 2021
