Wed Sep 29 17:37:33 CDT 2021
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
$DATA ../../../../data/spa/S2/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1735.07751076603        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1155E+02  5.2985E+01 -1.6613E+01  1.2472E+02  2.4544E+01  4.7439E+01  2.6746E+01  6.7266E+00  5.1904E+01  2.0080E+00
             6.3579E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1749.58794568478        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0296E+00  9.6599E-01  9.9586E-01  9.8329E-01  1.0030E+00  1.0020E+00  8.2705E-01  9.7205E-01  7.6814E-01  1.0399E+00
             8.6323E-01
 PARAMETER:  1.2915E-01  6.5400E-02  9.5855E-02  8.3152E-02  1.0301E-01  1.0196E-01 -8.9895E-02  7.1652E-02 -1.6379E-01  1.3916E-01
            -4.7072E-02
 GRADIENT:   3.4578E+01 -2.4058E+01 -2.1513E+01 -1.9968E+01  3.4423E+01 -2.0038E-01  3.0329E+00  4.4222E+00 -5.9117E+00  1.5560E+00
             8.3827E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1751.00162622301        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0188E+00  7.9747E-01  1.1940E+00  1.1143E+00  9.9046E-01  9.7927E-01  6.5567E-01  1.0307E+00  8.0571E-01  1.0497E+00
             8.4629E-01
 PARAMETER:  1.1862E-01 -1.2632E-01  2.7730E-01  2.0824E-01  9.0413E-02  7.9049E-02 -3.2209E-01  1.3027E-01 -1.1603E-01  1.4852E-01
            -6.6891E-02
 GRADIENT:   1.6810E+01  1.4334E+01  1.0653E+00  3.2466E+01  9.0068E-02 -8.1056E+00  7.9118E-01 -5.5633E-01  9.1685E+00 -3.3295E+00
            -4.8187E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1751.52009198873        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      544
 NPARAMETR:  1.0105E+00  7.7022E-01  1.1720E+00  1.1231E+00  9.6758E-01  9.9915E-01  7.0576E-01  9.7978E-01  7.5988E-01  1.0532E+00
             8.4629E-01
 PARAMETER:  1.1041E-01 -1.6108E-01  2.5868E-01  2.1614E-01  6.7045E-02  9.9152E-02 -2.4849E-01  7.9576E-02 -1.7460E-01  1.5186E-01
            -6.6890E-02
 GRADIENT:  -2.0090E+00  9.9254E+00  3.6302E+00  1.2388E+01 -6.4948E+00  1.7884E-01 -2.9716E-01 -8.7409E-01 -7.7509E-01  1.9413E-01
            -7.5473E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1752.19434881212        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0053E+00  4.7923E-01  1.4576E+00  1.3159E+00  9.6847E-01  1.0011E+00  8.0838E-01  1.2200E+00  6.6404E-01  1.0849E+00
             8.5198E-01
 PARAMETER:  1.0533E-01 -6.3558E-01  4.7678E-01  3.7452E-01  6.7965E-02  1.0110E-01 -1.1273E-01  2.9884E-01 -3.0941E-01  1.8151E-01
            -6.0189E-02
 GRADIENT:  -3.7258E+00  9.8254E+00  3.8353E+00  2.0397E+01 -1.0508E+01  3.0273E+00 -2.2436E-01  1.9807E-01 -2.5565E+00  1.3202E+00
             1.9521E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1752.72618384853        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0052E+00  3.0056E-01  1.6858E+00  1.4411E+00  9.8792E-01  9.8315E-01  8.2980E-01  1.4138E+00  6.2387E-01  1.0980E+00
             8.4324E-01
 PARAMETER:  1.0517E-01 -1.1021E+00  6.2225E-01  4.6543E-01  8.7847E-02  8.3005E-02 -8.6566E-02  4.4631E-01 -3.7182E-01  1.9351E-01
            -7.0508E-02
 GRADIENT:   2.2758E+00  9.4546E+00  1.5149E+00  4.3248E+01 -3.5609E+00 -2.7818E+00 -5.0336E-02 -1.2781E+00  1.7957E-01 -1.0943E+00
            -2.4424E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1753.59979101414        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1078
 NPARAMETR:  1.0024E+00  1.3512E-01  1.9966E+00  1.5539E+00  1.0193E+00  9.8165E-01  8.2386E-01  1.7090E+00  5.7749E-01  1.1252E+00
             8.4379E-01
 PARAMETER:  1.0241E-01 -1.9016E+00  7.9142E-01  5.4078E-01  1.1915E-01  8.1481E-02 -9.3757E-02  6.3588E-01 -4.4906E-01  2.1795E-01
            -6.9847E-02
 GRADIENT:   1.3264E+00  4.7514E+00 -9.3917E-01  3.9187E+01 -2.4712E+00 -2.0480E+00 -2.1125E-02  5.9755E-01 -1.7447E+00 -5.7588E-01
            -1.6423E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1754.43560633821        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1259             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0015E+00  3.8113E-02  2.1723E+00  1.5990E+00  1.0364E+00  9.8567E-01  7.5819E-01  1.8383E+00  5.5741E-01  1.1408E+00
             8.4669E-01
 PARAMETER:  1.0146E-01 -3.1672E+00  8.7577E-01  5.6937E-01  1.3571E-01  8.5564E-02 -1.7682E-01  7.0885E-01 -4.8445E-01  2.3171E-01
            -6.6425E-02
 GRADIENT:   6.3259E+02  2.6083E+00  1.1402E+01  1.5510E+03  1.7092E+01  6.2793E+01  1.0202E-02  4.4972E+00  3.9214E+01  2.9493E+00
             1.6290E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1754.52203053618        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1437
 NPARAMETR:  1.0003E+00  2.8706E-02  2.1753E+00  1.6151E+00  1.0292E+00  9.8469E-01  7.4887E-01  1.8327E+00  5.5294E-01  1.1403E+00
             8.4579E-01
 PARAMETER:  1.0032E-01 -3.4506E+00  8.7715E-01  5.7943E-01  1.2874E-01  8.4571E-02 -1.8919E-01  7.0581E-01 -4.9250E-01  2.3129E-01
            -6.7485E-02
 GRADIENT:   1.6508E-01  4.9170E-01  8.3865E-01 -8.7577E+00 -1.0487E+00  1.3302E-02  6.2734E-04 -2.8506E-01  4.6473E-01  3.7851E-01
            -8.4213E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1754.64066260788        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1617             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0009E+00  1.4347E-02  2.1740E+00  1.6167E+00  1.0289E+00  9.8473E-01  7.3653E-01  1.8377E+00  5.5097E-01  1.1367E+00
             8.4583E-01
 PARAMETER:  1.0094E-01 -4.1442E+00  8.7658E-01  5.8036E-01  1.2847E-01  8.4609E-02 -2.0581E-01  7.0854E-01 -4.9608E-01  2.2815E-01
            -6.7438E-02
 GRADIENT:   6.2973E+02  5.3910E-01  1.1760E+01  1.6285E+03  1.4405E+01  6.2608E+01  1.4815E-03  4.3263E+00  4.0029E+01  3.1186E+00
             1.1201E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1754.65733464908        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1795
 NPARAMETR:  1.0002E+00  1.0372E-02  2.1752E+00  1.6257E+00  1.0261E+00  9.8446E-01  7.3200E-01  1.8382E+00  5.4913E-01  1.1360E+00
             8.4570E-01
 PARAMETER:  1.0022E-01 -4.4687E+00  8.7712E-01  5.8593E-01  1.2573E-01  8.4335E-02 -2.1198E-01  7.0878E-01 -4.9942E-01  2.2750E-01
            -6.7595E-02
 GRADIENT:   7.5615E-01  1.7983E-01 -4.6245E-02 -1.3841E+01  6.3598E-01  1.0335E-01  1.0569E-04 -3.4834E-02  8.8417E-01  1.2176E-01
             3.1146E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1754.68882917369        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1973
 NPARAMETR:  1.0006E+00  1.0000E-02  2.1645E+00  1.6219E+00  1.0228E+00  9.8455E-01  7.2941E-01  1.8362E+00  5.4785E-01  1.1329E+00
             8.4518E-01
 PARAMETER:  1.0056E-01 -4.7659E+00  8.7220E-01  5.8360E-01  1.2258E-01  8.4426E-02 -2.1552E-01  7.0771E-01 -5.0175E-01  2.2477E-01
            -6.8205E-02
 GRADIENT:   1.5554E+00  0.0000E+00  2.6146E-01 -2.8790E+01  2.7983E-01  8.1767E-02  9.4826E-05  4.3062E-01  1.1796E-01  1.7940E-01
            -1.2963E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1973
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4169E-05 -3.1169E-04 -4.1595E-02 -1.1330E-02 -4.9711E-02
 SE:             2.9907E-02  1.6530E-04  1.8781E-02  2.8984E-02  2.0661E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9775E-01  5.9342E-02  2.6781E-02  6.9587E-01  1.6126E-02

 ETASHRINKSD(%)  1.0000E-10  9.9446E+01  3.7080E+01  2.8996E+00  3.0784E+01
 ETASHRINKVR(%)  1.0000E-10  9.9997E+01  6.0411E+01  5.7152E+00  5.2092E+01
 EBVSHRINKSD(%)  3.0011E-01  9.9477E+01  4.1386E+01  3.3417E+00  2.5449E+01
 EBVSHRINKVR(%)  5.9932E-01  9.9997E+01  6.5644E+01  6.5717E+00  4.4422E+01
 RELATIVEINF(%)  9.7832E+01  1.1805E-04  1.0739E+01  4.3401E+00  1.3796E+01
 EPSSHRINKSD(%)  4.5935E+01
 EPSSHRINKVR(%)  7.0770E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1754.6888291736948     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1019.5380026099566     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1754.689       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  2.16E+00  1.62E+00  1.02E+00  9.85E-01  7.29E-01  1.84E+00  5.48E-01  1.13E+00  8.45E-01
 


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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.54E+00  0.00E+00  2.67E+01
 
 TH 4
+       -1.51E+01  0.00E+00 -2.47E+01  1.32E+03
 
 TH 5
+       -1.21E+00  0.00E+00 -7.62E+01 -8.08E+01  4.87E+02
 
 TH 6
+        9.27E-01  0.00E+00 -1.97E-02 -2.93E+00  6.39E-01  2.04E+02
 
 TH 7
+        1.54E-01  0.00E+00  7.42E-03  8.40E-03 -1.51E-02 -2.90E-01  1.92E-01
 
 TH 8
+        2.36E-02  0.00E+00 -1.24E+01 -5.61E+00 -9.67E+00  1.11E-01  2.17E-02  1.83E+01
 
 TH 9
+        2.10E+00  0.00E+00  4.11E+00 -3.51E+00  1.39E+00 -2.21E+00 -1.58E-02  8.72E-01  5.87E+02
 
 TH10
+        8.96E-01  0.00E+00 -4.21E-01 -3.99E+00 -6.64E+01  3.85E-01  1.67E-02  6.26E+00  1.17E+00  6.12E+01
 
 TH11
+       -7.83E+00  0.00E+00 -2.63E+00 -1.90E+01 -2.78E+00  2.86E+00  2.39E-02  4.33E+00  2.19E+01  8.18E+00  2.90E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.757
Stop Time:
Wed Sep 29 17:38:05 CDT 2021
