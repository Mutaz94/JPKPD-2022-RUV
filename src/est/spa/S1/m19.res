Sat Sep 25 09:45:05 CDT 2021
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
$DATA ../../../../data/spa/S1/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1732.19360238882        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.5268E+01 -4.5145E+01 -2.2444E+01 -1.0010E+01  4.8002E+01  1.2239E+01  7.2188E+00  4.1384E+00  4.4221E+01 -5.9058E+00
             3.1406E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1741.29780231884        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0269E+00  1.0545E+00  9.8131E-01  9.8376E-01  9.8313E-01  9.6034E-01  9.6524E-01  9.8518E-01  8.0515E-01  1.0210E+00
             9.1143E-01
 PARAMETER:  1.2655E-01  1.5307E-01  8.1130E-02  8.3629E-02  8.2987E-02  5.9536E-02  6.4618E-02  8.5071E-02 -1.1673E-01  1.2079E-01
             7.2599E-03
 GRADIENT:   6.6316E+01  8.8089E-01  2.9216E+00  4.2014E+00  5.2291E+00 -1.8148E+00 -3.1135E+00  1.7285E+00  6.0189E+00 -2.0745E+00
            -7.7235E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1742.65270594751        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0208E+00  9.7290E-01  8.2403E-01  1.0279E+00  8.7112E-01  9.6668E-01  1.1511E+00  7.2348E-01  7.0609E-01  8.9161E-01
             9.2894E-01
 PARAMETER:  1.2060E-01  7.2523E-02 -9.3551E-02  1.2756E-01 -3.7978E-02  6.6109E-02  2.4070E-01 -2.2368E-01 -2.4801E-01 -1.4728E-02
             2.6289E-02
 GRADIENT:   4.2138E+01  6.1455E+00 -9.8143E+00  2.7829E+01  2.1612E+01  2.8780E-01  1.9290E+00  1.8655E+00 -1.0583E-01 -4.3017E+00
             1.2005E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1742.88599376472        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0124E+00  9.9213E-01  7.2152E-01  1.0020E+00  8.1951E-01  9.6825E-01  1.1111E+00  5.2302E-01  7.2420E-01  8.8149E-01
             9.2034E-01
 PARAMETER:  1.1234E-01  9.2100E-02 -2.2640E-01  1.0201E-01 -9.9050E-02  6.7738E-02  2.0533E-01 -5.4813E-01 -2.2268E-01 -2.6147E-02
             1.6990E-02
 GRADIENT:   1.6387E+01 -3.7074E-01 -9.7851E+00  1.3883E+01  1.2028E+01  4.0319E-02  7.8228E-01  1.8194E+00  8.7694E-01  2.0752E+00
            -6.8372E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1742.90969319430        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0100E+00  1.0050E+00  6.7485E-01  9.8727E-01  7.9402E-01  9.6939E-01  1.0954E+00  4.1801E-01  7.3008E-01  8.5104E-01
             9.2126E-01
 PARAMETER:  1.0991E-01  1.0495E-01 -2.9327E-01  8.7190E-02 -1.3064E-01  6.8914E-02  1.9112E-01 -7.7225E-01 -2.1460E-01 -6.1292E-02
             1.7986E-02
 GRADIENT:   7.8780E+00 -1.3956E-01 -5.4008E+00  7.5484E+00  5.8131E+00  5.3337E-03  4.4176E-01  1.1837E+00  5.5284E-01  1.4896E+00
            -1.9705E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1742.91334093240        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0090E+00  1.0102E+00  6.4834E-01  9.7970E-01  7.7903E-01  9.7000E-01  1.0886E+00  3.4056E-01  7.3296E-01  8.3334E-01
             9.2151E-01
 PARAMETER:  1.0894E-01  1.1015E-01 -3.3334E-01  7.9487E-02 -1.4970E-01  6.9543E-02  1.8492E-01 -9.7717E-01 -2.1066E-01 -8.2312E-02
             1.8263E-02
 GRADIENT:   4.3090E+00 -2.7908E-01 -3.4158E+00  4.3089E+00  3.3690E+00 -1.0409E-02  2.7560E-01  7.9898E-01  3.7415E-01  1.0544E+00
            -3.3757E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1742.91487672065        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      438
 NPARAMETR:  1.0085E+00  1.0128E+00  6.3269E-01  9.7542E-01  7.6997E-01  9.7037E-01  1.0851E+00  2.8273E-01  7.3462E-01  8.2284E-01
             9.2159E-01
 PARAMETER:  1.0845E-01  1.1272E-01 -3.5778E-01  7.5115E-02 -1.6141E-01  6.9920E-02  1.8170E-01 -1.1633E+00 -2.0840E-01 -9.4992E-02
             1.8345E-02
 GRADIENT:   2.4731E+00 -3.3047E-01 -2.2664E+00  2.5485E+00  2.0567E+00 -1.3894E-02  1.8241E-01  5.5620E-01  2.6530E-01  7.6186E-01
             2.9331E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1742.91593288571        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  1.0082E+00  1.0147E+00  6.2098E-01  9.7221E-01  7.6317E-01  9.7064E-01  1.0826E+00  2.2910E-01  7.3591E-01  8.1500E-01
             9.2161E-01
 PARAMETER:  1.0813E-01  1.1460E-01 -3.7646E-01  7.1820E-02 -1.7027E-01  7.0200E-02  1.7933E-01 -1.3736E+00 -2.0664E-01 -1.0457E-01
             1.8370E-02
 GRADIENT:   1.2288E+00 -3.4609E-01 -1.4240E+00  1.3179E+00  1.1534E+00 -1.5149E-02  1.1413E-01  3.6861E-01  1.8072E-01  5.2702E-01
             6.1118E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1742.91850518514        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      578
 NPARAMETR:  1.0080E+00  1.0160E+00  6.1278E-01  9.6996E-01  7.5842E-01  9.7083E-01  1.0808E+00  1.8117E-01  7.3686E-01  8.0953E-01
             9.2160E-01
 PARAMETER:  1.0794E-01  1.1592E-01 -3.8975E-01  6.9495E-02 -1.7651E-01  7.0393E-02  1.7768E-01 -1.6083E+00 -2.0536E-01 -1.1130E-01
             1.8358E-02
 GRADIENT:   4.6685E-01 -3.2901E-01 -8.5074E-01  5.4463E-01  5.8576E-01 -1.4476E-02  6.8149E-02  2.3256E-01  1.1945E-01  3.4856E-01
             7.0903E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1742.92556552558        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  1.0078E+00  1.0172E+00  6.0610E-01  9.6804E-01  7.5460E-01  9.7099E-01  1.0792E+00  1.2458E-01  7.3768E-01  8.0505E-01
             9.2156E-01
 PARAMETER:  1.0781E-01  1.1709E-01 -4.0070E-01  6.7513E-02 -1.8157E-01  7.0557E-02  1.7626E-01 -1.9828E+00 -2.0424E-01 -1.1685E-01
             1.8318E-02
 GRADIENT:  -5.9748E-02 -2.5945E-01 -2.5981E-01 -1.2410E-01  7.3716E-02 -7.7026E-03  1.4469E-02  1.1055E-01  4.1231E-02  1.3510E-01
             5.3772E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1743.70242344186        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  1.0277E+00  1.0272E+00  6.3424E-01  9.6804E-01  7.7588E-01  9.8478E-01  1.0710E+00  8.0028E-02  7.5060E-01  8.5308E-01
             9.2821E-01
 PARAMETER:  1.2734E-01  1.2683E-01 -3.5532E-01  6.7515E-02 -1.5376E-01  8.4663E-02  1.6856E-01 -2.4254E+00 -1.8689E-01 -5.8903E-02
             2.5508E-02
 GRADIENT:  -9.2935E+00 -2.1213E+00  2.7681E+00 -1.0279E+01 -9.3936E+00  1.4519E+00 -4.4143E-02  3.2199E-02  8.5428E-01  2.6322E+00
             3.1384E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1743.80447366049        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      962
 NPARAMETR:  1.0315E+00  1.0104E+00  6.4874E-01  9.8292E-01  7.8042E-01  9.8059E-01  1.0925E+00  3.4432E-02  7.4171E-01  8.5440E-01
             9.2244E-01
 PARAMETER:  1.3099E-01  1.1037E-01 -3.3272E-01  8.2773E-02 -1.4792E-01  8.0400E-02  1.8847E-01 -3.2688E+00 -1.9880E-01 -5.7351E-02
             1.9266E-02
 GRADIENT:  -5.9689E-01  6.7677E-02 -2.0748E-01  4.4782E-01  1.2678E-01 -8.3956E-02 -4.0342E-02  3.1611E-03  1.8488E-02  1.1204E-01
            -2.1893E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1743.80611856700        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1137
 NPARAMETR:  1.0318E+00  1.0123E+00  6.4773E-01  9.8155E-01  7.8061E-01  9.8082E-01  1.0911E+00  1.0000E-02  7.4237E-01  8.5348E-01
             9.2255E-01
 PARAMETER:  1.3128E-01  1.1225E-01 -3.3428E-01  8.1376E-02 -1.4768E-01  8.0636E-02  1.8720E-01 -4.7750E+00 -1.9790E-01 -5.8437E-02
             1.9388E-02
 GRADIENT:   3.8452E-02  1.3190E-02  2.1691E-02  2.6500E-03 -4.6396E-02  1.2247E-03  2.6838E-04  0.0000E+00  5.1560E-03 -1.2530E-03
             3.4919E-03

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1743.80612670051        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  1.0318E+00  1.0133E+00  6.4756E-01  9.8095E-01  7.8099E-01  9.8083E-01  1.0901E+00  1.0000E-02  7.4271E-01  8.5365E-01
             9.2255E-01
 PARAMETER:  1.3126E-01  1.1326E-01 -3.3455E-01  8.0762E-02 -1.4720E-01  8.0642E-02  1.8631E-01 -4.8527E+00 -1.9744E-01 -5.8231E-02
             1.9383E-02
 GRADIENT:  -4.3149E-03  1.0893E-03 -9.9706E-04  3.1900E-03  8.5339E-04  3.0111E-04 -1.4729E-04  0.0000E+00 -6.7495E-04  8.3717E-04
            -4.7568E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1229
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.4058E-05 -2.7029E-03 -4.5617E-04 -1.6422E-03 -1.0701E-02
 SE:             2.9853E-02  2.2257E-02  1.9815E-04  2.3606E-02  2.3683E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9909E-01  9.0334E-01  2.1327E-02  9.4454E-01  6.5138E-01

 ETASHRINKSD(%)  1.0000E-10  2.5436E+01  9.9336E+01  2.0917E+01  2.0659E+01
 ETASHRINKVR(%)  1.0000E-10  4.4401E+01  9.9996E+01  3.7458E+01  3.7050E+01
 EBVSHRINKSD(%)  3.8229E-01  2.5334E+01  9.9415E+01  2.1467E+01  1.9643E+01
 EBVSHRINKVR(%)  7.6312E-01  4.4250E+01  9.9997E+01  3.8326E+01  3.5427E+01
 RELATIVEINF(%)  9.8991E+01  3.1879E+00  3.3352E-04  3.7333E+00  5.3423E+00
 EPSSHRINKSD(%)  4.4555E+01
 EPSSHRINKVR(%)  6.9259E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1743.8061267005107     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1008.6553001367726     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1743.806       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.01E+00  6.48E-01  9.81E-01  7.81E-01  9.81E-01  1.09E+00  1.00E-02  7.43E-01  8.54E-01  9.23E-01
 


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
+        1.08E+03
 
 TH 2
+       -6.65E+00  5.26E+02
 
 TH 3
+        1.88E+01  2.64E+02  1.01E+03
 
 TH 4
+       -1.67E+01  4.69E+02 -5.24E+02  1.30E+03
 
 TH 5
+       -3.04E+00 -4.45E+02 -1.11E+03  5.12E+02  1.61E+03
 
 TH 6
+        5.40E-01 -1.95E+00  5.75E+00 -3.72E+00 -1.80E+00  2.03E+02
 
 TH 7
+        9.35E-01  3.52E+01 -2.17E+01 -1.12E+01  2.18E+00 -3.52E-01  5.93E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.53E+00 -2.64E+01 -3.90E+01  3.21E+01  7.14E+00 -4.38E-01  2.40E+01  0.00E+00  1.55E+02
 
 TH10
+       -1.45E+00 -1.54E+01 -8.87E+01 -2.06E+01 -5.63E+01  4.15E-01  1.58E+01  0.00E+00  1.63E+01  1.13E+02
 
 TH11
+       -7.72E+00 -1.37E+01 -3.71E+01 -8.99E-01  1.00E+01  2.20E+00  6.19E+00  0.00E+00  1.35E+01  2.38E+01  2.48E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.551
Stop Time:
Sat Sep 25 09:45:23 CDT 2021
