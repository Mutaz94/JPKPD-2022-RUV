Wed Sep 29 17:23:41 CDT 2021
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
$DATA ../../../../data/spa/S2/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1714.20467440892        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9842E+02  6.4995E+01 -2.4685E+01  1.4916E+02  7.1735E+01  5.2552E+01  1.2445E+01  1.6038E+00  1.3096E+01 -4.7784E+00
             8.7687E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1720.07728502039        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0538E+00  8.8517E-01  9.2586E-01  1.0223E+00  8.8646E-01  9.8072E-01  8.3260E-01  9.8782E-01  9.4226E-01  9.6335E-01
             9.3216E-01
 PARAMETER:  1.5245E-01 -2.1981E-02  2.2971E-02  1.2206E-01 -2.0515E-02  8.0529E-02 -8.3200E-02  8.7741E-02  4.0526E-02  6.2657E-02
             2.9754E-02
 GRADIENT:   8.8556E+00 -1.3600E+01 -1.7309E+01 -1.2955E+01  1.8777E+01 -6.4226E+00 -2.1723E-01  5.6181E+00 -1.0033E+01  6.2372E+00
            -1.9897E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1721.13305365793        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0546E+00  7.8710E-01  9.9704E-01  1.1070E+00  8.6847E-01  9.8162E-01  6.2167E-01  9.1006E-01  9.8774E-01  9.8100E-01
             9.8715E-01
 PARAMETER:  1.5320E-01 -1.3940E-01  9.7037E-02  2.0164E-01 -4.1024E-02  8.1446E-02 -3.7534E-01  5.7563E-03  8.7664E-02  8.0814E-02
             8.7064E-02
 GRADIENT:   1.1752E+01  1.3640E+01 -3.0884E+00  3.1023E+01 -5.7667E+00 -5.4731E+00 -1.1218E-01  6.9537E-01  5.1178E+00  4.1477E+00
             3.5094E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1721.66614838185        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0468E+00  6.8452E-01  1.0691E+00  1.1590E+00  8.6280E-01  9.9593E-01  7.2673E-01  9.3389E-01  9.1502E-01  9.5852E-01
             9.8089E-01
 PARAMETER:  1.4576E-01 -2.7904E-01  1.6677E-01  2.4755E-01 -4.7572E-02  9.5926E-02 -2.1920E-01  3.1601E-02  1.1194E-02  5.7632E-02
             8.0702E-02
 GRADIENT:  -1.6449E+00  3.2523E+00  3.0290E+00  2.2647E+00 -3.2331E+00  9.8009E-01  2.0282E-01 -8.6008E-01 -6.2157E-01 -7.0557E-01
             4.7337E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1721.69854409169        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      720
 NPARAMETR:  1.0467E+00  6.0604E-01  1.1139E+00  1.2104E+00  8.5636E-01  9.9228E-01  7.0808E-01  9.7262E-01  8.9231E-01  9.6513E-01
             9.8362E-01
 PARAMETER:  1.4567E-01 -4.0082E-01  2.0788E-01  2.9098E-01 -5.5061E-02  9.2251E-02 -2.4519E-01  7.2241E-02 -1.3942E-02  6.4509E-02
             8.3487E-02
 GRADIENT:   6.9079E-01  2.8502E+00  1.2401E+00  5.6628E+00 -2.3405E+00  1.5404E-02  4.2674E-01 -5.5773E-02  1.3198E+00  2.4349E-01
             1.5482E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1721.70352751299        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  1.0460E+00  5.5184E-01  1.1429E+00  1.2447E+00  8.5141E-01  9.9086E-01  6.6785E-01  9.9307E-01  8.7462E-01  9.6855E-01
             9.8308E-01
 PARAMETER:  1.4497E-01 -4.9450E-01  2.3357E-01  3.1888E-01 -6.0857E-02  9.0819E-02 -3.0369E-01  9.3046E-02 -3.3969E-02  6.8049E-02
             8.2932E-02
 GRADIENT:   1.0265E+00  2.5086E+00  8.9623E-01  5.9772E+00 -2.0966E+00 -2.2222E-01  3.6357E-01  6.1155E-02  1.3420E+00  3.9938E-01
             1.3545E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1721.70566464523        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  1.0481E+00  5.3113E-01  1.1472E+00  1.2452E+00  8.5202E-01  9.9193E-01  4.5376E-01  9.9546E-01  8.6434E-01  9.6485E-01
             9.7728E-01
 PARAMETER:  1.4700E-01 -5.3274E-01  2.3736E-01  3.1928E-01 -6.0148E-02  9.1897E-02 -6.9018E-01  9.5446E-02 -4.5787E-02  6.4212E-02
             7.7022E-02
 GRADIENT:   6.6867E+00 -4.3925E+00 -3.0487E+00 -1.8400E+01  7.9101E+00  3.0153E-01  3.0565E-02 -1.9985E-01 -4.2415E+00 -1.7122E+00
            -1.2902E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1721.79925364242        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  1.0448E+00  5.2915E-01  1.1538E+00  1.2545E+00  8.5040E-01  9.9115E-01  3.4699E-01  1.0014E+00  8.7929E-01  9.7560E-01
             9.7997E-01
 PARAMETER:  1.4383E-01 -5.3649E-01  2.4306E-01  3.2673E-01 -6.2053E-02  9.1112E-02 -9.5846E-01  1.0136E-01 -2.8636E-02  7.5296E-02
             7.9768E-02
 GRADIENT:  -7.0434E-01  7.5520E-01 -3.0508E-01  2.9271E-01 -7.2806E-01 -2.3463E-02  6.5289E-02  3.0048E-03  6.8529E-01  1.8079E-01
             1.8293E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1721.80324513170        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1431
 NPARAMETR:  1.0443E+00  5.1535E-01  1.1669E+00  1.2639E+00  8.5109E-01  9.9092E-01  2.0072E-01  1.0119E+00  8.7838E-01  9.8039E-01
             9.8009E-01
 PARAMETER:  1.4339E-01 -5.6291E-01  2.5433E-01  3.3419E-01 -6.1236E-02  9.0876E-02 -1.5058E+00  1.1180E-01 -2.9681E-02  8.0196E-02
             7.9890E-02
 GRADIENT:  -1.1335E+00  1.1533E+00  2.2448E-01  1.7717E+00 -1.8167E+00 -3.2873E-02  3.1854E-02  1.0297E-02  1.2508E+00  5.0128E-01
             2.8295E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1721.82268810252        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1611
 NPARAMETR:  1.0469E+00  5.0932E-01  1.1711E+00  1.2654E+00  8.5184E-01  9.9115E-01  6.9958E-02  1.0155E+00  8.7224E-01  9.7758E-01
             9.7954E-01
 PARAMETER:  1.4583E-01 -5.7468E-01  2.5795E-01  3.3536E-01 -6.0357E-02  9.1114E-02 -2.5599E+00  1.1535E-01 -3.6693E-02  7.7329E-02
             7.9323E-02
 GRADIENT:   4.7497E+00  1.0044E-01 -2.7009E-01 -3.2185E+00  3.7279E-01  1.0902E-01  5.7152E-03 -1.2774E-01 -8.5347E-01 -2.5740E-01
            -1.3281E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1721.82528861251        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1796             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0467E+00  5.0913E-01  1.1720E+00  1.2650E+00  8.5196E-01  9.9117E-01  1.0266E-02  1.0175E+00  8.7429E-01  9.7834E-01
             9.7955E-01
 PARAMETER:  1.4567E-01 -5.7504E-01  2.5869E-01  3.3507E-01 -6.0215E-02  9.1127E-02 -4.4789E+00  1.1735E-01 -3.4339E-02  7.8100E-02
             7.9335E-02
 GRADIENT:   7.6614E+02  8.5691E+01  6.1328E+00  5.0224E+02  9.0310E+00  4.8155E+01  6.6929E-03  2.5878E-01  1.2587E+01  7.7601E-01
             7.5266E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1721.82569721697        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1978             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0467E+00  5.0977E-01  1.1721E+00  1.2645E+00  8.5234E-01  9.9118E-01  1.0000E-02  1.0182E+00  8.7474E-01  9.7834E-01
             9.7958E-01
 PARAMETER:  1.4568E-01 -5.7379E-01  2.5882E-01  3.3470E-01 -5.9775E-02  9.1137E-02 -4.6127E+00  1.1804E-01 -3.3826E-02  7.8102E-02
             7.9366E-02
 GRADIENT:   7.6620E+02  8.5607E+01  6.0109E+00  5.0108E+02  9.2710E+00  4.8157E+01  0.0000E+00  2.8189E-01  1.2649E+01  7.5237E-01
             7.7444E-01

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1721.82583631608        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     2116
 NPARAMETR:  1.0468E+00  5.1014E-01  1.1718E+00  1.2639E+00  8.5301E-01  9.9119E-01  1.0000E-02  1.0190E+00  8.7521E-01  9.7856E-01
             9.7962E-01
 PARAMETER:  1.4567E-01 -5.7240E-01  2.5926E-01  3.3470E-01 -5.9989E-02  9.1137E-02 -4.6127E+00  1.1842E-01 -3.3796E-02  7.9165E-02
             7.9427E-02
 GRADIENT:  -3.0967E-02  1.1663E-01  2.7009E-01  9.4676E-01 -5.0098E-01 -4.4111E-03  0.0000E+00 -1.6746E-02 -1.2473E-01  4.7279E-02
             4.3871E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2116
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0976E-04 -2.6008E-04 -3.0424E-02 -4.3708E-03 -2.8953E-02
 SE:             2.9873E-02  1.0822E-04  1.6083E-02  2.9242E-02  2.2741E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8639E-01  1.6255E-02  5.8540E-02  8.8118E-01  2.0295E-01

 ETASHRINKSD(%)  1.0000E-10  9.9637E+01  4.6119E+01  2.0369E+00  2.3816E+01
 ETASHRINKVR(%)  1.0000E-10  9.9999E+01  7.0969E+01  4.0323E+00  4.1960E+01
 EBVSHRINKSD(%)  4.0485E-01  9.9673E+01  5.0039E+01  2.3883E+00  2.0986E+01
 EBVSHRINKVR(%)  8.0805E-01  9.9999E+01  7.5039E+01  4.7197E+00  3.7568E+01
 RELATIVEINF(%)  9.8151E+01  7.0705E-05  5.3088E+00  8.0068E+00  9.5659E+00
 EPSSHRINKSD(%)  4.5377E+01
 EPSSHRINKVR(%)  7.0164E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1721.8258363160846     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -986.67500975234645     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.06
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1721.826       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  5.10E-01  1.17E+00  1.26E+00  8.52E-01  9.91E-01  1.00E-02  1.02E+00  8.75E-01  9.79E-01  9.80E-01
 


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
+        1.02E+03
 
 TH 2
+       -2.47E+01  4.88E+02
 
 TH 3
+        5.47E+00  1.07E+02  2.07E+02
 
 TH 4
+       -1.15E+01  5.41E+02 -3.05E+01  8.66E+02
 
 TH 5
+        4.09E+00 -3.27E+02 -3.73E+02 -5.32E+01  1.01E+03
 
 TH 6
+       -2.35E-01 -4.34E+00  1.03E+00 -2.14E+00  3.74E-01  2.00E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -2.16E-01 -1.35E+01 -3.82E+01 -3.21E+00  7.17E+00 -1.20E-01  0.00E+00  2.82E+01
 
 TH 9
+        2.26E+00 -1.07E+02  6.09E+00 -4.14E-01 -2.27E+00  2.47E-01  0.00E+00  1.19E-01  2.37E+02
 
 TH10
+        9.23E-01  6.11E+00 -1.27E+01 -8.68E-01 -8.05E+01  2.70E-01  0.00E+00  2.12E+01  1.47E-01  8.59E+01
 
 TH11
+       -6.47E+00 -1.27E+01 -1.31E+01 -8.65E+00 -7.70E-01  1.45E+00  0.00E+00  9.40E+00  9.54E+00  1.34E+01  2.15E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.771
Stop Time:
Wed Sep 29 17:24:15 CDT 2021
