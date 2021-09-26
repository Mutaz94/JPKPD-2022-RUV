Sat Sep 25 02:17:48 CDT 2021
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
$DATA ../../../../data/int/SL3/dat45.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      983
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      883
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1116.98309075205        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0499E+02  6.0278E+01  1.1941E+02  8.8226E+00  1.9450E+02 -6.6139E+01 -1.1282E+02 -6.7400E+02 -2.4800E+02 -5.0829E+01
            -4.2219E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2818.97037313632        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6315E-01  1.1925E+00  1.2603E+00  8.7040E-01  1.1712E+00  1.2732E+00  8.6652E-01  1.0766E+00  1.1328E+00  8.6952E-01
             2.0130E+00
 PARAMETER:  6.2455E-02  2.7602E-01  3.3137E-01 -3.8806E-02  2.5802E-01  3.4150E-01 -4.3271E-02  1.7383E-01  2.2467E-01 -3.9817E-02
             7.9962E-01
 GRADIENT:  -1.3159E+01 -5.9570E+01 -1.4259E+01 -2.8412E+01  6.5167E+01  3.9188E+01 -5.2242E+00 -1.0708E+01 -1.6430E+01 -4.0727E+01
            -3.4598E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2835.38106634086        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0063E+00  1.3601E+00  2.3863E+00  8.5064E-01  1.5823E+00  1.2124E+00  4.4891E-01  7.8933E-01  1.2030E+00  1.3395E+00
             2.1060E+00
 PARAMETER:  1.0630E-01  4.0753E-01  9.6973E-01 -6.1763E-02  5.5891E-01  2.9259E-01 -7.0092E-01 -1.3657E-01  2.8479E-01  3.9231E-01
             8.4481E-01
 GRADIENT:   4.7333E+01  3.1394E+01 -1.6710E+01  5.8725E+01  7.3546E+01  2.0279E+01 -1.0869E+01 -4.0314E+00 -2.7014E+01 -2.0457E+01
            -2.5250E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2855.25717359563        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.8284E-01  1.5056E+00  1.5664E+00  7.1306E-01  1.4192E+00  1.1442E+00  5.6026E-01  1.6081E-01  1.4311E+00  1.3375E+00
             2.3181E+00
 PARAMETER:  8.2695E-02  5.0918E-01  5.4877E-01 -2.3819E-01  4.5007E-01  2.3473E-01 -4.7935E-01 -1.7275E+00  4.5842E-01  3.9083E-01
             9.4076E-01
 GRADIENT:   5.7852E+00  4.0152E+00  2.2513E-01  4.4637E+00 -4.0466E+00 -1.7266E+00 -9.9754E-01 -6.8968E-02  1.8394E+00  1.6393E-01
             3.8923E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2855.85079480953        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      328
 NPARAMETR:  9.7928E-01  1.5634E+00  1.6324E+00  6.6864E-01  1.4913E+00  1.1495E+00  5.7134E-01  2.2107E-01  1.4697E+00  1.3705E+00
             2.3148E+00
 PARAMETER:  7.9060E-02  5.4689E-01  5.9006E-01 -3.0252E-01  4.9966E-01  2.3935E-01 -4.5978E-01 -1.4093E+00  4.8507E-01  4.1519E-01
             9.3933E-01
 GRADIENT:  -9.1287E+00 -2.2806E+01 -1.1957E+00 -6.9084E+00  4.8670E+00 -4.1202E+00 -7.4580E-01 -1.1084E-01 -1.0206E+00  3.0357E-02
             1.9786E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2859.00498278232        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  9.8594E-01  1.7682E+00  1.6575E+00  5.3153E-01  1.6163E+00  1.1555E+00  5.6908E-01  2.3401E+00  1.7095E+00  1.4226E+00
             2.2879E+00
 PARAMETER:  8.5837E-02  6.6997E-01  6.0530E-01 -5.3200E-01  5.8013E-01  2.4454E-01 -4.6374E-01  9.5021E-01  6.3623E-01  4.5251E-01
             9.2765E-01
 GRADIENT:   2.1149E+00 -2.2415E+01 -7.0011E+00 -1.2510E+01  4.3836E-01 -2.7934E+00 -3.3902E-01  9.5182E-01  5.0299E-01  8.4553E+00
             2.3901E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2864.82287848366        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      689
 NPARAMETR:  9.7961E-01  1.4139E+00  3.7036E+00  8.1642E-01  1.6211E+00  1.1767E+00  5.8656E-01  3.2781E+00  1.2193E+00  1.2321E+00
             2.2329E+00
 PARAMETER:  7.9403E-02  4.4633E-01  1.4093E+00 -1.0283E-01  5.8313E-01  2.6273E-01 -4.3347E-01  1.2873E+00  2.9828E-01  3.0869E-01
             9.0332E-01
 GRADIENT:  -7.8060E+00  4.2277E+01  2.0396E+00  1.4783E+01  3.0552E+01  3.4172E+00 -7.0190E-01 -3.2807E+00 -8.5554E+00 -1.4221E+01
            -3.4044E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2865.33042774918        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  9.7706E-01  1.4117E+00  3.7094E+00  8.0083E-01  1.6194E+00  1.1488E+00  5.3165E-01  3.2775E+00  1.2977E+00  1.2322E+00
             2.2339E+00
 PARAMETER:  7.6795E-02  4.4480E-01  1.4109E+00 -1.2211E-01  5.8204E-01  2.3871E-01 -5.3178E-01  1.2871E+00  3.6058E-01  3.0881E-01
             9.0374E-01
 GRADIENT:  -2.4202E+00  3.2515E+01  3.3823E+00  1.6270E+00  3.3359E+01 -1.2550E+00  9.4760E-01 -2.9250E+00 -1.2091E+00 -1.4564E+01
            -2.9280E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2865.47391417902        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  9.8435E-01  1.4117E+00  3.7094E+00  7.9643E-01  1.6194E+00  1.1659E+00  4.8212E-01  3.2775E+00  1.3389E+00  1.2322E+00
             2.2339E+00
 PARAMETER:  8.4229E-02  4.4480E-01  1.4109E+00 -1.2761E-01  5.8204E-01  2.5350E-01 -6.2957E-01  1.2871E+00  3.9185E-01  3.0881E-01
             9.0374E-01
 GRADIENT:   2.9836E-01  1.9198E+01  2.9530E+00 -1.1934E+00  2.9392E+01  1.1344E-01  1.0559E-02 -3.6760E+00 -1.4118E-01 -1.5393E+01
            -3.0868E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2866.06026143846        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1115
 NPARAMETR:  9.8371E-01  1.4069E+00  3.6698E+00  7.9759E-01  1.6077E+00  1.1641E+00  4.8054E-01  3.3059E+00  1.3384E+00  1.2513E+00
             2.2408E+00
 PARAMETER:  8.3572E-02  4.4140E-01  1.4001E+00 -1.2616E-01  5.7480E-01  2.5195E-01 -6.3284E-01  1.2957E+00  3.9148E-01  3.2420E-01
             9.0683E-01
 GRADIENT:   9.4573E+00  2.7053E+01  2.8545E+00 -1.3562E+00  2.6510E+01  4.5692E+00  5.9747E-01 -1.7168E+00  1.0200E+00 -1.0232E+01
            -2.0448E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2866.07633362229        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1293
 NPARAMETR:  9.8364E-01  1.4069E+00  3.6694E+00  7.9769E-01  1.6075E+00  1.1654E+00  4.7539E-01  3.3055E+00  1.3407E+00  1.2518E+00
             2.2415E+00
 PARAMETER:  8.3507E-02  4.4137E-01  1.4000E+00 -1.2604E-01  5.7471E-01  2.5308E-01 -6.4362E-01  1.2956E+00  3.9318E-01  3.2462E-01
             9.0715E-01
 GRADIENT:  -9.0736E-01  1.6526E+01  2.2791E+00 -2.5161E+00  2.2742E+01  4.5435E-02  1.1845E-01 -2.2231E+00  1.3579E-01 -1.0500E+01
            -2.1373E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2866.83019007052        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1433
 NPARAMETR:  9.8377E-01  1.3920E+00  3.5460E+00  8.0001E-01  1.5707E+00  1.1642E+00  4.6912E-01  3.3033E+00  1.3413E+00  1.3114E+00
             2.2624E+00
 PARAMETER:  8.3634E-02  4.3074E-01  1.3658E+00 -1.2313E-01  5.5151E-01  2.5204E-01 -6.5689E-01  1.2949E+00  3.9368E-01  3.7107E-01
             9.1643E-01
 GRADIENT:   8.7274E+00  1.6129E+01  1.7036E+00 -7.8614E+00  7.2998E+00  4.5277E+00  1.1054E+00  3.6472E-01  2.1602E+00  3.3302E+00
             5.2385E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2866.83427987906        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     1638
 NPARAMETR:  9.8379E-01  1.3918E+00  3.5472E+00  8.0011E-01  1.5709E+00  1.1644E+00  4.6875E-01  3.3023E+00  1.3339E+00  1.3115E+00
             2.2629E+00
 PARAMETER:  8.3660E-02  4.3063E-01  1.3662E+00 -1.2301E-01  5.5165E-01  2.5217E-01 -6.5769E-01  1.2946E+00  3.8808E-01  3.7116E-01
             9.1666E-01
 GRADIENT:  -7.9990E-01  6.1182E+00  1.1381E+00 -9.9648E+00  3.9518E+00 -1.0965E-02  5.2638E-01 -1.8044E-01  2.6657E-02  2.9537E+00
             3.9573E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2866.85316952102        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1782
 NPARAMETR:  9.8385E-01  1.3914E+00  3.5449E+00  8.0058E-01  1.5707E+00  1.1622E+00  4.5324E-01  3.3011E+00  1.3348E+00  1.3105E+00
             2.2632E+00
 PARAMETER:  8.3715E-02  4.3030E-01  1.3655E+00 -1.2242E-01  5.5152E-01  2.5034E-01 -6.9133E-01  1.2942E+00  3.8881E-01  3.7038E-01
             9.1678E-01
 GRADIENT:   8.7436E+00  1.9302E+01  1.6267E+00 -7.3887E+00  7.4959E+00  3.8446E+00  2.8625E-01  2.1963E-01 -4.0388E-01  3.1569E+00
             5.5098E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2866.88223553101        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1963            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8391E-01  1.3905E+00  3.5382E+00  8.0142E-01  1.5699E+00  1.1644E+00  4.4383E-01  3.3026E+00  1.3467E+00  1.3083E+00
             2.2627E+00
 PARAMETER:  8.3775E-02  4.2964E-01  1.3636E+00 -1.2137E-01  5.5099E-01  2.5223E-01 -7.1232E-01  1.2947E+00  3.9764E-01  3.6875E-01
             9.1655E-01
 GRADIENT:   8.8252E+00  1.9849E+01  1.4860E+00 -5.9428E+00  7.5729E+00  4.6026E+00  3.3856E-01  2.7849E-01  9.9691E-01  2.8617E+00
             5.1420E+00

0ITERATION NO.:   74    OBJECTIVE VALUE:  -2866.88359351811        NO. OF FUNC. EVALS.: 152
 CUMULATIVE NO. OF FUNC. EVALS.:     2115
 NPARAMETR:  9.8393E-01  1.3904E+00  3.5390E+00  8.0150E-01  1.5700E+00  1.1644E+00  4.4355E-01  3.3020E+00  1.3461E+00  1.3082E+00
             2.2630E+00
 PARAMETER:  8.3781E-02  4.2964E-01  1.3636E+00 -1.2125E-01  5.5099E-01  2.5223E-01 -7.1234E-01  1.2947E+00  3.9763E-01  3.6856E-01
             9.1656E-01
 GRADIENT:  -1.3947E+05  3.2464E+04 -2.0448E+04  2.3005E+05 -5.0622E+04 -1.0141E-02  1.6748E-02  2.1495E+04  1.5458E-01 -7.5683E+04
            -3.0436E+04
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         4.0         2.6         3.3         2.5         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2115
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2461E-03 -4.2203E-02 -3.3780E-02  1.8459E-02 -3.0342E-02
 SE:             2.9656E-02  1.2035E-02  1.6933E-02  2.6596E-02  2.3881E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3963E-01  4.5370E-04  4.6052E-02  4.8764E-01  2.0389E-01

 ETASHRINKSD(%)  6.4873E-01  5.9682E+01  4.3272E+01  1.0901E+01  1.9996E+01
 ETASHRINKVR(%)  1.2932E+00  8.3744E+01  6.7819E+01  2.0614E+01  3.5993E+01
 EBVSHRINKSD(%)  9.6082E-01  6.1946E+01  4.7368E+01  1.0619E+01  1.5772E+01
 EBVSHRINKVR(%)  1.9124E+00  8.5519E+01  7.2299E+01  2.0110E+01  2.9057E+01
 RELATIVEINF(%)  9.8054E+01  1.9870E+00  1.3344E+01  1.2036E+01  4.0223E+01
 EPSSHRINKSD(%)  1.7735E+01
 EPSSHRINKVR(%)  3.2325E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          883
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1622.8454496394520     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2866.8835935181082     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1244.0381438786562     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    58.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2866.884       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.39E+00  3.54E+00  8.02E-01  1.57E+00  1.16E+00  4.44E-01  3.30E+00  1.35E+00  1.31E+00  2.26E+00
 


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
+        7.20E+08
 
 TH 2
+        3.85E+02  1.95E+07
 
 TH 3
+       -4.87E+01  1.30E+03  2.99E+05
 
 TH 4
+        4.77E+04 -7.37E+03  9.68E+02  7.38E+08
 
 TH 5
+       -2.74E+02  7.10E+03  1.67E+06  5.42E+03  9.32E+06
 
 TH 6
+       -2.75E+03  4.52E+02 -5.63E+01  2.79E+03 -3.15E+02  1.42E+02
 
 TH 7
+       -1.58E+03  2.16E+02 -3.19E+01  2.27E+08 -1.78E+02  2.68E-01  4.20E+01
 
 TH 8
+        5.53E+01 -1.46E+03  2.28E+02 -1.10E+03  3.01E+01  6.35E+01  3.71E+01  3.80E+05
 
 TH 9
+       -4.12E+02  3.25E+01 -8.35E+00  4.58E+02 -4.67E+01 -2.05E-01  3.57E+01  1.06E+01  6.65E+01
 
 TH10
+        1.47E+08  5.33E+01 -9.76E+00 -1.49E+08 -6.57E+01 -5.63E+02 -3.20E+02  1.07E+01 -8.51E+01  3.00E+07
 
 TH11
+       -1.24E+02  2.97E+03 -4.74E+02  2.25E+03 -6.20E+01 -1.29E+02 -7.00E+01 -2.59E+02 -1.42E+01 -1.32E+01  1.62E+06
 
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
 #CPUT: Total CPU Time in Seconds,       72.192
Stop Time:
Sat Sep 25 02:19:02 CDT 2021
