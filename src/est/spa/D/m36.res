Wed Sep 29 19:58:57 CDT 2021
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
$DATA ../../../../data/spa/D/dat36.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13644.2897781826        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3304E+02  1.8630E+02 -2.0247E+01  1.3831E+02  2.0295E+02 -1.1418E+03 -5.9122E+02 -3.5774E+01 -9.3152E+02 -4.2119E+02
            -2.7173E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -660.767317703104        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.5653E+00  1.3256E+00  1.0291E+00  1.6772E+00  1.0340E+00  1.5490E+00  1.1706E+00  9.7992E-01  1.1982E+00  1.0468E+00
             1.4940E+01
 PARAMETER:  5.4811E-01  3.8188E-01  1.2870E-01  6.1712E-01  1.3344E-01  5.3760E-01  2.5750E-01  7.9720E-02  2.8084E-01  1.4577E-01
             2.8040E+00
 GRADIENT:   6.7419E+01  3.3826E+01  3.1909E+00  5.6013E+01 -1.8785E+01  2.3880E+01 -8.6803E-01  3.2317E+00  7.7503E+00  4.2472E+00
             1.1739E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -685.648019020396        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3349E+00  7.8496E-01  1.7404E+00  1.7837E+00  1.6896E+00  1.3137E+00  2.2487E+00  2.5265E-01  7.9902E-01  2.6697E+00
             1.3985E+01
 PARAMETER:  3.8886E-01 -1.4212E-01  6.5414E-01  6.7870E-01  6.2448E-01  3.7283E-01  9.1033E-01 -1.2758E+00 -1.2437E-01  1.0820E+00
             2.7380E+00
 GRADIENT:  -5.6725E+00  2.0975E+01  6.8124E+00  4.8661E+01 -1.7388E+01 -1.4521E+00  6.3444E+00  7.3436E-02  5.0368E+00  4.6605E+00
             1.3842E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -710.806412889644        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1064E+00  4.1113E-01  1.8214E+00  1.5746E+00  3.2233E+00  1.3358E+00  2.2809E+00  8.1009E-02  6.5267E-01  4.9949E+00
             1.0572E+01
 PARAMETER:  2.0114E-01 -7.8884E-01  6.9958E-01  5.5399E-01  1.2704E+00  3.8951E-01  9.2456E-01 -2.4132E+00 -3.2668E-01  1.7084E+00
             2.4582E+00
 GRADIENT:  -2.5138E+01  2.2499E+00  8.9994E+00 -6.2907E+00 -6.8268E+00  1.8398E+01  2.3997E+00  5.5380E-04  2.7576E+00 -9.6450E-01
             6.6655E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -726.156645865711        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0384E+00  3.7899E-01  5.7640E-01  1.3467E+00  6.8349E+00  1.1862E+00  6.5449E-01  1.0000E-02  3.0512E-01  7.9539E+00
             1.0419E+01
 PARAMETER:  1.3772E-01 -8.7025E-01 -4.5095E-01  3.9768E-01  2.0220E+00  2.7077E-01 -3.2389E-01 -6.4251E+00 -1.0870E+00  2.1737E+00
             2.4436E+00
 GRADIENT:  -1.9895E+01  1.5064E+01 -9.5087E+00  1.1338E+01 -9.9513E+00 -1.7910E+01  2.2967E-01  0.0000E+00  1.7595E+00  1.0237E+01
            -1.1854E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -739.776330300651        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.5205E-01  7.7325E-02  2.6814E-01  1.2057E+00  7.0179E+00  1.1905E+00  2.7884E-01  1.0000E-02  4.0589E-02  2.6626E+00
             1.0732E+01
 PARAMETER:  5.0859E-02 -2.4597E+00 -1.2163E+00  2.8707E-01  2.0485E+00  2.7437E-01 -1.1771E+00 -1.1086E+01 -3.1043E+00  1.0793E+00
             2.4732E+00
 GRADIENT:   2.8138E+01 -1.4051E+00 -6.0873E+00  6.1865E+01  1.7503E+01 -1.4037E+01  8.2889E-04  0.0000E+00 -1.6724E-03 -7.1136E+00
            -2.4802E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -747.517544363666        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  7.0977E-01  2.5017E-02  9.1698E-02  7.7633E-01  4.6925E+00  1.2671E+00  6.8867E-02  1.0000E-02  1.0000E-02  1.0706E+00
             1.0629E+01
 PARAMETER: -2.4281E-01 -3.5882E+00 -2.2893E+00 -1.5318E-01  1.6460E+00  3.3669E-01 -2.5756E+00 -1.6700E+01 -4.7970E+00  1.6821E-01
             2.4636E+00
 GRADIENT:  -3.0979E+01  3.5555E+00 -3.6548E+01  1.0896E+02 -2.5632E+01 -6.9450E+00  3.5969E-04  0.0000E+00  0.0000E+00  8.4139E+00
            -3.1025E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -747.640217179937        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  6.6493E-01  2.0613E-02  7.1553E-02  6.7914E-01  4.7370E+00  1.2710E+00  4.4933E-02  1.0000E-02  1.0000E-02  9.1265E-01
             1.0601E+01
 PARAMETER: -3.0807E-01 -3.7818E+00 -2.5373E+00 -2.8693E-01  1.6554E+00  3.3979E-01 -3.0026E+00 -1.8043E+01 -5.2134E+00  8.5990E-03
             2.4609E+00
 GRADIENT:  -1.6084E+01  3.1280E+00 -5.1348E+01  1.2732E+02 -2.2854E+01 -3.3411E+00  1.8544E-04  0.0000E+00  0.0000E+00  7.9967E+00
            -8.6004E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -747.642707326748        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      633
 NPARAMETR:  6.6021E-01  2.0230E-02  6.9842E-02  6.6959E-01  4.7411E+00  1.2707E+00  4.3010E-02  1.0000E-02  1.0000E-02  8.9877E-01
             1.0597E+01
 PARAMETER: -3.1520E-01 -3.8006E+00 -2.5615E+00 -3.0109E-01  1.6563E+00  3.3960E-01 -3.0463E+00 -1.8177E+01 -5.2553E+00 -6.7335E-03
             2.4606E+00
 GRADIENT:  -3.1387E+01  2.8713E+00 -7.1814E+01  1.2391E+02 -3.0495E+01 -5.6430E+00  1.5625E-04  0.0000E+00  0.0000E+00  7.9440E+00
            -2.9361E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -754.928847873691        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  6.6769E-01  1.4207E-02  6.6450E-02  6.1279E-01  4.8472E+00  1.3188E+00  3.0974E-02  1.0000E-02  1.0000E-02  5.0650E-01
             1.0834E+01
 PARAMETER: -3.0393E-01 -4.1540E+00 -2.6113E+00 -3.8974E-01  1.6784E+00  3.7669E-01 -3.3746E+00 -1.9151E+01 -5.8650E+00 -5.8023E-01
             2.4827E+00
 GRADIENT:   1.3959E+01  2.8621E-01 -2.0398E+01  2.8709E+01 -6.8950E+00  1.1987E+01  8.0988E-06  0.0000E+00  0.0000E+00  9.7264E-01
            -9.9523E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -757.557399725342        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  5.5476E-01  1.0000E-02  4.2573E-02  4.4289E-01  5.3561E+00  1.2182E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.7794E-02
             1.0825E+01
 PARAMETER: -4.8921E-01 -6.2410E+00 -3.0565E+00 -7.1443E-01  1.7782E+00  2.9737E-01 -4.5372E+00 -2.3934E+01 -8.6904E+00 -3.1756E+00
             2.4818E+00
 GRADIENT:  -1.5202E+01  0.0000E+00  9.2044E-01  5.1807E+00 -2.1550E+00 -2.0527E+00  2.2875E-07  0.0000E+00  0.0000E+00  7.9354E-03
             1.7037E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -757.649090269266        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1164
 NPARAMETR:  5.6931E-01  1.0000E-02  4.4451E-02  4.5500E-01  5.4169E+00  1.2275E+00  1.0747E-02  1.0000E-02  1.0000E-02  2.2612E-02
             1.0848E+01
 PARAMETER: -4.6332E-01 -5.9767E+00 -3.0134E+00 -6.8747E-01  1.7895E+00  3.0495E-01 -4.4331E+00 -2.3401E+01 -8.3662E+00 -3.6893E+00
             2.4839E+00
 GRADIENT:   8.0268E-02  0.0000E+00 -1.0394E+00  1.2250E+00  1.2268E-01  2.7934E-02  1.2740E-06  0.0000E+00  0.0000E+00  2.2474E-03
            -3.7792E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -757.649768654870        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1329
 NPARAMETR:  5.6919E-01  1.0000E-02  4.4471E-02  4.5463E-01  5.4278E+00  1.2268E+00  1.0747E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.0851E+01
 PARAMETER: -4.6354E-01 -5.9767E+00 -3.0129E+00 -6.8827E-01  1.7915E+00  3.0444E-01 -4.4331E+00 -2.3401E+01 -8.3662E+00 -4.6821E+00
             2.4843E+00
 GRADIENT:  -1.0394E-01  0.0000E+00  1.6090E-01 -4.3526E-01  1.0158E-02 -2.1763E-02  1.2473E-06  0.0000E+00  0.0000E+00  0.0000E+00
             9.8689E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1329
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1525E-03  1.9959E-06  1.2299E-04 -2.5734E-04  7.5377E-06
 SE:             2.8209E-02  4.0353E-06  2.4621E-04  3.5215E-04  6.7229E-05
 N:                     100         100         100         100         100

 P VAL.:         9.1102E-01  6.2087E-01  6.1741E-01  4.6493E-01  9.1073E-01

 ETASHRINKSD(%)  5.4956E+00  9.9986E+01  9.9175E+01  9.8820E+01  9.9775E+01
 ETASHRINKVR(%)  1.0689E+01  1.0000E+02  9.9993E+01  9.9986E+01  9.9999E+01
 EBVSHRINKSD(%)  5.7453E+00  9.9984E+01  9.9133E+01  9.8786E+01  9.9735E+01
 EBVSHRINKVR(%)  1.1161E+01  1.0000E+02  9.9992E+01  9.9985E+01  9.9999E+01
 RELATIVEINF(%)  3.5185E+00  9.1847E-08  4.4589E-05  8.6016E-05  3.1091E-05
 EPSSHRINKSD(%)  6.6098E+00
 EPSSHRINKVR(%)  1.2783E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -757.64976865487040     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -22.498942091132221     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -757.650       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.69E-01  1.00E-02  4.45E-02  4.55E-01  5.43E+00  1.23E+00  1.07E-02  1.00E-02  1.00E-02  1.00E-02  1.09E+01
 


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
+        1.98E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.63E+03  0.00E+00  2.47E+05
 
 TH 4
+       -4.25E+02  0.00E+00 -3.23E+04  5.00E+03
 
 TH 5
+        1.17E+01  0.00E+00 -4.26E+02  5.76E+01  2.02E+00
 
 TH 6
+       -6.79E+00  0.00E+00  3.47E+02 -6.36E+01  6.83E-01  1.01E+02
 
 TH 7
+        2.99E-03  0.00E+00  3.57E-02 -4.01E-03 -2.52E-04 -1.06E-02 -7.70E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.55E+01  0.00E+00  2.27E+02 -1.82E+01 -4.82E-01  1.02E+00 -5.71E-05  0.00E+00  0.00E+00  0.00E+00  3.84E+00
 
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
 #CPUT: Total CPU Time in Seconds,       23.033
Stop Time:
Wed Sep 29 19:59:22 CDT 2021
