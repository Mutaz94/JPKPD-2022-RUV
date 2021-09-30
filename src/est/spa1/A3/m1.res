Wed Sep 29 23:52:59 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -78.0434722687312        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6219E+02  5.5958E+01  2.7594E+02  6.7988E+00  1.7724E+02  4.7928E+01 -1.0184E+02 -4.0044E+02 -1.0447E+02 -1.3544E+02
            -3.2518E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1474.44432836879        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0461E+00  1.0420E+00  8.9164E-01  1.1112E+00  9.9622E-01  8.3479E-01  1.0374E+00  1.0946E+00  9.9679E-01  9.3356E-01
             4.7073E+00
 PARAMETER:  1.4507E-01  1.4115E-01 -1.4694E-02  2.0541E-01  9.6208E-02 -8.0575E-02  1.3669E-01  1.9038E-01  9.6787E-02  3.1247E-02
             1.6491E+00
 GRADIENT:   7.5537E+01 -5.8896E+00 -2.3169E+01  2.1248E+01  1.6973E-02 -1.4060E+01  7.4822E+00  9.3747E+00  1.6794E+01  1.9529E+01
             2.7058E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1498.83091113219        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0221E+00  8.6063E-01  3.2848E-01  1.1114E+00  4.3509E-01  9.2377E-01  1.0393E+00  9.2880E-01  1.0061E+00  2.2976E-01
             4.1301E+00
 PARAMETER:  1.2188E-01 -5.0086E-02 -1.0133E+00  2.0564E-01 -7.3220E-01  2.0703E-02  1.3852E-01  2.6142E-02  1.0608E-01 -1.3707E+00
             1.5183E+00
 GRADIENT:   2.9291E+00  1.1079E+02  6.8896E+01  4.6208E+01 -1.4696E+02  1.8481E+00 -3.9285E+00  1.0167E+01 -5.4641E+00  1.6720E+00
             2.0990E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1539.89914194426        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      256
 NPARAMETR:  9.5646E-01  7.2885E-01  2.4043E-01  1.0899E+00  3.7667E-01  9.5727E-01  1.0069E+00  3.2333E-01  1.0469E+00  1.6266E-01
             3.0508E+00
 PARAMETER:  5.5483E-02 -2.1629E-01 -1.3253E+00  1.8608E-01 -8.7639E-01  5.6334E-02  1.0691E-01 -1.0291E+00  1.4581E-01 -1.7161E+00
             1.2154E+00
 GRADIENT:  -1.1699E+02  3.9057E+01 -2.6433E+00  4.4792E+01 -1.9566E+01 -3.0616E+00 -2.2080E+01 -1.7573E+00 -3.7187E+01 -1.4949E+00
            -1.1602E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1554.86796990285        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      434
 NPARAMETR:  1.0403E+00  5.4850E-01  1.7608E-01  1.0808E+00  2.7767E-01  8.7153E-01  1.3076E+00  1.1099E+00  1.2042E+00  7.1973E-02
             2.7233E+00
 PARAMETER:  1.3952E-01 -5.0057E-01 -1.6368E+00  1.7767E-01 -1.1813E+00 -3.7509E-02  3.6823E-01  2.0426E-01  2.8581E-01 -2.5315E+00
             1.1018E+00
 GRADIENT:   1.0639E+02  5.3726E+01 -2.9280E+01  2.1482E+01 -1.7483E+01 -3.9288E+01 -3.1090E-01  5.5178E+00 -4.1868E+01 -2.7558E-01
             9.0049E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1564.57275617783        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      621             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0362E+00  4.9690E-01  1.6883E-01  1.0834E+00  2.6090E-01  9.8558E-01  1.3171E+00  1.2106E+00  1.5383E+00  5.9325E-02
             2.6553E+00
 PARAMETER:  1.3560E-01 -5.9936E-01 -1.6789E+00  1.8010E-01 -1.2436E+00  8.5479E-02  3.7543E-01  2.9115E-01  5.3066E-01 -2.7247E+00
             1.0766E+00
 GRADIENT:   1.4401E+02  4.7006E+01  2.6776E+01  2.0515E+01  1.1548E+02  1.1604E+01  3.3973E+00  2.0457E+00  2.4417E+01  1.2768E-03
             4.4038E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1564.96744449700        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  1.0362E+00  4.9654E-01  1.6900E-01  1.0834E+00  2.6099E-01  9.7730E-01  1.3171E+00  1.2106E+00  1.4715E+00  8.8133E-02
             2.6533E+00
 PARAMETER:  1.3557E-01 -6.0009E-01 -1.6779E+00  1.8007E-01 -1.2433E+00  7.7043E-02  3.7546E-01  2.9111E-01  4.8631E-01 -2.3289E+00
             1.0758E+00
 GRADIENT:   8.6611E+01  3.6524E+01 -6.6252E+00  6.7102E+00 -3.7916E+00  5.6178E+00  1.6319E+00  1.6235E+00 -1.1937E-01 -2.0763E-01
             3.0979E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1565.40752032200        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      955
 NPARAMETR:  1.0361E+00  4.9566E-01  1.6944E-01  1.0833E+00  2.6122E-01  9.6110E-01  1.3171E+00  1.2106E+00  1.4583E+00  2.2477E-01
             2.6473E+00
 PARAMETER:  1.3547E-01 -6.0187E-01 -1.6752E+00  1.8003E-01 -1.2424E+00  6.0322E-02  3.7545E-01  2.9108E-01  4.7730E-01 -1.3927E+00
             1.0735E+00
 GRADIENT:   8.8810E+01  3.5592E+01 -8.1499E+00  5.3858E+00  2.7939E+00 -2.6611E-01  4.6568E+00  3.8477E+00  9.2616E-01  1.9526E-01
             3.1499E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1567.07304841531        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  1.0339E+00  4.7865E-01  1.7918E-01  1.0826E+00  2.6323E-01  9.6142E-01  1.3160E+00  1.2097E+00  1.4345E+00  2.4779E-01
             2.5285E+00
 PARAMETER:  1.3338E-01 -6.3679E-01 -1.6194E+00  1.7938E-01 -1.2347E+00  6.0654E-02  3.7462E-01  2.9036E-01  4.6078E-01 -1.2952E+00
             1.0276E+00
 GRADIENT:   8.7713E+01  2.4750E+01  1.7123E+01 -4.0261E+00 -2.7314E+00  7.3932E-02  9.6946E-01  4.0436E-02 -3.5120E-01  1.5526E-01
            -2.2617E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1569.28554056631        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  9.8374E-01  4.5308E-01  1.7368E-01  1.0880E+00  2.6159E-01  9.6087E-01  1.3073E+00  1.2097E+00  1.4413E+00  2.4023E-01
             2.5325E+00
 PARAMETER:  8.3602E-02 -6.9168E-01 -1.6506E+00  1.8432E-01 -1.2410E+00  6.0081E-02  3.6797E-01  2.9038E-01  4.6553E-01 -1.3261E+00
             1.0292E+00
 GRADIENT:  -2.5160E+01 -1.0026E+01 -1.4281E+01 -4.7442E-01  6.8561E+01  2.6359E+00 -2.1302E+00  6.0394E-02 -1.5046E+00 -7.4851E-02
             7.5828E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1569.41500802958        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1436
 NPARAMETR:  9.9616E-01  4.5274E-01  1.7407E-01  1.0878E+00  2.6181E-01  9.5179E-01  1.3545E+00  1.2093E+00  1.4479E+00  2.3170E-01
             2.5353E+00
 PARAMETER:  9.6150E-02 -6.9244E-01 -1.6483E+00  1.8412E-01 -1.2401E+00  5.0584E-02  4.0345E-01  2.9005E-01  4.7008E-01 -1.3623E+00
             1.0303E+00
 GRADIENT:   3.7900E+00 -8.8061E+00 -1.4117E+01 -1.4187E+00  6.9137E+01 -3.2684E-01  9.9483E-01 -1.2290E-01  7.5232E-01  8.0581E-02
             8.8102E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1570.01506706997        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1614
 NPARAMETR:  9.9513E-01  4.5250E-01  1.7662E-01  1.0874E+00  2.5816E-01  9.5298E-01  1.3335E+00  1.2087E+00  1.4505E+00  2.4161E-01
             2.5364E+00
 PARAMETER:  9.5117E-02 -6.9297E-01 -1.6337E+00  1.8380E-01 -1.2542E+00  5.1840E-02  3.8783E-01  2.8955E-01  4.7192E-01 -1.3204E+00
             1.0308E+00
 GRADIENT:   6.8714E-02  8.4064E+00  1.2820E+01 -1.0443E+00  1.3765E+01  9.4478E-02  3.8775E-02 -2.4305E-01 -1.9701E-01  4.0242E-03
             7.4418E+00

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1570.04933066310        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1712
 NPARAMETR:  9.9496E-01  4.4959E-01  1.7565E-01  1.0885E+00  2.5953E-01  9.5282E-01  1.3385E+00  1.2104E+00  1.4521E+00  2.4128E-01
             2.5278E+00
 PARAMETER:  9.5126E-02 -7.0022E-01 -1.6422E+00  1.8434E-01 -1.2511E+00  5.1715E-02  3.8765E-01  2.9045E-01  4.7270E-01 -1.3198E+00
             1.0255E+00
 GRADIENT:   1.1456E+00 -1.2947E+00 -2.1240E+03 -5.2919E-01 -2.7492E+03  4.8194E-02 -4.9969E-01 -1.7811E-01 -2.0148E-01  2.3542E-02
            -3.4758E+03
 NUMSIGDIG:         2.3         2.5         2.3         2.2         2.3         2.9         1.6         2.3         2.8         2.4
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1712
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.7574E-04  1.2291E-02 -9.5129E-03 -7.2392E-03  4.6585E-03
 SE:             2.9172E-02  2.0065E-02  2.1138E-02  2.7465E-02  8.3267E-03
 N:                     100         100         100         100         100

 P VAL.:         9.7879E-01  5.4018E-01  6.5268E-01  7.9210E-01  5.7585E-01

 ETASHRINKSD(%)  2.2696E+00  3.2778E+01  2.9185E+01  7.9894E+00  7.2105E+01
 ETASHRINKVR(%)  4.4878E+00  5.4812E+01  4.9852E+01  1.5340E+01  9.2218E+01
 EBVSHRINKSD(%)  2.4112E+00  3.2884E+01  2.8541E+01  7.3326E+00  7.3638E+01
 EBVSHRINKVR(%)  4.7644E+00  5.4954E+01  4.8936E+01  1.4128E+01  9.3050E+01
 RELATIVEINF(%)  9.5044E+01  9.8194E+00  1.1573E+01  7.7121E+01  6.7379E-01
 EPSSHRINKSD(%)  3.2776E+01
 EPSSHRINKVR(%)  5.4809E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1570.0493306630981     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -651.11079745842540     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.85
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1570.049       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  4.49E-01  1.75E-01  1.09E+00  2.59E-01  9.53E-01  1.33E+00  1.21E+00  1.45E+00  2.42E-01  2.52E+00
 


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
+        1.18E+03
 
 TH 2
+       -2.18E+01  1.42E+03
 
 TH 3
+       -5.69E+02 -9.62E+05  1.06E+06
 
 TH 4
+       -5.50E+00  9.25E+01  1.25E+03  3.81E+02
 
 TH 5
+       -3.12E+02 -3.36E+03 -2.00E+04  1.07E+03  8.49E+05
 
 TH 6
+        2.91E+00 -9.36E+00 -9.30E+01 -4.21E+00 -8.09E+01  1.97E+02
 
 TH 7
+        6.39E-01  4.82E+01 -4.26E+02 -2.70E+00 -4.05E+02  6.37E-02  2.58E+01
 
 TH 8
+       -2.50E+06  7.91E+05 -8.62E+05 -1.28E+00  1.61E+03  3.25E+00  1.70E+00  7.08E+05
 
 TH 9
+        1.12E+01 -2.66E+01 -6.93E+02 -1.08E+01 -5.02E+02  1.64E+00  4.55E+00 -2.37E+00  6.00E+01
 
 TH10
+       -1.16E+00  9.49E+00 -7.64E+02 -4.84E+00 -6.37E+02 -9.79E-02  1.24E+01  1.09E+01  9.69E+00  3.70E+01
 
 TH11
+       -7.35E+01 -1.10E+05  1.19E+05  1.66E+02  3.88E+03 -1.03E+01 -3.88E+01  2.14E+02 -7.80E+01 -8.06E+01  1.36E+04
 
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
 #CPUT: Total CPU Time in Seconds,       37.509
Stop Time:
Wed Sep 29 23:53:39 CDT 2021
