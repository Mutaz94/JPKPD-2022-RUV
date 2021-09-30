Wed Sep 29 13:29:16 CDT 2021
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
$DATA ../../../../data/spa/A3/dat37.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -513.750568114158        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0661E+02  4.0711E+01  1.5858E+02 -6.0684E+01  2.5669E+01  3.4000E+01 -2.3107E+01 -6.6809E+01 -2.7529E+01 -7.2411E+01
            -2.0677E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1227.89258145150        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9136E-01  1.0467E+00  8.1039E-01  1.0558E+00  9.2584E-01  9.3314E-01  1.0538E+00  1.0250E+00  9.4015E-01  1.0810E+00
             1.8732E+00
 PARAMETER:  9.1321E-02  1.4560E-01 -1.1024E-01  1.5433E-01  2.2941E-02  3.0798E-02  1.5236E-01  1.2466E-01  3.8284E-02  1.7792E-01
             7.2765E-01
 GRADIENT:   7.7337E+01  4.2040E+01  3.0963E+01  4.0989E+01 -1.2751E+00 -1.1614E+01 -9.8678E+00 -5.8616E+00 -6.5685E+00 -5.9638E+00
            -4.4212E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1246.29911632647        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      199
 NPARAMETR:  9.8902E-01  9.6842E-01  4.1217E-01  1.0660E+00  6.1102E-01  1.0002E+00  1.5446E+00  4.1379E-01  7.3115E-01  7.6353E-01
             1.9244E+00
 PARAMETER:  8.8957E-02  6.7907E-02 -7.8631E-01  1.6391E-01 -3.9263E-01  1.0023E-01  5.3476E-01 -7.8241E-01 -2.1314E-01 -1.6980E-01
             7.5462E-01
 GRADIENT:  -4.1185E+01  5.1785E+01 -3.1460E+01  8.5177E+01  5.8216E+01  3.6120E+00  2.6209E+01 -8.9844E-01 -2.5084E+01 -6.8666E-01
            -3.8724E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1335.16347397740        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0211E+00  9.3284E-01  4.0179E-01  1.0528E+00  5.7250E-01  9.6727E-01  1.1446E+00  1.6523E-01  7.7257E-01  5.3078E-01
             3.0724E+00
 PARAMETER:  1.2089E-01  3.0480E-02 -8.1183E-01  1.5141E-01 -4.5775E-01  6.6727E-02  2.3505E-01 -1.7004E+00 -1.5804E-01 -5.3341E-01
             1.2225E+00
 GRADIENT:  -2.2273E+00  2.4360E+01  1.6644E+01  3.6430E+00 -3.3700E+01  3.1956E+00 -9.4680E-01  5.7862E-03  3.9246E+00  4.2209E+00
            -4.6730E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1337.22229980001        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  1.0240E+00  8.1317E-01  3.8354E-01  1.0903E+00  5.3041E-01  9.5141E-01  1.3174E+00  3.7873E-01  7.0656E-01  2.5511E-01
             3.1369E+00
 PARAMETER:  1.2371E-01 -1.0682E-01 -8.5830E-01  1.8649E-01 -5.3411E-01  5.0189E-02  3.7563E-01 -8.7092E-01 -2.4735E-01 -1.2661E+00
             1.2432E+00
 GRADIENT:   5.8464E+00 -6.7617E+00 -3.9168E+00 -1.5849E+00  9.4143E+00 -1.8923E+00  9.1015E-01 -5.9962E-01  2.0621E+00  3.5473E-01
             4.5645E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1337.79151615022        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  1.0223E+00  7.6989E-01  3.7492E-01  1.1094E+00  5.0303E-01  9.6034E-01  1.3712E+00  6.3986E-01  6.6289E-01  1.8230E-01
             3.0805E+00
 PARAMETER:  1.2204E-01 -1.6150E-01 -8.8104E-01  2.0379E-01 -5.8710E-01  5.9528E-02  4.1569E-01 -3.4651E-01 -3.1115E-01 -1.6021E+00
             1.2251E+00
 GRADIENT:   1.2676E+00 -2.5701E-01 -5.3385E-01  4.4076E+00  2.2276E+00 -3.6569E-02 -1.9067E-01 -2.1461E-01 -1.9192E+00  4.0278E-01
             8.7166E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1338.22178358585        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0176E+00  9.7818E-01  3.2275E-01  9.9000E-01  5.4387E-01  9.6318E-01  1.1264E+00  7.3774E-01  6.8004E-01  1.1207E-01
             3.0830E+00
 PARAMETER:  1.1746E-01  7.7938E-02 -1.0309E+00  8.9952E-02 -5.0904E-01  6.2481E-02  2.1906E-01 -2.0417E-01 -2.8560E-01 -2.0886E+00
             1.2259E+00
 GRADIENT:  -3.9995E+00  9.7049E+00  1.1691E-01  1.6921E+01  1.1742E-01  1.2882E+00  8.6906E-01 -8.8034E-02 -4.0746E+00  2.9435E-01
            -1.0276E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1352.60695700559        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  1.0001E+00  1.6258E+00  1.5337E-01  5.9737E-01  7.3622E-01  9.1741E-01  7.4373E-01  2.6923E+00  7.1905E-01  1.0000E-02
             2.9841E+00
 PARAMETER:  1.0013E-01  5.8603E-01 -1.7749E+00 -4.1522E-01 -2.0623E-01  1.3798E-02 -1.9608E-01  1.0904E+00 -2.2982E-01 -6.9108E+00
             1.1933E+00
 GRADIENT:  -1.0691E+01  5.1791E+01  4.6561E+00  4.4239E+01 -1.6710E+01 -1.2823E+01 -1.2993E+01 -9.1908E-01 -7.7959E+00  0.0000E+00
             3.6695E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1357.53595138973        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1270
 NPARAMETR:  9.9295E-01  1.7024E+00  1.2533E-01  5.4039E-01  7.7571E-01  9.5674E-01  8.0133E-01  2.7726E+00  7.7358E-01  1.0000E-02
             2.7317E+00
 PARAMETER:  9.2926E-02  6.3205E-01 -1.9768E+00 -5.1545E-01 -1.5398E-01  5.5771E-02 -1.2149E-01  1.1198E+00 -1.5672E-01 -7.9980E+00
             1.1049E+00
 GRADIENT:  -1.1189E+01  3.3088E+01 -8.5396E+00  3.4360E+01  2.8979E+00  1.0705E+00  7.0257E-01 -1.8097E+00 -1.2534E+01  0.0000E+00
            -4.2137E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1360.82504222196        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1436
 NPARAMETR:  9.9424E-01  1.6942E+00  1.2842E-01  5.2765E-01  7.7743E-01  9.5418E-01  7.9489E-01  2.7675E+00  1.2367E+00  1.0000E-02
             2.7400E+00
 PARAMETER:  9.4220E-02  6.2720E-01 -1.9524E+00 -5.3932E-01 -1.5176E-01  5.3101E-02 -1.2955E-01  1.1180E+00  3.1241E-01 -8.0244E+00
             1.1080E+00
 GRADIENT:   4.7345E+01  5.8648E+01  1.1268E+01  2.4102E+01  6.6494E+00  5.5908E+00  1.0979E+01  6.7850E+00  6.6307E+00  0.0000E+00
             3.6757E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1361.88793419903        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1594
 NPARAMETR:  9.9213E-01  1.6946E+00  1.2298E-01  5.2030E-01  7.7312E-01  9.5205E-01  7.7866E-01  2.6864E+00  1.2453E+00  1.0000E-02
             2.6751E+00
 PARAMETER:  9.2099E-02  6.2746E-01 -1.9957E+00 -5.5334E-01 -1.5732E-01  5.0860E-02 -1.5018E-01  1.0882E+00  3.1934E-01 -8.0244E+00
             1.0840E+00
 GRADIENT:  -2.6774E+00 -1.4394E+01  3.6968E+00  7.5930E+00  1.9363E+00  3.6850E-01  5.4933E+00  2.4577E+00  3.0249E+00  0.0000E+00
             1.5009E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1363.40581482115        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1768
 NPARAMETR:  9.8759E-01  1.7466E+00  9.7280E-02  4.7727E-01  7.7761E-01  9.5378E-01  7.3946E-01  2.2434E+00  1.4486E+00  1.0000E-02
             2.5485E+00
 PARAMETER:  8.7372E-02  6.5725E-01 -2.2331E+00 -6.4053E-01 -1.5132E-01  5.1678E-02 -2.0211E-01  9.0672E-01  4.7125E-01 -8.0244E+00
             1.0368E+00
 GRADIENT:  -8.7516E+02 -4.3515E+00 -7.4548E+01 -2.7178E+02  1.1439E+03 -6.4944E-01 -4.3496E+02 -1.0058E+02  3.7051E+02  0.0000E+00
             1.5897E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1768
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8667E-03 -1.9446E-02  1.6100E-03  1.1043E-02 -3.3870E-04
 SE:             2.9165E-02  2.5249E-02  1.4625E-02  2.1488E-02  2.2768E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4897E-01  4.4120E-01  9.1234E-01  6.0731E-01  1.3685E-01

 ETASHRINKSD(%)  2.2939E+00  1.5412E+01  5.1005E+01  2.8013E+01  9.9237E+01
 ETASHRINKVR(%)  4.5352E+00  2.8448E+01  7.5995E+01  4.8178E+01  9.9994E+01
 EBVSHRINKSD(%)  2.6346E+00  1.6557E+01  4.8174E+01  2.7848E+01  9.9164E+01
 EBVSHRINKVR(%)  5.1998E+00  3.0372E+01  7.3140E+01  4.7941E+01  9.9993E+01
 RELATIVEINF(%)  9.2313E+01  9.4257E+00  1.5201E+01  1.1157E+01  9.6861E-04
 EPSSHRINKSD(%)  3.3905E+01
 EPSSHRINKVR(%)  5.6314E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1363.4058148211536     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -628.25498825741545     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1363.406       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  1.75E+00  9.70E-02  4.77E-01  7.78E-01  9.53E-01  7.39E-01  2.24E+00  1.45E+00  1.00E-02  2.55E+00
 


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
+        4.49E+05
 
 TH 2
+       -4.29E+01  4.71E+02
 
 TH 3
+       -3.52E+02  4.36E+02  9.53E+04
 
 TH 4
+        1.46E+05  1.64E+02 -3.60E+02  4.87E+04
 
 TH 5
+        2.20E+02 -5.37E+02  7.88E+03 -1.25E+03  3.15E+05
 
 TH 6
+       -3.54E+02 -9.10E+00 -1.52E+02 -1.32E+02  2.95E+02  1.96E+02
 
 TH 7
+        2.96E+05 -4.18E+00 -4.26E+02  9.63E+04  7.57E+02 -2.37E+02  1.96E+05
 
 TH 8
+        2.25E+04 -1.35E+00  5.23E+02  9.24E+01 -1.19E+03 -1.61E+01  1.48E+04  1.14E+03
 
 TH 9
+       -6.49E+04 -1.27E+01  5.90E+01 -2.11E+04  8.75E+01  5.13E+01 -4.28E+04 -3.25E+03  9.41E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.94E+00 -1.10E+01  4.14E+02 -7.14E+01 -9.75E+01  1.72E+01  4.97E+01 -5.13E+01  1.42E+01  0.00E+00  6.64E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.582
Stop Time:
Wed Sep 29 13:29:49 CDT 2021
