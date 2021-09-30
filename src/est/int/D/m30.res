Wed Sep 29 08:28:02 CDT 2021
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
$DATA ../../../../data/int/D/dat30.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E22.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   21087.8698796031        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5852E+02  3.3107E+02 -2.2362E+01  6.6388E+00  1.3251E+02 -1.1960E+03 -7.2855E+02 -9.3031E+01 -1.4146E+03 -7.0626E+02
            -4.5512E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1076.52597523024        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1955E+00  1.9970E+00  9.5780E-01  3.3550E+00  7.8438E-01  5.2694E+00  3.1727E+00  1.0313E+00  3.3560E+00  1.7572E+00
             1.2468E+01
 PARAMETER:  2.7857E-01  7.9165E-01  5.6888E-02  1.3104E+00 -1.4286E-01  1.7619E+00  1.2546E+00  1.3079E-01  1.3108E+00  6.6370E-01
             2.6231E+00
 GRADIENT:  -1.1559E+01  7.7348E+01 -2.2532E+01  1.5346E+02 -4.0356E+01  2.1263E+02  5.1465E+00  5.6064E+00  7.4255E+00 -2.7715E+00
             5.6974E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1216.41563640145        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  6.1218E-01  2.7526E+00  3.1760E+03  3.1765E+00  3.8714E+00  4.1924E+00  7.3351E+00  2.0164E+00  2.9187E+00  2.8006E+00
             1.1427E+01
 PARAMETER: -3.9073E-01  1.1125E+00  8.1634E+00  1.2558E+00  1.4536E+00  1.5333E+00  2.0927E+00  8.0133E-01  1.1711E+00  1.1298E+00
             2.5360E+00
 GRADIENT:  -8.4813E+01  5.5197E+01 -4.0360E-02  1.3450E+02  3.7882E+01  1.9394E+02  5.2203E+01  2.5986E-05  2.5999E+01  8.0609E+01
             5.8655E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1420.44456469817        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0812E+00  1.1887E+00  2.5812E+02  1.4799E+00  2.4160E+00  2.0821E+00  4.0783E+00  3.0791E+01  1.4442E+00  1.4977E+00
             1.0256E+01
 PARAMETER:  1.7805E-01  2.7285E-01  5.6534E+00  4.9197E-01  9.8209E-01  8.3338E-01  1.5057E+00  3.5272E+00  4.6758E-01  5.0395E-01
             2.4279E+00
 GRADIENT:  -7.5839E+01 -1.5498E+01 -4.4206E-01  2.9527E+01 -6.4092E+00  2.5300E+01 -4.6310E+01  2.4373E+01  2.9086E+00  4.2585E+01
             5.8073E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1538.59837082543        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.2659E+00  1.1332E+00  9.6398E+01  1.0083E+00  2.1437E+00  1.5618E+00  4.3958E+00  1.5287E+01  1.1069E+00  8.4027E-01
             8.2190E+00
 PARAMETER:  3.3577E-01  2.2505E-01  4.6685E+00  1.0827E-01  8.6255E-01  5.4581E-01  1.5807E+00  2.8270E+00  2.0152E-01 -7.4027E-02
             2.2064E+00
 GRADIENT:   5.9795E+01 -3.3446E+01  1.3431E+00 -5.2820E+01 -1.5796E+01 -1.0451E+02  5.7918E+01  1.4134E+01  6.1517E+00  1.8410E+01
             3.3070E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1609.11907281958        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      492
 NPARAMETR:  1.2255E+00  1.3645E+00  1.1051E+01  1.1480E+00  1.9849E+00  2.2949E+00  4.9429E+00  5.3236E+00  1.1471E+00  3.4212E-01
             6.4014E+00
 PARAMETER:  3.0331E-01  4.1082E-01  2.5026E+00  2.3800E-01  7.8556E-01  9.3068E-01  1.6980E+00  1.7721E+00  2.3724E-01 -9.7258E-01
             1.9565E+00
 GRADIENT:   8.0316E+00  9.3091E+00 -1.1194E+01  3.3191E+01 -3.4283E+01  6.8351E-01 -2.1266E+01  2.9143E+01 -6.2454E+00  2.2016E+00
            -1.7849E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1625.37251235363        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      668
 NPARAMETR:  1.2297E+00  9.8818E-01  5.0819E+01  1.3343E+00  2.2420E+00  2.2934E+00  5.8226E+00  4.4028E+00  1.2255E+00  4.0357E-01
             6.9069E+00
 PARAMETER:  3.0679E-01  8.8109E-02  4.0283E+00  3.8837E-01  9.0738E-01  9.3005E-01  1.8617E+00  1.5822E+00  3.0334E-01 -8.0740E-01
             2.0325E+00
 GRADIENT:   3.4590E+00  3.2245E+00  8.9260E-01  8.8554E-01 -4.2570E+00  2.6200E+00 -2.5270E+00 -1.2381E-01 -4.3145E+00  2.0205E+00
            -6.3803E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1628.71421353986        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.2321E+00  7.0458E-01  2.7785E+01  1.5337E+00  2.1857E+00  2.2849E+00  6.3490E+00  3.0239E+00  1.5167E+00  2.8854E-01
             6.9334E+00
 PARAMETER:  3.0874E-01 -2.5016E-01  3.4245E+00  5.2768E-01  8.8194E-01  9.2634E-01  1.9483E+00  1.2065E+00  5.1654E-01 -1.1429E+00
             2.0363E+00
 GRADIENT:   3.7233E+00  3.4844E-01  1.2190E+00 -1.9926E+00  1.2162E+01  9.5870E-01 -1.5713E-02 -1.2588E+00 -6.5408E-01  1.2836E+00
             4.7121E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1630.18729333679        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  1.2185E+00  6.0539E-01  1.7490E+01  1.5921E+00  2.0343E+00  2.2769E+00  6.5648E+00  2.2964E+00  1.6311E+00  2.2336E-01
             6.8939E+00
 PARAMETER:  2.9761E-01 -4.0188E-01  2.9616E+00  5.6507E-01  8.1014E-01  9.2282E-01  1.9817E+00  9.3133E-01  5.8923E-01 -1.3990E+00
             2.0306E+00
 GRADIENT:  -4.6431E-02 -4.4513E-01 -6.4731E-01  8.4320E-01  1.0174E+00 -8.9847E-01  3.7123E-01  6.5381E-01  1.0502E+00  6.9961E-01
            -2.1264E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1630.25740683876        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  1.2176E+00  6.3249E-01  1.6206E+01  1.5718E+00  2.0160E+00  2.2881E+00  6.5132E+00  1.9048E+00  1.5985E+00  1.7464E-01
             6.8954E+00
 PARAMETER:  2.9688E-01 -3.5809E-01  2.8854E+00  5.5225E-01  8.0110E-01  9.2773E-01  1.9738E+00  7.4437E-01  5.6907E-01 -1.6450E+00
             2.0309E+00
 GRADIENT:  -2.7520E-01  3.0950E-01  1.8618E-01  8.5020E-01 -2.1731E+00  9.1945E-01  5.9270E-01  7.2698E-02 -7.6567E-01  3.8911E-01
            -2.2591E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1630.26373156807        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1377
 NPARAMETR:  1.2177E+00  6.3531E-01  1.5482E+01  1.5704E+00  2.0089E+00  2.2845E+00  6.4900E+00  1.6358E+00  1.6114E+00  1.4415E-01
             6.9010E+00
 PARAMETER:  2.9698E-01 -3.5364E-01  2.8397E+00  5.5136E-01  7.9761E-01  9.2616E-01  1.9703E+00  5.9211E-01  5.7710E-01 -1.8369E+00
             2.0317E+00
 GRADIENT:  -2.8658E-01  3.0938E-01 -1.5676E-02  5.0075E-01 -7.7106E-01  2.6921E-01  3.1611E-01 -2.4660E-02 -1.6497E-01  2.5608E-01
            -5.2341E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1630.31123469562        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1556
 NPARAMETR:  1.2181E+00  6.2317E-01  1.4386E+01  1.5752E+00  1.9939E+00  2.2818E+00  6.4897E+00  8.9098E-01  1.6288E+00  6.5197E-02
             6.9083E+00
 PARAMETER:  2.9733E-01 -3.7293E-01  2.7662E+00  5.5437E-01  7.9010E-01  9.2495E-01  1.9702E+00 -1.5439E-02  5.8786E-01 -2.6303E+00
             2.0327E+00
 GRADIENT:  -1.6728E-01  2.5901E-02 -2.4010E-01  6.9866E-02  7.4762E-01 -2.7910E-01 -1.3416E-01 -5.6287E-02  3.8773E-01  4.6740E-02
             9.7760E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1631.15689623364        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1738
 NPARAMETR:  1.2241E+00  5.6266E-01  1.5152E+01  1.6109E+00  1.9888E+00  2.2889E+00  6.8462E+00  1.3267E+00  1.6214E+00  1.9767E-02
             6.9102E+00
 PARAMETER:  3.0224E-01 -4.7508E-01  2.8181E+00  5.7678E-01  7.8755E-01  9.2808E-01  2.0237E+00  3.8268E-01  5.8327E-01 -3.8238E+00
             2.0330E+00
 GRADIENT:   1.8125E+00  3.1788E-01 -2.1133E-01 -2.3577E-01 -1.7068E+00  8.6311E-01  6.6311E+00 -2.4243E-02 -7.8415E-01  4.2705E-03
             1.1059E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1631.78527385903        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1904             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2207E+00  4.8895E-01  1.5829E+01  1.6549E+00  1.9867E+00  2.3022E+00  7.2815E+00  1.4795E+00  1.6077E+00  1.0000E-02
             6.8978E+00
 PARAMETER:  2.9940E-01 -6.1549E-01  2.8619E+00  6.0376E-01  7.8646E-01  9.3387E-01  2.0853E+00  4.9171E-01  5.7479E-01 -4.6472E+00
             2.0312E+00
 GRADIENT:   3.6539E+01  7.7010E+00  1.5122E-01  4.3990E+01  3.6358E+00  6.2204E+01  2.6543E+02  8.2539E-02  6.9920E+00  0.0000E+00
             3.3448E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1631.82037252197        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2079
 NPARAMETR:  1.2236E+00  4.8176E-01  1.6971E+01  1.6789E+00  2.0039E+00  2.3014E+00  7.2450E+00  1.3782E+00  1.6587E+00  1.0000E-02
             6.9179E+00
 PARAMETER:  3.0181E-01 -6.3032E-01  2.9315E+00  6.1811E-01  7.9509E-01  9.3350E-01  2.0803E+00  4.2075E-01  6.0601E-01 -4.6472E+00
             2.0341E+00
 GRADIENT:   1.5104E+00  2.9656E-01 -1.5753E-01  2.3866E+00 -1.9231E+00  2.8850E+00  1.1432E+01 -3.0798E-02 -2.7343E-01  0.0000E+00
             1.8200E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1632.15980215662        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     2248             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2196E+00  4.2035E-01  1.7456E+01  1.7057E+00  2.0018E+00  2.3019E+00  7.5927E+00  1.5930E+00  1.6471E+00  1.0000E-02
             6.9010E+00
 PARAMETER:  2.9849E-01 -7.6666E-01  2.9597E+00  6.3396E-01  7.9405E-01  9.3373E-01  2.1272E+00  5.6560E-01  5.9904E-01 -4.6472E+00
             2.0317E+00
 GRADIENT:   3.6008E+01  8.0309E+00  5.7685E-03  4.8797E+01  4.2683E+00  6.2092E+01  2.8599E+02  1.3185E-01  9.8244E+00  0.0000E+00
             3.3672E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1632.29200487937        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2430
 NPARAMETR:  1.2239E+00  4.1994E-01  1.7977E+01  1.7252E+00  2.0057E+00  2.3083E+00  7.7821E+00  1.4963E+00  1.6762E+00  1.0000E-02
             6.9135E+00
 PARAMETER:  3.0202E-01 -7.6765E-01  2.9891E+00  6.4533E-01  7.9601E-01  9.3653E-01  2.1518E+00  5.0302E-01  6.1656E-01 -4.6472E+00
             2.0335E+00
 GRADIENT:   1.6778E+00  8.5782E-01 -3.3518E-01  2.3003E+00 -2.0949E+00  3.9475E+00  2.1001E+01 -8.7106E-03  4.7313E-01  0.0000E+00
            -3.4299E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1632.41600396384        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2613
 NPARAMETR:  1.2176E+00  3.7603E-01  1.8928E+01  1.7308E+00  2.0108E+00  2.2984E+00  8.0466E+00  1.5590E+00  1.6707E+00  1.0000E-02
             6.9071E+00
 PARAMETER:  2.9689E-01 -8.7808E-01  3.0407E+00  6.4857E-01  7.9855E-01  9.3223E-01  2.1852E+00  5.4404E-01  6.1324E-01 -4.6472E+00
             2.0325E+00
 GRADIENT:  -1.4669E-01 -4.6236E-02  6.7926E-02 -3.3797E+00 -1.8982E+00  2.4709E+00  2.4980E+01 -4.7261E-02 -1.9320E-01  0.0000E+00
            -7.0002E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1632.44221237971        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2796             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2198E+00  3.7011E-01  1.9397E+01  1.7427E+00  2.0162E+00  2.3010E+00  8.0972E+00  1.5906E+00  1.6775E+00  1.0000E-02
             6.9113E+00
 PARAMETER:  2.9867E-01 -8.9396E-01  3.0651E+00  6.5543E-01  8.0121E-01  9.3334E-01  2.1915E+00  5.6409E-01  6.1732E-01 -4.6472E+00
             2.0332E+00
 GRADIENT:   3.6042E+01  9.2972E+00  6.0721E-01  4.9181E+01  3.1891E+00  6.1769E+01  3.1542E+02 -5.2308E-02  1.2876E+01  0.0000E+00
             3.6321E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1632.45079141780        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2975             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2196E+00  3.6585E-01  1.9564E+01  1.7474E+00  2.0235E+00  2.3007E+00  8.1253E+00  1.8018E+00  1.6813E+00  1.0000E-02
             6.9136E+00
 PARAMETER:  2.9849E-01 -9.0553E-01  3.0737E+00  6.5813E-01  8.0482E-01  9.3323E-01  2.1950E+00  6.8881E-01  6.1959E-01 -4.6472E+00
             2.0335E+00
 GRADIENT:   3.5791E+01  9.2127E+00 -5.8540E-02  4.9649E+01  6.1545E+00  6.1632E+01  3.1718E+02  1.2337E-01  1.3223E+01  0.0000E+00
             3.7367E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1632.45947862848        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3151
 NPARAMETR:  1.2214E+00  3.6879E-01  2.0023E+01  1.7548E+00  2.0241E+00  2.3034E+00  8.1412E+00  1.6985E+00  1.6831E+00  1.0000E-02
             6.9141E+00
 PARAMETER:  2.9998E-01 -8.9754E-01  3.0969E+00  6.6236E-01  8.0514E-01  9.3440E-01  2.1969E+00  6.2972E-01  6.2063E-01 -4.6472E+00
             2.0336E+00
 GRADIENT:   9.2605E-01  3.8152E-01  2.0989E-02 -1.9101E-01 -1.9239E+00  3.2916E+00  2.6080E+01 -5.8235E-02  6.5315E-02  0.0000E+00
            -2.1984E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1632.46961199469        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3334
 NPARAMETR:  1.2207E+00  3.6370E-01  2.0293E+01  1.7576E+00  2.0274E+00  2.3021E+00  8.1767E+00  1.7504E+00  1.6844E+00  1.0000E-02
             6.9136E+00
 PARAMETER:  2.9940E-01 -9.1142E-01  3.1103E+00  6.6393E-01  8.0678E-01  9.3381E-01  2.2013E+00  6.5983E-01  6.2143E-01 -4.6472E+00
             2.0335E+00
 GRADIENT:   7.0416E-01  3.0055E-01 -2.7587E-02 -5.7285E-01 -1.4348E+00  3.0924E+00  2.6589E+01 -4.6186E-02  8.5714E-02  0.0000E+00
            -2.6324E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1632.47513485223        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3517
 NPARAMETR:  1.2199E+00  3.5882E-01  2.0731E+01  1.7598E+00  2.0304E+00  2.3005E+00  8.2079E+00  1.7633E+00  1.6860E+00  1.0000E-02
             6.9131E+00
 PARAMETER:  2.9878E-01 -9.2494E-01  3.1316E+00  6.6521E-01  8.0826E-01  9.3314E-01  2.2051E+00  6.6717E-01  6.2237E-01 -4.6472E+00
             2.0334E+00
 GRADIENT:   4.8932E-01  2.4086E-01  1.8075E-01 -1.2071E+00 -2.0128E+00  2.8699E+00  2.7016E+01 -8.5715E-02  6.9455E-02  0.0000E+00
            -4.1430E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1632.48150142380        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3701
 NPARAMETR:  1.2199E+00  3.5571E-01  2.0921E+01  1.7631E+00  2.0334E+00  2.3004E+00  8.2247E+00  1.8211E+00  1.6875E+00  1.0000E-02
             6.9132E+00
 PARAMETER:  2.9875E-01 -9.3363E-01  3.1407E+00  6.6705E-01  8.0972E-01  9.3310E-01  2.2071E+00  6.9943E-01  6.2322E-01 -4.6472E+00
             2.0334E+00
 GRADIENT:   4.6173E-01  1.7892E-01  5.7152E-02 -1.0271E+00 -1.4160E+00  2.8653E+00  2.7136E+01 -5.7709E-02  9.2513E-02  0.0000E+00
            -4.3355E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1632.48576377642        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3889             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2198E+00  3.5114E-01  2.1093E+01  1.7667E+00  2.0383E+00  2.3003E+00  8.2414E+00  1.9533E+00  1.6898E+00  1.0000E-02
             6.9140E+00
 PARAMETER:  2.9868E-01 -9.4658E-01  3.1489E+00  6.6912E-01  8.1212E-01  9.3306E-01  2.2092E+00  7.6954E-01  6.2460E-01 -4.6472E+00
             2.0335E+00
 GRADIENT:   3.5861E+01  9.4025E+00  1.7345E-01  5.2062E+01  5.6120E+00  6.1551E+01  3.2306E+02  8.9069E-02  1.3472E+01  0.0000E+00
             3.6586E+01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1632.48944636414        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4070             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2199E+00  3.5024E-01  2.1360E+01  1.7684E+00  2.0399E+00  2.3003E+00  8.2496E+00  1.9494E+00  1.6911E+00  1.0000E-02
             6.9143E+00
 PARAMETER:  2.9875E-01 -9.4915E-01  3.1615E+00  6.7006E-01  8.1291E-01  9.3303E-01  2.2102E+00  7.6752E-01  6.2535E-01 -4.6472E+00
             2.0336E+00
 GRADIENT:   3.5914E+01  9.4316E+00  3.2916E-01  5.2194E+01  5.0311E+00  6.1546E+01  3.2330E+02  4.1722E-02  1.3526E+01  0.0000E+00
             3.6535E+01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1632.49078673551        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4251
 NPARAMETR:  1.2198E+00  3.4838E-01  2.1529E+01  1.7707E+00  2.0427E+00  2.3002E+00  8.2615E+00  1.9967E+00  1.6923E+00  1.0000E-02
             6.9145E+00
 PARAMETER:  2.9871E-01 -9.5447E-01  3.1694E+00  6.7140E-01  8.1425E-01  9.3300E-01  2.2116E+00  7.9148E-01  6.2610E-01 -4.6472E+00
             2.0336E+00
 GRADIENT:   3.9073E-01  1.1854E-02 -2.8957E-01 -8.6752E-01  3.3653E-01  2.8535E+00  2.7403E+01  5.8130E-02  2.1339E-01  0.0000E+00
            -3.9916E-02

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1632.49350161113        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     4441             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2199E+00  3.4662E-01  2.2030E+01  1.7740E+00  2.0443E+00  2.3001E+00  8.2797E+00  1.9559E+00  1.6939E+00  1.0000E-02
             6.9144E+00
 PARAMETER:  2.9880E-01 -9.5953E-01  3.1924E+00  6.7324E-01  8.1505E-01  9.3297E-01  2.2138E+00  7.7083E-01  6.2704E-01 -4.6472E+00
             2.0336E+00
 GRADIENT:   3.5966E+01  9.5098E+00  6.1195E-01  5.2840E+01  3.9287E+00  6.1543E+01  3.2456E+02 -4.3331E-02  1.3591E+01  0.0000E+00
             3.6159E+01

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1632.49534953579        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4622
 NPARAMETR:  1.2199E+00  3.4550E-01  2.2121E+01  1.7751E+00  2.0461E+00  2.3001E+00  8.2851E+00  1.9873E+00  1.6949E+00  1.0000E-02
             6.9146E+00
 PARAMETER:  2.9878E-01 -9.6277E-01  3.1965E+00  6.7388E-01  8.1594E-01  9.3295E-01  2.2145E+00  7.8677E-01  6.2762E-01 -4.6472E+00
             2.0336E+00
 GRADIENT:   4.4266E-01  7.9696E-02  3.5765E-02 -8.5515E-01 -8.5631E-01  2.8597E+00  2.7690E+01 -3.9513E-02  1.8410E-01  0.0000E+00
            -3.6114E-01

0ITERATION NO.:  145    OBJECTIVE VALUE:  -1632.49745647559        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4810             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2199E+00  3.4316E-01  2.2360E+01  1.7780E+00  2.0498E+00  2.3000E+00  8.2977E+00  2.0635E+00  1.6967E+00  1.0000E-02
             6.9152E+00
 PARAMETER:  2.9877E-01 -9.6955E-01  3.2073E+00  6.7549E-01  8.1775E-01  9.3291E-01  2.2160E+00  8.2438E-01  6.2868E-01 -4.6472E+00
             2.0337E+00
 GRADIENT:   3.5888E+01  9.4632E+00  3.6733E-01  5.3184E+01  5.1951E+00  6.1482E+01  3.2571E+02  4.6771E-02  1.3754E+01  0.0000E+00
             3.6444E+01

0ITERATION NO.:  150    OBJECTIVE VALUE:  -1632.49809622123        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     4993
 NPARAMETR:  1.2199E+00  3.4255E-01  2.2475E+01  1.7789E+00  2.0507E+00  2.3000E+00  8.3019E+00  2.0619E+00  1.6973E+00  1.0000E-02
             6.9153E+00
 PARAMETER:  2.9878E-01 -9.7134E-01  3.2124E+00  6.7598E-01  8.1820E-01  9.3290E-01  2.2165E+00  8.2362E-01  6.2904E-01 -4.6472E+00
             2.0337E+00
 GRADIENT:   4.1834E-01  3.8768E-02 -8.0776E-02 -7.7718E-01 -1.7937E-01  2.8536E+00  2.7849E+01  4.1092E-03  2.1890E-01  0.0000E+00
            -2.1132E-01

0ITERATION NO.:  151    OBJECTIVE VALUE:  -1632.49809622123        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     5021
 NPARAMETR:  1.2200E+00  3.4106E-01  2.2907E+01  1.7816E+00  2.0525E+00  2.2999E+00  8.3156E+00  2.0450E+00  1.6987E+00  1.0000E-02
             6.9154E+00
 PARAMETER:  2.9878E-01 -9.7134E-01  3.2124E+00  6.7598E-01  8.1820E-01  9.3290E-01  2.2165E+00  8.2362E-01  6.2904E-01 -4.6472E+00
             2.0337E+00
 GRADIENT:  -5.7902E-03  1.7195E-02 -7.6140E-02 -1.5459E-01 -1.4691E-01  2.0904E-03 -6.3653E-02  4.1311E-03 -3.4090E-02  0.0000E+00
            -5.4241E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     5021
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.3640E-03  6.2770E-02 -7.1199E-03 -8.5560E-02  6.5992E-05
 SE:             2.9249E-02  2.1032E-02  3.8822E-03  1.8167E-02  1.7634E-04
 N:                     100         100         100         100         100

 P VAL.:         8.5449E-01  2.8409E-03  6.6655E-02  2.4841E-06  7.0823E-01

 ETASHRINKSD(%)  2.0136E+00  2.9540E+01  8.6994E+01  3.9138E+01  9.9409E+01
 ETASHRINKVR(%)  3.9867E+00  5.0354E+01  9.8308E+01  6.2959E+01  9.9997E+01
 EBVSHRINKSD(%)  2.6495E+00  3.1795E+01  8.6515E+01  2.8755E+01  9.9371E+01
 EBVSHRINKVR(%)  5.2288E+00  5.3480E+01  9.8182E+01  4.9242E+01  9.9996E+01
 RELATIVEINF(%)  9.4647E+01  2.6112E+01  4.6710E-01  2.9209E+01  1.0279E-03
 EPSSHRINKSD(%)  8.5521E+00
 EPSSHRINKVR(%)  1.6373E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1632.4980962212260     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       21.591263547184781     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   191.05
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    18.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1632.498       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.22E+00  3.43E-01  2.25E+01  1.78E+00  2.05E+00  2.30E+00  8.30E+00  2.06E+00  1.70E+00  1.00E-02  6.92E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        5.62E+01
 
 TH 2
+        2.06E+01  1.12E+01
 
 TH 3
+       -1.49E-01  8.13E-02  5.55E-03
 
 TH 4
+        3.86E+00  9.78E+00  3.07E-01  1.99E+01
 
 TH 5
+        1.59E+01 -9.44E+00 -6.21E-01 -3.45E+01  6.95E+01
 
 TH 6
+        8.08E+00  1.38E+00 -8.15E-02 -3.15E+00  9.03E+00  1.86E+00
 
 TH 7
+        2.09E+00  3.14E-01 -2.27E-02 -9.17E-01  2.52E+00  5.01E-01  1.35E-01
 
 TH 8
+        3.36E-02 -8.37E-02 -3.73E-03 -2.23E-01  4.19E-01  4.74E-02  1.35E-02  2.62E-03
 
 TH 9
+        1.82E+00 -1.83E+00 -9.96E-02 -5.72E+00  1.12E+01  1.37E+00  3.84E-01  6.83E-02  1.80E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.00E+00 -1.31E+00 -1.77E-02 -1.75E+00  2.04E+00  9.29E-03  1.66E-02  1.84E-02  3.84E-01  0.00E+00  6.83E-01
 
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
+        1.32E+02
 
 TH 2
+       -6.66E-01  8.24E+01
 
 TH 3
+        3.16E-02  1.92E-01  1.94E-02
 
 TH 4
+       -1.83E+00  2.60E+01 -2.06E-01  7.75E+01
 
 TH 5
+       -2.23E+00 -1.10E+01 -1.03E+00 -1.14E+01  9.89E+01
 
 TH 6
+       -3.49E-01 -3.18E-01  9.15E-03  2.21E-01  1.06E-01  3.37E+01
 
 TH 7
+        3.72E-02  6.38E+00 -2.66E-04 -2.77E+00  1.34E+00 -3.81E-02  1.38E+00
 
 TH 8
+       -5.95E-02 -3.23E-01 -4.23E-02  1.51E-01  8.30E-01 -9.78E-03 -1.84E-03  2.92E-01
 
 TH 9
+        6.72E-01 -5.29E+00 -1.15E-01 -1.08E+01  1.53E+00 -6.34E-01 -6.31E-02 -2.04E-01  3.42E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.02E+00 -1.71E+00 -1.30E-02 -7.19E+00  7.19E-01  9.84E-01  2.17E-01  1.69E-01  1.10E+00  0.00E+00  2.30E+01
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.34E+02
 
 TH 2
+        6.31E+01  7.62E+01
 
 TH 3
+        1.53E-01  1.79E-01  1.03E-02
 
 TH 4
+        5.68E+01  2.51E+01  1.32E-01  8.30E+01
 
 TH 5
+       -3.21E+01 -2.39E+01 -8.36E-01 -3.19E+01  9.28E+01
 
 TH 6
+        3.65E+01  2.39E+01 -1.92E-02 -3.36E+00 -1.81E+00  4.08E+01
 
 TH 7
+        1.02E+00  8.08E+00 -1.55E-02 -3.55E+00  6.91E-01  2.31E+00  2.22E+00
 
 TH 8
+        9.00E-02 -3.40E-02 -8.81E-03 -2.59E-02  5.69E-01  4.42E-02  7.93E-03  2.49E-02
 
 TH 9
+       -1.69E+00 -2.81E+00 -1.49E-01 -9.91E+00  1.53E+01  4.55E+00  3.03E-01  2.21E-01  2.49E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.54E+01  1.72E+00 -9.54E-02 -6.56E+01  2.25E+01 -1.53E+01  9.49E+00 -2.64E-01 -2.29E+00  0.00E+00  1.03E+03
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,      209.582
Stop Time:
Wed Sep 29 08:31:33 CDT 2021
