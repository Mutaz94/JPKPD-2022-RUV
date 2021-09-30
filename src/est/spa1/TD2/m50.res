Thu Sep 30 02:10:13 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat50.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2109.39597324047        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2707E+02 -8.2754E+00 -2.2762E+01  1.1365E+01  3.1499E+00  3.9954E+01  9.6995E+00  1.9620E+01  1.3348E+01  2.4015E+01
            -5.8719E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2118.22587535739        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9486E-01  1.0239E+00  1.1434E+00  1.0369E+00  1.0643E+00  1.0064E+00  9.5734E-01  8.8909E-01  9.5771E-01  8.8576E-01
             1.0812E+00
 PARAMETER:  9.4847E-02  1.2364E-01  2.3398E-01  1.3624E-01  1.6233E-01  1.0634E-01  5.6398E-02 -1.7558E-02  5.6785E-02 -2.1307E-02
             1.7804E-01
 GRADIENT:  -6.1140E+00  8.7883E+00 -8.0985E-02 -4.9753E+00  1.4808E+00  8.8016E-01  2.6512E+00  9.3889E+00 -4.3814E+00 -6.8800E+00
            -2.4925E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2123.44638554289        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.9697E-01  8.0821E-01  1.3427E+00  1.2050E+00  1.0734E+00  9.8141E-01  7.9793E-01  4.8382E-01  9.6318E-01  1.1084E+00
             1.0772E+00
 PARAMETER:  9.6964E-02 -1.1293E-01  3.9468E-01  2.8651E-01  1.7082E-01  8.1239E-02 -1.2573E-01 -6.2605E-01  6.2485E-02  2.0293E-01
             1.7432E-01
 GRADIENT:   5.0987E+00  2.8056E+01 -3.8465E+00  5.1067E+01 -9.2085E+00 -7.5954E+00  6.6446E-01  2.3808E+00  5.2804E+00  7.8360E+00
            -1.1468E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2126.18049206069        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.8833E-01  6.3908E-01  1.4157E+00  1.2907E+00  1.0390E+00  9.9985E-01  8.4591E-01  2.6879E-01  8.7965E-01  1.0800E+00
             1.0818E+00
 PARAMETER:  8.8259E-02 -3.4773E-01  4.4759E-01  3.5515E-01  1.3829E-01  9.9848E-02 -6.7346E-02 -1.2138E+00 -2.8226E-02  1.7694E-01
             1.7864E-01
 GRADIENT:  -9.6758E+00  1.2112E+01  6.1212E+00  1.4495E+01 -9.7464E+00  7.6231E-01  2.5402E-01  4.4051E-01 -1.8736E+00 -4.9889E-01
             6.1012E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2127.74440917520        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.9108E-01  3.3931E-01  1.4664E+00  1.4765E+00  9.6742E-01  9.9138E-01  6.0693E-01  3.5985E-02  8.0033E-01  1.1016E+00
             1.0767E+00
 PARAMETER:  9.1035E-02 -9.8085E-01  4.8278E-01  4.8965E-01  6.6878E-02  9.1347E-02 -3.9933E-01 -3.2247E+00 -1.2274E-01  1.9675E-01
             1.7387E-01
 GRADIENT:   6.2628E+00  5.2685E+00 -1.8051E-01  2.0675E+01 -4.6416E+00 -6.6469E-01 -3.1603E-02  5.8003E-03 -6.5369E-01  1.5937E+00
            -1.1381E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2127.91241246902        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.8438E-01  2.3830E-01  1.5453E+00  1.5385E+00  9.6603E-01  9.9192E-01  4.4250E-01  1.0007E-02  7.6783E-01  1.1117E+00
             1.0787E+00
 PARAMETER:  8.4255E-02 -1.3342E+00  5.3523E-01  5.3081E-01  6.5444E-02  9.1886E-02 -7.1531E-01 -4.5045E+00 -1.6419E-01  2.0586E-01
             1.7575E-01
 GRADIENT:  -5.3843E+00  3.0630E+00  4.6715E+00  1.1059E+01 -6.5007E+00  1.5079E-01  1.9066E-03  4.1699E-04 -1.7120E+00 -4.4390E-01
            -2.8637E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2128.02321141375        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  9.8636E-01  2.0431E-01  1.5480E+00  1.5512E+00  9.6371E-01  9.9133E-01  2.4553E-01  1.0000E-02  7.6237E-01  1.1173E+00
             1.0785E+00
 PARAMETER:  8.6268E-02 -1.4881E+00  5.3699E-01  5.3904E-01  6.3040E-02  9.1294E-02 -1.3043E+00 -4.9652E+00 -1.7133E-01  2.1089E-01
             1.7560E-01
 GRADIENT:   6.7683E-01  2.7254E-01 -1.1727E-02 -8.7639E+00  2.1562E+00  2.1089E-01  5.1700E-03  0.0000E+00  3.4918E-01 -7.9581E-02
             5.8950E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2128.02994726430        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  9.8643E-01  1.9770E-01  1.5467E+00  1.5544E+00  9.6039E-01  9.9091E-01  7.1891E-02  1.0000E-02  7.6010E-01  1.1170E+00
             1.0778E+00
 PARAMETER:  8.6341E-02 -1.5210E+00  5.3611E-01  5.4106E-01  5.9580E-02  9.0868E-02 -2.5326E+00 -4.9652E+00 -1.7430E-01  2.1063E-01
             1.7492E-01
 GRADIENT:   1.1159E+00  2.1033E-01  6.2455E-01 -1.0620E+01  8.8626E-01  9.3363E-02  6.9779E-04  0.0000E+00  6.9632E-02  1.1164E-01
             4.0957E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2128.03177299813        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1438
 NPARAMETR:  9.8644E-01  1.9553E-01  1.5445E+00  1.5554E+00  9.5919E-01  9.9087E-01  2.9835E-02  1.0000E-02  7.5966E-01  1.1160E+00
             1.0778E+00
 PARAMETER:  8.6344E-02 -1.5321E+00  5.3473E-01  5.4175E-01  5.8335E-02  9.0831E-02 -3.4121E+00 -4.9652E+00 -1.7488E-01  2.0977E-01
             1.7488E-01
 GRADIENT:   1.2214E+00  1.5129E-01  2.9848E-01 -1.0919E+01  1.2644E+00  9.5510E-02  1.7542E-04  0.0000E+00  1.5370E-01  7.7184E-02
             3.6408E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2128.03831352014        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:     1542
 NPARAMETR:  9.8706E-01  1.9530E-01  1.5224E+00  1.5554E+00  9.4794E-01  9.9129E-01  1.0000E-02  1.0000E-02  7.5873E-01  1.1067E+00
             1.0777E+00
 PARAMETER:  8.6972E-02 -1.5332E+00  5.2026E-01  5.4171E-01  4.6541E-02  9.1252E-02 -2.6377E+01 -4.9652E+00 -1.7611E-01  2.0136E-01
             1.7479E-01
 GRADIENT:   3.6215E+02  2.0983E+01  8.2900E+00  8.2530E+02  3.1874E+00  3.5092E+01  0.0000E+00  0.0000E+00  1.4674E+01  1.6793E+00
             1.5777E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2128.03947170863        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     1657
 NPARAMETR:  9.8654E-01  1.9423E-01  1.5197E+00  1.5560E+00  9.4693E-01  9.9099E-01  1.0000E-02  1.0000E-02  7.5894E-01  1.1053E+00
             1.0777E+00
 PARAMETER:  8.6445E-02 -1.5387E+00  5.1852E-01  5.4211E-01  4.5466E-02  9.0951E-02 -3.0798E+01 -4.9652E+00 -1.7583E-01  2.0011E-01
             1.7482E-01
 GRADIENT:   1.5328E+00  5.7161E-01  7.4179E-01 -7.4052E+00 -1.7152E+00  1.6559E-01  0.0000E+00  0.0000E+00 -1.7273E-01  1.9192E-01
            -1.2556E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2128.04135227688        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1847
 NPARAMETR:  9.8633E-01  1.9156E-01  1.5208E+00  1.5565E+00  9.4691E-01  9.9072E-01  1.0000E-02  1.0000E-02  7.5841E-01  1.1053E+00
             1.0777E+00
 PARAMETER:  8.6238E-02 -1.5526E+00  5.1921E-01  5.4246E-01  4.5447E-02  9.0678E-02 -3.2252E+01 -4.9652E+00 -1.7653E-01  2.0010E-01
             1.7483E-01
 GRADIENT:   1.1919E+00  3.1209E-01  6.3451E-01 -1.0211E+01 -1.1534E+00  7.9885E-02  0.0000E+00  0.0000E+00 -2.2098E-02  9.9894E-02
            -8.1350E-02

0ITERATION NO.:   56    OBJECTIVE VALUE:  -2128.04135227688        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1877
 NPARAMETR:  9.8637E-01  1.9093E-01  1.5193E+00  1.5562E+00  9.4786E-01  9.9075E-01  1.0000E-02  1.0000E-02  7.5854E-01  1.1048E+00
             1.0778E+00
 PARAMETER:  8.6238E-02 -1.5526E+00  5.1921E-01  5.4246E-01  4.5447E-02  9.0678E-02 -3.2252E+01 -4.9652E+00 -1.7653E-01  2.0010E-01
             1.7483E-01
 GRADIENT:  -1.0263E-01  1.1238E-01  6.5812E-01  9.3994E-01 -8.0592E-01 -1.4969E-02  0.0000E+00  0.0000E+00 -7.3849E-02  9.5842E-02
            -8.0105E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1877
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.9668E-04 -6.9329E-05 -1.7275E-04 -4.0904E-03 -2.4785E-02
 SE:             2.9853E-02  3.5395E-05  1.4941E-04  2.9443E-02  2.4730E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7871E-01  5.0146E-02  2.4757E-01  8.8951E-01  3.1622E-01

 ETASHRINKSD(%)  1.0000E-10  9.9881E+01  9.9499E+01  1.3622E+00  1.7153E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  9.9997E+01  2.7059E+00  3.1363E+01
 EBVSHRINKSD(%)  3.8563E-01  9.9889E+01  9.9489E+01  1.6088E+00  1.3990E+01
 EBVSHRINKVR(%)  7.6977E-01  1.0000E+02  9.9997E+01  3.1918E+00  2.6023E+01
 RELATIVEINF(%)  9.7630E+01  7.2501E-06  4.2048E-04  6.6190E+00  1.0133E+01
 EPSSHRINKSD(%)  3.1048E+01
 EPSSHRINKVR(%)  5.2456E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2128.0413522768836     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1209.1028190722109     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.19
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2128.041       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.92E-01  1.52E+00  1.56E+00  9.47E-01  9.91E-01  1.00E-02  1.00E-02  7.58E-01  1.11E+00  1.08E+00
 


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
+        1.15E+03
 
 TH 2
+       -2.76E+01  3.91E+02
 
 TH 3
+       -5.15E+00  5.55E+01  1.29E+02
 
 TH 4
+       -6.57E+00  4.71E+02 -4.57E+01  7.72E+02
 
 TH 5
+        7.54E+00 -2.27E+02 -2.74E+02 -1.88E+01  7.59E+02
 
 TH 6
+        9.78E-01 -4.76E+00 -7.45E-01 -2.12E+00  1.71E-01  2.00E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.32E+00 -9.93E+01  5.19E+00 -1.43E+00 -9.74E-01 -3.53E-01  0.00E+00  0.00E+00  3.26E+02
 
 TH10
+        2.10E+00  1.25E+01 -8.68E+00  2.24E+00 -6.03E+01  9.90E-01  0.00E+00  0.00E+00 -8.83E-01  8.30E+01
 
 TH11
+       -9.04E+00 -1.78E+01 -1.80E+01 -1.29E+01  1.62E+01  1.09E+00  0.00E+00  0.00E+00  7.61E+00  2.95E+01  3.71E+02
 
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
 #CPUT: Total CPU Time in Seconds,       34.863
Stop Time:
Thu Sep 30 02:10:50 CDT 2021
