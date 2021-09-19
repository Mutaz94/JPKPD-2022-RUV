Sat Sep 18 14:37:36 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat35.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1704.27221509408        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.8987E+01 -2.7577E+01 -4.3636E+01  1.0944E+01  6.6397E+01  1.1623E+01 -2.1221E+00  1.1498E+01  1.4190E+00 -1.6256E+01
             1.1635E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1712.34359906033        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:       91
 NPARAMETR:  1.0562E+00  9.9508E-01  1.1557E+00  9.8890E-01  1.0217E+00  9.4159E-01  1.0048E+00  8.5050E-01  1.0054E+00  1.1650E+00
             9.6540E-01
 PARAMETER:  1.5472E-01  9.5073E-02  2.4475E-01  8.8841E-02  1.2145E-01  3.9820E-02  1.0481E-01 -6.1936E-02  1.0539E-01  2.5270E-01
             6.4786E-02
 GRADIENT:   7.6683E+01 -1.4703E+01  4.5815E+00 -3.5674E+01 -1.2089E+01 -4.9635E+00  2.1636E+00  2.7447E+00 -1.7087E+00  2.3966E-01
            -1.5803E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1713.25558694834        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      189
 NPARAMETR:  1.0548E+00  9.5555E-01  1.2255E+00  1.0523E+00  1.0475E+00  9.6449E-01  9.2064E-01  7.2753E-01  1.0167E+00  1.2439E+00
             9.7088E-01
 PARAMETER:  1.5332E-01  5.4536E-02  3.0335E-01  1.5101E-01  1.4636E-01  6.3846E-02  1.7316E-02 -2.1810E-01  1.1654E-01  3.1828E-01
             7.0452E-02
 GRADIENT:  -4.5753E+00  1.0608E+01  1.9188E+00  1.4830E+01  3.2424E-01 -5.2087E-01  5.9244E-02 -2.2600E+00  4.6622E-01  3.4515E+00
            -5.5504E-02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1713.74794110818        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  1.0564E+00  7.9394E-01  1.4033E+00  1.1535E+00  1.0343E+00  9.5952E-01  9.8321E-01  9.8950E-01  9.4617E-01  1.2360E+00
             9.6738E-01
 PARAMETER:  1.5485E-01 -1.3074E-01  4.3882E-01  2.4277E-01  1.3377E-01  5.8679E-02  8.3064E-02  8.9443E-02  4.4667E-02  3.1190E-01
             6.6836E-02
 GRADIENT:   3.7478E+00  7.2993E+00  3.1583E+00  5.0051E+00 -1.1506E+01 -1.7657E+00  6.5608E-01  1.5511E+00 -1.0866E+00  1.7733E+00
            -1.0305E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1714.24499778087        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      542
 NPARAMETR:  1.0519E+00  6.0796E-01  1.5292E+00  1.2753E+00  1.0223E+00  9.6298E-01  8.4938E-01  1.0267E+00  9.0198E-01  1.2403E+00
             9.7007E-01
 PARAMETER:  1.5056E-01 -3.9764E-01  5.2473E-01  3.4315E-01  1.2201E-01  6.2278E-02 -6.3244E-02  1.2637E-01 -3.1662E-03  3.1533E-01
             6.9609E-02
 GRADIENT:  -1.5159E+00  6.0524E+00  9.0860E-01  1.4500E+01 -1.5691E+00  5.8398E-01 -1.5578E-01 -1.7913E-01 -1.0737E+00 -6.8125E-01
            -3.7091E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1714.58377270010        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  1.0494E+00  4.0539E-01  1.6424E+00  1.3967E+00  9.9634E-01  9.6074E-01  5.5578E-01  1.1084E+00  8.4571E-01  1.2370E+00
             9.7205E-01
 PARAMETER:  1.4819E-01 -8.0292E-01  5.9615E-01  4.3413E-01  9.6333E-02  5.9947E-02 -4.8739E-01  2.0290E-01 -6.7581E-02  3.1273E-01
             7.1648E-02
 GRADIENT:  -8.3792E-01  1.2037E+00  2.7919E-01  3.0615E+00  1.0045E+00  6.9644E-01 -5.0919E-02 -3.2845E-01 -5.8044E-01 -7.4879E-01
             4.7060E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1714.82125934782        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  1.0456E+00  1.8371E-01  1.7296E+00  1.5414E+00  9.5532E-01  9.5777E-01  1.7896E-01  1.1985E+00  7.6792E-01  1.2208E+00
             9.7217E-01
 PARAMETER:  1.4456E-01 -1.5944E+00  6.4792E-01  5.3267E-01  5.4294E-02  5.6854E-02 -1.6206E+00  2.8107E-01 -1.6406E-01  2.9950E-01
             7.1774E-02
 GRADIENT:  -1.9207E+00  1.3353E+00 -3.5603E-01  1.3336E+01 -2.0117E+00  4.9595E-01  7.9136E-04  1.0662E-01 -1.2147E+00  5.0593E-01
             5.9595E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1714.90369630562        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  1.0456E+00  1.2057E-01  1.8255E+00  1.5794E+00  9.6428E-01  9.5580E-01  8.3102E-02  1.2825E+00  7.4813E-01  1.2270E+00
             9.7018E-01
 PARAMETER:  1.4458E-01 -2.0156E+00  7.0187E-01  5.5706E-01  6.3627E-02  5.4792E-02 -2.3877E+00  3.4884E-01 -1.9018E-01  3.0455E-01
             6.9729E-02
 GRADIENT:   4.7837E-01  1.7622E-01  2.0865E-01  6.2044E-01 -6.3144E-01 -2.6700E-02  2.4533E-04  2.6183E-01 -2.0037E-01  7.0983E-02
            -2.8885E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1714.94735177121        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  1.0442E+00  5.4766E-02  1.8385E+00  1.6221E+00  9.4996E-01  9.5471E-01  1.8271E-02  1.2940E+00  7.2714E-01  1.2186E+00
             9.7015E-01
 PARAMETER:  1.4325E-01 -2.8047E+00  7.0898E-01  5.8369E-01  4.8668E-02  5.3651E-02 -3.9024E+00  3.5773E-01 -2.1863E-01  2.9774E-01
             6.9696E-02
 GRADIENT:  -2.0813E-01  1.5746E-01 -1.8618E-02  4.5754E+00 -1.1686E+00 -1.6831E-01  2.9738E-06  1.0019E-01 -3.8608E-01  2.1265E-01
            -8.5186E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1714.97170078875        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1419
 NPARAMETR:  1.0438E+00  2.0178E-02  1.9008E+00  1.6447E+00  9.5767E-01  9.5492E-01  1.0000E-02  1.3413E+00  7.1655E-01  1.2232E+00
             9.7012E-01
 PARAMETER:  1.4283E-01 -3.8031E+00  7.4226E-01  5.9756E-01  5.6747E-02  5.3869E-02 -5.9464E+00  3.9367E-01 -2.3331E-01  3.0148E-01
             6.9663E-02
 GRADIENT:  -5.7799E-02  1.9910E-02  9.3949E-02  5.2207E-01  3.3507E-01  5.9695E-02  0.0000E+00 -6.8971E-02 -1.4465E-01 -2.2106E-01
            -1.4138E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1714.97870761247        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1602             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0436E+00  1.0000E-02  1.8928E+00  1.6502E+00  9.5324E-01  9.5461E-01  1.0000E-02  1.3370E+00  7.1399E-01  1.2207E+00
             9.7010E-01
 PARAMETER:  1.4268E-01 -4.5107E+00  7.3807E-01  6.0091E-01  5.2116E-02  5.3550E-02 -7.4065E+00  3.9042E-01 -2.3688E-01  2.9942E-01
             6.9648E-02
 GRADIENT:   7.0456E+01  0.0000E+00  8.4363E-01  1.4101E+02  1.0901E+00  5.1402E+00  0.0000E+00  1.0053E-01  2.2686E+00  2.8050E-01
             8.6628E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1714.97878123874        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1695
 NPARAMETR:  1.0436E+00  1.0000E-02  1.8932E+00  1.6503E+00  9.5298E-01  9.5462E-01  1.0000E-02  1.3370E+00  7.1397E-01  1.2211E+00
             9.7010E-01
 PARAMETER:  1.4268E-01 -4.5107E+00  7.3826E-01  6.0098E-01  5.1843E-02  5.3556E-02 -7.4065E+00  3.9040E-01 -2.3691E-01  2.9977E-01
             6.9643E-02
 GRADIENT:   5.0241E-02  0.0000E+00  1.3383E-02 -8.0229E-02 -1.4555E-02 -4.5997E-03  0.0000E+00 -3.3957E-03  4.3941E-03  2.1426E-03
             4.9495E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1695
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.2578E-04 -3.7351E-06 -3.0331E-02 -7.6499E-03 -3.9369E-02
 SE:             2.9863E-02  1.7830E-06  1.5536E-02  2.9216E-02  2.2245E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9397E-01  3.6190E-02  5.0902E-02  7.9345E-01  7.6757E-02

 ETASHRINKSD(%)  1.0000E-10  9.9994E+01  4.7953E+01  2.1221E+00  2.5478E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  7.2911E+01  4.1992E+00  4.4464E+01
 EBVSHRINKSD(%)  4.1130E-01  9.9994E+01  5.2152E+01  2.3927E+00  2.0774E+01
 EBVSHRINKVR(%)  8.2092E-01  1.0000E+02  7.7106E+01  4.7281E+00  3.7232E+01
 RELATIVEINF(%)  9.7218E+01  1.6973E-08  5.9257E+00  5.7605E+00  1.2220E+01
 EPSSHRINKSD(%)  4.4562E+01
 EPSSHRINKVR(%)  6.9266E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.9787812387351     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -979.82795467499693     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.32
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.979       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  1.89E+00  1.65E+00  9.53E-01  9.55E-01  1.00E-02  1.34E+00  7.14E-01  1.22E+00  9.70E-01
 


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
+        1.11E+03
 
 TH 2
+        0.00E+00  1.51E+03
 
 TH 3
+       -5.13E-01  0.00E+00  4.38E+01
 
 TH 4
+       -1.22E+01  0.00E+00 -2.32E+01  7.64E+02
 
 TH 5
+       -1.35E+01  0.00E+00 -1.14E+02 -4.59E+01  5.00E+02
 
 TH 6
+       -9.99E+00  0.00E+00  5.91E-01 -1.70E+00 -2.47E+01  2.22E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -4.62E-01  0.00E+00 -1.63E+01 -2.86E+00  2.67E-02 -2.50E+00  0.00E+00  2.05E+01
 
 TH 9
+        5.42E+00  0.00E+00  4.29E+00 -3.50E+00 -3.50E+01  1.37E+01  0.00E+00 -1.94E+00  3.51E+02
 
 TH10
+       -2.24E+00  0.00E+00 -2.52E+00 -1.59E+00 -6.89E+01  9.07E-01  0.00E+00  9.26E+00  1.23E+00  5.83E+01
 
 TH11
+       -1.92E+01  0.00E+00 -3.94E+00 -1.12E+01 -2.72E+01  2.77E+00  0.00E+00 -1.12E+00 -4.83E+00  1.42E+01  2.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.265
Stop Time:
Sat Sep 18 14:38:03 CDT 2021
