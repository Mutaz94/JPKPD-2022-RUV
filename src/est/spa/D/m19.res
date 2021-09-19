Sat Sep 18 15:10:55 CDT 2021
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
$DATA ../../../../data/spa/D/dat19.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   959.211738184041        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.7527E+02 -1.4319E+02 -8.4250E+01 -3.0000E+02  2.0941E+02 -7.2791E+02 -3.6866E+02 -3.0325E+01 -7.4154E+02 -2.1979E+02
            -3.3458E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1031.26349767064        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4520E+00  1.2800E+00  5.5119E+00  2.0168E+00  2.2886E+00  2.4210E+00  3.2311E+00  1.0240E+00  5.6998E+00  1.5194E+00
             3.9547E+00
 PARAMETER:  4.7295E-01  3.4690E-01  1.8069E+00  8.0152E-01  9.2793E-01  9.8419E-01  1.2728E+00  1.2368E-01  1.8404E+00  5.1834E-01
             1.4749E+00
 GRADIENT:   8.2594E+01 -2.2483E+01 -5.2249E+00  4.4306E+01  9.1995E+00  4.9778E+01  3.1210E+01  2.3678E-01  1.1146E+02  6.1967E+00
             6.1442E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1062.00785608511        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.1694E+00  2.1809E+00  7.9387E+00  7.5130E-01  1.3147E+01  2.2287E+00  2.7564E+00  1.1702E+00  4.4379E+00  1.3986E+01
             4.9674E+00
 PARAMETER:  2.5649E-01  8.7973E-01  2.1717E+00 -1.8595E-01  2.6762E+00  9.0142E-01  1.1139E+00  2.5716E-01  1.5902E+00  2.7380E+00
             1.7029E+00
 GRADIENT:  -1.1369E+01  8.8471E+00  1.0170E+00  6.8496E+00 -4.0068E-01  3.7428E+01  5.0630E+01 -1.5908E-02  3.6467E+01  2.5733E+00
             1.8203E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1084.77010994614        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1455E+00  2.2572E+00  6.5107E+00  5.8872E-01  6.4734E+00  2.0879E+00  1.5607E+00  5.4256E+00  4.9940E+00  8.3872E+00
             4.0404E+00
 PARAMETER:  2.3580E-01  9.1411E-01  1.9734E+00 -4.2981E-01  1.9677E+00  8.3618E-01  5.4512E-01  1.7911E+00  1.7082E+00  2.2267E+00
             1.4963E+00
 GRADIENT:  -1.5995E+01  2.2689E+01 -1.3462E+00  1.3602E+01 -2.5233E+00  1.7902E+01  2.3047E+01  6.4748E+00  2.7550E+01  3.0768E+01
             5.6043E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1110.58816696891        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.1654E+00  1.3817E+00  5.5614E+00  1.0447E+00  2.8205E+00  1.9935E+00  9.9313E-01  2.2487E+00  3.2183E+00  2.6606E+00
             3.3277E+00
 PARAMETER:  2.5308E-01  4.2332E-01  1.8159E+00  1.4370E-01  1.1369E+00  7.8988E-01  9.3104E-02  9.1036E-01  1.2689E+00  1.0786E+00
             1.3023E+00
 GRADIENT:  -7.5157E-01  6.5504E+00 -3.1773E-01  3.5887E+00 -7.6235E+00  1.9181E+00  5.0292E+00  8.2293E-01  2.6549E+00  2.9435E+01
            -7.8294E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1132.66077282221        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.1697E+00  1.0163E+00  6.9618E+00  1.2476E+00  1.6030E+00  2.0451E+00  3.1885E-01  6.3042E-01  3.0249E+00  8.6878E-01
             3.5465E+00
 PARAMETER:  2.5673E-01  1.1617E-01  2.0404E+00  3.2125E-01  5.7190E-01  8.1547E-01 -1.0430E+00 -3.6137E-01  1.2069E+00 -4.0663E-02
             1.3660E+00
 GRADIENT:  -1.3275E+00 -1.2899E+01 -2.5022E+00  2.6553E+00 -3.1966E+00  8.1893E+00  9.2541E-01  1.9954E-01  1.5232E+01  6.8266E+00
             1.1273E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1136.77423801266        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.1720E+00  1.0787E+00  1.7012E+01  1.1794E+00  1.6041E+00  2.0207E+00  1.5177E-01  2.1129E-01  3.0530E+00  3.0381E-01
             3.5012E+00
 PARAMETER:  2.5867E-01  1.7577E-01  2.9339E+00  2.6501E-01  5.7257E-01  8.0346E-01 -1.7854E+00 -1.4545E+00  1.2161E+00 -1.0914E+00
             1.3531E+00
 GRADIENT:  -5.3922E-01  5.0555E-02  6.1589E-01  8.8190E-01 -1.6711E+00  2.7623E+00  4.6378E-02  3.6291E-03  1.0932E+00  9.0016E-01
             4.3473E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1137.21379781459        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  1.1960E+00  1.1151E+00  2.1503E+01  1.1496E+00  1.6175E+00  2.0699E+00  1.1719E-01  1.4487E-01  3.1794E+00  2.1633E-01
             3.4844E+00
 PARAMETER:  2.7895E-01  2.0895E-01  3.1682E+00  2.3940E-01  5.8085E-01  8.2749E-01 -2.0440E+00 -1.8319E+00  1.2567E+00 -1.4309E+00
             1.3483E+00
 GRADIENT:  -3.0060E-01  5.4912E-01  5.6694E-01 -8.1651E-02 -8.7738E-01 -4.4285E-01 -3.5515E-02  1.0034E-03 -8.1815E-01  4.5052E-01
            -7.4839E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1137.61726683132        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      803
 NPARAMETR:  1.1962E+00  1.1013E+00  1.1661E+01  1.1558E+00  1.5950E+00  2.0695E+00  5.2739E-02  3.4456E-02  3.1615E+00  5.8079E-02
             3.4920E+00
 PARAMETER:  2.7916E-01  1.9651E-01  2.5562E+00  2.4478E-01  5.6689E-01  8.2733E-01 -2.8424E+00 -3.2681E+00  1.2510E+00 -2.7459E+00
             1.3505E+00
 GRADIENT:   4.1389E-01 -5.8024E-02 -1.1902E-01 -1.6787E-01  8.3853E-02  9.7358E-02 -7.6732E-03  2.2168E-04 -4.4014E-02  3.0555E-02
            -5.0854E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1137.63160611497        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  1.1955E+00  1.0989E+00  1.2202E+01  1.1590E+00  1.5963E+00  2.0693E+00  1.5907E-02  1.0000E-02  3.1580E+00  1.0000E-02
             3.4936E+00
 PARAMETER:  2.7856E-01  1.9435E-01  2.6016E+00  2.4758E-01  5.6768E-01  8.2719E-01 -4.0410E+00 -5.2153E+00  1.2499E+00 -4.5490E+00
             1.3509E+00
 GRADIENT:   4.3242E-02 -8.0402E-02  2.0676E-03 -3.1331E-02 -1.5025E-02 -7.9146E-03 -7.1670E-04  0.0000E+00  2.3816E-02  0.0000E+00
             4.1883E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1137.63164429174        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1156
 NPARAMETR:  1.1954E+00  1.0993E+00  1.2198E+01  1.1589E+00  1.5965E+00  2.0696E+00  1.6426E-02  1.0000E-02  3.1584E+00  1.0000E-02
             3.4934E+00
 PARAMETER:  2.7851E-01  1.9466E-01  2.6013E+00  2.4747E-01  5.6779E-01  8.2733E-01 -4.0089E+00 -5.1648E+00  1.2501E+00 -4.5094E+00
             1.3509E+00
 GRADIENT:   1.1594E-02  6.4613E-03 -6.6925E-05 -1.4141E-03 -1.3398E-03  4.6994E-02 -7.8033E-04  0.0000E+00 -1.2120E-02  3.2020E-05
            -4.5272E-04

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1137.63164429174        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  1.1954E+00  1.0993E+00  1.2198E+01  1.1589E+00  1.5965E+00  2.0696E+00  1.6426E-02  1.0000E-02  3.1584E+00  1.0000E-02
             3.4934E+00
 PARAMETER:  2.7851E-01  1.9466E-01  2.6013E+00  2.4747E-01  5.6779E-01  8.2733E-01 -4.0089E+00 -5.1648E+00  1.2501E+00 -4.5094E+00
             1.3509E+00
 GRADIENT:   1.1594E-02  6.4613E-03 -6.6925E-05 -1.4141E-03 -1.3398E-03  4.6994E-02 -7.8033E-04  0.0000E+00 -1.2120E-02  3.2020E-05
            -4.5272E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1178
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0243E-03 -1.0984E-03  6.6964E-06  3.0711E-03  2.9579E-05
 SE:             2.9551E-02  2.7853E-04  5.8303E-06  2.8243E-02  8.3223E-05
 N:                     100         100         100         100         100

 P VAL.:         9.1849E-01  8.0297E-05  2.5074E-01  9.1341E-01  7.2228E-01

 ETASHRINKSD(%)  9.9912E-01  9.9067E+01  9.9980E+01  5.3838E+00  9.9721E+01
 ETASHRINKVR(%)  1.9883E+00  9.9991E+01  1.0000E+02  1.0478E+01  9.9999E+01
 EBVSHRINKSD(%)  1.0285E+00  9.9403E+01  9.9966E+01  3.8525E+00  9.9656E+01
 EBVSHRINKVR(%)  2.0464E+00  9.9996E+01  1.0000E+02  7.5565E+00  9.9999E+01
 RELATIVEINF(%)  9.7551E+01  1.5707E-03  8.9464E-06  4.1407E+01  8.5021E-04
 EPSSHRINKSD(%)  2.4240E+01
 EPSSHRINKVR(%)  4.2605E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1137.6316442917416     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -402.48081772800344     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1137.632       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.20E+00  1.10E+00  1.22E+01  1.16E+00  1.60E+00  2.07E+00  1.64E-02  1.00E-02  3.16E+00  1.00E-02  3.49E+00
 


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
+        1.77E+02
 
 TH 2
+       -4.19E+00  1.90E+02
 
 TH 3
+       -3.84E-02 -3.06E-02  8.69E-03
 
 TH 4
+       -3.33E+00  6.92E+01  5.44E-02  7.64E+01
 
 TH 5
+       -2.41E+00 -2.44E+01 -1.77E-01 -1.16E+01  5.02E+01
 
 TH 6
+        6.25E-03  1.99E+00 -2.50E-02 -7.55E-02 -1.64E-01  4.46E+01
 
 TH 7
+       -3.86E+00  3.39E-01 -1.04E-02 -3.09E-01 -6.74E-02 -1.82E-01 -1.87E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.52E-01 -3.42E+01 -1.41E-02  2.46E+00  3.45E+00 -9.34E-01  3.03E-01  0.00E+00  1.62E+01
 
 TH10
+        5.37E-01 -5.82E+00  3.30E-03 -5.77E+00  1.52E+00  7.58E-01  1.77E+01  0.00E+00  1.16E-01  2.72E+02
 
 TH11
+       -3.57E+00 -1.27E+01  5.67E-02 -4.21E+00 -9.11E-01  1.26E+00  9.82E-02  0.00E+00  2.32E+00 -3.73E-02  3.26E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.229
Stop Time:
Sat Sep 18 15:11:16 CDT 2021
