Sat Sep 25 09:57:58 CDT 2021
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
$DATA ../../../../data/spa/S1/dat55.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1669.43061491956        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2422E+01 -9.3765E+01 -1.7139E+01 -9.8101E+01  7.0161E+01  5.4942E+01 -9.2961E+00 -2.9024E+00 -1.4919E+00 -1.6634E+01
             7.3334E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1678.21578485091        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0066E+00  1.0559E+00  9.2855E-01  1.0271E+00  9.2227E-01  7.8249E-01  1.0478E+00  1.0365E+00  9.5935E-01  1.0741E+00
             9.6858E-01
 PARAMETER:  1.0654E-01  1.5437E-01  2.5866E-02  1.2676E-01  1.9086E-02 -1.4527E-01  1.4665E-01  1.3582E-01  5.8497E-02  1.7152E-01
             6.8078E-02
 GRADIENT:   6.8729E+01 -2.0806E+00  3.6461E+00 -1.2537E+01 -1.8556E+01 -3.5078E+01 -1.0454E+00  3.5076E+00 -1.8080E+00  1.0345E+01
            -4.5119E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1679.44664368444        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0002E+00  1.0187E+00  6.7753E-01  1.0413E+00  7.7962E-01  8.0677E-01  1.2051E+00  7.0116E-01  8.7587E-01  8.6218E-01
             9.7754E-01
 PARAMETER:  1.0016E-01  1.1850E-01 -2.8930E-01  1.4047E-01 -1.4894E-01 -1.1472E-01  2.8653E-01 -2.5502E-01 -3.2541E-02 -4.8293E-02
             7.7281E-02
 GRADIENT:   3.3309E+01  3.8638E+00 -1.8881E+01  2.1653E+01  1.5914E+01 -2.3070E+01  4.8794E+00  4.2898E+00 -7.0441E+00  6.2161E+00
             1.5301E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1680.87457731080        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.8879E-01  1.0348E+00  6.7310E-01  1.0222E+00  7.8089E-01  8.5551E-01  1.1448E+00  5.9035E-01  9.2973E-01  8.5032E-01
             9.7658E-01
 PARAMETER:  8.8724E-02  1.3419E-01 -2.9586E-01  1.2197E-01 -1.4732E-01 -5.6054E-02  2.3525E-01 -4.2705E-01  2.7137E-02 -6.2142E-02
             7.6297E-02
 GRADIENT:   7.5161E-02 -6.2840E-01 -4.5145E+00  4.3680E+00  3.7685E+00  5.1077E-01  2.8130E-01  1.7312E+00  2.3709E-01  1.7056E+00
            -4.7590E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1680.88002874985        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.8887E-01  1.0269E+00  6.2647E-01  1.0191E+00  7.4773E-01  8.5616E-01  1.1540E+00  4.8402E-01  9.2433E-01  8.1184E-01
             9.7656E-01
 PARAMETER:  8.8811E-02  1.2657E-01 -3.6765E-01  1.1888E-01 -1.9072E-01 -5.5299E-02  2.4324E-01 -6.2562E-01  2.1313E-02 -1.0845E-01
             7.6285E-02
 GRADIENT:  -9.3853E-01 -5.5175E-01 -3.2789E+00  2.8128E+00  2.3621E+00  5.7627E-01  3.5429E-01  1.3094E+00  2.6349E-01  1.3176E+00
             2.4258E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1680.88054279520        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.8899E-01  1.0253E+00  6.1611E-01  1.0181E+00  7.4040E-01  8.5615E-01  1.1554E+00  4.5688E-01  9.2322E-01  8.0322E-01
             9.7656E-01
 PARAMETER:  8.8925E-02  1.2503E-01 -3.8434E-01  1.1796E-01 -2.0056E-01 -5.5310E-02  2.4443E-01 -6.8332E-01  2.0116E-02 -1.1912E-01
             7.6277E-02
 GRADIENT:  -8.7379E-01 -5.1946E-01 -2.9648E+00  2.5021E+00  2.1077E+00  5.2649E-01  3.2847E-01  1.1892E+00  2.4379E-01  1.1970E+00
             2.6858E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1680.88062842029        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  9.8902E-01  1.0249E+00  6.1302E-01  1.0178E+00  7.3824E-01  8.5614E-01  1.1557E+00  4.4844E-01  9.2292E-01  8.0068E-01
             9.7655E-01
 PARAMETER:  8.8961E-02  1.2462E-01 -3.8935E-01  1.1766E-01 -2.0348E-01 -5.5318E-02  2.4473E-01 -7.0199E-01  1.9783E-02 -1.2229E-01
             7.6273E-02
 GRADIENT:  -8.4685E-01 -5.0480E-01 -2.8696E+00  2.4181E+00  2.0376E+00  5.0968E-01  3.1849E-01  1.1514E+00  2.3649E-01  1.1590E+00
             2.6448E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1681.31046515975        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  9.9835E-01  1.0124E+00  6.5888E-01  1.0343E+00  7.6123E-01  8.6277E-01  1.1826E+00  4.3724E-01  9.2372E-01  8.3674E-01
             9.7910E-01
 PARAMETER:  9.8347E-02  1.1236E-01 -3.1721E-01  1.3373E-01 -1.7282E-01 -4.7604E-02  2.6767E-01 -7.2727E-01  2.0658E-02 -7.8243E-02
             7.8874E-02
 GRADIENT:  -1.1514E+01  6.1062E-02  2.0139E+00 -3.2198E+00 -2.6894E+00  1.5702E+00  3.2156E-01  9.0553E-02  7.3175E-02 -4.1066E-01
             2.7850E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1681.39937418475        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      764
 NPARAMETR:  1.0010E+00  1.1209E+00  5.8234E-01  9.6490E-01  7.6165E-01  8.6085E-01  1.0849E+00  2.1574E-01  9.6742E-01  8.2753E-01
             9.7908E-01
 PARAMETER:  1.0101E-01  2.1413E-01 -4.4069E-01  6.4268E-02 -1.7227E-01 -4.9831E-02  1.8149E-01 -1.4337E+00  6.6878E-02 -8.9312E-02
             7.8854E-02
 GRADIENT:  -5.4204E+00  5.1010E+00  8.8270E-01  3.2447E+00 -5.0737E+00  4.5866E-01  2.8821E-01  9.5767E-02 -2.0636E-01  6.0727E-01
             5.6496E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1681.49065657790        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      939
 NPARAMETR:  1.0023E+00  1.2662E+00  5.5687E-01  8.7539E-01  8.2096E-01  8.5989E-01  9.8107E-01  9.8403E-02  1.0467E+00  8.6587E-01
             9.7835E-01
 PARAMETER:  1.0232E-01  3.3603E-01 -4.8542E-01 -3.3087E-02 -9.7278E-02 -5.0956E-02  8.0892E-02 -2.2187E+00  1.4561E-01 -4.4022E-02
             7.8109E-02
 GRADIENT:  -1.3524E+00  2.9514E+00  9.1700E-01  1.4827E+00 -2.1065E+00  7.9911E-02 -5.0579E-02  1.6135E-02 -1.8425E-01 -2.0129E-01
             2.0210E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1681.51953876050        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1116
 NPARAMETR:  1.0028E+00  1.3784E+00  5.3616E-01  8.0538E-01  8.7166E-01  8.5941E-01  9.1558E-01  4.2053E-02  1.1189E+00  9.0446E-01
             9.7835E-01
 PARAMETER:  1.0281E-01  4.2095E-01 -5.2333E-01 -1.1644E-01 -3.7360E-02 -5.1514E-02  1.1804E-02 -3.0688E+00  2.1237E-01 -4.2265E-04
             7.8110E-02
 GRADIENT:   3.8528E-01  3.1678E-01  2.7082E-02  4.2675E-01 -4.5998E-02 -6.6480E-02 -5.4112E-02  3.4900E-03 -7.6803E-03 -3.4536E-02
            -7.7339E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1681.52801725969        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1299
 NPARAMETR:  1.0026E+00  1.5005E+00  5.0914E-01  7.2927E-01  9.2757E-01  8.5935E-01  8.5738E-01  1.2865E-02  1.2079E+00  9.4677E-01
             9.7906E-01
 PARAMETER:  1.0257E-01  5.0578E-01 -5.7503E-01 -2.1571E-01  2.4811E-02 -5.1578E-02 -5.3869E-02 -4.2533E+00  2.8890E-01  4.5298E-02
             7.8837E-02
 GRADIENT:  -2.6494E-03  2.7089E-01 -1.6583E-02  2.2474E-01  7.0418E-02 -1.8071E-02 -1.9085E-02  3.4260E-04 -5.7744E-02  5.4811E-03
            -4.9666E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1681.52815974511        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1466
 NPARAMETR:  1.0026E+00  1.4991E+00  5.0916E-01  7.2987E-01  9.2672E-01  8.5939E-01  8.5798E-01  1.0000E-02  1.2072E+00  9.4581E-01
             9.7915E-01
 PARAMETER:  1.0257E-01  5.0487E-01 -5.7499E-01 -2.1489E-01  2.3894E-02 -5.1530E-02 -5.3177E-02 -5.1929E+00  2.8832E-01  4.4288E-02
             7.8929E-02
 GRADIENT:   1.7458E-02 -1.5725E-02  7.9045E-03 -3.0940E-02 -2.0384E-02  1.4773E-03 -2.0756E-03  0.0000E+00 -2.3491E-03  1.1852E-03
             1.5552E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1466
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3494E-04 -2.5269E-02 -3.4885E-04  1.8426E-02 -2.9154E-02
 SE:             2.9802E-02  2.3027E-02  1.3432E-04  2.4306E-02  2.2576E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9639E-01  2.7248E-01  9.3973E-03  4.4840E-01  1.9658E-01

 ETASHRINKSD(%)  1.6095E-01  2.2858E+01  9.9550E+01  1.8572E+01  2.4367E+01
 ETASHRINKVR(%)  3.2164E-01  4.0491E+01  9.9998E+01  3.3695E+01  4.2796E+01
 EBVSHRINKSD(%)  5.5111E-01  2.2704E+01  9.9614E+01  1.9171E+01  2.3425E+01
 EBVSHRINKVR(%)  1.0992E+00  4.0253E+01  9.9999E+01  3.4666E+01  4.1363E+01
 RELATIVEINF(%)  9.8877E+01  3.5074E+00  1.5427E-04  4.2061E+00  7.6563E+00
 EPSSHRINKSD(%)  4.5296E+01
 EPSSHRINKVR(%)  7.0075E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1681.5281597451103     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -946.37733318137214     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1681.528       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.50E+00  5.09E-01  7.30E-01  9.27E-01  8.59E-01  8.58E-01  1.00E-02  1.21E+00  9.46E-01  9.79E-01
 


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
+        1.48E+03
 
 TH 2
+       -7.17E+00  4.20E+02
 
 TH 3
+        6.28E+00  2.43E+02  7.07E+02
 
 TH 4
+       -1.60E+01  3.01E+02 -3.92E+02  9.39E+02
 
 TH 5
+       -2.43E+00 -3.06E+02 -6.16E+02  3.36E+02  8.51E+02
 
 TH 6
+        1.90E+00 -1.26E+00  2.31E+00 -4.26E+00 -3.96E+00  2.66E+02
 
 TH 7
+       -2.64E+00  1.19E+01 -4.05E+01 -5.51E-01  3.77E-01 -2.69E+00  1.05E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.55E+00 -2.31E+01 -3.71E+01  5.03E+01 -7.34E+00  4.59E-01  1.83E+01  0.00E+00  6.81E+01
 
 TH10
+       -1.39E+00 -1.67E+01 -5.41E+01 -1.00E+01 -6.12E+01  7.47E-01  2.02E+01  0.00E+00  7.30E+00  8.76E+01
 
 TH11
+       -9.36E+00 -1.35E+01 -2.58E+01 -6.31E-02 -1.91E+00  2.52E+00  7.62E+00  0.00E+00  6.06E+00  1.72E+01  2.20E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.528
Stop Time:
Sat Sep 25 09:58:21 CDT 2021
