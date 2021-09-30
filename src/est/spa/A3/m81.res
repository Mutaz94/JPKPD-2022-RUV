Wed Sep 29 13:50:04 CDT 2021
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
$DATA ../../../../data/spa/A3/dat81.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -106.993916370224        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2931E+02  7.7061E+01  6.8912E+01  2.8784E+01  2.2328E+02  1.2637E+01 -6.1042E+01 -2.9860E+01 -1.2293E+02 -1.7773E+02
            -2.6637E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1173.53804082209        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0673E+00  8.6660E-01  9.1806E-01  1.1050E+00  8.5751E-01  1.0407E+00  9.3986E-01  9.4515E-01  1.2949E+00  8.6818E-01
             2.3598E+00
 PARAMETER:  1.6515E-01 -4.3183E-02  1.4504E-02  1.9984E-01 -5.3725E-02  1.3985E-01  3.7979E-02  4.3584E-02  3.5844E-01 -4.1354E-02
             9.5859E-01
 GRADIENT:   2.2291E+02 -6.6187E+00  9.9324E-01  1.7047E+01  9.9753E+01 -6.5601E+00 -3.6664E-02 -3.5612E+00  4.4216E+00 -4.3094E+01
            -3.5005E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1204.63826245252        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0691E+00  4.4408E-01  4.0296E-01  1.3741E+00  3.5425E-01  9.2898E-01  2.5582E-01  7.4762E-01  1.2623E+00  2.4944E-01
             2.6038E+00
 PARAMETER:  1.6685E-01 -7.1176E-01 -8.0891E-01  4.1779E-01 -9.3774E-01  2.6336E-02 -1.2633E+00 -1.9086E-01  3.3296E-01 -1.2885E+00
             1.0570E+00
 GRADIENT:   1.6813E+02  1.4828E+02  2.1986E+02  2.1693E+02 -2.7885E+02 -6.4909E+01 -1.9193E+00 -4.2148E+01  1.1840E+01 -1.9098E+01
            -2.4692E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1266.15985128832        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      255
 NPARAMETR:  1.0225E+00  3.9612E-01  2.6223E-01  1.2209E+00  2.9783E-01  1.0383E+00  2.6309E-02  3.6170E-01  1.1765E+00  1.4594E-01
             3.7272E+00
 PARAMETER:  1.2224E-01 -8.2604E-01 -1.2385E+00  2.9959E-01 -1.1112E+00  1.3763E-01 -3.5378E+00 -9.1695E-01  2.6256E-01 -1.8245E+00
             1.4157E+00
 GRADIENT:  -4.5290E+01  1.3568E+01 -1.6520E+01  7.4029E+01 -1.2452E+01 -1.7492E+01 -8.3582E-03 -1.6338E+00  9.2017E-01 -2.3433E+00
             1.2395E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1274.86193897018        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  1.0495E+00  4.0311E-01  2.3913E-01  1.1406E+00  2.8355E-01  1.1464E+00  5.4469E-02  3.7384E-01  1.2283E+00  5.3303E-01
             3.2221E+00
 PARAMETER:  1.4828E-01 -8.0854E-01 -1.3307E+00  2.3152E-01 -1.1604E+00  2.3660E-01 -2.8101E+00 -8.8394E-01  3.0567E-01 -5.2917E-01
             1.2700E+00
 GRADIENT:   1.5033E+01 -1.5928E+01 -3.7593E+01  1.9788E+01  6.9549E+01  1.6277E+01 -6.0231E-03  6.3638E-01  4.4772E+00 -2.1905E+00
            -1.3941E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1276.91006305849        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  1.0412E+00  3.6855E-01  2.1501E-01  1.0907E+00  2.5452E-01  1.0940E+00  5.3940E-02  3.6625E-01  1.2336E+00  5.7307E-01
             3.2258E+00
 PARAMETER:  1.4036E-01 -8.9819E-01 -1.4371E+00  1.8686E-01 -1.2684E+00  1.8981E-01 -2.8199E+00 -9.0444E-01  3.0990E-01 -4.5675E-01
             1.2712E+00
 GRADIENT:  -1.1799E+00 -2.0026E+00 -2.7329E+00 -2.5416E+00  4.8630E+00  4.8039E-01 -1.7335E-03  7.5687E-01 -1.9676E-01 -2.4546E-01
            -2.8544E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1277.17717763411        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  1.0410E+00  3.6829E-01  2.3100E-01  1.1180E+00  2.6399E-01  1.0874E+00  1.1415E-01  1.7358E-01  1.2087E+00  5.8767E-01
             3.2822E+00
 PARAMETER:  1.4018E-01 -8.9889E-01 -1.3653E+00  2.1153E-01 -1.2318E+00  1.8382E-01 -2.0702E+00 -1.6511E+00  2.8953E-01 -4.3160E-01
             1.2885E+00
 GRADIENT:  -4.0234E+00  1.7641E+00  5.0447E-01  1.8466E+00 -1.9782E+00 -8.7079E-01 -5.7236E-03  1.6187E-01  1.8939E-01  4.4387E-01
             5.5472E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1277.40077987836        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  1.0429E+00  2.8465E-01  2.7386E-01  1.1966E+00  2.7578E-01  1.0804E+00  4.1106E-01  1.9139E-02  1.1381E+00  5.9405E-01
             3.2864E+00
 PARAMETER:  1.4199E-01 -1.1565E+00 -1.1951E+00  2.7951E-01 -1.1882E+00  1.7733E-01 -7.8901E-01 -3.8560E+00  2.2934E-01 -4.2079E-01
             1.2898E+00
 GRADIENT:   1.3351E-01  1.2414E+00  2.3757E+00  1.2352E+00 -4.2228E+00 -4.0539E-01 -1.2800E-02  2.0273E-03 -5.7226E-02 -3.4209E-01
             2.5319E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1277.43280176170        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1116
 NPARAMETR:  1.0415E+00  2.6445E-01  2.8306E-01  1.2122E+00  2.7899E-01  1.0789E+00  4.9705E-01  1.0000E-02  1.1234E+00  5.9540E-01
             3.2903E+00
 PARAMETER:  1.4066E-01 -1.2301E+00 -1.1621E+00  2.9246E-01 -1.1766E+00  1.7597E-01 -5.9907E-01 -4.6530E+00  2.1632E-01 -4.1852E-01
             1.2910E+00
 GRADIENT:   4.3519E+01  3.7311E+00  1.6260E+01  2.7670E+01  5.2385E+01  6.1344E+00  3.1680E-02  0.0000E+00  2.8662E+00  3.0946E-01
             1.0325E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1277.46818805518        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1278
 NPARAMETR:  1.0417E+00  2.3639E-01  2.8277E-01  1.2212E+00  2.7539E-01  1.0801E+00  2.3945E-01  1.0000E-02  1.1221E+00  6.1128E-01
             3.2735E+00
 PARAMETER:  1.4084E-01 -1.3423E+00 -1.1631E+00  2.9983E-01 -1.1896E+00  1.7710E-01 -1.3294E+00 -4.6530E+00  2.1518E-01 -3.9220E-01
             1.2859E+00
 GRADIENT:   1.8371E+00 -1.5530E-01 -3.3993E+00  1.8456E-01  3.3827E+00  6.1137E-01  4.6578E-03  0.0000E+00  4.1194E-01  5.3107E-01
            -8.0263E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1277.56153288241        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1454
 NPARAMETR:  1.0371E+00  1.6456E-01  2.9722E-01  1.2624E+00  2.7593E-01  1.0727E+00  1.0000E-02  1.0000E-02  1.0812E+00  5.9227E-01
             3.3302E+00
 PARAMETER:  1.3640E-01 -1.7045E+00 -1.1133E+00  3.3303E-01 -1.1876E+00  1.7016E-01 -6.3762E+00 -4.6530E+00  1.7806E-01 -4.2378E-01
             1.3030E+00
 GRADIENT:   3.2871E-01  9.5062E-01  3.7803E+00  2.8998E+00 -1.1274E+01 -1.1989E-02  0.0000E+00  0.0000E+00 -1.2846E+00 -1.3935E+00
             4.2917E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1277.78430435993        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1630
 NPARAMETR:  1.0310E+00  1.0853E-01  3.2232E-01  1.3063E+00  2.8765E-01  1.0680E+00  1.0000E-02  1.0000E-02  1.0543E+00  6.2231E-01
             3.2976E+00
 PARAMETER:  1.3051E-01 -2.1207E+00 -1.0322E+00  3.6723E-01 -1.1460E+00  1.6579E-01 -1.4501E+01 -4.6530E+00  1.5287E-01 -3.7432E-01
             1.2932E+00
 GRADIENT:  -8.2790E-02  2.5170E-01  1.8775E+00  5.5644E-01 -1.9431E+00  4.7250E-03  0.0000E+00  0.0000E+00  1.1945E-01  3.4253E-01
            -6.8752E-01

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1277.80582661444        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1723
 NPARAMETR:  1.0290E+00  8.1953E-02  3.2064E-01  1.3133E+00  2.8473E-01  1.0667E+00  1.0000E-02  1.0000E-02  1.0472E+00  6.1783E-01
             3.3053E+00
 PARAMETER:  1.2862E-01 -2.4016E+00 -1.0374E+00  3.7258E-01 -1.1562E+00  1.6461E-01 -1.9433E+01 -4.6530E+00  1.4614E-01 -3.8155E-01
             1.2955E+00
 GRADIENT:  -2.2962E-01  3.1015E-02  7.2736E-02  5.9450E-01 -7.1075E-01 -9.6846E-03  0.0000E+00  0.0000E+00 -1.4005E-01 -1.2548E-01
             1.3730E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1723
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.4042E-04 -1.5268E-05  9.4102E-05 -1.2710E-02  6.4870E-04
 SE:             2.8818E-02  1.0283E-05  2.2007E-04  2.6256E-02  1.9020E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9334E-01  1.3762E-01  6.6894E-01  6.2835E-01  9.7279E-01

 ETASHRINKSD(%)  3.4548E+00  9.9966E+01  9.9263E+01  1.2038E+01  3.6281E+01
 ETASHRINKVR(%)  6.7902E+00  1.0000E+02  9.9995E+01  2.2627E+01  5.9399E+01
 EBVSHRINKSD(%)  3.2488E+00  9.9966E+01  9.9256E+01  1.0808E+01  3.5965E+01
 EBVSHRINKVR(%)  6.3920E+00  1.0000E+02  9.9994E+01  2.0448E+01  5.8995E+01
 RELATIVEINF(%)  7.8603E+01  1.3936E-06  2.0599E-04  1.7197E+01  1.1522E+00
 EPSSHRINKSD(%)  2.8988E+01
 EPSSHRINKVR(%)  4.9573E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1277.8058266144351     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -542.65500005069691     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1277.806       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  8.20E-02  3.21E-01  1.31E+00  2.85E-01  1.07E+00  1.00E-02  1.00E-02  1.05E+00  6.18E-01  3.31E+00
 


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
+        8.65E+02
 
 TH 2
+       -9.48E+01  2.03E+02
 
 TH 3
+       -4.42E+01  5.69E+02  5.92E+03
 
 TH 4
+       -3.60E+01  1.81E+02 -3.34E+02  4.68E+02
 
 TH 5
+        2.03E+02 -1.26E+03 -8.53E+03 -3.11E+02  1.46E+04
 
 TH 6
+       -1.18E+00 -1.07E+01  2.26E+01 -1.06E+01 -4.67E-01  1.51E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.26E+00 -3.68E+01  7.35E+01 -1.31E+01  4.96E+01  5.97E-01  0.00E+00  0.00E+00  1.13E+02
 
 TH10
+       -9.54E+00 -3.87E+00 -1.19E+02  9.35E-01  1.73E+02  2.00E+00  0.00E+00  0.00E+00 -1.13E+00  8.40E+01
 
 TH11
+       -1.33E+01 -2.22E+00 -1.16E+01 -6.47E+00  4.77E+00  2.73E+00  0.00E+00  0.00E+00  7.39E+00  2.17E+01  2.89E+01
 
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
 #CPUT: Total CPU Time in Seconds,       27.604
Stop Time:
Wed Sep 29 13:50:33 CDT 2021
