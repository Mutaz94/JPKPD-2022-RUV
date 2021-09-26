Sat Sep 25 09:31:47 CDT 2021
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
$DATA ../../../../data/spa/A3/dat82.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m82.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -333.778993844099        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.9611E+01  1.5220E+01  3.5703E+01 -2.6350E+01  2.3569E+02  5.4137E+01 -7.7869E+01 -1.7612E+01 -1.2506E+02 -1.8429E+02
            -2.2438E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1279.57878644784        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0104E+00  9.6145E-01  1.0522E+00  1.0682E+00  9.0251E-01  6.4757E-01  1.0320E+00  9.3944E-01  1.1268E+00  9.0277E-01
             3.6365E+00
 PARAMETER:  1.1036E-01  6.0683E-02  1.5092E-01  1.6601E-01 -2.5737E-03 -3.3453E-01  1.3148E-01  3.7534E-02  2.1936E-01 -2.2922E-03
             1.3910E+00
 GRADIENT:   1.6425E+01 -1.1128E+01 -7.5101E+00 -1.4506E+01  5.3234E+00 -5.7098E+01  4.0305E+00  3.8399E+00  7.8506E+00  1.2989E+01
             7.1011E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1290.78146840954        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0195E+00  6.8729E-01  5.6679E-01  1.1931E+00  5.5131E-01  7.3429E-01  1.2163E+00  2.4665E-01  1.0148E+00  4.2213E-01
             3.6259E+00
 PARAMETER:  1.1928E-01 -2.7500E-01 -4.6777E-01  2.7655E-01 -4.9545E-01 -2.0885E-01  2.9585E-01 -1.2998E+00  1.1471E-01 -7.6245E-01
             1.3881E+00
 GRADIENT:   2.0813E+01  4.3256E+00 -1.3309E+00  5.7384E+00  5.5703E-01 -2.1865E+01 -2.5691E+00  4.3442E-01  4.0236E+00  1.5362E+00
             1.7667E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1297.50257015491        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0077E+00  5.2911E-01  2.5939E-01  1.1603E+00  3.1492E-01  8.8152E-01  1.7288E+00  5.3575E-02  9.9812E-01  2.9630E-01
             3.1664E+00
 PARAMETER:  1.0768E-01 -5.3656E-01 -1.2494E+00  2.4868E-01 -1.0554E+00 -2.6110E-02  6.4746E-01 -2.8267E+00  9.8117E-02 -1.1164E+00
             1.2526E+00
 GRADIENT:  -1.8901E+01  3.1277E+01 -8.9847E+00  8.2553E+01  1.6631E-01  2.2469E+01  1.0702E+01 -5.2111E-03 -1.1732E+01 -1.3935E+00
             5.3472E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1306.07704607798        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0030E+00  4.0674E-01  1.3470E-01  9.5863E-01  2.1088E-01  8.2843E-01  1.2778E+00  1.0000E-02  1.2950E+00  4.8929E-01
             2.7883E+00
 PARAMETER:  1.0302E-01 -7.9959E-01 -1.9047E+00  5.7755E-02 -1.4565E+00 -8.8219E-02  3.4517E-01 -5.9486E+00  3.5849E-01 -6.1479E-01
             1.1254E+00
 GRADIENT:   7.4998E+00  6.7347E+00 -4.5263E+00  9.1781E+00  7.1177E+00 -7.2199E+00  1.5890E+00  0.0000E+00 -3.4297E+00 -3.4162E+00
            -2.2474E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1306.84882365449        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0084E+00  3.3084E-01  1.0144E-01  8.9787E-01  1.7203E-01  8.1946E-01  1.0918E+00  1.0000E-02  1.4863E+00  7.2639E-01
             2.8212E+00
 PARAMETER:  1.0837E-01 -1.0061E+00 -2.1883E+00 -7.7277E-03 -1.6601E+00 -9.9107E-02  1.8785E-01 -7.9377E+00  4.9627E-01 -2.1967E-01
             1.1372E+00
 GRADIENT:   1.7117E+01  1.2260E+01 -9.8224E+00  3.1790E+01  2.8016E+00 -8.8610E+00  1.0815E+01  0.0000E+00 -2.4399E+00  6.5959E+00
             2.5169E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1307.39172012106        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  1.0086E+00  3.0212E-01  9.0594E-02  8.4703E-01  1.5886E-01  8.2806E-01  9.4258E-01  1.0000E-02  1.5478E+00  7.6216E-01
             2.9195E+00
 PARAMETER:  1.0859E-01 -1.0969E+00 -2.3014E+00 -6.6015E-02 -1.7397E+00 -8.8668E-02  4.0861E-02 -8.6309E+00  5.3683E-01 -1.7160E-01
             1.1714E+00
 GRADIENT:   9.2867E+00  1.0098E+01 -3.6233E+00  2.4580E+01 -1.1955E+01 -3.3023E+00  7.6197E+00  0.0000E+00 -1.2358E+00  7.4307E+00
             1.7489E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1307.81982301479        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  1.0068E+00  2.8283E-01  8.2073E-02  7.8983E-01  1.5006E-01  8.3623E-01  7.6562E-01  1.0000E-02  1.6030E+00  7.7551E-01
             2.9455E+00
 PARAMETER:  1.0678E-01 -1.1629E+00 -2.4001E+00 -1.3594E-01 -1.7967E+00 -7.8848E-02 -1.6708E-01 -9.1540E+00  5.7185E-01 -1.5424E-01
             1.1803E+00
 GRADIENT:   2.6648E+00  5.1233E+00  8.7985E-01  9.2019E+00 -1.1802E+01  1.0694E+00  4.4776E+00  0.0000E+00 -1.7077E-01  5.0206E+00
             1.6744E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1309.44347703640        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      652
 NPARAMETR:  1.0056E+00  3.0375E-01  9.1008E-02  8.0329E-01  1.6237E-01  8.4665E-01  6.6552E-01  1.0000E-02  1.5462E+00  7.6274E-01
             2.9162E+00
 PARAMETER:  1.0563E-01 -1.0915E+00 -2.2968E+00 -1.1905E-01 -1.7179E+00 -6.6466E-02 -3.0718E-01 -8.6278E+00  5.3581E-01 -1.7083E-01
             1.1703E+00
 GRADIENT:   2.5617E+00 -9.4769E+00 -4.0212E+00 -1.6647E+01  2.2582E+01  4.9527E+00  3.0238E+00  0.0000E+00  2.4237E+00  2.8300E+00
             1.4963E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1311.01003145895        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  1.0026E+00  3.3545E-01  9.7614E-02  8.5413E-01  1.7212E-01  8.3733E-01  3.6028E-01  1.0000E-02  1.5009E+00  7.6246E-01
             2.8289E+00
 PARAMETER:  1.0258E-01 -9.9229E-01 -2.2267E+00 -5.7676E-02 -1.6595E+00 -7.7534E-02 -9.2089E-01 -7.7537E+00  5.0606E-01 -1.7121E-01
             1.1399E+00
 GRADIENT:   3.4153E+00 -3.0768E+00 -4.4582E+00  5.1316E+00  9.6032E+00  4.1046E-02  4.8428E-01  0.0000E+00 -4.2147E-01 -2.0234E+00
            -2.6342E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1311.33810466722        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  1.0019E+00  3.2957E-01  9.4766E-02  8.3598E-01  1.6849E-01  8.3524E-01  8.1647E-02  1.0000E-02  1.5125E+00  7.8577E-01
             2.8453E+00
 PARAMETER:  1.0193E-01 -1.0100E+00 -2.2563E+00 -7.9155E-02 -1.6809E+00 -8.0039E-02 -2.4054E+00 -7.4649E+00  5.1379E-01 -1.4109E-01
             1.1457E+00
 GRADIENT:   4.6329E-01  5.6345E-01  7.1570E-01 -4.1263E-01 -1.6405E+00 -1.5446E-01  2.6864E-02  0.0000E+00 -3.7732E-01  9.2307E-01
            -3.2584E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1311.35447540510        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1177
 NPARAMETR:  1.0016E+00  3.2970E-01  9.4530E-02  8.3600E-01  1.6849E-01  8.3584E-01  1.1537E-02  1.0000E-02  1.5177E+00  7.8250E-01
             2.8480E+00
 PARAMETER:  1.0162E-01 -1.0096E+00 -2.2588E+00 -7.9130E-02 -1.6809E+00 -7.9320E-02 -4.3622E+00 -6.7674E+00  5.1717E-01 -1.4526E-01
             1.1466E+00
 GRADIENT:  -4.2819E-02  9.9834E-04  5.4290E-02 -5.5425E-02 -2.5017E-02  2.9285E-02  5.0847E-04  0.0000E+00  5.0238E-02 -1.4116E-02
             6.7238E-02

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1311.35454864864        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  1.0016E+00  3.2976E-01  9.4529E-02  8.3605E-01  1.6851E-01  8.3580E-01  1.0000E-02  1.0000E-02  1.5175E+00  7.8252E-01
             2.8479E+00
 PARAMETER:  1.0162E-01 -1.0094E+00 -2.2589E+00 -7.9064E-02 -1.6808E+00 -7.9368E-02 -4.6244E+00 -6.6737E+00  5.1709E-01 -1.4523E-01
             1.1466E+00
 GRADIENT:  -6.4328E-07 -9.7441E-03 -7.7713E-03 -1.1068E-02  5.1184E-02  1.0638E-02  0.0000E+00  0.0000E+00  2.5154E-02 -1.1491E-02
             4.0982E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1234
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.5597E-04 -1.3443E-05  3.3079E-04 -1.5088E-02  7.5432E-03
 SE:             2.8521E-02  1.4071E-04  1.4377E-04  2.4950E-02  2.5883E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7606E-01  9.2388E-01  2.1401E-02  5.4535E-01  7.7072E-01

 ETASHRINKSD(%)  4.4502E+00  9.9529E+01  9.9518E+01  1.6413E+01  1.3290E+01
 ETASHRINKVR(%)  8.7023E+00  9.9998E+01  9.9998E+01  3.0132E+01  2.4814E+01
 EBVSHRINKSD(%)  4.1812E+00  9.9510E+01  9.9537E+01  1.2308E+01  1.4402E+01
 EBVSHRINKVR(%)  8.1876E+00  9.9998E+01  9.9998E+01  2.3102E+01  2.6730E+01
 RELATIVEINF(%)  8.6159E+01  3.8626E-04  3.2984E-04  2.7733E+01  4.9680E+00
 EPSSHRINKSD(%)  3.4096E+01
 EPSSHRINKVR(%)  5.6567E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1311.3545486486407     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -576.20372208490255     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.15
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
 





 #OBJV:********************************************    -1311.355       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  3.30E-01  9.45E-02  8.36E-01  1.69E-01  8.36E-01  1.00E-02  1.00E-02  1.52E+00  7.83E-01  2.85E+00
 


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
+        1.46E+03
 
 TH 2
+        1.08E+02  2.44E+03
 
 TH 3
+       -7.80E+02  4.76E+03  2.82E+04
 
 TH 4
+       -1.04E+01 -5.02E+01 -1.46E+03  5.29E+02
 
 TH 5
+        5.32E+02 -9.19E+03 -2.67E+04 -4.37E+02  4.71E+04
 
 TH 6
+       -2.89E+00 -1.62E+01  7.13E+01 -1.67E+01  1.11E+01  2.33E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.29E+01 -1.91E+01  3.61E+02 -1.95E+01  5.62E+01  1.27E-01  0.00E+00  0.00E+00  4.50E+01
 
 TH10
+       -7.21E+00 -3.04E+01 -8.39E+00  1.33E+01  1.12E+02  1.09E+01  0.00E+00  0.00E+00  4.83E+00  1.72E+02
 
 TH11
+       -2.48E+01  9.56E-01  2.81E+01 -3.11E+00  5.36E+01  3.52E+00  0.00E+00  0.00E+00  5.92E+00  1.30E+01  3.44E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.455
Stop Time:
Sat Sep 25 09:32:08 CDT 2021
