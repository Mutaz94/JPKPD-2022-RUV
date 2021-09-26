Sat Sep 25 08:58:52 CDT 2021
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
$DATA ../../../../data/spa/A2/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1112.96098494114        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0800E+02 -6.5070E+01  1.8517E+01 -1.0826E+02  2.2091E+01 -1.6925E+01  2.1115E+00 -7.7183E+00  8.8996E+00 -4.9968E+00
            -9.5330E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1412.20690775785        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6031E-01  1.1070E+00  1.0569E+00  1.0207E+00  1.0558E+00  9.8480E-01  9.1911E-01  9.8299E-01  7.6356E-01  8.1688E-01
             1.9259E+00
 PARAMETER:  5.9505E-02  2.0164E-01  1.5534E-01  1.2046E-01  1.5428E-01  8.4684E-02  1.5652E-02  8.2846E-02 -1.6976E-01 -1.0226E-01
             7.5541E-01
 GRADIENT:   7.7396E+01 -3.4162E+00  1.0920E+01 -1.5006E+01  2.8752E+00 -4.8276E+00 -8.1408E+00 -5.5862E+00 -9.2061E+00  7.9623E-01
            -1.1900E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1419.20128575619        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.5536E-01  8.6214E-01  5.9866E-01  1.1535E+00  6.6156E-01  9.7621E-01  1.2540E+00  8.2847E-01  7.7321E-01  2.3759E-01
             1.9811E+00
 PARAMETER:  5.4333E-02 -4.8335E-02 -4.1306E-01  2.4280E-01 -3.1316E-01  7.5920E-02  3.2632E-01 -8.8178E-02 -1.5721E-01 -1.3372E+00
             7.8366E-01
 GRADIENT:   5.4376E+01  2.0355E+01 -2.1643E-01  5.8681E+01 -7.8705E+00 -1.0691E+01 -1.1777E+00  5.4416E+00  1.8034E+01  8.6652E-01
            -7.5313E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1429.83839926498        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.2748E-01  9.9482E-01  2.8318E-01  9.5150E-01  5.0525E-01  9.9638E-01  9.9909E-01  3.9653E-01  6.4950E-01  1.5180E-01
             2.2615E+00
 PARAMETER:  2.4716E-02  9.4804E-02 -1.1617E+00  5.0284E-02 -5.8270E-01  9.6371E-02  9.9088E-02 -8.2501E-01 -3.3155E-01 -1.7852E+00
             9.1601E-01
 GRADIENT:  -8.8138E+00  4.9570E+00  1.7638E+00  1.5276E+00 -7.8042E+00  7.2528E-01  2.6373E+00  9.4741E-01  2.8205E+00  1.5028E+00
             1.5746E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1430.64441178587        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.2531E-01  1.0796E+00  2.4390E-01  8.7830E-01  5.1438E-01  9.9275E-01  9.2012E-01  3.5510E-01  6.7795E-01  7.1509E-02
             2.1844E+00
 PARAMETER:  2.2369E-02  1.7663E-01 -1.3110E+00 -2.9770E-02 -5.6479E-01  9.2725E-02  1.6745E-02 -9.3536E-01 -2.8868E-01 -2.5379E+00
             8.8136E-01
 GRADIENT:  -2.2740E+00 -7.2750E-01  6.5098E+00 -1.6401E+01 -1.3854E+01 -2.7125E-01  1.3992E-01 -1.9890E-02 -1.5632E-01  3.4589E-01
            -2.5977E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1433.10729239190        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.2599E-01  1.4439E+00  1.8293E-01  6.8666E-01  6.4163E-01  9.8840E-01  7.5412E-01  2.9822E-01  8.0577E-01  1.0000E-02
             2.2337E+00
 PARAMETER:  2.3108E-02  4.6738E-01 -1.5986E+00 -2.7591E-01 -3.4375E-01  8.8333E-02 -1.8220E-01 -1.1099E+00 -1.1596E-01 -6.0191E+00
             9.0367E-01
 GRADIENT:   7.7872E+00  2.4149E+01 -1.2125E+00  2.9296E+01 -3.3161E+00  1.1613E-01 -2.4874E+00 -5.1036E-01 -5.6780E-01  0.0000E+00
             4.2957E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1446.05479795218        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  9.0836E-01  1.9615E+00  6.8511E-02  3.7584E-01  8.3188E-01  9.6462E-01  6.4878E-01  1.3303E-01  1.4756E+00  1.0000E-02
             2.0173E+00
 PARAMETER:  3.8883E-03  7.7369E-01 -2.5808E+00 -8.7860E-01 -8.4072E-02  6.3978E-02 -3.3266E-01 -1.9172E+00  4.8905E-01 -1.7355E+01
             8.0174E-01
 GRADIENT:  -1.0526E+01  1.1337E+02  8.1479E+00  2.4623E+01 -9.3039E+01 -1.2560E+01 -7.2505E+00 -4.0359E-01  2.1014E+00  0.0000E+00
            -1.4001E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1447.90512687850        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      589             RESET HESSIAN, TYPE I
 NPARAMETR:  9.0763E-01  1.9868E+00  6.4630E-02  3.6061E-01  9.1813E-01  1.0299E+00  6.8276E-01  1.2623E-01  1.4968E+00  1.0000E-02
             2.1089E+00
 PARAMETER:  3.0856E-03  7.8651E-01 -2.6391E+00 -9.1997E-01  1.4589E-02  1.2943E-01 -2.8162E-01 -1.9696E+00  5.0335E-01 -1.8195E+01
             8.4617E-01
 GRADIENT:  -1.8068E+01 -3.2666E+01 -2.2078E+00  5.6223E+01  7.2564E+01  1.1859E+01  1.4004E+01 -2.7530E-01 -3.7014E+00  0.0000E+00
             1.3119E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1450.07622540476        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  9.0764E-01  1.9866E+00  6.4635E-02  3.6057E-01  8.8792E-01  9.9305E-01  6.4811E-01  1.4661E-01  1.5578E+00  1.0000E-02
             2.0511E+00
 PARAMETER:  3.0895E-03  7.8645E-01 -2.6390E+00 -9.2007E-01 -1.8871E-02  9.3023E-02 -3.3370E-01 -1.8200E+00  5.4329E-01 -1.8195E+01
             8.1840E-01
 GRADIENT:  -1.5718E+01  2.0954E+01  3.6721E+00  4.3109E+01  7.5995E+00 -1.6247E+00 -4.8375E-01 -4.2489E-01  9.4475E-01  0.0000E+00
             4.5088E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1471.44413528117        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      745
 NPARAMETR:  9.0765E-01  1.9853E+00  6.4618E-02  3.6011E-01  9.0318E-01  9.8203E-01  6.5344E-01  2.1922E+00  1.6407E+00  1.0000E-02
             1.7752E+00
 PARAMETER:  3.1011E-03  7.8577E-01 -2.6393E+00 -9.2134E-01 -1.8380E-03  8.1865E-02 -3.2550E-01  8.8489E-01  5.9513E-01 -1.8195E+01
             6.7394E-01
 GRADIENT:  -6.1961E+00  4.9632E+00 -3.3587E+01  4.4252E+01  2.6143E+00 -1.8378E+00 -1.4344E+00  1.0056E+00 -2.2525E-01  0.0000E+00
            -2.2170E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1471.48421658833        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      821
 NPARAMETR:  9.0765E-01  1.9853E+00  6.4637E-02  3.6008E-01  9.0122E-01  9.8557E-01  6.5665E-01  2.1417E+00  1.6331E+00  1.0000E-02
             1.7841E+00
 PARAMETER:  3.1020E-03  7.8576E-01 -2.6390E+00 -9.2143E-01 -4.0104E-03  8.5463E-02 -3.2060E-01  8.6161E-01  5.9049E-01 -1.8195E+01
             6.7894E-01
 GRADIENT:  -6.2009E+00  8.0422E+00 -3.2742E+01  4.3772E+01  9.0846E-02 -5.1791E-01 -9.8681E-02 -4.9219E-02 -2.1386E-01  0.0000E+00
            -3.2563E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1476.26629379775        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  9.0826E-01  1.9575E+00  7.9137E-02  3.3582E-01  9.1580E-01  9.7385E-01  6.6157E-01  2.3740E+00  1.6038E+00  1.0000E-02
             1.7348E+00
 PARAMETER:  3.7805E-03  7.7166E-01 -2.4366E+00 -9.9117E-01  1.2038E-02  7.3502E-02 -3.1315E-01  9.6459E-01  5.7239E-01 -1.8195E+01
             6.5087E-01
 GRADIENT:  -7.7174E+00 -5.8036E+01 -9.4911E+00 -1.1512E+01  6.7506E+00 -3.7848E+00 -8.5950E-01 -4.1186E-01 -1.1701E+00  0.0000E+00
            -1.5160E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1476.32074374015        NO. OF FUNC. EVALS.: 112
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  9.0827E-01  1.9571E+00  7.9405E-02  3.3543E-01  9.1397E-01  9.8551E-01  6.6424E-01  2.3989E+00  1.6205E+00  1.0000E-02
             1.7357E+00
 PARAMETER:  3.7918E-03  7.7144E-01 -2.4332E+00 -9.9233E-01  1.0040E-02  8.5405E-02 -3.0911E-01  9.7500E-01  5.8273E-01 -1.8195E+01
             6.5141E-01
 GRADIENT:  -1.7859E+01 -8.0411E+01 -1.0468E+01 -1.4829E+01 -2.5095E-01 -2.9456E-01 -2.0844E-02  6.7401E-03 -8.3084E-02  0.0000E+00
            -1.8192E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1478.14411622474        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1197
 NPARAMETR:  9.0848E-01  2.0234E+00  8.7499E-02  3.3204E-01  9.4801E-01  9.7850E-01  6.5303E-01  2.6377E+00  1.6242E+00  1.0000E-02
             1.7307E+00
 PARAMETER:  4.0193E-03  8.0478E-01 -2.3361E+00 -1.0025E+00  4.6611E-02  7.8262E-02 -3.2613E-01  1.0699E+00  5.8502E-01 -1.8195E+01
             6.4853E-01
 GRADIENT:  -2.4524E+01  5.7343E+00 -2.4715E+00 -5.7721E+00 -3.7948E+00 -3.7094E+00 -1.8533E+00 -6.3624E-01 -2.0151E+00  0.0000E+00
            -4.2611E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1478.45159928591        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1352
 NPARAMETR:  9.1929E-01  2.0190E+00  9.1192E-02  3.3547E-01  9.5500E-01  9.8629E-01  6.5497E-01  2.7108E+00  1.6404E+00  1.0000E-02
             1.7402E+00
 PARAMETER:  1.5842E-02  8.0261E-01 -2.2948E+00 -9.9222E-01  5.3954E-02  8.6196E-02 -3.2316E-01  1.0972E+00  5.9496E-01 -1.8195E+01
             6.5399E-01
 GRADIENT:   1.1720E+01  1.7100E+01 -6.0487E-01 -1.1648E+00  9.5211E+00  7.6828E-01 -1.9197E-02 -1.0260E-01  1.0537E+00  0.0000E+00
             6.0230E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1478.81126132923        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  9.1015E-01  1.9159E+00  1.0539E-01  3.8163E-01  9.0087E-01  9.7414E-01  6.6523E-01  2.7230E+00  1.4391E+00  1.0000E-02
             1.6965E+00
 PARAMETER:  5.8571E-03  7.5017E-01 -2.1501E+00 -8.6331E-01 -4.3911E-03  7.3804E-02 -3.0762E-01  1.1017E+00  4.6401E-01 -1.8195E+01
             6.2854E-01
 GRADIENT:  -6.6752E+00 -2.4617E+01 -2.5619E+00 -8.2155E+00  1.3943E+01 -3.3863E+00 -2.3887E+00 -4.5545E-01 -1.2400E-01  0.0000E+00
             5.5984E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1478.92624720625        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     1572             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1737E-01  1.9139E+00  1.0580E-01  3.8256E-01  8.9991E-01  9.8579E-01  6.7266E-01  2.7231E+00  1.4358E+00  1.0000E-02
             1.6958E+00
 PARAMETER:  1.3753E-02  7.4912E-01 -2.1462E+00 -8.6088E-01 -5.4618E-03  8.5692E-02 -2.9652E-01  1.1018E+00  4.6175E-01 -1.8195E+01
             6.2818E-01
 GRADIENT:   1.1470E+01 -2.5269E+01 -2.5912E+00 -8.6769E+00  1.3206E+01  1.2669E+00  3.2832E-01 -3.6297E-01  1.2388E-01  0.0000E+00
             8.7419E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1479.27118233833        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1756
 NPARAMETR:  9.1778E-01  1.9372E+00  1.0795E-01  3.8356E-01  8.9968E-01  9.8371E-01  6.7050E-01  2.7213E+00  1.4367E+00  1.0000E-02
             1.6971E+00
 PARAMETER:  1.4206E-02  7.6126E-01 -2.1261E+00 -8.5827E-01 -5.7153E-03  8.3573E-02 -2.9973E-01  1.1011E+00  4.6237E-01 -1.8195E+01
             6.2894E-01
 GRADIENT:  -9.4321E-02  2.9803E+00 -1.5853E-01 -7.8206E+00 -1.0911E+01 -9.1928E-01 -3.2848E-01 -2.0829E+00 -1.8884E-01  0.0000E+00
            -3.3961E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1479.28787419152        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1948             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1788E-01  1.9367E+00  1.0792E-01  3.8389E-01  9.0019E-01  9.8615E-01  6.7140E-01  2.7249E+00  1.4369E+00  1.0000E-02
             1.6971E+00
 PARAMETER:  1.4317E-02  7.6096E-01 -2.1264E+00 -8.5740E-01 -5.1543E-03  8.6049E-02 -2.9840E-01  1.1024E+00  4.6250E-01 -1.8195E+01
             6.2892E-01
 GRADIENT:   1.1405E+01  2.5422E+01  8.4750E-01 -5.0381E+00 -8.4988E+00  1.1695E+00  4.0384E-01 -1.0544E+00  9.7476E-02  0.0000E+00
             5.7446E-01

0ITERATION NO.:   92    OBJECTIVE VALUE:  -1479.28787419152        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     2017
 NPARAMETR:  9.1779E-01  1.9356E+00  1.0808E-01  3.8366E-01  9.0012E-01  9.8610E-01  6.7126E-01  2.7227E+00  1.4374E+00  1.0000E-02
             1.6979E+00
 PARAMETER:  1.4317E-02  7.6096E-01 -2.1264E+00 -8.5740E-01 -5.1543E-03  8.6049E-02 -2.9840E-01  1.1024E+00  4.6250E-01 -1.8195E+01
             6.2892E-01
 GRADIENT:   7.7717E-02  8.7649E+02 -1.5701E+02  7.8242E+02  6.6412E+03  1.1709E-02  3.7073E-02  6.0492E+02 -7.1854E+02  0.0000E+00
            -1.0597E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2017
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.7530E-04 -1.7497E-02  1.2196E-04  2.8809E-02 -5.6541E-04
 SE:             2.9553E-02  2.6768E-02  1.6798E-02  2.1229E-02  3.0139E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8177E-01  5.1334E-01  9.9421E-01  1.7477E-01  6.0651E-02

 ETASHRINKSD(%)  9.9444E-01  1.0323E+01  4.3725E+01  2.8879E+01  9.8990E+01
 ETASHRINKVR(%)  1.9790E+00  1.9580E+01  6.8331E+01  4.9419E+01  9.9990E+01
 EBVSHRINKSD(%)  1.2179E+00  1.0849E+01  4.2787E+01  2.9631E+01  9.9034E+01
 EBVSHRINKVR(%)  2.4209E+00  2.0520E+01  6.7266E+01  5.0482E+01  9.9991E+01
 RELATIVEINF(%)  9.6571E+01  1.8437E+01  1.9473E+01  1.2777E+01  2.0252E-03
 EPSSHRINKSD(%)  3.9703E+01
 EPSSHRINKVR(%)  6.3643E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1479.2878741915160     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -744.13704762777786     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.03
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.46
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1479.288       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.18E-01  1.94E+00  1.08E-01  3.84E-01  9.00E-01  9.86E-01  6.71E-01  2.72E+00  1.44E+00  1.00E-02  1.70E+00
 


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
+        1.97E+07
 
 TH 2
+       -9.99E+02  7.75E+04
 
 TH 3
+        4.79E+03 -2.34E+03  3.17E+06
 
 TH 4
+       -2.07E+03  1.15E+03  2.76E+03  1.58E+06
 
 TH 5
+       -1.95E+04  4.95E+03 -2.86E+04  1.31E+04  2.06E+07
 
 TH 6
+        1.83E+07 -5.65E+02  2.79E+03 -1.17E+03 -1.11E+04  1.94E+02
 
 TH 7
+        9.02E+06 -3.36E+02  1.51E+03 -4.57E+02 -6.98E+03  8.40E+06  4.13E+06
 
 TH 8
+       -2.68E+02  1.34E+02 -6.81E+01  8.93E+02  1.57E+03 -1.53E+02 -6.97E+01  1.86E+04
 
 TH 9
+        2.74E+03 -5.53E+02  2.76E+03 -1.08E+03 -2.78E+06  1.56E+03  1.00E+03 -1.51E+02  3.77E+05
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.20E+03 -1.06E+05  4.28E+03 -9.69E+02 -6.86E+03  6.98E+02  4.21E+02 -1.70E+02  6.93E+02  0.00E+00  1.47E+05
 
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
 #CPUT: Total CPU Time in Seconds,       32.560
Stop Time:
Sat Sep 25 08:59:26 CDT 2021
