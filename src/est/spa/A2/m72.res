Wed Sep 29 13:00:45 CDT 2021
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
$DATA ../../../../data/spa/A2/dat72.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1086.16595844537        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3275E+02  3.3359E+01  1.3065E+02 -7.9170E+01 -4.4961E+01  5.5495E+01 -5.2390E+00 -4.3067E+01 -2.2042E+01 -2.9120E+01
            -1.0900E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1455.80607308592        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0681E+00  1.0548E+00  8.7498E-01  1.0993E+00  1.0026E+00  8.5888E-01  9.2976E-01  1.0270E+00  8.8410E-01  9.0505E-01
             2.6772E+00
 PARAMETER:  1.6586E-01  1.5335E-01 -3.3556E-02  1.9467E-01  1.0257E-01 -5.2132E-02  2.7171E-02  1.2662E-01 -2.3180E-02  2.3695E-04
             1.0848E+00
 GRADIENT:   1.1703E+02  3.8101E+01  1.9882E+00  5.7382E+01 -1.2312E+01 -2.1503E+01  1.4574E+00  2.3006E+00 -2.2010E+00  1.0548E+01
             4.5758E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1459.55299089001        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0610E+00  1.0493E+00  5.2363E-01  1.0573E+00  6.9996E-01  8.8882E-01  9.6549E-01  9.3391E-01  9.0947E-01  4.3434E-01
             2.5830E+00
 PARAMETER:  1.5925E-01  1.4808E-01 -5.4696E-01  1.5572E-01 -2.5673E-01 -1.7857E-02  6.4876E-02  3.1620E-02  5.1094E-03 -7.3393E-01
             1.0490E+00
 GRADIENT:   9.1627E+01  7.6822E+01  3.5088E+01  5.7111E+01 -6.4870E+01 -1.1167E+01  1.0575E+00  4.6581E+00  5.9667E+00  4.4543E+00
             4.1324E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1462.36808545916        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      318
 NPARAMETR:  1.0518E+00  9.9116E-01  4.7657E-01  1.0657E+00  6.4887E-01  9.1753E-01  1.0883E+00  9.3026E-01  8.6628E-01  4.1072E-01
             2.4323E+00
 PARAMETER:  1.5049E-01  9.1118E-02 -6.4114E-01  1.6362E-01 -3.3253E-01  1.3932E-02  1.8462E-01  2.7705E-02 -4.3544E-02 -7.8983E-01
             9.8883E-01
 GRADIENT:  -1.9500E+01  5.2397E+01  2.0693E+01  3.2453E+01 -4.6623E+01 -5.6129E+00  7.7330E+00  6.5603E+00  3.4027E+00  4.8871E+00
             1.1946E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1478.55496901011        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      495
 NPARAMETR:  1.0768E+00  6.0866E-01  2.1008E-01  1.0490E+00  3.1857E-01  1.0005E+00  1.2711E+00  1.2917E+00  8.6321E-01  6.0901E-02
             1.9178E+00
 PARAMETER:  1.7396E-01 -3.9650E-01 -1.4603E+00  1.4786E-01 -1.0439E+00  1.0052E-01  3.3986E-01  3.5596E-01 -4.7098E-02 -2.6985E+00
             7.5119E-01
 GRADIENT:   6.7159E+01  5.6311E+01  7.2396E+01  6.4853E+00 -1.2002E+02  1.7804E+01  1.2188E+01 -1.2418E+01  1.4285E+00  1.8401E-01
             3.0512E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1491.34506074350        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  1.0314E+00  6.4594E-01  1.8131E-01  9.9640E-01  3.2894E-01  9.3567E-01  1.0947E+00  1.4443E+00  9.1680E-01  5.6395E-02
             1.6239E+00
 PARAMETER:  1.3096E-01 -3.3705E-01 -1.6076E+00  9.6393E-02 -1.0119E+00  3.3513E-02  1.9044E-01  4.6762E-01  1.3134E-02 -2.7754E+00
             5.8485E-01
 GRADIENT:  -6.9418E-01 -3.8604E+00  6.7301E+00 -5.9581E+00 -9.2716E+00 -2.7028E+00  3.7221E+00 -3.1139E+00  1.9664E-01  1.7211E-01
            -8.1564E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1499.26898264754        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.0314E+00  1.1506E+00  1.9578E-01  8.1777E-01  5.3379E-01  9.3265E-01  8.3531E-01  1.9456E+00  9.1629E-01  6.6482E-02
             1.6814E+00
 PARAMETER:  1.3096E-01  2.4032E-01 -1.5308E+00 -1.0117E-01 -5.2776E-01  3.0275E-02 -7.9948E-02  7.6556E-01  1.2580E-02 -2.6108E+00
             6.1963E-01
 GRADIENT:   1.2915E+01  5.5745E+00  1.1676E+01  3.1187E+01 -1.8640E+00  4.1461E+00 -4.1539E+00  1.1007E+00 -5.5833E+00  1.2136E-01
            -1.7404E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1502.28347955195        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  1.0306E+00  1.3184E+00  1.8176E-01  7.2063E-01  6.0010E-01  9.2696E-01  8.2097E-01  2.1213E+00  1.0178E+00  3.6857E-02
             1.6968E+00
 PARAMETER:  1.3010E-01  3.7640E-01 -1.6051E+00 -2.2763E-01 -4.1065E-01  2.4152E-02 -9.7267E-02  8.5201E-01  1.1768E-01 -3.2007E+00
             6.2873E-01
 GRADIENT:   1.4512E+01  2.0224E+01  1.5214E+01 -2.1369E+00 -2.7895E+01  2.9638E+00  6.8534E+00  3.8304E+00  4.4918E-01  2.4690E-02
            -1.7272E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1502.43008168514        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1212
 NPARAMETR:  1.0308E+00  1.3192E+00  1.8212E-01  7.1900E-01  6.0046E-01  9.1949E-01  7.9322E-01  2.1174E+00  1.0269E+00  1.4329E-02
             1.7010E+00
 PARAMETER:  1.3034E-01  3.7702E-01 -1.6031E+00 -2.2990E-01 -4.1005E-01  1.6065E-02 -1.3166E-01  8.5021E-01  1.2658E-01 -4.1454E+00
             6.3122E-01
 GRADIENT:   1.5313E+01  1.9371E+01  1.6846E+01 -3.4978E+00 -3.0069E+01 -2.6274E-02 -1.3356E-02  3.3358E+00  4.3070E-01  3.1323E-03
            -1.7078E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1503.18405254989        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1392
 NPARAMETR:  1.0309E+00  1.3085E+00  1.7291E-01  7.1871E-01  6.0367E-01  9.2432E-01  7.9244E-01  2.0989E+00  1.0203E+00  1.0000E-02
             1.7345E+00
 PARAMETER:  1.3038E-01  3.6890E-01 -1.6550E+00 -2.3030E-01 -4.0472E-01  2.1301E-02 -1.3263E-01  8.4140E-01  1.2014E-01 -5.6741E+01
             6.5070E-01
 GRADIENT:   1.5992E+01 -2.6327E+01  3.4473E-01  8.0447E+00  1.6044E+01  2.0566E+00  2.1053E-01  4.1059E+00  1.7167E-01  0.0000E+00
            -9.5201E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1503.75213705738        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1581
 NPARAMETR:  1.0249E+00  1.3195E+00  1.6827E-01  7.1278E-01  6.0322E-01  9.1922E-01  7.8749E-01  2.0066E+00  1.0184E+00  1.0000E-02
             1.7803E+00
 PARAMETER:  1.2457E-01  3.7723E-01 -1.6822E+00 -2.3858E-01 -4.0547E-01  1.5767E-02 -1.3891E-01  7.9645E-01  1.1825E-01 -5.5044E+01
             6.7676E-01
 GRADIENT:   6.8820E-01 -1.6654E+01 -1.6334E-01  1.0200E+01  9.9147E+00  3.8623E-02  2.7956E-01  5.6023E-01 -4.9348E-01  0.0000E+00
            -1.7022E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1505.18299356106        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1758
 NPARAMETR:  1.0210E+00  1.5331E+00  1.3664E-01  5.9185E-01  6.9231E-01  9.1893E-01  7.2951E-01  1.9933E+00  1.1996E+00  1.0000E-02
             1.8588E+00
 PARAMETER:  1.2081E-01  5.2731E-01 -1.8904E+00 -4.2450E-01 -2.6772E-01  1.5449E-02 -2.1538E-01  7.8981E-01  2.8198E-01 -5.5044E+01
             7.1991E-01
 GRADIENT:  -8.7270E+00 -1.8999E+00 -2.6802E+00 -6.2327E-01 -4.7213E+00  3.9252E-01 -3.4441E-01 -1.4600E+00  7.1982E-01  0.0000E+00
             1.3159E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1505.32714862791        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1920
 NPARAMETR:  1.0256E+00  1.5527E+00  1.4048E-01  5.8672E-01  7.0732E-01  9.1760E-01  7.3008E-01  2.0972E+00  1.1992E+00  1.0000E-02
             1.8559E+00
 PARAMETER:  1.2525E-01  5.3997E-01 -1.8627E+00 -4.3321E-01 -2.4627E-01  1.4008E-02 -2.1460E-01  8.4061E-01  2.8163E-01 -5.5044E+01
             7.1837E-01
 GRADIENT:  -1.3506E-01 -6.7631E-01  2.3445E-01 -3.0957E-01  4.9211E-01  1.3216E-03  2.8948E-01 -1.8882E-01  1.7017E-01  0.0000E+00
             1.6610E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1920
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.1183E-04 -1.3298E-02 -6.3207E-03  9.0756E-03 -5.3288E-04
 SE:             2.9459E-02  2.5675E-02  1.8386E-02  2.3212E-02  2.9691E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8614E-01  6.0449E-01  7.3101E-01  6.9581E-01  7.2691E-02

 ETASHRINKSD(%)  1.3094E+00  1.3987E+01  3.8404E+01  2.2235E+01  9.9005E+01
 ETASHRINKVR(%)  2.6016E+00  2.6017E+01  6.2060E+01  3.9527E+01  9.9990E+01
 EBVSHRINKSD(%)  1.5754E+00  1.3771E+01  3.8045E+01  2.2013E+01  9.9016E+01
 EBVSHRINKVR(%)  3.1259E+00  2.5645E+01  6.1615E+01  3.9180E+01  9.9990E+01
 RELATIVEINF(%)  9.5408E+01  8.7523E+00  2.1016E+01  1.2377E+01  1.1974E-03
 EPSSHRINKSD(%)  4.1431E+01
 EPSSHRINKVR(%)  6.5697E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1505.3271486279136     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -770.17632206417545     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.56
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1505.327       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.55E+00  1.40E-01  5.87E-01  7.07E-01  9.18E-01  7.30E-01  2.10E+00  1.20E+00  1.00E-02  1.86E+00
 


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
+        1.22E+03
 
 TH 2
+       -2.55E+01  4.66E+04
 
 TH 3
+       -1.31E+02  1.50E+05  4.92E+05
 
 TH 4
+       -4.91E+01  1.52E+05 -1.27E+03  5.00E+05
 
 TH 5
+       -5.52E+01  2.21E+05 -1.26E+03  7.31E+05  2.23E+03
 
 TH 6
+        5.54E-01 -6.24E+00  7.43E+00 -1.21E+01  3.40E-01  2.22E+02
 
 TH 7
+        6.93E-01  4.95E+00 -3.67E+01 -6.41E+00  1.96E+01 -1.55E-01  2.10E+02
 
 TH 8
+        2.50E+01  2.16E+04  1.88E+03  7.11E+04  1.04E+05  2.27E+01  9.67E+00  1.01E+04
 
 TH 9
+        2.77E+00 -1.73E+01  6.11E+01  3.72E+01 -6.35E+00 -1.38E-01  1.73E+01  8.32E+00  4.80E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.46E+01 -2.90E+04 -2.59E+03 -9.56E+04 -1.40E+05 -2.61E+01  1.67E+00 -1.36E+04  5.78E+00  0.00E+00  1.84E+04
 
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
 #CPUT: Total CPU Time in Seconds,       32.407
Stop Time:
Wed Sep 29 13:01:19 CDT 2021
