Sun Oct 24 01:44:17 CDT 2021
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
$DATA ../../../../data/SD4/B/dat78.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1664.11147683765        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5393E+02 -6.6903E+01 -3.6571E+01 -4.7751E+01  1.2368E+01  3.2457E+01 -1.5892E+01  1.5554E+01 -7.4549E+00  1.2363E+01
             1.3613E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1676.34579065550        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7559E-01  1.1531E+00  1.2479E+00  9.9185E-01  1.1552E+00  9.9436E-01  1.1267E+00  8.9381E-01  1.0525E+00  9.4383E-01
             9.4505E-01
 PARAMETER:  7.5286E-02  2.4247E-01  3.2144E-01  9.1814E-02  2.4425E-01  9.4342E-02  2.1930E-01 -1.2262E-02  1.5116E-01  4.2187E-02
             4.3488E-02
 GRADIENT:  -3.3571E-02  1.7890E+00  1.3702E+01 -5.3021E+00 -2.4297E-01 -9.2579E-01 -3.0517E+00 -1.9392E+00  7.4105E-01 -2.3631E+01
            -2.1458E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1677.55664926321        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.7201E-01  1.1901E+00  1.2967E+00  9.7928E-01  1.2388E+00  9.9037E-01  1.1499E+00  7.5707E-01  1.1042E+00  1.1305E+00
             9.8514E-01
 PARAMETER:  7.1608E-02  2.7406E-01  3.5980E-01  7.9067E-02  3.1416E-01  9.0323E-02  2.3968E-01 -1.7830E-01  1.9915E-01  2.2265E-01
             8.5027E-02
 GRADIENT:  -9.0081E+00 -1.4743E-01  2.2110E+00  9.6198E+00  2.2877E+01 -2.5440E+00  4.2512E+00 -2.5772E+00  8.3808E+00 -4.5428E+00
            -1.2086E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1678.75891344262        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.7799E-01  1.2700E+00  1.0501E+00  9.1970E-01  1.1474E+00  1.0008E+00  1.1077E+00  5.6626E-01  1.0693E+00  1.0668E+00
             9.7999E-01
 PARAMETER:  7.7743E-02  3.3899E-01  1.4887E-01  1.6289E-02  2.3753E-01  1.0077E-01  2.0229E-01 -4.6870E-01  1.6698E-01  1.6467E-01
             7.9785E-02
 GRADIENT:   1.7633E+00  4.4740E+00  4.3551E-01  6.0651E+00 -2.8835E+00  1.3525E+00 -9.8831E-03  5.5714E-02  2.2015E-01  4.4149E-01
             9.0436E-03

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1678.92965886080        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.7810E-01  1.5179E+00  9.3117E-01  7.5606E-01  1.2319E+00  9.9657E-01  9.7514E-01  5.4308E-01  1.2174E+00  1.1032E+00
             9.8077E-01
 PARAMETER:  7.7853E-02  5.1734E-01  2.8687E-02 -1.7963E-01  3.0854E-01  9.6567E-02  7.4823E-02 -5.1050E-01  2.9672E-01  1.9824E-01
             8.0581E-02
 GRADIENT:  -1.7164E-01  4.6572E+00  3.2838E+00  1.3031E+00 -4.2210E+00 -6.4957E-01 -3.1296E-01 -5.2960E-02 -8.1972E-01 -4.5106E-01
            -6.5132E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1678.93676677221        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.7844E-01  1.5955E+00  8.7416E-01  7.0700E-01  1.2516E+00  9.9753E-01  9.4466E-01  4.9227E-01  1.2711E+00  1.1112E+00
             9.8127E-01
 PARAMETER:  7.8202E-02  5.6717E-01 -3.4491E-02 -2.4673E-01  3.2442E-01  9.7523E-02  4.3065E-02 -6.0873E-01  3.3986E-01  2.0541E-01
             8.1096E-02
 GRADIENT:  -1.3075E-01  7.2927E+00  2.8767E+00  3.1335E+00 -4.4343E+00 -3.9929E-01 -4.6606E-01  8.7593E-04 -7.9536E-01 -3.6518E-01
            -5.4354E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1678.98826423013        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  9.7865E-01  1.6088E+00  8.5118E-01  6.9293E-01  1.2584E+00  9.9860E-01  9.4138E-01  4.3878E-01  1.2886E+00  1.1148E+00
             9.8181E-01
 PARAMETER:  7.8423E-02  5.7550E-01 -6.1130E-02 -2.6683E-01  3.2980E-01  9.8594E-02  3.9591E-02 -7.2375E-01  3.5356E-01  2.0868E-01
             8.1646E-02
 GRADIENT:   1.6749E-01 -1.9844E+00  1.7269E-01  1.0858E+00  1.0832E+00 -1.6054E-02  1.5243E-02  6.1487E-02  1.5736E-01  5.2741E-02
             1.4757E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1679.00175494349        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  9.7829E-01  1.6408E+00  8.0128E-01  6.7411E-01  1.2512E+00  9.9892E-01  9.3688E-01  2.2976E-01  1.3012E+00  1.1064E+00
             9.8187E-01
 PARAMETER:  7.8050E-02  5.9518E-01 -1.2154E-01 -2.9436E-01  3.2412E-01  9.8917E-02  3.4803E-02 -1.3707E+00  3.6327E-01  2.0115E-01
             8.1703E-02
 GRADIENT:  -1.2404E+00  2.5515E+00 -1.2374E-01  3.7414E+00 -2.9576E-01 -1.0538E-02  3.5804E-01  1.2933E-02 -3.0007E-02  2.4429E-01
             1.4309E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1679.02084460521        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  9.7920E-01  1.6446E+00  7.9264E-01  6.6666E-01  1.2511E+00  9.9905E-01  9.3317E-01  1.6116E-01  1.3059E+00  1.1037E+00
             9.8162E-01
 PARAMETER:  7.8979E-02  5.9751E-01 -1.3239E-01 -3.0548E-01  3.2398E-01  9.9050E-02  3.0836E-02 -1.7254E+00  3.6691E-01  1.9868E-01
             8.1448E-02
 GRADIENT:   8.1021E-01 -3.1146E+00  3.8749E-01 -3.1452E-01 -6.7294E-01  5.5932E-02  4.6759E-02  3.5804E-03 -3.1980E-01 -8.9451E-03
            -8.3099E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1679.02284337196        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1607
 NPARAMETR:  9.7956E-01  1.6438E+00  7.9144E-01  6.6638E-01  1.2510E+00  9.9914E-01  9.3268E-01  8.5847E-02  1.3092E+00  1.1039E+00
             9.8184E-01
 PARAMETER:  7.9352E-02  5.9702E-01 -1.3390E-01 -3.0590E-01  3.2394E-01  9.9135E-02  3.0312E-02 -2.3552E+00  3.6938E-01  1.9887E-01
             8.1674E-02
 GRADIENT:   1.6423E+00 -4.3909E+00  2.9695E-01 -6.1061E-01 -1.6348E-01  8.7667E-02 -1.2820E-02  2.6687E-04 -1.0055E-02  1.5419E-02
            -7.3786E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1679.02315068375        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1787
 NPARAMETR:  9.7957E-01  1.6438E+00  7.9050E-01  6.6639E-01  1.2504E+00  9.9914E-01  9.3283E-01  4.5998E-02  1.3089E+00  1.1036E+00
             9.8190E-01
 PARAMETER:  7.9357E-02  5.9703E-01 -1.3509E-01 -3.0587E-01  3.2344E-01  9.9141E-02  3.0471E-02 -2.9792E+00  3.6922E-01  1.9857E-01
             8.1738E-02
 GRADIENT:   1.6423E+00 -4.2529E+00  3.2385E-01 -5.6833E-01 -2.8166E-01  8.7410E-02 -7.6270E-03  7.3856E-05 -5.5277E-03  3.7845E-02
            -5.9950E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1679.02343954129        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1967
 NPARAMETR:  9.7957E-01  1.6439E+00  7.8968E-01  6.6653E-01  1.2501E+00  9.9914E-01  9.3293E-01  3.8918E-02  1.3089E+00  1.1032E+00
             9.8192E-01
 PARAMETER:  7.9357E-02  5.9704E-01 -1.3612E-01 -3.0568E-01  3.2324E-01  9.9136E-02  3.0577E-02 -3.1463E+00  3.6921E-01  1.9819E-01
             8.1751E-02
 GRADIENT:   1.6263E+00 -4.2096E+00  1.4341E-01 -3.0699E-01 -4.3043E-02  8.2438E-02  9.5080E-04  1.1618E-04  5.2069E-02  3.6308E-02
            -1.8477E-02

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1679.02352600403        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     2063
 NPARAMETR:  9.7957E-01  1.6437E+00  7.8921E-01  6.6664E-01  1.2501E+00  9.9913E-01  9.3298E-01  1.0000E-02  1.3090E+00  1.1030E+00
             9.8196E-01
 PARAMETER:  7.9355E-02  5.9697E-01 -1.3672E-01 -3.0550E-01  3.2320E-01  9.9133E-02  3.0626E-02 -5.4581E+00  3.6928E-01  1.9801E-01
             8.1792E-02
 GRADIENT:   1.6137E+00 -4.3051E+00 -2.8144E-03 -1.4171E-01  2.3651E-01  7.8612E-02 -2.3253E-03  0.0000E+00  1.0198E-01  2.7762E-02
             1.4710E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2063
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2738E-04 -2.3819E-02 -3.0090E-04  1.8783E-02 -3.6834E-02
 SE:             2.9795E-02  2.3206E-02  1.1085E-04  2.2201E-02  2.2627E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9659E-01  3.0471E-01  6.6387E-03  3.9753E-01  1.0356E-01

 ETASHRINKSD(%)  1.8348E-01  2.2257E+01  9.9629E+01  2.5623E+01  2.4196E+01
 ETASHRINKVR(%)  3.6663E-01  3.9560E+01  9.9999E+01  4.4681E+01  4.2537E+01
 EBVSHRINKSD(%)  4.4337E-01  2.0848E+01  9.9689E+01  2.8135E+01  2.1860E+01
 EBVSHRINKVR(%)  8.8477E-01  3.7349E+01  9.9999E+01  4.8354E+01  3.8941E+01
 RELATIVEINF(%)  9.8938E+01  3.4387E+00  1.2165E-04  2.8300E+00  1.3582E+01
 EPSSHRINKSD(%)  4.3248E+01
 EPSSHRINKVR(%)  6.7792E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1679.0235260040279     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -943.87269944028969     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1679.024       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.64E+00  7.89E-01  6.67E-01  1.25E+00  9.99E-01  9.33E-01  1.00E-02  1.31E+00  1.10E+00  9.82E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       64.176
Stop Time:
Sun Oct 24 01:44:30 CDT 2021
