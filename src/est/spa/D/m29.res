Wed Sep 29 19:55:49 CDT 2021
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
$DATA ../../../../data/spa/D/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   10229.9468848100        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0112E+02  2.1115E+01 -8.4636E+01 -1.2150E+01  2.4358E+02 -1.1314E+03 -6.1112E+02 -1.6683E+01 -9.4101E+02 -3.5629E+02
            -2.0476E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -685.533161309008        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.3368E+00  1.2218E+00  1.0461E+00  1.5448E+00  1.0035E+00  1.6404E+00  1.2410E+00  9.8422E-01  1.3697E+00  1.1289E+00
             1.4468E+01
 PARAMETER:  3.9031E-01  3.0029E-01  1.4502E-01  5.3490E-01  1.0347E-01  5.9493E-01  3.1590E-01  8.4090E-02  4.1459E-01  2.2128E-01
             2.7719E+00
 GRADIENT:  -2.7055E+01  1.3648E+01 -7.0896E-02  1.8988E+01 -1.7084E+01  2.3522E+01  1.9106E+00  3.4083E+00  1.6044E+01  6.9862E+00
             2.1821E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -700.368011124881        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.3644E+00  9.3751E-01  2.4836E+00  2.0440E+00  2.2962E+00  1.5366E+00  3.4773E+00  5.0529E-01  1.5300E+00  4.2582E+00
             1.2755E+01
 PARAMETER:  4.1074E-01  3.5470E-02  1.0097E+00  8.1491E-01  9.3126E-01  5.2955E-01  1.3463E+00 -5.8263E-01  5.2530E-01  1.5488E+00
             2.6460E+00
 GRADIENT:   4.6487E+00  2.0047E+01  5.3230E+00  5.7757E+01 -1.2811E+01 -2.5030E+01  9.7650E+00  7.0987E-02  2.6425E+01  9.6824E+00
             1.6015E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -731.906835154091        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.1554E+00  1.0306E+00  1.8392E+00  1.3541E+00  2.6894E+00  1.8448E+00  1.4370E+00  1.7121E-01  1.4804E+00  5.0392E+00
             1.0148E+01
 PARAMETER:  2.4445E-01  1.3015E-01  7.0935E-01  4.0314E-01  1.0893E+00  7.1239E-01  4.6256E-01 -1.6649E+00  4.9230E-01  1.7172E+00
             2.4173E+00
 GRADIENT:  -1.5572E+01 -1.1445E+01  6.3433E+00 -3.5421E+01 -1.3729E+01  5.2407E+01  2.0412E+00  6.5460E-03  7.4397E+00  6.3057E+00
             8.7804E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -743.680882013159        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.1524E+00  8.4964E-01  1.3723E+00  1.5070E+00  4.5541E+00  1.4709E+00  2.1822E+00  2.0161E-01  1.0224E+00  6.8917E+00
             9.3044E+00
 PARAMETER:  2.4181E-01 -6.2948E-02  4.1647E-01  5.1012E-01  1.6160E+00  4.8587E-01  8.8035E-01 -1.5014E+00  1.2218E-01  2.0303E+00
             2.3305E+00
 GRADIENT:  -2.4238E-01  1.9566E+01  5.6006E+00  3.5545E+01 -3.8307E+00 -1.6868E+01  8.0007E-01  2.6099E-03 -4.6165E+00  2.7237E-01
            -2.1493E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -757.440744301168        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.1032E+00  3.7676E-01  8.0376E-01  1.5736E+00  5.9210E+00  1.5042E+00  8.6029E-01  1.2371E-02  8.9398E-01  7.5023E+00
             9.5578E+00
 PARAMETER:  1.9822E-01 -8.7615E-01 -1.1845E-01  5.5338E-01  1.8785E+00  5.0823E-01 -5.0485E-02 -4.2924E+00 -1.2073E-02  2.1152E+00
             2.3574E+00
 GRADIENT:   3.0313E-01  8.0916E+00  7.0141E-01  6.8985E+00 -4.9487E+00  2.5021E+00  1.9192E-01 -2.6097E-04 -1.1567E+00  7.0559E+00
             1.2843E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -770.703075588625        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      473
 NPARAMETR:  9.4610E-01  3.0118E-02  2.5463E-01  1.2832E+00  8.7043E+00  1.2883E+00  4.6374E-02  1.0000E-02  3.4328E-01  2.6749E+00
             9.2300E+00
 PARAMETER:  4.4594E-02 -3.4026E+00 -1.2679E+00  3.4935E-01  2.2638E+00  3.5335E-01 -2.9710E+00 -1.0396E+01 -9.6921E-01  1.0839E+00
             2.3225E+00
 GRADIENT:   4.1163E+01  1.5178E+00 -6.0954E+01  2.2407E+02  6.5797E+00 -6.6178E+01 -6.4296E-06  0.0000E+00 -1.5102E+01 -6.7650E+00
            -1.2502E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -797.244873628769        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      572
 NPARAMETR:  7.9135E-01  1.0000E-02  1.3311E-01  9.4438E-01  8.5843E+00  1.3619E+00  5.4800E-02  1.0000E-02  3.8189E-01  2.6907E+00
             9.6677E+00
 PARAMETER: -1.3402E-01 -4.9831E+00 -1.9166E+00  4.2771E-02  2.2499E+00  4.0888E-01 -2.8041E+00 -1.0396E+01 -8.6262E-01  1.0898E+00
             2.3688E+00
 GRADIENT:   4.4638E+00  0.0000E+00 -4.3837E+01  1.5199E+02  1.1236E+00 -2.7042E+01  6.4124E-07  0.0000E+00 -8.6564E+00 -9.2708E-02
            -4.8345E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -805.474013156706        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      672
 NPARAMETR:  6.4799E-01  1.3110E-02  6.5771E-02  6.3167E-01  8.5631E+00  1.3164E+00  5.2611E-02  1.0000E-02  4.3506E-01  2.6902E+00
             9.9826E+00
 PARAMETER: -3.3388E-01 -4.2344E+00 -2.6216E+00 -3.5939E-01  2.2475E+00  3.7492E-01 -2.8448E+00 -1.0396E+01 -7.3227E-01  1.0896E+00
             2.4008E+00
 GRADIENT:  -2.0353E+01 -3.4166E-01 -6.2950E+01  1.1615E+02  1.1573E+01 -2.8034E+01 -1.8621E-06  0.0000E+00 -8.0530E-01 -5.1603E+00
             3.8369E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -810.257912726251        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      856
 NPARAMETR:  6.4550E-01  2.1130E-02  6.3235E-02  5.9089E-01  8.5968E+00  1.3642E+00  1.0993E-01  1.0000E-02  4.4627E-01  2.6883E+00
             9.8918E+00
 PARAMETER: -3.3772E-01 -3.7571E+00 -2.6609E+00 -4.2613E-01  2.2514E+00  4.1056E-01 -2.1079E+00 -1.0396E+01 -7.0684E-01  1.0889E+00
             2.3917E+00
 GRADIENT:  -6.4739E+00 -4.4361E-01 -1.7877E+01  3.5734E+01  7.2242E+00 -1.1568E+01 -8.3460E-06  0.0000E+00 -1.3788E+00 -2.9565E+00
             4.4360E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -811.793128671381        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      957
 NPARAMETR:  6.2718E-01  2.0464E-02  6.3577E-02  5.8609E-01  8.7322E+00  1.4020E+00  1.2840E-01  1.0000E-02  4.4970E-01  2.6652E+00
             9.3025E+00
 PARAMETER: -3.6653E-01 -3.7891E+00 -2.6555E+00 -4.3429E-01  2.2670E+00  4.3792E-01 -1.9526E+00 -1.0396E+01 -6.9918E-01  1.0803E+00
             2.3303E+00
 GRADIENT:   1.9641E+00 -5.5392E-01  2.1288E+01  3.7777E+01  1.1021E+01  5.8647E-01 -1.3758E-05  0.0000E+00 -4.3589E+00 -3.6474E+00
             2.4297E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -812.337994303528        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  6.4274E-01  2.0686E-02  6.4529E-02  5.8325E-01  8.6206E+00  1.4223E+00  2.2447E+00  1.0000E-02  4.5206E-01  2.6781E+00
             9.2023E+00
 PARAMETER: -3.4202E-01 -3.7783E+00 -2.6406E+00 -4.3914E-01  2.2542E+00  4.5230E-01  9.0856E-01 -1.0396E+01 -6.9395E-01  1.0851E+00
             2.3195E+00
 GRADIENT:   1.1155E+00 -6.3609E-01  6.4074E+00 -4.6770E+00  8.6154E+00  7.2678E-01  5.5169E-04  0.0000E+00 -6.1605E+00 -3.5921E+00
            -3.2304E+00

0ITERATION NO.:   57    OBJECTIVE VALUE:  -812.458864472546        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1201
 NPARAMETR:  6.3978E-01  2.0811E-02  6.4581E-02  5.8136E-01  8.5342E+00  1.4157E+00  6.1046E-01  1.0000E-02  4.5499E-01  2.6873E+00
             9.2017E+00
 PARAMETER: -3.4782E-01 -3.7691E+00 -2.6411E+00 -4.4204E-01  2.2423E+00  4.4687E-01 -3.8964E-01 -1.0396E+01 -6.8819E-01  1.0894E+00
             2.3178E+00
 GRADIENT:  -3.3553E+00  1.3017E+01 -2.6963E+01  1.0640E+02 -7.5996E+01 -6.7110E-01  1.9509E-04  0.0000E+00 -4.3345E+01  4.4179E+01
            -2.4677E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1201
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.5278E-03  3.3468E-04  7.8581E-05 -1.2747E-02 -1.4653E-02
 SE:             2.8717E-02  1.8309E-04  2.5815E-04  1.6716E-02  7.1881E-03
 N:                     100         100         100         100         100

 P VAL.:         8.2018E-01  6.7566E-02  7.6082E-01  4.4570E-01  4.1501E-02

 ETASHRINKSD(%)  3.7935E+00  9.9387E+01  9.9135E+01  4.4000E+01  7.5919E+01
 ETASHRINKVR(%)  7.4431E+00  9.9996E+01  9.9993E+01  6.8640E+01  9.4201E+01
 EBVSHRINKSD(%)  4.0752E+00  9.9372E+01  9.9213E+01  4.7899E+01  8.1096E+01
 EBVSHRINKVR(%)  7.9843E+00  9.9996E+01  9.9994E+01  7.2855E+01  9.6426E+01
 RELATIVEINF(%)  9.9546E+00  4.2764E-04  1.0478E-04  4.6769E-01  5.8801E-01
 EPSSHRINKSD(%)  1.0459E+01
 EPSSHRINKVR(%)  1.9825E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -812.45886447254611     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -77.308037908807933     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -812.459       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.39E-01  2.09E-02  6.45E-02  5.82E-01  8.52E+00  1.41E+00  6.13E-01  1.00E-02  4.55E-01  2.69E+00  9.19E+00
 


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
+        1.24E+03
 
 TH 2
+       -4.53E+02  2.06E+05
 
 TH 3
+       -2.46E+03 -2.66E+05  2.78E+05
 
 TH 4
+       -2.67E+02 -6.59E+02 -7.34E+04  2.23E+04
 
 TH 5
+        8.76E-01 -7.10E+02  1.06E+03 -2.20E+02  1.30E+01
 
 TH 6
+        5.95E+00 -1.63E+01  1.82E+02 -5.75E+01  2.44E-01  8.24E+01
 
 TH 7
+        3.51E-03  5.86E-01 -3.88E-01  3.52E-02 -6.66E-03  9.78E-04  6.58E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.84E+04  1.68E+02  2.40E+04  8.28E+01  2.57E+00 -6.45E+03 -3.51E-02  0.00E+00  1.30E+04
 
 TH10
+       -4.67E+00  5.58E+03 -5.02E+03 -2.65E+01 -2.20E+01 -7.22E-01  1.04E-02  0.00E+00  6.29E+00  1.52E+02
 
 TH11
+       -1.54E+01 -7.54E+02  8.22E+02 -1.46E+01  2.81E+00  1.71E+00 -9.28E-04  0.00E+00  8.25E+00 -2.05E+01  7.81E+00
 
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
 #CPUT: Total CPU Time in Seconds,       24.997
Stop Time:
Wed Sep 29 19:56:16 CDT 2021
