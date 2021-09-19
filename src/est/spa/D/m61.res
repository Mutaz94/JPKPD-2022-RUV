Sat Sep 18 15:29:48 CDT 2021
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
$DATA ../../../../data/spa/D/dat61.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12984.1078553514        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -9.2910E+01  2.5270E+02 -3.3398E+01  1.4742E+02  3.8635E+02 -2.5551E+03 -6.6556E+02 -2.6411E+01 -1.1443E+03 -8.8675E+02
            -2.4082E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -624.534120205788        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.5935E+00  1.0175E+00  8.5121E-01  1.7019E+00  1.3088E+00  2.1600E+00  1.1871E+00  9.6331E-01  1.1868E+00  1.2784E+00
             1.4315E+01
 PARAMETER:  5.6592E-01  1.1738E-01 -6.1097E-02  6.3172E-01  3.6912E-01  8.7011E-01  2.7153E-01  6.2617E-02  2.7130E-01  3.4558E-01
             2.7613E+00
 GRADIENT:   5.3695E+01  1.2358E+01 -5.9070E+00  8.6766E+00 -5.6963E+00  3.0968E+01  5.9750E-01  6.2588E+00  7.3644E+00  2.6874E+00
             8.2470E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -640.422868239447        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.5101E+00  7.2572E-01  1.0519E+00  2.0635E+00  2.7320E+00  2.1513E+00  1.4519E+00  5.1929E-01  2.1592E+00  7.6946E+00
             1.1490E+01
 PARAMETER:  5.1216E-01 -2.2059E-01  1.5059E-01  8.2440E-01  1.1051E+00  8.6608E-01  4.7289E-01 -5.5530E-01  8.6975E-01  2.1405E+00
             2.5415E+00
 GRADIENT:   4.5780E+01  1.0621E+01 -1.5021E+00  3.9517E+01 -2.8295E+01 -6.5980E+01  4.1608E+00 -8.3738E-01  2.8785E+01  2.6288E+01
             4.3985E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -668.485385873699        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.3281E+00  8.4061E-01  8.4780E-01  1.6879E+00  5.4866E+00  2.6705E+00  7.0373E-01  7.8338E-02  2.0695E+00  1.0353E+01
             1.0859E+01
 PARAMETER:  3.8378E-01 -7.3626E-02 -6.5116E-02  6.2350E-01  1.8023E+00  1.0823E+00 -2.5136E-01 -2.4467E+00  8.2732E-01  2.4372E+00
             2.4850E+00
 GRADIENT:   1.4317E+00  1.4706E+01  5.6428E+00  9.6780E+00 -4.2329E+00  2.5049E+01  2.1157E+00 -3.2214E-03 -2.3088E+00 -3.8200E+00
             6.3855E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -732.150970706745        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.2212E-01  1.4656E-01  1.5913E-01  1.1046E+00  1.8140E+02  1.8812E+00  1.0000E-02  1.0000E-02  1.4007E+00  2.1406E+01
             7.7648E+00
 PARAMETER:  1.8916E-02 -1.8203E+00 -1.7381E+00  1.9945E-01  5.3007E+00  7.3189E-01 -5.2689E+00 -6.2476E+00  4.3694E-01  3.1637E+00
             2.1496E+00
 GRADIENT:   2.9216E+01  4.2826E+01 -5.3166E+01  8.0958E+01 -3.9746E-01 -7.2110E+01  0.0000E+00  0.0000E+00  4.3255E+01  1.8557E-01
            -1.7741E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -803.858759905592        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  5.4279E-01  2.2243E-02  3.5807E-02  3.8064E-01  6.8297E+02  1.9227E+00  1.0000E-02  1.0000E-02  4.9314E-01  2.2557E+01
             9.2213E+00
 PARAMETER: -5.1104E-01 -3.7057E+00 -3.2296E+00 -8.6591E-01  6.6264E+00  7.5373E-01 -1.4759E+01 -1.4471E+01 -6.0697E-01  3.2160E+00
             2.3215E+00
 GRADIENT:   5.1113E+00  1.7799E+00 -8.0014E+00  2.5417E+01  7.9332E-06  2.0411E+00  0.0000E+00  0.0000E+00 -7.2773E+00  2.2401E-04
            -1.7606E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -805.718689390485        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  4.9841E-01  1.4999E-02  2.8703E-02  3.1489E-01  9.7997E+02  1.9544E+00  1.0000E-02  1.0000E-02  5.8073E-01  2.3000E+01
             9.3447E+00
 PARAMETER: -5.9633E-01 -4.0998E+00 -3.4507E+00 -1.0555E+00  6.9875E+00  7.7009E-01 -1.7406E+01 -1.6009E+01 -4.4347E-01  3.2355E+00
             2.3348E+00
 GRADIENT:  -1.3778E+00  6.2713E-01  2.0547E+01 -2.6586E+01 -2.5735E-04  9.9103E+00  0.0000E+00  0.0000E+00 -4.5331E+00 -1.8519E-06
             1.3786E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -806.704277391955        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  4.8793E-01  1.3206E-02  2.7215E-02  3.0523E-01  1.2817E+03  1.9053E+00  1.0000E-02  1.0000E-02  6.8341E-01  2.4084E+01
             9.0200E+00
 PARAMETER: -6.1759E-01 -4.2271E+00 -3.5040E+00 -1.0867E+00  7.2560E+00  7.4462E-01 -1.8071E+01 -1.6198E+01 -2.8067E-01  3.2815E+00
             2.2994E+00
 GRADIENT:  -2.9577E-02  2.0796E-01 -4.6930E-02 -3.1719E-01  3.0170E-04 -6.5801E-02  0.0000E+00  0.0000E+00  9.0823E-02 -7.7149E-06
            -6.8076E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -806.716352150667        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      891
 NPARAMETR:  4.8884E-01  1.1502E-02  2.7419E-02  3.0705E-01  1.2490E+03  1.9070E+00  1.0000E-02  1.0000E-02  6.8215E-01  2.3227E+01
             9.0243E+00
 PARAMETER: -6.1573E-01 -4.3653E+00 -3.4965E+00 -1.0807E+00  7.2301E+00  7.4553E-01 -1.8370E+01 -1.6115E+01 -2.8250E-01  3.2453E+00
             2.2999E+00
 GRADIENT:   1.7134E-02  2.1585E-03  4.0355E-02 -7.1063E-02  4.5951E-04 -4.3446E-02  0.0000E+00  0.0000E+00  5.6451E-03 -1.0430E-05
             1.6309E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -806.724949154580        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  4.8771E-01  1.1301E-02  2.7378E-02  3.0695E-01  5.1127E+01  1.9051E+00  1.0000E-02  1.0000E-02  6.8192E-01  2.3201E+01
             9.0175E+00
 PARAMETER: -6.1804E-01 -4.3829E+00 -3.4980E+00 -1.0811E+00  4.0343E+00  7.4454E-01 -1.8370E+01 -1.6115E+01 -2.8284E-01  3.2442E+00
             2.2992E+00
 GRADIENT:   3.7246E+00  2.0231E-02  6.0202E+00  4.0889E+00  7.0945E-03  1.5195E+00  0.0000E+00  0.0000E+00 -7.2729E-03  1.0795E-02
             1.5302E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -806.728584803089        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1166
 NPARAMETR:  4.8927E-01  1.1191E-02  2.7392E-02  3.0695E-01  4.5425E+01  1.9079E+00  1.0000E-02  1.0000E-02  6.8189E-01  2.2310E+01
             9.0243E+00
 PARAMETER: -6.1484E-01 -4.3927E+00 -3.4975E+00 -1.0811E+00  3.9161E+00  7.4600E-01 -1.8370E+01 -1.6115E+01 -2.8289E-01  3.2051E+00
             2.2999E+00
 GRADIENT:   2.4110E-01  1.9511E-03 -7.5561E-02 -1.0485E-01  4.1892E-03  6.0126E-02  0.0000E+00  0.0000E+00  8.0681E-03  1.6520E-02
             2.1189E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -806.745070605696        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1346
 NPARAMETR:  4.8980E-01  1.0328E-02  2.7403E-02  3.0721E-01  1.9379E+01  1.9081E+00  1.0000E-02  1.0000E-02  6.8038E-01  1.3151E+01
             9.0298E+00
 PARAMETER: -6.1375E-01 -4.4729E+00 -3.4971E+00 -1.0802E+00  3.0642E+00  7.4610E-01 -1.8370E+01 -1.6115E+01 -2.8511E-01  2.6765E+00
             2.3005E+00
 GRADIENT:   2.2642E-01  1.1414E-02 -4.4608E-02 -1.7221E-01 -9.5105E-02 -6.0713E-02  0.0000E+00  0.0000E+00 -5.5213E-02  1.2777E-01
             2.8758E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -806.748513910719        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1524
 NPARAMETR:  4.8984E-01  1.0424E-02  2.7407E-02  3.0732E-01  1.7114E+01  1.9106E+00  1.0000E-02  1.0000E-02  6.8187E-01  1.1481E+01
             9.0158E+00
 PARAMETER: -6.1368E-01 -4.4636E+00 -3.4970E+00 -1.0799E+00  2.9399E+00  7.4741E-01 -1.8370E+01 -1.6115E+01 -2.8292E-01  2.5407E+00
             2.2990E+00
 GRADIENT:   2.6995E-01  1.0335E-01 -3.9279E-01  2.2954E-01 -5.7799E-01  2.8426E-01  0.0000E+00  0.0000E+00 -5.9126E-02  5.1369E-01
            -8.2448E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -806.749716963344        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1701
 NPARAMETR:  4.8981E-01  1.0177E-02  2.7416E-02  3.0735E-01  1.6675E+01  1.9084E+00  1.0000E-02  1.0000E-02  6.8216E-01  1.1181E+01
             9.0227E+00
 PARAMETER: -6.1373E-01 -4.4877E+00 -3.4966E+00 -1.0798E+00  2.9139E+00  7.4627E-01 -1.8370E+01 -1.6115E+01 -2.8249E-01  2.5142E+00
             2.2997E+00
 GRADIENT:   1.5630E-01  1.0237E-01  3.8806E-02 -2.2923E-01 -6.9894E-01 -1.3595E-01  0.0000E+00  0.0000E+00  3.3533E-02  5.9579E-01
            -2.7446E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -806.749947002523        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1871
 NPARAMETR:  4.8955E-01  1.0097E-02  2.7515E-02  3.0706E-01  1.6602E+01  1.9098E+00  1.0000E-02  1.0000E-02  6.8174E-01  1.1060E+01
             9.0449E+00
 PARAMETER: -6.1366E-01 -4.4911E+00 -3.4965E+00 -1.0796E+00  2.9066E+00  7.4624E-01 -1.8370E+01 -1.6115E+01 -2.8340E-01  2.5058E+00
             2.2999E+00
 GRADIENT:   2.2631E+03  1.5465E+02 -3.9526E+02  1.2834E+03 -2.3893E+02 -9.3063E+02  0.0000E+00  0.0000E+00 -2.4496E+03  2.7762E+02
            -6.0410E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1871
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4402E-03  1.5486E-06  1.5263E-04 -1.9353E-02 -9.9374E-04
 SE:             2.9221E-02  8.3776E-07  2.5594E-04  2.1289E-02  7.5617E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6069E-01  6.4535E-02  5.5093E-01  3.6334E-01  1.8879E-01

 ETASHRINKSD(%)  2.1062E+00  9.9997E+01  9.9143E+01  2.8678E+01  9.7467E+01
 ETASHRINKVR(%)  4.1680E+00  1.0000E+02  9.9993E+01  4.9131E+01  9.9936E+01
 EBVSHRINKSD(%)  2.0450E+00  9.9996E+01  9.9182E+01  3.0087E+01  9.7746E+01
 EBVSHRINKVR(%)  4.0481E+00  1.0000E+02  9.9993E+01  5.1121E+01  9.9949E+01
 RELATIVEINF(%)  3.6633E+00  1.7535E-08  5.2830E-05  3.0881E-01  3.3922E-03
 EPSSHRINKSD(%)  1.3143E+01
 EPSSHRINKVR(%)  2.4558E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -806.74994700252341     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -71.599120438785235     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -806.750       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.90E-01  1.01E-02  2.74E-02  3.07E-01  1.66E+01  1.91E+00  1.00E-02  1.00E-02  6.82E-01  1.11E+01  9.02E+00
 


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
+        3.84E+06
 
 TH 2
+       -5.61E+03  1.67E+08
 
 TH 3
+       -4.04E+03  7.56E+05  3.83E+07
 
 TH 4
+       -1.54E+04 -1.00E+05 -4.76E+04  3.15E+06
 
 TH 5
+        4.82E+00  2.07E+01 -7.17E+02  9.51E+01  1.50E+02
 
 TH 6
+        6.09E+01  3.62E+02 -1.46E+02  3.13E+01 -2.97E-01  1.71E+05
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.06E+03  2.68E+04 -1.20E+04  3.62E+03 -2.53E+01 -8.63E+02  0.00E+00  0.00E+00  9.30E+06
 
 TH10
+       -8.23E+00 -3.85E+01  1.24E+03 -1.64E+02 -7.17E-01  5.26E-01  0.00E+00  0.00E+00  4.38E+01  4.51E+02
 
 TH11
+       -1.78E+00  6.27E+01 -1.37E+03  1.94E+02  3.09E-01  2.19E-01  0.00E+00  0.00E+00 -5.26E+01 -5.29E-01  8.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.238
Stop Time:
Sat Sep 18 15:30:27 CDT 2021
