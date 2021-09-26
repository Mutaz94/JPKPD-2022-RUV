Sat Sep 25 02:26:18 CDT 2021
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
$DATA ../../../../data/int/SL3/dat57.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      985
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      885
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -507.904438893077        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.3426E+01 -6.4334E+01  6.6842E+01  4.8628E+01  1.7987E+02  3.6088E+01 -1.2700E+02 -1.9327E+02 -9.4078E+01 -5.1014E+01
            -6.2909E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2764.06821094398        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0738E+00  1.3205E+00  1.0494E+00  8.8570E-01  1.1172E+00  7.6619E-01  1.0969E+00  9.4854E-01  8.1718E-01  1.0780E+00
             2.7948E+00
 PARAMETER:  1.7124E-01  3.7802E-01  1.4821E-01 -2.1380E-02  2.1080E-01 -1.6633E-01  1.9249E-01  4.7168E-02 -1.0190E-01  1.7512E-01
             1.1278E+00
 GRADIENT:   6.6899E+01  6.9041E+00 -1.0544E+01  3.8588E-01 -7.1266E+00 -5.9681E+01  1.8183E+01  3.3757E+00 -3.9464E+00 -1.9974E+01
             5.5241E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2773.87821094630        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0674E+00  1.4815E+00  2.0774E+00  8.1938E-01  1.4502E+00  8.0306E-01  1.0412E+00  8.3158E-01  6.7107E-01  1.5048E+00
             2.8111E+00
 PARAMETER:  1.6521E-01  4.9309E-01  8.3111E-01 -9.9207E-02  4.7169E-01 -1.1932E-01  1.4038E-01 -8.4427E-02 -2.9889E-01  5.0867E-01
             1.1336E+00
 GRADIENT:   4.3171E+01  2.4980E+01 -2.9681E+00  3.2325E+01 -1.6281E+01 -3.7806E+01  1.8165E+01 -1.7801E+00 -9.6837E-01 -5.3454E-01
             7.4589E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2792.79766674489        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0378E+00  1.6377E+00  3.8300E+00  7.0139E-01  1.8466E+00  8.7540E-01  9.2233E-01  3.3035E+00  6.5642E-01  1.7699E+00
             2.6526E+00
 PARAMETER:  1.3714E-01  5.9329E-01  1.4429E+00 -2.5468E-01  7.1333E-01 -3.3073E-02  1.9152E-02  1.2950E+00 -3.2095E-01  6.7094E-01
             1.0755E+00
 GRADIENT:  -3.4215E+01 -1.3820E+01 -7.4807E+00  2.8015E+00  3.3835E+01 -5.8244E+00  8.0638E+00 -1.6370E+00  2.6962E+00  5.2474E+00
             7.2806E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2795.53972032571        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0494E+00  1.7595E+00  4.1740E+00  6.3311E-01  1.8059E+00  8.9162E-01  8.7552E-01  3.8851E+00  5.4575E-01  1.7780E+00
             2.6441E+00
 PARAMETER:  1.4826E-01  6.6503E-01  1.5289E+00 -3.5711E-01  6.9103E-01 -1.4715E-02 -3.2938E-02  1.4571E+00 -5.0559E-01  6.7551E-01
             1.0723E+00
 GRADIENT:  -2.0492E+00  2.2500E+01 -6.8013E-01  1.3026E+01  6.4091E-01  1.4295E+00  1.4970E+00  1.5270E+00  3.8979E-01  2.4521E+00
            -2.3637E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2796.54219505149        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      483
 NPARAMETR:  1.0628E+00  1.8884E+00  5.3276E+00  5.4995E-01  1.9210E+00  8.9713E-01  8.3899E-01  4.4661E+00  3.8338E-01  1.8666E+00
             2.6539E+00
 PARAMETER:  1.6093E-01  7.3574E-01  1.7729E+00 -4.9793E-01  7.5287E-01 -8.5504E-03 -7.5557E-02  1.5965E+00 -8.5874E-01  7.2412E-01
             1.0760E+00
 GRADIENT:   2.0759E+01  1.4388E+01 -1.1841E+00  1.1662E+01  2.6412E+00  2.9864E+00 -5.1087E-01  6.9268E-01 -1.7924E-01  3.4393E-01
             2.9402E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2796.59829069098        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      645
 NPARAMETR:  1.0619E+00  1.8892E+00  5.3596E+00  5.4810E-01  1.9214E+00  8.9642E-01  8.3836E-01  4.4797E+00  3.8318E-01  1.8679E+00
             2.6534E+00
 PARAMETER:  1.6003E-01  7.3615E-01  1.7789E+00 -5.0130E-01  7.5305E-01 -9.3496E-03 -7.6304E-02  1.5996E+00 -8.5924E-01  7.2480E-01
             1.0758E+00
 GRADIENT:   2.9128E+01  2.9441E+01 -1.0940E+00  1.2771E+01  5.2375E+00  3.2618E+00 -2.0099E-01  8.8734E-01 -1.0936E-01  1.2000E+00
             4.5745E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2796.60573390296        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      771
 NPARAMETR:  1.0618E+00  1.8882E+00  5.3548E+00  5.4813E-01  1.9204E+00  8.9634E-01  8.3825E-01  4.4833E+00  3.8374E-01  1.8670E+00
             2.6518E+00
 PARAMETER:  1.5998E-01  7.3563E-01  1.7780E+00 -5.0124E-01  7.5252E-01 -9.4381E-03 -7.6436E-02  1.6004E+00 -8.5778E-01  7.2435E-01
             1.0752E+00
 GRADIENT:   2.9102E+01  2.8273E+01 -1.0551E+00  1.2269E+01  5.0849E+00  3.2332E+00 -2.5917E-01  9.2361E-01 -1.0609E-01  1.1083E+00
             3.2849E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2796.61139557057        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  1.0618E+00  1.8840E+00  5.3560E+00  5.4801E-01  1.9195E+00  8.9634E-01  8.3074E-01  4.4826E+00  3.8375E-01  1.8669E+00
             2.6502E+00
 PARAMETER:  1.5993E-01  7.3342E-01  1.7782E+00 -5.0146E-01  7.5209E-01 -9.4390E-03 -8.5443E-02  1.6002E+00 -8.5776E-01  7.2426E-01
             1.0747E+00
 GRADIENT:   2.9135E+01  2.2037E+01 -1.0298E+00  1.0588E+01  5.0625E+00  3.2557E+00 -2.8590E+00  9.1505E-01 -1.9482E-01  1.0874E+00
             1.9478E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2796.61504945226        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      913
 NPARAMETR:  1.0617E+00  1.8781E+00  5.3577E+00  5.4785E-01  1.9184E+00  8.9634E-01  8.3047E-01  4.4816E+00  3.8376E-01  1.8666E+00
             2.6481E+00
 PARAMETER:  1.5985E-01  7.3027E-01  1.7785E+00 -5.0176E-01  7.5147E-01 -9.4404E-03 -8.5763E-02  1.6000E+00 -8.5773E-01  7.2412E-01
             1.0738E+00
 GRADIENT:   2.9300E+01  1.3281E+01 -9.6648E-01  6.8648E+00  5.0506E+00  3.3076E+00 -3.0407E+00  8.9358E-01 -1.6689E-01  1.0414E+00
             5.6512E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2796.61514089655        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1007
 NPARAMETR:  1.0617E+00  1.8767E+00  5.3581E+00  5.4781E-01  1.9181E+00  8.9634E-01  8.3137E-01  4.4814E+00  3.8377E-01  1.8665E+00
             2.6476E+00
 PARAMETER:  1.5983E-01  7.2952E-01  1.7786E+00 -5.0183E-01  7.5132E-01 -9.4408E-03 -8.4682E-02  1.5999E+00 -8.5772E-01  7.2409E-01
             1.0736E+00
 GRADIENT:   1.8681E+01 -6.3808E+00 -1.0341E+00  3.4749E+00  2.1216E+00  2.7985E+00 -2.9271E+00  7.3283E-01 -1.9954E-01  1.1089E-01
            -1.6164E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2799.09666454294        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1146
 NPARAMETR:  1.0459E+00  1.8788E+00  5.1140E+01  5.4434E-01  1.9946E+00  8.8349E-01  8.3458E-01  4.4030E+00  2.8842E-01  1.9367E+00
             2.6559E+00
 PARAMETER:  1.4490E-01  7.3064E-01  4.0346E+00 -5.0819E-01  7.9046E-01 -2.3870E-02 -8.0829E-02  1.5823E+00 -1.1433E+00  7.6098E-01
             1.0768E+00
 GRADIENT:  -1.2489E+01  9.1061E-01 -3.5426E-01 -3.6895E+00 -3.4088E+00 -2.0807E+00 -2.3097E+00 -1.0274E-01  2.6657E-01  4.3974E-01
             4.1552E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2799.80137817429        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1217
 NPARAMETR:  1.0504E+00  1.8465E+00  9.3444E+02  5.6806E-01  2.0106E+00  8.8771E-01  8.6261E-01  4.4031E+00  9.5823E-02  1.9373E+00
             2.6505E+00
 PARAMETER:  1.4915E-01  7.1331E-01  6.9399E+00 -4.6553E-01  7.9845E-01 -1.9115E-02 -4.7797E-02  1.5823E+00 -2.2453E+00  7.6128E-01
             1.0748E+00
 GRADIENT:  -4.4985E-01  4.9999E+00 -2.7861E-02 -3.1370E+00 -1.8767E+00 -3.3491E-01 -8.1066E-01 -4.3828E-04  4.2138E-02 -3.6501E-01
             1.7387E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2799.91225921645        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  1.0519E+00  1.6571E+00  4.6742E+06  7.0355E-01  1.9900E+00  8.8224E-01  9.6784E-01  4.3870E+00  1.0000E-02  1.9156E+00
             2.6513E+00
 PARAMETER:  1.5059E-01  6.0509E-01  1.5458E+01 -2.5162E-01  7.8815E-01 -2.5293E-02  6.7313E-02  1.5787E+00 -5.4861E+00  7.5006E-01
             1.0751E+00
 GRADIENT:  -3.5294E+02  5.6581E+01 -6.1073E-07  3.7410E+01 -2.4243E+01  7.4956E+01 -4.6816E+00  3.3244E-04  0.0000E+00  7.2292E+00
             1.3232E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2800.07144250954        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1491
 NPARAMETR:  1.0556E+00  1.6616E+00  4.8538E+06  6.9967E-01  2.0293E+00  8.8253E-01  9.6779E-01  4.3856E+00  1.0000E-02  1.9425E+00
             2.6528E+00
 PARAMETER:  1.5408E-01  6.0780E-01  1.5495E+01 -2.5715E-01  8.0767E-01 -2.4967E-02  6.7262E-02  1.5783E+00 -5.4861E+00  7.6397E-01
             1.0756E+00
 GRADIENT:  -3.5629E+01  3.6876E+01  4.2198E-04  1.7234E+01  2.4173E-01 -9.6781E+00  3.7985E-01  3.9874E-03  0.0000E+00  6.0179E-01
             1.3607E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2800.07505557344        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1661             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0551E+00  1.6609E+00  5.5432E+06  6.9993E-01  2.0331E+00  8.8280E-01  9.6728E-01  4.3741E+00  1.0000E-02  1.9439E+00
             2.6533E+00
 PARAMETER:  1.5368E-01  6.0734E-01  1.5628E+01 -2.5677E-01  8.0959E-01 -2.4655E-02  6.6731E-02  1.5757E+00 -5.4861E+00  7.6468E-01
             1.0758E+00
 GRADIENT:  -1.7647E+01  1.7726E+01  2.7323E-04  6.1024E+00  2.1371E+00  1.6308E+00 -4.5969E-01  3.0568E-03  0.0000E+00  4.4604E-01
             2.8548E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2800.07848276072        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1832
 NPARAMETR:  1.0549E+00  1.6603E+00  5.2838E+06  7.0109E-01  2.0318E+00  8.8341E-01  9.6763E-01  4.3645E+00  1.0000E-02  1.9480E+00
             2.6531E+00
 PARAMETER:  1.5340E-01  6.0702E-01  1.5580E+01 -2.5512E-01  8.0891E-01 -2.3972E-02  6.7098E-02  1.5735E+00 -5.4861E+00  7.6682E-01
             1.0757E+00
 GRADIENT:  -6.2926E+01  4.8494E+01  3.2481E-04  2.6417E+01  6.1977E-02 -1.0144E+01 -2.5166E+00  3.9951E-03  0.0000E+00  1.2384E+00
             1.7354E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2800.08016791030        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2029
 NPARAMETR:  1.0546E+00  1.6596E+00  5.3390E+06  7.0108E-01  2.0326E+00  8.8366E-01  9.6752E-01  4.2676E+00  1.0000E-02  1.9501E+00
             2.6529E+00
 PARAMETER:  1.5319E-01  6.0655E-01  1.5591E+01 -2.5514E-01  8.0930E-01 -2.3683E-02  6.6986E-02  1.5511E+00 -5.4861E+00  7.6790E-01
             1.0757E+00
 GRADIENT:  -6.2283E+01  3.3822E+01  8.0892E-06  1.8768E+01 -2.9000E+00 -1.7620E+01  2.2743E+00  5.7463E-04  0.0000E+00  1.6949E+00
             2.0811E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2800.08163178934        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2223
 NPARAMETR:  1.0548E+00  1.6605E+00  5.3039E+06  7.0099E-01  2.0315E+00  8.8445E-01  9.6643E-01  4.3342E+00  1.0000E-02  1.9510E+00
             2.6538E+00
 PARAMETER:  1.5336E-01  6.0715E-01  1.5584E+01 -2.5526E-01  8.0878E-01 -2.2791E-02  6.5858E-02  1.5665E+00 -5.4861E+00  7.6832E-01
             1.0760E+00
 GRADIENT:  -2.0849E+01  8.6786E+00  3.3572E-06  5.8341E+00 -1.2813E+00 -5.8809E-01 -1.5252E+00  3.3435E-04  0.0000E+00  5.1789E-01
             8.2036E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2800.08311088495        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2403             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0547E+00  1.6602E+00  5.3231E+06  7.0028E-01  2.0325E+00  8.8502E-01  9.6649E-01  4.3385E+00  1.0000E-02  1.9486E+00
             2.6520E+00
 PARAMETER:  1.5326E-01  6.0691E-01  1.5588E+01 -2.5627E-01  8.0926E-01 -2.2151E-02  6.5913E-02  1.5675E+00 -5.4861E+00  7.6713E-01
             1.0753E+00
 GRADIENT:  -8.9455E+01  3.3774E+01  2.6486E-04  1.5850E+01 -7.5716E-02  1.1101E+01 -5.4602E-01  3.0113E-03  0.0000E+00  1.2345E+00
             4.0753E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2800.08383617849        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2593
 NPARAMETR:  1.0546E+00  1.6612E+00  5.4071E+06  7.0029E-01  2.0346E+00  8.8522E-01  9.6658E-01  4.3478E+00  1.0000E-02  1.9500E+00
             2.6532E+00
 PARAMETER:  1.5316E-01  6.0756E-01  1.5603E+01 -2.5626E-01  8.1028E-01 -2.1914E-02  6.6008E-02  1.5697E+00 -5.4861E+00  7.6783E-01
             1.0758E+00
 GRADIENT:  -1.2616E+02  4.3801E+01 -4.2342E-05  2.8063E+01 -3.2272E+00  5.7098E+00 -3.0103E+00  1.1107E-04  0.0000E+00  2.3107E-01
             2.7245E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2800.08405304065        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2779
 NPARAMETR:  1.0546E+00  1.6604E+00  5.4486E+06  7.0025E-01  2.0331E+00  8.8530E-01  9.6672E-01  4.3432E+00  1.0000E-02  1.9484E+00
             2.6533E+00
 PARAMETER:  1.5314E-01  6.0707E-01  1.5595E+01 -2.5632E-01  8.0956E-01 -2.1842E-02  6.6159E-02  1.5686E+00 -5.4861E+00  7.6676E-01
             1.0758E+00
 GRADIENT:  -1.8348E-01  3.6684E-01 -7.2464E-03  8.6774E-02  3.9322E-02 -9.7251E-01  1.7157E+00  1.7267E-02  0.0000E+00 -9.7514E-02
            -9.3792E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2779
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5153E-03 -5.8734E-03  1.9878E-08 -5.2578E-04 -1.2392E-02
 SE:             2.9459E-02  2.7934E-02  1.4913E-08  2.0440E-04  2.7312E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5898E-01  8.3347E-01  1.8254E-01  1.0102E-02  6.5003E-01

 ETASHRINKSD(%)  1.3099E+00  6.4166E+00  1.0000E+02  9.9315E+01  8.5018E+00
 ETASHRINKVR(%)  2.6027E+00  1.2421E+01  1.0000E+02  9.9995E+01  1.6281E+01
 EBVSHRINKSD(%)  1.9209E+00  6.4782E+00  1.0000E+02  9.9374E+01  6.7955E+00
 EBVSHRINKVR(%)  3.8049E+00  1.2537E+01  1.0000E+02  9.9996E+01  1.3129E+01
 RELATIVEINF(%)  9.6035E+01  6.0582E+00  0.0000E+00  2.7059E-04  5.8321E+01
 EPSSHRINKSD(%)  1.5561E+01
 EPSSHRINKVR(%)  2.8701E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          885
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1626.5212037722706     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2800.0840530406463     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1173.5628492683757     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    81.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2800.084       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.66E+00  5.36E+06  7.00E-01  2.03E+00  8.85E-01  9.67E-01  4.34E+00  1.00E-02  1.95E+00  2.65E+00
 


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
+        1.28E+03
 
 TH 2
+        3.58E+02  3.99E+02
 
 TH 3
+        1.46E-06 -2.01E-07  4.98E-15
 
 TH 4
+        3.99E+02  5.27E+02  2.54E-07  1.10E+03
 
 TH 5
+        1.17E+02  1.33E+01 -2.03E-07  2.16E+01  6.88E+01
 
 TH 6
+        2.74E+02 -3.60E+02  2.53E-06  4.15E+02 -1.52E+02  1.25E+04
 
 TH 7
+       -1.42E+03 -2.45E+02  3.36E-06 -1.39E+02  1.10E+02  4.46E+03 -7.93E+01
 
 TH 8
+       -7.36E+01 -1.43E-01 -1.92E-08  5.45E-02 -3.81E-01  3.80E-01  4.28E+01  5.58E-02
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        1.00E+02  1.42E+01 -3.93E-08  1.08E+02 -7.44E+00  6.81E+01 -6.64E+02  1.67E+00  0.00E+00  4.22E+01
 
 TH11
+        1.12E+01 -3.14E+01 -2.43E-07 -7.91E+01  8.97E+00  2.83E+02  8.67E+01  7.73E-01  0.00E+00  2.08E+01  1.61E+02
 
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
 #CPUT: Total CPU Time in Seconds,       96.484
Stop Time:
Sat Sep 25 02:27:56 CDT 2021
