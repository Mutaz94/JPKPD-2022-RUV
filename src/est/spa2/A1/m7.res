Thu Sep 30 04:48:59 CDT 2021
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
$DATA ../../../../data/spa2/A1/dat7.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1878.70560888261        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8128E+02  4.2742E+01  6.9118E+01  1.5206E+02  1.2499E+02  4.6786E+00 -5.9711E-01 -1.1480E+02  2.8442E+01 -1.6484E+01
            -9.9781E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2076.49466320796        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.6350E-01  1.0600E+00  7.7188E-01  9.4038E-01  8.4535E-01  1.1122E+00  9.1559E-01  1.2600E+00  7.9632E-01  9.2779E-01
             2.4048E+00
 PARAMETER:  6.2814E-02  1.5826E-01 -1.5892E-01  3.8527E-02 -6.8009E-02  2.0631E-01  1.1809E-02  3.3109E-01 -1.2776E-01  2.5049E-02
             9.7745E-01
 GRADIENT:   2.2571E+01  3.6004E+01  1.1027E+01 -1.8392E+01 -6.5397E+01  2.8099E+01  3.7989E+00  1.2068E+01 -2.0879E+00  1.2379E+01
             2.9012E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2087.80714876754        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.7427E-01  7.5203E-01  5.0851E-01  1.1170E+00  5.8181E-01  1.0753E+00  1.4932E+00  8.4830E-01  7.6516E-01  3.6053E-01
             2.2207E+00
 PARAMETER:  7.3934E-02 -1.8498E-01 -5.7626E-01  2.1067E-01 -4.4161E-01  1.7261E-01  5.0095E-01 -6.4523E-02 -1.6768E-01 -9.2019E-01
             8.9782E-01
 GRADIENT:   4.5342E+01  5.6835E+01 -2.8871E+01  1.5017E+02  4.3107E+01  9.4882E+00  2.6578E+01  9.2987E+00 -5.9949E+00 -8.8228E-01
             2.1967E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2110.56415215932        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      250
 NPARAMETR:  9.4973E-01  7.2361E-01  4.5326E-01  1.1132E+00  5.4147E-01  1.0991E+00  1.2994E+00  7.2178E-01  8.5211E-01  4.8866E-01
             1.8980E+00
 PARAMETER:  4.8421E-02 -2.2350E-01 -6.9129E-01  2.0725E-01 -5.1346E-01  1.9452E-01  3.6187E-01 -2.2603E-01 -6.0037E-02 -6.1608E-01
             7.4082E-01
 GRADIENT:  -7.4080E+01  4.6583E+01 -5.4011E+01  9.6119E+01 -7.7848E+00 -3.1130E+00 -3.1188E-02  5.7479E+00  5.0995E+00 -1.5110E+00
             8.6578E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2130.60218824095        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      427
 NPARAMETR:  9.7681E-01  1.1392E+00  9.3483E-01  9.2021E-01  1.0013E+00  1.0836E+00  8.7899E-01  1.2842E+00  8.7464E-01  9.7851E-01
             1.7489E+00
 PARAMETER:  7.6541E-02  2.3035E-01  3.2609E-02  1.6846E-02  1.0129E-01  1.8032E-01 -2.8977E-02  3.5011E-01 -3.3941E-02  7.8271E-02
             6.5899E-01
 GRADIENT:  -1.9051E+01  3.1003E+01  9.1933E+00  2.6145E+01  1.6238E+00 -3.9944E+00  5.4763E-02 -1.2539E+01 -5.4928E+00 -1.8431E+00
             1.3976E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2134.11174794884        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  9.9129E-01  1.1319E+00  1.0050E+00  9.1075E-01  1.0373E+00  1.1050E+00  8.3039E-01  1.8938E+00  9.0145E-01  9.4959E-01
             1.6849E+00
 PARAMETER:  9.1248E-02  2.2387E-01  1.0496E-01  6.5126E-03  1.3667E-01  1.9986E-01 -8.5862E-02  7.3859E-01 -3.7515E-03  4.8275E-02
             6.2172E-01
 GRADIENT:   1.0316E+01 -6.5398E+00 -6.8313E+00  7.4636E+00  1.4002E+01  3.2404E+00 -7.6793E-01 -1.9250E+00 -5.4667E+00 -2.4531E+00
            -1.9747E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2134.57027574970        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      764
 NPARAMETR:  9.8570E-01  1.1454E+00  1.0444E+00  9.0128E-01  1.0487E+00  1.0982E+00  8.2160E-01  2.0003E+00  9.2465E-01  9.6177E-01
             1.6809E+00
 PARAMETER:  8.5593E-02  2.3578E-01  1.4346E-01 -3.9441E-03  1.4753E-01  1.9366E-01 -9.6502E-02  7.9330E-01  2.1660E-02  6.1018E-02
             6.1930E-01
 GRADIENT:   1.5687E+02  5.0487E+01  7.9979E-01  2.7196E+01  1.2536E+01  3.6241E+01  2.0851E+00  6.3187E+00 -1.2641E+00 -1.3018E-01
             7.7360E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2134.59411625065        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      947
 NPARAMETR:  9.8566E-01  1.1457E+00  1.0461E+00  9.0118E-01  1.0484E+00  1.0965E+00  8.0679E-01  1.9965E+00  9.2872E-01  9.6878E-01
             1.6833E+00
 PARAMETER:  8.5559E-02  2.3600E-01  1.4509E-01 -4.0555E-03  1.4723E-01  1.9215E-01 -1.1470E-01  7.9142E-01  2.6047E-02  6.8286E-02
             6.2077E-01
 GRADIENT:   4.5369E-01 -3.9141E+00 -5.8318E-01  1.9282E+00 -1.1770E+00  3.7392E-01  1.1330E-01 -1.9799E+00 -3.2847E+00 -1.8617E-01
             1.8910E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2134.68244862246        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1123
 NPARAMETR:  9.8541E-01  1.1552E+00  1.0446E+00  8.9542E-01  1.0539E+00  1.0952E+00  7.7509E-01  1.9902E+00  9.5185E-01  9.7365E-01
             1.6873E+00
 PARAMETER:  8.5301E-02  2.4426E-01  1.4367E-01 -1.0465E-02  1.5252E-01  1.9097E-01 -1.5477E-01  7.8821E-01  5.0656E-02  7.3300E-02
             6.2311E-01
 GRADIENT:  -1.4711E-01 -3.9770E+00 -9.2493E-02  2.7508E+00 -2.2472E-01 -8.3721E-02 -5.1578E-01 -3.4754E+00 -1.6266E+00 -7.1864E-01
             4.6234E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2134.70660568275        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1307
 NPARAMETR:  9.8522E-01  1.1577E+00  1.0433E+00  8.9305E-01  1.0553E+00  1.0950E+00  7.7330E-01  1.9898E+00  9.5830E-01  9.7796E-01
             1.6874E+00
 PARAMETER:  8.5105E-02  2.4642E-01  1.4235E-01 -1.3112E-02  1.5378E-01  1.9079E-01 -1.5708E-01  7.8805E-01  5.7406E-02  7.7718E-02
             6.2320E-01
 GRADIENT:  -4.8164E-01 -5.2329E+00  2.9904E-02  1.6713E+00 -2.7580E-01 -1.4765E-01 -1.7984E-01 -3.6975E+00 -7.0204E-01 -3.4010E-01
             5.1850E+00

0ITERATION NO.:   46    OBJECTIVE VALUE:  -2134.70660568275        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     1334
 NPARAMETR:  9.8535E-01  1.1574E+00  1.0434E+00  8.9314E-01  1.0551E+00  1.0952E+00  7.7437E-01  1.9883E+00  9.5820E-01  9.7894E-01
             1.6884E+00
 PARAMETER:  8.5105E-02  2.4642E-01  1.4235E-01 -1.3112E-02  1.5378E-01  1.9079E-01 -1.5708E-01  7.8805E-01  5.7406E-02  7.7718E-02
             6.2320E-01
 GRADIENT:  -6.1186E-01  1.8925E+02 -3.4005E+02 -4.8182E+02  3.1482E+02 -3.1408E-01 -1.6367E-01  4.5610E+01  4.8387E+02 -3.0887E-01
            -5.8840E+01
 NUMSIGDIG:         2.2         2.3         2.3         2.3         2.3         2.4         1.4         2.3         2.3         1.3
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1334
 NO. OF SIG. DIGITS IN FINAL EST.:  1.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6340E-03 -2.4282E-02 -2.7177E-02  1.1816E-02 -3.3553E-02
 SE:             2.9779E-02  1.8208E-02  1.9197E-02  2.4820E-02  2.1527E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5624E-01  1.8234E-01  1.5687E-01  6.3403E-01  1.1907E-01

 ETASHRINKSD(%)  2.3566E-01  3.9001E+01  3.5686E+01  1.6851E+01  2.7882E+01
 ETASHRINKVR(%)  4.7077E-01  6.2791E+01  5.8638E+01  3.0862E+01  4.7990E+01
 EBVSHRINKSD(%)  7.2964E-01  3.9128E+01  3.9158E+01  1.8444E+01  2.7032E+01
 EBVSHRINKVR(%)  1.4540E+00  6.2946E+01  6.2983E+01  3.3486E+01  4.6756E+01
 RELATIVEINF(%)  9.8498E+01  4.7846E+00  1.6061E+01  1.1103E+01  1.2653E+01
 EPSSHRINKSD(%)  2.8630E+01
 EPSSHRINKVR(%)  4.9064E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2134.7066056827530     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1031.9803658371459     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.71
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2134.707       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.16E+00  1.04E+00  8.93E-01  1.06E+00  1.10E+00  7.73E-01  1.99E+00  9.58E-01  9.78E-01  1.69E+00
 


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
+        9.41E+02
 
 TH 2
+       -1.64E+01  1.50E+04
 
 TH 3
+        1.01E+01  4.35E+02  5.46E+04
 
 TH 4
+       -6.33E+00  1.23E+03  9.05E+04  1.51E+05
 
 TH 5
+        1.36E+02 -4.80E+02 -5.02E+04 -8.30E+04  4.64E+04
 
 TH 6
+        2.32E+00 -9.54E+00  9.69E+00  8.65E+00  3.57E+04  1.61E+02
 
 TH 7
+        1.08E+00  2.18E+01 -7.36E+01 -1.39E+02  6.15E+04  1.09E-04  4.10E+01
 
 TH 8
+       -9.38E-01 -1.11E+02  6.44E+01  1.54E+02 -5.79E+01 -1.77E+00  3.89E+00  2.96E+02
 
 TH 9
+       -8.88E+00 -4.80E+02 -8.46E+04 -1.41E+05  7.76E+04 -1.29E+01  1.64E+02 -1.07E+02  1.31E+05
 
 TH10
+       -4.44E-02  4.03E+01 -8.50E+01 -1.40E+02 -7.63E+04  5.05E-01  1.18E+01  1.42E+01  1.39E+02  6.67E+01
 
 TH11
+       -6.77E+00  1.38E+02 -1.24E+02 -2.37E+02  6.63E+01  4.52E+00  4.51E+00  2.92E+02  1.66E+02 -2.21E+00  8.61E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.463
Stop Time:
Thu Sep 30 04:49:37 CDT 2021
