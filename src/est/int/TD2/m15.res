Wed Sep 29 07:04:10 CDT 2021
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
$DATA ../../../../data/int/TD2/dat15.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
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

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3214.59636057553        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5011E+02 -8.5249E+01  1.1075E+02  9.9118E+01  1.7746E+02  6.0741E+01 -4.0356E+01 -1.2787E+02 -9.1243E+01  7.6633E+00
            -1.2129E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3455.66486574637        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      128
 NPARAMETR:  1.0213E+00  1.0433E+00  8.5964E-01  9.6871E-01  9.6508E-01  9.2284E-01  1.0345E+00  1.1300E+00  1.1784E+00  9.6265E-01
             1.3805E+00
 PARAMETER:  1.2106E-01  1.4239E-01 -5.1239E-02  6.8210E-02  6.4458E-02  1.9703E-02  1.3388E-01  2.2224E-01  2.6412E-01  6.1930E-02
             4.2246E-01
 GRADIENT:  -6.1929E+01 -1.2659E+02 -1.8239E+01 -2.0964E+01  8.6930E+01 -4.3780E+00 -1.9022E+01 -2.5111E+01 -1.9853E+01 -8.0855E-01
            -1.5373E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3467.12955740285        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0044E+00  1.4117E+00  1.4792E+00  7.6705E-01  1.3724E+00  8.9048E-01  8.9537E-01  2.6654E+00  1.4133E+00  1.1896E+00
             1.3286E+00
 PARAMETER:  1.0439E-01  4.4479E-01  4.9149E-01 -1.6521E-01  4.1660E-01 -1.5989E-02 -1.0513E-02  1.0804E+00  4.4595E-01  2.7359E-01
             3.8415E-01
 GRADIENT:  -1.1281E+02 -8.8478E+01  6.9957E+00 -4.7955E+01  7.3175E+01 -2.2561E+01  6.0904E+00 -9.6684E+00  1.1002E+01 -5.5633E+00
            -2.1427E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3473.07677154015        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0058E+00  1.4583E+00  1.4871E+00  7.6135E-01  1.4029E+00  9.0237E-01  8.3425E-01  2.5760E+00  1.4395E+00  1.3101E+00
             1.3643E+00
 PARAMETER:  1.0573E-01  4.7728E-01  4.9682E-01 -1.7267E-01  4.3856E-01 -2.7337E-03 -8.1224E-02  1.0462E+00  4.6432E-01  3.7008E-01
             4.1065E-01
 GRADIENT:   1.4939E+02  2.7076E+02  2.1807E+01  2.2056E+01  1.6443E+02  4.9027E+00  9.0731E+00 -6.5441E+00  3.6042E+01  2.5052E+01
            -1.4387E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3477.61300228581        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      626
 NPARAMETR:  1.0058E+00  1.4583E+00  1.4871E+00  7.6136E-01  1.4028E+00  9.2347E-01  8.1136E-01  2.7262E+00  1.4396E+00  1.3101E+00
             1.4101E+00
 PARAMETER:  1.0574E-01  4.7727E-01  4.9682E-01 -1.7264E-01  4.3844E-01  2.0382E-02 -1.0904E-01  1.1029E+00  4.6436E-01  3.7008E-01
             4.4364E-01
 GRADIENT:  -1.0210E+02 -9.1300E+01 -1.2945E+01 -2.3229E+01  3.7981E+01 -6.5300E+00  5.9897E+00  1.4640E+01  1.0110E+01  6.5737E+00
            -9.1762E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3479.27547341044        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      803
 NPARAMETR:  1.0058E+00  1.4584E+00  1.4871E+00  7.6137E-01  1.4026E+00  9.3926E-01  7.7234E-01  2.7493E+00  1.4395E+00  1.3101E+00
             1.4583E+00
 PARAMETER:  1.0574E-01  4.7735E-01  4.9682E-01 -1.7264E-01  4.3832E-01  3.7338E-02 -1.5833E-01  1.1113E+00  4.6427E-01  3.7007E-01
             4.7730E-01
 GRADIENT:  -1.0060E+02 -8.8006E+01  4.4891E+00 -2.8489E+01  5.8455E+01  1.6604E-01 -8.8231E-01 -5.4080E+00  1.2544E+01  9.7028E+00
             4.9131E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3481.61911822370        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      942
 NPARAMETR:  1.0508E+00  1.4623E+00  1.4830E+00  8.1341E-01  1.2855E+00  9.3123E-01  8.3025E-01  2.8456E+00  1.2945E+00  1.2063E+00
             1.4558E+00
 PARAMETER:  1.4954E-01  4.8000E-01  4.9405E-01 -1.0652E-01  3.5115E-01  2.8748E-02 -8.6028E-02  1.1458E+00  3.5814E-01  2.8760E-01
             4.7556E-01
 GRADIENT:   8.2732E+00  1.9099E+01  2.0220E+01  5.6419E+00 -1.0714E+01  9.4778E-01 -3.6446E-01 -4.5766E-01 -8.4933E-01  8.0748E-02
             2.4004E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3485.01483126969        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1127
 NPARAMETR:  1.0503E+00  1.4604E+00  1.4853E+00  7.8654E-01  1.3146E+00  9.2974E-01  8.1097E-01  2.8595E+00  1.3257E+00  1.1949E+00
             1.4552E+00
 PARAMETER:  1.4903E-01  4.7873E-01  4.9561E-01 -1.4011E-01  3.7352E-01  2.7145E-02 -1.0952E-01  1.1506E+00  3.8194E-01  2.7805E-01
             4.7516E-01
 GRADIENT:   7.8383E+00 -2.1186E+01  7.5326E+00 -2.0443E+01  7.3488E-01  5.0539E-01 -4.6023E-01  3.4818E+00  3.4623E-01 -9.7049E-02
            -7.4079E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3485.21935083042        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1308             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0476E+00  1.4740E+00  1.4819E+00  7.8690E-01  1.3194E+00  9.2852E-01  8.1219E-01  2.8130E+00  1.3292E+00  1.2020E+00
             1.4641E+00
 PARAMETER:  1.4646E-01  4.8795E-01  4.9330E-01 -1.3966E-01  3.7720E-01  2.5834E-02 -1.0802E-01  1.1343E+00  3.8455E-01  2.8396E-01
             4.8121E-01
 GRADIENT:   3.1014E+02  2.8615E+02  1.6601E+01  2.4255E+01  7.5311E+01  1.6254E+01  4.3537E+00  1.2499E+01  1.6589E+01  8.4385E+00
             1.2795E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3485.27068149937        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1465
 NPARAMETR:  1.0475E+00  1.4758E+00  1.4829E+00  7.9069E-01  1.3176E+00  9.2837E-01  8.0881E-01  2.8060E+00  1.3295E+00  1.2006E+00
             1.4630E+00
 PARAMETER:  1.4641E-01  4.8917E-01  4.9402E-01 -1.3485E-01  3.7584E-01  2.5671E-02 -1.1220E-01  1.1317E+00  3.8479E-01  2.8283E-01
             4.8046E-01
 GRADIENT:   3.1026E+02  2.9342E+02  1.6933E+01  2.7976E+01  7.3995E+01  1.6197E+01  3.9525E+00  1.2072E+01  1.6825E+01  8.1109E+00
             1.0953E+01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -3485.27068149937        NO. OF FUNC. EVALS.:  62
 CUMULATIVE NO. OF FUNC. EVALS.:     1527
 NPARAMETR:  1.0475E+00  1.4758E+00  1.4829E+00  7.9069E-01  1.3176E+00  9.2837E-01  8.0881E-01  2.8060E+00  1.3295E+00  1.2006E+00
             1.4630E+00
 PARAMETER:  1.4641E-01  4.8917E-01  4.9402E-01 -1.3485E-01  3.7584E-01  2.5671E-02 -1.1220E-01  1.1317E+00  3.8479E-01  2.8283E-01
             4.8046E-01
 GRADIENT:  -1.0066E+00  3.0165E+02 -2.8534E+02  1.0883E+03  2.0697E-01 -1.1604E-01 -3.0680E-02  1.2655E+02  3.0346E-01  6.1985E-02
            -2.9946E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1527
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3670E-03 -4.6981E-02 -4.5807E-02  3.3538E-02 -4.5631E-02
 SE:             2.9804E-02  2.0411E-02  2.0656E-02  2.5514E-02  2.4365E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6342E-01  2.1347E-02  2.6581E-02  1.8868E-01  6.1097E-02

 ETASHRINKSD(%)  1.5429E-01  3.1622E+01  3.0800E+01  1.4525E+01  1.8374E+01
 ETASHRINKVR(%)  3.0834E-01  5.3244E+01  5.2114E+01  2.6940E+01  3.3371E+01
 EBVSHRINKSD(%)  6.3010E-01  3.2640E+01  3.3078E+01  1.5576E+01  1.6529E+01
 EBVSHRINKVR(%)  1.2562E+00  5.4626E+01  5.5214E+01  2.8726E+01  3.0325E+01
 RELATIVEINF(%)  9.8734E+01  1.4872E+01  3.5735E+01  2.7350E+01  3.5817E+01
 EPSSHRINKSD(%)  2.0857E+01
 EPSSHRINKVR(%)  3.7363E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3485.2706814993703     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1831.1813217309596     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    50.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3485.271       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.48E+00  1.48E+00  7.91E-01  1.32E+00  9.28E-01  8.09E-01  2.81E+00  1.33E+00  1.20E+00  1.46E+00
 


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
+        1.58E+05
 
 TH 2
+       -2.14E+01  7.53E+03
 
 TH 3
+        1.69E+01 -4.37E+01  6.69E+03
 
 TH 4
+       -1.19E+02  4.85E+04 -3.43E+00  3.25E+05
 
 TH 5
+        4.82E+04  5.32E+01 -2.27E+02 -6.92E+04  3.13E+02
 
 TH 6
+        3.75E+00 -5.83E+00  4.41E+00 -3.17E+01 -6.77E-01  2.25E+02
 
 TH 7
+        2.21E-01 -3.34E+01  4.03E+01 -3.02E+02  1.29E+00 -5.19E-01  8.02E+01
 
 TH 8
+       -3.65E+00  9.41E+00  9.67E+01 -6.32E+00  3.32E+01 -1.03E+00 -9.76E+00  3.51E+02
 
 TH 9
+        1.43E+00  1.48E+01 -3.04E+01  2.63E+02 -2.38E-01 -3.84E-01  1.60E+01  1.01E+01  6.23E+01
 
 TH10
+        3.18E-01 -1.41E+02  1.31E+02 -8.24E+02 -1.18E+01  3.73E-01  1.66E+01 -2.66E+01 -6.48E-01  6.76E+01
 
 TH11
+        6.33E+00 -7.21E+03 -2.72E+02 -4.83E+04 -2.06E+02  7.59E+00  5.47E+01  6.53E+01 -2.89E+01  1.38E+02  7.69E+03
 
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
 #CPUT: Total CPU Time in Seconds,       66.377
Stop Time:
Wed Sep 29 07:05:18 CDT 2021
