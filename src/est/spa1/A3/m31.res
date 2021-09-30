Thu Sep 30 00:09:09 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat31.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -213.560843235597        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2322E+02  1.3286E+02  1.2378E+02  1.4156E+02  4.0551E+02  2.3856E+01 -1.2615E+02 -1.6416E+02 -1.0819E+02 -2.7161E+02
            -3.0019E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1493.33863337810        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.5127E-01  7.7261E-01  8.4501E-01  1.1050E+00  7.3139E-01  9.4019E-01  1.0101E+00  9.4551E-01  1.2709E+00  9.9194E-01
             2.8549E+00
 PARAMETER:  5.0040E-02 -1.5799E-01 -6.8403E-02  1.9985E-01 -2.1281E-01  3.8332E-02  1.1001E-01  4.3973E-02  3.3973E-01  9.1906E-02
             1.1490E+00
 GRADIENT:  -1.1884E+01  4.0061E+00 -9.1824E+00  2.0975E+01  8.7840E+01 -1.6301E+01  1.4699E+00  1.0912E+01  1.5388E+01 -1.0822E+01
            -5.5932E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1516.67274021328        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      199
 NPARAMETR:  9.5616E-01  3.9982E-01  4.3992E-01  1.4383E+00  3.5981E-01  9.3878E-01  7.8066E-01  5.1865E-01  1.6681E+00  5.7194E-01
             2.6446E+00
 PARAMETER:  5.5170E-02 -8.1674E-01 -7.2116E-01  4.6344E-01 -9.2219E-01  3.6828E-02 -1.4762E-01 -5.5652E-01  6.1167E-01 -4.5872E-01
             1.0725E+00
 GRADIENT:  -4.9801E+01  5.9384E+01  1.4389E+02  1.4570E+02 -1.5232E+02 -2.6086E+01 -5.2653E-01 -2.4808E+00  5.2697E+01 -1.5502E+01
            -7.4930E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1557.53173792445        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0166E+00  3.8200E-01  3.0762E-01  1.1898E+00  2.9470E-01  1.0853E+00  4.9604E-01  2.6914E-02  1.4450E+00  7.6877E-01
             2.7684E+00
 PARAMETER:  1.1647E-01 -8.6234E-01 -1.0789E+00  2.7377E-01 -1.1218E+00  1.8181E-01 -6.0111E-01 -3.5151E+00  4.6809E-01 -1.6297E-01
             1.1183E+00
 GRADIENT:   6.5308E+01  2.8845E+01  6.7996E+01  5.9404E+01 -8.2160E+01  2.4514E+01  1.4821E-01  9.2801E-03  2.6425E+01  2.9831E+01
             1.0878E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1576.75833130769        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  9.7907E-01  3.7441E-01  1.9379E-01  9.9939E-01  2.3524E-01  1.0098E+00  4.5197E-01  1.0000E-02  1.4830E+00  6.1770E-01
             2.6634E+00
 PARAMETER:  7.8845E-02 -8.8241E-01 -1.5410E+00  9.9392E-02 -1.3471E+00  1.0977E-01 -6.9413E-01 -6.4665E+00  4.9405E-01 -3.8175E-01
             1.0796E+00
 GRADIENT:  -3.1472E+00  3.4147E+00  4.6676E+00 -1.7272E-01 -7.8311E+00 -1.6750E+00  7.0292E-01  0.0000E+00  1.8226E+00  4.6886E-01
             3.6920E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1577.06862999148        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  9.8037E-01  3.7018E-01  1.8778E-01  9.9068E-01  2.3146E-01  1.0147E+00  1.0411E-01  1.0000E-02  1.4888E+00  6.2890E-01
             2.6484E+00
 PARAMETER:  8.0175E-02 -8.9377E-01 -1.5725E+00  9.0640E-02 -1.3633E+00  1.1460E-01 -2.1623E+00 -6.4757E+00  4.9794E-01 -3.6378E-01
             1.0740E+00
 GRADIENT:   1.7363E-02  7.4876E-02 -1.1964E+00 -8.2520E-02  1.1065E+00 -1.0708E-01  2.5730E-02  0.0000E+00 -6.5181E-01 -7.3252E-01
            -1.3987E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1577.08268960210        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  9.8044E-01  3.7045E-01  1.8897E-01  9.9239E-01  2.3200E-01  1.0149E+00  3.7247E-02  1.0000E-02  1.4897E+00  6.3110E-01
             2.6527E+00
 PARAMETER:  8.0249E-02 -8.9303E-01 -1.5661E+00  9.2362E-02 -1.3610E+00  1.1476E-01 -3.1902E+00 -6.4100E+00  4.9860E-01 -3.6029E-01
             1.0756E+00
 GRADIENT:  -1.2830E+01  2.6196E+00  1.8248E+00  2.1261E+00 -6.9099E+00 -1.7246E+00 -2.6658E-02  0.0000E+00 -2.3348E+00  2.1242E-01
             6.1342E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1577.08482451390        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1091            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8046E-01  3.7012E-01  1.8880E-01  9.9231E-01  2.3198E-01  1.0148E+00  2.4069E-02  1.0000E-02  1.4898E+00  6.3140E-01
             2.6527E+00
 PARAMETER:  8.0271E-02 -8.9394E-01 -1.5671E+00  9.2285E-02 -1.3611E+00  1.1464E-01 -3.6268E+00 -6.4100E+00  4.9866E-01 -3.5981E-01
             1.0756E+00
 GRADIENT:   4.5572E+01  1.1673E+01  2.7260E+01  1.0791E+01  1.4987E+02  4.5348E+00  4.3494E-03  0.0000E+00  1.2063E+01  1.2218E+00
             1.0738E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1577.08552601125        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1267
 NPARAMETR:  9.8052E-01  3.6968E-01  1.8859E-01  9.9185E-01  2.3177E-01  1.0148E+00  1.0000E-02  1.0000E-02  1.4901E+00  6.3158E-01
             2.6526E+00
 PARAMETER:  8.0333E-02 -8.9512E-01 -1.5682E+00  9.1817E-02 -1.3620E+00  1.1469E-01 -4.8967E+00 -6.4100E+00  4.9884E-01 -3.5953E-01
             1.0756E+00
 GRADIENT:   2.0614E-01  6.3159E-02  3.4535E-01 -1.1154E-01 -1.7262E-01  2.4332E-02  0.0000E+00  0.0000E+00  7.7359E-02  3.9816E-03
             2.8405E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1577.08557895357        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1442
 NPARAMETR:  9.8042E-01  3.6951E-01  1.8837E-01  9.9171E-01  2.3162E-01  1.0147E+00  1.0000E-02  1.0000E-02  1.4902E+00  6.3168E-01
             2.6523E+00
 PARAMETER:  8.0228E-02 -8.9557E-01 -1.5694E+00  9.1671E-02 -1.3627E+00  1.1464E-01 -5.5264E+00 -6.4100E+00  4.9890E-01 -3.5938E-01
             1.0754E+00
 GRADIENT:  -8.9491E-03 -2.6770E-02  3.8454E-02  4.5764E-02  6.8083E-02 -3.0803E-03  0.0000E+00  0.0000E+00 -1.3463E-02 -3.6456E-03
             1.0677E-03

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1577.08558073382        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1542
 NPARAMETR:  9.8100E-01  3.7012E-01  1.8828E-01  9.9102E-01  2.3168E-01  1.0151E+00  1.0000E-02  1.0000E-02  1.4978E+00  6.3197E-01
             2.6509E+00
 PARAMETER:  8.0224E-02 -8.9563E-01 -1.5697E+00  9.1571E-02 -1.3628E+00  1.1465E-01 -5.6933E+00 -6.4100E+00  4.9900E-01 -3.5932E-01
             1.0754E+00
 GRADIENT:  -1.2203E-01 -8.1322E-02  1.6766E-02  5.2653E-02 -1.5053E-01 -1.0695E-02  0.0000E+00  0.0000E+00 -1.4667E-01 -6.6287E-03
             4.7306E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1542
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4610E-03 -9.1612E-05  1.8331E-04 -6.9299E-03  3.8661E-03
 SE:             2.9110E-02  1.3534E-04  2.2809E-04  2.7345E-02  2.3978E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5997E-01  4.9847E-01  4.2159E-01  7.9994E-01  8.7190E-01

 ETASHRINKSD(%)  2.4785E+00  9.9547E+01  9.9236E+01  8.3898E+00  1.9672E+01
 ETASHRINKVR(%)  4.8956E+00  9.9998E+01  9.9994E+01  1.6076E+01  3.5475E+01
 EBVSHRINKSD(%)  2.2445E+00  9.9540E+01  9.9300E+01  5.9859E+00  2.0076E+01
 EBVSHRINKVR(%)  4.4386E+00  9.9998E+01  9.9995E+01  1.1613E+01  3.6121E+01
 RELATIVEINF(%)  9.5474E+01  3.6826E-04  3.9168E-04  4.8549E+01  2.6904E+00
 EPSSHRINKSD(%)  2.6965E+01
 EPSSHRINKVR(%)  4.6659E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1577.0855807338212     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -658.14704752914849     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.03
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
 





 #OBJV:********************************************    -1577.086       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  3.69E-01  1.88E-01  9.92E-01  2.32E-01  1.01E+00  1.00E-02  1.00E-02  1.49E+00  6.32E-01  2.65E+00
 


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
+        1.07E+03
 
 TH 2
+       -3.61E+01  1.77E+03
 
 TH 3
+        1.45E+02  3.40E+03  1.50E+04
 
 TH 4
+       -9.21E+00  1.09E+02 -4.90E+02  4.53E+02
 
 TH 5
+        1.63E+02 -6.31E+03 -2.00E+04 -3.19E+02  3.36E+04
 
 TH 6
+        3.89E+00 -1.43E+01  8.64E+01 -7.71E+00 -9.09E+00  1.75E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.20E+00 -3.57E+01  1.45E+02 -4.12E+00  5.60E+01 -1.00E+00  0.00E+00  0.00E+00  6.67E+01
 
 TH10
+       -4.80E+00 -2.12E+01  2.68E+01  8.55E+00  1.07E+02  4.17E+00  0.00E+00  0.00E+00  3.05E+00  2.05E+02
 
 TH11
+       -1.65E+01 -8.14E+00 -6.48E+01 -5.60E+00  4.90E+01  2.37E+00  0.00E+00  0.00E+00  3.72E+00  2.24E+01  6.60E+01
 
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
 #CPUT: Total CPU Time in Seconds,       31.560
Stop Time:
Thu Sep 30 00:09:52 CDT 2021
