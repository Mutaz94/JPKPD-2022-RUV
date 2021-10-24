Sun Oct 24 01:10:23 CDT 2021
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
$DATA ../../../../data/SD3/D2/dat10.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m10.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1991.55991206867        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6677E+02 -5.1408E+01 -5.6493E+01 -3.9906E+01  8.1753E+01 -7.2083E+01 -5.4554E+01 -3.2728E+00 -1.1953E+02  9.5038E+00
             5.0086E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2042.91315391570        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      131
 NPARAMETR:  1.0571E+00  9.9579E-01  1.1213E+00  1.0660E+00  9.5762E-01  1.1499E+00  1.1865E+00  1.0284E+00  1.4383E+00  9.3258E-01
             9.1781E-01
 PARAMETER:  1.5548E-01  9.5785E-02  2.1445E-01  1.6392E-01  5.6691E-02  2.3967E-01  2.7099E-01  1.2799E-01  4.6349E-01  3.0199E-02
             1.4233E-02
 GRADIENT:   5.3779E+01  2.9852E-02 -1.6845E+00 -1.4533E+01 -4.4615E+01 -7.0483E+01 -1.2200E+01 -1.1101E+01  1.0389E+01  4.7525E+00
            -1.1251E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2050.48432526066        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.0557E+00  9.1216E-01  1.5680E+00  1.1572E+00  1.2130E+00  1.3002E+00  1.9586E+00  1.6398E+00  1.2003E+00  9.6074E-01
             9.6515E-01
 PARAMETER:  1.5423E-01  8.0587E-03  5.4982E-01  2.4598E-01  2.9310E-01  3.6248E-01  7.7222E-01  5.9456E-01  2.8255E-01  5.9951E-02
             6.4524E-02
 GRADIENT:   4.2689E+01 -1.6686E+01 -3.8596E+01  1.6791E+01  1.0711E+02 -9.8798E+00  3.9472E+00  6.4471E+00  1.0821E+01 -1.9744E+01
             2.6293E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2056.63478123531        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  1.0351E+00  6.0023E-01  1.8716E+00  1.4226E+00  1.0962E+00  1.3125E+00  2.3880E+00  1.5722E+00  1.0914E+00  9.7377E-01
             9.4398E-01
 PARAMETER:  1.3449E-01 -4.1045E-01  7.2679E-01  4.5252E-01  1.9188E-01  3.7193E-01  9.7048E-01  5.5250E-01  1.8749E-01  7.3424E-02
             4.2344E-02
 GRADIENT:   2.0165E+01  1.4052E+01 -4.9369E+00  5.1741E+01  2.3659E+01 -3.3472E+00  1.4805E-01 -3.6018E+00 -2.1919E+00 -3.4552E+00
             1.2092E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2058.76261771885        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  1.0225E+00  6.0895E-01  1.6743E+00  1.3560E+00  1.0306E+00  1.3277E+00  2.4922E+00  1.4726E+00  1.0598E+00  9.2559E-01
             9.2488E-01
 PARAMETER:  1.2220E-01 -3.9602E-01  6.1542E-01  4.0452E-01  1.3010E-01  3.8344E-01  1.0132E+00  4.8704E-01  1.5808E-01  2.2681E-02
             2.1909E-02
 GRADIENT:   3.8086E+00  4.9413E+00  4.3160E+00  3.2771E+00 -4.9030E+00  1.7906E+00  1.0650E+00 -3.1330E-01  2.8720E+00 -1.0747E+00
            -2.0873E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2058.80942469218        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      806
 NPARAMETR:  1.0202E+00  6.0771E-01  1.6796E+00  1.3532E+00  1.0379E+00  1.3350E+00  2.4793E+00  1.4690E+00  1.0607E+00  9.4119E-01
             9.2832E-01
 PARAMETER:  1.2001E-01 -3.9805E-01  6.1857E-01  4.0244E-01  1.3719E-01  3.8893E-01  1.0080E+00  4.8455E-01  1.5889E-01  3.9393E-02
             2.5621E-02
 GRADIENT:   9.5412E-01  2.8201E+00  2.1952E+00  5.5260E-01 -3.6626E-02  4.1706E+00  4.1461E-01 -5.9338E-01  3.1307E+00  1.6420E-02
             1.1028E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2058.99003762796        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      989
 NPARAMETR:  1.0197E+00  5.8985E-01  1.6601E+00  1.3548E+00  1.0345E+00  1.3244E+00  2.5801E+00  1.4674E+00  1.0095E+00  9.4451E-01
             9.2694E-01
 PARAMETER:  1.1952E-01 -4.2789E-01  6.0690E-01  4.0364E-01  1.3394E-01  3.8094E-01  1.0478E+00  4.8349E-01  1.0948E-01  4.2906E-02
             2.4138E-02
 GRADIENT:   4.4909E-01  4.5137E-01  3.4807E-01 -3.1328E+00  3.2402E+00  6.1583E-01 -7.4223E-01  2.5149E-01  4.3376E-01 -1.6665E-01
            -2.3191E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2059.05620464620        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1166
 NPARAMETR:  1.0168E+00  5.2719E-01  1.6377E+00  1.3946E+00  1.0095E+00  1.3113E+00  2.8163E+00  1.4305E+00  9.7480E-01  9.4014E-01
             9.2714E-01
 PARAMETER:  1.1671E-01 -5.4019E-01  5.9329E-01  4.3259E-01  1.0949E-01  3.7099E-01  1.1354E+00  4.5800E-01  7.4478E-02  3.8274E-02
             2.4350E-02
 GRADIENT:  -2.3491E+00  1.9213E+00  5.7247E-03  2.2532E+00  6.2573E-01 -3.6681E+00  9.0337E-01 -5.3328E-01 -6.4172E-01  9.2227E-01
             4.3901E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2059.10270043247        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1327
 NPARAMETR:  1.0221E+00  4.8767E-01  1.6364E+00  1.4135E+00  1.0005E+00  1.3600E+00  2.8854E+00  1.4392E+00  9.9228E-01  9.1779E-01
             9.2660E-01
 PARAMETER:  1.2189E-01 -6.1812E-01  5.9251E-01  4.4607E-01  1.0049E-01  4.0752E-01  1.1596E+00  4.6410E-01  9.2250E-02  1.4209E-02
             2.3771E-02
 GRADIENT:   6.0828E+02  8.7633E+01  8.4446E+00  7.7584E+02  1.4820E+01  3.9445E+02  1.4473E+02  3.1004E+00  1.4468E+01  1.4683E-01
             1.0785E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2059.13598148892        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1511
 NPARAMETR:  1.0221E+00  4.8767E-01  1.6566E+00  1.4104E+00  1.0010E+00  1.3203E+00  2.9021E+00  1.4359E+00  9.7962E-01  9.2509E-01
             9.2656E-01
 PARAMETER:  1.2189E-01 -6.1811E-01  6.0476E-01  4.4385E-01  1.0099E-01  3.7783E-01  1.1654E+00  4.6178E-01  7.9412E-02  2.2140E-02
             2.3725E-02
 GRADIENT:   5.1540E+00 -4.7385E-01  1.6922E+00 -9.5212E+00 -4.5662E-01 -5.6722E-01  1.2495E+00 -9.9406E-01 -4.3783E-01 -5.8651E-02
            -1.4127E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2059.18563475360        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1694
 NPARAMETR:  1.0203E+00  4.9027E-01  1.6497E+00  1.4099E+00  9.9965E-01  1.3379E+00  2.9103E+00  1.4399E+00  9.8006E-01  9.2237E-01
             9.2658E-01
 PARAMETER:  1.2014E-01 -6.1280E-01  6.0059E-01  4.4351E-01  9.9651E-02  3.9109E-01  1.1683E+00  4.6460E-01  7.9861E-02  1.9189E-02
             2.3746E-02
 GRADIENT:   2.8178E+00  4.9190E-02  9.7056E-01 -7.8606E+00 -4.8052E-01  5.1206E+00  1.8143E+00 -3.2743E-01  5.8002E-02 -1.3226E-01
            -1.5145E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2059.18756961660        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1880             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0208E+00  4.9099E-01  1.6473E+00  1.4094E+00  9.9955E-01  1.3386E+00  2.9209E+00  1.4418E+00  9.7660E-01  9.2355E-01
             9.2655E-01
 PARAMETER:  1.2063E-01 -6.1134E-01  5.9913E-01  4.4313E-01  9.9550E-02  3.9166E-01  1.1719E+00  4.6590E-01  7.6326E-02  2.0467E-02
             2.3709E-02
 GRADIENT:   6.0069E+02  8.8627E+01  1.2996E+01  7.6788E+02  6.7584E+00  3.7414E+02  1.5251E+02  2.6258E+00  1.3165E+01  7.1232E-01
             9.9351E-01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -2059.18756961660        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1939
 NPARAMETR:  1.0208E+00  4.9099E-01  1.6473E+00  1.4094E+00  9.9955E-01  1.3386E+00  2.9209E+00  1.4418E+00  9.7660E-01  9.2355E-01
             9.2655E-01
 PARAMETER:  1.2063E-01 -6.1134E-01  5.9913E-01  4.4313E-01  9.9550E-02  3.9166E-01  1.1719E+00  4.6590E-01  7.6326E-02  2.0467E-02
             2.3709E-02
 GRADIENT:   3.8064E-01  2.0342E-02 -1.2579E+04  4.0279E-01 -5.9671E-01 -1.7493E-02  3.2028E-02 -7.8185E-02  3.1079E-02 -1.1390E-02
            -6.0881E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1939
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3837E-04  3.0193E-02 -4.8807E-02 -3.0174E-02 -2.9189E-02
 SE:             2.9977E-02  2.0258E-02  1.7851E-02  2.3372E-02  1.9533E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8833E-01  1.3610E-01  6.2541E-03  1.9669E-01  1.3507E-01

 ETASHRINKSD(%)  1.0000E-10  3.2134E+01  4.0198E+01  2.1702E+01  3.4563E+01
 ETASHRINKVR(%)  1.0000E-10  5.3943E+01  6.4237E+01  3.8694E+01  5.7180E+01
 EBVSHRINKSD(%)  1.7851E-01  3.6889E+01  4.3514E+01  1.7609E+01  3.0794E+01
 EBVSHRINKVR(%)  3.5669E-01  6.0169E+01  6.8093E+01  3.2117E+01  5.2105E+01
 RELATIVEINF(%)  9.9360E+01  7.5341E+00  8.3912E+00  1.3513E+01  1.4010E+01
 EPSSHRINKSD(%)  3.4892E+01
 EPSSHRINKVR(%)  5.7609E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2059.1875696165980     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1140.2490364119253     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2059.188       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  4.91E-01  1.65E+00  1.41E+00  1.00E+00  1.34E+00  2.92E+00  1.44E+00  9.77E-01  9.24E-01  9.27E-01
 


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
 #CPUT: Total CPU Time in Seconds,       71.751
Stop Time:
Sun Oct 24 01:10:38 CDT 2021
