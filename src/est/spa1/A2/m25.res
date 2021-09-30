Wed Sep 29 23:11:31 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1579.87271116226        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2517E+02 -3.4845E+00  1.7060E+01 -2.4243E+01  2.8414E+01  5.6646E+01 -5.4887E+00  1.3548E+01 -2.4906E+00 -7.0752E+00
            -1.0633E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1793.86143205421        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.1851E+00  1.0638E+00  1.1379E+00  1.0840E+00  1.0506E+00  1.1472E+00  9.6497E-01  8.3784E-01  9.7737E-01  7.8546E-01
             2.5913E+00
 PARAMETER:  2.6979E-01  1.6185E-01  2.2917E-01  1.8066E-01  1.4940E-01  2.3735E-01  6.4338E-02 -7.6927E-02  7.7109E-02 -1.4148E-01
             1.0521E+00
 GRADIENT:   4.2984E+02  3.2632E+01  4.6699E+00  4.5940E+01 -2.3241E+01  4.0582E+01  4.9534E+00  4.5607E+00  4.2103E+00  1.2983E+01
             2.1182E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1829.95460597395        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.0649E+00  9.2228E-01  4.2573E-01  1.0974E+00  6.0246E-01  9.4873E-01  1.2852E+00  9.2987E-02  1.1494E+00  5.0899E-01
             2.1573E+00
 PARAMETER:  1.6289E-01  1.9092E-02 -7.5394E-01  1.9292E-01 -4.0674E-01  4.7366E-02  3.5088E-01 -2.2753E+00  2.3928E-01 -5.7532E-01
             8.6888E-01
 GRADIENT:   2.6384E+02 -1.3682E+01 -9.4737E+01  1.1412E+02  1.0928E+02  9.3240E+00  3.3836E+01  2.2499E-01  5.2676E+01  1.1969E+01
             1.2286E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1863.57556686526        NO. OF FUNC. EVALS.: 113
 CUMULATIVE NO. OF FUNC. EVALS.:      278
 NPARAMETR:  9.8772E-01  7.6605E-01  6.8464E-01  1.1664E+00  6.9080E-01  8.9987E-01  1.2929E+00  1.4897E-01  8.9980E-01  5.6784E-01
             1.9165E+00
 PARAMETER:  8.7641E-02 -1.6651E-01 -2.7886E-01  2.5393E-01 -2.6990E-01 -5.5105E-03  3.5685E-01 -1.8040E+00 -5.5827E-03 -4.6592E-01
             7.5048E-01
 GRADIENT:  -3.7833E+01  1.7480E+00  4.0579E+00 -1.3995E+01 -7.1290E+00 -7.3387E+00 -1.1398E+00  3.2410E-01  3.6218E+00 -2.3634E+00
             1.8067E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1865.98101201677        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      453
 NPARAMETR:  1.0024E+00  5.4059E-01  7.4758E-01  1.3145E+00  6.5473E-01  9.1259E-01  1.6058E+00  5.2436E-02  8.3586E-01  7.2182E-01
             1.8485E+00
 PARAMETER:  1.0240E-01 -5.1509E-01 -1.9092E-01  3.7344E-01 -3.2353E-01  8.5349E-03  5.7364E-01 -2.8482E+00 -7.9292E-02 -2.2598E-01
             7.1435E-01
 GRADIENT:   1.3253E+01  6.0838E+00 -3.5194E+00  1.3253E+01 -5.3155E-02 -3.9481E-01 -1.8321E-01  4.8203E-02 -4.3832E-01  2.5511E+00
            -2.6430E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1867.23723354278        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      628
 NPARAMETR:  9.8619E-01  2.6114E-01  8.3080E-01  1.4826E+00  6.2189E-01  9.1005E-01  2.2922E+00  1.0000E-02  8.0015E-01  8.1616E-01
             1.8607E+00
 PARAMETER:  8.6089E-02 -1.2427E+00 -8.5365E-02  4.9377E-01 -3.7500E-01  5.7412E-03  9.2953E-01 -5.0939E+00 -1.2295E-01 -1.0314E-01
             7.2096E-01
 GRADIENT:  -1.0967E+01  8.4502E+00  1.9460E+01  2.1147E+01 -3.0274E+01  1.0042E+00  2.6174E+00  0.0000E+00  4.5877E-01  3.8205E+00
             4.5738E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1868.71303479008        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      804
 NPARAMETR:  9.8678E-01  1.0694E-01  7.8778E-01  1.5499E+00  5.7366E-01  9.0165E-01  3.0984E+00  1.0000E-02  7.7920E-01  7.9001E-01
             1.8427E+00
 PARAMETER:  8.6693E-02 -2.1355E+00 -1.3854E-01  5.3819E-01 -4.5571E-01 -3.5309E-03  1.2309E+00 -8.2360E+00 -1.4949E-01 -1.3571E-01
             7.1122E-01
 GRADIENT:   3.4139E+00  1.7707E+00  7.2830E+00  2.0550E+01 -1.3385E+01 -1.2955E+00 -6.6482E-01  0.0000E+00 -1.2432E+00 -3.7334E-01
            -2.9063E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1869.33911041490        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  9.8068E-01  1.7541E-02  7.9030E-01  1.5922E+00  5.6123E-01  9.0317E-01  6.4054E+00  1.0000E-02  7.6581E-01  8.0037E-01
             1.8459E+00
 PARAMETER:  8.0491E-02 -3.9432E+00 -1.3535E-01  5.6510E-01 -4.7763E-01 -1.8438E-03  1.9571E+00 -1.4996E+01 -1.6682E-01 -1.2268E-01
             7.1298E-01
 GRADIENT:  -4.6735E+00  2.3700E-01  2.4166E+00  1.4141E+01 -5.5595E+00  1.9054E-01 -8.6787E-02  0.0000E+00 -1.4002E+00  8.8679E-03
             3.0987E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1869.46254003940        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1168             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8208E-01  1.0000E-02  7.8296E-01  1.5856E+00  5.5744E-01  9.0253E-01  8.8926E+00  1.0000E-02  7.6719E-01  7.9998E-01
             1.8437E+00
 PARAMETER:  8.1922E-02 -4.6288E+00 -1.4467E-01  5.6099E-01 -4.8439E-01 -2.5546E-03  2.2852E+00 -1.7626E+01 -1.6502E-01 -1.2317E-01
             7.1177E-01
 GRADIENT:   1.1235E+02  0.0000E+00  2.5027E+00  2.9355E+02  1.9373E+01  8.0426E+00  8.2055E-02  0.0000E+00  6.0728E+00  3.6335E-01
             7.2732E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1869.46572722840        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1345
 NPARAMETR:  9.8187E-01  1.0000E-02  7.8310E-01  1.5891E+00  5.5742E-01  9.0237E-01  1.1358E+01  1.0000E-02  7.6788E-01  7.9920E-01
             1.8435E+00
 PARAMETER:  8.1707E-02 -4.6288E+00 -1.4450E-01  5.6318E-01 -4.8444E-01 -2.7303E-03  2.5299E+00 -1.7626E+01 -1.6412E-01 -1.2414E-01
             7.1167E-01
 GRADIENT:  -2.6162E-01  0.0000E+00 -2.1631E-01  2.8574E+00 -3.9162E-01 -1.4082E-02  1.4720E-02  0.0000E+00  1.1231E-01 -2.6654E-02
            -9.0955E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1869.46572722840        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1367
 NPARAMETR:  9.8187E-01  1.0000E-02  7.8310E-01  1.5891E+00  5.5742E-01  9.0237E-01  1.1358E+01  1.0000E-02  7.6788E-01  7.9920E-01
             1.8435E+00
 PARAMETER:  8.1707E-02 -4.6288E+00 -1.4450E-01  5.6318E-01 -4.8444E-01 -2.7303E-03  2.5299E+00 -1.7626E+01 -1.6412E-01 -1.2414E-01
             7.1167E-01
 GRADIENT:  -2.6162E-01  0.0000E+00 -2.1631E-01  2.8574E+00 -3.9162E-01 -1.4082E-02  1.4720E-02  0.0000E+00  1.1231E-01 -2.6654E-02
            -9.0955E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1367
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0349E-03  1.1668E-03 -3.2388E-05 -8.8724E-03 -1.4952E-02
 SE:             2.9500E-02  2.0259E-03  2.1071E-04  2.8603E-02  2.2925E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7201E-01  5.6465E-01  8.7784E-01  7.5642E-01  5.1425E-01

 ETASHRINKSD(%)  1.1704E+00  9.3213E+01  9.9294E+01  4.1763E+00  2.3200E+01
 ETASHRINKVR(%)  2.3271E+00  9.9539E+01  9.9995E+01  8.1782E+00  4.1018E+01
 EBVSHRINKSD(%)  1.2988E+00  9.3941E+01  9.9265E+01  4.1678E+00  2.2985E+01
 EBVSHRINKVR(%)  2.5807E+00  9.9633E+01  9.9995E+01  8.1620E+00  4.0687E+01
 RELATIVEINF(%)  8.6990E+01  1.6693E-02  3.7576E-04  6.1799E+00  3.3312E+00
 EPSSHRINKSD(%)  2.8235E+01
 EPSSHRINKVR(%)  4.8497E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1869.4657272283989     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -950.52719402372622     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1869.466       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.00E-02  7.83E-01  1.59E+00  5.57E-01  9.02E-01  1.14E+01  1.00E-02  7.68E-01  7.99E-01  1.84E+00
 


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
+        1.38E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.94E+01  0.00E+00  9.75E+02
 
 TH 4
+       -1.79E+01  0.00E+00 -1.31E+02  6.86E+02
 
 TH 5
+        3.78E+01  0.00E+00 -1.70E+03 -9.58E+01  3.32E+03
 
 TH 6
+        4.03E+00  0.00E+00  2.17E+00 -6.33E+00 -2.52E+00  2.31E+02
 
 TH 7
+        7.60E-03  0.00E+00 -2.38E-02 -6.00E-02  7.91E-02 -7.94E-04  2.71E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.91E+00  0.00E+00  2.84E+01 -1.16E+01 -1.63E+00 -1.92E+00 -3.05E-02  0.00E+00  2.87E+02
 
 TH10
+       -2.04E-01  0.00E+00 -8.25E+00  5.41E+00 -6.36E+01  1.11E+00 -3.15E-02  0.00E+00  2.14E+00  1.10E+02
 
 TH11
+       -1.55E+01  0.00E+00 -1.87E+01 -1.17E+01  2.17E+00  3.60E+00 -5.30E-03  0.00E+00  8.58E+00  2.94E+01  1.32E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       27.031
Stop Time:
Wed Sep 29 23:11:59 CDT 2021
