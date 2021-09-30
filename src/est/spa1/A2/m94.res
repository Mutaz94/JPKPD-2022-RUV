Wed Sep 29 23:49:01 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat94.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -909.735380648625        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2721E+02  6.4221E+01  4.4226E+01  5.8771E+01  2.7261E+02  4.7647E+01 -7.1223E+01 -3.3470E+01 -8.2971E+01 -1.3661E+02
            -2.0557E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1642.32729131521        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0175E+00  7.1658E-01  8.4622E-01  1.1350E+00  6.8045E-01  8.2976E-01  1.1795E+00  8.8992E-01  1.5461E+00  1.1032E+00
             2.3918E+00
 PARAMETER:  1.1739E-01 -2.3326E-01 -6.6981E-02  2.2662E-01 -2.8500E-01 -8.6619E-02  2.6507E-01 -1.6619E-02  5.3571E-01  1.9818E-01
             9.7205E-01
 GRADIENT:   1.0429E+02  1.5695E+01  2.4237E+01  1.3803E+01 -8.4138E+00 -4.7898E+01  7.1163E+00  1.2459E+01  7.4988E+01  2.0674E+01
            -3.5363E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1655.01595607151        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  1.0138E+00  4.5175E-01  5.0281E-01  1.2900E+00  3.9810E-01  7.8558E-01  1.3496E+00  6.6933E-01  1.3734E+00  9.3618E-01
             2.2504E+00
 PARAMETER:  1.1372E-01 -6.9463E-01 -5.8754E-01  3.5464E-01 -8.2105E-01 -1.4134E-01  3.9981E-01 -3.0147E-01  4.1732E-01  3.4053E-02
             9.1110E-01
 GRADIENT:   1.6273E+00  5.0725E+01  9.8977E+01  4.7164E+01 -1.5064E+02 -8.4494E+01  1.6491E+00  1.0107E+01  4.8377E+01  3.4932E+01
            -4.0427E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1672.27528584958        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      356
 NPARAMETR:  1.0121E+00  3.8343E-01  4.4333E-01  1.3065E+00  3.7215E-01  7.8585E-01  1.5021E+00  1.3076E-01  1.2608E+00  9.2567E-01
             2.2867E+00
 PARAMETER:  1.1205E-01 -8.5860E-01 -7.1344E-01  3.6734E-01 -8.8846E-01 -1.4099E-01  5.0690E-01 -1.9344E+00  3.3178E-01  2.2758E-02
             9.2710E-01
 GRADIENT:  -4.5959E+00  3.7719E+01  5.6325E+01  6.2760E+01 -9.9417E+01 -8.3311E+01  2.7823E+00  4.3093E-01  2.8764E+01  4.1107E+01
            -2.1155E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1697.58876950434        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  1.0095E+00  2.3295E-01  3.7969E-01  1.2685E+00  3.3010E-01  9.6368E-01  1.2928E+00  6.9639E-02  1.1252E+00  6.6345E-01
             2.3462E+00
 PARAMETER:  1.0949E-01 -1.3569E+00 -8.6841E-01  3.3786E-01 -1.0084E+00  6.3006E-02  3.5680E-01 -2.5644E+00  2.1800E-01 -3.1030E-01
             9.5281E-01
 GRADIENT:   3.9221E-01  4.7314E+00  9.6877E+00  1.3056E+00 -1.5468E+01  5.3165E+00 -4.3082E-01  2.1594E-02  1.5759E+00 -8.4675E-01
            -1.8400E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1699.60661440958        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.9704E-01  6.5467E-02  4.2206E-01  1.3691E+00  3.3720E-01  9.3625E-01  3.1706E+00  1.0000E-02  1.0426E+00  6.7697E-01
             2.3570E+00
 PARAMETER:  9.7037E-02 -2.6262E+00 -7.6261E-01  4.1418E-01 -9.8709E-01  3.4128E-02  1.2539E+00 -4.8453E+00  1.4174E-01 -2.9013E-01
             9.5739E-01
 GRADIENT:  -6.3847E+00  9.6366E-01  6.7325E+00  9.6519E+00 -1.1079E+01 -2.6694E+00  2.4287E-02  0.0000E+00 -2.8701E+00 -1.8492E+00
            -2.6821E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1699.93731203564        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.9626E-01  1.0831E-02  4.2712E-01  1.3883E+00  3.3554E-01  9.4127E-01  8.0789E+00  1.0000E-02  1.0340E+00  6.9031E-01
             2.3644E+00
 PARAMETER:  9.6256E-02 -4.4253E+00 -7.5069E-01  4.2806E-01 -9.9201E-01  3.9478E-02  2.1893E+00 -8.5241E+00  1.3344E-01 -2.7061E-01
             9.6052E-01
 GRADIENT:  -2.5362E-01  6.8895E-02  1.5320E+00  1.4149E+00 -2.5740E+00  1.5330E-01  2.0765E-02  0.0000E+00 -2.3906E-01  1.1601E-01
             3.9208E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1699.95098785095        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.9633E-01  1.0000E-02  4.2627E-01  1.3861E+00  3.3532E-01  9.4094E-01  3.6245E+00  1.0000E-02  1.0350E+00  6.9011E-01
             2.3634E+00
 PARAMETER:  9.6327E-02 -4.5850E+00 -7.5267E-01  4.2650E-01 -9.9266E-01  3.9126E-02  1.3877E+00 -8.8571E+00  1.3436E-01 -2.7090E-01
             9.6009E-01
 GRADIENT:   1.8500E-01  0.0000E+00  9.5662E-02 -1.2357E+00  1.2970E-01  3.2326E-02  3.4030E-03  0.0000E+00 -8.4161E-02  2.6548E-02
             4.6042E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1699.95110619423        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  9.9628E-01  1.0000E-02  4.2628E-01  1.3869E+00  3.3533E-01  9.4090E-01  2.9671E+00  1.0000E-02  1.0352E+00  6.8998E-01
             2.3634E+00
 PARAMETER:  9.6270E-02 -4.5850E+00 -7.5266E-01  4.2706E-01 -9.9265E-01  3.9078E-02  1.1876E+00 -8.8571E+00  1.3464E-01 -2.7109E-01
             9.6009E-01
 GRADIENT:   2.7162E-02  0.0000E+00 -4.5089E-02 -1.8466E-01  1.8809E-02  3.3346E-03  2.2557E-03  0.0000E+00 -1.5678E-02  3.8439E-03
             4.3927E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1699.95137784618        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1414
 NPARAMETR:  9.9625E-01  1.0000E-02  4.2646E-01  1.3874E+00  3.3542E-01  9.4088E-01  1.3466E+00  1.0000E-02  1.0354E+00  6.9002E-01
             2.3634E+00
 PARAMETER:  9.6244E-02 -4.5850E+00 -7.5225E-01  4.2742E-01 -9.9236E-01  3.9059E-02  3.9755E-01 -8.8571E+00  1.3477E-01 -2.7104E-01
             9.6008E-01
 GRADIENT:  -4.8598E-02  0.0000E+00 -8.1415E-02  3.3241E-01 -2.9662E-02 -8.7307E-03  4.5445E-04  0.0000E+00  2.4209E-02 -7.3720E-03
            -1.0738E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1699.95160556438        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1541
 NPARAMETR:  9.9627E-01  1.0000E-02  4.2679E-01  1.3873E+00  3.3561E-01  9.4090E-01  3.7412E-01  1.0000E-02  1.0352E+00  6.9026E-01
             2.3633E+00
 PARAMETER:  9.6268E-02 -4.5850E+00 -7.5146E-01  4.2736E-01 -9.9182E-01  3.9077E-02 -8.8318E-01 -8.8571E+00  1.3456E-01 -2.7069E-01
             9.6008E-01
 GRADIENT:   1.6304E-02  0.0000E+00  3.9281E-02 -1.1397E-01  3.0403E-02  3.2148E-03  3.5794E-05  0.0000E+00 -8.1239E-03  3.1975E-03
             3.6847E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1541
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.8157E-04 -7.1719E-05  6.9958E-05 -6.9484E-03 -5.3824E-03
 SE:             2.9268E-02  5.5080E-05  2.4252E-04  2.8245E-02  2.2390E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8142E-01  1.9288E-01  7.7299E-01  8.0568E-01  8.1003E-01

 ETASHRINKSD(%)  1.9486E+00  9.9815E+01  9.9188E+01  5.3749E+00  2.4990E+01
 ETASHRINKVR(%)  3.8592E+00  1.0000E+02  9.9993E+01  1.0461E+01  4.3736E+01
 EBVSHRINKSD(%)  1.9348E+00  9.9820E+01  9.9202E+01  4.5834E+00  2.4576E+01
 EBVSHRINKVR(%)  3.8321E+00  1.0000E+02  9.9994E+01  8.9567E+00  4.3112E+01
 RELATIVEINF(%)  8.5289E+01  4.3869E-05  4.1163E-04  2.3074E+01  2.8976E+00
 EPSSHRINKSD(%)  2.6956E+01
 EPSSHRINKVR(%)  4.6646E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1699.9516055643805     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -781.01307235970785     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.51
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1699.952       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.00E-02  4.27E-01  1.39E+00  3.36E-01  9.41E-01  3.74E-01  1.00E-02  1.04E+00  6.90E-01  2.36E+00
 


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
+        0.00E+00  1.93E+02
 
 TH 3
+       -3.31E+01  0.00E+00  3.79E+03
 
 TH 4
+       -1.58E+01  0.00E+00 -2.42E+02  4.92E+02
 
 TH 5
+        8.70E+01  0.00E+00 -6.22E+03 -2.14E+02  1.20E+04
 
 TH 6
+        3.20E+00  0.00E+00  7.57E+00 -6.22E+00 -3.65E+00  2.06E+02
 
 TH 7
+       -1.01E-02  0.00E+00 -1.58E-02 -2.73E-03  1.02E-02  1.37E-02 -4.83E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.07E+00  0.00E+00  5.41E+01 -6.07E+00  6.87E+00 -6.45E-01  2.28E-02  0.00E+00  1.52E+02
 
 TH10
+       -4.09E+00  0.00E+00 -5.30E+01  7.20E+00 -2.86E+01  1.97E+00 -7.62E-03  0.00E+00  3.43E-01  1.41E+02
 
 TH11
+       -1.62E+01  0.00E+00 -1.01E+01 -7.89E+00  1.83E+01  3.16E+00  1.65E-03  0.00E+00  5.49E+00  2.51E+01  8.29E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.27
 #CPUT: Total CPU Time in Seconds,       31.111
Stop Time:
Wed Sep 29 23:49:35 CDT 2021
