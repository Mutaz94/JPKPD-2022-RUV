Sat Sep 25 10:55:05 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat6.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1684.82668127736        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9380E+01 -9.1833E+01 -8.9478E+01 -1.8740E+01  1.3767E+02  6.3519E+00  1.6758E+00  1.2691E+01  1.2805E+00  5.3213E+00
             1.4030E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1697.47584913976        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      110
 NPARAMETR:  1.0063E+00  1.0460E+00  1.1917E+00  1.0134E+00  9.7055E-01  9.8918E-01  9.0819E-01  8.8441E-01  1.0549E+00  9.1314E-01
             1.0102E+00
 PARAMETER:  1.0631E-01  1.4500E-01  2.7539E-01  1.1331E-01  7.0110E-02  8.9118E-02  3.7037E-03 -2.2829E-02  1.5343E-01  9.1372E-03
             1.1015E-01
 GRADIENT:   1.1578E+01  5.5934E+00  9.2704E+00  1.5129E+00 -1.1415E+01 -2.7715E+00  4.8142E+00  1.9562E+00  3.2302E+00 -8.2704E+00
             1.4222E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1699.44830462215        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      287
 NPARAMETR:  1.0019E+00  1.1475E+00  1.1289E+00  9.5038E-01  9.8646E-01  9.8145E-01  5.4825E-01  6.6374E-01  1.1825E+00  1.0536E+00
             9.7922E-01
 PARAMETER:  1.0193E-01  2.3759E-01  2.2123E-01  4.9111E-02  8.6371E-02  8.1272E-02 -5.0102E-01 -3.0987E-01  2.6764E-01  1.5221E-01
             7.9002E-02
 GRADIENT:   4.7855E-01  2.3427E+01  9.4573E+00  1.8360E+01 -3.0190E+01 -6.3916E+00  1.7170E-01 -8.5961E-01 -2.5196E+00  3.4312E+00
             3.5567E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.84819858934        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  1.0019E+00  1.3102E+00  1.1496E+00  8.2385E-01  1.0771E+00  9.9667E-01  3.2361E-01  8.5456E-01  1.4092E+00  1.1073E+00
             9.6905E-01
 PARAMETER:  1.0192E-01  3.7020E-01  2.3938E-01 -9.3767E-02  1.7429E-01  9.6664E-02 -1.0282E+00 -5.7170E-02  4.4300E-01  2.0195E-01
             6.8564E-02
 GRADIENT:   1.0132E-01  4.6923E+00  4.1282E-01  3.7998E+00 -2.2321E+00 -3.1671E-01 -2.6283E-01 -3.4741E-01 -5.8699E-01  5.4393E-01
             4.4477E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1700.94070096957        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      642
 NPARAMETR:  1.0031E+00  1.4057E+00  1.0686E+00  7.5813E-01  1.0869E+00  9.9779E-01  4.4659E-01  8.2070E-01  1.4671E+00  1.0980E+00
             9.6868E-01
 PARAMETER:  1.0305E-01  4.4055E-01  1.6633E-01 -1.7690E-01  1.8330E-01  9.7785E-02 -7.0610E-01 -9.7594E-02  4.8327E-01  1.9345E-01
             6.8178E-02
 GRADIENT:   1.8673E+00 -2.2357E+00  3.1084E+00 -2.8466E+00 -8.3431E+00  2.0384E-02  5.0212E-01  1.4740E-02  1.7041E+00  9.6429E-01
             3.9591E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1701.29211982018        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  1.0015E+00  1.6007E+00  8.5536E-01  6.4436E-01  1.1149E+00  9.9825E-01  5.9770E-01  6.0493E-01  1.5789E+00  1.0739E+00
             9.6694E-01
 PARAMETER:  1.0149E-01  5.7042E-01 -5.6231E-02 -3.3949E-01  2.0874E-01  9.8244E-02 -4.1467E-01 -4.0264E-01  5.5673E-01  1.7129E-01
             6.6384E-02
 GRADIENT:  -3.6293E+00  8.8395E+00 -7.3911E-01  9.5216E+00  2.7371E+00 -1.3024E-01  1.5530E+00  6.1377E-01  2.5755E+00  2.2860E+00
             7.6614E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1701.68204703038        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      995
 NPARAMETR:  1.0029E+00  1.6845E+00  7.3873E-01  5.7638E-01  1.1115E+00  9.9740E-01  6.0233E-01  3.7429E-01  1.6297E+00  1.0366E+00
             9.6771E-01
 PARAMETER:  1.0293E-01  6.2148E-01 -2.0282E-01 -4.5098E-01  2.0571E-01  9.7399E-02 -4.0696E-01 -8.8273E-01  5.8842E-01  1.3596E-01
             6.7179E-02
 GRADIENT:  -1.0924E+00 -9.7056E-01 -6.9056E-02 -2.0250E+00 -1.7786E-01 -4.6785E-01 -1.4156E+00  3.1189E-01 -3.2088E+00 -5.3115E-01
            -1.9939E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.81918618723        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1171
 NPARAMETR:  1.0036E+00  1.7638E+00  6.7827E-01  5.2315E-01  1.1300E+00  9.9852E-01  6.0851E-01  1.9927E-01  1.7574E+00  1.0379E+00
             9.6805E-01
 PARAMETER:  1.0355E-01  6.6746E-01 -2.8820E-01 -5.4789E-01  2.2224E-01  9.8514E-02 -3.9675E-01 -1.5131E+00  6.6383E-01  1.3724E-01
             6.7528E-02
 GRADIENT:   5.6723E-02 -3.1851E-01 -2.3620E-01 -8.5506E-01 -7.2858E-01 -5.9233E-02  1.1039E-01  1.1389E-01  8.2972E-02 -2.8378E-01
             7.9052E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.88475001639        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1348
 NPARAMETR:  1.0035E+00  1.7342E+00  6.9031E-01  5.4402E-01  1.1196E+00  9.9872E-01  6.1423E-01  3.9356E-02  1.7084E+00  1.0359E+00
             9.6755E-01
 PARAMETER:  1.0345E-01  6.5057E-01 -2.7062E-01 -5.0877E-01  2.1298E-01  9.8718E-02 -3.8739E-01 -3.1351E+00  6.3558E-01  1.3523E-01
             6.7015E-02
 GRADIENT:  -1.4680E-01 -1.1512E-01 -1.9318E-01  7.9586E-02 -2.7792E-01 -3.0709E-03  2.2688E-02  3.8648E-03  5.8334E-02  2.6806E-01
             4.3328E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1701.88677646415        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1510
 NPARAMETR:  1.0035E+00  1.7354E+00  6.9071E-01  5.4320E-01  1.1206E+00  9.9873E-01  6.1405E-01  1.0000E-02  1.7103E+00  1.0353E+00
             9.6760E-01
 PARAMETER:  1.0352E-01  6.5124E-01 -2.7004E-01 -5.1027E-01  2.1383E-01  9.8725E-02 -3.8768E-01 -4.5824E+00  6.3665E-01  1.3474E-01
             6.7067E-02
 GRADIENT:   3.1474E-04 -2.1459E-03 -3.7005E-03  2.7490E-03  4.5400E-03  5.5146E-04 -1.2304E-03  0.0000E+00 -1.1548E-04  4.0617E-04
             1.4162E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1510
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6765E-04 -4.3491E-02 -2.2423E-04  2.6155E-02 -3.4188E-02
 SE:             2.9853E-02  1.9415E-02  9.8740E-05  2.4388E-02  2.4297E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9552E-01  2.5084E-02  2.3152E-02  2.8352E-01  1.5941E-01

 ETASHRINKSD(%)  1.0000E-10  3.4959E+01  9.9669E+01  1.8296E+01  1.8600E+01
 ETASHRINKVR(%)  1.0000E-10  5.7697E+01  9.9999E+01  3.3244E+01  3.3741E+01
 EBVSHRINKSD(%)  4.1486E-01  3.3888E+01  9.9710E+01  1.9015E+01  1.6831E+01
 EBVSHRINKVR(%)  8.2800E-01  5.6292E+01  9.9999E+01  3.4415E+01  3.0829E+01
 RELATIVEINF(%)  9.9127E+01  3.7054E+00  1.9332E-04  6.8659E+00  2.0255E+01
 EPSSHRINKSD(%)  4.4111E+01
 EPSSHRINKVR(%)  6.8764E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1701.8867764641484     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.73594990041022     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.17
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1701.887       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.74E+00  6.91E-01  5.43E-01  1.12E+00  9.99E-01  6.14E-01  1.00E-02  1.71E+00  1.04E+00  9.68E-01
 


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
+        1.09E+03
 
 TH 2
+       -7.28E+00  4.64E+02
 
 TH 3
+        6.09E+00  1.42E+02  2.17E+02
 
 TH 4
+       -1.12E+01  4.02E+02 -1.10E+02  8.55E+02
 
 TH 5
+       -3.34E+00 -2.20E+02 -2.26E+02  1.26E+02  5.69E+02
 
 TH 6
+        2.35E+00 -1.19E+00  1.27E+00 -1.56E+00 -1.80E+00  1.94E+02
 
 TH 7
+        2.59E+00 -2.32E+01  1.52E+01 -2.40E+01 -2.89E+01 -4.85E-01  1.09E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.05E+00 -2.52E+01 -1.66E+01  4.66E+01  2.28E+00 -2.84E-01  2.48E+01  0.00E+00  3.34E+01
 
 TH10
+       -2.00E+00 -1.60E+01 -3.47E+01 -3.07E+00 -5.28E+01  2.29E+00  2.26E+01  0.00E+00  1.94E+00  8.99E+01
 
 TH11
+       -7.55E+00 -2.15E+01 -2.22E+01  1.55E-01  8.89E+00  3.18E+00  1.14E+01  0.00E+00  4.12E+00  1.99E+01  2.33E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.163
Stop Time:
Sat Sep 25 10:55:31 CDT 2021
