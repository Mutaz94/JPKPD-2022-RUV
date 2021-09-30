Wed Sep 29 12:39:24 CDT 2021
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
$DATA ../../../../data/spa/A2/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -815.153015599594        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6230E+02  6.0643E+01  4.8987E+01  9.4795E-01  9.0405E+01  4.3952E+01 -6.6108E+01 -4.9921E+00 -1.2061E+02 -6.7382E+01
            -1.4662E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1373.66029990671        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0559E+00  9.9916E-01  1.0427E+00  1.1607E+00  9.9430E-01  8.4801E-01  1.0307E+00  9.1230E-01  1.1174E+00  8.3843E-01
             3.4110E+00
 PARAMETER:  1.5440E-01  9.9162E-02  1.4177E-01  2.4906E-01  9.4285E-02 -6.4867E-02  1.3027E-01  8.2184E-03  2.1103E-01 -7.6225E-02
             1.3270E+00
 GRADIENT:   7.6668E+00  4.1715E+01 -5.4571E+00  7.8237E+01 -2.1157E+01 -2.6020E+01  7.0645E+00  7.0811E+00  1.1584E+01  1.9352E+01
             6.1537E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1384.33630172321        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0588E+00  7.3855E-01  6.3635E-01  1.2969E+00  6.3212E-01  9.2970E-01  1.2988E+00  2.7971E-01  1.1100E+00  3.9498E-01
             3.1401E+00
 PARAMETER:  1.5711E-01 -2.0307E-01 -3.5201E-01  3.5996E-01 -3.5867E-01  2.7105E-02  3.6141E-01 -1.1740E+00  2.0434E-01 -8.2892E-01
             1.2443E+00
 GRADIENT:   3.1901E+01  5.9497E+01  5.1208E+00  1.3344E+02 -2.6428E+01 -3.9504E-01  8.1869E+00  1.3155E+00  1.1052E+01  5.9672E+00
             3.5164E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1395.59946254762        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0361E+00  6.2753E-01  2.8184E-01  1.1536E+00  3.5980E-01  9.8865E-01  9.8117E-01  7.5521E-02  1.1169E+00  2.0229E-01
             2.6815E+00
 PARAMETER:  1.3545E-01 -3.6596E-01 -1.1664E+00  2.4288E-01 -9.2220E-01  8.8584E-02  8.0987E-02 -2.4833E+00  2.1053E-01 -1.4980E+00
             1.0864E+00
 GRADIENT:   6.2584E+00  1.0975E+02  8.5497E+01  1.1393E+02 -9.8034E+01  1.0157E+01 -2.0753E+01 -7.4076E-02 -1.0907E+01 -2.5610E+00
             1.4978E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1399.04799027818        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      346
 NPARAMETR:  1.0525E+00  5.8805E-01  2.7378E-01  1.1735E+00  3.4501E-01  9.6976E-01  1.1107E+00  5.4305E-02  1.0977E+00  2.1509E-01
             2.6450E+00
 PARAMETER:  1.5119E-01 -4.3094E-01 -1.1954E+00  2.6001E-01 -9.6417E-01  6.9295E-02  2.0496E-01 -2.8131E+00  1.9321E-01 -1.4367E+00
             1.0727E+00
 GRADIENT:  -3.0293E+01  8.6964E+01  5.4735E+01  1.0183E+02 -1.2720E+02 -2.0520E-01 -1.5523E+01 -4.1285E-02 -1.8426E+01 -2.4418E+00
            -4.6900E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1410.46899369487        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  1.0580E+00  4.3238E-01  2.4715E-01  1.1207E+00  3.0089E-01  9.7371E-01  1.4242E+00  1.0000E-02  1.1393E+00  3.3274E-01
             2.4906E+00
 PARAMETER:  1.5636E-01 -7.3845E-01 -1.2977E+00  2.1396E-01 -1.1010E+00  7.3354E-02  4.5362E-01 -5.0287E+00  2.3043E-01 -1.0004E+00
             1.0125E+00
 GRADIENT:  -7.5658E+00  7.0806E+00  7.0900E+00  2.2649E+00 -1.6922E+01  1.1751E+00 -1.5931E+00  0.0000E+00  2.7157E-01 -6.9124E-01
            -3.6538E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1411.09769697248        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      697
 NPARAMETR:  1.0623E+00  3.2440E-01  2.8935E-01  1.2051E+00  3.0426E-01  9.6463E-01  1.7267E+00  1.0000E-02  1.1033E+00  4.6634E-01
             2.5214E+00
 PARAMETER:  1.6041E-01 -1.0258E+00 -1.1401E+00  2.8659E-01 -1.0899E+00  6.3994E-02  6.4621E-01 -6.0487E+00  1.9830E-01 -6.6284E-01
             1.0248E+00
 GRADIENT:  -1.7433E+00  4.6779E+00  5.6270E+00  3.8762E+00 -1.3050E+01  1.1542E+00  5.9194E-02  0.0000E+00 -7.7312E-01  9.0147E-01
             3.3295E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1412.93971730502        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      874
 NPARAMETR:  1.0578E+00  1.5692E-01  3.9396E-01  1.3728E+00  3.4147E-01  9.4695E-01  2.6074E+00  1.0000E-02  1.0268E+00  5.6053E-01
             2.6228E+00
 PARAMETER:  1.5621E-01 -1.7520E+00 -8.3151E-01  4.1684E-01 -9.7449E-01  4.5495E-02  1.0583E+00 -8.3161E+00  1.2646E-01 -4.7887E-01
             1.0643E+00
 GRADIENT:   7.5610E+00  1.1125E+01  4.6768E+01  2.4604E+01 -6.8348E+01  9.4081E-02  2.4480E+00  0.0000E+00 -8.5510E+00  7.7162E-01
             3.7754E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1416.79386045540        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1049
 NPARAMETR:  1.0448E+00  4.4542E-02  3.5528E-01  1.3649E+00  3.1393E-01  9.4559E-01  4.0981E+00  1.0000E-02  1.0465E+00  5.5916E-01
             2.5603E+00
 PARAMETER:  1.4383E-01 -3.0113E+00 -9.3485E-01  4.1106E-01 -1.0586E+00  4.4054E-02  1.5105E+00 -1.3828E+01  1.4546E-01 -4.8133E-01
             1.0401E+00
 GRADIENT:  -1.9413E+00  8.8877E-01  4.7805E+00  4.0783E+00 -1.0529E+01 -8.5781E-02 -2.4693E-01  0.0000E+00  4.5745E-03  1.2992E-01
            -9.0032E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1417.41465823913        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1225            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0429E+00  1.0000E-02  3.6735E-01  1.3844E+00  3.1982E-01  9.4408E-01  8.1799E+00  1.0000E-02  1.0338E+00  5.5862E-01
             2.5701E+00
 PARAMETER:  1.4199E-01 -4.6587E+00 -9.0145E-01  4.2530E-01 -1.0400E+00  4.2459E-02  2.2017E+00 -2.0843E+01  1.3323E-01 -4.8228E-01
             1.0439E+00
 GRADIENT:   7.6619E+01  0.0000E+00  1.4881E+01  9.2178E+01  6.4501E+01  4.0176E+00  2.0865E-02  0.0000E+00  5.0255E+00  1.0170E+00
             8.4082E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1417.42995399841        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1404
 NPARAMETR:  1.0426E+00  1.0000E-02  3.6710E-01  1.3865E+00  3.1976E-01  9.4383E-01  1.1849E+01  1.0000E-02  1.0325E+00  5.5647E-01
             2.5734E+00
 PARAMETER:  1.4171E-01 -4.6587E+00 -9.0212E-01  4.2681E-01 -1.0402E+00  4.2187E-02  2.5723E+00 -2.0843E+01  1.3200E-01 -4.8614E-01
             1.0452E+00
 GRADIENT:  -2.4573E-01  0.0000E+00 -2.6459E+00  1.6145E+00  3.1801E+00 -2.4712E-02 -2.1456E-02  0.0000E+00 -4.0709E-02  7.8325E-02
             4.2770E-01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1417.43209474238        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1461
 NPARAMETR:  1.0427E+00  1.0000E-02  3.6719E-01  1.3860E+00  3.1972E-01  9.4387E-01  1.2158E+01  1.0000E-02  1.0326E+00  5.5646E-01
             2.5731E+00
 PARAMETER:  1.4179E-01 -4.6587E+00 -9.0187E-01  4.2644E-01 -1.0403E+00  4.2229E-02  2.5980E+00 -2.0843E+01  1.3203E-01 -4.8616E-01
             1.0451E+00
 GRADIENT:  -4.4184E-02  0.0000E+00 -1.8186E+00  8.7697E-01  2.2910E+00 -1.9762E-04  2.4365E-02  0.0000E+00 -1.0747E-02  4.7560E-02
             3.5829E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1461
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.3846E-04  7.1799E-04  8.6954E-05 -9.5468E-03 -2.9862E-03
 SE:             2.9043E-02  1.8388E-03  2.5857E-04  2.7625E-02  1.9591E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9070E-01  6.9619E-01  7.3665E-01  7.2966E-01  8.7885E-01

 ETASHRINKSD(%)  2.7030E+00  9.3840E+01  9.9134E+01  7.4515E+00  3.4369E+01
 ETASHRINKVR(%)  5.3328E+00  9.9621E+01  9.9992E+01  1.4348E+01  5.6925E+01
 EBVSHRINKSD(%)  2.6655E+00  9.4502E+01  9.9132E+01  6.6785E+00  3.4040E+01
 EBVSHRINKVR(%)  5.2599E+00  9.9698E+01  9.9992E+01  1.2911E+01  5.6493E+01
 RELATIVEINF(%)  7.7404E+01  3.4236E-02  2.4345E-04  1.7702E+01  1.1806E+00
 EPSSHRINKSD(%)  3.1958E+01
 EPSSHRINKVR(%)  5.3702E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1417.4320947423817     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -682.28126817864347     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1417.432       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  3.67E-01  1.39E+00  3.20E-01  9.44E-01  1.22E+01  1.00E-02  1.03E+00  5.56E-01  2.57E+00
 


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
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.05E+01  0.00E+00  6.12E+03
 
 TH 4
+       -3.05E+01  0.00E+00 -2.99E+02  4.74E+02
 
 TH 5
+        1.83E+02  0.00E+00 -9.11E+03 -2.17E+02  1.53E+04
 
 TH 6
+       -1.46E-01  0.00E+00  1.72E+01 -8.18E+00 -2.96E+00  1.99E+02
 
 TH 7
+       -4.18E-01  0.00E+00 -1.55E+02 -9.04E+01  1.49E+02  3.93E-01  1.16E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.28E+00  0.00E+00  3.96E+01 -8.58E+00  3.00E+01 -1.47E+00  9.16E-01  0.00E+00  1.42E+02
 
 TH10
+       -9.38E+00  0.00E+00 -8.82E+01  8.86E-01  1.32E+02  3.94E+00  2.42E+00  0.00E+00  9.03E-02  1.15E+02
 
 TH11
+       -1.47E+01  0.00E+00 -4.90E+00 -5.69E+00 -9.41E+00  3.58E+00 -1.98E+01  0.00E+00  6.88E+00  3.16E+01  4.44E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.559
Stop Time:
Wed Sep 29 12:39:50 CDT 2021
