Sat Sep 25 12:56:01 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1587.36209020719        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4942E+02 -4.6741E+01 -3.1788E+01 -8.4810E+00  3.9894E+01 -6.9583E+00 -4.1300E+01 -1.8483E+00 -2.3254E+01 -5.9271E-01
            -2.4167E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1601.70183481518        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.3676E-01  1.1390E+00  1.0712E+00  9.1889E-01  1.0671E+00  1.0061E+00  1.3967E+00  1.0415E+00  1.0847E+00  9.5640E-01
             1.0473E+00
 PARAMETER:  3.4668E-02  2.3016E-01  1.6878E-01  1.5409E-02  1.6492E-01  1.0612E-01  4.3412E-01  1.4067E-01  1.8129E-01  5.5420E-02
             1.4619E-01
 GRADIENT:   2.9340E+00  3.3757E+00  6.5215E+00 -4.3898E+00  4.2043E+00  7.5507E-01  6.8961E+00 -4.9751E+00  9.0355E+00 -5.1012E-01
            -2.6501E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1602.11126726594        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:      258
 NPARAMETR:  9.3711E-01  1.1370E+00  1.0766E+00  9.2547E-01  1.0612E+00  1.0247E+00  1.4037E+00  1.1404E+00  1.0731E+00  9.1626E-01
             1.0534E+00
 PARAMETER:  3.5043E-02  2.2840E-01  1.7376E-01  2.2551E-02  1.5940E-01  1.2444E-01  4.3911E-01  2.3135E-01  1.7057E-01  1.2543E-02
             1.5203E-01
 GRADIENT:  -3.3947E+01 -2.8396E+00  3.3720E+00 -5.4456E+00  3.1281E+00  1.3005E+00  3.9157E+00 -2.3388E+00  7.0737E+00 -3.7446E+00
            -5.2514E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1602.53998341269        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  9.4181E-01  1.0194E+00  1.2589E+00  1.0155E+00  1.0722E+00  1.0216E+00  1.4941E+00  1.3172E+00  9.8954E-01  9.4101E-01
             1.0567E+00
 PARAMETER:  4.0046E-02  1.1923E-01  3.3026E-01  1.1543E-01  1.6972E-01  1.2133E-01  5.0151E-01  3.7548E-01  8.9481E-02  3.9203E-02
             1.5511E-01
 GRADIENT:  -2.1269E+01  3.2034E+00  3.1698E+00  2.8174E+00  1.7174E+00  8.4525E-01  1.1778E+00 -1.8861E+00  2.2122E+00 -3.7480E+00
             3.2661E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1602.62462053934        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      629
 NPARAMETR:  9.4401E-01  1.0159E+00  1.2686E+00  1.0177E+00  1.0742E+00  1.0191E+00  1.4953E+00  1.3322E+00  9.8536E-01  9.5054E-01
             1.0565E+00
 PARAMETER:  4.2387E-02  1.1580E-01  3.3791E-01  1.1751E-01  1.7160E-01  1.1891E-01  5.0229E-01  3.8681E-01  8.5249E-02  4.9279E-02
             1.5497E-01
 GRADIENT:  -1.6319E+01  2.7354E+00  2.8479E+00  2.1425E+00  6.8755E-01  1.5072E-02  1.0413E+00 -1.5287E+00  1.8087E+00 -2.7499E+00
             4.4456E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1602.69317788590        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      779
 NPARAMETR:  9.5490E-01  1.0159E+00  1.2688E+00  1.0177E+00  1.0741E+00  1.0173E+00  1.4955E+00  1.3324E+00  9.6378E-01  9.9010E-01
             1.0566E+00
 PARAMETER:  5.3855E-02  1.1578E-01  3.3808E-01  1.1754E-01  1.7150E-01  1.1711E-01  5.0243E-01  3.8700E-01  6.3111E-02  9.0052E-02
             1.5504E-01
 GRADIENT:   4.6171E+01  6.8703E+00  3.3207E+00  9.5049E+00 -4.6575E+00  6.3641E+00  3.7744E+00 -8.9171E-01 -6.0560E-02  2.2588E+00
             1.2638E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1602.71781167457        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:      953             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5075E-01  1.0152E+00  1.2688E+00  1.0176E+00  1.0747E+00  1.0147E+00  1.4954E+00  1.3324E+00  9.6410E-01  9.8613E-01
             1.0566E+00
 PARAMETER:  4.9499E-02  1.1507E-01  3.3807E-01  1.1745E-01  1.7202E-01  1.1459E-01  5.0242E-01  3.8701E-01  6.3441E-02  8.6032E-02
             1.5503E-01
 GRADIENT:   3.6905E+01  6.3314E+00  3.0403E+00  9.2366E+00 -3.2906E+00  5.3235E+00  3.6911E+00 -9.3860E-01 -2.6107E-02  1.6876E+00
             1.2175E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1602.73133811798        NO. OF FUNC. EVALS.: 207
 CUMULATIVE NO. OF FUNC. EVALS.:     1160             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5078E-01  1.0151E+00  1.2689E+00  1.0176E+00  1.0747E+00  1.0187E+00  1.4956E+00  1.3327E+00  9.6408E-01  9.7211E-01
             1.0566E+00
 PARAMETER:  4.9533E-02  1.1503E-01  3.3815E-01  1.1747E-01  1.7203E-01  1.1849E-01  5.0250E-01  3.8717E-01  6.3416E-02  7.1718E-02
             1.5503E-01
 GRADIENT:   3.6988E+01  6.3492E+00  3.2794E+00  9.5456E+00 -1.2895E+00  7.0730E+00  3.5533E+00 -1.2064E+00 -1.1843E-01 -1.0007E-01
             8.5405E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1602.73380013864        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1327
 NPARAMETR:  9.5079E-01  1.0151E+00  1.2678E+00  1.0176E+00  1.0747E+00  1.0184E+00  1.4956E+00  1.3326E+00  9.6408E-01  9.7305E-01
             1.0566E+00
 PARAMETER:  4.9534E-02  1.1503E-01  3.3728E-01  1.1747E-01  1.7203E-01  1.1823E-01  5.0251E-01  3.8717E-01  6.3414E-02  7.2682E-02
             1.5503E-01
 GRADIENT:   3.6980E+01  6.2470E+00  3.0074E+00  9.6957E+00 -1.0512E+00  6.9534E+00  3.5613E+00 -1.1175E+00 -9.0973E-02  3.3411E-02
             8.9640E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1602.73403453988        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:     1409
 NPARAMETR:  9.5079E-01  1.0151E+00  1.2677E+00  1.0176E+00  1.0747E+00  1.0183E+00  1.4956E+00  1.3326E+00  9.6408E-01  9.7305E-01
             1.0566E+00
 PARAMETER:  4.9534E-02  1.1503E-01  3.3718E-01  1.1747E-01  1.7203E-01  1.1802E-01  5.0251E-01  3.8717E-01  6.3414E-02  7.2682E-02
             1.5503E-01
 GRADIENT:  -1.1468E+06 -1.9939E+06 -6.8017E+05 -1.9524E+06  1.3332E+06 -1.6310E-01  5.3562E+00  5.9229E+05 -1.4253E+01 -5.6256E+00
            -1.4797E+06
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.2         8.5         3.3         8.8         9.2
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1409
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.7406E-04 -1.1478E-03 -3.9377E-02 -6.9168E-03 -3.4981E-02
 SE:             2.9847E-02  2.2077E-02  1.4914E-02  2.1884E-02  2.0349E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7397E-01  9.5854E-01  8.2823E-03  7.5195E-01  8.5610E-02

 ETASHRINKSD(%)  7.5799E-03  2.6039E+01  5.0037E+01  2.6687E+01  3.1828E+01
 ETASHRINKVR(%)  1.5159E-02  4.5298E+01  7.5037E+01  4.6252E+01  5.3526E+01
 EBVSHRINKSD(%)  4.7155E-01  2.6559E+01  5.4716E+01  2.7663E+01  2.8693E+01
 EBVSHRINKVR(%)  9.4088E-01  4.6064E+01  7.9494E+01  4.7674E+01  4.9154E+01
 RELATIVEINF(%)  9.8538E+01  3.1049E+00  3.0763E+00  3.0528E+00  1.2497E+01
 EPSSHRINKSD(%)  4.4726E+01
 EPSSHRINKVR(%)  6.9448E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1602.7340345398823     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -867.58320797614408     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1602.734       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.51E-01  1.02E+00  1.27E+00  1.02E+00  1.07E+00  1.02E+00  1.50E+00  1.33E+00  9.64E-01  9.73E-01  1.06E+00
 


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
+        6.34E+09
 
 TH 2
+       -6.43E+03  4.21E+09
 
 TH 3
+       -1.75E+03  3.40E+03  3.14E+08
 
 TH 4
+       -6.28E+03  1.23E+04  1.85E+04  4.01E+09
 
 TH 5
+        4.06E+03 -7.87E+03 -9.41E+03 -4.28E+04  1.68E+09
 
 TH 6
+       -5.02E+09 -9.77E+03 -2.67E+03 -9.55E+03  6.17E+03  3.97E+09
 
 TH 7
+       -8.02E+08 -7.30E+03 -2.00E+03 -7.16E+03  4.62E+03 -1.02E+00  1.02E+08
 
 TH 8
+        1.45E+03  2.83E+05  7.74E+04  2.77E+05 -1.79E+05  2.21E+03  1.66E+03  2.15E+08
 
 TH 9
+        3.64E+04  2.96E+04  8.07E+03  2.89E+04 -1.87E+04  5.14E+00  4.62E+03 -6.68E+03  1.23E+10
 
 TH10
+        1.51E+04  1.23E+04  3.36E+03  1.20E+04 -7.86E+03 -3.44E+00  7.84E+08 -2.78E+03  6.11E+09  1.21E+10
 
 TH11
+       -4.59E+03  1.43E+06  3.92E+05  1.40E+06 -9.06E+05 -6.97E+03 -4.66E+08 -3.25E+05  2.11E+04  8.80E+03  2.14E+09
 
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
 #CPUT: Total CPU Time in Seconds,       26.659
Stop Time:
Sat Sep 25 12:56:29 CDT 2021
