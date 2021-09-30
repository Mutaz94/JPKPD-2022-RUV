Thu Sep 30 06:28:57 CDT 2021
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
$DATA ../../../../data/spa2/A3/dat48.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE: -0.269930910677558        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6990E+02  1.3962E+02  3.4833E+02 -1.8596E+01  2.0504E+02  1.3907E+01 -9.8455E+01 -4.3028E+02 -4.8190E+01 -1.1987E+02
            -4.1228E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1603.40646318658        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9170E-01  9.7224E-01  7.7459E-01  1.1606E+00  8.7672E-01  9.2988E-01  1.0145E+00  1.1138E+00  8.8569E-01  9.0960E-01
             5.2775E+00
 PARAMETER:  9.1662E-02  7.1849E-02 -1.5542E-01  2.4897E-01 -3.1571E-02  2.7306E-02  1.1435E-01  2.0777E-01 -2.1388E-02  5.2460E-03
             1.7635E+00
 GRADIENT:  -9.1072E+01 -1.7372E+01 -3.8364E+01  4.2310E+01  3.5945E+01 -2.4301E+01  1.2779E+01  1.0855E+01  1.9529E+01  2.6111E+01
             4.5892E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1630.05163820186        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.7880E-01  6.2410E-01  3.5859E-01  1.3230E+00  4.1932E-01  9.7048E-01  1.2058E+00  5.1082E-01  8.4155E-01  4.1061E-01
             4.8202E+00
 PARAMETER:  7.8573E-02 -3.7145E-01 -9.2557E-01  3.7992E-01 -7.6913E-01  7.0035E-02  2.8717E-01 -5.7174E-01 -7.2516E-02 -7.9010E-01
             1.6728E+00
 GRADIENT:  -1.0974E+02  1.0410E+02 -4.2415E+00  1.7788E+02 -7.8851E+01 -2.5159E+01 -5.5763E-01  2.9359E+00 -7.0032E+00  2.4260E+00
             3.7971E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1706.24824114069        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  9.6792E-01  4.7419E-01  3.3925E-01  1.3495E+00  3.5510E-01  9.9651E-01  1.9207E+00  6.1721E-01  1.2902E+00  2.1067E-01
             3.5752E+00
 PARAMETER:  6.7390E-02 -6.4615E-01 -9.8101E-01  3.9974E-01 -9.3536E-01  9.6500E-02  7.5272E-01 -3.8254E-01  3.5480E-01 -1.4574E+00
             1.3740E+00
 GRADIENT:  -8.0326E+01  7.8271E+01  5.4690E+01  1.1097E+02 -1.2857E+02 -1.5693E+01  2.6753E+01 -5.7236E+00  4.9591E+01 -2.6841E+00
             1.8269E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1760.18640206860        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  9.9994E-01  3.4387E-01  2.4523E-01  1.1683E+00  2.7998E-01  1.0371E+00  1.4162E+00  1.0784E+00  1.0865E+00  2.4649E-01
             2.7039E+00
 PARAMETER:  9.9943E-02 -9.6749E-01 -1.3055E+00  2.5559E-01 -1.1730E+00  1.3647E-01  4.4799E-01  1.7547E-01  1.8298E-01 -1.3004E+00
             1.0947E+00
 GRADIENT:   1.4243E+01  2.2848E+01  3.2450E+01  9.2492E+00 -4.5227E+01 -4.6811E+00 -1.3992E+01 -7.8158E+00 -1.1682E+01 -2.5147E+00
            -3.1270E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1768.24599355405        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      630
 NPARAMETR:  9.9051E-01  2.6994E-01  1.7782E-01  1.1136E+00  2.2766E-01  1.0631E+00  1.6203E+00  1.3793E+00  1.4210E+00  1.4450E-01
             2.5562E+00
 PARAMETER:  9.0469E-02 -1.2096E+00 -1.6270E+00  2.0758E-01 -1.3799E+00  1.6120E-01  5.8264E-01  4.2156E-01  4.5136E-01 -1.8345E+00
             1.0385E+00
 GRADIENT:   4.2819E+00 -1.9963E+00  3.1231E+00  1.2855E+01 -1.4948E+00  4.3556E+00  1.3279E+00 -4.0606E+00  4.1553E+00  3.9574E-01
            -2.6087E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1768.62296814869        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  9.8839E-01  2.7220E-01  1.7610E-01  1.0950E+00  2.2716E-01  1.0475E+00  1.6007E+00  1.4128E+00  1.4011E+00  5.3324E-02
             2.5489E+00
 PARAMETER:  8.8325E-02 -1.2012E+00 -1.6367E+00  1.9077E-01 -1.3821E+00  1.4640E-01  5.7047E-01  4.4555E-01  4.3724E-01 -2.8314E+00
             1.0357E+00
 GRADIENT:   4.3627E+01  2.6906E+01  5.0834E+01  1.8422E+01  2.4041E+02  4.7053E+00  1.5234E+00  4.4421E-01  6.4946E+00  1.9914E-01
             7.3054E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1768.68140702969        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  9.8758E-01  2.7167E-01  1.7542E-01  1.0938E+00  2.2699E-01  1.0477E+00  1.6090E+00  1.4240E+00  1.4120E+00  3.0732E-02
             2.5508E+00
 PARAMETER:  8.7506E-02 -1.2032E+00 -1.6406E+00  1.8964E-01 -1.3828E+00  1.4663E-01  5.7559E-01  4.5349E-01  4.4499E-01 -3.3825E+00
             1.0364E+00
 GRADIENT:  -8.6474E-01  8.8832E-01 -1.0660E+00 -9.0380E-01 -3.7605E+00 -6.9413E-01 -8.0432E-01 -7.5855E-01 -1.5686E+00  3.0501E-02
            -2.8669E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1768.69813900564        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  9.8786E-01  2.7015E-01  1.7551E-01  1.0947E+00  2.2669E-01  1.0493E+00  1.6189E+00  1.4300E+00  1.4206E+00  1.7550E-02
             2.5556E+00
 PARAMETER:  8.7787E-02 -1.2088E+00 -1.6401E+00  1.9052E-01 -1.3842E+00  1.4817E-01  5.8178E-01  4.5769E-01  4.5105E-01 -3.9427E+00
             1.0383E+00
 GRADIENT:  -4.0386E-01 -5.0237E-01 -3.1932E-01 -1.1766E-01 -3.2142E+00 -5.0705E-02  4.4694E-02 -4.6048E-02 -1.7125E-01  1.1748E-02
            -5.1591E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1768.72222725132        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  9.8818E-01  2.7437E-01  1.7784E-01  1.0961E+00  2.2955E-01  1.0487E+00  1.6178E+00  1.4222E+00  1.4083E+00  1.0000E-02
             2.5502E+00
 PARAMETER:  8.8115E-02 -1.2025E+00 -1.6332E+00  1.9243E-01 -1.3793E+00  1.4836E-01  5.8034E-01  4.5674E-01  4.4644E-01 -5.6214E+00
             1.0392E+00
 GRADIENT:   1.2563E-03 -3.3177E-01 -3.1248E-01  4.3475E-02 -2.0971E+00  2.1172E-02 -8.2188E-03  5.5788E-02  7.2067E-02  0.0000E+00
             2.1258E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1249
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.3632E-04  1.5481E-02 -1.1006E-02 -9.3711E-03  7.9999E-04
 SE:             2.9318E-02  2.3341E-02  2.4475E-02  2.7661E-02  4.2707E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7724E-01  5.0717E-01  6.5294E-01  7.3477E-01  6.1041E-02

 ETASHRINKSD(%)  1.7799E+00  2.1806E+01  1.8006E+01  7.3324E+00  9.8569E+01
 ETASHRINKVR(%)  3.5280E+00  3.8857E+01  3.2770E+01  1.4127E+01  9.9980E+01
 EBVSHRINKSD(%)  1.9634E+00  2.1411E+01  1.6758E+01  6.7498E+00  9.8758E+01
 EBVSHRINKVR(%)  3.8883E+00  3.8237E+01  3.0708E+01  1.3044E+01  9.9985E+01
 RELATIVEINF(%)  9.6054E+01  3.2224E+01  2.5543E+01  7.1876E+01  3.6012E-03
 EPSSHRINKSD(%)  2.9117E+01
 EPSSHRINKVR(%)  4.9756E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1768.7222272513163     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -665.99598740570923     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.64
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1768.722       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  2.72E-01  1.77E-01  1.10E+00  2.28E-01  1.05E+00  1.62E+00  1.43E+00  1.41E+00  1.00E-02  2.56E+00
 


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
+        9.95E+02
 
 TH 2
+        1.10E+01  3.58E+03
 
 TH 3
+       -4.16E+01  2.07E+03  1.15E+04
 
 TH 4
+       -5.90E+00  9.33E+01 -1.19E+02  4.01E+02
 
 TH 5
+        6.43E+01 -7.06E+03 -1.46E+04 -7.99E+02  3.88E+04
 
 TH 6
+        3.60E+00 -7.74E+00  1.23E+01 -4.04E+00  1.95E+01  1.65E+02
 
 TH 7
+       -4.16E-01  8.11E+01 -1.19E+01 -1.35E+00 -2.70E+01  3.28E-01  3.31E+01
 
 TH 8
+        1.14E+00  2.17E+01 -2.17E+01 -2.20E+00  6.59E+01  2.06E+00  3.86E+00  4.38E+01
 
 TH 9
+        8.88E+00 -3.76E+01  6.12E+01 -1.21E+01  4.88E+02  1.60E+00  3.20E+00 -2.03E+00  6.50E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.77E+01 -2.65E+01 -8.32E+01 -5.17E-01  7.27E+01  1.03E+00  3.91E+00  9.22E+00  9.00E+00  0.00E+00  7.76E+01
 
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
 #CPUT: Total CPU Time in Seconds,       33.999
Stop Time:
Thu Sep 30 06:29:32 CDT 2021
