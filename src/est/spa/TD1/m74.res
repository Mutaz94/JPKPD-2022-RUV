Sat Sep 18 14:14:49 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat74.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1692.90779519264        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3907E+01 -5.4466E+01 -1.5806E+01 -7.1730E+01 -1.7983E+01  2.3580E+01 -1.0065E+01  1.4342E+01 -5.0253E+00  1.6566E+01
             2.2760E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1699.86566583657        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  9.9293E-01  1.0731E+00  1.1234E+00  1.0146E+00  1.0995E+00  9.4262E-01  1.1022E+00  8.8565E-01  1.0220E+00  9.5924E-01
             9.8367E-01
 PARAMETER:  9.2904E-02  1.7056E-01  2.1637E-01  1.1449E-01  1.9482E-01  4.0907E-02  1.9729E-01 -2.1439E-02  1.2171E-01  5.8390E-02
             8.3533E-02
 GRADIENT:  -2.8563E+00  2.3614E+00  7.8943E-01  1.6837E+00  4.6647E+00 -2.1669E+00 -7.6210E-01  3.5261E+00  1.7925E+00 -7.3741E+00
             9.2304E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1700.65997963824        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      290
 NPARAMETR:  9.9226E-01  8.9662E-01  1.1132E+00  1.1280E+00  1.0145E+00  9.5581E-01  1.2111E+00  5.6483E-01  9.9775E-01  9.9733E-01
             9.5571E-01
 PARAMETER:  9.2229E-02 -9.1274E-03  2.0727E-01  2.2048E-01  1.1438E-01  5.4801E-02  2.9150E-01 -4.7123E-01  9.7747E-02  9.7324E-02
             5.4700E-02
 GRADIENT:  -1.0838E+00  4.0569E+00  7.8405E+00  1.0215E+01 -6.5321E+00  3.7268E+00 -2.2422E+00 -1.4663E+00  1.1052E+01  6.6775E-01
            -2.4213E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1701.48735310344        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      466
 NPARAMETR:  9.9331E-01  8.3993E-01  1.0197E+00  1.1507E+00  9.4397E-01  9.4579E-01  1.3776E+00  4.8297E-01  8.9975E-01  9.2269E-01
             9.5667E-01
 PARAMETER:  9.3288E-02 -7.4433E-02  1.1947E-01  2.4039E-01  4.2342E-02  4.4262E-02  4.2031E-01 -6.2779E-01 -5.6405E-03  1.9543E-02
             5.5702E-02
 GRADIENT:   6.1529E-01  2.7271E+00  2.5727E+00  8.4597E-01 -2.7097E+00 -3.7366E-01  2.1057E-02 -3.5225E-01 -1.1442E+00 -1.8108E-01
            -2.8849E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1701.72987079001        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      643
 NPARAMETR:  9.8980E-01  6.0237E-01  1.1221E+00  1.3010E+00  9.0383E-01  9.4394E-01  1.7210E+00  6.1451E-01  8.4274E-01  9.2102E-01
             9.5564E-01
 PARAMETER:  8.9751E-02 -4.0688E-01  2.1522E-01  3.6313E-01 -1.1135E-03  4.2304E-02  6.4289E-01 -3.8693E-01 -7.1101E-02  1.7722E-02
             5.4622E-02
 GRADIENT:  -3.5966E-01  2.2737E+00 -3.3391E-02  4.9724E+00 -1.8844E+00  1.9609E-01  1.7215E-02  2.7676E-01  1.9956E-01  1.3778E-01
            -1.0689E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1701.89487774743        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  9.8687E-01  4.2226E-01  1.1838E+00  1.4079E+00  8.7742E-01  9.4220E-01  2.1546E+00  6.5892E-01  8.0706E-01  9.2808E-01
             9.5544E-01
 PARAMETER:  8.6779E-02 -7.6213E-01  2.6874E-01  4.4212E-01 -3.0769E-02  4.0458E-02  8.6759E-01 -3.1715E-01 -1.1435E-01  2.5367E-02
             5.4412E-02
 GRADIENT:  -3.5833E-01  2.8850E-01 -1.4388E+00  4.2839E-01  2.0592E+00  6.9917E-01  3.6277E-01 -1.2717E-01  1.7470E-01  6.0777E-02
            -3.2027E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1702.03680146277        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  9.8369E-01  2.4660E-01  1.2882E+00  1.5239E+00  8.6704E-01  9.3675E-01  2.9013E+00  7.8206E-01  7.7607E-01  9.3990E-01
             9.5634E-01
 PARAMETER:  8.3553E-02 -1.3000E+00  3.5325E-01  5.2130E-01 -4.2664E-02  3.4657E-02  1.1652E+00 -1.4583E-01 -1.5351E-01  3.8014E-02
             5.5357E-02
 GRADIENT:  -1.2974E-02  2.2717E+00  3.7136E+00  1.2896E+01 -6.7803E+00 -3.0437E-01  1.5799E-01 -1.9380E-01 -1.3214E+00  2.8872E-01
            -1.4897E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1702.46227539240        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1169
 NPARAMETR:  9.8006E-01  9.4444E-02  1.4300E+00  1.6201E+00  8.8270E-01  9.3412E-01  4.5965E+00  9.7175E-01  7.5891E-01  9.4897E-01
             9.5714E-01
 PARAMETER:  7.9855E-02 -2.2597E+00  4.5767E-01  5.8248E-01 -2.4775E-02  3.1853E-02  1.6253E+00  7.1342E-02 -1.7588E-01  4.7622E-02
             5.6198E-02
 GRADIENT:  -2.0140E+00  8.0170E-01  5.4891E-01  8.8488E+00 -3.6130E+00 -4.4167E-01  1.3244E-01  9.3865E-01  1.8756E+00 -5.3883E-01
             6.4970E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1702.64972236220        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1347
 NPARAMETR:  9.7995E-01  5.3191E-02  1.5607E+00  1.6573E+00  9.2029E-01  9.3481E-01  5.8527E+00  1.1048E+00  7.3753E-01  1.0050E+00
             9.5615E-01
 PARAMETER:  7.9743E-02 -2.8339E+00  5.4515E-01  6.0520E-01  1.6935E-02  3.2587E-02  1.8669E+00  1.9967E-01 -2.0445E-01  1.0501E-01
             5.5156E-02
 GRADIENT:  -8.5746E-01 -9.9962E-01  3.4419E-01  2.9282E+01 -8.1160E+00  2.3184E-01 -5.6639E+00  3.0243E+00 -6.2170E-01  5.3297E+00
             9.0560E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1703.66503905659        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1529
 NPARAMETR:  9.7876E-01  1.2867E-02  1.4727E+00  1.6732E+00  8.8514E-01  9.3476E-01  1.0473E+01  9.6713E-01  7.4051E-01  9.8953E-01
             9.5911E-01
 PARAMETER:  7.8533E-02 -4.2531E+00  4.8707E-01  6.1477E-01 -2.2008E-02  3.2539E-02  2.4488E+00  6.6582E-02 -2.0041E-01  8.9471E-02
             5.8251E-02
 GRADIENT:  -1.7863E+00 -1.6301E+00  1.6276E-01  3.4286E+01 -3.0173E+00  4.6902E-01 -4.8965E+00  5.0038E-01  3.8857E+00  6.1791E+00
             1.6751E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1704.36426332521        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1712
 NPARAMETR:  9.7914E-01  1.0000E-02  1.3683E+00  1.6525E+00  8.4288E-01  9.3402E-01  1.1966E+01  8.9496E-01  7.3711E-01  9.2255E-01
             9.5690E-01
 PARAMETER:  7.8922E-02 -4.7293E+00  4.1356E-01  6.0227E-01 -7.0933E-02  3.1745E-02  2.5821E+00 -1.0979E-02 -2.0501E-01  1.9384E-02
             5.5945E-02
 GRADIENT:  -1.7031E-01  0.0000E+00  1.8071E-02  2.6232E-01  1.2922E-01  1.4707E-01 -3.0960E-01 -1.5381E-02  3.0028E-03  1.0872E-01
            -2.9641E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1704.36444770111        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1811
 NPARAMETR:  9.7921E-01  1.0000E-02  1.3686E+00  1.6528E+00  8.4274E-01  9.3367E-01  1.1961E+01  8.9667E-01  7.3729E-01  9.2211E-01
             9.5695E-01
 PARAMETER:  7.8986E-02 -4.7293E+00  4.1368E-01  6.0229E-01 -7.1110E-02  3.1352E-02  2.5824E+00 -8.9653E-03 -2.0478E-01  1.8927E-02
             5.6052E-02
 GRADIENT:  -1.0558E-02  0.0000E+00 -2.0229E+03 -1.3838E+03 -1.0933E-01 -1.1161E-02  3.2067E+02  2.6595E-02 -6.3658E-02  7.7964E-02
             3.2903E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1811
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3384E-04  9.4566E-03 -2.2246E-02 -8.1007E-03 -2.9562E-02
 SE:             2.9834E-02  6.3199E-03  1.5209E-02  2.8937E-02  2.1643E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9642E-01  1.3457E-01  1.4356E-01  7.7952E-01  1.7197E-01

 ETASHRINKSD(%)  5.2505E-02  7.8828E+01  4.9047E+01  3.0586E+00  2.7494E+01
 ETASHRINKVR(%)  1.0498E-01  9.5517E+01  7.4038E+01  6.0236E+00  4.7428E+01
 EBVSHRINKSD(%)  4.3260E-01  8.4437E+01  5.0800E+01  3.0901E+00  2.4412E+01
 EBVSHRINKVR(%)  8.6333E-01  9.7578E+01  7.5794E+01  6.0847E+00  4.2864E+01
 RELATIVEINF(%)  9.9086E+01  1.9841E+00  3.9758E+00  7.0951E+01  9.3839E+00
 EPSSHRINKSD(%)  4.3936E+01
 EPSSHRINKVR(%)  6.8568E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1704.3644477011107     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -969.21362113737257     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1704.364       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.00E-02  1.37E+00  1.65E+00  8.43E-01  9.34E-01  1.20E+01  8.97E-01  7.37E-01  9.22E-01  9.57E-01
 


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
+        1.32E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        3.08E+01  0.00E+00  6.53E+05
 
 TH 4
+        4.54E+00  0.00E+00 -3.31E+02  2.11E+05
 
 TH 5
+        5.83E+01  0.00E+00 -4.39E+06  7.45E+02  4.10E+03
 
 TH 6
+        1.36E+01  0.00E+00 -2.96E+01 -1.46E+01 -6.83E+01  2.12E+02
 
 TH 7
+       -4.38E-01  0.00E+00  8.49E+00  3.21E+01 -2.33E+01  3.07E-01  2.18E+02
 
 TH 8
+        5.03E-01  0.00E+00 -4.30E+02 -1.77E+02 -7.44E+02  8.27E+00  5.00E+00  2.13E+02
 
 TH 9
+       -1.04E+02  0.00E+00  2.45E+06 -1.71E+02 -1.64E+07  6.12E+01  2.99E+00  1.07E+03  9.18E+06
 
 TH10
+       -3.88E+01  0.00E+00  4.01E+06 -5.95E+02 -2.61E+03  1.76E+01  1.69E+01  6.73E+02  1.50E+07  2.34E+03
 
 TH11
+       -9.00E+00  0.00E+00 -2.73E+02 -1.27E+02 -4.92E+02  9.09E+01  3.30E+00  1.45E+02  6.87E+02  4.65E+02  2.87E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.093
Stop Time:
Sat Sep 18 14:15:22 CDT 2021
