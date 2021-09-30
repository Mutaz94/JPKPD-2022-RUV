Thu Sep 30 09:32:55 CDT 2021
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
$DATA ../../../../data/spa2/D/dat69.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20925.3626249203        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6886E+02  3.2656E+02  1.2545E-01  2.4498E+02  4.1930E+02 -2.1011E+03 -1.1084E+03 -3.7329E+01 -1.7737E+03 -7.9875E+02
            -4.0531E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -657.777784029318        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3309E+00  1.1802E+00  8.7440E-01  1.8196E+00  1.0828E+00  3.0838E+00  2.2966E+00  9.4995E-01  2.2888E+00  1.1652E+00
             1.2511E+01
 PARAMETER:  3.8584E-01  2.6572E-01 -3.4220E-02  6.9864E-01  1.7954E-01  1.2262E+00  9.3143E-01  4.8656E-02  9.2803E-01  2.5286E-01
             2.6266E+00
 GRADIENT:  -1.4199E+01 -2.0731E+01 -5.9148E+01  6.5527E+01  6.1457E+01  9.8412E+01 -1.5768E+01  4.7675E+00 -1.9309E+01  1.4282E+01
             1.9076E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -742.601320970957        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      200
 NPARAMETR:  1.2279E+00  1.7484E+00  2.7270E+00  1.5260E+00  1.8956E+00  4.1648E+00  4.8819E+00  6.8309E-01  4.0522E+00  1.5079E+00
             1.0973E+01
 PARAMETER:  3.0531E-01  6.5872E-01  1.1032E+00  5.2264E-01  7.3953E-01  1.5267E+00  1.6855E+00 -2.8113E-01  1.4993E+00  5.1074E-01
             2.4954E+00
 GRADIENT:  -2.6702E+01 -1.0331E+01 -2.5229E+01  3.3272E+01  1.4479E+01  9.7480E+01  1.7246E+01  1.8695E-01  4.8682E+01  2.0449E+01
             1.8580E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -818.865521077236        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.4116E+00  7.3773E-01  1.5487E+01  1.5348E+00  1.8673E+00  2.6465E+00  7.4108E+00  1.5469E+00  1.7203E+00  7.0066E-01
             9.8449E+00
 PARAMETER:  4.4472E-01 -2.0418E-01  2.8400E+00  5.2837E-01  7.2450E-01  1.0732E+00  2.1029E+00  5.3628E-01  6.4252E-01 -2.5573E-01
             2.3870E+00
 GRADIENT:  -1.9399E+00 -1.9278E-02  9.8981E-02 -1.0367E+01  3.7297E+00 -2.9039E+01  1.5416E+01  6.5997E-02  1.7082E+01  5.2235E+00
             9.4178E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -829.929161125254        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  1.4390E+00  7.1711E-01  1.6167E+01  1.4141E+00  1.8390E+00  2.9488E+00  6.7029E+00  1.1424E+00  1.3731E+00  3.6357E-01
             8.9416E+00
 PARAMETER:  4.6393E-01 -2.3253E-01  2.8830E+00  4.4652E-01  7.0923E-01  1.1814E+00  2.0025E+00  2.3314E-01  4.1706E-01 -9.1178E-01
             2.2907E+00
 GRADIENT:   8.0905E+00 -4.5403E+00  5.5209E-01  4.5931E+00  3.4034E+00  4.8857E+00 -5.6287E+00  3.0538E-02 -3.3303E-01  1.4822E+00
            -7.9595E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -834.323760017756        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  1.3590E+00  1.4869E+00  4.1055E+00  1.0241E+00  1.6455E+00  2.8900E+00  5.0008E+00  1.2195E-01  1.0608E+00  1.1234E-01
             8.9958E+00
 PARAMETER:  4.0677E-01  4.9672E-01  1.5123E+00  1.2379E-01  5.9802E-01  1.1613E+00  1.7096E+00 -2.0042E+00  1.5906E-01 -2.0862E+00
             2.2968E+00
 GRADIENT:  -6.5685E+00  2.6873E+00  1.3090E+00  4.0013E-01  8.4569E-01 -1.7094E+00 -1.0570E+00  3.1428E-03  2.6005E-01  1.5964E-01
             1.4729E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -835.133514357947        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      904
 NPARAMETR:  1.3963E+00  1.4205E+00  2.2682E+00  1.0074E+00  1.4357E+00  2.9040E+00  5.0480E+00  3.1805E-02  9.2605E-01  7.8398E-02
             9.0276E+00
 PARAMETER:  4.3382E-01  4.5098E-01  9.1899E-01  1.0740E-01  4.6167E-01  1.1661E+00  1.7190E+00 -3.3481E+00  2.3178E-02 -2.4460E+00
             2.3003E+00
 GRADIENT:   2.8366E-01  5.1026E-01  2.7341E-01  4.1221E-01 -7.6267E-01  5.8723E-01  2.6491E-02  8.4119E-04 -7.4638E-01  8.9602E-02
             1.4896E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -835.228324569771        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1085             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3977E+00  1.3971E+00  2.3279E+00  1.0154E+00  1.4457E+00  2.9835E+00  5.0037E+00  1.7351E-02  9.8593E-01  3.0494E-02
             9.0082E+00
 PARAMETER:  4.3485E-01  4.3441E-01  9.4497E-01  1.1527E-01  4.6857E-01  1.1931E+00  1.7102E+00 -3.9541E+00  8.5830E-02 -3.3902E+00
             2.2981E+00
 GRADIENT:   3.6475E+01  6.6153E+00  1.8409E-01  2.6432E+00  1.2808E+00  9.2357E+01  6.5386E+01  2.5684E-04  1.3216E-01  1.4509E-02
             2.9045E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -835.262304029688        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1223             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3971E+00  1.4178E+00  2.3201E+00  1.0094E+00  1.4449E+00  2.9517E+00  5.0021E+00  1.0000E-02  9.8345E-01  1.0324E-02
             8.9969E+00
 PARAMETER:  4.3441E-01  4.4911E-01  9.4159E-01  1.0932E-01  4.6804E-01  1.1824E+00  1.7099E+00 -4.6605E+00  8.3308E-02 -4.4733E+00
             2.2969E+00
 GRADIENT:   3.6522E+01  7.6528E+00  3.4506E-01  2.0773E+00  7.7451E-01  8.7640E+01  6.5856E+01  0.0000E+00  2.0567E-01  1.7038E-03
             2.7317E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -835.265042140074        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1407
 NPARAMETR:  1.3969E+00  1.4204E+00  2.3080E+00  1.0044E+00  1.4448E+00  2.9521E+00  4.9869E+00  1.0000E-02  9.7949E-01  1.0000E-02
             8.9989E+00
 PARAMETER:  4.3425E-01  4.5091E-01  9.3636E-01  1.0435E-01  4.6795E-01  1.1825E+00  1.7068E+00 -4.6605E+00  7.9279E-02 -4.6683E+00
             2.2971E+00
 GRADIENT:   6.2647E-01 -2.9306E-01  1.8311E-01 -5.6968E-01 -3.1097E-01  6.5414E+00 -1.8306E+00  0.0000E+00  1.6830E-01  0.0000E+00
             4.2752E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -835.268256872391        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1605             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3968E+00  1.4275E+00  2.2908E+00  1.0040E+00  1.4452E+00  2.9523E+00  4.9702E+00  1.0000E-02  9.6891E-01  1.0000E-02
             8.9963E+00
 PARAMETER:  4.3420E-01  4.5590E-01  9.2888E-01  1.0402E-01  4.6824E-01  1.1826E+00  1.7035E+00 -4.6605E+00  6.8419E-02 -4.6683E+00
             2.2968E+00
 GRADIENT:   3.6431E+01  7.7189E+00  2.8283E-01  2.5377E+00  1.0396E+00  8.7772E+01  6.4582E+01  0.0000E+00 -1.2729E-02  0.0000E+00
             2.6811E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -835.269193917308        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1796
 NPARAMETR:  1.3968E+00  1.4320E+00  2.2829E+00  1.0023E+00  1.4449E+00  2.9522E+00  4.9622E+00  1.0000E-02  9.6506E-01  1.0000E-02
             8.9959E+00
 PARAMETER:  4.3420E-01  4.5904E-01  9.2546E-01  1.0225E-01  4.6801E-01  1.1826E+00  1.7018E+00 -4.6605E+00  6.4431E-02 -4.6683E+00
             2.2968E+00
 GRADIENT:   5.9799E-01  2.1045E-02  1.1897E-01  8.8362E-01 -1.5601E-01  6.5992E+00 -2.6380E+00  0.0000E+00 -1.6841E-01  0.0000E+00
            -9.5522E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -835.270509496477        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1986
 NPARAMETR:  1.3968E+00  1.4344E+00  2.2750E+00  1.0005E+00  1.4445E+00  2.9522E+00  4.9559E+00  1.0000E-02  9.6378E-01  1.0000E-02
             8.9964E+00
 PARAMETER:  4.3421E-01  4.6072E-01  9.2199E-01  1.0054E-01  4.6776E-01  1.1826E+00  1.7006E+00 -4.6605E+00  6.3112E-02 -4.6683E+00
             2.2968E+00
 GRADIENT:   5.9894E-01 -2.5978E-02  1.0671E-01  7.5276E-01 -1.1847E-01  6.6006E+00 -2.6398E+00  0.0000E+00 -1.4375E-01  0.0000E+00
            -8.2661E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -835.272005413104        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2183             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3968E+00  1.4369E+00  2.2631E+00  9.9693E-01  1.4445E+00  2.9520E+00  4.9523E+00  1.0000E-02  9.6846E-01  1.0000E-02
             8.9989E+00
 PARAMETER:  4.3421E-01  4.6252E-01  9.1673E-01  9.6927E-02  4.6779E-01  1.1825E+00  1.6999E+00 -4.6605E+00  6.7951E-02 -4.6683E+00
             2.2971E+00
 GRADIENT:   3.6421E+01  7.6786E+00  1.9691E-01  1.4390E+00  1.3094E+00  8.7736E+01  6.4700E+01  0.0000E+00  2.2472E-01  0.0000E+00
             2.7747E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -835.272608894831        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2375
 NPARAMETR:  1.3969E+00  1.4394E+00  2.2566E+00  9.9562E-01  1.4442E+00  2.9520E+00  4.9477E+00  1.0000E-02  9.6706E-01  1.0000E-02
             8.9991E+00
 PARAMETER:  4.3422E-01  4.6423E-01  9.1385E-01  9.5612E-02  4.6753E-01  1.1825E+00  1.6989E+00 -4.6605E+00  6.6507E-02 -4.6683E+00
             2.2971E+00
 GRADIENT:   5.8796E-01 -2.5097E-01  1.0402E-02 -3.0082E-01  1.7693E-01  6.5506E+00 -2.0795E+00  0.0000E+00  1.2709E-01  0.0000E+00
             2.2549E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -835.273675746716        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     2578             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3969E+00  1.4435E+00  2.2504E+00  9.9538E-01  1.4433E+00  2.9521E+00  4.9394E+00  1.0000E-02  9.6051E-01  1.0000E-02
             8.9975E+00
 PARAMETER:  4.3424E-01  4.6705E-01  9.1109E-01  9.5365E-02  4.6691E-01  1.1825E+00  1.6973E+00 -4.6605E+00  5.9705E-02 -4.6683E+00
             2.2969E+00
 GRADIENT:   3.6435E+01  7.9635E+00  2.4509E-01  2.0647E+00  1.1150E+00  8.7771E+01  6.3965E+01  0.0000E+00  4.7636E-02  0.0000E+00
             2.7074E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -835.274083873224        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2769
 NPARAMETR:  1.3969E+00  1.4453E+00  2.2451E+00  9.9442E-01  1.4429E+00  2.9520E+00  4.9358E+00  1.0000E-02  9.5906E-01  1.0000E-02
             8.9975E+00
 PARAMETER:  4.3424E-01  4.6835E-01  9.0875E-01  9.4402E-02  4.6665E-01  1.1825E+00  1.6965E+00 -4.6605E+00  5.8198E-02 -4.6683E+00
             2.2969E+00
 GRADIENT:   5.9674E-01 -5.3792E-02  6.2628E-02  3.7474E-01 -3.0241E-02  6.5825E+00 -2.4911E+00  0.0000E+00 -5.9918E-02  0.0000E+00
            -5.0896E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -835.274630175364        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2966             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3969E+00  1.4469E+00  2.2378E+00  9.9255E-01  1.4426E+00  2.9519E+00  4.9331E+00  1.0000E-02  9.6039E-01  1.0000E-02
             8.9987E+00
 PARAMETER:  4.3424E-01  4.6945E-01  9.0547E-01  9.2518E-02  4.6646E-01  1.1825E+00  1.6960E+00 -4.6605E+00  5.9579E-02 -4.6683E+00
             2.2971E+00
 GRADIENT:   3.6429E+01  7.9344E+00  2.0625E-01  1.6212E+00  1.2169E+00  8.7755E+01  6.4066E+01  0.0000E+00  1.4855E-01  0.0000E+00
             2.7483E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -835.274913056939        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3156
 NPARAMETR:  1.3969E+00  1.4485E+00  2.2333E+00  9.9176E-01  1.4423E+00  2.9519E+00  4.9301E+00  1.0000E-02  9.5940E-01  1.0000E-02
             8.9987E+00
 PARAMETER:  4.3424E-01  4.7053E-01  9.0349E-01  9.1730E-02  4.6623E-01  1.1825E+00  1.6954E+00 -4.6605E+00  5.8555E-02 -4.6683E+00
             2.2971E+00
 GRADIENT:   5.9212E-01 -1.4782E-01  2.2192E-02 -6.5223E-02  7.6554E-02  6.5640E+00 -2.2688E+00  0.0000E+00  4.6894E-02  0.0000E+00
            -8.0478E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -835.275217585901        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3353             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3969E+00  1.4514E+00  2.2277E+00  9.9127E-01  1.4416E+00  2.9519E+00  4.9242E+00  1.0000E-02  9.5503E-01  1.0000E-02
             8.9978E+00
 PARAMETER:  4.3426E-01  4.7255E-01  9.0095E-01  9.1232E-02  4.6573E-01  1.1825E+00  1.6942E+00 -4.6605E+00  5.3988E-02 -4.6683E+00
             2.2970E+00
 GRADIENT:   3.6437E+01  8.1154E+00  2.3413E-01  1.9906E+00  1.0932E+00  8.7779E+01  6.3600E+01  0.0000E+00  4.1638E-02  0.0000E+00
             2.7070E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -835.275415196706        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3544
 NPARAMETR:  1.3969E+00  1.4525E+00  2.2241E+00  9.9062E-01  1.4413E+00  2.9519E+00  4.9221E+00  1.0000E-02  9.5434E-01  1.0000E-02
             8.9979E+00
 PARAMETER:  4.3426E-01  4.7332E-01  8.9935E-01  9.0579E-02  4.6554E-01  1.1825E+00  1.6937E+00 -4.6605E+00  5.3261E-02 -4.6683E+00
             2.2970E+00
 GRADIENT:   5.9751E-01 -4.1747E-02  4.9047E-02  2.8484E-01 -4.0597E-02  6.5811E+00 -2.4902E+00  0.0000E+00 -5.3658E-02  0.0000E+00
            -4.7021E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -835.275675809463        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3741             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3969E+00  1.4534E+00  2.2192E+00  9.8948E-01  1.4411E+00  2.9519E+00  4.9207E+00  1.0000E-02  9.5577E-01  1.0000E-02
             8.9988E+00
 PARAMETER:  4.3426E-01  4.7389E-01  8.9712E-01  8.9421E-02  4.6542E-01  1.1824E+00  1.6935E+00 -4.6605E+00  5.4759E-02 -4.6683E+00
             2.2971E+00
 GRADIENT:   3.6429E+01  8.0765E+00  1.9576E-01  1.6500E+00  1.1920E+00  8.7761E+01  6.3719E+01  0.0000E+00  1.2886E-01  0.0000E+00
             2.7420E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -835.275794001997        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3931
 NPARAMETR:  1.3969E+00  1.4543E+00  2.2164E+00  9.8910E-01  1.4409E+00  2.9518E+00  4.9190E+00  1.0000E-02  9.5490E-01  1.0000E-02
             8.9988E+00
 PARAMETER:  4.3426E-01  4.7451E-01  8.9591E-01  8.9038E-02  4.6525E-01  1.1824E+00  1.6931E+00 -4.6605E+00  5.3854E-02 -4.6683E+00
             2.2971E+00
 GRADIENT:   5.9248E-01 -1.0350E-01  1.4849E-02  1.2347E-02  4.5008E-02  6.5680E+00 -2.3444E+00  0.0000E+00  1.9159E-02  0.0000E+00
            -1.7333E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -835.275894014846        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     4126             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3969E+00  1.4561E+00  2.2128E+00  9.8868E-01  1.4402E+00  2.9518E+00  4.9153E+00  1.0000E-02  9.5212E-01  1.0000E-02
             8.9983E+00
 PARAMETER:  4.3428E-01  4.7576E-01  8.9427E-01  8.8612E-02  4.6479E-01  1.1824E+00  1.6924E+00 -4.6605E+00  5.0939E-02 -4.6683E+00
             2.2970E+00
 GRADIENT:   3.6443E+01  8.1923E+00  2.2652E-01  1.8773E+00  1.0696E+00  8.7778E+01  6.3429E+01  0.0000E+00  5.5559E-02  0.0000E+00
             2.7153E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -835.275998816032        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     4317
 NPARAMETR:  1.3969E+00  1.4567E+00  2.2103E+00  9.8828E-01  1.4401E+00  2.9518E+00  4.9142E+00  1.0000E-02  9.5210E-01  1.0000E-02
             8.9984E+00
 PARAMETER:  4.3427E-01  4.7615E-01  8.9315E-01  8.8210E-02  4.6473E-01  1.1824E+00  1.6921E+00 -4.6605E+00  5.0912E-02 -4.6683E+00
             2.2970E+00
 GRADIENT:   5.9758E-01 -5.1410E-02  3.1969E-02  1.6853E-01 -3.2782E-02  6.5766E+00 -2.4518E+00  0.0000E+00 -2.9348E-02  0.0000E+00
            -3.6659E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -835.276078280239        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4502
 NPARAMETR:  1.3969E+00  1.4572E+00  2.2051E+00  9.8728E-01  1.4401E+00  2.9517E+00  4.9134E+00  1.0000E-02  9.5323E-01  1.0000E-02
             8.9995E+00
 PARAMETER:  4.3427E-01  4.7651E-01  8.9221E-01  8.7824E-02  4.6462E-01  1.1824E+00  1.6919E+00 -4.6605E+00  5.1104E-02 -4.6683E+00
             2.2971E+00
 GRADIENT:   2.3480E-03  3.1179E-04  2.4825E-02  1.1115E-01 -1.3605E-02  4.0305E-03 -5.3894E-03  0.0000E+00 -6.2794E-03  0.0000E+00
            -7.5389E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4502
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.2085E-03  1.8444E-02 -1.6310E-05 -5.4055E-02  1.4812E-05
 SE:             2.9238E-02  2.5681E-02  1.9343E-05  1.0997E-02  7.6154E-05
 N:                     100         100         100         100         100

 P VAL.:         8.0526E-01  4.7264E-01  3.9912E-01  8.8637E-07  8.4578E-01

 ETASHRINKSD(%)  2.0490E+00  1.3964E+01  9.9935E+01  6.3160E+01  9.9745E+01
 ETASHRINKVR(%)  4.0560E+00  2.5978E+01  1.0000E+02  8.6428E+01  9.9999E+01
 EBVSHRINKSD(%)  1.7558E+00  8.8859E+00  9.9906E+01  7.0157E+01  9.9637E+01
 EBVSHRINKVR(%)  3.4807E+00  1.6982E+01  1.0000E+02  9.1094E+01  9.9999E+01
 RELATIVEINF(%)  9.6078E+01  3.9301E+01  1.6958E-05  4.1481E+00  2.6721E-04
 EPSSHRINKSD(%)  9.9337E+00
 EPSSHRINKVR(%)  1.8881E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -835.27607828023872     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       267.45016156536838     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   120.12
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    11.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -835.276       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.40E+00  1.46E+00  2.21E+00  9.88E-01  1.44E+00  2.95E+00  4.91E+00  1.00E-02  9.52E-01  1.00E-02  9.00E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.67E+01
 
 TH 2
+        3.75E-01  5.53E+00
 
 TH 3
+       -4.00E-01  2.01E+00  5.50E+00
 
 TH 4
+       -1.64E+01  2.09E+01 -1.28E+01  1.77E+02
 
 TH 5
+        2.59E+00 -1.25E+01 -2.86E+01  5.52E+01  1.50E+02
 
 TH 6
+        4.89E+00 -2.71E-01 -3.18E-01 -2.78E+00  1.82E+00  6.86E-01
 
 TH 7
+        3.26E+00 -2.64E+00 -3.13E-01 -1.46E+01  2.77E+00  5.68E-01  1.68E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.75E+00 -4.20E+00 -1.15E+00 -1.98E+01  7.65E+00  7.47E-01  2.45E+00  0.00E+00  3.69E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -9.26E-01 -1.02E+00  4.33E-01 -6.92E+00 -1.77E+00 -2.58E-02  5.24E-01  0.00E+00  8.26E-01  0.00E+00  6.43E-01
 
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
+        6.21E+01
 
 TH 2
+       -8.93E-01  1.82E+01
 
 TH 3
+        4.71E-01  1.19E+00  2.94E+00
 
 TH 4
+       -3.36E+00  2.51E+01 -4.21E+00  1.54E+02
 
 TH 5
+       -2.34E+00 -8.93E+00 -1.19E+01  9.31E+00  6.53E+01
 
 TH 6
+       -1.90E-01 -2.24E-01  6.76E-02  1.28E+00 -1.04E+00  2.00E+01
 
 TH 7
+        2.39E-02  2.03E+00 -5.83E-01 -1.42E+01  2.15E+00 -4.90E-01  5.09E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -3.44E-02 -1.72E+00 -1.37E+00 -1.98E+01  4.06E+00 -2.48E-01  2.05E+00  0.00E+00  1.14E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.05E+00 -2.06E+00  6.72E-03 -9.18E+00  6.95E-01  5.43E-01  5.88E-01  0.00E+00  2.43E+00  0.00E+00  7.70E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        6.32E+01
 
 TH 2
+        3.02E+01  1.82E+01
 
 TH 3
+        1.48E+00  9.49E-01  1.08E+00
 
 TH 4
+        3.84E+01  2.64E+01  1.01E+00  1.64E+02
 
 TH 5
+       -1.44E+01 -7.78E+00 -5.33E+00 -1.53E+01  3.19E+01
 
 TH 6
+        1.33E+01  6.29E+00 -2.96E-01 -1.60E+01 -3.16E-01  1.72E+01
 
 TH 7
+        3.24E+00  1.19E+00 -4.48E-01 -1.56E+01  2.36E+00  5.55E+00  3.74E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.61E+00 -3.08E+00 -8.78E-01 -3.63E+01  7.17E+00  4.81E+00  3.45E+00  0.00E+00  1.36E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.21E+01 -1.02E+01 -2.67E+00 -7.23E+01  1.93E+01  1.11E+01  6.43E+00  0.00E+00  2.22E+01  0.00E+00  1.72E+02
 
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
 #CPUT: Total CPU Time in Seconds,      131.905
Stop Time:
Thu Sep 30 09:35:08 CDT 2021
