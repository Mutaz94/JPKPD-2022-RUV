Wed Sep 29 10:14:39 CDT 2021
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
$DATA ../../../../data/int/D/dat98.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   61993.1720949864        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3232E+03  9.9038E+02  2.3278E+01  1.0018E+03 -5.5642E+01 -3.4412E+03 -2.2121E+03 -3.8243E+01 -3.0071E+03 -6.3479E+02
            -1.2048E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -312.329124662175        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0403E+00  2.0918E+00  8.7105E-01  1.5422E+00  9.2667E-01  4.5460E+00  4.3587E+00  9.8910E-01  1.6633E+00  1.1837E+00
             1.3293E+01
 PARAMETER:  1.3954E-01  8.3803E-01 -3.8058E-02  5.3320E-01  2.3841E-02  1.6143E+00  1.5722E+00  8.9040E-02  6.0882E-01  2.6867E-01
             2.6873E+00
 GRADIENT:  -1.4306E+01  6.1081E+01 -2.4645E+01  1.4792E+02 -3.7335E+01  1.8133E+02  2.2633E+00  4.0360E+00 -8.3296E+00  2.2820E+01
            -7.7059E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -383.146767943239        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  9.5074E-01  4.1593E+00  1.8264E+01  2.2863E+00  2.1139E+00  2.5685E+00  1.1030E+01  7.2475E-01  2.1882E+00  9.5473E-01
             1.3463E+01
 PARAMETER:  4.9485E-02  1.5253E+00  3.0050E+00  9.2694E-01  8.4852E-01  1.0433E+00  2.5006E+00 -2.2193E-01  8.8307E-01  5.3675E-02
             2.7000E+00
 GRADIENT:  -4.1394E+01  4.8598E+01 -4.7523E+00  9.7533E+01 -1.1747E+01  6.2716E+01  4.3409E+01  5.2774E-02  5.8431E+00  1.2824E+01
            -1.0295E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -465.995183823392        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1002E+00  1.8303E+00  6.9901E+00  9.1833E-01  2.1809E+00  1.9283E+00  4.1237E+00  6.6087E-01  1.0718E+00  3.9097E-01
             1.4095E+01
 PARAMETER:  1.9552E-01  7.0447E-01  2.0445E+00  1.4803E-02  8.7974E-01  7.5664E-01  1.5168E+00 -3.1421E-01  1.6938E-01 -8.3912E-01
             2.7458E+00
 GRADIENT:   1.5050E+01  1.2797E+00 -1.4528E+00 -6.3236E+00  1.9999E+00  1.3652E+01  3.1100E+00  5.3976E-02  7.4545E+00  2.1941E+00
             7.6612E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -472.215882539647        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      345
 NPARAMETR:  1.0159E+00  1.7365E+00  1.0800E+01  1.0219E+00  2.3414E+00  1.9285E+00  4.8720E+00  5.8403E-01  9.6561E-01  2.2540E-01
             1.3847E+01
 PARAMETER:  1.1575E-01  6.5189E-01  2.4796E+00  1.2164E-01  9.5074E-01  7.5673E-01  1.6835E+00 -4.3781E-01  6.5001E-02 -1.3899E+00
             2.7280E+00
 GRADIENT:  -1.7067E+01  4.0494E+00 -1.9858E+00 -7.3007E+00  8.8194E+00  6.6653E+00 -4.0498E+00  2.6120E-02  6.5421E+00  6.5920E-01
             5.5477E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -476.289716000708        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  1.0579E+00  1.2529E+00  1.7020E+01  1.1683E+00  2.3458E+00  1.8825E+00  5.6428E+00  7.3853E-01  7.1027E-01  1.2000E-01
             1.3914E+01
 PARAMETER:  1.5633E-01  3.2550E-01  2.9344E+00  2.5558E-01  9.5263E-01  7.3259E-01  1.8304E+00 -2.0309E-01 -2.4211E-01 -2.0202E+00
             2.7329E+00
 GRADIENT:  -3.0421E-01 -1.9758E+00 -7.8276E-01  1.5730E+00  2.3169E+00 -3.1408E-02 -1.6641E+00  2.4942E-02  4.2140E-02  1.5753E-01
            -6.9724E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -476.481842176821        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      698
 NPARAMETR:  1.0604E+00  1.3422E+00  2.2608E+01  1.1395E+00  2.3917E+00  1.8843E+00  5.5407E+00  6.2170E-01  6.4760E-01  8.9853E-02
             1.3958E+01
 PARAMETER:  1.5866E-01  3.9429E-01  3.2183E+00  2.3060E-01  9.7199E-01  7.3355E-01  1.8121E+00 -3.7530E-01 -3.3448E-01 -2.3096E+00
             2.7360E+00
 GRADIENT:   3.6141E-02 -4.9076E-01 -3.1803E-01  1.5939E+00  1.4961E+00 -1.4119E-01 -1.6332E+00  8.8490E-03 -2.1766E-01  8.6020E-02
            -4.2471E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -476.692626803893        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      876
 NPARAMETR:  1.0635E+00  1.3121E+00  6.4915E+01  1.1720E+00  2.4574E+00  1.8850E+00  5.6694E+00  5.0986E-01  7.2671E-01  2.0568E-02
             1.3998E+01
 PARAMETER:  1.6155E-01  3.7166E-01  4.2731E+00  2.5875E-01  9.9911E-01  7.3394E-01  1.8351E+00 -5.7362E-01 -2.1922E-01 -3.7840E+00
             2.7389E+00
 GRADIENT:   6.3704E-01  3.0617E-01 -1.4922E-02  3.6350E-01 -2.7956E-01 -5.3411E-01  6.5640E-01  6.7809E-04  1.9603E-02  4.4718E-03
             2.7207E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -476.717581221480        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  1.0625E+00  1.2969E+00  8.0926E+01  1.1762E+00  2.4663E+00  1.8905E+00  5.7242E+00  4.9486E-01  7.2124E-01  1.0000E-02
             1.4007E+01
 PARAMETER:  1.6059E-01  3.6000E-01  4.4935E+00  2.6230E-01  1.0027E+00  7.3682E-01  1.8447E+00 -6.0349E-01 -2.2678E-01 -5.1939E+00
             2.7396E+00
 GRADIENT:   1.3371E-01  1.4540E-01 -4.5167E-03 -2.2877E-01 -3.0558E-01  4.0951E-01  1.9287E+00  4.1237E-04 -1.0394E-01  0.0000E+00
             1.3166E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -476.719277455564        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1230
 NPARAMETR:  1.0625E+00  1.2895E+00  1.1091E+02  1.1827E+00  2.4789E+00  1.8896E+00  5.7135E+00  4.7583E-01  7.4286E-01  1.0000E-02
             1.4005E+01
 PARAMETER:  1.6065E-01  3.5428E-01  4.8087E+00  2.6781E-01  1.0078E+00  7.3637E-01  1.8428E+00 -6.4269E-01 -1.9725E-01 -5.2075E+00
             2.7394E+00
 GRADIENT:   1.1007E-01  2.8878E-02 -3.7061E-03  9.5616E-03  3.8798E-02  1.9612E-01  9.9666E-01  2.0404E-04 -6.5887E-03  0.0000E+00
             1.3033E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -476.739106281900        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1415
 NPARAMETR:  1.0617E+00  1.2639E+00  1.2535E+02  1.1872E+00  2.4823E+00  1.8905E+00  5.7908E+00  4.4576E-01  7.4164E-01  1.0000E-02
             1.3978E+01
 PARAMETER:  1.5988E-01  3.3422E-01  4.9311E+00  2.7164E-01  1.0092E+00  7.3683E-01  1.8563E+00 -7.0797E-01 -1.9889E-01 -5.2075E+00
             2.7375E+00
 GRADIENT:   3.0663E-01 -2.5330E-01 -4.5305E-03 -1.1819E+00  2.7009E-01  2.6125E-01  2.9017E+00  1.4249E-04 -7.0564E-02  0.0000E+00
            -7.4378E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -476.745327094216        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1599             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0616E+00  1.2603E+00  1.5346E+02  1.1957E+00  2.4833E+00  1.8906E+00  5.8063E+00  4.3280E-01  7.5183E-01  1.0000E-02
             1.3978E+01
 PARAMETER:  1.5982E-01  3.3136E-01  5.1334E+00  2.7875E-01  1.0096E+00  7.3692E-01  1.8589E+00 -7.3748E-01 -1.8524E-01 -5.2075E+00
             2.7375E+00
 GRADIENT:   5.1533E+00  2.6842E+00  2.8003E-03  4.6576E+00  1.7335E+00  1.0070E+01  6.2032E+01  9.1138E-05 -1.8640E-01  0.0000E+00
             4.7358E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -476.748259778449        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1775
 NPARAMETR:  1.0617E+00  1.2602E+00  1.6610E+02  1.1970E+00  2.4855E+00  1.8912E+00  5.8150E+00  4.1693E-01  7.6182E-01  1.0000E-02
             1.3993E+01
 PARAMETER:  1.5990E-01  3.3128E-01  5.2126E+00  2.7983E-01  1.0105E+00  7.3719E-01  1.8604E+00 -7.7483E-01 -1.7204E-01 -5.2075E+00
             2.7386E+00
 GRADIENT:  -9.1898E-04  1.4408E-01 -3.9712E-04 -5.1318E-01 -1.8053E-01  3.7016E-01  3.0769E+00  7.1264E-05 -3.2718E-02  0.0000E+00
             4.2826E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -476.753082022933        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1962
 NPARAMETR:  1.0614E+00  1.2507E+00  2.5042E+02  1.2025E+00  2.4877E+00  1.8909E+00  5.8415E+00  4.5426E-01  7.7480E-01  1.0000E-02
             1.3982E+01
 PARAMETER:  1.5958E-01  3.2369E-01  5.6231E+00  2.8444E-01  1.0113E+00  7.3705E-01  1.8650E+00 -6.8908E-01 -1.5515E-01 -5.2075E+00
             2.7378E+00
 GRADIENT:   4.7688E-02  2.1997E-01  2.6431E-03 -5.7235E-01 -6.1322E-01  2.8699E-01  3.4735E+00  3.6752E-05 -6.2372E-03  0.0000E+00
            -5.8661E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -476.755322966178        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2148
 NPARAMETR:  1.0615E+00  1.2446E+00  2.7575E+02  1.2056E+00  2.4894E+00  1.8910E+00  5.8473E+00  4.5272E-01  7.7808E-01  1.0000E-02
             1.3983E+01
 PARAMETER:  1.5967E-01  3.1880E-01  5.7195E+00  2.8699E-01  1.0120E+00  7.3712E-01  1.8660E+00 -6.9248E-01 -1.5092E-01 -5.2075E+00
             2.7379E+00
 GRADIENT:   5.0430E-02  1.4247E-01  2.2450E-03 -2.9420E-01 -5.8882E-01  3.2618E-01  3.2649E+00  3.0205E-05 -6.4460E-02  0.0000E+00
            -5.8680E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -476.757806509893        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2332             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0617E+00  1.2388E+00  3.3716E+02  1.2090E+00  2.4932E+00  1.8909E+00  5.8599E+00  4.4830E-01  7.8508E-01  1.0000E-02
             1.3981E+01
 PARAMETER:  1.5986E-01  3.1413E-01  5.9206E+00  2.8977E-01  1.0136E+00  7.3707E-01  1.8681E+00 -7.0229E-01 -1.4197E-01 -5.2075E+00
             2.7377E+00
 GRADIENT:   5.1097E+00  2.4783E+00  1.1705E-03  3.8968E+00  1.8561E+00  1.0021E+01  6.3778E+01  2.0371E-05  6.3746E-02  0.0000E+00
             4.8549E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -476.759524277121        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2517             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0618E+00  1.2346E+00  3.8073E+02  1.2121E+00  2.4965E+00  1.8910E+00  5.8699E+00  4.3950E-01  7.9225E-01  1.0000E-02
             1.3980E+01
 PARAMETER:  1.5999E-01  3.1072E-01  6.0421E+00  2.9236E-01  1.0149E+00  7.3709E-01  1.8698E+00 -7.2212E-01 -1.3288E-01 -5.2075E+00
             2.7377E+00
 GRADIENT:   5.1744E+00  2.4524E+00 -9.3873E-05  3.8544E+00  2.1372E+00  1.0001E+01  6.4060E+01  1.4385E-05  1.1118E-01  0.0000E+00
             4.8525E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -476.760299591873        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2706
 NPARAMETR:  1.0614E+00  1.2305E+00  4.0121E+02  1.2136E+00  2.4975E+00  1.8912E+00  5.8715E+00  4.4586E-01  7.9438E-01  1.0000E-02
             1.3983E+01
 PARAMETER:  1.5963E-01  3.0742E-01  6.0945E+00  2.9359E-01  1.0153E+00  7.3720E-01  1.8701E+00 -7.0775E-01 -1.3019E-01 -5.2075E+00
             2.7379E+00
 GRADIENT:  -6.7506E-03  1.7615E-02 -5.5723E-04 -5.1363E-01  1.1708E-02  3.0040E-01  3.3555E+00  1.3533E-05  1.8449E-02  0.0000E+00
            -2.9697E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -476.761422435119        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2897             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0615E+00  1.2267E+00  5.0318E+02  1.2168E+00  2.4986E+00  1.8912E+00  5.8798E+00  4.3304E-01  7.9657E-01  1.0000E-02
             1.3983E+01
 PARAMETER:  1.5966E-01  3.0434E-01  6.3210E+00  2.9620E-01  1.0157E+00  7.3723E-01  1.8715E+00 -7.3693E-01 -1.2743E-01 -5.2075E+00
             2.7378E+00
 GRADIENT:   4.9377E+00  2.3476E+00  1.0001E-04  4.4719E+00  2.0520E+00  1.0065E+01  6.4044E+01  1.2132E-05 -2.6314E-03  0.0000E+00
             4.8556E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -476.761838498120        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3080
 NPARAMETR:  1.0615E+00  1.2253E+00  6.0757E+02  1.2179E+00  2.4989E+00  1.8913E+00  5.8829E+00  4.2458E-01  7.9903E-01  1.0000E-02
             1.3982E+01
 PARAMETER:  1.5973E-01  3.0317E-01  6.5095E+00  2.9716E-01  1.0159E+00  7.3725E-01  1.8721E+00 -7.5666E-01 -1.2436E-01 -5.2075E+00
             2.7378E+00
 GRADIENT:   4.0911E-02  9.6478E-02  2.3041E-04  9.3914E-02 -2.5360E-01  3.3416E-01  3.2448E+00  7.6990E-06 -9.2206E-02  0.0000E+00
            -8.4881E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -476.762082423533        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3268
 NPARAMETR:  1.0616E+00  1.2238E+00  6.7963E+02  1.2191E+00  2.4999E+00  1.8913E+00  5.8861E+00  4.2300E-01  7.9989E-01  1.0000E-02
             1.3982E+01
 PARAMETER:  1.5981E-01  3.0195E-01  6.6215E+00  2.9811E-01  1.0163E+00  7.3725E-01  1.8726E+00 -7.6039E-01 -1.2328E-01 -5.2075E+00
             2.7378E+00
 GRADIENT:   6.0083E-02  1.0576E-01  1.6824E-04  2.3916E-01 -2.3347E-01  3.3469E-01  3.2165E+00  5.6169E-06 -1.2235E-01  0.0000E+00
            -9.2195E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -476.762593773609        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3452
 NPARAMETR:  1.0616E+00  1.2221E+00  7.3835E+02  1.2197E+00  2.5009E+00  1.8913E+00  5.8901E+00  4.3078E-01  8.0229E-01  1.0000E-02
             1.3983E+01
 PARAMETER:  1.5980E-01  3.0056E-01  6.7044E+00  2.9859E-01  1.0167E+00  7.3725E-01  1.8733E+00 -7.4216E-01 -1.2028E-01 -5.2075E+00
             2.7379E+00
 GRADIENT:   4.8769E-02  6.5562E-02 -6.2826E-05 -3.3204E-02 -1.0718E-01  3.1752E-01  3.3378E+00  4.0742E-06 -6.6482E-02  0.0000E+00
            -6.1221E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -476.763004655888        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3646             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0616E+00  1.2203E+00  9.9087E+02  1.2205E+00  2.5017E+00  1.8913E+00  5.8934E+00  4.2410E-01  8.0628E-01  1.0000E-02
             1.3983E+01
 PARAMETER:  1.5975E-01  2.9914E-01  6.9986E+00  2.9925E-01  1.0170E+00  7.3726E-01  1.8738E+00 -7.5780E-01 -1.1533E-01 -5.2075E+00
             2.7379E+00
 GRADIENT:   4.9742E+00  2.2597E+00 -1.1690E-05  4.2042E+00  2.1315E+00  1.0048E+01  6.4486E+01 -8.5312E-07  8.2528E-02  0.0000E+00
             4.8886E+01

0ITERATION NO.:  112    OBJECTIVE VALUE:  -476.763004655888        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     3703
 NPARAMETR:  1.0616E+00  1.2203E+00  9.9087E+02  1.2205E+00  2.5017E+00  1.8913E+00  5.8934E+00  4.2410E-01  8.0628E-01  1.0000E-02
             1.3983E+01
 PARAMETER:  1.5975E-01  2.9914E-01  6.9986E+00  2.9925E-01  1.0170E+00  7.3726E-01  1.8738E+00 -7.5780E-01 -1.1533E-01 -5.2075E+00
             2.7379E+00
 GRADIENT:   3.3703E-02  3.6909E-02 -4.6407E-05 -3.0591E-01 -7.9313E-02  3.0594E-01  3.4304E+00 -9.1834E-07  2.2734E-03  0.0000E+00
            -3.9375E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3703
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4274E-02  2.4454E-02 -8.4085E-07 -5.8440E-02  5.3956E-06
 SE:             2.7675E-02  2.4802E-02  4.6322E-06  1.1009E-02  8.9707E-05
 N:                     100         100         100         100         100

 P VAL.:         3.8044E-01  3.2414E-01  8.5596E-01  1.1096E-07  9.5204E-01

 ETASHRINKSD(%)  7.2843E+00  1.6910E+01  9.9984E+01  6.3117E+01  9.9699E+01
 ETASHRINKVR(%)  1.4038E+01  3.0960E+01  1.0000E+02  8.6396E+01  9.9999E+01
 EBVSHRINKSD(%)  9.6651E+00  1.1075E+01  9.9980E+01  6.9487E+01  9.9626E+01
 EBVSHRINKVR(%)  1.8396E+01  2.0923E+01  1.0000E+02  9.0690E+01  9.9999E+01
 RELATIVEINF(%)  8.1293E+01  4.1974E+01  9.0604E-07  4.8919E+00  3.1114E-04
 EPSSHRINKSD(%)  3.1290E+00
 EPSSHRINKVR(%)  6.1601E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -476.76300465588844     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       1177.3263551125224     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   143.74
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    18.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -476.763       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.22E+00  9.91E+02  1.22E+00  2.50E+00  1.89E+00  5.89E+00  4.24E-01  8.06E-01  1.00E-02  1.40E+01
 


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
+        3.23E+02
 
 TH 2
+       -1.54E+01  4.54E+00
 
 TH 3
+       -5.03E-05  8.08E-06  2.15E-11
 
 TH 4
+       -1.53E+02  3.06E+01  6.47E-05  2.23E+02
 
 TH 5
+        2.39E+01 -3.97E+00 -1.02E-05 -3.12E+01  4.81E+00
 
 TH 6
+       -8.00E+01  3.34E+00  2.81E-05  5.36E+01 -1.25E+01  7.15E+01
 
 TH 7
+        1.28E+01 -2.18E+00 -4.97E-06 -1.64E+01  2.38E+00 -4.96E+00  1.23E+00
 
 TH 8
+        5.98E-01 -1.83E-02 -8.94E-08 -2.35E-01  4.14E-02 -1.85E-01  2.10E-02  1.16E-03
 
 TH 9
+        4.37E+01 -8.22E+00 -1.80E-05 -6.08E+01  8.66E+00 -1.64E+01  4.52E+00  6.93E-02  1.67E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.11E+00 -9.87E-01 -1.15E-06 -5.47E+00  5.84E-01  1.20E+00  3.34E-01 -1.17E-02  1.43E+00  0.00E+00  6.61E-01
 
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
+        2.12E+02
 
 TH 2
+       -1.65E+00  1.77E+01
 
 TH 3
+        6.55E-06  1.53E-05  5.23E-10
 
 TH 4
+       -7.42E+00  2.41E+01  1.90E-05  1.52E+02
 
 TH 5
+       -1.75E+00 -3.24E+00 -7.72E-05 -1.50E+01  2.37E+01
 
 TH 6
+        5.90E+00 -6.89E-01  1.79E-06  2.93E+00 -2.33E+00  5.09E+01
 
 TH 7
+        6.65E-01  2.18E+00 -2.60E-06 -1.07E+01  6.49E-01 -5.97E-01  3.33E+00
 
 TH 8
+        3.70E-01 -6.76E-02 -4.17E-06  9.43E-02 -1.61E-02 -1.60E-02  1.65E-03  4.93E-02
 
 TH 9
+        8.60E-01 -2.48E+00 -2.45E-05 -3.92E+01  4.28E+00 -9.48E-01  2.29E+00  1.54E-01  1.64E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.91E+00 -2.02E+00 -1.54E-06 -9.25E+00  8.35E-01 -8.93E-02  3.44E-01  3.70E-04  2.71E+00  0.00E+00  4.28E+00
 
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
+        2.31E+02
 
 TH 2
+        3.82E+01  1.73E+01
 
 TH 3
+        3.26E-05  1.05E-05  1.85E-10
 
 TH 4
+        8.02E+01  2.60E+01  6.96E-05  1.62E+02
 
 TH 5
+       -1.42E+01 -2.50E+00 -5.02E-05 -2.48E+01  1.55E+01
 
 TH 6
+        6.29E+01  5.89E+00  1.63E-06 -6.69E+00 -1.93E+00  5.67E+01
 
 TH 7
+       -4.10E+00  2.34E+00 -5.73E-06 -1.21E+01  2.72E+00  9.74E-01  3.70E+00
 
 TH 8
+        6.85E-04  1.86E-04  1.70E-09  1.39E-03 -5.22E-04 -5.13E-05 -7.80E-05  2.64E-07
 
 TH 9
+       -7.04E+00 -2.18E+00 -1.53E-05 -3.91E+01  5.84E+00  2.80E+00  2.20E+00 -4.21E-04  2.08E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.94E+01 -3.15E+00 -2.04E-05 -4.24E+01  8.35E+00 -1.05E+01  4.89E+00  2.25E-04  1.43E+01  0.00E+00  7.72E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,      161.885
Stop Time:
Wed Sep 29 10:17:32 CDT 2021
