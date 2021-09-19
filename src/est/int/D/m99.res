Sat Sep 18 07:51:49 CDT 2021
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
$DATA ../../../../data/int/D/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   55224.8425700822        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.6915E+02  8.2678E+02  2.9709E+01  8.4013E+02  2.7567E+01 -4.5623E+03 -2.3511E+03 -8.5026E+01 -2.9707E+03 -1.1311E+03
            -1.0614E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -461.225280342109        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2807E+00  1.9652E+00  8.4582E-01  1.8036E+00  1.0004E+00  3.6702E+00  3.4540E+00  9.9770E-01  1.4327E+00  1.3700E+00
             1.3125E+01
 PARAMETER:  3.4744E-01  7.7559E-01 -6.7448E-02  6.8980E-01  1.0037E-01  1.4003E+00  1.3395E+00  9.7700E-02  4.5954E-01  4.1479E-01
             2.6745E+00
 GRADIENT:  -8.7658E+00  5.2844E+01 -2.4759E+01  2.5789E+02 -7.5893E+00  3.4161E+01 -1.6136E+02  4.9152E+00 -6.1756E+01  3.2264E+01
            -1.1533E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -524.850683656980        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2638E+00  6.8493E+00  2.6616E+00  2.8965E+00  4.8277E+00  2.3156E+00  4.4848E+00  6.5574E-01  2.4699E+01  5.2696E-01
             1.2652E+01
 PARAMETER:  3.3412E-01  2.0241E+00  1.0789E+00  1.1635E+00  1.6744E+00  9.3966E-01  1.6007E+00 -3.2199E-01  3.3068E+00 -5.4063E-01
             2.6378E+00
 GRADIENT:   5.2102E+01  5.3396E+01 -5.4300E+00  1.6419E+01  4.8283E+01  3.1376E+01  3.4244E+01 -1.7886E+00 -1.1717E+01  1.0928E+00
            -6.8827E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -570.904939819840        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  7.6681E-01  3.3109E+00  2.1675E+01  2.3268E+00  3.0144E+00  1.4707E+00  4.4029E+00  4.0609E+00  6.7446E+00  1.1860E+00
             1.3340E+01
 PARAMETER: -1.6552E-01  1.2972E+00  3.1761E+00  9.4449E-01  1.2034E+00  4.8574E-01  1.5823E+00  1.5014E+00  2.0087E+00  2.7056E-01
             2.6908E+00
 GRADIENT:  -1.1990E+02  3.6887E+01 -1.8306E+01  3.9078E+01  8.0031E+00  2.4440E+01  2.6398E+01 -1.0409E+01  5.4643E+00  1.6237E+01
             3.3686E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -599.281744708963        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.1605E-01  2.2572E+00  1.8711E+01  1.9570E+00  3.0793E+00  1.5644E+00  3.1291E+00  1.0802E+01  4.7511E+00  1.0564E+00
             1.3603E+01
 PARAMETER:  1.2314E-02  9.1414E-01  3.0291E+00  7.7142E-01  1.2247E+00  5.4752E-01  1.2407E+00  2.4797E+00  1.6584E+00  1.5488E-01
             2.7103E+00
 GRADIENT:  -5.8766E+01  3.5612E+01 -1.2529E+01  3.7072E+01  1.0720E+01  2.3963E+01  1.6435E+01  4.8824E+00 -4.8434E+01  1.3005E+01
             4.7095E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -675.106323909993        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0036E+00  1.5742E+00  1.8328E+02  1.3504E+00  2.8570E+00  1.2122E+00  6.9275E-01  6.7241E+00  5.9662E+00  8.4794E-02
             1.4238E+01
 PARAMETER:  1.0362E-01  5.5378E-01  5.3110E+00  4.0039E-01  1.1498E+00  2.9246E-01 -2.6709E-01  2.0057E+00  1.8861E+00 -2.3675E+00
             2.7559E+00
 GRADIENT:   2.8163E+01 -3.6393E+01 -4.2189E+00  7.8916E+00  6.6798E+00  1.9928E+00  7.8084E+00  4.3581E+00  7.3326E+00  6.5876E-02
             9.1969E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -690.019738588605        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.5918E-01  1.6991E+00  1.8592E+03  9.0722E-01  2.8402E+00  1.2272E+00  2.7093E-01  3.6323E+00  6.5827E+00  3.2247E-02
             1.3354E+01
 PARAMETER:  5.8320E-02  6.3009E-01  7.6279E+00  2.6287E-03  1.1439E+00  3.0476E-01 -1.2059E+00  1.3899E+00  1.9844E+00 -3.3343E+00
             2.6918E+00
 GRADIENT:   2.8597E+00  1.4022E+01  2.1405E-01  2.7754E+00 -2.5955E+00  2.7570E+00 -1.5520E+00  5.2491E-01 -4.2450E-01  7.8708E-03
            -8.7896E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -690.028992852337        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  9.5731E-01  1.6828E+00  1.7827E+03  9.1085E-01  2.8322E+00  1.2250E+00  2.7620E-01  3.6406E+00  6.5263E+00  3.3434E-02
             1.3356E+01
 PARAMETER:  5.6371E-02  6.2046E-01  7.5859E+00  6.6263E-03  1.1410E+00  3.0293E-01 -1.1866E+00  1.3921E+00  1.9758E+00 -3.2982E+00
             2.6920E+00
 GRADIENT:   1.8295E+00  7.3864E+00  1.8708E-01  2.0621E+00 -1.9524E+00  2.1425E+00 -1.2549E+00  5.7391E-01 -5.4730E-01  8.8061E-03
            -3.2121E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -690.040236781053        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  9.5363E-01  1.6651E+00  1.4856E+03  9.1701E-01  2.8227E+00  1.2191E+00  2.9654E-01  3.5438E+00  6.4699E+00  3.9069E-02
             1.3321E+01
 PARAMETER:  5.2520E-02  6.0987E-01  7.4036E+00  1.3366E-02  1.1377E+00  2.9814E-01 -1.1156E+00  1.3652E+00  1.9672E+00 -3.1424E+00
             2.6894E+00
 GRADIENT:   5.9665E-01  4.1495E-01  8.0374E-02  1.5462E+00 -9.6326E-01  7.0726E-01 -9.5888E-01  7.2985E-01 -6.4217E-01  1.2550E-02
            -1.4123E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -690.100123684984        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  9.4713E-01  1.6474E+00  9.5238E+02  9.1586E-01  2.8055E+00  1.2079E+00  3.4361E-01  2.9016E+00  6.4289E+00  5.4697E-02
             1.3243E+01
 PARAMETER:  4.5676E-02  5.9920E-01  6.9590E+00  1.2106E-02  1.1316E+00  2.8891E-01 -9.6826E-01  1.1653E+00  1.9608E+00 -2.8059E+00
             2.6835E+00
 GRADIENT:  -1.2334E+00 -8.0970E+00 -1.9844E-01  8.9965E-01  6.7904E-01 -1.7605E+00 -4.8084E-01  8.5973E-01 -4.8421E-01  2.5772E-02
            -2.7850E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -691.210687203742        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  9.4944E-01  1.7054E+00  5.5857E+02  8.4421E-01  2.7924E+00  1.2139E+00  3.7131E-01  8.7483E-01  6.6414E+00  6.5760E-02
             1.3295E+01
 PARAMETER:  4.8120E-02  6.3380E-01  6.4254E+00 -6.9351E-02  1.1269E+00  2.9380E-01 -8.9072E-01 -3.3723E-02  1.9933E+00 -2.6217E+00
             2.6874E+00
 GRADIENT:  -1.8186E+00  1.1183E+00 -5.6645E-01  4.5889E-01  6.1250E-01  1.6386E-01 -1.3432E+00  1.8043E-01 -4.1671E-01  3.4280E-02
            -1.5942E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -692.080319776339        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      805
 NPARAMETR:  9.5334E-01  1.8157E+00  7.0637E+02  7.7525E-01  2.7745E+00  1.2185E+00  6.5141E-01  1.0000E-02  6.8653E+00  2.1764E-01
             1.3286E+01
 PARAMETER:  5.2216E-02  6.9645E-01  6.6601E+00 -1.5457E-01  1.1205E+00  2.9766E-01 -3.2861E-01 -5.7132E+00  2.0265E+00 -1.4249E+00
             2.6867E+00
 GRADIENT:   2.2198E-01 -2.2312E+00  1.6509E-01  1.5599E+00 -1.2617E+00 -7.5120E-01  1.3576E+00  0.0000E+00  6.2773E-01  4.0563E-01
             4.1662E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -692.389314672496        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  9.5451E-01  1.8476E+00  6.3955E+02  7.1247E-01  2.7712E+00  1.2237E+00  6.6433E-01  1.0000E-02  7.0025E+00  2.4213E-01
             1.3282E+01
 PARAMETER:  5.3441E-02  7.1388E-01  6.5608E+00 -2.3901E-01  1.1193E+00  3.0190E-01 -3.0898E-01 -6.9661E+00  2.0463E+00 -1.3183E+00
             2.6864E+00
 GRADIENT:   6.8945E-01 -2.4590E+00 -2.5712E-02 -1.1130E-02 -1.3308E+00  4.8442E-01  8.0094E-01  0.0000E+00 -1.3277E+01  4.9454E-01
            -1.1739E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -693.265309803795        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  9.5463E-01  1.9546E+00  6.4611E+02  5.7164E-01  2.7887E+00  1.2192E+00  6.8573E-01  1.0000E-02  7.8995E+00  2.4020E-01
             1.3379E+01
 PARAMETER:  5.3572E-02  7.7020E-01  6.5710E+00 -4.5925E-01  1.1256E+00  2.9817E-01 -2.7727E-01 -1.1251E+01  2.1668E+00 -1.3263E+00
             2.6937E+00
 GRADIENT:  -1.3535E+00 -2.5614E+00  4.2937E-02 -3.7055E-01  1.9502E-01  9.8906E-01 -8.9651E-01  0.0000E+00 -6.5441E-02  4.7455E-01
             1.0379E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -693.607289979290        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  9.5344E-01  2.0047E+00  6.3356E+02  5.3656E-01  2.7921E+00  1.2202E+00  7.7070E-01  1.0000E-02  8.0809E+00  5.5542E-02
             1.3289E+01
 PARAMETER:  5.2323E-02  7.9550E-01  6.5513E+00 -5.2258E-01  1.1268E+00  2.9899E-01 -1.6045E-01 -1.6658E+01  2.1895E+00 -2.7906E+00
             2.6869E+00
 GRADIENT:  -2.6897E-01 -8.6505E-02  2.7721E-03  4.2012E-02  7.8104E-02  1.0090E-01 -5.9825E-02  0.0000E+00  1.1791E-01  2.4353E-02
             5.5265E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -693.619840505650        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1442             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5352E-01  2.0082E+00  6.3241E+02  5.3298E-01  2.7912E+00  1.2200E+00  7.7608E-01  1.0000E-02  8.0985E+00  1.0000E-02
             1.3284E+01
 PARAMETER:  5.2407E-02  7.9726E-01  6.5495E+00 -5.2928E-01  1.1265E+00  2.9884E-01 -1.5350E-01 -2.1190E+01  2.1917E+00 -4.8309E+00
             2.6866E+00
 GRADIENT:   2.5089E-01  1.7123E+00  1.0351E-02  4.1397E-01  2.1695E-01  1.4632E-01  1.6036E-02  0.0000E+00  1.7715E+01  0.0000E+00
             4.8597E+00

0ITERATION NO.:   79    OBJECTIVE VALUE:  -693.619864679506        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:     1567
 NPARAMETR:  9.5364E-01  2.0087E+00  6.3198E+02  5.3239E-01  2.7913E+00  1.2200E+00  7.7611E-01  1.0000E-02  8.0983E+00  1.0000E-02
             1.3285E+01
 PARAMETER:  5.2531E-02  7.9747E-01  6.5489E+00 -5.3039E-01  1.1265E+00  2.9889E-01 -1.5347E-01 -2.1190E+01  2.1916E+00 -4.8309E+00
             2.6866E+00
 GRADIENT:   1.0785E-02  5.5711E-03 -6.1938E-04  7.7488E-03 -5.2500E-04  9.6069E-03 -1.5814E-02  0.0000E+00  1.2140E-01  0.0000E+00
             3.6560E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1567
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8622E-02 -7.0312E-02  1.6736E-06  4.3714E-02  4.3395E-05
 SE:             2.5729E-02  1.2252E-02  2.9004E-06  2.3448E-02  9.2241E-05
 N:                     100         100         100         100         100

 P VAL.:         2.6594E-01  9.5714E-09  5.6394E-01  6.2284E-02  6.3803E-01

 ETASHRINKSD(%)  1.3806E+01  5.8953E+01  9.9990E+01  2.1445E+01  9.9691E+01
 ETASHRINKVR(%)  2.5705E+01  8.3151E+01  1.0000E+02  3.8292E+01  9.9999E+01
 EBVSHRINKSD(%)  1.9073E+01  6.2931E+01  9.9966E+01  1.4551E+01  9.9644E+01
 EBVSHRINKVR(%)  3.4509E+01  8.6259E+01  1.0000E+02  2.6985E+01  9.9999E+01
 RELATIVEINF(%)  6.4856E+01  8.1950E+00  1.0885E-05  4.4943E+01  1.1605E-03
 EPSSHRINKSD(%)  2.6095E+00
 EPSSHRINKVR(%)  5.1509E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -693.61986467950567     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       960.46949508890509     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -693.620       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.54E-01  2.01E+00  6.32E+02  5.32E-01  2.79E+00  1.22E+00  7.76E-01  1.00E-02  8.10E+00  1.00E-02  1.33E+01
 


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
+        4.81E+02
 
 TH 2
+       -2.68E+01  9.61E+01
 
 TH 3
+       -3.01E-03 -3.19E-04  2.99E-06
 
 TH 4
+       -1.22E+02  2.27E+01  4.93E-04  7.14E+01
 
 TH 5
+       -4.32E+00 -1.17E+01 -1.49E-03 -5.16E+00  1.80E+01
 
 TH 6
+       -1.37E+01  1.36E+00 -1.60E-04  2.44E+00  5.99E-02  9.68E+01
 
 TH 7
+       -3.59E+01 -4.43E+01 -3.66E-03  8.76E+01  6.31E+00 -4.37E+01 -6.40E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.37E+00 -5.20E+00  8.81E-05  2.73E+00  4.95E-01 -5.20E-03  1.85E+00  0.00E+00  1.68E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.35E+01 -7.78E+00 -1.43E-05 -2.81E+00  8.98E-01  2.13E+00  4.89E+00  0.00E+00  2.92E-01  0.00E+00  5.43E+00
 
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
 #CPUT: Total CPU Time in Seconds,       63.728
Stop Time:
Sat Sep 18 07:52:54 CDT 2021
