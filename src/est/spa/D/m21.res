Sat Sep 18 15:11:48 CDT 2021
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
$DATA ../../../../data/spa/D/dat21.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1211.94257952598        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.7047E+02 -1.2667E+02 -5.0492E+01 -2.7956E+02  2.2298E+02 -1.0079E+03 -4.2027E+02 -1.8599E+01 -7.4925E+02 -4.1475E+02
            -3.4084E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1054.07441900566        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4840E+00  9.8086E-01  2.1511E+00  1.8933E+00  2.1097E+00  2.6059E+00  2.1372E+00  9.5200E-01  4.5777E+00  2.0877E+00
             4.2820E+00
 PARAMETER:  4.9477E-01  8.0679E-02  8.6597E-01  7.3831E-01  8.4653E-01  1.0578E+00  8.5950E-01  5.0813E-02  1.6212E+00  8.3605E-01
             1.5544E+00
 GRADIENT:   6.9035E+01 -2.0879E+01 -2.7951E+01  3.9600E+01 -2.5194E+01  5.7250E+01  1.6512E+01  9.9997E-01  1.0188E+02  1.4857E+01
             9.1452E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1088.92741039346        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2891E+00  1.6175E+00  1.7815E+00  1.0965E+00  6.2894E+00  2.7031E+00  2.0516E+00  7.7077E-01  4.6143E+00  5.1361E+00
             4.3468E+00
 PARAMETER:  3.5396E-01  5.8085E-01  6.7747E-01  1.9211E-01  1.9389E+00  1.0944E+00  8.1864E-01 -1.6037E-01  1.6292E+00  1.7363E+00
             1.5694E+00
 GRADIENT:   2.3376E+01 -3.2245E+00 -1.6380E+01  2.2957E+01 -6.1061E+00  8.0329E+01  2.9738E+01  8.6369E-01  6.4333E+01  2.0084E-01
             1.1744E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1132.80068268590        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1940E+00  1.3700E+00  2.7976E+00  8.5552E-01  5.7116E+00  2.1120E+00  5.8593E-01  1.1865E-01  3.7924E+00  4.8924E+00
             3.6328E+00
 PARAMETER:  2.7730E-01  4.1482E-01  1.1288E+00 -5.6051E-02  1.8425E+00  8.4763E-01 -4.3456E-01 -2.0316E+00  1.4330E+00  1.6877E+00
             1.3900E+00
 GRADIENT:   1.3274E+00 -1.5043E+01  9.1983E+00 -1.3258E+01 -4.9050E+00 -2.2991E+00  5.1934E+00 -1.7140E-02  9.1419E+00  5.8827E+00
             6.0505E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1145.78243819520        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.1438E+00  8.8525E-01  1.4877E+00  1.3347E+00  1.4278E+01  2.0572E+00  1.6190E-01  1.1616E+00  2.6197E+00  4.6943E+00
             3.6237E+00
 PARAMETER:  2.3438E-01 -2.1890E-02  4.9724E-01  3.8870E-01  2.7588E+00  8.2136E-01 -1.7208E+00  2.4977E-01  1.0631E+00  1.6463E+00
             1.3875E+00
 GRADIENT:  -4.6450E+00  1.8613E+01  5.3841E+00  1.8027E+01 -1.3877E+00  3.7730E+00  2.4927E-01 -7.4912E-01 -3.5770E+00 -5.3106E-02
             4.2312E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1206.84262694243        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  8.8454E-01  8.3294E-02  1.4530E-01  9.5287E-01  3.0586E+03  1.2819E+00  1.0000E-02  2.2087E+00  1.0591E+00  7.0628E+00
             3.5031E+00
 PARAMETER: -2.2687E-02 -2.3854E+00 -1.8289E+00  5.1720E-02  8.1257E+00  3.4835E-01 -1.0102E+01  8.9241E-01  1.5745E-01  2.0548E+00
             1.3536E+00
 GRADIENT:   1.0689E+02  4.9897E+01 -7.7272E+01  7.4495E+01 -1.0535E-02 -5.9343E+01  0.0000E+00  1.4749E+01  1.6502E+01 -1.3649E-06
            -4.4393E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1244.09402456745        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  5.9196E-01  2.1890E-02  5.6655E-02  5.0235E-01  6.6085E+04  1.2521E+00  1.0000E-02  1.4025E+00  6.0563E-01  9.8588E+00
             3.4305E+00
 PARAMETER: -4.2432E-01 -3.7217E+00 -2.7708E+00 -5.8846E-01  1.1199E+01  3.2483E-01 -1.5055E+01  4.3828E-01 -4.0148E-01  2.3884E+00
             1.3327E+00
 GRADIENT:   1.0656E+01  2.3452E+00 -9.5293E+00  5.3235E+01 -1.0301E-05 -4.0552E+01  0.0000E+00 -7.1901E+00  1.8873E-01  7.3304E-09
            -3.2754E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1247.87609983826        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      601
 NPARAMETR:  6.1470E-01  2.3886E-02  6.4529E-02  5.3073E-01  5.2184E+04  1.3953E+00  1.0000E-02  1.5090E+00  6.1816E-01  9.7403E+00
             3.5458E+00
 PARAMETER: -3.8663E-01 -3.6345E+00 -2.6406E+00 -5.3351E-01  1.0963E+01  4.3313E-01 -1.4648E+01  5.1143E-01 -3.8101E-01  2.3763E+00
             1.3658E+00
 GRADIENT:  -1.5846E+00  2.8859E+00  1.9698E+00 -1.2538E+01 -4.5887E-05 -6.9227E-02  0.0000E+00  1.0149E+00  5.5189E+00  1.4161E-08
            -2.7229E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1251.05343865730        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  6.9827E-01  2.2393E-02  1.0517E-01  6.9320E-01  2.0298E+04  1.4147E+00  1.0000E-02  2.0130E+00  3.4486E-01  8.3533E+00
             3.5477E+00
 PARAMETER: -2.5915E-01 -3.6990E+00 -2.1522E+00 -2.6643E-01  1.0018E+01  4.4695E-01 -1.3340E+01  7.9961E-01 -9.6461E-01  2.2227E+00
             1.3663E+00
 GRADIENT:  -3.8855E+00 -5.8647E-01  7.7843E+00 -9.6352E+00  7.2407E-05 -1.4527E+00  0.0000E+00 -1.3492E+00  1.4166E+00  2.1994E-08
            -2.2360E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1251.31497819238        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      953
 NPARAMETR:  6.9229E-01  1.6643E-02  9.9900E-02  6.7891E-01  2.5263E+04  1.4237E+00  1.0000E-02  2.0402E+00  2.0925E-01  8.8794E+00
             3.5642E+00
 PARAMETER: -2.6775E-01 -3.9957E+00 -2.2036E+00 -2.8727E-01  1.0237E+01  4.5323E-01 -1.3625E+01  8.1304E-01 -1.4642E+00  2.2837E+00
             1.3709E+00
 GRADIENT:   2.8015E-01 -2.8914E-01  1.2755E-01  3.2573E-01  8.3370E-05  3.5603E-04  0.0000E+00 -2.2177E-01  1.5737E-01  4.9781E-11
             2.4786E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1251.38034037575        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  7.0171E-01  2.2672E-02  1.0343E-01  6.9126E-01  1.5951E+04  1.4233E+00  1.0000E-02  2.0392E+00  2.2099E-01  9.4133E+00
             3.5610E+00
 PARAMETER: -2.5423E-01 -3.6866E+00 -2.1689E+00 -2.6924E-01  9.7773E+00  4.5298E-01 -1.2647E+01  8.1256E-01 -1.4097E+00  2.3421E+00
             1.3700E+00
 GRADIENT:   2.7730E+00 -2.2472E-01  1.0564E-01  4.8610E+00  1.7876E-04 -1.1126E-01  0.0000E+00 -5.1771E+00 -6.0271E-02  5.3637E-08
            -1.6284E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1251.49245485204        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  6.9092E-01  2.2274E-02  9.8835E-02  6.7334E-01  1.3225E+04  1.4248E+00  1.0000E-02  2.0665E+00  1.3273E-01  1.0894E+01
             3.5639E+00
 PARAMETER: -2.6973E-01 -3.7043E+00 -2.2143E+00 -2.9550E-01  9.5899E+00  4.5406E-01 -1.2004E+01  8.2586E-01 -1.9194E+00  2.4882E+00
             1.3709E+00
 GRADIENT:   6.7131E-02 -2.0288E-02  2.0815E-02 -1.9145E-01  1.4421E-04  1.0722E-02  0.0000E+00  1.5610E-01  6.7905E-03  1.7572E-07
             7.7981E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1251.49259836771        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1492
 NPARAMETR:  6.9083E-01  2.2281E-02  9.8793E-02  6.7320E-01  1.2073E+04  1.4249E+00  1.0000E-02  2.0666E+00  1.2743E-01  1.0989E+01
             3.5641E+00
 PARAMETER: -2.6986E-01 -3.7040E+00 -2.2147E+00 -2.9571E-01  9.4987E+00  4.5411E-01 -1.1939E+01  8.2591E-01 -1.9602E+00  2.4969E+00
             1.3709E+00
 GRADIENT:   1.7456E-02  4.5441E-04 -4.3438E-02  5.5625E-02  1.5864E-04  6.2049E-03  0.0000E+00  1.9314E-02  2.3208E-04  2.1618E-07
             6.5762E-05

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1251.51250834964        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1678
 NPARAMETR:  6.9019E-01  2.2205E-02  9.8992E-02  6.7305E-01  6.2596E+01  1.4239E+00  1.0000E-02  2.0628E+00  1.2395E-01  1.0970E+01
             3.5642E+00
 PARAMETER: -2.7079E-01 -3.7074E+00 -2.2127E+00 -2.9594E-01  4.2367E+00  4.5341E-01 -1.1939E+01  8.2406E-01 -1.9878E+00  2.4952E+00
             1.3709E+00
 GRADIENT:  -1.6416E+00  2.2176E-01  1.6443E+00 -2.5099E+00  1.3514E-02 -3.4232E-01  0.0000E+00 -2.0444E-01 -2.7344E-04  8.7452E-03
             9.4220E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1251.58695169266        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1857
 NPARAMETR:  6.9390E-01  2.0351E-02  9.8357E-02  6.7504E-01  1.3020E+01  1.4277E+00  1.0000E-02  2.0547E+00  1.1251E-01  3.4752E+00
             3.5635E+00
 PARAMETER: -2.6543E-01 -3.7946E+00 -2.2192E+00 -2.9298E-01  2.6665E+00  4.5606E-01 -1.1939E+01  8.2015E-01 -2.0847E+00  1.3457E+00
             1.3707E+00
 GRADIENT:   6.6299E-01  2.3560E-01 -1.4425E+00  1.6727E+00 -7.4693E-02  5.0998E-02  0.0000E+00 -5.1735E-02 -3.8659E-03  4.3099E-02
            -3.6378E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1251.62012403144        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2033
 NPARAMETR:  6.9427E-01  2.0129E-02  9.9360E-02  6.7708E-01  1.1879E+01  1.4269E+00  1.0000E-02  2.0584E+00  8.1957E-02  7.2959E-01
             3.5653E+00
 PARAMETER: -2.6489E-01 -3.8056E+00 -2.2090E+00 -2.8996E-01  2.5748E+00  4.5551E-01 -1.1939E+01  8.2195E-01 -2.4016E+00 -2.1528E-01
             1.3713E+00
 GRADIENT:  -6.9993E-01  1.3876E-01  7.6184E-01 -4.9905E-01 -4.4933E-02 -3.2993E-01  0.0000E+00 -1.1345E+00 -1.2347E-02  7.5661E-03
            -5.2141E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1251.62732487615        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2208
 NPARAMETR:  6.9593E-01  2.0100E-02  9.9865E-02  6.7948E-01  1.1927E+01  1.4282E+00  1.0000E-02  2.0712E+00  8.4416E-02  1.8943E-01
             3.5647E+00
 PARAMETER: -2.6251E-01 -3.8070E+00 -2.2039E+00 -2.8643E-01  2.5788E+00  4.5644E-01 -1.1939E+01  8.2813E-01 -2.3720E+00 -1.5637E+00
             1.3711E+00
 GRADIENT:   1.1959E-01  1.3207E-02 -4.2391E-02  1.8363E-01  7.8309E-03 -2.3465E-02  0.0000E+00 -5.8393E-02 -6.0807E-03  5.2455E-04
            -2.0183E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1251.62836971798        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2383
 NPARAMETR:  6.9615E-01  2.0116E-02  9.9976E-02  6.8004E-01  1.1851E+01  1.4280E+00  1.0000E-02  2.0677E+00  1.0975E-01  1.2536E-01
             3.5644E+00
 PARAMETER: -2.6219E-01 -3.8062E+00 -2.2028E+00 -2.8560E-01  2.5724E+00  4.5628E-01 -1.1939E+01  8.2642E-01 -2.1095E+00 -1.9766E+00
             1.3710E+00
 GRADIENT:   2.3502E-02 -5.0228E-04  6.9161E-03  5.8958E-02  4.6437E-04 -1.9525E-03  0.0000E+00 -3.2967E-03 -2.0196E-05  2.3566E-04
             4.4744E-04

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1251.62855547156        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2561
 NPARAMETR:  6.9581E-01  2.0098E-02  9.9804E-02  6.7940E-01  1.1846E+01  1.4280E+00  1.0000E-02  2.0665E+00  1.0836E-01  3.3490E-02
             3.5643E+00
 PARAMETER: -2.6269E-01 -3.8072E+00 -2.2045E+00 -2.8654E-01  2.5720E+00  4.5630E-01 -1.1939E+01  8.2585E-01 -2.1223E+00 -3.2965E+00
             1.3710E+00
 GRADIENT:  -1.4569E-02  1.9316E-03  1.4157E-02 -2.9565E-02 -1.1504E-03  1.9610E-03  0.0000E+00  8.7247E-03  6.8562E-05  1.6944E-05
             3.3190E-03

0ITERATION NO.:   94    OBJECTIVE VALUE:  -1251.62856426697        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     2688
 NPARAMETR:  6.9582E-01  2.0097E-02  9.9802E-02  6.7941E-01  1.1851E+01  1.4280E+00  1.0000E-02  2.0664E+00  1.0854E-01  1.0000E-02
             3.5643E+00
 PARAMETER: -2.6267E-01 -3.8072E+00 -2.2046E+00 -2.8653E-01  2.5724E+00  4.5629E-01 -1.1939E+01  8.2581E-01 -2.1206E+00 -4.5914E+00
             1.3710E+00
 GRADIENT:   4.7342E-04 -1.6108E-04 -5.3683E-04  6.8123E-04  5.4481E-05 -1.9608E-04  0.0000E+00 -5.0899E-05  2.8028E-06  0.0000E+00
             2.1048E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2688
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8113E-03 -2.1569E-05 -3.2710E-03 -6.6832E-03 -1.2154E-06
 SE:             2.9536E-02  7.3322E-06  2.7507E-02  3.3131E-03  2.4339E-06
 N:                     100         100         100         100         100

 P VAL.:         9.2417E-01  3.2641E-03  9.0534E-01  4.3675E-02  6.1751E-01

 ETASHRINKSD(%)  1.0506E+00  9.9975E+01  7.8484E+00  8.8901E+01  9.9992E+01
 ETASHRINKVR(%)  2.0902E+00  1.0000E+02  1.5081E+01  9.8768E+01  1.0000E+02
 EBVSHRINKSD(%)  1.2218E+00  9.9961E+01  6.1204E+00  8.9740E+01  9.9991E+01
 EBVSHRINKVR(%)  2.4286E+00  1.0000E+02  1.1866E+01  9.8947E+01  1.0000E+02
 RELATIVEINF(%)  1.0642E+01  9.2298E-07  4.9920E+00  2.7117E-02  3.7065E-08
 EPSSHRINKSD(%)  2.3624E+01
 EPSSHRINKVR(%)  4.1667E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1251.6285642669670     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -516.47773770322885     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.24
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1251.629       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.96E-01  2.01E-02  9.98E-02  6.79E-01  1.19E+01  1.43E+00  1.00E-02  2.07E+00  1.09E-01  1.00E-02  3.56E+00
 


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
+        1.10E+03
 
 TH 2
+       -3.29E+02  1.69E+04
 
 TH 3
+       -1.12E+02 -7.98E+02  2.29E+04
 
 TH 4
+       -5.25E+02  2.83E+01 -6.12E+03  2.23E+03
 
 TH 5
+        2.41E-01 -5.73E+00 -1.31E+00  5.77E-01  6.67E-03
 
 TH 6
+        5.70E+00  1.92E+01  3.51E+01 -1.62E+01  1.09E-02  9.09E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        4.15E+00 -5.28E+01 -6.46E+01 -6.77E+01 -3.87E-02 -8.05E-02  0.00E+00  3.42E+01
 
 TH 9
+        2.88E+00 -4.66E+01  3.87E+01 -3.38E+01 -1.33E-02  1.38E+00  0.00E+00  7.78E+00  3.81E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.02E+01 -4.00E+01  4.53E+01 -1.83E+01  4.29E-03  2.29E+00  0.00E+00  3.47E+00  1.45E+00  0.00E+00  3.15E+01
 
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
 #CPUT: Total CPU Time in Seconds,       40.863
Stop Time:
Sat Sep 18 15:12:30 CDT 2021
