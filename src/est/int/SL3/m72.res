Sat Sep 25 02:36:44 CDT 2021
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
$DATA ../../../../data/int/SL3/dat72.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      981
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      881
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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -712.939623906350        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.7566E+01 -5.5931E+00  9.4314E+01  9.5392E+01  2.8602E+01  1.1223E+01 -1.3035E+02 -2.1467E+02 -9.1225E+01 -3.5591E+01
            -5.7841E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2760.15723416892        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0691E+00  1.3474E+00  1.0014E+00  8.2858E-01  1.2290E+00  9.0775E-01  1.1790E+00  9.9695E-01  8.3034E-01  1.1374E+00
             2.5380E+00
 PARAMETER:  1.6686E-01  3.9819E-01  1.0139E-01 -8.8042E-02  3.0617E-01  3.2084E-03  2.6468E-01  9.6945E-02 -8.5919E-02  2.2871E-01
             1.0314E+00
 GRADIENT:   2.3867E+01  1.5226E+01 -8.2230E+00  3.8665E-01  6.1847E-01 -2.2943E+01  3.2535E+01 -2.5410E+00 -1.3532E+01 -2.4581E+01
            -6.4091E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2767.84555009193        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0897E+00  1.5888E+00  9.3382E-01  7.1724E-01  1.3780E+00  9.6984E-01  9.3017E-01  1.3573E+00  7.3255E-01  1.5308E+00
             2.5360E+00
 PARAMETER:  1.8590E-01  5.6296E-01  3.1523E-02 -2.3234E-01  4.2060E-01  6.9377E-02  2.7617E-02  4.0547E-01 -2.1122E-01  5.2580E-01
             1.0306E+00
 GRADIENT:   6.5637E+01  7.3048E+01 -1.1290E+01  6.5453E+01 -1.1891E+01  2.3181E+00  5.9983E+00 -3.3568E+00 -1.4428E+01  6.6156E+00
            -5.6201E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2781.55629105045        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.0604E+00  1.7804E+00  1.5437E+00  5.9202E-01  1.7378E+00  9.6174E-01  7.8481E-01  3.5557E+00  1.0295E+00  1.6260E+00
             2.5463E+00
 PARAMETER:  1.5866E-01  6.7684E-01  5.3417E-01 -4.2422E-01  6.5264E-01  6.0993E-02 -1.4231E-01  1.3686E+00  1.2904E-01  5.8609E-01
             1.0346E+00
 GRADIENT:   1.5109E+00  5.0594E+01 -1.0498E+01  3.7197E+01  7.2107E+00  5.2320E-01  1.2414E+00  3.2109E+00  3.0801E+00  3.4694E+00
            -4.2944E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2785.37489615926        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0553E+00  1.7017E+00  2.6685E+00  6.0857E-01  1.7918E+00  9.6456E-01  8.0996E-01  4.1476E+00  8.8333E-01  1.6408E+00
             2.5450E+00
 PARAMETER:  1.5383E-01  6.3162E-01  1.0815E+00 -3.9665E-01  6.8323E-01  6.3919E-02 -1.1077E-01  1.5225E+00 -2.4058E-02  5.9517E-01
             1.0341E+00
 GRADIENT:  -7.1892E+00 -1.6172E+01  2.2091E-01 -1.7369E+01 -6.6369E+00  1.9469E+00 -1.9689E+00 -7.0416E-01  5.1023E-01  2.4275E+00
             2.0035E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2786.54535127113        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.0605E+00  1.4550E+00  3.7202E+00  8.0248E-01  1.7213E+00  9.5774E-01  9.7174E-01  4.2115E+00  6.5067E-01  1.4981E+00
             2.5378E+00
 PARAMETER:  1.5872E-01  4.7498E-01  1.4138E+00 -1.2005E-01  6.4306E-01  5.6823E-02  7.1334E-02  1.5378E+00 -3.2975E-01  5.0420E-01
             1.0313E+00
 GRADIENT:   3.3581E+00  2.6805E+01  3.4705E-01  1.2746E+01 -5.0095E-02 -1.0421E+00  4.5132E-01  7.9057E-02 -1.9223E+00 -4.7240E-01
            -2.6059E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2786.59162467096        NO. OF FUNC. EVALS.:  86
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  1.0601E+00  1.4343E+00  3.7588E+00  8.1543E-01  1.7135E+00  9.5783E-01  9.7924E-01  4.1938E+00  6.5170E-01  1.4861E+00
             2.5371E+00
 PARAMETER:  1.5835E-01  4.6065E-01  1.4241E+00 -1.0404E-01  6.3853E-01  5.6916E-02  7.9021E-02  1.5336E+00 -3.2817E-01  4.9619E-01
             1.0310E+00
 GRADIENT:   2.5776E+00  2.4188E+01  1.5234E-01  1.1804E+01  5.9221E-01 -9.9955E-01  3.5384E-01  1.9832E-01 -1.8620E+00 -5.2325E-01
            -4.8978E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2786.59528657341        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.0600E+00  1.4344E+00  3.7601E+00  8.1541E-01  1.7138E+00  9.6208E-01  9.7807E-01  4.1921E+00  6.5165E-01  1.4863E+00
             2.5377E+00
 PARAMETER:  1.5831E-01  4.6076E-01  1.4245E+00 -1.0407E-01  6.3870E-01  6.1344E-02  7.7826E-02  1.5332E+00 -3.2825E-01  4.9631E-01
             1.0313E+00
 GRADIENT:   2.5623E+00  2.4197E+01  1.7814E-01  1.1916E+01  6.1879E-01  6.7579E-01  1.2727E-01  1.4792E-01 -1.9233E+00 -5.1119E-01
             8.5556E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2786.59529737945        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.0600E+00  1.4344E+00  3.7601E+00  8.1540E-01  1.7138E+00  9.6078E-01  9.7757E-01  4.1921E+00  6.5165E-01  1.4863E+00
             2.5378E+00
 PARAMETER:  1.5831E-01  4.6072E-01  1.4245E+00 -1.0407E-01  6.3870E-01  5.9990E-02  7.7318E-02  1.5332E+00 -3.2825E-01  4.9632E-01
             1.0313E+00
 GRADIENT:   2.5190E+00  2.4126E+01  1.7792E-01  1.1912E+01  6.2743E-01  1.6628E-01  2.9808E-02  1.4674E-01 -1.9478E+00 -5.1247E-01
             8.2496E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2786.59547343719        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.0600E+00  1.4342E+00  3.7602E+00  8.1540E-01  1.7138E+00  9.5905E-01  9.7693E-01  4.1919E+00  6.5165E-01  1.4864E+00
             2.5378E+00
 PARAMETER:  1.5830E-01  4.6058E-01  1.4245E+00 -1.0408E-01  6.3870E-01  5.8192E-02  7.6665E-02  1.5331E+00 -3.2825E-01  4.9634E-01
             1.0313E+00
 GRADIENT:   2.4592E+00  2.3895E+01  1.7787E-01  1.1796E+01  6.4943E-01 -5.1105E-01 -9.7537E-02  1.4470E-01 -1.9772E+00 -5.1132E-01
             1.1269E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2786.59601767034        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  1.0600E+00  1.4335E+00  3.7604E+00  8.1537E-01  1.7138E+00  9.5591E-01  9.7579E-01  4.1912E+00  6.5164E-01  1.4865E+00
             2.5380E+00
 PARAMETER:  1.5828E-01  4.6010E-01  1.4245E+00 -1.0412E-01  6.3871E-01  5.4904E-02  7.5490E-02  1.5330E+00 -3.2826E-01  4.9641E-01
             1.0314E+00
 GRADIENT:   2.3340E+00  2.3064E+01  1.7843E-01  1.1264E+01  7.2040E-01 -1.7522E+00 -3.3278E-01  1.3945E-01 -2.0234E+00 -5.0084E-01
             2.6301E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2786.59779961774        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      870
 NPARAMETR:  1.0599E+00  1.4310E+00  3.7611E+00  8.1526E-01  1.7138E+00  9.5052E-01  9.7397E-01  4.1886E+00  6.5163E-01  1.4869E+00
             2.5386E+00
 PARAMETER:  1.5821E-01  4.5835E-01  1.4247E+00 -1.0425E-01  6.3871E-01  4.9250E-02  7.3620E-02  1.5324E+00 -3.2829E-01  4.9666E-01
             1.0316E+00
 GRADIENT:   2.0643E+00  2.0038E+01  1.8214E-01  9.0862E+00  9.6250E-01 -3.8899E+00 -7.3690E-01  1.2420E-01 -2.0708E+00 -4.5040E-01
             8.7631E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2786.60730956665        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      940
 NPARAMETR:  1.0598E+00  1.4249E+00  3.7627E+00  8.1501E-01  1.7138E+00  9.4641E-01  9.7315E-01  4.1825E+00  6.5159E-01  1.4877E+00
             2.5399E+00
 PARAMETER:  1.5803E-01  4.5413E-01  1.4251E+00 -1.0455E-01  6.3872E-01  4.4916E-02  7.2786E-02  1.5309E+00 -3.2834E-01  4.9726E-01
             1.0321E+00
 GRADIENT:   1.6908E+00  1.2726E+01  1.9394E-01  3.4475E+00  1.5225E+00 -5.5022E+00 -1.0476E+00  9.2503E-02 -2.0065E+00 -3.1037E-01
             2.4085E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2786.62681434669        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  1.0595E+00  1.4152E+00  3.7652E+00  8.1461E-01  1.7138E+00  9.5175E-01  9.7690E-01  4.1727E+00  6.5153E-01  1.4892E+00
             2.5416E+00
 PARAMETER:  1.5775E-01  4.4725E-01  1.4258E+00 -1.0505E-01  6.3870E-01  5.0547E-02  7.6626E-02  1.5286E+00 -3.2843E-01  4.9823E-01
             1.0328E+00
 GRADIENT:   1.4791E+00  1.0265E+00  2.1932E-01 -6.1346E+00  2.3795E+00 -3.2401E+00 -6.6474E-01  4.4253E-02 -1.6664E+00 -6.0260E-02
             4.8008E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2786.65012796956        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1106
 NPARAMETR:  1.0594E+00  1.4150E+00  3.7651E+00  8.1460E-01  1.7137E+00  9.6177E-01  9.8299E-01  4.1724E+00  6.5153E-01  1.4892E+00
             2.5409E+00
 PARAMETER:  1.5774E-01  4.4710E-01  1.4258E+00 -1.0506E-01  6.3865E-01  6.1019E-02  8.2841E-02  1.5285E+00 -3.2843E-01  4.9826E-01
             1.0325E+00
 GRADIENT:  -1.0595E+01 -6.7799E+00 -1.9332E-01 -7.9883E+00 -5.1732E-01 -1.5725E-01  3.5902E-01 -6.9976E-01 -1.5029E+00 -5.1754E-01
             2.6579E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2786.68992444409        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1268
 NPARAMETR:  1.0602E+00  1.4213E+00  3.7652E+00  8.1501E-01  1.7140E+00  9.6236E-01  9.7989E-01  4.1827E+00  6.5677E-01  1.4890E+00
             2.5392E+00
 PARAMETER:  1.5850E-01  4.5160E-01  1.4258E+00 -1.0456E-01  6.3884E-01  6.1633E-02  7.9685E-02  1.5310E+00 -3.2043E-01  4.9809E-01
             1.0319E+00
 GRADIENT:   3.4785E+00  8.5661E+00  1.5856E-01 -5.2201E-01  1.8941E+00  9.1503E-01  5.3010E-01  2.1008E-01 -1.3892E+00 -4.1340E-02
             2.5766E+00

0ITERATION NO.:   78    OBJECTIVE VALUE:  -2786.69054532093        NO. OF FUNC. EVALS.:  87
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  1.0603E+00  1.4213E+00  3.7653E+00  8.1502E-01  1.7140E+00  9.6235E-01  9.7985E-01  4.1827E+00  6.5698E-01  1.4890E+00
             2.5392E+00
 PARAMETER:  1.5851E-01  4.5158E-01  1.4258E+00 -1.0454E-01  6.3883E-01  6.1625E-02  7.9639E-02  1.5310E+00 -3.2011E-01  4.9809E-01
             1.0319E+00
 GRADIENT:   1.3635E+03 -4.8120E+02 -7.6497E+01  2.0798E+03 -1.7138E+02  1.4994E-02 -2.1751E+03  1.4286E+02  6.7828E+02 -2.1882E+02
            -2.1151E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1355
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.3626E-03 -8.7888E-03 -4.2212E-02  1.5787E-03 -3.3949E-02
 SE:             2.9464E-02  2.4637E-02  1.8027E-02  1.5721E-02  2.3307E-02
 N:                     100         100         100         100         100

 P VAL.:         8.5558E-01  7.2130E-01  1.9203E-02  9.2001E-01  1.4522E-01

 ETASHRINKSD(%)  1.2903E+00  1.7462E+01  3.9606E+01  4.7331E+01  2.1920E+01
 ETASHRINKVR(%)  2.5640E+00  3.1875E+01  6.3525E+01  7.2260E+01  3.9034E+01
 EBVSHRINKSD(%)  1.4660E+00  1.7823E+01  4.4159E+01  5.0662E+01  1.7791E+01
 EBVSHRINKVR(%)  2.9104E+00  3.2470E+01  6.8818E+01  7.5658E+01  3.2416E+01
 RELATIVEINF(%)  9.7010E+01  5.8766E+00  1.4524E+01  2.0265E+00  3.9949E+01
 EPSSHRINKSD(%)  1.7463E+01
 EPSSHRINKVR(%)  3.1877E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          881
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1619.1696955066332     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2786.6905453209329     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1167.5208498142997     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.18
 Elapsed covariance  time in seconds:    15.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2786.691       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.42E+00  3.77E+00  8.15E-01  1.71E+00  9.62E-01  9.80E-01  4.18E+00  6.57E-01  1.49E+00  2.54E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.35E-04  1.32E-03  1.04E-02  1.69E-04  2.19E-03  2.02E-06  1.96E-04  1.22E-02  4.18E-04  1.48E-03  8.68E-03
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.12E-07
 
 TH 2
+       -4.39E-07  1.74E-06
 
 TH 3
+       -3.46E-06  1.36E-05  1.07E-04
 
 TH 4
+        5.67E-08 -2.21E-07 -1.75E-06  2.87E-08
 
 TH 5
+       -7.29E-07  2.85E-06  2.24E-05 -3.68E-07  4.79E-06
 
 TH 6
+        4.26E-11 -1.76E-10 -9.45E-10  1.67E-11 -3.56E-10  4.10E-12
 
 TH 7
+       -6.57E-08  2.57E-07  2.03E-06 -3.32E-08  4.27E-07 -2.13E-11  3.85E-08
 
 TH 8
+        4.08E-06 -1.60E-05 -1.26E-04  2.06E-06 -2.66E-05  1.45E-09 -2.39E-06  1.49E-04
 
 TH 9
+        1.40E-07 -5.49E-07 -4.33E-06  7.08E-08 -9.10E-07  4.41E-11 -8.21E-08  5.10E-06  1.75E-07
 
 TH10
+       -4.96E-07  1.94E-06  1.53E-05 -2.51E-07  3.22E-06 -1.63E-10  2.91E-07 -1.80E-05 -6.20E-07  2.20E-06
 
 TH11
+       -2.42E-06  9.56E-06  7.49E-05 -1.22E-06  1.54E-05 -1.37E-10  1.42E-06 -8.61E-05 -3.03E-06  1.08E-05  7.54E-05
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.35E-04
 
 TH 2
+       -9.93E-01  1.32E-03
 
 TH 3
+       -9.97E-01  9.92E-01  1.04E-02
 
 TH 4
+        9.98E-01 -9.91E-01 -9.98E-01  1.69E-04
 
 TH 5
+       -9.93E-01  9.87E-01  9.87E-01 -9.93E-01  2.19E-03
 
 TH 6
+        6.28E-02 -6.60E-02 -4.50E-02  4.87E-02 -8.04E-02  2.02E-06
 
 TH 7
+       -9.99E-01  9.94E-01  9.98E-01 -1.00E+00  9.94E-01 -5.36E-02  1.96E-04
 
 TH 8
+        9.97E-01 -9.91E-01 -9.96E-01  9.97E-01 -9.93E-01  5.88E-02 -9.97E-01  1.22E-02
 
 TH 9
+        9.99E-01 -9.94E-01 -9.98E-01  1.00E+00 -9.94E-01  5.20E-02 -1.00E+00  9.97E-01  4.18E-04
 
 TH10
+       -9.98E-01  9.94E-01  9.97E-01 -9.99E-01  9.92E-01 -5.42E-02  9.99E-01 -9.96E-01 -9.99E-01  1.48E-03
 
 TH11
+       -8.31E-01  8.35E-01  8.32E-01 -8.33E-01  8.09E-01 -7.78E-03  8.33E-01 -8.11E-01 -8.33E-01  8.38E-01  8.68E-03
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.91E+09
 
 TH 2
+       -3.79E+08  3.40E+09
 
 TH 3
+        3.31E+07 -2.57E+08  2.57E+07
 
 TH 4
+        1.39E+10 -1.27E+11  9.71E+09  4.79E+12
 
 TH 5
+        9.82E+06  8.81E+07  6.05E+06 -3.25E+09  4.96E+07
 
 TH 6
+       -5.39E+09 -3.27E+09  3.74E+08  1.45E+11  6.51E+08  2.85E+11
 
 TH 7
+        1.39E+10 -8.95E+10  6.45E+09  3.38E+12 -3.45E+09  8.90E+10  2.51E+12
 
 TH 8
+       -1.12E+07  2.31E+07 -7.22E+05 -8.70E+08  3.13E+06 -3.46E+07 -6.34E+08  2.15E+06
 
 TH 9
+       -2.29E+09  1.46E+10 -1.07E+09 -5.47E+11  4.49E+08 -7.20E+09 -3.57E+11  7.23E+07  8.13E+10
 
 TH10
+       -3.09E+07  3.24E+08 -1.10E+07 -1.25E+10  4.75E+07  5.23E+08 -1.16E+10  1.17E+06  1.93E+09  3.80E+08
 
 TH11
+        1.95E+06 -1.05E+07  8.29E+05  3.92E+08  3.77E+04  1.23E+07  2.72E+08 -2.11E+05 -4.25E+07 -1.49E+06  9.46E+04
 
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
 #CPUT: Total CPU Time in Seconds,       44.458
Stop Time:
Sat Sep 25 02:37:30 CDT 2021
