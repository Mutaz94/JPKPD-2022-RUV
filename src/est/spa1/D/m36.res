Thu Sep 30 02:58:55 CDT 2021
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
$DATA ../../../../data/spa1/D/dat36.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16630.9752908454        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6468E+02  2.5233E+02 -3.4098E+01  1.5367E+02  1.1232E+02 -1.1858E+03 -6.8655E+02 -3.2837E+01 -1.1222E+03 -3.2383E+02
            -3.3771E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -703.706097899515        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4187E+00  1.1593E+00  1.0573E+00  1.5835E+00  1.0990E+00  1.9716E+00  1.2916E+00  9.5089E-01  1.4459E+00  1.0397E+00
             1.4370E+01
 PARAMETER:  4.4972E-01  2.4785E-01  1.5576E-01  5.5966E-01  1.9444E-01  7.7885E-01  3.5590E-01  4.9642E-02  4.6874E-01  1.3892E-01
             2.7652E+00
 GRADIENT:   1.3697E+01  9.3520E+00 -5.8621E+00  6.8984E+00 -6.0021E+00  7.1380E+01 -1.0647E+00  4.3932E+00  2.7743E-01  5.1048E+00
             2.5202E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -738.738948287525        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3594E+00  7.8318E-01  2.3724E+00  2.1113E+00  6.3999E+00  1.6019E+00  5.8894E+00  2.0934E-01  1.3950E+00  4.3766E+00
             1.2931E+01
 PARAMETER:  4.0701E-01 -1.4439E-01  9.6390E-01  8.4729E-01  1.9563E+00  5.7120E-01  1.8732E+00 -1.4638E+00  4.3291E-01  1.5763E+00
             2.6596E+00
 GRADIENT:   3.9951E+01  2.4806E+01  5.2850E+00  6.7375E+01  1.2533E+00  1.4948E+00  2.4372E+01 -1.7869E-02  7.2802E+00 -2.5747E-01
             1.9744E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -748.894517922267        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.3187E+00  7.7401E-01  1.6630E+00  1.9433E+00  4.8557E+00  1.5502E+00  5.0405E+00  1.5851E-01  1.2955E+00  8.6522E+00
             1.2641E+01
 PARAMETER:  3.7662E-01 -1.5617E-01  6.0864E-01  7.6436E-01  1.6802E+00  5.3841E-01  1.7175E+00 -1.7419E+00  3.5889E-01  2.2578E+00
             2.6369E+00
 GRADIENT:   3.9656E+01  2.3314E+01  3.6942E+00  5.6745E+01 -6.8145E+00  9.8814E-01  9.7019E+00 -1.4124E-02  5.3143E+00  1.9326E+01
             1.8171E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -764.999486321594        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:      364             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3060E+00  7.5151E-01  1.6233E+00  1.8673E+00  4.7848E+00  1.5474E+00  5.1323E+00  1.8251E-01  1.2449E+00  8.7097E+00
             1.1133E+01
 PARAMETER:  3.6695E-01 -1.8567E-01  5.8449E-01  7.2452E-01  1.6654E+00  5.3658E-01  1.7356E+00 -1.6010E+00  3.1903E-01  2.2644E+00
             2.5099E+00
 GRADIENT:   7.3515E+01  2.2932E+01 -2.6713E+00  6.6035E+01 -1.0277E+01 -1.1437E+01  1.2779E+01 -2.5292E-02 -4.0712E-01  3.2582E+01
             9.3140E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -774.652957397130        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.3093E+00  7.5271E-01  1.6316E+00  1.8467E+00  4.9203E+00  1.5555E+00  5.2021E+00  2.2654E-01  1.2414E+00  8.5522E+00
             1.0054E+01
 PARAMETER:  3.6948E-01 -1.8407E-01  5.8955E-01  7.1338E-01  1.6934E+00  5.4179E-01  1.7491E+00 -1.3849E+00  3.1621E-01  2.2462E+00
             2.4080E+00
 GRADIENT:   8.9061E+01  2.3144E+01 -6.0780E+00  6.8560E+01 -1.1798E+01 -2.3392E+01 -8.6389E+00 -4.1228E-02 -9.0519E+00  3.4882E+01
            -1.6428E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -777.614736328209        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      615
 NPARAMETR:  1.3110E+00  7.5338E-01  1.8408E+00  1.8352E+00  4.8858E+00  1.5604E+00  5.2443E+00  2.4415E-01  1.3590E+00  8.3309E+00
             9.9685E+00
 PARAMETER:  3.7079E-01 -1.8318E-01  7.1017E-01  7.0716E-01  1.6863E+00  5.4496E-01  1.7571E+00 -1.3100E+00  4.0675E-01  2.2200E+00
             2.3994E+00
 GRADIENT:   1.1003E+02  2.2596E+01  1.3033E+00  6.5347E+01 -9.3272E+00 -1.4474E+01  1.7172E+01 -3.9782E-02  1.0131E+00  2.8168E+01
             1.5996E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -779.569444677185        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      764
 NPARAMETR:  1.3104E+00  6.8771E-01  1.8413E+00  1.8352E+00  4.8925E+00  1.5610E+00  5.2468E+00  2.4500E-01  1.3588E+00  8.3019E+00
             9.9749E+00
 PARAMETER:  3.7033E-01 -2.7438E-01  7.1049E-01  7.0718E-01  1.6877E+00  5.4534E-01  1.7576E+00 -1.3065E+00  4.0657E-01  2.2165E+00
             2.4001E+00
 GRADIENT:   8.7034E+01  1.7511E+01 -3.4619E-02  4.4619E+01 -1.0008E+01 -2.4179E+01 -6.9912E+00 -3.9307E-02  3.1230E+00  2.5989E+01
            -4.0723E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -785.208783778791        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      919
 NPARAMETR:  1.2944E+00  6.8947E-01  1.8546E+00  1.7811E+00  4.9920E+00  1.5738E+00  5.1413E+00  4.7953E-01  1.3526E+00  5.2817E+00
             9.4685E+00
 PARAMETER:  3.5807E-01 -2.7183E-01  7.1769E-01  6.7722E-01  1.7078E+00  5.5349E-01  1.7373E+00 -6.3495E-01  4.0204E-01  1.7642E+00
             2.3480E+00
 GRADIENT:   1.1037E+02  2.1073E+01  7.3562E+00  5.2090E+01  8.5041E-01 -6.0663E+00  1.7856E+01 -1.8514E-01 -5.1667E+00 -2.3997E-01
            -9.3142E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -786.507248094113        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  1.2811E+00  6.8901E-01  1.8532E+00  1.6949E+00  4.9938E+00  1.5810E+00  4.6707E+00  4.9486E+00  1.3521E+00  4.8437E+00
             8.9324E+00
 PARAMETER:  3.4774E-01 -2.7250E-01  7.1694E-01  6.2763E-01  1.7082E+00  5.5809E-01  1.6413E+00  1.6991E+00  4.0165E-01  1.6777E+00
             2.2897E+00
 GRADIENT:   1.1014E+02  4.0122E+01 -1.7218E+00  5.5965E+01 -5.2776E+00 -2.8469E+00  5.3535E+00  2.5197E+01  1.1666E+01 -1.8903E+00
            -4.2437E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -786.652967015075        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1062
 NPARAMETR:  1.2792E+00  6.8879E-01  1.8522E+00  1.6855E+00  4.9928E+00  1.5814E+00  4.6281E+00  4.3241E+00  1.3522E+00  5.1013E+00
             8.9785E+00
 PARAMETER:  3.4621E-01 -2.7282E-01  7.1639E-01  6.2206E-01  1.7080E+00  5.5831E-01  1.6321E+00  1.5642E+00  4.0171E-01  1.7295E+00
             2.2948E+00
 GRADIENT:   1.0687E+02  4.0010E+01  4.4184E+00  4.7876E+01 -4.6784E+00 -3.9202E+00  7.1749E+00  1.0010E+01  8.5525E+00 -1.9359E+00
            -3.5994E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -786.838192166162        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1132
 NPARAMETR:  1.2707E+00  6.8801E-01  1.8487E+00  1.6433E+00  4.9933E+00  1.5837E+00  4.4335E+00  3.8400E+00  1.3522E+00  5.9120E+00
             9.0546E+00
 PARAMETER:  3.3957E-01 -2.7395E-01  7.1447E-01  5.9668E-01  1.7081E+00  5.5974E-01  1.5892E+00  1.4455E+00  4.0173E-01  1.8770E+00
             2.3033E+00
 GRADIENT:   9.8193E+01  3.7348E+01  9.8264E+00  3.0006E+01 -4.8348E+00 -3.8418E+00  6.5984E+00 -1.9764E+00  6.9063E+00 -7.0704E-01
            -2.2086E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -786.987124273469        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  1.2562E+00  6.8682E-01  1.8436E+00  1.5724E+00  5.0116E+00  1.5886E+00  4.1254E+00  3.6535E+00  1.3516E+00  6.5157E+00
             9.1268E+00
 PARAMETER:  3.2813E-01 -2.7569E-01  7.1170E-01  5.5263E-01  1.7118E+00  5.6285E-01  1.5172E+00  1.3957E+00  4.0126E-01  1.9742E+00
             2.3112E+00
 GRADIENT:   8.7102E+01  3.2999E+01  1.3425E+01  4.8093E+00 -6.2351E+00 -1.9975E+00  3.5984E+00 -7.2946E+00  6.3329E+00  1.7986E+00
            -5.1347E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -787.538473247026        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1273
 NPARAMETR:  1.2338E+00  6.8504E-01  1.8357E+00  1.4760E+00  5.0773E+00  1.5965E+00  3.7293E+00  3.9416E+00  1.3495E+00  6.3408E+00
             9.1554E+00
 PARAMETER:  3.1012E-01 -2.7827E-01  7.0744E-01  4.8934E-01  1.7248E+00  5.6784E-01  1.4162E+00  1.4716E+00  3.9973E-01  1.9470E+00
             2.3143E+00
 GRADIENT:   7.4050E+01  2.6405E+01  1.2732E+01 -2.2939E+01 -6.5034E+00  4.2286E+00 -4.6555E+00 -1.9109E+00  9.2267E+00 -5.8062E-02
             1.2929E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -792.124443917121        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1409
 NPARAMETR:  1.1311E+00  6.8589E-01  1.8406E+00  1.4836E+00  1.4561E+01  1.5546E+00  3.8190E+00  4.0419E+00  1.3454E+00  7.1555E+00
             8.9989E+00
 PARAMETER:  2.2315E-01 -2.7704E-01  7.1008E-01  4.9445E-01  2.7784E+00  5.4123E-01  1.4400E+00  1.4967E+00  3.9666E-01  2.0679E+00
             2.2971E+00
 GRADIENT:   1.4656E+01  3.9289E+01  1.9750E+01 -7.1130E+00 -1.1795E+00  2.3215E+00  1.1396E+00 -3.3483E+00  8.4182E+00 -5.6312E-01
             7.1744E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -797.860677642295        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1586
 NPARAMETR:  1.1405E+00  6.7578E-01  1.7593E+00  1.5246E+00  1.7042E+04  1.5391E+00  4.2391E+00  5.0650E+00  1.3245E+00  7.1998E+01
             9.4377E+00
 PARAMETER:  2.3149E-01 -2.9189E-01  6.6493E-01  5.2176E-01  9.8434E+00  5.3117E-01  1.5443E+00  1.7224E+00  3.8105E-01  4.3766E+00
             2.3447E+00
 GRADIENT:  -1.5225E+00  3.1780E+01  3.1856E+00 -4.9240E+00 -1.1610E-03 -6.0538E+00 -4.0595E+00  6.0790E+00  1.8244E+01 -3.6462E-05
             2.7172E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -799.177589812792        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1764
 NPARAMETR:  1.1213E+00  6.6950E-01  1.7271E+00  1.5386E+00  2.2846E+05  1.5379E+00  4.5749E+00  5.3837E+00  1.3044E+00  1.6432E+02
             9.1867E+00
 PARAMETER:  2.1452E-01 -3.0122E-01  6.4645E-01  5.3090E-01  1.2439E+01  5.3041E-01  1.6206E+00  1.7834E+00  3.6572E-01  5.2018E+00
             2.3178E+00
 GRADIENT:  -1.3558E+00  4.1757E+01  4.5339E+00 -3.1293E+00 -9.2893E-05 -2.1103E+00  2.1034E+00  6.2606E+00  1.5264E+01 -1.0331E-06
             1.3618E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -806.292229660843        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1944
 NPARAMETR:  1.1132E+00  5.4376E-01  1.3750E+00  1.6825E+00  1.3291E+12  1.4940E+00  5.7164E+00  4.5522E+00  9.5314E-01  2.3044E+04
             9.5024E+00
 PARAMETER:  2.0727E-01 -5.0924E-01  4.1847E-01  6.2030E-01  2.8016E+01  5.0146E-01  1.8433E+00  1.6156E+00  5.2008E-02  1.0145E+01
             2.3515E+00
 GRADIENT:  -6.4094E+00  4.1653E+01  2.1743E+00  7.6202E+01 -6.8986E-12 -2.3147E+00 -1.8846E+00 -1.7050E+00 -9.1982E+00  0.0000E+00
             4.3381E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -815.182644659641        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2123
 NPARAMETR:  1.0907E+00  4.4943E-01  1.1451E+00  1.4875E+00  1.5007E+16  1.4638E+00  6.4546E+00  4.1469E+00  7.8178E-01  4.2961E+05
             9.2104E+00
 PARAMETER:  1.8679E-01 -6.9977E-01  2.3545E-01  4.9711E-01  3.7347E+01  4.8105E-01  1.9648E+00  1.5223E+00 -1.4618E-01  1.3071E+01
             2.3203E+00
 GRADIENT:   1.3962E+01  2.7820E+01  8.8421E+00 -1.3052E+01  0.0000E+00  1.8305E+00  1.7320E+01 -6.5912E+00 -6.2012E+00  0.0000E+00
            -4.7894E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -850.879047639697        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     2265
 NPARAMETR:  9.9276E-01  6.4475E-02  6.7597E-01  1.4052E+00  6.1868E+16  1.3121E+00  6.4306E+00  3.1924E+00  7.1876E-01  6.1299E+05
             8.8452E+00
 PARAMETER:  9.2738E-02 -2.6415E+00 -2.9160E-01  4.4019E-01  3.8764E+01  3.7165E-01  1.9611E+00  1.2608E+00 -2.3023E-01  1.3426E+01
             2.2799E+00
 GRADIENT:   1.1461E+01  1.3461E+00  3.0418E+01  2.7472E+00  0.0000E+00 -1.0734E+01 -3.0592E-02 -1.1514E+01 -7.9005E-01  0.0000E+00
             5.9634E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -865.573454425393        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     2336
 NPARAMETR:  8.1650E-01  2.0138E-02  1.7155E-01  1.0110E+00  6.2342E+16  1.3042E+00  6.4144E+00  2.2609E+00  5.5809E-01  6.0694E+05
             8.4066E+00
 PARAMETER: -1.0273E-01 -3.8051E+00 -1.6629E+00  1.1095E-01  3.8771E+01  3.6557E-01  1.9585E+00  9.1576E-01 -4.8324E-01  1.3416E+01
             2.2290E+00
 GRADIENT:   8.5622E+00  5.1101E-01 -2.2214E+01  1.0309E+02  0.0000E+00 -5.8890E+00 -9.7871E-03  1.1131E+01 -8.4395E-02  0.0000E+00
            -5.4224E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -866.273489519228        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     2408
 NPARAMETR:  7.4836E-01  1.3931E-02  1.1589E-01  8.4181E-01  6.2506E+16  1.3075E+00  6.4082E+00  2.0628E+00  4.5718E-01  6.0487E+05
             8.3497E+00
 PARAMETER: -1.8987E-01 -4.1737E+00 -2.0551E+00 -7.2203E-02  3.8774E+01  3.6812E-01  1.9576E+00  8.2408E-01 -6.8267E-01  1.3413E+01
             2.2222E+00
 GRADIENT:   2.0835E+01  2.3953E-01 -2.9865E+01  1.0498E+02  0.0000E+00 -2.1608E+00 -2.5425E-03  1.4944E+01  5.5336E-02  0.0000E+00
            -5.6871E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -874.025381084650        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2587
 NPARAMETR:  7.1232E-01  1.0415E-02  1.0147E-01  7.2486E-01  6.2640E+16  1.3461E+00  6.4034E+00  1.6860E+00  3.4716E-01  6.0318E+05
             9.0155E+00
 PARAMETER: -2.3923E-01 -4.4645E+00 -2.1880E+00 -2.2178E-01  3.8776E+01  3.9718E-01  1.9568E+00  6.2238E-01 -9.5797E-01  1.3410E+01
             2.2989E+00
 GRADIENT:   5.3054E-01 -2.6909E-02  4.6379E-01 -1.0541E+00  0.0000E+00  4.4945E+00 -1.0070E-02 -1.0322E+00 -6.0805E-01  0.0000E+00
             1.2916E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -874.187463092494        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2763
 NPARAMETR:  7.1993E-01  1.1707E-02  1.0660E-01  7.4450E-01  6.2593E+16  1.3247E+00  6.4049E+00  1.6982E+00  4.1199E-01  6.0376E+05
             8.8888E+00
 PARAMETER: -2.2860E-01 -4.3476E+00 -2.1386E+00 -1.9504E-01  3.8775E+01  3.8117E-01  1.9571E+00  6.2956E-01 -7.8675E-01  1.3411E+01
             2.2848E+00
 GRADIENT:   4.7449E-01 -3.5362E-02  4.8530E-01  4.6960E-01  0.0000E+00 -6.6273E-01 -1.1082E-02 -7.3347E-01 -1.8114E-01  0.0000E+00
            -2.7432E-01

0ITERATION NO.:  118    OBJECTIVE VALUE:  -874.221702612459        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     2856
 NPARAMETR:  7.1265E-01  1.1593E-02  1.0273E-01  7.2986E-01  6.2601E+16  1.3264E+00  6.4047E+00  1.6777E+00  4.1816E-01  6.0367E+05
             8.8874E+00
 PARAMETER: -2.3877E-01 -4.3574E+00 -2.1756E+00 -2.1490E-01  3.8776E+01  3.8249E-01  1.9570E+00  6.1743E-01 -7.7189E-01  1.3411E+01
             2.2846E+00
 GRADIENT:   8.6921E-02 -4.0325E-02 -1.0051E-01  1.9540E-01  0.0000E+00  7.9701E-02 -9.8501E-03 -2.9257E-02 -5.1903E-03  0.0000E+00
             2.8647E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2856
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1302E-03 -1.3699E-03 -3.0914E-03 -1.7705E-02  1.1495E-21
 SE:             2.8340E-02  7.1552E-04  2.2646E-02  1.0846E-02  1.0730E-21
 N:                     100         100         100         100         100

 P VAL.:         8.0135E-01  5.5549E-02  8.9142E-01  1.0260E-01  2.8407E-01

 ETASHRINKSD(%)  5.0579E+00  9.7603E+01  2.4133E+01  6.3665E+01  1.0000E+02
 ETASHRINKVR(%)  9.8599E+00  9.9943E+01  4.2442E+01  8.6798E+01  1.0000E+02
 EBVSHRINKSD(%)  5.3769E+00  9.7664E+01  2.3101E+01  6.4455E+01  1.0000E+02
 EBVSHRINKVR(%)  1.0465E+01  9.9945E+01  4.0866E+01  8.7365E+01  1.0000E+02
 RELATIVEINF(%)  1.5554E+01  1.5278E-03  8.5035E-01  1.0092E-01  0.0000E+00
 EPSSHRINKSD(%)  1.0574E+01
 EPSSHRINKVR(%)  2.0030E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -874.22170261245878     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       44.716830592213910     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    62.25
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -874.222       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         7.13E-01  1.16E-02  1.03E-01  7.30E-01  6.26E+16  1.33E+00  6.40E+00  1.68E+00  4.18E-01  6.04E+05  8.89E+00
 


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
+       -2.82E+02 -1.04E+02
 
 TH 3
+       -6.46E+02 -4.72E+02  2.13E+04
 
 TH 4
+       -3.62E+02  2.98E+02 -4.90E+03  1.48E+03
 
 TH 5
+       -1.27E-20 -1.09E-21 -1.25E-19  3.24E-19 -8.02E-39
 
 TH 6
+        5.43E+00  7.68E+00  4.90E+01 -2.74E+01  1.20E-19  8.90E+01
 
 TH 7
+       -8.62E-03 -2.50E-01 -5.04E-03  6.06E-03 -1.74E-21  3.12E-02 -6.12E-04
 
 TH 8
+        5.67E+00 -1.53E+01 -8.35E+01 -2.27E+01 -2.42E-20  3.47E+00  4.42E-04  2.42E+01
 
 TH 9
+        3.79E+00 -1.57E+01  7.15E+01 -4.47E+01  5.76E-20  3.86E+00 -1.69E-03  1.49E+01  2.49E+01
 
 TH10
+        3.88E-02 -1.07E-02 -1.20E-02 -1.59E-02  5.65E-28 -8.86E-08  4.43E-09  1.01E-07  6.07E-08  2.92E-08
 
 TH11
+       -1.79E+01  7.42E-01  6.91E+01 -1.28E+01  1.79E-21  1.15E+00  4.18E-04  1.51E+00  2.07E+00 -1.31E-03  6.70E+00
 
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
 #CPUT: Total CPU Time in Seconds,       75.260
Stop Time:
Thu Sep 30 03:00:12 CDT 2021
