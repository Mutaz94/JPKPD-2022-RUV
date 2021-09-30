Thu Sep 30 02:50:12 CDT 2021
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
$DATA ../../../../data/spa1/D/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   7769.01143183854        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.2568E+02  3.1009E+01 -1.1172E+02 -1.0860E+02  2.1565E+02 -1.1952E+03 -4.2731E+02 -1.7401E+01 -8.6762E+02 -3.6373E+02
            -1.6759E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -815.782591133503        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3021E+00  9.6584E-01  1.2096E+00  2.0065E+00  1.1433E+00  2.6153E+00  1.6359E+00  1.0018E+00  3.4858E+00  1.3617E+00
             1.1706E+01
 PARAMETER:  3.6396E-01  6.5245E-02  2.9027E-01  7.9638E-01  2.3394E-01  1.0614E+00  5.9219E-01  1.0184E-01  1.3487E+00  4.0870E-01
             2.5601E+00
 GRADIENT:  -9.6626E+00 -1.1938E+01 -2.6858E+01  2.9230E+01 -8.4027E+00  7.5967E+01  8.0472E+00  5.8827E+00  9.0058E+01  8.5731E+00
             3.7682E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -866.773246504081        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.3214E+00  7.3034E-01  1.0054E+01  2.9331E+00  1.1559E+01  1.9730E+00  7.4513E+00  1.0395E+00  3.6817E+00  1.1651E+01
             9.7605E+00
 PARAMETER:  3.7869E-01 -2.1425E-01  2.4079E+00  1.1761E+00  2.5475E+00  7.7956E-01  2.1084E+00  1.3871E-01  1.4034E+00  2.5554E+00
             2.3783E+00
 GRADIENT:   5.4550E+01  6.7551E+00 -8.6034E-01  7.1591E+01 -4.6517E+00  2.9270E+01  1.6404E+01  1.3509E-02  7.9234E+01  1.4644E+01
             2.8487E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -955.181146100174        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0127E+00  1.4133E+00  1.3891E+01  1.2847E+00  1.3395E+01  2.2342E+00  3.3965E+00  5.8721E+00  2.7095E+00  1.1995E+01
             7.9989E+00
 PARAMETER:  1.1265E-01  4.4595E-01  2.7312E+00  3.5054E-01  2.6949E+00  9.0386E-01  1.3227E+00  1.8702E+00  1.0968E+00  2.5845E+00
             2.1793E+00
 GRADIENT:  -7.4515E+01  7.8783E+00  6.3184E+00 -5.0851E+00 -2.1102E+00  5.1854E+01  3.3014E+01  2.4825E+01  3.0039E+01  3.1824E+00
             2.5891E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -982.491457767470        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:      358             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1068E+00  1.1809E+00  1.0321E+01  1.1451E+00  1.0060E+01  2.3439E+00  3.0454E+00  8.4448E+00  2.4443E+00  1.0126E+01
             6.8169E+00
 PARAMETER:  2.0146E-01  2.6625E-01  2.4342E+00  2.3547E-01  2.4086E+00  9.5182E-01  1.2136E+00  2.2335E+00  9.9374E-01  2.4151E+00
             2.0194E+00
 GRADIENT:  -6.6863E+00 -1.6276E+01  4.8236E+00 -5.0432E+01 -3.1821E+00  9.7436E+01  2.8707E+01  1.0633E+02  9.2480E-01  9.2244E+00
             1.4244E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -993.814554827456        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  1.1160E+00  1.1802E+00  1.0229E+01  1.1459E+00  1.0354E+01  2.3305E+00  3.0325E+00  4.8589E+00  2.4498E+00  9.7463E+00
             6.7502E+00
 PARAMETER:  2.0979E-01  2.6567E-01  2.4252E+00  2.3615E-01  2.4373E+00  9.4608E-01  1.2094E+00  1.6808E+00  9.9600E-01  2.3769E+00
             2.0096E+00
 GRADIENT:  -2.2997E+01 -1.3309E+01  8.0795E+00 -4.9353E+01 -1.7560E+00  3.5674E+01  8.3513E+00 -1.5356E+01 -4.9719E+00  1.2307E+00
             1.2147E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -997.991729081107        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.1414E+00  1.2129E+00  3.1205E+00  1.2493E+00  1.0527E+01  1.9441E+00  3.0438E+00  5.7672E+00  2.5121E+00  9.6038E+00
             6.6550E+00
 PARAMETER:  2.3228E-01  2.9303E-01  1.2380E+00  3.2261E-01  2.4540E+00  7.6480E-01  1.2131E+00  1.8522E+00  1.0211E+00  2.3622E+00
             1.9954E+00
 GRADIENT:   1.1614E+01  1.3803E+01  7.8344E-02 -1.1272E+01 -3.2393E+00  2.2068E+01  2.5257E+01  5.9439E+01  2.4272E+01  8.7509E+00
             1.0372E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1014.20970111522        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      664
 NPARAMETR:  1.1369E+00  8.2748E-01  1.3734E+00  1.4788E+00  1.4636E+01  1.6981E+00  3.0867E+00  2.2484E+00  1.8489E+00  9.0123E+00
             6.2262E+00
 PARAMETER:  2.2830E-01 -8.9369E-02  4.1727E-01  4.9124E-01  2.7835E+00  6.2952E-01  1.2271E+00  9.1021E-01  7.1457E-01  2.2986E+00
             1.9288E+00
 GRADIENT:   3.8439E+01  4.0169E+01  7.2171E+00  2.2766E+01 -1.3523E+00 -3.6234E+01  7.0510E+00 -9.5205E+00 -6.7375E+00 -4.3713E-01
             2.3076E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1046.79538995515        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      737
 NPARAMETR:  1.1241E+00  7.1765E-02  1.7003E-01  9.7014E-01  9.2476E+01  1.8942E+00  3.1592E+00  3.0772E-01  1.2664E+00  7.1838E+00
             5.3058E+00
 PARAMETER:  2.1702E-01 -2.5344E+00 -1.6718E+00  6.9683E-02  4.6269E+00  7.3881E-01  1.2503E+00 -1.0786E+00  3.3619E-01  2.0718E+00
             1.7688E+00
 GRADIENT:   2.4580E+02  4.0164E+00  1.7875E+01 -1.0467E+01  1.1306E-01 -3.1153E+00  1.2251E+00 -1.7466E+00  8.3366E+00  1.6726E-04
            -1.5171E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1046.85193312971        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      811
 NPARAMETR:  1.1231E+00  6.4081E-02  1.5183E-01  9.2230E-01  9.9802E+01  1.9171E+00  3.1603E+00  2.8944E-01  1.2352E+00  7.1252E+00
             5.3638E+00
 PARAMETER:  2.1612E-01 -2.6476E+00 -1.7850E+00  1.9112E-02  4.7032E+00  7.5082E-01  1.2507E+00 -1.1398E+00  3.1126E-01  2.0636E+00
             1.7797E+00
 GRADIENT:   2.5577E+02  3.4956E+00  1.5390E+01 -8.6730E+00  1.1357E-01 -2.7407E+00  1.0743E+00 -1.4729E+00  7.9273E+00  1.7379E-04
            -1.3800E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1061.15973720487        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      950
 NPARAMETR:  1.1218E+00  8.1640E-02  1.8844E-01  1.0100E+00  8.4699E+01  2.1487E+00  3.1450E+00  3.7117E-01  1.2729E+00  7.3083E+00
             5.8480E+00
 PARAMETER:  2.1491E-01 -2.4054E+00 -1.5690E+00  1.0991E-01  4.5391E+00  8.6487E-01  1.2458E+00 -8.9109E-01  3.4127E-01  2.0890E+00
             1.8661E+00
 GRADIENT:   1.4265E+02  3.6754E+00  6.1405E+00 -3.5134E+01  6.6648E-02  1.2689E+01  1.1750E+00 -2.2822E+00  3.2715E+00  3.0527E-05
            -4.4238E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1069.48912941709        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1128
 NPARAMETR:  1.1063E+00  6.9025E-02  1.7086E-01  1.0405E+00  6.1308E+01  1.7549E+00  3.0816E+00  1.2827E+00  1.0884E+00  7.7162E+00
             6.1812E+00
 PARAMETER:  2.0106E-01 -2.5733E+00 -1.6669E+00  1.3975E-01  4.2159E+00  6.6243E-01  1.2255E+00  3.4894E-01  1.8474E-01  2.1433E+00
             1.9215E+00
 GRADIENT:   2.0956E+02  1.6407E+00 -5.2466E+01  5.6500E+01  2.4746E-01 -4.5102E+01  8.1030E-01 -1.0712E+01  1.2858E+01 -4.6244E-03
            -9.2252E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1088.61421390637        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1305
 NPARAMETR:  1.0932E+00  6.1846E-02  3.2761E-01  1.2093E+00  5.7665E+01  1.7921E+00  3.0399E+00  3.1363E+00  7.9654E-01  7.8471E+00
             6.2104E+00
 PARAMETER:  1.8911E-01 -2.6831E+00 -1.0159E+00  2.9000E-01  4.1546E+00  6.8337E-01  1.2118E+00  1.2430E+00 -1.2747E-01  2.1601E+00
             1.9262E+00
 GRADIENT:   1.3473E+02  2.3123E+00 -3.0359E+00 -4.7020E+00  1.5619E-01  1.2026E+01  5.7901E-01  1.0313E+01  3.5544E+00 -5.0887E-03
             1.2383E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1094.04884952376        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1482
 NPARAMETR:  1.0332E+00  1.0000E-02  2.3382E-01  1.0557E+00  1.7686E+02  1.7744E+00  3.0139E+00  2.8040E+00  3.7594E-01  6.8391E+00
             6.2498E+00
 PARAMETER:  1.3266E-01 -4.8180E+00 -1.3532E+00  1.5425E-01  5.2754E+00  6.7344E-01  1.2032E+00  1.1310E+00 -8.7832E-01  2.0227E+00
             1.9326E+00
 GRADIENT:   1.5216E+02  0.0000E+00  2.5582E+00 -1.5263E+01  7.0801E-02  1.0193E+01 -8.4399E-03  1.1381E+00 -4.8543E+00 -7.7068E-04
             1.4142E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1129.85167693040        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1659
 NPARAMETR:  6.1362E-01  1.0000E-02  6.8066E-02  5.6622E-01  2.2879E+06  1.5285E+00  2.7650E+00  1.7162E+00  1.0000E-02  2.1030E+00
             5.9673E+00
 PARAMETER: -3.8837E-01 -2.3316E+01 -2.5873E+00 -4.6878E-01  1.4743E+01  5.2431E-01  1.1170E+00  6.4010E-01 -7.0343E+00  8.4335E-01
             1.8863E+00
 GRADIENT:   9.3464E+00  0.0000E+00  1.0945E+01 -2.9449E+01  2.1884E-07  3.2712E-02  3.7874E-03  4.1747E+00  0.0000E+00  0.0000E+00
             6.5577E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1130.05820225155        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1853             RESET HESSIAN, TYPE I
 NPARAMETR:  6.2254E-01  1.0000E-02  7.2493E-02  5.9598E-01  1.3412E+06  1.5329E+00  7.9409E-02  1.7401E+00  1.0000E-02  2.1343E+00
             5.9430E+00
 PARAMETER: -3.7395E-01 -2.2674E+01 -2.5243E+00 -4.1754E-01  1.4209E+01  5.2713E-01 -2.4331E+00  6.5394E-01 -6.8025E+00  8.5812E-01
             1.8822E+00
 GRADIENT:   6.2037E+01  0.0000E+00  7.4141E+01  2.5738E+01  9.4638E-07  2.0156E+01  7.0663E-06  5.6264E+00  0.0000E+00  0.0000E+00
             1.8429E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1130.06394087330        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     2025
 NPARAMETR:  6.2131E-01  1.0000E-02  7.2129E-02  5.9467E-01  1.5003E+05  1.5316E+00  8.9684E-02  1.7394E+00  1.0000E-02  2.0271E+00
             5.9427E+00
 PARAMETER: -3.7592E-01 -2.2674E+01 -2.5293E+00 -4.1974E-01  1.2019E+01  5.2629E-01 -2.3115E+00  6.5352E-01 -6.8025E+00  8.0659E-01
             1.8822E+00
 GRADIENT:   6.1647E+01  0.0000E+00  7.3690E+01  2.7240E+01  8.4184E-06  1.9836E+01  8.8176E-06  5.8677E+00  0.0000E+00  0.0000E+00
             1.8423E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1130.16677834610        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     2146
 NPARAMETR:  6.2053E-01  1.0000E-02  7.1843E-02  5.9349E-01  1.1182E+01  1.5304E+00  8.8603E-02  1.7346E+00  1.0000E-02  2.1681E+00
             5.9398E+00
 PARAMETER: -3.7717E-01 -2.2674E+01 -2.5333E+00 -4.2174E-01  2.5143E+00  5.2556E-01 -2.3236E+00  6.5075E-01 -6.8025E+00  8.7383E-01
             1.8817E+00
 GRADIENT:  -3.4348E+00  0.0000E+00  2.6219E-01 -1.5922E+00  1.4792E-02 -7.5318E-01  3.5318E-06  1.1920E+00  0.0000E+00 -1.4102E-03
             2.0580E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1130.17394140527        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2330
 NPARAMETR:  6.2486E-01  1.0000E-02  7.2068E-02  5.9738E-01  1.1060E+01  1.5348E+00  8.1514E-02  1.7303E+00  1.0000E-02  5.1399E+00
             5.9407E+00
 PARAMETER: -3.7023E-01 -2.2674E+01 -2.5301E+00 -4.1520E-01  2.5033E+00  5.2841E-01 -2.4070E+00  6.4828E-01 -6.8025E+00  1.7370E+00
             1.8818E+00
 GRADIENT:   9.4242E-02  0.0000E+00 -4.1667E+00  5.2242E+00  2.8230E-02  1.4312E-01  2.9061E-06 -1.1496E-01  0.0000E+00 -3.9969E-03
            -1.0366E+00

0ITERATION NO.:   94    OBJECTIVE VALUE:  -1130.18025733008        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     2465
 NPARAMETR:  6.2491E-01  1.0000E-02  7.2195E-02  5.9628E-01  1.0954E+01  1.5349E+00  8.1367E-02  1.7307E+00  1.0000E-02  5.2469E+00
             5.9431E+00
 PARAMETER: -3.7014E-01 -2.2674E+01 -2.5284E+00 -4.1704E-01  2.4937E+00  5.2850E-01 -2.4088E+00  6.4855E-01 -6.8025E+00  1.7576E+00
             1.8822E+00
 GRADIENT:   7.3317E-01  0.0000E+00 -1.3298E+00 -1.3513E-01 -1.2213E-02  2.6441E-01  2.8791E-06  2.0900E-01  0.0000E+00  2.1456E-03
            -1.3604E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2465
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.5902E-03 -2.1137E-05  1.2169E-03 -5.1517E-04 -1.9685E-03
 SE:             2.9247E-02  8.5272E-06  2.6370E-02  2.8743E-04  8.9972E-04
 N:                     100         100         100         100         100

 P VAL.:         8.4842E-01  1.3184E-02  9.6319E-01  7.3082E-02  2.8673E-02

 ETASHRINKSD(%)  2.0190E+00  9.9971E+01  1.1658E+01  9.9037E+01  9.6986E+01
 ETASHRINKVR(%)  3.9972E+00  1.0000E+02  2.1957E+01  9.9991E+01  9.9909E+01
 EBVSHRINKSD(%)  2.2998E+00  9.9971E+01  1.0306E+01  9.9076E+01  9.7655E+01
 EBVSHRINKVR(%)  4.5467E+00  1.0000E+02  1.9550E+01  9.9991E+01  9.9945E+01
 RELATIVEINF(%)  6.8835E+00  2.7355E-07  8.8916E-01  4.9163E-05  6.8676E-03
 EPSSHRINKSD(%)  1.4807E+01
 EPSSHRINKVR(%)  2.7421E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1130.1802573300756     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -211.24172412540292     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.92
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1130.180       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.25E-01  1.00E-02  7.22E-02  5.96E-01  1.10E+01  1.53E+00  8.14E-02  1.73E+00  1.00E-02  5.25E+00  5.94E+00
 


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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -4.39E+02  0.00E+00  5.62E+04
 
 TH 4
+       -5.07E+02  0.00E+00 -1.10E+04  2.67E+03
 
 TH 5
+        2.46E-01  0.00E+00 -1.93E+00  3.97E-01  1.03E-02
 
 TH 6
+        4.23E+00  0.00E+00  4.96E+01 -1.81E+01  1.08E-02  7.39E+01
 
 TH 7
+        2.11E-02  0.00E+00 -1.38E-02 -1.33E-03  6.99E-06 -1.94E-03 -1.81E-02
 
 TH 8
+        6.97E+00  0.00E+00 -8.34E+01 -6.29E+01 -3.85E-02  1.45E+00 -1.21E-03  4.03E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -5.75E-04  0.00E+00  9.86E-02 -3.93E-02 -3.48E-03  3.15E-03 -5.31E-04  5.15E-03  0.00E+00  4.74E-03
 
 TH11
+       -1.34E+01  0.00E+00  9.47E+01 -2.04E+01 -7.87E-03  1.38E+00 -1.46E-06  3.10E+00  0.00E+00  1.08E-04  1.57E+01
 
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
 #CPUT: Total CPU Time in Seconds,       52.982
Stop Time:
Thu Sep 30 02:51:07 CDT 2021
