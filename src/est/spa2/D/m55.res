Thu Sep 30 09:17:32 CDT 2021
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
$DATA ../../../../data/spa2/D/dat55.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25870.2655720976        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.4998E+02  4.6617E+02 -1.2709E+01  3.6911E+02  2.0946E+02 -2.1577E+03 -9.2914E+02 -6.4851E+01 -1.4574E+03 -6.7935E+02
            -5.0806E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -569.946811240841        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3662E+00  1.2689E+00  9.6961E-01  1.5557E+00  1.0382E+00  2.0769E+00  1.3704E+00  9.6896E-01  1.2202E+00  1.0361E+00
             1.4433E+01
 PARAMETER:  4.1203E-01  3.3812E-01  6.9135E-02  5.4193E-01  1.3749E-01  8.3088E-01  4.1507E-01  6.8464E-02  2.9899E-01  1.3547E-01
             2.7695E+00
 GRADIENT:   1.6419E+01  2.8042E-01 -1.7716E+01  3.9448E+01  2.1999E+01  5.9431E+01 -1.2341E+01  4.8551E+00 -1.4037E+01  1.4442E+01
             1.6499E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -646.523432181214        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3463E+00  7.4621E-01  1.2696E+01  2.4562E+00  2.5605E+00  1.9015E+00  6.0782E+00  1.8745E-01  1.4630E+00  3.0761E-01
             1.4392E+01
 PARAMETER:  3.9735E-01 -1.9274E-01  2.6413E+00  9.9863E-01  1.0402E+00  7.4265E-01  1.9047E+00 -1.5742E+00  4.8047E-01 -1.0789E+00
             2.7666E+00
 GRADIENT:   2.8631E+01  1.1584E+01 -1.2073E+00  8.5436E+01 -8.2430E-01  4.1488E+01  2.1954E+01  3.5630E-03 -2.8872E+00  3.5203E-01
             2.1490E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -675.354929179634        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0631E+00  9.8527E-01  1.2719E+01  1.3859E+00  2.5910E+00  1.7544E+00  4.5333E+00  1.1737E+00  1.2535E+00  1.1682E+00
             1.1339E+01
 PARAMETER:  1.6118E-01  8.5165E-02  2.6431E+00  4.2635E-01  1.0520E+00  6.6215E-01  1.6115E+00  2.6016E-01  3.2596E-01  2.5543E-01
             2.5283E+00
 GRADIENT:  -4.2386E+01 -2.0746E+01 -7.0837E-01 -4.4775E+01  9.7382E+00  8.2930E+00  7.8317E+00  4.0276E-02  1.4509E+01  7.5059E+00
             1.3667E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -683.789045236375        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.1104E+00  1.3873E+00  4.2050E+00  1.1971E+00  1.9954E+00  1.6456E+00  4.1125E+00  6.9427E-01  8.1193E-01  4.2730E-01
             1.0385E+01
 PARAMETER:  2.0470E-01  4.2739E-01  1.5363E+00  2.7989E-01  7.9086E-01  5.9811E-01  1.5140E+00 -2.6490E-01 -1.0835E-01 -7.5027E-01
             2.4403E+00
 GRADIENT:   5.7825E+00 -7.7940E-01 -2.6930E+00 -1.7503E+00  1.3178E+00 -8.0821E+00  4.2483E-01  1.7541E-01  2.6687E+00  1.4448E+00
            -2.7823E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -692.962134603225        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  1.1380E+00  8.1886E-01  1.5623E+01  1.7369E+00  2.3712E+00  1.7535E+00  6.3599E+00  2.8218E-01  9.1240E-01  1.6985E-01
             1.0923E+01
 PARAMETER:  2.2926E-01 -9.9845E-02  2.8488E+00  6.5209E-01  9.6340E-01  6.6164E-01  1.9500E+00 -1.1652E+00  8.3187E-03 -1.6729E+00
             2.4908E+00
 GRADIENT:  -4.4353E+00  1.0069E+01 -6.8663E-01  3.4035E+01 -9.0969E-01  2.1300E+00  5.9357E+00  3.1875E-03 -1.1521E+01  1.5745E-01
             5.0878E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -697.105834659852        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  1.1228E+00  4.4047E-01  5.2806E+01  1.9270E+00  2.4606E+00  1.7402E+00  6.8961E+00  8.1790E-01  1.2616E+00  1.4829E-01
             1.0670E+01
 PARAMETER:  2.1582E-01 -7.1992E-01  4.0666E+00  7.5594E-01  1.0004E+00  6.5398E-01  2.0310E+00 -1.0102E-01  3.3237E-01 -1.8086E+00
             2.4674E+00
 GRADIENT:  -4.4210E+00  2.1397E+00 -8.2650E-02  8.0650E+00  2.1271E+00  2.1640E+00 -4.0504E-01  2.8195E-03 -2.9238E-01  1.2921E-01
            -4.0095E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -698.033325823712        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      794
 NPARAMETR:  1.1359E+00  2.9557E-01  1.1161E+02  2.0286E+00  2.4647E+00  1.7288E+00  7.6250E+00  1.1049E+00  1.3554E+00  9.1902E-02
             1.0691E+01
 PARAMETER:  2.2744E-01 -1.1188E+00  4.8150E+00  8.0734E-01  1.0021E+00  6.4740E-01  2.1314E+00  1.9978E-01  4.0408E-01 -2.2870E+00
             2.4694E+00
 GRADIENT:   1.7077E+00  1.1088E-01 -4.3309E-03  4.2702E+00  1.4389E+00  7.3865E-01  1.1842E+00  1.4345E-03  3.0918E+00  4.9926E-02
            -4.8952E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -698.315628129575        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      955
 NPARAMETR:  1.1313E+00  2.9140E-01  9.1051E+01  2.0025E+00  2.4110E+00  1.7235E+00  7.8553E+00  1.0958E+00  1.2889E+00  3.0654E-02
             1.0719E+01
 PARAMETER:  2.2336E-01 -1.1331E+00  4.6114E+00  7.9441E-01  9.8005E-01  6.4436E-01  2.1612E+00  1.9145E-01  3.5376E-01 -3.3850E+00
             2.4720E+00
 GRADIENT:   3.8737E-01  1.3218E+00  1.3270E-02 -1.3385E+00  2.3383E-01  6.6994E-01  9.0550E+00  2.0900E-03 -6.6807E-01  5.6748E-03
             2.3988E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -698.447133796619        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1140             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1308E+00  2.6969E-01  6.1047E+01  2.0117E+00  2.3906E+00  1.7213E+00  8.0341E+00  1.0745E+00  1.3004E+00  1.0000E-02
             1.0705E+01
 PARAMETER:  2.2288E-01 -1.2105E+00  4.2116E+00  7.9899E-01  9.7155E-01  6.4306E-01  2.1837E+00  1.7185E-01  3.6269E-01 -3.0631E+01
             2.4708E+00
 GRADIENT:   1.0589E+01  4.5232E+00  5.3797E-03  1.8721E+01  1.1712E+00  9.7427E+00  1.1209E+02  4.7298E-03  1.3116E+00  0.0000E+00
             3.0701E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -698.489817910927        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1306
 NPARAMETR:  1.1303E+00  2.5285E-01  5.6770E+01  2.0204E+00  2.3831E+00  1.7214E+00  8.0402E+00  1.0543E+00  1.3018E+00  1.0000E-02
             1.0702E+01
 PARAMETER:  2.2252E-01 -1.2750E+00  4.1390E+00  8.0329E-01  9.6842E-01  6.4314E-01  2.1845E+00  1.5284E-01  3.6372E-01 -3.0937E+01
             2.4704E+00
 GRADIENT:  -2.0798E-02 -8.1068E-02  6.2105E-03 -1.3034E+00  1.3412E-01  3.4162E-01  6.9047E+00  5.4956E-03 -2.3428E-01  0.0000E+00
            -9.3508E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -698.613040683831        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1494             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1321E+00  2.3184E-01  5.3279E+01  2.0392E+00  2.3816E+00  1.7165E+00  8.4277E+00  1.0171E+00  1.3168E+00  1.0000E-02
             1.0714E+01
 PARAMETER:  2.2410E-01 -1.3617E+00  4.0755E+00  8.1257E-01  9.6777E-01  6.4031E-01  2.2315E+00  1.1693E-01  3.7518E-01 -3.0937E+01
             2.4715E+00
 GRADIENT:   1.1002E+01  4.4122E+00  7.2746E-05  1.8613E+01  1.2485E+00  8.8414E+00  1.2293E+02  5.9944E-03  1.8320E+00  0.0000E+00
             3.0238E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -698.651620162898        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1678             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1317E+00  2.2215E-01  5.1739E+01  2.0541E+00  2.3739E+00  1.7204E+00  8.5140E+00  1.0091E+00  1.3188E+00  1.0000E-02
             1.0722E+01
 PARAMETER:  2.2374E-01 -1.4044E+00  4.0462E+00  8.1986E-01  9.6455E-01  6.4253E-01  2.2417E+00  1.0905E-01  3.7675E-01 -3.0937E+01
             2.4723E+00
 GRADIENT:   1.0261E+01  4.2391E+00  9.0667E-03  2.2053E+01  6.3872E-01  9.6731E+00  1.2431E+02  6.4611E-03  1.2690E+00  0.0000E+00
             2.9350E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -698.662504550291        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1857
 NPARAMETR:  1.1326E+00  2.1553E-01  5.0582E+01  2.0601E+00  2.3789E+00  1.7209E+00  8.5223E+00  9.8879E-01  1.3254E+00  1.0000E-02
             1.0732E+01
 PARAMETER:  2.2456E-01 -1.4347E+00  4.0236E+00  8.2276E-01  9.6662E-01  6.4286E-01  2.2427E+00  8.8725E-02  3.8168E-01 -3.0937E+01
             2.4732E+00
 GRADIENT:   1.9517E-01  5.9728E-02 -1.5505E-03 -8.4961E-01  1.0243E-01  3.6155E-01  1.1501E+01  6.5881E-03  7.1740E-02  0.0000E+00
             7.4225E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -698.674748799944        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2044             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1323E+00  2.1121E-01  5.0506E+01  2.0654E+00  2.3752E+00  1.7199E+00  8.6326E+00  9.7661E-01  1.3259E+00  1.0000E-02
             1.0730E+01
 PARAMETER:  2.2424E-01 -1.4549E+00  4.0221E+00  8.2534E-01  9.6507E-01  6.4227E-01  2.2555E+00  7.6332E-02  3.8207E-01 -3.0937E+01
             2.4730E+00
 GRADIENT:   1.0242E+01  4.0677E+00  6.3676E-03  2.2666E+01  7.1535E-01  9.6598E+00  1.2723E+02  6.5059E-03  1.4434E+00  0.0000E+00
             2.9307E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -698.677566533079        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2227
 NPARAMETR:  1.1324E+00  2.0889E-01  5.0016E+01  2.0691E+00  2.3744E+00  1.7196E+00  8.6769E+00  9.6337E-01  1.3279E+00  1.0000E-02
             1.0732E+01
 PARAMETER:  2.2434E-01 -1.4660E+00  4.0123E+00  8.2710E-01  9.6473E-01  6.4208E-01  2.2607E+00  6.2678E-02  3.8362E-01 -3.0937E+01
             2.4732E+00
 GRADIENT:   4.6483E-02  4.1888E-01  1.7765E-03 -4.0205E-01 -1.3754E-01  1.0947E-01  1.3787E+01  6.4810E-03 -1.7129E-01  0.0000E+00
            -3.6769E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -698.682747061961        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     2396
 NPARAMETR:  1.1325E+00  2.0510E-01  4.9344E+01  2.0715E+00  2.3752E+00  1.7196E+00  8.6847E+00  9.1271E-01  1.3302E+00  1.0000E-02
             1.0734E+01
 PARAMETER:  2.2445E-01 -1.4843E+00  3.9988E+00  8.2830E-01  9.6510E-01  6.4207E-01  2.2616E+00  8.6594E-03  3.8532E-01 -3.0937E+01
             2.4734E+00
 GRADIENT:   1.6562E-02  1.5480E-01  2.2622E-04 -4.5285E-01 -4.4660E-02  1.8716E-01  1.3202E+01  6.0356E-03 -9.1472E-02  0.0000E+00
            -4.9860E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -698.684964271286        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2581
 NPARAMETR:  1.1330E+00  2.0282E-01  4.9500E+01  2.0749E+00  2.3763E+00  1.7194E+00  8.7126E+00  8.9479E-01  1.3329E+00  1.0000E-02
             1.0738E+01
 PARAMETER:  2.2486E-01 -1.4955E+00  4.0020E+00  8.2993E-01  9.6555E-01  6.4200E-01  2.2648E+00 -1.1169E-02  3.8735E-01 -3.0937E+01
             2.4738E+00
 GRADIENT:   1.2820E-01  1.4087E-01  2.6797E-04 -3.9720E-01 -3.3321E-02  2.0317E-01  1.3381E+01  5.7962E-03 -3.9777E-02  0.0000E+00
            -3.3600E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -698.686513282294        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2771             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1329E+00  1.9959E-01  4.9134E+01  2.0771E+00  2.3775E+00  1.7195E+00  8.7367E+00  8.8592E-01  1.3350E+00  1.0000E-02
             1.0740E+01
 PARAMETER:  2.2481E-01 -1.5115E+00  3.9946E+00  8.3099E-01  9.6605E-01  6.4204E-01  2.2675E+00 -2.1132E-02  3.8895E-01 -3.0937E+01
             2.4739E+00
 GRADIENT:   1.0246E+01  3.7665E+00  2.6893E-03  2.2908E+01  8.7555E-01  9.7700E+00  1.2975E+02  5.8135E-03  1.6972E+00  0.0000E+00
             2.9428E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -698.687182449264        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2957
 NPARAMETR:  1.1330E+00  1.9886E-01  4.9049E+01  2.0784E+00  2.3776E+00  1.7193E+00  8.7494E+00  8.6493E-01  1.3359E+00  1.0000E-02
             1.0741E+01
 PARAMETER:  2.2487E-01 -1.5152E+00  3.9928E+00  8.3158E-01  9.6608E-01  6.4194E-01  2.2690E+00 -4.5101E-02  3.8961E-01 -3.0937E+01
             2.4740E+00
 GRADIENT:   5.3661E-02  3.6569E-02 -1.9788E-03 -6.7723E-01  7.1810E-02  2.5580E-01  1.3497E+01  5.5623E-03  3.3089E-02  0.0000E+00
            -2.2727E-01

0ITERATION NO.:   98    OBJECTIVE VALUE:  -698.687422654858        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:     3059
 NPARAMETR:  1.1331E+00  1.9780E-01  4.8654E+01  2.0793E+00  2.3786E+00  1.7192E+00  8.7611E+00  8.6285E-01  1.3373E+00  1.0000E-02
             1.0742E+01
 PARAMETER:  2.2491E-01 -1.5151E+00  4.0034E+00  8.3215E-01  9.6538E-01  6.4186E-01  2.2699E+00 -4.6516E-02  3.8939E-01 -3.0937E+01
             2.4740E+00
 GRADIENT:  -8.7907E-03  5.5509E-02  2.2649E-03  4.9284E-02 -8.1733E-02 -1.0217E-03 -5.5897E-02  1.9111E-03 -9.0826E-02  0.0000E+00
            -1.8616E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3059
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8072E-02  6.9956E-02  2.1071E-04 -8.1948E-02 -5.6855E-06
 SE:             2.8277E-02  1.9558E-02  9.7790E-05  1.7285E-02  4.5499E-05
 N:                     100         100         100         100         100

 P VAL.:         5.2275E-01  3.4788E-04  3.1181E-02  2.1305E-06  9.0056E-01

 ETASHRINKSD(%)  5.2678E+00  3.4478E+01  9.9672E+01  4.2092E+01  9.9848E+01
 ETASHRINKVR(%)  1.0258E+01  5.7069E+01  9.9999E+01  6.6466E+01  1.0000E+02
 EBVSHRINKSD(%)  8.7708E+00  4.3329E+01  9.9365E+01  3.0143E+01  9.9769E+01
 EBVSHRINKVR(%)  1.6772E+01  6.7884E+01  9.9996E+01  5.1200E+01  9.9999E+01
 RELATIVEINF(%)  8.0623E+01  2.0636E+01  5.6840E-04  2.5327E+01  6.9032E-05
 EPSSHRINKSD(%)  7.9611E+00
 EPSSHRINKVR(%)  1.5288E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -698.68742265485810     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       404.03881719074900     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    82.45
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    13.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -698.687       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.13E+00  1.99E-01  4.96E+01  2.08E+00  2.38E+00  1.72E+00  8.76E+00  8.64E-01  1.34E+00  1.00E-02  1.07E+01
 


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
+        5.05E+02
 
 TH 2
+       -1.86E+02  3.16E+02
 
 TH 3
+       -8.31E-03 -7.40E-04  4.06E-07
 
 TH 4
+       -2.04E+02  1.32E+02  5.63E-03  1.47E+02
 
 TH 5
+        2.07E+01 -1.34E+01 -6.71E-04 -1.62E+01  1.84E+00
 
 TH 6
+       -5.22E+01 -2.87E+01  2.70E-03  1.53E+01 -2.62E+00  5.66E+01
 
 TH 7
+       -1.06E+00  1.39E+01 -3.72E-04  5.32E-01  2.03E-02 -2.78E+00  9.11E-01
 
 TH 8
+       -1.17E-01  7.56E-02  2.10E-06  5.83E-02 -6.53E-03  2.84E-02  1.81E-03  4.37E-05
 
 TH 9
+        5.99E+01 -3.06E+01 -2.11E-03 -4.73E+01  5.30E+00 -4.39E+00  6.36E-01 -1.54E-02  1.62E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.18E+00  4.63E-01 -2.01E-04 -4.37E+00  5.40E-01  1.91E+00  2.78E-01  1.69E-03  1.89E+00  0.00E+00  1.07E+00
 
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
+        2.37E+02
 
 TH 2
+       -2.98E+00  2.01E+02
 
 TH 3
+       -1.05E-03 -2.37E-03  3.70E-05
 
 TH 4
+       -1.10E+01  3.71E+01  2.88E-03  7.91E+01
 
 TH 5
+       -9.80E-01 -4.82E+00 -1.25E-02 -9.25E+00  9.84E+00
 
 TH 6
+       -3.07E-02 -9.69E+00  5.81E-04  1.70E+00 -7.18E-01  4.97E+01
 
 TH 7
+        4.76E-01  1.14E+01 -2.66E-04 -2.24E+00  2.08E-01 -4.59E-01  1.26E+00
 
 TH 8
+       -1.32E-02  2.97E-02 -6.83E-05 -5.53E-03  1.54E-02  4.61E-02  1.69E-03  9.28E-01
 
 TH 9
+       -7.44E-01  9.06E-01 -1.34E-03 -2.66E+01  2.43E+00 -2.37E+00  4.15E-01 -1.08E-01  3.05E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.42E+00  7.06E-01  3.37E-04 -7.56E+00  6.14E-01  1.36E+00  3.11E-01  2.08E-03  2.00E+00  0.00E+00  5.86E+00
 
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
+        2.47E+02
 
 TH 2
+        8.43E+01  1.83E+02
 
 TH 3
+        2.42E-03  2.68E-03  8.87E-06
 
 TH 4
+        9.43E+01  1.42E+01  1.89E-03  1.01E+02
 
 TH 5
+       -1.26E+01  5.55E+00 -4.74E-03 -1.02E+01  4.89E+00
 
 TH 6
+        3.47E+01  3.10E+01  7.98E-04 -4.17E-02 -1.86E+00  5.26E+01
 
 TH 7
+       -8.86E-01  1.16E+01 -1.71E-04 -4.00E+00  1.30E+00  7.91E-01  1.20E+00
 
 TH 8
+        4.60E-04  2.12E-03 -3.03E-07  1.16E-04  2.16E-04  2.81E-04  7.36E-05  5.13E-07
 
 TH 9
+       -2.71E+00 -6.25E+00  7.30E-04 -8.75E+00  4.89E-01  5.46E+00 -4.98E-01  5.03E-04  2.57E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.13E+01 -1.34E+01 -1.89E-03 -1.03E+01  2.57E+00 -6.48E+00  3.33E-01 -1.28E-03 -1.27E+00  0.00E+00  7.42E+01
 
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
 #CPUT: Total CPU Time in Seconds,       95.669
Stop Time:
Thu Sep 30 09:19:09 CDT 2021
