Sat Sep 25 00:53:40 CDT 2021
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
$DATA ../../../../data/int/SL2/dat11.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      997
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

 TOT. NO. OF OBS RECS:      897
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
 RAW OUTPUT FILE (FILE): m11.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1458.95145436801        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8962E+02 -3.6938E+01  7.8809E+01  9.8770E+01  1.7727E+02  2.0808E+01 -8.3374E+01 -1.8833E+02 -1.5022E+02 -2.1359E+00
            -4.3152E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2863.11955538822        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.3291E-01  1.2365E+00  1.0183E+00  8.4612E-01  1.0675E+00  8.5960E-01  1.0264E+00  9.8768E-01  1.2612E+00  8.3456E-01
             2.1384E+00
 PARAMETER:  3.0552E-02  3.1228E-01  1.1818E-01 -6.7097E-02  1.6534E-01 -5.1283E-02  1.2602E-01  8.7600E-02  3.3206E-01 -8.0851E-02
             8.6008E-01
 GRADIENT:  -2.4362E+01 -1.5067E+01 -8.0144E+00 -1.7616E+01  1.0317E+01 -2.8304E+01 -4.8908E-01 -2.5431E+00 -7.1139E+00 -1.1440E+01
            -1.7872E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2867.92631259731        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.4412E-01  1.4966E+00  1.1724E+00  7.4918E-01  1.2671E+00  8.9397E-01  8.1687E-01  9.0179E-01  1.4631E+00  1.1075E+00
             2.1550E+00
 PARAMETER:  4.2495E-02  5.0319E-01  2.5904E-01 -1.8878E-01  3.3670E-01 -1.2084E-02 -1.0228E-01 -3.3786E-03  4.8054E-01  2.0207E-01
             8.6779E-01
 GRADIENT:   8.7565E+00  4.6698E+01  7.4934E+00  3.4109E+01 -2.3865E+00 -1.1462E+01 -3.5450E+00 -6.0110E+00  3.6104E+00 -1.4327E+00
            -1.6696E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2876.04248037364        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.4286E-01  1.4977E+00  1.2628E+00  7.2507E-01  1.2994E+00  9.2339E-01  8.4240E-01  1.7036E+00  1.3977E+00  1.1131E+00
             2.2979E+00
 PARAMETER:  4.1159E-02  5.0393E-01  3.3332E-01 -2.2148E-01  3.6190E-01  2.0295E-02 -7.1502E-02  6.3273E-01  4.3481E-01  2.0718E-01
             9.3199E-01
 GRADIENT:   6.1297E-01  4.8902E+00 -3.5191E+00  7.3046E+00  5.9370E-01  1.4649E+00  1.7854E+00  2.0551E+00 -1.8010E-01  1.6488E+00
             1.0050E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2877.49285616467        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  9.4653E-01  1.6827E+00  1.5526E+00  6.0620E-01  1.4890E+00  9.2215E-01  7.9991E-01  2.3478E+00  1.5608E+00  1.2067E+00
             2.2838E+00
 PARAMETER:  4.5048E-02  6.2043E-01  5.3992E-01 -4.0055E-01  4.9808E-01  1.8957E-02 -1.2326E-01  9.5346E-01  5.4518E-01  2.8786E-01
             9.2582E-01
 GRADIENT:   3.2701E+00 -8.6199E+00 -3.1638E+00 -1.3510E+00  4.6927E+00 -1.0431E-01  1.5549E+00  8.8807E-01  1.8743E+00 -2.2254E+00
             1.0435E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2877.94913195633        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  9.4461E-01  1.7363E+00  1.8004E+00  5.7228E-01  1.5566E+00  9.2086E-01  7.8883E-01  2.8827E+00  1.5690E+00  1.2655E+00
             2.2760E+00
 PARAMETER:  4.3021E-02  6.5176E-01  6.8799E-01 -4.5812E-01  5.4247E-01  1.7553E-02 -1.3721E-01  1.1587E+00  5.5046E-01  3.3551E-01
             9.2243E-01
 GRADIENT:  -1.6659E+00 -8.2965E+00 -3.3033E+00 -6.4537E+00  2.3279E+00 -6.8577E-01  3.2839E-01  3.6488E-01  1.9346E-01  1.7745E-02
             1.6750E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2877.95418999602        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      696
 NPARAMETR:  9.4496E-01  1.7366E+00  1.8023E+00  5.7233E-01  1.5569E+00  9.2302E-01  7.8529E-01  2.8826E+00  1.5635E+00  1.2657E+00
             2.2771E+00
 PARAMETER:  4.3390E-02  6.5195E-01  6.8904E-01 -4.5805E-01  5.4270E-01  1.9893E-02 -1.4171E-01  1.1587E+00  5.4692E-01  3.3560E-01
             9.2288E-01
 GRADIENT:  -7.7277E-01 -7.7683E+00 -3.2321E+00 -6.5772E+00  2.3433E+00  2.0766E-01 -3.8258E-01  2.5630E-01 -5.0560E-01 -4.6793E-02
             9.5904E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2877.98030872018        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  9.4531E-01  1.7444E+00  1.8060E+00  5.7287E-01  1.5561E+00  9.2260E-01  7.8757E-01  2.8812E+00  1.5686E+00  1.2657E+00
             2.2756E+00
 PARAMETER:  4.3760E-02  6.5640E-01  6.9111E-01 -4.5709E-01  5.4220E-01  1.9442E-02 -1.3881E-01  1.1582E+00  5.5021E-01  3.3559E-01
             9.2226E-01
 GRADIENT:   7.6550E-02  1.7507E+00 -2.6470E+00 -2.9356E+00 -3.3032E-02 -1.6471E-02 -5.4671E-03  6.1694E-03  8.6458E-03 -2.1558E-01
            -8.6029E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2878.35479167397        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1034             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4373E-01  1.7163E+00  2.2105E+00  5.9419E-01  1.5596E+00  9.2145E-01  7.8595E-01  2.9607E+00  1.5341E+00  1.2723E+00
             2.2798E+00
 PARAMETER:  4.2082E-02  6.4020E-01  8.9323E-01 -4.2056E-01  5.4442E-01  1.8188E-02 -1.4087E-01  1.1854E+00  5.2792E-01  3.4085E-01
             9.2410E-01
 GRADIENT:   4.6917E+00  2.7280E+01 -7.0965E+00 -7.9413E-01 -2.3762E+01  7.3816E-01 -5.1800E-01  1.0051E+01  3.8741E+00 -6.9491E+00
            -2.3167E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2878.40574776364        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1215
 NPARAMETR:  9.4392E-01  1.7141E+00  2.2066E+00  5.9469E-01  1.5614E+00  9.2153E-01  7.8602E-01  2.9502E+00  1.5357E+00  1.2715E+00
             2.2841E+00
 PARAMETER:  4.2284E-02  6.3887E-01  8.9146E-01 -4.1971E-01  5.4560E-01  1.8280E-02 -1.4077E-01  1.1819E+00  5.2898E-01  3.4017E-01
             9.2598E-01
 GRADIENT:   4.7744E+00  2.5845E+01 -4.5824E+00 -2.3084E+00 -1.7941E+01  7.3792E-01 -3.4401E-01  5.3717E+00  3.2025E+00 -5.5875E+00
             5.9451E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2878.40774120257        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1286
 NPARAMETR:  9.4393E-01  1.7136E+00  2.2065E+00  5.9472E-01  1.5617E+00  9.1757E-01  7.8601E-01  2.9494E+00  1.5358E+00  1.2714E+00
             2.2841E+00
 PARAMETER:  4.2293E-02  6.3857E-01  8.9140E-01 -4.1966E-01  5.4578E-01  1.3972E-02 -1.4078E-01  1.1816E+00  5.2903E-01  3.4014E-01
             9.2596E-01
 GRADIENT:   4.7175E+00  2.5434E+01 -4.3228E+00 -2.5716E+00 -1.7212E+01 -8.9775E-01 -3.2549E-01  4.9215E+00  3.1310E+00 -5.4217E+00
             6.2462E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2878.41736968744        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1356
 NPARAMETR:  9.4399E-01  1.7100E+00  2.2057E+00  5.9494E-01  1.5636E+00  9.1005E-01  7.8593E-01  2.9436E+00  1.5362E+00  1.2712E+00
             2.2838E+00
 PARAMETER:  4.2356E-02  6.3649E-01  8.9106E-01 -4.1929E-01  5.4699E-01  5.7491E-03 -1.4088E-01  1.1796E+00  5.2933E-01  3.3995E-01
             9.2582E-01
 GRADIENT:   4.6530E+00  2.2164E+01 -2.9044E+00 -4.3594E+00 -1.2723E+01 -4.0651E+00 -2.1903E-01  2.4790E+00  2.7409E+00 -4.4336E+00
             7.8685E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2878.43250173121        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  9.4414E-01  1.7009E+00  2.2038E+00  5.9550E-01  1.5684E+00  9.0690E-01  7.8574E-01  2.9289E+00  1.5374E+00  1.2706E+00
             2.2828E+00
 PARAMETER:  4.2515E-02  6.3116E-01  8.9020E-01 -4.1835E-01  5.5006E-01  2.2779E-03 -1.4113E-01  1.1746E+00  5.3009E-01  3.3948E-01
             9.2539E-01
 GRADIENT:   4.9344E+00  1.2983E+01 -9.7287E-01 -8.3305E+00 -3.6277E+00 -5.4294E+00 -2.3521E-02 -1.0804E+00  2.2584E+00 -2.8350E+00
             9.8043E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2878.43666940146        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1496
 NPARAMETR:  9.4424E-01  1.6946E+00  2.2026E+00  5.9589E-01  1.5717E+00  9.1132E-01  7.8561E-01  2.9190E+00  1.5382E+00  1.2702E+00
             2.2819E+00
 PARAMETER:  4.2623E-02  6.2747E-01  8.8962E-01 -4.1770E-01  5.5217E-01  7.1347E-03 -1.4130E-01  1.1712E+00  5.3060E-01  3.3915E-01
             9.2500E-01
 GRADIENT:   5.1867E+00  6.2875E+00 -4.8332E-01 -1.0788E+01  1.4862E+00 -3.5763E+00  5.7787E-02 -2.3467E+00  2.1624E+00 -2.2414E+00
             1.0137E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2878.43791358094        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1566
 NPARAMETR:  9.4427E-01  1.6926E+00  2.2022E+00  5.9602E-01  1.5728E+00  9.1512E-01  7.8556E-01  2.9158E+00  1.5384E+00  1.2700E+00
             2.2814E+00
 PARAMETER:  4.2656E-02  6.2626E-01  8.8945E-01 -4.1748E-01  5.5285E-01  1.1295E-02 -1.4135E-01  1.1702E+00  5.3076E-01  3.3906E-01
             9.2480E-01
 GRADIENT:   5.4698E+00  4.0479E+00 -4.2304E-01 -1.1583E+01  3.0044E+00 -1.9529E+00  7.3830E-02 -2.6301E+00  2.1366E+00 -2.1090E+00
             1.0004E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2878.43813348310        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1636
 NPARAMETR:  9.4428E-01  1.6918E+00  2.2021E+00  5.9608E-01  1.5732E+00  9.1721E-01  7.8555E-01  2.9148E+00  1.5385E+00  1.2700E+00
             2.2811E+00
 PARAMETER:  4.2667E-02  6.2581E-01  8.8939E-01 -4.1738E-01  5.5309E-01  1.3579E-02 -1.4137E-01  1.1698E+00  5.3081E-01  3.3902E-01
             9.2466E-01
 GRADIENT:   5.5370E+00  3.2223E+00 -4.2530E-01 -1.1870E+01  3.5185E+00 -1.0821E+00  7.7086E-02 -2.7057E+00  2.1335E+00 -2.0753E+00
             9.7975E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2878.43816016319        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1709
 NPARAMETR:  9.4428E-01  1.6915E+00  2.2020E+00  5.9611E-01  1.5733E+00  9.1832E-01  7.8554E-01  2.9144E+00  1.5385E+00  1.2700E+00
             2.2808E+00
 PARAMETER:  4.2671E-02  6.2560E-01  8.8938E-01 -4.1734E-01  5.5319E-01  1.4790E-02 -1.4138E-01  1.1697E+00  5.3083E-01  3.3901E-01
             9.2455E-01
 GRADIENT:   5.5811E+00  2.8452E+00 -4.3956E-01 -1.1996E+01  3.7333E+00 -6.2063E-01  7.7468E-02 -2.7251E+00  2.1348E+00 -2.0660E+00
             9.5940E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2878.43818849108        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1779
 NPARAMETR:  9.4429E-01  1.6911E+00  2.2020E+00  5.9614E-01  1.5735E+00  9.1951E-01  7.8554E-01  2.9141E+00  1.5385E+00  1.2700E+00
             2.2804E+00
 PARAMETER:  4.2673E-02  6.2537E-01  8.8938E-01 -4.1728E-01  5.5328E-01  1.6084E-02 -1.4138E-01  1.1695E+00  5.3084E-01  3.3901E-01
             9.2437E-01
 GRADIENT:   5.6263E+00  2.4577E+00 -4.7194E-01 -1.2124E+01  3.9228E+00 -1.3000E-01  7.6421E-02 -2.7255E+00  2.1376E+00 -2.0674E+00
             9.2230E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2878.43822514366        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1849
 NPARAMETR:  9.4429E-01  1.6907E+00  2.2020E+00  5.9618E-01  1.5736E+00  9.2054E-01  7.8554E-01  2.9139E+00  1.5385E+00  1.2700E+00
             2.2799E+00
 PARAMETER:  4.2674E-02  6.2516E-01  8.8939E-01 -4.1721E-01  5.5336E-01  1.7209E-02 -1.4139E-01  1.1695E+00  5.3084E-01  3.3901E-01
             9.2414E-01
 GRADIENT:   5.6626E+00  2.1140E+00 -5.1994E-01 -1.2233E+01  4.0603E+00  2.9509E-01  7.4003E-02 -2.7042E+00  2.1436E+00 -2.0796E+00
             8.7193E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2878.48963290667        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     1966
 NPARAMETR:  9.4415E-01  1.6990E+00  2.2038E+00  5.9569E-01  1.5692E+00  9.2516E-01  7.8572E-01  2.9272E+00  1.5375E+00  1.2705E+00
             2.2806E+00
 PARAMETER:  4.2527E-02  6.3003E-01  8.9020E-01 -4.1804E-01  5.5057E-01  2.2212E-02 -1.4116E-01  1.1741E+00  5.3013E-01  3.3945E-01
             9.2443E-01
 GRADIENT:  -2.4811E+00 -8.3317E+00 -1.4799E+00 -1.1783E+01 -6.5521E+00  1.1533E+00 -2.7319E-01 -1.2501E+00  1.5171E+00 -3.2291E+00
             6.2599E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2878.89625001366        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     2080
 NPARAMETR:  9.4103E-01  1.6512E+00  2.2700E+00  6.2953E-01  1.5696E+00  9.1743E-01  8.0680E-01  2.9828E+00  1.4457E+00  1.2929E+00
             2.2660E+00
 PARAMETER:  3.9215E-02  6.0150E-01  9.1979E-01 -3.6278E-01  5.5081E-01  1.3823E-02 -1.1467E-01  1.1929E+00  4.6860E-01  3.5689E-01
             9.1800E-01
 GRADIENT:  -2.8213E+00  2.8070E+00 -4.9219E+00 -6.8161E+00  8.4632E+00 -1.0762E+00 -2.7979E-01 -1.0248E+00  1.5046E-01  2.7496E+00
            -3.3522E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2878.90922522562        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     2217
 NPARAMETR:  9.4100E-01  1.6495E+00  2.2722E+00  6.3048E-01  1.5694E+00  9.2121E-01  8.0953E-01  2.9875E+00  1.4435E+00  1.2926E+00
             2.2664E+00
 PARAMETER:  3.9186E-02  6.0045E-01  9.2076E-01 -3.6127E-01  5.5067E-01  1.7928E-02 -1.1130E-01  1.1944E+00  4.6710E-01  3.5667E-01
             9.1818E-01
 GRADIENT:  -2.7296E+00  1.9308E+00 -5.6366E+00 -6.7930E+00  7.8603E+00  4.9379E-01 -1.8089E-02 -2.1131E-01  4.3309E-01  2.4385E+00
            -3.1709E+00

0ITERATION NO.:  109    OBJECTIVE VALUE:  -2878.90975659211        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     2352
 NPARAMETR:  9.4100E-01  1.6494E+00  2.2722E+00  6.3049E-01  1.5694E+00  9.2185E-01  8.0997E-01  2.9878E+00  1.4437E+00  1.2926E+00
             2.2664E+00
 PARAMETER:  3.9186E-02  6.0045E-01  9.2077E-01 -3.6127E-01  5.5067E-01  1.8530E-02 -1.1076E-01  1.1946E+00  4.6720E-01  3.5667E-01
             9.1818E-01
 GRADIENT:  -2.3893E+03  7.7517E+02  5.0793E+02 -6.6700E+02 -8.6267E+02 -3.0496E-01  4.2944E+03  3.9986E+02 -5.0902E+02  6.6901E+02
            -2.6550E+02
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         1.8         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2352
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2581E-03 -3.3071E-02 -2.4366E-02  3.5382E-02 -3.7772E-02
 SE:             2.9523E-02  2.0603E-02  1.5473E-02  2.3338E-02  2.3712E-02
 N:                     100         100         100         100         100

 P VAL.:         8.3213E-01  1.0845E-01  1.1531E-01  1.2949E-01  1.1117E-01

 ETASHRINKSD(%)  1.0951E+00  3.0977E+01  4.8163E+01  2.1816E+01  2.0562E+01
 ETASHRINKVR(%)  2.1782E+00  5.2359E+01  7.3130E+01  3.8872E+01  3.6897E+01
 EBVSHRINKSD(%)  1.4565E+00  3.0208E+01  5.4265E+01  2.6288E+01  1.5063E+01
 EBVSHRINKVR(%)  2.8918E+00  5.1291E+01  7.9083E+01  4.5665E+01  2.7858E+01
 RELATIVEINF(%)  9.7074E+01  7.6672E+00  9.7832E+00  8.8959E+00  3.6712E+01
 EPSSHRINKSD(%)  1.7602E+01
 EPSSHRINKVR(%)  3.2106E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2878.9097565921084     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1230.3340280229256     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    59.26
 Elapsed covariance  time in seconds:    14.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2878.910       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.41E-01  1.65E+00  2.27E+00  6.30E-01  1.57E+00  9.22E-01  8.10E-01  2.99E+00  1.44E+00  1.29E+00  2.27E+00
 


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
 
         2.17E-02  2.49E-03  5.17E-03  5.89E-04  2.27E-03  6.84E-02  1.80E-04  2.14E-01  4.17E-02  1.06E-01  5.71E-03
 


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
+        4.72E-04
 
 TH 2
+       -3.08E-05  6.19E-06
 
 TH 3
+       -6.42E-05  1.28E-05  2.68E-05
 
 TH 4
+        7.07E-06 -1.46E-06 -3.04E-06  3.47E-07
 
 TH 5
+        2.84E-05 -5.65E-06 -1.18E-05  1.34E-06  5.17E-06
 
 TH 6
+        2.22E-04 -4.54E-06 -7.54E-06  6.96E-07  2.63E-06  4.68E-03
 
 TH 7
+       -1.68E-06  4.40E-07  9.17E-07 -1.05E-07 -4.02E-07  5.44E-08  3.24E-08
 
 TH 8
+        3.91E-03 -4.75E-04 -9.87E-04  1.11E-04  4.35E-04  1.35E-03 -3.09E-05  4.58E-02
 
 TH 9
+        7.56E-04 -9.36E-05 -1.94E-04  2.19E-05  8.57E-05  2.56E-04 -6.11E-06  8.93E-03  1.74E-03
 
 TH10
+        2.30E-03 -1.48E-04 -3.08E-04  3.40E-05  1.36E-04  1.09E-03 -8.04E-06  1.89E-02  3.66E-03  1.12E-02
 
 TH11
+        6.31E-05 -1.36E-05 -2.82E-05  3.22E-06  1.24E-05  1.81E-05 -9.70E-07  1.05E-03  2.07E-04  3.03E-04  3.26E-05
 
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
+        2.17E-02
 
 TH 2
+       -5.69E-01  2.49E-03
 
 TH 3
+       -5.72E-01  9.98E-01  5.17E-03
 
 TH 4
+        5.53E-01 -9.97E-01 -9.99E-01  5.89E-04
 
 TH 5
+        5.75E-01 -9.98E-01 -1.00E+00  9.99E-01  2.27E-03
 
 TH 6
+        1.49E-01 -2.66E-02 -2.13E-02  1.73E-02  1.69E-02  6.84E-02
 
 TH 7
+       -4.31E-01  9.83E-01  9.85E-01 -9.88E-01 -9.83E-01  4.42E-03  1.80E-04
 
 TH 8
+        8.42E-01 -8.93E-01 -8.91E-01  8.84E-01  8.95E-01  9.21E-02 -8.02E-01  2.14E-01
 
 TH 9
+        8.34E-01 -9.01E-01 -9.00E-01  8.92E-01  9.03E-01  8.95E-02 -8.14E-01  1.00E+00  4.17E-02
 
 TH10
+        1.00E+00 -5.62E-01 -5.64E-01  5.45E-01  5.68E-01  1.50E-01 -4.23E-01  8.37E-01  8.29E-01  1.06E-01
 
 TH11
+        5.09E-01 -9.56E-01 -9.53E-01  9.57E-01  9.55E-01  4.63E-02 -9.44E-01  8.59E-01  8.67E-01  5.02E-01  5.71E-03
 
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
+        3.91E+11
 
 TH 2
+       -6.44E+09  7.53E+08
 
 TH 3
+        5.77E+09 -7.41E+08  1.15E+09
 
 TH 4
+        7.76E+10 -9.48E+09  9.69E+09  1.29E+11
 
 TH 5
+       -2.36E+09  1.82E+08  7.58E+07 -2.44E+09  3.70E+08
 
 TH 6
+       -3.37E+06  8.45E+04 -5.33E+04 -9.28E+05  7.35E+04  2.64E+02
 
 TH 7
+        3.46E+11 -4.44E+09 -5.86E+09  5.67E+10 -7.12E+09 -2.79E+06  6.97E+11
 
 TH 8
+        4.36E+09 -4.00E+08  4.79E+08  5.29E+09 -4.30E+07 -5.29E+04 -2.01E+09  3.29E+08
 
 TH 9
+       -2.37E+10  2.13E+09 -2.54E+09 -2.82E+10  2.39E+08  2.86E+05  9.70E+09 -1.75E+09  9.28E+09
 
 TH10
+       -7.97E+10  1.31E+09 -1.18E+09 -1.58E+10  4.82E+08  6.89E+05 -7.07E+10 -8.91E+08  4.84E+09  1.63E+10
 
 TH11
+       -3.59E+07  1.24E+07 -2.70E+07 -1.60E+08 -4.72E+06 -7.11E+02  3.62E+08 -9.02E+06  4.75E+07  7.31E+06  1.23E+06
 
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
 #CPUT: Total CPU Time in Seconds,       73.918
Stop Time:
Sat Sep 25 00:54:55 CDT 2021
