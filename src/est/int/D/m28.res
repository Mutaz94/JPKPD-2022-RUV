Wed Sep 29 08:25:03 CDT 2021
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
$DATA ../../../../data/int/D/dat28.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20095.7209661675        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9683E+02  1.3872E+02 -1.2899E+01 -1.3634E+02  1.7074E+02 -1.0319E+03 -5.6603E+02 -1.1192E+02 -1.0637E+03 -2.7326E+02
            -4.4599E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1016.12697095164        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0485E+00  1.9939E+00  9.0917E-01  3.1263E+00  8.4335E-01  4.2966E+00  3.0706E+00  1.0299E+00  2.3907E+00  1.5391E+00
             1.3057E+01
 PARAMETER:  1.4737E-01  7.9009E-01  4.7812E-03  1.2398E+00 -7.0378E-02  1.5578E+00  1.2219E+00  1.2947E-01  9.7158E-01  5.3117E-01
             2.6693E+00
 GRADIENT:  -3.6596E+01  6.2946E+01 -2.5428E+01  2.0365E+02 -1.3763E+01  1.5493E+02 -1.9393E+01  5.2825E+00 -4.2274E+01  2.8650E+01
             6.0683E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1129.59902449705        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.5847E-01  2.1321E+00  3.7454E+01  3.1266E+00  2.5537E+00  4.4218E+00  7.3495E+00  1.1918E+00  2.8249E+00  2.2758E+00
             1.2195E+01
 PARAMETER:  5.7587E-02  8.5709E-01  3.7231E+00  1.2399E+00  1.0376E+00  1.5866E+00  2.0946E+00  2.7543E-01  1.1385E+00  9.2231E-01
             2.6010E+00
 GRADIENT:  -3.9324E+01  4.3199E+01 -1.6852E+00  9.9825E+01 -3.4321E-01  1.7117E+02  4.2864E+01  3.4037E-02  2.7943E+01  5.9133E+01
             6.4448E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1414.13623468546        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.2934E+00  4.1347E-01  1.2514E+02  1.8249E+00  3.3198E+00  2.0863E+00  4.1659E+00  6.7406E+00  2.2035E+00  6.5915E-01
             8.6039E+00
 PARAMETER:  3.5730E-01 -7.8317E-01  4.9295E+00  7.0152E-01  1.2999E+00  8.3538E-01  1.5269E+00  2.0081E+00  8.9004E-01 -3.1680E-01
             2.2522E+00
 GRADIENT:   3.0116E+01 -1.1227E+01 -3.8001E+00 -6.4533E+01  9.7553E+01 -1.0238E+01  5.9086E+00  2.7033E+00 -2.0312E+00  3.3987E+00
             4.0729E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1481.17729215016        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.1452E+00  1.9993E-01  1.0959E+02  1.9616E+00  2.3396E+00  2.0788E+00  3.5252E+00  4.2490E+00  2.0411E+00  6.4448E-01
             6.7761E+00
 PARAMETER:  2.3559E-01 -1.5098E+00  4.7968E+00  7.7374E-01  9.4998E-01  8.3178E-01  1.3599E+00  1.5467E+00  8.1348E-01 -3.3931E-01
             2.0134E+00
 GRADIENT:   2.1944E+00 -2.0265E+00  3.7675E-01 -3.2454E+00  3.7821E+00 -3.2109E+00  2.0261E+00 -1.9426E+00 -2.0069E+00 -3.5408E+00
             1.7176E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1488.19165468747        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      436
 NPARAMETR:  1.1969E+00  2.0575E-01  9.8916E+01  2.2050E+00  2.3622E+00  2.2749E+00  3.2676E+00  4.4549E+00  2.1240E+00  8.9413E-01
             6.8442E+00
 PARAMETER:  2.7975E-01 -1.4811E+00  4.6943E+00  8.9071E-01  9.5957E-01  9.2194E-01  1.2841E+00  1.5940E+00  8.5329E-01 -1.1908E-02
             2.0234E+00
 GRADIENT:  -1.3145E+00 -5.2395E-01 -1.0680E+00 -1.9025E+00  5.4874E+00 -3.1853E+00  1.3656E+00 -1.4661E-01 -4.3259E+00  5.1241E-01
            -1.0881E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1489.10004130403        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      613
 NPARAMETR:  1.2076E+00  1.5299E-01  5.3730E+02  2.2198E+00  2.4328E+00  2.2912E+00  2.5466E+00  9.1248E+00  2.1311E+00  8.6618E-01
             6.8548E+00
 PARAMETER:  2.8863E-01 -1.7774E+00  6.3865E+00  8.9742E-01  9.8906E-01  9.2908E-01  1.0348E+00  2.3110E+00  8.5664E-01 -4.3658E-02
             2.0249E+00
 GRADIENT:   2.0325E+00 -1.4859E+00  6.2190E-01 -8.2693E+00  7.0696E+00 -1.2236E+00  7.3366E-01 -2.1468E+00  6.5544E-01  1.3111E+00
             5.3073E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1491.39323689542        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  1.1826E+00  4.5878E-01  7.2985E+02  2.0974E+00  2.4534E+00  2.2321E+00  2.9774E-01  1.0091E+01  2.1845E+00  8.0485E-01
             6.9086E+00
 PARAMETER:  2.6774E-01 -6.7919E-01  6.6928E+00  8.4071E-01  9.9748E-01  9.0294E-01 -1.1115E+00  2.4116E+00  8.8141E-01 -1.1710E-01
             2.0328E+00
 GRADIENT:   1.6032E+01  5.3367E+00  8.7770E-01  6.5361E+01  8.9752E+00  2.4810E+01  1.4868E-01 -3.0912E-01 -9.0019E+00 -2.9427E-01
             4.9218E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1494.23164669755        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      808
 NPARAMETR:  1.1465E+00  8.5831E-01  9.3858E+02  1.6012E+00  2.5287E+00  2.1102E+00  1.0000E-02  9.1793E+00  2.5513E+00  7.6366E-01
             6.8091E+00
 PARAMETER:  2.3672E-01 -5.2794E-02  6.9444E+00  5.7074E-01  1.0277E+00  8.4677E-01 -4.5446E+00  2.3170E+00  1.0366E+00 -1.6963E-01
             2.0183E+00
 GRADIENT:   5.6249E-01  3.8002E+00 -1.4525E-01  1.5093E+01  2.1681E+01  1.5051E+00  7.5280E-05 -3.9404E-03 -1.0104E+01 -2.6581E+00
             1.9010E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1495.63001360316        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  1.1430E+00  9.9353E-01  1.0479E+03  1.4246E+00  2.5287E+00  2.1084E+00  1.0000E-02  8.8115E+00  2.7942E+00  9.6215E-01
             6.7053E+00
 PARAMETER:  2.3368E-01  9.3509E-02  7.0546E+00  4.5389E-01  1.0277E+00  8.4592E-01 -6.0360E+00  2.2761E+00  1.1276E+00  6.1413E-02
             2.0029E+00
 GRADIENT:   5.0159E-01  6.9182E-01 -1.3584E-01  1.7487E+00  2.9233E+01  1.0883E+00  0.0000E+00  1.2086E-02  3.8896E-01  1.2970E-01
             1.1224E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1499.39481557296        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1001
 NPARAMETR:  1.2104E+00  1.0395E+00  9.9901E+02  1.4215E+00  2.4595E+00  2.2990E+00  1.0000E-02  8.9857E+00  2.9473E+00  9.0613E-01
             6.7494E+00
 PARAMETER:  2.9094E-01  1.3871E-01  7.0068E+00  4.5172E-01  9.9998E-01  9.3248E-01 -5.8337E+00  2.2956E+00  1.1809E+00  1.4240E-03
             2.0095E+00
 GRADIENT:   4.4623E+00  3.1931E+00 -1.0954E-01 -5.7082E+00  4.8670E+00 -2.5930E-01  0.0000E+00  7.5711E-03 -2.1594E+01 -1.7531E+00
            -2.9961E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1508.21736377202        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1180
 NPARAMETR:  1.1805E+00  1.1889E+00  1.0503E+03  1.3373E+00  2.4436E+00  2.2285E+00  1.0000E-02  8.8239E+00  3.5503E+00  9.6417E-01
             6.8463E+00
 PARAMETER:  2.6592E-01  2.7299E-01  7.0568E+00  3.9067E-01  9.9349E-01  9.0132E-01 -7.0189E+00  2.2775E+00  1.3670E+00  6.3517E-02
             2.0237E+00
 GRADIENT:  -2.1320E-01  7.0186E-01 -1.7606E-01  2.3930E+00  6.0230E-02 -1.0981E+00  0.0000E+00  1.1847E-01 -4.9179E-01  1.3842E-02
            -4.1239E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1508.35767747996        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1318
 NPARAMETR:  1.1804E+00  1.2060E+00  2.4765E+03  1.2917E+00  2.4421E+00  2.2370E+00  1.0000E-02  3.2085E+00  3.6055E+00  9.6748E-01
             6.8552E+00
 PARAMETER:  2.6589E-01  2.8728E-01  7.9146E+00  3.5598E-01  9.9284E-01  9.0513E-01 -7.2987E+00  1.2658E+00  1.3825E+00  6.6939E-02
             2.0250E+00
 GRADIENT:  -2.3823E-01 -1.4616E+00 -2.1181E-02 -9.7846E-01 -3.0377E-01  3.2422E-01  0.0000E+00  4.1056E-03  4.1884E-01  1.5396E-01
             2.1415E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1508.37198134991        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1496
 NPARAMETR:  1.1818E+00  1.2011E+00  1.1735E+04  1.3075E+00  2.4493E+00  2.2351E+00  1.0000E-02  1.7360E+00  3.5855E+00  9.5542E-01
             6.8601E+00
 PARAMETER:  2.6705E-01  2.8324E-01  9.4703E+00  3.6815E-01  9.9581E-01  9.0428E-01 -7.2987E+00  6.5157E-01  1.3769E+00  5.4394E-02
             2.0257E+00
 GRADIENT:   1.4440E-01 -4.0293E-01 -3.6323E-03  1.6675E-01  2.1857E-01  5.0086E-02  0.0000E+00  5.4008E-05 -8.7178E-02 -5.4595E-02
             3.3426E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1508.38020883531        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1678             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1819E+00  1.2060E+00  5.6697E+04  1.3010E+00  2.4477E+00  2.2408E+00  1.0000E-02  1.1591E+00  3.6035E+00  9.6060E-01
             6.8578E+00
 PARAMETER:  2.6709E-01  2.8734E-01  1.1045E+01  3.6311E-01  9.9516E-01  9.0683E-01 -7.2987E+00  2.4764E-01  1.3819E+00  5.9798E-02
             2.0254E+00
 GRADIENT:   2.3170E+01  5.1998E+00 -6.5148E-04  1.0468E+01  5.5915E+00  3.6181E+01  0.0000E+00  1.2707E-07  3.7579E+01  8.4401E-02
             3.6581E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1508.38121274279        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1818
 NPARAMETR:  1.1818E+00  1.2076E+00  4.2192E+05  1.2999E+00  2.4478E+00  2.2409E+00  1.0000E-02  1.1540E+00  3.6058E+00  9.5811E-01
             6.8569E+00
 PARAMETER:  2.6706E-01  2.8860E-01  1.3053E+01  3.6229E-01  9.9520E-01  9.0687E-01 -7.2987E+00  2.4327E-01  1.3825E+00  5.7203E-02
             2.0253E+00
 GRADIENT:  -2.9945E+02  5.5207E+01 -1.0182E-04  4.9081E+01 -6.4658E+00  6.4279E+01  0.0000E+00  1.3797E-05 -2.3201E+01  7.7537E-01
             3.6543E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1508.38154997176        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     2021             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1818E+00  1.2083E+00  4.6658E+05  1.2983E+00  2.4482E+00  2.2409E+00  1.0000E-02  1.2014E+00  3.6097E+00  9.5896E-01
             6.8578E+00
 PARAMETER:  2.6707E-01  2.8923E-01  1.3153E+01  3.6104E-01  9.9534E-01  9.0687E-01 -7.2987E+00  2.8351E-01  1.3836E+00  5.8092E-02
             2.0254E+00
 GRADIENT:  -2.8108E+02  7.4492E+01 -7.8679E-05  6.4700E+01  1.1813E+00  1.1568E+02  0.0000E+00  2.3342E-05  5.8291E+00 -6.4985E-02
             6.9822E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1508.38193484842        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2219
 NPARAMETR:  1.1819E+00  1.2099E+00  3.5309E+06  1.2965E+00  2.4480E+00  2.2409E+00  1.0000E-02  1.1932E+00  3.6137E+00  9.6132E-01
             6.8575E+00
 PARAMETER:  2.6710E-01  2.9050E-01  1.5177E+01  3.5969E-01  9.9526E-01  9.0687E-01 -7.2987E+00  2.7663E-01  1.3847E+00  6.0551E-02
             2.0253E+00
 GRADIENT:  -2.1445E+02  4.5938E+01 -1.1808E-05  5.9576E+01 -5.2048E+00 -5.2537E+01  0.0000E+00 -6.7031E-05  8.2291E+00  2.6628E-01
             5.9481E+01

0ITERATION NO.:   86    OBJECTIVE VALUE:  -1508.38193484842        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     2247
 NPARAMETR:  1.1816E+00  1.2099E+00  3.0337E+06  1.2960E+00  2.4480E+00  2.2409E+00  1.0000E-02  1.1940E+00  3.6139E+00  9.6102E-01
             6.8575E+00
 PARAMETER:  2.6710E-01  2.9050E-01  1.5177E+01  3.5969E-01  9.9526E-01  9.0687E-01 -7.2987E+00  2.7663E-01  1.3847E+00  6.0551E-02
             2.0253E+00
 GRADIENT:   1.2846E-01 -5.1692E-02  1.9952E-03  9.8803E-02 -1.2926E-02 -6.3855E-03  0.0000E+00 -6.1518E-02 -1.9769E-02  2.9416E-01
            -1.0492E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2247
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.4733E-03 -7.2320E-04  3.3449E-08  2.0716E-03 -9.4059E-03
 SE:             2.9165E-02  1.7709E-04  2.1096E-08  2.8237E-02  1.5469E-02
 N:                     100         100         100         100         100

 P VAL.:         7.9777E-01  4.4321E-05  1.1284E-01  9.4152E-01  5.4314E-01

 ETASHRINKSD(%)  2.2929E+00  9.9407E+01  1.0000E+02  5.4023E+00  4.8178E+01
 ETASHRINKVR(%)  4.5333E+00  9.9996E+01  1.0000E+02  1.0513E+01  7.3145E+01
 EBVSHRINKSD(%)  2.2269E+00  9.9633E+01  1.0000E+02  3.6071E+00  4.9854E+01
 EBVSHRINKVR(%)  4.4043E+00  9.9999E+01  1.0000E+02  7.0841E+00  7.4854E+01
 RELATIVEINF(%)  9.5500E+01  6.9262E-04  1.5493E-10  4.7799E+01  2.0298E+01
 EPSSHRINKSD(%)  8.9035E+00
 EPSSHRINKVR(%)  1.7014E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1508.3819348484158     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       145.70742491999499     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    71.60
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    19.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1508.382       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.18E+00  1.21E+00  3.53E+06  1.30E+00  2.45E+00  2.24E+00  1.00E-02  1.19E+00  3.61E+00  9.61E-01  6.86E+00
 


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
+        1.57E+02
 
 TH 2
+       -9.00E+00  1.69E+02
 
 TH 3
+        6.86E-09 -1.24E-07  3.32E-16
 
 TH 4
+       -8.92E-01  4.85E+01  1.12E-08  5.21E+01
 
 TH 5
+       -1.12E+00 -1.27E+01  2.31E-09 -6.28E+00  4.56E+01
 
 TH 6
+        2.67E-01  1.14E+00  6.48E-09  9.61E-02 -4.19E-01  3.61E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.48E+01  7.11E+00  3.02E-08  8.15E-01  9.10E-02  5.96E-01  0.00E+00  6.53E+00
 
 TH 9
+       -4.05E-01 -2.83E+01 -3.62E-09  1.63E+00  1.49E+00  1.67E-01  0.00E+00 -1.05E+00  1.20E+01
 
 TH10
+        5.41E+00 -2.29E+01 -3.55E-06  1.78E+01  6.09E+00  1.74E+00  0.00E+00  1.49E+01 -8.78E-02  6.05E+01
 
 TH11
+       -5.61E+00 -9.79E+00 -9.65E-10 -4.44E+00 -1.14E+00  5.32E-01  0.00E+00  7.76E-02  1.01E+00  3.99E+00  2.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       91.143
Stop Time:
Wed Sep 29 08:26:46 CDT 2021
