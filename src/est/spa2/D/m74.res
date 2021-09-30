Thu Sep 30 09:43:20 CDT 2021
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
$DATA ../../../../data/spa2/D/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   21583.9856429247        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9899E+02  4.0841E+02 -1.0688E+01  2.1307E+02  2.9240E+02 -2.8935E+03 -1.1078E+03 -4.8643E+01 -1.8025E+03 -8.6037E+02
            -4.0945E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -601.220546662601        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.6253E+00  1.2016E+00  9.0849E-01  1.8020E+00  1.0760E+00  2.8787E+00  1.9039E+00  9.7558E-01  2.0570E+00  1.1205E+00
             1.2804E+01
 PARAMETER:  5.8570E-01  2.8362E-01  4.0321E-03  6.8892E-01  1.7323E-01  1.1573E+00  7.4390E-01  7.5280E-02  8.2125E-01  2.1377E-01
             2.6497E+00
 GRADIENT:   5.3062E+01 -1.7604E+00 -2.9349E+01  3.9360E+01  2.9247E+01  5.5787E+01 -1.1432E+01  2.7081E+00 -4.4972E+01  1.2922E+01
             1.4887E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -675.798229079242        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.5526E+00  1.8965E+00  3.1915E+00  1.5105E+00  2.4382E+00  4.0439E+00  5.6870E+00  7.7727E-01  4.6616E+00  8.9751E-01
             1.1068E+01
 PARAMETER:  5.3990E-01  7.4000E-01  1.2605E+00  5.1241E-01  9.9127E-01  1.4972E+00  1.8382E+00 -1.5197E-01  1.6394E+00 -8.1255E-03
             2.5040E+00
 GRADIENT:   4.2267E+01  7.7279E+00 -2.5731E+01  2.8491E+01  1.0123E+01  1.5369E+02  1.0193E+02  2.0706E-01  6.8633E+01  4.2642E+00
             1.9533E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -730.208697492643        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      343
 NPARAMETR:  1.2921E+00  5.4328E-01  1.6125E+01  2.2841E+00  2.0363E+00  2.9827E+00  8.0094E+00  2.6839E+00  2.7626E+00  1.1565E+00
             1.0714E+01
 PARAMETER:  3.5624E-01 -5.1012E-01  2.8804E+00  9.2596E-01  8.1113E-01  1.1928E+00  2.1806E+00  1.0873E+00  1.1162E+00  2.4541E-01
             2.4715E+00
 GRADIENT:  -3.6712E+00  3.1003E+00 -1.4802E-01  1.5300E+01 -2.1868E+01  2.9836E+01  9.5718E+00  5.6770E-01  3.7907E+01  1.0380E+01
             1.3036E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -765.268802086241        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.2586E+00  3.8297E-01  4.1192E+01  1.6985E+00  2.4765E+00  2.6186E+00  8.2275E+00  5.6678E+00  1.4026E+00  6.2807E-01
             9.5090E+00
 PARAMETER:  3.3002E-01 -8.5979E-01  3.8182E+00  6.2977E-01  1.0069E+00  1.0627E+00  2.2075E+00  1.8348E+00  4.3832E-01 -3.6510E-01
             2.3522E+00
 GRADIENT:   9.3241E+00 -7.4212E+00 -2.7703E-01 -1.0116E+01 -1.1454E-02  2.5929E+00  1.4647E+00  1.5507E-01  2.8454E+00  2.5036E+00
             1.8949E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -769.534975961844        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  1.2025E+00  1.0872E+00  1.4533E+01  1.2962E+00  2.3562E+00  2.5985E+00  6.1788E+00  1.1939E-01  1.0022E+00  1.7604E-01
             9.3089E+00
 PARAMETER:  2.8443E-01  1.8361E-01  2.7764E+00  3.5947E-01  9.5707E-01  1.0549E+00  1.9211E+00 -2.0254E+00  1.0217E-01 -1.6371E+00
             2.3310E+00
 GRADIENT:  -2.6259E+00 -3.2714E-01 -1.1062E-01 -1.5422E+00 -6.8039E-01 -1.6665E+00 -1.0398E+00  3.1865E-04  2.5234E+00  2.1262E-01
            -1.0166E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -769.940207963936        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      873
 NPARAMETR:  1.2200E+00  1.2043E+00  1.6879E+01  1.2143E+00  2.3845E+00  2.6295E+00  6.0850E+00  3.9068E-02  7.5407E-01  4.0059E-02
             9.4249E+00
 PARAMETER:  2.9885E-01  2.8591E-01  2.9261E+00  2.9414E-01  9.6899E-01  1.0668E+00  1.9058E+00 -3.1424E+00 -1.8227E-01 -3.1174E+00
             2.3434E+00
 GRADIENT:   4.8768E-01 -1.4264E-01  1.9218E-02 -1.7340E+00 -7.6846E-01  2.6858E+00  2.4901E+00  2.1723E-05  4.9002E-02  1.0943E-02
             9.9453E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -769.952567368425        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1049
 NPARAMETR:  1.2200E+00  1.1913E+00  2.0686E+01  1.2326E+00  2.4265E+00  2.6106E+00  6.0875E+00  2.5331E-02  7.8449E-01  1.7574E-02
             9.4304E+00
 PARAMETER:  2.9882E-01  2.7501E-01  3.1295E+00  3.0916E-01  9.8644E-01  1.0596E+00  1.9062E+00 -3.5757E+00 -1.4273E-01 -3.9414E+00
             2.3439E+00
 GRADIENT:   2.2744E-01  3.8576E-02 -8.3302E-03  2.4380E-01  1.1517E-01  2.9685E-01  9.8882E-01  5.9584E-06 -4.1231E-02  2.0630E-03
             5.0526E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -769.978194570106        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1226
 NPARAMETR:  1.2215E+00  1.1807E+00  2.1521E+01  1.2375E+00  2.4261E+00  2.6191E+00  6.1623E+00  2.5150E-02  7.8292E-01  1.0000E-02
             9.4303E+00
 PARAMETER:  3.0010E-01  2.6613E-01  3.1690E+00  3.1313E-01  9.8630E-01  1.0628E+00  1.9185E+00 -3.5829E+00 -1.4473E-01 -4.6635E+00
             2.3439E+00
 GRADIENT:   6.3914E-01  2.5862E-01 -5.2981E-06 -2.2514E-01 -6.5431E-02  1.3752E+00  3.0045E+00  5.4888E-06 -1.4028E-01  0.0000E+00
             5.8636E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -769.986983121854        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1407
 NPARAMETR:  1.2198E+00  1.1681E+00  2.2116E+01  1.2421E+00  2.4292E+00  2.6270E+00  6.1907E+00  2.2763E-02  8.0056E-01  1.0000E-02
             9.4248E+00
 PARAMETER:  2.9870E-01  2.5535E-01  3.1963E+00  3.1680E-01  9.8758E-01  1.0658E+00  1.9230E+00 -3.6826E+00 -1.2245E-01 -4.6632E+00
             2.3433E+00
 GRADIENT:   3.1214E-01  8.1712E-02 -4.2517E-03 -1.5891E+00  1.5299E-01  2.3562E+00  3.6131E+00  4.2217E-06  1.8489E-01  0.0000E+00
             1.0586E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -769.989103680235        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1582
 NPARAMETR:  1.2214E+00  1.1598E+00  2.2831E+01  1.2499E+00  2.4311E+00  2.6179E+00  6.1916E+00  2.2031E-02  8.0666E-01  1.0000E-02
             9.4306E+00
 PARAMETER:  2.9997E-01  2.4829E-01  3.2281E+00  3.2307E-01  9.8836E-01  1.0624E+00  1.9232E+00 -3.7153E+00 -1.1485E-01 -4.6632E+00
             2.3440E+00
 GRADIENT:   5.7845E-01  1.2741E-01 -1.1614E-03 -4.4532E-02 -4.3953E-02  1.2616E+00  2.7145E+00  3.7577E-06 -6.7934E-02  0.0000E+00
             7.4228E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -769.997509729453        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1767             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2194E+00  1.1468E+00  2.3231E+01  1.2528E+00  2.4328E+00  2.6273E+00  6.2269E+00  1.9948E-02  8.1309E-01  1.0000E-02
             9.4206E+00
 PARAMETER:  2.9838E-01  2.3696E-01  3.2455E+00  3.2537E-01  9.8906E-01  1.0660E+00  1.9289E+00 -3.8146E+00 -1.0691E-01 -4.6632E+00
             2.3429E+00
 GRADIENT:   1.8218E+01  2.2839E+00  8.8307E-04  5.3672E+00  1.0316E+00  4.3083E+01  9.5669E+01  3.1316E-06  1.3233E-01  0.0000E+00
             2.8476E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -770.000520683688        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1947
 NPARAMETR:  1.2193E+00  1.1428E+00  2.3650E+01  1.2570E+00  2.4330E+00  2.6274E+00  6.2432E+00  1.6934E-02  8.1810E-01  1.0000E-02
             9.4189E+00
 PARAMETER:  2.9829E-01  2.3346E-01  3.2634E+00  3.2871E-01  9.8913E-01  1.0660E+00  1.9315E+00 -3.9784E+00 -1.0077E-01 -4.6632E+00
             2.3427E+00
 GRADIENT:   2.3503E-01  7.8368E-02 -1.1663E-03 -5.1922E-01 -2.0183E-02  2.4385E+00  3.5898E+00  2.1086E-06 -3.1339E-02  0.0000E+00
            -5.3625E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -770.006024680275        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2134
 NPARAMETR:  1.2189E+00  1.1248E+00  2.4927E+01  1.2680E+00  2.4366E+00  2.6274E+00  6.2836E+00  1.1634E-02  8.3897E-01  1.0000E-02
             9.4166E+00
 PARAMETER:  2.9795E-01  2.1757E-01  3.3159E+00  3.3740E-01  9.9060E-01  1.0660E+00  1.9379E+00 -4.3538E+00 -7.5578E-02 -4.6632E+00
             2.3425E+00
 GRADIENT:   1.5255E-01  5.4327E-02 -1.9039E-03 -6.3336E-01 -1.0441E-02  2.4545E+00  3.8598E+00  8.9226E-07  7.4816E-02  0.0000E+00
            -6.1355E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -770.006511365658        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2318
 NPARAMETR:  1.2191E+00  1.1214E+00  2.5267E+01  1.2696E+00  2.4387E+00  2.6272E+00  6.2895E+00  1.3029E-02  8.4183E-01  1.0000E-02
             9.4178E+00
 PARAMETER:  2.9812E-01  2.1454E-01  3.3295E+00  3.3868E-01  9.9144E-01  1.0659E+00  1.9389E+00 -4.2406E+00 -7.2177E-02 -4.6632E+00
             2.3426E+00
 GRADIENT:   1.8580E-01  6.2190E-03 -3.1833E-03 -7.4927E-01  5.4647E-02  2.4398E+00  3.8873E+00  1.1036E-06  1.0015E-01  0.0000E+00
             1.4774E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -770.007120542771        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2500
 NPARAMETR:  1.2191E+00  1.1188E+00  2.5657E+01  1.2713E+00  2.4400E+00  2.6272E+00  6.2946E+00  1.3237E-02  8.4403E-01  1.0000E-02
             9.4177E+00
 PARAMETER:  2.9815E-01  2.1221E-01  3.3448E+00  3.4004E-01  9.9202E-01  1.0659E+00  1.9397E+00 -4.2248E+00 -6.9566E-02 -4.6632E+00
             2.3426E+00
 GRADIENT:   1.9095E-01  3.9052E-03 -2.9151E-03 -6.4749E-01  5.5375E-02  2.4395E+00  3.8688E+00  1.0956E-06  8.2928E-02  0.0000E+00
             1.0287E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -770.007587601334        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2687
 NPARAMETR:  1.2192E+00  1.1168E+00  2.5929E+01  1.2725E+00  2.4407E+00  2.6272E+00  6.2968E+00  1.2272E-02  8.4457E-01  1.0000E-02
             9.4176E+00
 PARAMETER:  2.9818E-01  2.1048E-01  3.3554E+00  3.4101E-01  9.9228E-01  1.0659E+00  1.9400E+00 -4.3005E+00 -6.8932E-02 -4.6632E+00
             2.3426E+00
 GRADIENT:   1.9675E-01 -2.4350E-03 -2.0234E-03 -4.3313E-01  2.5403E-02  2.4467E+00  3.7708E+00  9.1588E-07  3.1345E-02  0.0000E+00
            -1.7100E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -770.007851517874        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2870
 NPARAMETR:  1.2193E+00  1.1156E+00  2.6155E+01  1.2736E+00  2.4412E+00  2.6268E+00  6.2997E+00  1.0906E-02  8.4576E-01  1.0000E-02
             9.4178E+00
 PARAMETER:  2.9827E-01  2.0943E-01  3.3641E+00  3.4183E-01  9.9248E-01  1.0657E+00  1.9405E+00 -4.4185E+00 -6.7517E-02 -4.6632E+00
             2.3426E+00
 GRADIENT:   2.1833E-01  1.1972E-02 -1.4441E-03 -3.2712E-01  5.4114E-03  2.3992E+00  3.7567E+00  6.9702E-07  1.4645E-02  0.0000E+00
            -5.2493E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -770.008067836233        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     3049
 NPARAMETR:  1.2191E+00  1.1136E+00  2.6690E+01  1.2759E+00  2.4412E+00  2.6272E+00  6.3043E+00  1.0000E-02  8.4627E-01  1.0000E-02
             9.4156E+00
 PARAMETER:  2.9816E-01  2.0831E-01  3.3705E+00  3.4244E-01  9.9266E-01  1.0659E+00  1.9409E+00 -4.5073E+00 -6.6062E-02 -4.6632E+00
             2.3425E+00
 GRADIENT:   3.2566E-03  1.6322E-02 -1.2407E-03 -2.7327E-01  9.2075E-03 -3.0931E-03 -3.3036E-02  4.9897E-06  1.2507E-02  0.0000E+00
             9.9920E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3049
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0825E-02  2.9287E-02  1.7698E-06 -6.2877E-02  9.3245E-06
 SE:             2.8490E-02  2.4569E-02  1.1764E-06  1.1463E-02  4.8498E-05
 N:                     100         100         100         100         100

 P VAL.:         7.0398E-01  2.3324E-01  1.3248E-01  4.1418E-08  8.4753E-01

 ETASHRINKSD(%)  4.5565E+00  1.7692E+01  9.9996E+01  6.1597E+01  9.9838E+01
 ETASHRINKVR(%)  8.9054E+00  3.2254E+01  1.0000E+02  8.5252E+01  1.0000E+02
 EBVSHRINKSD(%)  3.8916E+00  1.2229E+01  9.9993E+01  6.7386E+01  9.9759E+01
 EBVSHRINKVR(%)  7.6318E+00  2.2962E+01  1.0000E+02  8.9363E+01  9.9999E+01
 RELATIVEINF(%)  9.1357E+01  4.2037E+01  8.2022E-08  5.5028E+00  9.6706E-05
 EPSSHRINKSD(%)  8.4443E+00
 EPSSHRINKVR(%)  1.6176E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -770.00806783623295     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       332.71817200937414     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    72.00
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -770.008       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.22E+00  1.11E+00  2.63E+01  1.27E+00  2.44E+00  2.63E+00  6.30E+00  1.00E-02  8.47E-01  1.00E-02  9.42E+00
 


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
+        9.81E+01
 
 TH 2
+       -1.07E+00  1.81E+01
 
 TH 3
+       -1.06E-03  1.54E-03  1.29E-04
 
 TH 4
+       -3.11E+00  2.35E+01  3.63E-03  1.42E+02
 
 TH 5
+       -6.96E-01 -2.36E+00 -2.75E-02 -9.37E+00  1.02E+01
 
 TH 6
+        5.52E-01 -2.17E-01  1.02E-03  1.53E+00 -4.78E-01  2.42E+01
 
 TH 7
+        1.77E-01  2.20E+00 -6.50E-04 -9.17E+00  3.89E-01 -2.27E-01  2.89E+00
 
 TH 8
+       -1.76E-02 -2.29E-01 -6.55E-04 -1.30E-01 -9.53E-03  3.40E-03  5.61E-03  2.48E-01
 
 TH 9
+        4.58E-01 -1.63E+00 -4.27E-03 -3.93E+01  2.49E+00 -3.32E-01  1.96E+00 -3.43E-01  1.91E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.14E+00 -1.97E+00  6.88E-04 -1.00E+01  4.38E-01  9.29E-01  4.73E-01 -3.50E-03  3.28E+00  0.00E+00  7.22E+00
 
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
 #CPUT: Total CPU Time in Seconds,       82.836
Stop Time:
Thu Sep 30 09:44:45 CDT 2021
