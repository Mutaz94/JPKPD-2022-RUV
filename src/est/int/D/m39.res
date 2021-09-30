Wed Sep 29 08:42:25 CDT 2021
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
$DATA ../../../../data/int/D/dat39.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   28775.3714926839        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.6523E+02  4.7385E+02  2.6131E+01  2.5835E+02  2.5037E+01 -1.7066E+03 -8.8941E+02 -7.0058E+01 -1.4815E+03 -4.5726E+02
            -5.9967E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -861.045195662433        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1265E+00  1.5399E+00  9.0749E-01  1.9539E+00  9.7985E-01  3.0595E+00  2.4593E+00  9.8691E-01  1.9499E+00  1.2631E+00
             1.3699E+01
 PARAMETER:  2.1915E-01  5.3171E-01  2.9279E-03  7.6981E-01  7.9645E-02  1.2183E+00  9.9990E-01  8.6823E-02  7.6776E-01  3.3359E-01
             2.7173E+00
 GRADIENT:  -4.8288E+01  7.4326E+00 -2.5151E+01  1.4127E+02  3.3818E+00  9.2690E+01 -7.7507E+01  4.5110E+00  8.3301E+00  2.6531E+01
             5.4812E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -947.166122718725        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0548E+00  1.1137E+00  3.2721E+01  2.7808E+00  2.8519E+00  3.7148E+00  6.4281E+00  6.7036E-01  2.6186E+00  1.6847E+00
             1.2628E+01
 PARAMETER:  1.5336E-01  2.0765E-01  3.5880E+00  1.1228E+00  1.1480E+00  1.4123E+00  1.9607E+00 -2.9994E-01  1.0627E+00  6.2156E-01
             2.6359E+00
 GRADIENT:  -3.2356E+01  2.2304E+01 -3.5795E+00  1.0953E+02  1.0935E+01  1.4213E+02  3.3804E+01  5.1548E-02  3.8981E+01  3.9224E+01
             5.1656E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1140.02007407647        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0761E+00  1.2280E+00  8.8365E+00  7.7598E-01  2.8398E+00  1.8536E+00  3.8996E+00  1.4380E-01  1.2497E+00  1.2127E-01
             9.1128E+00
 PARAMETER:  1.7334E-01  3.0536E-01  2.2789E+00 -1.5363E-01  1.1437E+00  7.1711E-01  1.4609E+00 -1.8393E+00  3.2291E-01 -2.0097E+00
             2.3097E+00
 GRADIENT:  -2.8099E+01 -5.6483E+01 -4.4973E+00 -6.8775E+01  7.0989E+01 -3.3350E+01  5.4956E+01  1.7905E-03  1.0604E+01  2.7584E-01
             2.6259E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1173.17160273371        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1084E+00  1.7741E+00  2.1652E+00  6.6074E-01  2.0339E+00  2.0667E+00  2.9450E+00  3.4234E-02  1.1011E+00  1.4323E-01
             7.9572E+00
 PARAMETER:  2.0295E-01  6.7329E-01  8.7250E-01 -3.1440E-01  8.0995E-01  8.2593E-01  1.1801E+00 -3.2745E+00  1.9634E-01 -1.8433E+00
             2.1741E+00
 GRADIENT:   8.7437E+00 -6.9582E+00  1.8246E-01 -1.2778E+01  1.0947E+00  4.7044E+00 -1.9079E+00  4.4682E-04  6.9662E+00  4.8778E-01
            -2.5464E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1178.75117799058        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      393
 NPARAMETR:  1.0822E+00  1.5602E+00  2.5721E+00  8.2642E-01  1.9409E+00  2.0263E+00  3.3100E+00  3.7308E-02  5.6893E-01  2.0830E-01
             8.0605E+00
 PARAMETER:  1.7897E-01  5.4481E-01  1.0447E+00 -9.0655E-02  7.6317E-01  8.0623E-01  1.2970E+00 -3.1886E+00 -4.6399E-01 -1.4688E+00
             2.1870E+00
 GRADIENT:  -2.2198E+01 -1.2883E+01 -4.4491E-01  2.2782E+00 -5.0462E+00 -2.7131E+01 -3.5385E+01  9.3763E-04  5.3201E-01  9.0170E-01
            -3.3975E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1184.35194087504        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      571
 NPARAMETR:  1.1514E+00  1.3894E+00  6.2858E+00  1.0126E+00  2.2331E+00  2.1752E+00  4.0932E+00  9.1096E-02  8.5323E-01  1.2595E-01
             8.2220E+00
 PARAMETER:  2.4097E-01  4.2886E-01  1.9383E+00  1.1255E-01  9.0339E-01  8.7713E-01  1.5093E+00 -2.2958E+00 -5.8721E-02 -1.9719E+00
             2.2068E+00
 GRADIENT:   5.5485E+00  9.5165E-01  7.8116E-01 -3.4737E+00 -2.4327E-01 -5.7792E-01 -2.9567E+00  2.4917E-03  8.0906E-01  2.5875E-01
            -1.1375E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1185.07017910591        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      747
 NPARAMETR:  1.1372E+00  1.0992E+00  9.0638E+00  1.1812E+00  2.2788E+00  2.1718E+00  4.6318E+00  1.1061E-01  1.0027E+00  7.6720E-02
             8.2272E+00
 PARAMETER:  2.2857E-01  1.9457E-01  2.3043E+00  2.6651E-01  9.2366E-01  8.7558E-01  1.6329E+00 -2.1018E+00  1.0267E-01 -2.4676E+00
             2.2074E+00
 GRADIENT:   4.6636E-01  7.3900E-02  4.0863E-01 -9.4448E-01 -5.1954E-01 -6.8034E-01 -2.2137E-01  2.5051E-03  3.2447E-01  8.4770E-02
            -5.1682E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1185.16810478968        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      927
 NPARAMETR:  1.1372E+00  1.0692E+00  8.8744E+00  1.1974E+00  2.2658E+00  2.1819E+00  4.7421E+00  5.3784E-02  9.9326E-01  2.5741E-02
             8.2324E+00
 PARAMETER:  2.2856E-01  1.6687E-01  2.2832E+00  2.8016E-01  9.1794E-01  8.8019E-01  1.6565E+00 -2.8228E+00  9.3241E-02 -3.5597E+00
             2.2081E+00
 GRADIENT:   4.3911E-01  3.9321E-01  1.4588E-01 -3.3953E-01 -2.7950E-02  9.1407E-01  2.0734E+00  6.5971E-04 -1.8317E-01  9.4555E-03
             1.9188E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1185.17374502730        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1102
 NPARAMETR:  1.1364E+00  1.0571E+00  8.7617E+00  1.2046E+00  2.2586E+00  2.1768E+00  4.7344E+00  2.1433E-02  1.0064E+00  1.1821E-02
             8.2319E+00
 PARAMETER:  2.2782E-01  1.5554E-01  2.2704E+00  2.8611E-01  9.1473E-01  8.7786E-01  1.6549E+00 -3.7428E+00  1.0635E-01 -4.3378E+00
             2.2080E+00
 GRADIENT:   8.9825E-02  3.9785E-02 -2.7567E-02  6.4505E-01  7.0676E-02  1.4333E-01  7.1112E-01  1.0999E-04  6.6858E-03  1.9895E-03
             1.9474E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1185.18225300089        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1279
 NPARAMETR:  1.1374E+00  1.0533E+00  8.7907E+00  1.2070E+00  2.2583E+00  2.1809E+00  4.7658E+00  1.5164E-02  1.0035E+00  1.0000E-02
             8.2351E+00
 PARAMETER:  2.2875E-01  1.5195E-01  2.2737E+00  2.8813E-01  9.1463E-01  8.7976E-01  1.6615E+00 -4.0888E+00  1.0350E-01 -5.4476E+00
             2.2084E+00
 GRADIENT:   4.6371E-01  2.5061E-01 -2.3691E-02  3.1758E-01  5.0494E-02  7.9210E-01  1.7718E+00  5.5228E-05 -6.7203E-02  0.0000E+00
             7.2202E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1185.18763259187        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1463             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1358E+00  1.0425E+00  8.8209E+00  1.2092E+00  2.2576E+00  2.1794E+00  4.7966E+00  1.0000E-02  1.0048E+00  1.0000E-02
             8.2278E+00
 PARAMETER:  2.2735E-01  1.4164E-01  2.2771E+00  2.8997E-01  9.1432E-01  8.7907E-01  1.6679E+00 -4.5144E+00  1.0476E-01 -5.4476E+00
             2.2075E+00
 GRADIENT:   1.7925E+01  1.8845E+00  2.4309E-01  7.9081E+00  4.1332E+00  2.5193E+01  6.3805E+01  2.1522E-05  2.6438E-01  0.0000E+00
             3.9243E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1185.18899020373        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1646             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1363E+00  1.0379E+00  8.8512E+00  1.2119E+00  2.2574E+00  2.1800E+00  4.8042E+00  1.0000E-02  1.0093E+00  1.0000E-02
             8.2305E+00
 PARAMETER:  2.2778E-01  1.3721E-01  2.2806E+00  2.9215E-01  9.1422E-01  8.7934E-01  1.6695E+00 -4.5655E+00  1.0930E-01 -5.4476E+00
             2.2079E+00
 GRADIENT:   1.8097E+01  1.7373E+00  2.2554E-01  7.8526E+00  4.1118E+00  2.5293E+01  6.3991E+01  0.0000E+00  3.8924E-01  0.0000E+00
             3.9983E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1185.19004766431        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1824
 NPARAMETR:  1.1365E+00  1.0382E+00  8.8827E+00  1.2140E+00  2.2572E+00  2.1803E+00  4.8057E+00  1.0000E-02  1.0078E+00  1.0000E-02
             8.2309E+00
 PARAMETER:  2.2797E-01  1.3753E-01  2.2841E+00  2.9389E-01  9.1413E-01  8.7947E-01  1.6698E+00 -4.5655E+00  1.0774E-01 -5.4476E+00
             2.2079E+00
 GRADIENT:   2.3042E-01  1.5236E-01 -1.6990E-02 -8.4433E-02 -1.6331E-01  6.8306E-01  2.3928E+00  0.0000E+00 -7.1737E-02  0.0000E+00
            -1.4165E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1185.19176669412        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2011
 NPARAMETR:  1.1359E+00  1.0330E+00  8.9405E+00  1.2162E+00  2.2580E+00  2.1792E+00  4.8222E+00  1.0000E-02  1.0106E+00  1.0000E-02
             8.2279E+00
 PARAMETER:  2.2739E-01  1.3249E-01  2.2906E+00  2.9577E-01  9.1447E-01  8.7895E-01  1.6732E+00 -4.5655E+00  1.1056E-01 -5.4476E+00
             2.2075E+00
 GRADIENT:   4.3796E-02  1.2351E-01 -1.5191E-02 -4.5384E-01 -1.7385E-01  4.9123E-01  2.7769E+00  0.0000E+00 -4.7253E-03  0.0000E+00
            -6.9280E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1185.19315761229        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2205             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1362E+00  1.0254E+00  8.9868E+00  1.2197E+00  2.2599E+00  2.1800E+00  4.8275E+00  1.0000E-02  1.0135E+00  1.0000E-02
             8.2302E+00
 PARAMETER:  2.2770E-01  1.2512E-01  2.2958E+00  2.9857E-01  9.1530E-01  8.7935E-01  1.6743E+00 -4.5655E+00  1.1338E-01 -5.4476E+00
             2.2078E+00
 GRADIENT:   1.8062E+01  1.5114E+00  2.0719E-01  8.6433E+00  4.1018E+00  2.5309E+01  6.4365E+01  0.0000E+00  2.5621E-01  0.0000E+00
             3.9636E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1185.19385126320        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2393             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1363E+00  1.0237E+00  9.0291E+00  1.2209E+00  2.2606E+00  2.1799E+00  4.8348E+00  1.0315E-02  1.0165E+00  1.0000E-02
             8.2308E+00
 PARAMETER:  2.2775E-01  1.2343E-01  2.3005E+00  2.9961E-01  9.1564E-01  8.7929E-01  1.6758E+00 -4.4741E+00  1.1634E-01 -5.4476E+00
             2.2079E+00
 GRADIENT:   1.8089E+01  1.5113E+00  2.1403E-01  8.4760E+00  4.0552E+00  2.5286E+01  6.4689E+01  2.9048E-05  3.5836E-01  0.0000E+00
             3.9862E+01

0ITERATION NO.:   83    OBJECTIVE VALUE:  -1185.19405880332        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     2491
 NPARAMETR:  1.1353E+00  1.0254E+00  9.1700E+00  1.2249E+00  2.2598E+00  2.1785E+00  4.8512E+00  1.0000E-02  1.0146E+00  1.0000E-02
             8.2239E+00
 PARAMETER:  2.2799E-01  1.2440E-01  2.3014E+00  3.0001E-01  9.1563E-01  8.7946E-01  1.6756E+00 -4.6997E+00  1.1561E-01 -5.4476E+00
             2.2080E+00
 GRADIENT:   9.5571E-02 -1.1429E-02 -5.9207E-02 -2.7954E-01  4.2023E-02  5.8987E-02 -1.8034E-01  0.0000E+00  1.1114E-02  0.0000E+00
             4.2350E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2491
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0501E-02  2.9471E-02 -1.7483E-05 -6.5855E-02  3.4767E-05
 SE:             2.9045E-02  2.4858E-02  1.9169E-05  1.3315E-02  1.3447E-04
 N:                     100         100         100         100         100

 P VAL.:         7.1769E-01  2.3579E-01  3.6174E-01  7.5925E-07  7.9598E-01

 ETASHRINKSD(%)  2.6940E+00  1.6722E+01  9.9936E+01  5.5392E+01  9.9550E+01
 ETASHRINKVR(%)  5.3154E+00  3.0648E+01  1.0000E+02  8.0101E+01  9.9998E+01
 EBVSHRINKSD(%)  3.1860E+00  1.2150E+01  9.9922E+01  5.8949E+01  9.9487E+01
 EBVSHRINKVR(%)  6.2705E+00  2.2824E+01  1.0000E+02  8.3148E+01  9.9997E+01
 RELATIVEINF(%)  9.3578E+01  3.2716E+01  9.3388E-06  7.1923E+00  4.1804E-04
 EPSSHRINKSD(%)  7.0721E+00
 EPSSHRINKVR(%)  1.3644E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1185.1940588033247     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       468.89530096508611     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    82.04
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    16.69
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1185.194       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  1.02E+00  9.04E+00  1.22E+00  2.26E+00  2.18E+00  4.83E+00  1.00E-02  1.02E+00  1.00E-02  8.23E+00
 


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
+        1.44E+02
 
 TH 2
+       -4.93E+00  8.95E+00
 
 TH 3
+        5.58E-01 -8.38E-02  2.64E-03
 
 TH 4
+       -5.85E+01  3.53E+01 -4.72E-01  1.50E+02
 
 TH 5
+       -8.20E+00 -1.13E+00 -2.14E-02 -2.00E+00  6.94E-01
 
 TH 6
+        1.24E+01 -2.57E+00  6.39E-02 -1.32E+01 -3.65E-01  1.60E+00
 
 TH 7
+        7.37E+00 -2.91E+00  4.81E-02 -1.30E+01  5.52E-03  1.28E+00  1.18E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.54E+01 -8.09E+00  1.15E-01 -3.49E+01  3.34E-01  3.18E+00  3.07E+00  0.00E+00  8.16E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.83E-01 -2.01E+00  1.41E-02 -7.62E+00  3.19E-01  5.15E-01  6.01E-01  0.00E+00  1.76E+00  0.00E+00  8.54E-01
 
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
+        1.64E+02
 
 TH 2
+       -1.07E+00  3.55E+01
 
 TH 3
+        9.97E-02  2.76E-01  1.18E-01
 
 TH 4
+       -4.22E+00  4.09E+01 -6.85E-01  1.55E+02
 
 TH 5
+       -1.72E+00 -8.70E+00 -2.20E+00  6.38E-01  5.35E+01
 
 TH 6
+       -1.43E-01 -1.84E-01  1.62E-02  7.31E-01 -9.59E-02  3.45E+01
 
 TH 7
+        4.52E-01  4.08E+00 -6.60E-02 -1.20E+01  1.58E+00 -1.09E-01  4.95E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.94E-02 -3.51E+00 -3.42E-01 -3.05E+01  4.57E+00 -6.22E-01  2.79E+00  0.00E+00  2.24E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.03E+00 -3.45E+00 -9.95E-03 -1.07E+01  5.11E-01  1.60E+00  4.66E-01  0.00E+00  3.18E+00  0.00E+00  1.54E+01
 
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
+        1.70E+02
 
 TH 2
+        6.08E+01  3.57E+01
 
 TH 3
+        9.01E-01  3.86E-01  6.74E-02
 
 TH 4
+        9.27E+01  4.19E+01  6.19E-01  1.63E+02
 
 TH 5
+       -2.89E+01 -1.17E+01 -1.62E+00 -2.31E+01  4.22E+01
 
 TH 6
+        6.43E+00  4.82E+00 -1.45E-01 -2.85E+01  2.32E+00  3.02E+01
 
 TH 7
+        3.14E-01  4.61E+00 -1.26E-01 -1.38E+01  3.41E+00  6.28E+00  5.80E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.47E+01 -5.32E+00 -2.42E-01 -4.32E+01  9.62E+00  9.70E+00  3.94E+00  0.00E+00  2.17E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.16E+01 -1.56E+01 -3.33E-01 -8.44E+01  1.29E+01  1.56E+01  1.13E+01  0.00E+00  3.39E+01  0.00E+00  5.68E+02
 
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
 #CPUT: Total CPU Time in Seconds,       98.833
Stop Time:
Wed Sep 29 08:44:06 CDT 2021
