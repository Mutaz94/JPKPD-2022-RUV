Sun Oct 24 02:26:56 CDT 2021
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
$DATA ../../../../data/SD4/A3/dat26.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   445.960475108299        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6343E+02  1.2123E+02  7.0025E+01  6.3668E+01  2.3112E+02  6.4415E+01 -1.0504E+02 -3.6119E+01 -2.0562E+02 -1.6077E+02
            -3.6576E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1194.38193399541        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0447E+00  9.5641E-01  9.6192E-01  1.0884E+00  9.4058E-01  7.2539E-01  1.0222E+00  9.9674E-01  1.1464E+00  9.6015E-01
             5.2578E+00
 PARAMETER:  1.4374E-01  5.5433E-02  6.1173E-02  1.8468E-01  3.8745E-02 -2.2104E-01  1.2194E-01  9.6733E-02  2.3664E-01  5.9338E-02
             1.7597E+00
 GRADIENT:   7.0024E+01 -2.7389E+01 -1.7504E+01 -2.5318E+01  5.2776E+00 -3.3932E+00  8.6357E+00  6.5476E+00  1.8473E+01  2.1682E+01
             1.5477E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1208.97851679125        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0180E+00  5.4838E-01  4.4855E-01  1.3523E+00  4.4217E-01  7.9552E-01  8.3148E-01  3.8335E-01  1.1787E+00  2.5215E-01
             4.9217E+00
 PARAMETER:  1.1781E-01 -5.0080E-01 -7.0173E-01  4.0177E-01 -7.1607E-01 -1.2876E-01 -8.4547E-02 -8.5880E-01  2.6438E-01 -1.2777E+00
             1.6937E+00
 GRADIENT:  -5.1190E+01  3.9850E+01  5.7292E+00  1.1333E+02 -2.8407E+01  5.6128E+00 -2.4444E-01  3.2935E+00  1.0185E+01  2.6955E+00
             1.3648E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1225.07220505874        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  9.7220E-01  5.1652E-01  1.9961E-01  1.1611E+00  2.8414E-01  8.4032E-01  6.3188E-01  3.2527E-01  1.4315E+00  1.3546E-01
             3.8009E+00
 PARAMETER:  7.1809E-02 -5.6065E-01 -1.5114E+00  2.4940E-01 -1.1583E+00 -7.3969E-02 -3.5905E-01 -1.0231E+00  4.5873E-01 -1.8991E+00
             1.4352E+00
 GRADIENT:  -9.0409E+01  5.9465E+01  2.1823E+01  1.1306E+02 -2.7414E+01 -2.2720E-01 -8.5323E+00 -1.0645E+00 -2.7015E+00 -1.2358E+00
             3.3882E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1227.65278498611        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      344
 NPARAMETR:  9.7405E-01  5.3330E-01  1.9635E-01  1.1226E+00  2.8501E-01  8.4800E-01  6.9519E-01  3.9073E-01  1.3964E+00  1.2672E-01
             3.7086E+00
 PARAMETER:  7.3705E-02 -5.2867E-01 -1.5279E+00  2.1565E-01 -1.1552E+00 -6.4876E-02 -2.6358E-01 -8.3975E-01  4.3393E-01 -1.9658E+00
             1.4106E+00
 GRADIENT:  -9.9342E+01  5.6098E+01  1.8802E+01  8.5418E+01 -7.2374E+01  1.6895E+00 -9.7346E+00 -1.7393E+00 -1.1717E+01 -1.0914E+00
            -1.8914E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1235.56770112808        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  9.9446E-01  5.3481E-01  1.6532E-01  9.9314E-01  2.6376E-01  8.3067E-01  9.1423E-01  1.0711E+00  1.4493E+00  7.2037E-02
             3.2544E+00
 PARAMETER:  9.4446E-02 -5.2584E-01 -1.6999E+00  9.3118E-02 -1.2327E+00 -8.5523E-02  1.0331E-02  1.6869E-01  4.7110E-01 -2.5306E+00
             1.2800E+00
 GRADIENT:  -2.8405E+00  3.0014E+01  2.2704E+01  6.0901E+00 -5.4339E+01 -1.5477E+00 -5.0041E+00 -3.9633E+00 -1.8041E+01 -2.5505E-01
            -4.1418E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1238.07616449978        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      696
 NPARAMETR:  9.9511E-01  5.0346E-01  1.5183E-01  9.8432E-01  2.5198E-01  8.3002E-01  9.1317E-01  1.1010E+00  1.5998E+00  4.4636E-02
             3.4497E+00
 PARAMETER:  9.5097E-02 -5.8625E-01 -1.7850E+00  8.4196E-02 -1.2784E+00 -8.6308E-02  9.1660E-03  1.9620E-01  5.6991E-01 -3.0092E+00
             1.3383E+00
 GRADIENT:   6.0083E-01 -2.4372E-01  7.2195E-01  1.6402E+00 -1.6579E+00 -1.3528E-01  3.4871E-01  4.6664E-01  1.5413E+00 -7.7559E-02
             1.7941E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1238.12106300883        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      851
 NPARAMETR:  9.9403E-01  5.1481E-01  1.5163E-01  9.7966E-01  2.5520E-01  8.3026E-01  9.0373E-01  1.0976E+00  1.5929E+00  5.4882E-02
             3.4426E+00
 PARAMETER:  9.4014E-02 -5.6395E-01 -1.7863E+00  7.9449E-02 -1.2657E+00 -8.6021E-02 -1.2197E-03  1.9311E-01  5.6555E-01 -2.8026E+00
             1.3362E+00
 GRADIENT:   2.6329E+01  6.9121E+00  1.3313E+01  4.9572E+00  6.4076E+01  1.4470E+00  4.9710E-01  2.5326E-01  8.2387E+00 -6.5120E-02
             1.1651E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1238.63039813902        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1033
 NPARAMETR:  9.9513E-01  5.1700E-01  1.5270E-01  9.8267E-01  2.5722E-01  8.2983E-01  8.5744E-01  1.0950E+00  1.6045E+00  2.4637E-01
             3.4173E+00
 PARAMETER:  9.5113E-02 -5.5972E-01 -1.7793E+00  8.2518E-02 -1.2578E+00 -8.6536E-02 -5.3807E-02  1.9080E-01  5.7279E-01 -1.3009E+00
             1.3288E+00
 GRADIENT:   6.3656E+00 -8.5007E+00 -4.3721E+00  1.4446E+00  2.6765E+01  2.7040E-01  1.9297E+00  5.8122E-01  1.0037E+00 -8.2835E-02
             4.6332E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1239.54523747516        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1210
 NPARAMETR:  9.9150E-01  5.0539E-01  1.3995E-01  9.6395E-01  2.4231E-01  8.2957E-01  6.2212E-01  1.0410E+00  1.7386E+00  4.3237E-01
             3.3562E+00
 PARAMETER:  9.1464E-02 -5.8243E-01 -1.8665E+00  6.3287E-02 -1.3175E+00 -8.6845E-02 -3.7462E-01  1.4016E-01  6.5309E-01 -7.3847E-01
             1.3108E+00
 GRADIENT:   6.5708E+00  1.6080E+01  9.2366E+00 -1.0379E-03 -5.2150E+00  2.8811E-01 -3.7820E-01  4.0135E-01  3.7568E+00  1.0925E+00
             3.1046E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1245.22977280102        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1394
 NPARAMETR:  9.7950E-01  3.4312E-01  9.0243E-02  8.7062E-01  1.7340E-01  8.4465E-01  2.5951E-01  4.8176E-01  2.0355E+00  6.5040E-01
             3.2964E+00
 PARAMETER:  7.9292E-02 -9.6967E-01 -2.3053E+00 -3.8548E-02 -1.6522E+00 -6.8835E-02 -1.2490E+00 -6.3031E-01  8.1075E-01 -3.3017E-01
             1.2928E+00
 GRADIENT:  -1.7807E+01  3.0550E+00 -1.0871E+01  9.2476E-01  5.3202E-01  2.0942E+00  7.8532E-01  7.9831E-01 -4.7571E+00 -5.8726E+00
            -1.2264E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1245.94735789256        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1571
 NPARAMETR:  9.8636E-01  3.2925E-01  8.9423E-02  8.7722E-01  1.6981E-01  8.3540E-01  1.3870E-01  2.9875E-01  2.1432E+00  7.0016E-01
             3.2896E+00
 PARAMETER:  8.6264E-02 -1.0109E+00 -2.3144E+00 -3.1001E-02 -1.6731E+00 -7.9849E-02 -1.8755E+00 -1.1081E+00  8.6229E-01 -2.5644E-01
             1.2908E+00
 GRADIENT:   1.9275E+00 -3.0810E+00  4.3851E-01  2.8358E+00  5.1573E+00  1.0165E+00  2.4993E-01 -7.3098E-02 -1.8543E-01  2.2951E+00
             1.9175E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1246.15883425714        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1746
 NPARAMETR:  9.8241E-01  3.3158E-01  8.8342E-02  8.6886E-01  1.6910E-01  8.2861E-01  6.2355E-02  1.7692E-01  2.2059E+00  6.9559E-01
             3.2716E+00
 PARAMETER:  8.2252E-02 -1.0039E+00 -2.3265E+00 -4.0573E-02 -1.6773E+00 -8.8007E-02 -2.6749E+00 -1.6320E+00  8.9111E-01 -2.6300E-01
             1.2853E+00
 GRADIENT:  -2.3531E+00  1.3823E+00  3.7735E-01 -5.9908E-01 -3.1108E+00 -9.2297E-01  5.0314E-02  9.6139E-02  1.6610E+00 -5.9236E-01
            -1.0918E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1246.22395730588        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1921
 NPARAMETR:  9.8254E-01  3.3103E-01  8.8184E-02  8.7246E-01  1.6905E-01  8.2953E-01  1.1170E-02  5.2847E-02  2.2146E+00  6.9874E-01
             3.2735E+00
 PARAMETER:  8.2384E-02 -1.0055E+00 -2.3283E+00 -3.6437E-02 -1.6775E+00 -8.6892E-02 -4.3945E+00 -2.8404E+00  8.9508E-01 -2.5848E-01
             1.2859E+00
 GRADIENT:  -1.0936E+00 -1.7720E-03 -4.0770E-01  7.8616E-01 -1.2025E-01 -4.6705E-01  1.6383E-03  8.6185E-03  1.3951E+00 -8.6024E-02
            -6.4026E-01

0ITERATION NO.:   69    OBJECTIVE VALUE:  -1246.22855047572        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     2058
 NPARAMETR:  9.8263E-01  3.3143E-01  8.8143E-02  8.7196E-01  1.6911E-01  8.2952E-01  1.0000E-02  2.2343E-02  2.2153E+00  6.9837E-01
             3.2740E+00
 PARAMETER:  8.2474E-02 -1.0043E+00 -2.3285E+00 -3.7001E-02 -1.6772E+00 -8.6914E-02 -5.5700E+00 -3.6646E+00  8.9527E-01 -2.5902E-01
             1.2861E+00
 GRADIENT:  -1.1428E-01 -1.3410E-01  9.5548E+00  1.3884E-01  3.3163E-01 -1.5188E-01  0.0000E+00  1.4093E-03 -7.6258E+00 -9.0954E-02
             4.1279E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2058
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9836E-03 -1.8711E-05  9.1662E-04 -2.1264E-02  1.5012E-02
 SE:             2.8153E-02  1.3136E-04  3.3298E-04  2.4569E-02  2.4176E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1560E-01  8.8673E-01  5.9100E-03  3.8677E-01  5.3465E-01

 ETASHRINKSD(%)  5.6841E+00  9.9560E+01  9.8884E+01  1.7691E+01  1.9006E+01
 ETASHRINKVR(%)  1.1045E+01  9.9998E+01  9.9988E+01  3.2252E+01  3.4400E+01
 EBVSHRINKSD(%)  6.0390E+00  9.9485E+01  9.8873E+01  1.0210E+01  2.0913E+01
 EBVSHRINKVR(%)  1.1713E+01  9.9997E+01  9.9987E+01  1.9378E+01  3.7453E+01
 RELATIVEINF(%)  8.1342E+01  4.5810E-04  4.0625E-03  5.8109E+01  6.5766E+00
 EPSSHRINKSD(%)  3.1371E+01
 EPSSHRINKVR(%)  5.2900E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1246.2285504757199     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -511.07772391198171     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1246.229       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  3.31E-01  8.82E-02  8.72E-01  1.69E-01  8.30E-01  1.00E-02  2.32E-02  2.22E+00  6.98E-01  3.27E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       65.260
Stop Time:
Sun Oct 24 02:27:10 CDT 2021
