Sat Sep 25 09:19:06 CDT 2021
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
$DATA ../../../../data/spa/A3/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25.1395980985696        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1681E+02  4.4685E+01  9.3944E+01 -5.7406E+01  1.4999E+02 -2.9358E+00 -1.3709E+02 -4.1728E+01 -2.1994E+02 -1.8274E+02
            -2.7064E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1174.55211380130        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.3666E-01  9.7507E-01  1.0998E+00  1.0769E+00  1.0927E+00  8.6536E-01  1.1301E+00  9.4453E-01  1.2880E+00  8.2508E-01
             3.7863E+00
 PARAMETER:  3.4568E-02  7.4750E-02  1.9511E-01  1.7408E-01  1.8867E-01 -4.4610E-02  2.2234E-01  4.2936E-02  3.5306E-01 -9.2272E-02
             1.4314E+00
 GRADIENT:  -4.2286E+01 -3.1304E+01 -2.1551E+01 -3.1981E+01  2.8423E+01 -2.3275E+01 -5.6279E-02  3.9927E+00  1.0411E+00  7.8139E+00
            -2.8149E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1181.94901383545        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.4579E-01  7.0601E-01  1.2279E+00  1.3316E+00  9.4624E-01  9.4498E-01  1.0312E+00  5.8491E-01  1.1334E+00  4.0805E-01
             4.0825E+00
 PARAMETER:  4.4261E-02 -2.4812E-01  3.0527E-01  3.8636E-01  4.4742E-02  4.3411E-02  1.3071E-01 -4.3630E-01  2.2525E-01 -7.9637E-01
             1.5067E+00
 GRADIENT:  -2.7816E+01  1.9122E+01  9.4650E+00  3.5522E+01 -2.3143E+01  6.5086E+00 -1.1902E+00  1.7375E+00 -3.6197E+00  1.5710E+00
             1.6339E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1184.32095654730        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.5505E-01  7.4266E-01  1.1544E+00  1.2635E+00  9.5934E-01  9.2036E-01  1.2036E+00  1.5560E-01  1.1289E+00  2.7848E-01
             3.9952E+00
 PARAMETER:  5.4010E-02 -1.9752E-01  2.4354E-01  3.3385E-01  5.8489E-02  1.7007E-02  2.8535E-01 -1.7605E+00  2.2126E-01 -1.1784E+00
             1.4851E+00
 GRADIENT:   3.1989E+00  7.4995E-01  4.4714E-01 -5.3234E+00 -5.4126E-01 -6.2077E-01 -2.8202E-01  1.4062E-01 -3.4985E+00  6.7270E-01
            -3.1773E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1185.44136010598        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.5021E-01  4.6147E-01  1.1061E+00  1.4333E+00  8.4543E-01  9.1780E-01  1.6738E+00  1.0000E-02  1.0549E+00  2.5801E-02
             4.0207E+00
 PARAMETER:  4.8931E-02 -6.7333E-01  2.0085E-01  4.5999E-01 -6.7908E-02  1.4229E-02  6.1510E-01 -5.2774E+00  1.5347E-01 -3.5574E+00
             1.4914E+00
 GRADIENT:  -1.1763E+00  4.6451E-01 -6.4077E-01  2.8583E-01  1.5596E+00 -5.2331E-01  3.2907E-01  0.0000E+00  6.3013E-01  4.7458E-03
             3.2286E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1185.65095038377        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.4806E-01  3.1821E-01  1.0701E+00  1.5111E+00  7.8702E-01  9.1843E-01  1.9733E+00  1.0000E-02  1.0190E+00  1.0000E-02
             4.0010E+00
 PARAMETER:  4.6659E-02 -1.0450E+00  1.6777E-01  5.1282E-01 -1.3950E-01  1.4912E-02  7.7972E-01 -8.8071E+00  1.1883E-01 -5.8983E+00
             1.4865E+00
 GRADIENT:   3.4806E-01  1.0131E-01 -8.2948E-01  7.1080E-01  1.0630E+00 -1.8149E-02 -7.5070E-03  0.0000E+00  3.2815E-02  0.0000E+00
            -3.1510E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1185.70944114034        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  9.4653E-01  2.5514E-01  1.0687E+00  1.5453E+00  7.6824E-01  9.1763E-01  2.1393E+00  1.0000E-02  1.0034E+00  1.0000E-02
             4.0008E+00
 PARAMETER:  4.5051E-02 -1.2660E+00  1.6642E-01  5.3524E-01 -1.6365E-01  1.4043E-02  8.6050E-01 -1.1129E+01  1.0340E-01 -7.4400E+00
             1.4865E+00
 GRADIENT:  -2.1636E-01  1.1217E-01  2.6562E-01  5.2978E-01 -4.7562E-01 -9.8626E-03 -3.3430E-02  0.0000E+00 -1.3084E-01  0.0000E+00
            -8.4717E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1185.77257452440        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      587
 NPARAMETR:  9.4690E-01  1.7132E-01  1.0960E+00  1.6059E+00  7.6118E-01  9.1834E-01  2.5366E+00  1.0000E-02  9.7715E-01  1.0000E-02
             4.0138E+00
 PARAMETER:  4.5439E-02 -1.6642E+00  1.9165E-01  5.7367E-01 -1.7289E-01  1.4816E-02  1.0308E+00 -1.5501E+01  7.6882E-02 -1.0357E+01
             1.4897E+00
 GRADIENT:   8.2023E-01  4.7776E-01  4.8101E-01  5.3340E+00 -1.1621E+00  7.5287E-02 -3.6375E-02  0.0000E+00 -3.5337E-01  0.0000E+00
             5.7557E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1185.81132515305        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  9.4437E-01  1.1174E-01  1.0747E+00  1.6307E+00  7.3711E-01  9.1642E-01  2.9376E+00  1.0000E-02  9.6829E-01  1.0000E-02
             4.0042E+00
 PARAMETER:  4.2762E-02 -2.0916E+00  1.7200E-01  5.8903E-01 -2.0502E-01  1.2716E-02  1.1776E+00 -2.0626E+01  6.7778E-02 -1.3756E+01
             1.4874E+00
 GRADIENT:  -9.6685E-01  5.9850E-02  1.3248E-01 -7.5851E-02 -4.4538E-01 -2.3361E-01 -6.7899E-03  0.0000E+00 -1.4857E-02  0.0000E+00
            -4.4575E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1185.82969229164        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  9.4345E-01  5.7615E-02  1.0837E+00  1.6630E+00  7.2891E-01  9.1616E-01  3.5915E+00  1.0000E-02  9.5498E-01  1.0000E-02
             4.0054E+00
 PARAMETER:  4.1788E-02 -2.7540E+00  1.8037E-01  6.0861E-01 -2.1620E-01  1.2440E-02  1.3786E+00 -2.8775E+01  5.3935E-02 -1.9161E+01
             1.4877E+00
 GRADIENT:  -3.8184E-01  4.0329E-02  2.8184E-01  8.9418E-01 -4.5968E-01 -8.6383E-02  9.8924E-03  0.0000E+00  1.9673E-01  0.0000E+00
            -1.2711E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1185.83570715726        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  9.4283E-01  2.7553E-02  1.0812E+00  1.6794E+00  7.2104E-01  9.1597E-01  4.3301E+00  1.0000E-02  9.4678E-01  1.0000E-02
             4.0055E+00
 PARAMETER:  4.1135E-02 -3.4917E+00  1.7804E-01  6.1844E-01 -2.2705E-01  1.2230E-02  1.5656E+00 -3.8097E+01  4.5311E-02 -2.5336E+01
             1.4877E+00
 GRADIENT:  -2.0568E-01  1.8346E-02  2.2711E-01  9.6705E-01 -4.5302E-01 -2.1300E-02  7.3286E-03  0.0000E+00  1.5994E-02  0.0000E+00
            -5.7318E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1185.83737109613        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  9.4254E-01  1.3311E-02  1.0793E+00  1.6865E+00  7.1726E-01  9.1582E-01  5.0161E+00  1.0000E-02  9.4295E-01  1.0000E-02
             4.0054E+00
 PARAMETER:  4.0824E-02 -4.2192E+00  1.7630E-01  6.2263E-01 -2.3232E-01  1.2061E-02  1.7126E+00 -4.7413E+01  4.1254E-02 -3.1497E+01
             1.4876E+00
 GRADIENT:  -5.3496E-03  2.5912E-03  4.5148E-02  1.9776E-01 -1.0456E-01  9.2546E-03  3.9557E-03  0.0000E+00  5.7965E-03  0.0000E+00
             2.4086E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1185.83745081538        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1465
 NPARAMETR:  9.4246E-01  1.0040E-02  1.0795E+00  1.6881E+00  7.1668E-01  9.1573E-01  5.2634E+00  1.0000E-02  9.4197E-01  1.0000E-02
             4.0052E+00
 PARAMETER:  4.0733E-02 -4.5012E+00  1.7647E-01  6.2363E-01 -2.3313E-01  1.1967E-02  1.7608E+00 -5.1039E+01  4.0220E-02 -3.3892E+01
             1.4876E+00
 GRADIENT:   4.1297E-03  8.7282E-04  1.7101E-03  1.5706E-03 -5.5580E-04 -1.8217E-04  2.8812E-03  0.0000E+00  9.9561E-03  0.0000E+00
             6.6729E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1185.83927928968        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     1637
 NPARAMETR:  9.4235E-01  1.0000E-02  1.0795E+00  1.6875E+00  7.1666E-01  9.1566E-01  7.4643E-01  1.0000E-02  9.4242E-01  1.0000E-02
             4.0047E+00
 PARAMETER:  4.0623E-02 -4.5087E+00  1.7646E-01  6.2328E-01 -2.3316E-01  1.1888E-02 -1.9246E-01 -5.1090E+01  4.0693E-02 -3.3926E+01
             1.4875E+00
 GRADIENT:  -1.4603E-01 -5.2283E-04  1.2436E-02 -5.6674E-01  5.4042E-02 -1.1708E-02  1.0587E-04  0.0000E+00  3.9354E-02  0.0000E+00
            -3.4499E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1185.83985355642        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1816
 NPARAMETR:  9.4261E-01  1.1920E-02  1.0776E+00  1.6872E+00  7.1639E-01  9.1590E-01  1.9576E-01  1.0000E-02  9.4316E-01  1.0000E-02
             4.0058E+00
 PARAMETER:  4.0899E-02 -4.3295E+00  1.7472E-01  6.2308E-01 -2.3353E-01  1.2155E-02 -1.5309E+00 -5.1090E+01  4.1477E-02 -3.3926E+01
             1.4878E+00
 GRADIENT:   2.2762E-01 -1.9440E-03 -9.7501E-02  3.9482E-01  2.8476E-02  3.1482E-02  1.0149E-05  0.0000E+00  1.4953E-02  0.0000E+00
             1.1088E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1185.84399629070        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1995
 NPARAMETR:  9.4279E-01  3.1412E-02  1.0765E+00  1.6757E+00  7.1994E-01  9.1611E-01  1.0000E-02  1.0000E-02  9.4914E-01  1.0000E-02
             4.0049E+00
 PARAMETER:  4.1087E-02 -3.3606E+00  1.7374E-01  6.1623E-01 -2.2858E-01  1.2384E-02 -2.4357E+01 -5.1090E+01  4.7801E-02 -3.3926E+01
             1.4875E+00
 GRADIENT:  -3.7685E-01  5.8960E-03  1.2004E-01  3.0752E-01 -4.5659E-01 -6.2255E-03  0.0000E+00  0.0000E+00 -2.5720E-01  0.0000E+00
            -1.8059E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1185.84834790040        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2173
 NPARAMETR:  9.4368E-01  6.5461E-02  1.0844E+00  1.6574E+00  7.3172E-01  9.1653E-01  1.0000E-02  1.0000E-02  9.6142E-01  1.0000E-02
             4.0061E+00
 PARAMETER:  4.2032E-02 -2.6263E+00  1.8101E-01  6.0525E-01 -2.1235E-01  1.2840E-02 -5.4675E+01 -5.1090E+01  6.0661E-02 -3.3926E+01
             1.4878E+00
 GRADIENT:  -7.0853E-02  7.1381E-03  3.8648E-03  1.6426E-01 -1.8215E-01 -3.0013E-02  0.0000E+00  0.0000E+00 -8.9942E-04  0.0000E+00
            -5.8621E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1185.84871182441        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2348
 NPARAMETR:  9.4392E-01  7.4599E-02  1.0898E+00  1.6527E+00  7.3651E-01  9.1672E-01  1.0000E-02  1.0000E-02  9.6409E-01  1.0000E-02
             4.0067E+00
 PARAMETER:  4.2281E-02 -2.4956E+00  1.8598E-01  6.0239E-01 -2.0583E-01  1.3046E-02 -6.1831E+01 -5.1090E+01  6.3428E-02 -3.3926E+01
             1.4880E+00
 GRADIENT:   1.3652E-02 -2.9762E-03 -4.1381E-02 -1.4158E-01  2.2213E-02  3.7544E-04  0.0000E+00  0.0000E+00  1.9347E-02  0.0000E+00
             3.7575E-03

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1185.84876996230        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2510
 NPARAMETR:  9.4385E-01  7.1174E-02  1.0916E+00  1.6550E+00  7.3651E-01  9.1667E-01  1.0000E-02  1.0000E-02  9.6251E-01  1.0000E-02
             4.0068E+00
 PARAMETER:  4.2211E-02 -2.5426E+00  1.8763E-01  6.0378E-01 -2.0584E-01  1.2995E-02 -5.7842E+01 -5.1090E+01  6.1790E-02 -3.3926E+01
             1.4880E+00
 GRADIENT:   2.6105E-03 -7.7753E-04 -2.0699E-03 -1.9843E-02  8.4730E-03  1.3896E-03  0.0000E+00  0.0000E+00  1.6890E-03  0.0000E+00
             3.2123E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2510
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9222E-04 -2.3068E-05  1.0790E-04 -1.4978E-02 -2.6514E-05
 SE:             2.7996E-02  9.3948E-06  1.0434E-04  2.5734E-02  1.7063E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9452E-01  1.4074E-02  3.0110E-01  5.6056E-01  8.7652E-01

 ETASHRINKSD(%)  6.2112E+00  9.9969E+01  9.9650E+01  1.3787E+01  9.9428E+01
 ETASHRINKVR(%)  1.2037E+01  1.0000E+02  9.9999E+01  2.5674E+01  9.9997E+01
 EBVSHRINKSD(%)  5.8365E+00  9.9971E+01  9.9599E+01  1.3538E+01  9.9418E+01
 EBVSHRINKVR(%)  1.1332E+01  1.0000E+02  9.9998E+01  2.5243E+01  9.9997E+01
 RELATIVEINF(%)  8.0667E+01  3.9559E-07  1.1525E-04  5.8878E+00  1.3849E-04
 EPSSHRINKSD(%)  2.0670E+01
 EPSSHRINKVR(%)  3.7068E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1185.8487699623038     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -450.69794339856560     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1185.849       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.44E-01  7.12E-02  1.09E+00  1.65E+00  7.37E-01  9.17E-01  1.00E-02  1.00E-02  9.63E-01  1.00E-02  4.01E+00
 


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
+        1.32E+03
 
 TH 2
+       -8.59E+01  2.03E+02
 
 TH 3
+       -6.42E+00  7.02E+01  1.31E+02
 
 TH 4
+       -7.13E+01  2.15E+02  9.07E+00  3.23E+02
 
 TH 5
+        2.68E+01 -2.28E+02 -2.79E+02 -1.11E+02  6.61E+02
 
 TH 6
+        2.64E+00 -1.18E+01  9.05E+00 -1.60E+01 -1.40E+00  1.86E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.60E+01 -5.31E+01  1.20E+01 -5.99E+00  2.14E+01  1.20E+01  0.00E+00  0.00E+00  1.12E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.29E+01 -8.74E+00 -7.55E-01 -8.44E+00  3.87E+00  3.89E+00  0.00E+00  0.00E+00  7.69E+00  0.00E+00  2.75E+01
 
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
 #CPUT: Total CPU Time in Seconds,       39.268
Stop Time:
Sat Sep 25 09:19:46 CDT 2021
