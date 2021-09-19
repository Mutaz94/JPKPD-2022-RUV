Sat Sep 18 11:34:11 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat23.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1688.02945519678        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0775E+00  3.1015E+00 -4.2005E+01  5.5606E+01 -3.9471E+00  1.6614E+01  4.8935E+00  1.5615E+01  2.6839E+01  1.5746E+01
            -5.5478E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1697.87631267429        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0128E+00  1.1498E+00  1.3370E+00  9.4392E-01  1.2038E+00  9.3653E-01  9.7431E-01  8.6046E-01  8.2941E-01  9.0063E-01
             1.2013E+00
 PARAMETER:  1.1269E-01  2.3955E-01  3.9047E-01  4.2285E-02  2.8549E-01  3.4427E-02  7.3972E-02 -5.0292E-02 -8.7042E-02 -4.6649E-03
             2.8336E-01
 GRADIENT:   1.4416E+01  6.5904E+01 -9.6085E+00  8.7798E+01  3.5421E+01 -7.8926E+00 -4.7380E+00  5.8894E-01 -1.1900E+01 -2.3207E+01
             9.6125E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1702.50309836548        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0144E+00  1.0117E+00  1.5177E+00  1.0208E+00  1.2022E+00  9.5184E-01  6.1485E-01  4.4427E-01  9.9503E-01  1.0796E+00
             1.1758E+00
 PARAMETER:  1.1426E-01  1.1167E-01  5.1718E-01  1.2060E-01  2.8417E-01  5.0638E-02 -3.8638E-01 -7.1131E-01  9.5016E-02  1.7659E-01
             2.6191E-01
 GRADIENT:   2.4968E+01  4.2164E+01 -8.7819E+00  7.3522E+01  1.0908E+01 -7.1596E-01 -1.7372E+00 -5.5349E-02 -2.7270E+00 -6.8667E+00
             7.9120E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1704.00500141610        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0032E+00  9.9005E-01  1.7234E+00  9.9671E-01  1.2425E+00  9.5441E-01  5.1150E-01  5.2710E-01  1.0208E+00  1.1666E+00
             1.1501E+00
 PARAMETER:  1.0324E-01  8.9996E-02  6.4432E-01  9.6704E-02  3.1709E-01  5.3342E-02 -5.7040E-01 -5.4036E-01  1.2055E-01  2.5412E-01
             2.3987E-01
 GRADIENT:  -5.8086E-01  2.0915E-01 -4.9181E-01  2.7781E+00  1.7841E+00  6.0933E-01  2.1978E-01 -1.6230E-01  6.2325E-01  7.1612E-01
            -1.8150E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1704.09012735091        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0035E+00  1.0111E+00  1.8052E+00  9.8116E-01  1.2609E+00  9.5379E-01  4.4925E-01  8.8744E-01  1.0480E+00  1.1625E+00
             1.1454E+00
 PARAMETER:  1.0348E-01  1.1105E-01  6.9065E-01  8.0977E-02  3.3185E-01  5.2686E-02 -7.0017E-01 -1.9416E-02  1.4688E-01  2.5061E-01
             2.3573E-01
 GRADIENT:   2.4813E-01  6.0174E-02  7.4205E-01 -1.6247E-01 -9.0330E-01  3.3539E-01  4.7933E-02 -9.9045E-02  9.8521E-01  1.4166E-01
            -7.1128E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1704.45630565888        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  1.0196E+00  1.0159E+00  1.7965E+00  9.8346E-01  1.2640E+00  9.6173E-01  5.0741E-01  9.0682E-01  1.0336E+00  1.1557E+00
             1.1482E+00
 PARAMETER:  1.1943E-01  1.1574E-01  6.8584E-01  8.3325E-02  3.3425E-01  6.0982E-02 -5.7844E-01  2.1889E-03  1.3307E-01  2.4474E-01
             2.3819E-01
 GRADIENT:   5.8115E-01 -4.1765E-01 -2.6733E-01  1.4328E-01  6.4707E-01 -3.4245E-01 -2.0346E-01 -7.6986E-03 -5.5487E-01 -5.1095E-01
            -1.0923E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1704.54095251042        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      617
 NPARAMETR:  1.0198E+00  1.1697E+00  1.6665E+00  8.7969E-01  1.2934E+00  9.6270E-01  5.0693E-01  8.2477E-01  1.1305E+00  1.1647E+00
             1.1507E+00
 PARAMETER:  1.1964E-01  2.5676E-01  6.1074E-01 -2.8187E-02  3.5730E-01  6.1985E-02 -5.7938E-01 -9.2657E-02  2.2268E-01  2.5248E-01
             2.4040E-01
 GRADIENT:  -1.3273E+00 -1.4151E+00 -2.7451E-01 -8.1291E-01  1.4464E+00 -4.2092E-01 -3.0532E-01 -4.1292E-02 -4.2764E-01 -2.5944E-01
             2.4473E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1704.59754424324        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      792
 NPARAMETR:  1.0213E+00  1.2972E+00  1.5280E+00  7.9766E-01  1.3098E+00  9.6514E-01  5.3033E-01  7.3612E-01  1.2131E+00  1.1651E+00
             1.1508E+00
 PARAMETER:  1.2107E-01  3.6022E-01  5.2394E-01 -1.2607E-01  3.6987E-01  6.4513E-02 -5.3426E-01 -2.0636E-01  2.9316E-01  2.5277E-01
             2.4043E-01
 GRADIENT:   1.3391E-01  4.2378E+00  3.8778E-02  3.8496E+00 -7.4819E-01  1.5850E-01 -2.4022E-01 -1.7214E-02 -4.3059E-01  2.6578E-02
            -3.3280E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1704.63677373005        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  1.0224E+00  1.4445E+00  1.3836E+00  6.9665E-01  1.3367E+00  9.6596E-01  5.3467E-01  6.7240E-01  1.3400E+00  1.1673E+00
             1.1524E+00
 PARAMETER:  1.2213E-01  4.6777E-01  4.2466E-01 -2.6147E-01  3.9021E-01  6.5362E-02 -5.2611E-01 -2.9691E-01  3.9267E-01  2.5469E-01
             2.4183E-01
 GRADIENT:   8.1976E-01  1.7801E+00  2.8520E-02  1.1109E+00 -7.0933E-01  1.0816E-01 -1.0830E-01  4.6883E-02 -9.0164E-02  7.0580E-02
            -1.3970E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1704.64192276861        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1142
 NPARAMETR:  1.0223E+00  1.4878E+00  1.3297E+00  6.6634E-01  1.3442E+00  9.6597E-01  5.3977E-01  5.4550E-01  1.3800E+00  1.1680E+00
             1.1542E+00
 PARAMETER:  1.2202E-01  4.9730E-01  3.8493E-01 -3.0595E-01  3.9580E-01  6.5377E-02 -5.1662E-01 -5.0606E-01  4.2206E-01  2.5532E-01
             2.4341E-01
 GRADIENT:   5.1118E-02 -5.7640E-01 -8.7604E-02 -3.2444E-01  2.2143E-01  2.0237E-02  1.7058E-02  9.7933E-03 -2.0535E-02  6.2579E-03
             2.2469E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1704.64224106014        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1318
 NPARAMETR:  1.0222E+00  1.4895E+00  1.3284E+00  6.6540E-01  1.3444E+00  9.6591E-01  5.3934E-01  5.0946E-01  1.3821E+00  1.1686E+00
             1.1547E+00
 PARAMETER:  1.2200E-01  4.9846E-01  3.8397E-01 -3.0737E-01  3.9595E-01  6.5310E-02 -5.1741E-01 -5.7441E-01  4.2359E-01  2.5577E-01
             2.4388E-01
 GRADIENT:  -3.7983E-03 -1.0073E-02  1.8509E-02 -7.9299E-03 -3.3269E-02 -2.8743E-03  5.7732E-03 -1.2090E-04  2.6148E-02  6.1442E-03
             7.7518E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1704.64225983327        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1498
 NPARAMETR:  1.0223E+00  1.4889E+00  1.3282E+00  6.6582E-01  1.3441E+00  9.6592E-01  5.3962E-01  5.0987E-01  1.3811E+00  1.1684E+00
             1.1547E+00
 PARAMETER:  1.2201E-01  4.9804E-01  3.8386E-01 -3.0673E-01  3.9576E-01  6.5326E-02 -5.1689E-01 -5.7361E-01  4.2288E-01  2.5562E-01
             2.4383E-01
 GRADIENT:   2.2824E-02 -3.2090E-02  1.5206E-03 -1.9822E-02 -7.3678E-03  2.7500E-03  4.0876E-03  8.2405E-05  4.5811E-03  6.8436E-04
            -1.4766E-03

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1704.64225983327        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1520
 NPARAMETR:  1.0223E+00  1.4889E+00  1.3282E+00  6.6582E-01  1.3441E+00  9.6592E-01  5.3962E-01  5.0987E-01  1.3811E+00  1.1684E+00
             1.1547E+00
 PARAMETER:  1.2201E-01  4.9804E-01  3.8386E-01 -3.0673E-01  3.9576E-01  6.5326E-02 -5.1689E-01 -5.7361E-01  4.2288E-01  2.5562E-01
             2.4383E-01
 GRADIENT:   2.2824E-02 -3.2090E-02  1.5206E-03 -1.9822E-02 -7.3678E-03  2.7500E-03  4.0876E-03  8.2405E-05  4.5811E-03  6.8436E-04
            -1.4766E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1520
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.0325E-04 -3.2866E-02 -4.2658E-03  1.0794E-02 -3.5169E-02
 SE:             2.9798E-02  1.5058E-02  3.2765E-03  2.4967E-02  2.3435E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8653E-01  2.9066E-02  1.9295E-01  6.6550E-01  1.3343E-01

 ETASHRINKSD(%)  1.7323E-01  4.9553E+01  8.9023E+01  1.6358E+01  2.1490E+01
 ETASHRINKVR(%)  3.4616E-01  7.4551E+01  9.8795E+01  3.0041E+01  3.8362E+01
 EBVSHRINKSD(%)  5.6884E-01  4.9218E+01  8.9487E+01  1.6503E+01  1.9639E+01
 EBVSHRINKVR(%)  1.1344E+00  7.4212E+01  9.8895E+01  3.0283E+01  3.5421E+01
 RELATIVEINF(%)  9.8689E+01  1.6091E+00  3.9453E-01  4.7702E+00  2.2342E+01
 EPSSHRINKSD(%)  3.9979E+01
 EPSSHRINKVR(%)  6.3975E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1704.6422598332747     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -969.49143326953651     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1704.642       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.49E+00  1.33E+00  6.66E-01  1.34E+00  9.66E-01  5.40E-01  5.10E-01  1.38E+00  1.17E+00  1.15E+00
 


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
+        1.12E+03
 
 TH 2
+       -1.51E+01  4.63E+02
 
 TH 3
+        3.09E+00  2.76E+01  2.55E+01
 
 TH 4
+       -1.78E+01  5.67E+02 -1.61E+01  9.23E+02
 
 TH 5
+       -4.74E+00 -1.16E+02 -6.24E+01  4.73E+00  2.95E+02
 
 TH 6
+       -4.47E+00 -1.45E+00  6.16E-01 -3.06E+00  2.35E-02  2.03E+02
 
 TH 7
+        1.88E+00 -4.55E+01  1.05E+01 -3.04E+01 -1.10E+01  1.46E+00  6.27E+01
 
 TH 8
+        9.96E-01 -8.18E-01 -1.11E+00 -2.17E-01  6.48E-01  1.01E+00  4.93E-02  8.79E-02
 
 TH 9
+        1.49E+00 -2.89E+01 -1.34E+00  4.07E+01  1.36E+00 -1.12E+00  3.56E+01 -9.18E-02  5.40E+01
 
 TH10
+       -3.79E-01  1.32E+00 -2.76E+00 -3.81E+00 -4.09E+01 -4.94E-01  6.42E+00  8.86E-01  3.71E-01  6.46E+01
 
 TH11
+       -8.92E+00 -2.83E+01 -1.33E+01 -1.26E+01  3.23E+00  3.07E+00  6.76E+00  3.12E+00  5.36E+00  2.21E+01  1.79E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.044
Stop Time:
Sat Sep 18 11:34:36 CDT 2021
