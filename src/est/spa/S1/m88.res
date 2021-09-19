Sat Sep 18 11:21:01 CDT 2021
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
$DATA ../../../../data/spa/S1/dat88.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1645.49514836112        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2546E+01 -8.1782E+01 -4.6077E+01 -5.0092E+01  7.4157E+01  9.8315E+00 -1.0530E+01 -3.9621E+00  6.4165E-01  2.3547E+01
            -5.7238E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1657.45893625268        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0062E+00  1.1808E+00  1.1325E+00  9.7824E-01  1.0207E+00  9.6724E-01  1.1446E+00  1.1653E+00  9.7274E-01  6.1227E-01
             1.1778E+00
 PARAMETER:  1.0614E-01  2.6618E-01  2.2444E-01  7.7996E-02  1.2052E-01  6.6695E-02  2.3509E-01  2.5299E-01  7.2362E-02 -3.9058E-01
             2.6369E-01
 GRADIENT:   3.2483E+01  6.0281E+01  2.6042E+01  4.2232E+01 -1.9201E+00 -3.5493E+00  1.4472E+00 -1.8622E+01 -3.3905E+00 -9.0982E+00
             7.8749E-02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1660.25474384316        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0027E+00  1.0837E+00  9.1177E-01  1.0197E+00  8.8584E-01  9.8029E-01  1.2655E+00  1.3597E+00  9.0135E-01  4.3050E-01
             1.1501E+00
 PARAMETER:  1.0273E-01  1.8039E-01  7.6329E-03  1.1952E-01 -2.1220E-02  8.0098E-02  3.3544E-01  4.0726E-01 -3.8662E-03 -7.4281E-01
             2.3989E-01
 GRADIENT:   2.3313E+01  3.9740E+01  2.4398E+00  4.9911E+01 -1.9265E+01  1.2578E+00  4.6208E+00  3.6425E+00  1.2394E+00 -3.8554E-01
             6.9831E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1661.33622457171        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.9294E-01  1.1001E+00  9.2188E-01  9.6969E-01  9.1933E-01  9.7753E-01  1.1853E+00  1.2551E+00  9.3205E-01  5.2935E-01
             1.1464E+00
 PARAMETER:  9.2918E-02  1.9537E-01  1.8657E-02  6.9221E-02  1.5888E-02  7.7271E-02  2.6996E-01  3.2720E-01  2.9627E-02 -5.3610E-01
             2.3665E-01
 GRADIENT:   1.1606E+00 -2.1928E+00  1.1885E+00 -5.3727E+00 -1.4654E+00  1.2896E-01  5.5109E-01 -4.2263E-01  3.8100E-01 -1.8083E-01
             2.6077E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1661.68270973112        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      390
 NPARAMETR:  1.0065E+00  1.1138E+00  1.0358E+00  9.7750E-01  9.7104E-01  9.8384E-01  1.1500E+00  1.3902E+00  9.5114E-01  6.0326E-01
             1.1483E+00
 PARAMETER:  1.0652E-01  2.0778E-01  1.3518E-01  7.7242E-02  7.0611E-02  8.3710E-02  2.3975E-01  4.2946E-01  4.9910E-02 -4.0540E-01
             2.3826E-01
 GRADIENT:   1.2970E+00 -8.3680E-01 -1.3934E-02 -3.0288E-01  1.4266E-01  3.3463E-02  1.7998E-01 -4.7364E-02  8.0847E-03  1.4015E-01
             2.0532E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1661.84939837780        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      575
 NPARAMETR:  1.0092E+00  1.4011E+00  6.4974E-01  7.7972E-01  9.2685E-01  9.8527E-01  1.0055E+00  1.1578E+00  1.0629E+00  5.1225E-01
             1.1417E+00
 PARAMETER:  1.0920E-01  4.3726E-01 -3.3119E-01 -1.4882E-01  2.4031E-02  8.5164E-02  1.0552E-01  2.4650E-01  1.6096E-01 -5.6893E-01
             2.3253E-01
 GRADIENT:   1.8721E+00  1.0221E+01  1.7911E+00  5.6107E+00 -8.0339E+00 -5.5696E-01  5.7370E-01  2.3057E-01  3.3039E-01  8.0713E-01
            -9.7465E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1662.32520400592        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  1.0069E+00  1.6026E+00  3.4846E-01  6.1508E-01  8.4737E-01  9.8710E-01  9.2147E-01  8.0928E-01  1.1696E+00  3.5306E-01
             1.1391E+00
 PARAMETER:  1.0692E-01  5.7161E-01 -9.5423E-01 -3.8600E-01 -6.5620E-02  8.7013E-02  1.8213E-02 -1.1161E-01  2.5670E-01 -9.4112E-01
             2.3023E-01
 GRADIENT:  -8.3994E-01  2.4093E+01  1.0773E+01 -6.3937E+00 -3.3790E+01  7.8986E-02  2.6289E+00  6.0093E-01  2.7066E+00  1.3977E+00
             2.3568E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.03136248636        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      932
 NPARAMETR:  9.9019E-01  1.7744E+00  1.5249E-01  4.6431E-01  7.8915E-01  9.6678E-01  8.2737E-01  5.7958E-01  1.3154E+00  1.9466E-01
             1.0300E+00
 PARAMETER:  9.0146E-02  6.7347E-01 -1.7806E+00 -6.6720E-01 -1.3680E-01  6.6217E-02 -8.9500E-02 -4.4545E-01  3.7416E-01 -1.5365E+00
             1.2958E-01
 GRADIENT:  -1.4063E+01  5.5138E+01  8.9063E+00 -6.0600E+00 -5.2406E+01 -8.1057E+00 -5.3890E-01 -3.4413E+00 -3.0080E+00  4.2538E-01
            -1.5925E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1666.77639958547        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1109
 NPARAMETR:  9.9706E-01  1.7879E+00  1.5056E-01  4.4844E-01  8.1879E-01  9.8874E-01  7.9982E-01  7.0004E-01  1.3675E+00  1.6306E-01
             1.0714E+00
 PARAMETER:  9.7054E-02  6.8101E-01 -1.7934E+00 -7.0199E-01 -9.9922E-02  8.8674E-02 -1.2336E-01 -2.5662E-01  4.1296E-01 -1.7136E+00
             1.6899E-01
 GRADIENT:  -1.0837E-02  5.1804E-01 -1.2539E+00  7.6924E+00  8.9365E+00  7.7540E-01 -7.8618E+00 -3.4801E+00 -2.7211E+00 -1.7111E-01
             4.4923E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1666.87912469736        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  9.9628E-01  1.7879E+00  1.5055E-01  4.4809E-01  8.1686E-01  9.8579E-01  8.1828E-01  7.0266E-01  1.3681E+00  1.9674E-01
             1.0715E+00
 PARAMETER:  9.6277E-02  6.8106E-01 -1.7935E+00 -7.0277E-01 -1.0229E-01  8.5683E-02 -1.0056E-01 -2.5288E-01  4.1341E-01 -1.5259E+00
             1.6909E-01
 GRADIENT:  -1.6350E+00  4.3091E+00 -2.2914E-02  3.6323E+00 -2.6193E+00 -3.6947E-01  7.6322E-01 -3.2968E+00 -1.0941E+00  2.6255E-01
             2.8159E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1666.91922966081        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1468
 NPARAMETR:  9.9620E-01  1.7879E+00  1.5055E-01  4.4808E-01  8.1607E-01  9.9112E-01  8.1862E-01  7.0741E-01  1.3681E+00  1.5719E-01
             1.0715E+00
 PARAMETER:  9.6194E-02  6.8106E-01 -1.7935E+00 -7.0277E-01 -1.0326E-01  9.1082E-02 -1.0014E-01 -2.4614E-01  4.1341E-01 -1.7503E+00
             1.6909E-01
 GRADIENT:  -1.6499E+00  6.6833E+00 -8.1237E-02  3.7574E+00 -5.5151E-01  1.7018E+00 -3.1579E-01 -3.3963E+00 -2.1056E+00 -8.8521E-03
             1.7484E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1666.92878378939        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1661
 NPARAMETR:  9.9623E-01  1.7873E+00  1.5069E-01  4.4793E-01  8.1625E-01  9.8574E-01  8.1960E-01  7.0804E-01  1.3678E+00  1.5694E-01
             1.0716E+00
 PARAMETER:  9.6221E-02  6.8072E-01 -1.7926E+00 -7.0312E-01 -1.0303E-01  8.5632E-02 -9.8935E-02 -2.4526E-01  4.1321E-01 -1.7519E+00
             1.6918E-01
 GRADIENT:  -1.6437E+00  5.4808E+00 -5.0532E-02  3.3489E+00  3.9353E-02 -4.2006E-01  6.4916E-02 -3.3866E+00 -2.1607E+00 -5.6712E-03
             1.7197E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1668.03824047133        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     1776
 NPARAMETR:  9.9151E-01  1.7736E+00  1.5284E-01  4.4606E-01  8.1761E-01  9.8500E-01  8.1931E-01  1.1576E+00  1.3817E+00  1.3927E-01
             1.0533E+00
 PARAMETER:  9.1479E-02  6.7300E-01 -1.7784E+00 -7.0729E-01 -1.0137E-01  8.4890E-02 -9.9296E-02  2.4636E-01  4.2328E-01 -1.8713E+00
             1.5190E-01
 GRADIENT:   2.6767E+01  4.6980E+01 -4.5799E+00  6.7124E+00 -5.0451E+00  1.9371E+00  7.3876E-01  1.3532E+00 -7.5321E-01  1.5426E-01
             9.2362E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1668.54882667195        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1935
 NPARAMETR:  9.9891E-01  1.7891E+00  1.6617E-01  4.5293E-01  8.2962E-01  9.9268E-01  8.2453E-01  1.2621E+00  1.4047E+00  1.3144E-01
             1.0221E+00
 PARAMETER:  9.8913E-02  6.8173E-01 -1.6948E+00 -6.9203E-01 -8.6788E-02  9.2654E-02 -9.2942E-02  3.3280E-01  4.3985E-01 -1.9292E+00
             1.2187E-01
 GRADIENT:   3.5810E+00  5.5118E+00  3.5489E-01  1.2667E+00 -5.8613E+00  6.7274E-01  5.7128E-01 -3.6224E-01  9.7581E-01 -6.0431E-02
            -2.3800E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1668.60175777359        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2111
 NPARAMETR:  9.9802E-01  1.7908E+00  1.6993E-01  4.5289E-01  8.3657E-01  9.9092E-01  8.2079E-01  1.3129E+00  1.3930E+00  1.7140E-01
             1.0293E+00
 PARAMETER:  9.8020E-02  6.8269E-01 -1.6723E+00 -6.9211E-01 -7.8448E-02  9.0880E-02 -9.7492E-02  3.7222E-01  4.3148E-01 -1.6638E+00
             1.2885E-01
 GRADIENT:   3.9241E-01 -7.5818E-01  2.1238E-01  9.5291E-01 -1.6255E+00  5.8869E-03  2.0310E-01  2.2787E-02 -1.6834E-01  4.5030E-03
             4.0557E-01

0ITERATION NO.:   72    OBJECTIVE VALUE:  -1668.60213182907        NO. OF FUNC. EVALS.:  68
 CUMULATIVE NO. OF FUNC. EVALS.:     2179
 NPARAMETR:  9.9802E-01  1.7908E+00  1.6996E-01  4.5286E-01  8.3667E-01  9.9091E-01  8.2072E-01  1.3130E+00  1.3931E+00  1.7154E-01
             1.0293E+00
 PARAMETER:  9.8013E-02  6.8276E-01 -1.6723E+00 -6.9225E-01 -7.8330E-02  9.0872E-02 -9.7587E-02  3.7234E-01  4.3152E-01 -1.6613E+00
             1.2887E-01
 GRADIENT:  -4.0215E+05  5.8898E+04 -2.4025E+04 -5.8106E+04 -2.0107E+05  1.5953E-03 -4.0215E+05  5.3972E+04 -4.6604E+04  4.2882E-03
             1.5603E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         4.4         3.3         3.3         3.3         2.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2179
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2486E-04 -8.2394E-03 -9.5273E-03  1.3672E-02 -1.5806E-02
 SE:             2.9844E-02  2.8357E-02  1.2140E-02  2.4764E-02  6.1866E-03
 N:                     100         100         100         100         100

 P VAL.:         9.9666E-01  7.7139E-01  4.3259E-01  5.8088E-01  1.0622E-02

 ETASHRINKSD(%)  1.8306E-02  5.0018E+00  5.9328E+01  1.7037E+01  7.9274E+01
 ETASHRINKVR(%)  3.6608E-02  9.7534E+00  8.3458E+01  3.1172E+01  9.5704E+01
 EBVSHRINKSD(%)  4.7179E-01  5.2413E+00  5.9278E+01  1.5774E+01  8.0292E+01
 EBVSHRINKVR(%)  9.4135E-01  1.0208E+01  8.3417E+01  2.9060E+01  9.6116E+01
 RELATIVEINF(%)  9.8415E+01  1.9023E+01  6.0309E+00  1.1310E+01  6.5209E-01
 EPSSHRINKSD(%)  4.4427E+01
 EPSSHRINKVR(%)  6.9116E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1668.6021318290673     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -933.45130526532910     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.28
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.22
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1668.602       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  1.79E+00  1.70E-01  4.53E-01  8.37E-01  9.91E-01  8.21E-01  1.31E+00  1.39E+00  1.72E-01  1.03E+00
 


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
+        1.01E+09
 
 TH 2
+        2.89E+02  6.72E+06
 
 TH 3
+       -1.33E+03  5.28E+04  1.24E+08
 
 TH 4
+        1.53E+05 -1.23E+04  5.19E+04  1.02E+08
 
 TH 5
+       -4.93E+04  3.47E+03 -1.91E+04  1.84E+05  1.44E+09
 
 TH 6
+       -2.45E+03  1.99E+02 -8.60E+02 -7.88E+02 -2.93E+03  2.00E+02
 
 TH 7
+       -8.02E+03  6.65E+02 -2.90E+03 -2.56E+03 -9.54E+03 -2.98E+03  1.49E+09
 
 TH 8
+        2.35E+05 -1.91E+04  8.22E+04 -3.12E+04 -2.46E+08  5.00E+02  1.64E+03  4.20E+07
 
 TH 9
+       -5.95E+02 -3.87E+03  1.66E+04  2.55E+04 -8.24E+03 -4.08E+02 -1.33E+03  3.90E+04  2.78E+07
 
 TH10
+        1.18E+04 -9.72E+02  4.19E+03  3.74E+03  1.40E+04 -1.97E-01  1.44E+04 -2.41E+03  1.98E+03  1.82E+01
 
 TH11
+        2.94E+04 -2.41E+03  1.03E+04  9.36E+03  3.50E+04  1.85E+03  6.04E+03 -6.00E+03  4.91E+03 -8.90E+03  5.71E+08
 
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
 #CPUT: Total CPU Time in Seconds,       32.564
Stop Time:
Sat Sep 18 11:21:35 CDT 2021
