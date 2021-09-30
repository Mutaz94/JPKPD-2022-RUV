Wed Sep 29 08:19:52 CDT 2021
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
$DATA ../../../../data/int/D/dat25.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12663.1015807402        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5643E+02  9.8046E+01 -8.3082E+01 -7.6000E+01  2.1770E+02 -1.2504E+03 -5.7744E+02 -3.3848E+01 -9.6541E+02 -3.9642E+02
            -2.9733E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1160.65147241937        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.6092E+00  1.3934E+00  9.8931E-01  2.8835E+00  1.1508E+00  5.4935E+00  5.0834E+00  9.8256E-01  5.4192E+00  2.3184E+00
             1.0636E+01
 PARAMETER:  5.7573E-01  4.3172E-01  8.9253E-02  1.1590E+00  2.4050E-01  1.8036E+00  1.7260E+00  8.2403E-02  1.7899E+00  9.4086E-01
             2.4642E+00
 GRADIENT:   4.0807E+01 -1.7665E+01 -5.0295E+01  7.5421E+01 -6.9230E+00  2.1875E+02  1.2298E+02  3.2401E+00  1.2240E+02  6.2951E+01
             7.3744E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1289.15273218597        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      220
 NPARAMETR:  1.2456E+00  6.1041E+00  1.7637E+00  6.2518E+00  6.1006E+00  1.8349E+00  2.2323E+00  8.2392E-01  3.5130E+01  1.6950E+00
             1.0940E+01
 PARAMETER:  3.1963E-01  1.9090E+00  6.6739E-01  1.9329E+00  1.9084E+00  7.0698E-01  9.0303E-01 -9.3678E-02  3.6590E+00  6.2766E-01
             2.4924E+00
 GRADIENT:  -1.7715E+01  8.8285E+01  3.8614E+00  1.3181E+01  4.8898E+01  4.4072E+00 -9.4736E+01 -2.0653E+00 -2.0606E+01  1.7877E+01
             6.4169E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1501.13223781465        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  1.0597E+00  4.2787E+00  1.4588E+00  5.8889E+00  3.3661E+00  1.8350E+00  2.7282E+00  4.9044E+00  3.5189E+01  1.4505E+00
             9.0375E+00
 PARAMETER:  1.5796E-01  1.5536E+00  4.7759E-01  1.8731E+00  1.3138E+00  7.0705E-01  1.1036E+00  1.6901E+00  3.6607E+00  4.7191E-01
             2.3014E+00
 GRADIENT:  -3.4845E+00 -5.0496E+01  1.5747E+01  1.2377E+01 -6.2448E+00  5.9110E+01  3.1738E+01 -3.3191E+01 -9.8867E+00  3.9975E+01
             6.5707E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1628.66004895981        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  8.8737E-01  4.3116E+00  1.5145E+00  6.1572E+00  3.5170E+00  1.6364E+00  3.6345E+00  1.0539E+01  6.3732E+01  1.3196E+00
             6.2762E+00
 PARAMETER: -1.9497E-02  1.5613E+00  5.1508E-01  1.9176E+00  1.3576E+00  5.9249E-01  1.3905E+00  2.4551E+00  4.2547E+00  3.7732E-01
             1.9368E+00
 GRADIENT:  -2.3040E+01 -2.4111E+01  1.7189E+01  3.5474E+00  4.0959E+01  4.8611E+01  9.5595E+01 -9.3428E+01 -3.4777E+01  3.2761E+01
             1.8868E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1636.47342461482        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  9.0831E-01  4.3031E+00  1.5078E+00  6.1396E+00  3.5063E+00  1.6373E+00  2.8840E+00  1.0509E+01  7.0180E+01  1.0302E+00
             6.1788E+00
 PARAMETER:  3.8300E-03  1.5593E+00  5.1062E-01  1.9148E+00  1.3546E+00  5.9305E-01  1.1592E+00  2.4522E+00  4.3511E+00  1.2979E-01
             1.9211E+00
 GRADIENT:  -1.2905E+01 -1.5118E+01  1.8067E+01  3.3169E+00  5.1204E+01  4.8837E+01  4.5778E+01 -7.0117E+01 -1.8444E+01  1.9411E+01
             1.1814E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1704.37226049702        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      982            RESET HESSIAN, TYPE II
 NPARAMETR:  9.2214E-01  4.4096E+00  7.1994E-01  2.6847E-01  3.3067E+00  1.3652E+00  2.6514E+00  1.0647E+01  5.2766E+01  7.3072E-01
             5.9382E+00
 PARAMETER:  1.8946E-02  1.5838E+00 -2.2858E-01 -1.2150E+00  1.2960E+00  4.1129E-01  1.0751E+00  2.4653E+00  4.0659E+00 -2.1372E-01
             1.8814E+00
 GRADIENT:  -9.6407E+00  2.5290E+02  1.7462E+01  4.6125E+00  8.8289E+01  1.3762E+01  8.0410E+01  1.7899E+02  8.2509E+01  6.9302E+00
             1.1088E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1733.17181288737        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1109
 NPARAMETR:  1.0013E+00  3.7121E+00  4.9780E-01  3.2681E-02  3.2970E+00  1.5845E+00  2.4689E+00  5.0984E-01  2.6419E+01  4.8451E-01
             5.8523E+00
 PARAMETER:  1.0131E-01  1.4116E+00 -5.9756E-01 -3.3210E+00  1.2930E+00  5.6029E-01  1.0038E+00 -5.7366E-01  3.3741E+00 -6.2461E-01
             1.8668E+00
 GRADIENT:  -4.7234E+01  1.6380E+02  1.9701E+01  3.9569E+00  1.7772E+02 -4.5020E+01  3.0422E+01  1.8649E-01  2.4482E+01  2.4595E+00
             3.8396E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1736.81012083213        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1267
 NPARAMETR:  1.0028E+00  3.6073E+00  5.0165E-01  3.1155E-02  3.2741E+00  1.5971E+00  2.4999E+00  4.9704E-01  2.5438E+01  5.8772E-02
             5.7447E+00
 PARAMETER:  1.0275E-01  1.3830E+00 -5.8985E-01 -3.3688E+00  1.2860E+00  5.6817E-01  1.0162E+00 -5.9908E-01  3.3362E+00 -2.7341E+00
             1.8483E+00
 GRADIENT:  -5.4178E+01  1.2208E+01  1.5175E+01  6.4009E-01  1.8635E+02 -6.2517E+01 -2.4190E+01 -1.8854E-01 -4.6869E+00  3.0192E-02
            -4.2232E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1738.98878681163        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1427
 NPARAMETR:  1.0047E+00  3.6099E+00  4.9786E-01  3.1286E-02  3.2318E+00  1.5951E+00  2.4927E+00  4.5698E-01  2.5605E+01  2.1518E-02
             5.7740E+00
 PARAMETER:  1.0469E-01  1.3837E+00 -5.9744E-01 -3.3646E+00  1.2730E+00  5.6691E-01  1.0134E+00 -6.8312E-01  3.3428E+00 -3.7389E+00
             1.8534E+00
 GRADIENT:  -4.1004E+01  1.5495E+02  2.0178E+01  5.2889E+00  1.7726E+02 -4.1705E+01  3.6589E+01 -9.3515E-02  2.8027E+01  5.1400E-03
             7.3685E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1740.48210184146        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1574
 NPARAMETR:  1.0093E+00  3.6090E+00  4.9385E-01  3.1270E-02  3.2287E+00  1.6095E+00  2.4924E+00  4.5826E-01  2.5595E+01  1.9802E-02
             5.7756E+00
 PARAMETER:  1.0923E-01  1.3834E+00 -6.0552E-01 -3.3651E+00  1.2721E+00  5.7591E-01  1.0133E+00 -6.8032E-01  3.3424E+00 -3.8220E+00
             1.8536E+00
 GRADIENT:  -4.8688E+01  4.3555E+00  1.4551E+01  4.9423E+00  1.6630E+02 -5.7936E+01 -3.0241E+01 -5.3464E-03  2.1594E+00  3.4594E-03
            -2.6908E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1741.15895792024        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1754
 NPARAMETR:  1.0101E+00  3.6197E+00  4.9536E-01  3.1223E-02  3.2295E+00  1.6417E+00  2.4918E+00  4.3558E-01  2.5513E+01  1.0000E-02
             5.8253E+00
 PARAMETER:  1.1002E-01  1.3864E+00 -6.0247E-01 -3.3666E+00  1.2723E+00  5.9574E-01  1.0130E+00 -7.3109E-01  3.3392E+00 -4.8060E+00
             1.8622E+00
 GRADIENT:  -4.6994E+01  1.0544E+01  1.6339E+01  1.6207E+00  1.7278E+02 -4.9226E+01 -2.7277E+01 -1.2551E-01 -3.9851E+00  0.0000E+00
            -4.7700E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1741.86889421508        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1912             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0094E+00  3.6276E+00  4.9317E-01  3.1161E-02  3.2203E+00  1.6362E+00  2.4911E+00  4.6602E-01  2.5467E+01  1.0000E-02
             5.8724E+00
 PARAMETER:  1.0933E-01  1.3886E+00 -6.0691E-01 -3.3686E+00  1.2695E+00  5.9237E-01  1.0127E+00 -6.6353E-01  3.3374E+00 -4.8060E+00
             1.8703E+00
 GRADIENT:  -3.6901E+01  1.5281E+02  1.8877E+01  4.5031E+00  1.7314E+02 -2.9057E+01  3.4942E+01  2.5358E-01  2.7146E+01  0.0000E+00
             4.7592E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1742.62675492337        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     2045
 NPARAMETR:  1.0092E+00  3.6321E+00  4.8878E-01  3.1110E-02  3.2033E+00  1.6331E+00  2.4856E+00  2.7338E-01  2.5401E+01  1.0000E-02
             5.8681E+00
 PARAMETER:  1.0913E-01  1.3898E+00 -6.1584E-01 -3.3702E+00  1.2642E+00  5.9046E-01  1.0105E+00 -1.1969E+00  3.3348E+00 -4.8060E+00
             1.8695E+00
 GRADIENT:  -3.5332E+01  1.5300E+02  2.8315E+01  5.8212E+00  1.6695E+02 -3.0159E+01  3.2996E+01 -4.0546E-01  2.8946E+01  0.0000E+00
             4.4689E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1747.85719256371        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2224             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0106E+00  3.6543E+00  4.8869E-01  3.1067E-02  3.1986E+00  1.6348E+00  2.4773E+00  1.4949E+00  2.5347E+01  1.0000E-02
             5.9070E+00
 PARAMETER:  1.1051E-01  1.3959E+00 -6.1602E-01 -3.3716E+00  1.2627E+00  5.9151E-01  1.0072E+00  5.0205E-01  3.3327E+00 -4.8060E+00
             1.8761E+00
 GRADIENT:  -2.8776E+01  1.5644E+02  9.6275E+00  3.6241E+00  1.6788E+02 -2.2083E+01  3.5413E+01 -3.7651E+00  3.2791E+01  0.0000E+00
             8.6126E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1751.80528376364        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     2388             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0105E+00  3.6642E+00  4.8902E-01  3.1044E-02  3.1954E+00  1.6560E+00  2.4845E+00  3.5205E+00  2.5352E+01  1.0000E-02
             5.8943E+00
 PARAMETER:  1.1047E-01  1.3986E+00 -6.1535E-01 -3.3724E+00  1.2617E+00  6.0442E-01  1.0101E+00  1.3586E+00  3.3329E+00 -4.8060E+00
             1.8740E+00
 GRADIENT:  -3.0468E+01  1.6011E+02  3.6096E+00  3.3616E+00  1.6989E+02 -1.6529E+01  4.2218E+01  1.8729E+01  3.7155E+01  0.0000E+00
             1.0552E+02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1753.84320765962        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     2535
 NPARAMETR:  1.0125E+00  3.6694E+00  4.8860E-01  3.0989E-02  3.1962E+00  1.7285E+00  2.4748E+00  3.3039E+00  2.5339E+01  1.0000E-02
             5.7716E+00
 PARAMETER:  1.1242E-01  1.4000E+00 -6.1621E-01 -3.3741E+00  1.2620E+00  6.4727E-01  1.0061E+00  1.2951E+00  3.3323E+00 -4.8060E+00
             1.8529E+00
 GRADIENT:  -3.5701E+01  1.7209E+01 -4.0044E+00  1.0319E+00  1.7095E+02 -2.0595E+01 -2.2673E+01 -2.0878E+00 -1.7107E+01  0.0000E+00
             2.7139E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1754.13777570895        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2733
 NPARAMETR:  1.0135E+00  3.6724E+00  4.8904E-01  3.0948E-02  3.1940E+00  1.7272E+00  2.4743E+00  3.5282E+00  2.5324E+01  1.0000E-02
             5.7730E+00
 PARAMETER:  1.1340E-01  1.4009E+00 -6.1532E-01 -3.3754E+00  1.2613E+00  6.4652E-01  1.0060E+00  1.3608E+00  3.3318E+00 -4.8060E+00
             1.8532E+00
 GRADIENT:  -3.5295E+01  1.7655E+01 -3.4418E+00  9.4048E-01  1.7068E+02 -2.0828E+01 -2.2886E+01  5.1247E-01 -1.6833E+01  0.0000E+00
             2.7877E+01

0ITERATION NO.:   86    OBJECTIVE VALUE:  -1754.13777570895        NO. OF FUNC. EVALS.:  33
 CUMULATIVE NO. OF FUNC. EVALS.:     2766
 NPARAMETR:  1.0134E+00  3.6742E+00  4.8898E-01  3.0925E-02  3.1941E+00  1.7265E+00  2.4740E+00  3.5766E+00  2.5308E+01  1.0000E-02
             5.7703E+00
 PARAMETER:  1.1340E-01  1.4009E+00 -6.1532E-01 -3.3754E+00  1.2613E+00  6.4652E-01  1.0060E+00  1.3608E+00  3.3318E+00 -4.8060E+00
             1.8532E+00
 GRADIENT:   2.5522E+03 -3.6883E+02  1.9615E+02  7.4594E+01 -3.3483E+01  4.2456E+02  1.4752E+02 -1.5537E+00  6.9020E+01  0.0000E+00
             9.4928E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2766
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1205E-02 -3.4932E-02  8.5938E-03  1.6801E-02 -9.1160E-04
 SE:             2.9745E-02  3.0278E-02  8.8602E-03  1.0741E-02  1.3912E-04
 N:                     100         100         100         100         100

 P VAL.:         2.9414E-01  2.4862E-01  3.3208E-01  1.1776E-01  5.6888E-11

 ETASHRINKSD(%)  3.5060E-01  1.0000E-10  7.0317E+01  6.4018E+01  9.9534E+01
 ETASHRINKVR(%)  6.9997E-01  1.0000E-10  9.1189E+01  8.7053E+01  9.9998E+01
 EBVSHRINKSD(%)  5.8772E+00  1.7295E+00  7.3391E+01  6.9072E+01  9.9430E+01
 EBVSHRINKVR(%)  1.1409E+01  3.4290E+00  9.2920E+01  9.0434E+01  9.9997E+01
 RELATIVEINF(%)  8.8459E+01  8.6582E+01  7.0729E+00  8.6421E+00  3.1812E-03
 EPSSHRINKSD(%)  1.0719E+01
 EPSSHRINKVR(%)  2.0288E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1754.1377757089540     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -100.04841594054324     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   143.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    19.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1754.138       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  3.67E+00  4.89E-01  3.09E-02  3.19E+00  1.73E+00  2.47E+00  3.53E+00  2.53E+01  1.00E-02  5.77E+00
 


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
+        5.54E+05
 
 TH 2
+        1.25E+04  1.09E+03
 
 TH 3
+       -1.23E+02 -2.94E+01  1.28E+05
 
 TH 4
+        2.80E+03 -1.75E+04  1.19E+03  1.92E+06
 
 TH 5
+       -7.35E+01  8.22E+02 -3.84E+01 -2.29E+04  1.36E+03
 
 TH 6
+       -2.35E+02  1.31E+03 -2.01E+01  2.56E+03 -6.69E+01  5.73E+03
 
 TH 7
+       -1.74E+00  8.80E+02  8.88E+00  2.80E+04  1.63E+01 -3.40E-01  2.54E+03
 
 TH 8
+       -1.17E+02 -3.71E-01 -4.46E+01  1.46E+02 -3.86E+00 -1.20E+01  1.62E-01  3.23E+02
 
 TH 9
+        3.53E+00 -1.34E+01  1.86E+00  1.12E+03 -3.04E+01  3.20E+00  2.01E+01  1.93E-01  3.04E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        5.97E+03  9.98E+01  2.29E+03 -2.28E+02 -4.03E+00  4.26E+00  2.79E+02 -1.28E+00 -4.49E-01  0.00E+00  2.49E+02
 
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
 #CPUT: Total CPU Time in Seconds,      163.187
Stop Time:
Wed Sep 29 08:22:37 CDT 2021
