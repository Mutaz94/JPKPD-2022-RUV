Sat Sep 18 10:28:16 CDT 2021
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
$DATA ../../../../data/spa/A3/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   896.736670725003        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.4539E+00  7.8748E+01  1.4090E+02 -8.8812E+01  7.4324E+01  1.2051E+01 -6.2260E+01 -4.6403E+01 -1.4462E+02 -1.5378E+02
            -4.7450E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -825.606202421660        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.7405E+00  1.2098E+00  8.8303E-01  1.8521E+00  6.3864E-01  8.9501E-01  1.0322E+00  1.0139E+00  9.6759E-01  1.3803E+00
             1.5675E+01
 PARAMETER:  6.5415E-01  2.9047E-01 -2.4392E-02  7.1635E-01 -3.4842E-01 -1.0916E-02  1.3165E-01  1.1384E-01  6.7050E-02  4.2226E-01
             2.8521E+00
 GRADIENT:   2.5309E+02  1.1716E+01  1.2151E+01  6.3925E+00 -4.1464E+01 -1.0635E+01  7.2721E+00  5.5484E+00  2.5719E+01  2.2527E+01
             2.6153E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1042.30613369028        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.9441E-01  3.2465E-01  2.4010E-01  1.8253E+00  2.6604E-01  8.9799E-01  9.8140E-01  2.3507E-01  2.1758E-01  2.1520E-01
             8.4378E+00
 PARAMETER:  9.4398E-02 -1.0250E+00 -1.3267E+00  7.0177E-01 -1.2241E+00 -7.5995E-03  8.1223E-02 -1.3479E+00 -1.4252E+00 -1.4362E+00
             2.2327E+00
 GRADIENT:  -3.1075E+02  1.8218E+01  8.1604E-01  2.1919E+02 -3.4017E+01 -2.4379E+01  1.8870E+00  1.4525E+00  2.0280E+00  2.3869E+00
             2.3667E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1220.81237979796        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.1680E+00  5.8786E-01  1.9329E-01  1.1046E+00  3.3879E-01  1.0797E+00  6.3149E-02  1.0835E-02  3.7476E-01  1.3323E-01
             5.4482E+00
 PARAMETER:  2.5532E-01 -4.3126E-01 -1.5436E+00  1.9951E-01 -9.8236E-01  1.7664E-01 -2.6623E+00 -4.4249E+00 -8.8147E-01 -1.9156E+00
             1.7953E+00
 GRADIENT:   1.5115E+02 -9.4259E+01 -8.7614E+01  1.9379E+01  1.6047E+02  1.7949E+01  7.5072E-02  2.6232E-03  2.4992E+00  1.3559E+00
             1.0460E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1252.60769638517        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0679E+00  8.9153E-01  4.1619E-01  1.1083E+00  5.7957E-01  9.0671E-01  1.0000E-02  1.0000E-02  5.7805E-01  1.3493E-01
             4.6126E+00
 PARAMETER:  1.6569E-01 -1.4815E-02 -7.7661E-01  2.0278E-01 -4.4546E-01  2.0692E-03 -5.1936E+00 -4.6099E+00 -4.4810E-01 -1.9030E+00
             1.6288E+00
 GRADIENT:  -1.9306E+01  2.7156E+01  6.9811E+00  2.9532E+01 -3.1553E+01 -1.1558E+01  0.0000E+00  0.0000E+00 -1.5729E+00  6.0488E-01
            -1.1433E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1257.45762746808        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0783E+00  1.2391E+00  9.8171E-01  1.0035E+00  1.2282E+00  9.1254E-01  1.0000E-02  1.0000E-02  7.5413E-01  5.2398E-02
             4.6612E+00
 PARAMETER:  1.7542E-01  3.1437E-01  8.1541E-02  1.0349E-01  3.0555E-01  8.4820E-03 -8.6092E+00 -5.5798E+00 -1.8219E-01 -2.8489E+00
             1.6393E+00
 GRADIENT:   1.2851E+01  1.0910E+01 -1.5507E+00  1.9352E+01  2.1813E+00  2.5086E+00  0.0000E+00  0.0000E+00  1.5667E+00  3.1630E-02
             1.1293E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1257.81971288630        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  1.0669E+00  1.2116E+00  1.0360E+00  9.9772E-01  1.2266E+00  9.0273E-01  1.0000E-02  1.0000E-02  6.9858E-01  5.1162E-02
             4.6073E+00
 PARAMETER:  1.6478E-01  2.9191E-01  1.3536E-01  9.7720E-02  3.0428E-01 -2.3308E-03 -8.5876E+00 -5.1980E+00 -2.5871E-01 -2.8728E+00
             1.6276E+00
 GRADIENT:   1.0637E-01 -5.4952E-01 -1.7325E-01  1.2970E+00  6.4999E-01  8.5249E-02  0.0000E+00  0.0000E+00 -1.2279E-01  3.0403E-02
            -3.4072E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1258.06595027006        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.0677E+00  1.4261E+00  8.9662E-01  8.6525E-01  1.2775E+00  9.0375E-01  1.0000E-02  1.0000E-02  8.3961E-01  2.4450E-02
             4.5925E+00
 PARAMETER:  1.6548E-01  4.5497E-01 -9.1188E-03 -4.4740E-02  3.4488E-01 -1.1984E-03 -9.8690E+00 -4.7906E+00 -7.4819E-02 -3.6111E+00
             1.6244E+00
 GRADIENT:  -4.7773E+00  1.1190E+01  6.0067E-01  9.5170E+00 -2.2313E+00  1.7680E-01  0.0000E+00  0.0000E+00 -2.0676E-01  6.6502E-03
            -2.7526E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1258.32491888306        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      595
 NPARAMETR:  1.0672E+00  1.7121E+00  8.1762E-01  6.8239E-01  1.4478E+00  9.0695E-01  1.0000E-02  1.6163E-02  1.0751E+00  1.0000E-02
             4.5790E+00
 PARAMETER:  1.6502E-01  6.3774E-01 -1.0135E-01 -2.8216E-01  4.7007E-01  2.3294E-03 -1.2220E+01 -4.0250E+00  1.7238E-01 -4.9350E+00
             1.6215E+00
 GRADIENT:  -7.3254E+00  1.8760E+01  9.9039E-01  9.3833E+00 -2.7071E+00  2.0674E+00  0.0000E+00  5.4396E-04 -2.3036E-01  0.0000E+00
            -4.9514E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1258.36299155073        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  1.0684E+00  1.6436E+00  8.1533E-01  7.1713E-01  1.3907E+00  9.0102E-01  1.0000E-02  1.3975E-02  1.0116E+00  1.0000E-02
             4.5930E+00
 PARAMETER:  1.6614E-01  5.9689E-01 -1.0416E-01 -2.3250E-01  4.2977E-01 -4.2317E-03 -1.1661E+01 -4.1705E+00  1.1151E-01 -4.6565E+00
             1.6245E+00
 GRADIENT:  -1.3543E+00  3.3265E+00  2.6129E-01  1.8081E+00 -5.8991E-01  3.7600E-01  0.0000E+00  4.8829E-04 -3.2154E-02  0.0000E+00
            -9.6145E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1258.36384002417        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      737
 NPARAMETR:  1.0686E+00  1.6404E+00  8.0787E-01  7.1756E-01  1.3835E+00  9.0009E-01  1.0000E-02  1.3830E-02  1.0094E+00  1.0000E-02
             4.5955E+00
 PARAMETER:  1.6638E-01  5.9491E-01 -1.1335E-01 -2.3190E-01  4.2464E-01 -5.2650E-03 -1.1639E+01 -4.1809E+00  1.0938E-01 -4.6543E+00
             1.6251E+00
 GRADIENT:  -3.7593E-01  9.3584E-01  7.7874E-02  5.2879E-01 -1.7555E-01  1.0346E-01  0.0000E+00  4.8937E-04 -8.7414E-03  0.0000E+00
            -2.7281E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1258.49903465186        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.0697E+00  1.7846E+00  7.1552E-01  6.2380E-01  1.4378E+00  8.9525E-01  1.0000E-02  1.8787E-02  1.1565E+00  1.0000E-02
             4.6033E+00
 PARAMETER:  1.6739E-01  6.7921E-01 -2.3474E-01 -3.7192E-01  4.6309E-01 -1.0653E-02 -1.2942E+01 -3.8746E+00  2.4538E-01 -5.4473E+00
             1.6268E+00
 GRADIENT:  -1.6393E-01  6.3756E+00  1.9852E-01  2.2161E+00 -1.6306E+00 -1.4518E+00  0.0000E+00  7.2268E-04 -2.1555E-01  0.0000E+00
             3.5096E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1258.65576043632        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  1.0713E+00  1.9452E+00  6.6468E-01  5.1847E-01  1.5633E+00  9.0545E-01  1.0000E-02  2.8194E-02  1.4170E+00  1.0000E-02
             4.5873E+00
 PARAMETER:  1.6892E-01  7.6534E-01 -3.0845E-01 -5.5688E-01  5.4683E-01  6.7543E-04 -1.4800E+01 -3.4687E+00  4.4857E-01 -6.4820E+00
             1.6233E+00
 GRADIENT:   7.1400E-01 -2.2206E+00 -2.3965E-01  3.4858E-01  6.6875E-01  2.6142E+00  0.0000E+00  1.1399E-03  1.6781E-01  0.0000E+00
            -1.6525E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1258.87066639736        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  1.0721E+00  2.1709E+00  5.9550E-01  3.6830E-01  1.7299E+00  8.9308E-01  1.0000E-02  9.7695E-02  1.9152E+00  1.0000E-02
             4.6067E+00
 PARAMETER:  1.6964E-01  8.7513E-01 -4.1835E-01 -8.9885E-01  6.4807E-01 -1.3080E-02 -1.7953E+01 -2.2259E+00  7.4985E-01 -8.3264E+00
             1.6275E+00
 GRADIENT:   1.1758E+00  4.6243E-01 -7.5744E-02  6.3851E-01  3.5746E-01 -9.6975E-01  0.0000E+00  7.7366E-03 -7.2463E-02  0.0000E+00
             2.0339E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1258.96300775121        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1326
 NPARAMETR:  1.0717E+00  2.2990E+00  5.7085E-01  2.7820E-01  1.8242E+00  8.9596E-01  1.0000E-02  3.3595E-01  2.4650E+00  1.0000E-02
             4.5804E+00
 PARAMETER:  1.6925E-01  9.3247E-01 -4.6062E-01 -1.1794E+00  7.0113E-01 -9.8602E-03 -2.0390E+01 -9.9080E-01  1.0022E+00 -9.7170E+00
             1.6218E+00
 GRADIENT:   8.1841E-01 -2.0205E+00 -3.6638E-01  2.6118E-01  8.0881E-01  2.2859E-01  0.0000E+00  7.8388E-02  5.6430E-02  0.0000E+00
            -3.0429E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1259.08077030139        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1503
 NPARAMETR:  1.0703E+00  2.4074E+00  5.6034E-01  2.0135E-01  1.9034E+00  8.9774E-01  1.0000E-02  1.3393E-01  3.1351E+00  1.0000E-02
             4.5700E+00
 PARAMETER:  1.6796E-01  9.7855E-01 -4.7920E-01 -1.5027E+00  7.4364E-01 -7.8797E-03 -2.3403E+01 -1.9104E+00  1.2427E+00 -1.1049E+01
             1.6195E+00
 GRADIENT:  -2.8236E+00 -2.4408E+00 -6.4230E-01  5.7306E-02  1.5179E+00  4.9766E-01  0.0000E+00  1.2293E-02  1.1613E-01  0.0000E+00
            -3.7009E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1259.13483352428        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1678
 NPARAMETR:  1.0713E+00  2.4372E+00  6.1810E-01  1.7707E-01  1.8993E+00  8.9539E-01  1.0000E-02  4.0675E-02  3.3891E+00  1.0000E-02
             4.5729E+00
 PARAMETER:  1.6888E-01  9.9083E-01 -3.8111E-01 -1.6312E+00  7.4148E-01 -1.0501E-02 -2.4541E+01 -3.1021E+00  1.3206E+00 -1.1290E+01
             1.6202E+00
 GRADIENT:  -1.2564E-01 -2.5935E+00 -1.1552E-01 -4.6109E-01  2.7651E-01 -2.2305E-01  0.0000E+00  1.0692E-03 -1.9166E-01  0.0000E+00
             8.6193E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1259.14025829736        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1856
 NPARAMETR:  1.0714E+00  2.4277E+00  6.3171E-01  1.8485E-01  1.8948E+00  8.9608E-01  1.0000E-02  3.1885E-02  3.3161E+00  1.0000E-02
             4.5681E+00
 PARAMETER:  1.6896E-01  9.8695E-01 -3.5933E-01 -1.5882E+00  7.3909E-01 -9.7236E-03 -2.4200E+01 -3.3456E+00  1.2988E+00 -1.1052E+01
             1.6191E+00
 GRADIENT:  -2.4597E-01  8.4756E-02  2.6439E-03 -7.7250E-03 -2.3504E-02 -1.4436E-02  0.0000E+00  6.5244E-04 -1.9416E-02  0.0000E+00
            -6.1611E-03

0ITERATION NO.:   89    OBJECTIVE VALUE:  -1259.14057492771        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1983
 NPARAMETR:  1.0715E+00  2.4274E+00  6.3140E-01  1.8514E-01  1.8949E+00  8.9613E-01  1.0000E-02  1.0000E-02  3.3139E+00  1.0000E-02
             4.5680E+00
 PARAMETER:  1.6906E-01  9.8682E-01 -3.5982E-01 -1.5867E+00  7.3918E-01 -9.6655E-03 -2.4389E+01 -4.8541E+00  1.2981E+00 -1.0910E+01
             1.6191E+00
 GRADIENT:  -3.2670E-02  4.8674E-02 -1.0766E-04  3.2444E-03 -5.3709E-03  2.4637E-03  0.0000E+00  0.0000E+00 -3.3146E-03  0.0000E+00
            -8.7197E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1983
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2729E-02 -1.0365E-03  2.0902E-05  3.3080E-02  4.3549E-05
 SE:             2.8617E-02  3.9063E-04  6.8970E-06  1.7490E-02  5.5535E-05
 N:                     100         100         100         100         100

 P VAL.:         6.5646E-01  7.9699E-03  2.4411E-03  5.8575E-02  4.3294E-01

 ETASHRINKSD(%)  4.1278E+00  9.8691E+01  9.9977E+01  4.1407E+01  9.9814E+01
 ETASHRINKVR(%)  8.0852E+00  9.9983E+01  1.0000E+02  6.5669E+01  1.0000E+02
 EBVSHRINKSD(%)  4.9995E+00  9.8613E+01  9.9942E+01  4.4786E+01  9.9763E+01
 EBVSHRINKVR(%)  9.7491E+00  9.9981E+01  1.0000E+02  6.9514E+01  9.9999E+01
 RELATIVEINF(%)  8.4211E+01  3.5871E-03  2.8742E-05  7.0187E+00  3.3131E-04
 EPSSHRINKSD(%)  1.4284E+01
 EPSSHRINKVR(%)  2.6528E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1259.1405749277142     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -523.98974836397599     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.19
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1259.141       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  2.43E+00  6.31E-01  1.85E-01  1.89E+00  8.96E-01  1.00E-02  1.00E-02  3.31E+00  1.00E-02  4.57E+00
 


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
+        1.13E+03
 
 TH 2
+       -4.31E+01  3.90E+02
 
 TH 3
+       -2.54E-01  4.45E+00  5.77E-02
 
 TH 4
+       -2.25E+02  4.51E+02  5.98E+00  9.31E+02
 
 TH 5
+       -1.08E+01 -5.39E+01 -5.72E-01 -4.80E+01  8.17E+00
 
 TH 6
+        4.67E+01 -3.12E+01 -1.44E+00 -3.07E+01 -5.36E-01  2.59E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.37E+01 -4.18E+00  2.06E-02  2.84E+01  1.76E+00  5.83E-01  0.00E+00  0.00E+00  2.76E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.82E+01 -3.68E+01 -3.95E-01 -6.32E+00  6.39E+00  1.16E+01  0.00E+00  0.00E+00  3.35E+00  0.00E+00  7.28E+00
 
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
+        1.08E+03
 
 TH 2
+       -1.43E+02  3.89E+02
 
 TH 3
+       -9.97E-01  4.83E+00  5.82E+00
 
 TH 4
+       -2.16E+02  4.57E+02  5.76E+00  9.30E+02
 
 TH 5
+        8.60E+00 -5.24E+01 -1.98E+00 -4.88E+01  1.88E+01
 
 TH 6
+        5.07E+00 -2.69E+01 -1.38E+00 -3.58E+01  4.54E-01  1.96E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.34E+00 -3.70E+00  5.82E-04  2.78E+01  9.66E-01 -2.73E-01  0.00E+00  0.00E+00  2.94E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.04E+01 -2.73E+01  9.35E-01 -1.02E+01  2.92E+00  5.60E+00  0.00E+00  0.00E+00  2.19E+00  0.00E+00  2.49E+01
 
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
+        1.09E+03
 
 TH 2
+       -2.03E+02  3.60E+02
 
 TH 3
+        9.58E-01  1.35E+00  9.39E-01
 
 TH 4
+       -2.18E+02  4.49E+02  1.24E+01  9.43E+02
 
 TH 5
+        2.24E+01 -3.25E+01 -6.90E-01 -3.73E+01  1.06E+01
 
 TH 6
+       -4.53E+01 -6.57E+00 -9.99E-01 -3.21E+01 -1.53E+00  1.45E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.52E+00 -2.17E+00  1.04E+00  2.67E+01  2.18E-01 -1.78E+00  0.00E+00  0.00E+00  2.53E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.02E+02 -7.11E+01  2.40E-01 -7.19E+01  7.90E+00  1.10E+01  0.00E+00  0.00E+00  1.46E+00  0.00E+00  7.57E+01
 
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
 #CPUT: Total CPU Time in Seconds,       30.469
Stop Time:
Sat Sep 18 10:28:48 CDT 2021
