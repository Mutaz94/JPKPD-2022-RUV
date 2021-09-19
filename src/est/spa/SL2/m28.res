Sat Sep 18 12:10:36 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat28.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1635.72246263813        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.6752E+01 -1.1307E+02 -5.2171E+01 -1.1310E+02  9.2660E+01  1.2723E+01  1.1764E-01  2.9176E+00 -1.7155E+01  2.0408E+01
            -2.0088E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1647.01736516276        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7311E-01  1.0059E+00  1.1169E+00  1.0814E+00  9.4460E-01  9.2829E-01  8.9156E-01  1.0138E+00  1.1047E+00  7.0335E-01
             1.1365E+00
 PARAMETER:  7.2739E-02  1.0591E-01  2.1057E-01  1.7824E-01  4.3010E-02  2.5585E-02 -1.4782E-02  1.1375E-01  1.9956E-01 -2.5191E-01
             2.2797E-01
 GRADIENT:   3.5431E+00  7.5512E+00  2.1191E+01  3.3944E-01 -3.7560E+00 -1.6297E+01  3.2792E+00 -9.4025E+00  1.2004E+01 -1.4492E+01
             2.3188E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1650.00950186289        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.7314E-01  1.0788E+00  1.3276E+00  1.0367E+00  1.0540E+00  9.4595E-01  6.1825E-01  1.5191E+00  1.1599E+00  8.7065E-01
             1.0791E+00
 PARAMETER:  7.2777E-02  1.7583E-01  3.8340E-01  1.3607E-01  1.5264E-01  4.4430E-02 -3.8086E-01  5.1811E-01  2.4831E-01 -3.8513E-02
             1.7608E-01
 GRADIENT:   9.9466E+00  1.8023E+00  1.6976E+00 -1.6426E+00  1.7338E+00 -7.8529E+00  2.2258E+00  3.5038E+00  2.9448E+00 -1.9115E+00
             9.0122E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1650.59319544026        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.6821E-01  1.0043E+00  1.2075E+00  1.0802E+00  9.8854E-01  9.6478E-01  6.1848E-01  1.2629E+00  1.1102E+00  8.5564E-01
             1.0538E+00
 PARAMETER:  6.7698E-02  1.0426E-01  2.8858E-01  1.7719E-01  8.8473E-02  6.4142E-02 -3.8048E-01  3.3340E-01  2.0456E-01 -5.5901E-02
             1.5243E-01
 GRADIENT:   4.1246E-01  9.9928E-01  8.6423E-01 -2.6891E-01 -1.7690E+00 -1.8784E-01  6.1596E-01 -1.2501E-02  1.8663E-01  5.3617E-01
            -1.4132E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1650.89399732661        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  9.6804E-01  8.2650E-01  1.3803E+00  1.2002E+00  9.7603E-01  9.6164E-01  2.2626E-01  1.3374E+00  1.0711E+00  8.6726E-01
             1.0506E+00
 PARAMETER:  6.7513E-02 -9.0551E-02  4.2232E-01  2.8245E-01  7.5742E-02  6.0887E-02 -1.3860E+00  3.9071E-01  1.6869E-01 -4.2422E-02
             1.4936E-01
 GRADIENT:  -3.0094E+01  6.8825E+00  7.3933E+00  1.4180E+00 -1.5682E+01 -3.7364E+00 -1.4245E-01 -7.4464E-01  2.4738E+00 -5.5484E-02
            -1.7799E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1651.62885699206        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      504
 NPARAMETR:  9.7641E-01  6.3305E-01  1.3933E+00  1.3257E+00  9.2700E-01  9.6701E-01  9.9950E-02  1.2697E+00  9.6167E-01  8.3955E-01
             1.0550E+00
 PARAMETER:  7.6130E-02 -3.5721E-01  4.3167E-01  3.8196E-01  2.4201E-02  6.6457E-02 -2.2031E+00  3.3880E-01  6.0916E-02 -7.4887E-02
             1.5358E-01
 GRADIENT:  -5.6104E+00  2.2539E+00  1.8317E+00  1.7366E+00 -3.0872E+00 -7.2720E-01 -8.1771E-03 -2.8271E-02  4.4723E-01  1.7482E-02
            -7.0563E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1651.89673970984        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      680
 NPARAMETR:  9.7672E-01  4.2001E-01  1.3341E+00  1.4588E+00  8.4547E-01  9.6810E-01  2.2619E-02  1.1762E+00  8.6226E-01  7.8694E-01
             1.0552E+00
 PARAMETER:  7.6441E-02 -7.6747E-01  3.8828E-01  4.7762E-01 -6.7859E-02  6.7584E-02 -3.6890E+00  2.6226E-01 -4.8195E-02 -1.3960E-01
             1.5369E-01
 GRADIENT:   5.2326E-01  2.3058E+00  1.7003E+00  7.0802E+00 -3.6349E+00  1.9618E-01 -1.9275E-04 -3.9161E-01 -1.0672E+00 -9.6041E-02
            -3.8229E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1651.96915628746        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      856
 NPARAMETR:  9.7512E-01  2.9764E-01  1.3278E+00  1.5318E+00  8.1143E-01  9.6713E-01  1.0000E-02  1.1793E+00  8.1530E-01  7.6198E-01
             1.0560E+00
 PARAMETER:  7.4801E-02 -1.1119E+00  3.8350E-01  5.2642E-01 -1.0896E-01  6.6579E-02 -5.1363E+00  2.6491E-01 -1.0420E-01 -1.7183E-01
             1.5451E-01
 GRADIENT:   9.3268E-01  9.8318E-01  1.0098E+00  4.8283E+00 -1.9873E+00  1.5193E-01  0.0000E+00 -5.8921E-02 -8.7387E-01 -7.9998E-03
            -8.8053E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1652.01349740764        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1032
 NPARAMETR:  9.7331E-01  1.9034E-01  1.2961E+00  1.5890E+00  7.7417E-01  9.6627E-01  1.0000E-02  1.1644E+00  7.8143E-01  7.3771E-01
             1.0571E+00
 PARAMETER:  7.2947E-02 -1.5590E+00  3.5937E-01  5.6308E-01 -1.5596E-01  6.5685E-02 -7.0944E+00  2.5223E-01 -1.4663E-01 -2.0420E-01
             1.5552E-01
 GRADIENT:   8.6115E-01 -2.6179E-01  7.0738E-01 -5.2124E+00 -1.1339E+00  1.2724E-01  0.0000E+00  5.7052E-02  4.7560E-01  3.8205E-01
             4.0527E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1652.04879882536        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1211
 NPARAMETR:  9.7127E-01  1.2079E-01  1.2972E+00  1.6346E+00  7.5834E-01  9.6493E-01  1.0000E-02  1.1788E+00  7.5706E-01  7.1971E-01
             1.0562E+00
 PARAMETER:  7.0848E-02 -2.0137E+00  3.6019E-01  5.9143E-01 -1.7662E-01  6.4301E-02 -9.1485E+00  2.6452E-01 -1.7831E-01 -2.2891E-01
             1.5471E-01
 GRADIENT:  -1.0986E+00  3.1905E-01 -2.5314E-02  4.8259E+00 -2.4132E-01 -2.1231E-01  0.0000E+00  2.4238E-02 -4.5110E-01 -2.8168E-01
            -2.3471E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1652.07459180056        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1387
 NPARAMETR:  9.7081E-01  6.9079E-02  1.2764E+00  1.6632E+00  7.3937E-01  9.6526E-01  1.0000E-02  1.1676E+00  7.4138E-01  7.1002E-01
             1.0569E+00
 PARAMETER:  7.0374E-02 -2.5725E+00  3.4407E-01  6.0876E-01 -2.0195E-01  6.4643E-02 -1.1707E+01  2.5493E-01 -1.9925E-01 -2.4246E-01
             1.5538E-01
 GRADIENT:   1.8352E-01  2.0654E-01  4.4603E-03  5.3293E+00 -1.1347E+00  9.4405E-02  0.0000E+00  8.7514E-02 -6.3310E-01  2.1105E-01
             4.9519E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1652.09641876950        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1562
 NPARAMETR:  9.7005E-01  3.3775E-02  1.2774E+00  1.6822E+00  7.3174E-01  9.6476E-01  1.0000E-02  1.1737E+00  7.3177E-01  7.0147E-01
             1.0568E+00
 PARAMETER:  6.9596E-02 -3.2880E+00  3.4486E-01  6.2010E-01 -2.1233E-01  6.4126E-02 -1.5051E+01  2.6012E-01 -2.1229E-01 -2.5457E-01
             1.5529E-01
 GRADIENT:   1.1328E-01  4.9809E-02  3.9281E-01  1.7239E+00 -7.2613E-01  1.6348E-02  0.0000E+00 -7.1199E-02 -2.4528E-01 -9.3967E-02
            -6.6815E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1652.10621329434        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1738
 NPARAMETR:  9.6953E-01  1.4138E-02  1.2742E+00  1.6932E+00  7.2634E-01  9.6447E-01  1.0000E-02  1.1741E+00  7.2642E-01  6.9825E-01
             1.0569E+00
 PARAMETER:  6.9054E-02 -4.1589E+00  3.4231E-01  6.2662E-01 -2.1974E-01  6.3823E-02 -1.9139E+01  2.6048E-01 -2.1963E-01 -2.5918E-01
             1.5538E-01
 GRADIENT:  -1.5375E-01  1.9760E-02  3.9100E-01  1.5977E+00 -8.5080E-01 -2.3743E-02  0.0000E+00 -2.0196E-02 -1.4015E-01  2.3368E-02
            -1.5387E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1652.10879502345        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1913
 NPARAMETR:  9.6950E-01  1.0000E-02  1.2734E+00  1.6949E+00  7.2547E-01  9.6450E-01  1.0000E-02  1.1743E+00  7.2553E-01  6.9754E-01
             1.0570E+00
 PARAMETER:  6.9030E-02 -4.5308E+00  3.4169E-01  6.2761E-01 -2.2094E-01  6.3853E-02 -2.0891E+01  2.6067E-01 -2.2086E-01 -2.6020E-01
             1.5540E-01
 GRADIENT:   1.7186E-02  0.0000E+00  6.6074E-03  7.1744E-03 -2.3288E-02 -1.9281E-06  0.0000E+00  5.6668E-04 -1.1557E-03  5.3969E-03
             2.6050E-03

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1652.10879502345        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1935
 NPARAMETR:  9.6950E-01  1.0000E-02  1.2734E+00  1.6949E+00  7.2547E-01  9.6450E-01  1.0000E-02  1.1743E+00  7.2553E-01  6.9754E-01
             1.0570E+00
 PARAMETER:  6.9030E-02 -4.5308E+00  3.4169E-01  6.2761E-01 -2.2094E-01  6.3853E-02 -2.0891E+01  2.6067E-01 -2.2086E-01 -2.6020E-01
             1.5540E-01
 GRADIENT:   1.7186E-02  0.0000E+00  6.6074E-03  7.1744E-03 -2.3288E-02 -1.9281E-06  0.0000E+00  5.6668E-04 -1.1557E-03  5.3969E-03
             2.6050E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1935
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5194E-04 -2.7166E-06 -2.1903E-02 -4.7546E-03 -2.9106E-02
 SE:             2.9830E-02  1.8175E-06  1.9918E-02  2.9213E-02  1.8616E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9594E-01  1.3498E-01  2.7149E-01  8.7071E-01  1.1793E-01

 ETASHRINKSD(%)  6.5396E-02  9.9994E+01  3.3273E+01  2.1344E+00  3.7634E+01
 ETASHRINKVR(%)  1.3075E-01  1.0000E+02  5.5475E+01  4.2233E+00  6.1105E+01
 EBVSHRINKSD(%)  4.7205E-01  9.9994E+01  3.4334E+01  2.6277E+00  3.6731E+01
 EBVSHRINKVR(%)  9.4188E-01  1.0000E+02  5.6879E+01  5.1864E+00  5.9970E+01
 RELATIVEINF(%)  9.6075E+01  1.7227E-08  5.9122E+00  7.3188E+00  3.2106E+00
 EPSSHRINKSD(%)  4.4339E+01
 EPSSHRINKVR(%)  6.9019E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1652.1087950234464     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -916.95796845970824     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.38
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1652.109       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.70E-01  1.00E-02  1.27E+00  1.69E+00  7.25E-01  9.64E-01  1.00E-02  1.17E+00  7.26E-01  6.98E-01  1.06E+00
 


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
+        1.26E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.11E+00  0.00E+00  2.20E+02
 
 TH 4
+       -1.08E+01  0.00E+00 -2.40E+01  6.97E+02
 
 TH 5
+        1.68E+00  0.00E+00 -5.46E+02 -1.22E+02  1.76E+03
 
 TH 6
+        5.31E+00  0.00E+00  8.44E-01 -2.25E+00 -3.25E-01  2.15E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -3.17E+00  0.00E+00 -3.22E+01 -2.92E+00 -2.00E+01  4.32E-01  0.00E+00  4.18E+01
 
 TH 9
+        3.69E+00  0.00E+00  5.28E+00 -1.16E+00  1.03E+00  3.32E+00  0.00E+00  1.10E-01  3.38E+02
 
 TH10
+       -1.89E+00  0.00E+00  2.44E+00 -8.37E-01 -1.20E+02  1.78E+00  0.00E+00  2.45E+01  3.39E+00  9.15E+01
 
 TH11
+       -3.53E+00  0.00E+00 -2.70E+00 -7.09E+00 -1.57E+01  5.81E+00  0.00E+00  1.05E+01  1.26E+01  1.74E+01  1.88E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.844
Stop Time:
Sat Sep 18 12:11:04 CDT 2021
