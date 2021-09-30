Thu Sep 30 08:33:25 CDT 2021
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
$DATA ../../../../data/spa2/D/dat9.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1967.77088612678        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2014E+02 -3.5893E+01 -2.9067E+01 -8.0150E+01  1.3157E+02 -2.7259E+02 -2.1281E+02 -2.7725E+01 -3.0052E+02 -9.7600E+01
             3.7000E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2197.44000195844        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      170
 NPARAMETR:  8.9196E-01  9.6107E-01  1.0947E+00  1.0916E+00  1.0059E+00  1.7553E+00  2.2165E+00  1.1698E+00  1.9793E+00  1.3592E+00
             9.4232E-01
 PARAMETER: -1.4337E-02  6.0289E-02  1.9052E-01  1.8764E-01  1.0591E-01  6.6264E-01  8.9595E-01  2.5681E-01  7.8273E-01  4.0691E-01
             4.0594E-02
 GRADIENT:  -6.1589E+01 -1.9840E+01 -4.5879E+01  1.2167E+01  8.8142E+00  3.0070E+01  2.7438E+00 -5.0576E+00  5.5304E+01  2.9724E+01
             4.8256E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2216.21870377873        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  9.1727E-01  9.1763E-01  1.4811E+00  1.0929E+00  1.1156E+00  1.6975E+00  3.3758E+00  1.7068E+00  1.5060E+00  1.3855E+00
             9.1360E-01
 PARAMETER:  1.3645E-02  1.4041E-02  4.9277E-01  1.8882E-01  2.0940E-01  6.2913E-01  1.3166E+00  6.3459E-01  5.0947E-01  4.2608E-01
             9.6399E-03
 GRADIENT:  -4.3372E+01 -2.0400E+00 -1.1647E+01 -5.0057E+00 -2.0964E+01  2.0208E+01  3.3644E+01  1.2068E+01  4.9125E+01  2.3514E+01
            -3.1113E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2228.47722711467        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.6751E-01  8.3053E-01  1.6028E+00  1.1616E+00  1.1071E+00  1.6335E+00  3.3107E+00  1.6809E+00  1.2370E+00  1.2694E+00
             9.4290E-01
 PARAMETER:  6.6969E-02 -8.5693E-02  5.7174E-01  2.4980E-01  2.0176E-01  5.9072E-01  1.2971E+00  6.1932E-01  3.1268E-01  3.3854E-01
             4.1203E-02
 GRADIENT:  -2.4791E+00 -1.3682E+00 -4.4000E+00  5.9471E+00 -2.4278E+01  8.9120E-01  8.9951E+00  9.1386E+00  1.9993E+01  6.3275E+00
            -3.2571E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2233.21615380575        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.7244E-01  1.0507E+00  1.5445E+00  1.0228E+00  1.2385E+00  1.6341E+00  2.9216E+00  1.6373E+00  9.6636E-01  1.3553E+00
             9.5681E-01
 PARAMETER:  7.2053E-02  1.4946E-01  5.3469E-01  1.2254E-01  3.1390E-01  5.9107E-01  1.1721E+00  5.9307E-01  6.5782E-02  4.0406E-01
             5.5849E-02
 GRADIENT:   9.1962E-01  2.6886E+00  1.9077E+00  3.7039E+00 -2.9790E+00  7.4907E-01 -2.8642E+00  8.9306E-01 -1.0168E+00 -1.0088E+00
            -2.1094E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2233.35992636570        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      858
 NPARAMETR:  9.7200E-01  1.0326E+00  1.5296E+00  1.0170E+00  1.2401E+00  1.7282E+00  2.8874E+00  1.6095E+00  9.7683E-01  1.3621E+00
             9.5916E-01
 PARAMETER:  7.1600E-02  1.3204E-01  5.2499E-01  1.1684E-01  3.1516E-01  6.4708E-01  1.1604E+00  5.7595E-01  7.6560E-02  4.0902E-01
             5.8300E-02
 GRADIENT:   4.0137E+02  6.9013E+01  1.8494E+01  9.5758E+01  6.2354E+01  7.2138E+02  6.2123E+02  5.4193E+00  5.8550E+00  1.5428E+01
             1.9570E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2233.47166427521        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      999
 NPARAMETR:  9.7954E-01  1.0506E+00  1.5187E+00  1.0162E+00  1.2404E+00  1.6862E+00  2.8858E+00  1.6039E+00  9.6638E-01  1.3606E+00
             9.5862E-01
 PARAMETER:  7.9323E-02  1.4936E-01  5.1787E-01  1.1604E-01  3.1540E-01  6.2246E-01  1.1598E+00  5.7246E-01  6.5802E-02  4.0795E-01
             5.7735E-02
 GRADIENT:   4.0774E+02  8.2570E+01  1.7481E+01  9.9811E+01  6.2642E+01  6.7862E+02  6.1259E+02  5.2032E+00  4.6161E+00  1.4884E+01
             7.7257E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2233.50763716898        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1142             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7327E-01  1.0551E+00  1.5139E+00  1.0111E+00  1.2406E+00  1.6865E+00  2.8854E+00  1.6019E+00  9.8624E-01  1.3613E+00
             9.5886E-01
 PARAMETER:  7.2903E-02  1.5367E-01  5.1472E-01  1.1107E-01  3.1558E-01  6.2265E-01  1.1597E+00  5.7117E-01  8.6149E-02  4.0847E-01
             5.7991E-02
 GRADIENT:   4.0254E+02  8.5118E+01  1.7448E+01  9.0985E+01  6.2175E+01  6.8351E+02  6.1302E+02  5.6252E+00  6.5906E+00  1.5225E+01
             1.0980E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2233.51216972890        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1322             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7314E-01  1.0572E+00  1.5101E+00  1.0098E+00  1.2416E+00  1.6877E+00  2.8864E+00  1.5979E+00  9.8054E-01  1.3628E+00
             9.5906E-01
 PARAMETER:  7.2772E-02  1.5567E-01  5.1217E-01  1.0978E-01  3.1639E-01  6.2339E-01  1.1600E+00  5.6868E-01  8.0349E-02  4.0952E-01
             5.8196E-02
 GRADIENT:   4.0226E+02  8.6523E+01  1.7092E+01  9.0013E+01  6.2746E+01  6.8450E+02  6.1367E+02  5.3267E+00  6.1277E+00  1.5296E+01
             1.1477E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2233.51342897719        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1509
 NPARAMETR:  9.7310E-01  1.0594E+00  1.5076E+00  1.0080E+00  1.2419E+00  1.6892E+00  2.8809E+00  1.5960E+00  9.8267E-01  1.3636E+00
             9.5916E-01
 PARAMETER:  7.2727E-02  1.5773E-01  5.1052E-01  1.0793E-01  3.1661E-01  6.2425E-01  1.1581E+00  5.6749E-01  8.2518E-02  4.1012E-01
             5.8305E-02
 GRADIENT:   1.4088E+00 -5.7374E-01  3.9059E-01 -1.0869E+00 -7.1562E-01  1.5973E+01 -2.4419E+00  1.3140E-01  2.5684E-01  1.1810E-01
            -1.4351E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2233.51582886950        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1699             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7310E-01  1.0632E+00  1.5044E+00  1.0079E+00  1.2425E+00  1.6892E+00  2.8747E+00  1.5935E+00  9.8046E-01  1.3637E+00
             9.5923E-01
 PARAMETER:  7.2734E-02  1.6129E-01  5.0839E-01  1.0785E-01  3.1711E-01  6.2425E-01  1.1559E+00  5.6592E-01  8.0270E-02  4.1022E-01
             5.8380E-02
 GRADIENT:   4.0206E+02  9.0496E+01  1.6680E+01  8.8901E+01  6.3126E+01  6.8560E+02  6.0854E+02  5.1995E+00  5.9121E+00  1.5378E+01
             1.1522E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2233.51641953731        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1889
 NPARAMETR:  9.7311E-01  1.0646E+00  1.5027E+00  1.0071E+00  1.2427E+00  1.6892E+00  2.8718E+00  1.5922E+00  9.8065E-01  1.3639E+00
             9.5929E-01
 PARAMETER:  7.2740E-02  1.6256E-01  5.0725E-01  1.0712E-01  3.1731E-01  6.2425E-01  1.1550E+00  5.6515E-01  8.0457E-02  4.1033E-01
             5.8436E-02
 GRADIENT:   1.4046E+00  3.3754E-02 -1.1505E-02  6.0803E-01 -1.4761E-01  1.5971E+01 -2.8243E+00  5.2883E-02 -9.3653E-02  4.6352E-02
            -4.5268E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2233.51679516083        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2079             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7312E-01  1.0653E+00  1.5014E+00  1.0054E+00  1.2427E+00  1.6892E+00  2.8699E+00  1.5910E+00  9.8328E-01  1.3644E+00
             9.5937E-01
 PARAMETER:  7.2750E-02  1.6325E-01  5.0639E-01  1.0543E-01  3.1731E-01  6.2425E-01  1.1543E+00  5.6435E-01  8.3134E-02  4.1069E-01
             5.8526E-02
 GRADIENT:   4.0196E+02  9.1649E+01  1.6724E+01  8.5316E+01  6.2816E+01  6.8541E+02  6.0777E+02  5.1720E+00  6.1401E+00  1.5463E+01
             1.2087E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2233.51735490528        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2266
 NPARAMETR:  9.7312E-01  1.0665E+00  1.5003E+00  1.0054E+00  1.2430E+00  1.6892E+00  2.8678E+00  1.5903E+00  9.8234E-01  1.3643E+00
             9.5939E-01
 PARAMETER:  7.2751E-02  1.6436E-01  5.0565E-01  1.0538E-01  3.1753E-01  6.2425E-01  1.1536E+00  5.6395E-01  8.2182E-02  4.1067E-01
             5.8544E-02
 GRADIENT:   1.4059E+00 -1.0886E-01  6.1994E-02 -1.1299E-01 -3.2828E-01  1.5967E+01 -2.7150E+00  5.7331E-02  3.3927E-02  8.4697E-02
            -1.4899E-02

0ITERATION NO.:   66    OBJECTIVE VALUE:  -2233.51735490528        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     2290
 NPARAMETR:  9.7312E-01  1.0665E+00  1.5003E+00  1.0054E+00  1.2430E+00  1.6892E+00  2.8678E+00  1.5903E+00  9.8234E-01  1.3643E+00
             9.5939E-01
 PARAMETER:  7.2751E-02  1.6436E-01  5.0565E-01  1.0538E-01  3.1753E-01  6.2425E-01  1.1536E+00  5.6395E-01  8.2182E-02  4.1067E-01
             5.8544E-02
 GRADIENT:  -1.0619E-03 -8.0789E-02  9.8326E-02 -5.8329E-03 -7.7410E-02 -1.4420E-04  1.0633E-01  3.1012E-02  2.6639E-02 -2.8815E-03
            -1.4340E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2290
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1849E-04  9.5576E-03 -4.6789E-02 -1.4021E-02 -2.2541E-02
 SE:             2.9939E-02  2.5979E-02  1.6623E-02  2.0706E-02  2.3536E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9151E-01  7.1295E-01  4.8823E-03  4.9831E-01  3.3820E-01

 ETASHRINKSD(%)  1.0000E-10  1.2966E+01  4.4311E+01  3.0634E+01  2.1151E+01
 ETASHRINKVR(%)  1.0000E-10  2.4252E+01  6.8988E+01  5.1884E+01  3.7829E+01
 EBVSHRINKSD(%)  1.1075E-01  1.0350E+01  4.9195E+01  3.5895E+01  1.7317E+01
 EBVSHRINKVR(%)  2.2137E-01  1.9628E+01  7.4189E+01  5.8905E+01  3.1636E+01
 RELATIVEINF(%)  9.9770E+01  3.9801E+01  1.3369E+01  1.5815E+01  4.1234E+01
 EPSSHRINKSD(%)  3.1225E+01
 EPSSHRINKVR(%)  5.2701E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2233.5173549052774     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1130.7911150596703     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    52.71
 Elapsed covariance  time in seconds:    10.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2233.517       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.73E-01  1.07E+00  1.50E+00  1.01E+00  1.24E+00  1.69E+00  2.87E+00  1.59E+00  9.82E-01  1.36E+00  9.59E-01
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.94E-02  1.76E-01  1.96E-01  9.82E-02  9.66E-02  1.35E-01  3.34E-01  2.69E-01  1.82E-01  1.27E-01  4.26E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.44E-03
 
 TH 2
+        2.95E-03  3.10E-02
 
 TH 3
+        1.67E-03 -1.80E-02  3.82E-02
 
 TH 4
+        6.20E-04 -1.30E-02  1.16E-02  9.64E-03
 
 TH 5
+        4.59E-04  9.29E-03 -5.06E-04 -5.71E-03  9.32E-03
 
 TH 6
+        1.99E-03  3.59E-04  3.50E-03  1.78E-03  9.55E-04  1.83E-02
 
 TH 7
+        2.00E-03 -3.96E-02  3.79E-02  1.99E-02 -5.84E-03  2.35E-02  1.12E-01
 
 TH 8
+        1.44E-03 -1.97E-03  1.89E-02 -1.06E-05  7.82E-03  7.25E-04  1.69E-02  7.23E-02
 
 TH 9
+        3.54E-04  3.40E-03 -1.23E-03  3.26E-03 -5.01E-03 -1.94E-04 -1.76E-02 -2.15E-02  3.31E-02
 
 TH10
+       -8.38E-04  5.29E-03 -1.47E-04 -5.02E-03  6.30E-03  2.04E-03  1.17E-03  7.80E-03 -5.96E-03  1.61E-02
 
 TH11
+       -7.22E-05  2.51E-03 -1.91E-03 -1.67E-03  7.43E-04  6.10E-05 -4.79E-03 -8.81E-07 -3.43E-04  6.57E-04  1.82E-03
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        4.94E-02
 
 TH 2
+        3.39E-01  1.76E-01
 
 TH 3
+        1.73E-01 -5.23E-01  1.96E-01
 
 TH 4
+        1.28E-01 -7.50E-01  6.02E-01  9.82E-02
 
 TH 5
+        9.62E-02  5.46E-01 -2.68E-02 -6.03E-01  9.66E-02
 
 TH 6
+        2.98E-01  1.51E-02  1.32E-01  1.34E-01  7.31E-02  1.35E-01
 
 TH 7
+        1.21E-01 -6.72E-01  5.79E-01  6.05E-01 -1.81E-01  5.19E-01  3.34E-01
 
 TH 8
+        1.09E-01 -4.15E-02  3.59E-01 -4.00E-04  3.01E-01  1.99E-02  1.87E-01  2.69E-01
 
 TH 9
+        3.94E-02  1.06E-01 -3.47E-02  1.82E-01 -2.85E-01 -7.89E-03 -2.89E-01 -4.39E-01  1.82E-01
 
 TH10
+       -1.34E-01  2.37E-01 -5.91E-03 -4.03E-01  5.14E-01  1.19E-01  2.76E-02  2.29E-01 -2.58E-01  1.27E-01
 
 TH11
+       -3.43E-02  3.35E-01 -2.29E-01 -4.00E-01  1.80E-01  1.06E-02 -3.36E-01 -7.68E-05 -4.42E-02  1.21E-01  4.26E-02
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.19E+03
 
 TH 2
+       -4.87E+02  3.50E+02
 
 TH 3
+       -7.40E+01  4.24E+01  7.89E+01
 
 TH 4
+       -3.28E+02  1.98E+02 -6.60E+01  6.25E+02
 
 TH 5
+        1.08E+02 -1.23E+02 -7.02E+01  1.36E+02  3.21E+02
 
 TH 6
+        5.33E+01 -9.17E+01  1.32E+01 -2.35E+01  8.64E+00  1.38E+02
 
 TH 7
+       -1.12E+02  9.41E+01 -7.33E+00 -3.89E+00 -2.43E+01 -6.90E+01  6.60E+01
 
 TH 8
+       -2.38E+00 -3.51E+00 -1.35E+01 -5.58E+00 -6.61E+00 -6.44E-01  3.21E+00  2.21E+01
 
 TH 9
+        3.96E+01 -2.60E+01 -1.95E+01 -5.79E+01  2.25E+01 -2.32E+01  2.16E+01  1.48E+01  6.19E+01
 
 TH10
+        9.47E+01 -3.43E+01 -1.08E+01  3.85E+01 -2.73E+01  4.63E+00 -2.02E+01 -3.65E+00 -3.01E-01  1.03E+02
 
 TH11
+       -2.53E+01  3.25E+01 -3.02E+01  1.30E+02  4.47E+01 -7.48E+01  5.16E+01  7.63E-01  2.39E+01 -4.32E+00  7.18E+02
 
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
 #CPUT: Total CPU Time in Seconds,       63.137
Stop Time:
Thu Sep 30 08:34:30 CDT 2021
