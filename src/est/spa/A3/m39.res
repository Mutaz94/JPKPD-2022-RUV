Sat Sep 25 09:14:38 CDT 2021
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
$DATA ../../../../data/spa/A3/dat39.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1390.08386530167        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1895E+02  8.0486E+01  8.3834E+01 -1.5671E+01  8.3119E+01  2.4579E+01 -7.2222E+01 -3.7110E+01 -2.0081E+02 -1.4541E+02
            -5.6497E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -817.836606319106        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4423E+00  1.2929E+00  1.0164E+00  1.5660E+00  9.5839E-01  7.0853E-01  8.7890E-01  9.7003E-01  8.3358E-01  9.6614E-01
             1.5932E+01
 PARAMETER:  4.6621E-01  3.5693E-01  1.1624E-01  5.4853E-01  5.7495E-02 -2.4456E-01 -2.9089E-02  6.9575E-02 -8.2027E-02  6.5552E-02
             2.8683E+00
 GRADIENT:   1.1656E+02 -3.0233E+01 -1.0843E+01 -3.8504E+01  6.2642E+00 -4.7174E+00  6.5875E+00  3.6330E+00  1.6097E+01  6.6699E+00
             3.8531E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1032.52096865587        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0176E+00  1.1541E+00  1.3496E+00  1.3781E+00  2.9429E+00  6.1102E-01  1.0000E-02  2.4565E-01  1.0222E+00  6.4730E-01
             8.7881E+00
 PARAMETER:  1.1745E-01  2.4331E-01  3.9980E-01  4.2072E-01  1.1794E+00 -3.9262E-01 -5.6434E+00 -1.3038E+00  1.2196E-01 -3.3495E-01
             2.2734E+00
 GRADIENT:  -2.7849E+02  2.9747E+01 -3.9557E+00  8.6706E+01 -3.0129E+00 -4.3364E+01  0.0000E+00  6.7025E-02  3.0913E+01  1.1639E-01
             3.4850E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1192.75329591086        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.4871E-01  1.4349E+00  1.2459E+00  8.9125E-01  1.3089E+00  9.6140E-01  1.0000E-02  1.5862E-01  7.4361E-01  2.6641E-01
             4.8856E+00
 PARAMETER:  4.7343E-02  4.6111E-01  3.1988E-01 -1.5134E-02  3.6918E-01  6.0632E-02 -4.8303E+00 -1.7412E+00 -1.9623E-01 -1.2227E+00
             1.6863E+00
 GRADIENT:  -1.7086E+02  1.0807E+02  1.0077E+01  7.7888E+01 -4.1301E+01  7.3201E+00  0.0000E+00  3.9988E-02 -3.7591E+00  5.2269E-01
             1.8498E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1203.64854777102        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0169E+00  1.3298E+00  1.2530E+00  9.1962E-01  1.3757E+00  8.8353E-01  1.1557E-02  1.0303E-01  8.0901E-01  2.3026E-01
             4.7320E+00
 PARAMETER:  1.1675E-01  3.8500E-01  3.2553E-01  1.6203E-02  4.1898E-01 -2.3827E-02 -4.3605E+00 -2.1727E+00 -1.1195E-01 -1.3686E+00
             1.6544E+00
 GRADIENT:   1.6948E+00  1.4900E+00  6.6548E-01  1.4005E+00 -9.0437E-01 -1.6198E+00 -5.4842E-04  2.3298E-02  3.8131E-01  4.6247E-01
            -8.6134E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1203.92997045796        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0188E+00  1.3906E+00  1.0519E+00  8.7859E-01  1.3214E+00  8.9089E-01  3.0515E-02  1.4761E-02  8.5482E-01  5.3425E-02
             4.7403E+00
 PARAMETER:  1.1862E-01  4.2972E-01  1.5058E-01 -2.9433E-02  3.7869E-01 -1.5539E-02 -3.3895E+00 -4.1158E+00 -5.6863E-02 -2.8295E+00
             1.6561E+00
 GRADIENT:   2.4622E+00  4.1100E+00 -1.7225E-01  3.3982E+00 -6.3046E-01  3.7878E-01 -6.2768E-03  6.6543E-04 -3.1712E-01  2.6833E-02
             1.6906E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1203.96114247769        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  1.0161E+00  1.4235E+00  1.0846E+00  8.5191E-01  1.3558E+00  8.8817E-01  3.9675E-02  1.0845E-02  8.7146E-01  2.9175E-02
             4.7334E+00
 PARAMETER:  1.1595E-01  4.5310E-01  1.8121E-01 -6.0275E-02  4.0439E-01 -1.8588E-02 -3.1270E+00 -4.4241E+00 -3.7588E-02 -3.4345E+00
             1.6547E+00
 GRADIENT:  -1.2393E+00 -8.5583E-01  9.9356E-02 -1.3001E+00 -1.3526E-02 -1.7021E-01 -1.0480E-02  3.0168E-04 -6.3957E-02  7.6747E-03
            -3.0794E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1203.96130052250        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      519
 NPARAMETR:  1.0165E+00  1.4170E+00  1.0738E+00  8.5645E-01  1.3468E+00  8.8854E-01  3.9800E-02  1.0718E-02  8.7120E-01  3.0622E-02
             4.7340E+00
 PARAMETER:  1.1640E-01  4.4857E-01  1.7124E-01 -5.4963E-02  3.9773E-01 -1.8171E-02 -3.1239E+00 -4.4359E+00 -3.7886E-02 -3.3860E+00
             1.6548E+00
 GRADIENT:  -4.9301E-01 -5.4100E-01  1.5592E-02 -8.9711E-01 -1.3217E-02 -6.6320E-02 -1.0243E-02  3.0844E-04 -3.6974E-02  8.5825E-03
            -9.4010E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1203.96180043535        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.0172E+00  1.4109E+00  1.0591E+00  8.6086E-01  1.3360E+00  8.8904E-01  4.4550E-02  1.0000E-02  8.7163E-01  2.9037E-02
             4.7349E+00
 PARAMETER:  1.1702E-01  4.4425E-01  1.5743E-01 -4.9826E-02  3.8966E-01 -1.7612E-02 -3.0111E+00 -4.5420E+00 -3.7393E-02 -3.4392E+00
             1.6550E+00
 GRADIENT:   5.0811E-01 -1.5534E-01 -1.0124E-01 -3.8634E-01 -4.9589E-03  6.9318E-02 -1.2439E-02  0.0000E+00 -3.7512E-03  7.8649E-03
             1.8867E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1204.06255094021        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  1.0211E+00  1.4311E+00  9.7890E-01  8.5166E-01  1.3047E+00  8.9160E-01  2.7421E-01  1.0000E-02  8.6128E-01  1.0000E-02
             4.7415E+00
 PARAMETER:  1.2088E-01  4.5842E-01  7.8675E-02 -6.0568E-02  3.6595E-01 -1.4732E-02 -1.1939E+00 -6.1291E+00 -4.9336E-02 -5.0825E+00
             1.6564E+00
 GRADIENT:   7.0432E+00 -3.2135E-01 -7.4573E-01  1.5203E+00  8.1431E-01  1.2670E+00  1.1141E-01  0.0000E+00  1.1221E+00  0.0000E+00
             4.1088E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1204.19680187608        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  1.0170E+00  1.5225E+00  1.0260E+00  7.9149E-01  1.3865E+00  8.8842E-01  3.1832E-01  1.0000E-02  8.3199E-01  1.0000E-02
             4.7315E+00
 PARAMETER:  1.1687E-01  5.2032E-01  1.2564E-01 -1.3384E-01  4.2679E-01 -1.8307E-02 -1.0447E+00 -6.4458E+00 -8.3939E-02 -5.8391E+00
             1.6542E+00
 GRADIENT:  -5.6890E-01  6.6060E-01 -8.8795E-02  5.7605E-01  5.3637E-02  1.6289E-01  1.1474E-01  0.0000E+00  1.8607E-01  0.0000E+00
             2.4986E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1204.23679966561        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      876
 NPARAMETR:  1.0181E+00  1.6355E+00  1.0971E+00  7.2512E-01  1.4985E+00  8.8687E-01  3.3235E-01  1.0000E-02  8.2561E-01  1.0000E-02
             4.7537E+00
 PARAMETER:  1.1797E-01  5.9198E-01  1.9267E-01 -2.2142E-01  5.0447E-01 -2.0054E-02 -1.0016E+00 -6.7969E+00 -9.1633E-02 -6.7237E+00
             1.6589E+00
 GRADIENT:  -1.8440E+00  4.9500E+00  3.0666E-01  3.5702E+00 -4.4755E-01 -7.0176E-01  6.0305E-02  0.0000E+00  3.0356E-02  0.0000E+00
             4.8401E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1204.30890722710        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1053
 NPARAMETR:  1.0181E+00  1.8374E+00  9.2440E-01  5.8959E-01  1.5431E+00  8.8756E-01  3.2021E-01  1.0000E-02  9.3360E-01  1.0000E-02
             4.7541E+00
 PARAMETER:  1.1794E-01  7.0837E-01  2.1394E-02 -4.2834E-01  5.3381E-01 -1.9279E-02 -1.0388E+00 -9.8399E+00  3.1288E-02 -9.6242E+00
             1.6590E+00
 GRADIENT:  -3.6496E+00  1.0740E+01  8.1375E-01  4.0262E+00 -3.1725E+00 -7.5150E-01 -2.1412E-01  0.0000E+00 -2.8668E-01  0.0000E+00
            -1.5188E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1204.37070710174        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1231
 NPARAMETR:  1.0195E+00  1.9179E+00  7.8342E-01  5.3391E-01  1.5564E+00  8.9022E-01  3.1838E-01  1.0000E-02  1.0271E+00  1.0000E-02
             4.7558E+00
 PARAMETER:  1.1933E-01  7.5125E-01 -1.4409E-01 -5.2752E-01  5.4237E-01 -1.6284E-02 -1.0445E+00 -1.1904E+01  1.2679E-01 -1.1215E+01
             1.6594E+00
 GRADIENT:  -7.0311E-01  3.3049E+00  1.4576E-02  1.4911E+00 -4.8995E-01  2.5671E-01  2.3643E-02  0.0000E+00 -1.7804E-03  0.0000E+00
            -4.1425E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1204.40838113380        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1407
 NPARAMETR:  1.0196E+00  2.0736E+00  6.9121E-01  4.2968E-01  1.6336E+00  8.8915E-01  3.1836E-01  1.0000E-02  1.1386E+00  1.0000E-02
             4.7616E+00
 PARAMETER:  1.1939E-01  8.2926E-01 -2.6931E-01 -7.4472E-01  5.9080E-01 -1.7490E-02 -1.0446E+00 -1.5010E+01  2.2980E-01 -1.4446E+01
             1.6606E+00
 GRADIENT:  -6.8248E-01  3.1594E+00  1.6473E-01  8.4456E-01 -7.0888E-01 -1.5121E-01 -6.1637E-02  0.0000E+00 -2.1913E-02  0.0000E+00
            -2.3708E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1204.43281937488        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1582
 NPARAMETR:  1.0199E+00  2.2264E+00  5.3221E-01  3.2763E-01  1.6937E+00  8.8965E-01  3.2137E-01  1.0000E-02  1.3353E+00  1.0000E-02
             4.7656E+00
 PARAMETER:  1.1969E-01  9.0041E-01 -5.3071E-01 -1.0159E+00  6.2694E-01 -1.6929E-02 -1.0352E+00 -1.9664E+01  3.8918E-01 -1.8787E+01
             1.6614E+00
 GRADIENT:  -9.1843E-01  2.0669E+00 -2.6595E-02  5.9467E-01  2.2044E-03 -1.0411E-02  7.9744E-02  0.0000E+00  1.0133E-01  0.0000E+00
             1.7106E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1204.44260329924        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1759
 NPARAMETR:  1.0201E+00  2.3025E+00  4.7161E-01  2.7596E-01  1.7270E+00  8.8992E-01  3.2456E-01  1.0000E-02  1.3868E+00  1.0000E-02
             4.7685E+00
 PARAMETER:  1.1985E-01  9.3398E-01 -6.5161E-01 -1.1875E+00  6.4636E-01 -1.6628E-02 -1.0253E+00 -2.2316E+01  4.2698E-01 -2.1490E+01
             1.6620E+00
 GRADIENT:  -2.8244E-01  1.8952E+00  1.9975E-02  2.9147E-01 -3.3012E-01  1.6591E-02 -1.0943E-02  0.0000E+00 -8.7123E-02  0.0000E+00
            -2.8018E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1204.45105089228        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1935
 NPARAMETR:  1.0203E+00  2.3791E+00  3.9042E-01  2.2365E-01  1.7576E+00  8.8977E-01  3.2244E-01  1.0000E-02  1.6068E+00  1.0000E-02
             4.7705E+00
 PARAMETER:  1.2005E-01  9.6674E-01 -8.4054E-01 -1.3977E+00  6.6396E-01 -1.6794E-02 -1.0319E+00 -2.5997E+01  5.7422E-01 -2.5013E+01
             1.6624E+00
 GRADIENT:   1.3168E-01  1.4871E-01  3.2418E-02  5.0707E-02 -1.1935E-01  3.5975E-02  5.8093E-02  0.0000E+00  7.3308E-02  0.0000E+00
             1.5979E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1204.45619919778        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2112
 NPARAMETR:  1.0203E+00  2.4389E+00  3.3356E-01  1.8357E-01  1.7850E+00  8.8979E-01  3.2264E-01  1.0000E-02  1.7198E+00  1.0000E-02
             4.7728E+00
 PARAMETER:  1.2005E-01  9.9155E-01 -9.9792E-01 -1.5952E+00  6.7944E-01 -1.6764E-02 -1.0312E+00 -2.9256E+01  6.4219E-01 -2.8251E+01
             1.6629E+00
 GRADIENT:   7.1805E-03  6.5159E-01  3.3337E-02  4.3798E-02 -1.9458E-01 -4.5178E-02 -3.9846E-02  0.0000E+00 -1.7840E-02  0.0000E+00
            -1.3860E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1204.46051548295        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2287
 NPARAMETR:  1.0204E+00  2.5031E+00  2.4857E-01  1.4074E-01  1.8075E+00  8.9001E-01  3.2182E-01  1.0000E-02  1.9924E+00  1.0000E-02
             4.7745E+00
 PARAMETER:  1.2017E-01  1.0175E+00 -1.2920E+00 -1.8608E+00  6.9195E-01 -1.6517E-02 -1.0338E+00 -3.4231E+01  7.8935E-01 -3.2793E+01
             1.6633E+00
 GRADIENT:  -2.4674E-01  1.3682E+00 -1.6246E-02  1.3421E-01 -1.1078E-01 -4.6874E-02 -6.7907E-02  0.0000E+00 -8.4663E-03  0.0000E+00
            -2.1688E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1204.46283471912        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2464
 NPARAMETR:  1.0204E+00  2.5436E+00  2.0326E-01  1.1366E-01  1.8261E+00  8.9024E-01  3.2314E-01  1.0000E-02  2.1787E+00  1.0000E-02
             4.7755E+00
 PARAMETER:  1.2018E-01  1.0336E+00 -1.4933E+00 -2.0746E+00  7.0221E-01 -1.6259E-02 -1.0297E+00 -3.8015E+01  8.7873E-01 -3.6413E+01
             1.6635E+00
 GRADIENT:  -3.5466E-01  1.2481E+00 -2.4183E-02  1.0471E-01  1.4407E-02  2.5202E-02 -1.9323E-02  0.0000E+00 -2.5349E-02  0.0000E+00
            -1.7023E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1204.46546156463        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2640
 NPARAMETR:  1.0204E+00  2.5803E+00  1.6786E-01  8.8344E-02  1.8417E+00  8.8989E-01  3.2338E-01  1.0000E-02  2.4635E+00  1.0000E-02
             4.7767E+00
 PARAMETER:  1.2023E-01  1.0479E+00 -1.6846E+00 -2.3265E+00  7.1071E-01 -1.6652E-02 -1.0289E+00 -4.2273E+01  1.0016E+00 -4.0707E+01
             1.6637E+00
 GRADIENT:  -9.6897E-02  6.1436E-01  1.1005E-02  2.3510E-02 -7.5346E-02 -5.5990E-02  6.4509E-03  0.0000E+00 -6.0586E-03  0.0000E+00
            -3.8172E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1204.46695057035        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2817
 NPARAMETR:  1.0206E+00  2.6142E+00  1.2109E-01  6.4628E-02  1.8488E+00  8.8962E-01  3.2366E-01  1.0000E-02  2.9294E+00  1.0000E-02
             4.7780E+00
 PARAMETER:  1.2043E-01  1.0609E+00 -2.0112E+00 -2.6391E+00  7.1453E-01 -1.6961E-02 -1.0281E+00 -4.8113E+01  1.1748E+00 -4.6160E+01
             1.6640E+00
 GRADIENT:   4.8836E-01 -4.4121E-01  6.9946E-03 -1.6866E-02 -1.9788E-01 -1.0155E-01  5.3491E-02  0.0000E+00  8.9628E-03  0.0000E+00
             2.4073E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1204.46811511157        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2992
 NPARAMETR:  1.0205E+00  2.6336E+00  9.6684E-02  5.2465E-02  1.8604E+00  8.9013E-01  3.2282E-01  1.0000E-02  3.2263E+00  1.0000E-02
             4.7781E+00
 PARAMETER:  1.2032E-01  1.0684E+00 -2.2363E+00 -2.8476E+00  7.2077E-01 -1.6393E-02 -1.0307E+00 -5.2018E+01  1.2714E+00 -4.9771E+01
             1.6640E+00
 GRADIENT:  -1.7573E-01  1.1434E+00 -5.7567E-03  3.9800E-02 -1.0744E-01 -2.0328E-02 -1.4526E-02  0.0000E+00 -4.0861E-03  0.0000E+00
            -8.6213E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1204.46958384683        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3169
 NPARAMETR:  1.0205E+00  2.6627E+00  6.3213E-02  3.2763E-02  1.8734E+00  8.9010E-01  3.2263E-01  1.0000E-02  4.1017E+00  1.0000E-02
             4.7789E+00
 PARAMETER:  1.2032E-01  1.0793E+00 -2.6612E+00 -3.3184E+00  7.2774E-01 -1.6421E-02 -1.0312E+00 -6.0447E+01  1.5114E+00 -5.7954E+01
             1.6642E+00
 GRADIENT:  -1.9418E-01  8.9676E-01  3.6669E-03  1.6596E-02 -6.6331E-02 -2.9800E-02 -9.3409E-03  0.0000E+00  3.9996E-03  0.0000E+00
            -5.7379E-02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1204.47019084112        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3346
 NPARAMETR:  1.0205E+00  2.6774E+00  4.4974E-02  2.2731E-02  1.8797E+00  8.9015E-01  3.2255E-01  1.0000E-02  4.7697E+00  1.0000E-02
             4.7792E+00
 PARAMETER:  1.2033E-01  1.0849E+00 -3.0017E+00 -3.6840E+00  7.3109E-01 -1.6367E-02 -1.0315E+00 -6.7004E+01  1.6623E+00 -6.4285E+01
             1.6643E+00
 GRADIENT:  -1.7059E-01  7.5750E-01  2.5132E-03  5.9546E-03 -4.3124E-02 -2.4448E-02 -2.1635E-02  0.0000E+00 -5.2274E-03  0.0000E+00
            -9.4327E-02

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1204.47071303096        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3524
 NPARAMETR:  1.0206E+00  2.6902E+00  2.7464E-02  1.4156E-02  1.8841E+00  8.9022E-01  3.2243E-01  1.0000E-02  5.9763E+00  1.0000E-02
             4.7796E+00
 PARAMETER:  1.2037E-01  1.0896E+00 -3.4949E+00 -4.1576E+00  7.3344E-01 -1.6287E-02 -1.0319E+00 -7.5853E+01  1.8878E+00 -7.2580E+01
             1.6643E+00
 GRADIENT:  -1.5612E-01  8.9378E-01 -5.3774E-04  5.4104E-03 -5.6252E-02 -1.0239E-02 -2.7256E-02  0.0000E+00 -6.0185E-03  0.0000E+00
            -1.1450E-01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1204.47110624893        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3699
 NPARAMETR:  1.0206E+00  2.6954E+00  1.9568E-02  1.0340E-02  1.8856E+00  8.9018E-01  3.2257E-01  1.0000E-02  7.4087E+00  1.0000E-02
             4.7797E+00
 PARAMETER:  1.2038E-01  1.0916E+00 -3.8339E+00 -4.4718E+00  7.3422E-01 -1.6331E-02 -1.0314E+00 -8.1908E+01  2.1027E+00 -7.8192E+01
             1.6644E+00
 GRADIENT:  -5.7936E-02  2.5254E-01 -3.8997E-04  4.1274E-03 -8.3230E-03  8.5511E-04  4.7590E-03  0.0000E+00  3.0489E-03  0.0000E+00
             3.2738E-03

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1204.47116555715        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3883
 NPARAMETR:  1.0206E+00  2.6958E+00  1.9468E-02  1.0000E-02  1.8857E+00  8.9018E-01  3.2252E-01  1.0000E-02  7.4024E+00  1.0000E-02
             4.7798E+00
 PARAMETER:  1.2038E-01  1.0917E+00 -3.8390E+00 -4.5174E+00  7.3427E-01 -1.6335E-02 -1.0316E+00 -8.2474E+01  2.1018E+00 -7.8721E+01
             1.6644E+00
 GRADIENT:  -1.6805E-02  1.0169E-01  8.0484E-04  0.0000E+00 -1.0161E-02  1.7120E-03  3.8632E-04  0.0000E+00  1.0504E-03  0.0000E+00
             8.3368E-03

0ITERATION NO.:  144    OBJECTIVE VALUE:  -1204.47117371332        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     4011
 NPARAMETR:  1.0206E+00  2.6957E+00  1.9240E-02  1.0000E-02  1.8857E+00  8.9017E-01  3.2252E-01  1.0000E-02  7.3606E+00  1.0000E-02
             4.7797E+00
 PARAMETER:  1.2038E-01  1.0917E+00 -3.8508E+00 -4.5174E+00  7.3430E-01 -1.6341E-02 -1.0316E+00 -8.2474E+01  2.0961E+00 -7.8721E+01
             1.6644E+00
 GRADIENT:   1.2930E-02 -2.3718E-02 -1.1168E-04  0.0000E+00  7.8272E-03  2.9678E-04  2.5976E-04  0.0000E+00 -1.1195E-04  0.0000E+00
             2.3097E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     4011
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.7195E-03 -1.2422E-02  3.9162E-07 -6.6170E-04 -1.1724E-05
 SE:             2.8265E-02  1.4775E-02  2.1313E-07  1.6870E-03  5.2525E-05
 N:                     100         100         100         100         100

 P VAL.:         8.6739E-01  4.0050E-01  6.6138E-02  6.9488E-01  8.2338E-01

 ETASHRINKSD(%)  5.3093E+00  5.0502E+01  9.9999E+01  9.4348E+01  9.9824E+01
 ETASHRINKVR(%)  1.0337E+01  7.5500E+01  1.0000E+02  9.9681E+01  1.0000E+02
 EBVSHRINKSD(%)  5.4771E+00  5.0840E+01  9.9999E+01  9.4381E+01  9.9759E+01
 EBVSHRINKVR(%)  1.0654E+01  7.5833E+01  1.0000E+02  9.9684E+01  9.9999E+01
 RELATIVEINF(%)  8.0447E+01  9.5856E-02  7.7058E-10  1.5344E-03  2.3459E-05
 EPSSHRINKSD(%)  1.3057E+01
 EPSSHRINKVR(%)  2.4409E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1204.4711737133196     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -469.32034714958138     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    51.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1204.471       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.70E+00  1.92E-02  1.00E-02  1.89E+00  8.90E-01  3.23E-01  1.00E-02  7.36E+00  1.00E-02  4.78E+00
 


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
+        1.20E+03
 
 TH 2
+       -1.49E+02  3.40E+02
 
 TH 3
+        7.70E-01  1.58E+01  8.30E+01
 
 TH 4
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 5
+        8.60E+00 -4.66E+01 -7.25E+00  0.00E+00  1.82E+01
 
 TH 6
+        1.13E+01 -2.07E+01  1.51E+00  0.00E+00  4.96E-01  2.00E+02
 
 TH 7
+       -3.86E+00 -3.84E+01 -1.37E+00  0.00E+00  9.11E+00  1.98E+01  1.18E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.54E-02 -1.36E-02  1.05E-01  0.00E+00  2.45E-03  1.65E-02  3.85E-02  0.00E+00  1.23E-03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.74E+01 -1.46E+01 -4.05E-01  0.00E+00  1.72E+00  3.12E+00  1.78E+01  0.00E+00  9.52E-03  0.00E+00  2.36E+01
 
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
 #CPUT: Total CPU Time in Seconds,       58.776
Stop Time:
Sat Sep 25 09:15:39 CDT 2021
