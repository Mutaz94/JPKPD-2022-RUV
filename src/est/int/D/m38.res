Sat Sep 25 05:43:39 CDT 2021
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
$DATA ../../../../data/int/D/dat38.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   18602.4764594706        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.8578E+01  2.2684E+02  8.9718E+01  4.3505E+01  2.6701E+02 -2.7063E+03 -1.0880E+03 -2.8907E+02 -1.4951E+03 -8.1976E+02
            -3.8698E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1035.05140137545        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.9359E+00  1.2770E+00  7.1704E-01  2.2970E+00  1.1261E+00  8.2228E+00  4.3362E+00  1.1576E+00  4.7446E+00  2.6191E+00
             9.8657E+00
 PARAMETER:  7.6056E-01  3.4451E-01 -2.3262E-01  9.3160E-01  2.1875E-01  2.2069E+00  1.5670E+00  2.4637E-01  1.6570E+00  1.0628E+00
             2.3891E+00
 GRADIENT:   1.6293E+01 -4.4551E+01 -4.6340E+01  7.2923E+01 -7.6224E+00  1.7815E+02 -1.5071E+00  5.1869E+00  8.2092E+01  7.2849E+01
             5.5045E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1112.18129025498        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  2.6930E+00  6.3421E+00  5.1368E+00  1.6393E+00  1.4099E+01  5.3889E+00  3.4121E+00  3.0264E+00  2.2255E+01  2.3682E+00
             9.8095E+00
 PARAMETER:  1.0906E+00  1.9472E+00  1.7364E+00  5.9424E-01  2.7461E+00  1.7843E+00  1.3273E+00  1.2074E+00  3.2026E+00  9.6214E-01
             2.3834E+00
 GRADIENT:   5.7977E+01  8.7042E+01  6.4429E+00  1.8015E+01  1.0440E+01  1.0517E+02 -1.3061E+01  3.1450E-01  8.3013E+00 -1.9701E+00
             4.1601E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1214.39416730966        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.3475E+00  3.5820E+00  1.0027E+01  2.4485E+00  6.5501E+00  3.4194E+00  4.0431E+00  6.2405E+01  1.2810E+01  2.1199E+00
             9.4091E+00
 PARAMETER:  3.9825E-01  1.3759E+00  2.4053E+00  9.9547E-01  1.9795E+00  1.3295E+00  1.4970E+00  4.2336E+00  2.6502E+00  8.5139E-01
             2.3417E+00
 GRADIENT:   1.5157E+01 -1.8550E-01 -2.0375E+00  3.0704E+01  4.1073E+01  8.4589E+01  7.5731E+01  1.8507E+00  5.1646E+01  4.5832E+01
             4.0774E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1370.64526752288        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.9966E-01  1.4826E+00  3.9683E+00  9.2269E-01  2.3052E+00  1.8889E+00  2.5977E+00  2.9434E+00  5.7739E+00  1.0857E+00
             7.4990E+00
 PARAMETER:  9.9659E-02  4.9378E-01  1.4783E+00  1.9533E-02  9.3518E-01  7.3601E-01  1.0546E+00  1.1796E+00  1.8533E+00  1.8227E-01
             2.1148E+00
 GRADIENT:  -7.1523E+01 -9.1188E+01 -5.9733E+01  1.3397E+01  6.1814E+01 -1.0500E+02 -7.8934E+00 -4.9574E+01  4.5427E+01  1.8448E+01
             2.2127E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1507.97132985057        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.2254E+00  2.1829E+00  9.1274E+01  7.3873E-01  2.4542E+00  2.5036E+00  1.8316E+00  5.6596E+00  5.8414E+00  9.9088E-01
             6.1259E+00
 PARAMETER:  3.0324E-01  8.8063E-01  4.6139E+00 -2.0283E-01  9.9779E-01  1.0177E+00  7.0517E-01  1.8334E+00  1.8650E+00  9.0838E-02
             1.9125E+00
 GRADIENT:   3.2860E+01  4.8309E+01 -5.4263E+00  1.7866E+01  2.5292E+01 -1.4882E+01  2.9991E+00  2.8219E+00 -5.7022E+00  1.3126E+01
            -1.6580E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1527.52259570099        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.1141E+00  1.6075E+00  1.4253E+03  7.3955E-01  2.3505E+00  2.5164E+00  9.8070E-01  3.1684E+00  5.5989E+00  9.1921E-01
             6.6198E+00
 PARAMETER:  2.0804E-01  5.7467E-01  7.3622E+00 -2.0171E-01  9.5464E-01  1.0228E+00  8.0513E-02  1.2532E+00  1.8226E+00  1.5756E-02
             1.9901E+00
 GRADIENT:  -5.0814E+00 -7.1916E+00 -2.7652E-01 -1.1226E+00 -8.4407E+00 -5.7982E+00  2.2501E+00  2.8781E-01  6.5842E+00  1.2656E+01
             6.8283E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1528.67339698633        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      564
 NPARAMETR:  1.1309E+00  1.5805E+00  4.5102E+03  7.5937E-01  2.3782E+00  2.5554E+00  7.8628E-01  2.6485E+00  5.7456E+00  8.8637E-01
             6.6175E+00
 PARAMETER:  2.2301E-01  5.5771E-01  8.5141E+00 -1.7526E-01  9.6634E-01  1.0382E+00 -1.4044E-01  1.0740E+00  1.8484E+00 -2.0624E-02
             1.9897E+00
 GRADIENT:  -2.4559E+00  1.2369E-01  1.0534E-01  1.9305E-01 -6.5100E+00 -5.9954E+00 -1.6361E+00  1.6299E-01 -8.1755E+00  1.1018E+01
            -9.9102E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1532.25804817757        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      746
 NPARAMETR:  1.1920E+00  2.4704E+00  8.0516E+02  2.4019E-01  2.3541E+00  2.6427E+00  1.4934E+00  2.5812E+00  8.4641E+00  8.1614E-01
             6.6410E+00
 PARAMETER:  2.7564E-01  1.0044E+00  6.7910E+00 -1.3263E+00  9.5617E-01  1.0718E+00  5.0104E-01  1.0483E+00  2.2358E+00 -1.0317E-01
             1.9933E+00
 GRADIENT:   1.3329E+01  3.1281E+01 -5.6752E-01  1.6600E+00 -1.4022E+01  5.8928E+00 -8.3953E+00 -1.6922E-01 -4.1604E+00  8.9294E+00
             2.3071E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1535.29237133907        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  1.1496E+00  2.4249E+00  6.9330E+03  1.9764E-01  2.4118E+00  2.5992E+00  1.5717E+00  1.8056E+00  8.8968E+00  6.2609E-01
             6.6325E+00
 PARAMETER:  2.3938E-01  9.8579E-01  8.9440E+00 -1.5213E+00  9.8039E-01  1.0552E+00  5.5216E-01  6.9089E-01  2.2857E+00 -3.6825E-01
             1.9920E+00
 GRADIENT:   2.3683E+00  2.5617E+00  2.3969E-01 -7.7897E-01 -2.7546E+00  2.3854E-02  2.2582E+00  4.1569E-02  1.3353E+00  4.1639E+00
             4.8301E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1536.08976562788        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  1.1361E+00  2.3345E+00  8.5444E+05  2.9784E-01  2.4472E+00  2.5719E+00  1.4997E+00  7.8687E-01  7.8409E+00  3.6566E-01
             6.6426E+00
 PARAMETER:  2.2760E-01  9.4779E-01  1.3758E+01 -1.1112E+00  9.9494E-01  1.0446E+00  5.0526E-01 -1.3969E-01  2.1594E+00 -9.0605E-01
             1.9935E+00
 GRADIENT:  -8.7340E+01  2.1673E+01  3.8607E-03  3.6251E+00 -3.9137E+00  2.3541E+00 -2.1519E+00 -1.5477E-04 -2.7596E+00  9.2902E-01
             7.1080E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1536.25808641518        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1273
 NPARAMETR:  1.1410E+00  2.2927E+00  8.4928E+05  2.9651E-01  2.4472E+00  2.5964E+00  1.4964E+00  7.8747E-01  7.9441E+00  3.6547E-01
             6.6566E+00
 PARAMETER:  2.3188E-01  9.2972E-01  1.3752E+01 -1.1157E+00  9.9495E-01  1.0541E+00  5.0309E-01 -1.3893E-01  2.1724E+00 -9.0656E-01
             1.9956E+00
 GRADIENT:  -2.3936E+02  7.7632E+01  3.7711E-03  1.6380E+01 -4.9274E+00  5.8536E+01 -2.2549E+01  1.3494E-03  4.0769E+01  1.0461E+00
             3.5167E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1536.27561032696        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1473             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1310E+00  2.2989E+00  1.0160E+06  2.9200E-01  2.4425E+00  2.5945E+00  1.4921E+00  7.8734E-01  7.9988E+00  3.6213E-01
             6.6513E+00
 PARAMETER:  2.2312E-01  9.3243E-01  1.3931E+01 -1.1310E+00  9.9301E-01  1.0534E+00  5.0017E-01 -1.3909E-01  2.1793E+00 -9.1574E-01
             1.9948E+00
 GRADIENT:  -1.1233E+03 -2.7744E+02  2.6426E-03 -3.8979E+02 -1.3587E+00  3.2236E+02 -2.2712E+01 -4.4413E-03 -7.5834E+02  3.9660E-01
             3.6014E+02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1536.28514996074        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1573
 NPARAMETR:  1.1395E+00  2.2952E+00  1.0160E+06  2.9157E-01  2.4431E+00  2.5953E+00  1.4921E+00  7.8734E-01  7.9761E+00  3.6212E-01
             6.6513E+00
 PARAMETER:  2.3063E-01  9.3083E-01  1.3931E+01 -1.1325E+00  9.9327E-01  1.0537E+00  5.0020E-01 -1.3909E-01  2.1765E+00 -9.1579E-01
             1.9948E+00
 GRADIENT:  -6.6257E+02  1.1012E+02  3.0454E-03  6.2644E+01 -6.7523E+00  4.5515E+01 -7.1299E+00  1.0774E-03  9.5479E+01  1.0892E+00
             9.1704E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1536.30799294467        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1740
 NPARAMETR:  1.1397E+00  2.2958E+00  1.0107E+06  2.8649E-01  2.4482E+00  2.5969E+00  1.4907E+00  7.8653E-01  8.0028E+00  3.6021E-01
             6.6484E+00
 PARAMETER:  2.3073E-01  9.3108E-01  1.3926E+01 -1.1501E+00  9.9535E-01  1.0543E+00  4.9924E-01 -1.4012E-01  2.1798E+00 -9.2106E-01
             1.9944E+00
 GRADIENT:  -2.3375E+02  2.8220E+01  3.2627E-03  1.4047E+01 -3.2122E+00 -5.3196E+01 -8.6208E-01  1.2013E-04  1.9637E+01  8.8124E-01
             3.4507E+01

0ITERATION NO.:   71    OBJECTIVE VALUE:  -1536.30799294467        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     1765
 NPARAMETR:  1.1397E+00  2.2958E+00  1.0107E+06  2.8649E-01  2.4482E+00  2.5969E+00  1.4907E+00  7.8653E-01  8.0028E+00  3.6021E-01
             6.6484E+00
 PARAMETER:  2.3073E-01  9.3108E-01  1.3926E+01 -1.1501E+00  9.9535E-01  1.0543E+00  4.9924E-01 -1.4012E-01  2.1798E+00 -9.2106E-01
             1.9944E+00
 GRADIENT:   8.9608E-02 -4.8813E-01  5.0088E-02  1.8486E-01 -2.8890E-01  6.0240E-01 -4.9609E-01  4.0351E+00 -7.3476E-02  8.2338E-01
            -1.1007E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1765
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4160E-03 -9.2577E-02 -1.0198E-06  6.0690E-02 -1.1646E-03
 SE:             2.9153E-02  1.6969E-02  7.1735E-07  2.1700E-02  6.2238E-03
 N:                     100         100         100         100         100

 P VAL.:         8.5262E-01  4.8914E-08  1.5516E-01  5.1611E-03  8.5156E-01

 ETASHRINKSD(%)  2.3344E+00  4.3151E+01  9.9998E+01  2.7303E+01  7.9149E+01
 ETASHRINKVR(%)  4.6144E+00  6.7682E+01  1.0000E+02  4.7151E+01  9.5653E+01
 EBVSHRINKSD(%)  2.2806E+00  3.5698E+01  9.9995E+01  2.6407E+01  7.8330E+01
 EBVSHRINKVR(%)  4.5093E+00  5.8652E+01  1.0000E+02  4.5841E+01  9.5304E+01
 RELATIVEINF(%)  9.5421E+01  2.2364E+01  2.3636E-07  2.9520E+01  4.5296E+00
 EPSSHRINKSD(%)  8.6881E+00
 EPSSHRINKVR(%)  1.6621E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1536.3079929446683     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       117.78136682374247     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    63.51
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    20.92
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1536.308       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  2.30E+00  1.01E+06  2.86E-01  2.45E+00  2.60E+00  1.49E+00  7.87E-01  8.00E+00  3.60E-01  6.65E+00
 


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
+        1.46E+03
 
 TH 2
+        1.48E+02  6.89E+01
 
 TH 3
+       -4.20E-05  1.86E-05  1.13E-12
 
 TH 4
+        3.41E+03 -2.27E+02 -2.04E-05  1.46E+03
 
 TH 5
+        1.04E+02 -2.06E+00 -6.26E-06  1.29E+02  6.71E+01
 
 TH 6
+        1.26E+02  4.85E+01 -5.93E-07  3.48E-02  3.88E+01  4.52E+01
 
 TH 7
+       -1.18E+02 -2.85E+01  1.79E-05  8.10E+02 -3.77E+01 -2.21E+01  2.76E+02
 
 TH 8
+        7.17E+02 -5.25E+02  1.49E-04  1.10E+04  3.75E+01 -1.07E+02  1.51E+03  1.38E+04
 
 TH 9
+        3.03E+01 -1.21E+00  2.47E-07 -6.88E+00 -1.44E+00  4.45E+00 -1.36E+01  5.16E+00  2.22E+00
 
 TH10
+       -1.21E+03 -1.80E+02 -1.72E-05 -3.52E+02  1.13E+02 -3.10E+01  6.92E+02 -4.83E+03 -7.21E+01  1.77E+03
 
 TH11
+        2.75E+01 -1.79E-01  3.37E-07  1.48E+01 -8.99E-02  4.91E+00  8.17E+00 -8.13E+01 -5.95E-01  2.51E+01  2.48E+01
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       84.557
Stop Time:
Sat Sep 25 05:45:05 CDT 2021
