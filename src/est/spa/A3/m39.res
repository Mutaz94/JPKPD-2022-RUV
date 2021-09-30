Wed Sep 29 13:30:22 CDT 2021
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
 GRADIENT:   4.6295E+02  1.1048E+02  8.5822E+01  3.8492E+01  1.0594E+02  4.9552E+01 -6.6611E+01 -3.6838E+01 -1.8765E+02 -1.1646E+02
            -5.6439E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -811.644214940940        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3674E+00  1.2141E+00  9.4704E-01  1.4945E+00  9.8798E-01  5.4783E-01  9.2337E-01  9.8733E-01  9.1439E-01  9.8764E-01
             1.5967E+01
 PARAMETER:  4.1292E-01  2.9403E-01  4.5587E-02  5.0178E-01  8.7909E-02 -5.0179E-01  2.0279E-02  8.7248E-02  1.0497E-02  8.7566E-02
             2.8705E+00
 GRADIENT:   1.0502E+02 -4.8467E+01 -1.9630E+01 -5.9209E+01  1.6844E+01 -1.1700E+01  6.8072E+00  3.9161E+00  1.9549E+01  6.3961E+00
             4.1673E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1058.36989470334        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.4051E-01  1.2293E+00  1.6673E+00  1.3293E+00  2.2635E+00  7.1816E-01  1.0000E-02  1.6793E-01  1.4462E+00  1.1279E+00
             7.4090E+00
 PARAMETER:  3.8672E-02  3.0646E-01  6.1118E-01  3.8462E-01  9.1693E-01 -2.3106E-01 -6.8567E+00 -1.6842E+00  4.6892E-01  2.2040E-01
             2.1027E+00
 GRADIENT:  -4.2002E+02  4.9657E+01  1.6663E+00  1.2477E+02  2.2835E-01 -6.5831E+01  0.0000E+00  2.1171E-02  4.3329E+01  3.8987E-01
             3.2293E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1197.26496570341        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.0235E+00  1.0445E+00  1.5609E+00  1.1058E+00  1.4686E+00  8.0325E-01  1.0870E-01  4.4171E-02  3.3028E-01  2.4268E-01
             4.7198E+00
 PARAMETER:  1.2326E-01  1.4357E-01  5.4527E-01  2.0059E-01  4.8430E-01 -1.1909E-01 -2.1192E+00 -3.0197E+00 -1.0078E+00 -1.3160E+00
             1.6518E+00
 GRADIENT:   5.7869E+01 -8.4796E+00 -6.4817E+00 -5.7034E+00  9.2903E+00 -3.9821E+01 -2.7288E-01  3.1927E-03 -5.8150E+00  3.2365E-01
            -3.2212E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1203.27006523425        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0062E+00  1.1297E+00  1.7503E+00  1.0516E+00  1.4428E+00  8.8426E-01  4.0221E-02  3.9932E-02  6.5469E-01  2.3503E-01
             4.6896E+00
 PARAMETER:  1.0619E-01  2.2195E-01  6.5981E-01  1.5029E-01  4.6656E-01 -2.3008E-02 -3.1134E+00 -3.1206E+00 -3.2359E-01 -1.3480E+00
             1.6454E+00
 GRADIENT:   9.3063E-01  2.1367E+00  4.5393E-01  1.9228E+00  9.9785E-02 -2.8451E-01 -4.2535E-04  1.8570E-03 -4.5199E-01  4.2215E-01
             8.3085E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1203.60084693932        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      482
 NPARAMETR:  1.0157E+00  1.1647E+00  1.6582E+00  1.0347E+00  1.4411E+00  8.9478E-01  6.5102E-02  2.9728E-02  6.7171E-01  1.2041E-01
             4.7229E+00
 PARAMETER:  1.1559E-01  2.5245E-01  6.0575E-01  1.3408E-01  4.6541E-01 -1.1180E-02 -2.6318E+00 -3.4157E+00 -2.9793E-01 -2.0168E+00
             1.6524E+00
 GRADIENT:   2.1185E+00  6.3685E-01  3.5568E-01  1.4062E+00 -1.5502E-03  2.0962E+00 -2.0167E-02  1.1054E-03 -4.2761E-01  1.0435E-01
            -7.0664E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1204.02845247634        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      659
 NPARAMETR:  1.0149E+00  1.3272E+00  1.2877E+00  9.2405E-01  1.3832E+00  8.8475E-01  3.4484E-01  1.0000E-02  7.0988E-01  1.0000E-02
             4.7307E+00
 PARAMETER:  1.1481E-01  3.8310E-01  3.5283E-01  2.1015E-02  4.2436E-01 -2.2447E-02 -9.6468E-01 -4.6795E+00 -2.4266E-01 -5.2142E+00
             1.6541E+00
 GRADIENT:  -4.3361E+00  1.9438E+00  4.3831E-01  2.7341E+00 -6.7861E-01 -1.1110E+00  1.4506E-02  0.0000E+00 -3.3751E-02  0.0000E+00
            -2.1186E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1204.16558034767        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      836
 NPARAMETR:  1.0189E+00  1.4914E+00  1.0316E+00  8.1511E-01  1.3682E+00  8.8833E-01  3.0806E-01  1.0000E-02  8.0127E-01  1.0000E-02
             4.7546E+00
 PARAMETER:  1.1875E-01  4.9972E-01  1.3111E-01 -1.0444E-01  4.1352E-01 -1.8411E-02 -1.0775E+00 -5.7506E+00 -1.2156E-01 -7.5393E+00
             1.6591E+00
 GRADIENT:   1.2205E-01  3.1649E+00 -1.1839E-02  2.9906E+00 -9.7127E-01 -3.5462E-01 -9.2423E-02  0.0000E+00 -3.3431E-01  0.0000E+00
             1.6136E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1204.22818234618        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  1.0201E+00  1.6175E+00  9.0396E-01  7.2887E-01  1.3851E+00  8.9159E-01  3.1597E-01  1.0000E-02  8.7795E-01  1.0000E-02
             4.7429E+00
 PARAMETER:  1.1991E-01  5.8091E-01 -9.6898E-04 -2.1626E-01  4.2579E-01 -1.4753E-02 -1.0521E+00 -6.4340E+00 -3.0161E-02 -9.4166E+00
             1.6566E+00
 GRADIENT:   2.4429E+00  7.1942E-01 -2.2569E-01  9.6368E-01 -8.8353E-01  8.4233E-01  4.0077E-02  0.0000E+00 -8.2987E-02  0.0000E+00
            -2.4289E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1204.33565646127        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1191
 NPARAMETR:  1.0184E+00  1.8252E+00  9.2124E-01  5.9677E-01  1.5547E+00  8.8787E-01  3.1479E-01  1.0000E-02  9.4527E-01  1.0000E-02
             4.7548E+00
 PARAMETER:  1.1819E-01  7.0169E-01  1.7969E-02 -4.1622E-01  5.4130E-01 -1.8931E-02 -1.0559E+00 -6.4814E+00  4.3717E-02 -1.1780E+01
             1.6591E+00
 GRADIENT:  -2.2220E+00  3.1691E+00  1.5772E-01  1.7753E+00 -1.7393E-01 -5.3015E-01 -7.4234E-02  0.0000E+00 -5.7721E-02  0.0000E+00
            -4.3813E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1204.35442176487        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1368
 NPARAMETR:  1.0195E+00  1.9666E+00  8.9599E-01  5.0349E-01  1.6335E+00  8.9053E-01  3.2789E-01  1.0000E-02  9.9067E-01  1.0000E-02
             4.7599E+00
 PARAMETER:  1.1935E-01  7.7629E-01 -9.8212E-03 -5.8619E-01  5.9072E-01 -1.5941E-02 -1.0151E+00 -6.7701E+00  9.0631E-02 -1.4042E+01
             1.6602E+00
 GRADIENT:   4.5209E-01  4.4133E+00  4.6721E-01  1.4704E+00 -1.1544E+00  4.5121E-01 -1.1695E-02  0.0000E+00 -1.6656E-02  0.0000E+00
            -8.2252E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1204.40548890179        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1547
 NPARAMETR:  1.0193E+00  2.1324E+00  6.4957E-01  3.9219E-01  1.6620E+00  8.8822E-01  3.2582E-01  1.0000E-02  1.1260E+00  1.0000E-02
             4.7659E+00
 PARAMETER:  1.1912E-01  8.5723E-01 -3.3145E-01 -8.3600E-01  6.0803E-01 -1.8536E-02 -1.0214E+00 -8.3465E+00  2.1867E-01 -1.8588E+01
             1.6615E+00
 GRADIENT:  -2.0973E+00  6.3091E+00  1.2449E-01  1.5750E+00 -1.0971E+00 -6.9495E-01 -1.5870E-01  0.0000E+00 -2.9581E-01  0.0000E+00
            -8.7116E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1204.42180907993        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1725
 NPARAMETR:  1.0202E+00  2.2028E+00  5.2669E-01  3.4460E-01  1.6726E+00  8.8994E-01  3.1934E-01  1.0000E-02  1.3094E+00  1.0000E-02
             4.7659E+00
 PARAMETER:  1.2000E-01  8.8971E-01 -5.4115E-01 -9.6539E-01  6.1438E-01 -1.6596E-02 -1.0415E+00 -9.4033E+00  3.6954E-01 -2.1166E+01
             1.6615E+00
 GRADIENT:  -8.6786E-01  3.2921E+00 -1.8256E-01  1.1729E+00 -3.5615E-03 -5.0134E-02  9.4595E-03  0.0000E+00  3.9326E-02  0.0000E+00
            -9.8711E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1204.43043470244        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1903
 NPARAMETR:  1.0201E+00  2.2579E+00  5.2224E-01  3.0758E-01  1.7123E+00  8.8992E-01  3.1686E-01  1.0000E-02  1.4067E+00  1.0000E-02
             4.7681E+00
 PARAMETER:  1.1994E-01  9.1444E-01 -5.4963E-01 -1.0790E+00  6.3785E-01 -1.6623E-02 -1.0493E+00 -9.6104E+00  4.4125E-01 -2.2726E+01
             1.6619E+00
 GRADIENT:  -7.2222E-01  4.0812E+00  1.4469E-01  9.6209E-01 -5.6300E-01  4.1480E-02 -2.1397E-02  0.0000E+00  2.0908E-01  0.0000E+00
            -8.2241E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1204.43917923708        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2079
 NPARAMETR:  1.0201E+00  2.3296E+00  4.8280E-01  2.5882E-01  1.7533E+00  8.8983E-01  3.2092E-01  1.0000E-02  1.4694E+00  1.0000E-02
             4.7710E+00
 PARAMETER:  1.1989E-01  9.4569E-01 -6.2816E-01 -1.2516E+00  6.6152E-01 -1.6723E-02 -1.0366E+00 -1.0131E+01  4.8485E-01 -2.5340E+01
             1.6625E+00
 GRADIENT:  -3.7791E-01  2.6335E+00  1.7559E-01  4.7979E-01 -5.1790E-01  1.2734E-02 -3.4632E-02  0.0000E+00  8.7077E-02  0.0000E+00
            -1.3709E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1204.44032459001        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2254
 NPARAMETR:  1.0201E+00  2.3700E+00  4.4929E-01  2.3249E-01  1.7721E+00  8.8985E-01  3.2289E-01  1.0000E-02  1.5129E+00  1.0000E-02
             4.7724E+00
 PARAMETER:  1.1993E-01  9.6288E-01 -7.0010E-01 -1.3589E+00  6.7214E-01 -1.6702E-02 -1.0304E+00 -1.0547E+01  5.1403E-01 -2.7060E+01
             1.6628E+00
 GRADIENT:  -6.8021E-01  4.3341E+00  1.6653E-01  6.3124E-01 -7.0059E-01 -5.3934E-02 -6.1663E-02  0.0000E+00  2.1075E-02  0.0000E+00
            -3.8879E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1204.44332363345        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2429
 NPARAMETR:  1.0202E+00  2.4190E+00  3.9567E-01  1.9978E-01  1.7914E+00  8.8988E-01  3.2539E-01  1.0000E-02  1.5668E+00  1.0000E-02
             4.7740E+00
 PARAMETER:  1.1999E-01  9.8336E-01 -8.2717E-01 -1.5105E+00  6.8302E-01 -1.6667E-02 -1.0227E+00 -1.1231E+01  5.4904E-01 -2.9597E+01
             1.6632E+00
 GRADIENT:  -7.8278E-01  4.9328E+00  1.1491E-01  5.9664E-01 -7.2646E-01 -1.0940E-01 -7.2975E-02  0.0000E+00 -6.6741E-02  0.0000E+00
            -5.6779E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1204.44782056553        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2604
 NPARAMETR:  1.0203E+00  2.4691E+00  3.2352E-01  1.6562E-01  1.8069E+00  8.8997E-01  3.2595E-01  1.0000E-02  1.6585E+00  1.0000E-02
             4.7754E+00
 PARAMETER:  1.2010E-01  1.0038E+00 -1.0285E+00 -1.6981E+00  6.9160E-01 -1.6570E-02 -1.0210E+00 -1.2257E+01  6.0590E-01 -3.2924E+01
             1.6635E+00
 GRADIENT:  -6.0272E-01  4.0351E+00  4.4880E-02  4.2249E-01 -5.6372E-01 -1.3827E-01 -1.0693E-01  0.0000E+00 -1.2996E-01  0.0000E+00
            -6.7126E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1204.45660433272        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2788
 NPARAMETR:  1.0205E+00  2.4968E+00  2.7608E-01  1.4473E-01  1.8162E+00  8.9018E-01  3.2698E-01  1.0000E-02  1.9574E+00  1.0000E-02
             4.7772E+00
 PARAMETER:  1.2026E-01  1.0150E+00 -1.1871E+00 -1.8329E+00  6.9675E-01 -1.6336E-02 -1.0179E+00 -1.3051E+01  7.7164E-01 -3.5318E+01
             1.6639E+00
 GRADIENT:   3.4731E-01 -2.7161E+00  3.1734E-02 -7.0735E-02  3.4408E-01  2.4214E-01  3.1692E-01  0.0000E+00  7.7662E-02  0.0000E+00
             1.1424E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1204.46060412070        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     2958
 NPARAMETR:  1.0209E+00  2.4996E+00  2.5607E-01  1.4241E-01  1.8095E+00  8.9012E-01  3.2492E-01  1.0000E-02  1.9854E+00  1.0000E-02
             4.7771E+00
 PARAMETER:  1.1974E-01  1.0161E+00 -1.2520E+00 -1.8499E+00  6.9211E-01 -1.6581E-02 -1.0340E+00 -1.3051E+01  7.7802E-01 -3.5318E+01
             1.6631E+00
 GRADIENT:  -8.0907E-01 -8.0367E-02  4.8790E-03 -7.3916E-03 -4.5083E-02 -1.9012E-02 -6.2815E-02  0.0000E+00 -6.0487E-03  0.0000E+00
            -2.5160E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2958
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0594E-03 -1.2198E-02  5.8615E-06 -2.5890E-03 -1.0179E-05
 SE:             2.8254E-02  1.3682E-02  3.2258E-06  6.3935E-03  5.5592E-05
 N:                     100         100         100         100         100

 P VAL.:         8.5789E-01  3.7264E-01  6.9205E-02  6.8552E-01  8.5472E-01

 ETASHRINKSD(%)  5.3443E+00  5.4164E+01  9.9989E+01  7.8581E+01  9.9814E+01
 ETASHRINKVR(%)  1.0403E+01  7.8991E+01  1.0000E+02  9.5412E+01  1.0000E+02
 EBVSHRINKSD(%)  5.5227E+00  5.4526E+01  9.9982E+01  7.8683E+01  9.9747E+01
 EBVSHRINKVR(%)  1.0740E+01  7.9321E+01  1.0000E+02  9.5456E+01  9.9999E+01
 RELATIVEINF(%)  8.0552E+01  9.2474E-02  1.8095E-07  2.5268E-02  2.5941E-05
 EPSSHRINKSD(%)  1.3133E+01
 EPSSHRINKVR(%)  2.4541E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1204.4606041207021     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -469.30977755696392     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.83
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1204.461       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.50E+00  2.59E-01  1.42E-01  1.81E+00  8.90E-01  3.22E-01  1.00E-02  1.97E+00  1.00E-02  4.77E+00
 


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
+        1.23E+03
 
 TH 2
+       -1.13E+01  3.22E+02
 
 TH 3
+        6.22E+00  1.09E+01  4.62E-01
 
 TH 4
+       -3.58E+01  4.31E+02  1.44E+01  5.78E+02
 
 TH 5
+       -8.94E+00 -4.26E+01 -1.56E+00 -5.68E+01  5.79E+00
 
 TH 6
+       -1.76E+01 -2.37E+00  3.42E+00 -5.05E+00 -3.75E+00  2.23E+02
 
 TH 7
+       -1.97E+01 -3.29E+01 -7.87E-01 -4.39E+01  4.01E+00  2.75E+01  6.98E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.47E-01 -2.60E-01  8.84E-04 -3.45E-01  2.45E-02  8.27E-01  1.37E-01  0.00E+00  3.55E-03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.35E+01 -1.55E+01 -5.13E-01 -2.04E+01  2.09E+00  9.12E+00  3.04E+00  0.00E+00  5.71E-02  0.00E+00  1.55E+00
 
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
+        1.20E+03
 
 TH 2
+       -1.48E+02  3.41E+02
 
 TH 3
+        1.59E+00  1.80E+01  1.04E+01
 
 TH 4
+       -2.17E+02  4.38E+02  6.68E+00  6.11E+02
 
 TH 5
+        9.30E+00 -5.19E+01 -8.99E+00 -4.72E+01  2.20E+01
 
 TH 6
+        8.29E+00 -2.10E+01  1.14E+00 -2.87E+01  6.30E-01  1.97E+02
 
 TH 7
+       -1.73E+00 -3.55E+01 -1.79E+00 -2.94E+01  9.03E+00  1.65E+01  9.08E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.35E-01 -7.31E-01  4.34E-01  3.28E-01  2.03E-01  5.01E-01  2.00E+00  0.00E+00  2.96E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.74E+01 -1.48E+01 -4.77E-01 -1.55E+01  1.94E+00  3.11E+00  1.54E+01  0.00E+00  5.12E-01  0.00E+00  2.36E+01
 
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
+        1.20E+03
 
 TH 2
+       -2.86E+02  3.59E+02
 
 TH 3
+       -1.70E+00  9.21E+00  1.95E+00
 
 TH 4
+       -3.92E+02  4.74E+02  7.88E+00  6.45E+02
 
 TH 5
+        1.86E+01 -3.78E+01 -4.41E+00 -4.00E+01  1.18E+01
 
 TH 6
+        4.90E+00 -2.37E+01  5.07E+00 -4.36E+01 -9.51E+00  1.62E+02
 
 TH 7
+        4.15E+01 -7.94E+01 -3.96E-01 -8.85E+01  1.03E+01  1.50E+01  8.76E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.53E+00 -2.33E+00  1.59E-02 -2.58E+00  2.25E-01  7.29E-01  2.77E+00  0.00E+00  9.27E-02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.28E+02 -4.21E+01  6.74E-01 -5.88E+01  9.65E-01  2.04E+01  1.14E+01  0.00E+00  4.13E-01  0.00E+00  7.88E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       47.920
Stop Time:
Wed Sep 29 13:31:12 CDT 2021
