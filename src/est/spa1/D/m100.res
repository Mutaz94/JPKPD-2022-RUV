Thu Sep 30 03:53:27 CDT 2021
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
$DATA ../../../../data/spa1/D/dat100.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   38251.3947174877        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3246E+03  8.5795E+02 -3.2565E+01  9.3535E+02  8.9301E+01 -3.0289E+03 -1.4790E+03 -4.7610E+01 -2.1076E+03 -4.9090E+02
            -7.2502E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -306.171137296739        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.9854E-01  9.6804E-01  9.4981E-01  1.2035E+00  1.3624E+00  1.8011E+00  1.0701E+00  9.5834E-01  1.0340E+00  9.1791E-01
             1.4693E+01
 PARAMETER:  9.8542E-02  6.7513E-02  4.8503E-02  2.8520E-01  4.0922E-01  6.8842E-01  1.6771E-01  5.7447E-02  1.3339E-01  1.4349E-02
             2.7874E+00
 GRADIENT:  -5.5259E+01  1.2438E+01 -2.6418E+00  7.5582E+00 -6.4793E+00  3.6383E+01 -3.0533E+00  3.1041E+00 -6.2356E+00  1.9578E+00
            -5.1024E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -321.212031038709        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0569E+00  7.7812E-01  1.1460E+00  1.3839E+00  4.1477E+00  1.4428E+00  3.5648E-01  5.2947E-01  1.4426E+00  3.2497E-01
             1.5496E+01
 PARAMETER:  1.5538E-01 -1.5087E-01  2.3631E-01  4.2494E-01  1.5225E+00  4.6657E-01 -9.3149E-01 -5.3588E-01  4.6648E-01 -1.0240E+00
             2.8406E+00
 GRADIENT:  -4.9295E+01  1.2699E+01  6.0939E+00  1.8890E+01 -6.7770E+00 -3.6720E+01  1.8472E-01  1.5944E-01  1.2595E+01  2.0905E-02
             3.4420E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -365.151866563279        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0795E+00  2.9293E-01  3.1023E-01  1.4997E+00  1.4308E+01  1.9486E+00  4.2269E-02  1.0000E-02  1.3439E+00  1.2031E-01
             1.6527E+01
 PARAMETER:  1.7650E-01 -1.1278E+00 -1.0704E+00  5.0527E-01  2.7608E+00  7.6710E-01 -3.0637E+00 -7.0432E+00  3.9560E-01 -2.0177E+00
             2.9050E+00
 GRADIENT:   1.3730E+01  3.4849E+01 -2.1406E+01  5.9127E+01 -2.7601E+00 -1.9793E-01  1.4867E-02  0.0000E+00  8.6078E+00  1.2127E-03
             2.5993E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -413.372174715409        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  5.2490E-01  3.6385E-02  4.1728E-02  5.0209E-01  3.8410E+03  1.6804E+00  1.0000E-02  1.0000E-02  2.4406E-01  9.4594E-02
             1.7466E+01
 PARAMETER: -5.4454E-01 -3.2136E+00 -3.0766E+00 -5.8897E-01  8.3535E+00  6.1904E-01 -7.0048E+00 -1.4975E+01 -1.3104E+00 -2.2582E+00
             2.9603E+00
 GRADIENT:  -4.4026E+01  1.1669E+01 -7.8907E+00  8.6915E+01 -1.9963E-03 -6.6838E-01  0.0000E+00  0.0000E+00  3.3399E+00  2.9905E-09
             9.7520E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -422.685245422628        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      473
 NPARAMETR:  5.6174E-01  2.5598E-02  3.5384E-02  4.1917E-01  4.7633E+03  1.6826E+00  1.0000E-02  1.0000E-02  1.5876E-01  1.6299E-01
             1.6234E+01
 PARAMETER: -4.7672E-01 -3.5652E+00 -3.2415E+00 -7.6949E-01  8.5687E+00  6.2033E-01 -5.4660E+00 -1.6267E+01 -1.7404E+00 -1.7140E+00
             2.8871E+00
 GRADIENT:   5.1253E+01 -1.6036E-01 -4.3682E+01  4.0948E+01  6.8776E-04  1.1182E+00  0.0000E+00  0.0000E+00  6.5994E-01  6.6327E-10
            -1.0853E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -428.239325922826        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      649
 NPARAMETR:  4.2662E-01  1.0000E-02  2.1951E-02  2.8428E-01  1.5660E+04  1.6078E+00  1.0000E-02  1.0000E-02  5.1595E-02  1.5223E-01
             1.5739E+01
 PARAMETER: -7.5185E-01 -4.9482E+00 -3.7190E+00 -1.1578E+00  9.7589E+00  5.7489E-01 -5.8616E+00 -1.9134E+01 -2.8643E+00 -1.7823E+00
             2.8561E+00
 GRADIENT:  -1.4498E+00  0.0000E+00 -2.6213E+00  3.5283E+00  2.5790E-05  2.3308E+00  0.0000E+00  0.0000E+00  8.5630E-02 -3.1893E-12
             3.1840E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -428.297616362565        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      842
 NPARAMETR:  4.2724E-01  1.0000E-02  2.1933E-02  2.8356E-01  5.6579E+03  1.5959E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6246E-01
             1.5737E+01
 PARAMETER: -7.5040E-01 -4.9482E+00 -3.7198E+00 -1.1603E+00  8.7408E+00  5.6745E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7173E+00
             2.8560E+00
 GRADIENT:   1.9431E-01  0.0000E+00 -7.9129E-01  5.2320E-01  7.2589E-05  3.4328E-01  0.0000E+00  0.0000E+00  0.0000E+00 -2.9790E-11
            -1.7443E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -428.304155218922        NO. OF FUNC. EVALS.: 213
 CUMULATIVE NO. OF FUNC. EVALS.:     1055             RESET HESSIAN, TYPE I
 NPARAMETR:  4.2712E-01  1.0000E-02  2.1850E-02  2.8303E-01  3.3839E+03  1.5940E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5580E-01
             1.5736E+01
 PARAMETER: -7.5069E-01 -4.9482E+00 -3.7235E+00 -1.1622E+00  8.2268E+00  5.6625E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7592E+00
             2.8559E+00
 GRADIENT:   3.7564E+01  0.0000E+00  5.3107E+01  1.8251E+01  1.3451E-04  5.3374E+00  0.0000E+00  0.0000E+00  0.0000E+00 -3.8774E-11
             2.7785E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -428.306375436092        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  4.2671E-01  1.0000E-02  2.1850E-02  2.8267E-01  2.4450E+03  1.5940E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5462E-01
             1.5741E+01
 PARAMETER: -7.5165E-01 -4.9482E+00 -3.7236E+00 -1.1635E+00  7.9018E+00  5.6622E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7668E+00
             2.8563E+00
 GRADIENT:   2.9286E-01  0.0000E+00 -5.8815E-01  1.7347E-01  1.6487E-04  7.6258E-02  0.0000E+00  0.0000E+00  0.0000E+00 -1.4156E-10
             4.5082E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -428.310540945942        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     1445
 NPARAMETR:  4.2664E-01  1.0000E-02  2.1790E-02  2.8242E-01  1.6574E+03  1.5937E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6016E-01
             1.5734E+01
 PARAMETER: -7.5182E-01 -4.9482E+00 -3.7263E+00 -1.1644E+00  7.5130E+00  5.6607E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7316E+00
             2.8558E+00
 GRADIENT:   1.1187E+00  0.0000E+00 -2.5804E+00  2.2998E+00  2.8238E-04 -1.7123E-02  0.0000E+00  0.0000E+00  0.0000E+00 -3.3155E-10
            -7.8925E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -428.314686335781        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     1646             RESET HESSIAN, TYPE I
 NPARAMETR:  4.2633E-01  1.0000E-02  2.1750E-02  2.8202E-01  1.0662E+03  1.5935E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5630E-01
             1.5734E+01
 PARAMETER: -7.5253E-01 -4.9482E+00 -3.7281E+00 -1.1658E+00  7.0719E+00  5.6595E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7560E+00
             2.8558E+00
 GRADIENT:   3.7644E+01  0.0000E+00  5.3312E+01  1.8207E+01  4.2200E-04  5.3389E+00  0.0000E+00  0.0000E+00  0.0000E+00 -4.4673E-10
             2.7779E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -428.316794399947        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1834
 NPARAMETR:  4.2595E-01  1.0000E-02  2.1750E-02  2.8167E-01  8.0637E+02  1.5935E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5358E-01
             1.5740E+01
 PARAMETER: -7.5343E-01 -4.9482E+00 -3.7281E+00 -1.1670E+00  6.7925E+00  5.6591E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7735E+00
             2.8562E+00
 GRADIENT:   2.8035E-01  0.0000E+00 -6.1853E-01  1.4714E-01  4.9162E-04  7.5248E-02  0.0000E+00  0.0000E+00  0.0000E+00 -1.2852E-09
            -9.7832E-04

0ITERATION NO.:   65    OBJECTIVE VALUE:  -428.320878813195        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     2033
 NPARAMETR:  4.2590E-01  1.0000E-02  2.1691E-02  2.8143E-01  5.7728E+02  1.5932E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5451E-01
             1.5733E+01
 PARAMETER: -7.5356E-01 -4.9482E+00 -3.7308E+00 -1.1679E+00  6.4583E+00  5.6577E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7675E+00
             2.8558E+00
 GRADIENT:   1.1198E+00  0.0000E+00 -2.6368E+00  2.2988E+00  7.9997E-04 -1.9429E-02  0.0000E+00  0.0000E+00  0.0000E+00 -2.5407E-09
            -7.9952E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -428.325018908116        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     2236             RESET HESSIAN, TYPE I
 NPARAMETR:  4.2561E-01  1.0000E-02  2.1653E-02  2.8104E-01  3.9510E+02  1.5931E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5734E-01
             1.5733E+01
 PARAMETER: -7.5423E-01 -4.9482E+00 -3.7326E+00 -1.1693E+00  6.0791E+00  5.6565E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7493E+00
             2.8557E+00
 GRADIENT:   3.7755E+01  0.0000E+00  5.3441E+01  1.8244E+01  1.1297E-03  5.3343E+00  0.0000E+00  0.0000E+00  0.0000E+00 -3.2624E-09
             2.7770E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -428.327170595557        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2423
 NPARAMETR:  4.2524E-01  1.0000E-02  2.1653E-02  2.8070E-01  3.1076E+02  1.5930E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6391E-01
             1.5738E+01
 PARAMETER: -7.5511E-01 -4.9482E+00 -3.7326E+00 -1.1705E+00  5.8390E+00  5.6562E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.7085E+00
             2.8561E+00
 GRADIENT:   2.6412E-01  0.0000E+00 -6.3566E-01  1.0556E-01  1.2478E-03  7.4731E-02  0.0000E+00  0.0000E+00  0.0000E+00 -9.7918E-09
             3.4561E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -428.331312837657        NO. OF FUNC. EVALS.: 202
 CUMULATIVE NO. OF FUNC. EVALS.:     2625
 NPARAMETR:  4.2520E-01  1.0000E-02  2.1596E-02  2.8047E-01  2.3315E+02  1.5928E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6680E-01
             1.5732E+01
 PARAMETER: -7.5519E-01 -4.9482E+00 -3.7353E+00 -1.1713E+00  5.5517E+00  5.6549E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6909E+00
             2.8557E+00
 GRADIENT:   1.1150E+00  0.0000E+00 -2.6847E+00  2.2892E+00  1.9434E-03 -2.1390E-02  0.0000E+00  0.0000E+00  0.0000E+00 -1.7991E-08
            -8.0615E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -428.335603306403        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     2828             RESET HESSIAN, TYPE I
 NPARAMETR:  4.2494E-01  1.0000E-02  2.1558E-02  2.8010E-01  1.6830E+02  1.5926E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6975E-01
             1.5731E+01
 PARAMETER: -7.5581E-01 -4.9482E+00 -3.7370E+00 -1.1726E+00  5.2258E+00  5.6538E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6734E+00
             2.8557E+00
 GRADIENT:   3.7852E+01  0.0000E+00  5.3584E+01  1.8271E+01  2.6137E-03  5.3303E+00  0.0000E+00  0.0000E+00  0.0000E+00 -2.1223E-08
             2.7764E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -428.337902041527        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3015
 NPARAMETR:  4.2458E-01  1.0000E-02  2.1560E-02  2.7977E-01  1.3691E+02  1.5926E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6860E-01
             1.5737E+01
 PARAMETER: -7.5666E-01 -4.9482E+00 -3.7369E+00 -1.1738E+00  5.0193E+00  5.6536E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6802E+00
             2.8560E+00
 GRADIENT:   2.4474E-01  0.0000E+00 -6.5548E-01  6.6985E-02  2.7421E-03  7.3700E-02  0.0000E+00  0.0000E+00  0.0000E+00 -5.2790E-08
             6.4529E-03

0ITERATION NO.:   95    OBJECTIVE VALUE:  -428.342254154177        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     3214
 NPARAMETR:  4.2457E-01  1.0000E-02  2.1503E-02  2.7956E-01  1.0695E+02  1.5924E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7123E-01
             1.5731E+01
 PARAMETER: -7.5668E-01 -4.9482E+00 -3.7396E+00 -1.1746E+00  4.7724E+00  5.6523E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6647E+00
             2.8556E+00
 GRADIENT:   1.1029E+00  0.0000E+00 -2.7257E+00  2.2713E+00  4.1088E-03 -2.3175E-02  0.0000E+00  0.0000E+00  0.0000E+00 -8.8990E-08
            -8.1012E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -428.346870070466        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     3413             RESET HESSIAN, TYPE I
 NPARAMETR:  4.2434E-01  1.0000E-02  2.1467E-02  2.7920E-01  8.0815E+01  1.5922E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7051E-01
             1.5730E+01
 PARAMETER: -7.5723E-01 -4.9482E+00 -3.7413E+00 -1.1758E+00  4.4922E+00  5.6514E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6690E+00
             2.8556E+00
 GRADIENT:   3.7930E+01  0.0000E+00  5.3747E+01  1.8289E+01  5.2904E-03  5.3269E+00  0.0000E+00  0.0000E+00  0.0000E+00 -9.0596E-08
             2.7760E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -428.349419770881        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3600
 NPARAMETR:  4.2400E-01  1.0000E-02  2.1469E-02  2.7888E-01  6.7675E+01  1.5922E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6966E-01
             1.5736E+01
 PARAMETER: -7.5803E-01 -4.9482E+00 -3.7412E+00 -1.1770E+00  4.3147E+00  5.6513E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6739E+00
             2.8559E+00
 GRADIENT:   2.2101E-01  0.0000E+00 -6.7970E-01  3.0929E-02  5.2504E-03  7.1995E-02  0.0000E+00  0.0000E+00  0.0000E+00 -2.1515E-07
             7.6837E-03

0ITERATION NO.:  110    OBJECTIVE VALUE:  -428.354114993663        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     3799
 NPARAMETR:  4.2402E-01  1.0000E-02  2.1413E-02  2.7868E-01  5.4730E+01  1.5920E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6910E-01
             1.5729E+01
 PARAMETER: -7.5797E-01 -4.9482E+00 -3.7438E+00 -1.1777E+00  4.1024E+00  5.6502E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6772E+00
             2.8555E+00
 GRADIENT:   1.0815E+00  0.0000E+00 -2.7597E+00  2.2421E+00  7.5882E-03 -2.4705E-02  0.0000E+00  0.0000E+00  0.0000E+00 -3.2476E-07
            -8.1089E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -428.359168264014        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     4000             RESET HESSIAN, TYPE I
 NPARAMETR:  4.2383E-01  1.0000E-02  2.1378E-02  2.7834E-01  4.3015E+01  1.5919E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7322E-01
             1.5729E+01
 PARAMETER: -7.5842E-01 -4.9482E+00 -3.7454E+00 -1.1789E+00  3.8615E+00  5.6494E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6532E+00
             2.8555E+00
 GRADIENT:   3.7985E+01  0.0000E+00  5.3939E+01  1.8293E+01  9.3470E-03  5.3242E+00  0.0000E+00  0.0000E+00  0.0000E+00 -3.1703E-07
             2.7760E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -428.362006270510        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4188
 NPARAMETR:  4.2352E-01  1.0000E-02  2.1380E-02  2.7803E-01  3.6929E+01  1.5919E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6642E-01
             1.5735E+01
 PARAMETER: -7.5916E-01 -4.9482E+00 -3.7453E+00 -1.1800E+00  3.7090E+00  5.6495E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.6932E+00
             2.8559E+00
 GRADIENT:   1.9491E-01  0.0000E+00 -7.1713E-01  2.5501E-03  8.6021E-03  6.9353E-02  0.0000E+00  0.0000E+00  0.0000E+00 -6.7334E-07
             4.4550E-03

0ITERATION NO.:  125    OBJECTIVE VALUE:  -428.367108638726        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4367
 NPARAMETR:  4.2301E-01  1.0000E-02  2.1298E-02  2.7776E-01  3.0690E+01  1.5910E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.7938E-01
             1.5713E+01
 PARAMETER: -7.6035E-01 -4.9482E+00 -3.7492E+00 -1.1810E+00  3.5239E+00  5.6434E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -1.1752E+00
             2.8545E+00
 GRADIENT:   3.7718E+01  0.0000E+00  5.3245E+01  1.9657E+01  1.2242E-02  5.1436E+00  0.0000E+00  0.0000E+00  0.0000E+00 -1.8922E-06
             2.7136E+01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -428.368560095719        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     4443
 NPARAMETR:  4.1956E-01  1.0000E-02  2.1095E-02  2.7599E-01  3.0651E+01  1.5885E+00  1.0000E-02  1.0000E-02  1.0000E-02  6.4457E-01
             1.5652E+01
 PARAMETER: -7.6856E-01 -4.9482E+00 -3.7587E+00 -1.1874E+00  3.5227E+00  5.6276E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00 -3.3917E-01
             2.8506E+00
 GRADIENT:   3.5850E+01  0.0000E+00  5.2715E+01  2.2100E+01  9.6475E-03  4.7895E+00  0.0000E+00  0.0000E+00  0.0000E+00 -1.3178E-05
             2.5242E+01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -428.370202727405        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:     4561
 NPARAMETR:  4.2079E-01  1.0000E-02  2.1192E-02  2.7686E-01  2.9967E+01  1.5892E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7296E+01
             1.5669E+01
 PARAMETER: -7.6562E-01 -4.9482E+00 -3.7541E+00 -1.1842E+00  3.5001E+00  5.6326E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00  2.9505E+00
             2.8517E+00
 GRADIENT:  -1.1124E+00  0.0000E+00 -4.2079E+00  5.0515E+00  1.1530E-02 -4.1840E-01  0.0000E+00  0.0000E+00  0.0000E+00  6.3341E-05
            -2.7485E+00

0ITERATION NO.:  140    OBJECTIVE VALUE:  -428.382725097440        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:     4670
 NPARAMETR:  4.2095E-01  1.0000E-02  2.1186E-02  2.7624E-01  2.3456E+01  1.5903E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5685E+01
             1.5690E+01
 PARAMETER: -7.6524E-01 -4.9482E+00 -3.7544E+00 -1.1865E+00  3.2551E+00  5.6395E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00  2.8527E+00
             2.8530E+00
 GRADIENT:   3.5834E+01  0.0000E+00  5.6106E+01  1.7247E+01  1.2405E-02  5.2845E+00  0.0000E+00  0.0000E+00  0.0000E+00  2.2897E-02
             2.7268E+01

0ITERATION NO.:  145    OBJECTIVE VALUE:  -428.390519782377        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:     4781
 NPARAMETR:  4.2240E-01  1.0000E-02  2.1657E-02  2.7224E-01  2.2010E+01  1.5886E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.4017E+01
             1.6123E+01
 PARAMETER: -7.6282E-01 -4.9482E+00 -3.7569E+00 -1.1900E+00  3.1599E+00  5.6476E-01 -5.8616E+00 -1.9134E+01 -4.6505E+00  2.7679E+00
             2.8562E+00
 GRADIENT:  -3.9819E-01  0.0000E+00 -4.7672E+01  1.6021E+02 -3.0627E+01  2.7801E-01  0.0000E+00  0.0000E+00  0.0000E+00  3.5452E+01
            -6.6830E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4781
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8112E-03  2.1896E-06  5.4994E-05 -1.7629E-04 -5.0221E-04
 SE:             2.8628E-02  1.0193E-06  2.6189E-04  3.0715E-04  6.8684E-04
 N:                     100         100         100         100         100

 P VAL.:         9.2177E-01  3.1708E-02  8.3368E-01  5.6601E-01  4.6466E-01

 ETASHRINKSD(%)  4.0931E+00  9.9997E+01  9.9123E+01  9.8971E+01  9.7699E+01
 ETASHRINKVR(%)  8.0187E+00  1.0000E+02  9.9992E+01  9.9989E+01  9.9947E+01
 EBVSHRINKSD(%)  4.2595E+00  9.9995E+01  9.9102E+01  9.8920E+01  9.7841E+01
 EBVSHRINKVR(%)  8.3376E+00  1.0000E+02  9.9992E+01  9.9988E+01  9.9953E+01
 RELATIVEINF(%)  8.9102E-01  3.4545E-08  2.0454E-05  2.2791E-05  7.1580E-04
 EPSSHRINKSD(%)  3.5414E+00
 EPSSHRINKVR(%)  6.9575E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -428.39051978237734     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       490.54801342229536     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    86.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -428.391       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.22E-01  1.00E-02  2.11E-02  2.75E-01  2.13E+01  1.59E+00  1.00E-02  1.00E-02  1.00E-02  1.44E+01  1.57E+01
 


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
+        2.22E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.72E+04  0.00E+00  2.05E+06
 
 TH 4
+        4.58E+04  0.00E+00 -1.13E+05  5.62E+04
 
 TH 5
+        3.34E-01  0.00E+00 -3.50E+00  1.18E+00  1.06E+00
 
 TH 6
+       -9.44E+00  0.00E+00  5.01E+02 -1.65E+04 -8.26E-02  5.69E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -4.16E-01  0.00E+00 -3.87E+01 -1.82E+00 -4.60E-02  1.41E-01  0.00E+00  0.00E+00  0.00E+00  3.12E+00
 
 TH11
+       -2.26E+01  0.00E+00  3.61E+02 -1.82E+01 -4.01E-03  1.40E-01  0.00E+00  0.00E+00  0.00E+00  5.01E-03  4.03E+00
 
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
 #CPUT: Total CPU Time in Seconds,       94.826
Stop Time:
Thu Sep 30 03:55:04 CDT 2021
