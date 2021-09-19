Sat Sep 18 14:54:20 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1671.88498545771        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.8129E+01 -1.0914E+02 -6.4230E+01 -8.2039E+01  3.5725E+01  3.5461E+00 -1.8246E+01  1.8186E+01  6.8450E+00  5.6811E+00
             2.0044E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1691.67836184981        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6372E-01  1.3451E+00  1.9672E+00  8.4554E-01  1.4017E+00  9.7153E-01  1.2271E+00  7.1513E-01  8.6159E-01  1.1466E+00
             9.5239E-01
 PARAMETER:  6.3045E-02  3.9645E-01  7.7662E-01 -6.7781E-02  4.3766E-01  7.1120E-02  3.0464E-01 -2.3528E-01 -4.8979E-02  2.3683E-01
             5.1217E-02
 GRADIENT:   2.2649E+01  5.9100E+00  1.2242E+01 -3.8867E+01 -2.1153E+00 -6.7290E+00  2.0693E+01 -1.0265E+00 -3.8732E+00 -3.5038E+01
            -2.1840E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1696.54969454800        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8076E-01  1.2671E+00  2.2670E+00  9.0591E-01  1.4593E+00  1.0219E+00  1.0076E+00  4.1915E-01  1.1180E+00  1.2494E+00
             9.7313E-01
 PARAMETER:  8.0572E-02  3.3669E-01  9.1844E-01  1.1869E-03  4.7795E-01  1.2163E-01  1.0753E-01 -7.6953E-01  2.1157E-01  3.2270E-01
             7.2765E-02
 GRADIENT:   5.9790E+01 -8.9790E+00  2.6473E+00 -2.2457E+01  2.5895E+01  1.3468E+01  7.7639E+00 -1.9477E-01  7.1190E+00 -1.9143E+01
            -4.8985E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.62386958366        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  9.7052E-01  1.1188E+00  1.8675E+00  1.0217E+00  1.3203E+00  9.9720E-01  1.1209E+00  1.6616E-01  9.7881E-01  1.2730E+00
             9.7250E-01
 PARAMETER:  7.0076E-02  2.1225E-01  7.2462E-01  1.2143E-01  3.7789E-01  9.7197E-02  2.1413E-01 -1.6948E+00  7.8580E-02  3.4136E-01
             7.2113E-02
 GRADIENT:  -4.0574E+00 -9.2477E-01  5.1707E-01 -3.4533E+00 -3.9460E-01  3.0141E-01  2.4565E+00 -1.6527E-03  6.2031E-01 -1.1674E+00
             1.8994E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1698.71492884868        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  9.7240E-01  1.1340E+00  1.9468E+00  1.0146E+00  1.3408E+00  9.9658E-01  1.0212E+00  1.7078E-01  1.0208E+00  1.2947E+00
             9.7115E-01
 PARAMETER:  7.2009E-02  2.2571E-01  7.6619E-01  1.1451E-01  3.9329E-01  9.6577E-02  1.2096E-01 -1.6674E+00  1.2062E-01  3.5828E-01
             7.0722E-02
 GRADIENT:  -4.8200E-02 -2.2622E-01  2.8091E-02 -1.8003E-01  2.1146E-02 -1.3097E-02  5.7406E-03 -6.3814E-03  4.0127E-02 -6.8701E-02
             4.1678E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.72806406738        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      631
 NPARAMETR:  9.7323E-01  1.2588E+00  1.8282E+00  9.3300E-01  1.3628E+00  9.9733E-01  9.7069E-01  2.1027E-01  1.0723E+00  1.3001E+00
             9.6936E-01
 PARAMETER:  7.2866E-02  3.3017E-01  7.0331E-01  3.0648E-02  4.0955E-01  9.7323E-02  7.0252E-02 -1.4594E+00  1.6982E-01  3.6247E-01
             6.8879E-02
 GRADIENT:   4.4051E-01  1.7322E+00  2.0175E-03  1.8945E+00 -6.1633E-01  2.7324E-02 -2.0477E-01 -3.3108E-03 -4.0137E-01  1.2148E-01
            -4.4376E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.76820334830        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  9.7390E-01  1.4853E+00  1.6512E+00  7.8397E-01  1.4177E+00  9.9829E-01  8.8788E-01  2.8660E-01  1.2024E+00  1.3181E+00
             9.7032E-01
 PARAMETER:  7.3556E-02  4.9562E-01  6.0152E-01 -1.4338E-01  4.4903E-01  9.8292E-02 -1.8919E-02 -1.1497E+00  2.8434E-01  3.7618E-01
             6.9868E-02
 GRADIENT:  -3.5988E-01  4.4931E+00 -2.8693E-01  4.0416E+00  4.7925E-01 -1.1288E-02 -2.8680E-01  2.2987E-02 -3.0182E-01 -9.7130E-02
             2.0202E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1698.81202665398        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  9.7523E-01  1.6806E+00  1.4969E+00  6.5350E-01  1.4645E+00  9.9832E-01  8.3331E-01  3.2862E-01  1.3528E+00  1.3398E+00
             9.7032E-01
 PARAMETER:  7.4916E-02  6.1914E-01  5.0341E-01 -3.2542E-01  4.8154E-01  9.8317E-02 -8.2353E-02 -1.0129E+00  4.0217E-01  3.9254E-01
             6.9869E-02
 GRADIENT:   9.8739E-01  6.6187E+00  2.3340E-01  3.6351E+00 -1.8058E+00 -2.8520E-01 -2.7865E-01  6.2993E-02 -2.3073E-01  2.0254E-01
            -3.5605E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1698.87407213022        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1169
 NPARAMETR:  9.7471E-01  1.8296E+00  1.3203E+00  5.5016E-01  1.5042E+00  1.0004E+00  8.0527E-01  2.4908E-01  1.4934E+00  1.3554E+00
             9.7057E-01
 PARAMETER:  7.4383E-02  7.0412E-01  3.7789E-01 -4.9754E-01  5.0829E-01  1.0044E-01 -1.1658E-01 -1.2900E+00  5.0107E-01  4.0412E-01
             7.0128E-02
 GRADIENT:  -1.3240E+00  1.3750E+00 -8.5071E-01  1.7180E+00  2.0539E+00  3.4689E-01 -4.7082E-01  6.1956E-02 -9.5549E-04  3.8890E-02
             5.2423E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1698.92531084138        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1345
 NPARAMETR:  9.7595E-01  2.0215E+00  1.1706E+00  4.2098E-01  1.5576E+00  9.9996E-01  7.7116E-01  1.2088E-01  1.7490E+00  1.3854E+00
             9.7159E-01
 PARAMETER:  7.5660E-02  8.0383E-01  2.5751E-01 -7.6517E-01  5.4316E-01  9.9958E-02 -1.5986E-01 -2.0129E+00  6.5902E-01  4.2602E-01
             7.1176E-02
 GRADIENT:   1.7297E-01  2.8576E+00 -1.3771E-01  8.8281E-01  3.9740E-01 -5.4979E-02 -4.2889E-01  1.6526E-02 -2.7582E-01  9.7301E-02
             3.7137E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1698.93485664819        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1526
 NPARAMETR:  9.7588E-01  2.0274E+00  1.1706E+00  4.1511E-01  1.5601E+00  1.0000E+00  7.7036E-01  3.9451E-02  1.7710E+00  1.3855E+00
             9.7165E-01
 PARAMETER:  7.5586E-02  8.0678E-01  2.5753E-01 -7.7921E-01  5.4473E-01  1.0003E-01 -1.6089E-01 -3.1327E+00  6.7156E-01  4.2603E-01
             7.1236E-02
 GRADIENT:   1.5654E-03 -3.3221E-01 -5.1730E-02 -1.1520E-02  3.9041E-01 -2.9651E-02  7.7085E-02  1.7883E-03  1.0295E-01 -1.1387E-01
            -1.5379E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1698.93585001867        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1702
 NPARAMETR:  9.7588E-01  2.0271E+00  1.1703E+00  4.1536E-01  1.5591E+00  1.0001E+00  7.7032E-01  1.2188E-02  1.7681E+00  1.3858E+00
             9.7167E-01
 PARAMETER:  7.5583E-02  8.0659E-01  2.5730E-01 -7.7862E-01  5.4413E-01  1.0006E-01 -1.6095E-01 -4.3073E+00  6.6993E-01  4.2625E-01
             7.1257E-02
 GRADIENT:  -2.3724E-03 -1.8382E-01 -4.0031E-03 -9.1841E-02 -4.9366E-02 -1.1482E-02 -5.8885E-03  1.7002E-04 -2.9163E-02  8.1446E-03
            -1.0082E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1698.93593257796        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1877
 NPARAMETR:  9.7588E-01  2.0268E+00  1.1720E+00  4.1565E-01  1.5593E+00  1.0001E+00  7.7017E-01  1.0000E-02  1.7687E+00  1.3859E+00
             9.7166E-01
 PARAMETER:  7.5584E-02  8.0647E-01  2.5875E-01 -7.7791E-01  5.4426E-01  1.0008E-01 -1.6114E-01 -4.6193E+00  6.7023E-01  4.2634E-01
             7.1251E-02
 GRADIENT:   1.5896E-03  2.0565E-02 -5.5899E-05 -6.4712E-03  1.0954E-02 -6.5710E-03 -6.2459E-03  0.0000E+00  1.3233E-03  4.3276E-03
            -1.5427E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1698.93708666453        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2059
 NPARAMETR:  9.7584E-01  2.0124E+00  1.1832E+00  4.2532E-01  1.5551E+00  1.0009E+00  7.7258E-01  1.0000E-02  1.7457E+00  1.3831E+00
             9.7173E-01
 PARAMETER:  7.5546E-02  7.9934E-01  2.6826E-01 -7.5492E-01  5.4151E-01  1.0090E-01 -1.5802E-01 -1.2226E+01  6.5713E-01  4.2432E-01
             7.1323E-02
 GRADIENT:   1.4687E-02 -6.4876E-02  1.1168E-02 -1.5881E-02  1.8026E-03  3.3137E-01  6.6326E-02  0.0000E+00  1.7639E-02 -3.7268E-02
             1.0124E-01

0ITERATION NO.:   68    OBJECTIVE VALUE:  -1698.93725403737        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:     2161
 NPARAMETR:  9.7584E-01  2.0124E+00  1.1829E+00  4.2533E-01  1.5551E+00  1.0000E+00  7.7229E-01  1.0000E-02  1.7455E+00  1.3832E+00
             9.7153E-01
 PARAMETER:  7.5546E-02  7.9932E-01  2.6800E-01 -7.5489E-01  5.4151E-01  1.0001E-01 -1.5839E-01 -1.2226E+01  6.5705E-01  4.2441E-01
             7.1115E-02
 GRADIENT:   1.6787E-02 -8.0119E-02 -1.4552E-03 -1.4482E-02  2.2565E-02 -1.6443E-02 -1.0162E-02  0.0000E+00 -1.1093E-03 -2.2792E-02
             4.9300E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2161
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.6545E-04 -2.4066E-02 -6.5018E-05  1.6732E-02 -4.0842E-02
 SE:             2.9766E-02  2.2819E-02  3.3565E-05  1.9123E-02  2.3756E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8216E-01  2.9158E-01  5.2735E-02  3.8158E-01  8.5577E-02

 ETASHRINKSD(%)  2.8182E-01  2.3554E+01  9.9888E+01  3.5936E+01  2.0414E+01
 ETASHRINKVR(%)  5.6284E-01  4.1561E+01  1.0000E+02  5.8958E+01  3.6660E+01
 EBVSHRINKSD(%)  4.2270E-01  2.0394E+01  9.9889E+01  4.0566E+01  1.6863E+01
 EBVSHRINKVR(%)  8.4362E-01  3.6629E+01  1.0000E+02  6.4676E+01  3.0882E+01
 RELATIVEINF(%)  9.9058E+01  2.7521E+00  4.9207E-05  1.5300E+00  3.2548E+01
 EPSSHRINKSD(%)  4.1061E+01
 EPSSHRINKVR(%)  6.5262E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1698.9372540373688     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -963.78642747363062     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1698.937       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  2.01E+00  1.18E+00  4.25E-01  1.56E+00  1.00E+00  7.72E-01  1.00E-02  1.75E+00  1.38E+00  9.72E-01
 


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
+        1.15E+03
 
 TH 2
+       -7.24E+00  3.01E+02
 
 TH 3
+       -2.69E-01  8.59E+00  1.13E+01
 
 TH 4
+       -6.44E+00  4.06E+02 -2.44E+01  7.31E+02
 
 TH 5
+       -2.80E+00 -3.53E+01 -1.95E+01  4.22E+01  1.67E+02
 
 TH 6
+        7.72E+00 -2.23E+00  1.30E-02 -1.91E+00 -1.80E+00  1.90E+02
 
 TH 7
+        2.15E+00 -6.15E+00  1.10E+01 -1.28E+01 -9.60E+00  4.26E+00  1.47E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.20E-01 -6.15E+00 -5.09E+00  3.26E+01  1.91E+00  7.31E-03  1.91E+01  0.00E+00  1.65E+01
 
 TH10
+        7.04E-01 -5.19E+00 -3.15E+00  1.49E+00 -3.03E+01 -1.06E+00 -3.98E+00  0.00E+00  1.67E+00  5.55E+01
 
 TH11
+       -7.79E+00 -1.42E+01 -9.91E+00  8.70E-01  5.58E+00 -6.83E-01  1.60E+01  0.00E+00  1.98E+00  1.56E+01  2.57E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.470
Stop Time:
Sat Sep 18 14:54:57 CDT 2021
