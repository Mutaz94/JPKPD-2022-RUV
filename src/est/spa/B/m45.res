Wed Sep 29 11:15:30 CDT 2021
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
$DATA ../../../../data/spa/B/dat45.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1155.71503148910        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1172E+02 -9.8301E+00  1.9540E+01  3.2148E+01  7.7494E+01 -4.8844E+00  7.5045E+00 -1.9303E+02 -1.7050E+01  1.9792E+01
            -7.2545E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1641.28493995746        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      129
 NPARAMETR:  1.1380E+00  1.0585E+00  1.0911E+00  9.8621E-01  1.0013E+00  1.2428E+00  9.7887E-01  1.0364E+00  9.8931E-01  8.7304E-01
             9.4110E-01
 PARAMETER:  2.2928E-01  1.5688E-01  1.8718E-01  8.6117E-02  1.0130E-01  3.1737E-01  7.8645E-02  1.3571E-01  8.9252E-02 -3.5779E-02
             3.9289E-02
 GRADIENT:   2.4756E+02  1.2120E+01  3.0525E+00  7.9627E+00 -6.6863E+00 -1.5746E+01  5.1207E-02  3.4790E+00 -9.1199E+00  8.9920E+00
            -1.7345E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1645.95674702441        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  1.1201E+00  1.0102E+00  8.6847E-01  9.9679E-01  8.9533E-01  1.2942E+00  1.0691E+00  8.5748E-01  1.0278E+00  6.5928E-01
             9.8040E-01
 PARAMETER:  2.1342E-01  1.1019E-01 -4.1028E-02  9.6781E-02 -1.0567E-02  3.5790E-01  1.6685E-01 -5.3763E-02  1.2746E-01 -3.1660E-01
             8.0208E-02
 GRADIENT:   2.0385E+02 -1.6363E+01 -2.2826E+01  1.3047E+01  3.6245E+01  8.1009E+00  2.1180E+00  3.8609E+00  9.3683E+00 -2.2289E+00
             6.3334E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1658.77659051823        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      486
 NPARAMETR:  1.0099E+00  1.1391E+00  7.6879E-01  9.1004E-01  8.9630E-01  1.3169E+00  1.0059E+00  3.8145E-01  1.0362E+00  7.3850E-01
             9.2511E-01
 PARAMETER:  1.0988E-01  2.3023E-01 -1.6293E-01  5.7312E-03 -9.4754E-03  3.7526E-01  1.0591E-01 -8.6378E-01  1.3557E-01 -2.0313E-01
             2.2152E-02
 GRADIENT:   5.8190E+01  6.3449E+00  2.8107E+00  2.8522E+00 -4.7539E+00  4.1780E+01  9.4087E-01 -3.7157E-01 -1.7162E+00 -2.5496E-01
            -2.5034E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.93880253974        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      662
 NPARAMETR:  9.6628E-01  1.2206E+00  7.2507E-01  8.5442E-01  9.1585E-01  1.1698E+00  9.5142E-01  1.5653E-01  1.0909E+00  7.4046E-01
             9.8238E-01
 PARAMETER:  6.5702E-02  2.9936E-01 -2.2149E-01 -5.7337E-02  1.2097E-02  2.5682E-01  5.0198E-02 -1.7545E+00  1.8704E-01 -2.0048E-01
             8.2218E-02
 GRADIENT:  -9.9374E-02  1.9352E+00  7.7380E-01  1.6161E+00 -2.9470E+00 -1.9160E-01  4.2172E-01 -7.0253E-02  2.3404E-01  4.4473E-01
             3.6060E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1662.95442709169        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      838
 NPARAMETR:  9.6668E-01  1.2944E+00  7.0246E-01  8.0894E-01  9.4176E-01  1.1686E+00  9.1051E-01  1.5654E-01  1.1345E+00  7.4575E-01
             9.8284E-01
 PARAMETER:  6.6111E-02  3.5802E-01 -2.5316E-01 -1.1203E-01  3.9990E-02  2.5581E-01  6.2545E-03 -1.7544E+00  2.2623E-01 -1.9337E-01
             8.2693E-02
 GRADIENT:   3.3679E-01  3.9819E+00  7.5311E-01  3.0908E+00 -1.8043E+00 -6.4899E-01 -1.5128E-02 -4.6707E-02 -4.9473E-01 -2.3601E-02
            -3.4050E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1662.97260364289        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1018
 NPARAMETR:  9.6742E-01  1.3069E+00  6.9685E-01  7.9724E-01  9.4762E-01  1.1763E+00  9.0283E-01  1.5939E-01  1.1476E+00  7.4718E-01
             9.8357E-01
 PARAMETER:  6.6879E-02  3.6763E-01 -2.6118E-01 -1.2660E-01  4.6200E-02  2.6237E-01 -2.2211E-03 -1.7364E+00  2.3769E-01 -1.9145E-01
             8.3436E-02
 GRADIENT:   1.5879E+00 -2.0484E+00 -4.1286E-01 -5.6959E-01  8.5091E-01  2.0774E+00  1.1551E-02 -3.4881E-02 -2.3512E-02  1.2825E-01
             9.4220E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1662.97591348445        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1205             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6764E-01  1.3074E+00  6.9808E-01  7.9755E-01  9.4680E-01  1.1771E+00  9.0277E-01  1.7653E-01  1.1477E+00  7.4587E-01
             9.8330E-01
 PARAMETER:  6.7108E-02  3.6801E-01 -2.5942E-01 -1.2622E-01  4.5331E-02  2.6302E-01 -2.2853E-03 -1.6343E+00  2.3779E-01 -1.9321E-01
             8.3163E-02
 GRADIENT:   4.6951E+02  2.6710E+02  5.2630E+00  6.6364E+01  6.3108E+00  2.3491E+02  5.1729E+00  2.3321E-02  1.4857E+01  8.6175E-01
             6.9053E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1662.97909986452        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1388             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6763E-01  1.3071E+00  6.9696E-01  7.9819E-01  9.4794E-01  1.1770E+00  9.0243E-01  1.9131E-01  1.1486E+00  7.4534E-01
             9.8341E-01
 PARAMETER:  6.7093E-02  3.6785E-01 -2.6103E-01 -1.2541E-01  4.6535E-02  2.6300E-01 -2.6667E-03 -1.5539E+00  2.3855E-01 -1.9391E-01
             8.3266E-02
 GRADIENT:   4.6929E+02  2.6520E+02  2.8907E+00  6.7907E+01  1.0978E+01  2.3482E+02  5.1240E+00  3.7850E-02  1.5236E+01  8.5547E-01
             8.7653E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1662.98166415355        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1562
 NPARAMETR:  9.6817E-01  1.3119E+00  6.9762E-01  7.9572E-01  9.4962E-01  1.1718E+00  9.0131E-01  2.2257E-01  1.1495E+00  7.4505E-01
             9.8350E-01
 PARAMETER:  6.7650E-02  3.7149E-01 -2.6009E-01 -1.2851E-01  4.8309E-02  2.5852E-01 -3.9084E-03 -1.4025E+00  2.3932E-01 -1.9430E-01
             8.3359E-02
 GRADIENT:   2.8150E+00  2.6035E-01 -4.1621E-01  9.7059E-01  3.4611E-01  4.8451E-01  1.5191E-01 -4.7859E-02 -1.1050E-01 -1.0268E-02
             5.7633E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1662.98407132509        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1737
 NPARAMETR:  9.6718E-01  1.3197E+00  6.9895E-01  7.9090E-01  9.5454E-01  1.1678E+00  8.9612E-01  3.0919E-01  1.1547E+00  7.4371E-01
             9.8326E-01
 PARAMETER:  6.6628E-02  3.7742E-01 -2.5817E-01 -1.3459E-01  5.3478E-02  2.5513E-01 -9.6757E-03 -1.0738E+00  2.4386E-01 -1.9611E-01
             8.3117E-02
 GRADIENT:   1.1868E+00 -6.4821E-01 -1.4141E+00  9.4946E-01  1.7207E+00 -8.9040E-01  1.5194E-01 -1.4581E-02 -1.5277E-01  1.6707E-02
             1.3150E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1662.98582505117        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1913
 NPARAMETR:  9.6588E-01  1.3234E+00  7.0323E-01  7.8849E-01  9.5797E-01  1.1678E+00  8.9186E-01  3.5771E-01  1.1594E+00  7.4376E-01
             9.8301E-01
 PARAMETER:  6.5287E-02  3.8021E-01 -2.5207E-01 -1.3763E-01  5.7064E-02  2.5510E-01 -1.4448E-02 -9.2803E-01  2.4786E-01 -1.9604E-01
             8.2863E-02
 GRADIENT:  -9.1891E-01 -5.7505E-01 -8.4065E-01  1.8254E-02  4.9408E-01 -8.6037E-01  1.8349E-02  1.0194E-02 -7.5933E-02 -6.8335E-02
             2.1333E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1663.00644606581        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2090
 NPARAMETR:  9.6667E-01  1.3164E+00  7.2651E-01  7.9501E-01  9.6687E-01  1.1690E+00  8.9046E-01  4.1393E-01  1.1606E+00  7.5658E-01
             9.8357E-01
 PARAMETER:  6.6105E-02  3.7494E-01 -2.1951E-01 -1.2940E-01  6.6304E-02  2.5612E-01 -1.6014E-02 -7.8207E-01  2.4891E-01 -1.7894E-01
             8.3431E-02
 GRADIENT:   6.3900E-01 -4.3593E-02 -6.0982E-01  5.0147E-01  1.1331E+00 -3.6958E-01  1.1402E-01  8.0285E-03  9.0430E-02  2.4618E-01
             1.9034E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1663.01802653274        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2275             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6745E-01  1.3146E+00  7.2788E-01  7.9541E-01  9.6595E-01  1.1767E+00  8.9009E-01  4.1422E-01  1.1600E+00  7.5393E-01
             9.8308E-01
 PARAMETER:  6.6911E-02  3.7357E-01 -2.1762E-01 -1.2890E-01  6.5356E-02  2.6269E-01 -1.6436E-02 -7.8137E-01  2.4843E-01 -1.8245E-01
             8.2940E-02
 GRADIENT:   4.7117E+02  2.7488E+02  3.8156E+00  6.6915E+01  7.9835E+00  2.3561E+02  4.9902E+00  1.4486E-01  1.5699E+01  6.3737E-01
             6.1430E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1663.01911741831        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2453
 NPARAMETR:  9.6749E-01  1.3140E+00  7.2757E-01  7.9597E-01  9.6594E-01  1.1762E+00  8.9067E-01  4.2405E-01  1.1600E+00  7.5442E-01
             9.8331E-01
 PARAMETER:  6.6947E-02  3.7305E-01 -2.1805E-01 -1.2820E-01  6.5344E-02  2.6230E-01 -1.5785E-02 -7.5791E-01  2.4842E-01 -1.8181E-01
             8.3167E-02
 GRADIENT:   2.0237E+00 -7.7468E-01 -4.4328E-01 -3.2389E-01  8.3545E-01  2.1950E+00  1.4064E-02  1.6998E-02  1.9659E-01  1.1128E-01
             9.5635E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1663.02033444443        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2634
 NPARAMETR:  9.6744E-01  1.3125E+00  7.2787E-01  7.9662E-01  9.6560E-01  1.1767E+00  8.9117E-01  4.2628E-01  1.1593E+00  7.5422E-01
             9.8331E-01
 PARAMETER:  6.6903E-02  3.7196E-01 -2.1763E-01 -1.2738E-01  6.4994E-02  2.6268E-01 -1.5219E-02 -7.5267E-01  2.4779E-01 -1.8207E-01
             8.3168E-02
 GRADIENT:   1.9572E+00 -1.3741E+00 -6.9391E-01 -5.0382E-01  1.2763E+00  2.3571E+00  1.0501E-02  2.7100E-02  2.3454E-01  1.4578E-01
             1.4311E-01

0ITERATION NO.:   77    OBJECTIVE VALUE:  -1663.02056376326        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     2693
 NPARAMETR:  9.6750E-01  1.3130E+00  7.2808E-01  7.9673E-01  9.6549E-01  1.1760E+00  8.9119E-01  4.2308E-01  1.1591E+00  7.5404E-01
             9.8325E-01
 PARAMETER:  6.6964E-02  3.7230E-01 -2.1734E-01 -1.2723E-01  6.4881E-02  2.6215E-01 -1.5196E-02 -7.6019E-01  2.4765E-01 -1.8230E-01
             8.3105E-02
 GRADIENT:   9.2191E-02  9.5224E-01 -2.4811E-01 -1.4092E-01  5.4078E-01 -2.1163E-01 -6.4202E-03  2.7908E-03  9.9070E-02  5.3384E-02
             4.0956E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2693
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1712E-04 -1.8133E-02 -1.3136E-02  1.0884E-02 -2.6435E-02
 SE:             2.9901E-02  2.1893E-02  6.1845E-03  2.4462E-02  2.1503E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9687E-01  4.0753E-01  3.3670E-02  6.5638E-01  2.1892E-01

 ETASHRINKSD(%)  1.0000E-10  2.6656E+01  7.9281E+01  1.8048E+01  2.7964E+01
 ETASHRINKVR(%)  1.0000E-10  4.6206E+01  9.5707E+01  3.2839E+01  4.8108E+01
 EBVSHRINKSD(%)  3.0993E-01  2.6401E+01  8.1181E+01  1.8303E+01  2.7602E+01
 EBVSHRINKVR(%)  6.1889E-01  4.5832E+01  9.6459E+01  3.3257E+01  4.7585E+01
 RELATIVEINF(%)  9.9271E+01  2.5953E+00  3.3224E-01  3.9759E+00  6.4511E+00
 EPSSHRINKSD(%)  4.4012E+01
 EPSSHRINKVR(%)  6.8653E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1663.0205637632598     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.86973719952164     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.07
 Elapsed covariance  time in seconds:     6.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1663.021       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  1.31E+00  7.28E-01  7.97E-01  9.65E-01  1.18E+00  8.91E-01  4.23E-01  1.16E+00  7.54E-01  9.83E-01
 


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
 
         3.45E-02  6.52E-01  4.02E-01  4.33E-01  1.82E-01  7.05E-02  2.75E-01  9.99E-01  4.28E-01  1.50E-01  6.62E-02
 


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
+        1.19E-03
 
 TH 2
+        4.02E-03  4.26E-01
 
 TH 3
+       -2.41E-03 -2.43E-01  1.61E-01
 
 TH 4
+       -2.70E-03 -2.81E-01  1.63E-01  1.87E-01
 
 TH 5
+        8.69E-04  1.04E-01 -4.70E-02 -6.75E-02  3.32E-02
 
 TH 6
+        5.41E-04  8.25E-04 -7.70E-04 -9.66E-04  1.51E-05  4.97E-03
 
 TH 7
+       -2.00E-03 -1.71E-01  9.74E-02  1.12E-01 -4.10E-02  5.94E-04  7.57E-02
 
 TH 8
+       -4.33E-03 -5.35E-01  3.37E-01  3.58E-01 -1.18E-01 -2.17E-03  2.23E-01  9.99E-01
 
 TH 9
+        2.31E-03  2.68E-01 -1.48E-01 -1.77E-01  6.81E-02  6.26E-04 -1.08E-01 -3.44E-01  1.84E-01
 
 TH10
+       -7.32E-04  3.25E-02 -6.34E-03 -2.04E-02  1.65E-02 -6.22E-04 -1.25E-02 -5.78E-02  2.26E-02  2.25E-02
 
 TH11
+       -2.33E-05  1.51E-02 -7.61E-03 -9.74E-03  4.30E-03 -6.33E-05 -5.46E-03 -2.47E-02  9.59E-03  2.22E-03  4.38E-03
 
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
+        3.45E-02
 
 TH 2
+        1.78E-01  6.52E-01
 
 TH 3
+       -1.74E-01 -9.28E-01  4.02E-01
 
 TH 4
+       -1.80E-01 -9.96E-01  9.39E-01  4.33E-01
 
 TH 5
+        1.38E-01  8.72E-01 -6.41E-01 -8.56E-01  1.82E-01
 
 TH 6
+        2.22E-01  1.79E-02 -2.72E-02 -3.16E-02  1.17E-03  7.05E-02
 
 TH 7
+       -2.11E-01 -9.50E-01  8.81E-01  9.42E-01 -8.17E-01  3.06E-02  2.75E-01
 
 TH 8
+       -1.25E-01 -8.20E-01  8.39E-01  8.27E-01 -6.49E-01 -3.07E-02  8.10E-01  9.99E-01
 
 TH 9
+        1.56E-01  9.59E-01 -8.62E-01 -9.56E-01  8.72E-01  2.07E-02 -9.18E-01 -8.04E-01  4.28E-01
 
 TH10
+       -1.41E-01  3.31E-01 -1.05E-01 -3.14E-01  6.03E-01 -5.87E-02 -3.02E-01 -3.85E-01  3.51E-01  1.50E-01
 
 TH11
+       -1.02E-02  3.50E-01 -2.86E-01 -3.40E-01  3.56E-01 -1.36E-02 -3.00E-01 -3.73E-01  3.38E-01  2.24E-01  6.62E-02
 
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
+        1.03E+03
 
 TH 2
+        2.90E+01  5.77E+02
 
 TH 3
+       -4.01E+00  2.61E+02  4.70E+02
 
 TH 4
+        8.82E+00  3.82E+02 -2.80E+02  9.46E+02
 
 TH 5
+       -1.26E+02 -4.44E+02 -6.10E+02  3.10E+02  1.19E+03
 
 TH 6
+       -1.10E+02 -1.78E+00 -3.92E+01  6.50E+01  4.41E+01  2.27E+02
 
 TH 7
+        6.35E+01  1.39E+02  6.68E+01  3.48E+01 -8.56E+01 -3.98E+01  1.70E+02
 
 TH 8
+        7.43E+00 -1.56E+01 -2.37E+01  1.92E+00  5.17E+00  3.89E+00 -9.80E+00  6.62E+00
 
 TH 9
+        3.39E+01 -4.25E+01 -4.44E+01  4.17E+01 -6.32E+00  2.74E-01  1.90E+00  7.29E+00  8.48E+01
 
 TH10
+        1.05E+02  9.40E-01 -2.66E+01 -1.04E+01 -1.51E+02  7.52E+00 -1.48E+01  1.71E+01  2.59E+01  1.48E+02
 
 TH11
+        3.41E+01 -7.10E+01 -4.64E+01 -3.14E+01  8.41E+00  9.80E+00 -4.05E+01  1.21E+01  1.32E+01  2.27E+01  2.92E+02
 
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
 #CPUT: Total CPU Time in Seconds,       42.183
Stop Time:
Wed Sep 29 11:16:16 CDT 2021
