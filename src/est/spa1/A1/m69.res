Wed Sep 29 22:44:20 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat69.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1559.01834686266        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4124E+02  9.9446E+01 -2.7423E+01  1.9786E+02  1.7626E+02  4.3206E+01 -3.3660E+01  7.4097E+00 -3.0925E+01 -7.0457E+01
            -9.8005E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1819.45255520488        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0284E+00  8.4301E-01  9.9183E-01  1.0110E+00  8.0075E-01  9.7461E-01  1.1465E+00  8.5476E-01  1.1747E+00  1.0985E+00
             1.9817E+00
 PARAMETER:  1.2805E-01 -7.0772E-02  9.1798E-02  1.1097E-01 -1.2221E-01  7.4278E-02  2.3670E-01 -5.6930E-02  2.6105E-01  1.9392E-01
             7.8395E-01
 GRADIENT:   1.5547E+02  1.5622E+01  2.5254E+01 -1.2603E+01 -5.4368E+01  1.7689E+00  4.0525E+00  9.6583E+00  2.8950E+01  1.4577E+01
             7.3279E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1828.32240150948        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      264
 NPARAMETR:  1.0550E+00  6.0222E-01  7.5686E-01  1.2091E+00  6.6083E-01  1.0128E+00  2.2882E+00  2.6405E-01  8.9812E-01  8.5373E-01
             1.9188E+00
 PARAMETER:  1.5355E-01 -4.0713E-01 -1.7858E-01  2.8985E-01 -3.1426E-01  1.1276E-01  9.2775E-01 -1.2316E+00 -7.4489E-03 -5.8141E-02
             7.5168E-01
 GRADIENT:   5.1286E+01  3.7149E+01 -3.2763E+01  6.4930E+01  4.3935E+01  2.3740E+00  2.8642E+01  1.5601E+00 -8.6056E+00 -7.8524E-01
             4.6580E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1839.78691802214        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.0286E+00  5.1279E-01  6.4831E-01  1.1889E+00  5.6757E-01  9.8589E-01  1.7197E+00  6.5356E-02  9.4732E-01  7.6168E-01
             1.7752E+00
 PARAMETER:  1.2824E-01 -5.6788E-01 -3.3338E-01  2.7306E-01 -4.6640E-01  8.5786E-02  6.4213E-01 -2.6279E+00  4.5885E-02 -1.7223E-01
             6.7389E-01
 GRADIENT:   4.8130E-01  1.1022E+01  3.9207E+00  1.7604E+01 -3.5055E+00 -6.9374E+00 -3.4373E+00  1.0903E-01 -1.2272E+00 -4.5240E+00
            -1.9052E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1842.01271119395        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  1.0240E+00  3.1737E-01  6.8565E-01  1.2837E+00  5.4353E-01  9.9732E-01  2.4363E+00  1.0000E-02  9.0543E-01  7.9383E-01
             1.8120E+00
 PARAMETER:  1.2373E-01 -1.0477E+00 -2.7738E-01  3.4976E-01 -5.0967E-01  9.7314E-02  9.9046E-01 -5.0256E+00  6.5908E-04 -1.3089E-01
             6.9441E-01
 GRADIENT:   1.2544E+00  1.0963E+00  9.4733E-01  1.3565E+00 -1.8625E+00 -5.9084E-01  4.3236E-01  0.0000E+00  6.1573E-02 -2.0933E-01
            -1.0147E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1842.56319517648        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  1.0201E+00  2.1002E-01  8.2041E-01  1.3771E+00  5.8537E-01  9.9057E-01  2.5965E+00  1.0000E-02  8.8389E-01  9.1283E-01
             1.8161E+00
 PARAMETER:  1.1993E-01 -1.4606E+00 -9.7957E-02  4.2000E-01 -4.3552E-01  9.0526E-02  1.0542E+00 -7.2517E+00 -2.3417E-02  8.7938E-03
             6.9670E-01
 GRADIENT:   9.3068E-01  2.7890E+00  1.8758E+00  1.0723E+01 -5.1484E+00 -2.0192E+00  5.4876E-01  0.0000E+00 -2.4811E-01 -9.0960E-01
            -1.7313E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1843.30128898694        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      973
 NPARAMETR:  1.0134E+00  5.7678E-02  8.7789E-01  1.4684E+00  5.7982E-01  9.9668E-01  2.8553E+00  1.0000E-02  8.5547E-01  9.5534E-01
             1.8264E+00
 PARAMETER:  1.1333E-01 -2.7529E+00 -3.0231E-02  4.8415E-01 -4.4504E-01  9.6678E-02  1.1492E+00 -1.6126E+01 -5.6105E-02  5.4309E-02
             7.0234E-01
 GRADIENT:  -3.9548E+00  1.0078E+00  4.6193E+00  1.4547E+01 -1.0879E+01  1.4233E+00 -4.0051E-01  0.0000E+00 -6.0528E-01  1.6740E+00
             1.7807E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1843.77464787623        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1151
 NPARAMETR:  1.0132E+00  1.0000E-02  8.9172E-01  1.4878E+00  5.8074E-01  9.9199E-01  3.0751E+00  1.0000E-02  8.4634E-01  9.5638E-01
             1.8239E+00
 PARAMETER:  1.1314E-01 -4.5434E+00 -1.4603E-02  4.9733E-01 -4.4345E-01  9.1955E-02  1.2233E+00 -2.9697E+01 -6.6830E-02  5.5398E-02
             7.0099E-01
 GRADIENT:  -7.1189E-01  1.0479E-02 -6.2473E-01 -1.8323E+00  1.0899E+00  1.7072E-02 -1.9828E-02  0.0000E+00  9.4058E-01  2.7411E-01
             5.6725E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1843.78019383870        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1342
 NPARAMETR:  1.0138E+00  1.0000E-02  8.9107E-01  1.4860E+00  5.8061E-01  9.9198E-01  3.2350E+00  1.0000E-02  8.4449E-01  9.5534E-01
             1.8231E+00
 PARAMETER:  1.1374E-01 -4.5588E+00 -1.5327E-02  4.9609E-01 -4.4367E-01  9.1944E-02  1.2740E+00 -2.9808E+01 -6.9024E-02  5.4315E-02
             7.0052E-01
 GRADIENT:   6.8176E-01  0.0000E+00 -6.1496E-01 -5.1790E+00  1.7021E+00  2.7273E-02 -2.1825E-02  0.0000E+00  1.8232E-01  7.0125E-02
             9.0014E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1843.78269006603        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1535
 NPARAMETR:  1.0138E+00  1.0000E-02  8.9114E-01  1.4857E+00  5.7987E-01  9.9198E-01  3.4513E+00  1.0000E-02  8.4409E-01  9.5476E-01
             1.8229E+00
 PARAMETER:  1.1372E-01 -4.5588E+00 -1.5249E-02  4.9587E-01 -4.4495E-01  9.1945E-02  1.3388E+00 -2.9808E+01 -6.9493E-02  5.3700E-02
             7.0044E-01
 GRADIENT:   6.1845E-01  0.0000E+00  8.4469E-01 -5.6733E+00 -4.1799E-01  3.1033E-02 -2.4328E-02  0.0000E+00  3.6783E-02  4.0646E-02
            -4.8901E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1843.78492096369        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1725
 NPARAMETR:  1.0138E+00  1.0000E-02  8.8943E-01  1.4856E+00  5.7978E-01  9.9199E-01  3.6931E+00  1.0000E-02  8.4408E-01  9.5426E-01
             1.8230E+00
 PARAMETER:  1.1373E-01 -4.5588E+00 -1.7170E-02  4.9579E-01 -4.4511E-01  9.1955E-02  1.4065E+00 -2.9808E+01 -6.9506E-02  5.3180E-02
             7.0050E-01
 GRADIENT:   6.8344E-01  0.0000E+00 -4.6981E-01 -5.3434E+00  1.3598E+00  3.5504E-02 -2.7070E-02  0.0000E+00  8.8867E-03  3.0733E-02
             7.1670E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1843.78731975939        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1915
 NPARAMETR:  1.0134E+00  1.0000E-02  8.8916E-01  1.4858E+00  5.7923E-01  9.9182E-01  3.9471E+00  1.0000E-02  8.4401E-01  9.5383E-01
             1.8230E+00
 PARAMETER:  1.1336E-01 -4.5588E+00 -1.7474E-02  4.9595E-01 -4.4606E-01  9.1790E-02  1.4730E+00 -2.9808E+01 -6.9588E-02  5.2733E-02
             7.0047E-01
 GRADIENT:  -1.6936E-01  0.0000E+00  2.5281E-01 -4.6997E+00  1.1246E-01 -2.9446E-02 -3.0047E-02  0.0000E+00 -1.2476E-02  3.0135E-02
            -5.5436E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1843.79170536149        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2107
 NPARAMETR:  1.0129E+00  1.0000E-02  8.8788E-01  1.4860E+00  5.7869E-01  9.9159E-01  4.5706E+00  1.0000E-02  8.4385E-01  9.5318E-01
             1.8230E+00
 PARAMETER:  1.1281E-01 -4.5588E+00 -1.8918E-02  4.9609E-01 -4.4699E-01  9.1555E-02  1.6196E+00 -2.9808E+01 -6.9782E-02  5.2052E-02
             7.0046E-01
 GRADIENT:  -1.3793E+00  0.0000E+00  6.7078E-02 -3.7580E+00  9.1132E-02 -1.2358E-01 -3.7294E-02  0.0000E+00 -7.5255E-02  3.0060E-02
            -8.8939E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1843.81466317237        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2290
 NPARAMETR:  1.0098E+00  1.0000E-02  8.8264E-01  1.4886E+00  5.7684E-01  9.9076E-01  1.0874E+01  1.0000E-02  8.4331E-01  9.4887E-01
             1.8235E+00
 PARAMETER:  1.0975E-01 -4.5588E+00 -2.4839E-02  4.9783E-01 -4.5019E-01  9.0716E-02  2.4864E+00 -2.9808E+01 -7.0423E-02  4.7514E-02
             7.0077E-01
 GRADIENT:  -8.2524E+00  0.0000E+00 -1.6606E+00  3.0322E+00  1.1442E+00 -5.0681E-01  3.4739E-01  0.0000E+00 -1.2064E-01 -3.8788E-01
            -6.8062E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1843.84545217496        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     2454
 NPARAMETR:  1.0136E+00  1.0000E-02  8.8387E-01  1.4861E+00  5.7697E-01  9.9183E-01  1.0658E+01  1.0000E-02  8.4314E-01  9.5061E-01
             1.8233E+00
 PARAMETER:  1.1348E-01 -4.5588E+00 -2.3450E-02  4.9616E-01 -4.4997E-01  9.1801E-02  2.4663E+00 -2.9808E+01 -7.0620E-02  4.9352E-02
             7.0065E-01
 GRADIENT:   1.2934E-01  0.0000E+00 -3.4686E-01 -2.1104E+00  5.6420E-02 -2.1832E-02 -2.1006E-02  0.0000E+00 -4.3589E-02  9.2924E-03
            -7.2817E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2454
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.8616E-04  2.4583E-04 -6.2178E-05 -6.7153E-03 -1.5828E-02
 SE:             2.9610E-02  1.8990E-03  1.6991E-04  2.8777E-02  2.4094E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8151E-01  8.9700E-01  7.1441E-01  8.1548E-01  5.1121E-01

 ETASHRINKSD(%)  8.0258E-01  9.3638E+01  9.9431E+01  3.5951E+00  1.9284E+01
 ETASHRINKVR(%)  1.5987E+00  9.9595E+01  9.9997E+01  7.0609E+00  3.4849E+01
 EBVSHRINKSD(%)  1.0261E+00  9.4449E+01  9.9400E+01  3.4299E+00  1.8106E+01
 EBVSHRINKVR(%)  2.0417E+00  9.9692E+01  9.9996E+01  6.7421E+00  3.2934E+01
 RELATIVEINF(%)  9.2522E+01  1.8701E-02  3.4365E-04  7.2577E+00  6.7465E+00
 EPSSHRINKSD(%)  2.9326E+01
 EPSSHRINKVR(%)  5.0051E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1843.8454521749645     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -924.90691897029183     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1843.845       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.00E-02  8.84E-01  1.49E+00  5.77E-01  9.92E-01  1.07E+01  1.00E-02  8.43E-01  9.51E-01  1.82E+00
 


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
+        1.08E+03
 
 TH 2
+        0.00E+00  2.04E+03
 
 TH 3
+       -1.09E+01  0.00E+00  5.11E+02
 
 TH 4
+       -1.49E+01  0.00E+00 -1.11E+02  6.62E+02
 
 TH 5
+        1.92E+01  0.00E+00 -1.02E+03 -6.67E+01  2.49E+03
 
 TH 6
+        1.83E+00  0.00E+00  2.62E-01 -4.23E+00 -1.13E+00  1.94E+02
 
 TH 7
+       -2.04E-01  0.00E+00  4.11E+00 -5.35E+01 -8.81E+00  1.55E-01  2.90E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.40E+00  0.00E+00  1.47E+01 -5.30E+00 -3.65E+00 -6.27E-01 -7.81E-01  0.00E+00  2.44E+02
 
 TH10
+        5.84E-01  0.00E+00 -2.35E+01  4.49E+00 -6.90E+01  4.65E-01  3.29E+00  0.00E+00 -1.11E+00  1.07E+02
 
 TH11
+       -1.18E+01  0.00E+00 -5.38E+00 -1.14E+01  1.14E+01  2.32E+00 -3.16E+01  0.00E+00  7.45E+00  1.90E+01  1.34E+02
 
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
 #CPUT: Total CPU Time in Seconds,       52.558
Stop Time:
Wed Sep 29 22:45:14 CDT 2021
