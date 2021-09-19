Sat Sep 18 14:53:27 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.99709950280        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.3241E+00  4.8312E+00  1.1212E+01  6.7115E+00 -5.2864E+01 -7.5674E+00  3.3832E+00 -1.3512E+00  4.1570E+01  1.4923E+01
            -1.4694E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1700.07700374474        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0059E+00  1.0639E+00  9.9618E-01  9.5685E-01  1.0641E+00  9.9812E-01  1.0548E+00  1.0694E+00  7.9249E-01  9.2445E-01
             1.0257E+00
 PARAMETER:  1.0587E-01  1.6196E-01  9.6170E-02  5.5891E-02  1.6218E-01  9.8114E-02  1.5337E-01  1.6707E-01 -1.3257E-01  2.1446E-02
             1.2536E-01
 GRADIENT:   4.6348E+00 -1.9757E-01 -1.1607E+00 -9.9052E-01 -8.7600E-01 -7.9363E+00 -1.8609E+00 -1.8356E+00  3.8658E+00  4.6331E+00
            -6.4189E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1700.93422264516        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0155E+00  9.6712E-01  1.0138E+00  1.0325E+00  1.0184E+00  1.0308E+00  1.1723E+00  1.1088E+00  7.4189E-01  8.5146E-01
             1.0256E+00
 PARAMETER:  1.1537E-01  6.6564E-02  1.1375E-01  1.3201E-01  1.1828E-01  1.3035E-01  2.5899E-01  2.0325E-01 -1.9856E-01 -6.0807E-02
             1.2531E-01
 GRADIENT:   3.1442E+01  1.5094E+01 -6.6716E+00  3.6266E+01  2.5839E+00  7.0329E+00 -5.2662E-01  1.9179E+00  4.2868E+00  1.2420E+00
            -7.0533E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1701.08821634824        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      256
 NPARAMETR:  1.0049E+00  9.4207E-01  7.0265E-01  1.0140E+00  8.3060E-01  1.0178E+00  1.2882E+00  7.0273E-01  6.6195E-01  6.5386E-01
             1.0446E+00
 PARAMETER:  1.0485E-01  4.0323E-02 -2.5289E-01  1.1392E-01 -8.5605E-02  1.1761E-01  3.5328E-01 -2.5278E-01 -3.1257E-01 -3.2487E-01
             1.4360E-01
 GRADIENT:  -4.6581E+01  7.1437E+00 -1.3945E+01  2.1410E+01  2.1228E+01 -6.1516E+00  2.9758E+00  1.3304E+00 -2.9780E+00  1.1146E+00
            -7.1656E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1702.49395740765        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      437
 NPARAMETR:  1.0305E+00  7.7913E-01  1.8969E+00  1.1779E+00  1.2701E+00  1.0402E+00  1.3122E+00  1.8941E+00  6.7377E-01  1.1125E+00
             1.0467E+00
 PARAMETER:  1.3007E-01 -1.4958E-01  7.4021E-01  2.6370E-01  3.3910E-01  1.3941E-01  3.7167E-01  7.3875E-01 -2.9487E-01  2.0662E-01
             1.4562E-01
 GRADIENT:   2.0374E+01  6.6485E+00 -7.2905E+00  1.8944E+01  7.7538E+00  5.2701E+00 -1.9106E+00  2.6771E+00 -2.9545E+00  1.3617E+00
            -3.5032E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1703.23800246237        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      613
 NPARAMETR:  1.0179E+00  6.3770E-01  2.3691E+00  1.2659E+00  1.2928E+00  1.0251E+00  1.5217E+00  2.1435E+00  6.4684E-01  1.1150E+00
             1.0495E+00
 PARAMETER:  1.1778E-01 -3.4989E-01  9.6253E-01  3.3579E-01  3.5685E-01  1.2474E-01  5.1986E-01  8.6243E-01 -3.3566E-01  2.0889E-01
             1.4834E-01
 GRADIENT:  -2.4233E+00  2.5917E-01  1.7826E+00 -6.3792E+00 -2.8319E+00  5.5583E-02 -4.0716E-01  2.3460E-01 -9.1739E-01 -8.3213E-01
            -8.7597E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1703.39575740326        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  1.0218E+00  5.2769E-01  2.7375E+00  1.3469E+00  1.3134E+00  1.0250E+00  1.6535E+00  2.3221E+00  6.5211E-01  1.1456E+00
             1.0484E+00
 PARAMETER:  1.2153E-01 -5.3925E-01  1.1071E+00  3.9779E-01  3.7261E-01  1.2474E-01  6.0291E-01  9.4248E-01 -3.2755E-01  2.3595E-01
             1.4722E-01
 GRADIENT:   6.8669E+00  2.5608E+00  2.6565E+00  2.6288E+00 -4.4238E+00  2.7980E-01  3.1950E-01 -6.1648E-01  1.6627E+00 -8.4015E-02
            -3.4611E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1703.43533609703        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      965
 NPARAMETR:  1.0196E+00  4.3994E-01  3.2413E+00  1.4110E+00  1.3482E+00  1.0245E+00  1.7748E+00  2.5841E+00  6.4344E-01  1.1761E+00
             1.0507E+00
 PARAMETER:  1.1938E-01 -7.2111E-01  1.2760E+00  4.4433E-01  3.9876E-01  1.2416E-01  6.7367E-01  1.0494E+00 -3.4092E-01  2.6222E-01
             1.4947E-01
 GRADIENT:   2.8528E+00  2.3143E+00  2.7956E+00  3.5303E+00 -4.8770E+00  2.5835E-01 -8.3898E-02 -6.4068E-01  7.5651E-01 -1.8964E-01
            -2.2105E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1703.64172752277        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1142
 NPARAMETR:  1.0211E+00  2.7488E-01  4.3132E+00  1.5320E+00  1.3956E+00  1.0238E+00  2.2151E+00  3.0277E+00  6.2657E-01  1.2192E+00
             1.0540E+00
 PARAMETER:  1.2092E-01 -1.1914E+00  1.5617E+00  5.2655E-01  4.3331E-01  1.2354E-01  8.9530E-01  1.2078E+00 -3.6750E-01  2.9819E-01
             1.5257E-01
 GRADIENT:   6.8220E+00  2.7607E+00  3.2502E+00  1.0683E+01 -5.9639E+00  1.8982E-01  6.1228E-01 -3.0222E+00  7.5250E-01  2.1006E-01
             2.2844E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1704.15045528282        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1318
 NPARAMETR:  1.0179E+00  1.3144E-01  5.8034E+00  1.6405E+00  1.4537E+00  1.0251E+00  2.5941E+00  3.6281E+00  6.1933E-01  1.2504E+00
             1.0572E+00
 PARAMETER:  1.1774E-01 -1.9292E+00  1.8584E+00  5.9501E-01  4.7413E-01  1.2475E-01  1.0532E+00  1.3887E+00 -3.7911E-01  3.2345E-01
             1.5565E-01
 GRADIENT:  -5.8344E-01  1.6565E+00 -7.1932E-02  2.1059E+01 -1.3006E+00  3.1281E-01  1.2408E-01  4.1061E-01  4.8409E-01 -3.0870E-01
            -4.8372E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1704.39044462900        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1493
 NPARAMETR:  1.0178E+00  6.0135E-02  6.7761E+00  1.6876E+00  1.4773E+00  1.0241E+00  2.6921E+00  3.9218E+00  6.0981E-01  1.2620E+00
             1.0591E+00
 PARAMETER:  1.1762E-01 -2.7112E+00  2.0134E+00  6.2333E-01  4.9019E-01  1.2386E-01  1.0903E+00  1.4665E+00 -3.9460E-01  3.3268E-01
             1.5743E-01
 GRADIENT:  -9.0664E-01  3.5502E-01  1.7404E-01  7.5166E+00  6.0534E-02 -3.3360E-02  1.0756E-01 -1.5355E-01  8.8931E-02 -1.1520E-01
             3.0990E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1704.43679660080        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1659
 NPARAMETR:  1.0179E+00  5.7006E-02  6.7693E+00  1.6871E+00  1.4769E+00  1.0240E+00  1.8717E+00  3.9246E+00  6.1170E-01  1.2618E+00
             1.0591E+00
 PARAMETER:  1.1776E-01 -2.7646E+00  2.0124E+00  6.2302E-01  4.8994E-01  1.2367E-01  7.2682E-01  1.4673E+00 -3.9152E-01  3.3251E-01
             1.5738E-01
 GRADIENT:  -4.5375E-01  1.0149E-01  9.2981E-02 -2.1653E-01  1.5865E-01 -6.9979E-02  7.5453E-02  2.6945E-01  6.1284E-01 -1.9346E-01
             2.4097E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1704.48139655184        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1835
 NPARAMETR:  1.0182E+00  5.6444E-02  6.7693E+00  1.6871E+00  1.4769E+00  1.0244E+00  1.4961E-01  3.9246E+00  6.1173E-01  1.2618E+00
             1.0587E+00
 PARAMETER:  1.1801E-01 -2.7745E+00  2.0124E+00  6.2302E-01  4.8994E-01  1.2408E-01 -1.7997E+00  1.4673E+00 -3.9146E-01  3.3251E-01
             1.5707E-01
 GRADIENT:   1.1544E-01  2.6221E-02  9.4200E-02 -5.8971E-01  2.2356E-01  1.0644E-01  6.6271E-04  2.3017E-01 -8.7946E-02 -2.2078E-01
             6.8061E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1704.48184177058        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2010
 NPARAMETR:  1.0181E+00  5.5931E-02  6.7692E+00  1.6871E+00  1.4769E+00  1.0241E+00  3.0088E-02  3.9247E+00  6.1175E-01  1.2618E+00
             1.0587E+00
 PARAMETER:  1.1796E-01 -2.7836E+00  2.0124E+00  6.2302E-01  4.8994E-01  1.2382E-01 -3.4036E+00  1.4673E+00 -3.9144E-01  3.3251E-01
             1.5701E-01
 GRADIENT:   4.0163E-02 -1.4612E-03  1.0208E-01 -1.6942E+00  2.8318E-01  1.3582E-02  2.6784E-05  2.4245E-01 -2.3374E-03 -2.2626E-01
             5.5938E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1704.48186121717        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2186
 NPARAMETR:  1.0181E+00  5.5956E-02  6.7692E+00  1.6871E+00  1.4769E+00  1.0241E+00  1.0000E-02  3.9247E+00  6.1176E-01  1.2618E+00
             1.0586E+00
 PARAMETER:  1.1794E-01 -2.7832E+00  2.0124E+00  6.2302E-01  4.8994E-01  1.2380E-01 -4.8744E+00  1.4673E+00 -3.9142E-01  3.3251E-01
             1.5691E-01
 GRADIENT:  -1.3621E-02  1.4462E-04  9.8251E-02 -1.6261E+00  2.6915E-01  2.2673E-03  0.0000E+00  2.5421E-01 -4.0354E-04 -2.3448E-01
             1.8339E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1704.48324523765        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     2353
 NPARAMETR:  1.0181E+00  5.5606E-02  6.7310E+00  1.6875E+00  1.4764E+00  1.0241E+00  1.0000E-02  3.9168E+00  6.1176E-01  1.2629E+00
             1.0584E+00
 PARAMETER:  1.1791E-01 -2.7895E+00  2.0067E+00  6.2325E-01  4.8963E-01  1.2378E-01 -5.0961E+00  1.4653E+00 -3.9141E-01  3.3344E-01
             1.5671E-01
 GRADIENT:  -1.2471E-02  1.9100E-02  1.4718E-02 -6.2348E-01  4.1566E-01 -1.7489E-04  0.0000E+00  3.0599E-01  3.4878E-02 -8.1575E-02
            -1.1553E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1704.58592476955        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2533
 NPARAMETR:  1.0153E+00  2.8429E-02  4.9367E+00  1.6933E+00  1.3932E+00  1.0226E+00  1.0000E-02  3.3566E+00  6.0755E-01  1.2081E+00
             1.0527E+00
 PARAMETER:  1.1519E-01 -3.4604E+00  1.6967E+00  6.2667E-01  4.3159E-01  1.2236E-01 -5.0961E+00  1.3109E+00 -3.9832E-01  2.8908E-01
             1.5137E-01
 GRADIENT:  -6.7047E-01  1.8875E-02  3.0271E-01  9.8091E-02 -1.2555E-01 -1.2025E-01  0.0000E+00 -2.1436E-01 -2.7366E-01  2.7663E-01
            -1.1940E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1704.59437012525        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2709
 NPARAMETR:  1.0153E+00  2.4186E-02  4.6677E+00  1.6935E+00  1.3761E+00  1.0227E+00  1.0000E-02  3.2695E+00  6.0747E-01  1.1941E+00
             1.0514E+00
 PARAMETER:  1.1514E-01 -3.6220E+00  1.6407E+00  6.2679E-01  4.1924E-01  1.2244E-01 -5.0961E+00  1.2846E+00 -3.9845E-01  2.7741E-01
             1.5008E-01
 GRADIENT:   1.5703E-01  1.5311E-02  2.7755E-02  2.0475E-01 -2.3489E-01  2.2057E-02  0.0000E+00  2.1386E-03  6.8404E-02 -7.1343E-02
            -1.5675E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1704.59802091142        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2887
 NPARAMETR:  1.0151E+00  1.1774E-02  4.6967E+00  1.7017E+00  1.3762E+00  1.0224E+00  1.0000E-02  3.2821E+00  6.0423E-01  1.1936E+00
             1.0523E+00
 PARAMETER:  1.1496E-01 -4.3418E+00  1.6469E+00  6.3160E-01  4.1931E-01  1.2220E-01 -5.0961E+00  1.2885E+00 -4.0380E-01  2.7699E-01
             1.5101E-01
 GRADIENT:  -1.1245E-03 -6.5114E-04 -5.8617E-02 -7.6862E-01  3.9406E-01 -7.2732E-02  0.0000E+00 -1.7864E-02 -7.1502E-02  5.8867E-02
             3.1407E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1704.59860826115        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3062
 NPARAMETR:  1.0150E+00  1.0000E-02  4.6899E+00  1.7029E+00  1.3746E+00  1.0227E+00  1.0000E-02  3.2795E+00  6.0397E-01  1.1920E+00
             1.0516E+00
 PARAMETER:  1.1492E-01 -4.5091E+00  1.6454E+00  6.3232E-01  4.1817E-01  1.2240E-01 -5.0961E+00  1.2877E+00 -4.0424E-01  2.7564E-01
             1.5030E-01
 GRADIENT:  -9.8847E-03  2.5507E-04  4.3003E-03 -1.5101E-01  6.4121E-03  6.0949E-03  0.0000E+00 -1.1801E-03  1.1365E-02 -1.7778E-03
            -8.7473E-04

0ITERATION NO.:   96    OBJECTIVE VALUE:  -1704.59860826115        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     3084
 NPARAMETR:  1.0150E+00  1.0000E-02  4.6899E+00  1.7029E+00  1.3746E+00  1.0227E+00  1.0000E-02  3.2795E+00  6.0397E-01  1.1920E+00
             1.0516E+00
 PARAMETER:  1.1492E-01 -4.5091E+00  1.6454E+00  6.3232E-01  4.1817E-01  1.2240E-01 -5.0961E+00  1.2877E+00 -4.0424E-01  2.7564E-01
             1.5030E-01
 GRADIENT:  -9.8847E-03  2.5507E-04  4.3003E-03 -1.5101E-01  6.4121E-03  6.0949E-03  0.0000E+00 -1.1801E-03  1.1365E-02 -1.7778E-03
            -8.7473E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3084
 NO. OF SIG. DIGITS IN FINAL EST.:  4.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.7517E-03 -4.1910E-06 -4.9477E-02 -1.4352E-02 -6.1899E-02
 SE:             2.9777E-02  1.9564E-06  1.8178E-02  2.8618E-02  1.9742E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5309E-01  3.2177E-02  6.4919E-03  6.1601E-01  1.7167E-03

 ETASHRINKSD(%)  2.4206E-01  9.9993E+01  3.9102E+01  4.1264E+00  3.3861E+01
 ETASHRINKVR(%)  4.8353E-01  1.0000E+02  6.2914E+01  8.0826E+00  5.6256E+01
 EBVSHRINKSD(%)  4.5532E-01  9.9994E+01  4.7681E+01  4.0126E+00  2.8727E+01
 EBVSHRINKVR(%)  9.0858E-01  1.0000E+02  7.2628E+01  7.8642E+00  4.9202E+01
 RELATIVEINF(%)  9.7906E+01  2.1787E-08  1.3734E+01  5.1460E+00  2.1798E+01
 EPSSHRINKSD(%)  4.4090E+01
 EPSSHRINKVR(%)  6.8741E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1704.5986082611487     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -969.44778169741051     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.67
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1704.599       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E-02  4.69E+00  1.70E+00  1.37E+00  1.02E+00  1.00E-02  3.28E+00  6.04E-01  1.19E+00  1.05E+00
 


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
+        1.02E+03
 
 TH 2
+       -4.45E+01  7.80E+02
 
 TH 3
+       -1.81E+00 -1.64E+00  1.47E+00
 
 TH 4
+       -2.00E+01  4.50E+01 -5.24E+00  9.77E+02
 
 TH 5
+       -1.12E+00  1.66E+01 -9.05E+00 -4.67E+01  2.40E+02
 
 TH 6
+        2.71E+01  9.71E+01 -2.57E-01 -4.20E+00 -1.04E+01  2.01E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -9.82E-01  1.86E+00 -2.10E+00 -3.54E+00 -9.32E+00 -1.02E+00  0.00E+00  7.41E+00
 
 TH 9
+        7.95E+00 -2.91E-01  1.36E+00  6.10E+00 -2.83E+00 -2.14E+01  0.00E+00 -7.53E-01  4.77E+02
 
 TH10
+       -7.29E-01  2.14E+01  2.58E-01 -7.61E+00 -4.11E+01 -1.66E+01  0.00E+00  1.72E+00 -1.07E+01  5.61E+01
 
 TH11
+        9.06E+00 -3.78E+00 -1.83E+00 -2.18E+01  6.62E+00  1.83E+01  0.00E+00 -1.02E+00  2.43E+01  9.37E+00  2.10E+02
 
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
 #CPUT: Total CPU Time in Seconds,       51.139
Stop Time:
Sat Sep 18 14:54:20 CDT 2021
