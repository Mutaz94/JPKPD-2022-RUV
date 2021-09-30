Wed Sep 29 15:15:07 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1738.03879270502        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0916E+02  4.9364E+01 -2.2004E+01  1.2445E+02  3.0538E+01  4.7547E+01  2.6070E+01  9.3468E+00  5.1187E+01 -2.1230E-01
             6.8656E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1753.73014464573        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0275E+00  9.6864E-01  1.0165E+00  9.8196E-01  9.8721E-01  1.0002E+00  8.3522E-01  9.4744E-01  7.8424E-01  1.0684E+00
             8.7446E-01
 PARAMETER:  1.2712E-01  6.8137E-02  1.1632E-01  8.1791E-02  8.7132E-02  1.0020E-01 -8.0057E-02  4.6012E-02 -1.4304E-01  1.6612E-01
            -3.4145E-02
 GRADIENT:   2.7986E+01 -1.4714E+01 -7.2168E+00 -2.7215E+01  5.1674E-01 -3.3416E-01  4.1968E+00  6.0773E+00 -2.9393E+00  5.8601E+00
             2.0715E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1755.08951535729        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0287E+00  8.7940E-01  9.4600E-01  1.0448E+00  9.1274E-01  9.8173E-01  7.0360E-01  6.5927E-01  8.2577E-01  1.0186E+00
             8.5318E-01
 PARAMETER:  1.2829E-01 -2.8517E-02  4.4486E-02  1.4378E-01  8.6956E-03  8.1563E-02 -2.5155E-01 -3.1662E-01 -9.1436E-02  1.1845E-01
            -5.8780E-02
 GRADIENT:   3.1328E+01  1.0985E+00 -6.0748E+00  1.0931E+01  2.9972E+00 -8.1683E+00 -1.7233E-01  6.6967E-01  7.4964E+00  6.5379E-01
             9.6343E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1755.85241557901        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      543
 NPARAMETR:  1.0109E+00  7.8058E-01  1.0732E+00  1.1107E+00  9.2948E-01  9.9357E-01  7.9193E-01  7.7642E-01  7.5649E-01  1.0416E+00
             8.3138E-01
 PARAMETER:  1.1087E-01 -1.4771E-01  1.7068E-01  2.0496E-01  2.6875E-02  9.3545E-02 -1.3328E-01 -1.5307E-01 -1.7907E-01  1.4072E-01
            -8.4672E-02
 GRADIENT:  -4.7233E+00  5.6924E+00  2.2663E+00  4.7679E+00 -3.7355E+00 -2.2354E+00 -2.3360E-01 -2.0583E-01 -8.0705E-01 -7.1071E-01
            -1.2305E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1756.40592187019        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0106E+00  4.7842E-01  1.3081E+00  1.3112E+00  9.2012E-01  1.0044E+00  1.0843E+00  9.7803E-01  6.5969E-01  1.0712E+00
             8.3665E-01
 PARAMETER:  1.1052E-01 -6.3726E-01  3.6855E-01  3.7094E-01  1.6746E-02  1.0442E-01  1.8095E-01  7.7790E-02 -3.1599E-01  1.6879E-01
            -7.8346E-02
 GRADIENT:   5.0764E+00  8.1776E+00  1.8147E+00  1.9623E+01 -4.4791E+00  4.1721E+00  3.9589E-01 -2.0211E-01  3.4603E-02  5.6992E-01
             1.1799E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1756.83717366782        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0065E+00  2.9101E-01  1.4615E+00  1.4376E+00  9.2018E-01  9.9783E-01  1.2624E+00  1.1294E+00  6.1584E-01  1.0979E+00
             8.3614E-01
 PARAMETER:  1.0648E-01 -1.1344E+00  4.7948E-01  4.6299E-01  1.6816E-02  9.7831E-02  3.3298E-01  2.2168E-01 -3.8477E-01  1.9343E-01
            -7.8960E-02
 GRADIENT:   3.1155E+00  8.2048E+00 -1.8393E+00  4.0217E+01 -3.1553E+00  3.0609E+00 -9.0159E-02 -2.3988E-01 -1.5306E+00  2.3532E+00
             8.6125E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1757.31222348467        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  1.0015E+00  1.5435E-01  1.6328E+00  1.5305E+00  9.3492E-01  9.8380E-01  1.2777E+00  1.3107E+00  5.9044E-01  1.1158E+00
             8.3252E-01
 PARAMETER:  1.0151E-01 -1.7686E+00  5.9027E-01  5.2560E-01  3.2704E-02  8.3665E-02  3.4510E-01  3.7056E-01 -4.2689E-01  2.0955E-01
            -8.3295E-02
 GRADIENT:  -2.9948E+00  5.3702E+00 -3.3999E+00  4.7331E+01 -3.5137E+00 -1.4100E+00 -6.4003E-02  5.3421E-01 -1.1572E+00  2.5316E+00
            -5.9610E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1757.93719303682        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  1.0005E+00  6.7223E-02  1.8479E+00  1.5873E+00  9.6814E-01  9.8206E-01  1.0903E+00  1.4997E+00  5.7051E-01  1.1299E+00
             8.3207E-01
 PARAMETER:  1.0048E-01 -2.5997E+00  7.1405E-01  5.6203E-01  6.7624E-02  8.1901E-02  1.8647E-01  5.0528E-01 -4.6123E-01  2.2217E-01
            -8.3839E-02
 GRADIENT:  -2.3085E+00  1.8700E+00 -1.2714E+00  2.1608E+01 -1.7817E+00 -1.3366E+00 -1.7476E-03  3.0274E-01 -5.5751E-01  7.3045E-01
            -5.7296E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1758.10590882234        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1427
 NPARAMETR:  1.0009E+00  3.0551E-02  2.0084E+00  1.6141E+00  9.9409E-01  9.8383E-01  8.8047E-01  1.6256E+00  5.5897E-01  1.1422E+00
             8.3274E-01
 PARAMETER:  1.0092E-01 -3.3884E+00  7.9734E-01  5.7880E-01  9.4071E-02  8.3693E-02 -2.7303E-02  5.8588E-01 -4.8167E-01  2.3299E-01
            -8.3035E-02
 GRADIENT:  -3.6095E-01  7.8667E-01  7.3879E-01  1.1302E+01 -2.0748E+00 -2.7842E-01 -2.7547E-05 -9.6693E-02 -1.2856E+00 -2.0229E-01
            -2.9369E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1758.31302358569        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1611
 NPARAMETR:  1.0017E+00  1.0000E-02  2.0499E+00  1.6152E+00  1.0023E+00  9.8437E-01  7.7904E-01  1.6653E+00  5.5827E-01  1.1515E+00
             8.3253E-01
 PARAMETER:  1.0169E-01 -4.7326E+00  8.1778E-01  5.7947E-01  1.0228E-01  8.4246E-02 -1.4969E-01  6.1002E-01 -4.8291E-01  2.4108E-01
            -8.3288E-02
 GRADIENT:   2.4715E+00  0.0000E+00  2.6688E-01 -4.5511E+01  3.6207E+00  1.5881E-01  2.2608E-04  6.3263E-01  1.9125E+00  2.9811E-01
             3.3785E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1758.34526735344        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1781             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0024E+00  1.0000E-02  2.0303E+00  1.6181E+00  9.9824E-01  9.8448E-01  7.7533E-01  1.6538E+00  5.5570E-01  1.1438E+00
             8.3230E-01
 PARAMETER:  1.0235E-01 -4.7326E+00  8.0816E-01  5.8128E-01  9.8239E-02  8.4362E-02 -1.5447E-01  6.0309E-01 -4.8753E-01  2.3432E-01
            -8.3563E-02
 GRADIENT:   6.5809E+02  0.0000E+00  1.1494E+01  1.6911E+03  1.2877E+01  6.3901E+01  6.1664E-04  3.7195E+00  4.0307E+01  3.3546E+00
             7.8057E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1758.35253773359        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1961
 NPARAMETR:  1.0015E+00  1.0000E-02  2.0290E+00  1.6192E+00  9.9363E-01  9.8436E-01  7.7560E-01  1.6419E+00  5.5577E-01  1.1446E+00
             8.3239E-01
 PARAMETER:  1.0149E-01 -4.7326E+00  8.0754E-01  5.8196E-01  9.3611E-02  8.4238E-02 -1.5411E-01  5.9586E-01 -4.8740E-01  2.3510E-01
            -8.3452E-02
 GRADIENT:   1.8718E+00  0.0000E+00  1.5183E+00 -2.4449E+01 -8.2664E-01  1.1459E-01  1.2405E-04  6.9860E-02  2.0728E-01  2.3884E-01
            -6.5390E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1758.36648938635        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2152
 NPARAMETR:  1.0017E+00  1.0000E-02  2.0164E+00  1.6171E+00  9.9032E-01  9.8446E-01  7.7001E-01  1.6293E+00  5.5601E-01  1.1431E+00
             8.3238E-01
 PARAMETER:  1.0168E-01 -4.7326E+00  8.0134E-01  5.8065E-01  9.0268E-02  8.4340E-02 -1.6135E-01  5.8813E-01 -4.8697E-01  2.3375E-01
            -8.3470E-02
 GRADIENT:   2.5206E+00  0.0000E+00  2.0331E+00 -3.0953E+01 -1.2173E+00  1.5682E-01  1.3690E-04 -7.3258E-02  2.8119E-01  3.4725E-01
            -4.0389E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1758.37679525851        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2350             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0016E+00  1.0000E-02  1.9881E+00  1.6164E+00  9.8942E-01  9.8447E-01  7.7130E-01  1.6215E+00  5.5614E-01  1.1376E+00
             8.3241E-01
 PARAMETER:  1.0164E-01 -4.7326E+00  7.8720E-01  5.8022E-01  8.9368E-02  8.4350E-02 -1.5968E-01  5.8336E-01 -4.8674E-01  2.2888E-01
            -8.3429E-02
 GRADIENT:   6.5108E+02  0.0000E+00  1.1005E+01  1.6842E+03  1.3383E+01  6.3873E+01  6.2403E-04  3.5931E+00  4.0239E+01  3.2577E+00
             8.5738E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1758.38336128258        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2536
 NPARAMETR:  1.0014E+00  1.0000E-02  1.9898E+00  1.6176E+00  9.8489E-01  9.8440E-01  7.6864E-01  1.6102E+00  5.5618E-01  1.1387E+00
             8.3240E-01
 PARAMETER:  1.0141E-01 -4.7326E+00  7.8803E-01  5.8094E-01  8.4778E-02  8.4277E-02 -1.6314E-01  5.7639E-01 -4.8666E-01  2.2989E-01
            -8.3446E-02
 GRADIENT:   1.9440E+00  0.0000E+00  1.3963E+00 -2.4658E+01 -1.0204E+00  1.1382E-01  1.1952E-04 -5.6753E-03  2.1564E-01  2.7005E-01
            -6.7093E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1758.39361191162        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2725
 NPARAMETR:  1.0016E+00  1.0000E-02  1.9803E+00  1.6156E+00  9.8240E-01  9.8450E-01  7.6813E-01  1.6025E+00  5.5639E-01  1.1372E+00
             8.3241E-01
 PARAMETER:  1.0161E-01 -4.7326E+00  7.8323E-01  5.7973E-01  8.2246E-02  8.4378E-02 -1.6379E-01  5.7156E-01 -4.8629E-01  2.2857E-01
            -8.3426E-02
 GRADIENT:   2.4026E+00  0.0000E+00  1.6896E+00 -3.0980E+01 -1.1832E+00  1.6262E-01  1.3343E-04  7.9606E-03  2.8380E-01  3.3715E-01
            -1.0991E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1758.40286618649        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2921             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0016E+00  1.0000E-02  1.9614E+00  1.6151E+00  9.8148E-01  9.8450E-01  7.6931E-01  1.5956E+00  5.5648E-01  1.1331E+00
             8.3240E-01
 PARAMETER:  1.0157E-01 -4.7326E+00  7.7367E-01  5.7940E-01  8.1310E-02  8.4382E-02 -1.6226E-01  5.6723E-01 -4.8612E-01  2.2493E-01
            -8.3440E-02
 GRADIENT:   6.5007E+02  0.0000E+00  1.1945E+01  1.6793E+03  1.1547E+01  6.3796E+01  6.2665E-04  3.1905E+00  4.0213E+01  3.3602E+00
             8.1898E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1758.40853845921        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3113
 NPARAMETR:  1.0016E+00  1.0000E-02  1.9533E+00  1.6148E+00  9.7937E-01  9.8451E-01  7.6105E-01  1.5885E+00  5.5657E-01  1.1318E+00
             8.3240E-01
 PARAMETER:  1.0156E-01 -4.7326E+00  7.6953E-01  5.7920E-01  7.9155E-02  8.4387E-02 -1.7306E-01  5.6279E-01 -4.8596E-01  2.2384E-01
            -8.3439E-02
 GRADIENT:   2.4313E+00  0.0000E+00 -2.2100E-01 -3.0703E+01  2.0587E+00  1.4792E-01  1.2797E-04  3.7885E-01  2.3384E-01 -1.3809E-01
             1.2900E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1758.41493528932        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     3309             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0015E+00  1.0000E-02  1.9464E+00  1.6142E+00  9.7489E-01  9.8452E-01  7.6574E-01  1.5775E+00  5.5673E-01  1.1314E+00
             8.3239E-01
 PARAMETER:  1.0153E-01 -4.7326E+00  7.6596E-01  5.7886E-01  7.4570E-02  8.4402E-02 -1.6692E-01  5.5584E-01 -4.8567E-01  2.2348E-01
            -8.3456E-02
 GRADIENT:   6.4952E+02  0.0000E+00  1.3602E+01  1.6764E+03  8.1938E+00  6.3765E+01  6.3558E-04  2.7957E+00  4.0218E+01  3.8019E+00
             7.8115E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1758.41947323453        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3501
 NPARAMETR:  1.0015E+00  1.0000E-02  1.9387E+00  1.6140E+00  9.7382E-01  9.8453E-01  7.6241E-01  1.5726E+00  5.5679E-01  1.1302E+00
             8.3241E-01
 PARAMETER:  1.0152E-01 -4.7326E+00  7.6203E-01  5.7870E-01  7.3472E-02  8.4407E-02 -1.7127E-01  5.5276E-01 -4.8557E-01  2.2239E-01
            -8.3435E-02
 GRADIENT:   2.4391E+00  0.0000E+00  8.7312E-01 -3.0851E+01 -3.4660E-01  1.6338E-01  1.2941E-04  1.6229E-01  2.6901E-01  2.3217E-01
            -5.3324E-03

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1758.42294557110        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     3699             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0015E+00  1.0000E-02  1.9247E+00  1.6135E+00  9.7313E-01  9.8454E-01  7.5703E-01  1.5655E+00  5.5689E-01  1.1263E+00
             8.3243E-01
 PARAMETER:  1.0151E-01 -4.7326E+00  7.5474E-01  5.7839E-01  7.2761E-02  8.4422E-02 -1.7835E-01  5.4822E-01 -4.8539E-01  2.1892E-01
            -8.3405E-02
 GRADIENT:   6.4899E+02  0.0000E+00  1.1729E+01  1.6729E+03  1.1861E+01  6.3736E+01  6.5533E-04  2.9576E+00  4.0141E+01  3.1179E+00
             8.0702E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1758.42673777714        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3892
 NPARAMETR:  1.0015E+00  1.0000E-02  1.9176E+00  1.6131E+00  9.7139E-01  9.8456E-01  7.5592E-01  1.5590E+00  5.5698E-01  1.1254E+00
             8.3242E-01
 PARAMETER:  1.0149E-01 -4.7326E+00  7.5107E-01  5.7817E-01  7.0975E-02  8.4435E-02 -1.7982E-01  5.4405E-01 -4.8523E-01  2.1817E-01
            -8.3416E-02
 GRADIENT:   2.2645E+00  0.0000E+00 -5.3728E-01 -3.1117E+01  2.4514E+00  1.5391E-01  1.2628E-04  2.7590E-01  2.3607E-01 -2.6675E-01
            -2.1995E-03

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1758.43150042632        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     4088             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0015E+00  1.0000E-02  1.9140E+00  1.6128E+00  9.6787E-01  9.8456E-01  7.4645E-01  1.5521E+00  5.5708E-01  1.1258E+00
             8.3243E-01
 PARAMETER:  1.0147E-01 -4.7326E+00  7.4921E-01  5.7799E-01  6.7345E-02  8.4437E-02 -1.9242E-01  5.3963E-01 -4.8505E-01  2.1852E-01
            -8.3405E-02
 GRADIENT:   6.4837E+02  0.0000E+00  1.3204E+01  1.6709E+03  8.7587E+00  6.3688E+01  6.7276E-04  2.6703E+00  4.0150E+01  3.6199E+00
             7.8871E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1758.43313385777        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     4283
 NPARAMETR:  1.0015E+00  1.0000E-02  1.9092E+00  1.6126E+00  9.6620E-01  9.8456E-01  7.4647E-01  1.5475E+00  5.5714E-01  1.1250E+00
             8.3243E-01
 PARAMETER:  1.0147E-01 -4.7326E+00  7.4668E-01  5.7784E-01  6.5611E-02  8.4444E-02 -1.9240E-01  5.3665E-01 -4.8493E-01  2.1774E-01
            -8.3404E-02
 GRADIENT:   2.4316E+00  0.0000E+00  1.1112E+00 -3.1001E+01 -1.1390E+00  1.7297E-01  1.2100E-04  4.9950E-02  2.7897E-01  2.7503E-01
            -1.8845E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1758.43647034694        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     4480             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0015E+00  1.0000E-02  1.8994E+00  1.6123E+00  9.6612E-01  9.8457E-01  7.4449E-01  1.5435E+00  5.5719E-01  1.1226E+00
             8.3246E-01
 PARAMETER:  1.0145E-01 -4.7326E+00  7.4151E-01  5.7765E-01  6.5529E-02  8.4452E-02 -1.9505E-01  5.3403E-01 -4.8484E-01  2.1568E-01
            -8.3364E-02
 GRADIENT:   6.4788E+02  0.0000E+00  1.2156E+01  1.6685E+03  1.0684E+01  6.3649E+01  6.7516E-04  2.7684E+00  4.0109E+01  3.2771E+00
             8.0982E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1758.43697861753        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     4672
 NPARAMETR:  1.0014E+00  1.0000E-02  1.8938E+00  1.6121E+00  9.6551E-01  9.8458E-01  7.4301E-01  1.5399E+00  5.5724E-01  1.1212E+00
             8.3248E-01
 PARAMETER:  1.0145E-01 -4.7326E+00  7.3857E-01  5.7751E-01  6.4898E-02  8.4455E-02 -1.9705E-01  5.3172E-01 -4.8475E-01  2.1443E-01
            -8.3346E-02
 GRADIENT:   2.3732E+00  0.0000E+00 -6.1317E-01 -3.1119E+01  2.1841E+00  1.5094E-01  1.2001E-04  2.6660E-01  2.5144E-01 -2.4594E-01
             1.3688E-02

0ITERATION NO.:  128    OBJECTIVE VALUE:  -1758.43902386859        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     4777
 NPARAMETR:  1.0014E+00  1.0000E-02  1.8865E+00  1.6117E+00  9.6425E-01  9.8459E-01  7.3989E-01  1.5340E+00  5.5732E-01  1.1200E+00
             8.3247E-01
 PARAMETER:  1.0144E-01 -4.7326E+00  7.3942E-01  5.7749E-01  6.3109E-02  8.4456E-02 -1.9955E-01  5.3052E-01 -4.8470E-01  2.1547E-01
            -8.3362E-02
 GRADIENT:   7.4852E-03  0.0000E+00  6.0106E-01  4.7695E-01 -2.0782E-01 -2.5980E-03  5.7441E-04  1.1040E-01 -1.1597E-02  1.3859E-01
            -3.9337E-04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4777
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1416E-04 -3.0414E-04 -3.7837E-02 -9.5486E-03 -4.3518E-02
 SE:             2.9909E-02  1.6524E-04  1.8070E-02  2.9036E-02  2.1344E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9695E-01  6.5690E-02  3.6271E-02  7.4227E-01  4.1461E-02

 ETASHRINKSD(%)  1.0000E-10  9.9446E+01  3.9462E+01  2.7251E+00  2.8494E+01
 ETASHRINKVR(%)  1.0000E-10  9.9997E+01  6.3351E+01  5.3759E+00  4.8870E+01
 EBVSHRINKSD(%)  2.8886E-01  9.9475E+01  4.3909E+01  3.1099E+00  2.3519E+01
 EBVSHRINKVR(%)  5.7688E-01  9.9997E+01  6.8538E+01  6.1230E+00  4.1507E+01
 RELATIVEINF(%)  9.7668E+01  1.1426E-04  8.4191E+00  4.2536E+00  1.2148E+01
 EPSSHRINKSD(%)  4.5859E+01
 EPSSHRINKVR(%)  7.0687E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1758.4390238685905     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1023.2881973048524     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    65.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1758.439       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.90E+00  1.61E+00  9.64E-01  9.85E-01  7.41E-01  1.54E+00  5.57E-01  1.12E+00  8.32E-01
 


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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.50E+00  0.00E+00  4.49E+01
 
 TH 4
+       -1.38E+01  0.00E+00 -3.30E+01  1.30E+03
 
 TH 5
+       -3.88E-01  0.00E+00 -1.16E+02 -8.18E+01  5.79E+02
 
 TH 6
+       -2.50E-02  0.00E+00  5.15E-02 -2.57E+00  3.68E-01  2.04E+02
 
 TH 7
+       -1.44E-01  0.00E+00  1.29E-02  9.66E-03 -1.27E-01  2.18E-01  1.47E-01
 
 TH 8
+        6.03E-02  0.00E+00 -1.71E+01 -5.42E+00 -6.87E+00  8.28E-02 -2.61E-03  2.19E+01
 
 TH 9
+        2.07E+00  0.00E+00  5.38E+00 -1.29E+00  9.10E-01 -1.89E+00  2.26E-02  6.24E-01  5.72E+02
 
 TH10
+        7.69E-01  0.00E+00 -1.67E+00 -3.60E+00 -6.89E+01  4.74E-01 -1.64E-02  9.04E+00  1.63E+00  6.40E+01
 
 TH11
+       -7.22E+00  0.00E+00 -4.49E+00 -1.79E+01 -1.22E+00  2.17E+00  1.74E-01  6.10E+00  2.02E+01  8.43E+00  2.99E+02
 
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
 #CPUT: Total CPU Time in Seconds,       71.531
Stop Time:
Wed Sep 29 15:16:20 CDT 2021
