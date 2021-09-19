Sat Sep 18 07:17:15 CDT 2021
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
$DATA ../../../../data/int/D/dat64.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22905.2986983775        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.7823E+02  3.2666E+02 -4.4184E+01  2.4957E+02  5.2078E+02 -3.5733E+03 -1.3458E+03 -9.8926E+01 -1.9903E+03 -1.3834E+03
            -4.5544E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1113.71216195643        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  4.5429E+00  1.8937E+00  1.1315E+00  2.1966E+00  4.3646E-01  7.2789E+00  3.8488E+00  9.8717E-01  4.1051E+00  3.4678E+00
             8.7014E+00
 PARAMETER:  1.6136E+00  7.3853E-01  2.2350E-01  8.8690E-01 -7.2906E-01  2.0850E+00  1.4478E+00  8.7090E-02  1.5122E+00  1.3435E+00
             2.2635E+00
 GRADIENT:   6.0719E+01  8.1246E+01  2.6172E+00  6.8792E+01 -9.2413E+01  1.2097E+02  1.6033E+01  4.2488E+00  3.6904E+01  3.0940E+00
             2.3915E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1274.62227417001        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  3.8284E+00  1.1570E+00  8.8129E+01  2.7435E+00  1.7650E+00  6.0839E+00  9.8541E+00  2.2362E-01  2.9773E+00  2.4849E+00
             8.9692E+00
 PARAMETER:  1.4425E+00  2.4579E-01  4.5788E+00  1.1092E+00  6.6813E-01  1.9056E+00  2.3879E+00 -1.3978E+00  1.1910E+00  1.0102E+00
             2.2938E+00
 GRADIENT:   6.9324E+01  1.9186E+01  2.3661E+00  7.9598E+01 -6.0369E+01  9.1291E+01  2.8032E+01  3.9247E-04 -6.5873E+00  8.5027E+01
             3.3259E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1402.86127527585        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.8119E+00  1.4136E+00  3.4596E+00  9.8433E-01  1.6184E+00  4.2439E+00  6.3095E+00  2.0736E+00  1.9921E+00  8.9701E-01
             9.0575E+00
 PARAMETER:  6.9438E-01  4.4615E-01  1.3412E+00  8.4204E-02  5.8146E-01  1.5455E+00  1.9421E+00  8.2931E-01  7.8920E-01 -8.6913E-03
             2.3036E+00
 GRADIENT:   3.3484E+01 -3.8009E+00 -1.7841E+01 -1.5270E+01 -1.2211E+01  5.2854E+01  4.0887E+01  6.0833E+00  2.9076E+01  1.8124E+01
             4.0593E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1479.20036283813        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.2549E+00  6.2609E-01  2.1380E+01  1.2082E+00  1.9871E+00  3.7850E+00  7.1217E+00  4.1782E+00  1.9022E+00  6.8725E-01
             7.0325E+00
 PARAMETER:  3.2703E-01 -3.6826E-01  3.1624E+00  2.8914E-01  7.8666E-01  1.4311E+00  2.0631E+00  1.5299E+00  7.4301E-01 -2.7505E-01
             2.0505E+00
 GRADIENT:  -6.2196E+00 -1.2790E+01  2.0246E+00 -1.4391E+01  2.8297E+00  1.7173E+01  1.4396E+01 -1.5350E+00  1.7482E+01  8.1163E+00
             2.3584E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1481.03845479432        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.3229E+00  8.5397E-01  7.0555E+00  1.0800E+00  1.8223E+00  3.5892E+00  6.3788E+00  2.4852E+00  1.5128E+00  6.3684E-01
             7.0111E+00
 PARAMETER:  3.7982E-01 -5.7861E-02  2.0538E+00  1.7692E-01  7.0008E-01  1.3779E+00  1.9530E+00  1.0103E+00  5.1395E-01 -3.5123E-01
             2.0475E+00
 GRADIENT:   2.1101E+00 -1.2077E+01 -5.3172E+00 -3.9943E+00  6.7147E+00 -3.7031E+00  4.6157E+00  1.2983E+00  5.0167E+00  7.3018E+00
             4.1043E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1481.33485309170        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.3037E+00  1.0983E+00  6.8022E+00  1.0269E+00  1.8502E+00  3.5838E+00  6.0417E+00  2.3634E+00  1.4100E+00  6.6085E-01
             7.0216E+00
 PARAMETER:  3.6520E-01  1.9376E-01  2.0173E+00  1.2654E-01  7.1528E-01  1.3764E+00  1.8987E+00  9.6008E-01  4.4362E-01 -3.1423E-01
             2.0490E+00
 GRADIENT:  -7.2681E-01 -5.5137E+00 -3.7976E+00  1.9421E+00  7.0494E+00 -4.5308E+00  1.5245E+00  4.2733E-01  3.9345E-01  7.7490E+00
             1.0293E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1483.49556495159        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  1.2932E+00  1.3780E+00  1.0976E+01  9.5566E-01  1.9674E+00  3.6494E+00  5.7688E+00  3.4541E+00  1.3033E+00  3.6127E-01
             7.0573E+00
 PARAMETER:  3.5709E-01  4.2063E-01  2.4957E+00  5.4651E-02  7.7670E-01  1.3946E+00  1.8525E+00  1.3396E+00  3.6491E-01 -9.1814E-01
             2.0541E+00
 GRADIENT:  -1.9149E+00  1.6944E+00  1.8350E+00 -2.0607E+00 -4.7182E-01  2.8946E+00 -6.4118E-01 -6.1144E-01 -3.5162E+00  1.7937E+00
             4.4801E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1487.55300067773        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      642
 NPARAMETR:  1.3742E+00  1.0984E+00  1.1586E+01  1.0920E+00  1.9362E+00  3.7249E+00  6.8883E+00  3.5097E+00  1.4978E+00  1.9983E-01
             7.0458E+00
 PARAMETER:  4.1785E-01  1.9382E-01  2.5498E+00  1.8798E-01  7.6073E-01  1.4150E+00  2.0298E+00  1.3555E+00  5.0401E-01 -1.5103E+00
             2.0524E+00
 GRADIENT:   4.0725E+00  1.9654E+00 -2.0229E-01 -4.8455E+00  1.6466E+00 -7.4120E+00 -1.4248E+00  1.2621E+00  1.9539E+00  5.4324E-01
            -1.8752E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1488.38816615424        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      817
 NPARAMETR:  1.3410E+00  8.4704E-01  1.2195E+01  1.1786E+00  1.9229E+00  3.7934E+00  7.5745E+00  3.3863E+00  1.4571E+00  1.0653E-01
             7.0543E+00
 PARAMETER:  3.9342E-01 -6.6004E-02  2.6010E+00  2.6430E-01  7.5384E-01  1.4333E+00  2.1248E+00  1.3197E+00  4.7648E-01 -2.1393E+00
             2.0536E+00
 GRADIENT:   2.3475E-01 -2.2121E-01 -1.8659E-01 -8.6044E-02  8.3620E-01 -1.5041E-01  4.1650E-01  1.7804E-01 -1.3066E-01  1.4373E-01
            -1.7903E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1488.46565770472        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      992
 NPARAMETR:  1.3389E+00  8.6425E-01  1.2118E+01  1.1739E+00  1.9206E+00  3.7941E+00  7.5164E+00  3.3686E+00  1.4624E+00  1.2332E-02
             7.0611E+00
 PARAMETER:  3.9186E-01 -4.5888E-02  2.5947E+00  2.6033E-01  7.5265E-01  1.4335E+00  2.1171E+00  1.3145E+00  4.8009E-01 -4.2956E+00
             2.0546E+00
 GRADIENT:  -2.2705E-02 -3.0396E-02  7.3334E-02 -2.1746E-01 -2.8688E-01 -3.0151E-02  3.7876E-03 -4.3400E-02  8.9204E-02  1.9079E-03
             2.2089E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1488.46654136260        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1180             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3392E+00  8.6609E-01  1.2081E+01  1.1732E+00  1.9210E+00  3.7973E+00  7.5065E+00  3.3669E+00  1.4617E+00  1.0000E-02
             7.0609E+00
 PARAMETER:  3.9205E-01 -4.3771E-02  2.5916E+00  2.5972E-01  7.5286E-01  1.4343E+00  2.1158E+00  1.3140E+00  4.7962E-01 -4.7269E+00
             2.0546E+00
 GRADIENT:   4.3524E+00  1.9291E-01  8.2082E-02  1.1654E+00  5.2250E-01  1.9379E+01  3.6040E+01  4.0276E-02  3.4028E-01  0.0000E+00
             3.8426E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1488.46659776817        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  1.3392E+00  8.6639E-01  1.2076E+01  1.1733E+00  1.9210E+00  3.7966E+00  7.5052E+00  3.3672E+00  1.4608E+00  1.0000E-02
             7.0607E+00
 PARAMETER:  3.9204E-01 -4.3420E-02  2.5913E+00  2.5985E-01  7.5287E-01  1.4341E+00  2.1156E+00  1.3141E+00  4.7896E-01 -4.7269E+00
             2.0545E+00
 GRADIENT:   3.8657E-03 -2.7096E-02  1.5121E-03  6.8147E-04  1.2503E-02  2.1414E-01 -1.7437E-01  5.1679E-03  2.1194E-03  0.0000E+00
            -1.2486E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1488.46667215009        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1551            RESET HESSIAN, TYPE II
 NPARAMETR:  1.3392E+00  8.6759E-01  1.2069E+01  1.1730E+00  1.9211E+00  3.7974E+00  7.5011E+00  3.3667E+00  1.4610E+00  1.0000E-02
             7.0608E+00
 PARAMETER:  3.9206E-01 -4.2034E-02  2.5906E+00  2.5955E-01  7.5289E-01  1.4343E+00  2.1151E+00  1.3139E+00  4.7915E-01 -4.7269E+00
             2.0546E+00
 GRADIENT:   4.3534E+00  2.1674E-01  6.8657E-02  1.2963E+00  5.5107E-01  1.9382E+01  3.5915E+01  5.4248E-02  2.9209E-01  0.0000E+00
             3.7311E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1488.46669978831        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1680             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3392E+00  8.6818E-01  1.2067E+01  1.1727E+00  1.9210E+00  3.7974E+00  7.4984E+00  3.3664E+00  1.4610E+00  1.0000E-02
             7.0607E+00
 PARAMETER:  3.9205E-01 -4.1360E-02  2.5904E+00  2.5932E-01  7.5287E-01  1.4343E+00  2.1147E+00  1.3138E+00  4.7913E-01 -4.7269E+00
             2.0545E+00
 GRADIENT:   4.3535E+00  2.1715E-01  7.5967E-02  1.2812E+00  5.2451E-01  1.9382E+01  3.5868E+01  4.9775E-02  2.9199E-01  0.0000E+00
             3.7245E+00

0ITERATION NO.:   72    OBJECTIVE VALUE:  -1488.46669978831        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     1741
 NPARAMETR:  1.3391E+00  8.6825E-01  1.2061E+01  1.1728E+00  1.9210E+00  3.7974E+00  7.4962E+00  3.3658E+00  1.4603E+00  1.0000E-02
             7.0605E+00
 PARAMETER:  3.9205E-01 -4.1360E-02  2.5904E+00  2.5932E-01  7.5287E-01  1.4343E+00  2.1147E+00  1.3138E+00  4.7913E-01 -4.7269E+00
             2.0545E+00
 GRADIENT:   2.0792E-03 -6.3083E-03  5.4687E-03 -1.2059E-02  7.8834E-03 -5.4992E-05  2.5667E-02  2.5314E-03  1.5375E-02  0.0000E+00
             2.8773E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1741
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.9401E-03  4.0591E-02 -2.5185E-02 -7.1789E-02 -3.5381E-05
 SE:             2.9595E-02  2.3807E-02  9.9081E-03  1.5890E-02  1.6357E-04
 N:                     100         100         100         100         100

 P VAL.:         8.6743E-01  8.8190E-02  1.1025E-02  6.2545E-06  8.2875E-01

 ETASHRINKSD(%)  8.5240E-01  2.0245E+01  6.6807E+01  4.6765E+01  9.9452E+01
 ETASHRINKVR(%)  1.6975E+00  3.6391E+01  8.8982E+01  7.1661E+01  9.9997E+01
 EBVSHRINKSD(%)  7.0349E-01  1.3737E+01  7.4215E+01  5.0637E+01  9.9385E+01
 EBVSHRINKVR(%)  1.4020E+00  2.5587E+01  9.3351E+01  7.5633E+01  9.9996E+01
 RELATIVEINF(%)  9.8555E+01  4.1026E+01  2.9761E+00  1.3433E+01  1.7730E-03
 EPSSHRINKSD(%)  8.7961E+00
 EPSSHRINKVR(%)  1.6818E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1488.4666997883076     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       165.62265998010321     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    59.22
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1488.467       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.34E+00  8.68E-01  1.21E+01  1.17E+00  1.92E+00  3.80E+00  7.50E+00  3.37E+00  1.46E+00  1.00E-02  7.06E+00
 


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
+        4.42E+01
 
 TH 2
+        1.39E+00  6.43E+01
 
 TH 3
+        8.38E-03  1.13E-01  9.50E-02
 
 TH 4
+       -3.63E+00  1.84E+01 -4.67E-01  1.43E+02
 
 TH 5
+       -6.07E-01 -6.42E+00 -1.85E+00  1.38E+00  9.95E+01
 
 TH 6
+        3.46E-01  1.77E-01  3.58E-03  1.17E-01 -6.29E-01  1.32E+01
 
 TH 7
+        5.24E-02  2.60E+00 -5.19E-02 -7.05E+00  1.29E+00 -2.92E-02  1.97E+00
 
 TH 8
+       -4.45E-02  3.27E-01 -2.66E-01  6.05E-01  1.11E+00  2.18E-03  1.16E-01  1.69E+00
 
 TH 9
+       -1.65E-01 -2.41E+00 -6.23E-02 -1.90E+01  5.85E+00 -1.39E-01  1.10E+00 -3.55E-01  1.72E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.95E+00 -1.30E+00 -1.78E-02 -8.73E+00 -1.14E+00  6.09E-01  3.88E-01  1.25E-01  3.84E+00  0.00E+00  2.11E+01
 
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
 #CPUT: Total CPU Time in Seconds,       78.196
Stop Time:
Sat Sep 18 07:18:35 CDT 2021
