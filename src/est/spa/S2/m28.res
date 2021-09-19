Sat Sep 18 13:21:09 CDT 2021
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
$DATA ../../../../data/spa/S2/dat28.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1649.12880266761        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.4938E+01 -1.1226E+02 -5.3642E+01 -1.1149E+02  8.7850E+01  1.3603E+01  9.0735E-02  5.4051E+00 -1.5993E+01  1.8780E+01
             9.6861E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1661.61667154885        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7361E-01  1.0020E+00  1.0828E+00  1.0933E+00  9.4681E-01  9.7967E-01  8.9444E-01  9.9903E-01  1.0849E+00  7.9148E-01
             1.0220E+00
 PARAMETER:  7.3259E-02  1.0203E-01  1.7958E-01  1.8923E-01  4.5348E-02  7.9458E-02 -1.1558E-02  9.9026E-02  1.8148E-01 -1.3385E-01
             1.2171E-01
 GRADIENT:   1.3052E+01  8.0313E+00  2.0639E+00  2.0346E+01  1.0477E+01  6.3084E+00  3.1976E+00 -3.1417E+00  8.2772E+00 -7.4236E+00
             4.4865E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1662.15519731325        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.7119E-01  1.0481E+00  1.2355E+00  1.0653E+00  1.0246E+00  9.6503E-01  6.4642E-01  1.2311E+00  1.1529E+00  9.0177E-01
             1.0046E+00
 PARAMETER:  7.0770E-02  1.4694E-01  3.1150E-01  1.6330E-01  1.2429E-01  6.4401E-02 -3.3631E-01  3.0794E-01  2.4226E-01 -3.3921E-03
             1.0460E-01
 GRADIENT:   9.5167E+00  1.1264E+01  2.2938E+00  1.8075E+01  3.7932E+00  7.9957E-01  1.3611E+00 -1.3594E+00  6.6361E+00 -1.3385E+00
            -1.3256E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1662.24464363985        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.6824E-01  1.0079E+00  1.2098E+00  1.0814E+00  9.9796E-01  9.6330E-01  6.5860E-01  1.1798E+00  1.1130E+00  8.9507E-01
             1.0041E+00
 PARAMETER:  6.7727E-02  1.0792E-01  2.9047E-01  1.7826E-01  9.7959E-02  6.2615E-02 -3.1764E-01  2.6537E-01  2.0703E-01 -1.0850E-02
             1.0407E-01
 GRADIENT:   2.4890E+00  4.3587E+00  2.3448E+00  6.5932E+00 -1.0482E+00 -4.0425E-03  8.2422E-01 -1.0248E+00  2.3363E+00  4.9164E-01
            -1.1783E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.66681554256        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      363
 NPARAMETR:  9.8813E-01  9.2306E-01  1.2498E+00  1.1412E+00  9.7956E-01  9.7160E-01  6.4521E-01  1.2009E+00  1.0636E+00  8.8254E-01
             1.0073E+00
 PARAMETER:  8.8061E-02  1.9937E-02  3.2296E-01  2.3212E-01  7.9344E-02  7.1188E-02 -3.3818E-01  2.8310E-01  1.6164E-01 -2.4954E-02
             1.0723E-01
 GRADIENT:   1.4455E+01  4.1026E+00  6.3819E-01  2.4721E+00 -2.3831E+00  7.4703E-01  5.7531E-02  4.2242E-01 -8.2459E-01  6.8708E-02
            -2.5172E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1663.07542779767        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      539
 NPARAMETR:  9.8223E-01  6.9141E-01  1.3500E+00  1.2945E+00  9.3534E-01  9.6744E-01  4.9764E-01  1.1872E+00  9.6709E-01  8.7201E-01
             1.0073E+00
 PARAMETER:  8.2065E-02 -2.6903E-01  4.0009E-01  3.5811E-01  3.3160E-02  6.6897E-02 -5.9787E-01  2.7161E-01  6.6538E-02 -3.6952E-02
             1.0731E-01
 GRADIENT:   5.3029E+00  7.0320E+00  2.3991E+00  1.2046E+01 -5.0043E+00 -1.9663E-01 -3.4005E-01 -3.3337E-01 -1.8185E+00 -1.6429E-01
            -4.5792E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1663.33468136549        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      714
 NPARAMETR:  9.7560E-01  4.7719E-01  1.3795E+00  1.4268E+00  8.8123E-01  9.6616E-01  3.6546E-01  1.1677E+00  8.8610E-01  8.4649E-01
             1.0076E+00
 PARAMETER:  7.5294E-02 -6.3984E-01  4.2171E-01  4.5544E-01 -2.6432E-02  6.5576E-02 -9.0659E-01  2.5506E-01 -2.0921E-02 -6.6660E-02
             1.0760E-01
 GRADIENT:  -4.9555E+00  3.3034E+00  3.0337E+00  9.9837E+00 -5.5072E+00 -1.0713E-01 -8.7940E-02 -5.2647E-01  3.5843E-01  2.2101E-01
            -2.0606E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1663.43217645020        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  9.7582E-01  3.6046E-01  1.3727E+00  1.4931E+00  8.4914E-01  9.6551E-01  3.2912E-01  1.1681E+00  8.4036E-01  8.2339E-01
             1.0081E+00
 PARAMETER:  7.5524E-02 -9.2038E-01  4.1678E-01  5.0084E-01 -6.3527E-02  6.4902E-02 -1.0113E+00  2.5539E-01 -7.3923E-02 -9.4322E-02
             1.0811E-01
 GRADIENT:  -6.2964E-01  1.1844E-01  1.5259E-01 -6.3275E-02 -3.0432E-01 -8.4637E-03 -3.3227E-02 -6.2041E-03  1.4113E-01  3.4448E-02
            -5.1930E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1663.44051372753        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  9.7519E-01  2.9888E-01  1.3742E+00  1.5311E+00  8.3325E-01  9.6509E-01  3.1943E-01  1.1751E+00  8.1713E-01  8.1132E-01
             1.0083E+00
 PARAMETER:  7.4882E-02 -1.1077E+00  4.1785E-01  5.2598E-01 -8.2425E-02  6.4471E-02 -1.0412E+00  2.6135E-01 -1.0196E-01 -1.0909E-01
             1.0828E-01
 GRADIENT:   3.9271E-02  1.6876E-01  9.3606E-02  8.8516E-01 -2.7948E-01  8.1626E-03 -1.9835E-02 -1.2850E-02 -1.4037E-01 -7.1375E-03
            -2.0789E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1663.46382573936        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  9.7483E-01  2.8850E-01  1.3755E+00  1.5382E+00  8.3049E-01  9.6495E-01  8.4016E-01  1.1767E+00  8.0665E-01  8.0531E-01
             1.0080E+00
 PARAMETER:  7.4504E-02 -1.1431E+00  4.1880E-01  5.3061E-01 -8.5743E-02  6.4325E-02 -7.4164E-02  2.6275E-01 -1.1487E-01 -1.1652E-01
             1.0796E-01
 GRADIENT:  -4.4162E-01  4.0209E-02  5.1137E-01  4.1643E-01  1.5006E-01 -4.5814E-02 -6.7517E-03 -8.0480E-02 -3.3542E-01 -5.0038E-02
            -1.5012E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1663.46691633854        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1427
 NPARAMETR:  9.7576E-01  2.9110E-01  1.3626E+00  1.5369E+00  8.2689E-01  9.6549E-01  8.8918E-01  1.1655E+00  8.0807E-01  8.0269E-01
             1.0084E+00
 PARAMETER:  7.5461E-02 -1.1341E+00  4.0942E-01  5.2978E-01 -9.0087E-02  6.4882E-02 -1.7450E-02  2.5318E-01 -1.1310E-01 -1.1979E-01
             1.0832E-01
 GRADIENT:   1.6571E+00  1.9632E-01  1.9772E-01  2.4318E+00  1.4517E-01  1.4616E-01  3.1360E-03 -5.7624E-02  1.8457E-01 -2.0993E-02
             5.2600E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1663.46955362613        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  9.7518E-01  2.8534E-01  1.3470E+00  1.5386E+00  8.1999E-01  9.6519E-01  9.3951E-01  1.1534E+00  8.0551E-01  7.9754E-01
             1.0084E+00
 PARAMETER:  7.4863E-02 -1.1541E+00  3.9789E-01  5.3088E-01 -9.8462E-02  6.4571E-02  3.7602E-02  2.4269E-01 -1.1628E-01 -1.2622E-01
             1.0832E-01
 GRADIENT:   4.0942E-01  1.6857E-02 -2.2223E-03  7.1121E-01  6.2517E-02  3.2429E-02 -6.7755E-05 -1.0584E-02  4.6244E-02 -3.4441E-02
            -1.0772E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1663.47478456742        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1786
 NPARAMETR:  9.7559E-01  3.2081E-01  1.3403E+00  1.5174E+00  8.2654E-01  9.6539E-01  9.0685E-01  1.1434E+00  8.1651E-01  8.0170E-01
             1.0083E+00
 PARAMETER:  7.5290E-02 -1.0369E+00  3.9286E-01  5.1698E-01 -9.0509E-02  6.4782E-02  2.2209E-03  2.3397E-01 -1.0271E-01 -1.2102E-01
             1.0822E-01
 GRADIENT:   5.7015E-02  2.7042E-01  2.7586E-01  1.6590E+00 -5.6089E-01 -6.7051E-03 -6.5496E-03 -3.4489E-02 -1.3605E-01 -6.4938E-02
            -2.7594E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1663.47611151864        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1970             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7579E-01  3.3536E-01  1.3396E+00  1.5077E+00  8.3042E-01  9.6553E-01  8.7868E-01  1.1413E+00  8.2193E-01  8.0505E-01
             1.0082E+00
 PARAMETER:  7.5494E-02 -9.9256E-01  3.9235E-01  5.1061E-01 -8.5823E-02  6.4919E-02 -2.9330E-02  2.3220E-01 -9.6097E-02 -1.1685E-01
             1.0816E-01
 GRADIENT:   3.6241E+01  3.5109E+00  8.9048E-01  5.9537E+01  7.1462E-01  3.0134E+00  2.0181E-02  8.1167E-02  8.3156E-01  7.4070E-02
             8.4767E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1663.47611487311        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2145
 NPARAMETR:  9.7577E-01  3.3531E-01  1.3395E+00  1.5078E+00  8.3040E-01  9.6551E-01  8.8277E-01  1.1413E+00  8.2182E-01  8.0495E-01
             1.0082E+00
 PARAMETER:  7.5468E-02 -9.9269E-01  3.9230E-01  5.1064E-01 -8.5849E-02  6.4904E-02 -2.4695E-02  2.3220E-01 -9.6228E-02 -1.1697E-01
             1.0816E-01
 GRADIENT:  -3.3821E-02  3.4902E-04  5.0194E-03 -5.1799E-02  7.0234E-03 -3.3088E-03  2.3118E-04  2.1825E-03  7.9640E-04  2.0068E-03
             9.4674E-04

0ITERATION NO.:   72    OBJECTIVE VALUE:  -1663.47611508957        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     2204
 NPARAMETR:  9.7577E-01  3.3531E-01  1.3395E+00  1.5078E+00  8.3039E-01  9.6551E-01  8.8286E-01  1.1413E+00  8.2182E-01  8.0494E-01
             1.0082E+00
 PARAMETER:  7.5468E-02 -9.9271E-01  3.9229E-01  5.1064E-01 -8.5854E-02  6.4904E-02 -2.4584E-02  2.3219E-01 -9.6236E-02 -1.1698E-01
             1.0816E-01
 GRADIENT:  -4.5762E-02 -1.4259E-03  6.2520E-03  1.8709E-02  6.3300E-03 -2.5035E-03 -5.5961E-04  1.6671E-03 -1.3596E-03  1.1756E-03
             3.4044E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2204
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0552E-04 -7.7470E-03 -2.7379E-02 -2.9952E-03 -3.2194E-02
 SE:             2.9846E-02  5.1814E-03  1.8236E-02  2.8911E-02  2.0330E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9183E-01  1.3487E-01  1.3325E-01  9.1749E-01  1.1329E-01

 ETASHRINKSD(%)  1.2269E-02  8.2642E+01  3.8908E+01  3.1430E+00  3.1892E+01
 ETASHRINKVR(%)  2.4536E-02  9.6987E+01  6.2678E+01  6.1872E+00  5.3613E+01
 EBVSHRINKSD(%)  4.3672E-01  8.3435E+01  4.1068E+01  3.4536E+00  3.0249E+01
 EBVSHRINKVR(%)  8.7153E-01  9.7256E+01  6.5271E+01  6.7879E+00  5.1348E+01
 RELATIVEINF(%)  9.7260E+01  1.2510E-01  6.5618E+00  5.7742E+00  4.8287E+00
 EPSSHRINKSD(%)  4.4806E+01
 EPSSHRINKVR(%)  6.9536E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1663.4761150895688     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -928.32528852583062     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.97
 Elapsed covariance  time in seconds:     5.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1663.476       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  3.35E-01  1.34E+00  1.51E+00  8.30E-01  9.66E-01  8.83E-01  1.14E+00  8.22E-01  8.05E-01  1.01E+00
 


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
 
         2.95E-02  4.46E-01  4.45E-01  2.54E-01  2.55E-01  5.84E-02  5.91E-01  3.93E-01  1.56E-01  2.28E-01  5.98E-02
 


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
+        8.70E-04
 
 TH 2
+        3.74E-03  1.99E-01
 
 TH 3
+        2.75E-03  1.40E-01  1.98E-01
 
 TH 4
+       -1.92E-03 -1.12E-01 -7.41E-02  6.47E-02
 
 TH 5
+        1.86E-03  1.01E-01  1.06E-01 -5.52E-02  6.48E-02
 
 TH 6
+       -3.56E-04 -3.64E-03 -4.50E-03  1.83E-03 -2.79E-03  3.41E-03
 
 TH 7
+       -4.80E-03 -2.40E-01 -2.38E-01  1.32E-01 -1.48E-01  8.84E-03  3.49E-01
 
 TH 8
+        2.01E-03  1.17E-01  1.59E-01 -6.20E-02  8.70E-02 -5.03E-03 -2.01E-01  1.55E-01
 
 TH 9
+        1.33E-03  6.44E-02  4.70E-02 -3.62E-02  3.37E-02 -1.83E-03 -8.10E-02  3.93E-02  2.43E-02
 
 TH10
+        1.94E-03  7.80E-02  8.13E-02 -4.28E-02  4.96E-02 -9.63E-04 -1.11E-01  6.22E-02  2.58E-02  5.21E-02
 
 TH11
+       -1.78E-04 -2.02E-03 -3.12E-03  9.38E-04 -1.40E-03 -2.93E-04  1.83E-03 -1.45E-03 -2.77E-04 -2.55E-03  3.57E-03
 
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
+        2.95E-02
 
 TH 2
+        2.85E-01  4.46E-01
 
 TH 3
+        2.09E-01  7.04E-01  4.45E-01
 
 TH 4
+       -2.56E-01 -9.87E-01 -6.55E-01  2.54E-01
 
 TH 5
+        2.48E-01  8.91E-01  9.40E-01 -8.53E-01  2.55E-01
 
 TH 6
+       -2.07E-01 -1.40E-01 -1.73E-01  1.23E-01 -1.88E-01  5.84E-02
 
 TH 7
+       -2.76E-01 -9.13E-01 -9.07E-01  8.78E-01 -9.86E-01  2.56E-01  5.91E-01
 
 TH 8
+        1.73E-01  6.66E-01  9.10E-01 -6.20E-01  8.69E-01 -2.19E-01 -8.67E-01  3.93E-01
 
 TH 9
+        2.90E-01  9.27E-01  6.78E-01 -9.14E-01  8.49E-01 -2.02E-01 -8.80E-01  6.41E-01  1.56E-01
 
 TH10
+        2.88E-01  7.67E-01  8.01E-01 -7.38E-01  8.53E-01 -7.23E-02 -8.23E-01  6.92E-01  7.25E-01  2.28E-01
 
 TH11
+       -1.01E-01 -7.57E-02 -1.17E-01  6.17E-02 -9.17E-02 -8.40E-02  5.19E-02 -6.19E-02 -2.97E-02 -1.87E-01  5.98E-02
 
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
+        1.43E+03
 
 TH 2
+       -1.54E+02  4.28E+02
 
 TH 3
+       -6.70E+01  1.23E+02  1.99E+02
 
 TH 4
+       -1.81E+02  4.58E+02  1.46E+01  7.22E+02
 
 TH 5
+        3.37E+02 -3.03E+02 -4.34E+02 -1.23E+02  1.87E+03
 
 TH 6
+        1.09E+02 -8.48E+01 -5.17E+01  2.12E+00 -1.05E+02  4.52E+02
 
 TH 7
+        6.46E+01  7.48E+01  2.22E+01  4.36E+00  3.44E+02 -1.63E+02  2.58E+02
 
 TH 8
+        2.41E+01  8.23E+00 -2.51E+01 -5.61E+00  4.91E+01 -1.57E+01  4.64E+01  4.89E+01
 
 TH 9
+       -2.95E+01 -3.82E+01  3.29E+01  2.34E+01  3.79E+00  1.85E-01  7.44E+01  1.57E+01  3.39E+02
 
 TH10
+       -6.02E+01  7.65E+00 -5.01E+00  1.45E+01 -9.60E+01 -1.66E+01 -1.62E+01  7.31E+00 -6.12E+00  8.18E+01
 
 TH11
+        4.58E+01  6.30E+01  3.92E+01  3.22E+01 -1.59E+00 -2.66E+01  6.57E+01  5.85E+00 -8.71E+00  2.34E+01  3.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.392
Stop Time:
Sat Sep 18 13:21:41 CDT 2021
