Sat Sep 25 12:27:55 CDT 2021
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
$DATA ../../../../data/spa/S2/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1676.90669777208        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0630E+01 -4.8589E+01 -1.5505E+01 -5.9155E+01  9.7814E+00  2.4430E+01 -6.4542E+00  8.9848E+00 -7.2491E+00  1.4624E+01
             4.1263E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1680.46415773621        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.9137E-01  1.0327E+00  1.0432E+00  1.0259E+00  1.0172E+00  9.4683E-01  1.0424E+00  9.4630E-01  1.0322E+00  9.2622E-01
             1.0028E+00
 PARAMETER:  9.1330E-02  1.3213E-01  1.4233E-01  1.2556E-01  1.1701E-01  4.5367E-02  1.4153E-01  4.4804E-02  1.3172E-01  2.3360E-02
             1.0276E-01
 GRADIENT:  -1.1607E+01 -1.3041E-01  2.0254E+00 -1.1239E+00  4.3476E+00  3.3163E-01 -2.0538E+00  2.7439E+00  1.3388E+00 -4.3364E-01
             1.2937E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1681.74663182455        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      332
 NPARAMETR:  9.9337E-01  8.5484E-01  8.4522E-01  1.1386E+00  8.3543E-01  9.5084E-01  1.2974E+00  5.8772E-01  9.4085E-01  7.6440E-01
             9.7118E-01
 PARAMETER:  9.3346E-02 -5.6844E-02 -6.8161E-02  2.2983E-01 -7.9805E-02  4.9589E-02  3.6035E-01 -4.3150E-01  3.9026E-02 -1.6867E-01
             7.0760E-02
 GRADIENT:  -8.0037E+00  9.5188E+00 -1.3566E+01  3.4779E+01  1.7023E+01  1.4434E+00  5.0347E-01  7.1505E-01  2.6437E+00 -1.7216E+00
            -1.1114E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1682.44269957309        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  9.9613E-01  7.5775E-01  8.6260E-01  1.1828E+00  8.0293E-01  9.4644E-01  1.4204E+00  5.2840E-01  8.9792E-01  7.6740E-01
             9.9686E-01
 PARAMETER:  9.6120E-02 -1.7741E-01 -4.7806E-02  2.6789E-01 -1.1948E-01  4.4950E-02  4.5094E-01 -5.3789E-01 -7.6699E-03 -1.6474E-01
             9.6854E-02
 GRADIENT:   7.2326E-01  5.3453E+00 -5.7103E-02  8.5954E+00 -6.4956E-01  1.8002E-01  9.6492E-01 -5.6970E-01  1.5585E-01 -4.3861E-01
            -6.1137E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1682.86965410157        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      683
 NPARAMETR:  9.9223E-01  5.1754E-01  9.6186E-01  1.3310E+00  7.7276E-01  9.4334E-01  1.8220E+00  6.3247E-01  8.3707E-01  7.8166E-01
             9.9580E-01
 PARAMETER:  9.2198E-02 -5.5867E-01  6.1113E-02  3.8593E-01 -1.5779E-01  4.1671E-02  6.9996E-01 -3.5812E-01 -7.7852E-02 -1.4634E-01
             9.5793E-02
 GRADIENT:   4.3062E-01  3.8036E+00  4.5630E+00  8.6226E+00 -5.1205E+00  2.5427E-01  2.2660E-02 -7.1798E-01  6.1446E-02 -6.6256E-01
            -1.1533E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1683.10599551095        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      858
 NPARAMETR:  9.8871E-01  3.3465E-01  9.9811E-01  1.4371E+00  7.3972E-01  9.3961E-01  2.3850E+00  7.0426E-01  7.9802E-01  7.8862E-01
             9.9717E-01
 PARAMETER:  8.8642E-02 -9.9467E-01  9.8111E-02  4.6266E-01 -2.0149E-01  3.7713E-02  9.6920E-01 -2.5061E-01 -1.2562E-01 -1.3747E-01
             9.7170E-02
 GRADIENT:   3.1366E-01  2.3115E+00  1.5680E+00  7.6662E+00 -5.7804E+00 -8.2393E-02  3.5151E-01  7.0267E-01 -6.1711E-01  1.1303E+00
             3.5319E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1683.28161569077        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1036
 NPARAMETR:  9.8565E-01  2.1181E-01  1.0575E+00  1.5106E+00  7.3937E-01  9.3687E-01  3.0280E+00  7.7291E-01  7.7853E-01  7.9991E-01
             9.9704E-01
 PARAMETER:  8.5547E-02 -1.4521E+00  1.5591E-01  5.1250E-01 -2.0196E-01  3.4791E-02  1.2079E+00 -1.5760E-01 -1.5035E-01 -1.2326E-01
             9.7034E-02
 GRADIENT:  -7.7137E-02  6.2809E-01 -1.9232E-01  5.1742E+00 -4.2961E-02 -2.2751E-01 -2.2180E-01  3.9861E-02 -6.0708E-02 -6.1782E-02
             1.0413E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1683.40721611809        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1212
 NPARAMETR:  9.8306E-01  1.0700E-01  1.1142E+00  1.5760E+00  7.3826E-01  9.3519E-01  4.2618E+00  8.5122E-01  7.6097E-01  8.0542E-01
             9.9661E-01
 PARAMETER:  8.2919E-02 -2.1349E+00  2.0816E-01  5.5490E-01 -2.0345E-01  3.2995E-02  1.5497E+00 -6.1079E-02 -1.7316E-01 -1.1640E-01
             9.6609E-02
 GRADIENT:  -3.5222E-01  5.1223E-01  2.1248E+00  8.9144E+00 -3.6281E+00 -4.1892E-02 -1.4832E-01 -4.2125E-02 -1.7054E-01  1.2728E-02
            -1.2166E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1683.61478089354        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1391
 NPARAMETR:  9.8088E-01  3.7978E-02  1.0338E+00  1.6079E+00  6.8438E-01  9.3372E-01  6.9146E+00  7.9475E-01  7.5361E-01  7.7250E-01
             9.9613E-01
 PARAMETER:  8.0700E-02 -3.1708E+00  1.3329E-01  5.7496E-01 -2.7925E-01  3.1423E-02  2.0336E+00 -1.2973E-01 -1.8288E-01 -1.5812E-01
             9.6125E-02
 GRADIENT:  -2.0671E+00  4.8229E+00  2.9223E+00  1.8932E+01 -1.3578E+01 -6.9122E-01  8.3085E+00  1.4900E-02 -6.5284E+00  9.6928E-01
            -1.4980E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1684.29017396612        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1572
 NPARAMETR:  9.8108E-01  1.3645E-02  1.0244E+00  1.6024E+00  6.8163E-01  9.3393E-01  1.0743E+01  7.8703E-01  7.4967E-01  7.4284E-01
             9.9735E-01
 PARAMETER:  8.0903E-02 -4.1944E+00  1.2410E-01  5.7148E-01 -2.8328E-01  3.1648E-02  2.4743E+00 -1.3948E-01 -1.8812E-01 -1.9727E-01
             9.7342E-02
 GRADIENT:  -1.7733E-01 -2.0621E+00  2.2492E+00  5.5897E+00 -3.3241E+00  1.3429E-01 -4.6147E+00  4.9496E-01  2.8884E+00  9.3115E-01
             2.1917E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1684.32939420221        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1753
 NPARAMETR:  9.8121E-01  1.6410E-02  1.0286E+00  1.6008E+00  6.8385E-01  9.3402E-01  1.0046E+01  7.9634E-01  7.5052E-01  7.4541E-01
             9.9769E-01
 PARAMETER:  8.1030E-02 -4.0098E+00  1.2821E-01  5.7052E-01 -2.8002E-01  3.1739E-02  2.4072E+00 -1.2773E-01 -1.8698E-01 -1.9382E-01
             9.7684E-02
 GRADIENT:   2.8511E-02 -1.4874E+00  1.5621E+00  1.8253E+00 -2.8668E+00  1.1916E-01 -3.6228E+00  7.8797E-01  2.4110E+00  1.0994E+00
             5.4770E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1684.39611077918        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1930
 NPARAMETR:  9.8110E-01  1.2583E-02  1.0145E+00  1.6006E+00  6.7750E-01  9.3387E-01  1.1381E+01  7.7497E-01  7.4906E-01  7.3893E-01
             9.9753E-01
 PARAMETER:  8.0921E-02 -4.2754E+00  1.1442E-01  5.7037E-01 -2.8935E-01  3.1584E-02  2.5319E+00 -1.5493E-01 -1.8893E-01 -2.0255E-01
             9.7530E-02
 GRADIENT:   1.3100E-01  2.0454E+00 -2.0185E+00 -5.4243E+00  3.2203E+00 -1.6943E-01  4.2091E+00 -4.8331E-01 -2.6524E+00 -9.5216E-01
            -1.5147E-01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1684.39715404395        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1993
 NPARAMETR:  9.8109E-01  1.2475E-02  1.0148E+00  1.6010E+00  6.7732E-01  9.3390E-01  1.1431E+01  7.7677E-01  7.4868E-01  7.3999E-01
             9.9738E-01
 PARAMETER:  8.0903E-02 -4.2837E+00  1.1470E-01  5.7056E-01 -2.8964E-01  3.1618E-02  2.5365E+00 -1.5246E-01 -1.8944E-01 -2.0095E-01
             9.7406E-02
 GRADIENT:  -2.7939E-02  2.2994E+02  4.2974E-01 -1.7338E+03 -3.8657E-01  1.9237E-02  3.8792E+02  1.2503E-01  1.3712E-01  2.2266E-01
             6.4688E-02
 NUMSIGDIG:         4.0         3.3         2.5         3.3         3.3         3.4         3.3         2.1         3.6         2.2
                    2.6

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1993
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2004E-04  1.5096E-02 -1.8682E-02 -8.3078E-03 -1.9960E-02
 SE:             2.9818E-02  7.3457E-03  1.6363E-02  2.8796E-02  2.1069E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9411E-01  3.9869E-02  2.5355E-01  7.7296E-01  3.4346E-01

 ETASHRINKSD(%)  1.0686E-01  7.5391E+01  4.5183E+01  3.5298E+00  2.9416E+01
 ETASHRINKVR(%)  2.1362E-01  9.3944E+01  6.9951E+01  6.9350E+00  5.0179E+01
 EBVSHRINKSD(%)  4.6396E-01  8.0811E+01  4.6145E+01  3.5236E+00  2.7335E+01
 EBVSHRINKVR(%)  9.2576E-01  9.6318E+01  7.0997E+01  6.9230E+00  4.7198E+01
 RELATIVEINF(%)  9.9002E+01  2.9647E+00  3.1617E+00  6.6947E+01  5.7382E+00
 EPSSHRINKSD(%)  4.4233E+01
 EPSSHRINKVR(%)  6.8900E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1684.3971540439481     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -949.24632748020997     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.39
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1684.397       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.25E-02  1.01E+00  1.60E+00  6.77E-01  9.34E-01  1.14E+01  7.77E-01  7.49E-01  7.40E-01  9.97E-01
 


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
+        1.34E+03
 
 TH 2
+       -6.59E+02  8.64E+07
 
 TH 3
+       -7.06E+01  7.23E+03  1.16E+03
 
 TH 4
+        3.06E+01  1.51E+04 -5.70E+02  2.97E+05
 
 TH 5
+        1.52E+02 -1.92E+04 -3.06E+03  1.22E+03  6.44E+06
 
 TH 6
+        1.99E+01  8.11E+02  5.01E+01 -5.69E+01 -1.57E+02  2.43E+02
 
 TH 7
+       -1.10E+00 -6.14E+02  1.28E+01  2.67E+01 -3.40E+01  1.52E+00  2.94E+02
 
 TH 8
+       -6.81E+00  1.73E+03  1.47E+02 -1.32E+02 -5.32E+02 -8.70E+00  3.12E+00  1.17E+02
 
 TH 9
+       -7.81E+01  1.19E+04  1.44E+03 -8.72E+02 -8.91E+06  1.46E+02  2.12E+01  4.13E+02  3.59E+03
 
 TH10
+       -1.41E+01  2.17E+03  3.80E+02 -1.94E+02 -1.24E+03  1.62E+01  3.88E+00  1.59E+02  8.20E+02  3.24E+02
 
 TH11
+       -1.20E+01  4.33E+02  6.26E+01 -4.56E+01 -1.60E+02  2.26E+00  8.11E-01 -1.24E+01  9.95E+01  3.96E+01  2.36E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.246
Stop Time:
Sat Sep 25 12:28:29 CDT 2021
