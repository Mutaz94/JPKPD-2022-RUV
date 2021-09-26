Sat Sep 25 12:36:39 CDT 2021
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
$DATA ../../../../data/spa/S2/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1700.48627599293        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.1441E+01 -7.0681E+01 -3.1247E+01 -4.1125E+01  4.4952E+01  3.4041E+01 -1.5515E+01 -4.0758E+00  7.4813E+00  2.0528E+01
             3.4329E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1711.89507803460        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8943E-01  1.1467E+00  1.0032E+00  9.3758E-01  1.0082E+00  8.9418E-01  1.1989E+00  1.1402E+00  9.4313E-01  7.8910E-01
             9.2015E-01
 PARAMETER:  8.9371E-02  2.3693E-01  1.0322E-01  3.5544E-02  1.0816E-01 -1.1851E-02  2.8141E-01  2.3116E-01  4.1448E-02 -1.3686E-01
             1.6783E-02
 GRADIENT:   6.8753E+00  9.9805E+00  1.0920E+01 -7.3440E+00 -2.5932E+00 -8.1805E+00  4.3645E+00 -8.5785E+00 -8.4850E-01  3.5650E+00
            -1.7353E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1712.65780790610        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.8975E-01  1.0830E+00  9.6885E-01  9.7500E-01  9.5637E-01  9.1467E-01  1.1966E+00  1.3287E+00  9.2754E-01  6.7631E-01
             8.9429E-01
 PARAMETER:  8.9696E-02  1.7977E-01  6.8358E-02  7.4679E-02  5.5389E-02  1.0807E-02  2.7951E-01  3.8423E-01  2.4780E-02 -2.9111E-01
            -1.1726E-02
 GRADIENT:   1.2617E+01 -2.4489E+00 -2.1048E+00  2.5022E+00 -5.7027E+00  1.1667E+00 -2.4152E+00  1.5330E+00  3.4178E+00  3.4718E-01
            -1.3382E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1713.14530818916        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.8553E-01  1.0724E+00  1.1460E+00  9.9251E-01  1.0210E+00  9.1217E-01  1.2030E+00  1.5255E+00  9.1309E-01  7.3436E-01
             9.2411E-01
 PARAMETER:  8.5427E-02  1.6987E-01  2.3627E-01  9.2477E-02  1.2082E-01  8.0703E-03  2.8483E-01  5.2235E-01  9.0758E-03 -2.0875E-01
             2.1078E-02
 GRADIENT:  -1.1961E-01 -1.5355E-01 -1.6320E-01 -8.8757E-02  3.0111E-01  3.3184E-03  5.9194E-02  1.1360E-02  3.1452E-02 -2.0136E-02
             1.6682E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1713.74823130879        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      393
 NPARAMETR:  1.0047E+00  1.0988E+00  1.2166E+00  9.8699E-01  1.0530E+00  9.2295E-01  1.1800E+00  1.6455E+00  9.2841E-01  7.6154E-01
             9.2483E-01
 PARAMETER:  1.0465E-01  1.9419E-01  2.9602E-01  8.6909E-02  1.5161E-01  1.9823E-02  2.6554E-01  5.9804E-01  2.5715E-02 -1.7241E-01
             2.1859E-02
 GRADIENT:   1.1643E+00 -1.3053E-01 -6.4747E-02  1.2181E+00  2.9475E-01  6.8452E-02  7.4702E-03 -9.9101E-02 -4.0036E-02 -2.8313E-02
             4.3874E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1713.85162124465        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      568
 NPARAMETR:  1.0064E+00  1.2311E+00  1.0896E+00  9.0047E-01  1.0644E+00  9.2521E-01  1.0924E+00  1.6386E+00  9.8443E-01  7.5945E-01
             9.2449E-01
 PARAMETER:  1.0633E-01  3.0795E-01  1.8580E-01 -4.8382E-03  1.6243E-01  2.2265E-02  1.8834E-01  5.9383E-01  8.4306E-02 -1.7516E-01
             2.1484E-02
 GRADIENT:   3.0936E+00  2.4092E+00  4.6069E-01  2.6805E+00 -6.9161E-01  6.3852E-01  4.3033E-01 -1.7163E-01 -2.0053E-01  4.8743E-03
            -7.8042E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1713.85691322534        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  1.0064E+00  1.2360E+00  1.0844E+00  8.9703E-01  1.0648E+00  9.2371E-01  1.0890E+00  1.6374E+00  9.8692E-01  7.5936E-01
             9.2452E-01
 PARAMETER:  1.0637E-01  3.1188E-01  1.8100E-01 -8.6653E-03  1.6283E-01  2.0645E-02  1.8524E-01  5.9310E-01  8.6836E-02 -1.7528E-01
             2.1523E-02
 GRADIENT:   3.1192E+00  2.2001E+00  4.5489E-01  2.4725E+00 -5.7924E-01 -8.3114E-03  3.9892E-01 -1.9265E-01 -2.0060E-01 -9.6352E-03
            -6.9902E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1713.86726897947        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      941
 NPARAMETR:  1.0051E+00  1.2435E+00  1.0772E+00  8.9125E-01  1.0652E+00  9.2373E-01  1.0780E+00  1.6399E+00  9.9524E-01  7.5823E-01
             9.2420E-01
 PARAMETER:  1.0512E-01  3.1795E-01  1.7433E-01 -1.5135E-02  1.6313E-01  2.0660E-02  1.7515E-01  5.9466E-01  9.5233E-02 -1.7676E-01
             2.1177E-02
 GRADIENT:  -2.1728E-01  1.2003E+00  6.3042E-01  1.7080E+00 -1.3568E+00 -2.2020E-02 -1.4944E-01 -1.3993E-01  1.7503E-01 -9.8876E-02
            -1.6678E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1713.92066176386        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1124
 NPARAMETR:  1.0058E+00  1.3686E+00  9.6885E-01  8.1355E-01  1.0782E+00  9.2421E-01  1.0024E+00  1.6489E+00  1.0611E+00  7.5820E-01
             9.2404E-01
 PARAMETER:  1.0580E-01  4.1378E-01  6.8357E-02 -1.0634E-01  1.7532E-01  2.1189E-02  1.0243E-01  6.0010E-01  1.5929E-01 -1.7681E-01
             2.1001E-02
 GRADIENT:  -3.0434E-01  1.0144E+01  2.4510E+00  6.8903E+00 -5.5113E+00 -8.5713E-02 -9.5764E-01 -6.8703E-01 -8.1092E-02 -3.6017E-01
            -3.7001E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1713.97964710557        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1306
 NPARAMETR:  1.0059E+00  1.3823E+00  9.4960E-01  7.9906E-01  1.0813E+00  9.2450E-01  9.9479E-01  1.6565E+00  1.0753E+00  7.5962E-01
             9.2413E-01
 PARAMETER:  1.0589E-01  4.2375E-01  4.8282E-02 -1.2432E-01  1.7812E-01  2.1501E-02  9.4772E-02  6.0469E-01  1.7261E-01 -1.7494E-01
             2.1093E-02
 GRADIENT:  -1.8906E-01  2.5737E+00  1.2328E+00  2.1504E+00 -2.8259E+00  1.7964E-02 -5.3732E-01 -1.4193E-01  5.7135E-01 -1.4623E-01
            -1.4870E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1713.98347206129        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1498
 NPARAMETR:  1.0060E+00  1.3825E+00  9.4446E-01  7.9878E-01  1.0814E+00  9.2448E-01  9.9457E-01  1.6563E+00  1.0756E+00  7.5966E-01
             9.2415E-01
 PARAMETER:  1.0599E-01  4.2392E-01  4.2861E-02 -1.2467E-01  1.7821E-01  2.1476E-02  9.4556E-02  6.0457E-01  1.7286E-01 -1.7488E-01
             2.1118E-02
 GRADIENT:   4.9122E-03  1.3871E+00 -1.6463E-04  3.0585E+00 -3.5343E-01 -2.6813E-04 -5.4170E-01  1.8180E-01  6.9981E-01 -1.1322E-01
            -1.0075E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1713.99565754436        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1678
 NPARAMETR:  1.0061E+00  1.3822E+00  9.2348E-01  7.9595E-01  1.0753E+00  9.2462E-01  1.0027E+00  1.6092E+00  1.0666E+00  7.5618E-01
             9.2409E-01
 PARAMETER:  1.0612E-01  4.2371E-01  2.0397E-02 -1.2821E-01  1.7260E-01  2.1623E-02  1.0272E-01  5.7572E-01  1.6450E-01 -1.7947E-01
             2.1057E-02
 GRADIENT:   1.4915E-01 -7.5716E-01 -3.2660E-02  5.2849E-01  8.9869E-01  4.2427E-02  7.6155E-03 -1.4890E-01 -1.9355E-01  5.4098E-02
            -1.1077E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1713.99843810148        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1854            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0061E+00  1.3837E+00  9.0880E-01  7.9464E-01  1.0691E+00  9.2457E-01  1.0042E+00  1.5889E+00  1.0666E+00  7.4891E-01
             9.2426E-01
 PARAMETER:  1.0612E-01  4.2476E-01  4.3652E-03 -1.2987E-01  1.6684E-01  2.1570E-02  1.0422E-01  5.6301E-01  1.6451E-01 -1.8914E-01
             2.1239E-02
 GRADIENT:   5.3076E+01  3.1501E+01  3.3977E-01  7.5810E+00  7.3748E-01  5.3481E+00  7.3528E-01  1.3089E-01  8.4769E-01  7.8022E-02
             5.9822E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1714.02339224493        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     2018
 NPARAMETR:  1.0096E+00  1.4349E+00  8.7801E-01  7.6077E-01  1.0845E+00  9.2663E-01  9.7675E-01  1.6551E+00  1.0966E+00  7.5306E-01
             9.2515E-01
 PARAMETER:  1.0953E-01  4.6108E-01 -3.0099E-02 -1.7342E-01  1.8109E-01  2.3797E-02  7.6480E-02  6.0385E-01  1.9225E-01 -1.8361E-01
             2.2196E-02
 GRADIENT:   8.3342E+00 -2.1815E+00 -1.3768E+00  1.4537E+00  3.5933E+00  8.0532E-01  2.6298E-01  6.4447E-01 -1.1312E-01 -3.1981E-01
             2.6073E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1714.14051362788        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2194
 NPARAMETR:  1.0048E+00  1.6042E+00  7.7822E-01  6.5163E-01  1.1301E+00  9.2406E-01  8.9178E-01  1.7785E+00  1.2312E+00  7.8557E-01
             9.2392E-01
 PARAMETER:  1.0476E-01  5.7260E-01 -1.5074E-01 -3.2828E-01  2.2232E-01  2.1022E-02 -1.4540E-02  6.7580E-01  3.0802E-01 -1.4135E-01
             2.0872E-02
 GRADIENT:  -4.9901E+00  3.9980E+00  1.1246E+00  2.2311E+00  1.1480E+00 -3.5578E-01 -6.5704E-01 -3.1105E-01 -1.5842E-02 -1.7508E-01
            -3.7927E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1714.29128554900        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2371
 NPARAMETR:  1.0090E+00  1.8070E+00  5.6944E-01  5.1124E-01  1.1569E+00  9.2481E-01  8.2618E-01  1.8143E+00  1.4166E+00  7.8560E-01
             9.2475E-01
 PARAMETER:  1.0897E-01  6.9164E-01 -4.6310E-01 -5.7092E-01  2.4577E-01  2.1829E-02 -9.0940E-02  6.9570E-01  4.4827E-01 -1.4130E-01
             2.1768E-02
 GRADIENT:   4.8760E+00  5.9153E+00  2.1106E+00 -6.8935E-01  1.3125E-02  2.8217E-03 -7.5979E-01 -1.5518E-01 -1.1927E+00 -9.2210E-01
            -5.8815E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1714.37384980291        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2552
 NPARAMETR:  1.0076E+00  1.8348E+00  5.0608E-01  4.8453E-01  1.1468E+00  9.2434E-01  8.2143E-01  1.7424E+00  1.4529E+00  7.7181E-01
             9.2558E-01
 PARAMETER:  1.0753E-01  7.0695E-01 -5.8105E-01 -6.2458E-01  2.3693E-01  2.1323E-02 -9.6703E-02  6.5527E-01  4.7358E-01 -1.5902E-01
             2.2660E-02
 GRADIENT:   1.0566E+00 -3.8152E+00  1.2915E+00 -3.6742E+00 -6.5360E-01 -1.6144E-01  7.6431E-02  6.1165E-01 -7.1748E-01 -5.6658E-01
            -1.9233E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1714.37621601445        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     2718
 NPARAMETR:  1.0076E+00  1.8345E+00  5.0616E-01  4.8445E-01  1.1468E+00  9.2467E-01  8.2016E-01  1.7421E+00  1.4617E+00  7.7598E-01
             9.2515E-01
 PARAMETER:  1.0755E-01  7.0677E-01 -5.8091E-01 -6.2473E-01  2.3699E-01  2.1686E-02 -9.8254E-02  6.5510E-01  4.7962E-01 -1.5362E-01
             2.2205E-02
 GRADIENT:   1.1419E+00 -4.8965E+00  1.0437E+00 -3.3472E+00 -1.6267E+00 -3.0072E-02 -3.7346E-03  7.9888E-01  2.5681E-02  1.4590E-02
             8.1999E-03

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1714.44553509718        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     2877
 NPARAMETR:  1.0089E+00  1.8340E+00  4.9137E-01  4.8750E-01  1.1394E+00  9.2622E-01  8.2520E-01  1.5883E+00  1.4721E+00  7.7592E-01
             9.2511E-01
 PARAMETER:  1.0883E-01  7.0652E-01 -6.1056E-01 -6.1847E-01  2.3052E-01  2.3361E-02 -9.2125E-02  5.6268E-01  4.8666E-01 -1.5370E-01
             2.2155E-02
 GRADIENT:   4.0098E+00  1.3246E+00  7.0126E-01  3.2929E-01  1.2812E-01  4.8893E-01  7.5827E-01  2.8250E-01  9.7813E-01  2.9423E-01
             1.8022E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1714.52229733710        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3054
 NPARAMETR:  1.0076E+00  1.8387E+00  4.2148E-01  4.7874E-01  1.1057E+00  9.2532E-01  8.2693E-01  1.2044E+00  1.4513E+00  7.7592E-01
             9.2511E-01
 PARAMETER:  1.0757E-01  7.0908E-01 -7.6399E-01 -6.3660E-01  2.0044E-01  2.2383E-02 -9.0034E-02  2.8597E-01  4.7249E-01 -1.5370E-01
             2.2154E-02
 GRADIENT:   3.3031E-02  6.1708E-01  3.9727E-01 -5.3815E-02 -2.9486E+00 -6.9322E-03  1.1431E-01  3.6127E-02 -7.6818E-01  3.7454E+00
             1.0882E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1714.54953925239        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     3233
 NPARAMETR:  1.0078E+00  1.8654E+00  3.8779E-01  4.5962E-01  1.1026E+00  9.2497E-01  8.1963E-01  1.1046E+00  1.4786E+00  7.7593E-01
             9.2510E-01
 PARAMETER:  1.0772E-01  7.2348E-01 -8.4730E-01 -6.7736E-01  1.9763E-01  2.2007E-02 -9.8902E-02  1.9947E-01  4.9113E-01 -1.5370E-01
             2.2152E-02
 GRADIENT:   4.2229E-01  2.2888E+00  6.3080E-01 -2.4150E-01 -5.0462E+00 -1.3243E-01  7.0283E-02  1.0140E-01 -1.0369E+00  4.3981E+00
             1.1706E+00

0ITERATION NO.:  104    OBJECTIVE VALUE:  -1714.54979251514        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     3372
 NPARAMETR:  1.0077E+00  1.8654E+00  3.8775E-01  4.5964E-01  1.1025E+00  9.2512E-01  8.1957E-01  1.1020E+00  1.4786E+00  7.7592E-01
             9.2511E-01
 PARAMETER:  1.0772E-01  7.2348E-01 -8.4730E-01 -6.7725E-01  1.9763E-01  2.2064E-02 -9.8882E-02  1.9713E-01  4.9113E-01 -1.5370E-01
             2.2152E-02
 GRADIENT:   3.2557E+05 -4.4747E+00  4.1407E+04  1.0355E+05  3.5491E+05 -1.0385E-01  5.1620E-02 -3.5590E+05  1.4278E+05  2.2819E+05
            -7.0143E+05
 NUMSIGDIG:         3.3         7.9         3.3         3.3         3.3         2.3         2.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3372
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2195E-05 -2.5164E-02 -2.2351E-02  2.9109E-02 -4.0286E-02
 SE:             2.9866E-02  2.6130E-02  8.7388E-03  2.1754E-02  1.9194E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9914E-01  3.3553E-01  1.0538E-02  1.8086E-01  3.5829E-02

 ETASHRINKSD(%)  1.0000E-10  1.2462E+01  7.0724E+01  2.7121E+01  3.5696E+01
 ETASHRINKVR(%)  1.0000E-10  2.3372E+01  9.1429E+01  4.6887E+01  5.8650E+01
 EBVSHRINKSD(%)  4.2777E-01  1.2983E+01  7.2444E+01  2.7964E+01  3.2531E+01
 EBVSHRINKVR(%)  8.5371E-01  2.4280E+01  9.2407E+01  4.8108E+01  5.4479E+01
 RELATIVEINF(%)  9.9033E+01  5.7318E+00  4.4596E-01  2.8075E+00  1.0079E+01
 EPSSHRINKSD(%)  4.5660E+01
 EPSSHRINKVR(%)  7.0472E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.5497925151387     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -979.39896595140056     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.550       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.87E+00  3.88E-01  4.60E-01  1.10E+00  9.25E-01  8.20E-01  1.10E+00  1.48E+00  7.76E-01  9.25E-01
 


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
+        2.98E+09
 
 TH 2
+       -5.64E+00  1.93E+07
 
 TH 3
+        3.48E+05  1.49E+02  1.63E+08
 
 TH 4
+       -1.42E+05  2.80E+02 -4.73E+04  1.81E+08
 
 TH 5
+        2.86E+03 -1.84E+02  1.73E+05 -7.01E+04  3.69E+08
 
 TH 6
+        1.14E+04 -1.15E+00  3.76E+03 -6.09E+08  5.67E+03  2.05E+09
 
 TH 7
+       -3.15E+03  1.51E+01 -1.06E+03  6.87E+08 -1.56E+03 -2.31E+09  2.61E+09
 
 TH 8
+       -3.91E+05 -7.44E+00 -1.29E+05  7.08E+04 -1.95E+05 -5.69E+03  1.58E+03  3.72E+08
 
 TH 9
+        2.23E+08 -1.32E+01  7.35E+07 -2.11E+04  1.11E+08  1.70E+03 -4.59E+02 -5.84E+04  6.65E+07
 
 TH10
+        3.09E+03 -1.24E+01  1.00E+03 -1.29E+05  1.45E+03  1.04E+04 -2.85E+03 -6.77E+08  4.70E+02  1.23E+09
 
 TH11
+       -2.82E+04 -1.10E+01 -9.32E+03  1.66E+05 -1.41E+04 -1.33E+04  3.69E+03  4.59E+05 -2.61E+08 -3.62E+03  2.05E+09
 
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
 #CPUT: Total CPU Time in Seconds,       49.804
Stop Time:
Sat Sep 25 12:37:30 CDT 2021
