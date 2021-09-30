Wed Sep 29 11:18:39 CDT 2021
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
$DATA ../../../../data/spa/B/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1690.47054049652        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2429E+02 -1.0604E+01 -2.4312E+01  9.0600E+00  1.0184E+00  3.7620E+01  1.0251E+01  1.8434E+01  1.1257E+01  2.4757E+01
             7.5090E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1695.30432025768        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9359E-01  1.0360E+00  1.1812E+00  1.0221E+00  1.0768E+00  1.0245E+00  9.4149E-01  8.7297E-01  9.6188E-01  8.6307E-01
             9.8332E-01
 PARAMETER:  9.3571E-02  1.3534E-01  2.6649E-01  1.2188E-01  1.7401E-01  1.2423E-01  3.9710E-02 -3.5849E-02  6.1133E-02 -4.7262E-02
             8.3181E-02
 GRADIENT:  -4.6007E+00  1.4275E+01  1.8662E+01 -1.2270E+01 -4.3619E+00  5.3746E+00  2.6674E+00  3.3711E+00 -9.5719E+00 -1.7206E+01
            -1.4866E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1698.97572579520        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  9.9658E-01  9.1540E-01  1.2558E+00  1.1105E+00  1.0809E+00  9.8660E-01  7.8083E-01  5.6466E-01  9.6226E-01  9.9012E-01
             1.0297E+00
 PARAMETER:  9.6577E-02  1.1600E-02  3.2778E-01  2.0482E-01  1.7783E-01  8.6509E-02 -1.4740E-01 -4.7153E-01  6.1529E-02  9.0070E-02
             1.2926E-01
 GRADIENT:   3.4714E+00  1.6191E+01  5.3024E+00  1.4305E+01  7.2806E+00 -8.7867E+00 -7.9356E-01  4.6357E-01 -8.9008E+00 -8.1200E+00
             6.2816E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.69855151253        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  9.9423E-01  7.5728E-01  1.1347E+00  1.2057E+00  9.5285E-01  1.0060E+00  1.0116E+00  3.2994E-01  8.9884E-01  9.7574E-01
             9.8950E-01
 PARAMETER:  9.4216E-02 -1.7802E-01  2.2634E-01  2.8708E-01  5.1703E-02  1.0596E-01  1.1154E-01 -1.0088E+00 -6.6448E-03  7.5440E-02
             8.9441E-02
 GRADIENT:   6.4753E-02  1.4330E+01  6.7550E+00  1.6247E+01 -1.3229E+01 -8.2218E-01  5.5330E-01  3.8764E-01 -1.6903E+00  1.8531E+00
            -3.8553E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1701.68614489693        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  9.9120E-01  5.0188E-01  1.2332E+00  1.3657E+00  9.1797E-01  1.0019E+00  8.9033E-01  9.7738E-02  8.4017E-01  1.0199E+00
             1.0046E+00
 PARAMETER:  9.1163E-02 -5.8940E-01  3.0961E-01  4.1164E-01  1.4405E-02  1.0187E-01 -1.6164E-02 -2.2255E+00 -7.4152E-02  1.1970E-01
             1.0457E-01
 GRADIENT:   7.7528E-01  8.5211E+00  3.1534E+00  2.2301E+01 -4.4644E+00 -7.6416E-01 -1.5352E-02 -8.9244E-03 -1.5205E+00 -2.1879E-01
             7.7864E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1702.15684730351        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  9.8778E-01  3.1753E-01  1.2345E+00  1.4695E+00  8.6664E-01  1.0009E+00  7.3771E-01  1.4469E-02  7.8893E-01  1.0161E+00
             9.9828E-01
 PARAMETER:  8.7708E-02 -1.0472E+00  3.1064E-01  4.8490E-01 -4.3126E-02  1.0093E-01 -2.0421E-01 -4.1357E+00 -1.3707E-01  1.1601E-01
             9.8283E-02
 GRADIENT:  -5.0682E-01  3.0969E+00 -1.6853E+00  9.6822E+00  1.5903E+00  7.1009E-02 -8.5743E-02 -1.2150E-04 -2.3131E+00 -6.1664E-01
            -9.1021E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1702.17806490942        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1052
 NPARAMETR:  9.8664E-01  2.5369E-01  1.2470E+00  1.5105E+00  8.5146E-01  9.9962E-01  6.3979E-01  1.0000E-02  7.7284E-01  1.0203E+00
             9.9901E-01
 PARAMETER:  8.6548E-02 -1.2717E+00  3.2072E-01  5.1243E-01 -6.0801E-02  9.9625E-02 -3.4661E-01 -5.1248E+00 -1.5768E-01  1.2005E-01
             9.9007E-02
 GRADIENT:  -6.1132E-01  3.4681E+00  1.8407E+00  1.4568E+01 -4.6033E+00  1.4754E-02 -3.8754E-02  0.0000E+00 -1.5553E+00  2.5437E-01
            -6.7799E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1702.18260457762        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1227
 NPARAMETR:  9.8584E-01  2.0828E-01  1.2556E+00  1.5382E+00  8.4145E-01  9.9866E-01  5.5397E-01  1.0000E-02  7.6154E-01  1.0226E+00
             9.9976E-01
 PARAMETER:  8.5742E-02 -1.4689E+00  3.2760E-01  5.3063E-01 -7.2633E-02  9.8658E-02 -4.9065E-01 -6.0465E+00 -1.7242E-01  1.2232E-01
             9.9756E-02
 GRADIENT:  -5.4616E-01  3.0251E+00  3.7768E+00  1.4565E+01 -7.6930E+00 -2.8023E-02 -1.3077E-02  0.0000E+00 -7.5686E-01  7.4960E-01
            -3.7974E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1702.19378636513        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1402
 NPARAMETR:  9.8513E-01  1.6721E-01  1.2618E+00  1.5622E+00  8.3305E-01  9.9782E-01  4.6706E-01  1.0000E-02  7.5068E-01  1.0232E+00
             1.0003E+00
 PARAMETER:  8.5014E-02 -1.6885E+00  3.3253E-01  5.4611E-01 -8.2667E-02  9.7816E-02 -6.6129E-01 -7.1206E+00 -1.8677E-01  1.2294E-01
             1.0028E-01
 GRADIENT:  -4.1864E-01  2.2861E+00  4.1944E+00  1.2074E+01 -8.0918E+00 -4.7255E-02 -2.1182E-03  0.0000E+00 -2.3244E-01  8.8501E-01
            -1.4961E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1702.19925461457        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1577
 NPARAMETR:  9.8452E-01  1.3356E-01  1.2663E+00  1.5816E+00  8.2669E-01  9.9716E-01  3.8875E-01  1.0000E-02  7.4135E-01  1.0230E+00
             1.0006E+00
 PARAMETER:  8.4397E-02 -1.9132E+00  3.3609E-01  5.5846E-01 -9.0321E-02  9.7152E-02 -8.4482E-01 -8.2548E+00 -1.9928E-01  1.2277E-01
             1.0057E-01
 GRADIENT:  -3.1235E-01  1.6580E+00  3.7578E+00  9.3675E+00 -7.0837E+00 -4.9219E-02  7.0151E-04  0.0000E+00 -1.0683E-02  8.0750E-01
            -3.4238E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1702.27973436428        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1762
 NPARAMETR:  9.8454E-01  1.1445E-01  1.2662E+00  1.5844E+00  8.2578E-01  9.9706E-01  2.9524E-01  1.0000E-02  7.3768E-01  1.0208E+00
             1.0003E+00
 PARAMETER:  8.4416E-02 -2.0676E+00  3.3606E-01  5.6021E-01 -9.1425E-02  9.7054E-02 -1.1200E+00 -8.6864E+00 -2.0424E-01  1.2057E-01
             1.0032E-01
 GRADIENT:   7.6854E-01  1.5989E-01 -2.4604E-01 -1.2155E+01  1.6484E+00  9.1818E-02  1.7243E-03  0.0000E+00  9.9444E-01 -3.4314E-02
             1.9054E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1702.28408309667        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1956
 NPARAMETR:  9.8483E-01  1.1380E-01  1.2642E+00  1.5829E+00  8.2493E-01  9.9706E-01  9.6488E-02  1.0000E-02  7.3605E-01  1.0203E+00
             1.0001E+00
 PARAMETER:  8.4711E-02 -2.0733E+00  3.3447E-01  5.5926E-01 -9.2461E-02  9.7055E-02 -2.2383E+00 -8.6864E+00 -2.0646E-01  1.2011E-01
             1.0011E-01
 GRADIENT:   1.4730E+00  3.5431E-04 -3.9483E-01 -1.6261E+01  2.0932E+00  1.0790E-01  3.0393E-04  0.0000E+00  2.1409E-01 -5.6018E-03
             1.6302E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1702.28603659887        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2146             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8480E-01  1.1393E-01  1.2631E+00  1.5831E+00  8.2397E-01  9.9705E-01  7.2840E-02  1.0000E-02  7.3589E-01  1.0197E+00
             9.9993E-01
 PARAMETER:  8.4679E-02 -2.0721E+00  3.3355E-01  5.5939E-01 -9.3627E-02  9.7046E-02 -2.5195E+00 -8.6864E+00 -2.0668E-01  1.1953E-01
             9.9932E-02
 GRADIENT:   4.0500E+02  1.0077E+01  7.2990E+00  9.1794E+02  8.3898E+00  3.9327E+01  6.0472E-03  0.0000E+00  1.8936E+01  8.9906E-01
             8.5968E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1702.28734400678        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2339
 NPARAMETR:  9.8484E-01  1.1455E-01  1.2610E+00  1.5823E+00  8.2354E-01  9.9707E-01  5.8328E-02  1.0000E-02  7.3610E-01  1.0190E+00
             9.9992E-01
 PARAMETER:  8.4723E-02 -2.0667E+00  3.3194E-01  5.5887E-01 -9.4139E-02  9.7066E-02 -2.7417E+00 -8.6864E+00 -2.0639E-01  1.1881E-01
             9.9917E-02
 GRADIENT:   1.4550E+00  2.9896E-02 -3.4675E-01 -1.5930E+01  1.7556E+00  1.0599E-01  1.4494E-04  0.0000E+00  1.0290E-01 -6.7509E-03
             1.1656E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1702.28902629446        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2533             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8485E-01  1.1658E-01  1.2601E+00  1.5818E+00  8.2208E-01  9.9708E-01  4.0843E-02  1.0000E-02  7.3631E-01  1.0185E+00
             9.9965E-01
 PARAMETER:  8.4730E-02 -2.0492E+00  3.3117E-01  5.5853E-01 -9.5912E-02  9.7075E-02 -3.0980E+00 -8.6864E+00 -2.0611E-01  1.1830E-01
             9.9647E-02
 GRADIENT:   4.0519E+02  1.0594E+01  8.6694E+00  9.1601E+02  5.9044E+00  3.9344E+01  2.4549E-03  0.0000E+00  1.8701E+01  1.0194E+00
             6.9726E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1702.29050710440        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2719             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8487E-01  1.1680E-01  1.2585E+00  1.5814E+00  8.2169E-01  9.9709E-01  3.5864E-02  1.0000E-02  7.3643E-01  1.0179E+00
             9.9964E-01
 PARAMETER:  8.4753E-02 -2.0473E+00  3.2993E-01  5.5832E-01 -9.6393E-02  9.7090E-02 -3.2280E+00 -8.6864E+00 -2.0594E-01  1.1772E-01
             9.9644E-02
 GRADIENT:   4.0520E+02  1.0597E+01  8.2946E+00  9.1513E+02  6.3808E+00  3.9341E+01  1.9813E-03  0.0000E+00  1.8697E+01  9.8551E-01
             7.4008E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1702.29219115148        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2908
 NPARAMETR:  9.8487E-01  1.1776E-01  1.2569E+00  1.5806E+00  8.2138E-01  9.9710E-01  2.7068E-02  1.0000E-02  7.3678E-01  1.0172E+00
             9.9953E-01
 PARAMETER:  8.4756E-02 -2.0391E+00  3.2862E-01  5.5782E-01 -9.6772E-02  9.7098E-02 -3.5094E+00 -8.6864E+00 -2.0547E-01  1.1704E-01
             9.9528E-02
 GRADIENT:   1.3622E+00  1.9646E-01  8.7719E-01 -1.4204E+01 -5.8860E-01  9.1424E-02  4.6986E-05  0.0000E+00 -6.1301E-02  8.2934E-02
            -7.6115E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1702.29285301502        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3095
 NPARAMETR:  9.8465E-01  1.1786E-01  1.2557E+00  1.5804E+00  8.2101E-01  9.9700E-01  2.3562E-02  1.0000E-02  7.3699E-01  1.0168E+00
             9.9953E-01
 PARAMETER:  8.4530E-02 -2.0382E+00  3.2768E-01  5.5771E-01 -9.7224E-02  9.6995E-02 -3.6481E+00 -8.6864E+00 -2.0518E-01  1.1663E-01
             9.9527E-02
 GRADIENT:   8.5506E-01  1.8314E-01  6.5466E-01 -1.4232E+01 -3.3848E-01  4.9771E-02  3.8291E-05  0.0000E+00  1.9702E-02  7.6888E-02
            -3.6704E-02

0ITERATION NO.:   88    OBJECTIVE VALUE:  -1702.29324286968        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     3200
 NPARAMETR:  9.8489E-01  1.1919E-01  1.2550E+00  1.5799E+00  8.2035E-01  9.9712E-01  2.2561E-02  1.0000E-02  7.3705E-01  1.0164E+00
             9.9941E-01
 PARAMETER:  8.4791E-02 -2.0399E+00  3.2685E-01  5.5733E-01 -9.7026E-02  9.7132E-02 -3.6846E+00 -8.6864E+00 -2.0501E-01  1.1628E-01
             9.9606E-02
 GRADIENT:   2.1385E-02 -1.0876E-01 -1.7660E-01 -1.4920E-01  6.5260E-01  4.0618E-03  3.9839E-06  0.0000E+00  2.2709E-02 -3.8106E-03
             6.9322E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3200
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.4583E-05 -9.2906E-05 -1.9853E-04 -5.6203E-03 -2.3074E-02
 SE:             2.9828E-02  4.9980E-05  1.8067E-04  2.9208E-02  2.5292E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9907E-01  6.3045E-02  2.7183E-01  8.4741E-01  3.6161E-01

 ETASHRINKSD(%)  7.1815E-02  9.9833E+01  9.9395E+01  2.1489E+00  1.5270E+01
 ETASHRINKVR(%)  1.4358E-01  1.0000E+02  9.9996E+01  4.2517E+00  2.8208E+01
 EBVSHRINKSD(%)  4.0942E-01  9.9844E+01  9.9403E+01  2.2719E+00  1.2356E+01
 EBVSHRINKVR(%)  8.1717E-01  1.0000E+02  9.9996E+01  4.4921E+00  2.3185E+01
 RELATIVEINF(%)  9.6718E+01  1.1431E-05  3.9031E-04  5.8801E+00  5.4526E+00
 EPSSHRINKSD(%)  4.1674E+01
 EPSSHRINKVR(%)  6.5981E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1702.2932428696815     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -967.14241630594336     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1702.293       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.18E-01  1.25E+00  1.58E+00  8.21E-01  9.97E-01  2.27E-02  1.00E-02  7.37E-01  1.02E+00  1.00E+00
 


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
+       -3.02E+01  3.90E+02
 
 TH 3
+        1.13E+00  1.01E+02  2.56E+02
 
 TH 4
+       -9.96E+00  4.68E+02 -5.30E+01  7.83E+02
 
 TH 5
+        2.63E+00 -3.36E+02 -5.13E+02 -3.78E+01  1.24E+03
 
 TH 6
+        3.97E-01 -5.06E+00  3.29E-01 -2.64E+00 -8.61E-01  1.97E+02
 
 TH 7
+        3.39E-02 -3.59E-02  2.53E-02 -1.04E-02  6.52E-02  9.32E-02  7.16E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.92E+00 -9.81E+01  6.24E+00 -2.38E+00 -3.25E+00 -7.45E-01 -1.31E-02  0.00E+00  3.35E+02
 
 TH10
+        3.51E-01  1.27E+01 -8.47E+00 -1.84E+00 -7.69E+01  7.60E-01 -3.05E-03  0.00E+00  4.94E-01  1.04E+02
 
 TH11
+       -8.24E+00 -1.71E+01 -3.25E+01 -1.05E+01  2.38E+01  1.56E+00  4.41E-02  0.00E+00  1.03E+01  3.35E+01  2.29E+02
 
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
 #CPUT: Total CPU Time in Seconds,       45.873
Stop Time:
Wed Sep 29 11:19:26 CDT 2021
