Wed Sep 29 12:54:36 CDT 2021
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
$DATA ../../../../data/spa/A2/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1193.43494971493        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6381E+02 -2.4831E+01  7.7653E-01 -1.5829E+01  8.7075E+01  7.0892E+01 -4.2278E+00  1.7292E+00 -7.6685E+00 -1.8053E+01
            -9.5877E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1490.63833096100        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1610E+00  1.0663E+00  1.0825E+00  1.0358E+00  9.8292E-01  1.0187E+00  9.0640E-01  9.2309E-01  9.0158E-01  8.0792E-01
             2.3248E+00
 PARAMETER:  2.4932E-01  1.6423E-01  1.7926E-01  1.3515E-01  8.2777E-02  1.1852E-01  1.7231E-03  1.9970E-02 -3.6023E-03 -1.1329E-01
             9.4362E-01
 GRADIENT:   4.2610E+02  7.3899E+00  1.3118E+01 -9.3841E+00 -2.0809E+01  3.0011E+01  5.1526E+00  3.0180E+00 -1.2795E+00  1.0585E+01
            -2.2596E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1499.11103316690        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.1237E+00  1.2334E+00  6.1807E-01  9.1926E-01  8.1442E-01  9.9257E-01  8.5110E-01  5.2585E-01  1.1824E+00  4.8504E-01
             2.2160E+00
 PARAMETER:  2.1664E-01  3.0974E-01 -3.8115E-01  1.5818E-02 -1.0528E-01  9.2547E-02 -6.1227E-02 -5.4273E-01  2.6753E-01 -6.2353E-01
             8.9572E-01
 GRADIENT:   3.4688E+02  5.5629E+01  1.5862E+01  2.1876E+01 -3.1765E+01  3.1524E+01  7.1572E+00  1.3987E+00  2.8507E+01  2.5894E+00
            -2.0045E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1508.34203438392        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0346E+00  1.0749E+00  4.4089E-01  9.6874E-01  6.4039E-01  8.8101E-01  1.0097E+00  2.5993E-01  9.9274E-01  3.8516E-01
             2.1332E+00
 PARAMETER:  1.3403E-01  1.7225E-01 -7.1896E-01  6.8238E-02 -3.4568E-01 -2.6688E-02  1.0966E-01 -1.2474E+00  9.2713E-02 -8.5409E-01
             8.5762E-01
 GRADIENT:   9.4379E+01  2.4375E+01 -7.6680E+00  3.8126E+01  2.0499E+01  2.7103E+00  1.3881E+01  8.5870E-01  1.9257E+01  5.1541E+00
            -2.1816E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1512.62813869864        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  1.0507E+00  1.0013E+00  5.4028E-01  1.0183E+00  6.7428E-01  8.8293E-01  1.0165E+00  2.8279E-01  8.2113E-01  2.6411E-01
             2.3634E+00
 PARAMETER:  1.4945E-01  1.0127E-01 -5.1567E-01  1.1812E-01 -2.9411E-01 -2.4510E-02  1.1641E-01 -1.1630E+00 -9.7075E-02 -1.2314E+00
             9.6011E-01
 GRADIENT:   1.3008E+01  1.3168E+00  4.3900E+00 -8.0396E+00 -7.6990E+00 -1.2672E-01 -8.8834E-01  3.3816E-01 -2.3322E+00  8.9093E-01
             4.3168E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1513.01606616326        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0468E+00  9.5026E-01  4.7234E-01  1.0276E+00  6.1311E-01  8.8460E-01  1.0587E+00  9.9575E-02  8.1522E-01  1.4694E-01
             2.3340E+00
 PARAMETER:  1.4571E-01  4.8983E-02 -6.5006E-01  1.2718E-01 -3.8921E-01 -2.2622E-02  1.5702E-01 -2.2068E+00 -1.0430E-01 -1.8177E+00
             9.4760E-01
 GRADIENT:   3.0142E+00 -5.3646E+00 -1.3227E+00 -3.3658E+00  3.4430E+00 -4.6968E-01 -7.8292E-01  1.8139E-02 -6.3587E-01  1.1341E-01
             1.0592E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1513.20523779197        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      760
 NPARAMETR:  1.0413E+00  1.0960E+00  3.6909E-01  9.2161E-01  6.0155E-01  8.8866E-01  9.1478E-01  1.4720E-02  8.7086E-01  5.6793E-02
             2.3019E+00
 PARAMETER:  1.4049E-01  1.9171E-01 -8.9671E-01  1.8366E-02 -4.0824E-01 -1.8044E-02  1.0923E-02 -4.1186E+00 -3.8272E-02 -2.7683E+00
             9.3375E-01
 GRADIENT:  -6.6919E+00  4.6844E+00  9.2516E-01  3.1274E+00 -3.3922E+00  3.9452E-02  3.6122E-01  9.1998E-04 -1.9918E-01  6.9835E-02
             1.4598E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1513.37473833028        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  1.0392E+00  1.1847E+00  2.7576E-01  8.3484E-01  5.6841E-01  8.8980E-01  8.1962E-01  1.0000E-02  9.1249E-01  1.0000E-02
             2.2509E+00
 PARAMETER:  1.3842E-01  2.6945E-01 -1.1882E+00 -8.0509E-02 -4.6491E-01 -1.6755E-02 -9.8914E-02 -7.6860E+00  8.4229E-03 -4.6933E+00
             9.1132E-01
 GRADIENT:   1.0939E+00  2.8236E+01  1.4829E+01 -3.9203E+00 -3.3445E+01 -1.0093E-02 -8.6281E-01  0.0000E+00 -4.2712E-01  0.0000E+00
            -4.8750E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1515.18453475618        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1115
 NPARAMETR:  1.0271E+00  1.3091E+00  1.7270E-01  7.0537E-01  5.5168E-01  8.8854E-01  7.0760E-01  1.0000E-02  1.0329E+00  1.0000E-02
             2.1703E+00
 PARAMETER:  1.2676E-01  3.6932E-01 -1.6562E+00 -2.4903E-01 -4.9479E-01 -1.8180E-02 -2.4588E-01 -1.7055E+01  1.3239E-01 -1.0432E+01
             8.7486E-01
 GRADIENT:   5.6067E+00  2.6029E+01  1.1019E+01 -3.2098E+00 -3.8334E+01  2.3772E-03 -4.6885E+00  0.0000E+00  1.9588E+00  0.0000E+00
            -3.0025E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1523.80735916225        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1292
 NPARAMETR:  1.0068E+00  1.9306E+00  6.6529E-02  3.6977E-01  8.0411E-01  8.5964E-01  6.3011E-01  1.0000E-02  1.4567E+00  1.0000E-02
             2.1959E+00
 PARAMETER:  1.0681E-01  7.5782E-01 -2.6101E+00 -8.9488E-01 -1.1802E-01 -5.1239E-02 -3.6186E-01 -5.3969E+01  4.7616E-01 -3.4406E+01
             8.8659E-01
 GRADIENT:  -3.8758E+01  6.4512E+01 -9.1445E+00  2.5867E+01 -2.9843E+01 -9.1396E+00  9.4701E-01  0.0000E+00 -1.7183E+01  0.0000E+00
             3.3910E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1528.35762081277        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1481             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0257E+00  2.1004E+00  5.0266E-02  2.9454E-01  9.1688E-01  8.8363E-01  6.1790E-01  1.0000E-02  1.8576E+00  1.0000E-02
             2.1705E+00
 PARAMETER:  1.2538E-01  8.4212E-01 -2.8904E+00 -1.1223E+00  1.3227E-02 -2.3717E-02 -3.8143E-01 -6.8339E+01  7.1929E-01 -4.3326E+01
             8.7497E-01
 GRADIENT:   1.1147E+02  2.9955E+02  4.1964E+00  4.3703E+01  6.5578E+00  7.7942E+00  5.8732E+00  0.0000E+00 -9.5071E+00  0.0000E+00
             4.3375E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1530.15462385854        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1642
 NPARAMETR:  1.0245E+00  2.0644E+00  5.1161E-02  2.7413E-01  9.1543E-01  8.8051E-01  6.1273E-01  1.0000E-02  2.0397E+00  1.0000E-02
             2.1802E+00
 PARAMETER:  1.2417E-01  8.2482E-01 -2.8728E+00 -1.1942E+00  1.1633E-02 -2.7249E-02 -3.8983E-01 -6.8339E+01  8.1281E-01 -4.3326E+01
             8.7941E-01
 GRADIENT:   5.3023E+00 -2.3381E+01  7.5550E+00 -8.2682E-01  1.6352E+00  1.2966E+00  2.1098E+00  0.0000E+00  7.8788E-01  0.0000E+00
             9.1566E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1530.48619615564        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1831
 NPARAMETR:  1.0225E+00  2.0762E+00  4.9842E-02  2.7216E-01  9.1746E-01  8.7801E-01  6.0821E-01  1.0000E-02  2.0583E+00  1.0000E-02
             2.1607E+00
 PARAMETER:  1.2221E-01  8.3056E-01 -2.8989E+00 -1.2014E+00  1.3856E-02 -3.0101E-02 -3.9723E-01 -6.8339E+01  8.2186E-01 -4.3326E+01
             8.7042E-01
 GRADIENT:   4.6515E-01 -6.1850E+00  6.7873E+00  1.2185E+00 -4.5374E+00  4.7826E-02  1.2157E-01  0.0000E+00 -3.5966E-01  0.0000E+00
             3.2459E+00

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1530.48619615564        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1861
 NPARAMETR:  1.0212E+00  2.0763E+00  4.9803E-02  2.6891E-01  9.1838E-01  8.7785E-01  6.0765E-01  1.0000E-02  2.0414E+00  1.0000E-02
             2.1607E+00
 PARAMETER:  1.2221E-01  8.3056E-01 -2.8989E+00 -1.2014E+00  1.3856E-02 -3.0101E-02 -3.9723E-01 -6.8339E+01  8.2186E-01 -4.3326E+01
             8.7042E-01
 GRADIENT:   2.3249E+04 -8.6284E+00  5.2434E+01  2.3529E+03 -1.4208E+04  3.3405E-02  1.0202E-01  0.0000E+00  3.4003E+03  0.0000E+00
            -8.1123E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1861
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1861E-04 -1.4922E-02  9.2921E-05  9.9480E-03 -2.9066E-04
 SE:             2.9299E-02  2.6043E-02  7.5096E-05  2.0942E-02  2.5430E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9405E-01  5.6666E-01  2.1595E-01  6.3477E-01  2.5306E-01

 ETASHRINKSD(%)  1.8433E+00  1.2753E+01  9.9748E+01  2.9842E+01  9.9148E+01
 ETASHRINKVR(%)  3.6526E+00  2.3879E+01  9.9999E+01  5.0778E+01  9.9993E+01
 EBVSHRINKSD(%)  2.0379E+00  1.3741E+01  9.9783E+01  2.5853E+01  9.9061E+01
 EBVSHRINKVR(%)  4.0344E+00  2.5594E+01  1.0000E+02  4.5022E+01  9.9991E+01
 RELATIVEINF(%)  9.3342E+01  1.4694E+01  2.0433E-04  1.2870E+01  1.6650E-03
 EPSSHRINKSD(%)  3.2217E+01
 EPSSHRINKVR(%)  5.4055E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1530.4861961556423     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -795.33536959190417     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1530.486       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.08E+00  4.98E-02  2.72E-01  9.17E-01  8.78E-01  6.08E-01  1.00E-02  2.06E+00  1.00E-02  2.16E+00
 


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
+        9.10E+06
 
 TH 2
+       -5.31E+01  4.83E+04
 
 TH 3
+       -3.97E+02  5.71E+05  6.83E+06
 
 TH 4
+        1.75E+06  1.66E+02 -3.24E+03  6.59E+05
 
 TH 5
+       -6.20E+06 -4.69E+02 -6.51E+02  9.99E+02  8.44E+06
 
 TH 6
+        1.34E+03 -7.18E+00  4.43E+01  4.90E+02  8.82E+06  2.36E+02
 
 TH 7
+        2.35E+06 -1.34E+00 -2.08E+01  3.56E+02 -3.20E+06 -3.35E+06  1.22E+06
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.42E+05  2.47E+04  3.02E+05  1.30E+05  6.40E+01  9.92E+01  6.82E+01  0.00E+00  4.97E+04
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.02E+05 -2.19E+04 -2.67E+05 -1.15E+05  4.11E+05  2.30E+00 -1.56E+05  0.00E+00 -4.46E+04  0.00E+00  4.02E+04
 
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
 #CPUT: Total CPU Time in Seconds,       29.759
Stop Time:
Wed Sep 29 12:55:07 CDT 2021
