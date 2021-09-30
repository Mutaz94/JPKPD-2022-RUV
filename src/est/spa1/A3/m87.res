Thu Sep 30 00:38:29 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat87.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   31.1570003336899        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2724E+02 -2.2471E+01  1.7920E+02 -6.1193E+01  2.8462E+02  3.6041E+01 -6.0563E+01 -3.0636E+02 -1.1922E+02 -1.4314E+02
            -3.5971E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1473.49463630277        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  1.0028E+00  8.8151E-01  2.3384E+00  1.3080E+00  1.6342E+00  9.1121E-01  1.1731E+00  6.3150E-02  1.0349E+00  1.1317E+00
             2.8623E+00
 PARAMETER:  1.0279E-01 -2.6121E-02  9.4946E-01  3.6847E-01  5.9114E-01  7.0223E-03  2.5961E-01 -2.6622E+00  1.3432E-01  2.2376E-01
             1.1516E+00
 GRADIENT:   6.1270E+01  1.3484E+01 -1.9808E+01  1.1036E+02  9.5782E+01 -1.0631E+01  9.4263E+00 -2.0181E-03  3.2284E+01 -1.0785E+01
            -3.2382E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1506.51210103358        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      168
 NPARAMETR:  1.0041E+00  1.6313E+00  9.8475E-01  8.9587E-01  1.3509E+00  9.1666E-01  1.0092E+00  1.5056E-02  1.0902E+00  6.7262E-01
             3.2695E+00
 PARAMETER:  1.0414E-01  5.8940E-01  8.4631E-02 -9.9562E-03  4.0075E-01  1.2978E-02  1.0916E-01 -4.0960E+00  1.8635E-01 -2.9657E-01
             1.2846E+00
 GRADIENT:   1.2579E+01  1.5864E+02 -2.6553E+01  1.3071E+02  5.9892E+01 -1.0695E+01  2.1729E+01  2.3604E-04 -4.9691E+00  2.0764E+00
            -1.3044E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1536.34131841443        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      240
 NPARAMETR:  9.9756E-01  1.1667E+00  6.8754E-01  9.8260E-01  8.1885E-01  9.3667E-01  7.3689E-01  1.0000E-02  1.0217E+00  2.6237E-01
             3.5202E+00
 PARAMETER:  9.7558E-02  2.5419E-01 -2.7463E-01  8.2450E-02 -9.9859E-02  3.4573E-02 -2.0532E-01 -4.7939E+00  1.2145E-01 -1.2380E+00
             1.3585E+00
 GRADIENT:  -9.1354E+00 -3.1466E+00  4.3340E+00 -2.2445E+01 -1.7398E+00  8.5335E-01 -2.0136E+00  0.0000E+00 -3.5747E+00  7.2826E-01
             5.5080E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1537.90281940886        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0081E+00  9.8323E-01  4.9802E-01  1.1019E+00  6.1133E-01  9.6613E-01  8.3514E-01  1.0000E-02  9.7374E-01  2.0022E-01
             3.5204E+00
 PARAMETER:  1.0810E-01  8.3088E-02 -5.9711E-01  1.9702E-01 -3.9211E-01  6.5544E-02 -8.0150E-02 -5.0134E+00  7.3393E-02 -1.5083E+00
             1.3586E+00
 GRADIENT:  -2.4997E+01  2.3654E+01  5.5337E+00  3.0064E+01 -1.0811E+01  4.2782E+00 -3.2137E+00  0.0000E+00 -7.3544E+00  9.5418E-02
             9.0808E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1541.56521372525        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      559
 NPARAMETR:  1.0056E+00  7.0182E-01  3.4558E-01  1.1765E+00  4.1418E-01  9.4463E-01  1.0240E+00  1.0000E-02  1.0479E+00  1.3493E-01
             3.3227E+00
 PARAMETER:  1.0557E-01 -2.5407E-01 -9.6253E-01  2.6254E-01 -7.8144E-01  4.3038E-02  1.2370E-01 -4.9451E+00  1.4677E-01 -1.9030E+00
             1.3008E+00
 GRADIENT:  -2.5016E+01  1.3421E+01  1.5408E+01  2.8514E+01 -1.7523E+01 -8.6207E+00 -2.6532E+00  0.0000E+00 -3.3794E+00 -1.3415E+00
            -1.1517E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1543.46752905843        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  1.0173E+00  6.9794E-01  2.8262E-01  1.1218E+00  3.7456E-01  9.7839E-01  1.0539E+00  1.0000E-02  1.1179E+00  1.2661E-01
             3.2683E+00
 PARAMETER:  1.1712E-01 -2.5962E-01 -1.1636E+00  2.1491E-01 -8.8200E-01  7.8151E-02  1.5252E-01 -4.9851E+00  2.1148E-01 -1.9666E+00
             1.2843E+00
 GRADIENT:   4.3000E+00  5.5468E-02  4.0798E-01 -2.3391E-01  2.7309E+00  7.4237E-01  5.5295E-01  0.0000E+00  7.5939E-01 -1.3072E+00
            -3.9439E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1546.93541610679        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  9.8606E-01  5.7042E-01  2.4034E-01  1.1679E+00  3.0564E-01  9.7980E-01  9.8525E-01  1.0000E-02  1.1780E+00  5.2658E-01
             3.0737E+00
 PARAMETER:  8.5958E-02 -4.6138E-01 -1.3257E+00  2.5525E-01 -1.0853E+00  7.9594E-02  8.5139E-02 -5.2760E+00  2.6379E-01 -5.4136E-01
             1.2229E+00
 GRADIENT:  -6.3131E+01  2.9133E+01  8.2067E+00  4.0492E+01 -3.4361E+01 -5.9487E+00  5.7249E+00  0.0000E+00 -2.3354E+00  2.2327E+00
            -1.0427E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1549.97727126964        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  1.0115E+00  4.9051E-01  2.1018E-01  1.1250E+00  2.7268E-01  9.9202E-01  7.2130E-01  1.0000E-02  1.2341E+00  5.8037E-01
             3.0921E+00
 PARAMETER:  1.1147E-01 -6.1231E-01 -1.4598E+00  2.1778E-01 -1.1995E+00  9.1992E-02 -2.2670E-01 -5.5003E+00  3.1036E-01 -4.4409E-01
             1.2288E+00
 GRADIENT:  -5.4807E+00 -4.3493E+00 -1.0050E+01  1.6125E+00  1.5926E+01 -7.5278E-01  1.4517E+00  0.0000E+00  6.6221E-01  7.1607E-01
             6.2495E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1550.41010280619        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1265
 NPARAMETR:  1.0124E+00  4.4463E-01  1.8559E-01  1.1000E+00  2.4583E-01  9.9716E-01  3.8757E-01  1.0000E-02  1.2907E+00  6.4960E-01
             3.0339E+00
 PARAMETER:  1.1232E-01 -7.1052E-01 -1.5842E+00  1.9528E-01 -1.3031E+00  9.7157E-02 -8.4787E-01 -5.4834E+00  3.5517E-01 -3.3140E-01
             1.2098E+00
 GRADIENT:  -2.5044E+00 -5.2262E-01 -1.8217E+00  1.8292E+00  2.1779E+00 -5.7728E-01  3.1719E-01  0.0000E+00  4.3506E-01  1.6510E+00
            -1.2850E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1550.49978873843        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1440
 NPARAMETR:  1.0137E+00  4.3339E-01  1.7944E-01  1.0900E+00  2.3946E-01  9.9946E-01  1.5624E-01  1.0000E-02  1.3040E+00  6.5914E-01
             3.0288E+00
 PARAMETER:  1.1356E-01 -7.3612E-01 -1.6179E+00  1.8616E-01 -1.3294E+00  9.9459E-02 -1.7563E+00 -5.4003E+00  3.6545E-01 -3.1683E-01
             1.2082E+00
 GRADIENT:   2.3941E-01  6.2196E-02 -1.7756E-01 -1.6626E-01 -1.5360E-01 -6.3770E-02  1.8148E-02  0.0000E+00 -2.5186E-01 -1.5064E-01
            -2.4150E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1550.51072059813        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1623             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0136E+00  4.3107E-01  1.7807E-01  1.0882E+00  2.3807E-01  9.9986E-01  1.5623E-02  1.0000E-02  1.3092E+00  6.6326E-01
             3.0271E+00
 PARAMETER:  1.1347E-01 -7.4148E-01 -1.6256E+00  1.8450E-01 -1.3352E+00  9.9862E-02 -4.0590E+00 -5.2158E+00  3.6939E-01 -3.1059E-01
             1.2076E+00
 GRADIENT:   3.7299E+01  8.7361E+00  2.1247E+01  1.2782E+01  8.4950E+01  2.7738E+00  1.2806E-03  0.0000E+00  5.3379E+00  7.1034E-01
             1.1933E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1550.51083364046        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1792
 NPARAMETR:  1.0135E+00  4.3098E-01  1.7806E-01  1.0882E+00  2.3806E-01  9.9988E-01  1.0000E-02  1.0000E-02  1.3092E+00  6.6340E-01
             3.0272E+00
 PARAMETER:  1.1346E-01 -7.4169E-01 -1.6256E+00  1.8453E-01 -1.3352E+00  9.9884E-02 -4.8821E+00 -5.2158E+00  3.6945E-01 -3.1038E-01
             1.2076E+00
 GRADIENT:   8.8998E-02  2.0246E-01  2.0444E-01 -3.5127E-02 -6.2053E-01  6.0679E-03  0.0000E+00  0.0000E+00  4.8303E-02 -3.2729E-03
            -6.2834E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1792
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6352E-03 -1.1573E-04  1.9202E-04 -1.2170E-02  1.6210E-03
 SE:             2.8870E-02  1.5366E-04  2.0952E-04  2.6452E-02  2.3542E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5483E-01  4.5134E-01  3.5941E-01  6.4546E-01  9.4511E-01

 ETASHRINKSD(%)  3.2828E+00  9.9485E+01  9.9298E+01  1.1383E+01  2.1131E+01
 ETASHRINKVR(%)  6.4579E+00  9.9997E+01  9.9995E+01  2.1470E+01  3.7797E+01
 EBVSHRINKSD(%)  2.8758E+00  9.9491E+01  9.9347E+01  9.3186E+00  2.1628E+01
 EBVSHRINKVR(%)  5.6690E+00  9.9997E+01  9.9996E+01  1.7769E+01  3.8579E+01
 RELATIVEINF(%)  9.4124E+01  3.0685E-04  3.3720E-04  4.9074E+01  2.1454E+00
 EPSSHRINKSD(%)  2.5448E+01
 EPSSHRINKVR(%)  4.4420E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1550.5108336404569     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -631.57230043578420     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1550.511       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  4.31E-01  1.78E-01  1.09E+00  2.38E-01  1.00E+00  1.00E-02  1.00E-02  1.31E+00  6.63E-01  3.03E+00
 


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
+       -3.83E+01  1.55E+03
 
 TH 3
+       -9.45E+01  2.76E+03  1.42E+04
 
 TH 4
+       -1.22E+01  1.05E+02 -4.87E+02  4.54E+02
 
 TH 5
+        1.88E+02 -5.50E+03 -1.72E+04 -2.95E+02  2.74E+04
 
 TH 6
+        4.49E+00 -1.59E+01  7.46E+01 -1.22E+01 -9.06E+00  1.76E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.03E+01 -3.55E+01  1.78E+02 -1.10E+01  4.73E+01 -2.98E+00  0.00E+00  0.00E+00  7.69E+01
 
 TH10
+       -4.87E+00 -2.07E+01  6.65E+01  7.26E+00  1.04E+02  4.53E+00  0.00E+00  0.00E+00  5.58E+00  1.74E+02
 
 TH11
+       -1.77E+01 -1.26E+01 -7.61E+01 -5.48E+00  4.79E+01  3.03E+00  0.00E+00  0.00E+00  5.06E+00  1.73E+01  5.23E+01
 
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
 #CPUT: Total CPU Time in Seconds,       34.821
Stop Time:
Thu Sep 30 00:39:05 CDT 2021
