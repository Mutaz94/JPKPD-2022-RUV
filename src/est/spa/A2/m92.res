Wed Sep 29 13:09:18 CDT 2021
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
$DATA ../../../../data/spa/A2/dat92.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -926.706332254307        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3931E+02 -6.1761E+00  6.6471E-01  4.4180E+00  1.3380E+02  5.0996E+01 -4.4805E+01 -1.8797E+01 -7.3192E+01 -6.2767E+01
            -1.2252E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1336.72910864979        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1304E+00  1.0343E+00  1.0266E+00  1.0727E+00  9.5914E-01  1.1299E+00  1.0725E+00  9.8041E-01  1.1105E+00  1.0282E+00
             2.3935E+00
 PARAMETER:  2.2255E-01  1.3372E-01  1.2626E-01  1.7013E-01  5.8283E-02  2.2211E-01  1.6997E-01  8.0215E-02  2.0482E-01  1.2786E-01
             9.7275E-01
 GRADIENT:   3.7025E+02  7.8074E+00 -2.6764E+00  2.2962E+01 -1.0499E+01  4.2781E+01  4.4717E+00  3.4191E+00  1.4857E+01  1.0791E+01
            -8.6574E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1358.43238408370        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0428E+00  7.8213E-01  1.3997E+00  1.2633E+00  1.0368E+00  9.2080E-01  9.0600E-01  1.5724E-01  1.0208E+00  1.2111E+00
             2.4956E+00
 PARAMETER:  1.4194E-01 -1.4574E-01  4.3626E-01  3.3371E-01  1.3610E-01  1.7492E-02  1.2788E-03 -1.7500E+00  1.2057E-01  2.9153E-01
             1.0145E+00
 GRADIENT:   2.1772E+02  1.6309E+01 -2.1216E+00  7.6735E+01 -2.7608E+00 -3.5503E+00  1.7230E+00  4.6750E-02  1.0519E+01  1.3925E+01
            -7.2767E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1373.58744581808        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  9.7371E-01  1.0476E+00  1.7017E+00  1.0553E+00  1.2012E+00  9.1194E-01  4.0764E-01  2.2022E-01  1.1019E+00  6.9418E-01
             3.0068E+00
 PARAMETER:  7.3359E-02  1.4652E-01  6.3162E-01  1.5380E-01  2.8331E-01  7.8242E-03 -7.9738E-01 -1.4131E+00  1.9702E-01 -2.6502E-01
             1.2009E+00
 GRADIENT:  -6.2107E+00  3.9172E-01 -2.7695E-01 -1.5765E+00 -4.3523E+00  5.3626E+00  5.1664E-01  9.2864E-02  7.9235E-01  4.9150E+00
             1.5436E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1375.39920849835        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.7735E-01  1.0415E+00  2.1032E+00  1.0616E+00  1.2745E+00  8.9073E-01  2.0956E-01  4.2940E-01  1.1227E+00  2.4481E-01
             3.0170E+00
 PARAMETER:  7.7092E-02  1.4066E-01  8.4348E-01  1.5977E-01  3.4256E-01 -1.5711E-02 -1.4628E+00 -7.4536E-01  2.1570E-01 -1.3073E+00
             1.2043E+00
 GRADIENT:   4.8849E+00  1.9918E+00  1.1209E+00  2.2139E+00 -2.5321E+00 -1.4049E+00  1.9984E-01  1.9516E-01  2.4909E+00  2.8159E-01
             2.2524E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1376.03196156000        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      422
 NPARAMETR:  9.8475E-01  9.6854E-01  2.0756E+00  1.1156E+00  1.2609E+00  9.0534E-01  1.3080E-01  1.7935E-01  1.0729E+00  1.1653E-01
             3.0544E+00
 PARAMETER:  8.4633E-02  6.8038E-02  8.3026E-01  2.0939E-01  3.3185E-01  5.5028E-04 -1.9341E+00 -1.6184E+00  1.7037E-01 -2.0496E+00
             1.2166E+00
 GRADIENT:  -1.9293E+01 -1.7629E+00  7.3513E-03 -5.8381E+00 -1.4510E+00  8.5246E-01 -1.6044E-02  3.9840E-02  1.1845E-01  4.8460E-02
            -3.6891E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1376.38227678534        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      597
 NPARAMETR:  9.9068E-01  7.2970E-01  2.5026E+00  1.2922E+00  1.2529E+00  9.0159E-01  5.3028E-02  5.3975E-02  9.2258E-01  2.2881E-02
             3.0548E+00
 PARAMETER:  9.0640E-02 -2.1513E-01  1.0173E+00  3.5635E-01  3.2543E-01 -3.5983E-03 -2.8369E+00 -2.8192E+00  1.9419E-02 -3.6774E+00
             1.2167E+00
 GRADIENT:  -2.8445E+00  4.4722E+00  1.2548E+00  7.2383E+00 -4.5466E+00 -4.0254E-01 -1.9423E-04  2.6471E-03 -1.1368E+00  1.9153E-03
            -2.0654E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1376.49396482699        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      773
 NPARAMETR:  9.9121E-01  5.8363E-01  2.4050E+00  1.3892E+00  1.2118E+00  9.0377E-01  2.1227E-02  1.0000E-02  8.6084E-01  1.0000E-02
             3.0605E+00
 PARAMETER:  9.1168E-02 -4.3849E-01  9.7754E-01  4.2876E-01  2.9210E-01 -1.1858E-03 -3.7525E+00 -4.7445E+00 -4.9845E-02 -4.9759E+00
             1.2186E+00
 GRADIENT:   9.5852E-01  2.5244E+00 -1.1410E-01  7.2148E+00 -1.0572E-01  4.5390E-01  1.5399E-04  0.0000E+00 -5.8431E-01  0.0000E+00
            -2.0993E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1376.54180750223        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      949
 NPARAMETR:  9.8894E-01  4.6061E-01  2.2489E+00  1.4626E+00  1.1475E+00  9.0128E-01  1.0000E-02  1.0000E-02  8.1818E-01  1.0000E-02
             3.0573E+00
 PARAMETER:  8.8883E-02 -6.7520E-01  9.1046E-01  4.8021E-01  2.3762E-01 -3.9381E-03 -4.7620E+00 -6.7687E+00 -1.0068E-01 -6.3488E+00
             1.2175E+00
 GRADIENT:  -7.6979E-01  1.4338E+00  2.5478E-01  4.0708E+00 -2.3152E+00 -2.4816E-01  0.0000E+00  0.0000E+00 -4.4973E-01  0.0000E+00
            -5.4188E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1376.57185816538        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1125
 NPARAMETR:  9.8897E-01  3.6814E-01  2.3834E+00  1.5263E+00  1.1474E+00  9.0208E-01  1.0000E-02  1.0000E-02  7.8146E-01  1.0000E-02
             3.0583E+00
 PARAMETER:  8.8904E-02 -8.9930E-01  9.6855E-01  5.2283E-01  2.3754E-01 -3.0562E-03 -5.8372E+00 -8.3878E+00 -1.4659E-01 -7.7502E+00
             1.2179E+00
 GRADIENT:   6.9485E-01  9.5225E-01  2.2113E-01  3.7689E+00 -1.1303E+00  1.9274E-01  0.0000E+00  0.0000E+00 -3.6049E-01  0.0000E+00
            -2.6807E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1376.58056767182        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1300
 NPARAMETR:  9.8813E-01  3.0889E-01  2.4323E+00  1.5658E+00  1.1414E+00  9.0130E-01  1.0000E-02  1.0000E-02  7.6068E-01  1.0000E-02
             3.0585E+00
 PARAMETER:  8.8055E-02 -1.0748E+00  9.8883E-01  5.4843E-01  2.3226E-01 -3.9135E-03 -6.7150E+00 -9.7816E+00 -1.7354E-01 -8.8557E+00
             1.2179E+00
 GRADIENT:  -3.0781E-01  6.8610E-01  6.7305E-02  3.4072E+00 -3.0990E-01 -1.4563E-02  0.0000E+00  0.0000E+00 -2.9100E-01  0.0000E+00
            -1.8160E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1376.59151827902        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1475
 NPARAMETR:  9.8717E-01  2.3359E-01  2.4394E+00  1.6135E+00  1.1241E+00  9.0061E-01  1.0000E-02  1.0000E-02  7.3752E-01  1.0000E-02
             3.0583E+00
 PARAMETER:  8.7088E-02 -1.3542E+00  9.9173E-01  5.7840E-01  2.1699E-01 -4.6827E-03 -8.1657E+00 -1.2142E+01 -2.0446E-01 -1.0637E+01
             1.2179E+00
 GRADIENT:  -7.2822E-01  2.5790E-01 -7.5647E-02  1.1635E+00  3.9857E-01 -1.2948E-01  0.0000E+00  0.0000E+00 -8.7764E-02  0.0000E+00
            -3.2280E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1376.59201763867        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1650
 NPARAMETR:  9.8686E-01  2.0302E-01  2.4336E+00  1.6341E+00  1.1148E+00  9.0049E-01  1.0000E-02  1.0000E-02  7.2798E-01  1.0000E-02
             3.0582E+00
 PARAMETER:  8.6771E-02 -1.4945E+00  9.8937E-01  5.9110E-01  2.0867E-01 -4.8122E-03 -8.9068E+00 -1.3349E+01 -2.1748E-01 -1.1538E+01
             1.2178E+00
 GRADIENT:  -8.8627E-01  4.6274E-01  3.1273E-02  3.5047E+00 -2.7817E-01 -1.5234E-01  0.0000E+00  0.0000E+00 -2.7579E-01  0.0000E+00
            -1.7795E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1376.59224756307        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1825
 NPARAMETR:  9.8666E-01  1.8023E-01  2.4184E+00  1.6489E+00  1.1059E+00  9.0049E-01  1.0000E-02  1.0000E-02  7.2140E-01  1.0000E-02
             3.0581E+00
 PARAMETER:  8.6571E-02 -1.6135E+00  9.8312E-01  6.0008E-01  2.0068E-01 -4.8151E-03 -9.5446E+00 -1.4395E+01 -2.2656E-01 -1.2306E+01
             1.2178E+00
 GRADIENT:  -7.8160E-01  5.3328E-01  1.2106E-01  4.8130E+00 -8.1129E-01 -1.2788E-01  0.0000E+00  0.0000E+00 -3.8102E-01  0.0000E+00
            -2.6611E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1376.59382698633        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2000
 NPARAMETR:  9.8645E-01  1.5616E-01  2.3899E+00  1.6635E+00  1.0944E+00  9.0053E-01  1.0000E-02  1.0000E-02  7.1511E-01  1.0000E-02
             3.0580E+00
 PARAMETER:  8.6361E-02 -1.7569E+00  9.7124E-01  6.0891E-01  1.9020E-01 -4.7706E-03 -1.0323E+01 -1.5677E+01 -2.3532E-01 -1.3235E+01
             1.2178E+00
 GRADIENT:  -4.7343E-01  5.0372E-01  1.9881E-01  5.1980E+00 -1.2621E+00 -6.9485E-02  0.0000E+00  0.0000E+00 -4.1519E-01  0.0000E+00
            -3.0907E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1376.59746203296        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2175
 NPARAMETR:  9.8616E-01  1.2700E-01  2.3545E+00  1.6803E+00  1.0805E+00  9.0053E-01  1.0000E-02  1.0000E-02  7.0802E-01  1.0000E-02
             3.0578E+00
 PARAMETER:  8.6063E-02 -1.9636E+00  9.5632E-01  6.1897E-01  1.7742E-01 -4.7769E-03 -1.1464E+01 -1.7536E+01 -2.4528E-01 -1.4585E+01
             1.2177E+00
 GRADIENT:  -4.7316E-02  3.5354E-01  2.2724E-01  4.0213E+00 -1.4112E+00  3.5869E-03  0.0000E+00  0.0000E+00 -3.2928E-01  0.0000E+00
            -2.6745E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1376.60062328232        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2350
 NPARAMETR:  9.8581E-01  9.5369E-02  2.3493E+00  1.6996E+00  1.0718E+00  9.0038E-01  1.0000E-02  1.0000E-02  6.9977E-01  1.0000E-02
             3.0577E+00
 PARAMETER:  8.5708E-02 -2.2500E+00  9.5414E-01  6.3037E-01  1.6938E-01 -4.9394E-03 -1.3076E+01 -2.0087E+01 -2.5701E-01 -1.6476E+01
             1.2177E+00
 GRADIENT:   1.1406E-01  2.0666E-01  1.9321E-01  2.5976E+00 -1.1746E+00  2.8233E-02  0.0000E+00  0.0000E+00 -2.1759E-01  0.0000E+00
            -1.8996E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1376.60265445120        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2525
 NPARAMETR:  9.8550E-01  6.4256E-02  2.3762E+00  1.7198E+00  1.0695E+00  9.0018E-01  1.0000E-02  1.0000E-02  6.9103E-01  1.0000E-02
             3.0578E+00
 PARAMETER:  8.5398E-02 -2.6449E+00  9.6552E-01  6.4222E-01  1.6722E-01 -5.1574E-03 -1.5338E+01 -2.3600E+01 -2.6957E-01 -1.9109E+01
             1.2177E+00
 GRADIENT:   9.3991E-02  1.0572E-01  1.1644E-01  1.5345E+00 -6.8551E-01  2.3906E-02  0.0000E+00  0.0000E+00 -1.2688E-01  0.0000E+00
            -1.0918E-01

0ITERATION NO.:   88    OBJECTIVE VALUE:  -1376.60400133142        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:     2636
 NPARAMETR:  9.8535E-01  6.3866E-02  2.3806E+00  1.7196E+00  1.0685E+00  8.9997E-01  1.0000E-02  1.0000E-02  6.9066E-01  1.0000E-02
             3.0572E+00
 PARAMETER:  8.5408E-02 -2.6749E+00  9.6487E-01  6.4171E-01  1.6789E-01 -5.1801E-03 -1.5338E+01 -2.3600E+01 -2.6924E-01 -1.9109E+01
             1.2178E+00
 GRADIENT:   3.7567E-01 -2.4608E-02 -1.9477E-01 -1.2007E+00  4.5552E-01  6.4282E-02  0.0000E+00  0.0000E+00  1.5762E-01  0.0000E+00
             1.9067E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2636
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7412E-04 -1.8921E-05  6.0433E-05 -9.5959E-03 -1.9732E-05
 SE:             2.8777E-02  1.0061E-05  5.7240E-05  2.5632E-02  1.5541E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9517E-01  6.0024E-02  2.9107E-01  7.0813E-01  8.9896E-01

 ETASHRINKSD(%)  3.5936E+00  9.9966E+01  9.9808E+01  1.4129E+01  9.9479E+01
 ETASHRINKVR(%)  7.0581E+00  1.0000E+02  1.0000E+02  2.6262E+01  9.9997E+01
 EBVSHRINKSD(%)  3.5556E+00  9.9967E+01  9.9788E+01  1.3901E+01  9.9461E+01
 EBVSHRINKVR(%)  6.9847E+00  1.0000E+02  1.0000E+02  2.5870E+01  9.9997E+01
 RELATIVEINF(%)  8.6640E+01  2.6410E-07  2.7854E-05  2.4025E+00  1.2694E-04
 EPSSHRINKSD(%)  2.2130E+01
 EPSSHRINKVR(%)  3.9363E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1376.6040013314223     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -641.45317476768412     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.89
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1376.604       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  6.24E-02  2.37E+00  1.72E+00  1.07E+00  9.00E-01  1.00E-02  1.00E-02  6.91E-01  1.00E-02  3.06E+00
 


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
+        1.32E+03
 
 TH 2
+       -6.92E+01  3.06E+02
 
 TH 3
+       -4.99E-01  1.11E+01  8.02E+00
 
 TH 4
+       -8.26E+01  4.01E+02  1.97E+00  5.78E+02
 
 TH 5
+       -2.29E+00 -1.48E+02 -4.33E+01 -1.21E+02  2.74E+02
 
 TH 6
+       -1.91E+00 -1.08E+01  1.28E+00 -1.78E+01 -4.21E+00  2.11E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.49E+00 -7.20E+01  2.70E+00 -1.59E+01  1.58E+01 -1.93E+00  0.00E+00  0.00E+00  2.19E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.77E+01 -1.18E+01  4.16E-01 -1.05E+01  1.56E+00  5.19E+00  0.00E+00  0.00E+00  1.84E+01  0.00E+00  4.49E+01
 
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
 #CPUT: Total CPU Time in Seconds,       36.648
Stop Time:
Wed Sep 29 13:09:57 CDT 2021
