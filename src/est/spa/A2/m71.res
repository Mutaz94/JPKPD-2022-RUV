Sat Sep 25 08:49:32 CDT 2021
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
$DATA ../../../../data/spa/A2/dat71.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1345.91556421876        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5345E+02 -3.4686E+01  9.6717E+01 -1.4674E+02 -5.3739E+01  1.6204E+01 -1.6746E+01 -3.1348E+01 -1.2519E+01 -2.0223E+01
            -4.7326E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1473.95233517891        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.3036E-01  1.0510E+00  9.2582E-01  1.0798E+00  1.0486E+00  9.2282E-01  1.0492E+00  1.0186E+00  9.4056E-01  9.9371E-01
             1.6872E+00
 PARAMETER:  2.7811E-02  1.4973E-01  2.2928E-02  1.7682E-01  1.4746E-01  1.9676E-02  1.4800E-01  1.1840E-01  3.8719E-02  9.3687E-02
             6.2308E-01
 GRADIENT:  -5.1576E+01  5.7709E+00 -5.5806E+00  2.7079E+01  1.6997E+01 -6.8233E+00 -6.7488E+00 -1.0484E+00  1.6156E+00 -2.0505E+00
            -3.5140E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1475.47799723002        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.3576E-01  1.0318E+00  7.4257E-01  1.0769E+00  8.8235E-01  9.4137E-01  1.2408E+00  8.4545E-01  8.6472E-01  7.2352E-01
             1.7016E+00
 PARAMETER:  3.3604E-02  1.3133E-01 -1.9763E-01  1.7406E-01 -2.5163E-02  3.9582E-02  3.1573E-01 -6.7887E-02 -4.5347E-02 -2.2363E-01
             6.3159E-01
 GRADIENT:  -3.9913E+01  1.9565E+01  8.4517E+00  2.6900E+01 -2.9274E+00  5.4094E-01  1.8325E+00 -1.6062E+00 -4.5411E-02 -5.4704E+00
            -3.0069E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1479.94323265975        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.6633E-01  1.1084E+00  2.8066E-01  9.2221E-01  5.4976E-01  9.4013E-01  1.0962E+00  1.7811E-01  8.0403E-01  4.5654E-01
             1.8314E+00
 PARAMETER:  6.5751E-02  2.0289E-01 -1.1706E+00  1.9020E-02 -4.9827E-01  3.8259E-02  1.9182E-01 -1.6253E+00 -1.1812E-01 -6.8408E-01
             7.0506E-01
 GRADIENT:   4.5729E+01  7.0163E+01  3.5798E+01 -7.4990E-02 -8.0548E+01  1.8947E+00  5.5142E+00  1.8319E-01  7.8441E-01  6.1815E+00
             2.7634E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1487.96594549517        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.3173E-01  1.1531E+00  1.9465E-01  8.1476E-01  5.3827E-01  9.3165E-01  1.0193E+00  3.1897E-01  8.6169E-01  1.8163E-01
             1.6662E+00
 PARAMETER:  2.9284E-02  2.4247E-01 -1.5366E+00 -1.0486E-01 -5.1940E-01  2.9201E-02  1.1909E-01 -1.0427E+00 -4.8858E-02 -1.6058E+00
             6.1057E-01
 GRADIENT:  -9.4893E+00 -1.5453E+01 -5.5827E+00  5.7166E+00  1.6005E+01 -1.9700E-01  7.2890E+00 -9.6202E-02 -3.5053E+00  1.0378E+00
            -4.8709E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1494.36127149489        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.3363E-01  1.5717E+00  1.4623E-01  6.0025E-01  7.0350E-01  9.2252E-01  8.3457E-01  1.2505E+00  1.0596E+00  5.7613E-02
             1.6605E+00
 PARAMETER:  3.1321E-02  5.5215E-01 -1.8226E+00 -4.1041E-01 -2.5168E-01  1.9357E-02 -8.0837E-02  3.2354E-01  1.5793E-01 -2.7540E+00
             6.0713E-01
 GRADIENT:   2.5574E+00  4.4234E+01  5.6286E+00  1.3237E+01 -2.5900E+01 -3.5637E-01 -4.3476E+00 -1.5231E+00  8.0453E-01 -1.0159E-02
             9.5328E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1495.73061514289        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  9.3200E-01  1.6224E+00  1.3972E-01  5.5514E-01  7.4589E-01  9.2259E-01  8.3316E-01  1.5697E+00  1.1213E+00  4.5007E-02
             1.5816E+00
 PARAMETER:  2.9577E-02  5.8388E-01 -1.8681E+00 -4.8853E-01 -1.9317E-01  1.9426E-02 -8.2532E-02  5.5090E-01  2.1447E-01 -3.0009E+00
             5.5845E-01
 GRADIENT:   1.0999E+00  3.5994E+00  1.8851E+00 -9.9633E-01 -4.3938E+00  1.5506E-01  1.2665E-02 -5.4121E-01  2.3950E-01 -3.1670E-02
            -1.0162E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1496.40156351892        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.3864E-01  1.7600E+00  1.2724E-01  4.8978E-01  8.2461E-01  9.2486E-01  8.0558E-01  1.8262E+00  1.2244E+00  3.1672E-02
             1.6082E+00
 PARAMETER:  3.6674E-02  6.6533E-01 -1.9617E+00 -6.1380E-01 -9.2840E-02  2.1886E-02 -1.1619E-01  7.0221E-01  3.0245E-01 -3.3523E+00
             5.7509E-01
 GRADIENT:   6.3366E-01 -5.0805E-01 -4.8399E-01  8.0266E-01  7.3015E-01  1.8826E-01  7.7680E-02 -1.2502E-02 -3.5967E-02 -3.0126E-02
             7.5084E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1496.41366545195        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      784            RESET HESSIAN, TYPE II
 NPARAMETR:  9.3840E-01  1.7613E+00  1.2671E-01  4.8936E-01  8.2637E-01  9.2509E-01  8.0756E-01  1.8580E+00  1.2255E+00  5.2358E-02
             1.6093E+00
 PARAMETER:  3.6418E-02  6.6607E-01 -1.9658E+00 -6.1466E-01 -9.0714E-02  2.2132E-02 -1.1374E-01  7.1952E-01  3.0338E-01 -2.8496E+00
             5.7579E-01
 GRADIENT:   1.2854E+01  1.7749E+01 -7.7664E-01  6.0814E+00  2.8448E+00  1.4156E+00  1.2822E+00  8.2211E-01  4.1144E-01 -7.6231E-02
             1.5696E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1497.71396818435        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      859
 NPARAMETR:  9.4060E-01  1.7711E+00  1.2515E-01  4.8907E-01  8.2852E-01  9.2741E-01  8.0558E-01  1.8145E+00  1.2118E+00  4.2767E-01
             1.5857E+00
 PARAMETER:  3.8758E-02  6.7161E-01 -1.9783E+00 -6.1525E-01 -8.8110E-02  2.4644E-02 -1.1620E-01  6.9583E-01  2.9211E-01 -7.4940E-01
             5.6102E-01
 GRADIENT:   1.8090E+01  2.1653E+01 -1.1735E-01  7.9235E+00 -1.1932E+01  2.4946E+00  6.3424E+00  2.5246E+00  2.0304E+00  4.1376E-01
             1.2876E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1498.57639139965        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      931
 NPARAMETR:  9.3069E-01  1.8345E+00  1.1460E-01  4.4078E-01  8.8162E-01  9.1926E-01  7.6745E-01  1.8322E+00  1.2915E+00  5.2560E-01
             1.5211E+00
 PARAMETER:  2.8174E-02  7.0675E-01 -2.0663E+00 -7.1921E-01 -2.5995E-02  1.5811E-02 -1.6468E-01  7.0550E-01  3.5580E-01 -5.4321E-01
             5.1944E-01
 GRADIENT:  -8.7621E+00 -6.8989E+00  1.8082E+00 -3.6806E+00  3.9654E+00 -1.1109E+00 -6.8414E-01  7.4778E-01 -1.1139E+00 -1.6248E-01
            -3.1021E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1498.76130980092        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     1047
 NPARAMETR:  9.3102E-01  1.8618E+00  1.0473E-01  4.2209E-01  8.8846E-01  9.1980E-01  7.5625E-01  1.6290E+00  1.3543E+00  5.5812E-01
             1.5167E+00
 PARAMETER:  2.8522E-02  7.2152E-01 -2.1564E+00 -7.6254E-01 -1.8271E-02  1.6397E-02 -1.7938E-01  5.8795E-01  4.0331E-01 -4.8318E-01
             5.1652E-01
 GRADIENT:  -2.0619E+01 -2.8793E+01  1.1128E+00 -8.8145E+00 -3.8061E+00 -2.2639E+00 -1.8597E+00 -1.1789E+00 -4.9788E-01  5.4141E-02
            -2.2415E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1500.23265380254        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1230
 NPARAMETR:  9.2288E-01  2.1749E+00  5.3485E-02  2.6230E-01  1.1058E+00  9.1932E-01  7.1019E-01  1.0904E+00  1.9557E+00  8.9460E-01
             1.5419E+00
 PARAMETER:  1.9739E-02  8.7696E-01 -2.8284E+00 -1.2383E+00  2.0058E-01  1.5878E-02 -2.4222E-01  1.8653E-01  7.7076E-01 -1.1383E-02
             5.3300E-01
 GRADIENT:  -5.1404E+01  6.0989E-01  1.4603E+00  3.2103E+00 -6.9227E+00 -4.3057E+00 -1.4748E-01 -4.7842E-01 -1.7604E+00  6.0505E+00
             5.6330E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1502.76590756191        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1406
 NPARAMETR:  9.4187E-01  2.3432E+00  3.4751E-02  1.9479E-01  1.2403E+00  9.2994E-01  6.8456E-01  8.0827E-01  2.5955E+00  9.4076E-01
             1.5410E+00
 PARAMETER:  4.0113E-02  9.5152E-01 -3.2596E+00 -1.5358E+00  3.1535E-01  2.7367E-02 -2.7898E-01 -1.1286E-01  1.0538E+00  3.8934E-02
             5.3240E-01
 GRADIENT:  -1.6113E+00  6.6374E+01  1.2193E+00  6.6410E+00 -4.7671E+00  9.8828E-01 -5.3268E+00  8.1221E-01 -8.3797E+00 -4.3724E+00
            -8.1145E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1502.77478418320        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1602
 NPARAMETR:  9.4192E-01  2.3436E+00  3.4712E-02  1.9465E-01  1.2407E+00  9.2996E-01  6.8454E-01  8.0770E-01  2.5983E+00  9.4126E-01
             1.5408E+00
 PARAMETER:  4.0169E-02  9.5168E-01 -3.2607E+00 -1.5365E+00  3.1565E-01  2.7390E-02 -2.7901E-01 -1.1356E-01  1.0549E+00  3.9465E-02
             5.3232E-01
 GRADIENT:  -1.4665E+00  6.6468E+01  1.2223E+00  6.6562E+00 -4.7697E+00  9.9918E-01 -5.3212E+00  8.1597E-01 -8.3745E+00 -4.3691E+00
            -8.0889E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1503.42064105676        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1702
 NPARAMETR:  9.4197E-01  2.3406E+00  3.4770E-02  1.9449E-01  1.2509E+00  9.2733E-01  6.9077E-01  2.7405E-01  2.5970E+00  9.9830E-01
             1.5413E+00
 PARAMETER:  4.0220E-02  9.5043E-01 -3.2590E+00 -1.5374E+00  3.2386E-01  2.4550E-02 -2.6994E-01 -1.1944E+00  1.0543E+00  9.8300E-02
             5.3260E-01
 GRADIENT:   1.2341E+01  1.0683E+02  5.5418E+00  8.5613E+00  1.2740E+00  9.4594E-01 -1.6845E+00  1.3271E-01 -9.1030E+00 -3.4593E-01
             2.9429E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1503.96025997983        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1774
 NPARAMETR:  9.4200E-01  2.2973E+00  3.4791E-02  1.9411E-01  1.2293E+00  9.3636E-01  6.9223E-01  1.0000E-02  2.5976E+00  9.8003E-01
             1.5414E+00
 PARAMETER:  4.0251E-02  9.3173E-01 -3.2584E+00 -1.5393E+00  3.0646E-01  3.4249E-02 -2.6783E-01 -2.4516E+01  1.0546E+00  7.9823E-02
             5.3267E-01
 GRADIENT:   1.5695E+01  4.5427E+01  3.4632E+00  2.3616E+00 -2.8285E+00  5.0985E+00 -8.9739E-01  0.0000E+00 -7.7857E+00  1.5030E-01
             4.4102E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1504.28540695648        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1937
 NPARAMETR:  9.4179E-01  2.3006E+00  3.4492E-02  1.9339E-01  1.2353E+00  9.2702E-01  6.9608E-01  1.0000E-02  2.6867E+00  9.8432E-01
             1.5347E+00
 PARAMETER:  4.0026E-02  9.3315E-01 -3.2670E+00 -1.5430E+00  3.1127E-01  2.4218E-02 -2.6229E-01 -2.2501E+01  1.0883E+00  8.4195E-02
             5.2836E-01
 GRADIENT:   1.5823E+01  4.5325E+01  4.8667E+00  4.0185E+00 -9.2163E-01  1.2854E+00  1.0971E+00  0.0000E+00 -5.4876E+00 -4.5160E-02
             4.3942E+00

0ITERATION NO.:   88    OBJECTIVE VALUE:  -1504.32898208367        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     2028
 NPARAMETR:  9.4170E-01  2.3006E+00  3.4393E-02  1.9315E-01  1.2355E+00  9.2698E-01  6.9593E-01  1.0000E-02  2.6991E+00  9.8458E-01
             1.5334E+00
 PARAMETER:  3.9928E-02  9.3317E-01 -3.2699E+00 -1.5443E+00  3.1145E-01  2.4177E-02 -2.6246E-01 -2.2501E+01  1.0929E+00  8.4359E-02
             5.2746E-01
 GRADIENT:  -6.7395E+04  7.2182E+03 -2.0480E+03  2.1819E+03 -1.9097E+00  2.6958E-02  3.5365E-01  0.0000E+00  6.1442E+03 -1.3400E-01
            -1.2772E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2028
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.4835E-03 -2.2779E-02  1.9203E-04  1.7961E-02 -3.9627E-02
 SE:             2.9524E-02  2.7473E-02  6.6455E-05  1.9228E-02  2.0003E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5992E-01  4.0701E-01  3.8577E-03  3.5026E-01  4.7581E-02

 ETASHRINKSD(%)  1.0918E+00  7.9633E+00  9.9777E+01  3.5582E+01  3.2989E+01
 ETASHRINKVR(%)  2.1717E+00  1.5293E+01  1.0000E+02  5.8503E+01  5.5095E+01
 EBVSHRINKSD(%)  1.1641E+00  1.0076E+01  9.9802E+01  3.0515E+01  3.3967E+01
 EBVSHRINKVR(%)  2.3147E+00  1.9138E+01  1.0000E+02  5.1719E+01  5.6397E+01
 RELATIVEINF(%)  9.6336E+01  4.0083E+01  2.6512E-04  2.3191E+01  1.9922E+01
 EPSSHRINKSD(%)  4.1891E+01
 EPSSHRINKVR(%)  6.6233E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1504.3289820836706     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -769.17815551993237     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1504.329       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.42E-01  2.30E+00  3.44E-02  1.93E-01  1.24E+00  9.27E-01  6.96E-01  1.00E-02  2.70E+00  9.84E-01  1.53E+00
 


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
+        1.90E+08
 
 TH 2
+        6.63E+02  3.66E+05
 
 TH 3
+       -1.36E+04 -2.80E+02  1.32E+08
 
 TH 4
+        5.02E+03 -6.05E+03  1.16E+05  1.89E+07
 
 TH 5
+        4.65E+07  1.49E+02 -4.59E+03  2.15E+03  3.34E+02
 
 TH 6
+       -6.21E+03  2.68E+02 -5.29E+03  1.95E+03 -6.09E+00  2.09E+02
 
 TH 7
+       -9.80E+07 -6.53E+01  1.44E+03 -4.93E+02  6.80E+00 -4.26E+00  2.78E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.06E+02 -1.46E+03  2.79E+04  1.91E+06  1.81E+02  1.98E+02 -4.55E+01  0.00E+00  1.93E+05
 
 TH10
+        1.82E+08 -6.82E+01  1.41E+03 -4.60E+02 -4.22E+01 -2.46E+00  3.44E+00  0.00E+00 -4.36E+01  4.79E+01
 
 TH11
+        2.21E+07  2.91E+02 -5.66E+03  2.19E+03  5.41E+06 -7.21E+02  2.05E+02  0.00E+00  2.27E+02  1.91E+02  2.58E+06
 
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
 #CPUT: Total CPU Time in Seconds,       31.109
Stop Time:
Sat Sep 25 08:50:05 CDT 2021
