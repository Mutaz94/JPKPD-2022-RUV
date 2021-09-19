Sat Sep 18 07:06:15 CDT 2021
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
$DATA ../../../../data/int/D/dat51.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25271.1805406987        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3870E+02  4.4534E+02 -8.4362E+01  2.5246E+02  2.0549E+02 -2.0745E+03 -9.0704E+02 -8.6030E+01 -1.5294E+03 -7.3741E+02
            -5.2664E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -860.044346885769        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1619E+00  1.5866E+00  9.7467E-01  2.3637E+00  9.3078E-01  5.3931E+00  3.9211E+00  1.0008E+00  3.2238E+00  1.6474E+00
             1.2243E+01
 PARAMETER:  2.5004E-01  5.6157E-01  7.4348E-02  9.6023E-01  2.8273E-02  1.7851E+00  1.4664E+00  1.0080E-01  1.2706E+00  5.9922E-01
             2.6049E+00
 GRADIENT:  -1.0920E+01  1.3596E+01 -4.2917E+01  1.0370E+02 -9.4576E+00  1.5533E+02  1.7242E+01  4.3311E+00  4.2461E+01  4.0843E+01
             4.9898E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -927.487394365142        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.1422E+00  8.2425E-01  3.3462E+01  5.8632E+00  2.0350E+00  3.6043E+00  7.5665E+00  1.1562E+00  4.4792E+00  1.8642E+00
             1.1701E+01
 PARAMETER:  2.3298E-01 -9.3285E-02  3.6104E+00  1.8687E+00  8.1049E-01  1.3821E+00  2.1237E+00  2.4516E-01  1.5994E+00  7.2284E-01
             2.5596E+00
 GRADIENT:  -1.7568E+01  1.2427E+01 -6.5806E+00  1.1254E+02 -2.5488E+01  1.1602E+02  6.9792E+00 -7.2483E-01  1.1025E+01  4.4904E+01
             4.6956E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1187.53870478256        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0229E+00  7.1924E-01  4.2454E+01  1.5135E+00  2.2878E+00  1.9541E+00  7.0808E+00  2.3447E+00  1.6601E+00  8.8910E-01
             9.4387E+00
 PARAMETER:  1.2262E-01 -2.2956E-01  3.8484E+00  5.1443E-01  9.2760E-01  7.6995E-01  2.0574E+00  9.5217E-01  6.0690E-01 -1.7547E-02
             2.3448E+00
 GRADIENT:  -5.6840E+01  6.9668E+00  2.7286E-01 -4.3813E+01  1.1666E+01 -2.7717E+01  4.8111E+01 -7.2774E-02  3.0547E+01  8.5412E+00
             3.2704E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1225.79424768398        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0825E+00  3.3657E-01  7.5565E+01  1.5773E+00  2.2953E+00  2.2091E+00  7.0096E+00  1.5412E+00  1.3616E+00  6.7396E-01
             8.0261E+00
 PARAMETER:  1.7925E-01 -9.8896E-01  4.4250E+00  5.5572E-01  9.3089E-01  8.9258E-01  2.0473E+00  5.3259E-01  4.0865E-01 -2.9459E-01
             2.1827E+00
 GRADIENT:  -1.0792E+01 -6.8776E+00  7.6865E-01 -1.2089E+01  9.6127E+00 -4.7148E+00  6.7115E+00 -2.6279E-02  1.1123E+01  6.5739E-01
             4.8661E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1228.32990789815        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.1104E+00  9.0749E-01  1.6249E+01  1.2601E+00  2.1478E+00  2.2529E+00  5.0878E+00  9.0037E-01  1.0686E+00  6.8591E-01
             7.8255E+00
 PARAMETER:  2.0472E-01  2.9307E-03  2.8880E+00  3.3116E-01  8.6445E-01  9.1220E-01  1.7268E+00 -4.9458E-03  1.6633E-01 -2.7701E-01
             2.1574E+00
 GRADIENT:   2.6987E+00  1.9744E+00  1.3933E+00  3.1808E+00  4.8654E+00  1.3906E+00 -2.2782E+00 -9.8430E-02  2.2904E+00 -5.2653E-02
            -8.3319E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1228.58171973314        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  1.1175E+00  1.0697E+00  7.7053E+00  1.1345E+00  1.9535E+00  2.2387E+00  4.7525E+00  6.6207E-01  9.0883E-01  7.0134E-01
             7.8399E+00
 PARAMETER:  2.1112E-01  1.6734E-01  2.1419E+00  2.2622E-01  7.6960E-01  9.0591E-01  1.6587E+00 -3.1238E-01  4.4056E-03 -2.5476E-01
             2.1592E+00
 GRADIENT:   5.3187E+00  2.3830E+00  2.3751E+00 -4.1275E+00 -4.4275E+00 -5.4448E-01 -2.0089E+00 -1.4766E-01  1.1117E+00  6.4039E-01
            -6.4402E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1231.31162443050        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.0862E+00  8.2465E-01  6.2507E+00  1.2892E+00  1.8371E+00  2.2323E+00  5.2500E+00  2.1867E+00  9.6297E-01  5.9667E-01
             7.8725E+00
 PARAMETER:  1.8272E-01 -9.2796E-02  1.9327E+00  3.5402E-01  7.0816E-01  9.0302E-01  1.7582E+00  8.8237E-01  6.2264E-02 -4.1639E-01
             2.1634E+00
 GRADIENT:  -8.1272E+00  8.0010E-01 -4.4316E+00  2.0463E+01  1.3516E+01 -1.0938E+00 -2.0389E+00 -1.3416E+00 -6.7386E+00  2.8435E+00
             1.3247E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1234.04958442433        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      593
 NPARAMETR:  1.0946E+00  7.5622E-01  9.3658E+00  1.3237E+00  1.8431E+00  2.2286E+00  5.3516E+00  3.2370E+00  1.2128E+00  3.9425E-01
             7.7646E+00
 PARAMETER:  1.9039E-01 -1.7942E-01  2.3371E+00  3.8042E-01  7.1144E-01  9.0138E-01  1.7774E+00  1.2747E+00  2.9295E-01 -8.3077E-01
             2.1496E+00
 GRADIENT:  -2.9748E+00 -1.8786E-01  3.4032E+00 -3.9400E+00 -8.2695E+00 -2.3168E+00 -2.0931E-01 -2.4884E+00  3.9923E-01  1.6684E+00
            -7.9651E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1234.52954061037        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  1.1055E+00  7.6499E-01  8.8862E+00  1.3391E+00  1.8309E+00  2.2306E+00  5.3484E+00  3.3419E+00  1.2313E+00  3.0462E-01
             7.7724E+00
 PARAMETER:  2.0027E-01 -1.6789E-01  2.2845E+00  3.9202E-01  7.0479E-01  9.0229E-01  1.7768E+00  1.3065E+00  3.0803E-01 -1.0887E+00
             2.1506E+00
 GRADIENT:   8.0942E-01  8.9761E-01  5.5286E-01  3.2893E+00 -5.5034E+00 -2.0043E+00 -1.9369E-01  5.9847E-01  1.9563E-01  1.0447E+00
            -1.3686E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1234.65527073354        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  1.1041E+00  8.0787E-01  8.1425E+00  1.3009E+00  1.8365E+00  2.2468E+00  5.2204E+00  3.1646E+00  1.2042E+00  1.5309E-01
             7.7737E+00
 PARAMETER:  1.9902E-01 -1.1335E-01  2.1971E+00  3.6302E-01  7.0786E-01  9.0952E-01  1.7526E+00  1.2520E+00  2.8578E-01 -1.7767E+00
             2.1507E+00
 GRADIENT:   4.0114E-01  5.8904E-02 -6.4653E-03 -8.9921E-01  1.0693E-01  6.3595E-01  6.9186E-02 -6.5008E-01  4.4606E-01  2.4134E-01
            -1.1966E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1236.04359277845        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      876
 NPARAMETR:  1.1172E+00  6.3095E-01  1.1494E+01  1.4323E+00  1.8773E+00  2.2765E+00  6.0736E+00  3.6878E+00  1.2284E+00  1.0000E-02
             7.8062E+00
 PARAMETER:  2.1086E-01 -3.6053E-01  2.5418E+00  4.5930E-01  7.2985E-01  9.2266E-01  1.9039E+00  1.4050E+00  3.0569E-01 -5.4043E+00
             2.1549E+00
 GRADIENT:   3.0579E+00  2.6194E+00  1.5661E+00  3.5945E+00 -6.5065E+00 -2.3044E-01 -2.1522E+00 -2.2580E-01 -2.3963E+00  0.0000E+00
             1.6182E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1236.47094060230        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  1.1088E+00  4.6566E-01  1.1632E+01  1.5165E+00  1.8630E+00  2.2822E+00  6.7549E+00  3.5450E+00  1.3084E+00  1.0000E-02
             7.7903E+00
 PARAMETER:  2.0326E-01 -6.6430E-01  2.5538E+00  5.1639E-01  7.2219E-01  9.2515E-01  2.0103E+00  1.3655E+00  3.6883E-01 -5.0570E+00
             2.1529E+00
 GRADIENT:   1.9794E-01 -3.0377E-02  1.1943E-01 -2.4217E-01 -6.8132E-02  3.0042E-01  2.3428E-01 -6.0137E-02 -7.9093E-02  0.0000E+00
             3.3735E-01

0ITERATION NO.:   64    OBJECTIVE VALUE:  -1236.47174650011        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1181
 NPARAMETR:  1.1082E+00  4.7146E-01  1.1478E+01  1.5134E+00  1.8603E+00  2.2802E+00  6.7170E+00  3.5264E+00  1.3088E+00  1.0000E-02
             7.7889E+00
 PARAMETER:  2.0270E-01 -6.5192E-01  2.5404E+00  5.1438E-01  7.2076E-01  9.2428E-01  2.0046E+00  1.3603E+00  3.6910E-01 -4.9855E+00
             2.1527E+00
 GRADIENT:  -2.3420E-02 -1.3368E-02  4.8140E-03 -9.1273E-03  8.3183E-03 -1.5954E-02 -3.7555E-02 -2.8122E-03 -9.9574E-03  0.0000E+00
            -8.7857E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1181
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.3419E-03  5.3134E-02 -4.0620E-02 -6.8684E-02 -8.1772E-05
 SE:             2.9068E-02  2.1928E-02  1.3600E-02  1.7653E-02  1.4183E-04
 N:                     100         100         100         100         100

 P VAL.:         7.4792E-01  1.5387E-02  2.8189E-03  9.9946E-05  5.6424E-01

 ETASHRINKSD(%)  2.6181E+00  2.6540E+01  5.4439E+01  4.0860E+01  9.9525E+01
 ETASHRINKVR(%)  5.1676E+00  4.6036E+01  7.9242E+01  6.5025E+01  9.9998E+01
 EBVSHRINKSD(%)  2.8852E+00  2.4188E+01  6.1709E+01  3.8053E+01  9.9451E+01
 EBVSHRINKVR(%)  5.6872E+00  4.2525E+01  8.5338E+01  6.1625E+01  9.9997E+01
 RELATIVEINF(%)  9.4188E+01  2.8349E+01  5.6848E+00  1.9637E+01  1.2044E-03
 EPSSHRINKSD(%)  8.4604E+00
 EPSSHRINKVR(%)  1.6205E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1236.4717465001081     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       417.61761326830265     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.49
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1236.472       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.11E+00  4.71E-01  1.15E+01  1.51E+00  1.86E+00  2.28E+00  6.72E+00  3.53E+00  1.31E+00  1.00E-02  7.79E+00
 


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
+        1.65E+02
 
 TH 2
+       -7.29E-01  7.15E+01
 
 TH 3
+        7.36E-02  4.90E-01  1.72E-01
 
 TH 4
+       -3.36E+00  3.84E+01 -4.91E-01  1.39E+02
 
 TH 5
+       -1.60E+00 -1.42E+01 -2.29E+00 -4.45E+00  8.04E+01
 
 TH 6
+        4.18E+02 -3.60E+01  1.30E+00 -1.75E+01 -7.93E+04  5.07E+04
 
 TH 7
+        2.01E-01  7.14E+00 -5.84E-02 -5.83E+00  8.92E-01 -7.88E+03  2.24E+00
 
 TH 8
+       -1.79E-02 -1.21E+00 -4.75E-01  1.84E+00  1.21E+00 -4.49E-01  6.57E-02  2.46E+00
 
 TH 9
+       -4.84E-01 -3.91E+00 -3.48E-01 -2.74E+01  1.01E+01 -1.64E+00  1.20E+00  1.12E+00  2.53E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.33E+01  2.75E-01 -1.34E-01 -8.15E+00  5.24E+01 -6.39E+03  2.30E+00  2.76E-01  3.40E+00  0.00E+00  8.24E+02
 
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
 #CPUT: Total CPU Time in Seconds,       45.513
Stop Time:
Sat Sep 18 07:07:02 CDT 2021
