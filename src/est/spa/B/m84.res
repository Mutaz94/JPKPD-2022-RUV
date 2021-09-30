Wed Sep 29 11:34:26 CDT 2021
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
$DATA ../../../../data/spa/B/dat84.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1701.21070943331        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5852E+02 -1.1183E+01 -4.0633E+01  3.6134E+01  2.4713E+01  8.4334E+01 -8.0720E+00  1.4884E+01 -4.7136E+00  1.5350E+01
            -2.2969E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1714.43503690188        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0355E+00  1.0812E+00  1.2249E+00  9.9635E-01  1.1096E+00  8.0829E-01  1.0729E+00  8.9341E-01  1.0741E+00  9.2507E-01
             1.0924E+00
 PARAMETER:  1.3491E-01  1.7808E-01  3.0284E-01  9.6342E-02  2.0400E-01 -1.1283E-01  1.7034E-01 -1.2707E-02  1.7144E-01  2.2109E-02
             1.8839E-01
 GRADIENT:  -5.2478E+00  7.1111E+00 -2.5840E+00  1.1973E+01  1.1285E+01 -1.7490E+01  1.1167E+00  1.8381E+00  1.1439E+00 -1.3507E+01
             4.6107E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1715.07470917386        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0353E+00  1.0389E+00  1.3507E+00  1.0143E+00  1.1356E+00  8.1496E-01  1.0277E+00  7.7613E-01  1.1038E+00  1.0628E+00
             1.0927E+00
 PARAMETER:  1.3465E-01  1.3814E-01  4.0064E-01  1.1423E-01  2.2713E-01 -1.0461E-01  1.2731E-01 -1.5343E-01  1.9879E-01  1.6093E-01
             1.8867E-01
 GRADIENT:  -3.0657E+00  7.3146E-01  1.0472E+01 -5.6347E+00 -7.9291E+00 -1.3472E+01  1.5293E+00 -3.3398E+00  5.9183E+00 -2.0786E-01
             7.0281E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1715.93937381009        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0358E+00  9.0154E-01  1.3811E+00  1.1071E+00  1.0970E+00  8.4186E-01  1.1298E+00  8.8850E-01  1.0039E+00  1.0473E+00
             1.0698E+00
 PARAMETER:  1.3521E-01 -3.6517E-03  4.2290E-01  2.0179E-01  1.9257E-01 -7.2146E-02  2.2203E-01 -1.8221E-02  1.0392E-01  1.4624E-01
             1.6747E-01
 GRADIENT:   2.1948E+00  6.8185E-01 -1.0307E+00  3.0420E-01  3.1810E-01 -1.9573E-01 -2.6237E-01  1.1818E-01 -2.0088E-02  5.0502E-01
             6.8617E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1716.12779600476        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  1.0307E+00  6.4363E-01  1.6910E+00  1.2857E+00  1.0964E+00  8.3954E-01  1.3231E+00  1.0981E+00  9.0844E-01  1.0745E+00
             1.0648E+00
 PARAMETER:  1.3028E-01 -3.4063E-01  6.2532E-01  3.5127E-01  1.9206E-01 -7.4904E-02  3.7999E-01  1.9360E-01  3.9732E-03  1.7189E-01
             1.6283E-01
 GRADIENT:  -4.2471E+00  5.6223E+00  3.1715E+00  8.6208E+00 -4.0733E+00  1.2024E-01 -2.2392E-01 -5.4703E-01 -1.1156E+00 -8.9567E-01
            -1.5833E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1716.22335409246        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0300E+00  5.0869E-01  1.8242E+00  1.3782E+00  1.0899E+00  8.3940E-01  1.4986E+00  1.2070E+00  8.6540E-01  1.0853E+00
             1.0646E+00
 PARAMETER:  1.2953E-01 -5.7591E-01  7.0115E-01  4.2080E-01  1.8605E-01 -7.5070E-02  5.0450E-01  2.8816E-01 -4.4565E-02  1.8187E-01
             1.6261E-01
 GRADIENT:  -1.9464E+00  6.0115E+00  1.5864E+00  1.4390E+01 -4.1141E+00  7.7946E-01 -1.9600E-01  3.0550E-01 -1.4493E+00 -2.1367E-01
            -1.2026E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1716.38586741097        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0293E+00  3.5687E-01  1.9196E+00  1.4777E+00  1.0706E+00  8.3731E-01  1.8268E+00  1.2679E+00  8.2409E-01  1.0877E+00
             1.0659E+00
 PARAMETER:  1.2890E-01 -9.3039E-01  7.5214E-01  4.9047E-01  1.6820E-01 -7.7564E-02  7.0259E-01  3.3735E-01 -9.3479E-02  1.8404E-01
             1.6382E-01
 GRADIENT:   1.6162E+00  4.6638E+00  3.1272E-01  1.6081E+01 -3.2102E+00  6.5744E-01 -1.0892E-01  1.2235E-01 -1.3065E+00  4.5853E-01
            -1.1649E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1716.53706590179        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1257
 NPARAMETR:  1.0279E+00  2.1546E-01  2.0960E+00  1.5739E+00  1.0724E+00  8.3470E-01  2.3654E+00  1.4061E+00  7.9121E-01  1.1006E+00
             1.0668E+00
 PARAMETER:  1.2748E-01 -1.4350E+00  8.4005E-01  5.5358E-01  1.6989E-01 -8.0682E-02  9.6096E-01  4.4080E-01 -1.3419E-01  1.9583E-01
             1.6465E-01
 GRADIENT:   2.2887E+00  3.5868E+00 -4.7843E-01  1.8471E+01 -3.7068E+00  3.2771E-01  4.4332E-02  4.2166E-01 -8.4519E-01  1.0732E+00
             4.9071E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1716.78314632970        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1432
 NPARAMETR:  1.0263E+00  1.0595E-01  2.3210E+00  1.6480E+00  1.0910E+00  8.3224E-01  3.3634E+00  1.5674E+00  7.6877E-01  1.1175E+00
             1.0672E+00
 PARAMETER:  1.2600E-01 -2.1448E+00  9.4201E-01  5.9957E-01  1.8710E-01 -8.3633E-02  1.3129E+00  5.4942E-01 -1.6297E-01  2.1110E-01
             1.6501E-01
 GRADIENT:   1.0673E+00  1.8288E+00 -6.5892E-01  1.2843E+01 -2.8054E+00 -1.4312E-01  9.1830E-02  5.6992E-01  3.4231E-01  1.0904E+00
             6.3318E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1716.91318954416        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1608
 NPARAMETR:  1.0257E+00  4.9761E-02  2.5108E+00  1.6879E+00  1.1117E+00  8.3157E-01  4.9150E+00  1.6807E+00  7.5223E-01  1.1267E+00
             1.0670E+00
 PARAMETER:  1.2534E-01 -2.9005E+00  1.0206E+00  6.2346E-01  2.0585E-01 -8.4438E-02  1.6923E+00  6.1920E-01 -1.8472E-01  2.1929E-01
             1.6489E-01
 GRADIENT:   9.2205E-02  9.2959E-01  1.3329E-01  9.4175E+00 -1.7044E+00 -1.3602E-01  4.4911E-02 -6.9768E-02 -1.2325E+00  1.9466E-01
            -2.3790E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1717.10781014569        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1793             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0260E+00  1.3607E-02  2.5720E+00  1.6921E+00  1.1158E+00  8.3156E-01  9.0347E+00  1.7199E+00  7.4990E-01  1.1312E+00
             1.0663E+00
 PARAMETER:  1.2570E-01 -4.1972E+00  1.0447E+00  6.2597E-01  2.0961E-01 -8.4453E-02  2.3011E+00  6.4229E-01 -1.8782E-01  2.2324E-01
             1.6416E-01
 GRADIENT:   4.9291E+02  1.9764E+00  6.4560E+00  1.0923E+03  1.1296E+01  2.4023E+01  8.1022E-01  1.4172E+00  2.2419E+01  1.7227E+00
             1.6935E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1717.14456494557        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1977             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0256E+00  1.0000E-02  2.5672E+00  1.7016E+00  1.1132E+00  8.3137E-01  1.2625E+01  1.7187E+00  7.4745E-01  1.1287E+00
             1.0659E+00
 PARAMETER:  1.2530E-01 -4.5210E+00  1.0428E+00  6.3159E-01  2.0724E-01 -8.4675E-02  2.6357E+00  6.4160E-01 -1.9108E-01  2.2110E-01
             1.6384E-01
 GRADIENT:   4.9005E+02  9.0935E-01  6.2106E+00  1.1224E+03  9.8358E+00  2.3967E+01  1.0560E+00  1.4654E+00  2.2427E+01  1.6424E+00
             1.3360E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1717.15089111359        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     2138
 NPARAMETR:  1.0255E+00  1.0000E-02  2.5644E+00  1.7003E+00  1.1126E+00  8.3123E-01  1.1031E+01  1.7180E+00  7.4698E-01  1.1280E+00
             1.0658E+00
 PARAMETER:  1.2514E-01 -4.5566E+00  1.0417E+00  6.3083E-01  2.0666E-01 -8.4844E-02  2.5007E+00  6.4114E-01 -1.9172E-01  2.2048E-01
             1.6377E-01
 GRADIENT:   1.3033E+00  0.0000E+00  2.6357E-01 -2.2613E+01  6.9562E-01  1.9130E-03  2.3476E-02  1.0029E-01  2.5191E-01  1.3423E-01
             8.1965E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1717.15397008159        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2324             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0257E+00  1.0000E-02  2.5545E+00  1.7004E+00  1.1105E+00  8.3133E-01  1.0751E+01  1.7127E+00  7.4673E-01  1.1259E+00
             1.0659E+00
 PARAMETER:  1.2540E-01 -4.5566E+00  1.0378E+00  6.3085E-01  2.0485E-01 -8.4734E-02  2.4750E+00  6.3807E-01 -1.9205E-01  2.1854E-01
             1.6381E-01
 GRADIENT:   4.9094E+02  0.0000E+00  6.2576E+00  1.1186E+03  9.5375E+00  2.3948E+01  6.5388E-01  1.5325E+00  2.1959E+01  1.5115E+00
             1.3071E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1717.15814839461        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2509
 NPARAMETR:  1.0257E+00  1.0000E-02  2.5370E+00  1.6996E+00  1.1074E+00  8.3131E-01  1.1001E+01  1.7013E+00  7.4723E-01  1.1242E+00
             1.0657E+00
 PARAMETER:  1.2538E-01 -4.5566E+00  1.0310E+00  6.3042E-01  2.0205E-01 -8.4752E-02  2.4980E+00  6.3138E-01 -1.9139E-01  2.1704E-01
             1.6367E-01
 GRADIENT:   2.2616E+00  0.0000E+00  2.4557E-01 -2.2306E+01  3.7093E-01  3.9071E-02  2.2070E-02  1.1026E-01  2.1607E-01  8.2741E-02
             4.7624E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1717.16038392303        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2694
 NPARAMETR:  1.0252E+00  1.0000E-02  2.5247E+00  1.6997E+00  1.1045E+00  8.3116E-01  1.0425E+01  1.6921E+00  7.4720E-01  1.1230E+00
             1.0657E+00
 PARAMETER:  1.2489E-01 -4.5566E+00  1.0261E+00  6.3048E-01  1.9944E-01 -8.4930E-02  2.4442E+00  6.2597E-01 -1.9142E-01  2.1602E-01
             1.6362E-01
 GRADIENT:   7.9654E-01  0.0000E+00  4.5617E-01 -2.1145E+01 -3.7147E-01 -2.1032E-02  1.1258E-02  5.1515E-02  7.4901E-02  1.8569E-01
             8.9959E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1717.16172173688        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2885
 NPARAMETR:  1.0257E+00  1.0000E-02  2.5194E+00  1.6990E+00  1.1030E+00  8.3131E-01  1.0705E+01  1.6874E+00  7.4737E-01  1.1223E+00
             1.0657E+00
 PARAMETER:  1.2534E-01 -4.5566E+00  1.0240E+00  6.3005E-01  1.9804E-01 -8.4749E-02  2.4707E+00  6.2317E-01 -1.9119E-01  2.1542E-01
             1.6365E-01
 GRADIENT:   2.3168E+00  0.0000E+00  7.0154E-01 -2.2348E+01 -8.4754E-01  5.2755E-02  1.5875E-02 -1.3193E-03  1.5555E-01  2.4441E-01
             3.3370E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1717.16252817162        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3079             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0256E+00  1.0000E-02  2.5039E+00  1.6988E+00  1.1032E+00  8.3130E-01  1.0558E+01  1.6839E+00  7.4737E-01  1.1182E+00
             1.0656E+00
 PARAMETER:  1.2532E-01 -4.5566E+00  1.0179E+00  6.2989E-01  1.9823E-01 -8.4766E-02  2.4569E+00  6.2111E-01 -1.9120E-01  2.1174E-01
             1.6351E-01
 GRADIENT:   4.9063E+02  0.0000E+00  5.5937E+00  1.1154E+03  1.0744E+01  2.3922E+01  6.2590E-01  1.5315E+00  2.1884E+01  1.0374E+00
             1.2604E+00

0ITERATION NO.:   89    OBJECTIVE VALUE:  -1717.16434001740        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     3219
 NPARAMETR:  1.0256E+00  1.0000E-02  2.4968E+00  1.6985E+00  1.1019E+00  8.3129E-01  1.0572E+01  1.6774E+00  7.4749E-01  1.1182E+00
             1.0654E+00
 PARAMETER:  1.2531E-01 -4.5566E+00  1.0189E+00  6.2987E-01  1.9676E-01 -8.4764E-02  2.4594E+00  6.2006E-01 -1.9112E-01  2.1389E-01
             1.6356E-01
 GRADIENT:   2.0566E-02  0.0000E+00  3.3526E-01  2.1560E-01 -1.2822E-01  1.8266E-03  7.0467E-05  9.6683E-02 -1.2465E-02  1.2400E-01
             3.2604E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3219
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.1172E-04 -1.1171E-03 -2.9527E-02 -7.8101E-03 -4.5037E-02
 SE:             2.9761E-02  1.6232E-03  1.5888E-02  2.9166E-02  2.0810E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8628E-01  4.9132E-01  6.3112E-02  7.8887E-01  3.0452E-02

 ETASHRINKSD(%)  2.9754E-01  9.4562E+01  4.6772E+01  2.2889E+00  3.0283E+01
 ETASHRINKVR(%)  5.9419E-01  9.9704E+01  7.1668E+01  4.5255E+00  5.1396E+01
 EBVSHRINKSD(%)  6.5845E-01  9.4625E+01  5.0930E+01  2.5923E+00  2.6575E+01
 EBVSHRINKVR(%)  1.3126E+00  9.9711E+01  7.5921E+01  5.1174E+00  4.6088E+01
 RELATIVEINF(%)  9.3118E+01  5.4858E-03  6.7745E+00  1.9926E+00  8.9744E+00
 EPSSHRINKSD(%)  4.3230E+01
 EPSSHRINKVR(%)  6.7771E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1717.1643400173978     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -982.01351345365958     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.82
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     6.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1717.164       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.00E-02  2.51E+00  1.70E+00  1.10E+00  8.31E-01  1.06E+01  1.68E+00  7.47E-01  1.12E+00  1.07E+00
 


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
 
         2.59E-02  0.00E+00  1.11E+00  5.13E-02  2.12E-01  4.92E-02  2.44E+00  7.29E-01  5.84E-02  2.12E-01  6.45E-02
 


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
+        6.68E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.38E-03  0.00E+00  1.23E+00
 
 TH 4
+        1.15E-04  0.00E+00  3.65E-02  2.63E-03
 
 TH 5
+        5.52E-04  0.00E+00  2.27E-01  7.07E-03  4.49E-02
 
 TH 6
+       -2.34E-04  0.00E+00 -1.55E-03  7.82E-05 -6.00E-04  2.42E-03
 
 TH 7
+       -2.11E-02  0.00E+00 -4.93E-01 -2.84E-02 -1.01E-01  1.99E-02  5.97E+00
 
 TH 8
+        8.87E-04  0.00E+00  7.58E-01  2.34E-02  1.40E-01 -9.07E-04 -4.02E-01  5.31E-01
 
 TH 9
+       -2.76E-04  0.00E+00 -5.24E-03 -6.58E-04 -5.26E-04  1.06E-04  5.96E-02 -2.81E-03  3.41E-03
 
 TH10
+        6.43E-04  0.00E+00  1.76E-01  5.40E-03  3.54E-02 -3.25E-04 -5.68E-02  1.03E-01  4.00E-04  4.49E-02
 
 TH11
+        1.93E-04  0.00E+00  1.99E-02  3.02E-04  3.89E-03 -4.76E-04 -2.65E-02  1.26E-02 -6.36E-04  1.58E-03  4.16E-03
 
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
+        2.59E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        4.81E-02  0.00E+00  1.11E+00
 
 TH 4
+        8.66E-02  0.00E+00  6.41E-01  5.13E-02
 
 TH 5
+        1.01E-01  0.00E+00  9.67E-01  6.50E-01  2.12E-01
 
 TH 6
+       -1.84E-01  0.00E+00 -2.84E-02  3.10E-02 -5.75E-02  4.92E-02
 
 TH 7
+       -3.34E-01  0.00E+00 -1.82E-01 -2.27E-01 -1.95E-01  1.65E-01  2.44E+00
 
 TH 8
+        4.70E-02  0.00E+00  9.38E-01  6.26E-01  9.08E-01 -2.53E-02 -2.26E-01  7.29E-01
 
 TH 9
+       -1.83E-01  0.00E+00 -8.09E-02 -2.20E-01 -4.25E-02  3.69E-02  4.18E-01 -6.59E-02  5.84E-02
 
 TH10
+        1.17E-01  0.00E+00  7.48E-01  4.97E-01  7.89E-01 -3.11E-02 -1.10E-01  6.69E-01  3.23E-02  2.12E-01
 
 TH11
+        1.16E-01  0.00E+00  2.78E-01  9.13E-02  2.84E-01 -1.50E-01 -1.68E-01  2.67E-01 -1.69E-01  1.16E-01  6.45E-02
 
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
+        1.65E+03
 
 TH 2
+       -1.07E-13  6.40E-29
 
 TH 3
+        1.76E+01 -8.73E-15  9.01E+00
 
 TH 4
+       -2.62E+01  4.69E-14 -4.42E+00  6.91E+02
 
 TH 5
+       -1.32E+02  3.78E-14 -6.08E+01 -8.72E+01  4.32E+02
 
 TH 6
+        1.28E+02 -1.82E-13 -6.31E+00 -8.16E+01  6.82E+01  3.99E+02
 
 TH 7
+        2.44E-02 -2.41E-17  1.85E-03  1.25E-02 -1.16E-02 -4.33E-03  9.81E-06
 
 TH 8
+        2.71E+00 -4.64E-16  1.09E+00 -1.75E-01 -7.44E+00 -2.23E+00  2.79E-04  1.38E-01
 
 TH 9
+        1.33E+02 -1.34E-13  9.12E+00  1.07E+02 -6.23E+01 -2.25E+01  5.55E-02  1.41E+00  3.16E+02
 
 TH10
+        1.96E+01 -3.64E-15  9.67E+00  6.01E+00 -6.74E+01 -1.85E+01  2.21E-03  1.22E+00  1.14E+01  1.08E+01
 
 TH11
+       -4.39E+01  3.27E-14 -1.10E+00 -2.64E+01  1.01E+01 -6.35E+01 -9.65E-04  9.01E-02 -7.56E+00  1.85E-01  1.33E+01
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.51E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -3.23E+01  0.00E+00  1.58E+01
 
 TH 4
+        2.36E+01  0.00E+00 -2.07E+01  6.55E+02
 
 TH 5
+        1.21E+02  0.00E+00 -5.69E+01 -1.60E+01  4.13E+02
 
 TH 6
+       -9.47E+01  0.00E+00  5.11E+00  2.46E+01 -3.84E+01  1.91E+02
 
 TH 7
+       -4.08E-02  0.00E+00  9.35E-04 -3.28E-02  1.20E-02  2.28E-03  2.13E-05
 
 TH 8
+        2.02E+00  0.00E+00 -7.84E+00  4.98E+00 -8.67E+00  7.00E-01  2.83E-04  1.38E+01
 
 TH 9
+       -1.24E+02  0.00E+00 -1.85E+00 -1.16E+02  5.59E+01  2.38E+00  8.10E-02  3.07E+00  3.62E+02
 
 TH10
+        1.35E+01  0.00E+00  9.95E-01 -8.84E+00 -4.67E+01  1.01E+00  7.63E-03  5.28E+00  1.57E+01  5.31E+01
 
 TH11
+        2.52E+01  0.00E+00 -7.62E+00 -6.39E+01  3.93E+01 -2.06E+01  4.35E-03  4.48E+00  2.65E+00  1.50E+01  1.48E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,       50.052
Stop Time:
Wed Sep 29 11:35:18 CDT 2021
