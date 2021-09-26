Sat Sep 25 01:56:26 CDT 2021
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
$DATA ../../../../data/int/SL3/dat14.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      981
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      881
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   234.508038607567        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5845E+01 -1.2967E+02  8.8169E+01  1.2016E+02  2.1359E+02 -1.3183E+01 -1.6196E+02 -2.4395E+02 -1.1161E+02 -3.5303E+01
            -7.6117E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2377.98449079392        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1161E+00  1.5123E+00  8.9037E-01  7.9044E-01  1.0909E+00  9.9817E-01  1.1307E+00  1.0560E+00  7.7011E-01  9.6515E-01
             5.1693E+00
 PARAMETER:  2.0986E-01  5.1361E-01 -1.6118E-02 -1.3517E-01  1.8705E-01  9.8173E-02  2.2284E-01  1.5447E-01 -1.6122E-01  6.4531E-02
             1.7427E+00
 GRADIENT:   7.6853E+01 -1.6309E+01  2.8601E+00 -6.3963E+01 -3.8596E+01  2.0264E+00  3.9469E+01  3.6839E+00  9.4681E+00  2.7799E+00
             7.8361E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2452.86138623149        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0902E+00  1.8833E+00  1.2511E+00  6.4406E-01  1.2681E+00  1.1238E+00  7.1320E-01  4.1816E+00  4.1060E-01  1.7965E+00
             4.5623E+00
 PARAMETER:  1.8634E-01  7.3301E-01  3.2401E-01 -3.3996E-01  3.3755E-01  2.1671E-01 -2.3799E-01  1.5307E+00 -7.9013E-01  6.8582E-01
             1.6178E+00
 GRADIENT:   3.4674E+01  1.4683E+02 -7.5788E+00  1.0080E+02 -7.5588E+01  3.3852E+01 -1.7090E+01  2.3628E+01 -3.6880E+00  4.3628E+01
             6.5498E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2656.43928354738        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0366E+00  1.6449E+00  5.1845E+00  7.1378E-01  1.8987E+00  1.0720E+00  8.6689E-01  2.9511E+00  4.1372E-01  1.4977E+00
             2.9553E+00
 PARAMETER:  1.3596E-01  5.9767E-01  1.7457E+00 -2.3719E-01  7.4119E-01  1.6951E-01 -4.2840E-02  1.1822E+00 -7.8256E-01  5.0391E-01
             1.1836E+00
 GRADIENT:   7.3310E+00  1.6695E+00 -6.0329E+00  1.9459E+01  8.8028E+01  1.6321E+01 -1.4201E+01 -6.7288E+00 -3.9723E+00 -3.9677E+01
             4.6853E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2666.36048958048        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0299E+00  1.6592E+00  5.1884E+00  6.9314E-01  1.6630E+00  1.0203E+00  8.9282E-01  4.1851E+00  5.7456E-01  1.6125E+00
             2.8574E+00
 PARAMETER:  1.2942E-01  6.0636E-01  1.7464E+00 -2.6652E-01  6.0860E-01  1.2011E-01 -1.3365E-02  1.5315E+00 -4.5415E-01  5.7780E-01
             1.1499E+00
 GRADIENT:  -1.4753E+00  4.4781E+00 -7.9638E-01 -5.8854E-01 -5.7861E+00 -1.3037E+00  4.7937E+00  1.3914E+00 -1.1579E+00 -2.4771E-01
            -8.2469E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2666.79629580565        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      416
 NPARAMETR:  1.0357E+00  1.6675E+00  5.5394E+00  6.9799E-01  1.6958E+00  1.0304E+00  8.4723E-01  4.3036E+00  7.0051E-01  1.6228E+00
             2.8673E+00
 PARAMETER:  1.3507E-01  6.1132E-01  1.8119E+00 -2.5955E-01  6.2814E-01  1.2996E-01 -6.5784E-02  1.5595E+00 -2.5594E-01  5.8414E-01
             1.1534E+00
 GRADIENT:   1.8689E+00  3.6477E+00 -1.7340E-01  5.1606E+00  1.6826E+00  1.4384E+00  1.4299E+00  4.7090E-01 -2.9747E-01  1.0867E-01
            -1.6680E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2667.73693588003        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      591
 NPARAMETR:  1.0325E+00  1.9875E+00  3.9001E+00  4.7787E-01  1.7575E+00  1.0290E+00  7.4166E-01  4.6925E+00  8.3493E-01  1.6923E+00
             2.8656E+00
 PARAMETER:  1.3197E-01  7.8689E-01  1.4610E+00 -6.3841E-01  6.6388E-01  1.2860E-01 -1.9887E-01  1.6460E+00 -8.0412E-02  6.2607E-01
             1.1528E+00
 GRADIENT:  -4.5335E+00  3.3054E+00  6.4820E-01  1.5981E+00 -1.0925E-01  9.3173E-01  1.8393E+00 -6.7129E-01  2.5296E-02 -3.1167E+00
            -3.0911E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2668.75788708954        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      767
 NPARAMETR:  1.0341E+00  2.3252E+00  2.0334E+00  2.5167E-01  1.8235E+00  1.0300E+00  6.6123E-01  5.5533E+00  1.0160E+00  1.8015E+00
             2.8644E+00
 PARAMETER:  1.3352E-01  9.4382E-01  8.0971E-01 -1.2796E+00  7.0074E-01  1.2954E-01 -3.1365E-01  1.8144E+00  1.1583E-01  6.8861E-01
             1.1523E+00
 GRADIENT:  -1.8672E+00  1.9651E+01  2.0702E-01  3.7016E+00 -1.2014E+00  1.1875E+00  7.0977E-01 -8.7996E-01 -9.8590E-02 -4.9053E-01
            -4.4829E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2669.56537667154        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  1.0377E+00  2.5195E+00  8.7511E-01  1.1384E-01  1.8661E+00  1.0204E+00  6.1622E-01  8.1353E+00  1.2698E+00  1.8721E+00
             2.8774E+00
 PARAMETER:  1.3697E-01  1.0241E+00 -3.3401E-02 -2.0730E+00  7.2386E-01  1.2019E-01 -3.8415E-01  2.1962E+00  3.3887E-01  7.2705E-01
             1.1569E+00
 GRADIENT:   4.6777E+00  8.2206E+00 -2.8486E-01  1.8725E+00  3.9968E-01 -2.1465E+00 -1.0319E+00  3.6228E-01  4.3839E-02  3.0916E+00
             7.7256E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2669.91119407257        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1121
 NPARAMETR:  1.0346E+00  2.5916E+00  4.9504E-01  6.5748E-02  1.8766E+00  1.0273E+00  6.0531E-01  1.0868E+01  1.3813E+00  1.8632E+00
             2.8651E+00
 PARAMETER:  1.3397E-01  1.0523E+00 -6.0311E-01 -2.6219E+00  7.2948E-01  1.2694E-01 -4.0202E-01  2.4859E+00  4.2301E-01  7.2232E-01
             1.1526E+00
 GRADIENT:  -1.0821E+00  1.3989E+01 -3.0409E-01  1.6838E+00 -2.8535E-01  2.3174E-01 -1.6965E-01  7.5913E-01  5.3100E-02 -5.3868E-01
            -2.1462E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2670.22383492251        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1299
 NPARAMETR:  1.0343E+00  2.6414E+00  2.0832E-01  2.8323E-02  1.8881E+00  1.0277E+00  5.9420E-01  1.6053E+01  1.4949E+00  1.8755E+00
             2.8637E+00
 PARAMETER:  1.3376E-01  1.0713E+00 -1.4687E+00 -3.4641E+00  7.3558E-01  1.2732E-01 -4.2054E-01  2.8759E+00  5.0205E-01  7.2888E-01
             1.1521E+00
 GRADIENT:  -1.3364E+00  4.9698E+00 -1.5296E-01  6.0464E-01  1.3670E-01  3.7068E-01 -2.6020E-01  2.8999E-01  6.6497E-03 -5.2683E-01
            -2.1309E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2670.22711665484        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1464
 NPARAMETR:  1.0343E+00  2.6403E+00  2.0623E-01  2.8112E-02  1.8889E+00  1.0264E+00  5.9460E-01  1.6144E+01  1.4962E+00  1.8750E+00
             2.8653E+00
 PARAMETER:  1.3371E-01  1.0709E+00 -1.4788E+00 -3.4716E+00  7.3599E-01  1.2603E-01 -4.1987E-01  2.8816E+00  5.0296E-01  7.2859E-01
             1.1527E+00
 GRADIENT:  -1.4731E+00  2.0004E+00 -1.5748E-01  5.6458E-01  4.3867E-01 -8.2869E-02  6.1924E-05  3.0176E-01  7.7707E-03 -5.9209E-01
            -6.3158E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2670.29423153362        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1643
 NPARAMETR:  1.0335E+00  2.6398E+00  2.6969E-01  2.2992E-02  1.8890E+00  1.0255E+00  5.9197E-01  1.4846E+01  1.5431E+00  1.8868E+00
             2.8638E+00
 PARAMETER:  1.3299E-01  1.0707E+00 -1.2105E+00 -3.6726E+00  7.3603E-01  1.2515E-01 -4.2430E-01  2.7977E+00  5.3377E-01  7.3486E-01
             1.1521E+00
 GRADIENT:  -2.6142E+00 -1.3892E+01 -6.2310E-02  2.8578E-02 -5.3941E-01 -3.3434E-01 -4.8085E-01 -4.9977E-03 -1.2114E-02  2.0915E-01
            -7.9075E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2670.43118880721        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1821
 NPARAMETR:  1.0353E+00  2.6698E+00  7.1995E-01  1.0000E-02  1.8963E+00  1.0261E+00  5.8783E-01  1.6731E+01  3.3397E+00  1.8891E+00
             2.8653E+00
 PARAMETER:  1.3465E-01  1.0820E+00 -2.2857E-01 -4.5650E+00  7.3990E-01  1.2574E-01 -4.3132E-01  2.9173E+00  1.3059E+00  7.3609E-01
             1.1527E+00
 GRADIENT:   3.0306E-01  7.7853E+00 -1.2256E-02  0.0000E+00 -8.3789E-02 -2.4235E-01 -4.3407E-01 -7.0942E-04 -8.9297E-03  1.2676E-01
            -3.2655E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2670.45379533024        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2000
 NPARAMETR:  1.0352E+00  2.6640E+00  3.4253E+00  1.0000E-02  1.8945E+00  1.0259E+00  5.8821E-01  1.6200E+01  6.9097E+00  1.8888E+00
             2.8650E+00
 PARAMETER:  1.3458E-01  1.0798E+00  1.3312E+00 -5.1800E+00  7.3897E-01  1.2562E-01 -4.3067E-01  2.8850E+00  2.0329E+00  7.3592E-01
             1.1526E+00
 GRADIENT:   3.9755E-01 -4.2401E+00 -2.8738E-03  0.0000E+00 -2.9832E-01 -2.1429E-01  2.7622E-01 -2.9485E-05 -8.8460E-03  9.4630E-02
             1.2035E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2670.45913445621        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2179
 NPARAMETR:  1.0350E+00  2.6658E+00  2.2676E+01  1.0000E-02  1.8958E+00  1.0265E+00  5.8742E-01  1.6254E+01  7.9661E+00  1.8884E+00
             2.8651E+00
 PARAMETER:  1.3445E-01  1.0805E+00  3.2213E+00 -5.9307E+00  7.3963E-01  1.2613E-01 -4.3202E-01  2.8884E+00  2.1752E+00  7.3573E-01
             1.1526E+00
 GRADIENT:   7.6345E-02 -6.8918E-01 -4.6382E-04  0.0000E+00 -4.9126E-02 -4.5210E-02  5.3653E-02 -7.6304E-07 -2.8004E-04  2.3234E-02
            -8.1360E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2670.45962972541        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2353
 NPARAMETR:  1.0351E+00  2.6661E+00  4.7002E+02  1.0000E-02  1.8960E+00  1.0267E+00  5.8730E-01  1.6283E+01  7.9808E+00  1.8883E+00
             2.8651E+00
 PARAMETER:  1.3450E-01  1.0806E+00  6.2521E+00 -7.0775E+00  7.3978E-01  1.2644E-01 -4.3220E-01  2.8873E+00  2.1779E+00  7.3564E-01
             1.1526E+00
 GRADIENT:   9.6461E-02  8.4584E-02 -5.9205E-04  0.0000E+00  2.3120E-02  1.0092E-01  5.5825E-03 -3.4126E-04  2.5518E-03 -1.7079E-02
            -6.7103E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2353
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.2739E-04 -1.0483E-02  1.7681E-07 -7.9357E-04 -1.5278E-02
 SE:             2.9281E-02  2.7548E-02  1.1018E-07  2.0072E-03  2.7011E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7746E-01  7.0354E-01  1.0858E-01  6.9257E-01  5.7165E-01

 ETASHRINKSD(%)  1.9062E+00  7.7092E+00  1.0000E+02  9.3276E+01  9.5097E+00
 ETASHRINKVR(%)  3.7760E+00  1.4824E+01  1.0000E+02  9.9548E+01  1.8115E+01
 EBVSHRINKSD(%)  1.7375E+00  7.8994E+00  1.0000E+02  9.3290E+01  7.9151E+00
 EBVSHRINKVR(%)  3.4448E+00  1.5175E+01  1.0000E+02  9.9550E+01  1.5204E+01
 RELATIVEINF(%)  9.6423E+01  1.9442E-04  3.4943E-10  1.0321E-06  5.3680E+01
 EPSSHRINKSD(%)  1.5400E+01
 EPSSHRINKVR(%)  2.8429E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          881
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1619.1696955066332     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2670.4596297254143     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1051.2899342187811     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    55.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2670.460       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.67E+00  4.70E+02  1.00E-02  1.90E+00  1.03E+00  5.87E-01  1.62E+01  7.99E+00  1.89E+00  2.87E+00
 


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
+        1.05E+03
 
 TH 2
+       -1.15E+01  3.85E+02
 
 TH 3
+        1.16E-02 -1.61E-04  3.93E-07
 
 TH 4
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 5
+       -2.60E+01 -1.46E+01 -2.52E-04  0.00E+00  7.28E+01
 
 TH 6
+        9.99E+01 -9.69E+00 -6.69E-03  0.00E+00 -1.26E+00  3.18E+02
 
 TH 7
+       -4.09E+01 -6.91E+00  2.77E-03  0.00E+00 -2.24E+00  3.64E+01  4.55E+02
 
 TH 8
+        1.76E-01 -1.38E-02 -2.67E-06  0.00E+00  5.64E-02  1.85E-01 -3.23E-01  1.62E-03
 
 TH 9
+       -1.96E-01  3.11E-02  3.94E-05  0.00E+00 -8.17E-02  9.15E-01  4.81E-01 -2.88E-03  1.18E-03
 
 TH10
+        3.76E+00 -4.44E-01  4.03E-04  0.00E+00 -5.48E+00  2.11E+00 -5.04E+00 -1.95E-03 -5.34E-02  3.77E+01
 
 TH11
+       -1.49E+01 -1.54E+01 -2.95E-04  0.00E+00  1.12E+00  6.91E+00  1.20E+01 -6.75E-03 -6.01E-03  4.95E+00  1.42E+02
 
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
 #CPUT: Total CPU Time in Seconds,       66.747
Stop Time:
Sat Sep 25 01:57:35 CDT 2021
