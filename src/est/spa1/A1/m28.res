Wed Sep 29 21:58:26 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat28.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1800.82937143794        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2650E+02 -6.2428E+01 -1.3497E+01 -7.7768E+01  1.1044E+02  3.6252E+01 -3.7950E+00  7.2745E+00 -2.2735E+01 -1.3331E+01
            -5.4254E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1904.92853916156        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8281E-01  9.4102E-01  9.9066E-01  1.1354E+00  8.6553E-01  9.9495E-01  8.7364E-01  8.3612E-01  1.0992E+00  7.5944E-01
             1.8820E+00
 PARAMETER:  8.2657E-02  3.9209E-02  9.0613E-02  2.2694E-01 -4.4419E-02  9.4936E-02 -3.5090E-02 -7.8980E-02  1.9454E-01 -1.7517E-01
             7.3235E-01
 GRADIENT:   1.0005E+02  1.6536E+01  1.8044E+01  3.5358E+01 -2.1536E+01  1.6716E+01  6.0480E+00  5.2764E+00  1.9782E+01  4.1830E+00
             1.4391E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1912.66964849516        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      184
 NPARAMETR:  9.5905E-01  7.9838E-01  5.7280E-01  1.2142E+00  5.9769E-01  1.0489E+00  9.3378E-01  3.2816E-01  1.0187E+00  5.7842E-01
             1.7774E+00
 PARAMETER:  5.8187E-02 -1.2517E-01 -4.5722E-01  2.9405E-01 -4.1468E-01  1.4779E-01  3.1483E-02 -1.0143E+00  1.1848E-01 -4.4745E-01
             6.7515E-01
 GRADIENT:  -5.6881E+01  1.5215E+01 -3.3023E+01  6.3676E+01  2.1007E+01  2.2656E+01 -3.2980E+00  2.0159E+00  9.1763E+00 -1.3123E+00
             1.0845E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1927.26538546078        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      361
 NPARAMETR:  9.8005E-01  6.4017E-01  6.9502E-01  1.2717E+00  6.1059E-01  9.6927E-01  9.4819E-01  2.4170E-01  9.3153E-01  7.3716E-01
             1.5194E+00
 PARAMETER:  7.9847E-02 -3.4601E-01 -2.6381E-01  3.4034E-01 -3.9334E-01  6.8791E-02  4.6805E-02 -1.3201E+00  2.9077E-02 -2.0495E-01
             5.1830E-01
 GRADIENT:  -3.3069E+00  1.2012E+01  1.8813E+01 -9.5232E+00 -3.0484E+01 -4.2528E+00 -3.3754E+00  4.2639E-01 -5.2425E+00  7.6424E-01
            -4.0040E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1929.03445073063        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      536
 NPARAMETR:  9.7670E-01  4.5618E-01  7.0364E-01  1.3762E+00  5.7424E-01  9.7213E-01  1.3536E+00  1.1293E-01  8.8107E-01  7.3643E-01
             1.5323E+00
 PARAMETER:  7.6426E-02 -6.8487E-01 -2.5148E-01  4.1931E-01 -4.5472E-01  7.1732E-02  4.0274E-01 -2.0810E+00 -2.6619E-02 -2.0594E-01
             5.2676E-01
 GRADIENT:  -3.6647E+00  1.9344E+00 -7.0687E-01  6.1126E-01 -5.0152E-01 -2.1845E+00 -5.2207E-01 -3.1120E-02 -2.6822E-01 -4.2515E-02
             1.2317E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1929.37551799196        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  9.7432E-01  2.8450E-01  7.3119E-01  1.4733E+00  5.5296E-01  9.7574E-01  1.7413E+00  3.5860E-02  8.3769E-01  7.6219E-01
             1.5376E+00
 PARAMETER:  7.3979E-02 -1.1570E+00 -2.1308E-01  4.8748E-01 -4.9248E-01  7.5442E-02  6.5466E-01 -3.2281E+00 -7.7103E-02 -1.7156E-01
             5.3025E-01
 GRADIENT:   1.1704E+00  1.7175E+00  3.7619E+00  5.5615E+00 -4.9445E+00  5.3290E-01  6.0295E-02 -1.8118E-02 -3.5762E-01  1.3040E-01
            -9.1057E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1929.41333242960        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      892
 NPARAMETR:  9.7278E-01  2.4997E-01  7.2852E-01  1.4868E+00  5.4611E-01  9.7336E-01  1.8190E+00  2.9820E-02  8.3186E-01  7.6182E-01
             1.5408E+00
 PARAMETER:  7.2404E-02 -1.2864E+00 -2.1674E-01  4.9661E-01 -5.0493E-01  7.2996E-02  6.9830E-01 -3.4126E+00 -8.4097E-02 -1.7204E-01
             5.3231E-01
 GRADIENT:  -1.2672E-01  3.3145E-01  2.2877E-01 -9.0572E-01 -1.5633E-01 -1.2313E-01 -2.0950E-02 -1.2676E-02  4.7489E-02  1.1083E-02
             2.8074E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1929.91797900018        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  9.7060E-01  1.9280E-01  7.4268E-01  1.5200E+00  5.4055E-01  9.7469E-01  1.8250E+00  4.0608E-01  8.2406E-01  7.1997E-01
             1.5338E+00
 PARAMETER:  7.0156E-02 -1.5461E+00 -1.9749E-01  5.1869E-01 -5.1517E-01  7.4365E-02  7.0156E-01 -8.0122E-01 -9.3511E-02 -2.2854E-01
             5.2777E-01
 GRADIENT:  -9.9943E-01  7.0303E-01 -2.3453E+00 -1.0577E+00 -4.9570E+00  8.0611E-01 -1.9654E-01  3.3353E-01  1.1849E+00  3.2373E-01
            -3.9765E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1930.06746766161        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  9.6844E-01  1.1355E-01  7.7876E-01  1.5739E+00  5.4372E-01  9.6980E-01  2.8527E+00  4.7693E-01  7.9616E-01  7.0417E-01
             1.5376E+00
 PARAMETER:  6.7935E-02 -2.0755E+00 -1.5006E-01  5.5356E-01 -5.0932E-01  6.9332E-02  1.1483E+00 -6.4038E-01 -1.2796E-01 -2.5074E-01
             5.3024E-01
 GRADIENT:  -5.0602E-01  1.0130E+00  3.4903E+00  7.9230E+00 -8.7088E+00 -5.8242E-01  1.5230E-01 -2.7825E-01 -6.6893E-01 -2.4062E+00
             8.0896E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1930.18990773823        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1425
 NPARAMETR:  9.6701E-01  6.1515E-02  8.2271E-01  1.6090E+00  5.5758E-01  9.6951E-01  3.6021E+00  5.3422E-01  7.8310E-01  7.3079E-01
             1.5323E+00
 PARAMETER:  6.6450E-02 -2.6885E+00 -9.5154E-02  5.7563E-01 -4.8414E-01  6.9039E-02  1.3815E+00 -5.2694E-01 -1.4450E-01 -2.1363E-01
             5.2677E-01
 GRADIENT:  -3.3866E-02  3.8213E-01  1.5542E+00  3.7405E+00 -2.9783E+00 -3.6149E-01  7.9718E-02  5.3486E-03 -3.9374E-01  2.8887E-01
            -3.8673E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1930.20996747149        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1605
 NPARAMETR:  9.6689E-01  5.3720E-02  8.2626E-01  1.6102E+00  5.5846E-01  9.7020E-01  3.2676E+00  5.3861E-01  7.8278E-01  7.3133E-01
             1.5322E+00
 PARAMETER:  6.6333E-02 -2.8240E+00 -9.0842E-02  5.7636E-01 -4.8258E-01  6.9742E-02  1.2841E+00 -5.1876E-01 -1.4490E-01 -2.1289E-01
             5.2671E-01
 GRADIENT:   3.4057E-01  1.1548E-01  9.2885E-01 -3.9198E+00 -3.9938E-01 -1.1060E-02  9.8949E-03 -5.4935E-02 -2.3472E-02  4.2502E-02
            -6.6759E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1930.21030482165        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1781
 NPARAMETR:  9.6649E-01  4.4153E-02  8.2760E-01  1.6164E+00  5.5745E-01  9.7008E-01  3.2090E+00  5.4323E-01  7.8119E-01  7.3138E-01
             1.5325E+00
 PARAMETER:  6.5912E-02 -3.0201E+00 -8.9231E-02  5.8019E-01 -4.8439E-01  6.9623E-02  1.2660E+00 -5.1022E-01 -1.4694E-01 -2.1283E-01
             5.2689E-01
 GRADIENT:   3.7457E-02  1.1278E-01  6.2566E-01 -1.9121E+00 -5.0489E-01  3.5878E-03 -5.0929E-03 -5.0631E-02  8.2583E-02  1.0843E-01
             1.1425E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1930.21367216011        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1969
 NPARAMETR:  9.6658E-01  4.2092E-02  8.2771E-01  1.6157E+00  5.5701E-01  9.7007E-01  3.5279E+00  5.4506E-01  7.8038E-01  7.3073E-01
             1.5323E+00
 PARAMETER:  6.6007E-02 -3.0679E+00 -8.9095E-02  5.7976E-01 -4.8518E-01  6.9613E-02  1.3607E+00 -5.0685E-01 -1.4798E-01 -2.1371E-01
             5.2679E-01
 GRADIENT:   4.4658E-01  6.5605E-02  1.1724E+00 -5.6908E+00 -7.1954E-01  2.4895E-02  4.0456E-03 -2.8410E-02  4.9235E-02  9.9742E-02
            -8.4186E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1930.21448654316        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2152             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6653E-01  4.0783E-02  8.2735E-01  1.6163E+00  5.5662E-01  9.7004E-01  3.5877E+00  5.4602E-01  7.8006E-01  7.3008E-01
             1.5323E+00
 PARAMETER:  6.5961E-02 -3.0995E+00 -8.9532E-02  5.8016E-01 -4.8588E-01  6.9584E-02  1.3775E+00 -5.0510E-01 -1.4839E-01 -2.1460E-01
             5.2680E-01
 GRADIENT:   1.5156E+02  1.4897E+00  3.3179E+00  3.9800E+02  2.8290E+01  1.3099E+01  3.1160E-01  3.0755E-01  7.7117E+00  9.1114E-01
             5.2453E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1930.21508160241        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2326
 NPARAMETR:  9.6645E-01  3.9700E-02  8.2698E-01  1.6166E+00  5.5652E-01  9.7000E-01  3.6342E+00  5.4641E-01  7.7999E-01  7.2987E-01
             1.5324E+00
 PARAMETER:  6.5871E-02 -3.1264E+00 -8.9980E-02  5.8033E-01 -4.8606E-01  6.9539E-02  1.3904E+00 -5.0439E-01 -1.4848E-01 -2.1489E-01
             5.2682E-01
 GRADIENT:   3.3105E-01  4.1807E-02  2.4650E-01 -6.1001E+00  4.7180E-01  1.3873E-02  4.4305E-03  4.2072E-03  9.2010E-02  4.2031E-02
             2.3277E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1930.21574516065        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2505             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6657E-01  3.8492E-02  8.2632E-01  1.6169E+00  5.5608E-01  9.7005E-01  3.7217E+00  5.4669E-01  7.7972E-01  7.2936E-01
             1.5324E+00
 PARAMETER:  6.5994E-02 -3.1573E+00 -9.0773E-02  5.8052E-01 -4.8685E-01  6.9588E-02  1.4142E+00 -5.0387E-01 -1.4882E-01 -2.1559E-01
             5.2686E-01
 GRADIENT:   1.5178E+02  1.3756E+00  2.2353E+00  3.9791E+02  2.9846E+01  1.3109E+01  3.0770E-01  3.4099E-01  7.7911E+00  8.6993E-01
             5.2843E+00

0ITERATION NO.:   78    OBJECTIVE VALUE:  -1930.21582564730        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     2601
 NPARAMETR:  9.6652E-01  3.8311E-02  8.2633E-01  1.6174E+00  5.5605E-01  9.7003E-01  3.6695E+00  5.4661E-01  7.7969E-01  7.2934E-01
             1.5324E+00
 PARAMETER:  6.5948E-02 -3.1620E+00 -9.0767E-02  5.8080E-01 -4.8691E-01  6.9572E-02  1.4000E+00 -5.0402E-01 -1.4886E-01 -2.1562E-01
             5.2686E-01
 GRADIENT:   1.3196E-01  8.7853E-04 -9.8935E-02  2.9409E-01  9.3736E-01  1.1207E-02 -5.5909E-04  8.8820E-03  3.6215E-02  6.8899E-03
             4.2601E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2601
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4308E-04 -1.2985E-03 -5.6775E-03 -4.4260E-03 -1.3723E-02
 SE:             2.9701E-02  2.2869E-03  1.2409E-02  2.8940E-02  2.1840E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8541E-01  5.7017E-01  6.4728E-01  8.7845E-01  5.2978E-01

 ETASHRINKSD(%)  4.9932E-01  9.2338E+01  5.8429E+01  3.0470E+00  2.6832E+01
 ETASHRINKVR(%)  9.9615E-01  9.9413E+01  8.2719E+01  6.0012E+00  4.6465E+01
 EBVSHRINKSD(%)  7.8690E-01  9.2508E+01  5.8976E+01  3.2818E+00  2.6500E+01
 EBVSHRINKVR(%)  1.5676E+00  9.9439E+01  8.3170E+01  6.4559E+00  4.5978E+01
 RELATIVEINF(%)  9.2514E+01  3.3844E-02  1.7560E+00  8.8459E+00  3.9643E+00
 EPSSHRINKSD(%)  3.0550E+01
 EPSSHRINKVR(%)  5.1766E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1930.2158256472953     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1011.2772924426226     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.83
 Elapsed covariance  time in seconds:     6.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1930.216       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.67E-01  3.83E-02  8.26E-01  1.62E+00  5.56E-01  9.70E-01  3.67E+00  5.47E-01  7.80E-01  7.29E-01  1.53E+00
 


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
 
         2.91E-02  3.65E-01  2.10E-01  1.98E-01  1.12E-01  6.11E-02  1.72E+01  3.36E-01  1.00E-01  1.23E-01  9.26E-02
 


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
+        8.48E-04
 
 TH 2
+        2.68E-03  1.33E-01
 
 TH 3
+       -1.88E-03 -1.08E-02  4.40E-02
 
 TH 4
+       -1.75E-03 -6.82E-02  1.62E-02  3.91E-02
 
 TH 5
+       -5.07E-04  1.82E-02  1.91E-02 -4.23E-03  1.26E-02
 
 TH 6
+       -3.18E-04  1.95E-03 -3.13E-03 -1.98E-03 -1.27E-03  3.73E-03
 
 TH 7
+       -1.16E-01 -6.23E+00  4.56E-01  3.18E+00 -8.76E-01 -1.15E-01  2.97E+02
 
 TH 8
+       -2.74E-03 -2.13E-02  5.75E-02  2.57E-02  2.38E-02 -5.29E-03  1.07E+00  1.13E-01
 
 TH 9
+        8.06E-04  2.83E-02 -8.48E-03 -1.61E-02  9.46E-04  2.03E-04 -1.27E+00 -1.04E-02  1.01E-02
 
 TH10
+        1.78E-04  2.62E-03  1.95E-03 -1.45E-03  1.67E-03  7.34E-04 -1.14E-01 -9.83E-03 -5.58E-04  1.51E-02
 
 TH11
+        7.06E-04 -6.60E-03 -5.54E-03  1.57E-03 -4.07E-03  3.80E-04  2.75E-01 -1.02E-02 -1.17E-03 -6.30E-04  8.57E-03
 
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
+        2.91E-02
 
 TH 2
+        2.52E-01  3.65E-01
 
 TH 3
+       -3.08E-01 -1.41E-01  2.10E-01
 
 TH 4
+       -3.03E-01 -9.46E-01  3.90E-01  1.98E-01
 
 TH 5
+       -1.55E-01  4.45E-01  8.12E-01 -1.91E-01  1.12E-01
 
 TH 6
+       -1.78E-01  8.76E-02 -2.44E-01 -1.64E-01 -1.85E-01  6.11E-02
 
 TH 7
+       -2.31E-01 -9.90E-01  1.26E-01  9.32E-01 -4.53E-01 -1.09E-01  1.72E+01
 
 TH 8
+       -2.79E-01 -1.73E-01  8.14E-01  3.87E-01  6.29E-01 -2.57E-01  1.84E-01  3.36E-01
 
 TH 9
+        2.75E-01  7.70E-01 -4.02E-01 -8.08E-01  8.39E-02  3.30E-02 -7.32E-01 -3.07E-01  1.00E-01
 
 TH10
+        4.98E-02  5.83E-02  7.55E-02 -5.95E-02  1.21E-01  9.77E-02 -5.35E-02 -2.38E-01 -4.51E-02  1.23E-01
 
 TH11
+        2.62E-01 -1.95E-01 -2.85E-01  8.57E-02 -3.91E-01  6.71E-02  1.72E-01 -3.29E-01 -1.26E-01 -5.54E-02  9.26E-02
 
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
+        1.76E+03
 
 TH 2
+       -2.92E+02  9.59E+02
 
 TH 3
+       -1.74E+02  3.00E+02  1.00E+03
 
 TH 4
+       -7.01E+01  4.12E+02 -9.12E+01  7.07E+02
 
 TH 5
+        5.78E+02 -7.38E+02 -1.93E+03 -1.09E+02  4.32E+03
 
 TH 6
+        2.55E+02  1.19E+01 -7.62E+01  3.94E+01  2.49E+02  3.54E+02
 
 TH 7
+       -2.29E+00  1.23E+01  3.73E-01  1.06E+00  1.66E+00  1.14E+00  2.44E-01
 
 TH 8
+       -9.59E+00 -4.84E+01 -3.83E+01 -1.67E+01 -3.56E+01 -7.32E+00 -1.06E+00  4.42E+01
 
 TH 9
+        4.13E+01 -1.90E+02  5.63E+01  2.80E+01 -7.87E+00  4.77E+01 -2.77E+00 -2.91E+00  3.70E+02
 
 TH10
+       -5.86E+01 -2.51E+01  6.81E+00  1.84E+01 -1.38E+02 -3.10E+01 -1.21E+00  3.74E+01  1.79E+01  1.05E+02
 
 TH11
+       -1.43E+02  5.17E+01 -4.94E+01  3.31E+01  8.80E+01 -6.83E+00  8.65E-01  1.44E+01  1.32E+01  1.60E+01  1.65E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       46.623
Stop Time:
Wed Sep 29 21:59:15 CDT 2021
