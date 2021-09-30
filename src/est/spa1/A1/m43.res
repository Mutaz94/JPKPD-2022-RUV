Wed Sep 29 22:27:16 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1852.99304091481        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4077E+02  2.9798E+01  2.1517E+01  6.1183E+01  7.7388E+01  5.2655E+01  4.5027E+00  2.8307E+00  3.7704E+01 -4.4724E+01
            -4.5418E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1917.98202825912        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0077E+00  9.6017E-01  8.7397E-01  1.0271E+00  8.8020E-01  9.4911E-01  9.4327E-01  9.4352E-01  8.1700E-01  1.1478E+00
             1.8235E+00
 PARAMETER:  1.0769E-01  5.9353E-02 -3.4706E-02  1.2677E-01 -2.7609E-02  4.7770E-02  4.1592E-02  4.1860E-02 -1.0211E-01  2.3786E-01
             7.0077E-01
 GRADIENT:   1.7693E+02  2.4146E+00 -2.3073E+00  8.7948E-01 -8.1964E+00  9.1583E+00  3.8685E+00  1.2171E+01  3.4328E+00  1.7258E+01
             1.6610E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1929.85775196879        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      201
 NPARAMETR:  9.9165E-01  8.3352E-01  5.9925E-01  1.0829E+00  6.6146E-01  1.0455E+00  1.0019E+00  3.4401E-01  7.3202E-01  9.3172E-01
             1.7499E+00
 PARAMETER:  9.1614E-02 -8.2096E-02 -4.1207E-01  1.7964E-01 -3.1331E-01  1.4451E-01  1.0190E-01 -9.6708E-01 -2.1195E-01  2.9274E-02
             6.5958E-01
 GRADIENT:   3.7398E+00  1.3551E+01 -6.2252E+00  9.0731E+00 -6.4403E+00  3.0983E+01 -3.6782E+00  2.5699E+00 -1.4834E+01  1.8298E+01
             1.4202E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1950.50385949156        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  9.7515E-01  6.8153E-01  5.3872E-01  1.1394E+00  5.5973E-01  9.2426E-01  1.2699E+00  1.3088E-01  7.4875E-01  7.8281E-01
             1.4285E+00
 PARAMETER:  7.4835E-02 -2.8341E-01 -5.1856E-01  2.3049E-01 -4.8030E-01  2.1233E-02  3.3893E-01 -1.9335E+00 -1.8936E-01 -1.4486E-01
             4.5665E-01
 GRADIENT:  -2.8798E+01  2.0283E+01  2.2711E+01 -3.4555E+00 -4.1877E+01 -1.6353E+01 -2.6881E-02  3.3085E-01 -5.0848E+00  6.1170E+00
             6.5910E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1952.00440311477        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  9.8468E-01  5.6887E-01  5.3563E-01  1.1878E+00  5.3131E-01  9.5899E-01  1.4795E+00  8.1683E-02  7.4336E-01  7.1888E-01
             1.4362E+00
 PARAMETER:  8.4565E-02 -4.6411E-01 -5.2430E-01  2.7208E-01 -5.3240E-01  5.8122E-02  4.9174E-01 -2.4049E+00 -1.9658E-01 -2.3006E-01
             4.6198E-01
 GRADIENT:   2.2341E-01  2.8565E+00  9.5947E+00 -7.3089E+00 -1.5027E+01 -5.5770E-01 -7.2024E-02  9.7906E-02  1.4241E+00 -2.0517E+00
             9.6143E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1952.65558635256        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  9.8063E-01  4.0283E-01  6.9437E-01  1.3106E+00  5.8723E-01  9.6022E-01  2.0272E+00  1.3043E-01  6.8314E-01  8.7564E-01
             1.4487E+00
 PARAMETER:  8.0443E-02 -8.0924E-01 -2.6475E-01  3.7045E-01 -4.3233E-01  5.9412E-02  8.0664E-01 -1.9369E+00 -2.8105E-01 -3.2798E-02
             4.7065E-01
 GRADIENT:   4.3522E+00  3.8858E+00  7.6172E+00 -3.1649E+00 -1.6048E+01  1.9728E+00  9.4563E-01  3.7706E-02  3.3955E-01  4.8170E+00
             2.0681E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1952.86352582099        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  9.7655E-01  3.4491E-01  7.4607E-01  1.3545E+00  6.0839E-01  9.5389E-01  2.2472E+00  1.2498E-01  6.6580E-01  8.9962E-01
             1.4493E+00
 PARAMETER:  7.6269E-02 -9.6449E-01 -1.9293E-01  4.0344E-01 -3.9694E-01  5.2794E-02  9.0969E-01 -1.9796E+00 -3.0676E-01 -5.7864E-03
             4.7106E-01
 GRADIENT:  -8.5705E-01  2.9784E+00  1.9717E+00  1.0902E+01 -2.5577E+00  7.0600E-03 -4.9560E-02 -8.4094E-02 -6.8490E-01 -5.6817E-03
            -5.9391E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1952.91035664008        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  9.7496E-01  2.8232E-01  7.5707E-01  1.3856E+00  6.0148E-01  9.5240E-01  2.5442E+00  4.8911E-02  6.5787E-01  9.1114E-01
             1.4550E+00
 PARAMETER:  7.4644E-02 -1.1647E+00 -1.7829E-01  4.2613E-01 -4.0836E-01  5.1232E-02  1.0338E+00 -2.9178E+00 -3.1875E-01  6.9377E-03
             4.7502E-01
 GRADIENT:  -4.4745E-01  1.2043E+00  9.0173E-01  3.9611E+00 -2.2204E+00 -7.0072E-02  3.2588E-01 -2.4097E-02 -2.3778E-01  3.6385E-01
             8.7713E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1953.29105911210        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1274
 NPARAMETR:  9.7356E-01  2.5749E-01  7.6320E-01  1.4007E+00  5.9699E-01  9.5188E-01  2.6369E+00  4.1677E-01  6.5585E-01  8.8505E-01
             1.4503E+00
 PARAMETER:  7.3205E-02 -1.2568E+00 -1.7024E-01  4.3695E-01 -4.1585E-01  5.0686E-02  1.0696E+00 -7.7521E-01 -3.2182E-01 -2.2115E-02
             4.7177E-01
 GRADIENT:  -2.2594E+00  9.1292E-01 -4.0269E+00  6.6875E+00 -4.5529E-01 -1.5932E-01 -1.9031E+00  4.3171E-01  2.6923E+00  2.5242E+00
             2.1615E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1953.97253272243        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1452
 NPARAMETR:  9.7099E-01  1.4687E-01  8.0765E-01  1.4644E+00  5.9481E-01  9.5025E-01  3.8125E+00  5.2912E-01  6.2103E-01  8.7256E-01
             1.4487E+00
 PARAMETER:  7.0563E-02 -1.8182E+00 -1.1362E-01  4.8144E-01 -4.1951E-01  4.8968E-02  1.4383E+00 -5.3654E-01 -3.7637E-01 -3.6321E-02
             4.7064E-01
 GRADIENT:  -3.0319E-01  2.9336E+00  2.4860E+00  3.2530E+00 -6.9046E+00 -8.3061E-02  3.5662E+00 -7.4683E-01  6.7489E-02 -2.5519E+00
            -9.6539E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1954.90748946928        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1628
 NPARAMETR:  9.6820E-01  5.8645E-02  8.8616E-01  1.5253E+00  6.1721E-01  9.4850E-01  5.5330E+00  6.7562E-01  5.9361E-01  9.2035E-01
             1.4491E+00
 PARAMETER:  6.7681E-02 -2.7363E+00 -2.0855E-02  5.2219E-01 -3.8255E-01  4.7121E-02  1.8107E+00 -2.9213E-01 -4.2153E-01  1.6995E-02
             4.7096E-01
 GRADIENT:  -1.2205E+00  1.0958E+00 -8.4010E-01  1.5961E+01 -6.5842E+00 -7.0099E-02  1.5672E+00  1.7317E-01 -1.3990E+00  8.2885E-01
             1.0327E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1955.61271628223        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1807
 NPARAMETR:  9.6810E-01  3.5290E-02  9.4397E-01  1.5390E+00  6.4256E-01  9.4748E-01  6.7617E+00  7.1243E-01  5.8083E-01  9.5438E-01
             1.4433E+00
 PARAMETER:  6.7579E-02 -3.2442E+00  4.2338E-02  5.3115E-01 -3.4230E-01  4.6048E-02  2.0113E+00 -2.3908E-01 -4.4330E-01  5.3307E-02
             4.6691E-01
 GRADIENT:  -4.8035E-01 -3.3482E+00  3.8868E+00  1.5818E+01 -8.9752E+00 -1.7388E-01 -6.5574E+00  1.0171E-01  2.3767E+00  3.9762E+00
             4.0979E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1955.67520383302        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1989
 NPARAMETR:  9.6851E-01  3.5361E-02  9.4342E-01  1.5388E+00  6.4421E-01  9.4852E-01  6.7634E+00  7.3302E-01  5.8117E-01  9.3243E-01
             1.4446E+00
 PARAMETER:  6.8005E-02 -3.2422E+00  4.1760E-02  5.3101E-01 -3.3973E-01  4.7147E-02  2.0115E+00 -2.1058E-01 -4.4272E-01  3.0035E-02
             4.6784E-01
 GRADIENT:   7.6559E-01 -1.4006E+00  2.8174E-02  6.4182E+00  6.9162E-01  2.1140E-01 -2.0762E+00  1.7575E-01  6.2106E-01 -4.3602E-01
             5.3406E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1955.69962982264        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2151
 NPARAMETR:  9.6854E-01  3.5568E-02  9.4146E-01  1.5384E+00  6.4303E-01  9.4817E-01  6.7757E+00  7.2273E-01  5.8324E-01  9.3866E-01
             1.4443E+00
 PARAMETER:  6.8037E-02 -3.2363E+00  3.9679E-02  5.3072E-01 -3.4156E-01  4.6781E-02  2.0133E+00 -2.2472E-01 -4.3915E-01  3.6700E-02
             4.6763E-01
 GRADIENT:   1.0305E+00 -2.6026E-01  1.7724E-01  1.5066E+00  3.2491E-01  3.9532E-02  6.1049E-01  1.7910E-02 -2.6133E-02  8.1864E-02
             1.3357E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2151
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4777E-04  2.9518E-02 -1.4992E-02 -1.9649E-02 -1.4084E-02
 SE:             2.9755E-02  1.2221E-02  1.3449E-02  2.6818E-02  2.2019E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7727E-01  1.5722E-02  2.6495E-01  4.6375E-01  5.2241E-01

 ETASHRINKSD(%)  3.1702E-01  5.9057E+01  5.4945E+01  1.0155E+01  2.6232E+01
 ETASHRINKVR(%)  6.3304E-01  8.3237E+01  7.9701E+01  1.9279E+01  4.5583E+01
 EBVSHRINKSD(%)  6.9515E-01  7.2639E+01  5.6406E+01  7.7545E+00  2.3424E+01
 EBVSHRINKVR(%)  1.3855E+00  9.2514E+01  8.0996E+01  1.4908E+01  4.1361E+01
 RELATIVEINF(%)  9.8451E+01  4.7839E+00  2.6721E+00  4.1721E+01  8.6184E+00
 EPSSHRINKSD(%)  3.1420E+01
 EPSSHRINKVR(%)  5.2968E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1955.6996298226400     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1036.7610966179673     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1955.700       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.69E-01  3.56E-02  9.41E-01  1.54E+00  6.43E-01  9.48E-01  6.78E+00  7.23E-01  5.83E-01  9.39E-01  1.44E+00
 


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
+        1.30E+03
 
 TH 2
+        1.23E+02  1.86E+05
 
 TH 3
+       -2.34E+01 -8.17E+02  5.23E+02
 
 TH 4
+       -3.51E+01 -1.44E+04 -5.92E+01  5.16E+03
 
 TH 5
+        1.54E+02  2.52E+03 -1.32E+03 -3.92E+02  5.17E+04
 
 TH 6
+       -1.06E+00 -8.34E+00  4.65E+00 -4.17E+00 -3.67E+01  2.16E+02
 
 TH 7
+        1.52E+00  3.49E+02 -7.07E+00 -8.28E+01  2.29E+01  1.90E-01  1.88E+01
 
 TH 8
+       -3.20E+00 -1.75E+02 -2.37E+01  1.60E+01 -7.13E+01  6.43E-01 -1.79E+00  2.56E+01
 
 TH 9
+       -3.82E+01 -1.27E+03  2.64E+02  5.09E+01 -4.24E+04  5.71E+00 -6.93E+00  5.37E+01  3.65E+04
 
 TH10
+       -1.58E+01 -7.47E+02  3.53E+01  7.34E+01 -5.68E+02  3.28E+00 -7.16E+00  3.17E+01  1.78E+02  1.69E+02
 
 TH11
+       -1.94E+01 -2.76E+02  2.33E+01  3.48E+00 -8.43E+03  3.52E+00 -3.42E+00  1.64E+01  7.13E+03  4.00E+01  2.20E+02
 
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
 #CPUT: Total CPU Time in Seconds,       42.577
Stop Time:
Wed Sep 29 22:28:00 CDT 2021
