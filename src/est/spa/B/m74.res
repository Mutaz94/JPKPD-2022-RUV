Sat Sep 18 08:42:17 CDT 2021
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
$DATA ../../../../data/spa/B/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.60143315052        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4971E+01 -5.4029E+01 -1.6350E+01 -7.0355E+01 -1.6495E+01  2.3477E+01 -1.0435E+01  1.4260E+01 -5.6602E+00  1.6149E+01
             2.2676E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1698.46931307749        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  9.9215E-01  1.0718E+00  1.1243E+00  1.0142E+00  1.0983E+00  9.4289E-01  1.1055E+00  8.8500E-01  1.0260E+00  9.6129E-01
             9.8382E-01
 PARAMETER:  9.2115E-02  1.6934E-01  2.1717E-01  1.1414E-01  1.9373E-01  4.1196E-02  2.0030E-01 -2.2172E-02  1.2565E-01  6.0523E-02
             8.3689E-02
 GRADIENT:  -3.5400E+00  2.1586E+00  7.8651E-01  1.1972E+00  4.5848E+00 -2.1638E+00 -6.0532E-01  3.5745E+00  2.0201E+00 -7.4135E+00
             9.2769E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1699.28523557231        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      290
 NPARAMETR:  9.9143E-01  8.8845E-01  1.1183E+00  1.1335E+00  1.0121E+00  9.5670E-01  1.2173E+00  5.5955E-01  9.9806E-01  1.0018E+00
             9.5552E-01
 PARAMETER:  9.1395E-02 -1.8277E-02  2.1177E-01  2.2533E-01  1.1203E-01  5.5737E-02  2.9662E-01 -4.8063E-01  9.8063E-02  1.0175E-01
             5.4498E-02
 GRADIENT:  -1.6659E+00  4.9084E+00  8.1716E+00  1.1599E+01 -7.1441E+00  4.0087E+00 -2.3617E+00 -1.4695E+00  1.1047E+01  7.7541E-01
            -2.7589E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.15322469015        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      466
 NPARAMETR:  9.9274E-01  8.2903E-01  1.0244E+00  1.1566E+00  9.4129E-01  9.4565E-01  1.3919E+00  4.7703E-01  8.9866E-01  9.2693E-01
             9.5703E-01
 PARAMETER:  9.2715E-02 -8.7495E-02  1.2409E-01  2.4547E-01  3.9495E-02  4.4117E-02  4.3066E-01 -6.4017E-01 -6.8502E-03  2.4124E-02
             5.6079E-02
 GRADIENT:   7.3875E-01  2.7312E+00  2.6806E+00  4.1497E-01 -2.8020E+00 -4.8906E-01  5.3641E-02 -3.6074E-01 -1.2859E+00 -1.6897E-01
            -2.7035E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1700.47551638235        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      642
 NPARAMETR:  9.8846E-01  5.4980E-01  1.1472E+00  1.3339E+00  8.9643E-01  9.4300E-01  1.8252E+00  6.4974E-01  8.3368E-01  9.2387E-01
             9.5578E-01
 PARAMETER:  8.8389E-02 -4.9820E-01  2.3731E-01  3.8814E-01 -9.3367E-03  4.1313E-02  7.0167E-01 -3.3119E-01 -8.1903E-02  2.0815E-02
             5.4778E-02
 GRADIENT:  -4.4409E-01  2.7313E+00 -3.8596E-01  6.6141E+00 -2.6163E+00  2.3468E-02  3.1641E-03  7.6704E-01  7.2776E-02  1.9490E-01
             1.2013E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1700.71546352777        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  9.8496E-01  3.5062E-01  1.2206E+00  1.4569E+00  8.6899E-01  9.4050E-01  2.3974E+00  6.9817E-01  7.9732E-01  9.3374E-01
             9.5582E-01
 PARAMETER:  8.4842E-02 -9.4804E-01  2.9936E-01  4.7634E-01 -4.0422E-02  3.8655E-02  9.7438E-01 -2.5930E-01 -1.2650E-01  3.1444E-02
             5.4818E-02
 GRADIENT:  -6.4796E-01  2.4163E+00  1.2214E+00  1.0324E+01 -2.6715E+00  3.6682E-01  3.2961E-01 -1.0583E-01 -1.5931E-01 -1.7007E-02
            -3.7823E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1700.99900114288        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      994
 NPARAMETR:  9.8228E-01  1.8988E-01  1.3008E+00  1.5539E+00  8.5752E-01  9.3730E-01  3.3019E+00  7.8276E-01  7.7018E-01  9.5249E-01
             9.5504E-01
 PARAMETER:  8.2119E-02 -1.5613E+00  3.6297E-01  5.4079E-01 -5.3712E-02  3.5244E-02  1.2945E+00 -1.4493E-01 -1.6113E-01  5.1322E-02
             5.3998E-02
 GRADIENT:   3.5464E-01  1.3236E+00  3.2578E+00  7.5037E+00 -5.9914E+00  1.6872E-01 -2.3139E-01 -3.3557E-01 -1.6282E+00  1.5335E+00
            -6.5429E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.44479304253        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  9.7914E-01  6.4236E-02  1.5170E+00  1.6444E+00  9.0564E-01  9.3559E-01  5.4671E+00  1.0666E+00  7.4751E-01  9.7126E-01
             9.5714E-01
 PARAMETER:  7.8915E-02 -2.6452E+00  5.1671E-01  5.9740E-01  8.8815E-04  3.3427E-02  1.7987E+00  1.6451E-01 -1.9101E-01  7.0841E-02
             5.6197E-02
 GRADIENT:  -1.9322E+00  3.2926E+00 -2.9378E+00 -1.7312E+00  7.8899E-01  2.0963E-01  8.2424E+00  6.3373E-01 -2.4765E+00 -4.0200E+00
             7.7864E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.77826051285        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1348
 NPARAMETR:  9.7960E-01  4.1813E-02  1.5521E+00  1.6588E+00  9.1453E-01  9.3429E-01  6.4491E+00  1.0876E+00  7.4146E-01  9.9042E-01
             9.5629E-01
 PARAMETER:  7.9393E-02 -3.0745E+00  5.3962E-01  6.0610E-01  1.0659E-02  3.2026E-02  1.9639E+00  1.8393E-01 -1.9913E-01  9.0377E-02
             5.5306E-02
 GRADIENT:   4.8102E-02 -2.0024E+00  1.3355E+00  2.4443E+01 -5.6015E+00 -3.8734E-02 -7.3588E+00  2.5625E+00  2.0089E+00  4.0988E+00
             6.4601E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1703.21036954038        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1530
 NPARAMETR:  9.7857E-01  1.0000E-02  1.3565E+00  1.6545E+00  8.3910E-01  9.3404E-01  1.2215E+01  8.6248E-01  7.5030E-01  9.4187E-01
             9.5588E-01
 PARAMETER:  7.8339E-02 -4.7554E+00  4.0493E-01  6.0352E-01 -7.5426E-02  3.1760E-02  2.6026E+00 -4.7944E-02 -1.8728E-01  4.0110E-02
             5.4878E-02
 GRADIENT:  -3.3662E-01  0.0000E+00  1.5727E-01  1.0507E+01 -1.0628E+00  7.4415E-03 -1.6964E+00  1.6112E-01  5.8068E+00  2.8758E+00
             4.3589E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1703.27542896836        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1708
 NPARAMETR:  9.7874E-01  1.0000E-02  1.3427E+00  1.6480E+00  8.3240E-01  9.3393E-01  1.2019E+01  8.6276E-01  7.4228E-01  9.2287E-01
             9.5656E-01
 PARAMETER:  7.8507E-02 -4.7522E+00  3.9466E-01  5.9959E-01 -8.3440E-02  3.1650E-02  2.5865E+00 -4.7616E-02 -1.9803E-01  1.9734E-02
             5.5591E-02
 GRADIENT:   8.0878E-02  0.0000E+00  5.0317E-01 -1.6482E+00 -6.1169E-01 -3.3887E-02 -1.3987E+00  2.2009E-01  1.6023E+00  9.6347E-01
             6.5465E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1703.27893149951        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1883
 NPARAMETR:  9.7873E-01  1.0000E-02  1.3441E+00  1.6497E+00  8.3282E-01  9.3406E-01  1.2015E+01  8.6500E-01  7.4004E-01  9.2051E-01
             9.5706E-01
 PARAMETER:  7.8500E-02 -4.7536E+00  3.9572E-01  6.0057E-01 -8.2932E-02  3.1784E-02  2.5861E+00 -4.5027E-02 -2.0106E-01  1.7173E-02
             5.6107E-02
 GRADIENT:   3.5340E-03  0.0000E+00 -7.5022E-02  2.9594E-03  9.5614E-02 -9.5941E-03  1.5459E-01 -2.3416E-02 -2.2198E-02 -3.3098E-02
            -9.8261E-03

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1703.27893149951        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1912
 NPARAMETR:  9.7872E-01  1.0000E-02  1.3446E+00  1.6506E+00  8.3278E-01  9.3403E-01  1.1984E+01  8.6504E-01  7.4010E-01  9.2055E-01
             9.5706E-01
 PARAMETER:  7.8500E-02 -4.7536E+00  3.9572E-01  6.0057E-01 -8.2932E-02  3.1784E-02  2.5861E+00 -4.5027E-02 -2.0106E-01  1.7173E-02
             5.6107E-02
 GRADIENT:   8.8390E-03  0.0000E+00 -2.1418E+03 -1.4056E+03  1.5390E-01  1.0933E-02  3.2440E+02 -2.6050E-02 -2.0594E-01 -8.6951E-02
             3.9824E-06

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1912
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2736E-04  9.4584E-03 -2.1540E-02 -8.0595E-03 -2.8735E-02
 SE:             2.9833E-02  6.3051E-03  1.4921E-02  2.8938E-02  2.1830E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9659E-01  1.3358E-01  1.4885E-01  7.8062E-01  1.8806E-01

 ETASHRINKSD(%)  5.7133E-02  7.8877E+01  5.0013E+01  3.0556E+00  2.6868E+01
 ETASHRINKVR(%)  1.1423E-01  9.5538E+01  7.5013E+01  6.0179E+00  4.6517E+01
 EBVSHRINKSD(%)  4.3225E-01  8.4425E+01  5.1722E+01  3.0664E+00  2.3834E+01
 EBVSHRINKVR(%)  8.6263E-01  9.7574E+01  7.6693E+01  6.0388E+00  4.1988E+01
 RELATIVEINF(%)  9.9086E+01  1.9848E+00  3.7190E+00  7.0824E+01  9.2595E+00
 EPSSHRINKSD(%)  4.3902E+01
 EPSSHRINKVR(%)  6.8530E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1703.2789314995148     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -968.12810493577661     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1703.279       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.00E-02  1.34E+00  1.65E+00  8.33E-01  9.34E-01  1.20E+01  8.65E-01  7.40E-01  9.21E-01  9.57E-01
 


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
+        1.39E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        3.59E+01  0.00E+00  7.49E+05
 
 TH 4
+        7.69E+00  0.00E+00 -3.43E+02  2.16E+05
 
 TH 5
+        7.28E+01  0.00E+00 -4.78E+06  7.72E+02  6.27E+03
 
 TH 6
+       -1.91E+01  0.00E+00 -2.30E+01 -1.57E+01 -8.16E+01  3.05E+02
 
 TH 7
+       -4.42E-01  0.00E+00  8.64E+00  3.13E+01 -2.38E+01  3.38E-01  2.19E+02
 
 TH 8
+       -1.34E+01  0.00E+00 -4.46E+02 -1.81E+02 -9.37E+02  1.43E+01  4.96E+00  3.72E+02
 
 TH 9
+       -9.95E+01  0.00E+00  2.68E+06 -1.82E+02 -1.71E+07  1.04E+02  3.21E+00  1.71E+03  9.57E+06
 
 TH10
+       -1.46E+02  0.00E+00  4.33E+06 -5.88E+02 -4.28E+03  4.10E+01  1.63E+01  8.27E+02  1.55E+07  3.46E+03
 
 TH11
+        2.78E+01  0.00E+00 -2.96E+02 -1.29E+02 -6.52E+02  5.87E+00  3.27E+00  2.30E+02  1.08E+03  5.83E+02  3.84E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.673
Stop Time:
Sat Sep 18 08:42:52 CDT 2021
