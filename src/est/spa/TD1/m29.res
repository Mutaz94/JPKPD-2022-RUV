Wed Sep 29 18:03:36 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1618.58105038033        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6889E+02 -2.6798E+01 -4.8426E+01  2.5362E+01  6.6876E+01  2.3533E+01 -2.0254E+01  7.1996E+00 -2.3807E+01  1.7495E+01
            -3.4422E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1629.49710611653        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8573E-01  1.0431E+00  1.0974E+00  1.0138E+00  1.0142E+00  1.0901E+00  1.0965E+00  9.7328E-01  1.1330E+00  9.1208E-01
             1.1122E+00
 PARAMETER:  8.5623E-02  1.4221E-01  1.9296E-01  1.1366E-01  1.1409E-01  1.8628E-01  1.9209E-01  7.2917E-02  2.2491E-01  7.9690E-03
             2.0632E-01
 GRADIENT:  -1.4197E+00 -1.0768E+01 -1.3427E+01  8.8421E-01  1.4776E+01  8.4176E+00 -6.6610E+00  2.1114E+00  4.1325E-01  4.8657E+00
             1.0801E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1630.33890257633        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8647E-01  9.8773E-01  1.0772E+00  1.0443E+00  9.7383E-01  1.0643E+00  1.2897E+00  9.4643E-01  1.0825E+00  8.4229E-01
             1.0898E+00
 PARAMETER:  8.6376E-02  8.7651E-02  1.7435E-01  1.4336E-01  7.3478E-02  1.6229E-01  3.5438E-01  4.4943E-02  1.7932E-01 -7.1625E-02
             1.8598E-01
 GRADIENT:   7.9685E-01 -5.2077E+00 -7.1663E+00 -3.7483E+00  6.7160E+00 -1.0580E+00  1.9332E+00  2.0260E+00  2.4447E+00  9.1262E-01
             2.4228E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1630.64044137183        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.8437E-01  8.5212E-01  1.2326E+00  1.1463E+00  9.7283E-01  1.0661E+00  1.3476E+00  1.0013E+00  1.0334E+00  8.6938E-01
             1.0855E+00
 PARAMETER:  8.4248E-02 -6.0028E-02  3.0910E-01  2.3658E-01  7.2452E-02  1.6401E-01  3.9832E-01  1.0133E-01  1.3289E-01 -3.9980E-02
             1.8208E-01
 GRADIENT:  -2.8528E-01  5.4045E+00  2.2218E+00  4.9618E+00 -4.8187E+00  1.4293E-01 -2.1928E-01 -2.0951E-02 -1.7380E-01 -1.0532E-02
            -7.5939E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1630.87883466369        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.8178E-01  5.8150E-01  1.4552E+00  1.3297E+00  9.5828E-01  1.0623E+00  1.5993E+00  1.1350E+00  9.5102E-01  8.8619E-01
             1.0905E+00
 PARAMETER:  8.1611E-02 -4.4214E-01  4.7512E-01  3.8494E-01  5.7388E-02  1.6044E-01  5.6959E-01  2.2663E-01  4.9778E-02 -2.0829E-02
             1.8660E-01
 GRADIENT:   6.2507E-01  6.1823E+00  2.9349E+00  1.1394E+01 -4.1448E+00 -2.3510E-01 -5.6439E-02 -3.3150E-01 -1.2796E+00 -8.2734E-01
             5.2084E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1631.03811989449        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.7924E-01  3.9657E-01  1.5939E+00  1.4541E+00  9.4616E-01  1.0603E+00  1.8851E+00  1.2278E+00  9.0401E-01  8.9928E-01
             1.0903E+00
 PARAMETER:  7.9025E-02 -8.2491E-01  5.6618E-01  4.7440E-01  4.4658E-02  1.5860E-01  7.3398E-01  3.0522E-01 -9.1723E-04 -6.1631E-03
             1.8647E-01
 GRADIENT:   1.6184E-01  5.6828E+00  3.1529E+00  1.6371E+01 -5.5868E+00 -3.1325E-01 -1.6000E-01 -4.8588E-01 -1.2440E+00 -2.7748E-01
             3.5378E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1631.19322264945        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.7711E-01  2.6216E-01  1.6842E+00  1.5415E+00  9.3520E-01  1.0597E+00  2.2565E+00  1.2949E+00  8.7401E-01  9.0746E-01
             1.0880E+00
 PARAMETER:  7.6843E-02 -1.2388E+00  6.2127E-01  5.3274E-01  3.3001E-02  1.5797E-01  9.1381E-01  3.5840E-01 -3.4668E-02  2.8891E-03
             1.8430E-01
 GRADIENT:  -5.8273E-01  3.9845E+00  2.8060E+00  1.5457E+01 -6.7847E+00 -1.3447E-01 -2.4905E-01 -3.3233E-01 -5.5169E-01  6.0297E-01
            -3.7807E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1631.42902815043        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1238             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7646E-01  1.3920E-01  1.7616E+00  1.5951E+00  9.3480E-01  1.0597E+00  3.0147E+00  1.3612E+00  8.5154E-01  9.0529E-01
             1.0881E+00
 PARAMETER:  7.6180E-02 -1.8718E+00  6.6622E-01  5.6694E-01  3.2573E-02  1.5800E-01  1.2035E+00  4.0838E-01 -6.0706E-02  5.0008E-04
             1.8444E-01
 GRADIENT:   3.6822E+02  1.9271E+01  4.4386E+00  8.7393E+02  1.5920E+01  6.9394E+01  4.8278E+00  1.0276E+00  1.6970E+01 -2.2710E-01
             1.4997E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1631.51496393300        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1421
 NPARAMETR:  9.7622E-01  1.3965E-01  1.7670E+00  1.6047E+00  9.3001E-01  1.0596E+00  3.1426E+00  1.3623E+00  8.4516E-01  9.0532E-01
             1.0883E+00
 PARAMETER:  7.5928E-02 -1.8686E+00  6.6929E-01  5.7295E-01  2.7437E-02  1.5790E-01  1.2451E+00  4.0917E-01 -6.8232E-02  5.3511E-04
             1.8465E-01
 GRADIENT:   1.1188E+00  5.5124E-01  2.3609E-01 -1.5486E+01  6.6794E-01  2.5961E-01 -4.4238E-03 -2.6115E-02 -2.7696E-01  8.8653E-02
            -3.2660E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1631.52679973458        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1600
 NPARAMETR:  9.7563E-01  1.2327E-01  1.7676E+00  1.6172E+00  9.2637E-01  1.0590E+00  3.2608E+00  1.3624E+00  8.4398E-01  9.0474E-01
             1.0883E+00
 PARAMETER:  7.5329E-02 -1.9934E+00  6.6964E-01  5.8067E-01  2.3516E-02  1.5728E-01  1.2820E+00  4.0922E-01 -6.9632E-02 -1.0984E-04
             1.8466E-01
 GRADIENT:   3.8081E-01  6.3376E-01 -2.4468E-01 -1.1923E+01  6.1251E-01  5.8493E-02 -9.4686E-02 -1.9775E-02  3.8885E-01  1.3913E-01
             2.6389E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1631.54923558963        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1783
 NPARAMETR:  9.7581E-01  1.1014E-01  1.7699E+00  1.6206E+00  9.2474E-01  1.0592E+00  3.5205E+00  1.3634E+00  8.4191E-01  9.0349E-01
             1.0884E+00
 PARAMETER:  7.5514E-02 -2.1060E+00  6.7095E-01  5.8280E-01  2.1761E-02  1.5756E-01  1.3586E+00  4.0999E-01 -7.2085E-02 -1.4889E-03
             1.8468E-01
 GRADIENT:   1.1963E+00  2.0210E-01 -7.3265E-01 -2.0909E+01  2.4264E+00  2.2729E-01 -3.0892E-03 -1.1159E-02  8.5497E-01 -3.5300E-03
             6.7766E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1631.55379191465        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1962
 NPARAMETR:  9.7617E-01  1.0794E-01  1.7712E+00  1.6227E+00  9.2266E-01  1.0596E+00  3.7279E+00  1.3648E+00  8.3975E-01  9.0223E-01
             1.0883E+00
 PARAMETER:  7.5885E-02 -2.1262E+00  6.7167E-01  5.8411E-01  1.9504E-02  1.5789E-01  1.4158E+00  4.1102E-01 -7.4652E-02 -2.8891E-03
             1.8458E-01
 GRADIENT:   1.9831E+00  4.7148E-01  3.8719E-01 -1.9615E+01 -1.9224E-01  3.6306E-01  2.4354E-01  5.0552E-03  5.1849E-01  2.0795E-01
             2.4830E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1631.56063742108        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2138
 NPARAMETR:  9.7415E-01  9.8966E-02  1.7718E+00  1.6270E+00  9.2114E-01  1.0582E+00  3.6859E+00  1.3655E+00  8.3825E-01  9.0157E-01
             1.0882E+00
 PARAMETER:  7.3806E-02 -2.2130E+00  6.7199E-01  5.8672E-01  1.7860E-02  1.5661E-01  1.4045E+00  4.1149E-01 -7.6444E-02 -3.6151E-03
             1.8454E-01
 GRADIENT:  -1.8339E+00  2.0894E-01 -6.1015E-02 -2.2159E+01  6.6347E-01 -1.0437E-01 -4.2653E-02 -3.1735E-03  2.3745E-01  8.6730E-02
             4.2568E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1631.56665272624        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     2302             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7541E-01  9.5154E-02  1.7715E+00  1.6287E+00  9.2063E-01  1.0595E+00  3.7934E+00  1.3665E+00  8.3728E-01  9.0021E-01
             1.0883E+00
 PARAMETER:  7.5105E-02 -2.2523E+00  6.7181E-01  5.8777E-01  1.7307E-02  1.5776E-01  1.4333E+00  4.1222E-01 -7.7599E-02 -5.1323E-03
             1.8464E-01
 GRADIENT:   3.6652E+02  1.4120E+01  7.5668E+00  9.5168E+02  7.2381E+00  6.9194E+01  4.2189E+00  9.5670E-01  1.5958E+01  4.9735E-01
             1.5529E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1631.57098837777        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2480
 NPARAMETR:  9.7489E-01  9.2171E-02  1.7724E+00  1.6315E+00  9.1962E-01  1.0586E+00  3.8249E+00  1.3662E+00  8.3661E-01  9.0032E-01
             1.0882E+00
 PARAMETER:  7.4571E-02 -2.2841E+00  6.7234E-01  5.8951E-01  1.6209E-02  1.5699E-01  1.4415E+00  4.1201E-01 -7.8392E-02 -5.0087E-03
             1.8453E-01
 GRADIENT:  -1.3280E-01  2.2494E-01 -1.6914E-01 -2.1833E+01  6.6324E-01  6.8537E-02 -3.8738E-02 -1.3150E-02  1.5276E-01  2.1460E-02
             1.5959E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1631.57740854654        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2662             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7577E-01  8.8572E-02  1.7738E+00  1.6347E+00  9.1853E-01  1.0594E+00  4.0396E+00  1.3668E+00  8.3538E-01  8.9984E-01
             1.0882E+00
 PARAMETER:  7.5476E-02 -2.3239E+00  6.7312E-01  5.9148E-01  1.5019E-02  1.5767E-01  1.4961E+00  4.1246E-01 -7.9870E-02 -5.5429E-03
             1.8450E-01
 GRADIENT:   3.6743E+02  1.3639E+01  8.2471E+00  9.6710E+02  5.4614E+00  6.9036E+01  4.4103E+00  8.9824E-01  1.6023E+01  6.5062E-01
             1.4562E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1631.58095286233        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2836
 NPARAMETR:  9.7518E-01  8.3861E-02  1.7743E+00  1.6369E+00  9.1824E-01  1.0588E+00  4.0537E+00  1.3676E+00  8.3497E-01  8.9952E-01
             1.0882E+00
 PARAMETER:  7.4871E-02 -2.3786E+00  6.7341E-01  5.9283E-01  1.4705E-02  1.5710E-01  1.4996E+00  4.1304E-01 -8.0357E-02 -5.8935E-03
             1.8454E-01
 GRADIENT:   7.0097E-01  2.3880E-01 -3.6072E-01 -2.1803E+01  8.7013E-01  1.3405E-01 -1.3337E-03 -1.3105E-02  2.5415E-01 -1.2055E-02
             1.3867E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1631.58477624443        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     2997
 NPARAMETR:  9.7610E-01  8.0262E-02  1.7766E+00  1.6407E+00  9.1721E-01  1.0596E+00  4.1485E+00  1.3684E+00  8.3354E-01  8.9956E-01
             1.0882E+00
 PARAMETER:  7.5805E-02 -2.4225E+00  6.7470E-01  5.9515E-01  1.3579E-02  1.5786E-01  1.5228E+00  4.1362E-01 -8.2075E-02 -5.8469E-03
             1.8448E-01
 GRADIENT:   2.6158E+00  3.4431E-01  6.3887E-02 -1.9243E+01 -3.5508E-01  4.4022E-01 -1.7101E-03 -6.0589E-02 -1.6178E-02  8.4785E-02
            -5.7590E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1631.58967149946        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     3178
 NPARAMETR:  9.7482E-01  6.8648E-02  1.7800E+00  1.6498E+00  9.1651E-01  1.0587E+00  4.3628E+00  1.3703E+00  8.3278E-01  8.9992E-01
             1.0883E+00
 PARAMETER:  7.4496E-02 -2.5788E+00  6.7660E-01  6.0065E-01  1.2818E-02  1.5707E-01  1.5731E+00  4.1501E-01 -8.2981E-02 -5.4443E-03
             1.8459E-01
 GRADIENT:   3.5745E-01  3.2748E-01 -8.8440E-01 -1.6718E+01  1.1455E+00  1.5677E-01 -5.4443E-02 -7.7450E-02  5.1660E-01 -1.7472E-02
            -7.6346E-03

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1631.60362892880        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3362
 NPARAMETR:  9.7520E-01  5.6916E-02  1.7869E+00  1.6522E+00  9.1491E-01  1.0591E+00  4.9536E+00  1.3738E+00  8.2965E-01  9.0010E-01
             1.0883E+00
 PARAMETER:  7.4890E-02 -2.7662E+00  6.8048E-01  6.0210E-01  1.1075E-02  1.5740E-01  1.7001E+00  4.1761E-01 -8.6753E-02 -5.2458E-03
             1.8463E-01
 GRADIENT:   1.5448E+00  1.3727E-01 -2.7530E-01 -2.6397E+01  7.7584E-01  3.4720E-01  2.3700E-02 -8.3637E-02  3.8466E-01  1.1447E-01
             6.1499E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1631.60914852300        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3540
 NPARAMETR:  9.7451E-01  4.6302E-02  1.7910E+00  1.6643E+00  9.1353E-01  1.0585E+00  5.2373E+00  1.3766E+00  8.2772E-01  8.9976E-01
             1.0883E+00
 PARAMETER:  7.4181E-02 -2.9726E+00  6.8275E-01  6.0938E-01  9.5576E-03  1.5685E-01  1.7558E+00  4.1959E-01 -8.9079E-02 -5.6219E-03
             1.8463E-01
 GRADIENT:   3.6672E-01  2.6279E-01 -6.9814E-01 -1.7184E+01  4.4895E-01  1.4870E-01 -4.9296E-02 -1.4314E-01  3.0607E-01  6.4039E-02
            -1.7242E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1631.62171874466        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3723
 NPARAMETR:  9.7465E-01  3.9520E-02  1.7945E+00  1.6661E+00  9.1306E-01  1.0587E+00  5.8562E+00  1.3796E+00  8.2652E-01  8.9934E-01
             1.0883E+00
 PARAMETER:  7.4319E-02 -3.1310E+00  6.8470E-01  6.1049E-01  9.0495E-03  1.5700E-01  1.8675E+00  4.2181E-01 -9.0527E-02 -6.0953E-03
             1.8465E-01
 GRADIENT:   8.8494E-01  1.8301E-01 -8.0984E-01 -2.2043E+01  1.1035E+00  2.3175E-01 -1.1523E-03 -1.0277E-01  5.1827E-01  3.4853E-02
             3.3538E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1631.62236610240        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     3902
 NPARAMETR:  9.7425E-01  3.0100E-02  1.7999E+00  1.6753E+00  9.1224E-01  1.0583E+00  6.6641E+00  1.3841E+00  8.2457E-01  8.9940E-01
             1.0884E+00
 PARAMETER:  7.3911E-02 -3.4032E+00  6.8773E-01  6.1602E-01  8.1493E-03  1.5670E-01  1.9967E+00  4.2502E-01 -9.2892E-02 -6.0221E-03
             1.8466E-01
 GRADIENT:   3.0135E-01  2.1954E-01 -1.0622E+00 -1.6703E+01  8.3137E-01  1.3039E-01 -4.8151E-03 -9.4595E-02  4.7972E-01  6.4090E-02
             8.3220E-03

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1631.63480328839        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     4084
 NPARAMETR:  9.7442E-01  2.0145E-02  1.8065E+00  1.6794E+00  9.1116E-01  1.0585E+00  7.4905E+00  1.3880E+00  8.2235E-01  8.9874E-01
             1.0884E+00
 PARAMETER:  7.4090E-02 -3.8048E+00  6.9139E-01  6.1846E-01  6.9668E-03  1.5682E-01  2.1136E+00  4.2785E-01 -9.5595E-02 -6.7594E-03
             1.8467E-01
 GRADIENT:   9.4991E-01  1.1338E-01 -7.3552E-01 -2.1306E+01  6.3830E-01  2.2319E-01 -3.3749E-02 -1.5901E-01  2.8546E-01  3.7054E-03
            -3.8893E-03

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1631.64242472433        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     4254             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7462E-01  1.7567E-02  1.8105E+00  1.6790E+00  9.1057E-01  1.0588E+00  8.7428E+00  1.3915E+00  8.2068E-01  8.9828E-01
             1.0883E+00
 PARAMETER:  7.4290E-02 -3.9417E+00  6.9361E-01  6.1820E-01  6.3186E-03  1.5710E-01  2.2682E+00  4.3038E-01 -9.7620E-02 -7.2748E-03
             1.8460E-01
 GRADIENT:   3.6650E+02  3.0037E+00  8.2720E+00  1.0665E+03  5.2361E+00  6.8594E+01  1.1905E+00  8.5181E-01  1.6334E+01  6.5010E-01
             1.4735E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1631.64514455448        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     4438             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7451E-01  1.5043E-02  1.8119E+00  1.6819E+00  9.1058E-01  1.0586E+00  9.6301E+00  1.3931E+00  8.2053E-01  8.9814E-01
             1.0883E+00
 PARAMETER:  7.4181E-02 -4.0969E+00  6.9435E-01  6.1990E-01  6.3264E-03  1.5697E-01  2.3649E+00  4.3151E-01 -9.7801E-02 -7.4340E-03
             1.8464E-01
 GRADIENT:   3.6632E+02  2.7186E+00  8.0131E+00  1.0737E+03  5.5070E+00  6.8499E+01  1.1145E+00  8.7412E-01  1.6576E+01  6.2032E-01
             1.4832E+00

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1631.64569078224        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4617
 NPARAMETR:  9.7436E-01  1.0732E-02  1.8140E+00  1.6859E+00  9.1035E-01  1.0584E+00  1.0658E+01  1.3942E+00  8.2007E-01  8.9826E-01
             1.0883E+00
 PARAMETER:  7.4029E-02 -4.4345E+00  6.9551E-01  6.2232E-01  6.0793E-03  1.5679E-01  2.4663E+00  4.3231E-01 -9.8364E-02 -7.2949E-03
             1.8466E-01
 GRADIENT:   1.1266E+00  8.6333E-02 -5.2643E-01 -2.1069E+01  5.0165E-02  2.2982E-01 -8.3447E-03 -1.3715E-01  2.4332E-01  4.2755E-02
            -7.0875E-03

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1631.65009694005        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     4806
 NPARAMETR:  9.7451E-01  1.0000E-02  1.8171E+00  1.6845E+00  9.1031E-01  1.0587E+00  1.1891E+01  1.3981E+00  8.1902E-01  8.9742E-01
             1.0883E+00
 PARAMETER:  7.4178E-02 -4.5680E+00  6.9725E-01  6.2146E-01  6.0259E-03  1.5705E-01  2.5758E+00  4.3508E-01 -9.9649E-02 -8.2302E-03
             1.8460E-01
 GRADIENT:   1.4089E+00  0.0000E+00 -7.4765E-02 -2.4769E+01 -5.1204E-01  3.4745E-01  1.5904E-02 -8.7610E-02 -6.8138E-02  5.0655E-02
            -5.6240E-03

0ITERATION NO.:  138    OBJECTIVE VALUE:  -1631.65019169366        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     4907
 NPARAMETR:  9.7451E-01  1.0000E-02  1.8222E+00  1.6847E+00  9.1050E-01  1.0587E+00  1.1889E+01  1.4007E+00  8.1908E-01  8.9742E-01
             1.0881E+00
 PARAMETER:  7.4150E-02 -4.5680E+00  6.9747E-01  6.2165E-01  6.8790E-03  1.5698E-01  2.5845E+00  4.3798E-01 -9.8941E-02 -9.2302E-03
             1.8462E-01
 GRADIENT:  -3.4083E-02  0.0000E+00 -6.7159E-01  1.5661E-01  5.8479E-01 -2.3112E-02  3.6286E-03  5.1233E-02  1.5635E-01 -3.7688E-02
             4.9155E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4907
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7482E-05 -2.6135E-05 -2.6431E-02 -5.6360E-03 -3.8197E-02
 SE:             2.9810E-02  1.8021E-03  1.7891E-02  2.9257E-02  1.9877E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9926E-01  9.8843E-01  1.3960E-01  8.4724E-01  5.4648E-02

 ETASHRINKSD(%)  1.3323E-01  9.3963E+01  4.0062E+01  1.9857E+00  3.3410E+01
 ETASHRINKVR(%)  2.6629E-01  9.9635E+01  6.4074E+01  3.9321E+00  5.5658E+01
 EBVSHRINKSD(%)  4.5127E-01  9.4122E+01  4.3105E+01  2.2926E+00  3.1132E+01
 EBVSHRINKVR(%)  9.0050E-01  9.9654E+01  6.7629E+01  4.5325E+00  5.2573E+01
 RELATIVEINF(%)  9.5107E+01  7.7958E-03  6.6491E+00  2.6364E+00  5.6140E+00
 EPSSHRINKSD(%)  4.3962E+01
 EPSSHRINKVR(%)  6.8597E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1631.6501916936572     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -896.49936512991906     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    69.59
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     6.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1631.650       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  1.00E-02  1.82E+00  1.68E+00  9.11E-01  1.06E+00  1.20E+01  1.40E+00  8.20E-01  8.97E-01  1.09E+00
 


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
 
         3.10E-02  0.00E+00  3.42E-01  4.38E-02  8.51E-02  8.74E-02  6.21E-01  3.78E-01  5.65E-02  1.32E-01  6.46E-02
 


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
+        9.62E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.60E-03  0.00E+00  1.17E-01
 
 TH 4
+       -1.61E-04  0.00E+00  3.99E-03  1.92E-03
 
 TH 5
+       -6.58E-04  0.00E+00  2.59E-02  1.04E-03  7.23E-03
 
 TH 6
+       -3.00E-04  0.00E+00  1.70E-03  3.59E-04 -2.40E-05  7.63E-03
 
 TH 7
+       -3.00E-03  0.00E+00  9.31E-02  3.27E-04  2.36E-02 -3.78E-03  3.85E-01
 
 TH 8
+       -4.18E-03  0.00E+00  9.61E-02  1.38E-03  1.97E-02  3.88E-03  8.81E-02  1.43E-01
 
 TH 9
+        2.77E-04  0.00E+00 -2.83E-03 -7.14E-04 -5.36E-04 -3.98E-04  1.59E-02 -3.96E-03  3.19E-03
 
 TH10
+       -2.98E-04  0.00E+00  1.00E-02  8.57E-04  4.25E-03  1.64E-03  2.03E-03 -6.02E-03 -7.53E-04  1.74E-02
 
 TH11
+        5.21E-04  0.00E+00 -2.65E-03  4.37E-04 -6.40E-04  5.87E-04 -1.76E-03 -3.92E-03  4.43E-05 -6.14E-04  4.18E-03
 
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
+        3.10E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.45E-01  0.00E+00  3.42E-01
 
 TH 4
+       -1.18E-01  0.00E+00  2.66E-01  4.38E-02
 
 TH 5
+       -2.49E-01  0.00E+00  8.92E-01  2.80E-01  8.51E-02
 
 TH 6
+       -1.11E-01  0.00E+00  5.69E-02  9.37E-02 -3.23E-03  8.74E-02
 
 TH 7
+       -1.56E-01  0.00E+00  4.39E-01  1.20E-02  4.46E-01 -6.97E-02  6.21E-01
 
 TH 8
+       -3.56E-01  0.00E+00  7.43E-01  8.33E-02  6.13E-01  1.17E-01  3.75E-01  3.78E-01
 
 TH 9
+        1.58E-01  0.00E+00 -1.47E-01 -2.88E-01 -1.11E-01 -8.07E-02  4.55E-01 -1.85E-01  5.65E-02
 
 TH10
+       -7.27E-02  0.00E+00  2.22E-01  1.48E-01  3.78E-01  1.42E-01  2.47E-02 -1.20E-01 -1.01E-01  1.32E-01
 
 TH11
+        2.60E-01  0.00E+00 -1.20E-01  1.54E-01 -1.16E-01  1.04E-01 -4.37E-02 -1.60E-01  1.21E-02 -7.19E-02  6.46E-02
 
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
+        1.11E+03
 
 TH 2
+        5.81E-14  2.94E-29
 
 TH 3
+       -1.04E+01 -1.18E-14  4.36E+01
 
 TH 4
+       -5.21E+00  2.30E-14  3.20E-01  6.02E+02
 
 TH 5
+        6.92E+01  6.39E-14 -1.85E+02 -7.76E+01  8.45E+02
 
 TH 6
+        6.21E+01  2.24E-14 -1.37E+01 -9.35E+00  6.88E+01  1.41E+02
 
 TH 7
+       -1.06E-02 -8.35E-18  2.19E-03  4.82E-03 -7.28E-03  6.29E-04  4.54E-06
 
 TH 8
+       -9.20E+00 -5.25E-15 -1.04E+00 -9.71E+00 -3.98E+00 -2.00E+00  1.14E-04  2.31E+00
 
 TH 9
+       -8.91E+01 -6.28E-14  1.69E+01  1.08E+02 -5.28E+01  7.91E+00  3.90E-02 -2.37E+00  3.45E+02
 
 TH10
+       -2.07E+01 -1.42E-14  1.72E+01 -1.16E+01 -9.13E+01 -7.38E+00  7.13E-04  4.06E+00 -5.27E-01  1.57E+01
 
 TH11
+       -1.10E+02 -6.22E-14 -1.47E+01 -5.25E+01 -2.95E+01 -4.05E+01  5.46E-03  2.39E+01  1.68E+01  3.98E+01  2.58E+02
 
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
+        1.03E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.28E+01  0.00E+00  5.69E+01
 
 TH 4
+       -3.06E+01  0.00E+00 -7.21E+00  5.56E+02
 
 TH 5
+       -2.63E+01  0.00E+00 -1.80E+02 -5.37E+01  8.22E+02
 
 TH 6
+       -5.13E+01  0.00E+00  6.76E+00  2.51E+01 -9.26E+01  2.31E+02
 
 TH 7
+       -7.12E-03  0.00E+00  1.70E-03 -2.66E-02  9.59E-03 -2.37E-03  1.19E-05
 
 TH 8
+       -3.50E+01  0.00E+00 -1.25E+01 -1.76E+01 -8.67E+00  1.35E+01  7.01E-04  1.77E+01
 
 TH 9
+        6.99E+01  0.00E+00  7.27E+00 -9.41E+01  4.34E+01 -1.64E+01  4.20E-02 -8.76E+00  2.32E+02
 
 TH10
+       -2.39E+01  0.00E+00  9.41E-01 -1.50E+00 -7.97E+01  3.37E+01  5.19E-05  1.44E+01 -1.63E+01  5.53E+01
 
 TH11
+        6.46E+01  0.00E+00 -1.01E+01  3.60E+01 -1.75E+01  3.15E+01 -4.57E-04  1.02E+01  2.20E+00  1.33E+01  1.33E+02
 
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
 #CPUT: Total CPU Time in Seconds,       76.172
Stop Time:
Wed Sep 29 18:04:54 CDT 2021
