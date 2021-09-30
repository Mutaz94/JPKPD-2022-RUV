Wed Sep 29 10:13:04 CDT 2021
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
$DATA ../../../../data/int/D/dat97.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   52363.6848561350        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.9879E+02  8.5145E+02  5.0453E+02  7.0634E+02  4.7045E+02 -4.7032E+03 -2.0398E+03 -1.6194E+03 -2.8342E+03 -5.0006E+02
            -1.0067E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -495.902710317203        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  8.7377E-01  2.0855E+00  8.0812E-01  1.7116E+00  1.1512E+00  6.2516E+00  6.2716E+00  1.0033E+00  1.9868E+00  1.2671E+00
             1.2564E+01
 PARAMETER: -3.4937E-02  8.3500E-01 -1.1305E-01  6.3742E-01  2.4077E-01  1.9328E+00  1.9360E+00  1.0325E-01  7.8654E-01  3.3670E-01
             2.6308E+00
 GRADIENT:  -1.2402E+01  2.7745E+01 -7.6133E+01  1.6000E+02  4.3273E+01  2.1830E+02  1.3235E+02  2.7347E+00  9.6728E+00  1.8039E+01
             8.6837E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -596.417871055691        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  7.2008E-01  2.7428E+00  1.0038E+01  2.8134E+00  1.8906E+00  4.0457E+00  8.7656E+00  6.8941E-01  3.9483E+00  1.0731E+00
             1.2604E+01
 PARAMETER: -2.2839E-01  1.1090E+00  2.4064E+00  1.1344E+00  7.3687E-01  1.4976E+00  2.2708E+00 -2.7191E-01  1.4733E+00  1.7056E-01
             2.6340E+00
 GRADIENT:  -3.6949E+01  4.6049E+01 -1.7665E+01  9.7376E+01 -3.2676E+01  1.5130E+02  1.3847E+02  2.9287E-01  6.3296E+01  1.8586E+01
             1.3415E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -705.392542956381        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0141E+00  1.4406E+00  7.8367E+00  9.7713E-01  2.2391E+00  1.8681E+00  4.7201E+00  1.0000E-02  1.3348E+00  7.7973E-01
             1.3483E+01
 PARAMETER:  1.1400E-01  4.6509E-01  2.1588E+00  7.6860E-02  9.0608E-01  7.2492E-01  1.6518E+00 -5.6339E+00  3.8876E-01 -1.4881E-01
             2.7014E+00
 GRADIENT:   1.0801E+01 -2.1978E+01 -2.6757E+00 -4.1205E+01  1.7186E+00  1.5380E+01  9.0991E+00  0.0000E+00  1.7267E+01  1.0889E+01
             2.2995E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -729.481556792543        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.3238E-01  1.2767E+00  1.8246E+01  1.0316E+00  2.4464E+00  1.8200E+00  4.7854E+00  5.2141E-02  7.3127E-01  3.0627E-01
             1.1713E+01
 PARAMETER:  2.9986E-02  3.4431E-01  3.0040E+00  1.3110E-01  9.9461E-01  6.9885E-01  1.6656E+00 -2.8538E+00 -2.1297E-01 -1.0833E+00
             2.5607E+00
 GRADIENT:   8.9748E-01 -1.3229E+01 -4.8971E-01 -8.2961E+00  5.4748E+00  2.1709E+00  7.6737E+00  1.6760E-04  1.4398E+00  1.6682E+00
             1.1614E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -759.892731578593        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      481
 NPARAMETR:  9.5033E-01  1.0333E+00  3.1859E+01  1.3459E+00  2.4823E+00  1.8441E+00  7.3641E+00  1.0000E-02  8.4760E-01  2.8918E-01
             1.2156E+01
 PARAMETER:  4.9051E-02  1.3271E-01  3.5613E+00  3.9706E-01  1.0092E+00  7.1198E-01  2.0966E+00 -4.9733E+00 -6.5349E-02 -1.1407E+00
             2.5978E+00
 GRADIENT:  -1.0375E-01  6.9882E+00 -6.0206E-01  3.5992E+00 -3.5145E+00 -4.5683E-02  2.6155E+00  0.0000E+00 -2.1104E+00  1.2817E+00
            -2.2733E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -766.013540100000        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  9.6755E-01  6.0890E-01  1.8315E+02  1.5147E+00  2.5404E+00  1.8231E+00  9.3748E+00  1.0000E-02  1.0679E+00  5.8273E-02
             1.2055E+01
 PARAMETER:  6.7007E-02 -3.9609E-01  5.3103E+00  5.1524E-01  1.0323E+00  7.0052E-01  2.3380E+00 -7.1065E+00  1.6569E-01 -2.7426E+00
             2.5895E+00
 GRADIENT:   1.1526E+01  4.7116E+00 -7.7352E-02 -9.7502E+00 -5.5943E+00 -3.5400E+00  4.3353E+01  0.0000E+00  1.5169E+00  5.3647E-02
            -1.0307E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -767.002368486722        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  9.5209E-01  4.9561E-01  1.0143E+03  1.5671E+00  2.6131E+00  1.8461E+00  9.1296E+00  1.0000E-02  1.0902E+00  1.7266E-02
             1.2130E+01
 PARAMETER:  5.0908E-02 -6.0198E-01  7.0220E+00  5.4922E-01  1.0605E+00  7.1309E-01  2.3115E+00 -7.1065E+00  1.8637E-01 -3.9590E+00
             2.5957E+00
 GRADIENT:   1.3234E+00  1.0813E+00 -2.1789E-02  3.5749E+00 -5.4042E-01  8.3056E-01  3.0762E+01  0.0000E+00 -8.1487E-01  4.7238E-03
            -3.2357E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -767.445513384623        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1000             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4771E-01  4.3182E-01  1.6448E+04  1.5876E+00  2.6061E+00  1.8410E+00  9.6431E+00  1.0000E-02  1.0809E+00  1.0000E-02
             1.2094E+01
 PARAMETER:  4.6289E-02 -7.3975E-01  9.8079E+00  5.6220E-01  1.0579E+00  7.1033E-01  2.3662E+00 -7.1065E+00  1.7778E-01 -5.5421E+00
             2.5927E+00
 GRADIENT:   3.1762E+00  5.5169E+00 -1.1351E-03  2.1327E+01  5.3300E-01  1.0099E+01  3.3444E+02  0.0000E+00 -2.3256E+00  0.0000E+00
             4.1497E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -767.489839598511        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1176
 NPARAMETR:  9.5138E-01  4.2715E-01  3.5108E+04  1.5962E+00  2.6226E+00  1.8429E+00  9.6360E+00  1.0000E-02  1.1206E+00  1.0000E-02
             1.2146E+01
 PARAMETER:  5.0162E-02 -7.5062E-01  1.0566E+01  5.6764E-01  1.0642E+00  7.1135E-01  2.3655E+00 -7.1065E+00  2.1382E-01 -5.5421E+00
             2.5970E+00
 GRADIENT:  -3.8970E+00  4.7959E-01 -7.1875E-04  8.1300E-01  7.6041E-01  2.9684E+00  4.0647E+01  0.0000E+00 -2.0769E-01  0.0000E+00
             6.1005E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -767.643335091558        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1367
 NPARAMETR:  9.4533E-01  3.8221E-01  5.6847E+04  1.6085E+00  2.6126E+00  1.8409E+00  9.9590E+00  1.0000E-02  1.1191E+00  1.0000E-02
             1.2108E+01
 PARAMETER:  4.3775E-02 -8.6179E-01  1.1048E+01  5.7533E-01  1.0604E+00  7.1023E-01  2.3985E+00 -7.1065E+00  2.1254E-01 -5.5421E+00
             2.5939E+00
 GRADIENT:  -1.2935E+01 -3.4177E-01 -4.1051E-04  2.7687E+00 -1.7786E-01  4.6307E+00  4.6700E+01  0.0000E+00 -1.7216E+00  0.0000E+00
             4.2566E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -767.672824457030        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1560
 NPARAMETR:  9.4760E-01  3.7414E-01  7.1671E+04  1.6127E+00  2.6140E+00  1.8411E+00  1.0057E+01  1.0000E-02  1.1305E+00  1.0000E-02
             1.2106E+01
 PARAMETER:  4.6174E-02 -8.8313E-01  1.1280E+01  5.7792E-01  1.0609E+00  7.1039E-01  2.4082E+00 -7.1065E+00  2.2269E-01 -5.5421E+00
             2.5937E+00
 GRADIENT:  -3.4550E+01 -2.8451E+00 -3.6424E-04 -1.4400E+00  8.6734E-01  1.8199E+01  4.9301E+01  0.0000E+00 -4.7567E-02  0.0000E+00
             2.9566E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -767.697439260752        NO. OF FUNC. EVALS.: 202
 CUMULATIVE NO. OF FUNC. EVALS.:     1762             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4796E-01  3.5619E-01  1.2452E+05  1.6276E+00  2.6141E+00  1.8416E+00  1.0277E+01  1.0000E-02  1.1520E+00  1.0000E-02
             1.2103E+01
 PARAMETER:  4.6562E-02 -9.3228E-01  1.1832E+01  5.8709E-01  1.0609E+00  7.1065E-01  2.4300E+00 -7.1065E+00  2.4154E-01 -5.5421E+00
             2.5934E+00
 GRADIENT:  -2.7472E+02 -2.1378E+01 -3.1637E-04  2.4346E-01  1.1406E+01  1.3365E+02  3.8138E+02  0.0000E+00  1.3984E+01  0.0000E+00
             3.0251E+02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -767.701749003854        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1956
 NPARAMETR:  9.4836E-01  3.5103E-01  1.4879E+05  1.6296E+00  2.6138E+00  1.8417E+00  1.0278E+01  1.0000E-02  1.1538E+00  1.0000E-02
             1.2103E+01
 PARAMETER:  4.6976E-02 -9.4688E-01  1.2010E+01  5.8836E-01  1.0608E+00  7.1070E-01  2.4300E+00 -7.1065E+00  2.4305E-01 -5.5421E+00
             2.5934E+00
 GRADIENT:  -5.1931E+01 -2.2387E+00 -1.7710E-04  8.9228E-01  1.0666E+00  2.2580E+01  5.3254E+01  0.0000E+00  9.2286E-01  0.0000E+00
             2.7664E+01

0ITERATION NO.:   68    OBJECTIVE VALUE:  -767.702386277096        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:     2067
 NPARAMETR:  9.4830E-01  3.4992E-01  1.2967E+05  1.6302E+00  2.6157E+00  1.8424E+00  1.0342E+01  1.0000E-02  1.1533E+00  1.0000E-02
             1.2109E+01
 PARAMETER:  4.6928E-02 -9.5966E-01  1.1921E+01  5.8930E-01  1.0608E+00  7.1084E-01  2.4283E+00 -7.1065E+00  2.4137E-01 -5.5421E+00
             2.5934E+00
 GRADIENT:   4.2514E-03 -5.7126E-02  9.3345E-05  1.5628E-01 -8.6367E-02 -2.1570E-02 -5.3661E-01  0.0000E+00 -4.0732E-02  0.0000E+00
            -3.3842E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2067
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1318E-02  6.4138E-02  6.8669E-10 -9.2342E-02  1.2202E-05
 SE:             2.7457E-02  2.1316E-02  1.0116E-09  1.5594E-02  8.1673E-05
 N:                     100         100         100         100         100

 P VAL.:         4.3749E-01  2.6223E-03  4.9726E-01  3.1979E-09  8.8124E-01

 ETASHRINKSD(%)  8.0167E+00  2.8588E+01  1.0000E+02  4.7758E+01  9.9726E+01
 ETASHRINKVR(%)  1.5391E+01  4.9003E+01  1.0000E+02  7.2708E+01  9.9999E+01
 EBVSHRINKSD(%)  1.0327E+01  2.4346E+01  1.0000E+02  4.4148E+01  9.9606E+01
 EBVSHRINKVR(%)  1.9588E+01  4.2764E+01  1.0000E+02  6.8805E+01  9.9998E+01
 RELATIVEINF(%)  8.0207E+01  3.7531E+01  1.0000E-10  2.0165E+01  1.0000E-10
 EPSSHRINKSD(%)  3.9084E+00
 EPSSHRINKVR(%)  7.6641E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -767.70238627709648     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       886.38697349131428     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    74.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -767.702       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.48E-01  3.47E-01  1.36E+05  1.63E+00  2.61E+00  1.84E+00  1.03E+01  1.00E-02  1.15E+00  1.00E-02  1.21E+01
 


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
+        2.71E+02
 
 TH 2
+        6.36E+00  6.07E+01
 
 TH 3
+       -1.25E-08  4.56E-07  3.51E-15
 
 TH 4
+       -7.82E+00  2.20E+01  1.45E-07  1.47E+02
 
 TH 5
+        1.93E+00 -3.42E+00 -3.64E-09 -1.40E+01  2.43E+01
 
 TH 6
+       -3.16E+00 -1.53E+00 -1.20E-07 -1.51E+00  1.75E-03  3.99E+01
 
 TH 7
+        4.66E-01  4.19E+00  1.25E-08 -3.91E+00  2.30E-01  1.25E-01  9.49E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.34E+00 -3.91E+00 -8.66E-09 -3.85E+01  4.13E+00 -1.77E+00  6.10E-01  0.00E+00  3.44E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.14E+00 -2.26E+00 -8.36E-09 -9.71E+00  6.19E-01  1.86E+00  1.31E-01  0.00E+00  4.00E+00  0.00E+00  6.08E+00
 
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
 #CPUT: Total CPU Time in Seconds,       92.564
Stop Time:
Wed Sep 29 10:14:38 CDT 2021
