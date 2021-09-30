Wed Sep 29 19:04:40 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat50.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1688.48832179067        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2237E+02 -1.1697E+01 -2.3462E+01  7.0926E+00  1.7821E+00  3.7766E+01  9.4149E+00  1.8188E+01  1.0798E+01  2.1382E+01
             6.8928E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1691.91563673884        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.9652E-01  1.0395E+00  1.2009E+00  1.0210E+00  1.0983E+00  1.0278E+00  9.4073E-01  8.5939E-01  9.5476E-01  8.6475E-01
             9.7809E-01
 PARAMETER:  9.6515E-02  1.3870E-01  2.8309E-01  1.2082E-01  1.9379E-01  1.2740E-01  3.8901E-02 -5.1529E-02  5.3709E-02 -4.5318E-02
             7.7848E-02
 GRADIENT:   6.8117E-02  9.6591E+00  1.5764E+01 -1.2673E+01  1.4235E+01  6.8235E+00  1.1275E+00  2.2782E+00 -1.2232E+01 -2.6647E+01
            -1.9860E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1696.85747600918        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8166E-01  8.4053E-01  1.3776E+00  1.1805E+00  1.1290E+00  9.8961E-01  6.6922E-01  4.0561E-01  9.6850E-01  1.1594E+00
             1.0234E+00
 PARAMETER:  8.1486E-02 -7.3718E-02  4.2035E-01  2.6592E-01  2.2134E-01  8.9556E-02 -3.0164E-01 -8.0235E-01  6.7994E-02  2.4786E-01
             1.2309E-01
 GRADIENT:  -2.9926E+01  2.4994E+01 -5.7804E+00  5.0149E+01  1.6928E+01 -7.1829E+00 -6.1035E-01  2.2092E-01 -5.3160E+00  7.2920E-01
             6.0134E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.98616961749        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.9574E-01  6.9728E-01  1.2153E+00  1.2405E+00  9.7015E-01  1.0038E+00  8.5415E-01  2.6689E-01  9.0718E-01  1.0466E+00
             1.0000E+00
 PARAMETER:  9.5731E-02 -2.6057E-01  2.9497E-01  3.1548E-01  6.9696E-02  1.0380E-01 -5.7643E-02 -1.2209E+00  2.5867E-03  1.4552E-01
             1.0002E-01
 GRADIENT:   3.8340E+00  1.0080E+01  7.7447E+00  7.1698E+00 -1.5917E+01 -7.9910E-01  1.4921E-01  1.5299E-01 -7.3209E-01  1.3633E+00
            -5.1010E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1699.99916688757        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.8816E-01  3.5760E-01  1.3027E+00  1.4540E+00  9.0763E-01  1.0018E+00  8.2768E-01  5.6602E-02  8.0667E-01  1.0796E+00
             1.0019E+00
 PARAMETER:  8.8085E-02 -9.2835E-01  3.6446E-01  4.7433E-01  3.0851E-03  1.0184E-01 -8.9130E-02 -2.7717E+00 -1.1484E-01  1.7659E-01
             1.0193E-01
 GRADIENT:  -2.5136E+00  4.5114E+00  2.1955E+00  1.4259E+01 -5.3992E+00  5.3183E-01  1.1851E-01  1.5239E-03  6.5327E-01  2.1554E+00
             5.7361E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1700.06941613003        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.8800E-01  2.7345E-01  1.3176E+00  1.5068E+00  8.9032E-01  9.9905E-01  7.6306E-01  2.8168E-02  7.7606E-01  1.0695E+00
             9.9982E-01
 PARAMETER:  8.7930E-02 -1.1966E+00  3.7581E-01  5.1001E-01 -1.6176E-02  9.9049E-02 -1.7042E-01 -3.4696E+00 -1.5352E-01  1.6715E-01
             9.9824E-02
 GRADIENT:   1.5142E-01  4.0374E+00  1.9857E+00  1.7568E+01 -3.3415E+00 -9.8283E-03 -1.5970E-04 -4.7183E-04 -2.5430E+00 -3.7926E-01
            -1.3168E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1700.13569365949        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  9.8712E-01  2.0308E-01  1.3259E+00  1.5473E+00  8.7493E-01  9.9731E-01  6.0567E-01  1.2305E-02  7.5711E-01  1.0667E+00
             9.9975E-01
 PARAMETER:  8.7032E-02 -1.4941E+00  3.8208E-01  5.3654E-01 -3.3610E-02  9.7311E-02 -4.0142E-01 -4.2977E+00 -1.7824E-01  1.6461E-01
             9.9748E-02
 GRADIENT:   8.9416E-01  2.5286E+00  1.3814E+00  1.2164E+01 -1.8464E+00 -1.9591E-01  9.6529E-04 -1.1671E-04 -2.4092E+00 -9.2342E-01
            -1.3895E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1700.14824504896        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1232
 NPARAMETR:  9.8623E-01  1.5849E-01  1.3305E+00  1.5734E+00  8.6490E-01  9.9651E-01  4.6992E-01  1.0000E-02  7.4602E-01  1.0670E+00
             1.0003E+00
 PARAMETER:  8.6135E-02 -1.7420E+00  3.8556E-01  5.5322E-01 -4.5141E-02  9.6506E-02 -6.5520E-01 -5.0076E+00 -1.9301E-01  1.6483E-01
             1.0026E-01
 GRADIENT:   7.1466E-01  1.7807E+00  1.2068E+00  9.2379E+00 -1.6936E+00 -1.8914E-01  3.8289E-03  0.0000E+00 -1.6847E+00 -6.9535E-01
            -1.0293E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1700.22255776431        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1417
 NPARAMETR:  9.8610E-01  1.4463E-01  1.3307E+00  1.5737E+00  8.6287E-01  9.9673E-01  2.9595E-01  1.0000E-02  7.4626E-01  1.0705E+00
             1.0015E+00
 PARAMETER:  8.6006E-02 -1.8336E+00  3.8567E-01  5.5345E-01 -4.7491E-02  9.6726E-02 -1.1175E+00 -5.1156E+00 -1.9267E-01  1.6817E-01
             1.0146E-01
 GRADIENT:   1.1583E+00  3.1121E-01  2.0606E-01 -1.0233E+01  1.1009E+00  4.0086E-02  4.1308E-03  0.0000E+00  4.0972E-01  3.2424E-02
             1.4267E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1700.22978187746        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1608
 NPARAMETR:  9.8623E-01  1.4495E-01  1.3270E+00  1.5723E+00  8.6131E-01  9.9683E-01  5.3082E-02  1.0000E-02  7.4523E-01  1.0695E+00
             1.0010E+00
 PARAMETER:  8.6134E-02 -1.8314E+00  3.8289E-01  5.5252E-01 -4.9296E-02  9.6825E-02 -2.8359E+00 -5.1156E+00 -1.9406E-01  1.6723E-01
             1.0103E-01
 GRADIENT:   1.4363E+00  2.2878E-01  2.1045E-01 -1.2436E+01  9.9493E-01  8.3524E-02  2.4997E-04  0.0000E+00 -2.7204E-01  8.9546E-02
             2.4461E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1700.23308017205        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1800
 NPARAMETR:  9.8624E-01  1.4562E-01  1.3236E+00  1.5713E+00  8.6016E-01  9.9687E-01  2.3423E-02  1.0000E-02  7.4571E-01  1.0681E+00
             1.0009E+00
 PARAMETER:  8.6145E-02 -1.8268E+00  3.8032E-01  5.5190E-01 -5.0633E-02  9.6861E-02 -3.6541E+00 -5.1156E+00 -1.9342E-01  1.6584E-01
             1.0092E-01
 GRADIENT:   1.4346E+00  1.6961E-01 -4.1121E-02 -1.3142E+01  1.2700E+00  9.3113E-02  6.9915E-05  0.0000E+00 -1.6252E-01  3.3134E-02
             4.3528E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1700.23570355690        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1996             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8625E-01  1.4642E-01  1.3215E+00  1.5704E+00  8.5872E-01  9.9689E-01  1.7343E-02  1.0000E-02  7.4626E-01  1.0671E+00
             1.0008E+00
 PARAMETER:  8.6152E-02 -1.8213E+00  3.7880E-01  5.5132E-01 -5.2307E-02  9.6886E-02 -3.9546E+00 -5.1156E+00 -1.9268E-01  1.6498E-01
             1.0079E-01
 GRADIENT:   4.0655E+02  1.4717E+01  8.0164E+00  8.8447E+02  6.8644E+00  3.9399E+01  9.0630E-04  0.0000E+00  1.7438E+01  1.4127E+00
             7.9563E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1700.23673039461        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2183
 NPARAMETR:  9.8626E-01  1.4733E-01  1.3200E+00  1.5699E+00  8.5790E-01  9.9690E-01  1.1032E-02  1.0000E-02  7.4646E-01  1.0665E+00
             1.0007E+00
 PARAMETER:  8.6164E-02 -1.8151E+00  3.7765E-01  5.5100E-01 -5.3265E-02  9.6897E-02 -4.4069E+00 -5.1156E+00 -1.9241E-01  1.6436E-01
             1.0067E-01
 GRADIENT:   1.3967E+00  2.5169E-01  1.1468E+00 -1.3235E+01 -8.0715E-01  9.5962E-02  2.1141E-05  0.0000E+00 -6.4372E-02  1.4479E-01
            -7.3983E-02

0ITERATION NO.:   64    OBJECTIVE VALUE:  -1700.23726370922        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     2317
 NPARAMETR:  9.8628E-01  1.4661E-01  1.3182E+00  1.5695E+00  8.5817E-01  9.9692E-01  1.0000E-02  1.0000E-02  7.4664E-01  1.0659E+00
             1.0008E+00
 PARAMETER:  8.6180E-02 -1.8200E+00  3.7630E-01  5.5075E-01 -5.2955E-02  9.6912E-02 -4.6190E+00 -5.1156E+00 -1.9217E-01  1.6378E-01
             1.0078E-01
 GRADIENT:   1.4776E+00  3.9161E-02 -2.0364E-01 -1.5035E+01  1.3903E+00  1.0872E-01  0.0000E+00  0.0000E+00  1.1519E-01  3.1823E-03
             9.6327E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2317
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.3431E-05 -5.2520E-05 -1.9582E-04 -5.7958E-03 -2.4851E-02
 SE:             2.9827E-02  2.7135E-05  1.6807E-04  2.9201E-02  2.5174E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9884E-01  5.2926E-02  2.4399E-01  8.4267E-01  3.2356E-01

 ETASHRINKSD(%)  7.4883E-02  9.9909E+01  9.9437E+01  2.1737E+00  1.5664E+01
 ETASHRINKVR(%)  1.4971E-01  1.0000E+02  9.9997E+01  4.3001E+00  2.8875E+01
 EBVSHRINKSD(%)  4.1068E-01  9.9915E+01  9.9440E+01  2.2696E+00  1.2479E+01
 EBVSHRINKVR(%)  8.1967E-01  1.0000E+02  9.9997E+01  4.4877E+00  2.3401E+01
 RELATIVEINF(%)  9.7042E+01  3.4725E-06  3.7779E-04  5.8298E+00  6.1814E+00
 EPSSHRINKSD(%)  4.1546E+01
 EPSSHRINKVR(%)  6.5831E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1700.2372637092221     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -965.08643714548396     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.33
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1700.237       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.47E-01  1.32E+00  1.57E+00  8.58E-01  9.97E-01  1.00E-02  1.00E-02  7.47E-01  1.07E+00  1.00E+00
 


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
+        1.14E+03
 
 TH 2
+       -2.86E+01  3.86E+02
 
 TH 3
+        1.12E+00  8.37E+01  2.04E+02
 
 TH 4
+       -9.94E+00  4.68E+02 -4.63E+01  7.73E+02
 
 TH 5
+        2.29E+00 -2.91E+02 -4.14E+02 -3.44E+01  1.03E+03
 
 TH 6
+        2.29E-01 -4.81E+00  2.79E-01 -2.62E+00 -8.31E-01  1.97E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.87E+00 -9.83E+01  5.57E+00 -1.91E+00 -2.97E+00 -6.44E-01  0.00E+00  0.00E+00  3.27E+02
 
 TH10
+        3.67E-01  1.15E+01 -7.65E+00 -1.84E+00 -7.15E+01  6.71E-01  0.00E+00  0.00E+00  3.42E-01  9.50E+01
 
 TH11
+       -8.29E+00 -1.80E+01 -3.01E+01 -1.10E+01  2.19E+01  1.53E+00  0.00E+00  0.00E+00  1.03E+01  3.10E+01  2.30E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.178
Stop Time:
Wed Sep 29 19:05:20 CDT 2021
