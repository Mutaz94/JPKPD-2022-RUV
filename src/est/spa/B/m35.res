Wed Sep 29 11:09:35 CDT 2021
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
$DATA ../../../../data/spa/B/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1707.06389970836        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4552E+02  9.8732E+00 -4.1313E+01  8.1039E+01  6.7670E+01  6.5528E+01  8.1460E-01  1.1887E+01  7.1495E+00 -1.2028E+01
             1.2822E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1714.90705885663        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0553E+00  1.0258E+00  1.1564E+00  1.0205E+00  1.0360E+00  9.4136E-01  1.0116E+00  8.6010E-01  1.0051E+00  1.1379E+00
             9.5658E-01
 PARAMETER:  1.5384E-01  1.2551E-01  2.4535E-01  1.2027E-01  1.3541E-01  3.9568E-02  1.1156E-01 -5.0707E-02  1.0505E-01  2.2919E-01
             5.5609E-02
 GRADIENT:  -1.1425E+00  2.6284E+01 -9.0205E-02  3.1686E+01 -5.9177E+00 -1.0390E+01  1.2536E+00  2.4235E+00 -9.8146E-01 -4.1377E+00
            -6.7614E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1715.84309933220        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0562E+00  8.7922E-01  1.1769E+00  1.0936E+00  9.8318E-01  9.5919E-01  9.7511E-01  6.7486E-01  9.8405E-01  1.1529E+00
             9.7397E-01
 PARAMETER:  1.5466E-01 -2.8719E-02  2.6287E-01  1.8951E-01  8.3037E-02  5.8332E-02  7.4798E-02 -2.9325E-01  8.3924E-02  2.4232E-01
             7.3623E-02
 GRADIENT:   3.2465E+00  9.1354E+00  7.5722E+00  6.9061E+00 -1.4066E+01 -2.2205E+00 -6.1750E-01 -1.6218E+00  2.5772E+00  1.4841E-02
             5.6416E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1716.35019293005        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0516E+00  7.4672E-01  1.4235E+00  1.1884E+00  1.0373E+00  9.6155E-01  1.0248E+00  9.6022E-01  9.3332E-01  1.2049E+00
             9.7153E-01
 PARAMETER:  1.5030E-01 -1.9206E-01  4.5310E-01  2.7257E-01  1.3665E-01  6.0788E-02  1.2451E-01  5.9403E-02  3.0994E-02  2.8639E-01
             7.1114E-02
 GRADIENT:  -2.1614E+00  7.2154E+00  9.3118E-01  1.1653E+01 -3.5593E-02 -3.3840E-01  8.3290E-01 -6.7611E-03  1.5665E+00 -1.8487E+00
             4.9252E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1716.93399738118        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0503E+00  4.6051E-01  1.6331E+00  1.3669E+00  1.0067E+00  9.6052E-01  7.3840E-01  1.0916E+00  8.5634E-01  1.2319E+00
             9.7028E-01
 PARAMETER:  1.4908E-01 -6.7543E-01  5.9051E-01  4.1255E-01  1.0671E-01  5.9719E-02 -2.0327E-01  1.8761E-01 -5.5090E-02  3.0854E-01
             6.9834E-02
 GRADIENT:   3.3882E+00  3.9729E+00  3.5100E+00  3.4250E+00 -5.5375E+00  7.3236E-01 -1.1498E-01 -7.0487E-01 -1.2141E+00  2.2230E-01
            -3.9138E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1716.95217206370        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0486E+00  3.8933E-01  1.6889E+00  1.4183E+00  1.0008E+00  9.5897E-01  6.2025E-01  1.1421E+00  8.3785E-01  1.2326E+00
             9.7033E-01
 PARAMETER:  1.4745E-01 -8.4333E-01  6.2405E-01  4.4949E-01  1.0082E-01  5.8106E-02 -3.7763E-01  2.3290E-01 -7.6913E-02  3.0913E-01
             6.9884E-02
 GRADIENT:   1.5491E+00  5.4252E+00  3.2219E+00  1.3747E+01 -6.4222E+00  4.2456E-01 -6.2846E-02 -5.3660E-01 -1.5085E-01  3.5076E-01
            -4.2768E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1716.95536975215        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  1.0476E+00  3.5389E-01  1.7154E+00  1.4435E+00  9.9816E-01  9.5807E-01  5.5176E-01  1.1680E+00  8.2804E-01  1.2324E+00
             9.7040E-01
 PARAMETER:  1.4655E-01 -9.3876E-01  6.3963E-01  4.6705E-01  9.8157E-02  5.7171E-02 -4.9465E-01  2.5529E-01 -8.8699E-02  3.0896E-01
             6.9948E-02
 GRADIENT:   4.5236E-01  5.6025E+00  2.7234E+00  1.7665E+01 -6.2250E+00  2.1502E-01 -3.4505E-02 -3.8262E-01  4.0418E-01  3.9498E-01
            -3.8001E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1716.95984836290        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  1.0468E+00  3.2183E-01  1.7394E+00  1.4659E+00  9.9608E-01  9.5722E-01  4.8613E-01  1.1923E+00  8.1882E-01  1.2321E+00
             9.7045E-01
 PARAMETER:  1.4571E-01 -1.0337E+00  6.5352E-01  4.8247E-01  9.6076E-02  5.6283E-02 -6.2128E-01  2.7586E-01 -9.9888E-02  3.0876E-01
             7.0003E-02
 GRADIENT:  -5.9767E-01  5.4951E+00  2.1023E+00  2.0362E+01 -5.7211E+00  3.4296E-03 -1.3819E-02 -2.0090E-01  9.4964E-01  4.3092E-01
            -3.0651E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1716.96271130747        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1432
 NPARAMETR:  1.0459E+00  2.8915E-01  1.7619E+00  1.4887E+00  9.9361E-01  9.5640E-01  4.1857E-01  1.2152E+00  8.0869E-01  1.2311E+00
             9.7047E-01
 PARAMETER:  1.4485E-01 -1.1408E+00  6.6640E-01  4.9788E-01  9.3590E-02  5.5426E-02 -7.7092E-01  2.9492E-01 -1.1233E-01  3.0789E-01
             7.0027E-02
 GRADIENT:  -1.6503E+00  5.2818E+00  1.4106E+00  2.2851E+01 -5.0417E+00 -1.9463E-01 -2.1008E-03 -3.9010E-02  1.3307E+00  4.0885E-01
            -2.5162E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1716.96885087494        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1607
 NPARAMETR:  1.0448E+00  2.4847E-01  1.7917E+00  1.5167E+00  9.9127E-01  9.5538E-01  3.3407E-01  1.2458E+00  7.9578E-01  1.2302E+00
             9.7050E-01
 PARAMETER:  1.4378E-01 -1.2924E+00  6.8319E-01  5.1653E-01  9.1231E-02  5.4353E-02 -9.9640E-01  3.1979E-01 -1.2843E-01  3.0714E-01
             7.0052E-02
 GRADIENT:  -2.9531E+00  4.7834E+00  4.4257E-01  2.4868E+01 -3.9485E+00 -4.5007E-01  4.1808E-03  2.0145E-01  1.8194E+00  3.8687E-01
            -1.5388E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1716.97649877767        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1783
 NPARAMETR:  1.0437E+00  2.1026E-01  1.8204E+00  1.5427E+00  9.8941E-01  9.5446E-01  2.5723E-01  1.2746E+00  7.8323E-01  1.2293E+00
             9.7049E-01
 PARAMETER:  1.4281E-01 -1.4594E+00  6.9907E-01  5.3355E-01  8.9354E-02  5.3389E-02 -1.2578E+00  3.4262E-01 -1.4433E-01  3.0641E-01
             7.0049E-02
 GRADIENT:  -4.0815E+00  4.1815E+00 -4.7822E-01  2.5911E+01 -2.7666E+00 -6.7172E-01  4.2061E-03  4.0594E-01  2.1412E+00  3.3052E-01
            -7.0896E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1717.00078669074        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1958
 NPARAMETR:  1.0426E+00  1.6235E-01  1.8606E+00  1.5747E+00  9.8819E-01  9.5350E-01  1.6771E-01  1.3124E+00  7.6715E-01  1.2290E+00
             9.7047E-01
 PARAMETER:  1.4175E-01 -1.7180E+00  7.2089E-01  5.5409E-01  8.8122E-02  5.2384E-02 -1.6855E+00  3.7186E-01 -1.6507E-01  3.0617E-01
             7.0026E-02
 GRADIENT:  -5.0680E+00  3.2484E+00 -1.4752E+00  2.5182E+01 -1.2168E+00 -8.6852E-01  2.0009E-03  6.0031E-01  2.3721E+00  2.3928E-01
             2.7659E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1717.04357763980        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2133
 NPARAMETR:  1.0419E+00  1.1205E-01  1.9098E+00  1.6079E+00  9.8806E-01  9.5290E-01  8.7921E-02  1.3533E+00  7.4912E-01  1.2295E+00
             9.7039E-01
 PARAMETER:  1.4101E-01 -2.0888E+00  7.4699E-01  5.7496E-01  8.7991E-02  5.1750E-02 -2.3313E+00  4.0254E-01 -1.8885E-01  3.0658E-01
             6.9938E-02
 GRADIENT:  -5.1579E+00  2.2429E+00 -1.9751E+00  2.2747E+01 -4.9980E-02 -8.9889E-01  4.1474E-04  6.0047E-01  1.9046E+00  8.8445E-02
             3.7325E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1717.13297159794        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2308
 NPARAMETR:  1.0416E+00  6.3801E-02  1.9711E+00  1.6383E+00  9.9105E-01  9.5273E-01  3.1054E-02  1.3961E+00  7.3232E-01  1.2323E+00
             9.7027E-01
 PARAMETER:  1.4071E-01 -2.6520E+00  7.7858E-01  5.9365E-01  9.1005E-02  5.1577E-02 -3.3720E+00  4.3367E-01 -2.1153E-01  3.0884E-01
             6.9822E-02
 GRADIENT:  -4.2633E+00  1.1473E+00 -1.7614E+00  1.5168E+01  1.1590E+00 -7.5699E-01  2.8851E-05  2.6851E-01  1.5702E+00 -1.6407E-01
             2.5882E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1717.18104684991        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2483
 NPARAMETR:  1.0419E+00  3.7327E-02  2.0256E+00  1.6553E+00  9.9606E-01  9.5313E-01  1.1145E-02  1.4363E+00  7.2180E-01  1.2374E+00
             9.7020E-01
 PARAMETER:  1.4103E-01 -3.1880E+00  8.0587E-01  6.0400E-01  9.6054E-02  5.1996E-02 -4.3968E+00  4.6206E-01 -2.2601E-01  3.1300E-01
             6.9747E-02
 GRADIENT:  -2.6387E+00  6.3291E-01 -1.1166E+00  9.8514E+00  6.5341E-01 -4.7472E-01  1.8780E-06  2.0881E-01  7.2412E-01 -9.7343E-02
             4.0256E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1717.37786324869        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2672             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0443E+00  1.6277E-02  2.0430E+00  1.6501E+00  9.9403E-01  9.5433E-01  1.0000E-02  1.4455E+00  7.1381E-01  1.2391E+00
             9.6956E-01
 PARAMETER:  1.4338E-01 -4.0180E+00  8.1444E-01  6.0083E-01  9.4010E-02  5.3251E-02 -4.7238E+00  4.6842E-01 -2.3713E-01  3.1436E-01
             6.9091E-02
 GRADIENT:   7.1742E+02  4.7607E-01  1.0270E+01  1.3805E+03  8.2773E+00  5.2035E+01  0.0000E+00  1.4138E+00  2.2625E+01  4.1816E+00
             9.1368E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1717.39243897564        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2858
 NPARAMETR:  1.0432E+00  1.2447E-02  2.0362E+00  1.6575E+00  9.9172E-01  9.5411E-01  1.0000E-02  1.4410E+00  7.1315E-01  1.2358E+00
             9.6950E-01
 PARAMETER:  1.4228E-01 -4.2863E+00  8.1106E-01  6.0534E-01  9.1685E-02  5.3024E-02 -4.7238E+00  4.6535E-01 -2.3806E-01  3.1174E-01
             6.9025E-02
 GRADIENT:   1.7397E+00  4.3933E-02  1.6116E-01 -2.4923E+01  1.3913E+00  1.1576E-01  0.0000E+00  1.1404E-01  4.6861E-01  8.3080E-02
            -3.1745E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1717.40436295326        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3053             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0440E+00  1.0000E-02  2.0258E+00  1.6567E+00  9.8483E-01  9.5414E-01  1.0000E-02  1.4288E+00  7.1184E-01  1.2327E+00
             9.6962E-01
 PARAMETER:  1.4309E-01 -4.6868E+00  8.0599E-01  6.0483E-01  8.4715E-02  5.3060E-02 -4.7238E+00  4.5682E-01 -2.3990E-01  3.0923E-01
             6.9148E-02
 GRADIENT:   7.1482E+02  0.0000E+00  1.1563E+01  1.4047E+03  4.8076E+00  5.1915E+01  0.0000E+00  1.1005E+00  2.2917E+01  4.2175E+00
             7.4243E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1717.40994312506        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3248
 NPARAMETR:  1.0440E+00  1.0000E-02  2.0168E+00  1.6563E+00  9.8313E-01  9.5413E-01  1.0000E-02  1.4222E+00  7.1193E-01  1.2309E+00
             9.6958E-01
 PARAMETER:  1.4307E-01 -4.6868E+00  8.0151E-01  6.0461E-01  8.2991E-02  5.3040E-02 -4.7238E+00  4.5222E-01 -2.3978E-01  3.0777E-01
             6.9106E-02
 GRADIENT:   3.8392E+00  0.0000E+00  1.7605E+00 -2.9105E+01 -1.6896E+00  1.2721E-01  0.0000E+00 -1.0860E-01  1.0991E-01  3.9380E-01
            -8.2922E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1717.41796836260        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     3448             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0440E+00  1.0000E-02  2.0020E+00  1.6559E+00  9.8274E-01  9.5412E-01  1.0000E-02  1.4178E+00  7.1202E-01  1.2280E+00
             9.6968E-01
 PARAMETER:  1.4305E-01 -4.6868E+00  7.9416E-01  6.0434E-01  8.2585E-02  5.3034E-02 -4.7238E+00  4.4914E-01 -2.3966E-01  3.0539E-01
             6.9207E-02
 GRADIENT:   7.1408E+02  0.0000E+00  9.7636E+00  1.4020E+03  7.9580E+00  5.1852E+01  0.0000E+00  1.3808E+00  2.2832E+01  3.7670E+00
             8.1503E-01

0ITERATION NO.:   99    OBJECTIVE VALUE:  -1717.41949022771        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     3594
 NPARAMETR:  1.0440E+00  1.0000E-02  1.9849E+00  1.6553E+00  9.8171E-01  9.5410E-01  1.0000E-02  1.4117E+00  7.1212E-01  1.2244E+00
             9.6982E-01
 PARAMETER:  1.4304E-01 -4.6868E+00  7.9353E-01  6.0422E-01  8.0728E-02  5.3029E-02 -4.7238E+00  4.4632E-01 -2.3955E-01  3.0520E-01
             6.9168E-02
 GRADIENT:   1.1574E-02  0.0000E+00  9.5257E-01  3.9680E-01 -3.2016E-01  1.9259E-03  0.0000E+00  4.9071E-02 -7.1311E-03  1.9926E-01
            -3.4364E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3594
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.6409E-04 -3.7922E-06 -3.1781E-02 -7.8859E-03 -4.0795E-02
 SE:             2.9861E-02  1.7873E-06  1.5720E-02  2.9219E-02  2.2023E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9294E-01  3.3862E-02  4.3215E-02  7.8724E-01  6.3971E-02

 ETASHRINKSD(%)  1.0000E-10  9.9994E+01  4.7335E+01  2.1130E+00  2.6220E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  7.2264E+01  4.1813E+00  4.5564E+01
 EBVSHRINKSD(%)  4.1217E-01  9.9994E+01  5.1666E+01  2.4041E+00  2.1300E+01
 EBVSHRINKVR(%)  8.2263E-01  1.0000E+02  7.6639E+01  4.7503E+00  3.8064E+01
 RELATIVEINF(%)  9.7318E+01  1.6971E-08  6.3732E+00  5.6907E+00  1.2853E+01
 EPSSHRINKSD(%)  4.4539E+01
 EPSSHRINKVR(%)  6.9240E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1717.4194902277077     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -982.26866366396951     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    47.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1717.419       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  2.00E+00  1.66E+00  9.81E-01  9.54E-01  1.00E-02  1.41E+00  7.12E-01  1.23E+00  9.70E-01
 


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
+        1.11E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.98E+00  0.00E+00  3.59E+01
 
 TH 4
+       -1.15E+01  0.00E+00 -2.09E+01  7.62E+02
 
 TH 5
+       -1.12E+00  0.00E+00 -9.62E+01 -4.37E+01  4.89E+02
 
 TH 6
+       -3.29E-01  0.00E+00 -1.33E-01 -1.89E+00  3.19E-01  2.16E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        3.11E-01  0.00E+00 -1.45E+01 -3.18E+00 -4.28E+00 -6.81E-02  0.00E+00  1.96E+01
 
 TH 9
+        2.16E+00  0.00E+00  4.62E+00 -1.00E+00 -1.73E-01 -1.14E+00  0.00E+00 -2.65E-01  3.60E+02
 
 TH10
+        6.77E-01  0.00E+00 -2.12E+00 -2.05E+00 -5.80E+01  2.89E-01  0.00E+00  6.99E+00  2.37E-01  5.68E+01
 
 TH11
+       -8.81E+00  0.00E+00 -5.23E+00 -1.13E+01  3.35E+00  2.67E+00  0.00E+00  4.66E+00  1.21E+01  1.13E+01  2.29E+02
 
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
 #CPUT: Total CPU Time in Seconds,       54.002
Stop Time:
Wed Sep 29 11:10:31 CDT 2021
