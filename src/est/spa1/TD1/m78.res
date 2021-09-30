Thu Sep 30 01:38:03 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat78.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2126.80475444808        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5697E+02 -7.4132E+01 -5.2237E+01 -2.9515E+01  3.6245E+01  3.3047E+01 -1.5508E+01  1.8189E+01  1.1263E+01  3.8748E+00
             2.8070E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2143.78335910455        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.7500E-01  1.1965E+00  1.3219E+00  9.7543E-01  1.2058E+00  9.9128E-01  1.1392E+00  8.5752E-01  9.2278E-01  9.9990E-01
             9.2605E-01
 PARAMETER:  7.4684E-02  2.7944E-01  3.7907E-01  7.5124E-02  2.8717E-01  9.1243E-02  2.3030E-01 -5.3709E-02  1.9638E-02  9.9898E-02
             2.3173E-02
 GRADIENT:   2.5049E+00  6.4275E+00  4.5689E-01  1.1133E+01  1.3526E+01 -1.9559E+00 -7.1748E+00  2.1587E+00 -1.1265E-01 -3.7477E+01
            -4.7758E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2146.09451414470        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  9.7504E-01  1.1741E+00  1.9263E+00  9.8040E-01  1.3907E+00  1.0370E+00  1.2452E+00  6.5680E-01  8.9563E-01  1.5851E+00
             9.6696E-01
 PARAMETER:  7.4721E-02  2.6054E-01  7.5561E-01  8.0206E-02  4.2981E-01  1.3632E-01  3.1930E-01 -3.2038E-01 -1.0223E-02  5.6065E-01
             6.6401E-02
 GRADIENT:   4.3165E+00 -9.6757E+00  8.5442E+00 -2.8291E+01 -2.0606E+01  1.5836E+01  5.5274E+00 -8.7434E-01  1.8308E-01  2.4259E+01
            -1.1999E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2149.69744249925        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.7032E-01  9.7397E-01  1.6092E+00  1.1219E+00  1.2245E+00  9.9103E-01  1.3893E+00  5.2590E-01  8.4408E-01  1.2533E+00
             9.7539E-01
 PARAMETER:  6.9866E-02  7.3629E-02  5.7573E-01  2.1502E-01  3.0254E-01  9.0986E-02  4.2877E-01 -5.4265E-01 -6.9514E-02  3.2581E-01
             7.5084E-02
 GRADIENT:  -4.3464E+00  5.6169E+00  2.6918E-01  5.9770E+00 -1.2539E+00 -7.9095E-01  1.5596E+00 -4.5840E-01  4.4574E-01 -1.2134E+00
            -9.0368E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2150.15910066180        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.7036E-01  8.7900E-01  2.1802E+00  1.1967E+00  1.3306E+00  9.9108E-01  1.3875E+00  1.0997E+00  8.5920E-01  1.3509E+00
             9.7314E-01
 PARAMETER:  6.9913E-02 -2.8975E-02  8.7944E-01  2.7960E-01  3.8560E-01  9.1045E-02  4.2749E-01  1.9506E-01 -5.1752E-02  4.0075E-01
             7.2775E-02
 GRADIENT:  -2.8207E+00  4.4255E+00  1.7670E+00  2.7590E+00 -1.8052E+00 -8.3697E-01  1.0200E+00 -6.9268E-01  3.1417E-01  4.1697E-01
            -1.8806E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2150.32527996905        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.7094E-01  6.8712E-01  2.5940E+00  1.3340E+00  1.3331E+00  9.9281E-01  1.3971E+00  1.3045E+00  8.5405E-01  1.3674E+00
             9.7456E-01
 PARAMETER:  7.0505E-02 -2.7525E-01  1.0532E+00  3.8820E-01  3.8749E-01  9.2788E-02  4.3441E-01  3.6584E-01 -5.7766E-02  4.1294E-01
             7.4235E-02
 GRADIENT:   7.7611E-01  6.4119E+00  1.4959E+00  1.1021E+01 -5.3517E+00  8.7586E-02 -4.9093E-01  1.1622E-01 -9.1501E-01  5.3612E-01
            -1.6869E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2150.42159855800        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.7086E-01  5.2082E-01  2.8746E+00  1.4475E+00  1.3285E+00  9.9319E-01  1.3465E+00  1.4024E+00  8.4439E-01  1.3689E+00
             9.7711E-01
 PARAMETER:  7.0428E-02 -5.5236E-01  1.1559E+00  4.6983E-01  3.8405E-01  9.3163E-02  3.9753E-01  4.3819E-01 -6.9144E-02  4.1403E-01
             7.6845E-02
 GRADIENT:   2.6982E+00  5.5596E+00  9.8539E-01  1.3872E+01 -3.3961E+00  4.7120E-01 -2.3305E-01  2.2759E-02 -4.7235E-01 -3.6810E-01
             1.3738E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2150.44279301923        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  9.7029E-01  4.2933E-01  2.9732E+00  1.5101E+00  1.3177E+00  9.9270E-01  1.2904E+00  1.4066E+00  8.3248E-01  1.3665E+00
             9.7752E-01
 PARAMETER:  6.9843E-02 -7.4552E-01  1.1896E+00  5.1217E-01  3.7591E-01  9.2673E-02  3.5494E-01  4.4116E-01 -8.3344E-02  4.1228E-01
             7.7267E-02
 GRADIENT:   2.7585E+00  5.5873E+00  1.2625E+00  1.7788E+01 -2.8890E+00  4.6453E-01 -1.3000E-01 -5.9839E-01 -4.8039E-01 -5.9572E-01
             3.7338E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2150.50155496503        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  9.6919E-01  3.1934E-01  3.0443E+00  1.5830E+00  1.3005E+00  9.9170E-01  1.1779E+00  1.3942E+00  8.1273E-01  1.3626E+00
             9.7764E-01
 PARAMETER:  6.8701E-02 -1.0415E+00  1.2133E+00  5.5934E-01  3.6274E-01  9.1669E-02  2.6372E-01  4.3229E-01 -1.0736E-01  4.0937E-01
             7.7384E-02
 GRADIENT:   2.1266E+00  4.6548E+00  1.0678E+00  1.8944E+01 -1.6821E+00  3.5951E-01 -2.6941E-02 -1.0541E+00 -7.7026E-01 -6.1847E-01
             4.9125E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2150.57988223525        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1585
 NPARAMETR:  9.6781E-01  2.1489E-01  3.0797E+00  1.6510E+00  1.2812E+00  9.9054E-01  1.0042E+00  1.3739E+00  7.9027E-01  1.3588E+00
             9.7751E-01
 PARAMETER:  6.7278E-02 -1.4376E+00  1.2248E+00  6.0139E-01  3.4777E-01  9.0495E-02  1.0423E-01  4.1764E-01 -1.3538E-01  4.0659E-01
             7.7253E-02
 GRADIENT:   1.0381E+00  3.1721E+00  4.9552E-01  1.6826E+01 -2.8314E-01  2.3373E-01  3.2276E-02 -1.2143E+00 -4.9304E-01 -3.5621E-01
             6.0091E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2150.79892847826        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1749
 NPARAMETR:  9.6681E-01  1.7924E-01  3.0765E+00  1.6597E+00  1.2747E+00  9.8972E-01  2.8158E-01  1.4001E+00  7.8382E-01  1.3578E+00
             9.7666E-01
 PARAMETER:  6.6243E-02 -1.6190E+00  1.2238E+00  6.0662E-01  3.4269E-01  8.9669E-02 -1.1674E+00  4.3654E-01 -1.4357E-01  4.0584E-01
             7.6379E-02
 GRADIENT:  -2.8870E-01  4.6636E-01 -3.7301E-01 -1.2974E+01  1.7577E+00  8.3464E-02  8.5985E-03 -3.1640E-01 -5.8970E-01 -6.4101E-02
             7.9574E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2150.80918006685        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1931
 NPARAMETR:  9.6727E-01  1.7159E-01  3.0887E+00  1.6636E+00  1.2717E+00  9.8963E-01  8.6384E-02  1.4181E+00  7.8389E-01  1.3574E+00
             9.7557E-01
 PARAMETER:  6.6727E-02 -1.6627E+00  1.2278E+00  6.0900E-01  3.4038E-01  8.9581E-02 -2.3490E+00  4.4928E-01 -1.4348E-01  4.0557E-01
             7.5272E-02
 GRADIENT:   9.6983E-01  2.4428E-01 -5.5909E-02 -1.5873E+01  1.6283E-02  6.8057E-02  1.2208E-03  1.7429E-02  4.6589E-01  1.7702E-01
             2.4042E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2150.80999339609        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2119             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6754E-01  1.7129E-01  3.0864E+00  1.6635E+00  1.2711E+00  9.8966E-01  3.1752E-02  1.4165E+00  7.8299E-01  1.3561E+00
             9.7558E-01
 PARAMETER:  6.6999E-02 -1.6644E+00  1.2270E+00  6.0891E-01  3.3991E-01  8.9607E-02 -3.3498E+00  4.4822E-01 -1.4463E-01  4.0459E-01
             7.5276E-02
 GRADIENT:   4.2801E+02  2.0812E+01  5.7025E+00  1.1483E+03  1.6112E+01  3.6720E+01  4.0794E-03  4.3146E-01  1.7701E+01  5.5375E+00
             9.8807E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2150.81008315395        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2305
 NPARAMETR:  9.6753E-01  1.7109E-01  3.0842E+00  1.6635E+00  1.2708E+00  9.8965E-01  1.7825E-02  1.4154E+00  7.8298E-01  1.3559E+00
             9.7558E-01
 PARAMETER:  6.6994E-02 -1.6656E+00  1.2263E+00  6.0892E-01  3.3964E-01  8.9599E-02 -3.9271E+00  4.4739E-01 -1.4465E-01  4.0447E-01
             7.5279E-02
 GRADIENT:   1.5964E+00  2.1917E-01 -2.2315E-02 -1.6628E+01  2.1825E-02  8.6166E-02  8.6742E-05  7.1763E-03  7.9023E-02  4.0174E-02
            -2.4662E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2150.81012018825        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     2474
 NPARAMETR:  9.6753E-01  1.7099E-01  3.0829E+00  1.6635E+00  1.2706E+00  9.8965E-01  1.5954E-02  1.4148E+00  7.8296E-01  1.3558E+00
             9.7559E-01
 PARAMETER:  6.6991E-02 -1.6661E+00  1.2259E+00  6.0895E-01  3.3948E-01  8.9594E-02 -4.0381E+00  4.4698E-01 -1.4468E-01  4.0440E-01
             7.5284E-02
 GRADIENT:   1.0334E-02 -1.0841E-02  2.5279E-02 -2.4109E-01  9.8549E-02  1.3886E-03  1.2016E-05  8.7301E-03  2.6644E-02  2.2261E-03
             6.9211E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2474
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.7799E-04 -1.0809E-04 -2.1226E-02 -4.8621E-03 -3.9785E-02
 SE:             2.9823E-02  4.6016E-05  1.1280E-02  2.9465E-02  2.3097E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8186E-01  1.8828E-02  5.9874E-02  8.6893E-01  8.4980E-02

 ETASHRINKSD(%)  9.0159E-02  9.9846E+01  6.2210E+01  1.2895E+00  2.2621E+01
 ETASHRINKVR(%)  1.8024E-01  1.0000E+02  8.5719E+01  2.5625E+00  4.0124E+01
 EBVSHRINKSD(%)  3.6023E-01  9.9855E+01  6.5737E+01  1.4170E+00  1.8229E+01
 EBVSHRINKVR(%)  7.1917E-01  1.0000E+02  8.8261E+01  2.8139E+00  3.3136E+01
 RELATIVEINF(%)  9.8409E+01  1.1428E-05  3.4069E+00  5.4457E+00  1.8505E+01
 EPSSHRINKSD(%)  3.2117E+01
 EPSSHRINKVR(%)  5.3919E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2150.8101201882537     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1231.8715869835810     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.94
 Elapsed covariance  time in seconds:     7.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2150.810       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  1.71E-01  3.08E+00  1.66E+00  1.27E+00  9.90E-01  1.60E-02  1.41E+00  7.83E-01  1.36E+00  9.76E-01
 


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
 
         2.87E-02  6.99E-01  1.32E+00  4.53E-01  3.17E-01  7.71E-02  1.67E-01  8.39E-01  2.21E-01  2.28E-01  4.53E-02
 


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
+        8.23E-04
 
 TH 2
+        2.53E-03  4.89E-01
 
 TH 3
+       -4.61E-05  2.82E-01  1.75E+00
 
 TH 4
+       -1.77E-03 -3.15E-01 -1.40E-01  2.05E-01
 
 TH 5
+        5.83E-04  1.72E-01  3.38E-01 -1.04E-01  1.00E-01
 
 TH 6
+       -6.10E-04  1.67E-02  3.73E-02 -1.02E-02  9.36E-03  5.95E-03
 
 TH 7
+        6.51E-04  1.09E-01  1.29E-03 -7.16E-02  2.88E-02  2.00E-03  2.78E-02
 
 TH 8
+        7.51E-05  3.23E-01  9.60E-01 -1.88E-01  2.31E-01  2.45E-02  3.71E-02  7.04E-01
 
 TH 9
+        6.08E-04  1.50E-01  6.69E-02 -9.71E-02  4.99E-02  4.75E-03  3.52E-02  8.82E-02  4.88E-02
 
 TH10
+       -9.00E-04  7.73E-02  1.93E-01 -4.64E-02  5.18E-02  8.17E-03  9.09E-03  1.21E-01  2.18E-02  5.21E-02
 
 TH11
+        3.82E-04  1.10E-02 -4.03E-03 -7.39E-03  2.36E-03 -1.17E-04  2.46E-03 -1.32E-04  3.69E-03 -2.04E-04  2.05E-03
 
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
+        2.87E-02
 
 TH 2
+        1.26E-01  6.99E-01
 
 TH 3
+       -1.21E-03  3.04E-01  1.32E+00
 
 TH 4
+       -1.36E-01 -9.94E-01 -2.33E-01  4.53E-01
 
 TH 5
+        6.41E-02  7.77E-01  8.06E-01 -7.27E-01  3.17E-01
 
 TH 6
+       -2.76E-01  3.10E-01  3.65E-01 -2.91E-01  3.83E-01  7.71E-02
 
 TH 7
+        1.36E-01  9.34E-01  5.86E-03 -9.48E-01  5.45E-01  1.55E-01  1.67E-01
 
 TH 8
+        3.12E-03  5.51E-01  8.64E-01 -4.94E-01  8.70E-01  3.78E-01  2.65E-01  8.39E-01
 
 TH 9
+        9.59E-02  9.72E-01  2.29E-01 -9.70E-01  7.13E-01  2.79E-01  9.57E-01  4.76E-01  2.21E-01
 
 TH10
+       -1.37E-01  4.84E-01  6.39E-01 -4.49E-01  7.16E-01  4.64E-01  2.39E-01  6.30E-01  4.32E-01  2.28E-01
 
 TH11
+        2.94E-01  3.47E-01 -6.71E-02 -3.60E-01  1.65E-01 -3.34E-02  3.25E-01 -3.49E-03  3.69E-01 -1.98E-02  4.53E-02
 
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
+        2.07E+03
 
 TH 2
+        3.23E+03  2.31E+04
 
 TH 3
+       -1.47E+01  2.95E+01  9.02E+00
 
 TH 4
+        6.09E+02  3.87E+03 -1.15E+01  1.19E+03
 
 TH 5
+       -2.12E+02 -1.84E+03 -4.06E+01 -2.93E+02  4.20E+02
 
 TH 6
+       -1.57E+03 -1.24E+04 -2.73E+01 -1.85E+03  1.01E+03  7.05E+03
 
 TH 7
+       -1.83E+04 -1.29E+05 -1.63E+02 -1.94E+04  9.89E+03  7.05E+04  7.32E+05
 
 TH 8
+       -9.02E+02 -6.43E+03 -1.21E+01 -9.65E+02  4.89E+02  3.52E+03  3.65E+04  1.83E+03
 
 TH 9
+        7.85E+03  5.42E+04  7.06E+01  8.16E+03 -4.17E+03 -2.97E+04 -3.08E+05 -1.54E+04  1.30E+05
 
 TH10
+       -1.75E+03 -1.26E+04 -1.64E+01 -1.88E+03  9.35E+02  6.87E+03  7.16E+04  3.57E+03 -3.02E+04  7.05E+03
 
 TH11
+       -7.88E+03 -5.39E+04 -6.82E+01 -8.09E+03  4.10E+03  2.95E+04  3.06E+05  1.53E+04 -1.29E+05  2.99E+04  1.28E+05
 
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
 #CPUT: Total CPU Time in Seconds,       47.365
Stop Time:
Thu Sep 30 01:38:52 CDT 2021
