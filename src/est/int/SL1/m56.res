Sat Sep 25 00:22:21 CDT 2021
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
$DATA ../../../../data/int/SL1/dat56.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2763.89276881410        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3123E+01 -9.8149E+01  1.8104E+02  2.5765E+01  8.3713E+01 -1.7071E+01 -3.2806E+01 -1.7114E+02 -7.9498E+01  2.1117E+01
            -1.9845E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3265.00686298809        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0282E+00  1.3010E+00  9.8438E-01  8.6927E-01  1.1191E+00  1.0189E+00  9.4241E-01  1.2806E+00  1.1970E+00  7.3728E-01
             1.8012E+00
 PARAMETER:  1.2779E-01  3.6316E-01  8.4260E-02 -4.0099E-02  2.1250E-01  1.1872E-01  4.0684E-02  3.4732E-01  2.7983E-01 -2.0478E-01
             6.8845E-01
 GRADIENT:   4.4819E+01  2.1553E+01  1.4194E+01  1.0217E+00 -3.7590E+01 -1.0297E+01 -5.8540E+00 -7.6821E+00  5.3136E+00 -1.1588E+01
             6.1668E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3268.83838077737        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0372E+00  1.6198E+00  1.2808E+00  7.4520E-01  1.4539E+00  1.1513E+00  9.9354E-01  2.3376E+00  1.2977E+00  9.3845E-01
             1.7494E+00
 PARAMETER:  1.3656E-01  5.8233E-01  3.4751E-01 -1.9410E-01  4.7422E-01  2.4085E-01  9.3517E-02  9.4915E-01  3.6057E-01  3.6477E-02
             6.5925E-01
 GRADIENT:   5.6376E+01  8.5263E+01 -8.1869E+00  6.1452E+01  3.6146E+01  3.6236E+01  2.5060E+01 -4.0297E+00  5.9186E+00 -5.2176E+00
             2.2084E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3272.90524092288        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0241E+00  1.5608E+00  1.2470E+00  7.5607E-01  1.4075E+00  1.1072E+00  9.2632E-01  2.3321E+00  1.3339E+00  9.6960E-01
             1.7366E+00
 PARAMETER:  1.2384E-01  5.4517E-01  3.2078E-01 -1.7962E-01  4.4183E-01  2.0179E-01  2.3469E-02  9.4677E-01  3.8810E-01  6.9131E-02
             6.5191E-01
 GRADIENT:   3.5405E+01  5.3302E+01 -9.8816E+00  4.4328E+01  1.5780E+01  2.3144E+01  1.8014E+01 -1.3572E+00  8.4521E+00  3.4638E+00
             1.5572E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3275.69266888944        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0060E+00  1.5608E+00  1.2474E+00  7.1054E-01  1.4077E+00  1.0370E+00  8.0860E-01  2.3323E+00  1.3307E+00  9.3105E-01
             1.7368E+00
 PARAMETER:  1.0600E-01  5.4521E-01  3.2105E-01 -2.4174E-01  4.4193E-01  1.3634E-01 -1.1245E-01  9.4686E-01  3.8569E-01  2.8563E-02
             6.5203E-01
 GRADIENT:  -9.6144E-01  1.6464E+01 -4.2582E+00 -2.8474E+00  7.2809E-01 -2.3407E+00  1.3573E+00 -2.6314E+00 -1.9554E+00 -7.5207E+00
             1.3814E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3275.89125409838        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      504             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0048E+00  1.5608E+00  1.2474E+00  7.0755E-01  1.4077E+00  1.0449E+00  7.9022E-01  2.3323E+00  1.3848E+00  9.6060E-01
             1.7368E+00
 PARAMETER:  1.0479E-01  5.4520E-01  3.2105E-01 -2.4594E-01  4.4193E-01  1.4392E-01 -1.3545E-01  9.4686E-01  4.2558E-01  5.9803E-02
             6.5202E-01
 GRADIENT:  -3.1637E+00  1.1334E+01 -3.8306E+00 -2.6450E+00 -2.2429E+00  8.2167E-01  2.3669E+00 -2.4089E+00  3.2581E+00 -3.0445E+00
             1.6719E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3275.98086298444        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:      632
 NPARAMETR:  1.0142E+00  1.5613E+00  1.2476E+00  7.0746E-01  1.4080E+00  1.0493E+00  7.9028E-01  2.3311E+00  1.3852E+00  9.6055E-01
             1.7374E+00
 PARAMETER:  1.1415E-01  5.4550E-01  3.2123E-01 -2.4607E-01  4.4217E-01  1.4808E-01 -1.3537E-01  9.4635E-01  4.2581E-01  5.9749E-02
             6.5237E-01
 GRADIENT:   1.6996E+01  1.1461E+01 -3.7583E+00 -2.5786E+00 -2.2263E+00  2.6708E+00  2.3797E+00 -2.4409E+00  3.3070E+00 -3.0657E+00
             1.7113E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3276.21731105357        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      813            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0141E+00  1.5738E+00  1.2855E+00  7.1327E-01  1.4156E+00  1.0505E+00  7.7081E-01  2.3302E+00  1.3671E+00  9.8229E-01
             1.7376E+00
 PARAMETER:  1.1403E-01  5.5349E-01  3.5112E-01 -2.3789E-01  4.4754E-01  1.4928E-01 -1.6031E-01  9.4594E-01  4.1268E-01  8.2127E-02
             6.5252E-01
 GRADIENT:   1.6510E+01  3.2564E+01 -2.8236E-01  7.8189E+00 -5.5655E+00  3.1697E+00 -1.2834E-01 -5.0207E+00  3.6149E-02 -7.5538E-01
             1.6326E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3276.36360376291        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      969
 NPARAMETR:  1.0154E+00  1.5901E+00  1.3135E+00  7.0153E-01  1.4450E+00  1.0519E+00  7.6891E-01  2.3303E+00  1.3905E+00  1.0005E+00
             1.7376E+00
 PARAMETER:  1.1524E-01  5.6380E-01  3.7272E-01 -2.5449E-01  4.6808E-01  1.5058E-01 -1.6278E-01  9.4598E-01  4.2969E-01  1.0048E-01
             6.5248E-01
 GRADIENT:   2.1671E+00 -3.0636E-01  7.5238E-01  3.9098E+00  1.4048E+00  4.8902E-01  5.8848E-01 -7.2666E+00  2.6862E-01  1.8517E-01
             1.3908E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3276.45173186384        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1144
 NPARAMETR:  1.0143E+00  1.6290E+00  1.2738E+00  6.7309E-01  1.4664E+00  1.0506E+00  7.5057E-01  2.3307E+00  1.4327E+00  1.0183E+00
             1.7373E+00
 PARAMETER:  1.1418E-01  5.8799E-01  3.4203E-01 -2.9588E-01  4.8282E-01  1.4938E-01 -1.8692E-01  9.4616E-01  4.5956E-01  1.1812E-01
             6.5234E-01
 GRADIENT:   8.9701E-02  1.7033E-01 -4.1120E-02  8.6374E-02  1.0957E-03  2.0456E-02 -1.5998E-02 -6.3164E+00  2.7207E-03 -2.4812E-02
             1.0978E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3276.54580228470        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1326
 NPARAMETR:  1.0140E+00  1.6216E+00  1.2851E+00  6.7696E-01  1.4631E+00  1.0506E+00  7.5003E-01  2.3590E+00  1.4306E+00  1.0183E+00
             1.7205E+00
 PARAMETER:  1.1388E-01  5.8344E-01  3.5087E-01 -2.9015E-01  4.8055E-01  1.4938E-01 -1.8764E-01  9.5823E-01  4.5812E-01  1.1814E-01
             6.4259E-01
 GRADIENT:  -1.9735E-01 -3.3536E-01  8.5985E-02 -6.3666E-01 -9.9667E-02 -2.9231E-02  6.6152E-02 -6.8216E+00  4.9356E-02  2.0216E-02
            -8.2207E+00

0ITERATION NO.:   53    OBJECTIVE VALUE:  -3276.54668962809        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1422
 NPARAMETR:  1.0140E+00  1.6193E+00  1.2863E+00  6.7879E-01  1.4615E+00  1.0507E+00  7.5074E-01  2.3596E+00  1.4277E+00  1.0172E+00
             1.7201E+00
 PARAMETER:  1.1393E-01  5.8202E-01  3.5173E-01 -2.8744E-01  4.7948E-01  1.4943E-01 -1.8670E-01  9.5848E-01  4.5606E-01  1.1703E-01
             6.4239E-01
 GRADIENT:  -9.1312E-02 -1.6252E+03 -2.6780E+03  3.2646E+03 -1.9721E+03 -1.3906E-02 -5.0146E+03  9.9106E+02 -1.8790E+03  8.0214E+03
            -7.5208E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1422
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2327E-03 -4.2931E-02 -2.6477E-02  3.0396E-02 -4.0503E-02
 SE:             2.9782E-02  2.0983E-02  1.8774E-02  2.4898E-02  2.3516E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6698E-01  4.0763E-02  1.5845E-01  2.2217E-01  8.5002E-02

 ETASHRINKSD(%)  2.2635E-01  2.9703E+01  3.7104E+01  1.6587E+01  2.1220E+01
 ETASHRINKVR(%)  4.5219E-01  5.0583E+01  6.0440E+01  3.0423E+01  3.7937E+01
 EBVSHRINKSD(%)  6.4699E-01  2.9885E+01  4.4109E+01  1.8486E+01  1.9358E+01
 EBVSHRINKVR(%)  1.2898E+00  5.0838E+01  6.8761E+01  3.3554E+01  3.4969E+01
 RELATIVEINF(%)  9.8700E+01  1.1395E+01  2.1633E+01  1.8157E+01  2.4094E+01
 EPSSHRINKSD(%)  1.9003E+01
 EPSSHRINKVR(%)  3.4395E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3276.5466896280850     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1622.4573298596742     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3276.547       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.62E+00  1.29E+00  6.79E-01  1.46E+00  1.05E+00  7.51E-01  2.36E+00  1.43E+00  1.02E+00  1.72E+00
 


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
+        9.68E+02
 
 TH 2
+       -3.93E+02  2.68E+05
 
 TH 3
+       -1.13E+03  3.71E+03  1.15E+06
 
 TH 4
+        5.35E+03 -3.39E+03 -1.10E+04  6.18E+06
 
 TH 5
+       -5.11E+02  3.60E+05  4.70E+03 -4.78E+03  4.84E+05
 
 TH 6
+       -2.16E+00 -4.37E+02 -1.26E+03  6.01E+03 -5.72E+02  1.78E+02
 
 TH 7
+       -1.51E+04  2.95E+03  8.59E+03 -5.16E+04  3.88E+03 -1.75E+04  1.19E+07
 
 TH 8
+        7.91E+01 -7.55E+02 -4.98E+02  7.53E+02 -8.96E+02  8.86E+01 -6.01E+02  4.72E+04
 
 TH 9
+       -3.10E+06  2.52E+02  8.50E+02 -4.88E+03  3.38E+02 -2.29E+06  2.57E+06 -2.59E+01  5.52E+05
 
 TH10
+        8.18E+03 -6.58E+03 -1.95E+04  1.01E+07 -8.63E+03  9.18E+03 -7.12E+04  1.31E+03 -7.30E+03  1.66E+07
 
 TH11
+       -1.79E+02  1.61E+03  1.11E+03 -1.64E+03  1.94E+03 -1.90E+02  1.31E+03 -2.64E+03  7.30E+01 -2.83E+03  1.98E+05
 
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
 #CPUT: Total CPU Time in Seconds,       58.166
Stop Time:
Sat Sep 25 00:23:20 CDT 2021
