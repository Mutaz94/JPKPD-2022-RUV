Wed Sep 29 21:55:37 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1768.06596378602        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0554E+02  4.8369E+01  2.4161E+01  1.1982E+02  3.4757E+01  6.1818E+01  7.2661E+00 -4.8321E+01  2.1427E+01  1.9352E+01
            -6.9910E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1928.31519849597        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0464E+00  1.0624E+00  1.0507E+00  9.6603E-01  1.0142E+00  9.8856E-01  9.3940E-01  9.7506E-01  8.7565E-01  7.9590E-01
             1.8359E+00
 PARAMETER:  1.4539E-01  1.6058E-01  1.4947E-01  6.5442E-02  1.1415E-01  8.8495E-02  3.7483E-02  7.4748E-02 -3.2788E-02 -1.2828E-01
             7.0755E-01
 GRADIENT:   2.3827E+02  3.9435E+01  2.0022E+01  2.2619E+01 -3.2178E+01  3.4014E+01 -2.6511E-01 -3.1910E+00 -5.8563E+00  7.3206E+00
             6.8854E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1933.71781196228        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0522E+00  1.1408E+00  4.6721E-01  8.7331E-01  7.0379E-01  9.3245E-01  1.0210E+00  9.8331E-01  9.0836E-01  3.8945E-01
             1.7067E+00
 PARAMETER:  1.5084E-01  2.3171E-01 -6.6097E-01 -3.5462E-02 -2.5127E-01  3.0061E-02  1.2080E-01  8.3174E-02  3.8850E-03 -8.4302E-01
             6.3457E-01
 GRADIENT:   2.7141E+02  9.3877E+01  3.5808E+01  3.2717E+01 -7.3983E+01  8.2954E+00  1.1589E+01  8.6570E-03  9.6838E+00  6.3483E+00
             4.7265E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1935.82590339988        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0153E+00  1.1282E+00  4.3321E-01  8.6084E-01  6.8643E-01  9.2722E-01  9.9918E-01  1.0655E+00  8.7919E-01  3.7274E-01
             1.6327E+00
 PARAMETER:  1.1523E-01  2.2065E-01 -7.3653E-01 -4.9847E-02 -2.7626E-01  2.4430E-02  9.9177E-02  1.6340E-01 -2.8750E-02 -8.8687E-01
             5.9022E-01
 GRADIENT:   1.4376E+02  7.0052E+01  2.2423E+01  2.4034E+01 -4.8851E+01  7.6575E+00  9.9291E+00  1.6822E+00  1.3188E+00  6.5735E+00
             2.1165E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1956.92792362847        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      410
 NPARAMETR:  1.0418E+00  1.4600E+00  5.2306E-01  7.0215E-01  9.4826E-01  9.4340E-01  6.7442E-01  2.3747E+00  1.1470E+00  2.5087E-01
             1.5239E+00
 PARAMETER:  1.4093E-01  4.7845E-01 -5.4805E-01 -2.5361E-01  4.6871E-02  4.1734E-02 -2.9390E-01  9.6489E-01  2.3717E-01 -1.2828E+00
             5.2125E-01
 GRADIENT:   2.7288E+01  8.5024E-01 -4.4861E+00  4.6088E+01 -2.6165E+00  6.4896E-01 -8.4574E+00 -1.3233E+01 -6.7312E+00 -1.1371E+00
             2.2122E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1957.96965648377        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      576
 NPARAMETR:  1.0417E+00  1.4684E+00  5.2657E-01  6.9571E-01  9.5652E-01  9.3982E-01  7.3279E-01  2.4096E+00  1.1618E+00  2.5041E-01
             1.5228E+00
 PARAMETER:  1.4086E-01  4.8418E-01 -5.4136E-01 -2.6282E-01  5.5542E-02  3.7928E-02 -2.1089E-01  9.7946E-01  2.4994E-01 -1.2847E+00
             5.2058E-01
 GRADIENT:   2.7669E+01 -6.0715E+00 -4.3663E+00  4.4329E+01 -6.0016E-01 -8.8353E-01  7.4740E-01 -1.0994E+01  6.9326E-01 -1.0085E+00
             2.2948E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1958.24260134029        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      757             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0310E+00  1.4667E+00  5.2740E-01  6.6514E-01  9.5672E-01  9.4207E-01  7.2890E-01  2.4065E+00  1.1791E+00  2.5014E-01
             1.5204E+00
 PARAMETER:  1.3056E-01  4.8300E-01 -5.3979E-01 -3.0775E-01  5.5760E-02  4.0327E-02 -2.1622E-01  9.7818E-01  2.6473E-01 -1.2857E+00
             5.1898E-01
 GRADIENT:   2.5262E+02  1.2818E+02  8.1646E+00  5.2103E+01 -2.0570E-01  1.5939E+01  4.5835E+00  5.0307E+00  6.7222E+00 -3.2604E-01
             2.6443E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1958.92403646554        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      945            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0313E+00  1.4667E+00  5.2709E-01  6.6431E-01  9.5678E-01  9.4113E-01  7.3480E-01  2.4417E+00  1.1647E+00  2.5035E-01
             1.4758E+00
 PARAMETER:  1.3082E-01  4.8303E-01 -5.4039E-01 -3.0901E-01  5.5821E-02  3.9322E-02 -2.0815E-01  9.9270E-01  2.5243E-01 -1.2849E+00
             4.8921E-01
 GRADIENT:   2.6862E+02  1.3666E+02  1.0417E+01  5.6750E+01  1.1470E+01  1.5922E+01  5.0140E+00  3.0190E+00  2.2415E+00 -1.5410E+00
             1.2485E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1959.05395442539        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1122
 NPARAMETR:  1.0295E+00  1.4650E+00  5.2633E-01  6.6485E-01  9.4959E-01  9.4394E-01  6.9729E-01  2.4330E+00  1.2103E+00  2.6539E-01
             1.4684E+00
 PARAMETER:  1.2911E-01  4.8187E-01 -5.4183E-01 -3.0819E-01  4.8271E-02  4.2303E-02 -2.6056E-01  9.8914E-01  2.9090E-01 -1.2266E+00
             4.8418E-01
 GRADIENT:   1.4169E-01 -4.8073E+01  8.5915E+00  2.5836E+00 -2.1622E+00  5.5524E-01 -1.2406E+00 -1.8800E+01 -3.1799E-01 -2.2365E+00
             9.1175E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1960.16836967997        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1302
 NPARAMETR:  1.0332E+00  1.4661E+00  5.2598E-01  6.5965E-01  9.5936E-01  9.4127E-01  6.0314E-01  2.4354E+00  1.3459E+00  5.2482E-01
             1.4685E+00
 PARAMETER:  1.3265E-01  4.8263E-01 -5.4249E-01 -3.1605E-01  5.8514E-02  3.9470E-02 -4.0561E-01  9.9013E-01  3.9709E-01 -5.4469E-01
             4.8426E-01
 GRADIENT:   9.4887E+00 -6.7901E+01  4.3897E+00  4.0544E-01 -9.1805E+00 -6.6218E-01 -6.2817E-01 -2.3361E+01  8.1332E+00  2.6727E-01
             3.9671E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1960.54725865000        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1479
 NPARAMETR:  1.0307E+00  1.4666E+00  5.2584E-01  6.5997E-01  9.6378E-01  9.4304E-01  6.7189E-01  2.4359E+00  1.2364E+00  5.0819E-01
             1.4685E+00
 PARAMETER:  1.3024E-01  4.8292E-01 -5.4277E-01 -3.1556E-01  6.3108E-02  4.1354E-02 -2.9766E-01  9.9033E-01  3.1217E-01 -5.7690E-01
             4.8421E-01
 GRADIENT:   3.2684E+00 -6.9756E+01  3.4364E+00  9.2994E-02  5.8078E-01  1.5399E-01  1.6808E-01 -1.9821E+01 -8.4547E-01  2.9439E-01
             3.5576E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1964.99408029521        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1660
 NPARAMETR:  1.0306E+00  1.5626E+00  5.0438E-01  5.9522E-01  9.9428E-01  9.9032E-01  6.2579E-01  2.5771E+00  1.4938E+00  6.2816E-01
             1.4437E+00
 PARAMETER:  1.3010E-01  5.4638E-01 -5.8442E-01 -4.1883E-01  9.4266E-02  9.0276E-02 -3.6874E-01  1.0467E+00  5.0135E-01 -3.6496E-01
             4.6719E-01
 GRADIENT:   3.2544E+00 -6.0525E+01  9.4183E+00 -2.3833E+00 -4.1537E+01  1.7849E+01  6.3062E+00 -2.0426E+01  1.7518E+01  4.8043E+00
             4.0468E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1980.26991930118        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1839
 NPARAMETR:  1.0367E+00  1.9951E+00  4.2941E-01  3.8565E-01  1.2221E+00  1.0246E+00  6.0404E-01  3.2001E+00  2.0505E+00  8.4672E-01
             1.3522E+00
 PARAMETER:  1.3607E-01  7.9072E-01 -7.4535E-01 -8.5282E-01  3.0059E-01  1.2434E-01 -4.0412E-01  1.2632E+00  8.1810E-01 -6.6382E-02
             4.0173E-01
 GRADIENT:   1.5288E+01  1.0945E+02  6.4878E+00  3.7320E+01 -4.5971E+01  2.8824E+01  3.3888E-01 -1.4552E+01  1.0585E+01  4.7068E+00
            -7.0857E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1981.31654428866        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     1936
 NPARAMETR:  1.0359E+00  1.9877E+00  4.3078E-01  3.8656E-01  1.2206E+00  9.0306E-01  6.0512E-01  3.2185E+00  2.0572E+00  8.4710E-01
             1.3498E+00
 PARAMETER:  1.3528E-01  7.8700E-01 -7.4215E-01 -8.5047E-01  2.9937E-01 -1.9717E-03 -4.0233E-01  1.2689E+00  8.2135E-01 -6.5934E-02
             3.9996E-01
 GRADIENT:   3.4926E+02  7.2754E+02  1.0296E+01  9.5831E+01 -3.5834E+01  6.3744E-01  1.2743E+01  2.0867E+01  3.0908E+01  5.2778E+00
            -4.4788E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1981.95458867974        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:     2060             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0351E+00  1.9812E+00  4.3204E-01  3.8739E-01  1.2193E+00  9.3780E-01  6.0609E-01  3.2351E+00  2.0632E+00  8.4745E-01
             1.3476E+00
 PARAMETER:  1.3455E-01  7.8372E-01 -7.3923E-01 -8.4832E-01  2.9825E-01  3.5786E-02 -4.0072E-01  1.2741E+00  8.2427E-01 -6.5525E-02
             3.9836E-01
 GRADIENT:   3.4581E+02  7.1502E+02  1.0003E+01  9.4799E+01 -3.5117E+01  1.6269E+01  1.3846E+01  2.2068E+01  3.2304E+01  5.5795E+00
            -4.8466E+00

0ITERATION NO.:   74    OBJECTIVE VALUE:  -1982.51991980063        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     2174
 NPARAMETR:  1.0345E+00  1.9670E+00  4.3288E-01  3.8763E-01  1.2184E+00  9.3732E-01  6.0675E-01  3.2443E+00  2.0668E+00  8.4769E-01
             1.3462E+00
 PARAMETER:  1.3393E-01  7.7637E-01 -7.3749E-01 -8.4790E-01  2.9761E-01  3.4269E-02 -3.9975E-01  1.2765E+00  8.2578E-01 -6.5271E-02
             3.9737E-01
 GRADIENT:   7.9582E+02 -2.0386E+02 -2.7883E+02 -2.2058E+02  6.6335E+02 -3.8403E+00 -5.1976E+02 -9.9746E+01 -2.4115E+02 -2.0953E+03
             5.2173E+02
 NUMSIGDIG:         2.3         2.5         2.3         2.4         2.3         0.7         2.3         2.2         2.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2174
 NO. OF SIG. DIGITS IN FINAL EST.:  0.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.5140E-03 -7.1966E-02 -4.4822E-02  1.8596E-02 -4.5548E-02
 SE:             3.0099E-02  2.1058E-02  1.8971E-02  2.1674E-02  1.7999E-02
 N:                     100         100         100         100         100

 P VAL.:         9.0706E-01  6.3191E-04  1.8142E-02  3.9091E-01  1.1386E-02

 ETASHRINKSD(%)  1.0000E-10  2.9454E+01  3.6446E+01  2.7389E+01  3.9701E+01
 ETASHRINKVR(%)  1.0000E-10  5.0233E+01  5.9609E+01  4.7276E+01  6.3641E+01
 EBVSHRINKSD(%)  6.4408E-01  2.6149E+01  4.4959E+01  2.4277E+01  3.5812E+01
 EBVSHRINKVR(%)  1.2840E+00  4.5460E+01  6.9705E+01  4.2660E+01  5.8799E+01
 RELATIVEINF(%)  9.8677E+01  1.1761E+01  1.4953E+01  1.3563E+01  1.8097E+01
 EPSSHRINKSD(%)  3.4372E+01
 EPSSHRINKVR(%)  5.6930E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1982.5199198006267     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1063.5813865959540     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.45
 Elapsed covariance  time in seconds:     8.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1982.520       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.97E+00  4.33E-01  3.88E-01  1.22E+00  9.36E-01  6.07E-01  3.24E+00  2.07E+00  8.48E-01  1.35E+00
 


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
 
         3.30E-03  2.42E-02  1.69E-01  7.49E-03  8.42E-03  6.61E-02  1.17E-01  9.27E-02  3.88E-02  6.68E-02  1.24E-02
 


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
+        1.09E-05
 
 TH 2
+       -7.64E-05  5.88E-04
 
 TH 3
+       -3.70E-04  2.34E-03  2.86E-02
 
 TH 4
+       -2.46E-05  1.76E-04  8.60E-04  5.60E-05
 
 TH 5
+        2.76E-05 -1.96E-04 -8.79E-04 -6.25E-05  7.08E-05
 
 TH 6
+       -8.00E-05  6.28E-04 -4.79E-04  1.71E-04 -2.12E-04  4.37E-03
 
 TH 7
+       -2.58E-04  1.93E-03  5.18E-03  5.77E-04 -6.56E-04  2.99E-03  1.38E-02
 
 TH 8
+       -3.04E-04  2.15E-03  9.73E-03  6.87E-04 -7.79E-04  2.33E-03  7.23E-03  8.59E-03
 
 TH 9
+       -1.28E-04  8.97E-04  4.40E-03  2.89E-04 -3.25E-04  9.11E-04  2.98E-03  3.57E-03  1.50E-03
 
 TH10
+        1.82E-04 -1.25E-03 -9.21E-03 -4.16E-04  4.47E-04 -8.70E-04 -6.05E-03 -4.94E-03 -2.14E-03  4.46E-03
 
 TH11
+        4.07E-05 -2.87E-04 -1.36E-03 -9.22E-05  1.04E-04 -3.03E-04 -9.61E-04 -1.14E-03 -4.78E-04  6.73E-04  1.53E-04
 
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
+        3.30E-03
 
 TH 2
+       -9.56E-01  2.42E-02
 
 TH 3
+       -6.64E-01  5.71E-01  1.69E-01
 
 TH 4
+       -9.96E-01  9.68E-01  6.79E-01  7.49E-03
 
 TH 5
+        9.96E-01 -9.61E-01 -6.18E-01 -9.92E-01  8.42E-03
 
 TH 6
+       -3.67E-01  3.92E-01 -4.28E-02  3.46E-01 -3.81E-01  6.61E-02
 
 TH 7
+       -6.66E-01  6.78E-01  2.61E-01  6.57E-01 -6.64E-01  3.85E-01  1.17E-01
 
 TH 8
+       -9.94E-01  9.56E-01  6.21E-01  9.90E-01 -9.98E-01  3.81E-01  6.64E-01  9.27E-02
 
 TH 9
+       -9.98E-01  9.54E-01  6.71E-01  9.98E-01 -9.95E-01  3.55E-01  6.54E-01  9.94E-01  3.88E-02
 
 TH10
+        8.27E-01 -7.73E-01 -8.15E-01 -8.32E-01  7.95E-01 -1.97E-01 -7.71E-01 -7.98E-01 -8.25E-01  6.68E-02
 
 TH11
+        9.97E-01 -9.57E-01 -6.47E-01 -9.94E-01  9.95E-01 -3.70E-01 -6.60E-01 -9.94E-01 -9.96E-01  8.13E-01  1.24E-02
 
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
+        1.06E+08
 
 TH 2
+       -2.18E+06  3.44E+05
 
 TH 3
+        4.10E+07  2.73E+05  4.41E+07
 
 TH 4
+        2.43E+07 -4.34E+06 -7.97E+06  5.97E+07
 
 TH 5
+       -9.60E+06  1.58E+06  7.39E+06 -2.00E+07  1.59E+07
 
 TH 6
+        8.35E+04 -5.69E+03  2.00E+04  7.26E+04 -2.59E+04  4.39E+02
 
 TH 7
+        5.39E+07  3.54E+05  5.79E+07 -1.04E+07  9.71E+06  2.62E+04  7.61E+07
 
 TH 8
+       -1.93E+05 -8.05E+03 -2.15E+05  2.30E+05  3.00E+05 -7.26E+01 -2.82E+05  3.63E+04
 
 TH 9
+       -5.40E+06  7.81E+05 -1.41E+06 -1.04E+07  3.67E+06 -1.37E+04 -1.86E+06 -2.61E+04  2.26E+06
 
 TH10
+        1.54E+08  1.02E+06  1.65E+08 -2.98E+07  2.78E+07  7.48E+04  2.17E+08 -8.06E+05 -5.29E+06  6.20E+08
 
 TH11
+       -3.29E+06  7.13E+04  1.33E+06 -7.14E+05  7.39E+05 -1.18E+03  1.75E+06  2.18E+04  9.61E+04  5.01E+06  1.33E+06
 
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
 #CPUT: Total CPU Time in Seconds,       44.501
Stop Time:
Wed Sep 29 21:56:23 CDT 2021
