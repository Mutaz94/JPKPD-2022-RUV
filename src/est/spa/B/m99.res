Wed Sep 29 11:41:21 CDT 2021
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
$DATA ../../../../data/spa/B/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1683.00893202755        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1279E+02 -2.9926E+01 -2.4908E+01  1.8648E+01  5.7152E+01  7.5542E+01 -9.9298E+00  2.1846E+00  7.7235E+00  7.2952E+00
             8.6934E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1687.24392075388        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  1.0017E+00  1.0377E+00  1.0085E+00  1.0197E+00  9.7205E-01  9.1654E-01  1.0188E+00  9.9785E-01  9.9367E-01  9.8062E-01
             9.8948E-01
 PARAMETER:  1.0172E-01  1.3699E-01  1.0848E-01  1.1953E-01  7.1649E-02  1.2849E-02  1.1862E-01  9.7848E-02  9.3646E-02  8.0426E-02
             8.9427E-02
 GRADIENT:  -8.2889E+00 -3.2212E+00  3.8650E-01  4.0385E+00 -4.1768E+00 -1.4830E+00 -1.0007E+01 -1.7600E-01  1.7778E+00  8.5092E+00
             1.8919E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1688.63026083482        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  1.0035E+00  1.0467E+00  8.4145E-01  9.9529E-01  8.6659E-01  9.3369E-01  1.2900E+00  1.0744E+00  8.9348E-01  6.5198E-01
             9.5194E-01
 PARAMETER:  1.0346E-01  1.4565E-01 -7.2634E-02  9.5278E-02 -4.3185E-02  3.1394E-02  3.5465E-01  1.7172E-01 -1.2629E-02 -3.2774E-01
             5.0750E-02
 GRADIENT:  -6.0291E+00  1.2356E+01  1.8499E+01 -1.2867E+01 -2.4537E+01  4.9760E+00  5.0825E+00 -1.5941E+00 -3.2978E+00 -4.1445E+00
            -1.4114E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1689.90150870197        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      489
 NPARAMETR:  1.0056E+00  1.1740E+00  6.7538E-01  9.0730E-01  8.5846E-01  9.2429E-01  1.1341E+00  7.1544E-01  9.5616E-01  7.1696E-01
             9.8360E-01
 PARAMETER:  1.0559E-01  2.6042E-01 -2.9248E-01  2.7226E-03 -5.2611E-02  2.1271E-02  2.2582E-01 -2.3486E-01  5.5165E-02 -2.3273E-01
             8.3466E-02
 GRADIENT:  -5.3915E+00  1.2359E+00  1.1016E+00  2.0064E+00 -2.0981E+00  4.5774E-01 -8.5454E-01 -1.6909E-01 -7.4098E-01 -4.9757E-01
             5.3849E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1690.04004635839        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  1.0093E+00  1.3806E+00  5.5776E-01  7.7398E-01  8.9686E-01  9.2321E-01  1.0085E+00  5.7850E-01  1.0628E+00  7.3183E-01
             9.8197E-01
 PARAMETER:  1.0923E-01  4.2254E-01 -4.8383E-01 -1.5621E-01 -8.8579E-03  2.0102E-02  1.0851E-01 -4.4732E-01  1.6087E-01 -2.1221E-01
             8.1804E-02
 GRADIENT:   2.6498E+00  5.1421E+00  9.2347E-01  3.0369E+00 -2.9224E+00 -3.9703E-01  1.6249E-01  5.8477E-02  9.3376E-02  1.9590E-01
            -4.7366E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1690.05584605029        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.0082E+00  1.3950E+00  5.4494E-01  7.6083E-01  9.0014E-01  9.2406E-01  1.0000E+00  5.2415E-01  1.0737E+00  7.3366E-01
             9.8339E-01
 PARAMETER:  1.0817E-01  4.3289E-01 -5.0707E-01 -1.7334E-01 -5.2065E-03  2.1016E-02  1.0002E-01 -5.4598E-01  1.7108E-01 -2.0971E-01
             8.3248E-02
 GRADIENT:  -1.5390E-01  9.8299E-02  2.2779E-01  3.5411E-01 -1.1541E-01 -5.1897E-02  1.0692E-01  2.5880E-02  4.8618E-02 -9.0098E-03
             4.6281E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1690.05696192220        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  1.0083E+00  1.3971E+00  5.4207E-01  7.5927E-01  9.0036E-01  9.2420E-01  9.9894E-01  4.8352E-01  1.0756E+00  7.3810E-01
             9.8333E-01
 PARAMETER:  1.0823E-01  4.3438E-01 -5.1236E-01 -1.7540E-01 -4.9598E-03  2.1175E-02  9.8939E-02 -6.2666E-01  1.7288E-01 -2.0367E-01
             8.3187E-02
 GRADIENT:  -2.7580E-02 -1.6267E-01  1.0192E-01  5.8049E-01  3.7488E-01 -5.7294E-03  4.2496E-02 -7.3410E-03  3.2276E-02  5.1908E-02
            -4.2098E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1690.05699854125        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  1.0084E+00  1.4009E+00  5.3824E-01  7.5677E-01  9.0021E-01  9.2412E-01  9.9747E-01  4.6202E-01  1.0779E+00  7.3894E-01
             9.8326E-01
 PARAMETER:  1.0838E-01  4.3713E-01 -5.1944E-01 -1.7870E-01 -5.1277E-03  2.1089E-02  9.7467E-02 -6.7215E-01  1.7504E-01 -2.0254E-01
             8.3121E-02
 GRADIENT:   3.1592E-01  1.5460E-01  4.7438E-02  9.1242E-01  3.3556E-01 -5.1268E-02  1.1264E-01 -1.1661E-02  8.3241E-02  9.0875E-02
            -3.9331E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1690.05853449820        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1386
 NPARAMETR:  1.0091E+00  1.4018E+00  5.3604E-01  7.5445E-01  8.9999E-01  9.2476E-01  9.9571E-01  4.6191E-01  1.0787E+00  7.3829E-01
             9.8337E-01
 PARAMETER:  1.0904E-01  4.3777E-01 -5.2354E-01 -1.8176E-01 -5.3673E-03  2.1778E-02  9.5705E-02 -6.7238E-01  1.7576E-01 -2.0342E-01
             8.3229E-02
 GRADIENT:   2.0827E+00 -1.9922E+00  2.9928E-02 -8.1269E-01  3.8623E-01  2.2562E-01 -1.0025E-02  6.4906E-03  3.5892E-02  2.0309E-01
             5.1845E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1690.05898662402        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1566
 NPARAMETR:  1.0090E+00  1.4020E+00  5.3569E-01  7.5464E-01  8.9964E-01  9.2465E-01  9.9600E-01  4.7234E-01  1.0787E+00  7.3660E-01
             9.8344E-01
 PARAMETER:  1.0899E-01  4.3791E-01 -5.2420E-01 -1.8152E-01 -5.7618E-03  2.1660E-02  9.5995E-02 -6.5006E-01  1.7575E-01 -2.0571E-01
             8.3304E-02
 GRADIENT:   1.9376E+00 -1.6911E+00 -2.1766E-01 -2.7227E-01  5.7578E-01  1.8006E-01  9.5557E-02  2.5472E-02  1.4888E-01  2.1166E-01
             1.0917E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1690.05941585689        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1744
 NPARAMETR:  1.0090E+00  1.4021E+00  5.3555E-01  7.5473E-01  8.9940E-01  9.2464E-01  9.9592E-01  4.6988E-01  1.0784E+00  7.3561E-01
             9.8336E-01
 PARAMETER:  1.0899E-01  4.3798E-01 -5.2447E-01 -1.8140E-01 -6.0244E-03  2.1650E-02  9.5910E-02 -6.5527E-01  1.7544E-01 -2.0706E-01
             8.3216E-02
 GRADIENT:   1.9275E+00 -1.2648E+00 -3.0836E-02 -1.6143E-01  4.4059E-01  1.7427E-01  1.8999E-02  8.2849E-03  6.2765E-02  6.5798E-02
             1.7820E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1690.05957330430        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1930             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0090E+00  1.4020E+00  5.3539E-01  7.5476E-01  8.9898E-01  9.2464E-01  9.9594E-01  4.6613E-01  1.0779E+00  7.3483E-01
             9.8331E-01
 PARAMETER:  1.0899E-01  4.3790E-01 -5.2476E-01 -1.8136E-01 -6.4978E-03  2.1650E-02  9.5932E-02 -6.6329E-01  1.7503E-01 -2.0812E-01
             8.3172E-02
 GRADIENT:   4.7342E+02  2.9398E+02  9.4859E+00  7.7331E+01  9.5996E+00  4.5845E+01  5.6513E+00  2.0966E-01  7.7600E+00  9.6210E-01
             7.5383E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1690.05968104369        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2114
 NPARAMETR:  1.0090E+00  1.4019E+00  5.3522E-01  7.5481E-01  8.9878E-01  9.2464E-01  9.9602E-01  4.6551E-01  1.0778E+00  7.3459E-01
             9.8331E-01
 PARAMETER:  1.0899E-01  4.3780E-01 -5.2508E-01 -1.8129E-01 -6.7208E-03  2.1651E-02  9.6011E-02 -6.6463E-01  1.7489E-01 -2.0845E-01
             8.3169E-02
 GRADIENT:   1.9334E+00 -9.2885E-01  3.1431E-01 -4.1987E-01 -6.7155E-02  1.7368E-01 -5.4648E-02 -9.7925E-03 -3.8772E-02 -4.4306E-02
            -4.7684E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1690.05982069042        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2288
 NPARAMETR:  1.0090E+00  1.4015E+00  5.3476E-01  7.5496E-01  8.9861E-01  9.2464E-01  9.9632E-01  4.6900E-01  1.0779E+00  7.3466E-01
             9.8339E-01
 PARAMETER:  1.0899E-01  4.3771E-01 -5.2541E-01 -1.8122E-01 -6.8698E-03  2.1652E-02  9.6110E-02 -6.6378E-01  1.7483E-01 -2.0858E-01
             8.3183E-02
 GRADIENT:   1.7853E-03  2.3618E-01  2.3380E-01 -1.6027E-01  6.4961E-02 -3.3923E-04 -4.2038E-02 -4.2831E-03 -2.7853E-02 -2.3743E-02
            -2.8810E-02
 NUMSIGDIG:         5.2         3.5         3.0         3.1         3.5         5.1         2.7         1.9         3.0         2.9
                    3.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2288
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6866E-04 -1.3384E-02 -1.5354E-02  1.1056E-02 -2.4871E-02
 SE:             2.9835E-02  2.4725E-02  6.5873E-03  2.3662E-02  2.0125E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9549E-01  5.8830E-01  1.9761E-02  6.4031E-01  2.1654E-01

 ETASHRINKSD(%)  5.0092E-02  1.7167E+01  7.7932E+01  2.0730E+01  3.2577E+01
 ETASHRINKVR(%)  1.0016E-01  3.1387E+01  9.5130E+01  3.7162E+01  5.4542E+01
 EBVSHRINKSD(%)  4.8198E-01  1.7087E+01  8.0232E+01  2.1082E+01  3.2248E+01
 EBVSHRINKVR(%)  9.6164E-01  3.1255E+01  9.6092E+01  3.7719E+01  5.4096E+01
 RELATIVEINF(%)  9.9010E+01  4.7565E+00  3.0455E-01  4.0511E+00  5.4420E+00
 EPSSHRINKSD(%)  4.4896E+01
 EPSSHRINKVR(%)  6.9636E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1690.0598206904181     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -954.90899412667989     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.01
 Elapsed covariance  time in seconds:     5.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1690.060       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.40E+00  5.35E-01  7.55E-01  8.99E-01  9.25E-01  9.96E-01  4.66E-01  1.08E+00  7.34E-01  9.83E-01
 


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
 
         2.82E-02  4.65E-01  3.13E-01  3.15E-01  1.37E-01  7.66E-02  2.01E-01  1.87E+00  2.81E-01  2.40E-01  7.84E-02
 


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
+        7.93E-04
 
 TH 2
+       -1.22E-03  2.16E-01
 
 TH 3
+       -1.57E-04 -1.32E-01  9.77E-02
 
 TH 4
+        6.94E-04 -1.45E-01  9.14E-02  9.90E-02
 
 TH 5
+       -8.09E-04  5.08E-02 -2.14E-02 -3.29E-02  1.88E-02
 
 TH 6
+        2.42E-04 -1.55E-03  3.83E-04  1.14E-03 -2.60E-04  5.86E-03
 
 TH 7
+        4.26E-04 -8.65E-02  5.40E-02  5.78E-02 -1.94E-02  2.52E-05  4.03E-02
 
 TH 8
+       -1.00E-03 -6.61E-01  5.01E-01  4.56E-01 -1.20E-01 -1.75E-02  2.78E-01  3.51E+00
 
 TH 9
+       -9.53E-04  1.20E-01 -6.94E-02 -8.09E-02  3.22E-02  1.49E-04 -4.67E-02 -3.97E-01  7.88E-02
 
 TH10
+       -9.94E-04  7.69E-02 -4.57E-02 -5.09E-02  2.24E-02  2.47E-03 -3.12E-02 -3.64E-01  5.04E-02  5.76E-02
 
 TH11
+       -4.08E-05  3.70E-03  1.49E-04 -2.10E-03  2.63E-03  1.49E-04 -1.72E-03 -1.06E-02  2.80E-03  1.66E-03  6.15E-03
 
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
+        2.82E-02
 
 TH 2
+       -9.35E-02  4.65E-01
 
 TH 3
+       -1.79E-02 -9.07E-01  3.13E-01
 
 TH 4
+        7.83E-02 -9.91E-01  9.29E-01  3.15E-01
 
 TH 5
+       -2.10E-01  7.98E-01 -5.01E-01 -7.64E-01  1.37E-01
 
 TH 6
+        1.12E-01 -4.34E-02  1.60E-02  4.72E-02 -2.48E-02  7.66E-02
 
 TH 7
+        7.54E-02 -9.28E-01  8.61E-01  9.16E-01 -7.05E-01  1.64E-03  2.01E-01
 
 TH 8
+       -1.90E-02 -7.58E-01  8.55E-01  7.74E-01 -4.67E-01 -1.22E-01  7.39E-01  1.87E+00
 
 TH 9
+       -1.21E-01  9.20E-01 -7.91E-01 -9.16E-01  8.38E-01  6.92E-03 -8.29E-01 -7.55E-01  2.81E-01
 
 TH10
+       -1.47E-01  6.89E-01 -6.09E-01 -6.74E-01  6.80E-01  1.35E-01 -6.47E-01 -8.10E-01  7.48E-01  2.40E-01
 
 TH11
+       -1.85E-02  1.01E-01  6.10E-03 -8.50E-02  2.45E-01  2.48E-02 -1.09E-01 -7.23E-02  1.27E-01  8.81E-02  7.84E-02
 
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
+        1.41E+03
 
 TH 2
+        2.03E+01  4.76E+02
 
 TH 3
+        6.57E+01  2.58E+02  8.27E+02
 
 TH 4
+       -4.27E+01  2.84E+02 -6.03E+02  1.16E+03
 
 TH 5
+       -5.67E+01 -3.74E+02 -8.08E+02  5.19E+02  1.28E+03
 
 TH 6
+       -5.81E+01  9.04E+00  6.30E+00 -2.70E+01 -1.37E+01  1.88E+02
 
 TH 7
+       -1.10E+01  1.23E+02  3.91E+01  7.27E+00 -4.69E+01  1.63E+01  1.96E+02
 
 TH 8
+        7.01E+00 -1.13E+01 -3.38E+01  1.57E+01  6.85E+00  1.79E+00 -4.04E+00  4.19E+00
 
 TH 9
+        1.08E+01 -2.73E+01 -8.17E+01  1.19E+02 -4.77E+01 -7.74E+00 -2.81E+01  9.55E+00  1.48E+02
 
 TH10
+        6.52E+01 -1.52E+01 -2.76E+01 -3.52E+01 -1.18E+02 -7.76E+00 -4.45E+00  1.53E+01  1.10E+01  1.16E+02
 
 TH11
+       -6.82E+00 -4.51E+00 -3.82E+01  1.78E+00 -7.34E+01 -6.08E-01  9.08E+00  7.69E+00  1.78E+01  3.80E+01  1.96E+02
 
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
 #CPUT: Total CPU Time in Seconds,       34.596
Stop Time:
Wed Sep 29 11:41:57 CDT 2021
