Thu Sep 30 01:41:11 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat84.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2145.85036777148        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6683E+02 -8.7201E+00 -5.1075E+01  6.5260E+01  3.2643E+01  8.3409E+01 -6.6318E+00  1.4892E+01  1.6474E+00  1.2823E+01
            -3.2271E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2157.46665586545        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0241E+00  1.1241E+00  1.2644E+00  9.6400E-01  1.1271E+00  7.7396E-01  1.0618E+00  8.9578E-01  1.0543E+00  9.4241E-01
             1.0540E+00
 PARAMETER:  1.2386E-01  2.1699E-01  3.3457E-01  6.3340E-02  2.1962E-01 -1.5624E-01  1.5996E-01 -1.0062E-02  1.5290E-01  4.0685E-02
             1.5259E-01
 GRADIENT:  -4.3709E+01  1.6001E+01  7.0374E+00  1.0507E+01  1.7216E-01 -4.0326E+01  2.8635E+00 -2.5950E-01 -3.0400E+00 -2.1516E+01
             3.2811E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2159.50928151915        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0268E+00  1.2372E+00  1.4100E+00  8.8360E-01  1.2364E+00  8.0315E-01  9.1021E-01  8.2999E-01  1.1836E+00  1.1156E+00
             1.0577E+00
 PARAMETER:  1.2646E-01  3.1287E-01  4.4360E-01 -2.3752E-02  3.1222E-01 -1.1921E-01  5.9171E-03 -8.6337E-02  2.6859E-01  2.0943E-01
             1.5607E-01
 GRADIENT:  -3.2154E+01  7.6933E+00  1.9417E+01 -4.7551E+00 -2.5720E+00 -2.3057E+01  2.4513E+00 -6.9015E+00  5.1603E-01 -1.0539E+01
             7.5518E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2161.74236799975        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0373E+00  1.3201E+00  1.2375E+00  8.3115E-01  1.2264E+00  8.5043E-01  8.5644E-01  9.1105E-01  1.2291E+00  1.1467E+00
             1.0419E+00
 PARAMETER:  1.3660E-01  3.7768E-01  3.1311E-01 -8.4947E-02  3.0408E-01 -6.2017E-02 -5.4969E-02  6.8476E-03  3.0629E-01  2.3687E-01
             1.4104E-01
 GRADIENT:   8.9548E-01  4.1166E+00  7.6661E-01  4.5913E+00 -1.3304E+00  3.4715E-01 -6.2329E-01 -3.9917E-02 -2.9360E-01  2.5669E-02
            -1.5751E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2161.80588477671        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0385E+00  1.4706E+00  1.0551E+00  7.3441E-01  1.2411E+00  8.5030E-01  8.3163E-01  7.9292E-01  1.3269E+00  1.1327E+00
             1.0419E+00
 PARAMETER:  1.3779E-01  4.8566E-01  1.5364E-01 -2.0869E-01  3.1601E-01 -6.2164E-02 -8.4365E-02 -1.3203E-01  3.8285E-01  2.2458E-01
             1.4109E-01
 GRADIENT:   1.1677E+00  7.7813E+00 -1.8417E+00  9.4898E+00  2.1470E+00 -3.3743E-01 -4.4334E-01  6.1763E-01 -3.0172E-02 -1.6231E-01
            -1.5409E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2161.90995866220        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      905             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0396E+00  1.5838E+00  9.3523E-01  6.4549E-01  1.2562E+00  8.5182E-01  8.1121E-01  6.1542E-01  1.4226E+00  1.1296E+00
             1.0421E+00
 PARAMETER:  1.3880E-01  5.5985E-01  3.3036E-02 -3.3775E-01  3.2805E-01 -6.0383E-02 -1.0922E-01 -3.8545E-01  4.5252E-01  2.2184E-01
             1.4125E-01
 GRADIENT:   5.8381E+02  4.9430E+02  3.9285E+00  9.8711E+01  1.2670E+01  2.6474E+01  6.2721E+00 -9.3120E-02  2.3955E+01  1.1511E+00
             9.9657E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2161.94324381385        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  1.0387E+00  1.5911E+00  9.1914E-01  6.4624E-01  1.2564E+00  8.5187E-01  8.0906E-01  6.0010E-01  1.4311E+00  1.1311E+00
             1.0423E+00
 PARAMETER:  1.3793E-01  5.6440E-01  1.5682E-02 -3.3658E-01  3.2822E-01 -6.0322E-02 -1.1188E-01 -4.1067E-01  4.5843E-01  2.2322E-01
             1.4143E-01
 GRADIENT:  -5.9825E-01  1.9983E-01  6.7964E-01  1.6082E+00 -2.1464E-01 -6.0634E-03  8.4070E-03  3.2945E-02 -1.5683E-01  2.0561E-01
            -2.0859E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2161.94832253677        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  1.0380E+00  1.6246E+00  8.5944E-01  6.2633E-01  1.2488E+00  8.5233E-01  8.0617E-01  4.7790E-01  1.4594E+00  1.1219E+00
             1.0425E+00
 PARAMETER:  1.3727E-01  5.8525E-01 -5.1473E-02 -3.6788E-01  3.2218E-01 -5.9776E-02 -1.1546E-01 -6.3836E-01  4.7806E-01  2.1503E-01
             1.4163E-01
 GRADIENT:  -3.5261E+00  4.3273E+00  2.4054E-01  4.8612E+00 -1.9626E+00  9.4719E-03 -3.0446E-03  3.1306E-02  6.4010E-01  8.9236E-01
             2.0387E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2161.99452253575        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1445             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0401E+00  1.6359E+00  8.3241E-01  6.1251E-01  1.2486E+00  8.5254E-01  8.0521E-01  4.0994E-01  1.4694E+00  1.1123E+00
             1.0423E+00
 PARAMETER:  1.3935E-01  5.9218E-01 -8.3429E-02 -3.9019E-01  3.2200E-01 -5.9540E-02 -1.1665E-01 -7.9174E-01  4.8487E-01  2.0642E-01
             1.4146E-01
 GRADIENT:   5.8561E+02  5.5488E+02  4.0418E-01  1.0830E+02  1.7013E+01  2.6487E+01  6.6487E+00  9.6766E-02  2.6188E+01  1.6066E+00
             1.4715E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2161.99571275635        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1624             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0403E+00  1.6367E+00  8.3108E-01  6.1199E-01  1.2479E+00  8.5254E-01  8.0582E-01  3.9767E-01  1.4686E+00  1.1118E+00
             1.0423E+00
 PARAMETER:  1.3952E-01  5.9269E-01 -8.5025E-02 -3.9105E-01  3.2143E-01 -5.9535E-02 -1.1590E-01 -8.2214E-01  4.8432E-01  2.0600E-01
             1.4140E-01
 GRADIENT:   5.8687E+02  5.5634E+02  9.0033E-01  1.0825E+02  1.6363E+01  2.6477E+01  6.6548E+00  5.7571E-02  2.5969E+01  1.5812E+00
             1.3681E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2161.99630281952        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1804
 NPARAMETR:  1.0403E+00  1.6369E+00  8.2965E-01  6.1169E-01  1.2474E+00  8.5255E-01  8.0604E-01  3.9113E-01  1.4684E+00  1.1113E+00
             1.0423E+00
 PARAMETER:  1.3952E-01  5.9280E-01 -8.6753E-02 -3.9154E-01  3.2103E-01 -5.9521E-02 -1.1562E-01 -8.3870E-01  4.8416E-01  2.0552E-01
             1.4140E-01
 GRADIENT:   3.0952E+00 -4.6213E+00  4.2453E-01 -1.4697E-01 -1.0140E-01  4.9122E-02 -1.4288E-02 -1.0827E-02  5.2046E-02  2.8394E-02
            -5.0417E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2161.99713187147        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1987             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0403E+00  1.6372E+00  8.2679E-01  6.1144E-01  1.2469E+00  8.5256E-01  8.0630E-01  3.8942E-01  1.4686E+00  1.1107E+00
             1.0423E+00
 PARAMETER:  1.3953E-01  5.9298E-01 -9.0210E-02 -3.9194E-01  3.2069E-01 -5.9508E-02 -1.1530E-01 -8.4310E-01  4.8432E-01  2.0495E-01
             1.4142E-01
 GRADIENT:   5.8684E+02  5.5625E+02  4.4080E-01  1.0844E+02  1.6901E+01  2.6472E+01  6.6442E+00  8.3239E-02  2.6014E+01  1.6419E+00
             1.4624E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2161.99724375349        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2167             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0403E+00  1.6376E+00  8.2390E-01  6.1077E-01  1.2467E+00  8.5258E-01  8.0656E-01  3.8459E-01  1.4691E+00  1.1101E+00
             1.0423E+00
 PARAMETER:  1.3954E-01  5.9324E-01 -9.3703E-02 -3.9303E-01  3.2050E-01 -5.9487E-02 -1.1497E-01 -8.5558E-01  4.8466E-01  2.0446E-01
             1.4145E-01
 GRADIENT:   5.8684E+02  5.5602E+02  6.3249E-02  1.0848E+02  1.7457E+01  2.6470E+01  6.6536E+00  1.0249E-01  2.6083E+01  1.6873E+00
             1.5542E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2161.99885810145        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2350             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0403E+00  1.6390E+00  8.2300E-01  6.1038E-01  1.2461E+00  8.5259E-01  8.0650E-01  3.7611E-01  1.4694E+00  1.1094E+00
             1.0423E+00
 PARAMETER:  1.3955E-01  5.9406E-01 -9.4794E-02 -3.9367E-01  3.2002E-01 -5.9477E-02 -1.1505E-01 -8.7787E-01  4.8484E-01  2.0386E-01
             1.4140E-01
 GRADIENT:   5.8691E+02  5.5861E+02  5.1993E-01  1.0871E+02  1.6715E+01  2.6473E+01  6.6564E+00  7.4823E-02  2.5977E+01  1.6078E+00
             1.4333E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2161.99933440690        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2530
 NPARAMETR:  1.0403E+00  1.6391E+00  8.2165E-01  6.1002E-01  1.2458E+00  8.5260E-01  8.0657E-01  3.7198E-01  1.4694E+00  1.1090E+00
             1.0423E+00
 PARAMETER:  1.3955E-01  5.9415E-01 -9.6445E-02 -3.9427E-01  3.1977E-01 -5.9465E-02 -1.1496E-01 -8.8893E-01  4.8486E-01  2.0347E-01
             1.4141E-01
 GRADIENT:   3.0832E+00 -5.2385E+00 -7.4734E-02  8.1490E-02  4.0465E-01  4.5337E-02  4.1157E-03  1.8180E-02  1.6871E-01  8.1802E-02
             4.3536E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2161.99977396860        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2710
 NPARAMETR:  1.0401E+00  1.6405E+00  8.2045E-01  6.0973E-01  1.2453E+00  8.5259E-01  8.0668E-01  3.6508E-01  1.4692E+00  1.1085E+00
             1.0423E+00
 PARAMETER:  1.3933E-01  5.9499E-01 -9.7900E-02 -3.9473E-01  3.1940E-01 -5.9480E-02 -1.1483E-01 -9.0764E-01  4.8469E-01  2.0296E-01
             1.4140E-01
 GRADIENT:   2.3681E+00 -3.9457E+00  2.0629E-01  4.5448E-01 -6.9439E-02  3.3607E-02 -7.7915E-03  3.6533E-03  5.4277E-02  2.6442E-02
            -1.9362E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2162.00073018248        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2893             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0404E+00  1.6405E+00  8.1802E-01  6.0908E-01  1.2450E+00  8.5262E-01  8.0676E-01  3.6032E-01  1.4700E+00  1.1079E+00
             1.0423E+00
 PARAMETER:  1.3957E-01  5.9500E-01 -1.0087E-01 -3.9581E-01  3.1916E-01 -5.9436E-02 -1.1473E-01 -9.2076E-01  4.8529E-01  2.0250E-01
             1.4140E-01
 GRADIENT:   5.8697E+02  5.6011E+02  4.1646E-01  1.0884E+02  1.6795E+01  2.6473E+01  6.6532E+00  7.6090E-02  2.5945E+01  1.6048E+00
             1.4606E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2162.00119987756        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3073
 NPARAMETR:  1.0404E+00  1.6408E+00  8.1694E-01  6.0877E-01  1.2447E+00  8.5263E-01  8.0684E-01  3.5591E-01  1.4702E+00  1.1075E+00
             1.0423E+00
 PARAMETER:  1.3958E-01  5.9520E-01 -1.0219E-01 -3.9631E-01  3.1893E-01 -5.9426E-02 -1.1462E-01 -9.3308E-01  4.8541E-01  2.0215E-01
             1.4141E-01
 GRADIENT:   3.0838E+00 -5.3011E+00 -1.0106E-01  1.1170E-01  3.7626E-01  4.5537E-02  1.4341E-02  1.7298E-02  1.7567E-01  7.6957E-02
             5.0120E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2162.00171406426        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3255             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0404E+00  1.6417E+00  8.1616E-01  6.0822E-01  1.2441E+00  8.5265E-01  8.0689E-01  3.4380E-01  1.4700E+00  1.1068E+00
             1.0422E+00
 PARAMETER:  1.3959E-01  5.9572E-01 -1.0314E-01 -3.9723E-01  3.1838E-01 -5.9407E-02 -1.1457E-01 -9.6771E-01  4.8523E-01  2.0145E-01
             1.4135E-01
 GRADIENT:   5.8712E+02  5.6204E+02  1.0560E+00  1.0872E+02  1.5958E+01  2.6482E+01  6.6399E+00  3.9215E-02  2.5745E+01  1.4970E+00
             1.3321E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2162.00229965429        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3435
 NPARAMETR:  1.0404E+00  1.6419E+00  8.1493E-01  6.0805E-01  1.2440E+00  8.5265E-01  8.0695E-01  3.4304E-01  1.4703E+00  1.1066E+00
             1.0422E+00
 PARAMETER:  1.3959E-01  5.9584E-01 -1.0465E-01 -3.9750E-01  3.1832E-01 -5.9403E-02 -1.1450E-01 -9.6991E-01  4.8543E-01  2.0130E-01
             1.4137E-01
 GRADIENT:   3.0948E+00 -4.9137E+00  2.5204E-01 -3.5686E-02 -9.0522E-02  4.8781E-02 -9.3167E-05 -2.3189E-03  6.2310E-02  2.1514E-02
            -3.0478E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2162.00292779248        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3617             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0404E+00  1.6424E+00  8.1268E-01  6.0760E-01  1.2438E+00  8.5266E-01  8.0706E-01  3.3993E-01  1.4709E+00  1.1062E+00
             1.0423E+00
 PARAMETER:  1.3960E-01  5.9618E-01 -1.0741E-01 -3.9823E-01  3.1814E-01 -5.9390E-02 -1.1435E-01 -9.7901E-01  4.8587E-01  2.0091E-01
             1.4140E-01
 GRADIENT:   5.8705E+02  5.6236E+02  5.0649E-01  1.0901E+02  1.6645E+01  2.6475E+01  6.6585E+00  6.7391E-02  2.5894E+01  1.5716E+00
             1.4571E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2162.00328954781        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3795
 NPARAMETR:  1.0404E+00  1.6429E+00  8.1147E-01  6.0730E-01  1.2436E+00  8.5267E-01  8.0709E-01  3.3673E-01  1.4711E+00  1.1059E+00
             1.0423E+00
 PARAMETER:  1.3960E-01  5.9644E-01 -1.0890E-01 -3.9873E-01  3.1798E-01 -5.9382E-02 -1.1433E-01 -9.8849E-01  4.8599E-01  2.0063E-01
             1.4139E-01
 GRADIENT:   3.0838E+00 -5.3774E+00 -1.4361E-01  1.4826E-01  3.5492E-01  4.5357E-02  1.2850E-02  1.7350E-02  1.6941E-01  7.3109E-02
             5.0513E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -2162.00375987120        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3977             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0404E+00  1.6438E+00  8.1090E-01  6.0669E-01  1.2430E+00  8.5269E-01  8.0710E-01  3.2361E-01  1.4709E+00  1.1052E+00
             1.0422E+00
 PARAMETER:  1.3961E-01  5.9700E-01 -1.0961E-01 -3.9973E-01  3.1750E-01 -5.9364E-02 -1.1431E-01 -1.0282E+00  4.8590E-01  2.0003E-01
             1.4133E-01
 GRADIENT:   5.8721E+02  5.6454E+02  1.1181E+00  1.0894E+02  1.5865E+01  2.6484E+01  6.6439E+00  3.5342E-02  2.5697E+01  1.4733E+00
             1.3233E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -2162.00416754616        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     4143
 NPARAMETR:  1.0402E+00  1.6447E+00  8.0959E-01  6.0664E-01  1.2431E+00  8.5266E-01  8.0714E-01  3.2606E-01  1.4712E+00  1.1052E+00
             1.0422E+00
 PARAMETER:  1.3939E-01  5.9756E-01 -1.1123E-01 -3.9981E-01  3.1760E-01 -5.9390E-02 -1.1426E-01 -1.0207E+00  4.8610E-01  2.0002E-01
             1.4137E-01
 GRADIENT:  -6.8013E-01  6.5247E-01  1.0829E-01  6.3856E-01 -6.9071E-02 -1.5208E-02 -1.1363E-02  2.6883E-03 -7.1471E-02  8.2104E-03
            -1.5818E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4143
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2254E-04 -3.2394E-02 -8.4473E-03  2.4277E-02 -3.6459E-02
 SE:             2.9850E-02  2.1406E-02  3.5368E-03  2.3627E-02  2.3202E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8336E-01  1.3020E-01  1.6923E-02  3.0417E-01  1.1609E-01

 ETASHRINKSD(%)  1.0000E-10  2.8287E+01  8.8151E+01  2.0847E+01  2.2271E+01
 ETASHRINKVR(%)  1.0000E-10  4.8572E+01  9.8596E+01  3.7348E+01  3.9583E+01
 EBVSHRINKSD(%)  4.6987E-01  2.6677E+01  9.0014E+01  2.3164E+01  1.9676E+01
 EBVSHRINKVR(%)  9.3753E-01  4.6238E+01  9.9003E+01  4.0962E+01  3.5481E+01
 RELATIVEINF(%)  9.8891E+01  4.1430E+00  2.2032E-01  4.8825E+00  1.8111E+01
 EPSSHRINKSD(%)  3.2872E+01
 EPSSHRINKVR(%)  5.4938E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2162.0041675461571     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1243.0656343414844     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    69.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2162.004       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.64E+00  8.10E-01  6.07E-01  1.24E+00  8.53E-01  8.07E-01  3.26E-01  1.47E+00  1.11E+00  1.04E+00
 


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
+        1.40E+03
 
 TH 2
+       -7.20E+00  3.61E+02
 
 TH 3
+        5.85E+00  9.14E+01  1.81E+02
 
 TH 4
+       -6.87E+00  3.55E+02 -1.37E+02  8.35E+02
 
 TH 5
+        1.52E+00 -1.33E+02 -1.62E+02  1.53E+02  3.80E+02
 
 TH 6
+        1.42E+00 -1.03E+00  1.89E+00 -1.93E+00  3.28E-01  2.69E+02
 
 TH 7
+        1.89E+00  3.36E+00  9.92E+00 -1.78E+01 -1.80E+01 -5.39E-01  1.00E+02
 
 TH 8
+       -1.82E-01 -5.28E+00 -1.77E+01  2.97E+00  1.01E+00 -7.92E-02  2.08E+00  3.82E+00
 
 TH 9
+        1.62E-01 -1.82E+01 -2.02E+01  4.70E+01  4.38E+00 -4.15E-01  2.39E+01  1.54E+00  3.94E+01
 
 TH10
+        1.54E+00 -1.27E+01 -2.27E+01  4.41E+00 -4.73E+01  8.39E-01  2.26E+00  1.53E+00  5.30E+00  7.35E+01
 
 TH11
+       -1.04E+01 -1.74E+01 -2.02E+01 -3.68E+00  2.08E+00  2.44E+00  9.42E+00  2.21E+00  1.90E+00  1.81E+01  3.81E+02
 
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
 #CPUT: Total CPU Time in Seconds,       76.710
Stop Time:
Thu Sep 30 01:42:30 CDT 2021
