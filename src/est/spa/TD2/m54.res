Sat Sep 25 13:36:54 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat54.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1714.91392877482        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.4226E+01  2.1190E+01 -2.7534E+01  5.5270E+01 -8.1160E+00 -3.5437E+01 -4.2763E+00  1.3037E+01 -5.9372E+00  1.2191E+01
             1.2725E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1722.34291237293        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:       91
 NPARAMETR:  1.0292E+00  1.0417E+00  1.1262E+00  9.8153E-01  1.0732E+00  1.0411E+00  1.0279E+00  9.1449E-01  1.0444E+00  9.6817E-01
             9.8106E-01
 PARAMETER:  1.2883E-01  1.4084E-01  2.1882E-01  8.1361E-02  1.7066E-01  1.4027E-01  1.2755E-01  1.0611E-02  1.4347E-01  6.7657E-02
             8.0881E-02
 GRADIENT:   1.4482E+01  3.4543E+01 -3.9737E+00  5.1422E+01  1.6022E-02 -8.2622E+00 -8.6947E-02  4.4998E+00 -5.3302E-01 -6.2497E+00
            -3.7937E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1723.45021066957        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0279E+00  1.0090E+00  1.1310E+00  9.9528E-01  1.0779E+00  1.0451E+00  8.3348E-01  6.2959E-01  1.1182E+00  1.1072E+00
             9.6954E-01
 PARAMETER:  1.2757E-01  1.0897E-01  2.2307E-01  9.5264E-02  1.7499E-01  1.4408E-01 -8.2150E-02 -3.6269E-01  2.1173E-01  2.0187E-01
             6.9069E-02
 GRADIENT:   1.4075E+01  2.3826E+01 -7.4616E+00  4.7666E+01  1.9721E+00 -5.9657E+00  6.2860E-01 -9.7909E-05  6.1145E+00  1.5343E+00
            -1.8410E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1726.16937239237        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      238
 NPARAMETR:  1.0189E+00  1.1298E+00  1.4564E+00  8.6930E-01  1.2353E+00  1.0629E+00  3.1540E-01  1.1467E+00  1.2964E+00  1.2390E+00
             9.8380E-01
 PARAMETER:  1.1868E-01  2.2205E-01  4.7595E-01 -4.0065E-02  3.1133E-01  1.6096E-01 -1.0539E+00  2.3693E-01  3.5960E-01  3.1430E-01
             8.3663E-02
 GRADIENT:  -5.5183E+00  1.1773E+00  4.4596E+00 -1.5176E+01 -4.4775E+00  2.7747E+00  1.9346E+00  9.3075E-02 -1.8381E+00  9.0967E-01
             1.1101E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1726.56497675332        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      312
 NPARAMETR:  1.0219E+00  1.0008E+00  1.4256E+00  9.6131E-01  1.1811E+00  1.0544E+00  1.4631E-01  1.0244E+00  1.2262E+00  1.2185E+00
             9.7852E-01
 PARAMETER:  1.2170E-01  1.0080E-01  4.5462E-01  6.0538E-02  2.6647E-01  1.5294E-01 -1.8220E+00  1.2409E-01  3.0390E-01  2.9759E-01
             7.8288E-02
 GRADIENT:   2.6135E+00 -3.5132E+00  1.6214E-01 -9.3309E-01 -1.3662E+00 -8.8924E-01  4.4515E-01  3.8517E-02  3.3849E+00  5.7470E-01
            -1.6786E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1727.89125709527        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      433
 NPARAMETR:  1.0400E+00  1.0407E+00  1.4625E+00  9.3960E-01  1.2129E+00  1.0788E+00  4.1445E-02  1.0819E+00  1.2602E+00  1.2332E+00
             9.9067E-01
 PARAMETER:  1.3922E-01  1.3985E-01  4.8018E-01  3.7698E-02  2.9302E-01  1.7581E-01 -3.0834E+00  1.7871E-01  3.3126E-01  3.0964E-01
             9.0624E-02
 GRADIENT:  -3.0493E+01 -5.2068E+00 -1.4914E-01 -3.1309E+00  1.7547E+00 -4.3339E+00  1.5047E-02 -1.7723E-01 -1.2048E+00  1.7162E-01
             3.8615E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1728.44135018471        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      608
 NPARAMETR:  1.0577E+00  1.1784E+00  1.4147E+00  8.4608E-01  1.2518E+00  1.0897E+00  3.2831E-02  1.1456E+00  1.4070E+00  1.2522E+00
             9.8563E-01
 PARAMETER:  1.5609E-01  2.6415E-01  4.4693E-01 -6.7145E-02  3.2462E-01  1.8589E-01 -3.3164E+00  2.3597E-01  4.4148E-01  3.2492E-01
             8.5529E-02
 GRADIENT:  -1.8406E-01  1.1247E-01  2.0020E-01 -2.4882E-02 -3.1941E-01 -2.1746E-01 -3.4620E-03 -2.7994E-02  9.1334E-03 -3.6525E-02
             3.1963E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1728.45078002029        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      789
 NPARAMETR:  1.0581E+00  1.1916E+00  1.3994E+00  8.3806E-01  1.2546E+00  1.0913E+00  9.2551E-02  1.1405E+00  1.4185E+00  1.2551E+00
             9.8512E-01
 PARAMETER:  1.5643E-01  2.7531E-01  4.3607E-01 -7.6664E-02  3.2678E-01  1.8732E-01 -2.2800E+00  2.3143E-01  4.4959E-01  3.2723E-01
             8.5007E-02
 GRADIENT:   2.5729E-01  2.5078E-01 -5.3916E-01  1.0390E+00  4.9351E-01  3.0252E-01 -1.6422E-02  4.3927E-02  2.1296E-01  1.9239E-01
            -3.3114E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1728.46941693500        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      965
 NPARAMETR:  1.0580E+00  1.2407E+00  1.3878E+00  8.0479E-01  1.2694E+00  1.0895E+00  1.7048E-01  1.1666E+00  1.4629E+00  1.2585E+00
             9.8670E-01
 PARAMETER:  1.5637E-01  3.1570E-01  4.2771E-01 -1.1717E-01  3.3853E-01  1.8571E-01 -1.6691E+00  2.5407E-01  4.8044E-01  3.2994E-01
             8.6608E-02
 GRADIENT:  -3.4666E-01  5.0449E-01  4.0959E-01  3.9090E-01 -6.4846E-01 -4.1415E-01  5.1497E-03 -5.9247E-02  7.7766E-02 -6.0266E-02
             2.5859E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1728.47001878794        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1132
 NPARAMETR:  1.0581E+00  1.2416E+00  1.3818E+00  8.0402E-01  1.2695E+00  1.0905E+00  1.7148E-01  1.1667E+00  1.4634E+00  1.2585E+00
             9.8676E-01
 PARAMETER:  1.5647E-01  3.1637E-01  4.2336E-01 -1.1813E-01  3.3860E-01  1.8667E-01 -1.6633E+00  2.5420E-01  4.8077E-01  3.2989E-01
             8.6667E-02
 GRADIENT:  -2.0092E-01 -2.5181E-01 -3.2627E-01  1.9484E-01  3.3047E-01 -4.1134E-02  3.1534E-03  1.1469E-01 -1.5641E-02 -3.4865E-02
             1.8984E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1728.47036511011        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1310
 NPARAMETR:  1.0582E+00  1.2415E+00  1.3838E+00  8.0397E-01  1.2696E+00  1.0907E+00  1.7076E-01  1.1654E+00  1.4636E+00  1.2582E+00
             9.8653E-01
 PARAMETER:  1.5654E-01  3.1635E-01  4.2482E-01 -1.1819E-01  3.3870E-01  1.8677E-01 -1.6675E+00  2.5307E-01  4.8091E-01  3.2972E-01
             8.6434E-02
 GRADIENT:  -4.2709E-02 -7.8006E-02 -1.6420E-02  1.2295E-01  7.1642E-02  2.7301E-04  1.2817E-03  1.2555E-02  2.0185E-03 -1.2169E-01
             1.7944E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1728.47042415715        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     1401
 NPARAMETR:  1.0582E+00  1.2416E+00  1.3839E+00  8.0387E-01  1.2695E+00  1.0907E+00  1.7063E-01  1.1648E+00  1.4636E+00  1.2591E+00
             9.8648E-01
 PARAMETER:  1.5660E-01  3.1639E-01  4.2491E-01 -1.1832E-01  3.3863E-01  1.8681E-01 -1.6683E+00  2.5253E-01  4.8091E-01  3.3039E-01
             8.6387E-02
 GRADIENT:   6.2551E-02 -8.2386E-02  4.6249E-02  3.6610E-02 -1.2157E-01  1.6261E-02  1.1469E-03  1.9828E-03 -1.6087E-02  1.0407E-03
             1.8736E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1401
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0321E-04 -1.7998E-02 -2.1000E-02 -4.6516E-04 -3.1498E-02
 SE:             2.9892E-02  4.9902E-03  1.0625E-02  2.8901E-02  2.3771E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8657E-01  3.1014E-04  4.8097E-02  9.8716E-01  1.8515E-01

 ETASHRINKSD(%)  1.0000E-10  8.3282E+01  6.4405E+01  3.1794E+00  2.0365E+01
 ETASHRINKVR(%)  1.0000E-10  9.7205E+01  8.7330E+01  6.2577E+00  3.6583E+01
 EBVSHRINKSD(%)  3.5250E-01  8.5883E+01  6.7452E+01  3.2283E+00  1.7054E+01
 EBVSHRINKVR(%)  7.0376E-01  9.8007E+01  8.9406E+01  6.3524E+00  3.1200E+01
 RELATIVEINF(%)  9.9192E+01  2.7929E-01  5.0079E+00  1.4802E+01  3.1078E+01
 EPSSHRINKSD(%)  4.3026E+01
 EPSSHRINKVR(%)  6.7540E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1728.4704241571510     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -993.31959759341282     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.81
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1728.470       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.24E+00  1.38E+00  8.04E-01  1.27E+00  1.09E+00  1.71E-01  1.16E+00  1.46E+00  1.26E+00  9.86E-01
 


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
+        2.31E+09
 
 TH 2
+       -9.74E+08  4.11E+08
 
 TH 3
+        6.51E+08  2.38E+04  4.99E+01
 
 TH 4
+       -4.02E+09  1.70E+09  9.82E+04  7.53E+02
 
 TH 5
+       -8.90E+08  3.90E+03  2.51E+08  1.55E+09  3.43E+08
 
 TH 6
+       -8.96E+03  3.78E+03  4.39E-01  1.56E+04 -7.34E-02  1.53E+09
 
 TH 7
+        1.34E+09  2.85E+03 -3.79E+08  1.19E+04 -5.57E+03 -1.09E+09  7.83E+08
 
 TH 8
+        2.05E+05 -8.64E+04 -1.34E+01 -3.57E+05  5.01E+08  5.04E+03 -7.57E+08  7.33E+08
 
 TH 9
+        4.96E+04 -2.10E+04  1.04E+00 -8.63E+04  1.92E+00  2.11E+03 -3.16E+08 -6.07E-01  1.28E+08
 
 TH10
+        1.01E+00  3.88E+08  2.24E+04 -1.92E+00 -4.03E+01  3.57E+03  2.28E+00 -8.16E+04 -1.97E+04  5.88E+01
 
 TH11
+        3.88E+09 -5.12E+05  1.09E+09 -6.76E+09 -1.50E+09  3.16E+09 -1.15E+04 -2.19E+09  9.13E+08 -4.83E+05  6.52E+09
 
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
 #CPUT: Total CPU Time in Seconds,       23.062
Stop Time:
Sat Sep 25 13:37:21 CDT 2021
