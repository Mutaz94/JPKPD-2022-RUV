Sat Sep 18 04:02:11 CDT 2021
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
$DATA ../../../../data/int/SL2/dat98.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2141.38951118682        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0159E+02 -4.8715E+00  3.6416E+01  9.7089E+01  1.3095E+02  2.0586E+01 -7.1733E+01 -1.6288E+02 -9.7029E+01 -4.1432E+01
            -3.0901E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3042.18160253245        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9486E-01  1.2900E+00  1.3837E+00  8.3646E-01  1.2143E+00  9.7395E-01  1.0231E+00  8.4110E-01  1.2000E+00  1.1063E+00
             1.8665E+00
 PARAMETER:  9.4845E-02  3.5461E-01  4.2473E-01 -7.8579E-02  2.9414E-01  7.3603E-02  1.2279E-01 -7.3041E-02  2.8232E-01  2.0098E-01
             7.2408E-01
 GRADIENT:   4.7036E+01  1.0663E+01  6.1403E+00  1.0850E+01 -8.7603E+00  1.0964E+01  6.3476E+00 -6.7393E+00 -4.1665E+00 -1.6001E+01
            -2.0698E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3048.36401172228        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.8677E-01  1.4291E+00  1.2186E+00  7.6272E-01  1.2462E+00  9.6146E-01  7.6390E-01  4.8257E-01  1.2606E+00  1.3790E+00
             1.9360E+00
 PARAMETER:  8.6679E-02  4.5703E-01  2.9766E-01 -1.7086E-01  3.2013E-01  6.0698E-02 -1.6932E-01 -6.2864E-01  3.3161E-01  4.2133E-01
             7.6062E-01
 GRADIENT:   2.5331E+01  3.3569E+01 -3.8094E+00  2.1708E+01 -3.0207E+01  6.6760E+00 -7.1945E+00 -1.5380E+00 -4.9347E+00  1.5135E+01
            -1.3312E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3054.68158098587        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.7771E-01  1.4752E+00  1.3443E+00  7.1679E-01  1.3528E+00  9.4555E-01  8.2509E-01  4.2518E-01  1.3074E+00  1.3216E+00
             2.0353E+00
 PARAMETER:  7.7457E-02  4.8878E-01  3.9588E-01 -2.3298E-01  4.0221E-01  4.4013E-02 -9.2258E-02 -7.5524E-01  3.6802E-01  3.7888E-01
             8.1066E-01
 GRADIENT:   4.4148E-02  1.3544E+00  4.3457E-01  8.5411E-02 -1.9395E-02  6.6083E-01  8.2775E-01 -6.3239E-01 -8.5576E-02  1.1273E+00
            -6.5284E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3062.10274397501        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:      312
 NPARAMETR:  9.7951E-01  1.4578E+00  2.1946E+00  7.5086E-01  1.4853E+00  9.3245E-01  7.9822E-01  3.1797E+00  1.2389E+00  1.3461E+00
             2.0237E+00
 PARAMETER:  7.9302E-02  4.7696E-01  8.8601E-01 -1.8653E-01  4.9563E-01  3.0062E-02 -1.2538E-01  1.2568E+00  3.1424E-01  3.9722E-01
             8.0493E-01
 GRADIENT:   4.2055E+00  5.5872E+00 -8.3637E+00  5.6448E-01  1.2243E+01 -4.8668E+00  1.3512E-01  1.1077E+01  5.4364E+00  2.6330E-01
             4.1327E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3067.79343774601        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  9.7613E-01  1.1677E+00  3.4573E+00  9.5654E-01  1.4020E+00  9.5017E-01  1.0174E+00  3.1573E+00  9.5863E-01  1.2555E+00
             1.9592E+00
 PARAMETER:  7.5841E-02  2.5507E-01  1.3405E+00  5.5569E-02  4.3789E-01  4.8888E-02  1.1720E-01  1.2497E+00  5.7751E-02  3.2753E-01
             7.7253E-01
 GRADIENT:  -8.8667E-01  1.3055E+01  3.6800E+00 -6.6566E-01 -1.1736E+01  2.5666E+00  6.7637E-01 -2.2634E+00 -1.2652E+00  8.4441E-01
            -9.6740E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3067.93085646006        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      512             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8112E-01  1.1528E+00  3.4729E+00  9.6528E-01  1.4004E+00  9.4693E-01  1.0225E+00  3.1751E+00  9.6517E-01  1.2473E+00
             1.9609E+00
 PARAMETER:  8.0943E-02  2.4218E-01  1.3450E+00  6.4664E-02  4.3673E-01  4.5466E-02  1.2226E-01  1.2553E+00  6.4552E-02  3.2098E-01
             7.7342E-01
 GRADIENT:   1.1425E+01  9.0453E+00  2.7553E+00 -7.2024E-01 -8.4620E+00  1.3377E+00  1.2033E+00 -1.4862E+00  8.8984E-01  5.9338E-01
            -6.7508E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3067.93250452506        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      668
 NPARAMETR:  9.8113E-01  1.1528E+00  3.4729E+00  9.6528E-01  1.4004E+00  9.4678E-01  1.0225E+00  3.1751E+00  9.6063E-01  1.2473E+00
             1.9609E+00
 PARAMETER:  8.0954E-02  2.4218E-01  1.3450E+00  6.4664E-02  4.3673E-01  4.5312E-02  1.2226E-01  1.2553E+00  5.9832E-02  3.2098E-01
             7.7343E-01
 GRADIENT:   3.0416E-02  4.1937E+00  1.5066E+00 -3.3844E+00 -1.2943E+01  8.3593E-02  7.0642E-01 -2.5139E+00 -4.5850E-02  8.0706E-02
            -8.2286E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3067.96897302780        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  9.8141E-01  1.1526E+00  3.4616E+00  9.6529E-01  1.4023E+00  9.4745E-01  1.0225E+00  3.1820E+00  9.5773E-01  1.2473E+00
             1.9722E+00
 PARAMETER:  8.1230E-02  2.4206E-01  1.3417E+00  6.4675E-02  4.3812E-01  4.6018E-02  1.2226E-01  1.2575E+00  5.6809E-02  3.2095E-01
             7.7913E-01
 GRADIENT:   3.9404E-01  3.0696E+00  8.2682E-01 -3.3816E+00 -1.1500E+01  4.0010E-01  6.4356E-01 -1.7594E+00 -3.5254E-01  2.8739E-01
             4.2055E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3068.04054570154        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  9.8104E-01  1.1505E+00  3.4455E+00  9.6648E-01  1.4113E+00  9.4632E-01  1.0079E+00  3.2112E+00  9.6015E-01  1.2461E+00
             1.9699E+00
 PARAMETER:  8.0862E-02  2.4019E-01  1.3371E+00  6.5902E-02  4.4449E-01  4.4826E-02  1.0789E-01  1.2666E+00  5.9332E-02  3.1999E-01
             7.7799E-01
 GRADIENT:  -5.2903E-01  5.7487E-02 -1.0715E+00 -2.1503E+00 -3.6727E+00 -8.3038E-02 -4.1410E-01 -6.8973E-01 -5.9425E-01 -3.4323E-01
             1.8778E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3068.04404327839        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1171
 NPARAMETR:  9.8131E-01  1.1485E+00  3.4443E+00  9.6857E-01  1.4111E+00  9.4656E-01  1.0117E+00  3.2102E+00  9.6116E-01  1.2460E+00
             1.9703E+00
 PARAMETER:  8.1133E-02  2.3850E-01  1.3367E+00  6.8070E-02  4.4438E-01  4.5083E-02  1.1163E-01  1.2663E+00  6.0386E-02  3.1991E-01
             7.7819E-01
 GRADIENT:   1.1609E-01  5.0654E-01 -1.3418E+00 -7.0172E-01 -3.0609E+00  1.1175E-02 -9.8137E-02 -6.4560E-01 -5.2777E-02 -1.8378E-01
             2.4543E+00

0ITERATION NO.:   51    OBJECTIVE VALUE:  -3068.04404327839        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1201
 NPARAMETR:  9.8131E-01  1.1486E+00  3.4441E+00  9.6864E-01  1.4111E+00  9.4656E-01  1.0118E+00  3.2100E+00  9.6118E-01  1.2460E+00
             1.9704E+00
 PARAMETER:  8.1133E-02  2.3850E-01  1.3367E+00  6.8070E-02  4.4438E-01  4.5083E-02  1.1163E-01  1.2663E+00  6.0386E-02  3.1991E-01
             7.7819E-01
 GRADIENT:   1.1259E-01 -8.7948E+03  1.3755E+03 -6.6460E-01  4.7076E+03  1.1741E-02 -9.2556E-02  1.6569E+03 -4.6755E-02  6.5492E+03
            -2.6944E+03
 NUMSIGDIG:         3.3         3.3         3.3         2.1         3.3         3.6         1.9         3.3         2.5         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1201
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1424E-03 -2.6767E-02 -3.4430E-02  1.4760E-02 -3.6113E-02
 SE:             2.9642E-02  1.9618E-02  2.0035E-02  2.3846E-02  2.3204E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6926E-01  1.7244E-01  8.5698E-02  5.3595E-01  1.1964E-01

 ETASHRINKSD(%)  6.9496E-01  3.4276E+01  3.2882E+01  2.0111E+01  2.2263E+01
 ETASHRINKVR(%)  1.3851E+00  5.6804E+01  5.4952E+01  3.6178E+01  3.9569E+01
 EBVSHRINKSD(%)  1.0194E+00  3.4495E+01  3.5918E+01  2.2576E+01  1.9150E+01
 EBVSHRINKVR(%)  2.0284E+00  5.7091E+01  5.8935E+01  4.0055E+01  3.4632E+01
 RELATIVEINF(%)  9.7944E+01  7.3771E+00  2.0292E+01  1.0670E+01  3.9994E+01
 EPSSHRINKSD(%)  1.8771E+01
 EPSSHRINKVR(%)  3.4018E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3068.0440432783939     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1413.9546835099832     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.33
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3068.044       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.15E+00  3.44E+00  9.69E-01  1.41E+00  9.47E-01  1.01E+00  3.21E+00  9.61E-01  1.25E+00  1.97E+00
 


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
+        1.27E+03
 
 TH 2
+       -1.18E+03  6.99E+06
 
 TH 3
+        4.53E+01  3.02E+02  2.16E+04
 
 TH 4
+       -1.09E+01  6.76E+03 -2.53E+02  5.59E+07
 
 TH 5
+        2.86E+02  3.61E+05 -2.23E+04 -1.01E+06  1.17E+06
 
 TH 6
+       -2.49E+00 -4.07E+01  5.73E+00 -3.81E+00  1.28E+01  2.29E+02
 
 TH 7
+        1.11E+00 -2.76E+03  2.01E+02 -4.80E+07  8.68E+05 -1.30E+01  4.64E+01
 
 TH 8
+        4.81E+01  2.14E+02  5.41E+01 -1.02E+02  2.45E+04  2.47E+00  7.29E+01  3.18E+04
 
 TH 9
+        8.90E+00  1.29E+03 -1.68E+02  5.64E+07 -1.64E+02  6.86E+00  3.67E+01 -3.19E+01  9.25E+01
 
 TH10
+        4.94E+02  2.23E+03 -5.32E+02 -1.06E+03 -2.51E+05  2.23E+01  7.46E+02  3.24E+05 -3.34E+02  3.29E+06
 
 TH11
+       -1.17E+02 -3.49E+02 -9.91E+01  2.64E+01  1.72E+03 -2.96E+00 -9.59E+01 -1.64E+02  5.00E+01  1.35E+03  2.24E+05
 
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
 #CPUT: Total CPU Time in Seconds,       48.494
Stop Time:
Sat Sep 18 04:03:01 CDT 2021
