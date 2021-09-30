Wed Sep 29 20:22:16 CDT 2021
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
$DATA ../../../../data/spa/D/dat84.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20198.1940566683        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.6467E+02  5.0111E+02 -3.3851E+01  4.4194E+02  2.4424E+01 -2.1005E+03 -1.0298E+03 -3.9635E+01 -1.5330E+03 -3.9331E+02
            -3.8087E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -526.172765440094        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2480E+00  1.0476E+00  9.1367E-01  1.4378E+00  1.5716E+00  1.4758E+00  1.0761E+00  9.6260E-01  8.5569E-01  9.1114E-01
             1.5401E+01
 PARAMETER:  3.2156E-01  1.4649E-01  9.7125E-03  4.6314E-01  5.5212E-01  4.8923E-01  1.7334E-01  6.1884E-02 -5.5849E-02  6.9446E-03
             2.8344E+00
 GRADIENT:  -2.6064E+01  1.9192E+01 -6.4096E+00  3.6963E+01 -3.5245E+00  2.7407E+01  3.9581E-01  3.2673E+00  4.6251E+00  4.9512E-01
             6.6764E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -535.172779822199        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2299E+00  9.3847E-01  1.2552E+00  1.3299E+00  3.5799E+00  1.1966E+00  1.3169E+00  7.6481E-01  8.0639E-01  1.0169E+00
             1.5282E+01
 PARAMETER:  3.0693E-01  3.6499E-02  3.2730E-01  3.8514E-01  1.3753E+00  2.7948E-01  3.7531E-01 -1.6813E-01 -1.1518E-01  1.1673E-01
             2.8267E+00
 GRADIENT:   4.8346E+00  9.5732E+00  3.7155E+00  8.0801E+00 -3.9307E+00 -8.5193E+00  3.4820E+00  3.0617E-01  7.1568E+00  1.5345E-02
             6.4387E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -553.836547001999        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.1049E+00  4.6997E-01  7.0236E-01  1.4568E+00  6.8748E+00  1.2148E+00  1.5351E+00  7.1440E-02  5.5275E-01  8.6476E+00
             1.4149E+01
 PARAMETER:  1.9972E-01 -6.5509E-01 -2.5331E-01  4.7627E-01  2.0279E+00  2.9457E-01  5.2860E-01 -2.5389E+00 -4.9285E-01  2.2573E+00
             2.7496E+00
 GRADIENT:  -3.6315E+01  1.6666E+01 -1.4730E+00  5.0142E+01 -4.7975E+00 -1.3136E+01  1.3748E+00  7.1123E-03  7.1899E+00  9.5207E+00
             3.0718E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -585.536743618287        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  8.9318E-01  1.2624E-01  1.9996E-01  1.1188E+00  1.3054E+01  1.4339E+00  1.0898E+00  1.0000E-02  1.0581E-01  5.3447E+00
             1.3044E+01
 PARAMETER: -1.2963E-02 -1.9696E+00 -1.5096E+00  2.1225E-01  2.6691E+00  4.6037E-01  1.8602E-01 -6.8143E+00 -2.1461E+00  1.7761E+00
             2.6684E+00
 GRADIENT:  -1.5351E+01 -8.2044E+00 -9.6884E-01  1.0308E+02 -2.2477E-01 -4.2764E+00  8.8421E-01  0.0000E+00  5.3918E-01  1.9459E+01
            -1.9069E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -610.243137937610        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  5.1248E-01  2.1143E-02  2.8831E-02  3.5547E-01  1.4207E+01  1.5119E+00  3.1473E-01  1.0000E-02  1.0000E-02  4.2111E+00
             1.2540E+01
 PARAMETER: -5.6850E-01 -3.7564E+00 -3.4463E+00 -9.3431E-01  2.7537E+00  5.1333E-01 -1.0560E+00 -1.0321E+01 -7.4881E+00  1.5377E+00
             2.6289E+00
 GRADIENT:   6.0509E+01 -7.4207E+00 -8.6455E+01  1.6985E+02  1.0037E+01  1.6625E+01  1.0097E-01  0.0000E+00  0.0000E+00  6.1282E+00
            -3.7066E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -622.306210122900        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  3.1448E-01  1.0000E-02  1.0000E-02  1.3465E-01  1.4743E+01  1.2662E+00  5.6533E-02  1.0000E-02  1.0000E-02  3.6563E+00
             1.2142E+01
 PARAMETER: -1.0568E+00 -5.3372E+00 -4.7630E+00 -1.9051E+00  2.7908E+00  3.3604E-01 -2.7729E+00 -1.2209E+01 -1.2065E+01  1.3964E+00
             2.5967E+00
 GRADIENT:   6.9080E+01  0.0000E+00  0.0000E+00 -8.6739E+00  4.7655E+00 -2.0659E+01  4.7804E-03  0.0000E+00  0.0000E+00  7.2475E-01
             1.1588E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -623.784587676248        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      656
 NPARAMETR:  3.1716E-01  1.0000E-02  1.0000E-02  1.3656E-01  1.4672E+01  1.3713E+00  4.2185E-02  1.0000E-02  1.0000E-02  3.7256E+00
             1.2371E+01
 PARAMETER: -1.0484E+00 -5.3372E+00 -4.7630E+00 -1.8910E+00  2.7860E+00  4.1572E-01 -3.0657E+00 -1.2209E+01 -1.2065E+01  1.4152E+00
             2.6153E+00
 GRADIENT:   2.0545E+00  0.0000E+00  0.0000E+00 -1.5410E+01  1.7126E+00  1.9733E+00  2.8874E-03  0.0000E+00  0.0000E+00  5.7233E-01
            -5.2544E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -623.905277516418        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  3.1735E-01  1.0000E-02  1.0000E-02  1.3801E-01  1.4347E+01  1.3569E+00  1.7216E-02  1.0000E-02  1.0000E-02  3.5468E+00
             1.2435E+01
 PARAMETER: -1.0477E+00 -5.3372E+00 -4.7630E+00 -1.8805E+00  2.7636E+00  4.0521E-01 -3.9619E+00 -1.2209E+01 -1.2065E+01  1.3661E+00
             2.6205E+00
 GRADIENT:  -1.1912E-02  0.0000E+00  0.0000E+00  6.5782E+00  1.6615E+00 -1.3105E+00  5.0679E-04  0.0000E+00  0.0000E+00  4.7801E-01
            -2.9118E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -624.371785980117        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  3.1820E-01  1.0000E-02  1.0000E-02  1.3796E-01  1.0057E+01  1.3657E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.9144E+00
             1.2486E+01
 PARAMETER: -1.0451E+00 -5.3372E+00 -4.7630E+00 -1.8808E+00  2.4083E+00  4.1165E-01 -1.5052E+01 -1.2209E+01 -1.2065E+01  7.4939E-01
             2.6246E+00
 GRADIENT:   1.4068E+00  0.0000E+00  0.0000E+00  4.1835E+00 -3.2779E+00  2.4277E-01  0.0000E+00  0.0000E+00  0.0000E+00  2.1537E+00
            -2.8320E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -624.397369992749        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1187            RESET HESSIAN, TYPE II
 NPARAMETR:  3.1852E-01  1.0000E-02  1.0000E-02  1.3783E-01  9.7775E+00  1.3625E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7997E+00
             1.2523E+01
 PARAMETER: -1.0441E+00 -5.3372E+00 -4.7630E+00 -1.8817E+00  2.3801E+00  4.0934E-01 -1.6174E+01 -1.2209E+01 -1.2065E+01  6.8763E-01
             2.6276E+00
 GRADIENT:   6.6441E+01  0.0000E+00  0.0000E+00  2.9365E+01  1.1120E+00  3.6292E+00  0.0000E+00  0.0000E+00  0.0000E+00  5.7724E-01
             2.0803E+01

0ITERATION NO.:   53    OBJECTIVE VALUE:  -624.399274866182        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1279
 NPARAMETR:  3.1853E-01  1.0000E-02  1.0000E-02  1.3784E-01  9.7840E+00  1.3626E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7796E+00
             1.2524E+01
 PARAMETER: -1.0440E+00 -5.3372E+00 -4.7630E+00 -1.8817E+00  2.3807E+00  4.0936E-01 -1.6174E+01 -1.2209E+01 -1.2065E+01  6.7641E-01
             2.6277E+00
 GRADIENT:   1.1132E+00  0.0000E+00  0.0000E+00 -1.1723E+00  3.3694E-01  8.6787E-02  0.0000E+00  0.0000E+00  0.0000E+00  3.8843E-03
            -5.7988E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1279
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2826E-03  7.7380E-06  3.3461E-05 -8.2171E-05 -1.2436E-02
 SE:             2.8520E-02  3.0532E-05  2.5394E-04  3.1120E-04  6.9772E-03
 N:                     100         100         100         100         100

 P VAL.:         9.3621E-01  7.9993E-01  8.9517E-01  7.9174E-01  7.4689E-02

 ETASHRINKSD(%)  4.4556E+00  9.9898E+01  9.9149E+01  9.8957E+01  7.6625E+01
 ETASHRINKVR(%)  8.7126E+00  1.0000E+02  9.9993E+01  9.9989E+01  9.4536E+01
 EBVSHRINKSD(%)  4.5931E+00  9.9853E+01  9.9017E+01  9.8841E+01  8.3198E+01
 EBVSHRINKVR(%)  8.9752E+00  1.0000E+02  9.9990E+01  9.9987E+01  9.7177E+01
 RELATIVEINF(%)  1.7310E+00  1.8691E-05  1.8315E-05  2.3792E-05  1.6132E+00
 EPSSHRINKSD(%)  6.7671E+00
 EPSSHRINKVR(%)  1.3076E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -624.39927486618171     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       110.75155169755647     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.25
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -624.399       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.19E-01  1.00E-02  1.00E-02  1.38E-01  9.78E+00  1.36E+00  1.00E-02  1.00E-02  1.00E-02  1.78E+00  1.25E+01
 


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
+        5.28E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        0.00E+00  0.00E+00  0.00E+00
 
 TH 4
+       -1.39E+03  0.00E+00  0.00E+00  5.36E+04
 
 TH 5
+        8.69E-01  0.00E+00  0.00E+00 -5.14E+00  4.08E+00
 
 TH 6
+       -7.89E+00  0.00E+00  0.00E+00 -6.46E+01  1.82E-01  8.21E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -4.22E+00  0.00E+00  0.00E+00  2.85E+01 -2.21E+00  1.52E-01  0.00E+00  0.00E+00  0.00E+00  3.15E+00
 
 TH11
+       -3.55E+01  0.00E+00  0.00E+00 -2.01E+01  9.94E-02  1.27E+00  0.00E+00  0.00E+00  0.00E+00 -7.17E-02  2.37E+00
 
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
 #CPUT: Total CPU Time in Seconds,       23.059
Stop Time:
Wed Sep 29 20:22:41 CDT 2021
