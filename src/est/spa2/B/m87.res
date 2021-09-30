Thu Sep 30 04:30:57 CDT 2021
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
$DATA ../../../../data/spa2/B/dat87.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2174.83179797903        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1996E+02 -3.5357E+01  3.1740E+01 -2.9013E+01  1.2719E+02  4.3014E+01  6.5827E+00 -4.6237E+02 -1.1458E+02  1.1730E+01
            -7.7624E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2358.93126938147        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0145E+00  1.0023E+00  1.2414E+00  1.0192E+00  9.9857E-01  1.0170E+00  9.2535E-01  2.0046E+00  7.3452E-01  8.8544E-01
             9.5668E-01
 PARAMETER:  1.1438E-01  1.0228E-01  3.1627E-01  1.1898E-01  9.8567E-02  1.1687E-01  2.2412E-02  7.9546E-01 -2.0854E-01 -2.1665E-02
             5.5717E-02
 GRADIENT:   5.7565E+02  1.4189E+01  3.2138E+01 -4.4926E+01 -1.2321E+01  6.9854E+01 -2.3653E+01 -1.1285E+02 -7.5836E+01 -1.1639E+01
            -1.0069E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2379.68761487903        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      242
 NPARAMETR:  1.0169E+00  1.0099E+00  2.0378E+00  1.1408E+00  1.2113E+00  1.0153E+00  1.3376E+00  2.3122E+00  7.2050E-01  1.0262E+00
             1.0071E+00
 PARAMETER:  1.1674E-01  1.0983E-01  8.1187E-01  2.3174E-01  2.9171E-01  1.1516E-01  3.9086E-01  9.3819E-01 -2.2781E-01  1.2590E-01
             1.0704E-01
 GRADIENT:   9.7583E+00  5.8048E+01  3.0639E+01  6.9770E+01  3.2914E+01 -8.1139E+00 -1.8203E+00 -1.4092E+02 -3.4093E+01 -1.3423E+01
            -4.9667E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2403.56721826151        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      418
 NPARAMETR:  1.0113E+00  9.7138E-01  1.6262E+00  1.1236E+00  1.0988E+00  1.0630E+00  1.2905E+00  2.6376E+00  7.4725E-01  1.0181E+00
             1.0412E+00
 PARAMETER:  1.1120E-01  7.0958E-02  5.8622E-01  2.1652E-01  1.9426E-01  1.6108E-01  3.5501E-01  1.0699E+00 -1.9136E-01  1.1792E-01
             1.4034E-01
 GRADIENT:  -2.8033E+00  1.6949E+01  8.1850E+00  2.5913E+01 -2.0828E+00  1.0387E+01 -1.0930E+00 -6.9288E+01 -2.9235E+01  3.1658E+00
             8.8433E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2409.42073042831        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      608             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0122E+00  9.9040E-01  1.5706E+00  1.1056E+00  1.0916E+00  1.0453E+00  1.2009E+00  2.7530E+00  8.6532E-01  9.7896E-01
             1.0374E+00
 PARAMETER:  1.1214E-01  9.0351E-02  5.5146E-01  2.0038E-01  1.8768E-01  1.4428E-01  2.8304E-01  1.1127E+00 -4.4651E-02  7.8738E-02
             1.3676E-01
 GRADIENT:   4.7933E+02  6.3840E+01  2.4432E+01  1.9472E+02  3.2400E+01  8.7500E+01  4.1249E+01 -4.2346E+00 -2.0330E+00  1.5891E+00
             4.8912E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2409.54602319726        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  1.0122E+00  9.8940E-01  1.5706E+00  1.1034E+00  1.0916E+00  1.0314E+00  1.1870E+00  2.7537E+00  8.6859E-01  9.7897E-01
             1.0374E+00
 PARAMETER:  1.1215E-01  8.9339E-02  5.5144E-01  1.9835E-01  1.8766E-01  1.3095E-01  2.7143E-01  1.1130E+00 -4.0881E-02  7.8749E-02
             1.3674E-01
 GRADIENT:  -8.9704E-01  5.9431E+00  4.8242E+00  1.0140E+01 -4.6039E+00 -1.4301E+00  2.2270E+00 -4.8744E+01 -9.4835E+00 -3.2321E-01
             3.0899E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2409.81638866548        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  1.0122E+00  9.8351E-01  1.5706E+00  1.0951E+00  1.0916E+00  1.0282E+00  1.1597E+00  2.7544E+00  8.9274E-01  9.7897E-01
             1.0374E+00
 PARAMETER:  1.1215E-01  8.3370E-02  5.5143E-01  1.9084E-01  1.8766E-01  1.2784E-01  2.4818E-01  1.1132E+00 -1.3455E-02  7.8749E-02
             1.3674E-01
 GRADIENT:   4.7968E+02  4.8119E+01  2.5517E+01  1.6284E+02  3.5120E+01  7.1391E+01  3.4405E+01 -3.4161E+00  1.7268E+00  2.3593E+00
             4.8478E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2409.99989258378        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1070
 NPARAMETR:  1.0122E+00  9.9032E-01  1.5706E+00  1.0948E+00  1.0916E+00  1.0313E+00  1.1354E+00  2.7546E+00  9.1483E-01  9.7897E-01
             1.0374E+00
 PARAMETER:  1.1215E-01  9.0270E-02  5.5143E-01  1.9057E-01  1.8766E-01  1.3081E-01  2.2697E-01  1.1133E+00  1.0979E-02  7.8749E-02
             1.3674E-01
 GRADIENT:  -8.3035E-01 -4.2991E+00  6.4853E+00 -3.1575E-01 -3.3223E+00 -1.4862E+00  1.1124E+00 -4.7934E+01 -2.7597E+00  5.7044E-01
             2.8147E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2410.01278379065        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1214
 NPARAMETR:  1.0122E+00  9.8989E-01  1.5705E+00  1.0937E+00  1.0916E+00  1.0297E+00  1.1243E+00  2.7552E+00  9.1546E-01  9.7897E-01
             1.0374E+00
 PARAMETER:  1.1215E-01  8.9835E-02  5.5142E-01  1.8961E-01  1.8766E-01  1.2931E-01  2.1717E-01  1.1135E+00  1.1668E-02  7.8749E-02
             1.3674E-01
 GRADIENT:   4.7966E+02  5.0058E+01  2.6503E+01  1.6454E+02  3.4375E+01  7.2827E+01  2.8276E+01 -3.0589E+00  4.4938E+00  2.7394E+00
             4.3291E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2410.02826788842        NO. OF FUNC. EVALS.: 209
 CUMULATIVE NO. OF FUNC. EVALS.:     1423
 NPARAMETR:  1.0122E+00  9.9013E-01  1.5705E+00  1.0932E+00  1.0916E+00  1.0360E+00  1.1222E+00  2.7557E+00  9.1500E-01  9.7897E-01
             1.0374E+00
 PARAMETER:  1.1215E-01  9.0086E-02  5.5141E-01  1.8913E-01  1.8766E-01  1.3536E-01  2.1525E-01  1.1137E+00  1.1168E-02  7.8749E-02
             1.3674E-01
 GRADIENT:  -8.0250E-01 -6.2880E+00  6.6817E+00 -2.6557E+00 -3.1706E+00  3.4767E-01 -3.7327E-02 -4.7814E+01 -3.4949E+00  6.6303E-01
             2.6860E+00

0ITERATION NO.:   46    OBJECTIVE VALUE:  -2410.02826788842        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1445
 NPARAMETR:  1.0122E+00  9.9013E-01  1.5705E+00  1.0932E+00  1.0916E+00  1.0360E+00  1.1222E+00  2.7557E+00  9.1500E-01  9.7897E-01
             1.0374E+00
 PARAMETER:  1.1215E-01  9.0086E-02  5.5141E-01  1.8913E-01  1.8766E-01  1.3536E-01  2.1525E-01  1.1137E+00  1.1168E-02  7.8749E-02
             1.3674E-01
 GRADIENT:  -8.0250E-01 -6.2880E+00  6.6817E+00 -2.6557E+00 -3.1706E+00  3.4767E-01 -3.7327E-02 -4.7814E+01 -3.4949E+00  6.6303E-01
             2.6860E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1445
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0539E-03 -1.6147E-02 -4.6938E-02  8.8966E-03 -5.2158E-02
 SE:             2.9893E-02  1.9456E-02  2.7224E-02  2.5401E-02  2.0332E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4522E-01  4.0658E-01  8.4689E-02  7.2615E-01  1.0306E-02

 ETASHRINKSD(%)  1.0000E-10  3.4822E+01  8.7950E+00  1.4904E+01  3.1887E+01
 ETASHRINKVR(%)  1.0000E-10  5.7518E+01  1.6816E+01  2.7587E+01  5.3606E+01
 EBVSHRINKSD(%)  3.3262E-01  3.4859E+01  2.2431E+01  1.7556E+01  3.0355E+01
 EBVSHRINKVR(%)  6.6413E-01  5.7566E+01  3.9830E+01  3.2029E+01  5.1496E+01
 RELATIVEINF(%)  9.9314E+01  1.1733E+01  3.4699E+01  2.0882E+01  2.5756E+01
 EPSSHRINKSD(%)  3.1956E+01
 EPSSHRINKVR(%)  5.3700E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2410.0282678884150     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1307.3020280428079     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2410.028       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  9.90E-01  1.57E+00  1.09E+00  1.09E+00  1.04E+00  1.12E+00  2.76E+00  9.15E-01  9.79E-01  1.04E+00
 


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
+        8.81E+06
 
 TH 2
+       -5.64E+00  1.16E+07
 
 TH 3
+        1.17E+00  2.93E+01  1.51E+05
 
 TH 4
+       -3.23E+00  3.92E+02 -3.48E+01  1.33E+06
 
 TH 5
+       -2.44E+06 -1.54E+02 -8.20E+01  5.04E+01  2.70E+06
 
 TH 6
+        2.64E+00 -1.66E+00  3.51E-01 -8.73E-01 -1.15E-01  2.88E+06
 
 TH 7
+        5.40E-01  3.31E+01 -1.23E+00 -7.66E+00 -1.05E+01  4.80E+01  9.73E+05
 
 TH 8
+       -1.71E-01 -3.80E+00 -9.29E+00  5.43E+00 -6.89E+00 -2.06E-01 -3.62E-01  1.21E+04
 
 TH 9
+        1.12E+00 -6.26E+06  1.46E+00  3.00E+06  1.38E+01  1.27E+02 -5.13E+06  4.37E+00  6.78E+06
 
 TH10
+        5.10E+06 -1.10E+01 -4.63E-01  4.93E+00  2.73E+01 -2.38E-02 -4.42E+00  8.78E-01  5.18E+00  5.92E+06
 
 TH11
+       -3.52E+06 -1.98E+01 -8.75E+00  1.94E+06 -3.64E+01  2.53E+00  6.80E+00  3.76E+00  3.37E+00  1.94E+01  5.64E+06
 
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
 #CPUT: Total CPU Time in Seconds,       40.324
Stop Time:
Thu Sep 30 04:31:39 CDT 2021
