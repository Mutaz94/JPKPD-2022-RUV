Thu Sep 30 03:39:07 CDT 2021
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
$DATA ../../../../data/spa1/D/dat80.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26564.1128197719        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.6779E+02  5.2914E+02 -9.2166E+01  4.2348E+02  1.2175E+02 -2.1829E+03 -1.0752E+03 -4.4371E+01 -1.7854E+03 -4.1042E+02
            -5.1216E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -528.607737281815        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3878E+00  1.1768E+00  9.2384E-01  1.7983E+00  1.3103E+00  2.0275E+00  1.1565E+00  9.4039E-01  1.3214E+00  9.5735E-01
             1.4559E+01
 PARAMETER:  4.2770E-01  2.6281E-01  2.0787E-02  6.8686E-01  3.7027E-01  8.0678E-01  2.4539E-01  3.8535E-02  3.7868E-01  5.6415E-02
             2.7782E+00
 GRADIENT:   9.3191E+00  5.9235E+01 -1.0527E+01  1.0272E+02 -3.6412E+00  5.0763E+01 -7.4689E+00  5.6680E+00 -3.9505E+01  2.0743E+00
             8.7171E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -555.649869301162        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.3680E+00  1.2166E+00  1.2880E+00  1.7566E+00  4.4117E+00  1.8468E+00  5.1217E-01  3.6375E-01  2.7256E+00  8.3095E-01
             1.4388E+01
 PARAMETER:  4.1338E-01  2.9610E-01  3.5309E-01  6.6338E-01  1.5843E+00  7.1345E-01 -5.6909E-01 -9.1128E-01  1.1027E+00 -8.5188E-02
             2.7664E+00
 GRADIENT:  -6.2749E-01  1.9652E+01  1.1197E+00  4.2369E+01 -7.2779E+00 -1.9675E+01  1.6857E+00  3.6613E-01  2.2900E+01  2.8321E-01
             1.0985E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -601.925704398217        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.1454E+00  3.1889E-01  4.0167E-01  1.4802E+00  4.4128E+01  1.5392E+00  2.8778E-01  1.0000E-02  1.0589E+00  1.0777E+01
             1.2263E+01
 PARAMETER:  2.3576E-01 -1.0429E+00 -8.1213E-01  4.9216E-01  3.8871E+00  5.3123E-01 -1.1456E+00 -9.6729E+00  1.5727E-01  2.4774E+00
             2.6066E+00
 GRADIENT:   5.0572E+01  2.1336E+01  1.3205E+01  7.9375E+01 -6.8269E-02 -3.5409E+01  2.4647E-01  0.0000E+00 -3.2674E+01  6.9806E-02
            -6.0674E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -649.881175993835        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  6.5658E-01  4.0673E-02  6.0878E-02  5.6815E-01  1.2519E+04  1.5496E+00  1.0000E-02  1.0000E-02  7.8507E-01  1.2172E+01
             9.7455E+00
 PARAMETER: -3.2071E-01 -3.1022E+00 -2.6989E+00 -4.6537E-01  9.5350E+00  5.3797E-01 -6.1130E+00 -2.2227E+01 -1.4198E-01  2.5991E+00
             2.3768E+00
 GRADIENT:   7.7635E+01  8.1090E-01 -2.2895E+01  1.1083E+02  6.0276E-04 -4.1756E+00  0.0000E+00  0.0000E+00 -2.9497E+01  9.6706E-07
            -2.0490E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -673.456875945345        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      467
 NPARAMETR:  6.4708E-01  4.0820E-02  6.2609E-02  5.5777E-01  6.1029E+03  1.4849E+00  1.0000E-02  1.0000E-02  8.0992E-01  8.4904E+00
             1.1656E+01
 PARAMETER: -3.3529E-01 -3.0986E+00 -2.6708E+00 -4.8380E-01  8.8165E+00  4.9535E-01 -8.3507E+00 -2.1417E+01 -1.1082E-01  2.2389E+00
             2.5558E+00
 GRADIENT:  -9.0508E+00  2.7981E+00  2.1114E+01 -5.0872E+00 -1.3981E-04 -1.2994E+01  0.0000E+00  0.0000E+00 -5.2543E+00  1.1622E-06
            -2.1495E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -685.475575092994        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  3.7601E-01  1.0000E-02  1.4293E-02  1.9041E-01  2.3569E+05  1.3757E+00  1.0000E-02  1.0000E-02  4.7911E-01  4.6069E+00
             1.2820E+01
 PARAMETER: -8.7814E-01 -4.9562E+00 -4.1480E+00 -1.5586E+00  1.2470E+01  4.1895E-01 -1.6236E+01 -3.1137E+01 -6.3583E-01  1.6276E+00
             2.6510E+00
 GRADIENT:  -1.9546E+00  0.0000E+00 -1.5520E+01  1.5194E+01 -8.3978E-07 -2.6510E-02  0.0000E+00  0.0000E+00  6.4284E+00  6.9851E-12
             3.3641E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -688.230719606340        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      832             RESET HESSIAN, TYPE I
 NPARAMETR:  3.7887E-01  1.0000E-02  1.5184E-02  2.0090E-01  7.0391E+05  1.3882E+00  1.0000E-02  1.0000E-02  7.8396E-02  4.8916E+00
             1.2503E+01
 PARAMETER: -8.7057E-01 -4.8198E+00 -4.0875E+00 -1.5049E+00  1.3564E+01  4.2800E-01 -1.3532E+01 -3.3022E+01 -2.4460E+00  1.6875E+00
             2.6259E+00
 GRADIENT:   6.5240E+01  0.0000E+00  9.9354E+01  1.2292E+01  1.0546E-07  8.7767E+00  0.0000E+00  0.0000E+00  1.1331E-01  0.0000E+00
             2.2041E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -688.397869539160        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      970
 NPARAMETR:  3.7402E-01  1.0000E-02  1.5008E-02  2.0040E-01  3.5186E+06  1.3734E+00  1.0000E-02  1.0000E-02  4.0300E-02  4.9263E+00
             1.2500E+01
 PARAMETER: -8.8345E-01 -4.8198E+00 -4.0992E+00 -1.5074E+00  1.5174E+01  4.1727E-01 -1.3532E+01 -3.3022E+01 -3.1114E+00  1.6946E+00
             2.6258E+00
 GRADIENT:  -1.8909E+00  0.0000E+00 -5.5472E-01  2.0744E-01 -2.5238E-08  3.6292E-01  0.0000E+00  0.0000E+00  3.2330E-02  0.0000E+00
            -3.4789E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -688.407951093714        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1146
 NPARAMETR:  3.7595E-01  1.0000E-02  1.5046E-02  2.0072E-01  4.3336E+06  1.3694E+00  1.0000E-02  1.0000E-02  2.5720E-02  4.9221E+00
             1.2586E+01
 PARAMETER: -8.7830E-01 -4.8198E+00 -4.0966E+00 -1.5058E+00  1.5382E+01  4.1441E-01 -1.3532E+01 -3.3022E+01 -3.5605E+00  1.6937E+00
             2.6326E+00
 GRADIENT:   4.4578E-01  0.0000E+00  2.0104E-01 -1.8422E+00 -3.7258E-09 -3.3426E-01  0.0000E+00  0.0000E+00  1.3387E-02  0.0000E+00
             1.5742E+00

0ITERATION NO.:   49    OBJECTIVE VALUE:  -688.420740943897        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1283
 NPARAMETR:  3.7589E-01  1.0000E-02  1.5008E-02  2.0078E-01  7.0986E+06  1.3713E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.9516E+00
             1.2566E+01
 PARAMETER: -8.7845E-01 -4.8198E+00 -4.0992E+00 -1.5056E+00  1.5875E+01  4.1575E-01 -1.3532E+01 -3.3022E+01 -4.9662E+00  1.6997E+00
             2.6310E+00
 GRADIENT:   1.9252E+00  0.0000E+00 -5.0561E+00  3.8009E+00  1.8617E-08 -8.8240E-02  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -5.2981E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1283
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1069E-03  4.6865E-06  6.5461E-05 -1.7684E-04  9.4883E-11
 SE:             2.8702E-02  2.2536E-06  3.2813E-04  3.6694E-04  8.0603E-10
 N:                     100         100         100         100         100

 P VAL.:         9.6924E-01  3.7562E-02  8.4187E-01  6.2986E-01  9.0629E-01

 ETASHRINKSD(%)  3.8436E+00  9.9992E+01  9.8901E+01  9.8771E+01  1.0000E+02
 ETASHRINKVR(%)  7.5395E+00  1.0000E+02  9.9988E+01  9.9985E+01  1.0000E+02
 EBVSHRINKSD(%)  4.1791E+00  9.9989E+01  9.8925E+01  9.8750E+01  1.0000E+02
 EBVSHRINKVR(%)  8.1836E+00  1.0000E+02  9.9988E+01  9.9984E+01  1.0000E+02
 RELATIVEINF(%)  3.4263E+00  6.4040E-07  3.8484E-05  5.4276E-05  0.0000E+00
 EPSSHRINKSD(%)  4.8100E+00
 EPSSHRINKVR(%)  9.3886E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -688.42074094389750     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       230.51779226077520     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -688.421       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.76E-01  1.00E-02  1.50E-02  2.01E-01  7.10E+06  1.37E+00  1.00E-02  1.00E-02  1.00E-02  4.95E+00  1.26E+01
 


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
+        3.78E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -3.75E+04  0.00E+00  3.94E+06
 
 TH 4
+       -6.63E+01  0.00E+00 -3.25E+05  3.00E+04
 
 TH 5
+        3.21E-12  0.00E+00 -3.25E-11  2.63E-12 -3.65E-21
 
 TH 6
+       -3.50E+00  0.00E+00  1.22E+03 -1.46E+02 -7.18E-13  8.66E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -6.17E-06  0.00E+00 -1.02E-05 -4.63E-06 -2.37E-14  3.48E-06  0.00E+00  0.00E+00  0.00E+00 -9.04E-07
 
 TH11
+       -3.38E+01  0.00E+00  5.88E+02 -2.60E+01 -2.17E-14  9.76E-01  0.00E+00  0.00E+00  0.00E+00  2.25E-07  3.03E+00
 
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
 #CPUT: Total CPU Time in Seconds,       29.946
Stop Time:
Thu Sep 30 03:39:38 CDT 2021
