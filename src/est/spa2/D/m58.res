Thu Sep 30 09:22:02 CDT 2021
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
$DATA ../../../../data/spa2/D/dat58.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m58.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16252.3590801469        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9603E+02  3.9255E+02  1.6975E+01  2.6765E+02  4.3742E+02 -2.5958E+03 -9.6020E+02 -7.1093E+01 -1.5402E+03 -1.0945E+03
            -3.0832E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -572.690723930911        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1065E+00  9.3937E-01  6.4746E-01  2.0788E+00  1.3034E+00  4.0414E+00  2.7654E+00  9.5330E-01  2.6457E+00  1.8344E+00
             1.0917E+01
 PARAMETER:  2.0119E-01  3.7456E-02 -3.3471E-01  8.3179E-01  3.6499E-01  1.4966E+00  1.1172E+00  5.2177E-02  1.0729E+00  7.0669E-01
             2.4903E+00
 GRADIENT:  -1.5147E+01 -7.8351E+01 -9.2356E+01  9.7500E+01  9.7803E+01  1.6768E+02 -4.8876E+01  5.1863E+00  1.1768E+01  1.7468E+01
             1.9637E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -689.591680456722        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      241
 NPARAMETR:  1.8687E+00  5.9423E-01  3.2799E+00  1.6626E+00  3.3916E+00  8.2478E+00  1.1454E+01  9.3677E-01  1.0043E+00  1.5790E+00
             9.9040E+00
 PARAMETER:  7.2526E-01 -4.2048E-01  1.2878E+00  6.0838E-01  1.3213E+00  2.2099E+00  2.5383E+00  3.4686E-02  1.0428E-01  5.5682E-01
             2.3929E+00
 GRADIENT:   1.2594E+01  8.0963E+00 -1.5877E+01 -3.4394E+01  3.6333E+00  1.7601E+02  3.0637E+01  6.8216E-01  1.3681E+01  4.9724E+00
             2.2301E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -739.326755514310        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      417
 NPARAMETR:  2.9391E+00  4.9577E-01  8.0535E+00  1.6971E+00  8.5916E+00  5.4893E+00  1.2117E+01  9.1547E-01  1.4657E+00  6.7142E+00
             8.8709E+00
 PARAMETER:  1.1781E+00 -6.0163E-01  2.1861E+00  6.2894E-01  2.2508E+00  1.8028E+00  2.5946E+00  1.1680E-02  4.8230E-01  2.0042E+00
             2.2828E+00
 GRADIENT:   6.3342E+01  8.0723E+00  3.8119E-01 -5.2840E+01 -2.8683E+00  7.3537E+01  2.4265E+01  1.1554E-02  2.5360E+01  3.8873E+00
             1.2002E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -802.827529236175        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      594
 NPARAMETR:  1.4387E+00  3.3468E-01  7.4783E+00  2.0898E+00  4.8679E+00  3.2509E+00  1.0238E+01  1.1176E+00  1.5315E+00  2.6784E+00
             8.8697E+00
 PARAMETER:  4.6372E-01 -9.9459E-01  2.1120E+00  8.3706E-01  1.6827E+00  1.2789E+00  2.4261E+00  2.1114E-01  5.2624E-01  1.0852E+00
             2.2826E+00
 GRADIENT:   2.7024E+01  2.4274E+00 -4.8631E+00  6.8373E+01 -1.0748E+01  1.3121E+01 -9.4789E+00  1.2742E-01 -1.1250E+00  8.5906E+00
             6.0938E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -816.375984375242        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      772
 NPARAMETR:  1.2659E+00  2.5261E-01  1.0231E+01  1.7759E+00  4.8061E+00  3.1938E+00  1.1274E+01  1.4320E+00  1.2158E+00  2.3325E+00
             8.1231E+00
 PARAMETER:  3.3575E-01 -1.2759E+00  2.4254E+00  6.7433E-01  1.6699E+00  1.2612E+00  2.5225E+00  4.5909E-01  2.9538E-01  9.4695E-01
             2.1947E+00
 GRADIENT:   6.9433E+00 -3.8674E-01 -2.1253E+00  2.0818E+01 -2.9067E+00  8.2899E+00 -8.6226E+00  8.5389E-02 -1.4236E+00  8.6742E+00
            -1.2267E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -836.335381563606        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  1.2733E+00  2.4660E-01  5.4038E+02  1.4844E+00  2.1045E+00  3.2390E+00  1.0574E+01  3.9120E-02  1.2223E+00  5.0279E-01
             7.8247E+00
 PARAMETER:  3.4164E-01 -1.3000E+00  6.3923E+00  4.9499E-01  8.4409E-01  1.2753E+00  2.4584E+00 -3.1411E+00  3.0074E-01 -5.8759E-01
             2.1573E+00
 GRADIENT:   4.7763E+00 -4.0522E+00  5.3594E-02 -5.5726E+01 -1.1457E+01 -3.6311E+00 -3.2976E+01  9.6813E-08  6.3465E+00  2.5594E+00
             3.1389E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -840.682779770275        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:     1073             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2648E+00  2.7522E-01  3.8309E+02  1.5024E+00  2.1459E+00  3.2333E+00  9.4328E+00  4.8284E-02  1.2159E+00  5.4352E-02
             7.5715E+00
 PARAMETER:  3.3489E-01 -1.1902E+00  6.0483E+00  5.0705E-01  8.6358E-01  1.2735E+00  2.3442E+00 -2.9307E+00  2.9545E-01 -2.8123E+00
             2.1244E+00
 GRADIENT:   3.0757E+01  4.9550E+00  6.2361E-02 -4.8450E+00 -8.8084E+00  1.0599E+02  3.0676E+02  9.0766E-08  2.1865E+00  2.9970E-02
             8.0793E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -841.328431503855        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:     1198
 NPARAMETR:  1.2663E+00  2.7660E-01  3.2728E+02  1.4996E+00  2.3445E+00  3.3815E+00  9.1578E+00  5.0221E-02  1.1747E+00  1.0000E-02
             7.6866E+00
 PARAMETER:  3.3613E-01 -1.1852E+00  5.8908E+00  5.0519E-01  9.5205E-01  1.3183E+00  2.3146E+00 -2.8913E+00  2.6097E-01 -4.7622E+00
             2.1395E+00
 GRADIENT:   2.8711E+01  2.6430E+00  1.4895E-02 -1.3357E+01  8.4539E+00  1.2480E+02  2.9187E+02  6.3871E-08  2.0908E+00  0.0000E+00
             2.8213E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -841.577432936975        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1386
 NPARAMETR:  1.2651E+00  2.8038E-01  3.3212E+02  1.5026E+00  2.2971E+00  3.3911E+00  9.3669E+00  4.9821E-02  1.1715E+00  1.1278E-02
             7.6822E+00
 PARAMETER:  3.3518E-01 -1.1716E+00  5.9055E+00  5.0717E-01  9.3163E-01  1.3211E+00  2.3372E+00 -2.8993E+00  2.5825E-01 -4.3849E+00
             2.1389E+00
 GRADIENT:   3.5050E+00 -7.0419E+00  2.7644E-02 -2.9576E+01  3.3245E+00  1.4775E+01 -4.4026E+01  8.8069E-09  7.6434E-01  1.1335E-03
             1.5076E+00

0ITERATION NO.:   46    OBJECTIVE VALUE:  -841.577432936975        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1416
 NPARAMETR:  1.2653E+00  2.8023E-01  3.2820E+02  1.5022E+00  2.2981E+00  3.3888E+00  9.3556E+00  4.9856E-02  1.1716E+00  1.0794E-02
             7.6750E+00
 PARAMETER:  3.3518E-01 -1.1716E+00  5.9055E+00  5.0717E-01  9.3163E-01  1.3211E+00  2.3372E+00 -2.8993E+00  2.5825E-01 -4.3849E+00
             2.1389E+00
 GRADIENT:  -5.7014E+02  1.5792E+02  2.9064E-02  3.5139E+02 -2.0569E+02  1.5719E+02  4.3437E+01 -2.9259E-05 -7.4715E+02  1.1193E-03
             8.1844E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1416
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1220E-03  9.1562E-02  2.7971E-07 -7.3867E-02 -1.0863E-05
 SE:             2.9275E-02  2.0590E-02  8.4990E-07  1.5802E-02  7.3058E-05
 N:                     100         100         100         100         100

 P VAL.:         9.4221E-01  8.7150E-06  7.4208E-01  2.9508E-06  8.8180E-01

 ETASHRINKSD(%)  1.9259E+00  3.1022E+01  9.9997E+01  4.7060E+01  9.9755E+01
 ETASHRINKVR(%)  3.8147E+00  5.2421E+01  1.0000E+02  7.1974E+01  9.9999E+01
 EBVSHRINKSD(%)  1.1220E+00  2.4092E+01  9.9996E+01  4.3820E+01  9.9658E+01
 EBVSHRINKVR(%)  2.2315E+00  4.2380E+01  1.0000E+02  6.8438E+01  9.9999E+01
 RELATIVEINF(%)  9.7268E+01  3.7750E+01  3.7886E-08  1.9019E+01  2.1097E-04
 EPSSHRINKSD(%)  1.1219E+01
 EPSSHRINKVR(%)  2.1179E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -841.57743293697467     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       261.14880690863242     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.26
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -841.577       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.27E+00  2.80E-01  3.32E+02  1.50E+00  2.30E+00  3.39E+00  9.37E+00  4.98E-02  1.17E+00  1.13E-02  7.68E+00
 


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
+        2.66E+04
 
 TH 2
+        3.49E+02  4.48E+04
 
 TH 3
+       -3.58E-05  9.29E-05 -1.11E-07
 
 TH 4
+        1.48E+02  3.42E+01  9.01E-05  8.47E+03
 
 TH 5
+        8.59E+01 -1.16E+02 -4.20E-04 -6.35E+01  1.09E+03
 
 TH 6
+       -3.32E+02  4.27E+02  6.62E-06  1.85E+02 -5.40E+02  2.87E+02
 
 TH 7
+        5.29E+00  8.28E+00 -1.45E-06 -4.09E+00 -1.44E+00  6.35E+00  1.11E+01
 
 TH 8
+        4.52E-02 -2.51E-01 -2.15E-05  6.18E-03 -8.62E-03 -3.98E-02  1.62E-03  5.09E-01
 
 TH 9
+        1.45E-01 -2.63E-01 -1.05E-04 -2.82E+01  7.52E+03  1.13E-02  5.64E-01  2.93E-01  5.28E+04
 
 TH10
+       -7.90E-01  1.58E-01 -1.04E-04 -3.13E-02 -3.54E-02 -3.36E-02 -5.79E-03  1.98E-01 -3.47E-02  1.03E+01
 
 TH11
+        1.57E+02 -2.08E+02  1.02E-05 -9.91E+01 -1.23E+02  6.31E+01 -2.91E+00  1.20E-03 -8.58E+02 -9.78E-03  2.52E+01
 
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
 #CPUT: Total CPU Time in Seconds,       48.258
Stop Time:
Thu Sep 30 09:22:52 CDT 2021
