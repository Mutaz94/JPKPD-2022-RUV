Wed Sep 29 20:28:12 CDT 2021
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
$DATA ../../../../data/spa/D/dat99.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26821.4069739730        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.3732E+02  6.4615E+02  2.4293E+01  6.5240E+02  1.6676E+02 -3.6509E+03 -1.3080E+03 -6.7514E+01 -1.8498E+03 -8.7211E+02
            -4.8595E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -435.805906669306        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1892E+00  9.7921E-01  8.0359E-01  1.3425E+00  1.8049E+00  1.6411E+00  9.0957E-01  9.5679E-01  6.6438E-01  9.0051E-01
             1.5590E+01
 PARAMETER:  2.7330E-01  7.8988E-02 -1.1866E-01  3.9455E-01  6.9051E-01  5.9538E-01  5.2171E-03  5.5829E-02 -3.0890E-01 -4.7940E-03
             2.8466E+00
 GRADIENT:  -5.4583E+00  1.0188E+01 -3.8610E+00  1.5166E+01 -3.9443E+00  3.2018E+01  1.6345E+00  3.6558E+00  6.4388E+00  7.0515E-02
            -2.2548E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -449.225102956950        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.1908E+00  8.6320E-01  7.2775E-01  1.2909E+00  5.9445E+00  1.2721E+00  4.8553E-01  2.2234E-01  3.8043E-01  3.0959E+00
             1.6563E+01
 PARAMETER:  2.7461E-01 -4.7109E-02 -2.1779E-01  3.5536E-01  1.8825E+00  3.4066E-01 -6.2252E-01 -1.4035E+00 -8.6645E-01  1.2301E+00
             2.9072E+00
 GRADIENT:   1.0765E+01  1.4256E+01  3.7417E+00  1.8256E+01 -2.4878E+00 -2.9277E+01  3.8939E-01  1.3281E-01  2.9161E+00 -7.3325E-01
            -1.9593E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -464.236892554963        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1071E+00  3.9349E-01  4.9967E-01  1.4140E+00  2.0925E+01  1.3628E+00  6.1976E-01  2.8355E-02  1.2749E-01  7.6328E+00
             1.6634E+01
 PARAMETER:  2.0172E-01 -8.3271E-01 -5.9380E-01  4.4641E-01  3.1410E+00  4.0952E-01 -3.7843E-01 -3.4630E+00 -1.9597E+00  2.1325E+00
             2.9115E+00
 GRADIENT:  -1.1360E+01  7.6024E+00  1.0937E+01  1.2843E+01  8.1853E+00 -2.1424E+01  3.3456E-01  4.2658E-03  7.9112E-01 -1.6479E+01
             2.7950E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -476.511816734653        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0380E+00  1.8562E-01  4.2408E-01  1.4907E+00  2.6798E+01  1.6888E+00  7.2234E-01  1.0000E-02  7.6049E-02  1.4531E+01
             1.5736E+01
 PARAMETER:  1.3729E-01 -1.5841E+00 -7.5783E-01  4.9925E-01  3.3883E+00  6.2399E-01 -2.2526E-01 -4.7326E+00 -2.4764E+00  2.7763E+00
             2.8559E+00
 GRADIENT:  -3.9728E+01  5.1241E+00 -1.3234E+00 -1.4264E+01  1.7458E+01  2.7705E+01  4.0700E-02  0.0000E+00  4.7328E-01  2.9227E+01
             4.6665E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -477.013754350494        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  1.0373E+00  1.8335E-01  4.2430E-01  1.4917E+00  2.5235E+01  1.6945E+00  7.3405E-01  1.0000E-02  7.4819E-02  1.4462E+01
             1.5719E+01
 PARAMETER:  1.3666E-01 -1.5964E+00 -7.5731E-01  4.9994E-01  3.3282E+00  6.2736E-01 -2.0918E-01 -4.7100E+00 -2.4927E+00  2.7715E+00
             2.8549E+00
 GRADIENT:  -4.4215E+01  5.2810E+00 -3.2217E+00 -2.3211E+01  1.0252E+01  2.1118E+01  3.7812E-02  0.0000E+00  4.2890E-01 -2.3708E+01
             2.3509E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -477.973285597907        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0371E+00  1.8123E-01  4.2629E-01  1.4885E+00  2.5500E+01  1.5101E+00  5.2559E-01  1.0000E-02  2.3462E-02  1.4264E+01
             1.5749E+01
 PARAMETER:  1.3639E-01 -1.6080E+00 -7.5264E-01  4.9777E-01  3.3387E+00  5.1215E-01 -5.4323E-01 -4.7100E+00 -3.6524E+00  2.7577E+00
             2.8568E+00
 GRADIENT:  -5.2815E+01  5.3637E+00 -8.2296E+00  8.8365E+00  9.8655E+00 -7.4884E+00  1.9758E-02  0.0000E+00  3.9129E-02 -2.3136E+01
             1.1526E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -479.773123459037        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      698             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1078E+00  1.7993E-01  4.4319E-01  1.4861E+00  2.5477E+01  1.5413E+00  3.6932E-01  1.0000E-02  1.2802E-02  1.4406E+01
             1.5723E+01
 PARAMETER:  2.0238E-01 -1.6152E+00 -7.1376E-01  4.9617E-01  3.3378E+00  5.3265E-01 -8.9609E-01 -4.7100E+00 -4.2582E+00  2.7676E+00
             2.8551E+00
 GRADIENT:   6.0688E+00  2.9690E+00  3.1201E+00 -2.3716E+01  1.8850E+01  5.4704E+00  1.1847E-02  0.0000E+00  1.3955E-02  4.2239E+01
             2.4630E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -480.162834042171        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      859
 NPARAMETR:  1.1078E+00  1.7900E-01  4.3775E-01  1.4901E+00  2.5365E+01  1.5425E+00  7.7129E-02  1.0000E-02  1.0000E-02  1.4566E+01
             1.5778E+01
 PARAMETER:  2.0234E-01 -1.6204E+00 -7.2610E-01  4.9884E-01  3.3334E+00  5.3344E-01 -2.4623E+00 -4.7100E+00 -5.3943E+00  2.7787E+00
             2.8586E+00
 GRADIENT:   2.3652E+01  1.5161E+00 -1.3903E+01 -1.1530E+01  3.4499E+00 -1.6301E+01  5.6012E-04  0.0000E+00  0.0000E+00  1.6102E+01
            -2.1890E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -481.797296842545        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1052             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1066E+00  1.8176E-01  4.7133E-01  1.5019E+00  2.4853E+01  1.5435E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5041E+01
             1.6106E+01
 PARAMETER:  2.0127E-01 -1.6051E+00 -6.5220E-01  5.0674E-01  3.3130E+00  5.3406E-01 -7.2393E+00 -4.7100E+00 -8.0552E+00  2.8107E+00
             2.8792E+00
 GRADIENT:   2.1767E-02  2.7043E+00  7.9714E+00 -1.8662E+01  1.6660E+01 -4.7253E-01  0.0000E+00  0.0000E+00  0.0000E+00  6.7510E+01
             3.5992E+01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -481.797296842545        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1107
 NPARAMETR:  1.1066E+00  1.8176E-01  4.7133E-01  1.5019E+00  2.4853E+01  1.5435E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5041E+01
             1.6106E+01
 PARAMETER:  2.0127E-01 -1.6051E+00 -6.5220E-01  5.0674E-01  3.3130E+00  5.3406E-01 -7.2393E+00 -4.7100E+00 -8.0552E+00  2.8107E+00
             2.8792E+00
 GRADIENT:  -7.3593E+00  2.5166E+00  6.8310E+00 -2.8270E+01  6.9554E+00 -5.9099E+00  0.0000E+00  0.0000E+00  0.0000E+00 -9.9673E+00
             1.3349E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1107
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0816E-03 -2.1948E-05 -3.3506E-05 -4.5468E-04 -3.1420E-02
 SE:             2.7437E-02  6.8861E-06  6.8185E-05  2.0715E-04  7.3381E-03
 N:                     100         100         100         100         100

 P VAL.:         8.2458E-01  1.4361E-03  6.2315E-01  2.8168E-02  1.8551E-05

 ETASHRINKSD(%)  8.0841E+00  9.9977E+01  9.9772E+01  9.9306E+01  7.5416E+01
 ETASHRINKVR(%)  1.5515E+01  1.0000E+02  9.9999E+01  9.9995E+01  9.3956E+01
 EBVSHRINKSD(%)  7.4655E+00  9.9967E+01  9.9704E+01  9.9163E+01  7.5113E+01
 EBVSHRINKVR(%)  1.4374E+01  1.0000E+02  9.9999E+01  9.9993E+01  9.3806E+01
 RELATIVEINF(%)  2.9281E+01  1.1551E-06  1.8587E-04  6.2972E-04  5.3565E+00
 EPSSHRINKSD(%)  8.9517E-01
 EPSSHRINKVR(%)  1.7823E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -481.79729684254477     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       253.35352972119341     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -481.797       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.11E+00  1.82E-01  4.71E-01  1.50E+00  2.49E+01  1.54E+00  1.00E-02  1.00E-02  1.00E-02  1.50E+01  1.61E+01
 


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
+        2.02E+04
 
 TH 2
+       -3.28E+02  2.92E+04
 
 TH 3
+       -2.45E+01 -7.94E+02  2.72E+04
 
 TH 4
+       -2.27E+02 -9.55E+02  7.33E+02  5.53E+03
 
 TH 5
+       -6.23E+01 -5.21E-02  7.79E+01 -1.47E+01  4.43E-01
 
 TH 6
+       -8.00E+03  3.08E+03 -1.71E+02  4.74E+01 -1.16E+00  3.63E+03
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        7.67E+00 -8.63E+01  6.59E+01 -1.44E+00 -1.09E+00  4.71E-01  0.00E+00  0.00E+00  0.00E+00  3.98E+00
 
 TH11
+       -1.19E+01  2.75E+00 -6.48E+01 -8.93E-01  1.42E-01 -1.14E+00  0.00E+00  0.00E+00  0.00E+00 -1.24E+00  3.84E+00
 
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
 #CPUT: Total CPU Time in Seconds,       31.647
Stop Time:
Wed Sep 29 20:28:45 CDT 2021
