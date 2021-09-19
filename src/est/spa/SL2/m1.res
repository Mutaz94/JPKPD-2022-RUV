Sat Sep 18 12:01:34 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1657.53967399359        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1901E+01 -7.8989E+00  3.5330E+01 -4.7239E+01 -5.1892E+01  6.9517E+00 -4.3691E-01 -2.2682E+00  7.5075E+00 -2.3261E+01
            -6.1587E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1662.50797536286        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.9803E-01  1.0551E+00  1.0033E+00  1.0053E+00  1.0943E+00  9.9537E-01  9.9954E-01  9.6997E-01  9.4989E-01  1.2190E+00
             9.8513E-01
 PARAMETER:  9.8028E-02  1.5359E-01  1.0333E-01  1.0528E-01  1.9015E-01  9.5360E-02  9.9545E-02  6.9508E-02  4.8591E-02  2.9807E-01
             8.5021E-02
 GRADIENT:   3.7884E+01  1.0570E+01 -3.2220E+00  1.6013E+01  2.2490E+00  5.3668E+00  1.8414E+00 -3.4570E-01 -4.7946E+00 -2.7337E+00
            -8.8035E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1663.41931946165        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      254
 NPARAMETR:  1.0022E+00  1.1423E+00  1.0629E+00  9.5301E-01  1.1809E+00  9.8862E-01  7.8407E-01  1.1139E+00  1.0814E+00  1.3410E+00
             9.9030E-01
 PARAMETER:  1.0223E-01  2.3303E-01  1.6102E-01  5.1873E-02  2.6629E-01  8.8554E-02 -1.4326E-01  2.0786E-01  1.7827E-01  3.9345E-01
             9.0251E-02
 GRADIENT:   3.1317E+00  2.2477E+00 -2.8944E+00  1.2298E+01 -3.2602E-01 -2.0570E+00 -1.0738E+00 -1.7159E+00  2.3362E+00  1.6096E+00
            -4.9992E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.30465677981        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      430
 NPARAMETR:  1.0002E+00  1.2793E+00  1.1334E+00  8.6094E-01  1.2908E+00  9.9353E-01  7.2899E-01  1.4363E+00  1.1602E+00  1.3992E+00
             1.0007E+00
 PARAMETER:  1.0019E-01  3.4634E-01  2.2521E-01 -4.9734E-02  3.5523E-01  9.3505E-02 -2.1610E-01  4.6208E-01  2.4863E-01  4.3588E-01
             1.0075E-01
 GRADIENT:  -1.7654E+00 -9.2988E-01  1.0010E+00  1.4324E+00 -2.7126E-01  6.0216E-02 -5.9732E-02 -8.1030E-01 -2.7798E-01  2.2134E-01
             2.3120E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1665.02598247025        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  1.0036E+00  1.6103E+00  9.0751E-01  6.4498E-01  1.3767E+00  9.9474E-01  6.6470E-01  1.6339E+00  1.4406E+00  1.4208E+00
             9.9798E-01
 PARAMETER:  1.0361E-01  5.7645E-01  2.9470E-03 -3.3853E-01  4.1966E-01  9.4731E-02 -3.0842E-01  5.9097E-01  4.6505E-01  4.5119E-01
             9.7974E-02
 GRADIENT:   2.4084E+00  1.3509E+01  3.2702E+00  5.9482E+00 -6.0765E+00 -1.9829E-01 -7.6873E-02  3.4085E-01  9.5365E-03 -3.5619E-01
            -9.3774E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1665.60160008203        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      785
 NPARAMETR:  1.0038E+00  1.7983E+00  6.0348E-01  5.1314E-01  1.3683E+00  9.9724E-01  6.8248E-01  1.1756E+00  1.6335E+00  1.3732E+00
             1.0013E+00
 PARAMETER:  1.0376E-01  6.8686E-01 -4.0503E-01 -5.6720E-01  4.1354E-01  9.7239E-02 -2.8202E-01  2.6178E-01  5.9071E-01  4.1716E-01
             1.0126E-01
 GRADIENT:  -1.1315E-01  6.8715E+00 -1.2377E+00  5.6483E+00 -4.0138E-01  2.0580E-02  1.8800E+00  3.7023E-01  1.9800E-02  1.3943E-01
             2.4724E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.76243936808        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      961
 NPARAMETR:  1.0041E+00  1.8578E+00  5.3958E-01  4.6978E-01  1.3832E+00  9.9751E-01  6.6676E-01  1.0663E+00  1.7244E+00  1.3740E+00
             1.0011E+00
 PARAMETER:  1.0407E-01  7.1937E-01 -5.1696E-01 -6.5548E-01  4.2443E-01  9.7508E-02 -3.0532E-01  1.6415E-01  6.4489E-01  4.1774E-01
             1.0111E-01
 GRADIENT:   1.8017E-01  2.5798E+00 -1.8618E+00  4.4547E+00  1.8988E+00  2.5878E-02 -3.1308E-01  1.9440E-01 -4.4352E-01 -2.9555E-02
             1.8445E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.77148699761        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1160             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0040E+00  1.8609E+00  5.3632E-01  4.6737E-01  1.3844E+00  9.9748E-01  6.6641E-01  1.0176E+00  1.7304E+00  1.3740E+00
             1.0011E+00
 PARAMETER:  1.0402E-01  7.2108E-01 -5.2302E-01 -6.6063E-01  4.2523E-01  9.7481E-02 -3.0586E-01  1.1743E-01  6.4832E-01  4.1775E-01
             1.0110E-01
 GRADIENT:   4.5908E+01  9.1693E+01 -6.4383E-01  1.4265E+01  4.1203E+00  4.9637E+00  1.0321E+00 -1.2050E-01  1.7146E+00  2.5409E-01
            -7.2947E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1665.77257041837        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1317
 NPARAMETR:  1.0040E+00  1.8609E+00  5.3632E-01  4.6737E-01  1.3844E+00  9.9744E-01  6.6641E-01  1.0348E+00  1.7304E+00  1.3740E+00
             1.0013E+00
 PARAMETER:  1.0399E-01  7.2108E-01 -5.2302E-01 -6.6063E-01  4.2523E-01  9.7437E-02 -3.0586E-01  1.3421E-01  6.4832E-01  4.1775E-01
             1.0130E-01
 GRADIENT:   6.2138E-04  3.1676E+00 -1.3246E+00  4.1421E+00  2.0580E+00 -1.5317E-03 -4.9725E-01  1.6642E-05 -5.6561E-01 -1.7458E-01
            -8.8298E-05

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1665.79975969810        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1498
 NPARAMETR:  1.0040E+00  1.8605E+00  5.4042E-01  4.6338E-01  1.3828E+00  9.9739E-01  6.6751E-01  1.0535E+00  1.7404E+00  1.3741E+00
             1.0011E+00
 PARAMETER:  1.0394E-01  7.2083E-01 -5.1540E-01 -6.6921E-01  4.2411E-01  9.7382E-02 -3.0421E-01  1.5215E-01  6.5414E-01  4.1777E-01
             1.0107E-01
 GRADIENT:   9.1489E-02 -2.3068E+00  1.0652E-01  1.2668E-01 -1.1247E+00  1.3859E-02  2.1180E-01 -5.7694E-02 -4.3613E-01 -6.4023E-03
             2.6277E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1665.79999868800        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     1567
 NPARAMETR:  1.0039E+00  1.8605E+00  5.4045E-01  4.6334E-01  1.3828E+00  9.9738E-01  6.6726E-01  1.0540E+00  1.7406E+00  1.3741E+00
             1.0011E+00
 PARAMETER:  1.0394E-01  7.2086E-01 -5.1538E-01 -6.6926E-01  4.2417E-01  9.7381E-02 -3.0427E-01  1.5260E-01  6.5424E-01  4.1779E-01
             1.0106E-01
 GRADIENT:   7.0033E-02  9.8419E+04 -1.3764E+05  5.3001E+04  1.6728E+05  1.1210E-02  1.9530E-01  4.6487E+05  1.0839E+05 -1.6983E+05
            -7.0207E+05
 NUMSIGDIG:         3.5         3.3         3.3         3.3         3.3         3.5         2.1         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1567
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0954E-04 -4.6381E-02 -2.4026E-02  3.5927E-02 -5.1861E-02
 SE:             2.9879E-02  2.1578E-02  8.1382E-03  2.3002E-02  2.2495E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9173E-01  3.1594E-02  3.1547E-03  1.1832E-01  2.1142E-02

 ETASHRINKSD(%)  1.0000E-10  2.7712E+01  7.2736E+01  2.2939E+01  2.4639E+01
 ETASHRINKVR(%)  1.0000E-10  4.7745E+01  9.2567E+01  4.0616E+01  4.3207E+01
 EBVSHRINKSD(%)  4.3552E-01  2.5913E+01  7.5639E+01  2.5083E+01  2.0927E+01
 EBVSHRINKVR(%)  8.6914E-01  4.5111E+01  9.4065E+01  4.3874E+01  3.7474E+01
 RELATIVEINF(%)  9.9075E+01  4.2444E+00  9.2197E-01  4.6975E+00  2.4140E+01
 EPSSHRINKSD(%)  4.4960E+01
 EPSSHRINKVR(%)  6.9706E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.7999986879970     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.64917212425883     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.800       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.86E+00  5.40E-01  4.63E-01  1.38E+00  9.97E-01  6.67E-01  1.05E+00  1.74E+00  1.37E+00  1.00E+00
 


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
+        1.63E+09
 
 TH 2
+       -1.27E+08  9.86E+06
 
 TH 3
+       -3.47E+03  6.12E+03  2.29E+08
 
 TH 4
+        3.11E+03 -4.99E+03 -2.05E+08  1.84E+08
 
 TH 5
+        1.65E+03 -2.93E+03  4.31E+04 -2.53E+04  5.16E+07
 
 TH 6
+       -1.16E+00  1.03E+03 -4.95E+03  4.44E+03  2.35E+03  1.99E+02
 
 TH 7
+       -8.37E+08  6.51E+07  1.61E+03 -1.46E+03 -7.79E+02 -8.76E+08  4.30E+08
 
 TH 8
+       -1.06E+09  8.22E+07  1.33E+05 -1.19E+05 -6.31E+04 -1.11E+09  5.43E+08  6.85E+08
 
 TH 9
+        8.52E+02 -1.48E+03 -5.59E+07  5.02E+07 -2.75E+04  1.21E+03 -3.76E+02 -3.25E+04  1.37E+07
 
 TH10
+        2.96E+08 -2.30E+07 -1.68E+03  1.49E+03  7.51E+02 -2.40E+03 -1.52E+08 -1.92E+08  4.08E+02  5.38E+07
 
 TH11
+        1.68E+09 -1.31E+08  1.70E+04 -1.53E+04 -8.06E+03 -1.36E+04 -8.63E+08 -1.09E+09 -4.15E+03  3.05E+08  1.73E+09
 
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
 #CPUT: Total CPU Time in Seconds,       27.127
Stop Time:
Sat Sep 18 12:02:03 CDT 2021
