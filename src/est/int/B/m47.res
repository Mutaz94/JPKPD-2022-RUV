Fri Sep 24 19:32:00 CDT 2021
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
$DATA ../../../../data/int/B/dat47.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3496.33812187649        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.6270E+01 -9.1468E+01  1.3213E+02 -1.4023E+02  3.4404E+01  2.5928E+01  1.3139E+00 -4.9011E+02 -1.2657E+02 -3.8957E+00
            -1.4289E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3695.76405989703        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0344E+00  1.0115E+00  9.0874E-01  9.8165E-01  9.8346E-01  9.4042E-01  9.6738E-01  2.5276E+00  1.1359E+00  1.0124E+00
             9.4978E-01
 PARAMETER:  1.3386E-01  1.1148E-01  4.3002E-03  8.1475E-02  8.3323E-02  3.8573E-02  6.6834E-02  1.0273E+00  2.2744E-01  1.1233E-01
             4.8472E-02
 GRADIENT:   1.8612E+02 -7.0410E+01 -2.2336E+01 -1.1473E+02 -3.3851E+01 -3.5346E+00  1.2232E+00 -1.6073E+01  2.5117E+01 -5.3680E+00
            -2.0373E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3698.07767669111        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.0339E+00  1.0468E+00  9.0963E-01  9.6485E-01  1.0025E+00  9.4339E-01  8.8647E-01  2.5291E+00  1.1567E+00  1.0682E+00
             9.5424E-01
 PARAMETER:  1.3333E-01  1.4572E-01  5.2772E-03  6.4216E-02  1.0249E-01  4.1727E-02 -2.0506E-02  1.0279E+00  2.4556E-01  1.6596E-01
             5.3163E-02
 GRADIENT:   1.8315E+02 -6.2791E+01 -2.2849E+01 -1.1422E+02 -3.4789E+01 -2.0870E+00  4.4353E-02 -1.6499E+01  2.5187E+01 -4.9144E+00
            -1.9841E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3701.57315777167        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:      249
 NPARAMETR:  1.0298E+00  1.0562E+00  9.1900E-01  9.6688E-01  1.0115E+00  9.4144E-01  8.8195E-01  2.5482E+00  1.1437E+00  1.0763E+00
             9.5967E-01
 PARAMETER:  1.2940E-01  1.5466E-01  1.5530E-02  6.6319E-02  1.1144E-01  3.9651E-02 -2.5619E-02  1.0354E+00  2.3430E-01  1.7352E-01
             5.8834E-02
 GRADIENT:   1.7133E+02 -5.8668E+01 -2.1496E+01 -1.0691E+02 -3.2439E+01 -1.9941E+00  7.0376E-02 -1.5085E+01  2.3417E+01 -4.3129E+00
            -1.8621E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3702.56352662395        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.0290E+00  1.0607E+00  9.2226E-01  9.6640E-01  1.0155E+00  9.4124E-01  8.7917E-01  2.5534E+00  1.1413E+00  1.0799E+00
             9.6134E-01
 PARAMETER:  1.2863E-01  1.5895E-01  1.9070E-02  6.5821E-02  1.1537E-01  3.9442E-02 -2.8774E-02  1.0374E+00  2.3219E-01  1.7683E-01
             6.0569E-02
 GRADIENT:   1.1214E+02 -7.2584E+01 -2.1812E+01 -1.1332E+02 -3.7912E+01 -5.8142E+00 -4.5407E-01 -1.8270E+01  2.1024E+01 -5.2462E+00
            -1.8305E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3705.44351304785        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.0290E+00  1.0662E+00  9.2255E-01  9.7409E-01  1.0159E+00  9.4216E-01  8.7871E-01  2.5538E+00  1.1294E+00  1.0835E+00
             9.6916E-01
 PARAMETER:  1.2856E-01  1.6414E-01  1.9381E-02  7.3752E-02  1.1576E-01  4.0415E-02 -2.9298E-02  1.0376E+00  2.2173E-01  1.8019E-01
             6.8671E-02
 GRADIENT:   1.6736E+02 -4.7478E+01 -2.1478E+01 -9.2241E+01 -3.4217E+01 -1.5218E+00  6.9336E-01 -1.4821E+01  2.1968E+01 -3.5324E+00
            -1.6425E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3705.93826159685        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  1.0290E+00  1.0668E+00  9.2255E-01  9.7549E-01  1.0159E+00  9.4246E-01  8.7821E-01  2.5538E+00  1.1260E+00  1.0842E+00
             9.7068E-01
 PARAMETER:  1.2856E-01  1.6462E-01  1.9381E-02  7.5188E-02  1.1576E-01  4.0743E-02 -2.9872E-02  1.0376E+00  2.1867E-01  1.8088E-01
             7.0243E-02
 GRADIENT:   1.1148E+02 -6.2272E+01 -2.2283E+01 -9.8788E+01 -4.0620E+01 -5.2358E+00  1.7836E-01 -1.8327E+01  1.9635E+01 -4.4440E+00
            -1.6081E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3705.96681841251        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      974
 NPARAMETR:  1.0290E+00  1.0667E+00  9.2257E-01  9.7552E-01  1.0159E+00  9.5485E-01  8.7823E-01  2.5538E+00  1.1261E+00  1.0843E+00
             9.7066E-01
 PARAMETER:  1.2856E-01  1.6458E-01  1.9406E-02  7.5214E-02  1.1573E-01  5.3801E-02 -2.9847E-02  1.0376E+00  2.1872E-01  1.8092E-01
             7.0218E-02
 GRADIENT:   1.0864E+02 -6.2294E+01 -2.2268E+01 -9.8783E+01 -4.0617E+01 -9.0360E-03  1.7927E-01 -1.8332E+01  1.9649E+01 -4.4265E+00
            -1.6080E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3705.98675302749        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1164             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0289E+00  1.0667E+00  9.2261E-01  9.7561E-01  1.0158E+00  9.5424E-01  8.7826E-01  2.5546E+00  1.1261E+00  1.0844E+00
             9.7068E-01
 PARAMETER:  1.2853E-01  1.6458E-01  1.9450E-02  7.5303E-02  1.1571E-01  5.3164E-02 -2.9817E-02  1.0379E+00  2.1876E-01  1.8099E-01
             7.0237E-02
 GRADIENT:   1.6431E+02 -4.6368E+01 -2.1540E+01 -9.0294E+01 -3.4449E+01  3.5881E+00  7.0913E-01 -1.4776E+01  2.1478E+01 -3.3456E+00
            -1.6055E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3706.01680967191        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1349            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0289E+00  1.0668E+00  9.2261E-01  9.7561E-01  1.0158E+00  9.5464E-01  8.7825E-01  2.5546E+00  1.1259E+00  1.0844E+00
             9.7080E-01
 PARAMETER:  1.2853E-01  1.6468E-01  1.9450E-02  7.5303E-02  1.1571E-01  5.3576E-02 -2.9820E-02  1.0379E+00  2.1861E-01  1.8103E-01
             7.0367E-02
 GRADIENT:   1.6420E+02 -4.6265E+01 -2.1543E+01 -9.0264E+01 -3.4537E+01  3.7504E+00  7.1684E-01 -1.4776E+01  2.1440E+01 -3.3413E+00
            -1.6026E+02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -3706.01680967191        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  1.0289E+00  1.0668E+00  9.2262E-01  9.7562E-01  1.0158E+00  9.5473E-01  8.7826E-01  2.5546E+00  1.1260E+00  1.0844E+00
             9.7079E-01
 PARAMETER:  1.2853E-01  1.6468E-01  1.9450E-02  7.5303E-02  1.1571E-01  5.3576E-02 -2.9820E-02  1.0379E+00  2.1861E-01  1.8103E-01
             7.0367E-02
 GRADIENT:   1.0928E+02  1.6799E+05 -1.3840E+05 -1.3849E+05  2.3914E+05 -8.4484E-02 -2.7677E+05 -1.0983E+02 -1.2658E+05 -7.6447E+04
             2.7660E+05
 NUMSIGDIG:         6.9         3.3         3.3         3.3         3.3         2.3         3.3         6.0         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1424
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.6261E-02 -2.4563E-02 -4.3100E-03  7.0169E-02 -3.3602E-02
 SE:             2.9560E-02  2.0506E-02  2.7317E-02  2.5600E-02  2.4644E-02
 N:                     100         100         100         100         100

 P VAL.:         1.1759E-01  2.3098E-01  8.7463E-01  6.1254E-03  1.7273E-01

 ETASHRINKSD(%)  9.7124E-01  3.1303E+01  8.4852E+00  1.4237E+01  1.7440E+01
 ETASHRINKVR(%)  1.9330E+00  5.2807E+01  1.6251E+01  2.6448E+01  3.1838E+01
 EBVSHRINKSD(%)  2.6073E-01  3.0942E+01  1.3504E+01  7.4679E+00  1.9703E+01
 EBVSHRINKVR(%)  5.2078E-01  5.2309E+01  2.5184E+01  1.4378E+01  3.5524E+01
 RELATIVEINF(%)  9.9477E+01  2.2092E+01  6.7623E+01  5.7542E+01  3.1665E+01
 EPSSHRINKSD(%)  1.7883E+01
 EPSSHRINKVR(%)  3.2568E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3706.0168096719053     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2051.9274499034946     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.52
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3706.017       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.07E+00  9.23E-01  9.76E-01  1.02E+00  9.55E-01  8.78E-01  2.55E+00  1.13E+00  1.08E+00  9.71E-01
 


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
+        3.96E+08
 
 TH 2
+        5.24E+00  2.24E+08
 
 TH 3
+       -3.34E+04 -4.27E+08  1.63E+09
 
 TH 4
+       -2.10E+00  3.94E+04  7.69E+08  1.45E+09
 
 TH 5
+        4.45E+08 -3.28E+04 -6.38E+08 -6.03E+08  1.00E+09
 
 TH 6
+       -1.23E+02  1.15E+03 -2.18E+03 -2.06E+03  1.70E+03  2.12E+02
 
 TH 7
+        5.01E+00 -4.48E+08  3.16E+04  2.99E+04 -2.48E+04 -2.29E+03  8.97E+08
 
 TH 8
+        6.35E+01 -1.49E+07 -6.96E+00  2.68E+07 -1.77E+01 -3.36E-01  2.98E+07  1.97E+06
 
 TH 9
+        1.82E+00  1.55E+04  3.05E+08  2.48E+03 -1.96E+03 -8.16E+02  1.19E+04  1.06E+07  1.14E+08
 
 TH10
+        6.29E-01 -2.01E+08  3.82E+08  7.58E+03 -6.27E+03 -1.03E+03  4.01E+08  1.33E+07  1.43E+08  1.80E+08
 
 TH11
+       -1.03E+01  4.06E+08 -7.73E+08  1.39E+04 -1.16E+04  2.07E+03 -3.00E+04 -2.69E+07 -2.90E+08 -3.63E+08  7.34E+08
 
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
 #CPUT: Total CPU Time in Seconds,       58.480
Stop Time:
Fri Sep 24 19:33:02 CDT 2021
