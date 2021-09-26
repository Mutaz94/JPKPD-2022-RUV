Sat Sep 25 11:00:25 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1679.73600467666        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3354E-01  1.4161E+00 -4.4033E+01  5.4844E+01 -1.6507E-01  1.6068E+01  4.7366E+00  1.6078E+01  2.5735E+01  1.4393E+01
            -6.9577E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1694.78928128467        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0186E+00  1.0872E+00  1.2486E+00  9.1098E-01  1.1303E+00  9.5518E-01  9.7220E-01  8.9449E-01  8.6819E-01  9.4761E-01
             1.1602E+00
 PARAMETER:  1.1847E-01  1.8358E-01  3.2201E-01  6.7629E-03  2.2252E-01  5.4144E-02  7.1809E-02 -1.1505E-02 -4.1346E-02  4.6188E-02
             2.4857E-01
 GRADIENT:   3.6781E+01 -2.3529E+01 -1.8978E+00 -3.6673E+01 -6.9639E+00 -7.7700E-02 -2.0886E+00  4.1911E+00 -5.8397E+00 -8.4804E+00
            -3.7851E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1697.24231345158        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0087E+00  8.8439E-01  1.8030E+00  1.0721E+00  1.2037E+00  9.5130E-01  6.9063E-01  6.4971E-01  9.9900E-01  1.2567E+00
             1.1679E+00
 PARAMETER:  1.0868E-01 -2.2858E-02  6.8945E-01  1.6958E-01  2.8540E-01  5.0073E-02 -2.7014E-01 -3.3123E-01  9.8997E-02  3.2847E-01
             2.5520E-01
 GRADIENT:   1.3971E+01 -4.8925E-01  2.9238E+00  1.1271E+01 -2.1476E+01 -6.5933E-01  3.8584E+00  5.9122E-01  1.7092E+01  1.6826E+01
             5.4090E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.73733395051        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0039E+00  1.0301E+00  1.6986E+00  9.7136E-01  1.2466E+00  9.5322E-01  5.4561E-01  6.2305E-01  1.0390E+00  1.1525E+00
             1.1686E+00
 PARAMETER:  1.0390E-01  1.2966E-01  6.2980E-01  7.0943E-02  3.2045E-01  5.2086E-02 -5.0585E-01 -3.7313E-01  1.3828E-01  2.4197E-01
             2.5578E-01
 GRADIENT:  -1.9605E+00 -5.1135E-02 -1.8665E+00  1.0104E+00  2.4388E+00 -5.8865E-01  5.7117E-01  3.3321E-01  3.8291E-01 -8.2066E-01
            -1.6220E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1699.10645544916        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      343
 NPARAMETR:  1.0148E+00  1.0334E+00  1.7760E+00  9.7028E-01  1.2651E+00  9.6278E-01  4.7242E-01  5.7410E-01  1.0551E+00  1.1734E+00
             1.1731E+00
 PARAMETER:  1.1467E-01  1.3285E-01  6.7434E-01  6.9833E-02  3.3513E-01  6.2070E-02 -6.4988E-01 -4.5495E-01  1.5363E-01  2.5994E-01
             2.5961E-01
 GRADIENT:  -1.2218E+01 -2.3970E+00 -6.4410E-01 -2.3750E+00  8.6468E-01 -3.3271E-01 -3.4768E-01  1.4012E-01 -7.7108E-01 -4.0400E-01
             3.7855E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1699.34186791677        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  1.0189E+00  1.2204E+00  1.6445E+00  8.4587E-01  1.3007E+00  9.6287E-01  5.0537E-01  2.9175E-01  1.1722E+00  1.1835E+00
             1.1758E+00
 PARAMETER:  1.1876E-01  2.9914E-01  5.9742E-01 -6.7390E-02  3.6290E-01  6.2158E-02 -5.8247E-01 -1.1318E+00  2.5889E-01  2.6844E-01
             2.6195E-01
 GRADIENT:  -5.2774E+00 -1.5235E+00  1.4147E+00 -2.2140E+00 -3.1210E+00 -7.7614E-01  1.1133E-02  3.2270E-02 -1.3380E-01  2.5217E-01
             8.7378E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1699.47203074903        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      695
 NPARAMETR:  1.0219E+00  1.3846E+00  1.4594E+00  7.3934E-01  1.3285E+00  9.6573E-01  5.2411E-01  8.2280E-02  1.2950E+00  1.1805E+00
             1.1755E+00
 PARAMETER:  1.2170E-01  4.2540E-01  4.7802E-01 -2.0199E-01  3.8407E-01  6.5133E-02 -5.4605E-01 -2.3976E+00  3.5849E-01  2.6596E-01
             2.6173E-01
 GRADIENT:  -7.3345E-01  3.3000E+00 -3.4716E-01  3.0980E+00  2.8412E-01 -1.4301E-01 -5.3336E-03  4.9076E-03 -2.2699E-01  2.4727E-02
            -5.1492E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1699.49645558396        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      872
 NPARAMETR:  1.0227E+00  1.4652E+00  1.4040E+00  6.8236E-01  1.3477E+00  9.6657E-01  5.1705E-01  2.9558E-02  1.3788E+00  1.1839E+00
             1.1765E+00
 PARAMETER:  1.2241E-01  4.8197E-01  4.3932E-01 -2.8220E-01  3.9839E-01  6.6003E-02 -5.5961E-01 -3.4214E+00  4.2122E-01  2.6885E-01
             2.6254E-01
 GRADIENT:   9.4380E-02  7.3410E-02  6.4571E-02  3.9131E-02 -1.4937E-01  3.0649E-02  1.1964E-02  6.8378E-04  3.8029E-02  1.0876E-02
            -2.3417E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1699.49665214624        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1052
 NPARAMETR:  1.0226E+00  1.4697E+00  1.3989E+00  6.7927E-01  1.3485E+00  9.6658E-01  5.1761E-01  1.9059E-02  1.3832E+00  1.1839E+00
             1.1766E+00
 PARAMETER:  1.2238E-01  4.8505E-01  4.3565E-01 -2.8673E-01  3.9901E-01  6.6011E-02 -5.5854E-01 -3.8602E+00  4.2439E-01  2.6885E-01
             2.6260E-01
 GRADIENT:  -2.2350E-02 -3.0138E-02  2.3523E-02 -8.1158E-03 -7.5469E-02  2.1801E-02  2.6853E-02  2.8814E-04  4.9319E-02  2.5841E-02
             1.5967E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1699.49676910123        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1179
 NPARAMETR:  1.0226E+00  1.4698E+00  1.3983E+00  6.7923E-01  1.3485E+00  9.6653E-01  5.1748E-01  1.0000E-02  1.3830E+00  1.1838E+00
             1.1765E+00
 PARAMETER:  1.2239E-01  4.8510E-01  4.3529E-01 -2.8680E-01  3.9902E-01  6.5953E-02 -5.5879E-01 -4.6159E+00  4.2423E-01  2.6871E-01
             2.6258E-01
 GRADIENT:  -1.9642E-03 -1.3559E-02 -5.9539E-03 -5.2559E-03  1.9202E-02 -1.0981E-03  3.7450E-03  0.0000E+00 -5.5620E-03 -1.7778E-03
             2.2380E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1179
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.2891E-04 -3.2091E-02 -6.1363E-05  9.1956E-03 -3.4927E-02
 SE:             2.9788E-02  1.4370E-02  5.5757E-05  2.5212E-02  2.3486E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8583E-01  2.5536E-02  2.7109E-01  7.1531E-01  1.3697E-01

 ETASHRINKSD(%)  2.0708E-01  5.1859E+01  9.9813E+01  1.5536E+01  2.1320E+01
 ETASHRINKVR(%)  4.1374E-01  7.6824E+01  1.0000E+02  2.8658E+01  3.8095E+01
 EBVSHRINKSD(%)  5.8656E-01  5.1962E+01  9.9802E+01  1.5306E+01  1.9385E+01
 EBVSHRINKVR(%)  1.1697E+00  7.6923E+01  1.0000E+02  2.8270E+01  3.5012E+01
 RELATIVEINF(%)  9.8656E+01  1.5215E+00  1.4495E-04  5.1645E+00  2.2412E+01
 EPSSHRINKSD(%)  3.9547E+01
 EPSSHRINKVR(%)  6.3455E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1699.4967691012291     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -964.34594253749094     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.60
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.00
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1699.497       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.47E+00  1.40E+00  6.79E-01  1.35E+00  9.67E-01  5.17E-01  1.00E-02  1.38E+00  1.18E+00  1.18E+00
 


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
+        1.12E+03
 
 TH 2
+       -1.57E+01  4.62E+02
 
 TH 3
+        3.04E+00  2.33E+01  2.21E+01
 
 TH 4
+       -2.05E+01  5.66E+02 -1.31E+01  9.12E+02
 
 TH 5
+       -4.12E+00 -1.13E+02 -5.66E+01 -2.13E+00  2.88E+02
 
 TH 6
+        2.71E-02 -1.52E+00  1.01E+00 -3.65E+00 -6.50E-01  2.01E+02
 
 TH 7
+        1.52E+00 -4.79E+01  8.10E+00 -3.15E+01 -9.31E+00  1.08E+00  6.06E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.23E+00 -3.21E+01 -6.62E-01  3.92E+01  1.19E+00 -7.59E-01  3.48E+01  0.00E+00  5.63E+01
 
 TH10
+       -1.12E+00  1.84E+00 -1.88E+00 -6.00E+00 -4.06E+01  1.72E+00  6.82E+00  0.00E+00  2.80E-01  6.28E+01
 
 TH11
+       -9.79E+00 -2.77E+01 -1.05E+01 -1.53E+01  4.55E+00  1.79E+00  5.61E+00  0.00E+00  5.65E+00  2.27E+01  1.76E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.660
Stop Time:
Sat Sep 25 11:00:46 CDT 2021
