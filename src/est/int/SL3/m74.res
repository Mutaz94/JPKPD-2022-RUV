Sat Sep 25 02:37:57 CDT 2021
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
$DATA ../../../../data/int/SL3/dat74.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      981
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      881
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -200.537165511476        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0732E+01 -2.7731E+01  1.5011E+02 -8.9187E+00  8.6800E+01  1.9476E+01 -8.2166E+01 -2.7013E+02 -9.3555E+01 -4.4707E+01
            -6.7297E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2647.51429587137        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0488E+00  1.3381E+00  9.9298E-01  8.7812E-01  1.1770E+00  8.5072E-01  1.1189E+00  9.1755E-01  1.0728E+00  1.0075E+00
             2.9566E+00
 PARAMETER:  1.4766E-01  3.9125E-01  9.2956E-02 -2.9974E-02  2.6294E-01 -6.1669E-02  2.1237E-01  1.3953E-02  1.7026E-01  1.0744E-01
             1.1840E+00
 GRADIENT:   7.7986E+01  1.2087E+01 -9.0411E+00  2.3535E+00 -3.2306E+00 -3.9778E+01  2.0391E+01  6.3902E+00 -5.3560E+00 -1.7349E+01
             1.2600E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2655.99474006330        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0435E+00  1.4930E+00  1.1886E+00  8.2138E-01  1.3449E+00  9.0249E-01  8.5728E-01  2.2798E-01  1.1807E+00  1.2835E+00
             2.9931E+00
 PARAMETER:  1.4255E-01  5.0079E-01  2.7274E-01 -9.6775E-02  3.9630E-01 -2.5996E-03 -5.3986E-02 -1.3785E+00  2.6611E-01  3.4963E-01
             1.1963E+00
 GRADIENT:   5.4206E+01  4.2118E+01  7.2264E+00  4.1478E+01 -1.0747E+01 -1.4650E+01  7.5207E-01  1.6670E-01 -4.5404E+00 -2.4257E+00
             1.5809E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2670.17139157281        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0197E+00  1.7779E+00  1.0334E+00  6.0760E-01  1.5679E+00  9.2239E-01  7.0339E-01  1.9237E-01  1.5145E+00  1.4923E+00
             2.7513E+00
 PARAMETER:  1.1955E-01  6.7545E-01  1.3289E-01 -3.9825E-01  5.4971E-01  1.9208E-02 -2.5184E-01 -1.5483E+00  5.1507E-01  5.0033E-01
             1.1121E+00
 GRADIENT:   4.8351E+00  2.4370E+01 -3.2207E-01  2.2708E+01  4.5428E+00 -7.0611E+00 -3.9392E+00  8.8267E-02 -2.2100E+00  3.5830E+00
             1.4209E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2674.84757707415        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      327
 NPARAMETR:  1.0198E+00  2.1008E+00  8.2376E-01  3.7481E-01  1.8146E+00  9.3579E-01  6.9998E-01  2.9987E-02  2.0728E+00  1.5813E+00
             2.7418E+00
 PARAMETER:  1.1964E-01  8.4233E-01 -9.3870E-02 -8.8133E-01  6.9587E-01  3.3639E-02 -2.5670E-01 -3.4070E+00  8.2890E-01  5.5825E-01
             1.1086E+00
 GRADIENT:  -7.6595E-01 -1.3156E+01 -2.5044E-02  4.0819E+00  5.2737E+00 -2.2858E+00  3.4257E+00  1.3322E-03  3.5622E+00 -6.8324E+00
             7.0736E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2677.45512317943        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      502
 NPARAMETR:  1.0202E+00  2.4611E+00  4.8515E-01  1.5314E-01  2.0763E+00  9.4085E-01  6.8654E-01  1.0000E-02  3.0335E+00  1.7990E+00
             2.7230E+00
 PARAMETER:  1.2000E-01  1.0006E+00 -6.2330E-01 -1.7764E+00  8.3059E-01  3.9027E-02 -2.7608E-01 -8.2371E+00  1.2097E+00  6.8721E-01
             1.1017E+00
 GRADIENT:  -4.3171E-01  2.8583E+01 -7.3151E-01  4.6347E+00  4.3653E+00 -5.9162E-01 -2.2339E+00  0.0000E+00  8.5592E-01 -1.0088E+00
            -1.7045E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2677.95517143570        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  1.0201E+00  2.5271E+00  3.7339E-01  9.4140E-02  2.1361E+00  9.4240E-01  6.8869E-01  1.0000E-02  3.5988E+00  1.8584E+00
             2.7222E+00
 PARAMETER:  1.1993E-01  1.0271E+00 -8.8514E-01 -2.2630E+00  8.5900E-01  4.0675E-02 -2.7297E-01 -1.1303E+01  1.3806E+00  7.1972E-01
             1.1014E+00
 GRADIENT:  -2.5123E-01  5.0367E-01 -1.3802E+00 -5.9963E-01  4.8507E+00 -8.2115E-02  1.3322E+00  0.0000E+00 -3.4678E+00  4.4767E-01
             1.3418E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2678.10320061903        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      853
 NPARAMETR:  1.0201E+00  2.5214E+00  4.2569E-01  9.7781E-02  2.1195E+00  9.4253E-01  6.8830E-01  1.0000E-02  3.6448E+00  1.8539E+00
             2.7201E+00
 PARAMETER:  1.1986E-01  1.0248E+00 -7.5403E-01 -2.2250E+00  8.5116E-01  4.0808E-02 -2.7353E-01 -1.0877E+01  1.3933E+00  7.1731E-01
             1.1007E+00
 GRADIENT:  -1.4207E-01 -1.9717E-01 -9.5274E-02  1.5616E-01  3.6479E-03  4.8186E-02 -1.4096E-01  0.0000E+00  1.9723E-01  6.0969E-02
            -6.5467E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2678.10355580067        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1016
 NPARAMETR:  1.0201E+00  2.5216E+00  4.2733E-01  9.7583E-02  2.1197E+00  9.4242E-01  6.8829E-01  1.0000E-02  3.6444E+00  1.8539E+00
             2.7202E+00
 PARAMETER:  1.1993E-01  1.0249E+00 -7.5019E-01 -2.2271E+00  8.5125E-01  4.0696E-02 -2.7354E-01 -1.0881E+01  1.3932E+00  7.1728E-01
             1.1007E+00
 GRADIENT:  -5.0197E-03 -1.2435E-01 -1.2787E-03  8.5275E-04 -1.2525E-02  2.1310E-04 -1.0815E-03  0.0000E+00 -1.8308E-02 -1.6539E-02
            -3.7130E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1016
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3157E-03 -2.8482E-02 -1.7148E-05  3.8015E-02 -2.5089E-02
 SE:             2.9303E-02  2.5649E-02  1.7536E-05  1.5888E-02  2.6049E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6419E-01  2.6679E-01  3.2814E-01  1.6726E-02  3.3548E-01

 ETASHRINKSD(%)  1.8320E+00  1.4073E+01  9.9941E+01  4.6773E+01  1.2732E+01
 ETASHRINKVR(%)  3.6305E+00  2.6166E+01  1.0000E+02  7.1669E+01  2.3843E+01
 EBVSHRINKSD(%)  1.8871E+00  1.0422E+01  9.9928E+01  6.2168E+01  8.8855E+00
 EBVSHRINKVR(%)  3.7386E+00  1.9758E+01  1.0000E+02  8.5688E+01  1.6981E+01
 RELATIVEINF(%)  9.6188E+01  2.6214E+01  4.1160E-05  4.5053E+00  6.7340E+01
 EPSSHRINKSD(%)  1.6309E+01
 EPSSHRINKVR(%)  2.9959E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          881
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1619.1696955066332     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2678.1035558006747     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1058.9338602940416     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2678.104       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.52E+00  4.27E-01  9.76E-02  2.12E+00  9.42E-01  6.88E-01  1.00E-02  3.64E+00  1.85E+00  2.72E+00
 


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
+        1.16E+03
 
 TH 2
+       -1.52E+01  3.16E+02
 
 TH 3
+        1.09E+00  1.33E+01  4.31E+01
 
 TH 4
+       -9.88E+00  3.08E+02 -1.46E+02  2.47E+03
 
 TH 5
+       -2.66E+00 -1.16E+01 -5.43E+00  3.60E+01  5.99E+01
 
 TH 6
+        6.24E+00 -5.81E+00 -5.19E+00  8.91E+00 -1.05E+00  2.24E+02
 
 TH 7
+        5.32E+00  2.29E+01  1.75E+01 -2.64E+02  2.79E+00 -2.26E+00  3.50E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.39E-01 -7.54E+00 -3.63E+00  8.45E+01 -1.03E+00  4.40E-01 -1.29E+01  0.00E+00  5.82E+00
 
 TH10
+        8.22E-02 -4.49E+00 -4.13E+00  5.01E+01 -5.75E+00 -1.38E+00 -3.22E+00  0.00E+00  9.46E-01  3.83E+01
 
 TH11
+       -1.73E+01 -1.45E+01 -2.42E-01 -7.41E+00  4.42E-01  2.73E+00  1.01E+01  0.00E+00  7.32E-01  3.59E+00  1.55E+02
 
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
 #CPUT: Total CPU Time in Seconds,       38.626
Stop Time:
Sat Sep 25 02:38:37 CDT 2021
