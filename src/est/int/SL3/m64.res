Sat Sep 25 02:32:06 CDT 2021
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
$DATA ../../../../data/int/SL3/dat64.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      984
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

 TOT. NO. OF OBS RECS:      884
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   271.230689817861        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.7226E+01 -8.4715E+01 -1.7310E+01  2.1074E+02  1.9866E+02 -3.0186E+01 -1.5533E+02 -2.0982E+02 -2.4367E+02 -3.9836E+01
            -7.5880E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2344.45507464740        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1003E+00  1.5717E+00  1.0548E+00  8.2527E-01  1.1288E+00  1.1535E+00  1.1065E+00  1.0602E+00  1.0720E+00  1.0677E+00
             5.1504E+00
 PARAMETER:  1.9554E-01  5.5216E-01  1.5332E-01 -9.2049E-02  2.2116E-01  2.4276E-01  2.0121E-01  1.5842E-01  1.6952E-01  1.6546E-01
             1.7391E+00
 GRADIENT:   3.3316E+01  7.8837E+01  3.7653E+00  2.7820E+01 -4.6277E+01  3.0397E+01  3.3989E+01  2.4280E+00  1.4808E+01 -8.6331E-01
             7.6404E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2425.65405641480        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.8917E-01  2.1880E+00  5.2098E+00  6.3862E-01  2.0432E+00  1.4502E+00  1.0378E+00  3.8298E+00  2.0948E+00  1.5785E+00
             4.2791E+00
 PARAMETER:  8.9114E-02  8.8297E-01  1.7505E+00 -3.4844E-01  8.1452E-01  4.7167E-01  1.3706E-01  1.4428E+00  8.3946E-01  5.5646E-01
             1.5537E+00
 GRADIENT:  -6.9368E+01  1.9880E+02 -4.7431E+00  9.8307E+01  7.1234E+01  8.0213E+01  3.6684E+01  3.7027E+00  2.5768E+01 -3.3107E+01
             6.1346E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2641.90321852692        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0270E+00  1.7373E+00  4.1298E+00  5.9899E-01  1.7845E+00  1.0644E+00  8.6658E-01  3.9473E+00  1.0022E+00  1.8712E+00
             2.7906E+00
 PARAMETER:  1.2662E-01  6.5235E-01  1.5182E+00 -4.1251E-01  6.7915E-01  1.6239E-01 -4.3196E-02  1.4730E+00  1.0217E-01  7.2656E-01
             1.1263E+00
 GRADIENT:   7.0096E-01 -8.4452E+00 -5.9872E+00  6.6287E-01  2.8416E+01 -1.1308E+00  7.6571E+00 -4.1735E+00 -1.3126E+00  3.3419E+00
            -2.6409E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2644.34420450504        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      295
 NPARAMETR:  1.0286E+00  1.7051E+00  5.8185E+00  6.2322E-01  1.6974E+00  1.0629E+00  7.1101E-01  4.3730E+00  1.3299E+00  1.8329E+00
             2.8142E+00
 PARAMETER:  1.2821E-01  6.3362E-01  1.8610E+00 -3.7285E-01  6.2911E-01  1.6105E-01 -2.4106E-01  1.5755E+00  3.8512E-01  7.0588E-01
             1.1347E+00
 GRADIENT:   3.3854E+00  4.1255E-01  3.6352E-01 -1.5648E+00 -4.8438E+00 -1.5244E+00  4.6991E+00  1.1032E+00  1.5660E+00 -1.9548E+00
            -1.9176E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2645.28682012494        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      365
 NPARAMETR:  1.0259E+00  1.6176E+00  9.6681E+00  6.8074E-01  1.7518E+00  1.0703E+00  5.7356E-01  4.6850E+00  1.4020E+00  1.8977E+00
             2.8137E+00
 PARAMETER:  1.2554E-01  5.8096E-01  2.3688E+00 -2.8457E-01  6.6065E-01  1.6797E-01 -4.5589E-01  1.6444E+00  4.3787E-01  7.4063E-01
             1.1345E+00
 GRADIENT:  -2.0108E+00  2.8652E+00  6.6552E-01 -9.6036E-01  2.0448E-01  9.1061E-01  1.0676E+00 -6.9013E-01 -2.0083E-01  5.7947E-01
             9.2120E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2645.55499781415        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      543
 NPARAMETR:  1.0343E+00  1.7196E+00  8.5950E+00  6.2271E-01  1.7833E+00  1.0757E+00  5.7361E-01  4.7813E+00  1.4996E+00  1.9241E+00
             2.8152E+00
 PARAMETER:  1.3373E-01  6.4207E-01  2.2512E+00 -3.7368E-01  6.7849E-01  1.7297E-01 -4.5581E-01  1.6647E+00  5.0520E-01  7.5445E-01
             1.1350E+00
 GRADIENT:   4.7138E+00  9.3049E+00 -1.5897E+00  9.2602E+00  2.1048E+00  1.0950E+00 -5.2175E-01  9.8357E-01  6.4295E-01  5.4902E-01
            -2.4353E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2646.15755940236        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      719
 NPARAMETR:  1.0310E+00  1.9005E+00  7.8230E+00  4.8699E-01  1.8009E+00  1.0721E+00  5.7552E-01  5.1899E+00  1.7465E+00  1.9525E+00
             2.8141E+00
 PARAMETER:  1.3053E-01  7.4210E-01  2.1571E+00 -6.1950E-01  6.8830E-01  1.6960E-01 -4.5247E-01  1.7467E+00  6.5762E-01  7.6913E-01
             1.1347E+00
 GRADIENT:  -1.1336E+00  1.0821E-01  1.2056E+00 -1.9907E+00  4.9673E-01 -1.8717E-01 -7.6356E-01 -3.3526E+00 -1.0847E-01  2.6839E+00
             9.4077E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2646.19946024401        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  1.0312E+00  1.8922E+00  7.8426E+00  4.9401E-01  1.8046E+00  1.0725E+00  5.8266E-01  5.2117E+00  1.7184E+00  1.9451E+00
             2.8153E+00
 PARAMETER:  1.3060E-01  7.3771E-01  2.1594E+00 -6.0516E-01  6.9028E-01  1.6998E-01 -4.4013E-01  1.7510E+00  6.4134E-01  7.6526E-01
             1.1350E+00
 GRADIENT:  -9.7653E-01 -6.9248E+01 -4.6931E+01  1.5800E+02 -1.3871E+02 -4.7745E-02  2.1206E+02  5.9151E+01 -7.2407E+01 -1.2920E+02
            -8.5427E+01
 NUMSIGDIG:         2.1         3.3         3.3         3.3         3.3         2.8         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      862
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2560E-03 -4.2899E-02 -3.0380E-02  2.6839E-02 -2.4618E-02
 SE:             2.9402E-02  1.7497E-02  1.1299E-02  2.2599E-02  2.5885E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3884E-01  1.4216E-02  7.1732E-03  2.3498E-01  3.4159E-01

 ETASHRINKSD(%)  1.4992E+00  4.1382E+01  6.2146E+01  2.4291E+01  1.3281E+01
 ETASHRINKVR(%)  2.9760E+00  6.5639E+01  8.5671E+01  4.2681E+01  2.4799E+01
 EBVSHRINKSD(%)  1.5998E+00  3.8837E+01  7.2608E+01  2.8310E+01  8.5309E+00
 EBVSHRINKVR(%)  3.1741E+00  6.2591E+01  9.2497E+01  4.8605E+01  1.6334E+01
 RELATIVEINF(%)  9.6734E+01  3.5462E+00  5.1344E+00  4.8757E+00  6.7555E+01
 EPSSHRINKSD(%)  1.6420E+01
 EPSSHRINKVR(%)  3.0143E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          884
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1624.6833267058612     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2646.1994602440141     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1021.5161335381529     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2646.199       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.89E+00  7.84E+00  4.94E-01  1.80E+00  1.07E+00  5.83E-01  5.21E+00  1.72E+00  1.95E+00  2.82E+00
 


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
+        8.76E+02
 
 TH 2
+       -1.47E+02  1.39E+04
 
 TH 3
+       -9.21E+00  1.07E+02  9.32E+01
 
 TH 4
+        1.12E+03 -1.73E+03 -1.52E+02  2.78E+05
 
 TH 5
+       -2.60E+02  6.01E+02  4.20E+01 -6.68E+04  1.62E+04
 
 TH 6
+        4.69E+00 -2.92E+01 -1.70E+00  2.05E+02 -4.88E+01  1.60E+02
 
 TH 7
+        2.50E+03 -1.04E+03 -6.68E+01  8.89E+03 -2.00E+03  4.65E+02  3.63E+05
 
 TH 8
+        1.25E+01 -1.33E+02 -1.36E+01  2.00E+02 -5.80E+01  2.19E+00  8.86E+01  3.43E+02
 
 TH 9
+       -7.30E+02  1.39E+02  9.41E+00 -1.40E+03  3.23E+02 -1.36E+02 -8.40E+04 -1.07E+01  1.95E+04
 
 TH10
+       -1.47E+02  1.12E+03  7.30E+01 -2.51E+03  6.98E+02 -2.70E+01 -1.10E+03 -9.68E+01  1.69E+02  1.17E+04
 
 TH11
+       -9.42E+01  3.53E+02  2.52E+01 -1.47E+03  4.10E+02 -1.30E+01 -6.07E+02 -3.36E+01  9.96E+01  4.40E+02  2.64E+03
 
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
 #CPUT: Total CPU Time in Seconds,       38.389
Stop Time:
Sat Sep 25 02:32:46 CDT 2021
