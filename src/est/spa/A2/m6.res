Wed Sep 29 12:33:07 CDT 2021
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
$DATA ../../../../data/spa/A2/dat6.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -834.326360370397        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3911E+02 -6.1060E+00 -3.3351E+01  4.9190E+01  2.4346E+02  4.9618E+01 -3.4116E+01 -3.9230E+00 -6.1847E+01 -7.8452E+01
            -1.4585E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1363.02272488528        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.5812E-01  9.4362E-01  1.0250E+00  1.0130E+00  8.9139E-01  8.2604E-01  1.0075E+00  9.5567E-01  1.1129E+00  9.2327E-01
             2.4362E+00
 PARAMETER:  5.7218E-02  4.1967E-02  1.2472E-01  1.1296E-01 -1.4969E-02 -9.1116E-02  1.0746E-01  5.4658E-02  2.0695E-01  2.0169E-02
             9.9044E-01
 GRADIENT:  -7.8038E+01 -4.0361E+01 -2.0542E+01 -4.6828E+01  6.1124E+01 -5.1299E+01  1.1565E-01  5.7284E+00 -5.1976E+00 -5.8710E-02
            -6.1853E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1373.24636144700        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.6093E-01  6.4931E-01  6.6977E-01  1.2223E+00  5.7370E-01  8.5169E-01  7.6819E-01  5.4315E-01  1.0815E+00  5.1100E-01
             2.4384E+00
 PARAMETER:  6.0143E-02 -3.3184E-01 -3.0082E-01  3.0071E-01 -4.5565E-01 -6.0536E-02 -1.6371E-01 -5.1037E-01  1.7831E-01 -5.7138E-01
             9.9133E-01
 GRADIENT:  -8.0291E+01  3.6104E+01  1.6222E+01  7.5951E+01 -3.1634E+00 -4.2447E+01 -6.2548E+00  3.8445E-01  1.5197E+00 -1.1805E+01
            -7.2455E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1387.83213409288        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      316
 NPARAMETR:  1.0143E+00  6.0289E-01  6.3284E-01  1.2459E+00  5.5305E-01  9.2588E-01  1.1079E+00  1.1594E-01  1.0376E+00  5.8512E-01
             2.7023E+00
 PARAMETER:  1.1418E-01 -4.0602E-01 -3.5754E-01  3.1989E-01 -4.9230E-01  2.2985E-02  2.0244E-01 -2.0547E+00  1.3686E-01 -4.3593E-01
             1.0941E+00
 GRADIENT:  -3.4334E+00  4.7536E+00 -1.5211E+01  1.7471E+01  2.4845E+01 -7.4799E+00 -1.5921E+00  1.2902E-01  4.6161E+00 -1.4652E+00
            -4.4359E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1391.62710427480        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  1.0135E+00  3.8579E-01  4.8527E-01  1.2956E+00  4.1039E-01  9.4489E-01  1.6781E+00  1.0000E-02  9.3738E-01  5.8425E-01
             2.6112E+00
 PARAMETER:  1.1344E-01 -8.5246E-01 -6.2306E-01  3.5897E-01 -7.9064E-01  4.3316E-02  6.1768E-01 -5.2279E+00  3.5332E-02 -4.3743E-01
             1.0598E+00
 GRADIENT:  -7.8955E+00  9.4887E+00  5.0201E+00 -2.0811E+00 -1.3228E+01 -3.0266E+00  2.0827E+00  0.0000E+00 -5.4048E+00  2.1916E+00
            -3.9705E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1395.09174001998        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      668
 NPARAMETR:  1.0088E+00  1.7423E-01  5.2889E-01  1.4140E+00  4.0443E-01  9.4634E-01  2.2443E+00  1.0000E-02  9.1806E-01  5.7134E-01
             2.6598E+00
 PARAMETER:  1.0876E-01 -1.6474E+00 -5.3697E-01  4.4643E-01 -8.0527E-01  4.4850E-02  9.0838E-01 -9.7982E+00  1.4507E-02 -4.5977E-01
             1.0782E+00
 GRADIENT:  -9.2598E-01  1.5317E+00 -2.7587E+00  3.9419E+00  3.6813E+00 -1.4666E-01 -1.2664E+00  0.0000E+00 -6.3471E-01 -1.5251E+00
            -1.2338E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1396.40977019320        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  1.0049E+00  4.7216E-02  5.2101E-01  1.4648E+00  3.8298E-01  9.4190E-01  5.4856E+00  1.0000E-02  9.0641E-01  5.8460E-01
             2.6666E+00
 PARAMETER:  1.0487E-01 -2.9530E+00 -5.5199E-01  4.8174E-01 -8.5977E-01  4.0147E-02  1.8021E+00 -1.9116E+01  1.7337E-03 -4.3683E-01
             1.0808E+00
 GRADIENT:   1.4616E+00  5.5267E-01  3.0517E+00  2.4550E+00 -5.0098E+00 -9.9786E-01  1.6757E-01  0.0000E+00  2.9835E+00  2.4982E-01
             4.0519E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1396.68873862446        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1018
 NPARAMETR:  1.0026E+00  1.1069E-02  5.2648E-01  1.4809E+00  3.8246E-01  9.4556E-01  1.1348E+01  1.0000E-02  8.9002E-01  5.9323E-01
             2.6442E+00
 PARAMETER:  1.0257E-01 -4.4036E+00 -5.4155E-01  4.9264E-01 -8.6114E-01  4.4019E-02  2.5290E+00 -2.9450E+01 -1.6516E-02 -4.2217E-01
             1.0724E+00
 GRADIENT:   7.9778E-01  1.0135E-01 -9.5112E-01 -5.5424E-01  1.4242E+00  4.6285E-01  3.9606E-02  0.0000E+00 -3.9495E-01  2.2736E-01
            -6.3172E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1396.69635314556        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1081
 NPARAMETR:  1.0024E+00  1.0000E-02  5.2711E-01  1.4802E+00  3.8278E-01  9.4389E-01  1.1596E+01  1.0000E-02  8.9073E-01  5.9092E-01
             2.6488E+00
 PARAMETER:  1.0211E-01 -4.5056E+00 -5.3936E-01  4.9338E-01 -8.6075E-01  4.2896E-02  2.5764E+00 -3.0165E+01 -1.5613E-02 -4.2825E-01
             1.0743E+00
 GRADIENT:  -4.1958E-01  2.5293E-02  7.5269E-01  1.5948E+00 -6.6722E-01  1.3103E-01  1.7859E-02  0.0000E+00  1.5695E-02 -8.6939E-02
             6.0050E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1081
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1736E-04  1.1297E-03  7.0356E-05 -1.0829E-02 -6.3009E-03
 SE:             2.8985E-02  1.9709E-03  2.2363E-04  2.7072E-02  1.8859E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9677E-01  5.6651E-01  7.5306E-01  6.8914E-01  7.3829E-01

 ETASHRINKSD(%)  2.8950E+00  9.3397E+01  9.9251E+01  9.3069E+00  3.6822E+01
 ETASHRINKVR(%)  5.7061E+00  9.9564E+01  9.9994E+01  1.7748E+01  6.0085E+01
 EBVSHRINKSD(%)  2.7847E+00  9.4027E+01  9.9228E+01  8.7569E+00  3.6503E+01
 EBVSHRINKVR(%)  5.4918E+00  9.9643E+01  9.9994E+01  1.6747E+01  5.9681E+01
 RELATIVEINF(%)  8.1828E+01  1.9710E-02  2.4118E-04  8.5121E+00  1.3723E+00
 EPSSHRINKSD(%)  3.0647E+01
 EPSSHRINKVR(%)  5.1902E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1396.6963531455647     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -661.54552658182649     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1396.696       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  5.28E-01  1.48E+00  3.83E-01  9.44E-01  1.19E+01  1.00E-02  8.91E-01  5.90E-01  2.65E+00
 


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
+        1.17E+03
 
 TH 2
+       -4.06E+01  3.20E+05
 
 TH 3
+       -1.79E+01  2.12E+02  2.33E+03
 
 TH 4
+       -4.23E+01  1.11E+02 -1.51E+02  5.33E+02
 
 TH 5
+        1.18E+02 -5.27E+02 -4.27E+03 -2.26E+02  8.66E+03
 
 TH 6
+       -9.80E-01 -4.21E+00  1.17E+01 -1.09E+01  3.08E+00  1.97E+02
 
 TH 7
+        8.20E-03  1.34E+03 -1.14E-02 -4.48E-02  7.92E-02  1.49E-04  1.96E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.03E+01 -2.44E+01  3.73E+01 -1.31E+01  1.24E+01  2.12E+00 -1.30E-02  0.00E+00  1.77E+02
 
 TH10
+       -1.09E+01 -5.27E+00 -5.33E+01 -3.30E+00  5.37E+01 -1.51E+00 -6.60E-03  0.00E+00  1.23E+00  9.88E+01
 
 TH11
+       -1.48E+01 -5.58E+00 -6.30E+00 -6.36E+00 -7.33E+00  3.94E+00 -1.02E-02  0.00E+00  7.78E+00  2.74E+01  4.40E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.845
Stop Time:
Wed Sep 29 12:33:28 CDT 2021
