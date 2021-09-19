Sat Sep 18 11:25:50 CDT 2021
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
$DATA ../../../../data/spa/S1/dat100.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1666.66416882802        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2604E+02 -5.1700E+01 -3.0593E+01 -1.0794E+01  3.8581E+01 -1.1909E+01 -1.5589E+00 -5.1575E-01  3.2617E+01  5.8348E+00
             1.6511E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1674.43173340910        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6360E-01  1.1908E+00  9.9768E-01  8.7042E-01  1.0966E+00  1.0127E+00  1.1155E+00  1.1095E+00  6.7266E-01  9.7270E-01
             9.9331E-01
 PARAMETER:  6.2925E-02  2.7459E-01  9.7674E-02 -3.8782E-02  1.9224E-01  1.1263E-01  2.0932E-01  2.0393E-01 -2.9652E-01  7.2324E-02
             9.3284E-02
 GRADIENT:   4.3637E+01 -3.3293E+01 -1.5155E+00 -3.4112E+01  3.2443E+01 -2.3077E+00 -4.2606E+00 -5.7997E+00 -4.5529E+00 -1.1004E+01
             8.0305E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1675.49022121673        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.7486E-01  1.3163E+00  9.1793E-01  8.2705E-01  1.0951E+00  1.0426E+00  1.1476E+00  1.3037E+00  6.4157E-01  1.0494E+00
             9.7222E-01
 PARAMETER:  7.4543E-02  3.7484E-01  1.4369E-02 -8.9888E-02  1.9088E-01  1.4174E-01  2.3763E-01  3.6520E-01 -3.4384E-01  1.4822E-01
             7.1825E-02
 GRADIENT:   2.2444E+01  1.6324E+01 -7.9163E+00  1.7338E+01  2.7595E+00  2.0476E+00  1.0346E+01  1.6319E+00 -2.7136E+00  3.6006E+00
             3.0402E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1676.89270011431        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      412
 NPARAMETR:  9.6418E-01  1.3638E+00  9.7957E-01  7.9287E-01  1.1402E+00  1.0379E+00  1.0082E+00  1.4353E+00  7.8152E-01  1.0348E+00
             9.6334E-01
 PARAMETER:  6.3526E-02  4.1031E-01  7.9362E-02 -1.3209E-01  2.3123E-01  1.3722E-01  1.0818E-01  4.6135E-01 -1.4652E-01  1.3419E-01
             6.2655E-02
 GRADIENT:   6.6193E-01  1.0064E+00 -9.3373E-01  2.5130E+00  2.2933E+00  4.5555E-01  6.1161E-01  2.3469E-01 -9.3423E-01 -1.7466E+00
            -1.4786E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1677.01562241378        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  9.6430E-01  1.5198E+00  8.7729E-01  6.9229E-01  1.1795E+00  1.0378E+00  9.1767E-01  1.5076E+00  8.6537E-01  1.0693E+00
             9.6591E-01
 PARAMETER:  6.3652E-02  5.1859E-01 -3.0912E-02 -2.6776E-01  2.6510E-01  1.3713E-01  1.4087E-02  5.1050E-01 -4.4599E-02  1.6699E-01
             6.5320E-02
 GRADIENT:  -4.1249E-01  3.9227E+00  1.4450E+00  5.1477E-01 -2.3571E+00  9.6459E-02  1.5198E-01 -3.5076E-01 -5.2987E-01 -1.3546E-03
            -4.9524E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1677.04474457400        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      763
 NPARAMETR:  9.6510E-01  1.5990E+00  7.5625E-01  6.3722E-01  1.1745E+00  1.0381E+00  8.7634E-01  1.4500E+00  9.2300E-01  1.0543E+00
             9.6654E-01
 PARAMETER:  6.4480E-02  5.6936E-01 -1.7938E-01 -3.5063E-01  2.6087E-01  1.3741E-01 -3.2000E-02  4.7155E-01  1.9878E-02  1.5283E-01
             6.5972E-02
 GRADIENT:   1.0018E-01  5.8542E-01 -2.1990E-02  4.0128E-01 -1.7089E-01 -2.7041E-02  3.6539E-03  3.4165E-02  1.8522E-02  3.4733E-02
             3.4331E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1677.04487121658        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  9.6508E-01  1.6028E+00  7.5136E-01  6.3437E-01  1.1750E+00  1.0382E+00  8.7461E-01  1.4467E+00  9.2535E-01  1.0542E+00
             9.6649E-01
 PARAMETER:  6.4453E-02  5.7176E-01 -1.8587E-01 -3.5512E-01  2.6123E-01  1.3749E-01 -3.3981E-02  4.6930E-01  2.2418E-02  1.5275E-01
             6.5917E-02
 GRADIENT:   6.0344E-03  9.6699E-02 -1.8816E-02  1.0433E-01  3.6256E-02 -2.1210E-03 -1.0697E-02 -8.8263E-03  2.6389E-04 -3.5468E-03
             3.3551E-03

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1677.04487129586        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:      997
 NPARAMETR:  9.6508E-01  1.6032E+00  7.5109E-01  6.3414E-01  1.1750E+00  1.0382E+00  8.7447E-01  1.4470E+00  9.2554E-01  1.0542E+00
             9.6649E-01
 PARAMETER:  6.4453E-02  5.7197E-01 -1.8623E-01 -3.5549E-01  2.6131E-01  1.3749E-01 -3.4143E-02  4.6946E-01  2.2622E-02  1.5280E-01
             6.5916E-02
 GRADIENT:  -1.4404E-02 -1.1214E+01 -1.7093E-02  9.6477E-02 -2.9786E+05 -6.3824E-03 -1.0022E-02  1.6566E+05 -7.7835E+05 -3.6009E-03
            -2.5338E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      997
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5515E-04 -1.4341E-02 -3.5705E-02  1.2191E-02 -3.8778E-02
 SE:             2.9857E-02  2.4997E-02  1.2646E-02  1.8625E-02  2.1596E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9318E-01  5.6617E-01  4.7532E-03  5.1277E-01  7.2552E-02

 ETASHRINKSD(%)  1.0000E-10  1.6256E+01  5.7633E+01  3.7602E+01  2.7651E+01
 ETASHRINKVR(%)  1.0000E-10  2.9869E+01  8.2050E+01  6.1065E+01  4.7657E+01
 EBVSHRINKSD(%)  3.7727E-01  1.6102E+01  6.1797E+01  4.0488E+01  2.3859E+01
 EBVSHRINKVR(%)  7.5312E-01  2.9611E+01  8.5405E+01  6.4584E+01  4.2026E+01
 RELATIVEINF(%)  9.8941E+01  1.9293E+00  7.1988E-01  7.9822E-01  1.4380E+01
 EPSSHRINKSD(%)  4.4916E+01
 EPSSHRINKVR(%)  6.9658E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1677.0448712958628     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -941.89404473212460     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.58
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1677.045       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.65E-01  1.60E+00  7.51E-01  6.34E-01  1.18E+00  1.04E+00  8.74E-01  1.45E+00  9.26E-01  1.05E+00  9.66E-01
 


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
+        4.18E+09
 
 TH 2
+       -6.51E+00  4.10E+02
 
 TH 3
+        7.21E+00  5.53E+01  1.60E+02
 
 TH 4
+        1.27E+03 -9.41E+07 -6.17E+08  1.15E+03
 
 TH 5
+        1.31E+09  1.38E+08 -4.53E+08  2.81E+08  4.30E+02
 
 TH 6
+        2.60E+00 -1.84E+00  9.35E-01 -3.85E+00 -9.14E+02  1.82E+02
 
 TH 7
+        2.26E+00  2.34E+01 -1.22E+00 -2.55E+01  7.25E+08 -2.34E-01  1.35E+02
 
 TH 8
+       -2.97E+08 -3.12E+07 -2.05E+08  2.20E+01 -1.87E+08  4.12E+02  2.66E+03  7.24E+00
 
 TH 9
+        4.36E+09  2.29E+08 -1.68E+01 -9.33E+08  1.37E+09 -3.03E+03  2.40E+09  4.80E+00  4.54E+09
 
 TH10
+        1.25E+09  2.63E+08 -8.64E+08  5.36E+08 -6.46E+01  5.99E-01 -1.74E+00 -3.56E+08  2.61E+09  7.06E+01
 
 TH11
+       -4.17E+09 -2.20E+08 -1.03E+01  6.11E+00 -1.31E+09  2.68E+00  8.16E+00  3.10E+00 -4.35E+09 -2.50E+09  4.17E+09
 
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
 #CPUT: Total CPU Time in Seconds,       18.994
Stop Time:
Sat Sep 18 11:26:10 CDT 2021
