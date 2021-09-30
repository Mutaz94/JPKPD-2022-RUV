Thu Sep 30 00:21:06 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat53.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   691.382306825879        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4113E+02  1.2969E+00  2.0856E+02 -5.6466E+01  2.1259E+02  3.4516E+01 -2.6572E+01 -2.6631E+02 -1.4729E+01 -9.8732E+01
            -5.0723E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1425.55692640703        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8346E-01  1.0693E+00  9.5909E-01  1.1255E+00  1.0085E+00  8.1851E-01  8.4551E-01  9.6936E-01  7.4308E-01  8.6089E-01
             5.3398E+00
 PARAMETER:  8.3318E-02  1.6698E-01  5.8227E-02  2.1826E-01  1.0850E-01 -1.0027E-01 -6.7811E-02  6.8876E-02 -1.9695E-01 -4.9787E-02
             1.7752E+00
 GRADIENT:  -8.2257E+01 -7.2919E+00 -1.8645E+01  5.4881E+00  1.4845E+00 -4.1843E+01  1.1543E+01  6.6698E+00  1.6059E+01  1.6053E+01
             3.3446E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1467.03940496765        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.7555E-01  9.1351E-01  2.5828E-01  1.1079E+00  4.2422E-01  9.6373E-01  2.3977E-01  1.1793E-01  1.1441E+00  2.0759E-01
             4.5746E+00
 PARAMETER:  7.5251E-02  9.5412E-03 -1.2537E+00  2.0246E-01 -7.5751E-01  6.3051E-02 -1.3281E+00 -2.0376E+00  2.3461E-01 -1.4722E+00
             1.6205E+00
 GRADIENT:  -7.1958E+01  9.5443E+01  1.1414E+01  6.2587E+01 -8.3598E+01 -9.6332E+00 -6.2744E-01  1.6701E-01  2.7686E+01  5.8209E-01
             2.2094E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1509.07588783912        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      252
 NPARAMETR:  9.8013E-01  7.1675E-01  3.1359E-01  1.1547E+00  4.1516E-01  9.9156E-01  6.8023E-01  5.1593E-02  8.5296E-01  1.1091E-01
             3.4970E+00
 PARAMETER:  7.9930E-02 -2.3303E-01 -1.0597E+00  2.4383E-01 -7.7908E-01  9.1522E-02 -2.8532E-01 -2.8644E+00 -5.9044E-02 -2.0990E+00
             1.3519E+00
 GRADIENT:  -2.8982E+01  1.2318E+01 -9.7905E+00  2.2092E+01 -2.2682E+00 -8.9204E-01 -3.2042E+00 -8.7199E-03 -6.8294E+00 -1.9094E-01
            -3.1470E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1510.64702687598        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      430
 NPARAMETR:  9.9364E-01  6.6714E-01  3.6893E-01  1.1929E+00  4.3784E-01  9.8426E-01  7.3865E-01  5.7983E-02  8.3871E-01  8.4457E-02
             3.6006E+00
 PARAMETER:  9.3616E-02 -3.0476E-01 -8.9714E-01  2.7640E-01 -7.2590E-01  8.4135E-02 -2.0293E-01 -2.7476E+00 -7.5896E-02 -2.3715E+00
             1.3811E+00
 GRADIENT:   1.9776E-01  2.0608E+00  8.5250E-02  2.4320E+00 -6.1016E-01 -1.7596E-01 -5.8296E-01  6.5331E-03  6.7982E-01 -3.1326E-03
            -2.4076E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1510.73308515127        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.9314E-01  6.1219E-01  3.6999E-01  1.2144E+00  4.2225E-01  9.8385E-01  8.0575E-01  5.6380E-02  8.2262E-01  6.3530E-02
             3.6065E+00
 PARAMETER:  9.3117E-02 -3.9072E-01 -8.9427E-01  2.9425E-01 -7.6215E-01  8.3715E-02 -1.1598E-01 -2.7756E+00 -9.5259E-02 -2.6562E+00
             1.3827E+00
 GRADIENT:   7.9048E-02 -3.5471E-01 -2.1325E-01 -2.9734E-01  6.8914E-01 -5.6591E-02 -6.8950E-02  4.7665E-03  1.6117E-02 -2.2184E-02
            -4.6118E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1510.73952123575        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  9.9323E-01  6.1258E-01  3.6517E-01  1.2119E+00  4.1903E-01  9.8441E-01  8.1002E-01  2.0904E-02  8.2486E-01  7.0135E-02
             3.6055E+00
 PARAMETER:  9.3211E-02 -3.9007E-01 -9.0739E-01  2.9221E-01 -7.6981E-01  8.4289E-02 -1.1070E-01 -3.7678E+00 -9.2544E-02 -2.5573E+00
             1.3825E+00
 GRADIENT:   3.1529E-02 -3.6159E-01 -6.1575E-01 -1.9654E-01  1.0319E+00  9.5336E-04  3.9740E-02  5.4305E-04  3.3003E-02 -2.9560E-02
             1.1195E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1510.78035020114        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      961
 NPARAMETR:  9.9337E-01  6.0984E-01  3.6613E-01  1.2132E+00  4.1799E-01  9.8465E-01  7.7852E-01  1.0000E-02  8.2803E-01  1.7323E-01
             3.5951E+00
 PARAMETER:  9.3351E-02 -3.9456E-01 -9.0478E-01  2.9327E-01 -7.7230E-01  8.4527E-02 -1.5036E-01 -6.1495E+00 -8.8711E-02 -1.6531E+00
             1.3796E+00
 GRADIENT:   1.7450E-01  5.3515E-01  5.3967E-01 -3.9326E-01  9.0695E-01  2.7568E-02  2.7266E-03  0.0000E+00  5.8170E-03  1.0369E-02
             9.4633E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1510.78035020114        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      983
 NPARAMETR:  9.9337E-01  6.0984E-01  3.6613E-01  1.2132E+00  4.1799E-01  9.8465E-01  7.7852E-01  1.0000E-02  8.2803E-01  1.7323E-01
             3.5951E+00
 PARAMETER:  9.3351E-02 -3.9456E-01 -9.0478E-01  2.9327E-01 -7.7230E-01  8.4527E-02 -1.5036E-01 -6.1495E+00 -8.8711E-02 -1.6531E+00
             1.3796E+00
 GRADIENT:   1.7450E-01  5.3515E-01  5.3967E-01 -3.9326E-01  9.0695E-01  2.7568E-02  2.7266E-03  0.0000E+00  5.8170E-03  1.0369E-02
             9.4633E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      983
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7430E-03 -8.6335E-04  2.9655E-05 -1.4039E-02  2.5029E-04
 SE:             2.8822E-02  1.2516E-02  1.9674E-04  2.3975E-02  5.7401E-03
 N:                     100         100         100         100         100

 P VAL.:         9.5178E-01  9.4501E-01  8.8019E-01  5.5817E-01  9.6522E-01

 ETASHRINKSD(%)  3.4437E+00  5.8068E+01  9.9341E+01  1.9682E+01  8.0770E+01
 ETASHRINKVR(%)  6.7689E+00  8.2417E+01  9.9996E+01  3.5490E+01  9.6302E+01
 EBVSHRINKSD(%)  3.3837E+00  5.8596E+01  9.9324E+01  1.9655E+01  8.0856E+01
 EBVSHRINKVR(%)  6.6530E+00  8.2857E+01  9.9995E+01  3.5446E+01  9.6335E+01
 RELATIVEINF(%)  9.2574E+01  7.8757E-01  1.8314E-04  1.5749E+01  6.9084E-02
 EPSSHRINKSD(%)  1.8219E+01
 EPSSHRINKVR(%)  3.3119E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1510.7803502011432     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -591.84181699647047     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.17
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1510.780       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  6.10E-01  3.66E-01  1.21E+00  4.18E-01  9.85E-01  7.79E-01  1.00E-02  8.28E-01  1.73E-01  3.60E+00
 


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
+        1.09E+03
 
 TH 2
+       -7.36E+01  7.77E+02
 
 TH 3
+       -4.66E+01  1.33E+03  3.72E+03
 
 TH 4
+       -5.23E+01  3.26E+02 -2.28E+02  6.94E+02
 
 TH 5
+        1.47E+02 -2.20E+03 -5.12E+03 -1.75E+02  7.57E+03
 
 TH 6
+        2.03E+00 -8.20E+00  2.14E+01 -1.67E+01  3.23E+00  1.78E+02
 
 TH 7
+       -8.89E-01 -1.56E+01 -4.52E+01 -2.78E+00  7.17E+01  9.07E-01  1.16E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.50E+00 -2.74E+01  2.14E+01 -2.87E+00  5.97E+01 -7.38E-02  1.13E+01  0.00E+00  1.17E+02
 
 TH10
+       -2.92E+00 -2.19E+01 -6.89E+01 -3.31E+00  1.22E+02  2.21E-01  6.19E+00  0.00E+00  7.94E-01  1.33E+01
 
 TH11
+       -1.57E+01 -1.05E+01 -2.65E+01 -1.58E+01  2.66E+01  2.97E+00  6.41E+00  0.00E+00  1.14E+01  9.18E+00  4.50E+01
 
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
 #CPUT: Total CPU Time in Seconds,       23.480
Stop Time:
Thu Sep 30 00:21:31 CDT 2021
