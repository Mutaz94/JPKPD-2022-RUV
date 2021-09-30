Wed Sep 29 00:27:05 CDT 2021
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
$DATA ../../../../data/int/A3/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2012.52760179535        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2019E+02  1.9019E+02  1.9325E+02  2.8556E+01  7.4356E+01  5.2815E+01 -6.9167E+01 -1.1242E+02  3.0297E+01 -2.9681E+01
            -3.4225E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2999.53030791042        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0483E+00  9.4993E-01  8.0643E-01  1.0729E+00  9.3878E-01  9.2140E-01  1.0601E+00  9.8742E-01  7.5767E-01  1.0321E+00
             2.3701E+00
 PARAMETER:  1.4716E-01  4.8633E-02 -1.1513E-01  1.7039E-01  3.6826E-02  1.8141E-02  1.5836E-01  8.7344E-02 -1.7750E-01  1.3156E-01
             9.6292E-01
 GRADIENT:   2.0228E+02  2.4804E+01 -3.4153E+01  5.6242E+01  3.2421E+01 -1.1679E+01 -3.1816E-01  1.4021E+01 -2.8442E+00  4.8169E+00
             2.7211E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3009.24631939550        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      184
 NPARAMETR:  1.0491E+00  6.3157E-01  5.9242E-01  1.2775E+00  6.1661E-01  9.7539E-01  1.3400E+00  6.5431E-01  7.2579E-01  8.2163E-01
             2.3039E+00
 PARAMETER:  1.4797E-01 -3.5954E-01 -4.2353E-01  3.4494E-01 -3.8353E-01  7.5084E-02  3.9265E-01 -3.2417E-01 -2.2050E-01 -9.6460E-02
             9.3458E-01
 GRADIENT:   8.2988E+01  3.4409E+01 -1.1574E+01  2.2192E+02  1.5804E+01  1.4215E+00 -2.3133E+00  1.2126E+01 -2.8190E+01  1.2125E+01
             2.3468E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3034.05636085937        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      362
 NPARAMETR:  1.0203E+00  5.6892E-01  5.1245E-01  1.2686E+00  5.5157E-01  9.4088E-01  1.1703E+00  3.5574E-01  8.7057E-01  7.5225E-01
             2.0812E+00
 PARAMETER:  1.2007E-01 -4.6402E-01 -5.6856E-01  3.3789E-01 -4.9499E-01  3.9061E-02  2.5729E-01 -9.3355E-01 -3.8603E-02 -1.8469E-01
             8.3293E-01
 GRADIENT:   2.9798E+01  2.0873E+01 -7.9576E+00  1.3941E+02  3.1821E+01 -1.1358E+01 -3.8098E+01  2.3333E+00  1.8281E+01 -6.3587E+00
             5.9632E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3085.48543210108        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  1.0039E+00  2.3694E-01  1.8665E-01  1.1452E+00  2.2934E-01  1.0010E+00  1.6542E+00  3.9239E-01  9.4451E-01  6.2464E-01
             1.8549E+00
 PARAMETER:  1.0394E-01 -1.3400E+00 -1.5785E+00  2.3555E-01 -1.3725E+00  1.0099E-01  6.0334E-01 -8.3551E-01  4.2909E-02 -3.7059E-01
             7.1782E-01
 GRADIENT:  -7.7749E+00 -1.1720E+01  9.1217E+01  1.4097E+02 -1.5673E+02  7.8458E+00 -1.2860E+01 -5.2863E+00 -2.5668E+01 -4.8096E+00
             4.4172E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3114.52963866431        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  1.0027E+00  1.6744E-01  9.9958E-02  9.1135E-01  1.6568E-01  9.8961E-01  1.6347E+00  1.1264E+00  1.4535E+00  5.3645E-01
             1.6443E+00
 PARAMETER:  1.0274E-01 -1.6871E+00 -2.2030E+00  7.1749E-03 -1.6977E+00  8.9555E-02  5.9145E-01  2.1902E-01  4.7397E-01 -5.2278E-01
             5.9734E-01
 GRADIENT:   3.8504E+00 -6.1135E+00  8.1697E-02  2.1562E+01  4.2304E+01  2.5279E+00 -3.0594E+01 -3.8769E+00 -4.7470E+00  4.2416E+00
            -2.0347E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3117.06625596317        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      875
 NPARAMETR:  9.9505E-01  1.6328E-01  9.4158E-02  8.6124E-01  1.6129E-01  9.7352E-01  1.8759E+00  1.1774E+00  1.5286E+00  5.2663E-01
             1.6453E+00
 PARAMETER:  9.5042E-02 -1.7123E+00 -2.2628E+00 -4.9387E-02 -1.7246E+00  7.3165E-02  7.2909E-01  2.6329E-01  5.2438E-01 -5.4126E-01
             5.9795E-01
 GRADIENT:   9.3603E+01  1.6605E+02  1.9109E+02 -1.1397E+00  1.3530E+03  5.2154E+00  3.5152E+01  1.2509E+00  1.8784E+01  1.1624E+01
             3.2100E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3117.55872918376        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1031
 NPARAMETR:  1.0000E+00  1.6325E-01  9.4130E-02  8.8018E-01  1.6079E-01  9.8687E-01  1.7884E+00  1.1773E+00  1.5169E+00  5.2126E-01
             1.6458E+00
 PARAMETER:  1.0004E-01 -1.7125E+00 -2.2631E+00 -2.7634E-02 -1.7276E+00  8.6778E-02  6.8133E-01  2.6319E-01  5.1670E-01 -5.5150E-01
             5.9824E-01
 GRADIENT:  -2.2939E-01  3.6434E-01 -5.1439E+00  3.4007E+00  3.3796E+01  1.7385E+00 -2.7338E+00 -3.2408E+00 -2.1102E+00 -7.4844E-01
            -7.9409E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3117.84174460060        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1208            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9994E-01  1.6284E-01  9.3489E-02  8.7195E-01  1.5964E-01  9.8222E-01  1.8041E+00  1.1943E+00  1.5400E+00  5.2821E-01
             1.6496E+00
 PARAMETER:  9.9944E-02 -1.7150E+00 -2.2699E+00 -3.7023E-02 -1.7349E+00  8.2058E-02  6.9005E-01  2.7755E-01  5.3178E-01 -5.3825E-01
             6.0053E-01
 GRADIENT:   1.0486E+02  1.6650E+02  1.9898E+02  1.2315E+01  1.3327E+03  8.7671E+00  2.2364E+01  1.8253E+00  1.8500E+01  1.1929E+01
             9.7483E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3117.87455774124        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1374
 NPARAMETR:  9.9996E-01  1.6320E-01  9.4118E-02  8.7553E-01  1.5995E-01  9.8185E-01  1.8050E+00  1.2223E+00  1.5479E+00  5.1548E-01
             1.6527E+00
 PARAMETER:  1.0050E-01 -1.7177E+00 -2.2697E+00 -3.3154E-02 -1.7348E+00  8.1933E-02  6.9045E-01  3.0159E-01  5.3550E-01 -5.5710E-01
             6.0071E-01
 GRADIENT:   1.0597E+00 -1.1525E+03 -8.5611E+02 -3.0400E-01 -1.4724E+03  1.6806E-01 -6.1947E-02  6.5261E+03 -5.9995E-01  4.2420E-01
            -3.2998E+03
 NUMSIGDIG:         2.0         2.3         2.3         2.4         2.7         2.4         3.5         2.3         2.3         1.8
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1374
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.8120E-04  1.8102E-02  2.5956E-02 -5.8693E-04  3.1616E-02
 SE:             2.9717E-02  2.7381E-02  2.1906E-02  2.6901E-02  2.0424E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7634E-01  5.0853E-01  2.3606E-01  9.8259E-01  1.2164E-01

 ETASHRINKSD(%)  4.4505E-01  8.2711E+00  2.6613E+01  9.8783E+00  3.1576E+01
 ETASHRINKVR(%)  8.8812E-01  1.5858E+01  4.6143E+01  1.8781E+01  5.3181E+01
 EBVSHRINKSD(%)  8.9798E-01  6.8481E+00  2.2796E+01  7.0872E+00  3.4145E+01
 EBVSHRINKVR(%)  1.7879E+00  1.3227E+01  4.0395E+01  1.3672E+01  5.6632E+01
 RELATIVEINF(%)  9.8059E+01  5.1772E+01  1.4667E+01  5.6853E+01  9.4903E+00
 EPSSHRINKSD(%)  2.2548E+01
 EPSSHRINKVR(%)  4.0012E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3117.8745577412437     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1463.7851979728330     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3117.875       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.62E-01  9.35E-02  8.75E-01  1.60E-01  9.82E-01  1.80E+00  1.22E+00  1.55E+00  5.18E-01  1.65E+00
 


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
+        1.13E+03
 
 TH 2
+       -1.19E+02  6.46E+05
 
 TH 3
+       -2.06E+02 -2.70E+04  2.24E+06
 
 TH 4
+        1.00E+01  2.03E+03  2.33E+03  6.44E+06
 
 TH 5
+       -1.93E+02 -6.52E+05 -2.04E+06  1.52E+03  7.03E+05
 
 TH 6
+        1.17E+00 -7.83E+01 -1.15E+02 -1.10E-01 -1.45E+02  1.99E+02
 
 TH 7
+       -3.27E-01  3.66E+01 -5.06E+01  1.24E-01 -1.29E+02 -1.78E-01  4.40E+01
 
 TH 8
+        9.27E+01  4.09E+03  5.46E+03  1.52E+06  4.71E+03  5.88E+01 -1.07E+05  3.60E+05
 
 TH 9
+        8.14E+00  9.40E+02  1.47E+03  6.80E+05  1.12E+03 -5.26E-01 -4.78E+04  1.61E+05  7.19E+04
 
 TH10
+       -5.97E+00 -1.55E+03 -1.90E+03  2.74E+01  8.25E+05 -2.60E-01  7.56E+00  1.25E+03  2.24E+01  1.63E+02
 
 TH11
+       -5.08E+01  1.12E+03  1.36E+03  6.17E+02  1.72E+03 -1.96E+01 -2.70E+00  1.17E+03  2.70E+02 -4.49E+02  5.10E+04
 
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
 #CPUT: Total CPU Time in Seconds,       48.496
Stop Time:
Wed Sep 29 00:27:55 CDT 2021
