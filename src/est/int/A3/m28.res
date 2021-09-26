Fri Sep 24 22:12:33 CDT 2021
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
$DATA ../../../../data/int/A3/dat28.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1099.84896219469        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4730E+02  7.6583E+01  1.7865E+02 -1.5861E+02  1.6819E+02  4.8865E+00 -2.2955E+02 -1.0796E+02 -4.3276E+01 -2.0582E+02
            -4.9028E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2804.37254798001        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7035E-01  8.9102E-01  9.9276E-01  1.1952E+00  8.8645E-01  8.6752E-01  9.8999E-01  7.1892E-01  8.8833E-01  1.1392E+00
             2.2457E+00
 PARAMETER:  6.9906E-02 -1.5389E-02  9.2735E-02  2.7829E-01 -2.0531E-02 -4.2119E-02  8.9938E-02 -2.3000E-01 -1.8417E-02  2.3035E-01
             9.0902E-01
 GRADIENT:   2.9207E+01 -1.7153E+01  3.7141E+01  2.5612E+01  1.4296E+01 -4.5236E+01 -1.2224E+01  1.0527E+01 -2.2203E+01  1.1737E+01
            -2.3929E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2821.03180283240        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.8359E-01  7.3714E-01  6.9718E-01  1.2795E+00  6.6004E-01  9.0369E-01  1.0819E+00  1.1279E-01  9.9650E-01  8.1871E-01
             2.2497E+00
 PARAMETER:  8.3453E-02 -2.0497E-01 -2.6071E-01  3.4651E-01 -3.1546E-01 -1.2709E-03  1.7876E-01 -2.0822E+00  9.6498E-02 -1.0002E-01
             9.1078E-01
 GRADIENT:   6.4651E+01  4.2431E+01  3.2575E+00  5.4010E+01  2.1278E+00 -2.9280E+01 -1.0732E+01  5.5717E-01  7.0365E+00 -2.7447E+01
            -2.0938E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2836.23919798237        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.5536E-01  5.0902E-01  4.9833E-01  1.3312E+00  4.6023E-01  9.7003E-01  1.2783E+00  5.1893E-02  9.6342E-01  7.7177E-01
             2.4707E+00
 PARAMETER:  5.4336E-02 -5.7527E-01 -5.9649E-01  3.8605E-01 -6.7602E-01  6.9574E-02  3.4555E-01 -2.8586E+00  6.2730E-02 -1.5906E-01
             1.0045E+00
 GRADIENT:  -1.6243E+01  1.3582E+01  1.7623E+01  4.6731E+01 -7.2413E+00  4.4099E-01  1.2663E+00  1.6559E-01 -8.3020E+00  5.3287E+00
             5.6070E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2837.63083801795        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.4811E-01  3.4837E-01  3.2333E-01  1.3267E+00  3.1186E-01  9.9191E-01  1.4362E+00  1.5598E-02  1.0232E+00  6.4322E-01
             2.5178E+00
 PARAMETER:  4.6718E-02 -9.5450E-01 -1.0291E+00  3.8273E-01 -1.0652E+00  9.1877E-02  4.6202E-01 -4.0606E+00  1.2291E-01 -3.4127E-01
             1.0234E+00
 GRADIENT:  -3.4015E+01  5.9798E+00  5.3296E+01  8.8085E+01 -5.1224E+01  7.0701E+00  5.3314E+00  1.0737E-02 -2.2562E+01  4.1837E+00
             1.2886E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2838.97614319879        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.4812E-01  2.8006E-01  2.4526E-01  1.2630E+00  2.5212E-01  1.0047E+00  1.5006E+00  1.0000E-02  1.1080E+00  5.9286E-01
             2.4988E+00
 PARAMETER:  4.6721E-02 -1.1728E+00 -1.3054E+00  3.3347E-01 -1.2779E+00  1.0469E-01  5.0586E-01 -5.1730E+00  2.0256E-01 -4.2280E-01
             1.0158E+00
 GRADIENT:  -3.2246E+01 -8.6541E+00  5.5209E+01  6.9949E+01 -6.5459E+01  1.1298E+01  6.4202E+00  0.0000E+00 -2.0653E+01 -9.7437E-01
             1.2927E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2841.08711659346        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  9.5386E-01  2.4308E-01  1.9643E-01  1.1765E+00  2.1766E-01  9.9491E-01  1.5228E+00  1.0000E-02  1.2168E+00  5.9244E-01
             2.4108E+00
 PARAMETER:  5.2763E-02 -1.3144E+00 -1.5275E+00  2.6257E-01 -1.4248E+00  9.4900E-02  5.2056E-01 -6.2329E+00  2.9627E-01 -4.2350E-01
             9.7996E-01
 GRADIENT:  -1.6597E+01 -1.0754E+01  3.0045E+01  3.0108E+01 -3.5677E+01  7.8134E+00  3.4095E+00  0.0000E+00 -1.1066E+01 -2.5551E+00
             6.9836E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2846.81394700860        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      578
 NPARAMETR:  9.6122E-01  2.8659E-01  2.4608E-01  1.1970E+00  2.5896E-01  9.8350E-01  1.4295E+00  1.0000E-02  1.1806E+00  6.1599E-01
             2.3737E+00
 PARAMETER:  6.0451E-02 -1.1497E+00 -1.3021E+00  2.7979E-01 -1.2511E+00  8.3362E-02  4.5735E-01 -5.6089E+00  2.6599E-01 -3.8453E-01
             9.6446E-01
 GRADIENT:  -2.8173E+00 -1.5500E+01  1.2471E+01 -2.1550E+01 -5.9802E+00  4.1582E+00 -5.5439E+00  0.0000E+00  1.0720E+00  1.7962E+00
             3.2398E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2847.71610860488        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      753
 NPARAMETR:  9.6207E-01  3.0634E-01  2.6192E-01  1.2294E+00  2.7290E-01  9.7324E-01  1.4620E+00  1.0000E-02  1.1581E+00  6.1604E-01
             2.3405E+00
 PARAMETER:  6.1336E-02 -1.0831E+00 -1.2397E+00  3.0651E-01 -1.1986E+00  7.2876E-02  4.7982E-01 -5.2897E+00  2.4678E-01 -3.8444E-01
             9.5036E-01
 GRADIENT:  -5.0298E-02  3.9548E-02  1.2751E-02  2.3675E-02 -1.7354E-02 -6.2574E-03  1.9572E-02  0.0000E+00 -3.3790E-02  8.7844E-03
             1.4183E-01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -2847.71610860488        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      775
 NPARAMETR:  9.6207E-01  3.0634E-01  2.6192E-01  1.2294E+00  2.7290E-01  9.7324E-01  1.4620E+00  1.0000E-02  1.1581E+00  6.1604E-01
             2.3405E+00
 PARAMETER:  6.1336E-02 -1.0831E+00 -1.2397E+00  3.0651E-01 -1.1986E+00  7.2876E-02  4.7982E-01 -5.2897E+00  2.4678E-01 -3.8444E-01
             9.5036E-01
 GRADIENT:  -5.0298E-02  3.9548E-02  1.2751E-02  2.3675E-02 -1.7354E-02 -6.2574E-03  1.9572E-02  0.0000E+00 -3.3790E-02  8.7844E-03
             1.4183E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      775
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7429E-03  6.1350E-03 -5.3626E-07 -5.6921E-03  4.0338E-03
 SE:             2.9418E-02  2.4612E-02  2.8884E-04  2.8189E-02  2.3857E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5275E-01  8.0315E-01  9.9852E-01  8.3998E-01  8.6573E-01

 ETASHRINKSD(%)  1.4466E+00  1.7546E+01  9.9032E+01  5.5616E+00  2.0076E+01
 ETASHRINKVR(%)  2.8723E+00  3.2014E+01  9.9991E+01  1.0814E+01  3.6121E+01
 EBVSHRINKSD(%)  1.5457E+00  1.5880E+01  9.9011E+01  4.5456E+00  2.0401E+01
 EBVSHRINKVR(%)  3.0675E+00  2.9239E+01  9.9990E+01  8.8845E+00  3.6640E+01
 RELATIVEINF(%)  9.6909E+01  1.6144E+01  6.2226E-04  6.4212E+01  3.5111E+00
 EPSSHRINKSD(%)  1.8824E+01
 EPSSHRINKVR(%)  3.4104E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2847.7161086048827     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1193.6267488364720     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.71
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2847.716       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.62E-01  3.06E-01  2.62E-01  1.23E+00  2.73E-01  9.73E-01  1.46E+00  1.00E-02  1.16E+00  6.16E-01  2.34E+00
 


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
+        1.23E+03
 
 TH 2
+        1.13E+00  3.99E+03
 
 TH 3
+        2.12E+01 -4.33E+02  1.40E+04
 
 TH 4
+       -4.57E+00 -1.49E+01 -2.57E+02  4.99E+02
 
 TH 5
+       -8.09E+00 -3.58E+03 -1.64E+04 -1.86E+02  2.52E+04
 
 TH 6
+        6.97E+00 -2.94E+00  2.22E+01 -4.77E+00 -1.38E+01  1.98E+02
 
 TH 7
+       -3.66E-01  3.98E+01  2.46E+01 -2.10E-01 -4.38E+01 -4.96E-01  4.57E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.05E+00 -6.66E+00  1.34E+02 -3.74E+00  1.43E+01  1.01E+00 -4.59E-01  0.00E+00  1.18E+02
 
 TH10
+       -3.32E+00 -1.59E+01 -6.15E+01  1.03E+01 -1.06E+01  7.14E-02  8.13E+00  0.00E+00  6.76E-01  2.26E+02
 
 TH11
+       -1.58E+01 -1.53E+01 -1.09E+02 -8.21E+00  1.13E+02  2.21E+00  8.36E+00  0.00E+00  5.55E+00  2.22E+01  1.98E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.928
Stop Time:
Fri Sep 24 22:13:03 CDT 2021
