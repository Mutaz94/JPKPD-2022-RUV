Wed Sep 29 13:37:30 CDT 2021
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
$DATA ../../../../data/spa/A3/dat54.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   952.364333782445        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8034E+02  1.7294E+02  1.7029E+02  1.3350E+01  1.6768E+02  2.5793E+01 -7.5323E+01 -7.6474E+01 -2.2999E+02 -1.8326E+02
            -4.6394E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -737.445948109402        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0554E+00  7.5558E-01  6.4911E-01  1.1088E+00  6.1276E-01  9.6116E-01  9.7234E-01  9.9302E-01  1.3790E+00  1.1881E+00
             1.7895E+00
 PARAMETER:  1.5395E-01 -1.8027E-01 -3.3216E-01  2.0326E-01 -3.8979E-01  6.0390E-02  7.1951E-02  9.2999E-02  4.2136E-01  2.7235E-01
             6.8192E-01
 GRADIENT:   2.5034E+02  9.0746E+01  9.5372E+01  3.8423E+01 -3.3488E+01 -1.1079E+01 -9.9471E+00  8.0114E-01 -1.3810E+01  1.5694E+01
            -1.2660E+03

0ITERATION NO.:   10    OBJECTIVE VALUE:  -770.279915750245        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0595E+00  6.0199E-01  2.0577E-01  1.1504E+00  2.5286E-01  7.9883E-01  4.7004E-01  9.6506E-01  1.0434E+00  8.7201E-01
             1.7542E+00
 PARAMETER:  1.5784E-01 -4.0751E-01 -1.4810E+00  2.4010E-01 -1.2749E+00 -1.2461E-01 -6.5493E-01  6.4439E-02  1.4246E-01 -3.6949E-02
             6.6204E-01
 GRADIENT:   2.6251E+02  3.9481E+02  3.2047E+02  2.9748E+02 -2.5142E+02 -9.4976E+01 -3.0103E+01 -1.5944E+02 -2.5491E+02 -2.7708E+01
            -8.7800E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1048.70478033148        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0624E+00  6.0269E-01  2.0472E-01  1.0317E+00  3.5072E-01  9.8806E-01  1.0000E-02  1.7105E+00  1.7181E+00  1.0261E+00
             2.1109E+00
 PARAMETER:  1.6057E-01 -4.0636E-01 -1.4861E+00  1.3123E-01 -9.4776E-01  8.7989E-02 -5.3123E+00  6.3677E-01  6.4120E-01  1.2574E-01
             8.4713E-01
 GRADIENT:   9.8061E+01 -7.1518E+01  5.1067E+01  8.5497E+00  1.8359E+02  7.9188E+00  0.0000E+00 -6.5672E+01 -4.7674E+01  4.9864E+01
            -4.2628E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1208.03617825312        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      473
 NPARAMETR:  1.0316E+00  4.7808E-01  1.1993E-01  9.2601E-01  2.4435E-01  9.9391E-01  8.9487E-02  1.8688E+00  2.1007E+00  3.1227E-01
             3.4031E+00
 PARAMETER:  1.3114E-01 -6.3798E-01 -2.0209E+00  2.3135E-02 -1.3091E+00  9.3888E-02 -2.3137E+00  7.2532E-01  8.4225E-01 -1.0639E+00
             1.3247E+00
 GRADIENT:   1.5287E+01 -6.7505E+01 -2.0647E+01 -2.0129E+01  9.9232E+01  1.5172E+01  1.1959E-01  5.7766E+00  2.1395E+01 -1.0548E+00
            -2.0928E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1212.19450954854        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      649
 NPARAMETR:  1.0279E+00  5.2234E-01  1.2404E-01  9.5505E-01  2.5134E-01  9.5419E-01  2.5599E-01  1.8278E+00  1.9141E+00  3.1373E-01
             3.5243E+00
 PARAMETER:  1.2747E-01 -5.4943E-01 -1.9871E+00  5.4006E-02 -1.2810E+00  5.3107E-02 -1.2626E+00  7.0311E-01  7.4927E-01 -1.0592E+00
             1.3597E+00
 GRADIENT:  -1.7690E+00 -2.1231E+00 -1.4411E+00  5.1253E-01  3.4300E+00  1.0327E+00  6.3764E-01  2.5471E+00  3.6684E+00  6.6420E-03
            -9.5724E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1212.46651133539        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  1.0294E+00  5.1390E-01  1.2208E-01  9.4653E-01  2.4746E-01  9.5231E-01  1.6141E-01  1.7831E+00  1.8910E+00  3.5561E-01
             3.5919E+00
 PARAMETER:  1.2901E-01 -5.6573E-01 -2.0031E+00  4.5052E-02 -1.2965E+00  5.1138E-02 -1.7238E+00  6.7836E-01  7.3710E-01 -9.3392E-01
             1.3787E+00
 GRADIENT:  -2.0529E+00 -3.7935E+00 -4.3583E-01 -2.8968E+00  5.4517E+00 -1.1725E-01  2.9256E-01  6.6067E-01  1.5977E+00  1.2967E+00
             1.6482E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1212.64023596612        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1001
 NPARAMETR:  1.0306E+00  5.1703E-01  1.2271E-01  9.5266E-01  2.4908E-01  9.5112E-01  4.2848E-02  1.7877E+00  1.8738E+00  3.2777E-01
             3.6095E+00
 PARAMETER:  1.3011E-01 -5.5965E-01 -1.9979E+00  5.1504E-02 -1.2900E+00  4.9882E-02 -3.0501E+00  6.8094E-01  7.2795E-01 -1.0154E+00
             1.3836E+00
 GRADIENT:  -4.2141E-01 -8.3902E-01 -8.4378E-01  4.1763E-01  1.6208E+00 -3.1918E-01  1.4505E-02  9.9515E-02 -1.0345E-01  1.2176E-01
             4.4143E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1212.64809187659        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  1.0309E+00  5.1708E-01  1.2308E-01  9.5227E-01  2.4916E-01  9.5200E-01  1.0000E-02  1.7885E+00  1.8752E+00  3.2583E-01
             3.6072E+00
 PARAMETER:  1.3038E-01 -5.5956E-01 -1.9949E+00  5.1092E-02 -1.2897E+00  5.0811E-02 -4.6177E+00  6.8139E-01  7.2869E-01 -1.0214E+00
             1.3829E+00
 GRADIENT:   6.3393E-02  1.7875E-02  7.9900E-02 -4.4843E-02 -1.0509E-01  6.4847E-03  0.0000E+00  4.8738E-03  5.6102E-02  6.9796E-03
            -3.1661E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1163
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.1154E-03 -4.4945E-04  1.9416E-02 -5.1611E-03  2.0678E-02
 SE:             2.8272E-02  1.9582E-04  2.1118E-02  2.5813E-02  1.1192E-02
 N:                     100         100         100         100         100

 P VAL.:         8.2875E-01  2.1722E-02  3.5787E-01  8.4153E-01  6.4660E-02

 ETASHRINKSD(%)  5.2857E+00  9.9344E+01  2.9254E+01  1.3523E+01  6.2506E+01
 ETASHRINKVR(%)  1.0292E+01  9.9996E+01  4.9950E+01  2.5217E+01  8.5942E+01
 EBVSHRINKSD(%)  5.4122E+00  9.9321E+01  3.0327E+01  1.0802E+01  6.3442E+01
 EBVSHRINKVR(%)  1.0532E+01  9.9995E+01  5.1457E+01  2.0437E+01  8.6635E+01
 RELATIVEINF(%)  8.5490E+01  3.4166E-04  2.0125E+01  6.6420E+01  8.0963E-01
 EPSSHRINKSD(%)  3.1405E+01
 EPSSHRINKVR(%)  5.2947E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1212.6480918765892     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -477.49726531285103     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.83
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1212.648       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  5.17E-01  1.23E-01  9.52E-01  2.49E-01  9.52E-01  1.00E-02  1.79E+00  1.88E+00  3.26E-01  3.61E+00
 


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
+        1.04E+03
 
 TH 2
+       -1.56E+01  2.27E+03
 
 TH 3
+       -4.36E+02  2.52E+03  1.14E+04
 
 TH 4
+       -2.58E+01  1.44E+02 -3.00E+02  2.76E+02
 
 TH 5
+        3.39E+02 -6.94E+03 -1.05E+04 -2.82E+02  2.45E+04
 
 TH 6
+       -3.26E+00 -2.25E+01  2.37E+00 -6.08E+00  1.11E+02  1.75E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.60E-01  6.70E+00 -7.15E+01  3.87E+00 -5.21E+01  3.18E+00  0.00E+00  2.12E+01
 
 TH 9
+        1.52E+01 -4.56E+01  4.15E+01 -9.62E+00  1.41E+02  3.72E+00  0.00E+00  1.99E+00  2.50E+01
 
 TH10
+       -2.88E+00 -1.07E+02 -5.58E+01 -1.25E+00  4.57E+02 -8.73E-01  0.00E+00  4.11E+00  2.63E+00  6.10E+01
 
 TH11
+       -1.90E+01 -4.22E+01  3.18E+01 -1.64E+00  1.11E+02  2.08E+00  0.00E+00  2.29E+00  5.84E+00  1.47E+01  2.27E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.671
Stop Time:
Wed Sep 29 13:37:56 CDT 2021
