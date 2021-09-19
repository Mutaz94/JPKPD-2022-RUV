Sat Sep 18 09:15:46 CDT 2021
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
$DATA ../../../../data/spa/A1/dat46.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1463.78644143949        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2041E+02 -9.1336E+01 -5.8637E+00 -1.1350E+02  4.4142E+01  1.6154E+01 -1.1005E+00  1.0222E+00  1.7319E+01 -1.0211E+01
            -3.5844E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1548.22710706824        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.7345E-01  1.0734E+00  1.0266E+00  1.0325E+00  1.0055E+00  9.4031E-01  9.6556E-01  9.7230E-01  8.6255E-01  9.6527E-01
             1.5524E+00
 PARAMETER:  7.3092E-02  1.7079E-01  1.2624E-01  1.3197E-01  1.0545E-01  3.8453E-02  6.4955E-02  7.1905E-02 -4.7867E-02  6.4654E-02
             5.3983E-01
 GRADIENT:   3.0111E+01 -1.6660E+01  1.1510E+00 -2.4122E+01  2.6801E+00 -2.2472E+00 -6.1406E-02  1.6123E+00  4.9825E+00  1.8907E-01
            -2.7931E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1549.85537750135        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  9.7185E-01  9.4051E-01  7.9067E-01  1.1147E+00  8.2824E-01  9.5730E-01  1.1758E+00  5.3944E-01  7.4901E-01  8.4816E-01
             1.5702E+00
 PARAMETER:  7.1451E-02  3.8665E-02 -1.3488E-01  2.0859E-01 -8.8457E-02  5.6360E-02  2.6192E-01 -5.1723E-01 -1.8900E-01 -6.4684E-02
             5.5123E-01
 GRADIENT:   1.9416E+01 -4.3304E+00 -1.8584E+01  1.7116E+01  2.3743E+01  3.7723E+00  1.5136E+00  1.8659E+00  1.2123E+00  3.7876E+00
            -1.7511E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1550.93163488829        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  9.6589E-01  9.1147E-01  7.0097E-01  1.1169E+00  7.4708E-01  9.4780E-01  1.1954E+00  3.2172E-01  7.2812E-01  7.4109E-01
             1.6436E+00
 PARAMETER:  6.5297E-02  7.3049E-03 -2.5529E-01  2.1056E-01 -1.9158E-01  4.6385E-02  2.7846E-01 -1.0341E+00 -2.1729E-01 -1.9963E-01
             5.9686E-01
 GRADIENT:  -1.7783E-01  1.2505E+00 -1.4732E+00  3.8632E+00  5.5152E-01  7.0200E-02  3.9522E-01  5.1881E-01  4.1046E-01  7.6934E-01
             2.3875E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1551.27475589334        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      413
 NPARAMETR:  9.6909E-01  7.5995E-01  7.5901E-01  1.2135E+00  7.2369E-01  9.4839E-01  1.3802E+00  1.6399E-01  6.9438E-01  7.8215E-01
             1.6362E+00
 PARAMETER:  6.8600E-02 -1.7450E-01 -1.7574E-01  2.9349E-01 -2.2339E-01  4.7010E-02  4.2220E-01 -1.7079E+00 -2.6474E-01 -1.4571E-01
             5.9238E-01
 GRADIENT:  -1.6980E+00  2.8429E+00  1.7794E+00  4.1395E+00 -3.2999E+00  5.1946E-02 -3.1997E-01  6.1594E-02 -3.4545E-01  9.1774E-02
            -6.3500E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1551.30873414305        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  9.6928E-01  7.0953E-01  7.6285E-01  1.2395E+00  7.1004E-01  9.4747E-01  1.4600E+00  9.8138E-02  6.8490E-01  7.8557E-01
             1.6367E+00
 PARAMETER:  6.8798E-02 -2.4315E-01 -1.7069E-01  3.1468E-01 -2.4244E-01  4.6035E-02  4.7845E-01 -2.2214E+00 -2.7848E-01 -1.4135E-01
             5.9268E-01
 GRADIENT:   1.3509E-03 -2.9467E-01  5.6056E-02 -7.3335E-01 -1.5616E-01  4.6812E-03  4.4621E-02  2.2328E-02  9.6973E-03  1.1084E-01
             5.9272E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1551.32006436418        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      767
 NPARAMETR:  9.6958E-01  7.3015E-01  7.5421E-01  1.2268E+00  7.1209E-01  9.4794E-01  1.4283E+00  3.0343E-02  6.8837E-01  7.8188E-01
             1.6369E+00
 PARAMETER:  6.9103E-02 -2.1451E-01 -1.8209E-01  3.0441E-01 -2.3955E-01  4.6532E-02  4.5647E-01 -3.3952E+00 -2.7343E-01 -1.4605E-01
             5.9278E-01
 GRADIENT:   6.0272E-02 -1.1140E-01  1.3360E-01 -7.0772E-01 -1.4104E-01  4.1721E-02  4.9100E-02  2.0275E-03 -6.0647E-02  2.7130E-02
            -6.5089E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1551.32113205895        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      942
 NPARAMETR:  9.6953E-01  7.2752E-01  7.5521E-01  1.2287E+00  7.1178E-01  9.4779E-01  1.4313E+00  1.0000E-02  6.8833E-01  7.8257E-01
             1.6371E+00
 PARAMETER:  6.9051E-02 -2.1812E-01 -1.8076E-01  3.0599E-01 -2.3998E-01  4.6375E-02  4.5856E-01 -4.5168E+00 -2.7348E-01 -1.4517E-01
             5.9294E-01
 GRADIENT:  -8.1542E-03  8.2446E-03 -1.0526E-02  3.7033E-02  8.9163E-03 -3.5703E-03 -2.0844E-03  0.0000E+00  2.2388E-03  1.2903E-03
             3.3420E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1551.32113207486        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      999
 NPARAMETR:  9.6953E-01  7.2733E-01  7.5526E-01  1.2288E+00  7.1175E-01  9.4779E-01  1.4316E+00  1.0000E-02  6.8828E-01  7.8258E-01
             1.6371E+00
 PARAMETER:  6.9052E-02 -2.1838E-01 -1.8069E-01  3.0607E-01 -2.4003E-01  4.6380E-02  4.5879E-01 -4.5198E+00 -2.7356E-01 -1.4515E-01
             5.9294E-01
 GRADIENT:  -2.1483E-04  3.7379E-04  3.1846E-04  7.6714E-04 -5.8047E-04 -1.1019E-04 -5.3336E-05  0.0000E+00  5.0027E-05 -4.0832E-05
            -1.1719E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      999
 NO. OF SIG. DIGITS IN FINAL EST.:  4.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4249E-04  1.3543E-02 -2.9039E-04 -1.9445E-02 -8.4857E-03
 SE:             2.9558E-02  2.0084E-02  1.8395E-04  2.2627E-02  2.0821E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9076E-01  5.0010E-01  1.1442E-01  3.9011E-01  6.8360E-01

 ETASHRINKSD(%)  9.7541E-01  3.2717E+01  9.9384E+01  2.4198E+01  3.0247E+01
 ETASHRINKVR(%)  1.9413E+00  5.4730E+01  9.9996E+01  4.2541E+01  5.1345E+01
 EBVSHRINKSD(%)  1.1822E+00  3.3852E+01  9.9383E+01  2.3514E+01  2.8942E+01
 EBVSHRINKVR(%)  2.3504E+00  5.6244E+01  9.9996E+01  4.1499E+01  4.9508E+01
 RELATIVEINF(%)  9.6764E+01  3.1454E+00  2.1908E-04  5.4462E+00  2.4045E+00
 EPSSHRINKSD(%)  3.7492E+01
 EPSSHRINKVR(%)  6.0927E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1551.3211320748642     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -816.17030551112600     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.13
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1551.321       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.70E-01  7.27E-01  7.55E-01  1.23E+00  7.12E-01  9.48E-01  1.43E+00  1.00E-02  6.88E-01  7.83E-01  1.64E+00
 


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
+        1.29E+03
 
 TH 2
+       -2.23E+01  4.64E+02
 
 TH 3
+        2.60E+01  2.56E+02  7.82E+02
 
 TH 4
+       -3.29E+01  4.06E+02 -2.84E+02  9.21E+02
 
 TH 5
+       -3.48E+00 -5.03E+02 -1.12E+03  2.70E+02  1.80E+03
 
 TH 6
+        3.14E+00 -4.07E+00  9.78E+00 -9.99E+00 -6.12E+00  2.11E+02
 
 TH 7
+        2.03E+00  3.32E+01  2.76E+00 -7.88E+00 -1.12E+01 -9.20E-02  2.52E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.26E+00 -1.57E+01 -2.66E+01 -1.15E+01  2.50E+01 -1.57E+00  2.10E+01  0.00E+00  1.53E+02
 
 TH10
+       -3.78E+00  4.51E+00 -4.01E+01 -2.39E+01 -3.55E+01 -1.74E-01  6.10E+00  0.00E+00  4.94E+00  8.86E+01
 
 TH11
+       -1.21E+01 -1.14E+01 -3.14E+01 -1.14E+01  5.77E+00  3.38E+00  4.07E+00  0.00E+00  1.81E+01  2.56E+01  9.10E+01
 
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
 #CPUT: Total CPU Time in Seconds,       16.942
Stop Time:
Sat Sep 18 09:16:04 CDT 2021
