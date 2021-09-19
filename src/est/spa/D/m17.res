Sat Sep 18 15:10:06 CDT 2021
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
$DATA ../../../../data/spa/D/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -752.688618423151        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.1723E+02 -1.5213E+02 -3.4170E+01 -3.8258E+02  1.8823E+02 -7.3926E+02 -2.8255E+02 -2.5982E+01 -6.2487E+02 -2.4380E+02
            -1.0485E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1341.25272966788        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0422E+00  8.6559E-01  1.3355E+00  1.3334E+00  1.1414E+00  1.7277E+00  1.6147E+00  1.0826E+00  1.8948E+00  1.8328E+00
             1.0775E+00
 PARAMETER:  1.4135E-01 -4.4348E-02  3.8931E-01  3.8770E-01  2.3226E-01  6.4682E-01  5.7914E-01  1.7936E-01  7.3914E-01  7.0584E-01
             1.7467E-01
 GRADIENT:  -3.7439E+01  1.2918E+01 -5.7273E+01  2.4780E+01  5.2987E+01  1.0159E+01 -1.6769E+01 -1.5473E+01 -3.4195E+01 -2.9453E+01
            -9.1220E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1387.50828341610        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  1.2995E+00  1.7450E-01  7.2566E+00  1.7835E+00  2.0562E+00  1.8534E+00  6.3951E-01  6.4671E+00  1.5172E+00  2.1077E+00
             1.3379E+00
 PARAMETER:  3.6198E-01 -1.6459E+00  2.0819E+00  6.7858E-01  8.2086E-01  7.1703E-01 -3.4705E-01  1.9667E+00  5.1684E-01  8.4558E-01
             3.9108E-01
 GRADIENT:   1.7637E+02  1.3503E+01  9.5053E+00  5.1644E+01  1.7243E+01  4.5340E+01 -1.5729E-01  4.4017E+01 -1.8327E+01 -1.8958E+01
             3.5701E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1394.91910484613        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      238
 NPARAMETR:  1.3055E+00  1.5647E-01  6.8519E+00  1.6890E+00  1.9429E+00  1.9051E+00  4.8547E-01  5.1504E+00  1.5456E+00  2.5009E+00
             1.2266E+00
 PARAMETER:  3.6662E-01 -1.7549E+00  2.0245E+00  6.2416E-01  7.6420E-01  7.4452E-01 -6.2263E-01  1.7391E+00  5.3542E-01  1.0167E+00
             3.0423E-01
 GRADIENT:   1.9190E+02  4.0647E+00  1.3090E+01  3.3106E-01  1.4638E+00  7.0524E+01  1.1008E-03  2.1140E+01 -1.4184E+01  1.5095E+01
             1.4026E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1401.05734717907        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.2312E+00  1.6108E-01  6.8638E+00  1.7345E+00  1.9369E+00  1.8377E+00  5.2263E-01  5.0393E+00  1.5940E+00  2.4417E+00
             1.2002E+00
 PARAMETER:  3.0798E-01 -1.7259E+00  2.0263E+00  6.5074E-01  7.6111E-01  7.0852E-01 -5.4888E-01  1.7173E+00  5.6622E-01  9.9271E-01
             2.8250E-01
 GRADIENT:   5.4009E+01  4.1661E+00  1.2202E+01 -4.4694E+01 -3.1980E+00 -1.9643E+01  7.3030E-03 -2.2325E+01 -1.1135E+01  2.9536E+00
             2.3523E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1401.26712335927        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:      575
 NPARAMETR:  1.2303E+00  1.6076E-01  6.8729E+00  1.7360E+00  1.9380E+00  1.8397E+00  5.2249E-01  5.0686E+00  1.5941E+00  2.4425E+00
             1.1993E+00
 PARAMETER:  3.0726E-01 -1.7278E+00  2.0276E+00  6.5161E-01  7.6165E-01  7.0962E-01 -5.4915E-01  1.7231E+00  5.6634E-01  9.9303E-01
             2.8173E-01
 GRADIENT:   5.2674E+01  3.6456E+00  4.1754E+00 -4.7035E+01 -9.7114E+00 -1.8905E+01  8.3580E-03 -1.1596E+01 -9.9785E+00 -1.3617E+00
            -4.5147E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1401.84833118887        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  1.2302E+00  1.6084E-01  6.8771E+00  1.7364E+00  1.9376E+00  1.9191E+00  5.1835E-01  5.1475E+00  1.5944E+00  2.4418E+00
             1.2293E+00
 PARAMETER:  3.0716E-01 -1.7273E+00  2.0282E+00  6.5182E-01  7.6143E-01  7.5185E-01 -5.5710E-01  1.7385E+00  5.6652E-01  9.9273E-01
             3.0647E-01
 GRADIENT:   4.8132E+01  3.5801E+00  6.1287E+00 -4.7280E+01 -8.1531E+00 -1.9248E+00  1.2034E-02 -1.2387E+01 -1.1087E+01  1.5800E+00
             9.2119E+00

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1401.84833118887        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  1.2302E+00  1.6084E-01  6.8771E+00  1.7364E+00  1.9376E+00  1.9191E+00  5.1807E-01  5.1474E+00  1.5944E+00  2.4418E+00
             1.2293E+00
 PARAMETER:  3.0716E-01 -1.7273E+00  2.0282E+00  6.5182E-01  7.6143E-01  7.5185E-01 -5.5710E-01  1.7385E+00  5.6652E-01  9.9273E-01
             3.0647E-01
 GRADIENT:  -7.4805E+03  6.7248E+02  7.7283E+02  3.4956E+03 -2.0226E+03  3.0698E+03  1.1032E-02  4.2131E+02  4.0618E+03  2.3334E+03
            -7.5352E+03
 NUMSIGDIG:         3.3         3.3         3.6         3.3         3.3         3.3         1.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      783
 NO. OF SIG. DIGITS IN FINAL EST.:  1.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.6333E-02 -4.4540E-03 -6.9368E-02  1.4519E-02 -6.1722E-02
 SE:             2.9031E-02  1.1186E-03  1.8501E-02  2.8734E-02  2.1915E-02
 N:                     100         100         100         100         100

 P VAL.:         1.1049E-01  6.8395E-05  1.7738E-04  6.1336E-01  4.8557E-03

 ETASHRINKSD(%)  2.7435E+00  9.6253E+01  3.8018E+01  3.7361E+00  2.6583E+01
 ETASHRINKVR(%)  5.4117E+00  9.9860E+01  6.1582E+01  7.3325E+00  4.6100E+01
 EBVSHRINKSD(%)  6.4748E-01  9.7397E+01  5.2479E+01  3.6521E+00  2.3423E+01
 EBVSHRINKVR(%)  1.2908E+00  9.9932E+01  7.7417E+01  7.1708E+00  4.1360E+01
 RELATIVEINF(%)  9.8145E+01  2.5236E-02  1.9528E+01  3.5918E+01  4.7503E+01
 EPSSHRINKSD(%)  4.5071E+01
 EPSSHRINKVR(%)  6.9828E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1401.8483311888656     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -666.69750462512741     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.93
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1401.848       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.23E+00  1.61E-01  6.88E+00  1.74E+00  1.94E+00  1.92E+00  5.18E-01  5.15E+00  1.59E+00  2.44E+00  1.23E+00
 


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
+        4.05E+06
 
 TH 2
+        5.37E+03  7.48E+06
 
 TH 3
+       -4.18E+01 -1.58E+02  4.68E+03
 
 TH 4
+       -4.80E+02 -1.40E+03 -1.00E+02  4.50E+05
 
 TH 5
+        3.40E+05  4.64E+05  7.86E+01  1.14E+05  3.52E+05
 
 TH 6
+        2.18E+03  1.44E+06 -5.95E+01 -7.45E+02 -8.84E+04  2.77E+05
 
 TH 7
+       -5.86E+00  2.54E+00  1.39E-01  1.86E+00 -4.44E+05  2.30E+00 -7.60E-01
 
 TH 8
+       -6.79E+02  2.33E+05 -2.66E+03  2.18E+02 -1.95E+02  9.37E+01 -1.32E-01  1.14E+04
 
 TH 9
+        6.83E+03 -2.39E+03 -1.82E+02 -2.33E+03 -1.40E+05 -1.80E+03  8.16E-01  2.86E+02  7.06E+05
 
 TH10
+        2.05E+05  1.14E+06  5.14E+01  6.97E+04  5.25E+04  4.81E+02 -2.18E+00 -5.45E+00 -8.57E+04  3.23E+04
 
 TH11
+        4.19E+02  3.88E+02 -2.01E+01  3.09E+02  3.41E+05 -8.70E+01 -5.31E+00  1.71E+05 -2.59E+02  2.07E+05  4.07E+06
 
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
 #CPUT: Total CPU Time in Seconds,       22.597
Stop Time:
Sat Sep 18 15:10:30 CDT 2021
