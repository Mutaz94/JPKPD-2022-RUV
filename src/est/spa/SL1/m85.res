Sat Sep 25 10:47:08 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat85.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1648.60928719097        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.3791E+01 -7.5520E+01 -2.3372E+01 -8.0648E+01  9.2394E+00  1.9229E+01 -2.3539E+01  5.7682E+00 -1.7398E+01  1.9135E+01
            -1.5318E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1654.32777821041        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:       92
 NPARAMETR:  1.0056E+00  1.1651E+00  1.1156E+00  9.1439E-01  1.1154E+00  9.2415E-01  1.3136E+00  9.7881E-01  1.0937E+00  8.1403E-01
             1.0325E+00
 PARAMETER:  1.0560E-01  2.5283E-01  2.0942E-01  1.0498E-02  2.0925E-01  2.1120E-02  3.7280E-01  7.8584E-02  1.8955E-01 -1.0576E-01
             1.3199E-01
 GRADIENT:   1.0562E+02 -1.4672E+01  2.3311E+01 -4.6785E+01  1.6463E+01 -1.2903E+01  1.2005E+01 -8.0251E+00  4.1446E+00 -1.0827E+01
             2.4134E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1657.16060304579        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      168
 NPARAMETR:  9.9768E-01  1.0544E+00  8.1733E-01  9.7620E-01  9.2179E-01  9.3544E-01  1.4227E+00  7.0665E-01  9.8555E-01  5.4941E-01
             1.0294E+00
 PARAMETER:  9.7677E-02  1.5297E-01 -1.0171E-01  7.5912E-02  1.8560E-02  3.3266E-02  4.5257E-01 -2.4722E-01  8.5442E-02 -4.9891E-01
             1.2898E-01
 GRADIENT:   7.7459E+01 -2.0072E+01 -3.6110E+00 -1.4845E+01  4.7922E+01 -7.4788E+00  5.4317E+00 -1.0860E+00  2.8365E+00 -1.3984E+01
            -4.6887E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1660.12431413564        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      242
 NPARAMETR:  9.6724E-01  1.0894E+00  6.8643E-01  9.6694E-01  8.4845E-01  9.5051E-01  1.3357E+00  4.8161E-01  9.4357E-01  6.0860E-01
             1.0033E+00
 PARAMETER:  6.6686E-02  1.8566E-01 -2.7625E-01  6.6381E-02 -6.4349E-02  4.9244E-02  3.8943E-01 -6.3063E-01  4.1915E-02 -3.9659E-01
             1.0329E-01
 GRADIENT:  -8.1864E-01  9.0178E-01 -3.8300E+00  5.1305E+00  8.8467E+00 -9.8327E-01 -9.7468E-01  1.0464E+00 -4.6128E-01 -1.6834E+00
            -3.4707E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1660.45276612089        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      338
 NPARAMETR:  9.7170E-01  1.0267E+00  6.3165E-01  9.9012E-01  7.8570E-01  9.4994E-01  1.3980E+00  2.4758E-01  9.1398E-01  5.7349E-01
             1.0095E+00
 PARAMETER:  7.1296E-02  1.2633E-01 -3.5942E-01  9.0074E-02 -1.4118E-01  4.8648E-02  4.3502E-01 -1.2960E+00  1.0058E-02 -4.5602E-01
             1.0944E-01
 GRADIENT:  -2.9181E+01 -8.6143E+00 -1.9923E+00 -8.2425E+00  2.6294E-01 -5.9931E+00 -2.9910E+00  2.5384E-01  3.6694E-01 -1.3568E-01
            -2.0088E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1660.83152012489        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      514
 NPARAMETR:  9.8269E-01  1.0953E+00  6.3914E-01  9.6079E-01  8.2072E-01  9.6340E-01  1.3516E+00  1.7390E-01  9.4390E-01  6.1120E-01
             1.0115E+00
 PARAMETER:  8.2540E-02  1.9104E-01 -3.4763E-01  5.9998E-02 -9.7575E-02  6.2708E-02  4.0130E-01 -1.6493E+00  4.2262E-02 -3.9233E-01
             1.1139E-01
 GRADIENT:  -1.9903E+00  1.5548E-01 -6.5016E-01  4.9079E-01  2.6223E-01 -1.1729E-01  1.8054E-02  6.2594E-02  1.9485E-01  2.1854E-02
             3.3204E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1660.85260182104        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      689
 NPARAMETR:  9.8355E-01  1.0969E+00  6.3166E-01  9.5837E-01  8.1724E-01  9.6388E-01  1.3516E+00  5.4582E-02  9.4362E-01  6.1090E-01
             1.0111E+00
 PARAMETER:  8.3408E-02  1.9246E-01 -3.5940E-01  5.7483E-02 -1.0182E-01  6.3211E-02  4.0128E-01 -2.8081E+00  4.1969E-02 -3.9283E-01
             1.1109E-01
 GRADIENT:  -4.0880E-02 -1.0177E-01 -3.3956E-01 -8.9642E-02  2.1015E-01  5.8875E-02  2.8209E-01  4.4535E-03  6.5842E-02  1.2016E-01
             2.0934E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1660.85477108498        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      864
 NPARAMETR:  9.8356E-01  1.0999E+00  6.3077E-01  9.5657E-01  8.1812E-01  9.6377E-01  1.3470E+00  1.0000E-02  9.4507E-01  6.1105E-01
             1.0108E+00
 PARAMETER:  8.3427E-02  1.9522E-01 -3.6081E-01  5.5603E-02 -1.0074E-01  6.3096E-02  3.9787E-01 -4.6799E+00  4.3502E-02 -3.9258E-01
             1.1072E-01
 GRADIENT:  -1.6113E-02 -7.3103E-02 -8.3998E-02  4.8435E-03  1.1871E-01  7.0249E-03  7.6961E-03  0.0000E+00  6.8130E-03  1.1349E-02
            -1.3264E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1660.85477773970        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      956
 NPARAMETR:  9.8357E-01  1.0996E+00  6.3086E-01  9.5678E-01  8.1799E-01  9.6375E-01  1.3473E+00  1.0000E-02  9.4490E-01  6.1093E-01
             1.0108E+00
 PARAMETER:  8.3433E-02  1.9496E-01 -3.6068E-01  5.5815E-02 -1.0091E-01  6.3077E-02  3.9810E-01 -4.6613E+00  4.3324E-02 -3.9277E-01
             1.1074E-01
 GRADIENT:  -4.6139E-04 -2.2557E-03  1.4199E-03 -3.2108E-03  4.0827E-04  5.9151E-05 -1.2436E-03  0.0000E+00 -5.0454E-04 -6.5240E-04
            -3.7202E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      956
 NO. OF SIG. DIGITS IN FINAL EST.:  4.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6116E-04  1.3292E-03 -4.9799E-04 -2.5311E-03 -8.6122E-03
 SE:             2.9843E-02  2.4843E-02  2.0060E-04  2.4558E-02  1.8850E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9569E-01  9.5733E-01  1.3047E-02  9.1791E-01  6.4776E-01

 ETASHRINKSD(%)  2.3719E-02  1.6773E+01  9.9328E+01  1.7726E+01  3.6849E+01
 ETASHRINKVR(%)  4.7432E-02  3.0733E+01  9.9995E+01  3.2310E+01  6.0120E+01
 EBVSHRINKSD(%)  4.6105E-01  1.5697E+01  9.9403E+01  1.7975E+01  3.7661E+01
 EBVSHRINKVR(%)  9.1997E-01  2.8931E+01  9.9996E+01  3.2720E+01  6.1139E+01
 RELATIVEINF(%)  9.8971E+01  7.9560E+00  4.0043E-04  7.9564E+00  3.3180E+00
 EPSSHRINKSD(%)  4.3529E+01
 EPSSHRINKVR(%)  6.8110E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1660.8547777397041     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -925.70395117596593     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1660.855       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.10E+00  6.31E-01  9.57E-01  8.18E-01  9.64E-01  1.35E+00  1.00E-02  9.45E-01  6.11E-01  1.01E+00
 


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
+        1.22E+03
 
 TH 2
+       -5.10E+00  3.69E+02
 
 TH 3
+        1.63E+01  2.91E+02  1.10E+03
 
 TH 4
+       -1.28E+01  1.98E+02 -5.34E+02  9.31E+02
 
 TH 5
+       -8.05E+00 -4.35E+02 -1.21E+03  6.12E+02  1.73E+03
 
 TH 6
+        1.51E+00 -3.08E-01  2.95E+00 -3.02E+00 -2.05E+00  2.11E+02
 
 TH 7
+        9.51E-01  2.73E+01 -4.33E+01 -1.48E+01  2.34E+01  8.70E-02  5.73E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.40E-01 -1.96E+01 -4.07E+01  4.34E+01 -2.17E+01  3.01E-01  9.27E+00  0.00E+00  1.11E+02
 
 TH10
+       -1.08E+00 -1.38E+01 -7.63E+01 -3.53E+01 -6.62E+01  3.69E-01  1.63E+01  0.00E+00  1.89E+01  1.01E+02
 
 TH11
+       -7.68E+00 -1.07E+01 -3.77E+01 -7.97E+00 -4.63E+00  1.83E+00  5.72E+00  0.00E+00  1.25E+01  2.73E+01  2.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.074
Stop Time:
Sat Sep 25 10:47:25 CDT 2021
