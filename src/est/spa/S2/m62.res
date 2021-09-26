Sat Sep 25 12:24:08 CDT 2021
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
$DATA ../../../../data/spa/S2/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1685.15717690476        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5796E+00  4.8683E+00 -3.4052E+01  4.2806E+01  7.2023E+01  3.5830E+01 -7.0328E+00  4.2672E+00 -2.7588E+01  3.0428E+00
             1.4639E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1690.38542366147        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0188E+00  9.2532E-01  9.9061E-01  1.0308E+00  9.3593E-01  8.7491E-01  1.0100E+00  9.6237E-01  1.1957E+00  8.9316E-01
             9.4488E-01
 PARAMETER:  1.1863E-01  2.2389E-02  9.0565E-02  1.3032E-01  3.3780E-02 -3.3635E-02  1.0990E-01  6.1647E-02  2.7873E-01 -1.2988E-02
             4.3299E-02
 GRADIENT:   6.4961E+01 -7.0637E+00 -2.0584E+01  2.8554E+01  4.8637E+01 -1.3709E+01  5.0476E+00  2.1704E+00  2.1738E+01 -5.4560E+00
            -7.7111E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1692.28142394614        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0155E+00  9.3319E-01  7.2855E-01  1.0201E+00  7.8747E-01  9.0761E-01  1.3113E+00  5.6299E-01  1.0826E+00  7.8020E-01
             9.0552E-01
 PARAMETER:  1.1536E-01  3.0852E-02 -2.1670E-01  1.1985E-01 -1.3893E-01  3.0639E-03  3.7100E-01 -4.7449E-01  1.7937E-01 -1.4821E-01
             7.4906E-04
 GRADIENT:   5.2125E+01  1.8973E+01 -1.6216E+01  4.7431E+01  9.8623E+00  7.0747E-01  1.1613E+01  3.9966E+00  1.6692E+01  7.2002E+00
            -2.5078E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1694.25339190667        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0002E+00  8.8192E-01  6.8863E-01  1.0113E+00  7.4589E-01  9.1481E-01  1.2824E+00  4.1930E-01  9.9144E-01  7.3901E-01
             9.5556E-01
 PARAMETER:  1.0023E-01 -2.5652E-02 -2.7306E-01  1.1125E-01 -1.9317E-01  1.0965E-02  3.4875E-01 -7.6917E-01  9.1402E-02 -2.0244E-01
             5.4547E-02
 GRADIENT:  -3.4405E+00 -2.1712E+00 -3.8325E+00 -2.6608E+00 -6.5535E-02  2.7272E+00  1.2476E+00  1.9888E+00 -7.7184E-01  5.9122E+00
            -2.2402E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1694.64116704513        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0033E+00  8.9318E-01  5.9770E-01  9.9061E-01  6.9644E-01  9.0808E-01  1.2720E+00  1.7025E-01  9.9731E-01  6.3151E-01
             9.6680E-01
 PARAMETER:  1.0330E-01 -1.2962E-02 -4.1467E-01  9.0570E-02 -2.6177E-01  3.5729E-03  3.4060E-01 -1.6705E+00  9.7307E-02 -3.5965E-01
             6.6236E-02
 GRADIENT:   1.1667E+00 -1.0660E+00 -1.7447E+00 -1.5207E-01  1.0114E+00 -5.9062E-01  8.4568E-02  3.3029E-01  4.1201E-01 -3.2168E-01
             4.8414E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1695.46188021002        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      421
 NPARAMETR:  1.0135E+00  9.1340E-01  6.4708E-01  9.9178E-01  7.3602E-01  9.1805E-01  1.2606E+00  3.6885E-02  1.0091E+00  7.0020E-01
             9.7094E-01
 PARAMETER:  1.1343E-01  9.4178E-03 -3.3528E-01  9.1745E-02 -2.0650E-01  1.4493E-02  3.3157E-01 -3.1999E+00  1.0902E-01 -2.5640E-01
             7.0511E-02
 GRADIENT:  -2.0042E+01 -8.7661E-01  1.6282E+00 -4.8876E+00 -1.4733E+00 -1.7086E-01 -1.2898E+00  8.1845E-03 -1.3935E+00 -6.2343E-01
             2.5995E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1695.59024268656        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      596
 NPARAMETR:  1.0209E+00  8.7259E-01  6.5598E-01  1.0189E+00  7.2505E-01  9.1809E-01  1.3193E+00  3.1989E-02  9.9326E-01  7.0592E-01
             9.6517E-01
 PARAMETER:  1.2066E-01 -3.6293E-02 -3.2163E-01  1.1870E-01 -2.2152E-01  1.4539E-02  3.7713E-01 -3.3424E+00  9.3233E-02 -2.4826E-01
             6.4550E-02
 GRADIENT:  -6.3544E-01 -1.5719E-01 -1.1980E-01 -2.3161E-01 -4.3640E-02  8.9769E-02 -1.6880E-02  5.6728E-03 -5.5325E-02  1.7597E-01
             7.2511E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1695.59304519811        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      772
 NPARAMETR:  1.0211E+00  8.7441E-01  6.5531E-01  1.0179E+00  7.2536E-01  9.1787E-01  1.3177E+00  1.0000E-02  9.9439E-01  7.0471E-01
             9.6507E-01
 PARAMETER:  1.2091E-01 -3.4211E-02 -3.2265E-01  1.1771E-01 -2.2109E-01  1.4304E-02  3.7589E-01 -4.7616E+00  9.4370E-02 -2.4996E-01
             6.4442E-02
 GRADIENT:   1.7242E-03  8.2087E-03  5.8569E-03  1.3218E-02 -6.8307E-04 -9.5286E-03  1.5369E-03  0.0000E+00  5.4014E-03 -6.2843E-03
            -9.7277E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1695.59304519811        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      794
 NPARAMETR:  1.0211E+00  8.7441E-01  6.5531E-01  1.0179E+00  7.2536E-01  9.1787E-01  1.3177E+00  1.0000E-02  9.9439E-01  7.0471E-01
             9.6507E-01
 PARAMETER:  1.2091E-01 -3.4211E-02 -3.2265E-01  1.1771E-01 -2.2109E-01  1.4304E-02  3.7589E-01 -4.7616E+00  9.4370E-02 -2.4996E-01
             6.4442E-02
 GRADIENT:   1.7242E-03  8.2087E-03  5.8569E-03  1.3218E-02 -6.8307E-04 -9.5286E-03  1.5369E-03  0.0000E+00  5.4014E-03 -6.2843E-03
            -9.7277E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      794
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.8495E-05  2.2086E-03 -5.0088E-04 -3.1255E-03 -6.8806E-03
 SE:             2.9841E-02  2.1765E-02  2.1052E-04  2.5918E-02  2.1785E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9870E-01  9.1917E-01  1.7347E-02  9.0402E-01  7.5213E-01

 ETASHRINKSD(%)  2.8930E-02  2.7085E+01  9.9295E+01  1.3171E+01  2.7016E+01
 ETASHRINKVR(%)  5.7852E-02  4.6834E+01  9.9995E+01  2.4607E+01  4.6734E+01
 EBVSHRINKSD(%)  4.6720E-01  2.6546E+01  9.9354E+01  1.3197E+01  2.7059E+01
 EBVSHRINKVR(%)  9.3221E-01  4.6045E+01  9.9996E+01  2.4652E+01  4.6797E+01
 RELATIVEINF(%)  9.8838E+01  4.8561E+00  5.1302E-04  9.8928E+00  4.1661E+00
 EPSSHRINKSD(%)  4.4651E+01
 EPSSHRINKVR(%)  6.9365E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1695.5930451981055     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -960.44221863436735     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.00
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1695.593       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  8.74E-01  6.55E-01  1.02E+00  7.25E-01  9.18E-01  1.32E+00  1.00E-02  9.94E-01  7.05E-01  9.65E-01
 


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
+        1.25E+03
 
 TH 2
+       -7.18E+00  4.73E+02
 
 TH 3
+        2.09E+01  3.79E+02  1.19E+03
 
 TH 4
+       -1.12E+01  2.61E+02 -3.89E+02  8.25E+02
 
 TH 5
+       -4.20E+00 -6.39E+02 -1.43E+03  4.46E+02  2.24E+03
 
 TH 6
+       -1.18E+00 -7.95E-01  3.83E+00 -1.82E+00 -1.69E+00  2.32E+02
 
 TH 7
+        1.06E+00  3.47E+01 -2.04E+01 -7.46E+00 -3.66E+00  1.45E-01  3.83E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.31E+00 -2.73E+01 -2.62E+01  2.87E+01 -9.21E+00  1.21E+00  1.24E+01  0.00E+00  1.16E+02
 
 TH10
+       -2.39E+00 -1.31E+01 -1.24E+02 -2.97E+01 -3.77E+01  5.81E-02  1.71E+01  0.00E+00  1.33E+01  1.30E+02
 
 TH11
+       -7.83E+00 -1.04E+01 -3.80E+01 -7.61E+00  1.72E+01  3.36E+00  5.03E+00  0.00E+00  8.60E+00  2.47E+01  2.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       13.717
Stop Time:
Sat Sep 25 12:24:23 CDT 2021
