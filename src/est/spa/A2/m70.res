Wed Sep 29 13:00:05 CDT 2021
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
$DATA ../../../../data/spa/A2/dat70.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m70.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1156.94124576461        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8762E+02  5.7013E+01  3.7093E+01  4.3244E+01  5.6115E+01  5.6912E+01 -1.9302E+01 -5.8064E+00 -3.9861E+01 -2.1627E+01
            -9.5500E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1449.57450833332        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1933E+00  9.8883E-01  1.0609E+00  1.0500E+00  9.7814E-01  1.0720E+00  9.8649E-01  9.0380E-01  1.0737E+00  7.8016E-01
             2.4095E+00
 PARAMETER:  2.7676E-01  8.8764E-02  1.5911E-01  1.4877E-01  7.7900E-02  1.6950E-01  8.6399E-02 -1.1496E-03  1.7107E-01 -1.4825E-01
             9.7943E-01
 GRADIENT:   4.6981E+02  1.8187E+01  1.6418E+01  1.4453E+01 -2.2514E+01  2.1311E+01  4.4838E+00  3.0746E+00  5.6911E+00  8.7405E+00
             1.0688E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1458.48724464856        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.1647E+00  9.7534E-01  5.5079E-01  1.0209E+00  7.0140E-01  1.0790E+00  1.0297E+00  5.7914E-01  1.0907E+00  4.2827E-01
             2.2341E+00
 PARAMETER:  2.5246E-01  7.5027E-02 -4.9639E-01  1.2068E-01 -2.5467E-01  1.7603E-01  1.2923E-01 -4.4622E-01  1.8679E-01 -7.4799E-01
             9.0385E-01
 GRADIENT:   4.1892E+02  1.0223E+01 -1.5384E+01  4.6015E+01  2.5155E+01  3.3787E+01 -2.3635E+00  4.7703E+00  1.5641E+01  4.1705E+00
            -1.3952E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1470.87168945564        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0374E+00  1.0304E+00  4.4374E-01  9.5044E-01  6.4841E-01  9.4316E-01  1.0636E+00  3.5704E-01  1.0820E+00  4.0747E-01
             2.1234E+00
 PARAMETER:  1.3671E-01  1.2994E-01 -7.1251E-01  4.9175E-02 -3.3322E-01  4.1483E-02  1.6163E-01 -9.2991E-01  1.7884E-01 -7.9778E-01
             8.5301E-01
 GRADIENT:   1.1331E+02  3.0574E+01 -6.1024E+00  3.3887E+01  1.3492E+01  5.5669E+00  5.8210E+00  2.1252E+00  1.4028E+01  5.2445E+00
            -3.0241E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1476.39386112286        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      386
 NPARAMETR:  1.0380E+00  8.8287E-01  6.4381E-01  1.0728E+00  7.2158E-01  9.3395E-01  1.2668E+00  3.1432E-01  9.4639E-01  3.1397E-01
             2.4337E+00
 PARAMETER:  1.3728E-01 -2.4579E-02 -3.4034E-01  1.7029E-01 -2.2631E-01  3.1671E-02  3.3650E-01 -1.0573E+00  4.4896E-02 -1.0585E+00
             9.8942E-01
 GRADIENT:  -2.0746E+01  1.5038E+01  7.2447E-01  1.1418E+01 -2.0103E+00 -2.5121E+00  2.9018E+00  7.2279E-01  1.9982E+00  1.7203E+00
             1.5454E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1480.77747642245        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  1.0367E+00  5.1560E-01  4.1741E-01  1.1864E+00  4.5143E-01  9.5824E-01  1.8569E+00  7.5775E-02  8.3598E-01  1.2079E-01
             2.2355E+00
 PARAMETER:  1.3609E-01 -5.6242E-01 -7.7368E-01  2.7092E-01 -6.9533E-01  5.7344E-02  7.1891E-01 -2.4800E+00 -7.9150E-02 -2.0137E+00
             9.0445E-01
 GRADIENT:  -2.0773E+01  3.8034E+00  2.3983E+00  1.8171E+01 -3.8163E+00  2.5616E+00  1.5930E+00 -4.6421E-02  8.7458E-01 -7.3410E-01
            -8.3392E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1481.48047672585        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  1.0480E+00  5.4501E-01  3.6530E-01  1.1444E+00  4.1915E-01  9.5835E-01  1.6519E+00  5.8765E-02  8.4407E-01  2.4459E-01
             2.1827E+00
 PARAMETER:  1.4685E-01 -5.0694E-01 -9.0704E-01  2.3484E-01 -7.6953E-01  5.7457E-02  6.0193E-01 -2.7342E+00 -6.9521E-02 -1.3082E+00
             8.8056E-01
 GRADIENT:   2.4221E+00  1.0306E+01  5.2460E+00  3.5517E+00 -1.2352E+01  9.0867E-01 -4.1474E-01  3.9670E-02 -1.2810E+00 -5.4323E-01
            -3.0013E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1482.68425906203        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      912
 NPARAMETR:  1.0471E+00  4.1124E-01  3.9730E-01  1.2166E+00  4.0856E-01  9.5240E-01  2.0302E+00  1.0404E-02  8.2980E-01  3.8720E-01
             2.1679E+00
 PARAMETER:  1.4600E-01 -7.8857E-01 -8.2307E-01  2.9607E-01 -7.9513E-01  5.1227E-02  8.0814E-01 -4.4655E+00 -8.6567E-02 -8.4881E-01
             8.7374E-01
 GRADIENT:   3.8594E+00  1.6976E-01  1.2974E+00 -2.2595E+00 -9.6235E-01 -8.2134E-01 -9.7277E-01  1.7668E-03  1.0749E+00 -1.7526E-01
             5.7268E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1482.72789643677        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1039
 NPARAMETR:  1.0450E+00  3.8832E-01  4.0056E-01  1.2295E+00  4.0581E-01  9.5421E-01  2.1382E+00  1.0000E-02  8.2244E-01  4.0424E-01
             2.1610E+00
 PARAMETER:  1.4397E-01 -8.4594E-01 -8.1489E-01  3.0658E-01 -8.0188E-01  5.3126E-02  8.5996E-01 -4.8168E+00 -9.5481E-02 -8.0575E-01
             8.7059E-01
 GRADIENT:   7.8920E-02  3.9237E-01  1.4695E-01  8.2909E-01 -4.6836E-01 -2.8636E-02  1.2849E-01  0.0000E+00 -6.5434E-02 -3.0160E-02
            -4.6500E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1039
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6700E-03  3.9469E-02 -5.1557E-04 -2.5113E-02  1.7668E-02
 SE:             2.9365E-02  2.0282E-02  2.2787E-04  2.5606E-02  1.3773E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5465E-01  5.1657E-02  2.3661E-02  3.2672E-01  1.9956E-01

 ETASHRINKSD(%)  1.6243E+00  3.2051E+01  9.9237E+01  1.4217E+01  5.3858E+01
 ETASHRINKVR(%)  3.2222E+00  5.3830E+01  9.9994E+01  2.6412E+01  7.8709E+01
 EBVSHRINKSD(%)  1.9089E+00  3.4765E+01  9.9149E+01  1.3436E+01  5.1359E+01
 EBVSHRINKVR(%)  3.7814E+00  5.7444E+01  9.9993E+01  2.5066E+01  7.6340E+01
 RELATIVEINF(%)  9.5614E+01  9.8746E+00  2.3318E-04  2.7161E+01  6.8323E-01
 EPSSHRINKSD(%)  3.5041E+01
 EPSSHRINKVR(%)  5.7803E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1482.7278964367670     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -747.57706987302879     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.43
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1482.728       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  3.88E-01  4.01E-01  1.23E+00  4.06E-01  9.54E-01  2.14E+00  1.00E-02  8.22E-01  4.04E-01  2.16E+00
 


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
+        1.08E+03
 
 TH 2
+       -4.19E+01  7.40E+02
 
 TH 3
+        1.23E+01  1.05E+03  5.50E+03
 
 TH 4
+       -3.71E+01  2.25E+02 -7.09E+02  8.24E+02
 
 TH 5
+        6.01E+01 -1.86E+03 -7.17E+03  3.96E+02  1.01E+04
 
 TH 6
+       -7.25E-01 -6.82E+00  1.04E+01 -8.98E+00  1.02E+01  2.02E+02
 
 TH 7
+        1.70E+00  4.95E+01 -1.64E+01 -5.16E+00 -5.04E+00  3.38E-01  1.30E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.59E+00 -1.45E+01 -3.06E+01 -1.68E+01  7.80E+01 -1.32E+00  7.75E+00  0.00E+00  1.72E+02
 
 TH10
+       -3.62E+00  1.03E+01 -2.43E+02 -2.59E+01  2.74E+02 -1.11E+00  6.80E+00  0.00E+00 -4.13E+00  9.89E+01
 
 TH11
+       -1.26E+01 -3.79E+00 -7.92E+01 -1.32E+01  6.09E+01  3.99E+00  3.07E+00  0.00E+00  1.08E+01  2.61E+01  5.81E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.074
Stop Time:
Wed Sep 29 13:00:25 CDT 2021
