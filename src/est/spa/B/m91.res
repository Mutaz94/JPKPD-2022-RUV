Sat Sep 18 08:48:35 CDT 2021
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
$DATA ../../../../data/spa/B/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1605.74806092900        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8547E+02 -8.2834E+01 -5.7863E+01 -6.5209E+01  4.7278E+01 -1.8819E+01 -4.8090E+00  1.8421E+01 -8.6852E+00  2.5988E+01
            -4.8011E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1614.72325269373        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9253E-01  1.0582E+00  1.2873E+00  9.8574E-01  1.1233E+00  1.0937E+00  1.0100E+00  8.3867E-01  1.0792E+00  7.9436E-01
             1.2278E+00
 PARAMETER:  9.2504E-02  1.5656E-01  3.5256E-01  8.5634E-02  2.1630E-01  1.8960E-01  1.0999E-01 -7.5932E-02  1.7625E-01 -1.3022E-01
             3.0519E-01
 GRADIENT:   1.3390E+02 -4.9579E+01 -7.4863E+00 -4.9158E+01  6.8344E+01  2.0728E+01  2.9188E+00  1.2831E+00  3.0396E+00 -2.0219E+01
             2.8303E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1619.51816650537        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.7690E-01  9.4801E-01  1.1068E+00  1.0809E+00  1.0274E+00  1.0952E+00  8.3640E-01  1.9743E-01  1.0182E+00  1.0406E+00
             1.1490E+00
 PARAMETER:  7.6625E-02  4.6612E-02  2.0152E-01  1.7781E-01  1.2700E-01  1.9097E-01 -7.8652E-02 -1.5224E+00  1.1805E-01  1.3976E-01
             2.3891E-01
 GRADIENT:   1.0730E+02 -3.3395E+01 -4.0177E+01 -5.7011E+00  6.2745E+01  2.4783E+01 -3.9257E+00  1.0960E-01 -8.7090E+00  1.4854E+01
             1.6566E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1625.98038781429        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.2272E-01  9.9441E-01  9.9590E-01  1.0388E+00  9.5706E-01  1.0107E+00  1.0995E+00  2.4584E-01  1.0255E+00  8.5510E-01
             1.1033E+00
 PARAMETER:  1.9568E-02  9.4390E-02  9.5890E-02  1.3811E-01  5.6107E-02  1.1063E-01  1.9487E-01 -1.3031E+00  1.2515E-01 -5.6539E-02
             1.9832E-01
 GRADIENT:   2.4228E-01 -1.7447E+01 -8.4367E+00 -1.3957E+01  1.6867E+01 -4.4569E+00  8.1928E-01  4.6846E-01  5.0404E+00  5.5292E+00
            -9.3899E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1627.01315526137        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.2023E-01  9.1211E-01  9.0778E-01  1.0958E+00  8.6202E-01  1.0215E+00  1.2357E+00  1.7682E-01  9.3302E-01  7.2768E-01
             1.1107E+00
 PARAMETER:  1.6869E-02  8.0087E-03  3.2429E-03  1.9146E-01 -4.8480E-02  1.2131E-01  3.1161E-01 -1.6326E+00  3.0672E-02 -2.1789E-01
             2.0501E-01
 GRADIENT:  -6.4340E+00  3.7498E+00  2.0892E+00  2.7918E+00 -3.6376E+00 -3.0910E-01 -2.4491E-01  3.7273E-01 -1.6064E+00 -2.6857E-01
             1.6793E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1627.47867712991        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.2521E-01  8.5743E-01  8.9114E-01  1.1198E+00  8.3361E-01  1.0238E+00  1.2959E+00  7.0522E-02  9.2139E-01  7.1147E-01
             1.1086E+00
 PARAMETER:  2.2267E-02 -5.3821E-02 -1.5257E-02  2.1313E-01 -8.1984E-02  1.2351E-01  3.5920E-01 -2.5518E+00  1.8130E-02 -2.4042E-01
             2.0310E-01
 GRADIENT:   5.5013E+00 -3.1222E+00  9.3194E-01 -4.8394E+00 -1.8042E+00  9.8266E-01 -5.6962E-01  6.7204E-02  9.4048E-01  3.2876E-01
             9.7521E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1628.35265981097        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      555
 NPARAMETR:  9.4004E-01  6.5406E-01  9.0391E-01  1.2600E+00  7.6305E-01  1.0333E+00  1.6571E+00  1.0000E-02  8.3368E-01  6.9797E-01
             1.1051E+00
 PARAMETER:  3.8163E-02 -3.2455E-01 -1.0276E-03  3.3113E-01 -1.7043E-01  1.3278E-01  6.0505E-01 -5.3003E+00 -8.1905E-02 -2.5958E-01
             1.9992E-01
 GRADIENT:   1.1580E+01  8.4193E+00  4.7954E+00  1.5190E+01 -9.1253E+00  1.4147E+00 -2.9432E-01  0.0000E+00 -9.2800E-01  9.7827E-04
            -1.5059E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1628.61848777579        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  9.3364E-01  5.4459E-01  9.0029E-01  1.3108E+00  7.3161E-01  1.0277E+00  1.9125E+00  1.0000E-02  8.0403E-01  6.9005E-01
             1.1068E+00
 PARAMETER:  3.1338E-02 -5.0772E-01 -5.0340E-03  3.7067E-01 -2.1251E-01  1.2737E-01  7.4840E-01 -7.9268E+00 -1.1812E-01 -2.7099E-01
             2.0145E-01
 GRADIENT:  -1.4278E-02  3.3655E-02 -2.2672E-01  2.1607E-01  2.8504E-01 -6.8220E-02  4.3219E-02  0.0000E+00  1.0682E-01 -5.8795E-02
            -1.7341E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1628.61864530085        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  9.3362E-01  5.4229E-01  9.0126E-01  1.3122E+00  7.3137E-01  1.0279E+00  1.9177E+00  1.0000E-02  8.0314E-01  6.9112E-01
             1.1068E+00
 PARAMETER:  3.1313E-02 -5.1195E-01 -3.9594E-03  3.7169E-01 -2.1283E-01  1.2752E-01  7.5114E-01 -7.9843E+00 -1.1923E-01 -2.6945E-01
             2.0144E-01
 GRADIENT:   9.3943E-03  1.3217E-02  2.1234E-02  9.4579E-03 -5.7700E-02  6.4517E-03  3.3488E-03  0.0000E+00  2.1785E-02  5.2833E-03
            -7.7474E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      823
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3827E-04  2.2504E-02 -3.9990E-04 -1.9994E-02 -1.5264E-03
 SE:             2.9763E-02  1.9804E-02  2.1165E-04  2.5102E-02  2.1635E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8825E-01  2.5581E-01  5.8835E-02  4.2573E-01  9.4375E-01

 ETASHRINKSD(%)  2.9023E-01  3.3654E+01  9.9291E+01  1.5906E+01  2.7519E+01
 ETASHRINKVR(%)  5.7962E-01  5.5982E+01  9.9995E+01  2.9282E+01  4.7464E+01
 EBVSHRINKSD(%)  5.3963E-01  3.6022E+01  9.9258E+01  1.4717E+01  2.5233E+01
 EBVSHRINKVR(%)  1.0764E+00  5.9068E+01  9.9994E+01  2.7269E+01  4.4100E+01
 RELATIVEINF(%)  9.8434E+01  5.0507E+00  4.6924E-04  1.2361E+01  3.9254E+00
 EPSSHRINKSD(%)  4.1497E+01
 EPSSHRINKVR(%)  6.5774E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1628.6186453008497     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -893.46781873711154     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1628.619       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.34E-01  5.42E-01  9.01E-01  1.31E+00  7.31E-01  1.03E+00  1.92E+00  1.00E-02  8.03E-01  6.91E-01  1.11E+00
 


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
+        1.19E+03
 
 TH 2
+       -1.31E+01  4.56E+02
 
 TH 3
+        1.33E+01  2.67E+02  8.00E+02
 
 TH 4
+       -6.14E+00  3.09E+02 -2.19E+02  7.43E+02
 
 TH 5
+       -4.96E+00 -5.70E+02 -1.29E+03  2.73E+02  2.42E+03
 
 TH 6
+        1.09E+00 -2.43E+00  2.24E+00 -3.44E+00 -5.43E+00  1.83E+02
 
 TH 7
+        1.16E+00  3.74E+01  2.68E+00 -7.28E+00 -5.44E+00 -1.74E-01  1.61E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.73E+00 -2.21E+01 -2.39E+01  4.49E+00  2.01E+01 -1.65E+00  1.04E+01  0.00E+00  1.70E+02
 
 TH10
+       -1.15E+00  4.11E+00 -6.98E+01 -3.52E+01 -3.33E+01  2.30E+00  6.01E+00  0.00E+00  9.65E+00  1.33E+02
 
 TH11
+       -1.11E+01 -9.62E+00 -4.03E+01 -8.28E+00  2.18E+01  1.44E+00  2.27E+00  0.00E+00  1.11E+01  3.50E+01  1.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.818
Stop Time:
Sat Sep 18 08:48:52 CDT 2021
