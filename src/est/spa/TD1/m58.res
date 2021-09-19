Sat Sep 18 14:09:09 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat58.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m58.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1629.75961812473        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6020E+02 -1.5834E+01  2.8103E+01 -4.0887E+01 -1.7353E+01  1.3752E+01 -1.9229E+00 -3.9929E+00  1.2714E+01 -9.4709E+00
             1.8057E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1636.54469760856        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.4400E-01  1.0336E+00  8.8663E-01  9.9975E-01  9.7613E-01  9.6771E-01  1.0716E+00  9.9629E-01  8.8804E-01  1.0618E+00
             9.5692E-01
 PARAMETER:  4.2375E-02  1.3302E-01 -2.0327E-02  9.9752E-02  7.5838E-02  6.7178E-02  1.6919E-01  9.6282E-02 -1.8734E-02  1.5993E-01
             5.5968E-02
 GRADIENT:   3.3410E+01 -2.4152E-01 -2.7639E+00  2.9696E+00 -8.4415E-01  8.6003E+00 -4.7357E+00  4.9342E+00 -1.3936E+00  7.7908E+00
             5.2694E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1637.55480451484        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.4183E-01  9.0813E-01  6.9477E-01  1.0700E+00  7.9821E-01  9.8572E-01  1.3444E+00  6.8291E-01  7.3893E-01  8.4248E-01
             9.3329E-01
 PARAMETER:  4.0074E-02  3.6300E-03 -2.6418E-01  1.6763E-01 -1.2538E-01  8.5615E-02  3.9594E-01 -2.8139E-01 -2.0255E-01 -7.1401E-02
             3.0957E-02
 GRADIENT:   2.5881E+01  1.4584E+01 -2.1553E+01  4.4305E+01  2.3436E+01  1.4896E+01  2.9979E+00  4.3258E+00 -1.2016E+01  2.2872E+00
            -5.9081E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1639.05543211279        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.3126E-01  8.8360E-01  5.9131E-01  1.0477E+00  7.1024E-01  9.4475E-01  1.2821E+00  3.5429E-01  8.1584E-01  7.5593E-01
             9.4641E-01
 PARAMETER:  2.8783E-02 -2.3746E-02 -4.2542E-01  1.4660E-01 -2.4215E-01  4.3163E-02  3.4847E-01 -9.3763E-01 -1.0354E-01 -1.7980E-01
             4.4917E-02
 GRADIENT:  -7.0177E+00 -2.4617E-01  4.7644E-01 -1.6291E+00 -3.1217E+00 -2.9688E+00 -8.9636E-01  1.2769E+00  1.7440E+00  2.3637E+00
             1.3372E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1640.46407103390        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      411
 NPARAMETR:  9.4968E-01  7.2727E-01  6.4804E-01  1.1534E+00  6.8801E-01  9.6349E-01  1.5625E+00  2.2118E-01  7.5504E-01  7.9363E-01
             9.4704E-01
 PARAMETER:  4.8367E-02 -2.1845E-01 -3.3381E-01  2.4273E-01 -2.7395E-01  6.2806E-02  5.4629E-01 -1.4088E+00 -1.8099E-01 -1.3113E-01
             4.5588E-02
 GRADIENT:   4.4956E+00  5.5254E+00  2.2821E+00  8.8474E+00 -3.8492E+00  1.0876E+00 -9.2904E-02  5.0397E-02 -4.7667E-01 -1.2038E+00
            -9.8964E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1640.54716005849        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  9.4746E-01  6.8214E-01  6.5374E-01  1.1735E+00  6.7889E-01  9.6010E-01  1.6449E+00  1.6385E-01  7.4323E-01  8.0911E-01
             9.4963E-01
 PARAMETER:  4.6025E-02 -2.8253E-01 -3.2504E-01  2.5999E-01 -2.8730E-01  5.9283E-02  5.9769E-01 -1.7088E+00 -1.9675E-01 -1.1182E-01
             4.8322E-02
 GRADIENT:  -9.1059E-03 -6.8649E-02  3.8303E-02 -1.3254E-01 -7.7703E-02 -2.3281E-02  2.5080E-03  1.2442E-02  2.9798E-02  2.3446E-02
             2.1667E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1640.55234588686        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      767
 NPARAMETR:  9.4760E-01  6.9635E-01  6.4505E-01  1.1643E+00  6.7840E-01  9.6030E-01  1.6146E+00  9.0731E-02  7.4840E-01  8.0646E-01
             9.4931E-01
 PARAMETER:  4.6176E-02 -2.6190E-01 -3.3843E-01  2.5210E-01 -2.8802E-01  5.9488E-02  5.7908E-01 -2.2999E+00 -1.8982E-01 -1.1510E-01
             4.7981E-02
 GRADIENT:  -2.1542E-01 -9.9828E-02 -8.0772E-02 -9.7562E-02  1.9497E-01 -6.8011E-02 -3.7008E-02  5.2337E-03 -2.1814E-03 -1.0287E-01
            -4.1680E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1640.55500976474        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      942
 NPARAMETR:  9.4775E-01  6.9900E-01  6.4282E-01  1.1625E+00  6.7791E-01  9.6063E-01  1.6091E+00  1.6500E-02  7.4953E-01  8.0752E-01
             9.4934E-01
 PARAMETER:  4.6340E-02 -2.5811E-01 -3.4190E-01  2.5060E-01 -2.8874E-01  5.9838E-02  5.7570E-01 -4.0044E+00 -1.8831E-01 -1.1379E-01
             4.8014E-02
 GRADIENT:   5.1878E-02 -4.1250E-03  7.5411E-02 -3.4920E-02 -6.6254E-02  4.1965E-02 -1.2266E-02  1.7235E-04  1.2662E-02 -1.5524E-02
            -1.5548E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1640.55507548999        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  9.4773E-01  6.9929E-01  6.4267E-01  1.1624E+00  6.7794E-01  9.6055E-01  1.6087E+00  1.0000E-02  7.4957E-01  8.0751E-01
             9.4936E-01
 PARAMETER:  4.6318E-02 -2.5769E-01 -3.4212E-01  2.5045E-01 -2.8869E-01  5.9746E-02  5.7540E-01 -4.6064E+00 -1.8825E-01 -1.1380E-01
             4.8030E-02
 GRADIENT:  -1.1989E-02 -9.2921E-03 -1.0886E-02  8.4655E-03  2.1540E-02  3.7864E-03 -2.8447E-03  0.0000E+00  2.3687E-03 -1.2527E-03
            -3.0239E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1069
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4940E-04  1.4746E-02 -5.8339E-04 -1.4187E-02  1.5306E-03
 SE:             2.9860E-02  2.1438E-02  2.4278E-04  2.4858E-02  2.3073E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9334E-01  4.9155E-01  1.6264E-02  5.6818E-01  9.4711E-01

 ETASHRINKSD(%)  1.0000E-10  2.8181E+01  9.9187E+01  1.6723E+01  2.2702E+01
 ETASHRINKVR(%)  1.0000E-10  4.8420E+01  9.9993E+01  3.0650E+01  4.0249E+01
 EBVSHRINKSD(%)  4.1107E-01  2.8049E+01  9.9274E+01  1.6801E+01  2.1523E+01
 EBVSHRINKVR(%)  8.2044E-01  4.8231E+01  9.9995E+01  3.0780E+01  3.8414E+01
 RELATIVEINF(%)  9.8830E+01  5.9290E+00  5.3733E-04  9.9951E+00  4.3886E+00
 EPSSHRINKSD(%)  4.4510E+01
 EPSSHRINKVR(%)  6.9209E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1640.5550754899934     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -905.40424892625526     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1640.555       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.48E-01  6.99E-01  6.43E-01  1.16E+00  6.78E-01  9.61E-01  1.61E+00  1.00E-02  7.50E-01  8.08E-01  9.49E-01
 


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
+        1.33E+03
 
 TH 2
+       -1.04E+01  4.80E+02
 
 TH 3
+        2.58E+01  3.41E+02  1.56E+03
 
 TH 4
+       -9.86E+00  3.44E+02 -5.02E+02  1.03E+03
 
 TH 5
+       -6.11E+00 -5.55E+02 -1.70E+03  4.83E+02  2.28E+03
 
 TH 6
+       -4.93E-01 -3.41E+00  5.60E+00 -2.04E+00 -2.35E+00  2.13E+02
 
 TH 7
+        1.83E+00  4.00E+01 -1.83E+01 -1.15E+01  4.53E+00  1.26E-01  2.76E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.15E+00 -2.80E+01 -3.76E+01  2.09E+01  1.28E+01 -8.63E-01  1.25E+01  0.00E+00  1.78E+02
 
 TH10
+       -2.76E+00 -1.12E+01 -1.31E+02 -3.24E+01 -3.05E+01 -3.73E+00  7.89E+00  0.00E+00  1.74E+01  1.23E+02
 
 TH11
+       -8.94E+00 -1.10E+01 -4.77E+01 -1.06E+01  1.85E+01  4.68E+00  3.72E+00  0.00E+00  1.72E+01  2.18E+01  2.29E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.644
Stop Time:
Sat Sep 18 14:09:28 CDT 2021
