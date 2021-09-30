Wed Sep 29 19:23:48 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat89.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1642.65888332212        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9549E+02 -4.7955E+01 -2.2058E+01 -1.1287E+01  7.9399E+01  6.1832E+01 -3.7841E+00 -9.3521E-01  1.4150E+01 -1.7095E+01
            -2.1332E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1649.14929654285        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9679E-01  1.0452E+00  9.8405E-01  1.0202E+00  9.5031E-01  8.6162E-01  1.0158E+00  1.0141E+00  9.4618E-01  1.0570E+00
             1.0927E+00
 PARAMETER:  9.6789E-02  1.4422E-01  8.3919E-02  1.2004E-01  4.9034E-02 -4.8936E-02  1.1566E-01  1.1396E-01  4.4673E-02  1.5539E-01
             1.8870E-01
 GRADIENT:  -2.9749E+00 -6.2382E+00  2.5618E+00 -1.2489E+01 -9.1982E+00 -1.9547E+01 -2.0689E+00  7.6322E-01  2.5963E+00  3.2743E+00
             1.5677E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1649.88719736724        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.9580E-01  9.2696E-01  9.2340E-01  1.0996E+00  8.6839E-01  8.8028E-01  1.2735E+00  9.9365E-01  8.1895E-01  9.5114E-01
             1.0805E+00
 PARAMETER:  9.5792E-02  2.4151E-02  2.0308E-02  1.9497E-01 -4.1117E-02 -2.7515E-02  3.4174E-01  9.3628E-02 -9.9733E-02  4.9906E-02
             1.7742E-01
 GRADIENT:  -5.7914E+00  8.0906E+00 -5.3046E+00  5.1716E+00 -5.0081E+00 -1.0915E+01  1.6240E+00  5.6669E+00 -4.0383E+00  4.2091E+00
             1.0660E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1651.34750883372        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.9956E-01  8.1050E-01  7.8994E-01  1.1612E+00  7.5611E-01  9.1976E-01  1.4341E+00  6.6978E-01  7.8783E-01  8.1982E-01
             1.0458E+00
 PARAMETER:  9.9564E-02 -1.1010E-01 -1.3580E-01  2.4945E-01 -1.7957E-01  1.6354E-02  4.6055E-01 -3.0081E-01 -1.3848E-01 -9.8667E-02
             1.4483E-01
 GRADIENT:   2.9886E+00  1.1278E+01 -5.7909E+00  2.3367E+01  8.6408E+00  5.5737E+00  7.1766E-02  1.2691E+00 -2.8170E+00 -2.3386E+00
            -4.2610E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1652.46260456174        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.9918E-01  6.5458E-01  6.4307E-01  1.2108E+00  6.1707E-01  9.0611E-01  1.6827E+00  3.7271E-01  7.5457E-01  7.0559E-01
             1.0625E+00
 PARAMETER:  9.9175E-02 -3.2377E-01 -3.4151E-01  2.9125E-01 -3.8278E-01  1.4088E-03  6.2041E-01 -8.8697E-01 -1.8161E-01 -2.4872E-01
             1.6066E-01
 GRADIENT:  -5.7374E-01  1.8457E+00  1.4935E+00  2.1086E+00 -2.7156E+00 -5.9753E-01 -5.4085E-01  5.6671E-03  2.9037E-01 -8.0866E-01
             8.2656E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1652.50344214740        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.9813E-01  5.5484E-01  7.1882E-01  1.2769E+00  6.3058E-01  9.1031E-01  1.8904E+00  4.8647E-01  7.3110E-01  7.7063E-01
             1.0646E+00
 PARAMETER:  9.8132E-02 -4.8908E-01 -2.3014E-01  3.4444E-01 -3.6112E-01  6.0254E-03  7.3681E-01 -6.2058E-01 -2.1320E-01 -1.6054E-01
             1.6258E-01
 GRADIENT:   2.2805E+00  2.6708E-01  3.0507E+00 -6.1361E+00 -7.0662E+00  1.7163E+00 -1.5611E-01 -1.5914E-01  8.4678E-01  2.0725E+00
             2.0110E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.55457984269        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  9.9697E-01  5.5491E-01  7.2005E-01  1.2783E+00  6.3227E-01  9.0611E-01  1.9014E+00  5.1340E-01  7.2788E-01  7.4935E-01
             1.0602E+00
 PARAMETER:  9.6968E-02 -4.8895E-01 -2.2843E-01  3.4555E-01 -3.5843E-01  1.4048E-03  7.4261E-01 -5.6670E-01 -2.1761E-01 -1.8854E-01
             1.5848E-01
 GRADIENT:  -6.2760E-01  3.8406E-01 -5.2060E-01 -6.3157E-01 -9.5230E-01 -9.3031E-02 -1.4666E-01 -4.1730E-02 -7.1745E-02 -2.4603E-01
            -1.2110E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1652.56285624909        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1249             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9892E-01  5.5184E-01  7.2805E-01  1.2812E+00  6.3556E-01  9.0672E-01  1.9177E+00  5.2670E-01  7.2720E-01  7.5419E-01
             1.0601E+00
 PARAMETER:  9.8920E-02 -4.9450E-01 -2.1738E-01  3.4780E-01 -3.5324E-01  2.0811E-03  7.5110E-01 -5.4113E-01 -2.1855E-01 -1.8211E-01
             1.5837E-01
 GRADIENT:   3.3055E+02  4.5722E+01  9.8564E+00  2.9245E+02  3.8258E+01  2.2347E+01  2.7699E+01  6.3136E-01  9.6011E+00  1.1471E+00
             1.2366E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1652.56603308791        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1323
 NPARAMETR:  9.9792E-01  5.5111E-01  7.2803E-01  1.2807E+00  6.3552E-01  9.0630E-01  1.9160E+00  5.2656E-01  7.2708E-01  7.5413E-01
             1.0601E+00
 PARAMETER:  9.7920E-02 -4.9582E-01 -2.1741E-01  3.4740E-01 -3.5332E-01  1.6168E-03  7.5026E-01 -5.4140E-01 -2.1872E-01 -1.8219E-01
             1.5836E-01
 GRADIENT:   1.3365E+00  1.8001E-01  2.6549E-01 -1.3204E+00 -1.7478E+00  2.7079E-03  6.1648E-02  1.7279E-02 -1.3503E-02  1.6961E-02
            -1.1284E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1323
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.3235E-05  2.3907E-02 -2.3523E-02 -2.0921E-02 -1.5205E-03
 SE:             2.9806E-02  2.0216E-02  1.1571E-02  2.4604E-02  2.1268E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9750E-01  2.3699E-01  4.2051E-02  3.9516E-01  9.4300E-01

 ETASHRINKSD(%)  1.4450E-01  3.2273E+01  6.1237E+01  1.7572E+01  2.8751E+01
 ETASHRINKVR(%)  2.8879E-01  5.4130E+01  8.4974E+01  3.2057E+01  4.9235E+01
 EBVSHRINKSD(%)  5.7526E-01  3.4277E+01  6.2844E+01  1.6739E+01  2.6383E+01
 EBVSHRINKVR(%)  1.1472E+00  5.6805E+01  8.6194E+01  3.0675E+01  4.5805E+01
 RELATIVEINF(%)  9.8105E+01  5.5651E+00  1.2074E+00  1.1412E+01  3.9303E+00
 EPSSHRINKSD(%)  4.4405E+01
 EPSSHRINKVR(%)  6.9093E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1652.5660330879055     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -917.41520652416727     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1652.566       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  5.51E-01  7.28E-01  1.28E+00  6.36E-01  9.06E-01  1.92E+00  5.27E-01  7.27E-01  7.54E-01  1.06E+00
 


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
+        1.34E+03
 
 TH 2
+       -1.80E+01  4.62E+02
 
 TH 3
+        2.45E+01  2.67E+02  1.12E+03
 
 TH 4
+       -1.05E+01  3.61E+02 -3.28E+02  9.03E+02
 
 TH 5
+       -3.30E+00 -5.64E+02 -1.56E+03  3.29E+02  2.63E+03
 
 TH 6
+       -1.38E+00 -2.41E+00  3.75E+00 -2.87E+00 -2.09E+00  2.37E+02
 
 TH 7
+        1.62E+00  3.74E+01 -6.45E+00 -8.78E+00  3.31E-01 -2.09E-01  1.81E+01
 
 TH 8
+        2.40E-01 -8.38E+00 -8.90E+01 -4.67E+00  4.03E+01  3.99E-01  1.66E+00  3.43E+01
 
 TH 9
+        2.28E+00 -2.70E+01 -1.52E+01  1.16E+01  5.26E+00 -9.48E-01  6.26E+00  6.30E+00  2.09E+02
 
 TH10
+       -1.28E+00  6.62E+00 -7.33E+01 -4.21E+01 -4.43E+01  4.29E-02  9.08E+00  3.23E+01  6.76E+00  1.01E+02
 
 TH11
+       -9.35E+00 -4.86E+00 -1.62E+01 -1.02E+01  7.89E+00  3.50E+00  1.19E+00  8.14E+00  1.39E+01  1.68E+01  1.89E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.470
Stop Time:
Wed Sep 29 19:24:12 CDT 2021
