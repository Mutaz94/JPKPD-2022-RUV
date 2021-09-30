Wed Sep 29 14:10:53 CDT 2021
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
$DATA ../../../../data/spa/S1/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.90735207591        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7871E+02 -1.1627E+01 -4.5528E+01  5.6960E+01  9.2896E+01  3.6098E+01  9.4432E+00  8.5285E+00  1.0681E+01  3.9182E+00
            -4.8515E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1657.82638565304        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7267E-01  1.0201E+00  1.0556E+00  1.0053E+00  9.5918E-01  1.0088E+00  9.4402E-01  9.5504E-01  9.8028E-01  9.3355E-01
             1.0287E+00
 PARAMETER:  7.2288E-02  1.1987E-01  1.5414E-01  1.0529E-01  5.8324E-02  1.0871E-01  4.2395E-02  5.3994E-02  8.0080E-02  3.1236E-02
             1.2832E-01
 GRADIENT:   3.8382E-01  6.1362E+00  7.0761E+00 -2.5301E+00 -2.0543E+00  2.3862E+00  4.0048E+00  2.1366E+00 -3.3828E+00 -1.5716E+00
             2.7143E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1658.73857305695        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.7769E-01  9.5755E-01  9.9152E-01  1.0403E+00  9.0720E-01  9.9171E-01  7.5389E-01  7.1825E-01  1.0235E+00  9.8262E-01
             1.0199E+00
 PARAMETER:  7.7439E-02  5.6620E-02  9.1482E-02  1.3954E-01  2.6054E-03  9.1676E-02 -1.8251E-01 -2.3094E-01  1.2320E-01  8.2463E-02
             1.1971E-01
 GRADIENT:   1.1564E+01  4.3987E+00  3.6850E+00  3.3960E+00 -1.1880E+01 -4.5927E+00 -3.5040E-01  1.1299E+00  1.3575E+00  4.0737E+00
             7.1865E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1659.62875514021        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.7335E-01  1.1441E+00  7.7908E-01  9.1613E-01  8.9010E-01  1.0058E+00  8.6141E-01  4.4054E-01  1.0571E+00  8.9258E-01
             1.0187E+00
 PARAMETER:  7.2986E-02  2.3464E-01 -1.4964E-01  1.2405E-02 -1.6416E-02  1.0574E-01 -4.9186E-02 -7.1976E-01  1.5549E-01 -1.3634E-02
             1.1848E-01
 GRADIENT:  -2.7941E+00  7.7566E+00  1.2531E+00  6.2980E+00 -6.6949E+00  3.1061E-01  1.6611E-01  8.1277E-01 -2.1187E-01  6.1129E-01
             6.8157E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1659.96050429340        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.7552E-01  1.2448E+00  7.1890E-01  8.4649E-01  9.1257E-01  1.0048E+00  8.2158E-01  2.0725E-01  1.1182E+00  9.0489E-01
             1.0158E+00
 PARAMETER:  7.5211E-02  3.1901E-01 -2.3003E-01 -6.6661E-02  8.5132E-03  1.0474E-01 -9.6521E-02 -1.4738E+00  2.1173E-01  5.6660E-05
             1.1571E-01
 GRADIENT:   1.3024E+00  3.5732E+00  2.1939E-01  2.2723E+00 -2.9334E+00 -2.0614E-01  4.6505E-01  1.8477E-01  1.8133E-02  9.7844E-01
            -5.4543E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1660.00140298468        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.7509E-01  1.2973E+00  6.9918E-01  8.1063E-01  9.3185E-01  1.0052E+00  7.9445E-01  1.2571E-01  1.1584E+00  9.0973E-01
             1.0179E+00
 PARAMETER:  7.4776E-02  3.6028E-01 -2.5785E-01 -1.0995E-01  2.9419E-02  1.0521E-01 -1.3011E-01 -1.9738E+00  2.4702E-01  5.3969E-03
             1.1770E-01
 GRADIENT:   1.9991E-01 -2.9154E-01 -6.3521E-01  4.8709E-01  8.8369E-01 -2.3601E-02  1.3834E-01  6.5982E-02  2.5748E-01  1.0760E-01
            -5.3378E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1660.03124681077        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  9.7453E-01  1.2963E+00  6.9928E-01  8.1079E-01  9.3114E-01  1.0051E+00  7.9479E-01  2.8049E-02  1.1574E+00  9.1131E-01
             1.0182E+00
 PARAMETER:  7.4202E-02  3.5953E-01 -2.5771E-01 -1.0974E-01  2.8658E-02  1.0508E-01 -1.2968E-01 -3.4738E+00  2.4616E-01  7.1297E-03
             1.1807E-01
 GRADIENT:  -1.0416E+00 -5.7114E-02  1.5522E-01 -3.6896E-01 -2.5111E-01 -6.8439E-02  1.0034E-01  3.2523E-03  9.2268E-02  1.6966E-01
             7.1068E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1660.03571050923        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1201
 NPARAMETR:  9.7493E-01  1.2930E+00  6.9749E-01  8.1309E-01  9.2840E-01  1.0053E+00  7.9665E-01  1.2014E-02  1.1534E+00  9.0754E-01
             1.0180E+00
 PARAMETER:  7.4611E-02  3.5695E-01 -2.6026E-01 -1.0692E-01  2.5707E-02  1.0528E-01 -1.2734E-01 -4.3217E+00  2.4274E-01  2.9784E-03
             1.1786E-01
 GRADIENT:  -1.9534E-01  5.8674E-02 -2.3538E-01  1.9910E-01  2.3365E-01 -1.6371E-03 -2.6381E-02  6.2416E-04 -2.1357E-02 -2.2193E-02
             2.7614E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1660.03706397921        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  9.7506E-01  1.2911E+00  6.9794E-01  8.1402E-01  9.2735E-01  1.0053E+00  7.9848E-01  1.0000E-02  1.1523E+00  9.0647E-01
             1.0180E+00
 PARAMETER:  7.4747E-02  3.5551E-01 -2.5962E-01 -1.0577E-01  2.4572E-02  1.0532E-01 -1.2504E-01 -4.7477E+00  2.4176E-01  1.7999E-03
             1.1780E-01
 GRADIENT:   1.1679E-01  1.1339E-01  1.1381E-01 -2.6626E-01 -4.0466E-01  1.6888E-02  6.6053E-02  0.0000E+00  4.6747E-02 -6.5832E-02
            -2.2652E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1259
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1284E-04 -2.3855E-02 -3.2889E-04  1.1884E-02 -2.4990E-02
 SE:             2.9834E-02  1.9855E-02  1.4440E-04  2.5069E-02  2.3629E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8094E-01  2.2958E-01  2.2750E-02  6.3547E-01  2.9024E-01

 ETASHRINKSD(%)  5.2856E-02  3.3483E+01  9.9516E+01  1.6014E+01  2.0839E+01
 ETASHRINKVR(%)  1.0568E-01  5.5754E+01  9.9998E+01  2.9464E+01  3.7336E+01
 EBVSHRINKSD(%)  4.4941E-01  3.2985E+01  9.9558E+01  1.6371E+01  2.0031E+01
 EBVSHRINKVR(%)  8.9679E-01  5.5090E+01  9.9998E+01  3.0062E+01  3.6049E+01
 RELATIVEINF(%)  9.8978E+01  1.8753E+00  2.1230E-04  3.8600E+00  7.5143E+00
 EPSSHRINKSD(%)  4.3936E+01
 EPSSHRINKVR(%)  6.8569E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1660.0370639792138     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -924.88623741547565     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1660.037       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.29E+00  6.98E-01  8.14E-01  9.27E-01  1.01E+00  7.98E-01  1.00E-02  1.15E+00  9.06E-01  1.02E+00
 


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
+        1.15E+03
 
 TH 2
+       -6.44E+00  4.94E+02
 
 TH 3
+        1.08E+01  2.44E+02  4.80E+02
 
 TH 4
+       -1.32E+01  3.87E+02 -2.15E+02  8.87E+02
 
 TH 5
+       -2.41E+00 -3.87E+02 -5.65E+02  2.30E+02  1.00E+03
 
 TH 6
+        1.74E+00 -1.00E+00  2.91E+00 -3.97E+00 -1.28E+00  1.94E+02
 
 TH 7
+        8.64E-01  8.47E+00  1.80E+00 -1.01E+01 -1.76E+01 -5.21E-02  6.44E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.10E+00 -2.76E+01 -2.57E+01  4.02E+01 -6.52E-02 -3.62E-01  2.51E+01  0.00E+00  7.71E+01
 
 TH10
+       -7.31E-01 -1.17E+01 -5.28E+01 -1.11E+01 -5.69E+01 -4.02E-02  2.40E+01  0.00E+00  5.59E+00  1.01E+02
 
 TH11
+       -8.11E+00 -2.01E+01 -3.31E+01 -1.54E+00  4.47E+00  2.45E+00  8.99E+00  0.00E+00  7.79E+00  2.13E+01  2.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.716
Stop Time:
Wed Sep 29 14:11:17 CDT 2021
