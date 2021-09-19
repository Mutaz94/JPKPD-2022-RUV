Sat Sep 18 00:50:45 CDT 2021
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
$DATA ../../../../data/int/A2/dat52.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m52.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2639.42339914066        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.7402E+01 -1.5199E+01  1.2130E+02 -1.1609E+02  4.3434E+01  9.4157E+00 -1.1261E+02 -5.2453E+01 -2.0500E+01 -2.9890E+01
            -2.2498E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3226.84663048880        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8043E-01  1.0902E+00  1.0055E+00  1.0436E+00  1.0624E+00  9.9266E-01  1.0612E+00  7.5170E-01  1.0162E+00  8.0479E-01
             1.8002E+00
 PARAMETER:  8.0232E-02  1.8638E-01  1.0546E-01  1.4266E-01  1.6054E-01  9.2633E-02  1.5935E-01 -1.8542E-01  1.1605E-01 -1.1718E-01
             6.8790E-01
 GRADIENT:  -2.1252E+01 -7.4613E+00  9.9355E+00  1.5206E+01  3.3126E+01  8.3561E+00 -2.6773E+00  3.8937E+00  4.1750E+00 -1.2796E+01
            -2.9260E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3229.78229941971        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.9649E-01  8.4306E-01  7.7778E-01  1.1865E+00  8.0605E-01  9.6040E-01  1.0835E+00  2.8194E-01  9.4480E-01  8.4238E-01
             1.7864E+00
 PARAMETER:  9.6481E-02 -7.0716E-02 -1.5131E-01  2.7097E-01 -1.1561E-01  5.9596E-02  1.8020E-01 -1.1661E+00  4.3217E-02 -7.1526E-02
             6.8023E-01
 GRADIENT:   1.5722E+01 -1.2460E+01 -1.7445E+01  5.4592E+01  3.5908E+01 -4.4751E+00 -7.6689E-01  1.1206E+00 -1.3714E+01  6.3406E+00
            -3.1514E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3231.54208759087        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.9029E-01  9.0575E-01  8.1253E-01  1.1292E+00  8.4745E-01  9.7026E-01  1.0891E+00  2.8486E-01  9.9358E-01  8.1283E-01
             1.8154E+00
 PARAMETER:  9.0241E-02  1.0034E-03 -1.0761E-01  2.2147E-01 -6.5527E-02  6.9807E-02  1.8538E-01 -1.1558E+00  9.3559E-02 -1.0723E-01
             6.9628E-01
 GRADIENT:   6.6814E-01 -6.9371E-01 -1.7868E+00 -3.4150E+00  2.2060E-01 -1.3035E-01  2.5185E-01  5.6058E-01  4.1888E-01 -2.3688E-01
             1.3542E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3231.67003280946        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.9004E-01  8.8761E-01  7.9938E-01  1.1396E+00  8.3000E-01  9.7086E-01  1.0830E+00  1.0466E-01  9.9147E-01  8.2656E-01
             1.8151E+00
 PARAMETER:  8.9989E-02 -1.9219E-02 -1.2391E-01  2.3065E-01 -8.6329E-02  7.0427E-02  1.7971E-01 -2.1570E+00  9.1433E-02 -9.0477E-02
             6.9617E-01
 GRADIENT:   1.8093E-01 -1.5186E-01  1.4764E-01 -6.1613E-01 -5.2011E-01  8.2830E-02  9.4779E-02  2.9884E-02 -4.4479E-02  2.4709E-01
             4.5886E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3231.72972849824        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      418
 NPARAMETR:  9.9406E-01  8.9488E-01  8.0472E-01  1.1390E+00  8.3739E-01  9.7321E-01  1.0838E+00  1.8093E-02  9.9313E-01  8.3094E-01
             1.8172E+00
 PARAMETER:  9.4045E-02 -1.1061E-02 -1.1727E-01  2.3018E-01 -7.7462E-02  7.2844E-02  1.8051E-01 -3.9122E+00  9.3110E-02 -8.5193E-02
             6.9732E-01
 GRADIENT:  -2.3130E+00 -1.1089E+00 -3.1949E-01 -2.5010E+00 -1.8741E-01  1.7020E-01  1.8525E-01  1.7302E-04  2.0443E-01  7.4954E-02
            -7.1411E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3231.73503367069        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      594
 NPARAMETR:  9.9509E-01  9.0075E-01  8.0894E-01  1.1379E+00  8.4270E-01  9.7270E-01  1.0807E+00  1.5356E-02  9.9278E-01  8.3335E-01
             1.8179E+00
 PARAMETER:  9.5080E-02 -4.5223E-03 -1.1204E-01  2.2918E-01 -7.1143E-02  7.2324E-02  1.7763E-01 -4.0762E+00  9.2758E-02 -8.2304E-02
             6.9767E-01
 GRADIENT:  -1.0576E-02 -1.6159E-03  8.6420E-03  4.5621E-04 -7.0518E-03 -8.0704E-05 -8.4310E-04  2.1178E-05 -3.3847E-03 -6.2352E-04
            -3.6731E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3231.73503919608        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      771
 NPARAMETR:  9.9510E-01  9.0073E-01  8.0891E-01  1.1379E+00  8.4267E-01  9.7271E-01  1.0807E+00  1.0000E-02  9.9279E-01  8.3337E-01
             1.8179E+00
 PARAMETER:  9.5083E-02 -4.5507E-03 -1.1207E-01  2.2919E-01 -7.1174E-02  7.2326E-02  1.7762E-01 -4.9614E+00  9.2765E-02 -8.2280E-02
             6.9767E-01
 GRADIENT:  -2.7861E-03 -1.1913E-03  4.4023E-03 -3.1386E-04 -2.9591E-03  5.7265E-04 -1.3473E-04  0.0000E+00 -1.6471E-03 -1.0236E-03
            -3.7543E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -3231.73503933394        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:      834
 NPARAMETR:  9.9512E-01  9.0074E-01  8.0889E-01  1.1379E+00  8.4269E-01  9.7270E-01  1.0807E+00  1.0000E-02  9.9285E-01  8.3346E-01
             1.8179E+00
 PARAMETER:  9.5084E-02 -4.5437E-03 -1.1207E-01  2.2919E-01 -7.1168E-02  7.2325E-02  1.7762E-01 -4.9666E+00  9.2767E-02 -8.2275E-02
             6.9767E-01
 GRADIENT:  -5.7351E-03 -9.2014E-04  3.1925E-03  1.1216E-03 -1.3967E-03  4.9067E-06 -2.8892E-04  0.0000E+00 -1.6617E-03 -1.4748E-03
            -2.0391E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      834
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0170E-03 -1.3487E-02 -1.1417E-04  3.7652E-03 -1.3493E-02
 SE:             2.9636E-02  2.2210E-02  2.2196E-04  2.7512E-02  2.3394E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7263E-01  5.4368E-01  6.0700E-01  8.9114E-01  5.6410E-01

 ETASHRINKSD(%)  7.1531E-01  2.5594E+01  9.9256E+01  7.8318E+00  2.1628E+01
 ETASHRINKVR(%)  1.4255E+00  4.4637E+01  9.9994E+01  1.5050E+01  3.8578E+01
 EBVSHRINKSD(%)  9.1441E-01  2.5503E+01  9.9271E+01  8.0391E+00  2.2779E+01
 EBVSHRINKVR(%)  1.8205E+00  4.4501E+01  9.9995E+01  1.5432E+01  4.0370E+01
 RELATIVEINF(%)  9.8156E+01  1.1204E+01  1.5721E-03  3.8808E+01  7.3433E+00
 EPSSHRINKSD(%)  1.7841E+01
 EPSSHRINKVR(%)  3.2500E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3231.7350393339434     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1577.6456795655326     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3231.735       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  9.01E-01  8.09E-01  1.14E+00  8.43E-01  9.73E-01  1.08E+00  1.00E-02  9.93E-01  8.33E-01  1.82E+00
 


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
+        1.17E+03
 
 TH 2
+       -5.47E+00  6.54E+02
 
 TH 3
+       -5.39E-01  9.06E+01  8.15E+02
 
 TH 4
+       -7.88E+00  2.00E+02 -8.19E+01  7.36E+02
 
 TH 5
+       -2.27E+00 -6.16E+02 -7.09E+02  2.20E+02  1.32E+03
 
 TH 6
+        6.03E+00 -2.57E+00  2.64E+00 -4.81E+00 -2.83E+00  2.03E+02
 
 TH 7
+        3.37E-01  1.39E+01 -1.07E+01  5.26E+00 -4.05E+00 -1.14E+00  5.43E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.59E+00 -1.73E+01  2.18E+01  1.62E+01 -2.63E+00 -1.76E+00  1.19E+01  0.00E+00  1.44E+02
 
 TH10
+        7.16E-01 -5.85E+00 -3.68E+01 -5.58E+00 -9.57E+00  2.31E+00  2.77E+01  0.00E+00  2.73E+00  1.02E+02
 
 TH11
+       -1.30E+01 -2.10E+01 -1.56E+01 -1.03E+01 -5.97E+00  3.15E+00  4.05E+00  0.00E+00  7.30E+00  1.63E+01  3.47E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.009
Stop Time:
Sat Sep 18 00:51:17 CDT 2021
