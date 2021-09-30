Wed Sep 29 19:06:01 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat53.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1681.15747278463        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5824E+02 -4.2944E+01 -3.0875E+01  1.2134E+01  3.4403E+01  3.5732E+01  1.4562E+01  7.1465E+00  5.9057E+01  2.1268E+01
             5.5403E-02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1696.15405385861        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.7794E-01  1.1000E+00  1.0621E+00  1.0030E+00  1.0391E+00  1.1028E+00  9.3787E-01  9.8002E-01  7.8923E-01  9.1248E-01
             1.0119E+00
 PARAMETER:  7.7694E-02  1.9533E-01  1.6024E-01  1.0300E-01  1.3834E-01  1.9781E-01  3.5853E-02  7.9820E-02 -1.3670E-01  8.4101E-03
             1.1181E-01
 GRADIENT:  -8.0457E-01  1.3562E+01 -1.9485E+00  2.8072E+01  1.0937E+01  1.0597E+01  4.1370E+00 -2.6987E+00  1.2287E+01  1.3664E+00
            -3.2955E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1696.99736938123        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.8266E-01  1.0795E+00  9.7321E-01  1.0051E+00  9.6777E-01  1.0360E+00  9.7041E-01  1.0710E+00  7.0374E-01  8.1522E-01
             1.0232E+00
 PARAMETER:  8.2513E-02  1.7646E-01  7.2845E-02  1.0506E-01  6.7240E-02  1.3535E-01  6.9967E-02  1.6856E-01 -2.5135E-01 -1.0430E-01
             1.2298E-01
 GRADIENT:   7.3356E+00  1.3141E+01  2.3026E+00  1.4068E+01 -1.3937E+01 -1.5024E+01  5.9794E-01  2.3160E+00  1.4737E+00  2.3151E+00
             2.6827E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1697.67131733925        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.8110E-01  1.1539E+00  7.3759E-01  9.4063E-01  8.8664E-01  1.0804E+00  9.5208E-01  7.2536E-01  7.0128E-01  7.4260E-01
             1.0163E+00
 PARAMETER:  8.0916E-02  2.4315E-01 -2.0437E-01  3.8793E-02 -2.0312E-02  1.7734E-01  5.0894E-02 -2.2109E-01 -2.5484E-01 -1.9760E-01
             1.1618E-01
 GRADIENT:  -3.1485E-03  9.1031E+00  2.1984E+00  7.9157E+00 -5.9095E+00  1.5370E+00  1.0237E+00  3.9391E-01  1.1791E-01  4.3592E-01
            -1.3828E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1697.96722669842        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  9.8286E-01  1.2814E+00  5.3418E-01  8.3995E-01  8.2840E-01  1.0739E+00  8.6823E-01  4.3902E-01  7.3216E-01  6.7017E-01
             1.0169E+00
 PARAMETER:  8.2716E-02  3.4794E-01 -5.2703E-01 -7.4415E-02 -8.8261E-02  1.7129E-01 -4.1295E-02 -7.2321E-01 -2.1175E-01 -3.0023E-01
             1.1675E-01
 GRADIENT:   8.7023E-01  2.8923E+00 -6.5382E-01  4.3048E+00 -4.9693E-01 -1.6066E+00 -1.9829E-01  3.3781E-01 -1.4802E-01 -2.0678E-01
             1.5308E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.05268663045        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      875
 NPARAMETR:  9.8231E-01  1.4036E+00  4.7711E-01  7.6016E-01  8.6045E-01  1.0786E+00  8.0656E-01  2.8428E-01  7.8608E-01  6.9426E-01
             1.0169E+00
 PARAMETER:  8.2147E-02  4.3901E-01 -6.4002E-01 -1.7422E-01 -5.0299E-02  1.7562E-01 -1.1498E-01 -1.1578E+00 -1.4070E-01 -2.6490E-01
             1.1678E-01
 GRADIENT:  -2.7922E-02  1.0733E+00 -1.6024E-01  1.0405E+00 -4.9289E-01  1.5308E-01  3.7111E-01  8.9538E-02  5.1047E-02  3.3114E-01
             6.2491E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.06997417796        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.8329E-01  1.4081E+00  4.7426E-01  7.5724E-01  8.6172E-01  1.0793E+00  8.0322E-01  1.7854E-01  7.9032E-01  6.9824E-01
             1.0173E+00
 PARAMETER:  8.3154E-02  4.4227E-01 -6.4600E-01 -1.7808E-01 -4.8827E-02  1.7628E-01 -1.1913E-01 -1.6230E+00 -1.3531E-01 -2.5919E-01
             1.1715E-01
 GRADIENT:   1.8737E+00  1.4785E+00  4.1907E-01  1.1906E+00 -1.8446E-01  4.0393E-01 -3.3690E-02  1.2737E-02 -4.7469E-02 -2.0495E-03
             4.2403E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1698.07897599948        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.8187E-01  1.4188E+00  4.6203E-01  7.4821E-01  8.6034E-01  1.0778E+00  7.9830E-01  1.5769E-02  7.9599E-01  6.9679E-01
             1.0172E+00
 PARAMETER:  8.1705E-02  4.4982E-01 -6.7211E-01 -1.9007E-01 -5.0424E-02  1.7494E-01 -1.2528E-01 -4.0497E+00 -1.2817E-01 -2.6127E-01
             1.1704E-01
 GRADIENT:  -7.8949E-01 -1.1813E+00 -1.3470E-01 -3.4121E-02  9.1463E-01 -1.3842E-01  1.8600E-01  1.7942E-04  1.1882E-01  1.9178E-01
             6.9545E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1698.08170658513        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1380
 NPARAMETR:  9.8267E-01  1.4190E+00  4.6186E-01  7.4838E-01  8.5994E-01  1.0796E+00  7.9820E-01  1.0000E-02  7.9600E-01  6.9581E-01
             1.0172E+00
 PARAMETER:  8.2523E-02  4.4997E-01 -6.7250E-01 -1.8985E-01 -5.0895E-02  1.7656E-01 -1.2540E-01 -4.8557E+00 -1.2816E-01 -2.6268E-01
             1.1701E-01
 GRADIENT:   7.5889E-01 -1.7298E-01  1.2347E-01  3.1694E-01  3.7888E-01  5.1505E-01  1.3140E-01  0.0000E+00  1.1374E-01  7.4584E-02
             1.1645E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1380
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.0598E-04 -1.1491E-02 -3.1468E-04  7.4989E-03 -1.8191E-02
 SE:             2.9899E-02  2.4612E-02  1.5470E-04  2.2504E-02  2.1107E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8917E-01  6.4057E-01  4.1938E-02  7.3896E-01  3.8878E-01

 ETASHRINKSD(%)  1.0000E-10  1.7547E+01  9.9482E+01  2.4609E+01  2.9289E+01
 ETASHRINKVR(%)  1.0000E-10  3.2016E+01  9.9997E+01  4.3162E+01  4.9999E+01
 EBVSHRINKSD(%)  3.5876E-01  1.7639E+01  9.9525E+01  2.5422E+01  2.8936E+01
 EBVSHRINKVR(%)  7.1624E-01  3.2166E+01  9.9998E+01  4.4381E+01  4.9499E+01
 RELATIVEINF(%)  9.9240E+01  2.7900E+00  1.0625E-04  1.8949E+00  3.7388E+00
 EPSSHRINKSD(%)  4.3287E+01
 EPSSHRINKVR(%)  6.7837E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1698.0817065851274     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -962.93088002138927     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1698.082       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.42E+00  4.62E-01  7.48E-01  8.60E-01  1.08E+00  7.98E-01  1.00E-02  7.96E-01  6.96E-01  1.02E+00
 


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
+        9.80E+02
 
 TH 2
+       -7.58E+00  5.98E+02
 
 TH 3
+        4.85E+00  3.14E+02  1.24E+03
 
 TH 4
+       -1.97E+01  4.75E+02 -9.06E+02  1.75E+03
 
 TH 5
+       -5.21E+00 -4.73E+02 -1.25E+03  7.84E+02  1.57E+03
 
 TH 6
+        5.56E-01 -1.46E+00  1.99E+00 -4.83E+00 -1.36E+00  1.69E+02
 
 TH 7
+        3.18E-01  2.23E+01 -4.89E+01 -1.38E+01  9.48E+00 -2.83E-01  1.48E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.01E+00 -2.08E+01 -5.62E+01  6.44E+01 -4.71E+00 -2.88E-01  2.63E+01  0.00E+00  1.14E+02
 
 TH10
+       -7.58E-01 -1.77E+01 -6.05E+01 -1.84E+01 -8.57E+01 -7.97E-02  2.66E+01  0.00E+00  1.88E+01  1.19E+02
 
 TH11
+       -5.01E+00 -1.64E+01 -2.87E+01 -2.29E+00 -6.63E+00  1.97E+00  1.13E+01  0.00E+00  1.67E+01  2.42E+01  2.07E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       22.072
Stop Time:
Wed Sep 29 19:06:25 CDT 2021
