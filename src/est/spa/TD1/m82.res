Wed Sep 29 18:31:02 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat82.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m82.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1621.38373910727        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2874E+02 -2.9268E+01 -6.5045E+01  8.0220E+01  1.3535E+02  6.5604E+01 -1.4538E+01  4.0024E+00  9.3816E+00 -2.3775E+01
            -5.6218E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1633.48644469423        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0014E+00  1.0605E+00  1.0546E+00  9.7407E-01  9.5788E-01  9.3806E-01  1.0524E+00  9.8753E-01  9.9034E-01  1.0342E+00
             1.1497E+00
 PARAMETER:  1.0136E-01  1.5875E-01  1.5312E-01  7.3730E-02  5.6966E-02  3.6055E-02  1.5104E-01  8.7453E-02  9.0288E-02  1.3367E-01
             2.3953E-01
 GRADIENT:   7.2335E+00 -9.9239E+00 -7.9088E+00  3.9162E+00 -8.2696E-01  1.5570E-01 -9.0453E+00  3.5455E+00  3.5834E+00 -2.2782E+00
             7.8566E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1634.93384913799        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.9784E-01  1.1427E+00  1.0989E+00  9.2755E-01  1.0218E+00  9.4108E-01  1.1600E+00  9.1553E-01  9.3607E-01  1.1188E+00
             1.1325E+00
 PARAMETER:  9.7839E-02  2.3336E-01  1.9434E-01  2.4796E-02  1.2156E-01  3.9272E-02  2.4840E-01  1.1751E-02  3.3939E-02  2.1222E-01
             2.2446E-01
 GRADIENT:  -9.8068E-01  1.8396E+00 -1.1837E+00  3.9634E+00  4.6402E+00  1.3395E+00 -1.6808E-01  9.0358E-01 -1.7836E+00  5.5383E-01
             5.4396E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1635.34797517781        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.9998E-01  1.2690E+00  8.5969E-01  8.4092E-01  9.7527E-01  9.3854E-01  1.0709E+00  4.7957E-01  1.0112E+00  1.0714E+00
             1.1246E+00
 PARAMETER:  9.9977E-02  3.3825E-01 -5.1179E-02 -7.3259E-02  7.4960E-02  3.6569E-02  1.6852E-01 -6.3487E-01  1.1113E-01  1.6900E-01
             2.1743E-01
 GRADIENT:   4.3359E-01  7.1733E+00  1.2657E+00  6.8659E+00 -6.3935E+00 -5.0596E-01 -1.2689E+00  2.9358E-01  4.4499E-01  7.0484E-01
            -1.2048E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1635.57539741728        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  1.0006E+00  1.4648E+00  7.2279E-01  7.0932E-01  1.0198E+00  9.4081E-01  9.7727E-01  1.8521E-01  1.1149E+00  1.0783E+00
             1.1276E+00
 PARAMETER:  1.0059E-01  4.8169E-01 -2.2463E-01 -2.4345E-01  1.1964E-01  3.8982E-02  7.7009E-02 -1.5863E+00  2.0881E-01  1.7541E-01
             2.2007E-01
 GRADIENT:  -4.8439E-02  2.8399E+00 -4.2964E-01  2.9242E+00  5.9616E-01  1.7533E-01 -2.6266E-02  6.0796E-02 -4.6989E-01 -2.3183E-01
             1.7680E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1635.60297825846        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  1.0007E+00  1.4824E+00  7.1313E-01  6.9512E-01  1.0266E+00  9.4052E-01  9.6667E-01  6.2262E-02  1.1371E+00  1.0839E+00
             1.1271E+00
 PARAMETER:  1.0074E-01  4.9363E-01 -2.3809E-01 -2.6367E-01  1.2627E-01  3.8679E-02  6.6105E-02 -2.6764E+00  2.2851E-01  1.8056E-01
             2.1964E-01
 GRADIENT:   3.5992E-01 -1.1559E+00 -3.4712E-01  3.4096E-01  9.1102E-01  5.9963E-02  1.1145E-01  7.0469E-03  1.6025E-01 -1.9861E-02
             1.0248E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1635.60396414771        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1062
 NPARAMETR:  1.0005E+00  1.4852E+00  7.1173E-01  6.9356E-01  1.0264E+00  9.4030E-01  9.6515E-01  3.4170E-02  1.1373E+00  1.0841E+00
             1.1269E+00
 PARAMETER:  1.0055E-01  4.9557E-01 -2.4006E-01 -2.6592E-01  1.2606E-01  3.8448E-02  6.4525E-02 -3.2764E+00  2.2867E-01  1.8073E-01
             2.1947E-01
 GRADIENT:  -1.2438E-01 -1.6755E-01  1.5933E-01  1.8000E-01 -3.0423E-01 -2.9171E-02  9.4539E-04  2.0884E-03 -6.8730E-02 -5.5046E-05
            -3.7350E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1635.60613526248        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  1.0011E+00  1.4840E+00  7.1159E-01  6.9352E-01  1.0263E+00  9.4062E-01  9.6520E-01  1.0000E-02  1.1378E+00  1.0839E+00
             1.1269E+00
 PARAMETER:  1.0112E-01  4.9471E-01 -2.4026E-01 -2.6598E-01  1.2595E-01  3.8788E-02  6.4576E-02 -5.2686E+00  2.2909E-01  1.8055E-01
             2.1943E-01
 GRADIENT:   1.3186E+00 -1.5005E+00  2.7654E-02 -4.9028E-01  5.6805E-02  1.0531E-01 -1.3295E-02  0.0000E+00  2.0224E-02  1.8903E-02
            -1.1733E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1635.60613526248        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1266
 NPARAMETR:  1.0011E+00  1.4840E+00  7.1159E-01  6.9352E-01  1.0263E+00  9.4062E-01  9.6520E-01  1.0000E-02  1.1378E+00  1.0839E+00
             1.1269E+00
 PARAMETER:  1.0112E-01  4.9471E-01 -2.4026E-01 -2.6598E-01  1.2595E-01  3.8788E-02  6.4576E-02 -5.2686E+00  2.2909E-01  1.8055E-01
             2.1943E-01
 GRADIENT:   1.3186E+00 -1.5005E+00  2.7654E-02 -4.9028E-01  5.6805E-02  1.0531E-01 -1.3295E-02  0.0000E+00  2.0224E-02  1.8903E-02
            -1.1733E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1266
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.7395E-04 -1.7674E-02 -2.7909E-04  1.1805E-02 -2.7248E-02
 SE:             2.9778E-02  2.3248E-02  1.0500E-04  2.1608E-02  2.3338E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8998E-01  4.4711E-01  7.8598E-03  5.8485E-01  2.4300E-01

 ETASHRINKSD(%)  2.3919E-01  2.2116E+01  9.9648E+01  2.7612E+01  2.1814E+01
 ETASHRINKVR(%)  4.7781E-01  3.9341E+01  9.9999E+01  4.7599E+01  3.8869E+01
 EBVSHRINKSD(%)  6.0281E-01  2.1434E+01  9.9693E+01  2.9794E+01  1.9573E+01
 EBVSHRINKVR(%)  1.2020E+00  3.8275E+01  9.9999E+01  5.0712E+01  3.5315E+01
 RELATIVEINF(%)  9.8534E+01  2.7955E+00  9.4119E-05  2.1247E+00  1.0905E+01
 EPSSHRINKSD(%)  4.3065E+01
 EPSSHRINKVR(%)  6.7584E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1635.6061352624788     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -900.45530869874062     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1635.606       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.48E+00  7.12E-01  6.94E-01  1.03E+00  9.41E-01  9.65E-01  1.00E-02  1.14E+00  1.08E+00  1.13E+00
 


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
+        1.24E+03
 
 TH 2
+       -7.38E+00  3.48E+02
 
 TH 3
+        1.39E+01  9.11E+01  2.30E+02
 
 TH 4
+       -2.08E+01  3.56E+02 -2.10E+02  8.95E+02
 
 TH 5
+       -5.23E+00 -1.61E+02 -2.84E+02  2.50E+02  5.82E+02
 
 TH 6
+       -9.80E-02 -1.32E+00  2.49E+00 -4.44E+00 -1.33E+00  2.20E+02
 
 TH 7
+        1.02E+00  1.64E+01  9.20E-01 -1.86E+01 -6.39E+00 -8.42E-01  8.77E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.09E+00 -1.45E+01 -2.74E+01  4.08E+01  6.85E+00 -6.48E-01  2.14E+01  0.00E+00  4.73E+01
 
 TH10
+       -8.25E-01 -1.17E+01 -2.87E+01 -2.00E+00 -5.19E+01  3.94E-01  7.04E+00  0.00E+00  7.84E+00  7.24E+01
 
 TH11
+       -8.69E+00 -1.59E+01 -2.56E+01  4.19E+00  2.80E+00  3.39E+00  8.16E+00  0.00E+00  8.02E+00  1.65E+01  1.70E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.506
Stop Time:
Wed Sep 29 18:31:27 CDT 2021
