Wed Sep 29 20:44:41 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/28/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat92.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA (0.09 FIX)x3
$SIGMA  1  FIX;        [P]
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
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

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   23044.0978558630        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   1.0435E+03  7.7950E+02 -5.7032E+02 -1.2840E+03 -2.2212E+03 -5.2290E+02 -4.3555E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -444.636369584046        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  1.0922E+00  1.3099E+00  2.3394E+00  1.5637E+00  6.2489E-01  9.1971E-01  1.5951E+01
 PARAMETER:  1.8822E-01  3.6993E-01  9.4989E-01  5.4708E-01 -3.7017E-01  1.6301E-02  2.8695E+00
 GRADIENT:  -1.0166E+02  1.0710E+02 -9.5592E+00  2.7190E+01 -3.8496E+00  3.5319E-01  5.7451E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -461.672390240941        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      110
 NPARAMETR:  1.1873E+00  1.1477E+00  4.4419E+00  1.2880E+00  5.6074E-01  4.9626E+00  1.5178E+01
 PARAMETER:  2.7171E-01  2.3775E-01  1.5911E+00  3.5306E-01 -4.7850E-01  1.7019E+00  2.8199E+00
 GRADIENT:   3.3344E+01 -1.3882E+01 -2.4240E+00  1.0196E+01  4.1160E+00  4.5181E-02  2.3647E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -463.052511993801        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  1.1394E+00  1.1251E+00  6.5100E+00  1.2317E+00  5.0457E-01  9.1876E+00  1.4840E+01
 PARAMETER:  2.3053E-01  2.1789E-01  1.9733E+00  3.0841E-01 -5.8404E-01  2.3179E+00  2.7973E+00
 GRADIENT:   7.0386E+00  6.4116E+00 -1.9678E+00  3.4026E-01  2.7798E+00  1.4016E-01  8.8691E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -463.339247821206        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      216
 NPARAMETR:  1.1104E+00  1.0942E+00  8.5231E+00  1.2061E+00  4.3486E-01  1.1637E+01  1.4697E+01
 PARAMETER:  2.0471E-01  1.9003E-01  2.2428E+00  2.8736E-01 -7.3273E-01  2.5541E+00  2.7876E+00
 GRADIENT:  -2.8893E+00  2.1852E+00 -1.4712E+00 -1.3393E+00  1.7908E+00  5.9875E-01  6.3411E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -463.442525869192        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      268
 NPARAMETR:  1.1003E+00  1.0786E+00  1.0748E+01  1.1952E+00  3.7522E-01  1.3621E+01  1.4627E+01
 PARAMETER:  1.9562E-01  1.7569E-01  2.4747E+00  2.7830E-01 -8.8026E-01  2.7116E+00  2.7829E+00
 GRADIENT:  -1.8310E+00 -4.5093E-01 -8.6698E-01 -1.1980E+00  7.3329E-01  7.5895E-01  1.7816E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -463.940291192931        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  1.1131E+00  1.0811E+00  1.3203E+01  1.1767E+00  2.7017E-01  1.5022E+01  1.4979E+01
 PARAMETER:  2.0716E-01  1.7802E-01  2.6805E+00  2.6270E-01 -1.2087E+00  2.8095E+00  2.8067E+00
 GRADIENT:  -3.5741E+00 -3.4360E+00 -7.8156E-01 -4.9166E+00 -5.6980E-01  4.8069E-01 -1.1933E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -464.258515187510        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      466
 NPARAMETR:  1.1379E+00  1.1119E+00  1.2534E+01  1.2079E+00  3.2328E-01  1.4903E+01  1.5244E+01
 PARAMETER:  2.2918E-01  2.0608E-01  2.6285E+00  2.8892E-01 -1.0292E+00  2.8016E+00  2.8242E+00
 GRADIENT:   4.6621E-01  3.2822E+00 -7.9994E-01 -8.7957E-01 -6.2120E-01  5.8485E-01 -4.3893E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -464.290261797994        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      558
 NPARAMETR:  1.1348E+00  1.1088E+00  1.2721E+01  1.2096E+00  3.5370E-01  1.4858E+01  1.5296E+01
 PARAMETER:  2.2642E-01  2.0331E-01  2.6432E+00  2.9031E-01 -9.3931E-01  2.7986E+00  2.8276E+00
 GRADIENT:  -3.0697E+00 -3.7380E+00 -7.6639E-01 -1.0852E-01  1.2284E-01  4.5018E-01  4.3986E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -464.410668367306        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      674
 NPARAMETR:  1.1371E+00  1.1129E+00  1.8422E+01  1.2142E+00  3.7257E-01  1.4642E+01  1.5247E+01
 PARAMETER:  2.2848E-01  2.0701E-01  3.0135E+00  2.9409E-01 -8.8734E-01  2.7839E+00  2.8244E+00
 GRADIENT:  -1.2654E+00  1.7879E-01  2.0428E-02 -2.4941E-01  1.1755E-01 -7.9413E-01  2.9211E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -464.498944524330        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      790
 NPARAMETR:  1.1373E+00  1.1146E+00  2.1444E+01  1.2129E+00  3.5826E-01  1.7325E+01  1.5252E+01
 PARAMETER:  2.2866E-01  2.0847E-01  3.1654E+00  2.9301E-01 -9.2650E-01  2.9522E+00  2.8247E+00
 GRADIENT:  -9.5693E-01  3.4261E+00 -2.0842E-01 -7.8568E-01 -1.5648E-01  2.9040E-01 -2.3771E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -464.508014525874        NO. OF FUNC. EVALS.:  37
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  1.1374E+00  1.1123E+00  2.2267E+01  1.2147E+00  3.6195E-01  1.7549E+01  1.5256E+01
 PARAMETER:  2.2878E-01  2.0643E-01  3.2031E+00  2.9446E-01 -9.1625E-01  2.9650E+00  2.8249E+00
 GRADIENT:   3.0439E-01 -1.8779E-01 -9.2129E-02 -3.9872E-02  1.0925E-02  1.4776E-01 -7.1414E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      827
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2137E-02 -1.3267E-02 -8.1349E-03
 SE:             2.6733E-02  8.9150E-03  3.5577E-03
 N:                     100         100         100

 P VAL.:         6.4981E-01  1.3671E-01  2.2222E-02

 ETASHRINKSD(%)  1.0442E+01  7.0134E+01  8.8081E+01
 ETASHRINKVR(%)  1.9794E+01  9.1080E+01  9.8579E+01
 EBVSHRINKSD(%)  9.9271E+00  7.1172E+01  9.1730E+01
 EBVSHRINKVR(%)  1.8869E+01  9.1690E+01  9.9316E+01
 RELATIVEINF(%)  3.0889E+01  3.3639E+00  2.8720E-01
 EPSSHRINKSD(%)  3.3095E+00
 EPSSHRINKVR(%)  6.5095E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -464.50801452587359     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       270.64281203786459     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     7.47
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     2.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -464.508       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.14E+00  1.11E+00  2.23E+01  1.21E+00  3.62E-01  1.75E+01  1.53E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 


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
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        4.82E+02
 
 TH 2
+       -2.43E+02  5.42E+02
 
 TH 3
+       -1.41E+00  8.52E-02  1.09E-02
 
 TH 4
+        1.25E+01 -6.68E+01  8.50E-01  9.12E+01
 
 TH 5
+       -1.10E+01 -5.25E+01 -1.61E-02 -4.40E+00  2.80E+01
 
 TH 6
+        6.30E-01  3.62E-01 -2.11E-02 -2.25E-01 -3.79E-02  5.91E-02
 
 TH 7
+       -9.83E+00 -1.22E+01  4.49E-03  2.43E+00  3.44E+00 -3.79E-02  2.28E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,        9.709
Stop Time:
Wed Sep 29 20:44:52 CDT 2021
