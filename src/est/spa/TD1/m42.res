Sat Sep 18 14:03:20 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1704.99155913337        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.4032E+01  1.5117E+01 -2.4560E+01  5.4547E+01  6.0459E+01 -3.7374E-01  4.1228E+00  1.5831E+00 -2.1850E+00  2.1860E+00
            -3.5259E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1710.84614123646        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  1.0560E+00  9.2532E-01  9.5473E-01  1.0148E+00  8.9598E-01  1.0005E+00  9.2795E-01  9.9026E-01  1.0355E+00  8.8867E-01
             1.0280E+00
 PARAMETER:  1.5448E-01  2.2386E-02  5.3669E-02  1.1465E-01 -9.8410E-03  1.0055E-01  2.5222E-02  9.0210E-02  1.3485E-01 -1.8032E-02
             1.2762E-01
 GRADIENT:   5.5805E+00  2.0474E+00  4.5384E+00 -2.6881E+00 -9.8149E+00  7.6435E-02  2.8546E+00  2.2569E+00  7.4145E+00  1.5988E+00
             8.7516E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1711.05011574480        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  1.0613E+00  8.1514E-01  9.3612E-01  1.0873E+00  8.5018E-01  1.0038E+00  8.2803E-01  8.0219E-01  9.7160E-01  8.9214E-01
             9.9455E-01
 PARAMETER:  1.5954E-01 -1.0440E-01  3.3993E-02  1.8373E-01 -6.2312E-02  1.0379E-01 -8.8710E-02 -1.2041E-01  7.1192E-02 -1.4137E-02
             9.4533E-02
 GRADIENT:   1.7744E+01  7.8327E+00  1.7169E+00  1.2414E+01 -1.5672E+00  1.3047E+00 -2.3892E+00 -1.6807E+00 -3.4902E+00 -2.7053E+00
            -6.9809E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1711.55212608848        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  1.0507E+00  9.6766E-01  7.8482E-01  9.7541E-01  8.4023E-01  1.0048E+00  9.9238E-01  6.6291E-01  1.0012E+00  8.4779E-01
             1.0148E+00
 PARAMETER:  1.4950E-01  6.7129E-02 -1.4230E-01  7.5106E-02 -7.4075E-02  1.0482E-01  9.2354E-02 -3.1112E-01  1.0116E-01 -6.5120E-02
             1.1471E-01
 GRADIENT:  -8.9821E+00 -3.2150E+00 -9.9583E-01 -4.1472E+00  6.0482E-01  1.0375E+00  3.3418E-01  1.0643E+00  1.8951E+00  5.6185E-01
             3.6402E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1711.89018969997        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0571E+00  1.0336E+00  6.6303E-01  9.2562E-01  8.0662E-01  1.0031E+00  1.0025E+00  3.4345E-01  9.9956E-01  8.1619E-01
             1.0059E+00
 PARAMETER:  1.5556E-01  1.3303E-01 -3.1093E-01  2.2713E-02 -1.1491E-01  1.0309E-01  1.0250E-01 -9.6872E-01  9.9561E-02 -1.0311E-01
             1.0592E-01
 GRADIENT:   1.3703E+00 -2.9669E+00 -4.0935E+00  4.4156E-01  5.4073E+00 -1.7907E-01  6.5762E-01  4.3953E-01  3.2796E-02  1.3378E-01
             4.9059E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.03338806870        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0571E+00  1.0905E+00  6.2456E-01  8.8764E-01  8.1078E-01  1.0037E+00  9.5893E-01  9.8001E-02  1.0284E+00  8.2018E-01
             1.0051E+00
 PARAMETER:  1.5555E-01  1.8666E-01 -3.7070E-01 -1.9190E-02 -1.0976E-01  1.0373E-01  5.8068E-02 -2.2228E+00  1.2805E-01 -9.8229E-02
             1.0513E-01
 GRADIENT:   6.3873E-01 -1.5650E+00 -1.8589E+00  2.4585E-02  2.1250E+00 -8.8253E-02 -2.0539E-01  3.5177E-02  1.7236E-01  3.3278E-01
             2.7212E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1712.05846090250        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  1.0567E+00  1.0617E+00  6.3496E-01  9.0615E-01  8.0250E-01  1.0038E+00  9.8262E-01  1.2545E-02  1.0119E+00  8.1980E-01
             1.0048E+00
 PARAMETER:  1.5519E-01  1.5986E-01 -3.5419E-01  1.4504E-03 -1.2003E-01  1.0377E-01  8.2462E-02 -4.2784E+00  1.1178E-01 -9.8698E-02
             1.0475E-01
 GRADIENT:   3.8401E-02  8.6447E-02  9.9876E-02  5.7654E-02  1.8508E-02 -2.4971E-02 -4.3223E-02  4.2550E-04 -4.1272E-02 -7.7419E-02
            -6.8396E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1712.05859241473        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  1.0567E+00  1.0607E+00  6.3480E-01  9.0663E-01  8.0195E-01  1.0038E+00  9.8351E-01  1.0000E-02  1.0114E+00  8.1966E-01
             1.0049E+00
 PARAMETER:  1.5517E-01  1.5893E-01 -3.5445E-01  1.9815E-03 -1.2071E-01  1.0383E-01  8.3375E-02 -4.5614E+00  1.1134E-01 -9.8868E-02
             1.0485E-01
 GRADIENT:  -7.4260E-03 -3.9131E-02 -3.5346E-04 -2.9176E-02  2.7352E-02 -1.4145E-03 -4.8880E-03  0.0000E+00  3.2120E-03 -5.8168E-03
             1.2628E-04

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1712.05859241473        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  1.0567E+00  1.0607E+00  6.3480E-01  9.0663E-01  8.0195E-01  1.0038E+00  9.8351E-01  1.0000E-02  1.0114E+00  8.1966E-01
             1.0049E+00
 PARAMETER:  1.5517E-01  1.5893E-01 -3.5445E-01  1.9815E-03 -1.2071E-01  1.0383E-01  8.3375E-02 -4.5614E+00  1.1134E-01 -9.8868E-02
             1.0485E-01
 GRADIENT:  -7.4260E-03 -3.9131E-02 -3.5346E-04 -2.9176E-02  2.7352E-02 -1.4145E-03 -4.8880E-03  0.0000E+00  3.2120E-03 -5.8168E-03
             1.2628E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1255
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.9650E-07 -1.2121E-02 -4.0153E-04  5.3447E-03 -1.7253E-02
 SE:             2.9851E-02  2.0686E-02  1.7805E-04  2.5484E-02  2.3086E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9998E-01  5.5791E-01  2.4123E-02  8.3388E-01  4.5488E-01

 ETASHRINKSD(%)  1.0000E-10  3.0699E+01  9.9404E+01  1.4627E+01  2.2658E+01
 ETASHRINKVR(%)  1.0000E-10  5.1974E+01  9.9996E+01  2.7114E+01  4.0182E+01
 EBVSHRINKSD(%)  4.2775E-01  3.0362E+01  9.9462E+01  1.4839E+01  2.2257E+01
 EBVSHRINKVR(%)  8.5368E-01  5.1506E+01  9.9997E+01  2.7477E+01  3.9560E+01
 RELATIVEINF(%)  9.8979E+01  2.3694E+00  2.9194E-04  4.8026E+00  5.0358E+00
 EPSSHRINKSD(%)  4.4502E+01
 EPSSHRINKVR(%)  6.9199E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1712.0585924147254     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -976.90776585098718     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.35
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1712.059       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.06E+00  6.35E-01  9.07E-01  8.02E-01  1.00E+00  9.84E-01  1.00E-02  1.01E+00  8.20E-01  1.00E+00
 


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
+        9.79E+02
 
 TH 2
+       -6.60E+00  5.16E+02
 
 TH 3
+        1.37E+01  3.34E+02  8.54E+02
 
 TH 4
+       -1.38E+01  3.63E+02 -3.52E+02  9.63E+02
 
 TH 5
+       -2.09E+00 -5.36E+02 -9.88E+02  3.53E+02  1.54E+03
 
 TH 6
+        2.26E-01 -3.18E-01  3.31E+00 -5.40E+00 -5.01E-02  1.95E+02
 
 TH 7
+        1.11E+00  2.13E+01 -9.94E+00 -7.68E+00 -1.63E+01 -3.81E-01  5.10E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.79E+00 -2.62E+01 -3.28E+01  3.50E+01 -2.35E+00  1.08E-01  2.23E+01  0.00E+00  1.06E+02
 
 TH10
+       -1.31E+00 -1.02E+01 -7.96E+01 -1.85E+01 -5.15E+01 -7.00E-01  2.12E+01  0.00E+00  5.68E+00  1.18E+02
 
 TH11
+       -6.02E+00 -1.46E+01 -3.74E+01 -4.18E+00  7.79E+00  2.52E+00  6.93E+00  0.00E+00  8.55E+00  1.91E+01  2.08E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.161
Stop Time:
Sat Sep 18 14:03:43 CDT 2021
