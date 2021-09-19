Sat Sep 18 13:07:25 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat89.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1593.90620235962        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0433E+01 -6.8140E+01 -2.8157E+01 -2.7744E+01  8.4615E+01  4.1342E+01 -9.3888E+00 -4.8990E+00  1.9273E+01 -6.6679E+01
            -8.0006E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1612.34036317072        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0237E+00  1.1690E+00  1.0674E+00  9.4122E-01  1.0964E+00  8.3895E-01  1.0369E+00  1.0489E+00  7.9100E-01  1.5387E+00
             1.1718E+00
 PARAMETER:  1.2339E-01  2.5612E-01  1.6526E-01  3.9424E-02  1.9199E-01 -7.5606E-02  1.3625E-01  1.4775E-01 -1.3446E-01  5.3097E-01
             2.5851E-01
 GRADIENT:   1.0438E+02  5.7826E+00 -1.3815E+01  1.1828E+01  3.2421E+00 -2.5238E+01 -2.1820E+00 -6.5102E+00 -1.4361E+01  2.0819E+01
            -2.1923E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1616.64769895880        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0139E+00  1.1828E+00  2.0283E+00  9.7043E-01  1.4979E+00  8.1359E-01  1.0573E+00  2.4996E+00  8.2374E-01  1.7572E+00
             1.1861E+00
 PARAMETER:  1.1385E-01  2.6788E-01  8.0719E-01  6.9979E-02  5.0404E-01 -1.0630E-01  1.5568E-01  1.0161E+00 -9.3895E-02  6.6373E-01
             2.7070E-01
 GRADIENT:   7.8168E+01  1.8885E+01 -1.7724E+01  3.3188E+01  4.8289E+01 -3.7095E+01  6.4764E+00  2.3336E+00 -1.1892E+00  1.6291E+00
             3.6815E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1623.03294703477        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8365E-01  1.1632E+00  2.1881E+00  9.6602E-01  1.3266E+00  8.9604E-01  6.9442E-01  2.5453E+00  1.0463E+00  1.6014E+00
             1.1535E+00
 PARAMETER:  8.3512E-02  2.5114E-01  8.8304E-01  6.5434E-02  3.8265E-01 -9.7665E-03 -2.6467E-01  1.0342E+00  1.4526E-01  5.7086E-01
             2.4278E-01
 GRADIENT:  -1.3275E+01  1.4648E+01 -2.2708E-01  1.5092E+01 -3.5519E+00  3.9646E+00  2.0318E+00  6.2221E-01  4.4252E+00  4.6326E-01
            -1.2024E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1623.30667373095        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  9.8914E-01  1.0868E+00  2.3974E+00  1.0091E+00  1.3287E+00  8.8442E-01  6.7174E-01  2.6686E+00  9.9740E-01  1.5988E+00
             1.1543E+00
 PARAMETER:  8.9081E-02  1.8327E-01  9.7440E-01  1.0907E-01  3.8420E-01 -2.2825E-02 -2.9788E-01  1.0815E+00  9.7400E-02  5.6923E-01
             2.4352E-01
 GRADIENT:  -2.7737E+01 -2.8825E+00 -9.4344E-01 -3.5247E+00 -4.7760E+00 -2.5148E+00  4.0090E-01  2.0198E+00 -7.8350E-01 -2.0146E+00
            -1.9313E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1623.61814518138        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      502
 NPARAMETR:  9.9904E-01  1.2174E+00  2.2741E+00  9.2344E-01  1.3658E+00  8.8933E-01  6.0046E-01  2.6893E+00  1.0911E+00  1.6345E+00
             1.1602E+00
 PARAMETER:  9.9044E-02  2.9669E-01  9.2157E-01  2.0347E-02  4.1170E-01 -1.7285E-02 -4.1006E-01  1.0893E+00  1.8716E-01  5.9135E-01
             2.4862E-01
 GRADIENT:  -2.0194E+00  1.5892E+00  5.4685E-01  1.9777E+00  7.6807E-01 -2.0543E-01 -2.3723E-01 -1.0316E+00 -4.9192E-02  1.4621E-01
             6.8365E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1623.67147388584        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  1.0003E+00  1.3474E+00  2.1759E+00  8.3732E-01  1.3884E+00  8.9007E-01  6.3170E-01  2.8358E+00  1.1518E+00  1.6451E+00
             1.1595E+00
 PARAMETER:  1.0030E-01  3.9821E-01  8.7743E-01 -7.7551E-02  4.2817E-01 -1.6453E-02 -3.5934E-01  1.1423E+00  2.4130E-01  5.9780E-01
             2.4797E-01
 GRADIENT:  -1.9323E-01  2.1621E+00  1.0198E+00  1.4022E+00 -6.8400E-01 -8.9649E-02  2.0428E-01 -6.0010E-01  1.5509E-01  2.4663E-01
             2.5736E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1623.70021709142        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      854
 NPARAMETR:  1.0005E+00  1.4094E+00  1.9832E+00  7.9327E-01  1.3878E+00  8.9004E-01  6.2210E-01  2.7923E+00  1.1992E+00  1.6396E+00
             1.1580E+00
 PARAMETER:  1.0053E-01  4.4318E-01  7.8470E-01 -1.3159E-01  4.2772E-01 -1.6484E-02 -3.7465E-01  1.1269E+00  2.8166E-01  5.9445E-01
             2.4673E-01
 GRADIENT:  -6.9426E-02 -1.4177E-01 -3.8848E-01  5.9562E-01 -1.5150E-01 -1.9808E-01  6.3947E-02  4.2077E-01  1.6906E-01 -7.5613E-02
            -2.8193E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1623.70548859440        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1029
 NPARAMETR:  1.0008E+00  1.4606E+00  1.9168E+00  7.5853E-01  1.3960E+00  8.9062E-01  6.2285E-01  2.8205E+00  1.2325E+00  1.6434E+00
             1.1588E+00
 PARAMETER:  1.0076E-01  4.7886E-01  7.5063E-01 -1.7638E-01  4.3359E-01 -1.5842E-02 -3.7346E-01  1.1369E+00  3.0902E-01  5.9674E-01
             2.4739E-01
 GRADIENT:   5.4934E-04  1.1468E-01  1.8384E-02  6.7845E-02  4.7954E-04 -3.1746E-03 -4.1242E-03 -2.4338E-02 -2.1731E-03  3.8798E-03
             3.3547E-03

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1623.70548960350        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1086
 NPARAMETR:  1.0008E+00  1.4614E+00  1.9152E+00  7.5798E-01  1.3960E+00  8.9062E-01  6.2287E-01  2.8211E+00  1.2330E+00  1.6433E+00
             1.1588E+00
 PARAMETER:  1.0076E-01  4.7937E-01  7.4982E-01 -1.7710E-01  4.3363E-01 -1.5841E-02 -3.7343E-01  1.1371E+00  3.0946E-01  5.9672E-01
             2.4737E-01
 GRADIENT:  -3.6223E-03  2.8465E-02  1.4253E-03  2.0276E-02 -1.5205E-04 -3.8868E-03 -1.5863E-03 -2.9852E-03  1.3716E-04  4.4584E-04
            -1.7267E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1086
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4656E-04 -3.7129E-02 -4.2338E-02  1.3693E-02 -5.3103E-02
 SE:             2.9767E-02  1.5579E-02  1.3948E-02  2.4905E-02  2.2326E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9607E-01  1.7159E-02  2.4027E-03  5.8244E-01  1.7384E-02

 ETASHRINKSD(%)  2.7734E-01  4.7808E+01  5.3271E+01  1.6566E+01  2.5204E+01
 ETASHRINKVR(%)  5.5392E-01  7.2760E+01  7.8164E+01  3.0387E+01  4.4055E+01
 EBVSHRINKSD(%)  6.7529E-01  4.7122E+01  6.0545E+01  1.8861E+01  1.8987E+01
 EBVSHRINKVR(%)  1.3460E+00  7.2039E+01  8.4433E+01  3.4164E+01  3.4369E+01
 RELATIVEINF(%)  9.8454E+01  1.3232E+00  5.9723E+00  3.1724E+00  3.9627E+01
 EPSSHRINKSD(%)  4.4291E+01
 EPSSHRINKVR(%)  6.8965E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1623.7054896034952     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -888.55466303975697     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1623.705       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.46E+00  1.92E+00  7.58E-01  1.40E+00  8.91E-01  6.23E-01  2.82E+00  1.23E+00  1.64E+00  1.16E+00
 


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
+        1.37E+03
 
 TH 2
+       -2.01E+01  3.99E+02
 
 TH 3
+       -6.47E-01  1.31E+01  1.05E+01
 
 TH 4
+       -1.88E+01  5.23E+02 -1.27E+01  8.40E+02
 
 TH 5
+       -2.47E+00 -4.03E+01 -1.08E+01 -1.00E+01  1.48E+02
 
 TH 6
+       -6.10E+00 -1.54E+00  6.85E-01 -8.51E-01 -2.62E-01  2.49E+02
 
 TH 7
+        4.55E-01 -4.00E+01  5.37E-01 -1.69E+01 -2.47E+00 -1.42E+00  5.09E+01
 
 TH 8
+        5.42E-01 -1.05E+01 -6.64E+00  5.36E+00 -3.07E+00  1.46E-01  3.35E-01  6.28E+00
 
 TH 9
+       -9.84E-01 -1.86E+01 -4.93E-01  3.86E+01  1.12E+00 -3.14E+00  4.18E+01  1.15E+00  5.95E+01
 
 TH10
+        1.21E+00  7.49E-02 -2.11E-01 -1.30E+00 -2.49E+01 -2.97E-01  8.40E-01 -1.41E-01  1.19E-01  3.80E+01
 
 TH11
+       -1.34E+01 -1.56E+01  1.08E+00 -1.86E+01 -7.71E-01  4.14E+00  7.63E+00 -1.60E+00  7.99E+00  7.14E+00  1.62E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.117
Stop Time:
Sat Sep 18 13:07:45 CDT 2021
