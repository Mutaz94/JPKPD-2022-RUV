Sat Sep 25 13:09:39 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat86.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1697.25877960219        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.8857E+01 -5.6960E+01 -5.1814E+01 -8.7939E+00  8.2231E+01  9.0915E+00  7.0275E+00  7.9390E+00  1.9152E+01  2.0611E+00
            -1.8122E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1703.63641209298        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  1.0286E+00  1.0454E+00  1.0610E+00  9.9229E-01  9.8035E-01  9.7096E-01  9.5357E-01  9.6970E-01  9.2748E-01  9.6132E-01
             1.0612E+00
 PARAMETER:  1.2816E-01  1.4439E-01  1.5923E-01  9.2261E-02  8.0150E-02  7.0531E-02  5.2453E-02  6.9232E-02  2.4714E-02  6.0552E-02
             1.5941E-01
 GRADIENT:  -3.2761E+00 -7.4953E+00 -1.4460E+00 -1.0536E+01  3.8503E+00 -5.3740E+00  3.8055E+00  1.8606E+00 -4.7719E-01 -5.0102E-01
             3.6872E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1703.92768573031        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  1.0280E+00  1.0771E+00  1.0546E+00  9.7471E-01  9.9430E-01  9.7709E-01  8.3698E-01  8.7845E-01  9.8262E-01  9.9067E-01
             1.0545E+00
 PARAMETER:  1.2758E-01  1.7426E-01  1.5314E-01  7.4381E-02  9.4282E-02  7.6825E-02 -7.7951E-02 -2.9600E-02  8.2465E-02  9.0630E-02
             1.5304E-01
 GRADIENT:  -4.9603E+00 -6.7631E+00  6.8547E-01 -3.7565E+00  5.6886E+00 -2.9546E+00 -2.6395E-01 -9.8638E-01  1.8743E+00 -1.1078E+00
             5.9733E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1704.36950447903        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  1.0336E+00  1.2539E+00  8.1648E-01  8.6179E-01  9.5967E-01  9.9050E-01  8.2602E-01  6.8297E-01  1.0347E+00  9.3493E-01
             1.0484E+00
 PARAMETER:  1.3308E-01  3.2626E-01 -1.0275E-01 -4.8746E-02  5.8836E-02  9.0459E-02 -9.1139E-02 -2.8131E-01  1.3414E-01  3.2712E-02
             1.4723E-01
 GRADIENT:   2.8027E+00  9.4909E+00  1.9294E+00  7.9648E+00 -7.2321E+00  1.2699E+00  1.0371E+00 -5.4568E-03 -7.5174E-01  1.0698E+00
            -1.1320E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1704.85386497800        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  1.0326E+00  1.5565E+00  6.2395E-01  6.6386E-01  1.0267E+00  9.8791E-01  7.1363E-01  5.4513E-01  1.2530E+00  9.5678E-01
             1.0511E+00
 PARAMETER:  1.3205E-01  5.4244E-01 -3.7169E-01 -3.0969E-01  1.2632E-01  8.7833E-02 -2.3740E-01 -5.0673E-01  3.2552E-01  5.5815E-02
             1.4985E-01
 GRADIENT:  -1.9844E+00  1.6075E+01  3.3928E-01  9.5904E+00 -6.1498E+00 -5.4478E-01 -3.2126E-01  2.9831E-01  5.8928E-01  1.9088E+00
            -5.1294E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1705.01718798993        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      841
 NPARAMETR:  1.0334E+00  1.6675E+00  5.3294E-01  5.8045E-01  1.0512E+00  9.8930E-01  6.9500E-01  3.7075E-01  1.3490E+00  9.4953E-01
             1.0518E+00
 PARAMETER:  1.3283E-01  6.1133E-01 -5.2935E-01 -4.4394E-01  1.4996E-01  8.9246E-02 -2.6384E-01 -8.9224E-01  3.9936E-01  4.8216E-02
             1.5055E-01
 GRADIENT:  -5.2382E-01  7.0598E-01 -8.5384E-01  6.7363E-01 -1.9318E-01 -1.7395E-01  4.6105E-01  1.9964E-01  1.8354E-01  6.0974E-01
             1.4825E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1705.05775894944        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1018
 NPARAMETR:  1.0337E+00  1.6803E+00  5.1573E-01  5.6921E-01  1.0502E+00  9.8984E-01  6.9458E-01  2.0025E-01  1.3646E+00  9.4611E-01
             1.0518E+00
 PARAMETER:  1.3310E-01  6.1897E-01 -5.6218E-01 -4.6351E-01  1.4894E-01  8.9792E-02 -2.6444E-01 -1.5082E+00  4.1085E-01  4.4607E-02
             1.5052E-01
 GRADIENT:   5.9532E-02 -2.1980E+00  1.2983E-02 -1.5300E+00 -7.6431E-01  1.0243E-02  7.3095E-01  3.8882E-02  3.8367E-01  4.8891E-01
             1.7007E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1705.07679026184        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1193
 NPARAMETR:  1.0336E+00  1.6735E+00  5.1220E-01  5.7463E-01  1.0432E+00  9.8980E-01  6.9410E-01  6.4031E-02  1.3526E+00  9.3869E-01
             1.0516E+00
 PARAMETER:  1.3308E-01  6.1490E-01 -5.6905E-01 -4.5404E-01  1.4230E-01  8.9743E-02 -2.6515E-01 -2.6484E+00  4.0202E-01  3.6729E-02
             1.5032E-01
 GRADIENT:  -9.7234E-02  2.2821E-01 -5.5637E-02  1.8140E-01  1.5128E-02 -3.4269E-02 -6.0018E-02  2.8302E-03 -3.0383E-02  8.2398E-02
            -3.5777E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1705.07816743049        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1357
 NPARAMETR:  1.0336E+00  1.6739E+00  5.1136E-01  5.7411E-01  1.0430E+00  9.8983E-01  6.9431E-01  1.3669E-02  1.3532E+00  9.3787E-01
             1.0517E+00
 PARAMETER:  1.3309E-01  6.1517E-01 -5.7068E-01 -4.5493E-01  1.4207E-01  8.9777E-02 -2.6484E-01 -4.1927E+00  4.0246E-01  3.5858E-02
             1.5041E-01
 GRADIENT:   5.2582E+01  9.3067E+01  1.8933E+01  2.3595E+01 -8.3669E+00  2.8307E+00 -2.0189E+00  2.3153E-01  1.5916E+01 -5.2713E+00
            -6.5911E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1705.07857446925        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1515
 NPARAMETR:  1.0336E+00  1.6739E+00  5.1130E-01  5.7412E-01  1.0430E+00  9.8988E-01  6.9436E-01  1.0000E-02  1.3531E+00  9.3774E-01
             1.0517E+00
 PARAMETER:  1.3307E-01  6.1515E-01 -5.7080E-01 -4.5492E-01  1.4209E-01  8.9727E-02 -2.6476E-01 -4.8423E+00  4.0238E-01  3.5720E-02
             1.5042E-01
 GRADIENT:   1.9980E-02  2.6170E-02  2.1931E-03 -9.0195E-03  4.4197E-02 -7.1246E+00 -3.8383E-03  0.0000E+00 -2.0906E-02 -6.0501E-03
             5.5313E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1515
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5428E-04 -3.1775E-02 -2.7103E-04  2.3846E-02 -3.3486E-02
 SE:             2.9830E-02  2.2345E-02  1.1121E-04  2.3547E-02  2.2762E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9320E-01  1.5502E-01  1.4802E-02  3.1121E-01  1.4126E-01

 ETASHRINKSD(%)  6.5658E-02  2.5142E+01  9.9627E+01  2.1113E+01  2.3744E+01
 ETASHRINKVR(%)  1.3127E-01  4.3963E+01  9.9999E+01  3.7769E+01  4.1850E+01
 EBVSHRINKSD(%)  4.8778E-01  2.4572E+01  9.9683E+01  2.2571E+01  2.2379E+01
 EBVSHRINKVR(%)  9.7317E-01  4.3106E+01  9.9999E+01  4.0047E+01  3.9750E+01
 RELATIVEINF(%)  9.8985E+01  3.4154E+00  1.1395E-04  3.8877E+00  1.0768E+01
 EPSSHRINKSD(%)  4.4015E+01
 EPSSHRINKVR(%)  6.8656E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1705.0785744692544     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -969.92774790551618     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1705.079       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.67E+00  5.11E-01  5.74E-01  1.04E+00  9.90E-01  6.94E-01  1.00E-02  1.35E+00  9.38E-01  1.05E+00
 


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
+        6.05E+03
 
 TH 2
+       -7.31E+00  5.81E+02
 
 TH 3
+        8.45E+00  2.02E+02  1.58E+03
 
 TH 4
+       -1.48E+01  4.07E+02 -3.26E+02  2.51E+03
 
 TH 5
+       -5.49E+00 -2.71E+02 -4.30E+02  2.85E+02  5.02E+03
 
 TH 6
+        2.17E+02 -1.82E+00  4.09E+00 -3.07E+02 -1.40E+01  9.90E+03
 
 TH 7
+        6.87E+00  2.70E+00 -5.85E+00 -6.60E+00 -1.48E+01  6.06E-01  2.96E+03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.93E-01 -2.16E+01 -3.39E+01  5.70E+01  8.86E-01 -4.30E+00  2.80E+01  0.00E+00  3.67E+02
 
 TH10
+        8.98E+00 -1.85E+01 -4.04E+01  1.91E+01 -5.79E+01  3.40E+01  1.90E+01  0.00E+00  1.29E+01  1.09E+04
 
 TH11
+       -6.20E+04 -2.04E+01 -2.40E+01  3.26E+00 -9.10E-01  2.02E+02  8.66E+00  0.00E+00  1.04E+01  1.68E+01  5.85E+05
 
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
 #CPUT: Total CPU Time in Seconds,       24.557
Stop Time:
Sat Sep 25 13:10:06 CDT 2021
