Sat Sep 18 12:58:09 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat63.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1651.37134291568        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.9584E+01 -1.4147E+02 -6.8816E+01 -1.3397E+02  8.9731E+01 -2.9581E+01 -2.0497E+01  1.3109E+01 -2.2774E+01  5.7138E+00
            -3.2615E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.50102032969        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0846E+00  1.1338E+00  1.1613E+00  1.0859E+00  1.0011E+00  1.0628E+00  1.0289E+00  9.3202E-01  1.0202E+00  9.3344E-01
             1.0908E+00
 PARAMETER:  1.8125E-01  2.2555E-01  2.4956E-01  1.8242E-01  1.0113E-01  1.6086E-01  1.2850E-01  2.9597E-02  1.2004E-01  3.1123E-02
             1.8688E-01
 GRADIENT:   1.2191E+02  6.5385E+01  1.3570E+01  7.1347E+01 -3.6534E+01  4.6530E-01 -9.1327E+00  1.1912E+00 -1.1599E+01 -6.1595E+00
             1.0926E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1671.88804588973        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0565E+00  1.1425E+00  9.3799E-01  1.0455E+00  9.4257E-01  1.0502E+00  1.0203E+00  5.0470E-01  1.0466E+00  9.7714E-01
             1.0721E+00
 PARAMETER:  1.5497E-01  2.3325E-01  3.5989E-02  1.4453E-01  4.0854E-02  1.4899E-01  1.2015E-01 -5.8379E-01  1.4559E-01  7.6879E-02
             1.6957E-01
 GRADIENT:   6.0413E+01  2.9008E+01 -9.6162E+00  4.7318E+01 -6.9606E+00 -3.1281E+00 -9.7937E+00  1.2685E+00 -3.3887E+00  8.5173E+00
            -2.2754E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1673.19129447220        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0309E+00  1.0741E+00  9.4993E-01  1.0548E+00  9.2719E-01  1.0555E+00  1.2077E+00  5.0524E-01  9.9653E-01  8.9820E-01
             1.0755E+00
 PARAMETER:  1.3043E-01  1.7148E-01  4.8631E-02  1.5337E-01  2.4404E-02  1.5399E-01  2.8870E-01 -5.8273E-01  9.6523E-02 -7.3679E-03
             1.7280E-01
 GRADIENT:   3.8352E+00  2.6549E+00 -4.2185E+00  6.8075E+00  1.2075E+00 -1.5744E+00 -1.0281E+00  1.2398E+00 -1.4386E+00  2.1199E+00
            -8.7299E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1673.29731917236        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0278E+00  1.0484E+00  8.7145E-01  1.0542E+00  8.7874E-01  1.0625E+00  1.2713E+00  3.0010E-01  9.8669E-01  8.3095E-01
             1.0806E+00
 PARAMETER:  1.2746E-01  1.4729E-01 -3.7596E-02  1.5278E-01 -2.9268E-02  1.6058E-01  3.4006E-01 -1.1036E+00  8.6603E-02 -8.5187E-02
             1.7756E-01
 GRADIENT:  -3.9615E+00 -3.3880E+00 -1.7706E+00 -3.5453E+00  1.3709E+00  7.4951E-01  6.0889E-01  5.2537E-01  1.2650E+00  2.1590E-01
             1.0156E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1673.58501375690        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0303E+00  1.0527E+00  8.3735E-01  1.0525E+00  8.6187E-01  1.0604E+00  1.2722E+00  7.1620E-02  9.7614E-01  8.1730E-01
             1.0789E+00
 PARAMETER:  1.2981E-01  1.5140E-01 -7.7512E-02  1.5113E-01 -4.8649E-02  1.5862E-01  3.4073E-01 -2.5364E+00  7.5846E-02 -1.0175E-01
             1.7591E-01
 GRADIENT:   6.9011E-01  9.3488E-01  5.0971E-01  1.0050E+00 -4.2750E-01 -5.7141E-02 -1.7890E-01  2.5046E-02 -1.3622E-01  5.6378E-02
            -2.5980E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1674.22885421821        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  1.0537E+00  1.0271E+00  8.3947E-01  1.0720E+00  8.5239E-01  1.0777E+00  1.3162E+00  3.8479E-02  9.6774E-01  8.0750E-01
             1.0832E+00
 PARAMETER:  1.5235E-01  1.2679E-01 -7.4985E-02  1.6949E-01 -5.9712E-02  1.7487E-01  3.7473E-01 -3.1576E+00  6.7212E-02 -1.1381E-01
             1.7995E-01
 GRADIENT:  -1.6613E+00 -8.9594E-01 -8.1321E-01 -2.1304E+00  1.5979E+00  1.3729E+00 -3.2632E-02  6.4434E-03  7.1068E-01 -2.4544E-01
             1.0149E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1674.40225862774        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  1.0537E+00  8.6175E-01  8.2990E-01  1.1677E+00  7.8030E-01  1.0737E+00  1.5340E+00  2.3419E-02  8.9480E-01  7.5556E-01
             1.0824E+00
 PARAMETER:  1.5234E-01 -4.8789E-02 -8.6454E-02  2.5507E-01 -1.4808E-01  1.7109E-01  5.2789E-01 -3.6542E+00 -1.1151E-02 -1.8030E-01
             1.7914E-01
 GRADIENT:  -4.9280E-01  1.0672E+00  2.0654E-01  1.1546E+00 -6.5939E-01  8.9020E-04 -9.7061E-02  3.0777E-03 -1.8111E-01  2.6088E-02
            -6.7353E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1674.41092425824        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      864
 NPARAMETR:  1.0538E+00  8.2629E-01  8.2522E-01  1.1855E+00  7.6533E-01  1.0733E+00  1.5896E+00  1.3967E-02  8.8318E-01  7.4287E-01
             1.0829E+00
 PARAMETER:  1.5237E-01 -9.0804E-02 -9.2107E-02  2.7019E-01 -1.6745E-01  1.7077E-01  5.6349E-01 -4.1711E+00 -2.4230E-02 -1.9723E-01
             1.7968E-01
 GRADIENT:  -8.3495E-02 -5.8148E-02  1.6820E-01 -3.5320E-01 -3.0139E-01 -9.1856E-02 -4.9599E-02  1.1658E-03 -3.0028E-02 -2.5654E-02
             4.9167E-03

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1674.41134385575        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  1.0538E+00  8.2808E-01  8.2632E-01  1.1849E+00  7.6664E-01  1.0736E+00  1.5870E+00  1.0000E-02  8.8394E-01  7.4444E-01
             1.0829E+00
 PARAMETER:  1.5242E-01 -8.8651E-02 -9.0772E-02  2.6965E-01 -1.6574E-01  1.7101E-01  5.6185E-01 -4.6608E+00 -2.3367E-02 -1.9512E-01
             1.7965E-01
 GRADIENT:   1.1461E-03 -5.0058E-04  2.6060E-02 -3.3530E-02 -3.3956E-02  1.0967E-03 -3.3662E-03  0.0000E+00 -1.6174E-03 -7.3524E-03
             1.7489E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      991
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6044E-05  1.1810E-02 -4.7150E-04 -1.2745E-02 -4.2050E-03
 SE:             2.9828E-02  2.1759E-02  2.0082E-04  2.4616E-02  2.1677E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9743E-01  5.8729E-01  1.8882E-02  6.0464E-01  8.4619E-01

 ETASHRINKSD(%)  7.0767E-02  2.7106E+01  9.9327E+01  1.7533E+01  2.7378E+01
 ETASHRINKVR(%)  1.4148E-01  4.6865E+01  9.9995E+01  3.1993E+01  4.7261E+01
 EBVSHRINKSD(%)  4.4134E-01  2.7261E+01  9.9355E+01  1.7216E+01  2.6465E+01
 EBVSHRINKVR(%)  8.8074E-01  4.7091E+01  9.9996E+01  3.1468E+01  4.5925E+01
 RELATIVEINF(%)  9.8796E+01  5.4635E+00  4.6631E-04  8.9769E+00  4.4636E+00
 EPSSHRINKSD(%)  4.2868E+01
 EPSSHRINKVR(%)  6.7359E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1674.4113438557499     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.26051729201174     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1674.411       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  8.28E-01  8.26E-01  1.18E+00  7.67E-01  1.07E+00  1.59E+00  1.00E-02  8.84E-01  7.44E-01  1.08E+00
 


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
+        8.60E+02
 
 TH 2
+       -6.31E+00  3.57E+02
 
 TH 3
+        1.58E+01  2.27E+02  7.34E+02
 
 TH 4
+       -7.91E+00  2.46E+02 -2.77E+02  7.09E+02
 
 TH 5
+       -5.58E+00 -4.44E+02 -1.04E+03  3.58E+02  1.83E+03
 
 TH 6
+        1.84E-01 -9.85E-01  2.89E+00 -2.11E+00 -2.91E+00  1.69E+02
 
 TH 7
+        1.15E+00  3.13E+01 -5.33E+00 -1.03E+01 -6.01E-02 -4.32E-01  2.75E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.51E+00 -1.82E+01 -3.35E+01  2.17E+01  9.46E+00 -7.52E-01  1.12E+01  0.00E+00  1.33E+02
 
 TH10
+       -2.14E+00 -9.23E+00 -7.52E+01 -3.16E+01 -3.62E+01  6.85E-01  1.37E+01  0.00E+00  1.21E+01  1.09E+02
 
 TH11
+       -6.04E+00 -9.62E+00 -3.40E+01 -8.48E+00  1.12E+01  2.45E+00  3.01E+00  0.00E+00  1.01E+01  2.46E+01  1.89E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.030
Stop Time:
Sat Sep 18 12:58:27 CDT 2021
