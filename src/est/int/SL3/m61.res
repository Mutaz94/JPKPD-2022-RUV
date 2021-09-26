Sat Sep 25 02:30:01 CDT 2021
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
$DATA ../../../../data/int/SL3/dat61.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      987
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

 TOT. NO. OF OBS RECS:      887
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   493.010717469197        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.8998E+00 -2.1636E+01  1.4813E+02  7.9374E+01  1.0304E+02  1.5815E+01 -1.4941E+02 -2.4229E+02 -1.2760E+02 -2.6851E+01
            -8.1538E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2359.85373122740        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0961E+00  1.4339E+00  9.0367E-01  8.9898E-01  1.1790E+00  8.3774E-01  1.1258E+00  1.0366E+00  8.0684E-01  9.6218E-01
             5.2051E+00
 PARAMETER:  1.9178E-01  4.6042E-01 -1.2921E-03 -6.4901E-03  2.6464E-01 -7.7044E-02  2.1847E-01  1.3591E-01 -1.1463E-01  6.1447E-02
             1.7496E+00
 GRADIENT:   8.5251E+00  2.5456E+01 -1.9427E+01  4.3933E+01  1.5544E+01 -3.3710E+01  2.5205E+01  5.2866E+00  1.5549E+01  4.5175E+00
             7.6140E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2471.85001223232        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0771E+00  1.6177E+00  1.2104E+00  7.6211E-01  1.3915E+00  9.2716E-01  9.7037E-01  6.0743E+00  5.0978E-01  1.5033E+00
             4.2032E+00
 PARAMETER:  1.7428E-01  5.8098E-01  2.9093E-01 -1.7167E-01  4.3036E-01  2.4371E-02  6.9920E-02  1.9041E+00 -5.7378E-01  5.0769E-01
             1.5358E+00
 GRADIENT:   1.8562E+01  4.7803E+01 -2.1239E+01  9.7260E+01 -5.1404E+01 -2.0667E+00  3.2760E+00  6.4655E+01  1.2603E+00  2.7063E+01
             5.6264E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2632.29145873957        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0384E+00  1.4208E+00  4.5564E+00  8.0194E-01  1.8291E+00  9.4394E-01  1.0956E+00  3.1854E+00  4.1763E-01  1.5313E+00
             2.9832E+00
 PARAMETER:  1.3767E-01  4.5121E-01  1.6165E+00 -1.2072E-01  7.0381E-01  4.2310E-02  1.9133E-01  1.2586E+00 -7.7317E-01  5.2614E-01
             1.1930E+00
 GRADIENT:  -6.2346E+00 -4.5775E+01 -8.9486E+00 -2.2135E+01  1.3012E+01  1.1026E+00  9.1185E-01  6.1557E+00 -9.5165E-02 -1.2711E+01
             2.1779E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2640.34365318126        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0371E+00  1.2526E+00  3.8415E+01  9.6144E-01  2.1462E+00  9.3440E-01  1.1468E+00  3.2292E+00  5.4777E-01  1.8826E+00
             2.9468E+00
 PARAMETER:  1.3647E-01  3.2523E-01  3.7484E+00  6.0673E-02  8.6369E-01  3.2148E-02  2.3697E-01  1.2722E+00 -5.0190E-01  7.3264E-01
             1.1807E+00
 GRADIENT:  -1.2101E+01 -9.6181E+00 -1.5690E+00 -6.2742E+00  1.5237E+01 -4.0623E+00 -4.1509E+00 -4.1676E-01 -1.3102E+00  1.7309E+00
             8.6607E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2641.93574613977        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0420E+00  1.0615E+00  2.7749E+02  1.1072E+00  2.1296E+00  9.4461E-01  1.2401E+00  4.0110E+00  6.3404E-01  1.8901E+00
             2.9356E+00
 PARAMETER:  1.4116E-01  1.5967E-01  5.7258E+00  2.0186E-01  8.5595E-01  4.3013E-02  3.1518E-01  1.4890E+00 -3.5564E-01  7.3662E-01
             1.1769E+00
 GRADIENT:  -1.4337E-01 -1.3402E+00 -2.8943E-01  1.4140E+00  1.7921E+00 -1.5093E-01 -1.2206E+00 -2.5168E-02  3.0572E-01  1.6565E+00
             2.4816E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2642.23057576883        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      510
 NPARAMETR:  1.0467E+00  9.3251E-01  1.7495E+03  1.1965E+00  2.1408E+00  9.4668E-01  1.4032E+00  5.3905E+00  5.9579E-01  1.8890E+00
             2.9368E+00
 PARAMETER:  1.4568E-01  3.0122E-02  7.5671E+00  2.7940E-01  8.6119E-01  4.5205E-02  4.3877E-01  1.7846E+00 -4.1787E-01  7.3604E-01
             1.1773E+00
 GRADIENT:   2.6839E+00 -5.2290E-01 -5.6912E-02 -1.8668E+00  2.7059E-01  1.2757E-01  7.6148E-01 -1.6519E-03  1.6484E-01  2.3900E-01
             1.4759E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2642.29553975175        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      685
 NPARAMETR:  1.0457E+00  7.8663E-01  2.4671E+04  1.3011E+00  2.1425E+00  9.4717E-01  1.5392E+00  8.5232E+00  6.0387E-01  1.8893E+00
             2.9336E+00
 PARAMETER:  1.4465E-01 -1.3999E-01  1.0213E+01  3.6325E-01  8.6199E-01  4.5722E-02  5.3128E-01  2.2428E+00 -4.0440E-01  7.3622E-01
             1.1762E+00
 GRADIENT:  -2.4881E-01  7.5240E-01 -4.8711E-03  2.1324E+00  2.2760E-01  2.0658E-01 -2.7495E-01 -6.5630E-05 -8.9839E-02  1.7443E-01
            -1.1188E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2642.30153823878        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  1.0456E+00  7.3282E-01  8.0966E+04  1.3372E+00  2.1417E+00  9.4604E-01  1.6169E+00  1.0627E+01  5.9926E-01  1.8882E+00
             2.9348E+00
 PARAMETER:  1.4462E-01 -2.1085E-01  1.1402E+01  3.9055E-01  8.6158E-01  4.4533E-02  5.8050E-01  2.4634E+00 -4.1206E-01  7.3560E-01
             1.1766E+00
 GRADIENT:  -4.0846E+00  3.6787E-01 -1.5926E-03  1.1719E+00  3.2798E-03 -7.1122E-01 -1.2403E-01  1.5622E-05 -1.7587E-01  8.6727E-02
             2.6176E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -2642.30160360574        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  1.0457E+00  7.3313E-01  8.2833E+04  1.3374E+00  2.1410E+00  9.4614E-01  1.6178E+00  1.0678E+01  5.9995E-01  1.8878E+00
             2.9349E+00
 PARAMETER:  1.4471E-01 -2.1043E-01  1.1425E+01  3.9072E-01  8.6126E-01  4.4630E-02  5.8106E-01  2.4682E+00 -4.1090E-01  7.3543E-01
             1.1767E+00
 GRADIENT:   8.0672E-02  1.7902E-01 -2.3735E-03  2.2418E-01 -7.7915E-02 -1.9589E-01  9.1310E-02  8.3291E-03 -5.1378E-03 -1.0732E-02
            -3.7471E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1008
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0231E-03  1.3780E-02  9.0528E-06 -2.8189E-02 -1.6069E-02
 SE:             2.9252E-02  2.0831E-02  6.1588E-06  1.8484E-02  2.6649E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7210E-01  5.0827E-01  1.4159E-01  1.2725E-01  5.4652E-01

 ETASHRINKSD(%)  2.0035E+00  3.0215E+01  9.9979E+01  3.8075E+01  1.0722E+01
 ETASHRINKVR(%)  3.9668E+00  5.1300E+01  1.0000E+02  6.1653E+01  2.0294E+01
 EBVSHRINKSD(%)  2.0306E+00  3.1821E+01  9.9985E+01  3.6334E+01  8.8912E+00
 EBVSHRINKVR(%)  4.0201E+00  5.3517E+01  1.0000E+02  5.9466E+01  1.6992E+01
 RELATIVEINF(%)  9.5793E+01  2.2925E+00  1.2162E-06  1.9949E+00  4.7585E+01
 EPSSHRINKSD(%)  1.5186E+01
 EPSSHRINKVR(%)  2.8066E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          887
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1630.1969579050892     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2642.3016036057443     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1012.1046457006551     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.58
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2642.302       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  7.33E-01  8.28E+04  1.34E+00  2.14E+00  9.46E-01  1.62E+00  1.07E+01  6.00E-01  1.89E+00  2.93E+00
 


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
+        1.56E+03
 
 TH 2
+       -1.56E+02  9.84E+02
 
 TH 3
+        4.05E-05  2.24E-05 -1.24E-12
 
 TH 4
+       -6.14E+01  4.78E+02  1.43E-05  7.41E+02
 
 TH 5
+        7.62E+00 -6.40E+00 -5.86E-07 -6.87E+00  5.95E+01
 
 TH 6
+       -7.06E+02 -6.62E+02  1.45E-05 -1.08E+02 -1.31E+01  1.39E+03
 
 TH 7
+       -6.28E+01  5.73E+01  4.19E-06  1.70E+00 -1.79E+00 -1.33E+02  3.39E+01
 
 TH 8
+       -9.62E-02 -5.04E-01  4.30E-08 -1.78E-01 -4.13E-02  7.72E-01 -1.56E-01  1.18E-02
 
 TH 9
+       -1.87E+02  1.63E+02  1.68E-05 -4.42E+01 -7.66E+00 -3.40E+02  3.97E+01 -1.25E+00  2.19E+02
 
 TH10
+       -1.78E+01 -1.30E+01 -5.68E-07 -1.56E+01 -4.29E+00  2.11E+00 -3.56E+00 -6.10E-02  5.17E+00  4.13E+01
 
 TH11
+       -1.85E+01 -1.23E+01  3.82E-07 -2.46E+01  2.15E+00  2.00E+01  1.56E+00 -3.49E-02  7.11E+00  1.95E+00  1.38E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.505
Stop Time:
Sat Sep 25 02:30:40 CDT 2021
