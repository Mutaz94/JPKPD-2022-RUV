Sat Sep 18 12:45:40 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1485.40096171866        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.5910E+01 -8.7123E+01 -3.1847E+01 -4.5867E+01  1.0957E+02 -3.6513E+01 -5.2377E+01 -7.6789E+00 -3.2575E+01 -3.6470E+01
            -1.9322E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1515.17354372860        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.4724E-01  1.2585E+00  1.0855E+00  1.0011E+00  9.4399E-01  1.0650E+00  1.5267E+00  1.0088E+00  1.1297E+00  1.0906E+00
             1.5654E+00
 PARAMETER:  4.5795E-02  3.2992E-01  1.8201E-01  1.0114E-01  4.2360E-02  1.6298E-01  5.2312E-01  1.0879E-01  2.2193E-01  1.8671E-01
             5.4811E-01
 GRADIENT:  -6.7930E+01  1.0411E+02  3.6788E+01  5.8811E+01 -9.3297E+01 -9.1752E+00  1.5274E+01 -2.0428E+00  9.1932E+00  6.6686E+00
             4.8131E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1525.58617906578        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.5824E-01  1.1473E+00  1.2405E+00  9.7176E-01  1.0502E+00  1.0636E+00  1.2929E+00  1.6126E+00  9.4359E-01  9.0905E-01
             1.5003E+00
 PARAMETER:  5.7338E-02  2.3744E-01  3.1550E-01  7.1358E-02  1.4898E-01  1.6168E-01  3.5688E-01  5.7783E-01  4.1939E-02  4.6399E-03
             5.0569E-01
 GRADIENT:  -4.1067E+01  9.8198E+00  9.2127E+00  3.9906E+00  9.6425E+00 -7.1799E+00 -1.1334E+01 -7.4556E-01 -1.4800E+01 -1.5785E+01
             2.3355E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1532.30893685314        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.7922E-01  1.0958E+00  7.8177E-01  9.6964E-01  8.6450E-01  1.0833E+00  1.4464E+00  7.4347E-01  1.0193E+00  9.1013E-01
             1.3678E+00
 PARAMETER:  7.9006E-02  1.9148E-01 -1.4620E-01  6.9171E-02 -4.5610E-02  1.8000E-01  4.6910E-01 -1.9642E-01  1.1916E-01  5.8340E-03
             4.1319E-01
 GRADIENT:   2.6931E+00 -9.8190E-02 -4.2768E+00  5.5534E+00  4.4970E+00  5.4349E-01  2.6572E+00  1.6798E+00  1.9296E+00  1.2358E+00
             5.4359E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1532.64531380031        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.7952E-01  1.1424E+00  5.8382E-01  9.1200E-01  7.7169E-01  1.0846E+00  1.3905E+00  3.5012E-01  1.0039E+00  7.7075E-01
             1.3654E+00
 PARAMETER:  7.9309E-02  2.3315E-01 -4.3817E-01  7.8824E-03 -1.5917E-01  1.8121E-01  4.2967E-01 -9.4948E-01  1.0386E-01 -1.6039E-01
             4.1144E-01
 GRADIENT:   2.1800E-01 -1.2664E+00 -1.4449E+00 -1.1727E+00 -4.5719E-01  4.9985E-01 -5.4471E-01  6.8719E-01  6.4189E-02  7.8760E-01
             3.9643E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1532.77890909876        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.7970E-01  1.1128E+00  5.5131E-01  9.2318E-01  7.3757E-01  1.0839E+00  1.4249E+00  7.6755E-02  9.8981E-01  7.2255E-01
             1.3651E+00
 PARAMETER:  7.9493E-02  2.0686E-01 -4.9546E-01  2.0072E-02 -2.0440E-01  1.8059E-01  4.5408E-01 -2.4671E+00  8.9754E-02 -2.2497E-01
             4.1121E-01
 GRADIENT:   8.9781E-02 -1.8533E-01  4.3775E-01 -1.6635E-01 -1.1405E+00  5.4377E-02 -4.9828E-01  3.5262E-02  6.0111E-01 -1.5011E-01
            -6.0408E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1533.16091569376        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      494
 NPARAMETR:  9.8616E-01  1.1568E+00  6.2573E-01  9.1426E-01  8.0436E-01  1.0897E+00  1.3912E+00  1.0000E-02  1.0258E+00  8.3799E-01
             1.3662E+00
 PARAMETER:  8.6066E-02  2.4563E-01 -3.6883E-01  1.0360E-02 -1.1771E-01  1.8593E-01  4.3018E-01 -6.7626E+00  1.2550E-01 -7.6745E-02
             4.1202E-01
 GRADIENT:  -9.1199E+00 -1.8233E+00  2.0437E+00 -2.1510E+00 -3.5404E+00 -2.6173E+00 -9.8815E-01  0.0000E+00  9.3859E-02  1.3601E-01
            -3.7345E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1533.30168562893        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  9.9107E-01  1.2947E+00  6.0325E-01  8.3544E-01  8.5964E-01  1.0968E+00  1.2756E+00  1.0000E-02  1.0946E+00  8.8061E-01
             1.3671E+00
 PARAMETER:  9.1031E-02  3.5832E-01 -4.0542E-01 -7.9799E-02 -5.1239E-02  1.9241E-01  3.4344E-01 -8.8169E+00  1.9038E-01 -2.7143E-02
             4.1270E-01
 GRADIENT:  -2.4303E-01  3.3761E-01  2.4652E-02  3.2105E-01 -1.0456E-01 -1.0730E-01 -3.0486E-02  0.0000E+00 -4.2439E-02  1.2672E-04
            -4.1348E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1533.30184306785        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      761
 NPARAMETR:  9.9120E-01  1.2978E+00  6.0275E-01  8.3327E-01  8.6110E-01  1.0971E+00  1.2732E+00  1.0000E-02  1.0969E+00  8.8169E-01
             1.3672E+00
 PARAMETER:  9.1159E-02  3.6066E-01 -4.0626E-01 -8.2399E-02 -4.9547E-02  1.9265E-01  3.4150E-01 -8.8774E+00  1.9246E-01 -2.5912E-02
             4.1276E-01
 GRADIENT:  -1.8849E-03  2.4287E-02  2.5899E-02 -1.8861E-02 -4.7140E-02 -1.1927E-02  1.7429E-03  0.0000E+00  4.9361E-05  2.4029E-04
            -6.1874E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      761
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8251E-04 -7.3465E-03 -3.4298E-04  2.5849E-03 -1.9282E-02
 SE:             2.9709E-02  2.4347E-02  1.3066E-04  2.2636E-02  2.0703E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9510E-01  7.6285E-01  8.6636E-03  9.0908E-01  3.5166E-01

 ETASHRINKSD(%)  4.7028E-01  1.8433E+01  9.9562E+01  2.4168E+01  3.0642E+01
 ETASHRINKVR(%)  9.3834E-01  3.3469E+01  9.9998E+01  4.2495E+01  5.1895E+01
 EBVSHRINKSD(%)  7.0559E-01  1.7997E+01  9.9605E+01  2.4896E+01  3.0174E+01
 EBVSHRINKVR(%)  1.4062E+00  3.2755E+01  9.9998E+01  4.3594E+01  5.1243E+01
 RELATIVEINF(%)  9.8493E+01  5.9684E+00  1.6682E-04  4.7324E+00  6.0411E+00
 EPSSHRINKSD(%)  4.2142E+01
 EPSSHRINKVR(%)  6.6525E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1533.3018430678476     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -798.15101650410941     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1533.302       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  1.30E+00  6.03E-01  8.33E-01  8.61E-01  1.10E+00  1.27E+00  1.00E-02  1.10E+00  8.82E-01  1.37E+00
 


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
+        9.26E+02
 
 TH 2
+       -5.66E+00  2.83E+02
 
 TH 3
+        1.10E+01  1.50E+02  5.25E+02
 
 TH 4
+       -1.44E+01  2.11E+02 -3.42E+02  7.70E+02
 
 TH 5
+       -5.61E+00 -2.40E+02 -5.93E+02  3.74E+02  9.38E+02
 
 TH 6
+        9.36E-01 -1.64E+00  1.84E+00 -5.41E+00 -2.78E+00  1.58E+02
 
 TH 7
+        4.21E-01  2.07E+01 -2.03E+01 -1.63E+01  8.34E+00 -3.59E-01  5.89E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.62E+00 -1.44E+01 -3.38E+01  3.90E+01 -3.76E+00 -8.76E-02  1.06E+01  0.00E+00  6.43E+01
 
 TH10
+       -1.06E+00 -1.09E+01 -4.02E+01 -1.36E+01 -6.07E+01  9.44E-01  1.19E+01  0.00E+00  9.69E+00  6.50E+01
 
 TH11
+       -7.97E+00 -1.04E+01 -2.32E+01 -1.96E+00 -2.70E+00  3.88E+00  4.99E+00  0.00E+00  7.72E+00  1.63E+01  1.18E+02
 
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
 #CPUT: Total CPU Time in Seconds,       13.546
Stop Time:
Sat Sep 18 12:45:55 CDT 2021
