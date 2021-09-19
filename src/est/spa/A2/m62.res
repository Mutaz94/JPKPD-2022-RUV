Sat Sep 18 09:58:10 CDT 2021
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
$DATA ../../../../data/spa/A2/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -755.485808316535        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.2559E+01  4.7629E+01  7.5632E+01 -5.1563E+00  1.3748E+02  2.7653E+01 -1.5895E+01 -4.5870E+01 -5.9012E+01 -9.5328E+01
            -1.6272E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1375.11926490147        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0487E+00  9.5216E-01  9.0127E-01  1.0814E+00  8.7403E-01  7.9743E-01  8.9860E-01  1.0305E+00  9.6814E-01  9.1361E-01
             3.2766E+00
 PARAMETER:  1.4754E-01  5.0974E-02 -3.9516E-03  1.7822E-01 -3.4638E-02 -1.2636E-01 -6.9152E-03  1.3009E-01  6.7619E-02  9.6468E-03
             1.2868E+00
 GRADIENT:  -1.3802E+00  2.2202E+01 -3.0886E+00  3.6585E+01 -1.0285E+01 -2.7006E+01  7.6172E+00  4.9571E+00  7.7815E+00  1.7472E+01
             4.0701E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1383.02277896694        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0525E+00  8.1310E-01  5.5516E-01  1.1333E+00  6.1652E-01  8.7166E-01  6.4389E-01  6.7881E-01  1.0614E+00  5.3453E-01
             3.1683E+00
 PARAMETER:  1.5117E-01 -1.0690E-01 -4.8851E-01  2.2518E-01 -3.8367E-01 -3.7354E-02 -3.4023E-01 -2.8741E-01  1.5959E-01 -5.2636E-01
             1.2532E+00
 GRADIENT:  -1.5244E+00  2.4896E+01 -7.2169E-01  5.9368E+01 -1.2733E+01 -5.4136E+00  1.6061E+00  2.0401E+00  1.4974E+01  4.8432E+00
             3.7338E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1392.97236581229        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0364E+00  6.8890E-01  2.2476E-01  9.8982E-01  3.5498E-01  9.3020E-01  1.0324E+00  1.2434E+00  9.0089E-01  2.0660E-01
             2.5305E+00
 PARAMETER:  1.3574E-01 -2.7265E-01 -1.3927E+00  8.9763E-02 -9.3569E-01  2.7646E-02  1.3187E-01  3.1785E-01 -4.3709E-03 -1.4770E+00
             1.0284E+00
 GRADIENT:  -1.8237E+00  2.0415E+01  1.8287E+00  3.7396E+01 -1.4861E+01  7.4939E+00  3.6852E+00  3.3647E+00 -1.2249E+01  1.2115E+00
             4.4926E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1394.53617340926        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0316E+00  7.2059E-01  2.0358E-01  9.3445E-01  3.5392E-01  9.0871E-01  9.0196E-01  1.2533E+00  1.0140E+00  1.7984E-01
             2.5073E+00
 PARAMETER:  1.3107E-01 -2.2769E-01 -1.4917E+00  3.2202E-02 -9.3867E-01  4.2684E-03 -3.1854E-03  3.2579E-01  1.1392E-01 -1.6157E+00
             1.0192E+00
 GRADIENT:   1.4350E+00  3.4192E+00  3.1874E+00 -2.5600E+00 -8.8722E+00  3.2103E-01  3.6026E-02  2.6743E-02  5.3461E-01  5.4581E-01
             1.5752E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1397.44659224340        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.0218E+00  1.0322E+00  2.0146E-01  8.0919E-01  4.6922E-01  8.8632E-01  7.7382E-01  1.5152E+00  1.0825E+00  6.6852E-02
             2.5258E+00
 PARAMETER:  1.2154E-01  1.3171E-01 -1.5022E+00 -1.1173E-01 -6.5669E-01 -2.0674E-02 -1.5641E-01  5.1552E-01  1.7924E-01 -2.6053E+00
             1.0266E+00
 GRADIENT:  -6.3709E+00  4.8627E+00 -1.8255E-01  1.0871E+01 -3.8331E+00 -1.9770E+00  1.0136E+00  2.8131E-01 -2.6771E+00  1.5926E-01
            -5.8259E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1406.72816766252        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  1.0165E+00  1.3503E+00  1.7946E-01  6.4863E-01  6.0125E-01  8.7975E-01  6.8288E-01  2.2021E+00  1.3611E+00  1.6511E-02
             2.3517E+00
 PARAMETER:  1.1635E-01  4.0033E-01 -1.6178E+00 -3.3290E-01 -4.0874E-01 -2.8117E-02 -2.8143E-01  8.8939E-01  4.0831E-01 -4.0037E+00
             9.5513E-01
 GRADIENT:  -3.6844E+00  1.1819E+01  7.0706E-01  1.0272E+01 -4.8497E+00 -5.2153E-01 -1.1793E-01 -1.3784E+00 -8.6625E-01  5.8839E-03
            -1.1974E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1408.42985385056        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      565
 NPARAMETR:  1.0235E+00  1.6020E+00  1.5509E-01  5.1521E-01  7.1940E-01  8.8302E-01  6.4809E-01  2.4018E+00  1.5673E+00  1.0000E-02
             2.5147E+00
 PARAMETER:  1.2324E-01  5.7127E-01 -1.7638E+00 -5.6318E-01 -2.2934E-01 -2.4410E-02 -3.3372E-01  9.7620E-01  5.4933E-01 -5.6922E+00
             1.0222E+00
 GRADIENT:  -4.0139E+00  1.7019E+01  1.4199E+00  3.9364E+00 -9.6094E+00  9.0205E-01  5.2703E-01 -4.0530E-01 -8.4180E-01  0.0000E+00
             1.3364E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1408.58418281238        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  1.0254E+00  1.6757E+00  1.4294E-01  4.6729E-01  7.6071E-01  8.8039E-01  6.3424E-01  2.4553E+00  1.6859E+00  1.0000E-02
             2.5423E+00
 PARAMETER:  1.2506E-01  6.1624E-01 -1.8453E+00 -6.6081E-01 -1.7350E-01 -2.7386E-02 -3.5533E-01  9.9823E-01  6.2231E-01 -6.3950E+00
             1.0331E+00
 GRADIENT:  -1.6220E-01  5.9609E-01  4.2170E-01 -4.9928E-01 -6.6515E-01  2.1649E-02 -3.7341E-02 -3.2387E-02  2.8600E-02  0.0000E+00
            -1.0331E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1408.58543635630        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  1.0253E+00  1.6741E+00  1.4223E-01  4.6810E-01  7.5941E-01  8.8037E-01  6.3455E-01  2.4473E+00  1.6829E+00  1.0000E-02
             2.5427E+00
 PARAMETER:  1.2503E-01  6.1525E-01 -1.8503E+00 -6.5908E-01 -1.7521E-01 -2.7409E-02 -3.5483E-01  9.9498E-01  6.2053E-01 -6.4005E+00
             1.0332E+00
 GRADIENT:   1.8934E-03  3.0620E-02  1.7629E-02 -1.4085E-02 -1.8617E-02  2.4005E-03 -7.8047E-03 -7.2664E-03  1.8018E-03  0.0000E+00
            -7.6802E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      867
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6361E-03 -2.5220E-02 -6.7574E-03  2.1060E-02 -4.2014E-04
 SE:             2.9062E-02  2.2731E-02  1.5478E-02  2.2607E-02  2.4275E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5510E-01  2.6720E-01  6.6242E-01  3.5156E-01  8.3494E-02

 ETASHRINKSD(%)  2.6384E+00  2.3850E+01  4.8147E+01  2.4263E+01  9.9187E+01
 ETASHRINKVR(%)  5.2073E+00  4.2011E+01  7.3113E+01  4.2639E+01  9.9993E+01
 EBVSHRINKSD(%)  2.8611E+00  2.3142E+01  4.8100E+01  2.4126E+01  9.9125E+01
 EBVSHRINKVR(%)  5.6403E+00  4.0928E+01  7.3064E+01  4.2431E+01  9.9992E+01
 RELATIVEINF(%)  9.1377E+01  6.3448E+00  1.2009E+01  1.1659E+01  9.6848E-04
 EPSSHRINKSD(%)  3.4347E+01
 EPSSHRINKVR(%)  5.6897E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1408.5854363562983     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -673.43460979256008     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.68
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1408.585       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.67E+00  1.42E-01  4.68E-01  7.59E-01  8.80E-01  6.35E-01  2.45E+00  1.68E+00  1.00E-02  2.54E+00
 


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
+        1.29E+03
 
 TH 2
+       -5.90E+01  5.89E+02
 
 TH 3
+       -1.09E+02  4.49E+02  2.46E+03
 
 TH 4
+       -6.86E+01  2.42E+02 -8.49E+02  1.02E+03
 
 TH 5
+       -2.54E+01 -6.94E+02 -1.06E+03  3.17E+02  1.46E+03
 
 TH 6
+       -8.24E-01 -1.21E+01  1.95E+00 -1.10E+01  5.85E+00  2.23E+02
 
 TH 7
+        1.48E+00 -2.42E+01 -1.93E+01  1.71E+01  3.65E+01  2.27E+00  1.58E+02
 
 TH 8
+        1.96E+00 -1.17E+00 -4.65E+01  1.12E+00 -1.70E+01  1.19E+00  1.50E+00  6.65E+00
 
 TH 9
+        3.81E+00 -1.16E+01  6.27E+00  2.40E+01 -8.24E+00  1.54E+00  1.37E+01 -1.64E-02  2.03E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.66E+01 -1.50E+01  2.00E+01 -2.77E+00 -2.05E+01  3.04E+00  2.02E+01  1.35E+00  8.97E+00  0.00E+00  3.97E+01
 
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
 #CPUT: Total CPU Time in Seconds,       16.983
Stop Time:
Sat Sep 18 09:58:29 CDT 2021
