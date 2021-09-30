Wed Sep 29 13:40:04 CDT 2021
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
$DATA ../../../../data/spa/A3/dat60.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1189.78674624056        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3832E+02  9.1423E+01  1.0522E+02  1.6649E+01  2.6221E+02  4.1141E+01 -6.5725E+01 -6.9539E+01 -1.3573E+02 -1.8805E+02
            -5.1678E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -807.231134474175        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.2965E+00  1.0739E+00  8.4637E-01  1.4012E+00  6.6508E-01  6.7947E-01  1.0404E+00  1.1166E+00  1.1868E+00  1.2024E+00
             1.5779E+01
 PARAMETER:  3.5969E-01  1.7127E-01 -6.6800E-02  4.3734E-01 -3.0784E-01 -2.8645E-01  1.3961E-01  2.1029E-01  2.7126E-01  2.8432E-01
             2.8587E+00
 GRADIENT:  -6.2023E+01 -2.4964E+01 -7.2311E+00 -3.1033E+01 -5.2723E+00  1.1581E+01  9.0069E+00  5.2455E+00  3.1959E+01  1.9802E+01
             4.6235E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -916.901123148411        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0122E+00  2.7867E-01  5.9290E-02  1.3288E+00  8.1394E-02  6.7625E-01  4.7654E-01  1.9833E+00  1.7911E+00  1.4460E+00
             8.6125E+00
 PARAMETER:  1.1212E-01 -1.1777E+00 -2.7253E+00  3.8426E-01 -2.4085E+00 -2.9119E-01 -6.4120E-01  7.8475E-01  6.8284E-01  4.6877E-01
             2.2532E+00
 GRADIENT:  -1.2227E+02  8.7498E+01  1.0613E+02 -4.2110E+00 -2.5896E+02 -4.8544E+00  1.6959E+00  1.5638E+01  5.5887E+00  1.5616E+01
             2.0941E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1129.15053683707        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  8.5941E-01  1.7645E-01  2.5496E-02  1.7177E+00  9.0912E-02  7.4610E-01  1.4521E-02  4.7973E-01  5.8725E+00  8.2775E-01
             4.4605E+00
 PARAMETER: -5.1506E-02 -1.6347E+00 -3.5692E+00  6.4100E-01 -2.2979E+00 -1.9289E-01 -4.1321E+00 -6.3453E-01  1.8703E+00 -8.9040E-02
             1.5953E+00
 GRADIENT:  -1.0300E+02 -9.7609E+00  5.0744E+01 -3.3479E+01  8.1481E+01 -1.5581E+01  3.8122E-03  1.5944E+00 -8.6895E-02  3.3370E+01
            -7.0576E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1147.35879332055        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  9.1543E-01  1.8124E-01  2.0949E-02  2.1300E+00  9.0948E-02  7.6598E-01  1.0000E-02  2.3388E-01  5.8456E+00  6.9291E-01
             5.0061E+00
 PARAMETER:  1.1636E-02 -1.6079E+00 -3.7656E+00  8.5614E-01 -2.2975E+00 -1.6660E-01 -7.0176E+00 -1.3529E+00  1.8657E+00 -2.6686E-01
             1.7106E+00
 GRADIENT:  -1.4384E+01  3.5266E+00 -1.4892E+01 -9.7266E+00 -7.7811E+00 -7.9703E+00  0.0000E+00  1.8047E-01 -4.7978E+00  7.2519E+00
             7.2794E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1169.21891195996        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      544
 NPARAMETR:  7.7977E-01  2.2274E-01  4.5959E-02  7.9786E+00  1.0313E-01  7.0153E-01  1.0000E-02  8.9021E-02  6.2244E+00  1.7854E-01
             5.4194E+00
 PARAMETER: -1.4876E-01 -1.4017E+00 -2.9800E+00  2.1768E+00 -2.1718E+00 -2.5450E-01 -7.8748E+00 -2.3189E+00  1.9285E+00 -1.6230E+00
             1.7900E+00
 GRADIENT:   7.5638E+01  3.6749E+01  1.1959E+01 -1.4866E+00 -5.7843E+01  5.1503E-01  0.0000E+00  4.8014E-02  1.6742E+01 -2.1723E+00
             1.9312E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1175.07158024105        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      719
 NPARAMETR:  7.0013E-01  2.1158E-01  5.2899E-02  9.7455E+00  1.0538E-01  6.4052E-01  1.0000E-02  8.9750E-02  5.6176E+00  1.3517E-01
             5.2203E+00
 PARAMETER: -2.5649E-01 -1.4531E+00 -2.8394E+00  2.3768E+00 -2.1502E+00 -3.4547E-01 -7.4824E+00 -2.3107E+00  1.8259E+00 -1.9013E+00
             1.7525E+00
 GRADIENT:   2.9187E+00 -1.6429E+00  1.4290E+00 -8.4663E-01 -1.9727E+00  1.1478E+00  0.0000E+00  2.1655E-02  9.4374E-01 -1.0460E+00
            -2.2227E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1175.62478738083        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  6.9322E-01  2.1572E-01  5.3348E-02  9.8076E+00  1.0578E-01  6.1593E-01  1.0000E-02  4.7085E-02  5.5219E+00  2.6956E-01
             5.2102E+00
 PARAMETER: -2.6641E-01 -1.4338E+00 -2.8309E+00  2.3832E+00 -2.1464E+00 -3.8463E-01 -7.5837E+00 -2.9558E+00  1.8087E+00 -1.2110E+00
             1.7506E+00
 GRADIENT:   8.9294E+00 -1.8416E+00  2.7450E+01  6.9787E+00  1.2255E+02 -4.6828E-02  0.0000E+00  7.5055E-03  1.3183E+01  1.4540E+00
             2.5246E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1175.65762360111        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  6.9445E-01  2.1293E-01  5.2684E-02  9.6303E+00  1.0469E-01  6.2482E-01  1.0000E-02  1.9822E-02  5.4858E+00  3.0119E-01
             5.1618E+00
 PARAMETER: -2.6464E-01 -1.4468E+00 -2.8435E+00  2.3649E+00 -2.1568E+00 -3.7030E-01 -7.5837E+00 -3.8210E+00  1.8022E+00 -1.1000E+00
             1.7413E+00
 GRADIENT:   4.1696E+00 -5.9594E-01  2.7864E+01  5.3324E+00  1.1127E+02 -4.4915E-01  0.0000E+00  1.6802E-03  1.1361E+01  1.9831E+00
             1.9484E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1175.65819650437        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1007
 NPARAMETR:  6.9554E-01  2.1203E-01  5.2466E-02  9.6130E+00  1.0441E-01  6.2752E-01  1.0000E-02  1.6387E-02  5.4885E+00  3.0492E-01
             5.1565E+00
 PARAMETER: -2.6307E-01 -1.4510E+00 -2.8476E+00  2.3631E+00 -2.1594E+00 -3.6598E-01 -7.5837E+00 -4.0113E+00  1.8027E+00 -1.0877E+00
             1.7403E+00
 GRADIENT:  -5.6370E+00 -4.6101E+00  2.8283E+00 -2.5347E+00  6.8917E+00 -8.9805E-01  0.0000E+00  5.7866E-04 -4.9728E-01  1.4971E+00
             5.0255E+00

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1175.81182356778        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  6.9715E-01  2.1285E-01  5.1926E-02  9.9988E+00  1.0403E-01  6.5366E-01  1.0000E-02  1.0000E-02  5.5341E+00  2.9237E-01
             5.1160E+00
 PARAMETER: -2.6075E-01 -1.4472E+00 -2.8579E+00  2.4025E+00 -2.1631E+00 -3.2518E-01 -7.5837E+00 -4.6250E+00  1.8109E+00 -1.1298E+00
             1.7324E+00
 GRADIENT:  -4.3724E+00  1.0083E+00 -9.8721E-01 -5.5212E-01 -3.3075E+00  5.9900E-02  0.0000E+00  0.0000E+00 -1.1828E+00  5.8970E-02
            -2.2691E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1101
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.5818E-03 -1.8418E-05  2.8410E-04 -2.6089E-02  2.9349E-02
 SE:             1.7856E-02  1.0743E-04  1.7976E-04  2.2778E-02  1.3068E-02
 N:                     100         100         100         100         100

 P VAL.:         7.1242E-01  8.6388E-01  1.1400E-01  2.5207E-01  2.4713E-02

 ETASHRINKSD(%)  4.0181E+01  9.9640E+01  9.9398E+01  2.3689E+01  5.6220E+01
 ETASHRINKVR(%)  6.4217E+01  9.9999E+01  9.9996E+01  4.1767E+01  8.0833E+01
 EBVSHRINKSD(%)  4.2254E+01  9.9534E+01  9.9382E+01  1.5830E+01  5.8601E+01
 EBVSHRINKVR(%)  6.6655E+01  9.9998E+01  9.9996E+01  2.9153E+01  8.2861E+01
 RELATIVEINF(%)  7.0485E+00  8.9280E-04  8.7619E-04  2.6680E+01  4.9533E+00
 EPSSHRINKSD(%)  1.3937E+01
 EPSSHRINKVR(%)  2.5932E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1175.8118235677846     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -440.66099700404641     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1175.812       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.97E-01  2.13E-01  5.19E-02  1.00E+01  1.04E-01  6.54E-01  1.00E-02  1.00E-02  5.53E+00  2.92E-01  5.12E+00
 


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
+        1.76E+03
 
 TH 2
+        4.01E+02  5.56E+03
 
 TH 3
+        1.16E+04  5.83E+03  1.58E+05
 
 TH 4
+        8.38E+00 -2.06E+00 -3.71E+01  2.68E-01
 
 TH 5
+       -2.33E+03 -2.37E+04 -7.09E+04  6.89E+01  1.77E+05
 
 TH 6
+       -3.09E+01 -1.51E+01  1.07E+02 -1.41E-01  1.67E+02  5.32E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -9.50E+00  3.89E+00  6.52E+01 -2.63E-01  5.83E+01  1.81E+00  0.00E+00  0.00E+00  2.53E+00
 
 TH10
+        2.09E+01 -1.93E+02  4.12E+02  6.04E-04  1.51E+03  1.57E+00  0.00E+00  0.00E+00  2.11E+00  8.36E+01
 
 TH11
+       -2.98E+01 -3.39E+01 -2.81E+02 -5.73E-02  3.79E+02  9.89E+00  0.00E+00  0.00E+00  8.44E-01  1.47E+01  1.88E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.964
Stop Time:
Wed Sep 29 13:40:29 CDT 2021
