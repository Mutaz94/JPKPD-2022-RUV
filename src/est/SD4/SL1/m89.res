Sun Oct 24 03:07:49 CDT 2021
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
$DATA ../../../../data/SD4/SL1/dat89.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1670.20767839922        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0484E+02 -2.0375E+01 -1.7783E+01  2.2295E+01  5.2952E+01  6.4343E+01  1.7244E+00  1.3772E+00  2.4367E+01 -1.2700E+01
             1.1584E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1675.16995303707        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9691E-01  1.0447E+00  1.0020E+00  1.0059E+00  9.8022E-01  8.6603E-01  9.9979E-01  1.0010E+00  9.1878E-01  1.0493E+00
             9.6537E-01
 PARAMETER:  9.6905E-02  1.4372E-01  1.0197E-01  1.0592E-01  8.0017E-02 -4.3833E-02  9.9793E-02  1.0100E-01  1.5297E-02  1.4810E-01
             6.4760E-02
 GRADIENT:   4.0459E+00  2.8016E+00  2.6319E+00  1.0368E+00 -6.6399E+00 -1.4506E+01 -1.7014E+00 -2.1455E-01  2.5954E+00 -2.6014E+00
            -5.5417E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1675.75144740836        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.9325E-01  9.0562E-01  1.1537E+00  1.1070E+00  9.8575E-01  8.8496E-01  1.1723E+00  1.0872E+00  8.1473E-01  1.1076E+00
             9.7380E-01
 PARAMETER:  9.3224E-02  8.6543E-04  2.4294E-01  2.0165E-01  8.5649E-02 -2.2213E-02  2.5895E-01  1.8362E-01 -1.0490E-01  2.0221E-01
             7.3454E-02
 GRADIENT:  -2.1489E+00  1.3129E+01  4.2943E+00  1.1331E+01 -1.0997E+01 -5.1426E+00 -8.4929E-01 -9.4719E-01 -3.3197E+00  3.6220E+00
            -1.8424E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1676.30551733060        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.9261E-01  7.6132E-01  1.2774E+00  1.2023E+00  9.8463E-01  9.0563E-01  1.3348E+00  1.1619E+00  7.7777E-01  1.0775E+00
             9.7672E-01
 PARAMETER:  9.2587E-02 -1.7270E-01  3.4479E-01  2.8425E-01  8.4513E-02  8.7357E-04  3.8878E-01  2.5010E-01 -1.5133E-01  1.7462E-01
             7.6446E-02
 GRADIENT:   7.3062E-01  9.5957E+00  2.7899E+00  1.6559E+01 -9.1727E-02  4.3332E+00 -1.3544E+00 -1.7157E+00  6.8786E-01 -2.9669E+00
            -7.4166E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1678.44980542277        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.8795E-01  3.8971E-01  1.1485E+00  1.4068E+00  8.1397E-01  8.8297E-01  2.3600E+00  9.6572E-01  6.5940E-01  9.7827E-01
             9.7604E-01
 PARAMETER:  8.7881E-02 -8.4235E-01  2.3843E-01  4.4135E-01 -1.0584E-01 -2.4469E-02  9.5867E-01  6.5117E-02 -3.1643E-01  7.8035E-02
             7.5746E-02
 GRADIENT:  -1.6040E+00  4.5795E+00 -8.1283E-01  1.3360E+01 -2.4365E+00 -5.0431E+00 -4.9481E-01  1.6028E-01 -1.1626E+00  3.0984E+00
             1.2749E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1679.11635324771        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.8769E-01  2.6383E-01  9.2372E-01  1.4452E+00  6.8273E-01  8.9907E-01  3.0293E+00  7.4303E-01  6.5504E-01  8.2978E-01
             9.7644E-01
 PARAMETER:  8.7615E-02 -1.2325E+00  2.0657E-02  4.6823E-01 -2.8165E-01 -6.3932E-03  1.2083E+00 -1.9702E-01 -3.2306E-01 -8.6590E-02
             7.6157E-02
 GRADIENT:   5.3757E-01  3.3350E-02 -5.8889E+00  9.8590E+00  9.9825E+00  1.9687E+00 -2.8732E-01 -5.9473E-02  4.3360E-01 -1.5930E+00
            -1.1848E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1679.63656109657        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  9.9031E-01  3.0118E-01  6.9469E-01  1.3798E+00  5.6049E-01  8.9720E-01  2.8617E+00  4.8035E-01  6.6889E-01  7.0657E-01
             9.8129E-01
 PARAMETER:  9.0260E-02 -1.1000E+00 -2.6428E-01  4.2197E-01 -4.7895E-01 -8.4816E-03  1.1514E+00 -6.3324E-01 -3.0214E-01 -2.4734E-01
             8.1113E-02
 GRADIENT:  -1.3277E+00  3.0447E+00  2.2089E+01  1.0713E-01 -2.9798E+01  2.3198E-01  1.4163E+00 -2.6566E-01 -6.1951E-01 -1.5581E+00
             2.0102E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1680.31632416786        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  9.9498E-01  3.9527E-01  5.9244E-01  1.3048E+00  5.2647E-01  8.9909E-01  2.3818E+00  3.4225E-01  6.8771E-01  6.3828E-01
             9.8061E-01
 PARAMETER:  9.4968E-02 -8.2817E-01 -4.2350E-01  3.6607E-01 -5.4156E-01 -6.3672E-03  9.6784E-01 -9.7220E-01 -2.7439E-01 -3.4897E-01
             8.0424E-02
 GRADIENT:   7.8423E-01  8.3786E-01  3.4473E+00 -2.0944E-01 -4.8564E+00  1.2068E-01 -3.6006E-01  5.5536E-01  3.2811E-01 -1.2192E+00
            -3.7920E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1680.44089601632        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1423             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9583E-01  4.0039E-01  5.6462E-01  1.2926E+00  5.1197E-01  8.9929E-01  2.3597E+00  2.4287E-01  6.8815E-01  6.3771E-01
             9.8278E-01
 PARAMETER:  9.5821E-02 -8.1531E-01 -4.7161E-01  3.5668E-01 -5.6948E-01 -6.1452E-03  9.5855E-01 -1.3152E+00 -2.7375E-01 -3.4987E-01
             8.2629E-02
 GRADIENT:   3.7144E+02  5.9596E+01  3.5362E+01  3.8077E+02  9.9462E+01  2.4624E+01  5.9042E+01  5.0203E-01  1.3677E+01  2.3080E+00
             7.6135E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1680.44089601632        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:     1481
 NPARAMETR:  9.9583E-01  4.0039E-01  5.6462E-01  1.2926E+00  5.1197E-01  8.9929E-01  2.3597E+00  2.4287E-01  6.8815E-01  6.3771E-01
             9.8278E-01
 PARAMETER:  9.5821E-02 -8.1531E-01 -4.7161E-01  3.5668E-01 -5.6948E-01 -6.1452E-03  9.5855E-01 -1.3152E+00 -2.7375E-01 -3.4987E-01
             8.2629E-02
 GRADIENT:   1.2318E+00 -9.0946E-02  4.3161E-01 -3.4194E+00 -3.6318E-01  4.5750E-02  5.9060E-01  6.4420E-03  4.4037E-02 -2.3867E-02
            -3.8685E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1481
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0544E-03  3.9939E-02 -1.5276E-02 -2.7556E-02  1.9315E-02
 SE:             2.9843E-02  2.0534E-02  6.7790E-03  2.5599E-02  2.1698E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7181E-01  5.1779E-02  2.4230E-02  2.8173E-01  3.7336E-01

 ETASHRINKSD(%)  2.3861E-02  3.1207E+01  7.7289E+01  1.4241E+01  2.7310E+01
 ETASHRINKVR(%)  4.7716E-02  5.2676E+01  9.4842E+01  2.6453E+01  4.7161E+01
 EBVSHRINKSD(%)  5.1048E-01  3.3590E+01  7.7939E+01  1.2830E+01  2.4652E+01
 EBVSHRINKVR(%)  1.0184E+00  5.5897E+01  9.5133E+01  2.4013E+01  4.3227E+01
 RELATIVEINF(%)  9.8339E+01  9.8773E+00  2.9933E-01  2.3067E+01  3.1029E+00
 EPSSHRINKSD(%)  4.4659E+01
 EPSSHRINKVR(%)  6.9373E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1680.4408960163169     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -945.29006945257868     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1680.441       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  4.00E-01  5.65E-01  1.29E+00  5.12E-01  8.99E-01  2.36E+00  2.43E-01  6.88E-01  6.38E-01  9.83E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       45.593
Stop Time:
Sun Oct 24 03:07:59 CDT 2021
