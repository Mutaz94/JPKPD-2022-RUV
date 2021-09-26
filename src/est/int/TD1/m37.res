Sat Sep 25 03:59:03 CDT 2021
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
$DATA ../../../../data/int/TD1/dat37.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3402.17740759419        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6669E+01  2.9265E+01  1.1750E+02 -4.5868E+01  1.4777E+01  1.0467E+01 -1.0785E+02 -5.1215E+02 -1.4531E+02 -4.9717E+01
            -1.3264E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3646.44350755654        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0258E+00  9.6487E-01  9.8464E-01  1.0051E+00  9.9597E-01  9.8299E-01  1.4005E+00  2.6810E+00  1.0875E+00  1.1052E+00
             9.1414E-01
 PARAMETER:  1.2545E-01  6.4238E-02  8.4518E-02  1.0510E-01  9.5964E-02  8.2841E-02  4.3685E-01  1.0862E+00  1.8384E-01  2.0007E-01
             1.0227E-02
 GRADIENT:   1.1773E+02  7.8104E+00 -1.1388E+01 -1.3150E+01 -2.8936E+01  2.0094E+00 -2.1106E+01 -1.5973E+01  2.4247E+01  8.4692E-01
            -2.4007E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3658.30940371891        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0274E+00  1.1819E+00  1.1108E+00  9.4430E-01  1.1578E+00  1.0818E+00  9.6527E-01  3.2927E+00  1.0497E+00  1.2163E+00
             9.4751E-01
 PARAMETER:  1.2707E-01  2.6715E-01  2.0505E-01  4.2691E-02  2.4654E-01  1.7859E-01  6.4649E-02  1.2917E+00  1.4846E-01  2.9579E-01
             4.6087E-02
 GRADIENT:   1.1100E+02  6.4167E+01  5.5453E+00  2.1925E+01 -5.8492E+01  4.2783E+01 -3.9879E+01 -2.0345E+01  1.8671E+00 -1.3969E+01
            -1.8952E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3691.32923005679        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      240
 NPARAMETR:  1.0231E+00  1.1782E+00  1.2742E+00  9.4487E-01  1.1868E+00  9.2995E-01  9.0360E-01  3.4289E+00  1.1507E+00  1.1385E+00
             9.6046E-01
 PARAMETER:  1.2288E-01  2.6395E-01  3.4231E-01  4.3296E-02  2.7124E-01  2.7377E-02 -1.3658E-03  1.3322E+00  2.4040E-01  2.2975E-01
             5.9655E-02
 GRADIENT:   1.1143E+02  6.5357E+01  1.3442E+01  2.7527E+01 -4.0638E+01 -2.0511E+01 -2.9184E+01  3.3709E+01  2.1946E+01 -2.7561E+01
            -1.6988E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3691.56304827920        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  1.0229E+00  1.1778E+00  1.2739E+00  9.4494E-01  1.1868E+00  9.3035E-01  9.0638E-01  3.4256E+00  1.1490E+00  1.1382E+00
             9.6080E-01
 PARAMETER:  1.2266E-01  2.6362E-01  3.4206E-01  4.3367E-02  2.7128E-01  2.7807E-02  1.7072E-03  1.3313E+00  2.3891E-01  2.2943E-01
             6.0010E-02
 GRADIENT:   4.8608E+01  3.4310E+01  1.1425E+01  1.8409E+01 -5.3403E+01 -2.6117E+01 -2.9968E+01  2.4845E+01  1.9602E+01 -2.8729E+01
            -1.6907E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3694.35494406763        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      567
 NPARAMETR:  1.0229E+00  1.1776E+00  1.2737E+00  9.4498E-01  1.1869E+00  9.8997E-01  9.0635E-01  3.4274E+00  1.1491E+00  1.3348E+00
             9.6084E-01
 PARAMETER:  1.2262E-01  2.6352E-01  3.4192E-01  4.3407E-02  2.7139E-01  8.9916E-02  1.6670E-03  1.3318E+00  2.3900E-01  3.8877E-01
             6.0050E-02
 GRADIENT:   4.3020E+01  1.8172E+01  1.4714E+01  2.4460E+01 -4.5162E+01  3.0347E-02 -1.9085E+01  2.4296E+01  1.8787E+01 -1.3376E-02
            -1.5719E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3697.78339593669        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      755
 NPARAMETR:  1.0184E+00  1.1676E+00  1.2589E+00  9.4418E-01  1.2141E+00  9.8281E-01  9.0694E-01  2.9598E+00  1.1449E+00  1.3215E+00
             9.7119E-01
 PARAMETER:  1.1821E-01  2.5492E-01  3.3021E-01  4.2565E-02  2.9403E-01  8.2662E-02  2.3259E-03  1.1851E+00  2.3533E-01  3.7874E-01
             7.0765E-02
 GRADIENT:   3.3419E+01 -5.4703E+00  1.8323E+01  2.2777E+01 -1.5446E+01 -2.4380E+00 -2.0475E+01 -1.6327E+01  1.5911E+01  1.2741E+00
            -1.3532E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3702.36041060796        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      895
 NPARAMETR:  1.0108E+00  1.1667E+00  1.2193E+00  9.3827E-01  1.2185E+00  9.8763E-01  9.5544E-01  3.0095E+00  1.1169E+00  1.3096E+00
             9.8997E-01
 PARAMETER:  1.1079E-01  2.5416E-01  2.9828E-01  3.6285E-02  2.9761E-01  8.7556E-02  5.4422E-02  1.2018E+00  2.1053E-01  3.6971E-01
             8.9915E-02
 GRADIENT:   6.8876E+01  1.4443E+01  1.1629E+01  2.4961E+01 -2.2336E-01  5.5387E+00 -1.6931E+01 -2.0812E+00  1.6270E+01  5.2476E+00
            -9.0487E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3703.14080527418        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      980
 NPARAMETR:  1.0095E+00  1.1666E+00  1.2181E+00  9.3706E-01  1.2197E+00  9.8682E-01  9.7499E-01  3.0110E+00  1.1057E+00  1.3015E+00
             9.9220E-01
 PARAMETER:  1.0949E-01  2.5410E-01  2.9729E-01  3.4995E-02  2.9861E-01  8.6734E-02  7.4669E-02  1.2023E+00  2.0048E-01  3.6355E-01
             9.2169E-02
 GRADIENT:   6.5080E+01  1.3801E+01  1.0875E+01  2.3560E+01 -3.3356E-01  5.2395E+00 -1.6304E+01 -2.0015E+00  1.5661E+01  5.1270E+00
            -8.5363E+01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -3703.15687161664        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  1.0095E+00  1.1666E+00  1.2181E+00  9.3703E-01  1.2197E+00  9.8693E-01  9.7498E-01  3.0109E+00  1.1054E+00  1.3016E+00
             9.9233E-01
 PARAMETER:  1.0949E-01  2.5410E-01  2.9729E-01  3.4956E-02  2.9865E-01  8.6739E-02  7.4669E-02  1.2023E+00  2.0023E-01  3.6355E-01
             9.2302E-02
 GRADIENT:  -1.7828E+05  7.6807E+04  2.0220E+01 -1.9520E+05  6.5348E+04 -3.8607E-01  1.9519E+05  1.6198E+04 -9.7483E+04 -5.3694E+04
            -9.7693E+04
 NUMSIGDIG:         3.3         3.3         7.1         3.3         3.3         1.8         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1150
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.3800E-03 -4.6428E-02 -4.3322E-02  2.2582E-02 -4.6605E-02
 SE:             2.9928E-02  2.2771E-02  2.6387E-02  2.5338E-02  2.4930E-02
 N:                     100         100         100         100         100

 P VAL.:         8.5734E-01  4.1454E-02  1.0063E-01  3.7280E-01  6.1565E-02

 ETASHRINKSD(%)  1.0000E-10  2.3716E+01  1.1599E+01  1.5115E+01  1.6480E+01
 ETASHRINKVR(%)  1.0000E-10  4.1807E+01  2.1853E+01  2.7946E+01  3.0244E+01
 EBVSHRINKSD(%)  2.5959E-01  2.8795E+01  1.5173E+01  1.2375E+01  1.6398E+01
 EBVSHRINKVR(%)  5.1850E-01  4.9299E+01  2.8044E+01  2.3219E+01  3.0107E+01
 RELATIVEINF(%)  9.9480E+01  2.5292E+01  6.7559E+01  4.5147E+01  4.2671E+01
 EPSSHRINKSD(%)  2.0586E+01
 EPSSHRINKVR(%)  3.6934E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3703.1568716166439     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2049.0675118482332     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.07
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3703.157       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.17E+00  1.22E+00  9.37E-01  1.22E+00  9.87E-01  9.75E-01  3.01E+00  1.11E+00  1.30E+00  9.92E-01
 


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
+        3.99E+08
 
 TH 2
+        4.82E+02  5.55E+07
 
 TH 3
+       -2.03E+04  7.57E+03  7.44E+07
 
 TH 4
+        1.83E+04 -6.52E+03 -1.44E+08  1.11E+09
 
 TH 5
+        1.01E+04 -3.93E+03  6.12E+03 -5.45E+03  3.68E+07
 
 TH 6
+       -2.07E+03  7.57E+02  2.75E-01 -2.41E+03  6.22E+02  2.16E+02
 
 TH 7
+        3.53E+04 -1.32E+04  1.38E+08 -2.08E+04 -1.15E+04  2.33E+03  5.13E+08
 
 TH 8
+        3.95E+01 -1.66E+04  6.12E+02 -5.59E+02 -3.15E+02  6.29E+01 -1.08E+03  3.71E+05
 
 TH 9
+        5.94E+03 -2.23E+03 -6.09E+07  9.15E+03 -6.05E+07 -1.03E+03 -6.70E+03 -1.78E+02  9.96E+07
 
 TH10
+       -1.45E+03  5.20E+02 -4.73E+03  4.29E+03  2.36E+03 -4.79E+02  8.27E+03  4.30E+01  1.39E+03  2.18E+07
 
 TH11
+        3.57E+03 -1.39E+03 -1.08E+01  2.04E+04 -1.11E+03 -2.28E+03 -4.03E+03 -1.06E+02  1.79E+03  8.58E+02  4.96E+08
 
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
 #CPUT: Total CPU Time in Seconds,       49.165
Stop Time:
Sat Sep 25 03:59:53 CDT 2021
