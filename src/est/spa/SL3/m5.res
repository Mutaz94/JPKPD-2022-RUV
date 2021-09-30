Wed Sep 29 16:19:41 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat5.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1620.70682876488        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0037E+02 -2.0634E+01 -2.9824E+01  2.8959E+01  4.7690E+01  3.6959E+01  4.3972E-02  1.0114E+01  2.2924E+01 -2.1732E+01
            -1.3704E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1631.14206437390        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.4244E-01  1.0445E+00  1.1971E+00  1.0252E+00  1.0924E+00  9.7182E-01  1.0440E+00  7.8926E-01  7.9652E-01  1.3221E+00
             1.0113E+00
 PARAMETER:  4.0715E-02  1.4351E-01  2.7992E-01  1.2492E-01  1.8840E-01  7.1418E-02  1.4306E-01 -1.3665E-01 -1.2751E-01  3.7922E-01
             1.1123E-01
 GRADIENT:   3.6356E+02  6.2934E+01 -3.3086E+00  1.0615E+02  5.1011E-01  3.6955E+01 -8.4875E+00  3.1780E+00 -1.2349E+01  1.1522E+01
            -1.2276E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1633.15682424703        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  9.2885E-01  9.4429E-01  1.4664E+00  1.0960E+00  1.1600E+00  9.5509E-01  1.1630E+00  4.9354E-01  8.9562E-01  1.3723E+00
             1.0722E+00
 PARAMETER:  2.6192E-02  4.2681E-02  4.8281E-01  1.9163E-01  2.4843E-01  5.4054E-02  2.5102E-01 -6.0616E-01 -1.0242E-02  4.1651E-01
             1.6974E-01
 GRADIENT:  -4.0856E+01  1.7383E+01  2.2901E+00  2.4068E+01 -4.0848E+00 -3.5246E-01  5.2693E+00  2.8604E-01  7.2683E+00 -1.9214E+00
             7.3040E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1634.19244663401        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      402
 NPARAMETR:  9.4426E-01  9.1441E-01  1.4469E+00  1.1003E+00  1.1513E+00  9.5288E-01  1.0841E+00  3.0129E-01  8.8224E-01  1.3878E+00
             1.0556E+00
 PARAMETER:  4.2651E-02  1.0522E-02  4.6943E-01  1.9560E-01  2.4088E-01  5.1729E-02  1.8073E-01 -1.0997E+00 -2.5296E-02  4.2772E-01
             1.5416E-01
 GRADIENT:  -4.3241E-01  2.9089E+00 -1.6147E-01  3.6082E+00  3.9702E-01 -4.8901E-01  3.3690E-01  1.0530E-01  8.7370E-01 -3.0258E-01
            -2.5405E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1634.36485900976        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  9.4253E-01  6.8742E-01  1.6492E+00  1.2484E+00  1.1267E+00  9.5243E-01  1.0628E+00  1.0554E-02  8.4453E-01  1.4192E+00
             1.0631E+00
 PARAMETER:  4.0814E-02 -2.7481E-01  6.0026E-01  3.2186E-01  2.1929E-01  5.1264E-02  1.6090E-01 -4.4512E+00 -6.8975E-02  4.5008E-01
             1.6123E-01
 GRADIENT:   4.6338E-02  2.7810E+00  1.5313E+00  3.7064E+00 -3.1766E+00  4.2214E-02 -2.7604E-01  3.5908E-05 -7.2519E-02 -2.8292E-01
            -6.2821E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1634.38723247196        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      756
 NPARAMETR:  9.4173E-01  6.0515E-01  1.7071E+00  1.3014E+00  1.1183E+00  9.5174E-01  1.0757E+00  1.0000E-02  8.2574E-01  1.4272E+00
             1.0639E+00
 PARAMETER:  3.9964E-02 -4.0228E-01  6.3481E-01  3.6340E-01  2.1179E-01  5.0536E-02  1.7298E-01 -5.9586E+00 -9.1476E-02  4.5569E-01
             1.6194E-01
 GRADIENT:  -1.2383E-01  1.9011E+00  8.1276E-01  2.8223E+00 -1.6514E+00  9.6169E-03 -2.0546E-01  0.0000E+00 -1.1295E-01 -2.6441E-01
            -7.0759E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1634.38752987218        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      936
 NPARAMETR:  9.4157E-01  5.8747E-01  1.7198E+00  1.3130E+00  1.1163E+00  9.5160E-01  1.0787E+00  1.0000E-02  8.2168E-01  1.4289E+00
             1.0641E+00
 PARAMETER:  3.9793E-02 -4.3193E-01  6.4218E-01  3.7229E-01  2.0998E-01  5.0388E-02  1.7576E-01 -6.3132E+00 -9.6401E-02  4.5687E-01
             1.6213E-01
 GRADIENT:  -1.2893E-01  1.9171E+00  7.6112E-01  3.0896E+00 -1.5762E+00  6.1482E-03 -1.8848E-01  0.0000E+00 -1.0483E-01 -2.4105E-01
            -8.0326E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1634.39925338715        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1120
 NPARAMETR:  9.4136E-01  5.8071E-01  1.7168E+00  1.3157E+00  1.1163E+00  9.5153E-01  1.1193E+00  1.0000E-02  8.1776E-01  1.4306E+00
             1.0634E+00
 PARAMETER:  3.9575E-02 -4.4350E-01  6.4046E-01  3.7439E-01  2.1002E-01  5.0320E-02  2.1274E-01 -6.3651E+00 -1.0118E-01  4.5810E-01
             1.6151E-01
 GRADIENT:  -4.4843E-01  4.9444E-01 -1.2101E-01 -4.1207E-01  4.6071E-01 -6.3873E-03 -7.4607E-02  0.0000E+00  3.8654E-01  6.7135E-02
             3.9779E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1634.40344722557        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1295
 NPARAMETR:  9.4110E-01  5.6765E-01  1.7122E+00  1.3256E+00  1.1092E+00  9.5136E-01  1.2285E+00  1.0000E-02  8.0202E-01  1.4277E+00
             1.0617E+00
 PARAMETER:  3.9289E-02 -4.6625E-01  6.3778E-01  3.8186E-01  2.0361E-01  5.0134E-02  3.0582E-01 -6.3651E+00 -1.2062E-01  4.5607E-01
             1.5986E-01
 GRADIENT:  -8.3931E-01  1.5014E+00  5.2114E-01  1.9151E+00 -8.2503E-01 -5.2062E-02 -6.0461E-02  0.0000E+00 -3.0510E-01  4.5514E-01
            -5.6827E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1634.40980394332        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1477             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4159E-01  5.5706E-01  1.7065E+00  1.3292E+00  1.1062E+00  9.5143E-01  1.2739E+00  1.0000E-02  7.9760E-01  1.4239E+00
             1.0623E+00
 PARAMETER:  3.9816E-02 -4.8507E-01  6.3445E-01  3.8457E-01  2.0094E-01  5.0215E-02  3.4211E-01 -6.3651E+00 -1.2614E-01  4.5338E-01
             1.6042E-01
 GRADIENT:   3.3862E+02  4.3040E+01  3.9268E+00  3.5025E+02  1.1208E+01  2.6860E+01  3.7937E+00  0.0000E+00  7.4951E+00  7.0048E+00
             1.3186E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1634.41124740202        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1656
 NPARAMETR:  9.4154E-01  5.5758E-01  1.7047E+00  1.3304E+00  1.1046E+00  9.5144E-01  1.2646E+00  1.0000E-02  7.9764E-01  1.4216E+00
             1.0624E+00
 PARAMETER:  3.9764E-02 -4.8414E-01  6.3339E-01  3.8547E-01  1.9951E-01  5.0223E-02  3.3479E-01 -6.3651E+00 -1.2609E-01  4.5177E-01
             1.6055E-01
 GRADIENT:   5.4000E-01  4.7043E-01  1.2408E-01 -1.1249E+00  2.4090E-01  1.9707E-02  2.2894E-03  0.0000E+00 -6.0090E-02  8.0795E-02
            -4.6270E-02

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1634.41124740202        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1678
 NPARAMETR:  9.4154E-01  5.5758E-01  1.7047E+00  1.3304E+00  1.1046E+00  9.5144E-01  1.2646E+00  1.0000E-02  7.9764E-01  1.4216E+00
             1.0624E+00
 PARAMETER:  3.9764E-02 -4.8414E-01  6.3339E-01  3.8547E-01  1.9951E-01  5.0223E-02  3.3479E-01 -6.3651E+00 -1.2609E-01  4.5177E-01
             1.6055E-01
 GRADIENT:   5.4000E-01  4.7043E-01  1.2408E-01 -1.1249E+00  2.4090E-01  1.9707E-02  2.2894E-03  0.0000E+00 -6.0090E-02  8.0795E-02
            -4.6270E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1678
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.8064E-05 -9.9458E-03 -1.3516E-04 -6.0516E-03 -3.3237E-02
 SE:             2.9759E-02  1.1734E-02  9.2172E-05  2.6771E-02  2.4477E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9871E-01  3.9665E-01  1.4256E-01  8.2116E-01  1.7451E-01

 ETASHRINKSD(%)  3.0356E-01  6.0690E+01  9.9691E+01  1.0314E+01  1.7998E+01
 ETASHRINKVR(%)  6.0620E-01  8.4547E+01  9.9999E+01  1.9564E+01  3.2757E+01
 EBVSHRINKSD(%)  5.0907E-01  6.0559E+01  9.9687E+01  1.0328E+01  1.3733E+01
 EBVSHRINKVR(%)  1.0155E+00  8.4444E+01  9.9999E+01  1.9589E+01  2.5579E+01
 RELATIVEINF(%)  9.6221E+01  1.7498E-01  1.3764E-04  9.6992E-01  8.8013E+00
 EPSSHRINKSD(%)  4.0754E+01
 EPSSHRINKVR(%)  6.4899E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1634.4112474020217     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.26042083828349     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.47
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1634.411       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.42E-01  5.58E-01  1.70E+00  1.33E+00  1.10E+00  9.51E-01  1.26E+00  1.00E-02  7.98E-01  1.42E+00  1.06E+00
 


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
+       -1.90E+01  3.53E+02
 
 TH 3
+        1.69E+00  2.28E+01  3.88E+01
 
 TH 4
+       -1.07E+01  4.91E+02 -2.53E+01  7.95E+02
 
 TH 5
+       -3.90E-01 -1.12E+02 -9.81E+01  2.03E+00  3.43E+02
 
 TH 6
+        1.19E+00 -3.03E+00  7.03E-01 -3.35E+00 -1.36E+00  2.16E+02
 
 TH 7
+        4.35E-01 -2.51E+00  1.91E+00 -2.88E+00 -1.89E+00 -1.04E-01  3.17E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.74E+00 -2.20E+01  1.99E+00  9.60E+00 -9.59E-01 -1.04E+00  2.11E+01  0.00E+00  2.04E+02
 
 TH10
+        6.84E-01  3.97E+00 -4.42E+00 -4.35E+00 -3.64E+01 -3.47E-02  1.19E+00  0.00E+00 -2.00E-01  5.33E+01
 
 TH11
+       -1.15E+01 -1.71E+01 -1.52E+01 -1.36E+01  9.76E+00  2.61E+00  1.83E+00  0.00E+00  1.06E+01  1.70E+01  2.13E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.550
Stop Time:
Wed Sep 29 16:20:10 CDT 2021
