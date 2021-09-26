Sat Sep 25 02:20:56 CDT 2021
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
$DATA ../../../../data/int/SL3/dat49.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      976
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

 TOT. NO. OF OBS RECS:      876
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   94.7653744144389        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4537E+02  6.5627E+01  1.3886E+02 -5.2927E+01  1.4618E+02 -1.5206E+01 -1.3266E+02 -5.1942E+02 -1.7665E+02 -6.2775E+01
            -6.7399E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2612.94566042321        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.5352E-01  1.1717E+00  1.3955E+00  9.2737E-01  1.1782E+00  1.0486E+00  1.4095E+00  9.7027E-01  1.1325E+00  1.1444E+00
             2.6615E+00
 PARAMETER:  5.2401E-02  2.5844E-01  4.3327E-01  2.4597E-02  2.6396E-01  1.4749E-01  4.4321E-01  6.9823E-02  2.2441E-01  2.3492E-01
             1.0789E+00
 GRADIENT:  -2.6430E+01 -1.0371E+01 -7.0141E-01 -1.8264E+01 -1.2207E+01  7.0955E+00  2.1850E+01  1.3986E-01  1.1479E+01  3.7615E+00
             1.0205E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2616.66998782326        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.7598E-01  1.3356E+00  2.1979E+00  8.6569E-01  1.4974E+00  1.0437E+00  1.0346E+00  9.6757E-01  1.2303E+00  1.3685E+00
             2.6766E+00
 PARAMETER:  7.5689E-02  3.8937E-01  8.8752E-01 -4.4226E-02  5.0375E-01  1.4274E-01  1.3406E-01  6.7037E-02  3.0723E-01  4.1373E-01
             1.0846E+00
 GRADIENT:   1.8411E+01 -4.6475E+00 -7.6127E-01  7.8818E+00  1.5117E+01  6.2300E+00 -3.1361E+00 -2.7418E+00  6.5883E+00 -6.6331E+00
             1.8067E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2621.23043551127        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.6972E-01  1.3817E+00  2.6755E+00  8.4943E-01  1.5985E+00  1.0341E+00  1.0634E+00  2.6247E+00  1.1558E+00  1.3871E+00
             2.6229E+00
 PARAMETER:  6.9248E-02  4.2331E-01  1.0841E+00 -6.3184E-02  5.6904E-01  1.3350E-01  1.6145E-01  1.0650E+00  2.4477E-01  4.2719E-01
             1.0643E+00
 GRADIENT:   6.7312E+00  8.7522E-01 -1.0841E+01  1.0915E+01  2.5655E+01  2.0280E+00 -2.0574E+00  1.4232E+00  2.4138E+00 -2.3037E-01
            -9.6662E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2621.40603063796        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      369             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6981E-01  1.3762E+00  2.7077E+00  8.3494E-01  1.5978E+00  1.0316E+00  1.1102E+00  2.6427E+00  1.1151E+00  1.3874E+00
             2.6237E+00
 PARAMETER:  6.9349E-02  4.1931E-01  1.0961E+00 -8.0399E-02  5.6865E-01  1.3108E-01  2.0456E-01  1.0718E+00  2.0892E-01  4.2740E-01
             1.0646E+00
 GRADIENT:   7.3516E+00 -1.3633E+01 -8.7967E+00 -7.7058E+00  2.1723E+01  1.2016E+00  7.9468E-01  7.4015E-01  1.0749E+00  4.3339E-01
             2.3295E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2621.45806834089        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      485
 NPARAMETR:  9.6860E-01  1.3762E+00  2.7077E+00  8.4473E-01  1.5978E+00  1.0312E+00  1.1115E+00  2.6427E+00  1.1025E+00  1.3874E+00
             2.6237E+00
 PARAMETER:  6.8095E-02  4.1931E-01  1.0961E+00 -6.8742E-02  5.6865E-01  1.3077E-01  2.0574E-01  1.0718E+00  1.9759E-01  4.2740E-01
             1.0646E+00
 GRADIENT:  -2.1182E+00 -1.1446E+01 -1.0018E+01  4.8304E-01  2.0582E+01 -4.3125E-01 -1.8144E-01  5.5263E-01 -1.9800E-02  2.9074E-01
            -2.2780E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2622.06635487251        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      619
 NPARAMETR:  9.6966E-01  1.3722E+00  2.9864E+00  8.4433E-01  1.5760E+00  1.0283E+00  1.1089E+00  2.6413E+00  1.1016E+00  1.3830E+00
             2.6229E+00
 PARAMETER:  6.9187E-02  4.1638E-01  1.1941E+00 -6.9210E-02  5.5486E-01  1.2792E-01  2.0334E-01  1.0713E+00  1.9674E-01  4.2426E-01
             1.0643E+00
 GRADIENT:   7.0805E+00  2.9962E-01 -1.1118E-01 -1.0254E+01  9.1610E-02  2.6299E-02 -1.4961E-01 -2.7595E+00 -7.4774E-01 -1.5595E-01
             2.0828E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2622.11157071457        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      797
 NPARAMETR:  9.6966E-01  1.3827E+00  3.0408E+00  8.4433E-01  1.5902E+00  1.0323E+00  1.1114E+00  2.6413E+00  1.1016E+00  1.3935E+00
             2.6253E+00
 PARAMETER:  6.9187E-02  4.2404E-01  1.2121E+00 -6.9210E-02  5.6383E-01  1.3182E-01  2.0560E-01  1.0713E+00  1.9674E-01  4.3179E-01
             1.0652E+00
 GRADIENT:   1.7864E-01  9.7235E-03 -9.0281E-03 -5.7204E+00 -1.2756E-02  1.1675E-03 -1.7369E-04 -3.8377E+00 -9.4706E-01 -4.2859E-03
             7.4142E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2622.21041392203        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      919
 NPARAMETR:  9.6942E-01  1.3818E+00  3.0460E+00  8.4651E-01  1.5898E+00  1.0321E+00  1.1111E+00  2.7131E+00  1.1048E+00  1.3932E+00
             2.6251E+00
 PARAMETER:  6.8944E-02  4.2341E-01  1.2138E+00 -6.6635E-02  5.6361E-01  1.3157E-01  2.0534E-01  1.0981E+00  1.9965E-01  4.3158E-01
             1.0651E+00
 GRADIENT:  -3.2800E-01  6.3209E-01 -1.1267E+00 -3.6639E+00  1.1504E-01 -1.1062E-01  8.6729E-02 -2.3871E+00 -2.7377E-01  1.0141E-01
             6.4396E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2622.21088409794        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1105
 NPARAMETR:  9.6953E-01  1.3819E+00  3.0462E+00  8.4650E-01  1.5896E+00  1.0323E+00  1.1096E+00  2.7129E+00  1.1076E+00  1.3928E+00
             2.6253E+00
 PARAMETER:  6.9057E-02  4.2344E-01  1.2139E+00 -6.6642E-02  5.6350E-01  1.3177E-01  2.0398E-01  1.0980E+00  2.0215E-01  4.3131E-01
             1.0652E+00
 GRADIENT:  -9.9173E-02  6.1724E-01 -1.1113E+00 -3.6090E+00  6.5085E-02 -3.2113E-02  7.4259E-02 -2.3602E+00 -1.1137E-01  2.4112E-02
             8.5410E-01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -2622.21088409794        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     1132
 NPARAMETR:  9.6953E-01  1.3818E+00  3.0466E+00  8.4650E-01  1.5896E+00  1.0323E+00  1.1093E+00  2.7126E+00  1.1075E+00  1.3928E+00
             2.6255E+00
 PARAMETER:  6.9057E-02  4.2344E-01  1.2139E+00 -6.6642E-02  5.6350E-01  1.3177E-01  2.0398E-01  1.0980E+00  2.0215E-01  4.3131E-01
             1.0652E+00
 GRADIENT:   1.5975E+00  6.4930E+04 -2.2647E+04 -2.0096E+01  1.3032E+00  4.1731E+05  6.7482E-02  5.0035E+04  2.7202E+05 -4.2988E-01
            -5.1685E+04
 NUMSIGDIG:         9.1         3.3         3.3         8.0         8.5         3.3         2.3         3.3         3.3         9.1
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1132
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7443E-03 -2.0099E-02 -3.0078E-02  1.3601E-02 -3.0969E-02
 SE:             2.9386E-02  2.1688E-02  1.4731E-02  2.1383E-02  2.3587E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5267E-01  3.5407E-01  4.1167E-02  5.2472E-01  1.8918E-01

 ETASHRINKSD(%)  1.5533E+00  2.7344E+01  5.0649E+01  2.8365E+01  2.0982E+01
 ETASHRINKVR(%)  3.0824E+00  4.7211E+01  7.5645E+01  4.8684E+01  3.7561E+01
 EBVSHRINKSD(%)  1.5995E+00  2.7174E+01  5.5737E+01  3.1502E+01  1.8112E+01
 EBVSHRINKVR(%)  3.1735E+00  4.6963E+01  8.0408E+01  5.3080E+01  3.2944E+01
 RELATIVEINF(%)  9.6782E+01  5.1328E+00  6.0615E+00  4.5017E+00  3.0409E+01
 EPSSHRINKSD(%)  1.7285E+01
 EPSSHRINKVR(%)  3.1582E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          876
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1609.9803101745865     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2622.2108840979363     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1012.2305739233498     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.87
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2622.211       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.70E-01  1.38E+00  3.05E+00  8.47E-01  1.59E+00  1.03E+00  1.11E+00  2.71E+00  1.11E+00  1.39E+00  2.63E+00
 


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
+        2.92E+09
 
 TH 2
+        7.41E+02  2.60E+02
 
 TH 3
+       -1.18E+02  5.92E+00  1.01E+06
 
 TH 4
+       -1.48E+01 -8.00E+03  1.30E+03  3.84E+09
 
 TH 5
+       -1.58E+08 -2.62E+07 -6.57E+01 -1.81E+08  1.52E+02
 
 TH 6
+        4.85E+00  1.17E+03 -1.85E+02 -1.13E+00 -1.28E+00  7.43E+08
 
 TH 7
+        8.77E-01  2.54E+03 -4.00E+02 -3.11E+00 -1.33E+00 -4.47E+08  5.25E+01
 
 TH 8
+        1.47E+02 -5.75E+02  4.48E+02 -1.63E+03  5.16E+06  2.30E+02  4.96E+02  1.55E+06
 
 TH 9
+        2.13E+00  4.61E+02  1.66E+07  1.75E+01  3.75E+00  3.06E+03  6.63E+03  9.43E+01  5.48E+08
 
 TH10
+        4.63E-01  3.91E+07 -6.19E+06  8.18E+00 -2.55E+07  4.32E-01  4.92E+00 -7.69E+06  1.02E+08  7.62E+07
 
 TH11
+       -1.72E+02  5.96E+02 -4.81E+02  1.73E+03 -7.19E+01 -2.43E+02 -5.28E+02  1.65E+06 -9.24E+01  3.44E+01  1.76E+06
 
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
 #CPUT: Total CPU Time in Seconds,       44.379
Stop Time:
Sat Sep 25 02:21:41 CDT 2021
