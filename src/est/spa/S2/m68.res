Sat Sep 25 12:25:52 CDT 2021
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
$DATA ../../../../data/spa/S2/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1735.07751076603        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.1130E+00  1.8975E+01 -1.8031E+01  5.9269E+01  1.7919E+01  3.6969E+00  2.4405E+01  6.5893E+00  4.7491E+01  1.1824E+00
             6.2916E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1748.33059780394        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9457E-01  9.9809E-01  1.0023E+00  1.0071E+00  9.9592E-01  9.8921E-01  8.5690E-01  9.7145E-01  8.0327E-01  1.0188E+00
             8.5350E-01
 PARAMETER:  9.4551E-02  9.8086E-02  1.0225E-01  1.0703E-01  9.5914E-02  8.9151E-02 -5.4435E-02  7.1036E-02 -1.1907E-01  1.1859E-01
            -5.8409E-02
 GRADIENT:   1.3598E+01  5.2721E+01 -1.2094E+01  9.0829E+01  1.4191E+01  3.5799E-01  7.0816E+00  2.9477E+00  3.0921E+00 -1.4358E+00
             3.6334E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1749.46302229516        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      268
 NPARAMETR:  1.0115E+00  1.0679E+00  9.7260E-01  9.5719E-01  1.0166E+00  1.0020E+00  7.0271E-01  8.6079E-01  8.9700E-01  1.0825E+00
             8.4475E-01
 PARAMETER:  1.1148E-01  1.6568E-01  7.2216E-02  5.6251E-02  1.1646E-01  1.0202E-01 -2.5281E-01 -4.9906E-02 -8.6951E-03  1.7927E-01
            -6.8713E-02
 GRADIENT:  -9.7036E+00  4.3383E+01 -4.2584E+00  6.7223E+01  3.8070E+00 -4.8518E-01  2.9560E+00 -8.3036E-01  8.8843E+00  4.1760E+00
            -5.0910E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1751.05442502701        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  1.0155E+00  1.0166E+00  9.5593E-01  9.5675E-01  9.7904E-01  1.0018E+00  6.9256E-01  8.5960E-01  8.4372E-01  1.0223E+00
             8.4590E-01
 PARAMETER:  1.1542E-01  1.1650E-01  5.4933E-02  5.5789E-02  7.8819E-02  1.0181E-01 -2.6736E-01 -5.1286E-02 -6.9935E-02  1.2209E-01
            -6.7355E-02
 GRADIENT:  -4.5198E-01  1.7796E+00  1.0593E+00 -1.9397E+00 -9.4998E-01 -3.6200E-01  2.4186E-01 -5.0771E-02 -2.9152E-01 -4.8746E-01
            -4.8476E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1751.51293180096        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  1.0114E+00  7.5158E-01  1.2104E+00  1.1388E+00  9.7441E-01  9.9732E-01  6.9940E-01  1.0283E+00  7.5624E-01  1.0507E+00
             8.4795E-01
 PARAMETER:  1.1131E-01 -1.8558E-01  2.9093E-01  2.3001E-01  7.4075E-02  9.7317E-02 -2.5753E-01  1.2791E-01 -1.7940E-01  1.4945E-01
            -6.4931E-02
 GRADIENT:  -2.0636E+00  1.3096E+01  4.0717E+00  2.1479E+01 -7.1005E+00 -5.3002E-01 -2.6766E-01 -4.9376E-01 -2.6918E-01 -1.3478E+00
            -5.5554E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1752.09882799123        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  1.0099E+00  5.1016E-01  1.3977E+00  1.2945E+00  9.6080E-01  9.9754E-01  5.8143E-01  1.1667E+00  6.9519E-01  1.0888E+00
             8.4978E-01
 PARAMETER:  1.0987E-01 -5.7304E-01  4.3483E-01  3.5816E-01  6.0008E-02  9.7534E-02 -4.4226E-01  2.5419E-01 -2.6357E-01  1.8507E-01
            -6.2773E-02
 GRADIENT:   2.6144E+00  1.0121E+01  1.6964E+00  3.0363E+01 -9.9690E+00  1.2272E+00 -9.6212E-02  5.9133E-01  7.5606E-01  2.7435E+00
             1.1360E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1753.17666867917        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  1.0048E+00  2.7410E-01  1.7419E+00  1.4455E+00  9.9733E-01  9.8984E-01  3.6367E-01  1.4780E+00  6.2479E-01  1.1014E+00
             8.4807E-01
 PARAMETER:  1.0477E-01 -1.1943E+00  6.5495E-01  4.6843E-01  9.7329E-02  8.9787E-02 -9.1150E-01  4.9070E-01 -3.7035E-01  1.9659E-01
            -6.4787E-02
 GRADIENT:   2.5837E-02  2.4774E+00  1.1141E+00  7.8678E+00  4.9478E-01  9.3550E-02  1.2731E-02 -1.2242E-01  8.2718E-01 -1.6127E+00
             4.5935E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1754.01178108005        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1152
 NPARAMETR:  1.0020E+00  1.1889E-01  1.8309E+00  1.5470E+00  9.7332E-01  9.8623E-01  1.6677E-01  1.5567E+00  5.7295E-01  1.0965E+00
             8.4590E-01
 PARAMETER:  1.0202E-01 -2.0295E+00  7.0479E-01  5.3631E-01  7.2958E-02  8.6130E-02 -1.6911E+00  5.4258E-01 -4.5696E-01  1.9208E-01
            -6.7354E-02
 GRADIENT:  -5.9670E-01  2.1007E+00  1.1736E+00  2.2258E+01 -4.9462E+00 -2.5444E-01 -6.3593E-04 -4.8764E-01 -3.9774E+00  2.0426E-01
            -6.8705E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1754.53609511205        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1327
 NPARAMETR:  1.0006E+00  3.6812E-02  2.0196E+00  1.6020E+00  9.9505E-01  9.8359E-01  4.5771E-02  1.7282E+00  5.5389E-01  1.1156E+00
             8.4345E-01
 PARAMETER:  1.0056E-01 -3.2019E+00  8.0292E-01  5.7126E-01  9.5040E-02  8.3458E-02 -2.9841E+00  6.4710E-01 -4.9079E-01  2.0939E-01
            -7.0257E-02
 GRADIENT:  -1.1632E+00  4.5833E-01  9.5754E-01  1.1352E+01 -4.9862E+00 -6.7342E-01  1.8673E-06  7.0931E-01 -1.4545E+00  8.5705E-01
            -1.2028E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1754.70987529672        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1502
 NPARAMETR:  1.0007E+00  1.0000E-02  2.0720E+00  1.6205E+00  1.0030E+00  9.8474E-01  1.0000E-02  1.7679E+00  5.4892E-01  1.1181E+00
             8.4504E-01
 PARAMETER:  1.0071E-01 -4.5981E+00  8.2850E-01  5.8271E-01  1.0299E-01  8.4619E-02 -4.5923E+00  6.6977E-01 -4.9980E-01  2.1160E-01
            -6.8366E-02
 GRADIENT:   3.9341E-02  0.0000E+00 -2.0680E-01  1.0516E+01 -9.9247E-01 -6.1887E-03  0.0000E+00  4.5387E-01 -4.8552E-02  2.2818E-01
            -3.1276E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1754.72399567372        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1664
 NPARAMETR:  1.0006E+00  1.0000E-02  2.0244E+00  1.6162E+00  9.9317E-01  9.8479E-01  1.0000E-02  1.7244E+00  5.4942E-01  1.1096E+00
             8.4577E-01
 PARAMETER:  1.0057E-01 -4.5590E+00  8.0525E-01  5.8006E-01  9.3151E-02  8.4672E-02 -4.5473E+00  6.4486E-01 -4.9890E-01  2.0401E-01
            -6.7512E-02
 GRADIENT:  -2.0369E-02  0.0000E+00  1.0830E-02 -1.4168E-02  2.7011E-03  8.4935E-03  0.0000E+00 -7.9845E-03  1.2609E-02 -9.9974E-03
             1.2939E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1664
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.5118E-04 -4.1666E-06 -3.9851E-02 -1.0301E-02 -4.7493E-02
 SE:             2.9904E-02  2.2641E-06  1.8758E-02  2.8982E-02  2.0774E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9597E-01  6.5729E-02  3.3632E-02  7.2226E-01  2.2242E-02

 ETASHRINKSD(%)  1.0000E-10  9.9992E+01  3.7158E+01  2.9075E+00  3.0406E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.0509E+01  5.7304E+00  5.1566E+01
 EBVSHRINKSD(%)  2.9850E-01  9.9993E+01  4.1405E+01  3.3069E+00  2.5426E+01
 EBVSHRINKVR(%)  5.9611E-01  1.0000E+02  6.5667E+01  6.5045E+00  4.4387E+01
 RELATIVEINF(%)  9.7747E+01  2.2211E-08  9.9768E+00  4.3954E+00  1.2666E+01
 EPSSHRINKSD(%)  4.5902E+01
 EPSSHRINKVR(%)  7.0734E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1754.7239956737235     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1019.5731691099853     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.84
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1754.724       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  2.02E+00  1.62E+00  9.93E-01  9.85E-01  1.00E-02  1.72E+00  5.49E-01  1.11E+00  8.46E-01
 


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
+        1.13E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.35E+00  0.00E+00  3.43E+01
 
 TH 4
+       -1.62E+01  0.00E+00 -2.79E+01  1.32E+03
 
 TH 5
+        4.17E+00  0.00E+00 -9.19E+01 -8.24E+01  5.06E+02
 
 TH 6
+        2.70E+01  0.00E+00 -2.81E+00 -5.05E+00  1.00E+01  1.60E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -2.78E+00  0.00E+00 -1.42E+01 -5.60E+00 -9.00E+00  3.77E+00  0.00E+00  2.02E+01
 
 TH 9
+        1.41E+01  0.00E+00  4.37E+00 -5.64E-01  7.86E+00  8.09E+00  0.00E+00  1.80E+00  5.80E+02
 
 TH10
+       -4.16E+00  0.00E+00 -7.10E-01 -5.09E+00 -6.56E+01 -4.07E+00  0.00E+00  7.30E+00 -4.11E-01  5.92E+01
 
 TH11
+        1.64E+01  0.00E+00 -3.57E+00 -1.56E+01  6.01E+00 -3.30E+01  0.00E+00  5.19E+00  1.81E+01  4.08E+00  2.39E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.747
Stop Time:
Sat Sep 25 12:26:19 CDT 2021
