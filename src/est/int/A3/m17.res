Fri Sep 24 22:06:08 CDT 2021
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
$DATA ../../../../data/int/A3/dat17.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1205.51209805510        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3047E+01  1.1603E+02  2.1974E+02 -1.2680E+02  7.2134E+01  2.1964E+01 -1.7448E+02 -8.0126E+01 -2.2497E+01 -1.6347E+02
            -4.9137E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2832.17974434183        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0253E+00  9.5549E-01  8.4122E-01  1.0733E+00  9.3030E-01  8.7046E-01  1.4381E+00  8.3602E-01  9.2633E-01  1.3756E+00
             2.0630E+00
 PARAMETER:  1.2499E-01  5.4471E-02 -7.2904E-02  1.7076E-01  2.7753E-02 -3.8730E-02  4.6330E-01 -7.9098E-02  2.3479E-02  4.1889E-01
             8.2416E-01
 GRADIENT:   2.0833E+01  1.1221E+01  4.1038E+00 -3.5500E+01  1.2764E+00 -2.7255E+01 -9.0761E+00  2.0195E+01  1.0786E+00  1.5298E+01
            -3.5929E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2854.78729505862        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0178E+00  6.5790E-01  5.2741E-01  1.3082E+00  5.5627E-01  9.2545E-01  1.8572E+00  1.7421E-01  1.0267E+00  9.8104E-01
             2.0301E+00
 PARAMETER:  1.1764E-01 -3.1871E-01 -5.3977E-01  3.6867E-01 -4.8651E-01  2.2529E-02  7.1909E-01 -1.6475E+00  1.2636E-01  8.0861E-02
             8.0806E-01
 GRADIENT:   6.5297E-02  9.0940E+01 -3.2560E+01  1.3075E+02  1.7549E+01 -5.4869E+00  1.3849E+01  1.9721E+00  9.0088E+00  2.2090E+01
            -3.2975E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2899.03383098092        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0244E+00  3.6777E-01  3.2680E-01  1.3501E+00  3.5178E-01  9.9448E-01  1.9432E+00  3.9108E-02  1.0538E+00  6.6250E-01
             2.3621E+00
 PARAMETER:  1.2415E-01 -9.0030E-01 -1.0184E+00  4.0018E-01 -9.4474E-01  9.4460E-02  7.6435E-01 -3.1414E+00  1.5236E-01 -3.1173E-01
             9.5955E-01
 GRADIENT:   1.6921E+00  1.9349E+01 -5.8130E+00  1.6935E+02  4.4521E+01  2.0082E+01  1.6161E+01  4.9818E-02 -1.2651E+01 -4.0237E+00
             1.1215E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2934.73066550246        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0224E+00  1.7288E-01  1.1930E-01  1.0171E+00  1.5900E-01  9.4135E-01  1.7228E+00  1.0000E-02  1.4598E+00  7.0422E-01
             2.1501E+00
 PARAMETER:  1.2216E-01 -1.6551E+00 -2.0261E+00  1.1693E-01 -1.7388E+00  3.9555E-02  6.4394E-01 -6.8283E+00  4.7833E-01 -2.5067E-01
             8.6552E-01
 GRADIENT:   2.1857E+00 -5.2018E-01  1.5372E+02  7.8389E+01 -2.3222E+02  2.2167E+00 -9.1537E+00  0.0000E+00 -3.4067E+01 -5.2144E+01
             8.5292E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2949.03029566569        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.0137E+00  1.4932E-01  9.2349E-02  8.6876E-01  1.4209E-01  8.9021E-01  1.8392E+00  1.0000E-02  1.7866E+00  8.2649E-01
             1.9555E+00
 PARAMETER:  1.1359E-01 -1.8017E+00 -2.2822E+00 -4.0684E-02 -1.8513E+00 -1.6297E-02  7.0934E-01 -8.0967E+00  6.8029E-01 -9.0564E-02
             7.7063E-01
 GRADIENT:  -4.5240E+00 -7.6885E+00  4.2749E+01  1.6580E+01 -6.1489E+01 -1.5980E+01  5.1545E+00  0.0000E+00 -2.0407E+01 -3.3491E+00
            -8.0221E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2949.55006322709        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  1.0139E+00  1.4894E-01  9.1723E-02  8.6429E-01  1.4196E-01  8.9374E-01  1.8342E+00  1.0000E-02  1.7976E+00  8.2626E-01
             1.9619E+00
 PARAMETER:  1.1385E-01 -1.8042E+00 -2.2890E+00 -4.5852E-02 -1.8522E+00 -1.2342E-02  7.0662E-01 -8.1463E+00  6.8644E-01 -9.0848E-02
             7.7390E-01
 GRADIENT:  -3.7999E+00 -7.7252E+00  3.7387E+01  1.5030E+01 -5.4286E+01 -1.4269E+01  4.5289E+00  0.0000E+00 -1.9162E+01 -3.0020E+00
            -7.3294E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2959.50394237063        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      626
 NPARAMETR:  1.0214E+00  1.6571E-01  1.0288E-01  8.9803E-01  1.5926E-01  9.3354E-01  1.7839E+00  1.0000E-02  1.7627E+00  7.9177E-01
             2.0163E+00
 PARAMETER:  1.2121E-01 -1.6975E+00 -2.1742E+00 -7.5528E-03 -1.7372E+00  3.1228E-02  6.7883E-01 -7.7581E+00  6.6686E-01 -1.3348E-01
             8.0126E-01
 GRADIENT:   4.3898E+00 -6.1128E-01 -6.4118E+00 -2.8759E+00  4.3085E-01  9.7873E-01 -1.7620E+00  0.0000E+00 -6.9782E+00 -1.9002E+00
            -7.4660E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2959.66502041080        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  1.0192E+00  1.6693E-01  1.0420E-01  9.0909E-01  1.6028E-01  9.3104E-01  1.7950E+00  1.0000E-02  1.8035E+00  7.9352E-01
             2.0215E+00
 PARAMETER:  1.1905E-01 -1.6902E+00 -2.1614E+00  4.6895E-03 -1.7308E+00  2.8552E-02  6.8502E-01 -7.7060E+00  6.8974E-01 -1.3127E-01
             8.0382E-01
 GRADIENT:   1.8547E-03  2.4291E-02 -6.9035E-02  4.7048E-02  2.9814E-02 -5.4330E-03  7.1210E-03  0.0000E+00 -2.1815E-02 -6.1069E-03
             5.0087E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -2959.66502041080        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  1.0192E+00  1.6693E-01  1.0420E-01  9.0909E-01  1.6028E-01  9.3104E-01  1.7950E+00  1.0000E-02  1.8035E+00  7.9352E-01
             2.0215E+00
 PARAMETER:  1.1905E-01 -1.6902E+00 -2.1614E+00  4.6895E-03 -1.7308E+00  2.8552E-02  6.8502E-01 -7.7060E+00  6.8974E-01 -1.3127E-01
             8.0382E-01
 GRADIENT:   1.8547E-03  2.4291E-02 -6.9035E-02  4.7048E-02  2.9814E-02 -5.4330E-03  7.1210E-03  0.0000E+00 -2.1815E-02 -6.1069E-03
             5.0087E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      823
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3119E-03  1.2275E-02  3.5247E-04 -5.0294E-03  8.9339E-03
 SE:             2.9384E-02  2.6451E-02  2.5085E-04  2.6972E-02  2.8113E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3729E-01  6.4260E-01  1.5999E-01  8.5208E-01  7.5065E-01

 ETASHRINKSD(%)  1.5603E+00  1.1386E+01  9.9160E+01  9.6417E+00  5.8175E+00
 ETASHRINKVR(%)  3.0962E+00  2.1475E+01  9.9993E+01  1.8354E+01  1.1297E+01
 EBVSHRINKSD(%)  1.5778E+00  9.3413E+00  9.9208E+01  6.1951E+00  7.2995E+00
 EBVSHRINKVR(%)  3.1308E+00  1.7810E+01  9.9994E+01  1.2006E+01  1.4066E+01
 RELATIVEINF(%)  9.6821E+01  5.1949E+01  1.4424E-03  6.0570E+01  2.2076E+01
 EPSSHRINKSD(%)  2.0894E+01
 EPSSHRINKVR(%)  3.7422E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2959.6650204108023     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1305.5756606423915     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.44
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2959.665       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.67E-01  1.04E-01  9.09E-01  1.60E-01  9.31E-01  1.80E+00  1.00E-02  1.80E+00  7.94E-01  2.02E+00
 


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
+        1.20E+03
 
 TH 2
+        1.97E+00  1.04E+04
 
 TH 3
+        5.30E+01 -5.84E+03  6.47E+04
 
 TH 4
+        1.49E+01 -9.08E+01 -1.25E+03  3.74E+02
 
 TH 5
+       -2.66E+01 -3.00E+03 -4.44E+04 -2.82E+02  5.96E+04
 
 TH 6
+        2.12E+00 -7.62E+00  3.14E+01 -1.48E+01 -2.08E+01  2.09E+02
 
 TH 7
+       -1.51E-01  5.83E+01 -4.61E+00 -1.71E+00 -3.50E+01 -3.25E-01  3.93E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.29E+01  1.39E+01  2.24E+02 -1.19E+01  7.70E+01  1.63E-01  2.10E-01  0.00E+00  3.89E+01
 
 TH10
+       -2.49E+00  4.68E+00  2.73E+02  9.99E+00  5.19E+01 -4.14E+00  5.58E+00  0.00E+00  1.81E+00  2.42E+02
 
 TH11
+       -2.06E+01 -2.76E+01  7.22E+00  4.64E+00  7.75E+01  2.82E+00  5.29E+00  0.00E+00  7.18E+00  5.07E+00  2.57E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.791
Stop Time:
Fri Sep 24 22:06:50 CDT 2021
