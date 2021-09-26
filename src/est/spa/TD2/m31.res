Sat Sep 25 13:27:57 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1667.04431393661        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0742E+02 -5.7164E+01 -6.0636E+01 -3.2501E+00  9.2305E+01 -1.2186E+00  1.0341E+01  9.5841E+00  1.6814E+01  1.1828E+00
             1.2882E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1674.06128299063        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.6084E-01  1.0009E+00  1.1468E+00  1.0037E+00  9.8227E-01  9.9765E-01  8.3927E-01  9.1064E-01  9.3044E-01  9.8304E-01
             1.0194E+00
 PARAMETER:  6.0049E-02  1.0086E-01  2.3696E-01  1.0372E-01  8.2107E-02  9.7644E-02 -7.5226E-02  6.3962E-03  2.7897E-02  8.2893E-02
             1.1922E-01
 GRADIENT:   1.9393E+01 -1.1905E+01  6.5895E+00 -3.1215E+01 -6.4801E+00  5.9456E-01  2.0346E+00 -1.4677E-01 -9.5281E+00 -3.2575E+00
             1.5487E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1675.22610752572        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  9.6801E-01  9.4164E-01  1.1584E+00  1.0526E+00  9.5145E-01  1.0080E+00  6.0594E-01  8.4428E-01  1.0326E+00  1.0538E+00
             9.7606E-01
 PARAMETER:  6.7483E-02  3.9864E-02  2.4706E-01  1.5129E-01  5.0237E-02  1.0800E-01 -4.0097E-01 -6.9277E-02  1.3208E-01  1.5239E-01
             7.5773E-02
 GRADIENT:   4.0826E+01  5.9129E+00  1.0869E+01  6.9045E+00 -2.9312E+01  5.1398E+00  2.8470E+00 -8.6779E-01  1.0153E+01  5.8485E+00
             3.8371E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1676.05621814551        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.5308E-01  1.0471E+00  1.2469E+00  9.8457E-01  1.0393E+00  9.9935E-01  4.6373E-01  1.0875E+00  1.0987E+00  1.0833E+00
             9.7498E-01
 PARAMETER:  5.1939E-02  1.4605E-01  3.2065E-01  8.4454E-02  1.3859E-01  9.9351E-02 -6.6844E-01  1.8390E-01  1.9411E-01  1.8003E-01
             7.4658E-02
 GRADIENT:   6.6082E+00  7.5493E-01 -6.2639E-01  1.4659E+00 -3.6338E+00  8.9332E-01  2.2610E+00  8.5794E-01  4.5605E+00  3.2618E+00
            -1.2118E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1676.70247083966        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  9.4942E-01  1.0093E+00  1.3641E+00  1.0032E+00  1.0621E+00  9.9726E-01  8.9233E-02  1.1692E+00  1.1140E+00  1.0957E+00
             9.8139E-01
 PARAMETER:  4.8093E-02  1.0924E-01  4.1051E-01  1.0321E-01  1.6020E-01  9.7261E-02 -2.3165E+00  2.5635E-01  2.0800E-01  1.9137E-01
             8.1215E-02
 GRADIENT:  -1.4237E+00 -2.3640E+00  1.2444E+00 -2.4576E+00 -1.3341E+00 -8.5474E-02  1.4812E-01 -2.9400E-01  1.5207E+00 -2.3686E-02
             1.7570E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1677.38712834740        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      464
 NPARAMETR:  9.7152E-01  1.1116E+00  1.3695E+00  9.4205E-01  1.1041E+00  1.0098E+00  1.6548E-02  1.2604E+00  1.1981E+00  1.1212E+00
             9.8518E-01
 PARAMETER:  7.1103E-02  2.0581E-01  4.1446E-01  4.0306E-02  1.9900E-01  1.0980E-01 -4.0015E+00  3.3142E-01  2.8076E-01  2.1437E-01
             8.5066E-02
 GRADIENT:   4.6056E+00  2.1638E+00  1.9435E-01  3.3564E+00 -1.9619E-01  7.9599E-01  2.6604E-03 -2.2252E-01 -1.3025E-01  4.9177E-01
             2.1585E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1677.41278389731        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      639
 NPARAMETR:  9.6975E-01  1.1519E+00  1.3545E+00  9.1211E-01  1.1146E+00  1.0081E+00  1.6507E-02  1.2889E+00  1.2360E+00  1.1232E+00
             9.8546E-01
 PARAMETER:  6.9278E-02  2.4139E-01  4.0343E-01  8.0069E-03  2.0851E-01  1.0803E-01 -4.0040E+00  3.5376E-01  3.1189E-01  2.1615E-01
             8.5352E-02
 GRADIENT:   2.5614E-01  1.6044E-01  3.1960E-02  1.3336E-01 -7.9106E-02  5.3325E-02  2.6287E-03 -6.1027E-03 -1.4538E-02  4.0491E-02
             2.0845E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1677.41364388994        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      817
 NPARAMETR:  9.6963E-01  1.1517E+00  1.3542E+00  9.1210E-01  1.1145E+00  1.0079E+00  1.0000E-02  1.2887E+00  1.2360E+00  1.1229E+00
             9.8540E-01
 PARAMETER:  6.9163E-02  2.4126E-01  4.0323E-01  7.9973E-03  2.0843E-01  1.0789E-01 -4.5857E+00  3.5363E-01  3.1190E-01  2.1591E-01
             8.5293E-02
 GRADIENT:   8.6579E-03  8.9922E-03  2.7275E-03  4.0165E-04 -7.4446E-03  2.0778E-03  0.0000E+00 -2.4549E-04 -4.0003E-03 -2.6524E-04
             2.3670E-04

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1677.41364388994        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  9.6963E-01  1.1517E+00  1.3542E+00  9.1210E-01  1.1145E+00  1.0079E+00  1.0000E-02  1.2887E+00  1.2360E+00  1.1229E+00
             9.8540E-01
 PARAMETER:  6.9163E-02  2.4126E-01  4.0323E-01  7.9973E-03  2.0843E-01  1.0789E-01 -4.5857E+00  3.5363E-01  3.1190E-01  2.1591E-01
             8.5293E-02
 GRADIENT:   8.6579E-03  8.9922E-03  2.7275E-03  4.0165E-04 -7.4446E-03  2.0778E-03  0.0000E+00 -2.4549E-04 -4.0003E-03 -2.6524E-04
             2.3670E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      839
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.4177E-04 -9.7165E-04 -2.8083E-02 -2.6659E-03 -2.9045E-02
 SE:             2.9855E-02  2.7718E-04  1.2893E-02  2.9150E-02  2.3785E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8018E-01  4.5578E-04  2.9394E-02  9.2713E-01  2.2202E-01

 ETASHRINKSD(%)  1.0000E-10  9.9071E+01  5.6807E+01  2.3437E+00  2.0319E+01
 ETASHRINKVR(%)  1.0000E-10  9.9991E+01  8.1344E+01  4.6325E+00  3.6509E+01
 EBVSHRINKSD(%)  4.1163E-01  9.9198E+01  6.0401E+01  2.7903E+00  1.7734E+01
 EBVSHRINKVR(%)  8.2156E-01  9.9994E+01  8.4319E+01  5.5027E+00  3.2324E+01
 RELATIVEINF(%)  9.9053E+01  7.1887E-04  5.4393E+00  1.2323E+01  2.2996E+01
 EPSSHRINKSD(%)  4.4015E+01
 EPSSHRINKVR(%)  6.8657E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1677.4136438899384     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -942.26281732620021     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1677.414       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.70E-01  1.15E+00  1.35E+00  9.12E-01  1.11E+00  1.01E+00  1.00E-02  1.29E+00  1.24E+00  1.12E+00  9.85E-01
 


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
+        1.15E+03
 
 TH 2
+       -1.56E+01  5.97E+02
 
 TH 3
+        2.47E+00  6.08E+01  6.26E+01
 
 TH 4
+       -1.46E+01  5.74E+02 -7.37E+00  8.29E+02
 
 TH 5
+       -1.57E+00 -2.00E+02 -1.23E+02 -3.46E+01  4.89E+02
 
 TH 6
+        3.80E+00 -1.45E+00  1.22E+00 -7.71E-01 -1.24E+00  1.93E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.01E-01 -1.67E+01 -1.77E+01 -9.89E-01  1.35E+00  3.91E-03  0.00E+00  1.44E+01
 
 TH 9
+        2.59E+00 -1.06E+02  1.74E+00  7.08E+00  2.58E+00 -4.12E-01  0.00E+00 -9.19E-02  1.15E+02
 
 TH10
+        8.31E-01 -6.28E+00 -6.10E+00 -3.40E+00 -4.81E+01 -5.99E-01  0.00E+00  6.58E+00  6.95E-01  8.05E+01
 
 TH11
+       -9.46E+00 -3.21E+01 -9.81E+00 -1.30E+01  1.86E+00  2.49E+00  0.00E+00  4.54E+00  9.73E+00  1.23E+01  2.21E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.245
Stop Time:
Sat Sep 25 13:28:14 CDT 2021
