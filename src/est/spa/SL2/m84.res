Sat Sep 18 12:30:09 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat84.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1684.41017540432        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.1965E+01 -6.2761E+01 -3.6937E+01 -5.7208E+01  2.9317E+01  4.3050E+01 -1.2355E+01  1.1716E+01 -1.3311E+01  3.0080E+00
            -3.0914E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1693.23551530611        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0614E+00  1.0849E+00  1.2847E+00  9.7459E-01  1.1192E+00  8.3153E-01  1.1168E+00  8.5619E-01  1.0902E+00  1.0025E+00
             1.1312E+00
 PARAMETER:  1.5957E-01  1.8144E-01  3.5056E-01  7.4263E-02  2.1264E-01 -8.4485E-02  2.1048E-01 -5.5265E-02  1.8635E-01  1.0248E-01
             2.2330E-01
 GRADIENT:   1.0951E+02 -9.8234E+00  2.2960E+01 -4.0494E+01 -1.0432E+01 -2.0602E+01  6.0931E+00 -2.0496E+00  3.6023E+00 -1.5326E+01
             8.6502E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1695.18796937417        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0503E+00  8.7538E-01  1.3607E+00  1.1400E+00  1.0886E+00  8.5541E-01  1.0320E+00  4.5146E-01  1.0754E+00  1.1819E+00
             1.1450E+00
 PARAMETER:  1.4906E-01 -3.3098E-02  4.0799E-01  2.3105E-01  1.8487E-01 -5.6170E-02  1.3151E-01 -6.9527E-01  1.7267E-01  2.6708E-01
             2.3543E-01
 GRADIENT:   7.4124E+01 -1.4367E+00  9.3205E+00  1.5054E+01  6.7871E-01 -6.7313E+00  6.0284E-01 -8.5014E-01  1.6042E+01  2.2027E+00
             1.1287E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1697.03894130678        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0320E+00  9.7468E-01  1.1273E+00  1.0632E+00  1.0358E+00  8.6826E-01  1.1452E+00  4.1455E-01  1.0222E+00  1.0853E+00
             1.0916E+00
 PARAMETER:  1.3151E-01  7.4356E-02  2.1983E-01  1.6129E-01  1.3518E-01 -4.1261E-02  2.3560E-01 -7.8057E-01  1.2192E-01  1.8186E-01
             1.8765E-01
 GRADIENT:   1.5898E+01 -4.4681E+00 -1.1088E+00  1.5963E+00  5.4517E+00 -2.2508E+00  5.7827E-02  1.4469E-01  2.1547E+00  7.1077E-01
            -7.0464E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1697.34028054030        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      388
 NPARAMETR:  1.0444E+00  9.6551E-01  1.0904E+00  1.0753E+00  1.0095E+00  8.7855E-01  1.1929E+00  3.6822E-01  9.9866E-01  1.0603E+00
             1.0931E+00
 PARAMETER:  1.4343E-01  6.4906E-02  1.8650E-01  1.7262E-01  1.0950E-01 -2.9487E-02  2.7637E-01 -8.9908E-01  9.8656E-02  1.5857E-01
             1.8902E-01
 GRADIENT:   1.6877E+00  4.2085E-01 -1.0588E+00  1.1515E+00  3.3106E-01  6.7790E-02  1.4589E-01  1.9463E-01 -1.2040E-02  4.3200E-01
             3.3983E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1697.38264055030        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      563
 NPARAMETR:  1.0440E+00  9.3265E-01  1.0812E+00  1.0932E+00  9.9193E-01  8.7805E-01  1.2233E+00  2.0579E-01  9.8340E-01  1.0577E+00
             1.0941E+00
 PARAMETER:  1.4306E-01  3.0272E-02  1.7803E-01  1.8911E-01  9.1896E-02 -3.0052E-02  3.0158E-01 -1.4809E+00  8.3265E-02  1.5614E-01
             1.8995E-01
 GRADIENT:   1.1315E+00 -6.4868E-01  4.4761E-03 -1.0614E+00 -3.6395E-01 -3.1234E-02 -3.3635E-02  3.7927E-02 -2.9399E-01  7.4685E-02
            -3.0146E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1697.41318750305        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      739
 NPARAMETR:  1.0444E+00  1.0247E+00  1.0458E+00  1.0353E+00  1.0167E+00  8.7902E-01  1.1487E+00  6.5394E-02  1.0243E+00  1.0612E+00
             1.0949E+00
 PARAMETER:  1.4345E-01  1.2440E-01  1.4476E-01  1.3473E-01  1.1657E-01 -2.8949E-02  2.3861E-01 -2.6273E+00  1.2405E-01  1.5944E-01
             1.9068E-01
 GRADIENT:   1.2983E-01 -1.8466E-01 -2.3654E-01 -1.2172E-01  3.4924E-01  4.6374E-02 -1.8404E-03  5.4242E-03 -1.2959E-02  2.0821E-02
             3.1311E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1697.41609360619        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  1.0443E+00  1.0120E+00  1.0495E+00  1.0434E+00  1.0126E+00  8.7879E-01  1.1584E+00  1.3527E-02  1.0188E+00  1.0608E+00
             1.0950E+00
 PARAMETER:  1.4332E-01  1.1196E-01  1.4832E-01  1.4245E-01  1.1249E-01 -2.9208E-02  2.4703E-01 -4.2030E+00  1.1859E-01  1.5907E-01
             1.9075E-01
 GRADIENT:  -4.3626E-03 -8.5839E-02 -1.9251E-02 -9.3200E-02  4.6048E-02 -7.8751E-03  1.7706E-03  2.1383E-04  1.0514E-02 -1.3931E-03
            -2.4032E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1697.41615184076        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1007
 NPARAMETR:  1.0443E+00  1.0132E+00  1.0491E+00  1.0427E+00  1.0129E+00  8.7882E-01  1.1575E+00  1.0000E-02  1.0192E+00  1.0609E+00
             1.0950E+00
 PARAMETER:  1.4333E-01  1.1314E-01  1.4794E-01  1.4181E-01  1.1278E-01 -2.9179E-02  2.4625E-01 -4.5350E+00  1.1904E-01  1.5910E-01
             1.9076E-01
 GRADIENT:   2.4824E-04  4.4820E-04  3.1064E-04 -6.2709E-04 -9.7148E-04 -7.1824E-04  6.4914E-04  0.0000E+00  1.1967E-04  4.7613E-04
             4.3503E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1007
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.6693E-04 -9.6602E-03 -3.1160E-04  3.1265E-04 -2.5694E-02
 SE:             2.9763E-02  1.8657E-02  1.4491E-04  2.4643E-02  2.3592E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8480E-01  6.0460E-01  3.1537E-02  9.8988E-01  2.7611E-01

 ETASHRINKSD(%)  2.9040E-01  3.7498E+01  9.9515E+01  1.7442E+01  2.0964E+01
 ETASHRINKVR(%)  5.7996E-01  6.0935E+01  9.9998E+01  3.1842E+01  3.7532E+01
 EBVSHRINKSD(%)  6.4106E-01  3.7070E+01  9.9542E+01  1.7779E+01  1.8823E+01
 EBVSHRINKVR(%)  1.2780E+00  6.0398E+01  9.9998E+01  3.2397E+01  3.4102E+01
 RELATIVEINF(%)  9.7886E+01  1.5031E+00  3.2100E-04  3.1098E+00  7.9634E+00
 EPSSHRINKSD(%)  4.1995E+01
 EPSSHRINKVR(%)  6.6355E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1697.4161518407595     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -962.26532527702136     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1697.416       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.01E+00  1.05E+00  1.04E+00  1.01E+00  8.79E-01  1.16E+00  1.00E-02  1.02E+00  1.06E+00  1.10E+00
 


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
+        1.30E+03
 
 TH 2
+       -1.06E+01  3.35E+02
 
 TH 3
+        1.70E+01  1.12E+02  2.24E+02
 
 TH 4
+       -1.42E+01  3.38E+02 -1.03E+02  6.68E+02
 
 TH 5
+       -6.45E+00 -2.20E+02 -3.08E+02  1.24E+02  6.32E+02
 
 TH 6
+       -3.17E+00 -8.42E-01  3.97E+00 -3.16E+00 -3.24E-01  2.48E+02
 
 TH 7
+        2.61E+00  1.54E+01  9.68E+00 -6.98E+00 -1.17E+01  7.36E-01  2.95E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.52E+00 -2.18E+01 -1.58E+01  2.22E+01  5.42E+00 -1.39E+00  2.04E+01  0.00E+00  9.57E+01
 
 TH10
+       -1.55E+00 -1.81E+00 -2.57E+01 -1.10E+01 -4.87E+01 -7.52E-01  6.41E+00  0.00E+00  3.94E+00  7.75E+01
 
 TH11
+       -9.33E+00 -1.75E+01 -3.52E+01 -6.22E+00  8.41E+00  5.68E+00  5.11E+00  0.00E+00  6.06E+00  2.06E+01  1.89E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.420
Stop Time:
Sat Sep 18 12:30:28 CDT 2021
