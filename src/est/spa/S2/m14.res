Wed Sep 29 17:09:36 CDT 2021
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
$DATA ../../../../data/spa/S2/dat14.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1706.34787744647        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0378E+02 -6.7896E+01 -5.9680E+01  5.7216E+00  8.0967E+01  4.0288E+01  4.4211E+00  1.4116E+01  3.5135E+01  1.0654E+01
             1.0949E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1716.81516752710        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0086E+00  1.1281E+00  1.1417E+00  9.8255E-01  1.0313E+00  1.0378E+00  9.6314E-01  9.2888E-01  8.4615E-01  9.2388E-01
             9.6660E-01
 PARAMETER:  1.0861E-01  2.2052E-01  2.3248E-01  8.2397E-02  1.3081E-01  1.3710E-01  6.2442E-02  2.6228E-02 -6.7055E-02  2.0829E-02
             6.6025E-02
 GRADIENT:  -6.1278E+00  5.2981E+00  1.6395E+01 -1.3290E+01 -1.3173E+01  7.0855E+00 -1.5798E+00  8.5601E-01 -6.7139E+00 -1.1330E+01
            -1.6568E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1717.79185736782        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  1.0099E+00  1.0892E+00  1.1911E+00  1.0076E+00  1.0634E+00  1.0084E+00  9.7507E-01  4.3821E-01  9.3083E-01  1.1105E+00
             1.0030E+00
 PARAMETER:  1.0990E-01  1.8543E-01  2.7491E-01  1.0762E-01  1.6146E-01  1.0840E-01  7.4754E-02 -7.2506E-01  2.8319E-02  2.0483E-01
             1.0303E-01
 GRADIENT:  -2.6212E+00 -2.0987E+00  1.3371E+01 -9.7004E+00 -3.7942E+00 -3.8571E+00  6.0380E+00 -3.7540E-01  1.0634E+01  2.3150E+00
             1.5144E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1719.62397777710        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  1.0146E+00  1.1486E+00  9.3772E-01  9.6696E-01  9.7292E-01  1.0170E+00  9.6226E-01  3.4592E-01  8.8665E-01  9.7532E-01
             9.8975E-01
 PARAMETER:  1.1448E-01  2.3852E-01  3.5695E-02  6.6402E-02  7.2545E-02  1.1684E-01  6.1533E-02 -9.6156E-01 -2.0308E-02  7.5009E-02
             8.9694E-02
 GRADIENT:   3.1679E+00  3.7816E+00  1.4096E+00  3.7544E+00 -5.0120E+00 -1.2913E+00  3.6986E-01  5.7140E-01  1.2635E+00  1.3422E+00
             4.8353E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1719.90736673547        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  1.0144E+00  1.3717E+00  8.0836E-01  8.2074E-01  1.0261E+00  1.0210E+00  8.3909E-01  1.1008E-01  9.9778E-01  9.7587E-01
             9.9039E-01
 PARAMETER:  1.1426E-01  4.1607E-01 -1.1275E-01 -9.7548E-02  1.2578E-01  1.2080E-01 -7.5435E-02 -2.1065E+00  9.7776E-02  7.5577E-02
             9.0341E-02
 GRADIENT:   2.8050E-01  4.6885E-01 -1.1917E+00  2.7957E+00  2.5901E+00 -2.5735E-01  1.5762E-01  5.8306E-02  5.2776E-01 -2.3888E-02
             1.9223E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1719.93799364762        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      892
 NPARAMETR:  1.0149E+00  1.3722E+00  8.0714E-01  8.1966E-01  1.0246E+00  1.0221E+00  8.3870E-01  3.3590E-02  9.9545E-01  9.7665E-01
             9.9017E-01
 PARAMETER:  1.1479E-01  4.1641E-01 -1.1426E-01 -9.8864E-02  1.2429E-01  1.2182E-01 -7.5899E-02 -3.2935E+00  9.5441E-02  7.6376E-02
             9.0126E-02
 GRADIENT:   1.4342E+00  6.7754E-01  7.0099E-02  9.0934E-01  2.2534E-01  1.5504E-01  2.9278E-02  5.3239E-03  1.5438E-02  1.2347E-01
             4.3766E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1719.94063054078        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1070             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0156E+00  1.3729E+00  8.0402E-01  8.1758E-01  1.0236E+00  1.0222E+00  8.3807E-01  1.0000E-02  9.9603E-01  9.7408E-01
             9.8995E-01
 PARAMETER:  1.1550E-01  4.1690E-01 -1.1813E-01 -1.0141E-01  1.2333E-01  1.2198E-01 -7.6652E-02 -4.7298E+00  9.6019E-02  7.3733E-02
             8.9900E-02
 GRADIENT:   5.1347E+02  2.6653E+02  1.5691E+00  4.6447E+01  9.2729E+00  5.9184E+01  4.1043E+00  0.0000E+00  3.8758E+00  7.4148E-01
             7.6920E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1719.94113128826        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1257
 NPARAMETR:  1.0153E+00  1.3722E+00  8.0415E-01  8.1798E-01  1.0231E+00  1.0222E+00  8.3848E-01  1.0000E-02  9.9528E-01  9.7386E-01
             9.8995E-01
 PARAMETER:  1.1521E-01  4.1639E-01 -1.1797E-01 -1.0091E-01  1.2281E-01  1.2198E-01 -7.6164E-02 -4.7298E+00  9.5268E-02  7.3513E-02
             8.9903E-02
 GRADIENT:   2.3213E+00 -1.4251E+00  3.1448E-01 -1.4004E+00 -2.1509E-01  2.2015E-01 -5.1729E-02  0.0000E+00 -8.0560E-02  2.8262E-02
            -2.5660E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1719.94162876022        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1358
 NPARAMETR:  1.0153E+00  1.3714E+00  8.0457E-01  8.1845E-01  1.0224E+00  1.0222E+00  8.3904E-01  1.0000E-02  9.9468E-01  9.7338E-01
             9.8987E-01
 PARAMETER:  1.1519E-01  4.1649E-01 -1.1864E-01 -9.9903E-02  1.2290E-01  1.2196E-01 -7.5843E-02 -4.7298E+00  9.5668E-02  7.3404E-02
             8.9942E-02
 GRADIENT:  -2.1848E-02  6.3474E-01 -1.3964E-01  3.4460E-01  6.9353E-01 -5.9440E-03 -2.4380E-02  0.0000E+00  4.6102E-02  3.8108E-02
             3.2729E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1358
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.8304E-04 -1.8128E-02 -3.1510E-04  9.0905E-03 -2.5477E-02
 SE:             2.9840E-02  2.1334E-02  1.2889E-04  2.3217E-02  2.3669E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9511E-01  3.9547E-01  1.4499E-02  6.9539E-01  2.8174E-01

 ETASHRINKSD(%)  3.3406E-02  2.8529E+01  9.9568E+01  2.2222E+01  2.0707E+01
 ETASHRINKVR(%)  6.6802E-02  4.8919E+01  9.9998E+01  3.9505E+01  3.7127E+01
 EBVSHRINKSD(%)  3.9833E-01  2.8165E+01  9.9599E+01  2.3281E+01  1.8935E+01
 EBVSHRINKVR(%)  7.9506E-01  4.8397E+01  9.9998E+01  4.1142E+01  3.4285E+01
 RELATIVEINF(%)  9.8927E+01  1.6668E+00  1.4118E-04  2.0697E+00  8.6090E+00
 EPSSHRINKSD(%)  4.3282E+01
 EPSSHRINKVR(%)  6.7830E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1719.9416287602241     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -984.79080219648597     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.25
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1719.942       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.37E+00  8.04E-01  8.19E-01  1.02E+00  1.02E+00  8.39E-01  1.00E-02  9.96E-01  9.74E-01  9.90E-01
 


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
+        1.02E+03
 
 TH 2
+       -6.46E+00  4.52E+02
 
 TH 3
+        1.11E+01  1.45E+02  3.08E+02
 
 TH 4
+       -1.39E+01  4.42E+02 -2.06E+02  9.92E+02
 
 TH 5
+       -3.21E+00 -2.48E+02 -3.82E+02  2.35E+02  7.35E+02
 
 TH 6
+       -2.00E-02 -1.33E+00  2.85E+00 -3.76E+00 -1.26E+00  1.88E+02
 
 TH 7
+        8.72E-01  1.81E+01  4.22E+00 -1.52E+01 -8.81E+00  1.89E-01  8.41E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.21E+00 -2.16E+01 -2.66E+01  4.20E+01  5.37E+00 -1.87E-01  2.79E+01  0.00E+00  7.84E+01
 
 TH10
+       -3.36E-01 -1.11E+01 -3.82E+01 -9.59E+00 -5.70E+01 -4.99E-02  1.69E+01  0.00E+00  8.33E+00  8.80E+01
 
 TH11
+       -7.06E+00 -2.20E+01 -3.14E+01  2.76E+00  4.30E+00  1.46E+00  6.40E+00  0.00E+00  1.13E+01  2.15E+01  2.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.205
Stop Time:
Wed Sep 29 17:10:02 CDT 2021
