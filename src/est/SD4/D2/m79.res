Sun Oct 24 04:42:53 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1517.95390383054        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.1728E+02 -8.1521E+01 -4.1702E+01 -5.8288E+01  4.4105E+01 -1.9149E+01 -1.2586E+02 -1.4700E+01 -1.5236E+02 -1.0707E+01
             2.0207E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1604.67672236524        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      130
 NPARAMETR:  1.1112E+00  1.3265E+00  1.0832E+00  9.6800E-01  1.3048E+00  1.1593E+00  2.1267E+00  1.1125E+00  1.6165E+00  9.5340E-01
             9.1079E-01
 PARAMETER:  2.0543E-01  3.8252E-01  1.7990E-01  6.7473E-02  3.6604E-01  2.4784E-01  8.5457E-01  2.0664E-01  5.8026E-01  5.2279E-02
             6.5531E-03
 GRADIENT:   9.1342E+01 -6.8582E+00 -5.7503E+01  9.7184E+01  1.3547E+02 -9.9931E+00  2.0674E+01 -9.4326E+00  3.7213E+01 -1.5603E+01
            -9.1678E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1617.48520528928        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  1.0868E+00  1.0987E+00  1.5689E+00  1.1935E+00  1.1795E+00  1.1470E+00  2.8375E+00  1.9374E+00  9.9974E-01  7.6651E-01
             9.5736E-01
 PARAMETER:  1.8325E-01  1.9413E-01  5.5040E-01  2.7690E-01  2.6508E-01  2.3713E-01  1.1429E+00  7.6136E-01  9.9744E-02 -1.6591E-01
             5.6422E-02
 GRADIENT:   5.7850E+01  5.2317E+01 -2.0708E+01  1.0982E+02  3.0850E+01 -1.0780E+01  1.7725E+01  5.5663E+00  4.9963E+00 -1.6052E+01
             9.3113E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1630.50300121291        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  1.0482E+00  1.0061E+00  1.9301E+00  1.0890E+00  1.2764E+00  1.1984E+00  2.5511E+00  2.0883E+00  1.0029E+00  1.0216E+00
             9.2353E-01
 PARAMETER:  1.4706E-01  1.0608E-01  7.5758E-01  1.8525E-01  3.4405E-01  2.8101E-01  1.0365E+00  8.3635E-01  1.0286E-01  1.2138E-01
             2.0446E-02
 GRADIENT:  -3.3449E+00  4.8950E+00  3.2361E+00  5.6214E-01 -7.5070E+00  8.8556E+00  1.1407E+01 -1.2024E-01  1.1148E+01  2.0714E+00
            -1.2362E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1632.39879273206        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  1.0558E+00  1.0583E+00  1.6394E+00  1.0210E+00  1.2658E+00  1.1784E+00  2.4487E+00  2.0267E+00  7.5953E-01  9.9571E-01
             9.2226E-01
 PARAMETER:  1.5426E-01  1.5666E-01  5.9431E-01  1.2079E-01  3.3567E-01  2.6413E-01  9.9556E-01  8.0639E-01 -1.7505E-01  9.5698E-02
             1.9067E-02
 GRADIENT:   7.2167E+00 -1.5604E+00 -1.9973E+00  1.0340E+00  5.4478E+00  2.1420E+00  2.1469E+00 -1.9555E-01  2.5348E-01 -2.9556E-01
            -2.2161E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1632.44098312514        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  1.0544E+00  1.1349E+00  1.5801E+00  9.7717E-01  1.2690E+00  1.1742E+00  2.3184E+00  2.0481E+00  7.8199E-01  9.9240E-01
             9.2714E-01
 PARAMETER:  1.5294E-01  2.2652E-01  5.5751E-01  7.6908E-02  3.3821E-01  2.6059E-01  9.4090E-01  8.1690E-01 -1.4592E-01  9.2371E-02
             2.4354E-02
 GRADIENT:   4.0869E+00 -8.0878E-01 -2.9700E-01 -1.8093E+00 -3.7751E-01  7.6205E-01  1.6985E+00  4.2736E-02  3.9620E-01  3.2497E-01
             1.5575E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1632.44529187914        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  1.0548E+00  1.1364E+00  1.5881E+00  9.7714E-01  1.2705E+00  1.1748E+00  2.3177E+00  2.0658E+00  7.7701E-01  9.8986E-01
             9.2719E-01
 PARAMETER:  1.5333E-01  2.2789E-01  5.6256E-01  7.6877E-02  3.3942E-01  2.6106E-01  9.4059E-01  8.2552E-01 -1.5231E-01  8.9812E-02
             2.4398E-02
 GRADIENT:   4.7252E+00 -3.2886E-01 -2.1323E-01 -1.0387E+00 -2.7447E-01  9.5608E-01  1.4894E+00  8.5188E-02  2.2856E-01 -4.8232E-02
             6.5079E-02

0ITERATION NO.:   34    OBJECTIVE VALUE:  -1632.44684469507        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1147
 NPARAMETR:  1.0548E+00  1.1356E+00  1.5919E+00  9.7771E-01  1.2709E+00  1.1748E+00  2.3201E+00  2.0668E+00  7.7128E-01  9.9116E-01
             9.2708E-01
 PARAMETER:  1.5336E-01  2.2715E-01  5.6495E-01  7.7458E-02  3.3976E-01  2.6110E-01  9.4163E-01  8.2599E-01 -1.5971E-01  9.1118E-02
             2.4285E-02
 GRADIENT:   8.5451E-02  3.7560E-02  6.9134E-02 -4.6906E-01 -4.8603E-01  1.4532E-02 -1.9644E-01 -1.7444E-01 -1.4393E-02 -1.2540E-02
            -6.2680E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1147
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.6703E-04  1.1704E-02 -5.3744E-02 -2.6510E-02 -4.4472E-02
 SE:             2.9898E-02  2.6175E-02  1.7306E-02  1.5726E-02  1.9354E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7953E-01  6.5478E-01  1.8992E-03  9.1842E-02  2.1573E-02

 ETASHRINKSD(%)  1.0000E-10  1.2309E+01  4.2024E+01  4.7317E+01  3.5161E+01
 ETASHRINKVR(%)  1.0000E-10  2.3103E+01  6.6388E+01  7.2245E+01  5.7960E+01
 EBVSHRINKSD(%)  2.8788E-01  1.0233E+01  4.8574E+01  5.1296E+01  3.0431E+01
 EBVSHRINKVR(%)  5.7493E-01  1.9419E+01  7.3554E+01  7.6279E+01  5.1602E+01
 RELATIVEINF(%)  9.9280E+01  1.7011E+01  6.6747E+00  4.0366E+00  1.7491E+01
 EPSSHRINKSD(%)  4.5791E+01
 EPSSHRINKVR(%)  7.0614E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1632.4468446950737     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -897.29601813133547     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1632.447       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.14E+00  1.59E+00  9.78E-01  1.27E+00  1.17E+00  2.32E+00  2.07E+00  7.71E-01  9.91E-01  9.27E-01
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       41.021
Stop Time:
Sun Oct 24 04:43:02 CDT 2021
