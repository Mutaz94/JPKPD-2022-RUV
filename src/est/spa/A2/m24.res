Sat Sep 18 09:44:27 CDT 2021
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
$DATA ../../../../data/spa/A2/dat24.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1098.57649670645        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1364E+02 -2.1916E+01  1.0193E+02 -1.5316E+02  2.7912E+01  1.1386E+01 -3.6310E+01 -3.5015E+01 -6.4272E+01 -8.8771E+01
            -8.5922E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1387.38161546003        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.5328E-01  9.6525E-01  8.9240E-01  1.1406E+00  9.4409E-01  9.2844E-01  1.0482E+00  9.9124E-01  1.0618E+00  1.1660E+00
             2.1490E+00
 PARAMETER:  5.2156E-02  6.4632E-02 -1.3840E-02  2.3154E-01  4.2461E-02  2.5746E-02  1.4710E-01  9.1199E-02  1.5996E-01  2.5360E-01
             8.6501E-01
 GRADIENT:  -5.5726E+01  1.1198E+01 -8.1045E+00  2.1895E+01  1.5499E+01 -1.0240E+01 -2.9429E+00  5.1219E+00 -1.2284E+00 -6.3147E+00
            -1.9450E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1392.32252035440        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.6348E-01  8.3987E-01  6.0151E-01  1.2006E+00  6.9671E-01  9.4068E-01  1.4762E+00  2.5218E-01  9.3717E-01  1.0139E+00
             2.1725E+00
 PARAMETER:  6.2802E-02 -7.4507E-02 -4.0832E-01  2.8280E-01 -2.6139E-01  3.8851E-02  4.8949E-01 -1.2776E+00  3.5110E-02  1.1379E-01
             8.7587E-01
 GRADIENT:  -3.8105E+01  1.7057E+01 -3.5819E+01  5.6354E+01  4.0706E+01 -6.4144E+00  8.1349E+00  8.1850E-01  4.6923E-01  7.8098E+00
             1.3151E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1399.85666379513        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.9275E-01  6.1002E-01  3.5280E-01  1.2301E+00  4.0438E-01  9.5869E-01  1.3263E+00  8.0371E-02  8.9297E-01  6.2080E-01
             2.1902E+00
 PARAMETER:  9.2723E-02 -3.9427E-01 -9.4185E-01  3.0709E-01 -8.0541E-01  5.7814E-02  3.8241E-01 -2.4211E+00 -1.3203E-02 -3.7674E-01
             8.8400E-01
 GRADIENT:   2.6495E+01  4.3695E+01  3.6394E+01  5.7818E+01 -5.7102E+01 -1.3734E+00 -7.3745E+00  1.6639E-02 -5.2243E+00 -4.8419E+00
             6.2084E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1403.30578627010        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  9.8429E-01  5.1316E-01  2.5591E-01  1.1540E+00  3.2334E-01  9.7606E-01  1.3985E+00  1.7784E-02  9.3257E-01  5.3903E-01
             2.1041E+00
 PARAMETER:  8.4167E-02 -5.6716E-01 -1.2629E+00  2.4326E-01 -1.0290E+00  7.5773E-02  4.3541E-01 -3.9295E+00  3.0187E-02 -5.1799E-01
             8.4387E-01
 GRADIENT:   1.0122E+01  7.3848E+00  1.2292E+01 -8.8085E+00 -2.6937E+01  2.5143E+00 -1.4677E+00 -1.3215E-03 -9.8795E-01 -1.5040E+00
             1.3044E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1404.08565978576        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      502
 NPARAMETR:  9.8110E-01  4.4466E-01  3.3830E-01  1.2584E+00  3.6373E-01  9.5811E-01  1.7089E+00  2.7889E-02  8.8774E-01  6.3732E-01
             2.1200E+00
 PARAMETER:  8.0920E-02 -7.1045E-01 -9.8381E-01  3.2983E-01 -9.1133E-01  5.7212E-02  6.3583E-01 -3.4795E+00 -1.9073E-02 -3.5049E-01
             8.5140E-01
 GRADIENT:  -2.7124E+00  3.8943E+00  2.0861E+00  7.4022E+00 -5.1945E+00 -1.6602E+00 -1.2009E+00  2.0148E-03 -8.5615E-01 -1.1856E+00
            -9.6843E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1405.81938503378        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      679
 NPARAMETR:  9.7763E-01  2.2045E-01  4.2616E-01  1.4064E+00  3.7889E-01  9.5447E-01  2.9837E+00  1.0000E-02  8.4916E-01  7.8386E-01
             2.1506E+00
 PARAMETER:  7.7380E-02 -1.4121E+00 -7.5294E-01  4.4100E-01 -8.7052E-01  5.3401E-02  1.1932E+00 -5.4650E+00 -6.3503E-02 -1.4352E-01
             8.6574E-01
 GRADIENT:   5.0556E+00  1.3437E+00  2.3749E+01  9.9470E+00 -3.5228E+01  1.5133E+00 -3.5146E+00  0.0000E+00  1.7733E+00  4.2712E+00
             3.2264E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1406.42354252657        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      854
 NPARAMETR:  9.7291E-01  1.8769E-01  4.3756E-01  1.4192E+00  3.8974E-01  9.4973E-01  3.4698E+00  1.0000E-02  8.2984E-01  7.7201E-01
             2.1411E+00
 PARAMETER:  7.2541E-02 -1.5729E+00 -7.2654E-01  4.5008E-01 -8.4228E-01  4.8419E-02  1.3441E+00 -6.1239E+00 -8.6522E-02 -1.5876E-01
             8.6134E-01
 GRADIENT:  -5.3190E-02 -1.5998E-02 -7.0414E-01  1.2899E-01  1.0559E+00  4.0436E-01  7.0935E-02  0.0000E+00  2.8305E-01  1.0525E-01
            -1.3738E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1406.42457611917        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      946
 NPARAMETR:  9.7299E-01  1.8770E-01  4.3610E-01  1.4183E+00  3.8856E-01  9.4870E-01  3.4640E+00  1.0000E-02  8.2916E-01  7.7009E-01
             2.1422E+00
 PARAMETER:  7.2621E-02 -1.5729E+00 -7.2988E-01  4.4943E-01 -8.4530E-01  4.7335E-02  1.3424E+00 -6.1302E+00 -8.7337E-02 -1.6125E-01
             8.6182E-01
 GRADIENT:  -2.1455E-02 -5.5173E-03  2.3730E-02  2.3768E-02 -2.4593E-02 -5.4510E-04 -8.6967E-03  0.0000E+00  7.3978E-03 -1.1226E-02
            -2.0001E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      946
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7259E-03  4.0394E-02 -2.5704E-04 -2.4797E-02  1.1045E-02
 SE:             2.9278E-02  1.5667E-02  2.2088E-04  2.5896E-02  2.0867E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5299E-01  9.9310E-03  2.4456E-01  3.3827E-01  5.9657E-01

 ETASHRINKSD(%)  1.9161E+00  4.7513E+01  9.9260E+01  1.3246E+01  3.0094E+01
 ETASHRINKVR(%)  3.7955E+00  7.2451E+01  9.9995E+01  2.4738E+01  5.1132E+01
 EBVSHRINKSD(%)  1.9901E+00  5.8017E+01  9.9214E+01  1.1542E+01  2.4466E+01
 EBVSHRINKVR(%)  3.9407E+00  8.2374E+01  9.9994E+01  2.1751E+01  4.2946E+01
 RELATIVEINF(%)  9.4735E+01  7.1793E+00  2.5060E-04  2.9741E+01  2.2536E+00
 EPSSHRINKSD(%)  3.7039E+01
 EPSSHRINKVR(%)  6.0359E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1406.4245761191717     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -671.27374955543348     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.09
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1406.425       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.73E-01  1.88E-01  4.36E-01  1.42E+00  3.89E-01  9.49E-01  3.46E+00  1.00E-02  8.29E-01  7.70E-01  2.14E+00
 


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
+        1.25E+03
 
 TH 2
+       -1.34E+01  1.88E+03
 
 TH 3
+       -1.14E+01 -8.22E+01  3.96E+03
 
 TH 4
+       -3.24E+01  3.99E+01 -3.06E+02  6.69E+02
 
 TH 5
+        8.50E+01 -2.74E+02 -5.18E+03 -6.23E+01  7.55E+03
 
 TH 6
+        9.00E-02  1.69E+01 -1.94E+00 -1.22E+01  1.08E+01  1.97E+02
 
 TH 7
+        4.20E+00  1.10E+02 -3.19E+01 -1.35E+01  3.37E+01  1.51E+00  8.09E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.93E+00 -1.76E+02  4.84E+01  6.01E+00  7.63E+00 -1.31E+01 -5.68E+00  0.00E+00  2.03E+02
 
 TH10
+       -7.58E+00 -8.50E+01 -5.68E+01  7.50E+00 -3.72E+01 -4.10E+00 -7.13E+00  0.00E+00 -3.38E+00  1.31E+02
 
 TH11
+       -1.51E+01 -2.47E+01 -2.07E+01 -7.22E+00  1.42E+01  2.79E+00 -1.04E+00  0.00E+00  1.49E+01  2.20E+01  5.52E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.153
Stop Time:
Sat Sep 18 09:44:46 CDT 2021
