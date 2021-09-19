Sat Sep 18 08:23:36 CDT 2021
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
$DATA ../../../../data/spa/B/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.47321181186        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.7623E+01 -7.6209E+01 -1.4673E+01 -1.1049E+02 -4.0769E+01  3.0672E+01 -3.8873E+00  1.7335E+01  3.5232E+00  2.7517E+01
             6.9283E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1708.91178581975        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8867E-01  1.1210E+00  1.1969E+00  1.0103E+00  1.1829E+00  8.7561E-01  1.0297E+00  8.7697E-01  9.4971E-01  8.6083E-01
             1.0009E+00
 PARAMETER:  8.8603E-02  2.1420E-01  2.7971E-01  1.1024E-01  2.6798E-01 -3.2832E-02  1.2931E-01 -3.1282E-02  4.8399E-02 -4.9859E-02
             1.0092E-01
 GRADIENT:   3.0450E+01  3.9056E+00  2.1804E+00  7.8509E+00  3.7972E+01 -2.0413E+01 -3.4639E+00 -2.0402E+00 -7.4910E+00 -1.9950E+01
            -7.4436E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1710.93613277542        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8331E-01  1.0996E+00  9.2468E-01  1.0065E+00  1.0461E+00  9.2060E-01  1.0940E+00  3.9902E-01  8.9822E-01  8.7955E-01
             9.5546E-01
 PARAMETER:  8.3165E-02  1.9498E-01  2.1691E-02  1.0648E-01  1.4508E-01  1.7275E-02  1.8987E-01 -8.1875E-01 -7.3370E-03 -2.8345E-02
             5.4440E-02
 GRADIENT:   1.7397E+01 -7.1257E+00 -1.5176E+01  9.0367E+00  3.2978E+01 -8.6477E-02 -2.9133E+00  4.4901E-01 -1.2698E+01 -1.5636E+00
            -1.8538E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1711.90865311751        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.7812E-01  1.0982E+00  8.9780E-01  1.0038E+00  1.0080E+00  9.2378E-01  1.0916E+00  3.5416E-01  9.5222E-01  8.4733E-01
             9.9127E-01
 PARAMETER:  7.7876E-02  1.9366E-01 -7.8102E-03  1.0380E-01  1.0793E-01  2.0715E-02  1.8765E-01 -9.3800E-01  5.1046E-02 -6.5664E-02
             9.1232E-02
 GRADIENT:   3.1649E-02 -1.1931E+00 -2.3333E+00  5.5533E-01  9.9137E-01  7.7136E-01 -6.3730E-01  5.9880E-01 -3.8433E-01  7.2577E-01
            -8.5006E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1711.98905900306        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.7800E-01  1.0989E+00  8.5259E-01  9.9913E-01  9.8017E-01  9.2178E-01  1.1119E+00  1.8897E-01  9.5608E-01  8.1749E-01
             9.9877E-01
 PARAMETER:  7.7751E-02  1.9431E-01 -5.9474E-02  9.9128E-02  7.9976E-02  1.8553E-02  2.0606E-01 -1.5662E+00  5.5085E-02 -1.0151E-01
             9.8771E-02
 GRADIENT:  -2.1770E+00 -6.2710E-01 -1.3411E+00 -4.4929E-01 -1.7652E+00 -2.9422E-01  9.9006E-01  2.0373E-01  1.6698E+00  5.5553E-01
             2.9633E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.80386367089        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.9716E-01  1.2375E+00  8.2377E-01  9.2534E-01  1.0339E+00  9.3083E-01  1.0178E+00  2.9379E-02  1.0179E+00  8.4779E-01
             9.9574E-01
 PARAMETER:  9.7156E-02  3.1308E-01 -9.3864E-02  2.2406E-02  1.3332E-01  2.8326E-02  1.1768E-01 -3.4275E+00  1.1776E-01 -6.5117E-02
             9.5729E-02
 GRADIENT:   6.5921E+00  7.5517E+00  1.8904E+00  8.4164E+00 -5.9449E+00  4.1862E-01  1.1734E-01  3.8533E-03  2.4291E-01  2.4331E-01
             8.0196E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1713.25383819520        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  9.9498E-01  1.5201E+00  7.4894E-01  7.3943E-01  1.1708E+00  9.2923E-01  8.6817E-01  1.0000E-02  1.2066E+00  9.2071E-01
             9.9685E-01
 PARAMETER:  9.4970E-02  5.1879E-01 -1.8909E-01 -2.0187E-01  2.5770E-01  2.6604E-02 -4.1364E-02 -6.4039E+00  2.8780E-01  1.7389E-02
             9.6840E-02
 GRADIENT:  -5.5937E-01  2.0334E+00  5.7319E-01  6.4141E-01 -2.6309E+00 -4.5032E-01  7.4992E-02  0.0000E+00  4.4803E-01 -9.5256E-02
             4.5665E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1713.26651880221        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      795
 NPARAMETR:  9.9525E-01  1.5617E+00  7.3842E-01  7.1177E-01  1.1966E+00  9.3036E-01  8.5118E-01  1.0000E-02  1.2386E+00  9.3893E-01
             9.9631E-01
 PARAMETER:  9.5239E-02  5.4580E-01 -2.0324E-01 -2.4000E-01  2.7948E-01  2.7813E-02 -6.1135E-02 -6.9309E+00  3.1396E-01  3.6985E-02
             9.6303E-02
 GRADIENT:   2.5587E-02  5.3087E-02  4.1847E-02  2.0084E-03 -3.5723E-02  1.7116E-02 -1.1357E-03  0.0000E+00 -7.8777E-03 -3.7102E-03
            -1.6011E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1713.26651880221        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  9.9524E-01  1.5616E+00  7.3827E-01  7.1176E-01  1.1966E+00  9.3029E-01  8.5119E-01  1.0000E-02  1.2387E+00  9.3897E-01
             9.9638E-01
 PARAMETER:  9.5239E-02  5.4580E-01 -2.0324E-01 -2.4000E-01  2.7948E-01  2.7813E-02 -6.1135E-02 -6.9309E+00  3.1396E-01  3.6985E-02
             9.6303E-02
 GRADIENT:   1.0121E-02  8.3115E-02  2.1461E-02  3.6405E-03 -3.5020E-02  1.4714E-02 -1.6430E-03  0.0000E+00 -8.4842E-03 -3.7574E-03
            -1.6414E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      823
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.8170E-05 -2.3169E-02 -3.2604E-04  1.6565E-02 -3.4228E-02
 SE:             2.9814E-02  2.2653E-02  1.2681E-04  2.3229E-02  2.1851E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9764E-01  3.0641E-01  1.0141E-02  4.7577E-01  1.1724E-01

 ETASHRINKSD(%)  1.1823E-01  2.4110E+01  9.9575E+01  2.2179E+01  2.6797E+01
 ETASHRINKVR(%)  2.3632E-01  4.2407E+01  9.9998E+01  3.9438E+01  4.6413E+01
 EBVSHRINKSD(%)  4.8761E-01  2.3203E+01  9.9610E+01  2.3362E+01  2.5661E+01
 EBVSHRINKVR(%)  9.7285E-01  4.1022E+01  9.9998E+01  4.1266E+01  4.4737E+01
 RELATIVEINF(%)  9.8875E+01  3.3668E+00  2.1167E-04  3.7296E+00  9.4866E+00
 EPSSHRINKSD(%)  4.2860E+01
 EPSSHRINKVR(%)  6.7350E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1713.2665188022081     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -978.11569223846993     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1713.267       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  1.56E+00  7.38E-01  7.12E-01  1.20E+00  9.30E-01  8.51E-01  1.00E-02  1.24E+00  9.39E-01  9.96E-01
 


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
+        1.28E+03
 
 TH 2
+       -6.10E+00  3.89E+02
 
 TH 3
+        1.32E+01  1.58E+02  3.53E+02
 
 TH 4
+       -1.63E+01  3.24E+02 -2.20E+02  8.49E+02
 
 TH 5
+       -3.80E+00 -1.90E+02 -2.97E+02  2.15E+02  4.88E+02
 
 TH 6
+       -3.60E+00 -8.04E-01  2.02E+00 -3.34E+00 -1.11E+00  2.27E+02
 
 TH 7
+        4.83E-01  1.71E+01 -7.44E+00 -1.41E+01 -9.39E+00 -4.86E-01  1.10E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.75E+00 -2.13E+01 -2.95E+01  4.76E+01 -1.06E+00  1.17E-01  1.69E+01  0.00E+00  5.77E+01
 
 TH10
+        1.80E+00 -1.07E+01 -3.26E+01 -8.01E+00 -5.97E+01 -1.41E+00  1.19E+01  0.00E+00  5.99E+00  8.16E+01
 
 TH11
+       -7.37E+00 -1.96E+01 -3.11E+01  1.68E+00  1.06E+00 -1.32E+00  9.21E+00  0.00E+00  4.97E+00  1.58E+01  2.24E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.456
Stop Time:
Sat Sep 18 08:23:52 CDT 2021
