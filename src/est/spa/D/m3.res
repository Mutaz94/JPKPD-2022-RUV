Wed Sep 29 19:35:01 CDT 2021
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
$DATA ../../../../data/spa/D/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1633.43378498907        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9906E+02 -4.5644E+01 -4.8602E+01  1.6100E+01  5.6888E+01  1.3708E+01 -2.5723E+01  6.0679E+00 -1.7606E+01  1.8109E+01
            -5.0739E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1651.70423767562        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0170E+00  1.1770E+00  1.1970E+00  9.3314E-01  1.1120E+00  1.2510E+00  1.3284E+00  9.8114E-01  1.1507E+00  8.3023E-01
             1.0317E+00
 PARAMETER:  1.1684E-01  2.6298E-01  2.7979E-01  3.0797E-02  2.0612E-01  3.2396E-01  3.8400E-01  8.0964E-02  2.4034E-01 -8.6052E-02
             1.3120E-01
 GRADIENT:  -9.1996E+00 -9.7151E+00  1.3081E+01 -9.9197E+00  1.8242E+01  2.6056E+01  3.7235E+00 -8.3636E+00  1.0521E+01 -7.0619E+00
            -1.4021E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1652.93180434836        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0219E+00  1.1915E+00  1.0869E+00  9.2207E-01  1.0606E+00  1.2310E+00  1.3568E+00  1.1319E+00  1.1804E+00  7.3697E-01
             1.0213E+00
 PARAMETER:  1.2166E-01  2.7525E-01  1.8332E-01  1.8868E-02  1.5879E-01  3.0782E-01  4.0515E-01  2.2393E-01  2.6586E-01 -2.0520E-01
             1.2110E-01
 GRADIENT:  -3.0368E+00 -7.1810E+00  5.1434E+00 -4.6915E+00  4.7795E+00  1.9605E+01  7.6404E+00 -1.3421E+00  1.7442E+01 -5.4803E+00
            -1.8522E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1655.08916204979        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0286E+00  1.3311E+00  8.6608E-01  8.4464E-01  1.0300E+00  1.1732E+00  1.2312E+00  9.2525E-01  1.1020E+00  7.4849E-01
             1.0229E+00
 PARAMETER:  1.2815E-01  3.8601E-01 -4.3780E-02 -6.8848E-02  1.2957E-01  2.5972E-01  3.0802E-01  2.2314E-02  1.9715E-01 -1.8969E-01
             1.2266E-01
 GRADIENT:   4.3630E+00  1.0936E+01  2.6575E+00  9.3674E+00 -4.6279E+00 -7.3672E-02  4.9033E-01 -7.2927E-02 -3.0803E-02 -1.7760E-01
             3.1470E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1655.61139712421        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.0262E+00  1.6416E+00  5.7547E-01  6.3624E-01  1.0534E+00  1.1750E+00  1.0419E+00  5.6427E-01  1.2928E+00  7.3486E-01
             1.0212E+00
 PARAMETER:  1.2588E-01  5.9566E-01 -4.5257E-01 -3.5218E-01  1.5206E-01  2.6125E-01  1.4101E-01 -4.7222E-01  3.5684E-01 -2.0808E-01
             1.2096E-01
 GRADIENT:  -1.8991E+00  1.3630E+01  9.5185E-01  7.9156E+00 -8.3030E+00  1.9714E-01 -2.1152E+00  5.4711E-01 -9.6567E-01  2.2085E-01
            -9.5958E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1655.85062658459        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      904             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0317E+00  1.7649E+00  4.9589E-01  5.2548E-01  1.1145E+00  1.1882E+00  9.8792E-01  1.7870E-01  1.4670E+00  7.7546E-01
             1.0249E+00
 PARAMETER:  1.3117E-01  6.6809E-01 -6.0141E-01 -5.4345E-01  2.0841E-01  2.7241E-01  8.7846E-02 -1.6220E+00  4.8324E-01 -1.5430E-01
             1.2461E-01
 GRADIENT:   5.8381E+02  6.9923E+02  4.5835E+00  1.2649E+02  1.5137E+01  2.3849E+02  1.2027E+01  1.1488E-01  1.9301E+01  1.0497E+00
             2.3690E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1656.02782989353        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1085             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0296E+00  1.7807E+00  4.9500E-01  5.3087E-01  1.1150E+00  1.1849E+00  9.8220E-01  5.8253E-02  1.4720E+00  7.7643E-01
             1.0226E+00
 PARAMETER:  1.2917E-01  6.7701E-01 -6.0320E-01 -5.3323E-01  2.0887E-01  2.6969E-01  8.2040E-02 -2.7430E+00  4.8661E-01 -1.5304E-01
             1.2232E-01
 GRADIENT:   5.7395E+02  7.3175E+02  4.1974E+00  1.3439E+02  1.5546E+01  2.3714E+02  1.1030E+01  1.4449E-02  2.0157E+01  4.3063E-01
             9.4636E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1656.03004597970        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  1.0284E+00  1.7828E+00  4.9564E-01  5.3149E-01  1.1142E+00  1.1797E+00  9.8367E-01  2.0183E-02  1.4704E+00  7.7693E-01
             1.0226E+00
 PARAMETER:  1.2798E-01  6.7816E-01 -6.0190E-01 -5.3207E-01  2.0812E-01  2.6522E-01  8.3535E-02 -3.8029E+00  4.8552E-01 -1.5241E-01
             1.2233E-01
 GRADIENT:   1.5534E+00 -5.7602E+00 -1.2355E-01  2.2921E-02  1.9347E-01  1.9133E+00  1.4058E-01  9.5193E-04  1.3176E-01 -7.0054E-02
            -6.8396E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1656.03575734656        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1419
 NPARAMETR:  1.0295E+00  1.7761E+00  4.9788E-01  5.3135E-01  1.1125E+00  1.1823E+00  9.8400E-01  1.0000E-02  1.4660E+00  7.7749E-01
             1.0229E+00
 PARAMETER:  1.2865E-01  6.7647E-01 -6.0066E-01 -5.3033E-01  2.0765E-01  2.6665E-01  8.3580E-02 -4.8285E+00  4.8549E-01 -1.5322E-01
             1.2206E-01
 GRADIENT:  -3.5864E-01  1.9853E+00 -3.7659E-01  6.1286E-01  9.6461E-01 -1.7980E-01 -3.8439E-02  0.0000E+00  2.4743E-01 -4.1067E-02
            -1.2917E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1419
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1997E-04 -2.1906E-02 -2.7963E-04  2.4577E-02 -3.5158E-02
 SE:             2.9912E-02  2.6241E-02  1.0109E-04  2.2210E-02  1.9860E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9147E-01  4.0383E-01  5.6743E-03  2.6847E-01  7.6680E-02

 ETASHRINKSD(%)  1.0000E-10  1.2089E+01  9.9661E+01  2.5594E+01  3.3466E+01
 ETASHRINKVR(%)  1.0000E-10  2.2717E+01  9.9999E+01  4.4637E+01  5.5732E+01
 EBVSHRINKSD(%)  3.1640E-01  1.2196E+01  9.9703E+01  2.7039E+01  3.2739E+01
 EBVSHRINKVR(%)  6.3180E-01  2.2905E+01  9.9999E+01  4.6767E+01  5.4759E+01
 RELATIVEINF(%)  9.9341E+01  8.2699E+00  9.7011E-05  4.5020E+00  9.3496E+00
 EPSSHRINKSD(%)  4.4191E+01
 EPSSHRINKVR(%)  6.8854E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1656.0357573465606     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -920.88493078282238     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1656.036       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.78E+00  4.96E-01  5.32E-01  1.11E+00  1.18E+00  9.84E-01  1.00E-02  1.47E+00  7.76E-01  1.02E+00
 


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
+        7.47E+02
 
 TH 2
+       -2.93E+00  2.88E+02
 
 TH 3
+        6.00E+00  1.33E+02  4.39E+02
 
 TH 4
+       -9.97E+00  2.06E+02 -4.17E+02  1.01E+03
 
 TH 5
+       -4.63E+00 -1.68E+02 -4.07E+02  4.11E+02  6.85E+02
 
 TH 6
+        9.67E-02 -5.24E-01  1.23E+00 -2.01E+00 -6.24E-01  1.41E+02
 
 TH 7
+       -8.78E-02  1.27E+01 -2.91E+01 -1.79E+01  7.45E+00 -6.45E-01  1.26E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.43E-01 -1.32E+01 -3.67E+01  6.34E+01 -7.34E+00 -2.63E-01  8.71E+00  0.00E+00  3.66E+01
 
 TH10
+        2.82E-01 -1.27E+01 -2.91E+01 -1.20E+01 -8.45E+01  1.77E-01  1.45E+01  0.00E+00  9.41E+00  8.32E+01
 
 TH11
+       -4.56E+00 -1.36E+01 -2.20E+01  2.11E+00 -1.85E+01  2.06E+00  9.14E+00  0.00E+00  5.88E+00  2.29E+01  2.01E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.359
Stop Time:
Wed Sep 29 19:35:29 CDT 2021
