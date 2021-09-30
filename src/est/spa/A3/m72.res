Wed Sep 29 13:45:45 CDT 2021
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
$DATA ../../../../data/spa/A3/dat72.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   82.6175035893417        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8617E+02  1.1950E+02  1.9640E+02 -5.9091E+01  1.0020E+02  5.8583E+01 -6.8132E+01 -7.4114E+01 -1.6285E+02 -1.3799E+02
            -3.0472E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1028.39383995966        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0479E+00  8.5697E-01  7.6842E-01  1.1008E+00  7.7101E-01  9.4019E-01  1.0989E+00  9.5542E-01  1.3244E+00  1.1818E+00
             1.8192E+00
 PARAMETER:  1.4681E-01 -5.4356E-02 -1.6342E-01  1.9603E-01 -1.6005E-01  3.8329E-02  1.9429E-01  5.4391E-02  3.8097E-01  2.6700E-01
             6.9839E-01
 GRADIENT:   2.3553E+02  4.9171E+01  6.9493E+01  8.9498E+00 -1.5059E+01  1.2869E+01 -1.0241E+01 -3.7073E-01 -2.2574E-01 -8.2729E-01
            -7.6112E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1075.52859710150        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      201
 NPARAMETR:  1.0580E+00  7.3531E-01  3.1135E-01  1.1304E+00  4.0016E-01  8.1892E-01  1.2626E+00  2.8650E-01  1.2231E+00  7.5006E-01
             1.8051E+00
 PARAMETER:  1.5635E-01 -2.0746E-01 -1.0669E+00  2.2261E-01 -8.1588E-01 -9.9768E-02  3.3317E-01 -1.1500E+00  3.0142E-01 -1.8760E-01
             6.9060E-01
 GRADIENT:   8.5320E+01  1.4119E+02  1.1816E+02  3.8399E+01 -1.0878E+02 -6.0327E+01 -6.8736E+00 -1.6876E+00 -3.0556E+01  1.7673E+01
            -6.9573E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1251.64167625071        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0455E+00  5.0157E-01  1.4504E-01  1.2184E+00  2.4865E-01  9.5333E-01  6.7342E-01  1.0000E-02  1.9619E+00  5.8355E-01
             2.4313E+00
 PARAMETER:  1.4446E-01 -5.9000E-01 -1.8308E+00  2.9755E-01 -1.2917E+00  5.2206E-02 -2.9539E-01 -6.2483E+00  7.7390E-01 -4.3863E-01
             9.8842E-01
 GRADIENT:   5.4945E+01  6.2134E+01  3.2677E+01  1.0618E+02 -5.5541E+01 -6.4447E+00 -3.9495E+00  0.0000E+00  4.2810E+01  2.2560E+00
            -1.9926E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1298.39576048206        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0340E+00  4.4689E-01  1.4377E-01  9.8785E-01  2.4266E-01  9.5848E-01  8.5347E-01  1.0000E-02  1.4855E+00  4.2546E-01
             3.1661E+00
 PARAMETER:  1.3342E-01 -7.0544E-01 -1.8396E+00  8.7778E-02 -1.3161E+00  5.7593E-02 -5.8444E-02 -6.5405E+00  4.9574E-01 -7.5459E-01
             1.2525E+00
 GRADIENT:  -7.5907E-01 -1.3637E+01 -8.4349E+00 -1.8808E+01  2.3691E+01  1.8104E+00  2.4107E-02  0.0000E+00 -1.2631E+00  2.1119E+00
            -6.0680E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1299.38117564913        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  1.0386E+00  4.6774E-01  1.6686E-01  1.0488E+00  2.6352E-01  9.4630E-01  1.0402E+00  1.0039E-02  1.3991E+00  2.9905E-01
             3.2190E+00
 PARAMETER:  1.3789E-01 -6.5985E-01 -1.6906E+00  1.4764E-01 -1.2336E+00  4.4801E-02  1.3943E-01 -4.5012E+00  4.3581E-01 -1.1072E+00
             1.2691E+00
 GRADIENT:  -2.6902E+00 -6.4122E+00 -2.7229E+00  3.9509E+00  1.2159E+01 -1.5234E-01 -4.4008E-01 -8.4431E-04 -1.0821E+00  8.8873E-01
            -4.0938E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1299.63049347198        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  1.0374E+00  4.9487E-01  1.6001E-01  1.0277E+00  2.6527E-01  9.4400E-01  1.0131E+00  1.2458E-02  1.4200E+00  2.5375E-01
             3.2345E+00
 PARAMETER:  1.3671E-01 -6.0346E-01 -1.7325E+00  1.2733E-01 -1.2270E+00  4.2371E-02  1.1300E-01 -4.2854E+00  4.5065E-01 -1.2714E+00
             1.2739E+00
 GRADIENT:   7.3818E-01  2.6482E-01 -3.9006E-01  1.1777E+00  2.3604E-01 -1.5807E+00 -2.2351E-01 -3.1626E-03 -8.2236E-01  3.7384E-01
            -9.0529E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1299.67619555883        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1078
 NPARAMETR:  1.0366E+00  5.0254E-01  1.6033E-01  1.0246E+00  2.6791E-01  9.4771E-01  1.0406E+00  1.9427E-02  1.4207E+00  1.9836E-01
             3.2441E+00
 PARAMETER:  1.3594E-01 -5.8808E-01 -1.7305E+00  1.2435E-01 -1.2171E+00  4.6298E-02  1.3981E-01 -3.8411E+00  4.5114E-01 -1.5177E+00
             1.2768E+00
 GRADIENT:   1.5077E-02 -1.5412E-03  6.0493E-02  2.4234E-02  1.1443E-01  4.2213E-03  3.3383E-02 -8.2268E-03  1.2932E-04  3.9320E-03
             9.4597E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1299.67619555883        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1100
 NPARAMETR:  1.0366E+00  5.0254E-01  1.6033E-01  1.0246E+00  2.6791E-01  9.4771E-01  1.0406E+00  1.9427E-02  1.4207E+00  1.9836E-01
             3.2441E+00
 PARAMETER:  1.3594E-01 -5.8808E-01 -1.7305E+00  1.2435E-01 -1.2171E+00  4.6298E-02  1.3981E-01 -3.8411E+00  4.5114E-01 -1.5177E+00
             1.2768E+00
 GRADIENT:   1.5077E-02 -1.5412E-03  6.0493E-02  2.4234E-02  1.1443E-01  4.2213E-03  3.3383E-02 -8.2268E-03  1.2932E-04  3.9320E-03
             9.4597E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1100
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1211E-04  9.9317E-03 -1.4644E-04 -1.5139E-02  8.5676E-03
 SE:             2.8672E-02  1.9048E-02  3.6912E-04  2.6303E-02  7.2961E-03
 N:                     100         100         100         100         100

 P VAL.:         9.9131E-01  6.0209E-01  6.9158E-01  5.6491E-01  2.4029E-01

 ETASHRINKSD(%)  3.9446E+00  3.6187E+01  9.8763E+01  1.1881E+01  7.5557E+01
 ETASHRINKVR(%)  7.7336E+00  5.9279E+01  9.9985E+01  2.2351E+01  9.4025E+01
 EBVSHRINKSD(%)  3.6551E+00  3.5963E+01  9.8918E+01  9.9940E+00  7.6284E+01
 EBVSHRINKVR(%)  7.1766E+00  5.8993E+01  9.9988E+01  1.8989E+01  9.4375E+01
 RELATIVEINF(%)  8.7437E+01  2.6426E+00  1.5776E-03  5.0801E+01  2.0523E-01
 EPSSHRINKSD(%)  2.9987E+01
 EPSSHRINKVR(%)  5.0982E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1299.6761955588277     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -564.52536899508948     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1299.676       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  5.03E-01  1.60E-01  1.02E+00  2.68E-01  9.48E-01  1.04E+00  1.94E-02  1.42E+00  1.98E-01  3.24E+00
 


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
+        1.07E+03
 
 TH 2
+       -7.48E+01  1.62E+03
 
 TH 3
+       -5.14E+02  3.60E+03  1.42E+04
 
 TH 4
+       -2.65E+01  8.50E+01 -7.60E+02  4.26E+02
 
 TH 5
+        4.95E+02 -5.76E+03 -1.60E+04  8.97E+01  2.45E+04
 
 TH 6
+       -1.68E-01 -3.65E+00  5.15E+01 -1.21E+01  5.46E+01  1.89E+02
 
 TH 7
+       -3.23E+00  1.43E+01 -1.39E+02 -5.17E-01  1.25E+02  2.89E+00  2.69E+01
 
 TH 8
+       -1.05E-01 -7.08E-01 -3.83E+00 -1.31E-01  6.78E+00  4.30E-02  2.93E-01 -2.20E+01
 
 TH 9
+        8.86E+00 -2.04E+01  1.45E+02 -8.16E+00  9.39E+01 -4.44E-01  3.92E+00  1.09E-01  6.45E+01
 
 TH10
+       -1.57E+00 -2.11E+01 -2.15E+02 -3.46E+00  3.16E+02  2.08E+00  1.77E+01  6.27E-01  8.94E-01  2.96E+01
 
 TH11
+       -1.70E+01 -3.25E+00 -8.84E+01 -4.73E+00  6.65E+01  2.92E+00  1.08E+01  3.90E-01  3.39E+00  9.91E+00  3.03E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.732
Stop Time:
Wed Sep 29 13:46:08 CDT 2021
