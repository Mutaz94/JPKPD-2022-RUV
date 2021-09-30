Wed Sep 29 12:17:13 CDT 2021
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
$DATA ../../../../data/spa/A1/dat67.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1271.78280702962        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3157E+02 -4.0933E+01  2.7479E+00 -6.1906E+01  5.8952E+01  4.5250E+01 -2.7039E+01 -6.1213E+00 -5.6361E+01 -1.0216E+01
            -6.7433E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1457.33803700277        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0535E+00  1.1874E+00  1.2642E+00  1.0024E+00  1.1722E+00  1.1056E+00  9.7732E-01  8.9737E-01  8.7718E-01  6.9607E-01
             2.7860E+00
 PARAMETER:  1.5208E-01  2.7175E-01  3.3446E-01  1.0238E-01  2.5891E-01  2.0039E-01  7.7064E-02 -8.2908E-03 -3.1049E-02 -2.6230E-01
             1.1246E+00
 GRADIENT:   1.5710E+02 -2.9130E+00 -2.6371E+00 -1.5638E+01 -7.3886E+00  4.4564E+01 -1.2885E+00  2.3047E+00 -6.8071E+00  8.4529E+00
             1.0693E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1467.50112101282        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.0270E+00  1.0084E+00  1.0493E+00  1.1202E+00  9.6437E-01  9.8229E-01  1.2580E+00  6.3832E-01  9.3971E-01  3.4781E-01
             2.5875E+00
 PARAMETER:  1.2668E-01  1.0839E-01  1.4809E-01  2.1353E-01  6.3721E-02  8.2131E-02  3.2953E-01 -3.4891E-01  3.7818E-02 -9.5610E-01
             1.0507E+00
 GRADIENT:   1.2242E+02  2.2697E+01  1.2977E+01  2.5824E+01 -3.5487E+01  2.4895E+00  1.0567E+01  1.3301E+00  7.9949E+00  1.6903E+00
             6.8188E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1471.92479801652        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      268
 NPARAMETR:  9.9699E-01  1.0626E+00  1.1171E+00  1.0785E+00  1.0515E+00  9.8674E-01  1.0935E+00  5.5767E-01  9.7217E-01  4.5537E-01
             2.3868E+00
 PARAMETER:  9.6989E-02  1.6070E-01  2.1071E-01  1.7560E-01  1.5026E-01  8.6647E-02  1.8942E-01 -4.8398E-01  7.1771E-02 -6.8664E-01
             9.6996E-01
 GRADIENT:  -1.3084E+01 -1.3726E+00 -8.8920E-01 -2.2146E+00 -1.0464E+00 -2.8626E+00  5.6575E-01  3.9635E-01  1.5219E+00  4.9111E-01
             7.9243E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1472.19198167540        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.0002E+00  9.6744E-01  1.2810E+00  1.1481E+00  1.0810E+00  9.9182E-01  1.0824E+00  1.3667E-01  9.5669E-01  6.0809E-01
             2.3549E+00
 PARAMETER:  1.0018E-01  6.6902E-02  3.4762E-01  2.3811E-01  1.7793E-01  9.1790E-02  1.7915E-01 -1.8902E+00  5.5725E-02 -3.9744E-01
             9.5652E-01
 GRADIENT:  -2.6325E+00  1.5844E+00  6.9620E-01  1.1595E+00 -1.4487E+00 -8.0723E-01 -5.1146E-03  2.2713E-02  1.9670E-01  1.5254E-01
            -1.1981E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1472.24494483359        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      618
 NPARAMETR:  9.9905E-01  7.9256E-01  1.4054E+00  1.2615E+00  1.0627E+00  9.9232E-01  1.1818E+00  3.0588E-02  8.9902E-01  5.6894E-01
             2.3858E+00
 PARAMETER:  9.9048E-02 -1.3249E-01  4.4035E-01  3.3233E-01  1.6083E-01  9.2287E-02  2.6701E-01 -3.3871E+00 -6.4481E-03 -4.6399E-01
             9.6953E-01
 GRADIENT:  -1.5194E+00  6.2348E-01 -2.1445E-02  8.6243E-01 -8.0903E-02  6.6138E-02 -1.7383E-01  1.0539E-03 -4.8153E-01 -1.3537E-01
             5.4073E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1472.25517888767        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  9.9853E-01  7.0959E-01  1.4581E+00  1.3167E+00  1.0507E+00  9.9116E-01  1.2484E+00  1.0072E-02  8.8039E-01  6.1756E-01
             2.3718E+00
 PARAMETER:  9.8530E-02 -2.4307E-01  4.7717E-01  3.7510E-01  1.4943E-01  9.1122E-02  3.2188E-01 -4.4979E+00 -2.7385E-02 -3.8198E-01
             9.6367E-01
 GRADIENT:  -7.1852E-01  1.6661E+00  8.1835E-01  2.9121E+00 -1.5595E+00 -1.8336E-01  2.1457E-01  1.1411E-04  2.1748E-01  1.1406E-02
             3.7591E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1472.26087012218        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  9.9854E-01  6.7906E-01  1.4657E+00  1.3341E+00  1.0435E+00  9.9138E-01  1.2660E+00  1.0000E-02  8.7320E-01  6.2726E-01
             2.3680E+00
 PARAMETER:  9.8536E-02 -2.8704E-01  4.8231E-01  3.8823E-01  1.4255E-01  9.1341E-02  3.3587E-01 -4.9147E+00 -3.5591E-02 -3.6640E-01
             9.6205E-01
 GRADIENT:   1.4070E-01  4.6822E-01  3.3527E-01  2.4726E-01 -7.3407E-01 -1.1715E-02  4.7562E-02  0.0000E+00  7.7198E-03  2.6606E-02
            -2.6947E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1472.26087012218        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1001
 NPARAMETR:  9.9854E-01  6.7906E-01  1.4657E+00  1.3341E+00  1.0435E+00  9.9138E-01  1.2660E+00  1.0000E-02  8.7320E-01  6.2726E-01
             2.3680E+00
 PARAMETER:  9.8536E-02 -2.8704E-01  4.8231E-01  3.8823E-01  1.4255E-01  9.1341E-02  3.3587E-01 -4.9147E+00 -3.5591E-02 -3.6640E-01
             9.6205E-01
 GRADIENT:   1.4070E-01  4.6822E-01  3.3527E-01  2.4726E-01 -7.3407E-01 -1.1715E-02  4.7562E-02  0.0000E+00  7.7198E-03  2.6606E-02
            -2.6947E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1001
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.8333E-05 -1.3555E-03  8.6516E-05 -1.0366E-02 -1.5577E-02
 SE:             2.9240E-02  1.3514E-02  1.0109E-04  2.4097E-02  1.3045E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9895E-01  9.2010E-01  3.9208E-01  6.6707E-01  2.3245E-01

 ETASHRINKSD(%)  2.0439E+00  5.4727E+01  9.9661E+01  1.9272E+01  5.6298E+01
 ETASHRINKVR(%)  4.0461E+00  7.9504E+01  9.9999E+01  3.4829E+01  8.0901E+01
 EBVSHRINKSD(%)  2.0646E+00  5.5850E+01  9.9645E+01  1.8711E+01  5.7473E+01
 EBVSHRINKVR(%)  4.0866E+00  8.0508E+01  9.9999E+01  3.3921E+01  8.1914E+01
 RELATIVEINF(%)  9.1906E+01  2.1888E-01  8.8459E-05  9.7227E-01  1.1110E+00
 EPSSHRINKSD(%)  2.8136E+01
 EPSSHRINKVR(%)  4.8355E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1472.2608701221784     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -737.11004355844022     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1472.261       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  6.79E-01  1.47E+00  1.33E+00  1.04E+00  9.91E-01  1.27E+00  1.00E-02  8.73E-01  6.27E-01  2.37E+00
 


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
+        1.09E+03
 
 TH 2
+       -3.49E+01  2.89E+02
 
 TH 3
+        6.13E+00  5.61E+01  6.18E+01
 
 TH 4
+       -3.91E+01  3.52E+02 -1.81E+00  5.43E+02
 
 TH 5
+       -2.64E+00 -1.85E+02 -1.55E+02 -5.54E+01  4.41E+02
 
 TH 6
+        1.69E-01 -6.42E+00  2.49E+00 -1.08E+01 -4.19E+00  1.85E+02
 
 TH 7
+        4.66E-01  1.63E+00  2.93E+00 -2.94E+00 -3.70E-01  3.69E-01  6.96E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.32E+00 -1.38E+01  6.31E+00 -7.66E-01 -1.42E+00  9.69E-01  1.62E+01  0.00E+00  1.19E+02
 
 TH10
+       -4.28E+00 -1.02E+00 -1.71E-01 -6.92E+00 -4.93E+00  9.91E-02  1.63E+00  0.00E+00  3.53E-01  1.25E+01
 
 TH11
+       -1.27E+01 -2.65E+00 -1.01E+00 -1.04E+01 -1.22E+01  3.04E+00  3.03E+00  0.00E+00  8.75E+00  2.05E+01  6.11E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.458
Stop Time:
Wed Sep 29 12:17:33 CDT 2021
