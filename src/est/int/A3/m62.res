Wed Sep 29 00:17:03 CDT 2021
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
$DATA ../../../../data/int/A3/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -283.163917372744        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0681E+02  2.5290E+02  2.8515E+02  9.5068E+01  3.6031E+02  6.8360E+01 -1.4040E+02 -1.9400E+02 -4.5643E+01 -1.8469E+02
            -6.5466E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2712.99753494809        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0373E+00  9.0653E-01  8.2238E-01  1.0908E+00  7.9713E-01  7.5674E-01  8.7573E-01  9.7183E-01  1.0110E+00  1.0109E+00
             2.9026E+00
 PARAMETER:  1.3658E-01  1.8715E-03 -9.5547E-02  1.8690E-01 -1.2674E-01 -1.7873E-01 -3.2692E-02  7.1424E-02  1.1094E-01  1.1086E-01
             1.1656E+00
 GRADIENT:   1.1090E+02  5.0560E+01  1.8114E+01  4.2992E+01 -9.0959E+00 -3.7374E+01  1.5492E+01  1.2047E+01  2.5842E+00  6.3403E+00
             1.2362E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2727.93387683447        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0415E+00  7.3274E-01  5.7425E-01  1.1218E+00  5.7502E-01  8.0573E-01  4.2807E-01  6.9651E-01  1.0820E+00  1.0892E+00
             2.8328E+00
 PARAMETER:  1.4066E-01 -2.1096E-01 -4.5468E-01  2.1493E-01 -4.5334E-01 -1.1601E-01 -7.4847E-01 -2.6167E-01  1.7883E-01  1.8544E-01
             1.1413E+00
 GRADIENT:   5.8738E+01  1.2069E+02  1.7690E+01 -5.2253E+00 -1.0836E+02 -1.9202E+01  1.3554E+00  1.3088E+01  1.4308E+01  8.8587E+00
             8.2501E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2752.89788378665        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  1.0154E+00  4.2564E-01  3.7674E-01  1.1986E+00  3.6520E-01  9.0127E-01  5.3733E-02  1.7222E-01  1.0547E+00  9.9357E-01
             2.6505E+00
 PARAMETER:  1.1527E-01 -7.5417E-01 -8.7619E-01  2.8118E-01 -9.0732E-01 -3.9536E-03 -2.8237E+00 -1.6590E+00  1.5330E-01  9.3548E-02
             1.0748E+00
 GRADIENT:  -9.5373E+00 -9.4866E+00  6.6546E+01  5.0433E+01 -7.5870E+01  1.8841E+01 -1.2555E-01  6.2057E-01  5.9710E-01  6.9923E+00
            -3.2907E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2757.38118931578        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      578
 NPARAMETR:  1.0195E+00  4.3771E-01  3.6331E-01  1.1582E+00  3.7810E-01  8.5065E-01  3.8442E-02  1.0085E-01  1.0533E+00  9.6867E-01
             2.6837E+00
 PARAMETER:  1.1929E-01 -7.2619E-01 -9.1250E-01  2.4682E-01 -8.7260E-01 -6.1759E-02 -3.1586E+00 -2.1941E+00  1.5190E-01  6.8167E-02
             1.0872E+00
 GRADIENT:  -4.8055E-01  5.1380E-01  1.5176E-01 -2.9518E+00 -1.2209E+00 -1.0972E+00 -6.6300E-02  2.0364E-01  1.2698E-01 -1.0819E-01
             7.8784E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2757.46421560146        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  1.0196E+00  4.4063E-01  3.6651E-01  1.1608E+00  3.8099E-01  8.5302E-01  2.7831E-02  1.9995E-02  1.0516E+00  9.7102E-01
             2.6834E+00
 PARAMETER:  1.1945E-01 -7.1955E-01 -9.0373E-01  2.4913E-01 -8.6497E-01 -5.8968E-02 -3.4816E+00 -3.8123E+00  1.5031E-01  7.0595E-02
             1.0871E+00
 GRADIENT:   3.7869E-02  3.7535E-01  1.4635E-01 -1.7049E-01 -4.1888E-01 -3.4984E-02 -3.4004E-02  8.1456E-03  8.6715E-03 -4.8223E-03
            -2.7954E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2761.54105626341        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  1.0168E+00  4.3731E-01  3.6415E-01  1.1628E+00  3.8383E-01  8.5089E-01  7.6448E-01  1.0000E-02  1.0518E+00  9.6629E-01
             2.6783E+00
 PARAMETER:  1.1666E-01 -7.2712E-01 -9.1020E-01  2.5083E-01 -8.5755E-01 -6.1468E-02 -1.6856E-01 -5.7758E+00  1.5053E-01  6.5712E-02
             1.0852E+00
 GRADIENT:  -1.0902E+01 -1.4112E+01 -7.6429E-03 -1.6863E+00  3.2078E+01 -1.7589E+00  3.2161E+00  0.0000E+00 -3.0254E+00  2.3503E+01
             4.5642E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2772.90828391595        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1115
 NPARAMETR:  1.0288E+00  3.3797E-01  2.4103E-01  1.1109E+00  2.7441E-01  8.6688E-01  1.0923E+00  1.0000E-02  1.1897E+00  7.0481E-01
             2.5088E+00
 PARAMETER:  1.2844E-01 -9.8481E-01 -1.3228E+00  2.0518E-01 -1.1931E+00 -4.2854E-02  1.8829E-01 -5.7758E+00  2.7374E-01 -2.4982E-01
             1.0198E+00
 GRADIENT:   2.4658E+01  4.7909E+01  1.6980E+01  2.2802E+01 -4.4708E+01 -6.4861E-02  7.7717E-01  0.0000E+00 -2.0693E-01 -2.0258E+01
            -3.7016E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2777.70517551603        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  1.0213E+00  2.7584E-01  1.8740E-01  1.0267E+00  2.3227E-01  8.7027E-01  1.1513E+00  1.0000E-02  1.2600E+00  7.5631E-01
             2.5106E+00
 PARAMETER:  1.2112E-01 -1.1880E+00 -1.5745E+00  1.2634E-01 -1.3599E+00 -3.8950E-02  2.4092E-01 -5.7758E+00  3.3108E-01 -1.7931E-01
             1.0205E+00
 GRADIENT:   7.4482E-02 -1.7778E+00  3.0686E-01 -2.7105E-01  1.2230E+00 -6.1245E-02 -6.7368E-01  0.0000E+00 -1.0504E-01  4.0908E-02
            -2.8914E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -2777.71312789108        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1382
 NPARAMETR:  1.0213E+00  2.7576E-01  1.8686E-01  1.0263E+00  2.3185E-01  8.7049E-01  1.1577E+00  1.0000E-02  1.2613E+00  7.5577E-01
             2.5103E+00
 PARAMETER:  1.2110E-01 -1.1882E+00 -1.5774E+00  1.2593E-01 -1.3617E+00 -3.8699E-02  2.4642E-01 -5.7758E+00  3.3217E-01 -1.8002E-01
             1.0204E+00
 GRADIENT:  -5.5336E-02 -1.8736E-01 -7.8902E-02  9.4642E-02  2.9721E-01 -8.3626E-03  1.9151E-02  0.0000E+00 -9.2858E-03 -7.4756E-02
            -1.3176E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1382
 NO. OF SIG. DIGITS IN FINAL EST.:  2.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3479E-03  5.8035E-03  8.8808E-05 -5.1718E-03  2.6557E-03
 SE:             2.9288E-02  2.3067E-02  2.9664E-04  2.7816E-02  2.5940E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6329E-01  8.0135E-01  7.6465E-01  8.5250E-01  9.1846E-01

 ETASHRINKSD(%)  1.8804E+00  2.2724E+01  9.9006E+01  6.8120E+00  1.3098E+01
 ETASHRINKVR(%)  3.7254E+00  4.0284E+01  9.9990E+01  1.3160E+01  2.4481E+01
 EBVSHRINKSD(%)  1.9472E+00  2.2206E+01  9.9173E+01  5.3805E+00  1.3316E+01
 EBVSHRINKVR(%)  3.8565E+00  3.9481E+01  9.9993E+01  1.0471E+01  2.4859E+01
 RELATIVEINF(%)  9.6119E+01  1.2386E+01  5.9050E-04  4.5878E+01  6.1296E+00
 EPSSHRINKSD(%)  1.8765E+01
 EPSSHRINKVR(%)  3.4009E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2777.7131278910788     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1123.6237681226680     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2777.713       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.76E-01  1.87E-01  1.03E+00  2.32E-01  8.70E-01  1.16E+00  1.00E-02  1.26E+00  7.56E-01  2.51E+00
 


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
+        1.35E+03
 
 TH 2
+       -5.68E+00  6.67E+03
 
 TH 3
+        9.47E+01 -2.09E+03  1.83E+04
 
 TH 4
+       -1.19E+01 -8.99E+01 -6.32E+02  5.97E+02
 
 TH 5
+       -1.83E+01 -4.66E+03 -1.74E+04 -2.69E+02  2.74E+04
 
 TH 6
+        5.20E+00 -9.84E+00  5.24E+01 -6.90E+00  4.82E+01  2.40E+02
 
 TH 7
+       -2.02E+00  3.57E+01  1.24E+02 -2.79E+00 -4.97E+01 -1.29E-01  5.53E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.04E+00 -4.34E+00  1.64E+02 -4.38E+00  2.56E+01 -1.77E+00  1.02E+00  0.00E+00  9.54E+01
 
 TH10
+       -2.71E+00 -8.74E+00 -3.04E+00  1.41E+01 -2.29E+01 -2.02E-01  4.55E+00  0.00E+00  4.79E-01  2.01E+02
 
 TH11
+       -1.95E+01 -9.53E+00 -1.34E+02 -7.28E+00  9.12E+01  4.47E+00  1.22E+01  0.00E+00  5.48E+00  1.49E+01  1.70E+02
 
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
 #CPUT: Total CPU Time in Seconds,       49.406
Stop Time:
Wed Sep 29 00:17:54 CDT 2021
