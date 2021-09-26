Sat Sep 25 00:58:33 CDT 2021
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
$DATA ../../../../data/int/SL2/dat17.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      998
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

 TOT. NO. OF OBS RECS:      898
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -844.949962277226        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.6262E+01 -5.7919E+01  2.6031E+02 -1.2015E+01  4.1231E+01 -6.3943E+00 -1.0794E+02 -2.3759E+02 -9.1261E+01 -8.1322E+01
            -5.6309E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2819.03789983304        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0782E+00  1.3437E+00  8.3957E-01  8.6283E-01  1.1688E+00  1.0022E+00  1.1915E+00  9.6246E-01  9.3545E-01  1.2543E+00
             2.4570E+00
 PARAMETER:  1.7525E-01  3.9543E-01 -7.4867E-02 -4.7537E-02  2.5597E-01  1.0224E-01  2.7521E-01  6.1734E-02  3.3272E-02  3.2658E-01
             9.9895E-01
 GRADIENT:   6.6675E+01 -1.7017E+01 -3.3154E+00 -8.4298E-01 -2.2897E+01 -1.1316E+00  1.7473E+01  5.1301E+00  4.0601E-01 -1.1014E+01
            -1.0087E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2828.70316712537        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0890E+00  1.6698E+00  8.5295E-01  7.4775E-01  1.3890E+00  1.0020E+00  9.2378E-01  4.5952E-01  1.1843E+00  1.5278E+00
             2.4931E+00
 PARAMETER:  1.8527E-01  6.1271E-01 -5.9051E-02 -1.9068E-01  4.2861E-01  1.0198E-01  2.0714E-02 -6.7758E-01  2.6915E-01  5.2382E-01
             1.0135E+00
 GRADIENT:   8.5850E+01  8.5561E+01  1.1665E+01  7.3489E+01 -2.0590E+01 -2.8350E+00  2.0388E+00 -2.6372E-01  5.7403E+00  6.0106E+00
            -5.8642E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2841.62682163353        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0422E+00  1.8680E+00  6.4253E-01  5.3960E-01  1.6010E+00  1.0185E+00  8.3505E-01  2.4042E-01  1.2055E+00  1.6282E+00
             2.5301E+00
 PARAMETER:  1.4134E-01  7.2485E-01 -3.4235E-01 -5.1692E-01  5.7063E-01  1.1828E-01 -8.0261E-02 -1.3254E+00  2.8690E-01  5.8750E-01
             1.0283E+00
 GRADIENT:  -8.2905E+00 -1.1429E+01  3.9026E+00  3.5782E+00  2.0381E+00  6.1983E+00 -4.3940E+00 -3.8440E-02 -4.2440E+00  1.8072E+00
             2.7962E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2850.06282408474        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0459E+00  2.2619E+00  3.3830E-01  2.8637E-01  1.9054E+00  9.9965E-01  7.3510E-01  3.7859E-02  1.9662E+00  1.7971E+00
             2.4930E+00
 PARAMETER:  1.4492E-01  9.1623E-01 -9.8382E-01 -1.1505E+00  7.4472E-01  9.9648E-02 -2.0775E-01 -3.1739E+00  7.7609E-01  6.8617E-01
             1.0135E+00
 GRADIENT:  -1.9807E+00  2.1559E+01 -1.7519E+00  7.2416E+00  4.0800E+00 -1.4902E+00 -3.9233E+00  1.4589E-03  7.9636E-01  1.6119E+00
            -3.6710E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2851.34084790583        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      479
 NPARAMETR:  1.0530E+00  2.4000E+00  2.4516E-01  1.9818E-01  2.0298E+00  1.0051E+00  7.2659E-01  1.2117E-02  2.4040E+00  1.8516E+00
             2.4963E+00
 PARAMETER:  1.5161E-01  9.7548E-01 -1.3058E+00 -1.5186E+00  8.0791E-01  1.0514E-01 -2.1940E-01 -4.3132E+00  9.7712E-01  7.1607E-01
             1.0148E+00
 GRADIENT:   1.8113E+00 -1.8560E+00 -1.4736E+00  1.2980E+00  1.2535E+00 -7.9190E-02  5.7229E-01  1.7715E-04 -2.6632E-01  1.2595E-01
             4.8160E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2851.40940795747        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      655            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0520E+00  2.4081E+00  2.5999E-01  1.9200E-01  2.0395E+00  1.0052E+00  7.2214E-01  1.3069E-02  2.5114E+00  1.8481E+00
             2.4959E+00
 PARAMETER:  1.5066E-01  9.7882E-01 -1.2471E+00 -1.5502E+00  8.1272E-01  1.0516E-01 -2.2554E-01 -4.2375E+00  1.0208E+00  7.1418E-01
             1.0146E+00
 GRADIENT:   1.1349E+01  4.5864E+01  6.1741E-02  1.3438E+00  3.7581E+00  8.6555E-01  4.8489E-01  2.0682E-04  5.8639E-01  1.1694E+00
             1.8111E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2851.40943666139        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      830
 NPARAMETR:  1.0519E+00  2.4081E+00  2.5986E-01  1.9183E-01  2.0396E+00  1.0051E+00  7.2205E-01  1.0000E-02  2.5129E+00  1.8482E+00
             2.4958E+00
 PARAMETER:  1.5062E-01  9.7885E-01 -1.2476E+00 -1.5512E+00  8.1278E-01  1.0513E-01 -2.2567E-01 -4.5391E+00  1.0214E+00  7.1421E-01
             1.0146E+00
 GRADIENT:  -7.9693E-02 -2.4211E-01 -1.0140E-03 -3.1475E-02 -3.5203E-02 -9.2563E-03 -2.0399E-02  0.0000E+00  3.3184E-03 -2.4158E-03
            -4.9248E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2851.40945030750        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  1.0520E+00  2.4081E+00  2.5998E-01  1.9197E-01  2.0396E+00  1.0052E+00  7.2213E-01  1.0000E-02  2.5118E+00  1.8481E+00
             2.4959E+00
 PARAMETER:  1.5065E-01  9.7883E-01 -1.2471E+00 -1.5504E+00  8.1275E-01  1.0516E-01 -2.2555E-01 -4.5543E+00  1.0210E+00  7.1418E-01
             1.0146E+00
 GRADIENT:  -8.4622E-03 -5.2052E-02  9.2378E-04  5.6745E-04 -5.2071E-03 -1.1726E-03 -2.3293E-03  0.0000E+00  1.1918E-03 -2.0602E-03
            -2.5309E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      922
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5964E-03 -2.4995E-02 -4.3052E-05  3.4897E-02 -2.1231E-02
 SE:             2.9473E-02  2.6590E-02  3.1703E-05  1.7053E-02  2.6544E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5680E-01  3.4720E-01  1.7448E-01  4.0718E-02  4.2382E-01

 ETASHRINKSD(%)  1.2608E+00  1.0921E+01  9.9894E+01  4.2871E+01  1.1073E+01
 ETASHRINKVR(%)  2.5058E+00  2.0648E+01  1.0000E+02  6.7362E+01  2.0919E+01
 EBVSHRINKSD(%)  1.4364E+00  1.1129E+01  9.9888E+01  4.9990E+01  8.1459E+00
 EBVSHRINKVR(%)  2.8522E+00  2.1019E+01  1.0000E+02  7.4990E+01  1.5628E+01
 RELATIVEINF(%)  9.7107E+01  2.1770E+01  6.1349E-05  5.7570E+00  5.4379E+01
 EPSSHRINKSD(%)  1.6950E+01
 EPSSHRINKVR(%)  3.1027E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2851.4094503075021     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1200.9958446719099     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2851.409       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.41E+00  2.60E-01  1.92E-01  2.04E+00  1.01E+00  7.22E-01  1.00E-02  2.51E+00  1.85E+00  2.50E+00
 


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
+        9.62E+02
 
 TH 2
+       -8.80E+00  2.98E+02
 
 TH 3
+        3.08E+00  4.00E+01  1.79E+02
 
 TH 4
+       -2.24E+01  3.33E+02 -2.35E+02  1.37E+03
 
 TH 5
+       -1.55E+00 -1.83E+01 -1.96E+01  8.13E+01  6.63E+01
 
 TH 6
+        5.05E+00 -2.88E+00  1.78E+00 -1.06E+01 -8.42E-01  1.82E+02
 
 TH 7
+        3.59E+00 -9.45E-01 -9.28E-01 -2.22E+01 -5.79E-01 -1.49E+00  2.50E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.09E-01 -4.02E+00 -1.68E+01  5.41E+01 -2.80E-01 -1.05E-01  5.91E+00  0.00E+00  6.41E+00
 
 TH10
+        3.29E-01 -2.15E+00  1.00E+01  7.12E+00 -5.49E+00  2.09E-01  5.15E+00  0.00E+00  6.76E-01  3.90E+01
 
 TH11
+       -1.35E+01 -1.48E+01 -6.33E+00  2.38E+00  3.48E-01  3.64E+00  6.41E+00  0.00E+00  2.15E+00  2.83E+00  1.87E+02
 
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
 #CPUT: Total CPU Time in Seconds,       34.702
Stop Time:
Sat Sep 25 00:59:09 CDT 2021
