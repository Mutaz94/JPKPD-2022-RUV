Sat Sep 18 08:31:29 CDT 2021
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
$DATA ../../../../data/spa/B/dat46.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1706.51301285750        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2014E+02 -9.4185E+01 -4.8171E+01 -6.9120E+01  6.2322E+01  2.1083E+01  1.1196E+00  1.1533E+01  2.2011E+01  1.4304E+01
             6.9333E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1721.49987972635        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      110
 NPARAMETR:  9.5636E-01  1.0654E+00  1.0580E+00  1.0101E+00  9.8376E-01  9.0846E-01  9.7475E-01  9.2810E-01  8.5917E-01  9.1363E-01
             8.0106E-01
 PARAMETER:  5.5378E-02  1.6337E-01  1.5641E-01  1.1009E-01  8.3625E-02  3.9992E-03  7.4422E-02  2.5381E-02 -5.1790E-02  9.6677E-03
            -1.2182E-01
 GRADIENT:  -1.8608E+01 -3.3931E+00  1.7097E+01 -2.7843E+01 -1.2348E+01 -1.7723E+01 -6.8479E+00 -2.5123E-01 -1.1798E+01 -5.2162E+00
            -1.8299E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1723.64550126031        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      286
 NPARAMETR:  9.6072E-01  9.5321E-01  8.6182E-01  1.0888E+00  8.5981E-01  9.3571E-01  1.1865E+00  5.2630E-01  8.0752E-01  8.0008E-01
             7.9633E-01
 PARAMETER:  5.9924E-02  5.2079E-02 -4.8705E-02  1.8510E-01 -5.1043E-02  3.3547E-02  2.7100E-01 -5.4189E-01 -1.1378E-01 -1.2304E-01
            -1.2774E-01
 GRADIENT:  -8.7325E+00  8.2279E+00 -8.7075E+00  2.9656E+01  2.9421E+01 -5.7406E+00 -9.1220E-01 -4.1097E-02 -5.3772E+00 -1.1161E+01
            -1.8237E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1724.96429967017        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      462
 NPARAMETR:  9.6567E-01  9.1450E-01  8.0485E-01  1.0969E+00  8.0640E-01  9.5167E-01  1.2135E+00  3.5142E-01  8.1016E-01  8.2641E-01
             8.2806E-01
 PARAMETER:  6.5067E-02  1.0619E-02 -1.1710E-01  1.9245E-01 -1.1518E-01  5.0459E-02  2.9350E-01 -9.4577E-01 -1.1052E-01 -9.0663E-02
            -8.8668E-02
 GRADIENT:   2.7594E+00  2.5786E+00 -1.2125E+00  4.3037E+00 -1.0591E+00  1.0492E+00 -1.5111E-01  4.5930E-01 -1.5453E-01  1.3904E+00
             1.3259E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1725.04080354576        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  9.6422E-01  8.5606E-01  7.9335E-01  1.1263E+00  7.7906E-01  9.4910E-01  1.2860E+00  2.5312E-01  7.9196E-01  8.1433E-01
             8.2571E-01
 PARAMETER:  6.3564E-02 -5.5418E-02 -1.3149E-01  2.1892E-01 -1.4967E-01  4.7754E-02  3.5155E-01 -1.2739E+00 -1.3324E-01 -1.0539E-01
            -9.1514E-02
 GRADIENT:  -4.7740E-01 -1.4506E+00 -1.4423E+00 -2.3650E-01  1.9361E+00  1.3429E-01  6.1962E-02  8.1787E-02  3.3956E-01  3.2737E-01
             3.2980E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1725.07249937538        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  9.6497E-01  9.1407E-01  7.5954E-01  1.0893E+00  7.8288E-01  9.4959E-01  1.2228E+00  1.3603E-01  8.0747E-01  8.0736E-01
             8.2560E-01
 PARAMETER:  6.4338E-02  1.0147E-02 -1.7504E-01  1.8557E-01 -1.4477E-01  4.8270E-02  3.0115E-01 -1.8949E+00 -1.1385E-01 -1.1398E-01
            -9.1641E-02
 GRADIENT:  -5.3307E-02 -1.5364E-01 -5.2316E-01  1.7012E-01  4.6438E-01 -1.0869E-02  3.8438E-02  1.7358E-02  2.1270E-02  1.1926E-01
             2.7598E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1725.07641144249        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  9.6505E-01  9.1896E-01  7.5298E-01  1.0858E+00  7.8111E-01  9.4970E-01  1.2179E+00  6.1588E-02  8.0867E-01  8.0598E-01
             8.2567E-01
 PARAMETER:  6.4429E-02  1.5483E-02 -1.8371E-01  1.8230E-01 -1.4704E-01  4.8393E-02  2.9712E-01 -2.6873E+00 -1.1237E-01 -1.1569E-01
            -9.1563E-02
 GRADIENT:  -3.5119E-02  8.2391E-02  2.9178E-02  1.4482E-01 -7.6475E-02 -6.5037E-03 -6.9307E-03  1.3096E-03 -1.0663E-03 -2.0792E-03
            -2.1685E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1725.07706343073        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  9.6511E-01  9.2266E-01  7.5070E-01  1.0833E+00  7.8141E-01  9.4977E-01  1.2141E+00  1.3463E-02  8.0977E-01  8.0566E-01
             8.2575E-01
 PARAMETER:  6.4482E-02  1.9501E-02 -1.8675E-01  1.8001E-01 -1.4665E-01  4.8466E-02  2.9404E-01 -4.2078E+00 -1.1100E-01 -1.1610E-01
            -9.1464E-02
 GRADIENT:   8.9863E-03  5.5369E-03  6.3417E-03  9.8475E-03  4.2029E-03  1.1218E-03 -4.7530E-04  5.4663E-05 -4.2301E-03 -8.0670E-03
             2.1138E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1725.07707655697        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  9.6510E-01  9.2277E-01  7.5062E-01  1.0832E+00  7.8141E-01  9.4977E-01  1.2140E+00  1.0000E-02  8.0983E-01  8.0567E-01
             8.2574E-01
 PARAMETER:  6.4478E-02  1.9626E-02 -1.8685E-01  1.7993E-01 -1.4665E-01  4.8464E-02  2.9392E-01 -4.7031E+00 -1.1093E-01 -1.1608E-01
            -9.1474E-02
 GRADIENT:  -3.4419E-03 -2.2436E-04  4.1314E-04  5.9595E-04 -1.0662E-03 -3.5958E-04 -4.1312E-04  0.0000E+00  1.6315E-03  2.3332E-05
            -1.8626E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1255
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0845E-05  3.6459E-04 -4.9701E-04 -3.4715E-03 -8.5643E-03
 SE:             2.9872E-02  2.1465E-02  2.1558E-04  2.4771E-02  2.3633E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9944E-01  9.8645E-01  2.1142E-02  8.8854E-01  7.1706E-01

 ETASHRINKSD(%)  1.0000E-10  2.8090E+01  9.9278E+01  1.7015E+01  2.0826E+01
 ETASHRINKVR(%)  1.0000E-10  4.8289E+01  9.9995E+01  3.1135E+01  3.7315E+01
 EBVSHRINKSD(%)  3.2646E-01  2.8048E+01  9.9359E+01  1.7279E+01  1.9964E+01
 EBVSHRINKVR(%)  6.5186E-01  4.8229E+01  9.9996E+01  3.1572E+01  3.5942E+01
 RELATIVEINF(%)  9.9033E+01  3.3567E+00  4.4124E-04  5.3900E+00  5.0350E+00
 EPSSHRINKSD(%)  4.5120E+01
 EPSSHRINKVR(%)  6.9882E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1725.0770765569735     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -989.92624999323527     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1725.077       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.65E-01  9.23E-01  7.51E-01  1.08E+00  7.81E-01  9.50E-01  1.21E+00  1.00E-02  8.10E-01  8.06E-01  8.26E-01
 


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
+        1.31E+03
 
 TH 2
+       -5.84E+00  4.78E+02
 
 TH 3
+        1.91E+01  2.83E+02  8.76E+02
 
 TH 4
+       -9.48E+00  3.75E+02 -3.72E+02  9.93E+02
 
 TH 5
+       -2.94E+00 -4.95E+02 -1.10E+03  4.32E+02  1.78E+03
 
 TH 6
+        2.99E-01 -1.03E+00  4.71E+00 -2.33E+00 -1.88E+00  2.19E+02
 
 TH 7
+        1.58E+00  3.28E+01 -8.49E+00 -8.53E+00 -4.77E+00  2.32E-01  4.29E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.47E+00 -2.50E+01 -4.02E+01  2.68E+01  1.13E+01 -2.94E+00  2.15E+01  0.00E+00  1.55E+02
 
 TH10
+       -9.77E-01 -1.13E+01 -8.77E+01 -2.54E+01 -4.76E+01  3.54E-01  1.31E+01  0.00E+00  1.22E+01  1.33E+02
 
 TH11
+       -6.65E+00 -1.47E+01 -4.90E+01 -2.85E+00  1.09E+01  6.32E-01  5.04E+00  0.00E+00  1.02E+01  2.81E+01  3.06E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.276
Stop Time:
Sat Sep 18 08:31:51 CDT 2021
