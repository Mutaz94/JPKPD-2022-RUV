Sat Sep 25 01:50:24 CDT 2021
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
$DATA ../../../../data/int/SL3/dat3.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      983
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

 TOT. NO. OF OBS RECS:      883
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   257.528395861847        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5483E+00 -4.8350E+01  1.2914E+02  5.5418E+01  1.4008E+02  1.2283E+00 -1.4835E+02 -2.1738E+02 -1.3764E+02 -3.9147E+01
            -7.6806E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2354.44377545155        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1205E+00  1.5906E+00  9.2445E-01  8.4505E-01  1.2975E+00  9.0969E-01  1.1825E+00  1.0219E+00  9.8478E-01  9.3000E-01
             5.1663E+00
 PARAMETER:  2.1378E-01  5.6413E-01  2.1441E-02 -6.8360E-02  3.6046E-01  5.3478E-03  2.6765E-01  1.2164E-01  8.4665E-02  2.7425E-02
             1.7422E+00
 GRADIENT:   6.9549E+01  3.3105E+01 -2.5645E+01  4.2579E+01  4.4497E+01 -2.6378E+01  4.4275E+01  2.8463E+00  4.5247E+00 -4.7050E+00
             7.4883E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2482.79897267119        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0788E+00  2.1217E+00  6.0043E-01  6.4422E-01  1.4666E+00  9.1902E-01  8.9757E-01  1.8986E+00  2.9995E+00  1.3571E+00
             4.0101E+00
 PARAMETER:  1.7583E-01  8.5222E-01 -4.1011E-01 -3.3971E-01  4.8296E-01  1.5556E-02 -8.0676E-03  7.4111E-01  1.1984E+00  4.0533E-01
             1.4888E+00
 GRADIENT:   4.1315E+01  1.1109E+02 -1.5814E+01  8.5461E+01 -5.1612E+01 -2.8756E+01  3.4857E+01  1.1629E+01  5.6430E+01 -5.5525E+00
             5.2533E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2631.97660803314        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0262E+00  1.6614E+00  1.0556E+00  6.7823E-01  1.3945E+00  1.0013E+00  6.9036E-01  9.9222E-01  1.4540E+00  1.3937E+00
             2.8816E+00
 PARAMETER:  1.2591E-01  6.0767E-01  1.5407E-01 -2.8827E-01  4.3251E-01  1.0134E-01 -2.7055E-01  9.2193E-02  4.7431E-01  4.3195E-01
             1.1583E+00
 GRADIENT:  -2.4065E+01  5.4722E+00  2.2778E+00  1.7275E+01 -7.6630E+00  2.0154E-01 -6.9124E+00  7.8329E-01  6.9616E+00 -1.1284E+01
            -1.9119E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2639.55629633780        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0346E+00  1.9983E+00  7.3964E-01  4.4561E-01  1.6108E+00  9.9468E-01  7.0052E-01  5.9649E-01  1.7136E+00  1.5757E+00
             2.8755E+00
 PARAMETER:  1.3397E-01  7.9231E-01 -2.0159E-01 -7.0832E-01  5.7676E-01  9.4667E-02 -2.5593E-01 -4.1669E-01  6.3862E-01  5.5468E-01
             1.1562E+00
 GRADIENT:  -6.5394E+00  2.5262E+00 -5.6675E-01  3.9254E+00  1.2517E+00 -1.7209E+00 -8.8905E-02  1.8821E-01 -2.7178E+00 -7.3578E+00
            -7.3451E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2644.83421008526        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  1.0366E+00  2.4827E+00  3.0059E-01  1.2933E-01  1.9899E+00  9.8870E-01  7.2201E-01  3.7779E-02  3.3925E+00  1.9643E+00
             2.8378E+00
 PARAMETER:  1.3598E-01  1.0094E+00 -1.1020E+00 -1.9454E+00  7.8809E-01  8.8637E-02 -2.2572E-01 -3.1760E+00  1.3216E+00  7.7514E-01
             1.1430E+00
 GRADIENT:  -1.2742E-01  3.2473E+01 -1.5776E+00  4.4604E+00  6.9592E+00 -4.0201E+00  1.8501E+01  1.7915E-03  7.3741E+00  8.5003E+00
            -1.3723E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2645.95463332614        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      517
 NPARAMETR:  1.0416E+00  2.4812E+00  3.1399E-01  1.3812E-01  1.9673E+00  1.0016E+00  6.6944E-01  4.5012E-02  3.1788E+00  1.8956E+00
             2.8624E+00
 PARAMETER:  1.4079E-01  1.0087E+00 -1.0584E+00 -1.8796E+00  7.7667E-01  1.0159E-01 -3.0132E-01 -3.0008E+00  1.2565E+00  7.3954E-01
             1.1516E+00
 GRADIENT:   4.7645E-01  1.2842E+01 -1.1545E+00  2.7437E+00  2.8124E+00 -7.8534E-02  2.6064E-01  2.6377E-03  1.4527E+00  3.5098E-01
            -1.0331E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2646.04418964118        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      692
 NPARAMETR:  1.0413E+00  2.4741E+00  3.3467E-01  1.3541E-01  1.9602E+00  1.0014E+00  6.6749E-01  4.4197E-02  3.1719E+00  1.8937E+00
             2.8636E+00
 PARAMETER:  1.4046E-01  1.0059E+00 -9.9462E-01 -1.8994E+00  7.7304E-01  1.0140E-01 -3.0423E-01 -3.0191E+00  1.2543E+00  7.3851E-01
             1.1521E+00
 GRADIENT:   1.0952E-01 -4.3843E-01 -2.6909E-02  9.5991E-03 -9.5432E-02 -3.7024E-02 -1.7252E-03  2.5249E-03  1.1969E-01 -8.9159E-02
             5.5094E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2646.04443480855        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  1.0412E+00  2.4748E+00  3.3470E-01  1.3514E-01  1.9611E+00  1.0015E+00  6.6766E-01  4.3915E-02  3.1691E+00  1.8950E+00
             2.8630E+00
 PARAMETER:  1.4040E-01  1.0061E+00 -9.9453E-01 -1.9015E+00  7.7350E-01  1.0151E-01 -3.0398E-01 -3.0255E+00  1.2535E+00  7.3923E-01
             1.1519E+00
 GRADIENT:   8.3289E-03  3.2339E-02  2.3345E-03 -7.3007E-03 -3.2930E-02 -8.3526E-04 -9.6685E-03  2.4918E-03  3.4522E-03 -1.0078E-02
            -1.0320E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      818
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7815E-03 -3.1440E-02 -7.6568E-05  3.8794E-02 -2.5101E-02
 SE:             2.9357E-02  2.5110E-02  7.0884E-05  1.7075E-02  2.5905E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5161E-01  2.1055E-01  2.8006E-01  2.3085E-02  3.3255E-01

 ETASHRINKSD(%)  1.6495E+00  1.5877E+01  9.9763E+01  4.2798E+01  1.3216E+01
 ETASHRINKVR(%)  3.2718E+00  2.9233E+01  9.9999E+01  6.7280E+01  2.4686E+01
 EBVSHRINKSD(%)  1.8301E+00  1.3139E+01  9.9622E+01  5.4900E+01  9.0721E+00
 EBVSHRINKVR(%)  3.6268E+00  2.4552E+01  9.9999E+01  7.9660E+01  1.7321E+01
 RELATIVEINF(%)  9.6295E+01  2.2636E+01  1.0093E-03  5.7668E+00  6.1017E+01
 EPSSHRINKSD(%)  1.6256E+01
 EPSSHRINKVR(%)  2.9869E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          883
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1622.8454496394520     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2646.0444348085507     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1023.1989851690987     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.97
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2646.044       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.47E+00  3.35E-01  1.35E-01  1.96E+00  1.00E+00  6.68E-01  4.39E-02  3.17E+00  1.90E+00  2.86E+00
 


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
+        9.86E+02
 
 TH 2
+       -1.37E+01  3.16E+02
 
 TH 3
+        2.53E+00  1.86E+01  6.21E+01
 
 TH 4
+       -2.30E+01  4.08E+02 -1.09E+02  1.33E+03
 
 TH 5
+       -2.17E+00 -1.98E+01 -1.21E+01  6.79E+01  6.62E+01
 
 TH 6
+        5.22E+00 -4.41E+00  8.43E-01 -7.02E+00 -9.67E-01  1.81E+02
 
 TH 7
+        3.42E+00 -1.24E+01  1.26E+01  2.47E+01 -4.27E+00  2.36E-01  2.48E+02
 
 TH 8
+        1.24E+00 -2.16E-01  6.72E-01 -1.20E+00  2.66E-01 -6.78E+00 -5.58E-02 -4.27E-01
 
 TH 9
+        2.38E-01 -1.55E+00 -3.25E+00  3.41E+01  5.63E-01 -2.22E-01  6.38E+00  3.23E-02  3.75E+00
 
 TH10
+        7.40E-01 -5.75E+00 -2.80E+00  3.78E+01 -6.36E+00  3.01E-01 -1.90E+00 -1.47E-01  1.37E+00  3.48E+01
 
 TH11
+       -1.49E+01 -1.59E+01 -1.91E+00 -1.07E+01  5.02E-01  2.00E+00  1.02E+01 -2.86E-02  9.08E-01  3.39E+00  1.40E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.551
Stop Time:
Sat Sep 25 01:50:57 CDT 2021
