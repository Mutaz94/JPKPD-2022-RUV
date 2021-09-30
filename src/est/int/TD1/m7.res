Wed Sep 29 05:58:42 CDT 2021
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
$DATA ../../../../data/int/TD1/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3053.74120946052        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9520E+02  4.4911E+01  2.2047E+02  1.0840E+02  1.5504E+02  1.0294E+01 -6.1591E+01 -9.3362E+02 -2.2398E+02 -4.5100E+01
            -2.8026E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3456.89967603550        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.5290E-01  1.1416E+00  9.3960E-01  9.0553E-01  1.0297E+00  8.1769E-01  9.9934E-01  2.1275E+00  9.2667E-01  1.1407E+00
             1.4601E+00
 PARAMETER:  5.1758E-02  2.3241E-01  3.7696E-02  7.6387E-04  1.2923E-01 -1.0127E-01  9.9343E-02  8.5494E-01  2.3837E-02  2.3167E-01
             4.7848E-01
 GRADIENT:  -1.1009E+02 -2.7284E+01 -1.3614E+00 -5.2311E+00 -6.2730E+01 -1.5999E+02  1.1989E+01 -2.5057E+02 -2.3851E+00 -1.7447E+01
             4.6042E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3490.81468205403        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      339
 NPARAMETR:  9.1817E-01  1.3633E+00  9.4461E-01  7.1282E-01  1.1345E+00  8.5842E-01  9.5418E-01  2.5253E+00  5.6897E-01  1.3498E+00
             1.4587E+00
 PARAMETER:  1.4632E-02  4.0992E-01  4.3016E-02 -2.3853E-01  2.2623E-01 -5.2668E-02  5.3097E-02  1.0264E+00 -4.6393E-01  3.9996E-01
             4.7752E-01
 GRADIENT:  -2.0918E+02 -9.7913E+01  1.4260E+01 -1.7265E+02 -1.4587E+02 -1.3923E+02 -7.9007E+00 -2.0238E+02 -1.5988E+01 -4.2270E+01
             4.4629E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3519.80723021124        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      533             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1975E-01  1.2280E+00  9.7278E-01  7.9935E-01  1.0524E+00  8.9413E-01  1.0274E+00  2.8813E+00  4.2662E-01  1.4629E+00
             1.4412E+00
 PARAMETER:  1.6344E-02  3.0541E-01  7.2408E-02 -1.2396E-01  1.5109E-01 -1.1909E-02  1.2706E-01  1.1583E+00 -7.5185E-01  4.8043E-01
             4.6548E-01
 GRADIENT:   4.4448E+01  8.5452E+01  1.7002E+00 -9.1267E+01 -1.1188E+02 -7.2450E+01 -1.6376E+01 -6.6090E+01 -5.8980E+00 -2.0213E+00
             4.5547E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3524.35327402177        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:      681
 NPARAMETR:  9.1845E-01  1.2280E+00  9.7276E-01  8.0103E-01  1.0524E+00  9.2387E-01  1.0414E+00  2.8991E+00  4.2674E-01  1.4629E+00
             1.4412E+00
 PARAMETER:  1.4931E-02  3.0542E-01  7.2380E-02 -1.2186E-01  1.5105E-01  2.0817E-02  1.4052E-01  1.1644E+00 -7.5158E-01  4.8043E-01
             4.6544E-01
 GRADIENT:   5.2852E+01  9.1116E+01  1.1554E+00 -8.7572E+01 -1.1139E+02 -5.3457E+01 -1.1226E+01 -6.3342E+01 -5.2125E+00 -1.8211E+00
             4.5658E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3532.85920076286        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      850
 NPARAMETR:  9.1862E-01  1.2274E+00  9.7294E-01  8.0122E-01  1.0524E+00  1.1223E+00  1.0415E+00  2.8944E+00  4.2616E-01  1.4629E+00
             1.4410E+00
 PARAMETER:  1.5121E-02  3.0488E-01  7.2564E-02 -1.2162E-01  1.5105E-01  2.1539E-01  1.4065E-01  1.1628E+00 -7.5295E-01  4.8044E-01
             4.6531E-01
 GRADIENT:  -1.2063E+02 -1.0339E+02 -6.2730E-01 -1.3120E+02 -1.4564E+02 -2.6910E-01 -2.0880E+01 -1.3607E+02 -9.6070E+00 -1.9495E+01
             4.5113E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3533.91956100004        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1011
 NPARAMETR:  9.1880E-01  1.2272E+00  9.7307E-01  8.0143E-01  1.0526E+00  1.1096E+00  1.0413E+00  2.9042E+00  4.2579E-01  1.4630E+00
             1.4391E+00
 PARAMETER:  1.5318E-02  3.0473E-01  7.2698E-02 -1.2136E-01  1.5127E-01  2.0396E-01  1.4047E-01  1.1661E+00 -7.5381E-01  4.8050E-01
             4.6402E-01
 GRADIENT:   1.0993E+02  8.9425E+01  1.1818E+00 -8.9045E+01 -1.1145E+02  7.8222E+01 -1.1622E+01 -6.2821E+01 -5.2576E+00 -1.8587E+00
             4.5599E+02

0ITERATION NO.:   34    OBJECTIVE VALUE:  -3534.50144469311        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1149
 NPARAMETR:  9.1886E-01  1.2273E+00  9.7309E-01  8.0150E-01  1.0527E+00  1.1512E+00  1.0413E+00  2.9107E+00  4.2574E-01  1.4631E+00
             1.4382E+00
 PARAMETER:  1.5356E-02  3.0485E-01  7.2702E-02 -1.2130E-01  1.5137E-01  2.4328E-01  1.4048E-01  1.1687E+00 -7.5376E-01  4.8053E-01
             4.6336E-01
 GRADIENT:  -1.3323E+04  8.5559E+03 -1.3207E+04 -1.1021E+04 -1.4993E+02  9.5822E+00  1.8782E+04  3.5645E+03  3.4950E+03 -2.0848E+01
             4.2892E+02
 NUMSIGDIG:         2.3         2.3         2.3         2.3         4.7         0.7         2.3         2.3         2.3         5.0
                    3.7

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1149
 NO. OF SIG. DIGITS IN FINAL EST.:  0.7
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0011E-02  4.0426E-02 -3.5199E-02  3.1474E-02  7.1598E-02
 SE:             2.8487E-02  2.8398E-02  3.3013E-02  1.5416E-02  2.6537E-02
 N:                     100         100         100         100         100

 P VAL.:         3.5153E-02  1.5458E-01  2.8634E-01  4.1196E-02  6.9739E-03

 ETASHRINKSD(%)  4.5646E+00  4.8623E+00  1.0000E-10  4.8353E+01  1.1099E+01
 ETASHRINKVR(%)  8.9208E+00  9.4882E+00  1.0000E-10  7.3326E+01  2.0966E+01
 EBVSHRINKSD(%)  3.8824E-01  9.7864E+00  2.5265E+01  5.4774E+01  1.1104E+01
 EBVSHRINKVR(%)  7.7497E-01  1.8615E+01  4.4147E+01  7.9546E+01  2.0975E+01
 RELATIVEINF(%)  9.9218E+01  2.6604E+01  4.8412E+01  6.0329E+00  4.7188E+01
 EPSSHRINKSD(%)  3.8010E+01
 EPSSHRINKVR(%)  6.1573E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3534.5014446931118     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1880.4120849247010     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.10
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.46
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3534.501       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.19E-01  1.23E+00  9.73E-01  8.01E-01  1.05E+00  1.15E+00  1.04E+00  2.91E+00  4.26E-01  1.46E+00  1.44E+00
 


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
+        7.82E+06
 
 TH 2
+        1.92E+06  9.44E+05
 
 TH 3
+       -9.77E+01 -1.82E+06  6.97E+06
 
 TH 4
+       -1.09E+02 -1.82E+06  7.64E+02  1.16E+07
 
 TH 5
+       -4.51E+06 -1.11E+06  4.26E+06  7.09E+06  6.93E+06
 
 TH 6
+       -2.20E+01  2.97E+01 -1.23E+02 -1.26E+02 -9.79E-01  1.40E+02
 
 TH 7
+        6.58E+01 -4.27E+02  2.02E+03  2.02E+03 -2.83E+06  8.14E+01  3.09E+06
 
 TH 8
+        2.15E+05 -1.90E+04  7.29E+04  3.33E+05  2.03E+05  6.59E+00 -1.35E+05  9.19E+03
 
 TH 9
+        3.02E+01  5.50E+05  4.51E+02 -2.91E+02 -1.29E+06  3.71E+01 -5.62E+02 -1.01E+05  6.41E+05
 
 TH10
+       -5.84E-01 -1.17E+01  1.32E-01  9.13E+01  7.90E+01 -1.81E-01 -2.65E+00 -1.83E+00  2.17E+01  2.67E+05
 
 TH11
+       -5.93E+00 -9.05E+00 -1.89E+00  4.12E+01  6.69E+00  8.36E-01  6.52E+00 -1.90E+04  1.29E+01 -1.41E+05  1.50E+05
 
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
 #CPUT: Total CPU Time in Seconds,       52.668
Stop Time:
Wed Sep 29 05:59:36 CDT 2021
