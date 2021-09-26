Sat Sep 25 02:54:27 CDT 2021
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
$DATA ../../../../data/int/SL3/dat98.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      984
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

 TOT. NO. OF OBS RECS:      884
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   118.315616408969        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0090E+02 -4.1670E+01 -2.6294E+01  1.8105E+02  1.5733E+02  2.0297E+01 -1.5044E+02 -1.6320E+02 -1.7433E+02 -6.4555E+01
            -7.2952E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2311.83074209805        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0339E+00  1.3332E+00  1.0994E+00  8.7499E-01  1.0731E+00  8.4703E-01  1.1479E+00  1.0224E+00  9.7210E-01  1.0276E+00
             5.2077E+00
 PARAMETER:  1.3332E-01  3.8755E-01  1.9476E-01 -3.3547E-02  1.7053E-01 -6.6021E-02  2.3792E-01  1.2212E-01  7.1708E-02  1.2723E-01
             1.7501E+00
 GRADIENT:  -5.7569E+00 -3.5805E+00  6.7387E-01 -2.6328E+01 -3.3918E+01 -2.1602E+01  2.1373E+01  3.8919E+00  1.8276E+01  1.1413E+01
             8.0763E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2368.65089264694        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.7981E-01  1.8471E+00  9.4337E+00  7.1041E-01  2.3500E+00  7.5166E-01  1.2301E+00  1.7996E+00  1.3087E+00  1.6768E+00
             4.5522E+00
 PARAMETER:  7.9607E-02  7.1363E-01  2.3443E+00 -2.4191E-01  9.5443E-01 -1.8548E-01  3.0712E-01  6.8757E-01  3.6904E-01  6.1689E-01
             1.6156E+00
 GRADIENT:  -1.6964E+02  9.0140E+01 -5.0279E+00  4.6390E+01  1.2397E+02 -8.9245E+01  5.5530E+01 -3.8049E-01  1.9431E+01 -1.4794E+01
             6.9270E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2613.19050164463        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.8599E-01  1.4998E+00  4.5767E+00  7.6183E-01  1.6710E+00  9.5658E-01  9.2290E-01  3.0413E+00  1.0458E+00  1.5728E+00
             2.8331E+00
 PARAMETER:  8.5887E-02  5.0531E-01  1.6210E+00 -1.7204E-01  6.1343E-01  5.5613E-02  1.9762E-02  1.2123E+00  1.4481E-01  5.5283E-01
             1.1414E+00
 GRADIENT:  -1.0958E+01  1.2563E+01  2.4201E+00  5.5515E+00  1.2790E+01  7.5852E+00  1.8082E+00 -7.7823E+00  3.6630E+00 -3.3063E-01
             1.7380E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2614.60964905711        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.8950E-01  1.4077E+00  4.4621E+00  8.1489E-01  1.5841E+00  9.3596E-01  1.0007E+00  3.5640E+00  8.6712E-01  1.4984E+00
             2.8150E+00
 PARAMETER:  8.9442E-02  4.4197E-01  1.5956E+00 -1.0470E-01  5.6005E-01  3.3821E-02  1.0066E-01  1.3709E+00 -4.2583E-02  5.0441E-01
             1.1350E+00
 GRADIENT:  -2.1552E+00  3.7270E+00  1.8904E-01 -4.0273E-01 -7.6432E-01  2.7900E-02 -1.9423E+00 -4.0299E-01 -5.5020E-01  1.6638E-02
             5.7030E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2614.67271638312        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      456
 NPARAMETR:  9.9270E-01  1.4230E+00  4.6164E+00  8.0679E-01  1.6068E+00  9.3738E-01  1.0069E+00  3.6466E+00  8.7175E-01  1.5187E+00
             2.8097E+00
 PARAMETER:  9.2669E-02  4.5278E-01  1.6296E+00 -1.1470E-01  5.7424E-01  3.5331E-02  1.0691E-01  1.3938E+00 -3.7249E-02  5.1788E-01
             1.1331E+00
 GRADIENT:   1.5367E-01 -3.0995E-01 -1.1294E-01  1.3851E-01  5.1988E-01 -5.0601E-02 -1.8234E-02 -2.1473E-03  7.6782E-02  4.1834E-02
            -2.7195E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2614.67362067184        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      632            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9265E-01  1.4339E+00  4.5843E+00  7.9948E-01  1.6086E+00  9.3750E-01  1.0017E+00  3.6572E+00  8.7547E-01  1.5217E+00
             2.8101E+00
 PARAMETER:  9.2620E-02  4.6039E-01  1.6226E+00 -1.2379E-01  5.7534E-01  3.5465E-02  1.0166E-01  1.3967E+00 -3.2997E-02  5.1981E-01
             1.1332E+00
 GRADIENT:   5.6382E+00  6.3719E+00  3.0481E-01  1.0470E+00  2.2757E+00  5.7678E-01  1.1411E-01  3.2578E-01  5.0646E-02  4.3850E-01
             2.0002E+00

0ITERATION NO.:   33    OBJECTIVE VALUE:  -2614.67362112129        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  9.9265E-01  1.4339E+00  4.5842E+00  7.9948E-01  1.6086E+00  9.3751E-01  1.0017E+00  3.6572E+00  8.7556E-01  1.5217E+00
             2.8101E+00
 PARAMETER:  9.2620E-02  4.6039E-01  1.6226E+00 -1.2379E-01  5.7535E-01  3.5468E-02  1.0170E-01  1.3967E+00 -3.2897E-02  5.1981E-01
             1.1332E+00
 GRADIENT:   2.1984E-03 -1.5857E-03 -3.0380E-04 -1.0608E-03 -3.6125E-04 -1.2367E-03  1.8199E-03  6.4233E-04  1.5704E-03  1.5993E-04
             1.2319E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      725
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6858E-03 -1.2631E-02 -3.6871E-02  3.7065E-03 -2.9388E-02
 SE:             2.9283E-02  2.3113E-02  1.5031E-02  1.7744E-02  2.4091E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5409E-01  5.8473E-01  1.4165E-02  8.3454E-01  2.2252E-01

 ETASHRINKSD(%)  1.8970E+00  2.2567E+01  4.9645E+01  4.0554E+01  1.9292E+01
 ETASHRINKVR(%)  3.7581E+00  4.0041E+01  7.4644E+01  6.4662E+01  3.4862E+01
 EBVSHRINKSD(%)  1.9825E+00  2.3051E+01  5.6217E+01  4.2863E+01  1.5284E+01
 EBVSHRINKVR(%)  3.9257E+00  4.0788E+01  8.0831E+01  6.7353E+01  2.8232E+01
 RELATIVEINF(%)  9.5984E+01  2.9997E+00  5.5693E+00  1.5843E+00  4.0181E+01
 EPSSHRINKSD(%)  1.6685E+01
 EPSSHRINKVR(%)  3.0585E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          884
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1624.6833267058612     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2614.6736211212856     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -989.99029441542439     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2614.674       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  1.43E+00  4.58E+00  7.99E-01  1.61E+00  9.38E-01  1.00E+00  3.66E+00  8.76E-01  1.52E+00  2.81E+00
 


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
+        1.24E+03
 
 TH 2
+       -1.56E+01  3.29E+02
 
 TH 3
+        3.23E-01  1.47E+00  8.62E-01
 
 TH 4
+       -2.91E+01  4.38E+02 -6.59E+00  7.23E+02
 
 TH 5
+       -4.11E+00 -3.71E+01 -6.39E+00  2.46E+01  1.31E+02
 
 TH 6
+        1.04E+01 -5.19E+00  9.43E-02 -1.51E+01 -2.05E+00  2.03E+02
 
 TH 7
+        2.61E+00  1.41E+01  1.70E-01 -3.13E+01  6.61E-01 -4.84E+00  7.94E+01
 
 TH 8
+       -1.60E-01 -2.14E+00 -5.36E-01  4.20E+00  2.14E+00 -1.55E-02  1.39E-01  2.13E+00
 
 TH 9
+        7.19E+00 -7.87E+00 -5.19E-01 -4.18E+00  2.39E+00 -7.53E+00  3.42E+01  1.48E+00  3.64E+01
 
 TH10
+        1.17E+00 -4.40E+00 -1.52E+00  1.66E+01 -1.60E+01  4.67E-01  1.72E+00  1.33E+00  5.61E-01  4.06E+01
 
 TH11
+       -1.75E+01 -1.36E+01 -1.66E-01 -1.53E+01 -2.17E+00  2.57E+00  4.61E+00  1.08E+00  3.84E+00  5.32E+00  1.43E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.044
Stop Time:
Sat Sep 25 02:54:58 CDT 2021
