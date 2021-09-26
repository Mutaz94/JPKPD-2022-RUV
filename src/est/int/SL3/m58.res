Sat Sep 25 02:27:56 CDT 2021
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
$DATA ../../../../data/int/SL3/dat58.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      982
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

 TOT. NO. OF OBS RECS:      882
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
 RAW OUTPUT FILE (FILE): m58.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -399.819070036282        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7251E+02  4.3701E+01  7.8183E+01  7.5607E+01  4.6116E+01  2.0127E+01 -8.6557E+01 -2.0882E+02 -1.2048E+02 -6.3177E+01
            -6.1830E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2611.76462778456        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.4431E-01  1.3016E+00  1.1042E+00  8.6136E-01  1.2337E+00  7.9433E-01  1.0823E+00  9.0985E-01  9.2381E-01  1.2311E+00
             2.7468E+00
 PARAMETER:  4.2701E-02  3.6361E-01  1.9914E-01 -4.9240E-02  3.1003E-01 -1.3025E-01  1.7907E-01  5.5276E-03  2.0750E-02  3.0789E-01
             1.1104E+00
 GRADIENT:  -4.9349E+01  1.5003E+01 -1.1835E+01  2.8525E+01 -9.0220E+00 -6.0175E+01  1.7390E+01  4.2132E+00 -4.2815E+00 -1.4126E+01
             5.0477E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2620.78248775526        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.6381E-01  1.4062E+00  1.9813E+00  8.3685E-01  1.5507E+00  8.2516E-01  1.0150E+00  6.5221E-01  9.2527E-01  1.4780E+00
             2.7657E+00
 PARAMETER:  6.3143E-02  4.4086E-01  7.8376E-01 -7.8107E-02  5.3868E-01 -9.2179E-02  1.1484E-01 -3.2738E-01  2.2327E-02  4.9067E-01
             1.1173E+00
 GRADIENT:   1.8449E+01  4.6824E+01 -1.7927E+00  5.8902E+01 -3.4441E+00 -4.1772E+01  1.0599E+01 -7.4838E-01 -2.3596E+00 -1.6292E+01
             7.1974E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2639.16823229826        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.5549E-01  1.5448E+00  4.9638E+00  7.1042E-01  2.0431E+00  9.1012E-01  9.3644E-01  2.3893E+00  7.9223E-01  1.9578E+00
             2.6140E+00
 PARAMETER:  5.4474E-02  5.3492E-01  1.7022E+00 -2.4190E-01  8.1445E-01  5.8265E-03  3.4329E-02  9.7102E-01 -1.3291E-01  7.7184E-01
             1.0609E+00
 GRADIENT:   2.5559E+00 -6.7897E+00 -5.5489E+00  5.4212E+00  1.2368E+01 -2.9269E+00 -4.5217E+00 -2.6067E-02  1.1483E+00  1.1127E+00
            -1.3827E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2646.92860547785        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.5408E-01  1.3318E+00  6.6015E+01  8.5927E-01  2.2272E+00  9.1750E-01  1.2602E+00  5.4538E+00  1.2451E-01  2.0551E+00
             2.6168E+00
 PARAMETER:  5.2997E-02  3.8650E-01  4.2899E+00 -5.1673E-02  9.0077E-01  1.3892E-02  3.3127E-01  1.7963E+00 -1.9833E+00  8.2033E-01
             1.0620E+00
 GRADIENT:  -2.4857E+00  1.3779E+01 -4.5963E-01  5.7788E+00  4.5306E+00 -2.7678E-02  3.0631E+00 -3.3709E-01  2.8296E-01 -4.7842E+00
            -4.5672E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2647.76607489352        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.5509E-01  1.1101E+00  3.8808E+02  9.9401E-01  2.2265E+00  9.1528E-01  1.4811E+00  1.0735E+01  3.1781E-02  2.0944E+00
             2.6189E+00
 PARAMETER:  5.4054E-02  2.0448E-01  6.0612E+00  9.3988E-02  9.0042E-01  1.1470E-02  4.9276E-01  2.4735E+00 -3.3489E+00  8.3928E-01
             1.0628E+00
 GRADIENT:   5.6267E-01 -2.7159E+00 -1.0219E-01 -3.1314E+00  6.4540E-01 -4.3086E-01  1.1564E+00 -6.9178E-02  1.4694E-02 -1.4261E-01
            -1.6594E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2647.97481942848        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  9.5683E-01  1.1703E+00  5.2310E+03  9.5944E-01  2.2480E+00  9.1845E-01  1.4141E+00  4.7785E+01  1.0000E-02  2.1095E+00
             2.6217E+00
 PARAMETER:  5.5869E-02  2.5723E-01  8.6624E+00  5.8591E-02  9.1005E-01  1.4936E-02  4.4649E-01  3.9667E+00 -5.2197E+00  8.4644E-01
             1.0638E+00
 GRADIENT:  -1.3058E+00  1.5571E-01 -6.3070E-03  1.0577E-01  2.2897E-01 -2.6195E-02  9.0303E-02 -5.8411E-03  0.0000E+00  5.8825E-02
            -4.2909E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2647.98366791619        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      682             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5724E-01  1.1714E+00  7.9059E+04  9.5894E-01  2.2469E+00  9.1835E-01  1.4112E+00  2.2296E+02  1.0000E-02  2.1090E+00
             2.6223E+00
 PARAMETER:  5.6298E-02  2.5818E-01  1.1378E+01  5.8069E-02  9.0954E-01  1.4826E-02  4.4442E-01  5.5070E+00 -7.2627E+00  8.4620E-01
             1.0640E+00
 GRADIENT:   3.9475E+00  6.8035E+00 -1.6045E-04  7.0542E+00  2.8512E+00  4.5303E-01  2.6936E-01 -5.1476E-04  0.0000E+00  1.3046E+00
             1.4991E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2647.98394516410        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  9.5737E-01  1.1712E+00  7.8323E+04  9.5867E-01  2.2471E+00  9.1836E-01  1.4122E+00  2.2332E+02  1.0000E-02  2.1086E+00
             2.6223E+00
 PARAMETER:  5.6434E-02  2.5801E-01  1.1369E+01  5.7790E-02  9.0962E-01  1.4830E-02  4.4517E-01  5.5086E+00 -7.2627E+00  8.4603E-01
             1.0640E+00
 GRADIENT:  -3.7592E+00  4.7234E-01 -4.3075E-04  5.5616E-01 -7.1928E-02 -5.7298E-02 -2.3429E-02 -3.3852E-04  0.0000E+00  4.8687E-03
             1.3610E-01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -2647.98394516410        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:      895
 NPARAMETR:  9.5738E-01  1.1712E+00  7.8257E+04  9.5868E-01  2.2471E+00  9.1835E-01  1.4123E+00  2.2209E+02  1.0000E-02  2.1087E+00
             2.6223E+00
 PARAMETER:  5.6434E-02  2.5801E-01  1.1369E+01  5.7790E-02  9.0962E-01  1.4830E-02  4.4517E-01  5.5086E+00 -7.2627E+00  8.4603E-01
             1.0640E+00
 GRADIENT:  -2.1015E-01 -1.6744E-01  9.4167E-04 -2.6911E-01 -8.2427E-04  1.4646E-01 -1.9472E-02  1.3732E-03  0.0000E+00 -6.1499E-02
             1.0524E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      895
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3029E-03 -2.3768E-03  3.9644E-06 -9.9161E-04 -1.4189E-02
 SE:             2.9351E-02  2.7948E-02  7.2809E-05  3.0541E-04  2.6692E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3746E-01  9.3222E-01  9.5658E-01  1.1673E-03  5.9503E-01

 ETASHRINKSD(%)  1.6719E+00  6.3723E+00  9.9756E+01  9.8977E+01  1.0577E+01
 ETASHRINKVR(%)  3.3159E+00  1.2338E+01  9.9999E+01  9.9990E+01  2.0035E+01
 EBVSHRINKSD(%)  1.8490E+00  5.9061E+00  9.9819E+01  9.9110E+01  7.8085E+00
 EBVSHRINKVR(%)  3.6638E+00  1.1463E+01  1.0000E+02  9.9992E+01  1.5007E+01
 RELATIVEINF(%)  9.6194E+01  1.2019E+01  2.1043E-04  1.0663E-03  5.5406E+01
 EPSSHRINKSD(%)  1.5326E+01
 EPSSHRINKVR(%)  2.8303E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          882
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1621.0075725730426     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2647.9839451640960     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1026.9763725910534     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2647.984       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.17E+00  7.83E+04  9.59E-01  2.25E+00  9.18E-01  1.41E+00  2.23E+02  1.00E-02  2.11E+00  2.62E+00
 


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
+        2.59E+03
 
 TH 2
+       -4.62E+02  4.34E+02
 
 TH 3
+       -2.65E-05 -2.09E-05 -4.19E-13
 
 TH 4
+       -9.13E+02  4.11E+02 -5.58E-05  1.67E+03
 
 TH 5
+        2.79E+01 -1.50E+01  2.63E-06  2.14E+01  4.57E+01
 
 TH 6
+       -1.26E+03 -4.79E+01 -1.38E-05  3.00E+02  3.37E+01  1.83E+02
 
 TH 7
+        8.04E+01  5.95E+01 -4.81E-06  8.60E+00 -2.08E+00  5.35E+01  9.13E+01
 
 TH 8
+        7.11E-02 -3.21E-03  3.27E-09  8.68E-02  2.29E-03  8.92E-02  3.56E-03 -5.10E-07
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        4.61E+01  6.36E+00 -7.51E-08 -4.28E+01 -3.54E+00 -1.25E+01  2.39E-01 -1.62E-04  0.00E+00  2.96E+01
 
 TH11
+       -1.10E+01 -1.69E+01  5.67E-07 -8.41E+01  1.71E+00  3.15E+01  6.50E+00 -4.98E-04  0.00E+00  4.82E+00  1.72E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       31.424
Stop Time:
Sat Sep 25 02:28:29 CDT 2021
