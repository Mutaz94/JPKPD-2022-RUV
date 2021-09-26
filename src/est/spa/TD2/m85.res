Sat Sep 25 13:50:24 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat85.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1640.68569416267        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.5805E+01 -4.9114E+01 -1.1394E+01 -6.1860E+01 -5.5272E+00  3.2255E+01 -1.4977E+01  9.9572E+00 -1.4857E+01  9.8764E+00
            -3.0429E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1646.72738022641        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8265E-01  1.1169E+00  1.1171E+00  9.8843E-01  1.1152E+00  8.9442E-01  1.2027E+00  8.8127E-01  1.0801E+00  9.1341E-01
             1.0905E+00
 PARAMETER:  8.2500E-02  2.1060E-01  2.1072E-01  8.8362E-02  2.0906E-01 -1.1579E-02  2.8461E-01 -2.6394E-02  1.7707E-01  9.4332E-03
             1.8661E-01
 GRADIENT:   5.3862E+01  1.7941E+01  1.1577E+01  1.7889E+01  2.6060E+01 -8.3814E+00  4.6409E+00 -1.7028E+00  6.7590E+00 -1.3742E+01
            -1.2491E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1649.25056242045        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.7955E-01  9.1139E-01  8.4754E-01  1.1016E+00  9.0656E-01  9.2248E-01  1.5883E+00  2.9547E-01  9.2108E-01  8.4587E-01
             1.0924E+00
 PARAMETER:  7.9341E-02  7.2176E-03 -6.5412E-02  1.9679E-01  1.9020E-03  1.9309E-02  5.6266E-01 -1.1192E+00  1.7796E-02 -6.7395E-02
             1.8835E-01
 GRADIENT:   4.0848E+01  3.9201E+00 -3.9875E+01  5.3497E+01  5.8218E+01  3.6784E+00  1.6361E+01  8.8901E-01  4.9470E+00  6.5612E-01
             9.9323E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1651.00251070157        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.6599E-01  9.9259E-01  8.1559E-01  1.0296E+00  9.0192E-01  9.1728E-01  1.3466E+00  3.1147E-01  9.5143E-01  8.2479E-01
             1.0734E+00
 PARAMETER:  6.5400E-02  9.2566E-02 -1.0385E-01  1.2916E-01 -3.2267E-03  1.3654E-02  3.9756E-01 -1.0664E+00  5.0206E-02 -9.2631E-02
             1.7081E-01
 GRADIENT:   4.7987E+00 -1.6820E+00 -6.5597E+00  6.9409E+00  8.7662E+00  1.1917E+00  1.9392E+00  8.2272E-01  5.4231E-01  8.9942E-01
             1.7899E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1651.02094886479        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.6435E-01  1.0163E+00  7.9235E-01  1.0108E+00  8.9589E-01  9.1606E-01  1.3076E+00  2.4965E-01  9.6062E-01  8.1430E-01
             1.0704E+00
 PARAMETER:  6.3704E-02  1.1612E-01 -1.3275E-01  1.1073E-01 -9.9357E-03  1.2325E-02  3.6821E-01 -1.2877E+00  5.9821E-02 -1.0542E-01
             1.6803E-01
 GRADIENT:  -1.7939E-01 -1.2693E+00 -2.1429E+00  2.5931E-01  1.1859E+00  4.6081E-01  1.6835E-01  5.3113E-01 -8.7690E-02  7.9147E-01
             5.6543E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1651.11126416269        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.6359E-01  1.0293E+00  7.7018E-01  9.9910E-01  8.8745E-01  9.1470E-01  1.2892E+00  9.6276E-02  9.6621E-01  8.0322E-01
             1.0682E+00
 PARAMETER:  6.2913E-02  1.2891E-01 -1.6113E-01  9.9097E-02 -1.9399E-02  1.0844E-02  3.5403E-01 -2.2405E+00  6.5624E-02 -1.1913E-01
             1.6598E-01
 GRADIENT:  -2.7555E+00 -8.3318E-01  9.1046E-01 -3.6557E+00 -3.0573E+00 -2.7419E-01 -7.6159E-01  7.7121E-02 -3.1567E-01  3.6317E-01
            -3.6455E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1651.51633691722        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.7819E-01  1.0044E+00  7.9168E-01  1.0227E+00  8.8931E-01  9.2400E-01  1.3374E+00  1.0000E-02  9.5855E-01  8.1552E-01
             1.0708E+00
 PARAMETER:  7.7952E-02  1.0443E-01 -1.3360E-01  1.2246E-01 -1.7309E-02  2.0960E-02  3.9074E-01 -4.7676E+00  5.7668E-02 -1.0393E-01
             1.6839E-01
 GRADIENT:   2.3806E+00  6.2296E-01  6.2239E-01 -2.5781E-01 -5.4980E-01  3.9879E-01  1.4773E-01  0.0000E+00  3.9082E-01 -8.7521E-02
             2.7598E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1651.52744877974        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.7711E-01  9.6766E-01  7.9644E-01  1.0446E+00  8.7553E-01  9.2284E-01  1.3786E+00  1.0000E-02  9.3951E-01  8.1139E-01
             1.0706E+00
 PARAMETER:  7.6849E-02  6.7127E-02 -1.2760E-01  1.4360E-01 -3.2929E-02  1.9698E-02  4.2109E-01 -4.9737E+00  3.7606E-02 -1.0901E-01
             1.6825E-01
 GRADIENT:   8.3288E-03 -3.6962E-03 -3.3156E-03 -1.7222E-03  7.5852E-03  1.3223E-03  1.7424E-03  0.0000E+00  2.0631E-04 -3.2885E-04
             4.1333E-04

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1651.52744877974        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  9.7711E-01  9.6766E-01  7.9644E-01  1.0446E+00  8.7553E-01  9.2284E-01  1.3786E+00  1.0000E-02  9.3951E-01  8.1139E-01
             1.0706E+00
 PARAMETER:  7.6849E-02  6.7127E-02 -1.2760E-01  1.4360E-01 -3.2929E-02  1.9698E-02  4.2109E-01 -4.9737E+00  3.7606E-02 -1.0901E-01
             1.6825E-01
 GRADIENT:   8.3288E-03 -3.6962E-03 -3.3156E-03 -1.7222E-03  7.5852E-03  1.3223E-03  1.7424E-03  0.0000E+00  2.0631E-04 -3.2885E-04
             4.1333E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      725
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.5738E-05  2.8252E-03 -4.5103E-04 -6.3610E-03 -1.2837E-02
 SE:             2.9798E-02  2.2253E-02  1.8398E-04  2.4089E-02  2.1607E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9878E-01  8.9897E-01  1.4228E-02  7.9173E-01  5.5242E-01

 ETASHRINKSD(%)  1.7420E-01  2.5450E+01  9.9384E+01  1.9298E+01  2.7615E+01
 ETASHRINKVR(%)  3.4811E-01  4.4423E+01  9.9996E+01  3.4872E+01  4.7604E+01
 EBVSHRINKSD(%)  5.6074E-01  2.5012E+01  9.9423E+01  1.9308E+01  2.6724E+01
 EBVSHRINKVR(%)  1.1183E+00  4.3767E+01  9.9997E+01  3.4889E+01  4.6307E+01
 RELATIVEINF(%)  9.8517E+01  4.7057E+00  4.3380E-04  6.4746E+00  5.4257E+00
 EPSSHRINKSD(%)  4.2647E+01
 EPSSHRINKVR(%)  6.7106E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1651.5274487797410     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -916.37662221600283     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1651.527       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.77E-01  9.68E-01  7.96E-01  1.04E+00  8.76E-01  9.23E-01  1.38E+00  1.00E-02  9.40E-01  8.11E-01  1.07E+00
 


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
+       -9.37E+00  3.68E+02
 
 TH 3
+        1.83E+01  2.13E+02  6.47E+02
 
 TH 4
+       -1.31E+01  2.67E+02 -2.81E+02  7.64E+02
 
 TH 5
+       -6.85E+00 -3.38E+02 -7.44E+02  3.26E+02  1.18E+03
 
 TH 6
+       -6.16E-01 -1.35E+00  1.71E+00 -2.93E+00 -1.86E+00  2.28E+02
 
 TH 7
+        1.31E+00  3.09E+01 -8.86E+00 -1.23E+01  3.42E-01  4.47E-01  3.92E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.15E+00 -2.05E+01 -3.23E+01  2.96E+01  2.05E+00 -1.12E+00  1.25E+01  0.00E+00  1.05E+02
 
 TH10
+       -3.54E+00 -7.84E+00 -6.03E+01 -2.75E+01 -5.61E+01 -8.73E-01  1.13E+01  0.00E+00  1.05E+01  9.25E+01
 
 TH11
+       -9.22E+00 -1.29E+01 -3.99E+01 -5.09E+00  7.31E+00  2.48E+00  4.52E+00  0.00E+00  1.16E+01  2.66E+01  1.93E+02
 
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
 #CPUT: Total CPU Time in Seconds,       13.065
Stop Time:
Sat Sep 25 13:50:39 CDT 2021
