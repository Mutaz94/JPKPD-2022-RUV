Sat Sep 18 15:25:07 CDT 2021
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
$DATA ../../../../data/spa/D/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12014.3949245032        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.5572E+01  1.2195E+02 -4.0256E+01  5.3235E+01  2.1942E+02 -1.4446E+03 -5.6838E+02 -2.2060E+01 -9.4299E+02 -6.0501E+02
            -2.3743E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -654.054641692566        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4863E+00  1.1555E+00  9.4993E-01  1.6747E+00  1.2411E+00  1.7981E+00  1.1775E+00  9.6576E-01  1.1637E+00  1.2131E+00
             1.4563E+01
 PARAMETER:  4.9629E-01  2.4455E-01  4.8633E-02  6.1561E-01  3.1599E-01  6.8675E-01  2.6342E-01  6.5156E-02  2.5161E-01  2.9316E-01
             2.7785E+00
 GRADIENT:   2.9257E+01  1.6195E+01 -4.2218E+00  2.7700E+01 -8.2027E+00  4.3059E+01  9.0129E-01  4.3452E+00  1.0596E+01  3.5182E+00
             1.3587E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -672.726351036808        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.4011E+00  1.1998E+00  1.1067E+00  1.5962E+00  2.9013E+00  1.4578E+00  2.3893E+00  4.4818E-01  9.9222E-01  6.9582E+00
             1.3131E+01
 PARAMETER:  4.3728E-01  2.8219E-01  2.0138E-01  5.6760E-01  1.1651E+00  4.7692E-01  9.7100E-01 -7.0256E-01  9.2193E-02  2.0399E+00
             2.6749E+00
 GRADIENT:   5.0623E+01  2.8553E+01 -7.4170E+00  3.3833E+01 -2.0852E+01 -1.4725E+01  1.0634E+01  3.4424E-01  1.0132E+01  1.6145E+01
             8.3262E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -723.364853701122        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0804E+00  2.7392E-01  7.9433E-01  1.4722E+00  8.4087E+00  1.2857E+00  1.4760E+00  1.0426E-01  2.1704E-01  1.0339E+01
             1.0740E+01
 PARAMETER:  1.7730E-01 -1.1949E+00 -1.3026E-01  4.8674E-01  2.2293E+00  3.5129E-01  4.8930E-01 -2.1609E+00 -1.4277E+00  2.4360E+00
             2.4740E+00
 GRADIENT:   2.5910E+00 -3.3146E+00  1.9604E+01 -7.2150E+01  8.3657E+00  4.0540E+00  6.9544E-01 -5.2673E-03  1.3385E+00 -1.3591E+00
             2.5962E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -727.874527734935        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0832E+00  1.8611E-01  6.8104E-01  1.5730E+00  6.8246E+00  1.3089E+00  1.6411E+00  5.9250E-02  1.4808E-01  8.6223E+00
             1.0638E+01
 PARAMETER:  1.7993E-01 -1.5814E+00 -2.8413E-01  5.5300E-01  2.0205E+00  3.6922E-01  5.9534E-01 -2.7260E+00 -1.8100E+00  2.2544E+00
             2.4644E+00
 GRADIENT:  -1.1249E+00  4.7033E+00 -2.4458E+00  1.6672E+01  1.2251E+00 -3.5641E+00  2.9124E-01 -6.1205E-03  6.0887E-01 -2.1838E+00
            -6.7637E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -731.502135095105        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0667E+00  2.9723E-02  5.8602E-01  1.5851E+00  4.2596E+00  1.3338E+00  9.3056E-01  1.0000E-02  1.8999E-02  6.4463E+00
             1.0515E+01
 PARAMETER:  1.6460E-01 -3.4158E+00 -4.3439E-01  5.6066E-01  1.5492E+00  3.8803E-01  2.8030E-02 -4.5358E+00 -3.8634E+00  1.9635E+00
             2.4528E+00
 GRADIENT:   9.7674E+00  8.0394E-01 -1.8350E+00  1.3258E+01 -1.8173E+00 -2.2330E+00  2.4824E-03  0.0000E+00  1.2057E-02 -2.5950E+00
            -8.2702E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -731.695061510690        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.0433E+00  2.7709E-02  5.6316E-01  1.5571E+00  4.3516E+00  1.3503E+00  8.6381E-01  1.0000E-02  1.7094E-02  6.5486E+00
             1.0475E+01
 PARAMETER:  1.4238E-01 -3.4860E+00 -4.7419E-01  5.4280E-01  1.5706E+00  4.0029E-01 -4.6399E-02 -4.8041E+00 -3.9690E+00  1.9792E+00
             2.4490E+00
 GRADIENT:  -1.7091E+00  7.3109E-01 -1.4894E+00  1.6579E+00 -3.8943E+00  1.3812E+00  1.9969E-03  0.0000E+00  1.1170E-02  4.4614E+00
            -4.4600E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -731.739431316790        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  1.0430E+00  2.5978E-02  5.5329E-01  1.5507E+00  4.3979E+00  1.3482E+00  8.2688E-01  1.0000E-02  1.5680E-02  6.5553E+00
             1.0469E+01
 PARAMETER:  1.4215E-01 -3.5505E+00 -4.9187E-01  5.3868E-01  1.5811E+00  3.9879E-01 -9.0090E-02 -4.8932E+00 -4.0554E+00  1.9803E+00
             2.4485E+00
 GRADIENT:  -3.7294E+00  7.1434E-01 -2.2655E+00  2.3666E+00 -1.1941E+00  2.0223E+00  1.4922E-03  0.0000E+00  9.0896E-03 -1.8450E+00
            -8.4132E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -731.739569180855        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  1.0430E+00  2.5995E-02  5.5335E-01  1.5509E+00  4.3997E+00  1.3483E+00  8.2684E-01  1.0000E-02  1.5613E-02  6.5525E+00
             1.0475E+01
 PARAMETER:  1.4215E-01 -3.5507E+00 -4.9189E-01  5.3868E-01  1.5811E+00  3.9878E-01 -9.0198E-02 -4.8933E+00 -4.0556E+00  1.9803E+00
             2.4485E+00
 GRADIENT:   4.4526E+02 -1.4646E+01 -7.1421E+01 -1.4536E+02 -5.6953E+01 -1.5483E+02 -2.3294E-02  0.0000E+00  9.0657E-03  5.8463E+01
            -4.1692E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      627
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9682E-03 -2.0828E-04 -1.6398E-04 -3.9523E-04 -2.9783E-02
 SE:             2.7970E-02  1.3408E-04  7.7654E-05  4.5087E-04  9.3260E-03
 N:                     100         100         100         100         100

 P VAL.:         7.2155E-01  1.2033E-01  3.4709E-02  3.8071E-01  1.4054E-03

 ETASHRINKSD(%)  6.2973E+00  9.9551E+01  9.9740E+01  9.8490E+01  6.8757E+01
 ETASHRINKVR(%)  1.2198E+01  9.9998E+01  9.9999E+01  9.9977E+01  9.0239E+01
 EBVSHRINKSD(%)  6.2017E+00  9.9486E+01  9.9704E+01  9.8431E+01  5.2374E+01
 EBVSHRINKVR(%)  1.2019E+01  9.9997E+01  9.9999E+01  9.9975E+01  7.7317E+01
 RELATIVEINF(%)  2.9899E+01  1.0622E-04  1.2256E-04  6.0238E-04  7.1534E+00
 EPSSHRINKSD(%)  6.0197E+00
 EPSSHRINKVR(%)  1.1677E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -731.73956918085514     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       3.4112573828830364     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.89
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -731.740       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.60E-02  5.53E-01  1.55E+00  4.40E+00  1.35E+00  8.27E-01  1.00E-02  1.57E-02  6.56E+00  1.05E+01
 


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
+        8.09E+05
 
 TH 2
+       -6.77E+04  1.71E+06
 
 TH 3
+       -1.01E+05  3.17E+04  2.62E+05
 
 TH 4
+       -1.96E+04  6.98E+03  1.69E+04  3.37E+04
 
 TH 5
+       -1.71E+03  5.96E+02  1.47E+03  1.07E+03  5.68E+02
 
 TH 6
+       -5.10E+04  2.05E+04  2.20E+04  4.31E+03  3.85E+02  5.96E+04
 
 TH 7
+        8.07E+00 -4.10E+01 -9.79E+00  3.85E+00  9.25E-01  5.31E+00 -5.39E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.90E+03 -1.08E+04 -1.44E+03 -3.05E+02 -2.63E+01 -8.83E+02  5.07E+01  0.00E+00 -5.70E+01
 
 TH10
+        5.27E+02 -1.85E+02 -4.60E+02 -3.24E+02 -7.25E+01 -1.19E+02 -5.78E-01  0.00E+00  8.13E+00  2.60E+02
 
 TH11
+       -3.75E+02  1.24E+02  3.03E+02  2.12E+02  5.14E+01  8.44E+01 -6.74E-02  0.00E+00 -5.47E+00 -3.41E+01  5.47E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.857
Stop Time:
Sat Sep 18 15:25:27 CDT 2021
