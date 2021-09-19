Sat Sep 18 08:28:50 CDT 2021
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
$DATA ../../../../data/spa/B/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1708.21332936770        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9337E+01 -1.3220E+01 -2.5424E+01  1.7212E+01  1.9128E+01  5.7361E+01 -6.5797E+00  6.8119E+00 -2.4391E+00  7.5953E+00
             1.2507E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1713.58209437259        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0202E+00  1.0579E+00  1.0585E+00  9.5240E-01  1.0650E+00  8.1152E-01  1.0851E+00  9.5310E-01  1.0268E+00  9.6341E-01
             9.7935E-01
 PARAMETER:  1.1998E-01  1.5632E-01  1.5688E-01  5.1228E-02  1.6294E-01 -1.0885E-01  1.8169E-01  5.1966E-02  1.2647E-01  6.2721E-02
             7.9136E-02
 GRADIENT:   8.6108E+01 -1.7035E+01 -9.7169E+00 -2.0360E-01  3.8050E+01 -1.6168E+01 -1.7636E+00  3.9853E-01  3.3514E+00 -6.8522E+00
             2.7744E-02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1714.75869501316        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0134E+00  1.0758E+00  8.1533E-01  9.2830E-01  9.4324E-01  8.2316E-01  1.2777E+00  5.4287E-01  9.9005E-01  8.6883E-01
             9.7489E-01
 PARAMETER:  1.1328E-01  1.7305E-01 -1.0416E-01  2.5598E-02  4.1569E-02 -9.4611E-02  3.4509E-01 -5.1088E-01  8.9995E-02 -4.0603E-02
             7.4567E-02
 GRADIENT:   5.4783E+01 -5.1084E-01 -1.0470E+01  2.0923E+00  1.8840E+01 -1.0778E+01  1.4407E+01  1.5386E+00  6.1549E+00 -4.8727E-01
             1.2818E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1716.06371109399        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0148E+00  1.2091E+00  7.5938E-01  8.5405E-01  9.6928E-01  8.5837E-01  1.0623E+00  4.3282E-01  1.0439E+00  9.0857E-01
             9.7557E-01
 PARAMETER:  1.1469E-01  2.8988E-01 -1.7526E-01 -5.7771E-02  6.8795E-02 -5.2725E-02  1.6040E-01 -7.3744E-01  1.4299E-01  4.1149E-03
             7.5264E-02
 GRADIENT:   2.8742E+00  7.7967E-01 -1.1833E+00  3.4342E+00 -6.5002E-01  1.9792E+00 -3.3236E-01  5.8844E-01 -6.4049E-01  1.3652E+00
             8.0595E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1716.66342132368        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      473
 NPARAMETR:  1.0144E+00  1.4891E+00  6.2720E-01  6.7136E-01  1.0616E+00  8.5697E-01  8.9970E-01  1.1840E-01  1.2518E+00  9.4997E-01
             9.7777E-01
 PARAMETER:  1.1429E-01  4.9817E-01 -3.6649E-01 -2.9845E-01  1.5981E-01 -5.4352E-02 -5.6886E-03 -2.0337E+00  3.2461E-01  4.8675E-02
             7.7517E-02
 GRADIENT:  -5.6657E-01 -6.7320E-01 -1.2498E+00  1.9839E+00  2.2160E+00  1.2000E+00  1.1550E-01  4.2384E-02 -2.1394E-01  3.5085E-02
             1.2899E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1716.77744678251        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  1.0147E+00  1.6618E+00  5.5704E-01  5.5931E-01  1.1335E+00  8.5458E-01  8.2414E-01  3.1290E-02  1.4335E+00  9.9500E-01
             9.7674E-01
 PARAMETER:  1.1460E-01  6.0790E-01 -4.8513E-01 -4.8105E-01  2.2534E-01 -5.7148E-02 -9.3411E-02 -3.3645E+00  4.6014E-01  9.4988E-02
             7.6464E-02
 GRADIENT:   1.3523E-01  2.6229E+00  5.2557E-01  6.8054E-01 -1.2038E+00  1.2564E-01 -1.4617E-01  2.4641E-03 -2.4949E-02 -9.4831E-02
            -8.2559E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1716.78002198864        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      828
 NPARAMETR:  1.0147E+00  1.6756E+00  5.4898E-01  5.4888E-01  1.1405E+00  8.5432E-01  8.1934E-01  2.0993E-02  1.4509E+00  9.9903E-01
             9.7669E-01
 PARAMETER:  1.1459E-01  6.1619E-01 -4.9969E-01 -4.9987E-01  2.3143E-01 -5.7454E-02 -9.9250E-02 -3.7636E+00  4.7217E-01  9.9034E-02
             7.6419E-02
 GRADIENT:   8.6161E-02 -2.5572E-01 -8.3958E-02  1.2557E-02  2.3484E-01  7.9670E-03 -1.0178E-02  1.1377E-03  2.3432E-02 -3.9211E-03
             1.4648E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1716.78049097074        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1004
 NPARAMETR:  1.0147E+00  1.6756E+00  5.4897E-01  5.4898E-01  1.1402E+00  8.5430E-01  8.1939E-01  1.0000E-02  1.4505E+00  9.9883E-01
             9.7666E-01
 PARAMETER:  1.1456E-01  6.1619E-01 -4.9972E-01 -4.9970E-01  2.3119E-01 -5.7476E-02 -9.9193E-02 -4.6790E+00  4.7191E-01  9.8825E-02
             7.6383E-02
 GRADIENT:  -1.4552E-03  4.7295E-02 -2.6224E-03  1.6831E-02  1.2712E-02 -1.1958E-03 -1.5844E-02  0.0000E+00 -4.9908E-03 -1.0852E-02
            -7.8187E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1716.78049097074        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  1.0147E+00  1.6756E+00  5.4897E-01  5.4898E-01  1.1402E+00  8.5430E-01  8.1939E-01  1.0000E-02  1.4505E+00  9.9883E-01
             9.7666E-01
 PARAMETER:  1.1456E-01  6.1619E-01 -4.9972E-01 -4.9970E-01  2.3119E-01 -5.7476E-02 -9.9193E-02 -4.6790E+00  4.7191E-01  9.8825E-02
             7.6383E-02
 GRADIENT:  -1.4552E-03  4.7295E-02 -2.6224E-03  1.6831E-02  1.2712E-02 -1.1958E-03 -1.5844E-02  0.0000E+00 -4.9908E-03 -1.0852E-02
            -7.8187E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1026
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9448E-05 -3.0627E-02 -2.9104E-04  2.6830E-02 -3.6964E-02
 SE:             2.9806E-02  2.3643E-02  1.0272E-04  2.2902E-02  2.2489E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9948E-01  1.9519E-01  4.6073E-03  2.4139E-01  1.0024E-01

 ETASHRINKSD(%)  1.4637E-01  2.0793E+01  9.9656E+01  2.3276E+01  2.4660E+01
 ETASHRINKVR(%)  2.9252E-01  3.7262E+01  9.9999E+01  4.1134E+01  4.3240E+01
 EBVSHRINKSD(%)  5.5473E-01  2.0097E+01  9.9702E+01  2.5108E+01  2.2872E+01
 EBVSHRINKVR(%)  1.1064E+00  3.6155E+01  9.9999E+01  4.3912E+01  4.0513E+01
 RELATIVEINF(%)  9.8802E+01  4.3946E+00  1.0552E-04  3.8618E+00  1.2851E+01
 EPSSHRINKSD(%)  4.4560E+01
 EPSSHRINKVR(%)  6.9264E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1716.7804909707374     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -981.62966440699927     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1716.780       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.68E+00  5.49E-01  5.49E-01  1.14E+00  8.54E-01  8.19E-01  1.00E-02  1.45E+00  9.99E-01  9.77E-01
 


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
+        1.46E+03
 
 TH 2
+       -7.07E+00  3.98E+02
 
 TH 3
+        1.28E+01  1.54E+02  3.69E+02
 
 TH 4
+       -1.99E+01  3.30E+02 -2.94E+02  1.00E+03
 
 TH 5
+       -5.22E+00 -1.78E+02 -3.07E+02  2.64E+02  5.12E+02
 
 TH 6
+        1.41E-01 -1.24E+00  3.81E+00 -4.09E+00 -6.53E-01  2.65E+02
 
 TH 7
+        1.92E+00  9.05E+00 -1.44E+01 -1.37E+01 -1.05E+01  4.99E+00  1.27E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.91E+00 -1.92E+01 -3.13E+01  5.75E+01  2.49E-01 -4.06E-01  1.82E+01  0.00E+00  3.96E+01
 
 TH10
+        2.20E-01 -1.50E+01 -3.59E+01 -4.02E+00 -6.02E+01 -4.50E+00  1.46E+01  0.00E+00  4.95E+00  8.08E+01
 
 TH11
+       -1.16E+01 -1.75E+01 -2.69E+01  3.50E+00 -1.87E-01 -1.69E+00  1.63E+01  0.00E+00  3.47E+00  1.43E+01  2.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.339
Stop Time:
Sat Sep 18 08:29:10 CDT 2021
