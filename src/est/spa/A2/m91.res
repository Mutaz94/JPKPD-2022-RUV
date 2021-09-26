Sat Sep 25 08:56:30 CDT 2021
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
$DATA ../../../../data/spa/A2/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1201.55450543660        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8502E+02 -4.0508E+01 -4.0795E+01 -2.0465E+01  1.2101E+02 -3.0293E+01 -1.5879E+01  1.2752E+01 -3.1361E+01 -2.7223E+01
            -7.5773E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1432.26395099814        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.7183E-01  1.0040E+00  1.1471E+00  1.0446E+00  9.7323E-01  1.0583E+00  9.6840E-01  8.8731E-01  1.0710E+00  7.9520E-01
             1.9920E+00
 PARAMETER:  7.1430E-02  1.0396E-01  2.3728E-01  1.4363E-01  7.2866E-02  1.5665E-01  6.7892E-02 -1.9557E-02  1.6860E-01 -1.2916E-01
             7.8913E-01
 GRADIENT:   6.2404E+01  6.4277E+00  6.2102E-01  9.7075E+00  1.6151E+01  5.5406E+00  2.9811E+00  3.2985E+00  1.5732E+00 -1.1188E+00
            -4.1333E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1437.24366041641        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.5800E-01  7.5815E-01  8.1970E-01  1.1805E+00  7.5329E-01  1.0102E+00  1.1810E+00  2.0702E-01  9.6059E-01  7.4243E-01
             1.9895E+00
 PARAMETER:  5.7094E-02 -1.7688E-01 -9.8822E-02  2.6591E-01 -1.8330E-01  1.1015E-01  2.6639E-01 -1.4749E+00  5.9796E-02 -1.9783E-01
             7.8789E-01
 GRADIENT:   3.2698E+01 -8.7249E+00 -3.7272E+01  3.1489E+01  5.8391E+01 -1.1938E+01 -4.5275E+00  6.4582E-01  1.6394E+00  7.1806E+00
            -2.1061E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1443.66146498301        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.4241E-01  5.4026E-01  5.6610E-01  1.2468E+00  5.1818E-01  1.0517E+00  1.7589E+00  1.0816E-01  8.3200E-01  4.0957E-01
             2.0986E+00
 PARAMETER:  4.0687E-02 -5.1571E-01 -4.6898E-01  3.2060E-01 -5.5744E-01  1.5039E-01  6.6468E-01 -2.1241E+00 -8.3927E-02 -7.9266E-01
             8.4127E-01
 GRADIENT:  -8.8895E+00  1.2914E+01  7.2228E-01  2.0757E+01 -7.2679E+00  3.3957E+00  1.2905E+00  2.3981E-01 -1.8227E+00 -4.5843E-01
             1.2743E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1444.60587212736        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.4656E-01  4.0193E-01  4.8695E-01  1.2819E+00  4.3618E-01  1.0430E+00  2.0574E+00  1.8035E-02  8.2082E-01  4.8276E-01
             1.9085E+00
 PARAMETER:  4.5075E-02 -8.1149E-01 -6.1960E-01  3.4831E-01 -7.2969E-01  1.4211E-01  8.2143E-01 -3.9154E+00 -9.7456E-02 -6.2824E-01
             7.4631E-01
 GRADIENT:   4.6430E+00  7.9263E+00  7.3594E+00  1.4172E+01 -1.1589E+01 -5.5080E-01 -1.9560E+00  9.6063E-03 -2.8404E+00 -1.5360E+00
            -1.7339E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1444.88093558724        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.4408E-01  3.3116E-01  4.5247E-01  1.2934E+00  4.0111E-01  1.0434E+00  2.3260E+00  1.0000E-02  8.2116E-01  5.0441E-01
             1.9366E+00
 PARAMETER:  4.2455E-02 -1.0052E+00 -6.9304E-01  3.5725E-01 -8.1352E-01  1.4244E-01  9.4416E-01 -5.0408E+00 -9.7036E-02 -5.8437E-01
             7.6094E-01
 GRADIENT:  -9.8911E-01  4.0239E+00  3.6316E+00  2.7426E+00 -8.6900E+00 -6.3664E-02 -1.7452E-01  0.0000E+00 -2.9898E-01  7.4858E-01
            -7.9544E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1445.20720220183        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  9.4325E-01  2.8971E-01  4.6318E-01  1.3133E+00  4.0165E-01  1.0419E+00  2.5593E+00  1.0000E-02  8.1401E-01  5.1371E-01
             1.9414E+00
 PARAMETER:  4.1571E-02 -1.1389E+00 -6.6965E-01  3.7254E-01 -8.1218E-01  1.4100E-01  1.0397E+00 -5.6216E+00 -1.0578E-01 -5.6610E-01
             7.6342E-01
 GRADIENT:  -6.4506E-02  8.2466E-02  3.4634E-01 -1.8874E-01 -7.7708E-01  1.3142E-02  1.3086E-02  0.0000E+00  3.0095E-02  1.2189E-01
             3.6465E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1445.97518200370        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  9.4616E-01  3.1030E-01  5.6095E-01  1.3483E+00  4.6448E-01  1.0429E+00  2.5495E+00  1.0000E-02  8.0471E-01  5.3382E-01
             1.9857E+00
 PARAMETER:  4.4655E-02 -1.0702E+00 -4.7813E-01  3.9885E-01 -6.6683E-01  1.4200E-01  1.0359E+00 -4.9097E+00 -1.1728E-01 -5.2769E-01
             7.8599E-01
 GRADIENT:   2.1528E-01 -8.7462E-03  2.1450E-02 -6.2066E-01  5.0453E-02  4.1668E-02  4.9220E-02  0.0000E+00  9.7570E-02  5.4393E-02
             1.6184E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1445.97534835730        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  9.4606E-01  3.0959E-01  5.6098E-01  1.3490E+00  4.6435E-01  1.0428E+00  2.5524E+00  1.0000E-02  8.0435E-01  5.3373E-01
             1.9854E+00
 PARAMETER:  4.4555E-02 -1.0725E+00 -4.7807E-01  3.9935E-01 -6.6712E-01  1.4187E-01  1.0370E+00 -4.9170E+00 -1.1771E-01 -5.2787E-01
             7.8583E-01
 GRADIENT:   4.0702E-02  1.6714E-02  4.3475E-02  6.9172E-02 -5.9558E-02 -9.6279E-03 -3.4620E-03  0.0000E+00  6.1221E-03 -9.4334E-03
             3.0507E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      663
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4152E-03  3.8598E-02 -3.3857E-04 -2.7702E-02  1.0668E-02
 SE:             2.9410E-02  1.8633E-02  2.0887E-04  2.5194E-02  1.7223E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6162E-01  3.8310E-02  1.0503E-01  2.7153E-01  5.3567E-01

 ETASHRINKSD(%)  1.4724E+00  3.7578E+01  9.9300E+01  1.5597E+01  4.2300E+01
 ETASHRINKVR(%)  2.9230E+00  6.1035E+01  9.9995E+01  2.8761E+01  6.6707E+01
 EBVSHRINKSD(%)  1.5236E+00  4.4337E+01  9.9168E+01  1.3734E+01  3.8674E+01
 EBVSHRINKVR(%)  3.0239E+00  6.9016E+01  9.9993E+01  2.5582E+01  6.2391E+01
 RELATIVEINF(%)  9.6088E+01  7.7612E+00  2.7531E-04  2.6848E+01  1.4117E+00
 EPSSHRINKSD(%)  3.5785E+01
 EPSSHRINKVR(%)  5.8764E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1445.9753483573047     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -710.82452179356653     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.24
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1445.975       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.46E-01  3.10E-01  5.61E-01  1.35E+00  4.64E-01  1.04E+00  2.55E+00  1.00E-02  8.04E-01  5.34E-01  1.99E+00
 


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
+        1.11E+03
 
 TH 2
+       -3.52E+01  6.38E+02
 
 TH 3
+        1.79E+01  5.93E+02  2.59E+03
 
 TH 4
+       -2.37E+01  2.72E+02 -3.32E+02  7.15E+02
 
 TH 5
+        2.78E+01 -1.23E+03 -4.07E+03  1.84E+02  6.82E+03
 
 TH 6
+        4.45E+00 -6.19E+00  9.90E+00 -9.58E+00 -6.50E+00  1.71E+02
 
 TH 7
+        1.18E+00  4.14E+01  1.84E+00 -3.81E+00 -9.49E+00  2.17E-01  8.12E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.49E+00 -1.48E+01 -1.96E+01 -1.70E+01  6.08E+01 -1.56E+00  6.09E+00  0.00E+00  1.83E+02
 
 TH10
+       -4.71E+00  1.99E+01 -1.15E+02 -2.67E+01  1.26E+02 -5.06E-01  3.48E+00  0.00E+00 -6.01E+00  1.09E+02
 
 TH11
+       -1.28E+01 -2.61E+00 -4.64E+01 -1.24E+01  2.77E+01  3.16E+00  1.74E+00  0.00E+00  1.04E+01  3.04E+01  6.70E+01
 
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
 #CPUT: Total CPU Time in Seconds,       12.683
Stop Time:
Sat Sep 25 08:56:44 CDT 2021
