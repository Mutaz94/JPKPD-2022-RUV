Wed Sep 29 18:34:16 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat87.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1674.43736395736        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1968E+02 -6.5267E+01 -7.8223E+01  2.7782E+01  1.2957E+02  4.9011E+01 -1.7004E+00  9.1457E+00  1.1403E+01  2.9076E+00
            -1.8477E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1685.49918347602        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  1.0344E+00  1.0947E+00  1.1176E+00  1.0252E+00  9.7221E-01  1.0577E+00  9.9989E-01  9.6877E-01  9.8876E-01  9.5691E-01
             1.0620E+00
 PARAMETER:  1.3380E-01  1.9052E-01  2.1118E-01  1.2493E-01  7.1820E-02  1.5613E-01  9.9887E-02  6.8270E-02  8.8699E-02  5.5952E-02
             1.6019E-01
 GRADIENT:   3.5567E+01  2.7849E+01 -2.3760E+00  3.7195E+01 -1.3612E+01  4.6815E+00  1.5709E+00  2.0238E+00 -1.3018E+00  1.7804E+00
             3.4373E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1685.92357674927        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0313E+00  1.0897E+00  1.2847E+00  1.0272E+00  1.0355E+00  1.0256E+00  9.3292E-01  1.0325E+00  9.8431E-01  1.0299E+00
             1.0536E+00
 PARAMETER:  1.3086E-01  1.8591E-01  3.5053E-01  1.2688E-01  1.3489E-01  1.2530E-01  3.0561E-02  1.3197E-01  8.4184E-02  1.2945E-01
             1.5225E-01
 GRADIENT:   3.3072E+01  2.2884E+01  3.4607E+00  2.7422E+01 -4.3099E+00 -7.3097E+00 -2.2234E+00 -2.7986E+00 -6.6482E+00  2.5652E+00
            -9.1241E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1686.33992690042        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      545             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0287E+00  1.0821E+00  1.3555E+00  1.0314E+00  1.0573E+00  1.0401E+00  9.4725E-01  1.1606E+00  9.9110E-01  1.0314E+00
             1.0533E+00
 PARAMETER:  1.2828E-01  1.7893E-01  4.0420E-01  1.3090E-01  1.5572E-01  1.3931E-01  4.5805E-02  2.4894E-01  9.1058E-02  1.3093E-01
             1.5193E-01
 GRADIENT:   5.6123E+02  8.9058E+01  4.2996E+00  1.1488E+02  1.1721E+01  7.9432E+01  2.6180E+00 -1.3701E+00  3.4991E+00  1.9710E+00
             1.0064E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1686.36773760811        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  1.0287E+00  1.0821E+00  1.3504E+00  1.0009E+00  1.0573E+00  1.0453E+00  9.5696E-01  1.1592E+00  9.9060E-01  1.0321E+00
             1.0524E+00
 PARAMETER:  1.2828E-01  1.7893E-01  4.0041E-01  1.0094E-01  1.5572E-01  1.4432E-01  5.6004E-02  2.4770E-01  9.0558E-02  1.3158E-01
             1.5104E-01
 GRADIENT:   2.7820E+01 -1.1996E+01  3.2220E+00 -2.2025E+01 -8.2582E-02  8.2506E-01  3.7042E-01 -1.7666E+00 -4.5500E+00  1.4738E+00
            -7.3960E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1686.52638432341        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      905
 NPARAMETR:  1.0287E+00  1.0821E+00  1.3504E+00  1.0156E+00  1.0573E+00  1.0434E+00  9.5169E-01  1.1592E+00  9.9060E-01  1.0321E+00
             1.0524E+00
 PARAMETER:  1.2828E-01  1.7893E-01  4.0041E-01  1.1548E-01  1.5572E-01  1.4246E-01  5.0481E-02  2.4770E-01  9.0559E-02  1.3158E-01
             1.5104E-01
 GRADIENT:   2.7636E+01  1.8204E+00  1.3505E+00  5.9957E-02  1.8349E+00  4.7436E-03 -6.6985E-02 -1.6416E+00 -3.9785E+00  1.2957E+00
            -2.6293E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1686.72880458841        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  1.0169E+00  1.0784E+00  1.3465E+00  1.0143E+00  1.0558E+00  1.0384E+00  9.5082E-01  1.1829E+00  9.9953E-01  1.0275E+00
             1.0524E+00
 PARAMETER:  1.1678E-01  1.7550E-01  3.9749E-01  1.1421E-01  1.5432E-01  1.3773E-01  4.9565E-02  2.6797E-01  9.9529E-02  1.2709E-01
             1.5112E-01
 GRADIENT:   4.9099E+02  6.8788E+01  3.1805E+00  7.5511E+01  1.2744E+01  8.0911E+01  3.5480E+00 -1.5917E-01  4.9861E+00  2.1463E+00
             1.5489E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1686.76153977457        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1208             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0160E+00  1.0794E+00  1.3470E+00  1.0151E+00  1.0548E+00  1.0436E+00  9.3958E-01  1.1942E+00  1.0052E+00  1.0234E+00
             1.0524E+00
 PARAMETER:  1.1591E-01  1.7645E-01  3.9790E-01  1.1503E-01  1.5338E-01  1.4263E-01  3.7677E-02  2.7747E-01  1.0514E-01  1.2310E-01
             1.5107E-01
 GRADIENT:   4.8557E+02  7.0889E+01  3.0597E+00  7.8296E+01  1.1762E+01  8.6079E+01  3.1407E+00  1.5523E-01  5.9480E+00  1.6468E+00
             1.5007E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1686.76685933619        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1397
 NPARAMETR:  1.0161E+00  1.0799E+00  1.3493E+00  1.0151E+00  1.0542E+00  1.0439E+00  9.3792E-01  1.2066E+00  1.0057E+00  1.0227E+00
             1.0517E+00
 PARAMETER:  1.1599E-01  1.7689E-01  3.9959E-01  1.1503E-01  1.5275E-01  1.4298E-01  3.5907E-02  2.8784E-01  1.0567E-01  1.2248E-01
             1.5036E-01
 GRADIENT:   2.3161E+00 -9.1316E-01 -7.2064E-01 -1.2497E+00  7.0689E-01  5.5584E-01 -3.3933E-02  3.2600E-02 -1.4716E+00  8.0653E-01
             4.2706E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1686.76685933619        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     1421
 NPARAMETR:  1.0161E+00  1.0799E+00  1.3493E+00  1.0151E+00  1.0542E+00  1.0439E+00  9.3792E-01  1.2066E+00  1.0057E+00  1.0227E+00
             1.0517E+00
 PARAMETER:  1.1599E-01  1.7689E-01  3.9959E-01  1.1503E-01  1.5275E-01  1.4298E-01  3.5907E-02  2.8784E-01  1.0567E-01  1.2248E-01
             1.5036E-01
 GRADIENT:   1.8897E+05 -6.1958E+04  5.4877E+04 -1.9867E+01  1.4350E+05 -3.2983E-03 -2.4157E-02  7.6059E+04 -1.0372E+05  8.9474E+04
            -1.4611E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1421
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1641E-04 -1.9939E-02 -2.8378E-02  7.3049E-03 -3.5734E-02
 SE:             2.9832E-02  1.6593E-02  1.3028E-02  2.5213E-02  2.2309E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9689E-01  2.2950E-01  2.9385E-02  7.7202E-01  1.0921E-01

 ETASHRINKSD(%)  5.8907E-02  4.4410E+01  5.6355E+01  1.5534E+01  2.5261E+01
 ETASHRINKVR(%)  1.1778E-01  6.9098E+01  8.0951E+01  2.8655E+01  4.4141E+01
 EBVSHRINKSD(%)  4.3937E-01  4.3877E+01  6.0019E+01  1.6992E+01  2.2337E+01
 EBVSHRINKVR(%)  8.7681E-01  6.8502E+01  8.4015E+01  3.1097E+01  3.9684E+01
 RELATIVEINF(%)  9.8395E+01  7.0525E-01  2.0521E+00  1.7438E+00  1.1180E+01
 EPSSHRINKSD(%)  4.4098E+01
 EPSSHRINKVR(%)  6.8750E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1686.7668593361859     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -951.61603277244774     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.76
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1686.767       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.08E+00  1.35E+00  1.02E+00  1.05E+00  1.04E+00  9.38E-01  1.21E+00  1.01E+00  1.02E+00  1.05E+00
 


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
+        3.95E+07
 
 TH 2
+       -3.66E+02  3.00E+07
 
 TH 3
+       -1.22E+02 -4.94E+03  3.77E+06
 
 TH 4
+        3.98E+07  4.25E+02  8.70E+06  8.04E+07
 
 TH 5
+       -3.45E+03 -1.78E+07  6.31E+06 -2.91E+07  2.11E+07
 
 TH 6
+        3.12E+07 -5.80E+02  2.06E+02 -9.46E+02 -2.28E+07  1.80E+02
 
 TH 7
+        2.02E+02 -1.22E+02  4.76E+01 -2.08E+02  3.63E+07  1.28E-01  3.06E+01
 
 TH 8
+       -1.34E+07  1.90E+04  6.69E+03 -1.35E+07  1.17E+03 -1.06E+07  1.68E+07  4.53E+06
 
 TH 9
+       -6.37E+02 -1.75E+01 -8.98E+03  4.42E+07 -1.10E+03 -1.04E+03 -1.94E+02  1.49E+07  4.85E+07
 
 TH10
+        3.71E+07 -2.29E+07  7.60E+03 -3.75E+07 -6.85E+02  8.81E+02  1.93E+02 -2.78E+02 -1.41E+03  3.49E+07
 
 TH11
+       -2.95E+07  8.45E+04 -3.00E+04 -2.96E+07 -2.16E+07  2.33E+07 -1.46E+02 -1.00E+07  1.52E+05 -1.29E+05  2.20E+07
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       26.182
Stop Time:
Wed Sep 29 18:34:44 CDT 2021
